#!/usr/bin/env python3
"""
Build the Space Coast flood hot-spot dataset (Stage IV edition, 2002-present).

For every flood report the script finds the PEAK Stage IV rainfall accumulation
for each duration (1, 2, 3, 6, 12, 24, 48, 72 h) in the hours leading up to it.
The map can then answer "what has flooded when >= X in fell in N hours".

Sources (public, no keys):
  1. IEM Local Storm Reports, WFO MLB: FLASH FLOOD, FLOOD (+ HEAVY RAIN obs)
       https://mesonet.agron.iastate.edu/geojson/lsr.php
  2. NCEI Storm Events, Brevard Co. FL: Flash Flood, Flood
       https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/
  3. NCEP Stage IV QPE (IEM archive), ~4.7 km HRAP grid
       hourly : .../archive/data/YYYY/MM/DD/stage4/ST4.YYYYMMDDHH.01h.grib
       6-hourly: .../archive/data/YYYY/MM/DD/stage4/ST4.YYYYMMDDHH.06h.grib
     (validated by live directory probe 2026-11-24; HH = END of accumulation, UTC)

Peak definition
  Durations 1-24 h: rolling sums of 01h files.
  Durations 48/72 h: rolling sums of 06h files (RFC-QC'd 6-h product).
  Windows considered: those ENDING between (report_begin - 12 h) and (report_end + 2 h)
  (+6 h for the 6-hourly durations to reach the next synoptic time).
  Any missing hour inside a window invalidates that window.
  Two values per duration:
    pt = Stage IV grid cell nearest the report
    nb = max over grid cells within NB_KM of the report (report locations are
         often rounded to 0.01 deg / "2 N MIMS"-style offsets)

Cache: data/stage4/{01h,06h}_{YYYY}.npz (bbox subset, uint16 tenths of mm) and
data/stage4/grid.npz. Commit them: reruns only download hours not yet cached.
"""
import concurrent.futures as cf
import csv
import glob
import gzip
import io
import json
import math
import os
import re
import sys
import time
from datetime import date, datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import numpy as np
import pygrib
import requests

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
BBOX = dict(south=27.95, north=28.85, west=-81.00, east=-80.45)
WFO = "MLB"
COUNTY = "BREVARD"
START_YEAR = int(os.environ.get("START_YEAR", 2002))
END_YEAR = datetime.now(timezone.utc).year

DURS_1H = [1, 2, 3, 6, 12, 24]
DURS_6H = [48, 72]
DURS = DURS_1H + DURS_6H
LEAD_H = 12          # windows may end this many hours before the report begin time
LAG_H = 2            # ...and this many hours after the report end time
MAX_EVENT_H = 72     # cap on Storm Events end-begin span
NB_KM = 6.0          # neighborhood radius for "nb" values
MAX_FETCH = int(os.environ.get("MAX_FETCH", 25000))   # per run; reruns resume
WORKERS = int(os.environ.get("WORKERS", 6))

LSR_FLOOD_TYPES = {"FLASH FLOOD", "FLOOD"}
LSR_RAIN_TYPE = "HEAVY RAIN"
SE_FLOOD_TYPES = {"Flash Flood", "Flood"}

IEM = "https://mesonet.agron.iastate.edu"
SE_DIR = "https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "flood_hotspots.json")
S4DIR = os.path.join(ROOT, "data", "stage4")
TZ = ZoneInfo("America/New_York")
MISSING = np.uint16(65535)

S = requests.Session()
S.headers["User-Agent"] = "space-coast-flood-hotspots/2.0 (GitHub Actions)"
S.mount("https://", requests.adapters.HTTPAdapter(pool_connections=WORKERS, pool_maxsize=WORKERS))


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, params=None, timeout=180, allow_404=False):
    last = None
    for attempt in range(5):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404 and allow_404:
                return None
            last = f"HTTP {r.status_code}"
            if r.status_code not in (429, 500, 502, 503, 504):
                break
        except requests.RequestException as e:
            last = repr(e)
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"GET failed {url} {params or ''}: {last}")


def in_bbox(lat, lon):
    return BBOX["south"] <= lat <= BBOX["north"] and BBOX["west"] <= lon <= BBOX["east"]


def km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def hkey(dt):
    return int(dt.strftime("%Y%m%d%H"))


def hstr(dt):
    return dt.strftime("%Y-%m-%dT%HZ")


def floor_h(dt):
    return dt.replace(minute=0, second=0, microsecond=0)


def ceil_h(dt):
    f = floor_h(dt)
    return f if f == dt else f + timedelta(hours=1)


# --------------------------------------------------------------------------
# 1. Local Storm Reports
# --------------------------------------------------------------------------
def fetch_lsr():
    floods, rains, seen = [], [], set()
    for y in range(START_YEAR, END_YEAR + 1):
        params = dict(sts=f"{y}-01-01T00:00Z", ets=f"{y + 1}-01-01T00:00Z", wfos=WFO)
        try:
            fc = get(f"{IEM}/geojson/lsr.php", params).json()
        except Exception as e:
            log(f"  LSR {y}: FAILED {e}")
            continue
        nf = nr = 0
        for f in fc.get("features", []):
            p = f.get("properties", {})
            tt = (p.get("typetext") or "").upper()
            if tt not in LSR_FLOOD_TYPES and tt != LSR_RAIN_TYPE:
                continue
            lat, lon = fnum(p.get("lat")), fnum(p.get("lon"))
            if lat is None or lon is None or not in_bbox(lat, lon):
                continue
            valid = p.get("valid")
            key = (valid, round(lat, 3), round(lon, 3), tt)
            if key in seen:
                continue
            seen.add(key)
            t = datetime.fromisoformat(valid.replace("Z", "+00:00"))
            rec = dict(
                src="LSR", type=tt.title(),
                t=t.strftime("%Y-%m-%dT%H:%MZ"), t_end=t.strftime("%Y-%m-%dT%H:%MZ"),
                d=t.astimezone(TZ).date().isoformat(),
                lat=round(lat, 4), lon=round(lon, 4),
                city=p.get("city") or "", remark=(p.get("remark") or "").strip(),
                rsrc=p.get("source") or "",
            )
            if tt == LSR_RAIN_TYPE:
                mag = fnum(p.get("magf"))
                if mag is None or mag <= 0:
                    continue
                rec["inches"] = round(mag, 2)
                rains.append(rec)
                nr += 1
            else:
                floods.append(rec)
                nf += 1
        log(f"  LSR {y}: {nf} flood, {nr} heavy-rain in bbox")
        time.sleep(1)
    return floods, rains


# --------------------------------------------------------------------------
# 2. NCEI Storm Events
# --------------------------------------------------------------------------
def storm_events_files():
    html = get(SE_DIR).text
    best = {}
    for m in re.finditer(r"StormEvents_details-ftp_v1\.0_d(\d{4})_c(\d{8})\.csv\.gz", html):
        y, c = int(m.group(1)), m.group(2)
        if y >= START_YEAR and (y not in best or c > best[y][0]):
            best[y] = (c, m.group(0))
    return {y: SE_DIR + v[1] for y, v in sorted(best.items())}


def se_dt(ym, day, hhmm, tzs):
    """Storm Events local standard time -> UTC. tzs like 'EST-5'."""
    ym = str(ym).strip()
    hhmm = str(hhmm or "0").strip().zfill(4)
    m = re.search(r"([+-]\d+)\s*$", tzs or "")
    off = int(m.group(1)) if m else -5
    local = datetime(int(ym[:4]), int(ym[4:6]), int(str(day).strip()),
                     int(hhmm[:2]) % 24, int(hhmm[2:]) % 60)
    return (local - timedelta(hours=off)).replace(tzinfo=timezone.utc)


def fetch_storm_events():
    out = []
    for y, url in storm_events_files().items():
        try:
            raw = get(url, timeout=300).content
        except Exception as e:
            log(f"  SE {y}: FAILED {e}")
            continue
        n = 0
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
            for row in csv.DictReader(io.TextIOWrapper(gz, encoding="latin-1", newline="")):
                if row.get("STATE") != "FLORIDA" or row.get("EVENT_TYPE") not in SE_FLOOD_TYPES:
                    continue
                if COUNTY not in (row.get("CZ_NAME") or "").upper():
                    continue
                lat, lon = fnum(row.get("BEGIN_LAT")), fnum(row.get("BEGIN_LON"))
                if lat is None or lon is None or not in_bbox(lat, lon):
                    continue
                try:
                    tz = row.get("CZ_TIMEZONE", "")
                    tb = se_dt(row["BEGIN_YEARMONTH"], row["BEGIN_DAY"], row.get("BEGIN_TIME"), tz)
                    te = se_dt(row["END_YEARMONTH"], row["END_DAY"], row.get("END_TIME"), tz)
                except (KeyError, ValueError):
                    continue
                te = min(max(te, tb), tb + timedelta(hours=MAX_EVENT_H))
                out.append(dict(
                    src="StormEvents", type=row["EVENT_TYPE"],
                    t=tb.strftime("%Y-%m-%dT%H:%MZ"), t_end=te.strftime("%Y-%m-%dT%H:%MZ"),
                    d=tb.astimezone(TZ).date().isoformat(),
                    lat=round(lat, 4), lon=round(lon, 4),
                    city=(row.get("BEGIN_LOCATION") or "").strip(),
                    remark=(row.get("EVENT_NARRATIVE") or "").strip(),
                    episode=(row.get("EPISODE_NARRATIVE") or "").strip()[:600],
                    cause=(row.get("FLOOD_CAUSE") or "").strip(),
                    event_id=row.get("EVENT_ID"),
                ))
                n += 1
        log(f"  SE {y}: {n} flood events in bbox")
    return out


def merge(lsr, se):
    """Drop Storm Events duplicating an LSR (<=1.5 km, <=1 day); keep the narrative."""
    kept = []
    for e in se:
        ed = date.fromisoformat(e["d"])
        dup = next((r for r in lsr
                    if abs((date.fromisoformat(r["d"]) - ed).days) <= 1
                    and float(km(r["lat"], r["lon"], e["lat"], e["lon"])) <= 1.5), None)
        if dup is not None:
            if e["remark"] and "se_narrative" not in dup:
                dup["se_narrative"] = e["remark"]
        else:
            kept.append(e)
    return lsr + kept


# --------------------------------------------------------------------------
# 3. Stage IV cache
# --------------------------------------------------------------------------
class Stage4:
    def __init__(self):
        os.makedirs(S4DIR, exist_ok=True)
        self.lat = self.lon = None
        self.idx = {}                        # grid shape -> (iy, ix) into native grid
        self.data = {"01h": {}, "06h": {}}   # kind -> {YYYYMMDDHH: uint16[ny,nx]}
        self.dirty = {"01h": set(), "06h": set()}
        g = os.path.join(S4DIR, "grid.npz")
        if os.path.exists(g):
            z = np.load(g)
            self.lat, self.lon = z["lat"], z["lon"]
        for f in glob.glob(os.path.join(S4DIR, "*_*.npz")):
            kind = os.path.basename(f).split("_")[0]
            if kind not in self.data:
                continue
            z = np.load(f)
            for h, v in zip(z["h"], z["v"]):
                self.data[kind][int(h)] = v
        log(f"  Stage IV cache: {len(self.data['01h'])} hourly, {len(self.data['06h'])} 6-hourly grids")

    # ---- grid handling
    def _index_for(self, msg):
        lats, lons = msg.latlons()
        lons = ((lons + 180.0) % 360.0) - 180.0
        if self.lat is None:
            pad = 0.15
            m = ((lats >= BBOX["south"] - pad) & (lats <= BBOX["north"] + pad)
                 & (lons >= BBOX["west"] - pad) & (lons <= BBOX["east"] + pad))
            rows, cols = np.where(m.any(1))[0], np.where(m.any(0))[0]
            iy, ix = np.meshgrid(np.arange(rows.min(), rows.max() + 1),
                                 np.arange(cols.min(), cols.max() + 1), indexing="ij")
            self.lat, self.lon = lats[iy, ix].astype("f4"), lons[iy, ix].astype("f4")
            np.savez_compressed(os.path.join(S4DIR, "grid.npz"), lat=self.lat, lon=self.lon)
            log(f"  Stage IV reference subset {self.lat.shape} from native {lats.shape}")
            return iy, ix
        # a different native grid: map every reference cell to its nearest native cell
        pad = 0.4
        cand = np.where((lats >= BBOX["south"] - pad) & (lats <= BBOX["north"] + pad)
                        & (lons >= BBOX["west"] - pad) & (lons <= BBOX["east"] + pad))
        clat, clon = lats[cand], lons[cand]
        iy = np.empty(self.lat.shape, int)
        ix = np.empty(self.lat.shape, int)
        for j in np.ndindex(self.lat.shape):
            k = int(np.argmin(km(self.lat[j], self.lon[j], clat, clon)))
            iy[j], ix[j] = cand[0][k], cand[1][k]
        log(f"  Stage IV remapped native grid {lats.shape} onto reference subset")
        return iy, ix

    def decode(self, content):
        msg = pygrib.fromstring(content)
        vals = np.ma.filled(np.ma.asarray(msg.values, dtype="f4"), np.nan)
        if vals.shape not in self.idx:
            self.idx[vals.shape] = self._index_for(msg)
        iy, ix = self.idx[vals.shape]
        sub = vals[iy, ix]                                   # mm
        ok = np.isfinite(sub) & (sub >= 0) & (sub < 6000)
        return np.where(ok, np.round(sub * 10.0), MISSING).astype(np.uint16)

    # ---- fetching
    @staticmethod
    def url(kind, dt):
        return (f"{IEM}/archive/data/{dt:%Y/%m/%d}/stage4/"
                f"ST4.{dt:%Y%m%d%H}.{kind}.grib")

    def _fetch_one(self, kind, dt):
        try:
            r = get(self.url(kind, dt), timeout=90, allow_404=True)
        except Exception as e:
            return kind, dt, None, f"error {e}"
        if r is None:
            return kind, dt, None, "404"
        return kind, dt, r.content, None

    def ensure(self, needed):
        """needed: set of (kind, datetime). Downloads uncached hours."""
        todo = sorted(k for k in needed if hkey(k[1]) not in self.data[k[0]])
        recent = datetime.now(timezone.utc) - timedelta(days=5)
        log(f"  Stage IV: {len(needed)} grids needed, {len(todo)} not cached")
        if len(todo) > MAX_FETCH:
            log(f"  capping this run at MAX_FETCH={MAX_FETCH}; rerun to continue")
            todo = todo[:MAX_FETCH]
        if not todo:
            return
        # decode one serially first so the reference grid is set before threads run
        first, rest = todo[0], todo[1:]
        self._store(*self._fetch_one(*first), recent)
        n = miss = 0
        with cf.ThreadPoolExecutor(WORKERS) as ex:
            for kind, dt, content, err in ex.map(lambda k: self._fetch_one(*k), rest):
                if self._store(kind, dt, content, err, recent) is False:
                    miss += 1
                n += 1
                if n % 500 == 0:
                    log(f"    {n}/{len(rest)} fetched ({miss} missing)")
                    self.save()
        self.save()
        log(f"  Stage IV fetched {n + 1}, missing {miss}")

    def _store(self, kind, dt, content, err, recent):
        if content is None:
            if err == "404" and dt < recent:   # permanently missing: remember it
                self.data[kind][hkey(dt)] = None
                self.dirty[kind].add(dt.year)
            return False
        try:
            self.data[kind][hkey(dt)] = self.decode(content)
            self.dirty[kind].add(dt.year)
            return True
        except Exception as e:
            log(f"    decode fail {kind} {dt}: {e}")
            return False

    def save(self):
        if self.lat is None:
            return
        shape = self.lat.shape
        for kind, years in self.dirty.items():
            for y in sorted(years):
                items = sorted((h, v) for h, v in self.data[kind].items() if h // 1000000 == y)
                if not items:
                    continue
                h = np.array([i[0] for i in items], dtype=np.int64)
                v = np.stack([i[1] if i[1] is not None else np.full(shape, MISSING, np.uint16)
                              for i in items])
                np.savez_compressed(os.path.join(S4DIR, f"{kind}_{y}.npz"), h=h, v=v)
            years.clear()

    # ---- extraction
    def series(self, kind, times, cells):
        """(T, P) inches; NaN where missing."""
        out = np.full((len(times), len(cells[0])), np.nan, "f4")
        for i, dt in enumerate(times):
            a = self.data[kind].get(hkey(dt))
            if a is None:
                continue
            v = a[cells].astype("f4")
            out[i] = np.where(v == MISSING, np.nan, v / 10.0 / 25.4)
        return out


def windows_for(r):
    tb = datetime.fromisoformat(r["t"].replace("Z", "+00:00"))
    te = datetime.fromisoformat(r["t_end"].replace("Z", "+00:00"))
    end_lo = floor_h(tb - timedelta(hours=LEAD_H))
    end_hi1 = ceil_h(te + timedelta(hours=LAG_H))
    end_hi6 = end_hi1 + timedelta(hours=(-end_hi1.hour) % 6)
    end_lo6 = end_lo - timedelta(hours=end_lo.hour % 6)
    h1 = [end_lo - timedelta(hours=max(DURS_1H) - 1) + timedelta(hours=i)
          for i in range(int((end_hi1 - end_lo).total_seconds() // 3600) + max(DURS_1H))]
    n6 = int((end_hi6 - end_lo6).total_seconds() // 21600)
    h6 = [end_lo6 - timedelta(hours=max(DURS_6H) - 6) + timedelta(hours=6 * i)
          for i in range(n6 + max(DURS_6H) // 6)]
    return end_lo, end_hi1, end_lo6, end_hi6, h1, h6


def peak(series, times, step, dur, lo, hi):
    """Max dur-hour total over windows ending in [lo, hi]; per-cell then max over cells."""
    k = dur // step
    best, best_t = None, None
    for i in range(k - 1, len(times)):
        if not (lo <= times[i] <= hi):
            continue
        w = series[i - k + 1:i + 1]
        bad = np.isnan(w).any(axis=0)
        if bad.all():
            continue
        s = np.where(bad, -1.0, np.nansum(w, axis=0)).max()
        if s >= 0 and (best is None or s > best):
            best, best_t = float(s), times[i]
    return (round(best, 2) if best is not None else None), (hstr(best_t) if best_t else None)


def attach_stage4(reports):
    s4 = Stage4()
    plans = []
    needed = set()
    for r in reports:
        lo, hi1, lo6, hi6, h1, h6 = windows_for(r)
        plans.append((lo, hi1, lo6, hi6, h1, h6))
        needed.update(("01h", t) for t in h1)
        needed.update(("06h", t) for t in h6)
    s4.ensure(needed)
    if s4.lat is None:
        log("  no Stage IV grid available; skipping rainfall")
        return
    flat_lat, flat_lon = s4.lat.ravel(), s4.lon.ravel()
    for r, (lo, hi1, lo6, hi6, h1, h6) in zip(reports, plans):
        d = km(r["lat"], r["lon"], flat_lat, flat_lon)
        c0 = int(np.argmin(d))
        nb = np.where(d <= NB_KM)[0]
        if c0 not in nb:
            nb = np.append(nb, c0)
        cells_pt = np.unravel_index(np.array([c0]), s4.lat.shape)
        cells_nb = np.unravel_index(nb, s4.lat.shape)
        out = {"pt": [], "nb": [], "te": [], "via": []}
        s1p, s1n = s4.series("01h", h1, cells_pt), s4.series("01h", h1, cells_nb)
        s6p, s6n = s4.series("06h", h6, cells_pt), s4.series("06h", h6, cells_nb)
        for dur in DURS:
            via = "1h"
            if dur in DURS_1H:
                vp, tp = peak(s1p, h1, 1, dur, lo, hi1)
                vn, _ = peak(s1n, h1, 1, dur, lo, hi1)
                if vp is None and dur % 6 == 0:   # hourly gap: fall back to 6-h files
                    vp, tp = peak(s6p, h6, 6, dur, lo6, hi6)
                    vn, _ = peak(s6n, h6, 6, dur, lo6, hi6)
                    via = "6h"
            else:
                vp, tp = peak(s6p, h6, 6, dur, lo6, hi6)
                vn, _ = peak(s6n, h6, 6, dur, lo6, hi6)
                via = "6h"
            out["pt"].append(vp)
            out["nb"].append(vn)
            out["te"].append(tp)
            out["via"].append(via if vp is not None else None)
        r["s4"] = out
        r["cell_km"] = round(float(d[c0]), 1)


# --------------------------------------------------------------------------
def main():
    log(f"Years {START_YEAR}-{END_YEAR}")
    log("Fetching Local Storm Reports...")
    lsr_flood, lsr_rain = fetch_lsr()
    log("Fetching NCEI Storm Events...")
    try:
        se = fetch_storm_events()
    except Exception as e:
        log(f"  Storm Events skipped: {e}")
        se = []
    reports = merge(lsr_flood, se)
    reports.sort(key=lambda r: r["t"])
    log(f"{len(reports)} flood reports after merge ({len(lsr_flood)} LSR, {len(se)} SE raw)")
    log("Stage IV peak accumulations...")
    attach_stage4(reports)
    for i, r in enumerate(reports):
        r["id"] = i
    lsr_rain.sort(key=lambda r: r["t"])
    out = dict(
        generated=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
        bbox=BBOX, years=[START_YEAR, END_YEAR],
        durations_h=DURS, nb_km=NB_KM, lead_h=LEAD_H, lag_h=LAG_H,
        reports=reports, rain_obs=lsr_rain,
    )
    with open(OUT, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    n_ok = sum(1 for r in reports if r.get("s4") and any(v is not None for v in r["s4"]["pt"]))
    log(f"Wrote {OUT}: {len(reports)} flood reports ({n_ok} with Stage IV), "
        f"{len(lsr_rain)} heavy-rain obs")


if __name__ == "__main__":
    main()

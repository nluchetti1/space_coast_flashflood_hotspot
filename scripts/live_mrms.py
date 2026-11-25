#!/usr/bin/env python3
"""
Real-time MRMS QPE subset for the Space Coast flood hot-spot dashboard.

Writes live/mrms_latest.json: the latest MRMS accumulation grids (0.01 deg)
over the dashboard bbox for 1, 3, 6, 12, 24, 48, 72 h, from two products:
  ro = RadarOnly_QPE_xxH            (low latency, radar only)
  ms = MultiSensor_QPE_xxH_Pass2    (gauge-corrected, posted ~57 min after valid time)

Sources (validated by live probe 2026-11-25):
  primary : https://mrms.ncep.noaa.gov/data/2D/<PRODUCT>/MRMS_<PRODUCT>.latest.grib2.gz
  fallback: https://noaa-mrms-pds.s3.amazonaws.com/CONUS/<PRODUCT>_00.00/YYYYMMDD/
            MRMS_<PRODUCT>_00.00_YYYYMMDD-HHMMSS.grib2.gz
Grid: regular lat/lon, 7000 x 3500, 0.01 deg, first point 54.995N 230.005E, rows run
north->south (jScansPositively=0), PNG packing (template 5.41). Values in mm; <0 = missing.
"""
import gzip
import json
import os
import re
import sys
import time
from datetime import datetime, timedelta, timezone

import numpy as np
import pygrib
import requests

BBOX = dict(south=27.95, north=28.85, west=-81.00, east=-80.45)
PAD = 0.07                      # room for the 6 km neighborhood at the bbox edges
DURS = [1, 3, 6, 12, 24, 48, 72]
PRODUCTS = {
    "ro": lambda h: f"RadarOnly_QPE_{h:02d}H",
    "ms": lambda h: f"MultiSensor_QPE_{h:02d}H_Pass2",
}
NCEP = "https://mrms.ncep.noaa.gov/data/2D"
AWS = "https://noaa-mrms-pds.s3.amazonaws.com"

OUT_DIR = os.environ.get("LIVE_OUT", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "live"))
S = requests.Session()
S.headers["User-Agent"] = "space-coast-flood-hotspots-live/1.0 (GitHub Actions)"


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get(url, timeout=90):
    for attempt in range(3):
        try:
            r = S.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except requests.RequestException as e:
            log(f"  {url}: {e}")
        time.sleep(3 * (attempt + 1))
    return None


def fetch_latest(product):
    r = get(f"{NCEP}/{product}/MRMS_{product}.latest.grib2.gz")
    if r is not None:
        return r.content, "ncep"
    # fallback: newest key in today's (or yesterday's) AWS folder
    now = datetime.now(timezone.utc)
    for d in (now, now - timedelta(days=1)):
        prefix = f"CONUS/{product}_00.00/{d:%Y%m%d}/"
        lst = get(f"{AWS}/?list-type=2&prefix={prefix}")
        if lst is None:
            continue
        keys = re.findall(r"<Key>([^<]+)</Key>", lst.text)
        while "<IsTruncated>true</IsTruncated>" in lst.text:
            tok = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", lst.text)
            if not tok:
                break
            lst = get(f"{AWS}/?list-type=2&prefix={prefix}&continuation-token={requests.utils.quote(tok.group(1))}")
            if lst is None:
                break
            keys += re.findall(r"<Key>([^<]+)</Key>", lst.text)
        if keys:
            r = get(f"{AWS}/{sorted(keys)[-1]}")
            if r is not None:
                return r.content, "aws"
    return None, None


def subset(raw_gz):
    msg = pygrib.fromstring(gzip.decompress(raw_gz))
    lat0 = msg["latitudeOfFirstGridPointInDegrees"]
    lon0 = msg["longitudeOfFirstGridPointInDegrees"]
    lon0 = lon0 - 360 if lon0 > 180 else lon0
    dlon = msg["iDirectionIncrementInDegrees"]
    dlat = msg["jDirectionIncrementInDegrees"]
    if msg["jScansPositively"] != 0:
        raise ValueError("unexpected scan order")
    n, s = BBOX["north"] + PAD, BBOX["south"] - PAD
    w, e = BBOX["west"] - PAD, BBOX["east"] + PAD
    r0, r1 = int(round((lat0 - n) / dlat)), int(round((lat0 - s) / dlat))
    c0, c1 = int(round((w - lon0) / dlon)), int(round((e - lon0) / dlon))
    vals = np.ma.filled(np.ma.asarray(msg.values, dtype="f4"), -1.0)[r0:r1 + 1, c0:c1 + 1]
    inch100 = np.where(vals >= 0, np.round(vals / 25.4 * 100.0), -1).astype(np.int32)
    valid = f"{msg['dataDate']:08d}{msg['dataTime']:04d}"
    valid = datetime.strptime(valid, "%Y%m%d%H%M").strftime("%Y-%m-%dT%H:%MZ")
    grid = dict(lat_n=round(lat0 - r0 * dlat, 4), lon_w=round(lon0 + c0 * dlon, 4),
                dlat=dlat, dlon=dlon, ny=inch100.shape[0], nx=inch100.shape[1])
    return grid, valid, inch100


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    out = dict(generated=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ"),
               durations_h=DURS, units="hundredths of inch; -1 missing",
               grid=None, valid={}, src={}, ro={}, ms={})
    for key, name in PRODUCTS.items():
        out["valid"][key] = {}
        for h in DURS:
            prod = name(h)
            raw, where = fetch_latest(prod)
            if raw is None:
                log(f"  {prod}: unavailable")
                continue
            try:
                grid, valid, arr = subset(raw)
            except Exception as e:
                log(f"  {prod}: decode failed {e}")
                continue
            if out["grid"] is None:
                out["grid"] = grid
            elif (grid["ny"], grid["nx"]) != (out["grid"]["ny"], out["grid"]["nx"]):
                log(f"  {prod}: grid mismatch, skipped")
                continue
            out[key][str(h)] = arr.ravel().tolist()
            out["valid"][key][str(h)] = valid
            out["src"][prod] = where
            log(f"  {prod}: valid {valid} via {where}, max {arr.max() / 100:.2f} in")
    if out["grid"] is None:
        log("No MRMS products decoded; not writing.")
        sys.exit(1)
    path = os.path.join(OUT_DIR, "mrms_latest.json")
    with open(path, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    log(f"Wrote {path} ({os.path.getsize(path) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()

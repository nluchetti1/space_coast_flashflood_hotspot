# Space Coast Flood Hot Spots

An interactive map of where flooding has been reported around Cape Canaveral SFS, KSC, Patrick SFB and Merritt Island, and how much rain, over how long, it took to cause it.

Choose a duration (1, 2, 3, 6, 12, 24, 48 or 72 h) and a rainfall amount, for example 4" in 3 h. The map then shows every spot that has flooded with that much rain or less in that window. Click any spot to see its trigger for every duration, or click anywhere on the map to check the hot spots within 2 km.

## How it works

| Step | Source | What |
|---|---|---|
| 1 | IEM Local Storm Reports (WFO MLB, 2002–present) | `FLASH FLOOD` / `FLOOD` reports with lat/lon and remarks; `HEAVY RAIN` reports as observed totals |
| 2 | NCEI Storm Events (Brevard, 2002–present) | `Flash Flood` / `Flood` events with begin/end times and narratives. Events within 1.5 km and 1 day of an LSR are merged into it |
| 3 | NCEP Stage IV (IEM archive) | Peak accumulation for each duration before each report |

**Peak accumulation.** For each report and duration, the build takes the largest rolling total among windows that end between 12 h before the report's begin time and 2 h after its end time.
- **1–24 h:** sums of hourly (`01h`) files.
- **48/72 h:** sums of 6-hourly (`06h`) files, which get more RFC QC.
- **Gaps:** if an hourly gap kills a 6, 12 or 24 h window, that value falls back to 6-hourly files and is marked †.
- **Two values per duration:**
  - the nearest grid cell
  - the largest value within 6 km, because report locations are often coarse ("2 N MIMS")

**Hot spots.** Reports are clustered on the page (radius 250 m to 2 km). A spot's **trigger** for a duration is the lowest peak that has ever produced a flood report there.

Bounding box (edit `BBOX` in `scripts/build_data.py`): 27.95–28.85 N, 81.00–80.45 W.

## Setup

1. Push this folder to a new GitHub repo.
2. **Settings → Pages →** Deploy from branch `main`, folder `/ (root)`.
3. **Actions → Build flood data → Run workflow.**
   - The first run downloads every Stage IV hour around each report, roughly 50–100 files per report minus overlap. That's likely several thousand to about 20,000 files, at 0.2–1 MB each.
   - Files are subset to the bounding box and cached in `data/stage4/*.npz`, which the workflow commits.
   - If a run hits `MAX_FETCH` or the time limit, run it again and it resumes from the cache.
4. After that, it runs monthly and only fetches new hours.

To run locally: `pip install requests numpy pygrib && python scripts/build_data.py`, then `python -m http.server`.

## Caveats

- **Reporting bias.** Flooding is reported where people see it, so the secured interiors of CCSFS and KSC are under-sampled.
- **Trigger is an upper bound.** There's no record of heavy rain that *didn't* flood a spot, so the true threshold may be lower.
- **Stage IV resolution.** Stage IV is about 4.7 km, so a single convective core on the Cape can be smeared out or missed.
- **Hourly QC.** Hourly Stage IV gets less QC than the 6-hourly product.
- **Short-duration peaks.** The peak 1–3 h burst before a report isn't necessarily what caused the flooding. After a long soak, the 24–72 h totals are the better guide.
- **Report times are approximate.** Local storm reports carry the time someone reported the flooding, and Storm Events times are often estimates. The ±12 h / +2 h search window allows for that.
- **Mixed flood types.** Tropical events can mix lagoon or surge water with rainfall flooding. `COASTAL FLOOD` and `STORM SURGE` reports are excluded.

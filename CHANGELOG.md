# Changelog

## 2026-05-17 / 2026-05-18 — Local refactor and pipeline run (pbb2291 + Claude)

### Fixed: `correct_temps()` in `src/helper/Functions_TempCorrection.py`

The temperature correction function had three bugs introduced during a prior refactor (adding 3D distance calculation and sky radiation term):

1. **Double merge** — function called `pixel_df.merge(imgtraj_df)` then immediately called `merge_imgtraj_pixel_dfs()` which merged again, creating `_x`/`_y` duplicate columns for `Altitude`, `AirTemp`, `RH`, `STRD`, `DateTime`.

2. **Column rename clobbered drone coordinates** — `merge_imgtraj_pixel_dfs()` renamed `Northing` → `Northing_pix` and `Easting` → `Easting_pix`, destroying the drone position columns needed for distance calculation, causing `KeyError: 'Northing'`.

3. **`phi`/`tau` computed on `df`, read from `pixel_df`** — signal and transmissivity were set on an intermediate `df` variable but `Temperature_eVeg` was computed from `pixel_df`, causing `KeyError: 'phi'`.

**Fix:** Eliminated the redundant `merge_imgtraj_pixel_dfs` call. Renamed drone Northing/Easting to `Northing_uav`/`Easting_uav` before merging to avoid collision with pixel coordinates. All computation now operates on a single `df`.

**Distance calculation:** Switched from 3D distance (which requires `Northing_pix`/`Easting_pix` pixel coordinates not present in the backup pixelpickles) to simple vertical distance `Altitude - Elevation`. The 3D version is preserved in comments for future use when pixel coordinates are available.

### Added: temperature correction run across all 6 sites

Applied `correct_temps()` to all six pixelpickle files, adding `Temperature_eVeg` column (corrected surface temperature in °C). The correction shifts temperatures ~2–4°C warmer than raw values, with ground pixels receiving the largest correction (emissivity 0.91) and woody vegetation the smallest (0.97).

Pixelpickle files (not tracked in git) located at `out/pixelpickles/`:
- BuffaloCamp_pixels.obj — 4.85M pixels
- Hlangwine_pixels.obj — 30.4M pixels
- Letaba_pixels.obj — 24.6M pixels
- Makhohlola_pixels.obj — 1.9M pixels
- Nkuhlu_pixels.obj — 4.7M pixels
- RoanCamp_pixels.obj — 61.5M pixels

### Fixed: `src/3_process_figs_stats/Figures_Stats_AllSites.ipynb`

Updated the cross-site statistics and figures notebook to run locally:

- **Temperature column:** switched all analysis from raw `Temperature` to `Temperature_eVeg`
- **Pickle path:** `./data/out/pixelpickles/` → `../../out/pixelpickles/`
- **Directory creation:** added `parents=True, exist_ok=True` to `figdir.mkdir()`; added automatic `stats/` dir creation
- **pandas 2.0 compatibility:** replaced `.iteritems()` with `.items()` (3 occurrences)
- **Image shapefile section:** commented out — requires raw thermal rasters from Harvard cluster

### Added: `src/3_process_figs_stats/Figures_Stats_AllSites_run.py`

Standalone Python script generated from the notebook via `nbconvert`, for running Stage 3 without Jupyter.

### Generated outputs

**Stats** (`src/3_process_figs_stats/stats/`): 89 CSV files including per-site descriptive stats, per-image Welch t-tests and Cohen's D effect sizes for temperature, height, and PAI, and all-sites summaries.

**Figures** (`src/3_process_figs_stats/figs/final_120722/AllSites/`): 55 PNG files including publication figures (Fig 2: KDE by exclosure and topo position; Fig 3: height × temperature scatter grids; effect size boxplots; vegetation-stratified PAI and temperature plots).

Key results: inside exclosures are consistently cooler than outside at all 6 sites. RoanCamp shows the largest effect (median ~3°C cooler inside). All sites show significantly taller vegetation inside exclosures.

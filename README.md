# ThermalLandscapes

Analysis of drone-based thermal camera imagery to study surface temperature differences in African savanna exclosures. Six sites in Kruger National Park (South Africa) with livestock/wildlife exclosure treatments are compared to assess how large herbivore exclusion affects landscape surface temperatures.

## Study Design

Each site contains exclosures that fence out different combinations of animals (Inside) versus unfenced control areas (Outside). Drone thermal imagery flown at ~100m altitude captures surface temperatures across these treatments. The goal is to quantify temperature differences between exclosure treatments and relate them to vegetation structure (height, plant area index).

**Sites:** BuffaloCamp, Hlangwine, Letaba, Makhohlola, Nkuhlu, RoanCamp

## Repository Structure

```
src/
  1_correct_images_for_weather/   # Stage 1: build image trajectories with weather data
  2_process_sites/                # Stage 2: per-site pixel extraction (runs on server)
      RectifiedImageProcessing_*.ipynb
  3_process_figs_stats/           # Stage 3: cross-site statistics and publication figures
      Figures_Stats_AllSites.ipynb
      Figures_Stats_AllSites_run.py   # standalone script version
  helper/
      Functions_TempCorrection.py     # radiative transfer temperature correction
      Functions_ThermalLandscapes.py  # utility functions

data/
  in/
    ImageTrajectories/    # per-flight GPS + weather CSV files
    Polygons_AD_current/  # exclosure boundary shapefiles

out/
  pixelpickles/   # per-site pixel DataFrames (gitignored, large files)
  figs/           # output figures
  raster/         # output rasters
```

## Pipeline

### Stage 1 — Weather/Trajectory Processing
Runs on the Harvard cluster. Merges drone flight trajectories with meteorological data (air temperature, relative humidity, downwelling thermal radiation). Outputs: `data/in/ImageTrajectories/*_ImageTrajectory_withWeather.csv`.

### Stage 2 — Per-Site Pixel Extraction
Runs on the Harvard cluster (requires raw rectified thermal `.tif` files and LiDAR products). For each site:
1. Finds thermal images that overlap exclosure polygons
2. Extracts pixel-level temperature, elevation, canopy height, and plant area index (PAI)
3. Classifies pixels into vegetation types (ground/grass/shrub/tree)
4. Saves raw pixel DataFrame to `out/pixelpickles/{site}_pixels.obj`
5. Applies temperature correction (see below)

### Temperature Correction
`Functions_TempCorrection.correct_temps()` applies a radiative transfer correction to raw camera temperatures accounting for:
- **Emissivity** by vegetation type (ground: 0.91, grass: 0.96, woody: 0.97)
- **Atmospheric transmissivity** (tau) based on drone-to-pixel distance, air temperature, and relative humidity
- **Sky thermal radiation** (downwelling longwave, STRD from weather data)

Output column: `Temperature_eVeg` (corrected surface temperature in °C).

### Stage 3 — Statistics and Figures
Runs locally from `src/3_process_figs_stats/`. Requires corrected pixelpickles for all 6 sites.

```bash
cd src/3_process_figs_stats
python Figures_Stats_AllSites_run.py
```

Outputs (relative to `src/3_process_figs_stats/`):
- `stats/` — per-site and all-sites CSVs: descriptive stats, Welch t-tests, Cohen's D effect sizes for temperature, height, and PAI
- `figs/final_120722/AllSites/` — publication figures

## Setup

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt
```

A `.env` file is used on the Harvard cluster for path configuration. Local runs use hardcoded relative paths in the Stage 3 script.

## Notes

- **Pixelpickles are not tracked in git** (too large). Place `.obj` files in `out/pixelpickles/` before running Stage 3.
- **Stage 2 cannot run locally** — it requires raw thermal rasters and LiDAR products on the Harvard cluster.
- The image shapefile section in `Figures_Stats_AllSites.ipynb` (cell ~30) also requires server raster access and is skipped locally.
- Site `Nkuhlu` may appear as `Nkuhlu` in pickle filenames depending on which backup copy is used.

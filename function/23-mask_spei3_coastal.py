"""
23-mask_spei3_coastal.py
========================
Mask global SPEI-3 rasters to the US coastline + 80 km ocean buffer and
save the cropped, masked NetCDF files.

Strategy
--------
The coastal mask is derived from the already-processed SPEI-48 files in
spatial_SPEI/spei_subsetting/ — those files were previously masked to the
exact tl_2019_us_coastline + 80 km buffer, so we simply read the valid-cell
footprint from one of those files and reuse it.  This avoids re-reading the
large shapefile and is fully consistent with the SPEI-48 pipeline.

If no SPEI-48 reference file is available, the script falls back to building
the mask directly from tl_2019_us_coastline.shp using shapely + pyproj.

Run
---
    python3 23-mask_spei3_coastal.py

Progress is printed to stdout.  Safe to re-run; already-processed files are
skipped.

Inputs
------
  SPEI-3 source files:
    <MALONE_LAB>/Research/Natural_CH4_CO2/Drought/ECMWF_DroughtIndices_Global/
    SPEI3_genlogistic_global_era5_moda_ref1991to2020_YYYYMM.nc

  SPEI-48 reference (for mask):
    <MALONE_LAB>/Research/WUE_CUE/spatial_SPEI/spei_subsetting/
    *_US_OCEAN_COAST_80km.nc   (any one file is enough)

  Coastline shapefile (fallback only):
    <MALONE_LAB>/Research/WUE_CUE/spatial_SPEI/coastline/
    tl_2019_us_coastline.shp

Output
------
  <MALONE_LAB>/Research/WUE_CUE/spatial_SPEI/SPEI-3/
  SPEI3_genlogistic_global_era5_moda_ref1991to2020_YYYYMM_US_OCEAN_COAST_80km.nc
"""

import os, sys, glob, time, logging
import numpy as np
import netCDF4 as nc
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# *** EDIT THESE IF YOUR MOUNT POINT IS DIFFERENT ***
# On macOS with Finder:   /Volumes/MaloneLab/...
# On Windows:             \\corellia.environment.yale.edu\MaloneLab\...
# ---------------------------------------------------------------------------
MALONE_ROOT = "/Volumes/MaloneLab"   # adjust if needed

INPUT_DIR   = f"{MALONE_ROOT}/Research/Natural_CH4_CO2/Drought/ECMWF_DroughtIndices_Global"
SPEI48_DIR  = f"{MALONE_ROOT}/Research/WUE_CUE/spatial_SPEI/spei_subsetting"
COAST_SHP   = f"{MALONE_ROOT}/Research/WUE_CUE/spatial_SPEI/coastline/tl_2019_us_coastline.shp"
OUTPUT_DIR  = f"{MALONE_ROOT}/Research/WUE_CUE/spatial_SPEI/SPEI-3"

BUFFER_KM   = 80
FILL        = np.float32(9.969209968386869e+36)

# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger()

# ---------------------------------------------------------------------------
def build_mask_from_spei48(spei48_dir, global_lat, global_lon):
    """
    Derive the coastal mask and cropped grid from an existing SPEI-48 file.
    Returns (mask2d, lat_crop, lon_crop, lat_idx, lon_idx).
    """
    ref_files = sorted(glob.glob(os.path.join(spei48_dir, "*_US_OCEAN_COAST_80km.nc")))
    if not ref_files:
        return None

    log.info(f"Deriving mask from SPEI-48 reference: {os.path.basename(ref_files[0])}")
    f    = nc.Dataset(ref_files[0])
    lat_c = np.array(f.variables["lat"][:])
    lon_c = np.array(f.variables["lon"][:])

    # Find the SPEI variable (name may vary)
    spei_var_name = next((v for v in f.variables if v.startswith("SPEI")), None)
    if spei_var_name is None:
        f.close(); return None

    data  = np.array(f.variables[spei_var_name][0, :, :])
    fv    = getattr(f.variables[spei_var_name], "_FillValue", 9.969e36)
    f.close()

    mask2d = np.isfinite(data) & (data != fv)
    log.info(f"  Mask: {mask2d.sum():,} valid coastal cells over {mask2d.shape}")

    lat_idx = np.array([np.argmin(np.abs(global_lat - v)) for v in lat_c])
    lon_idx = np.array([np.argmin(np.abs(global_lon - v)) for v in lon_c])
    return mask2d, lat_c, lon_c, lat_idx, lon_idx


def build_mask_from_shapefile(shp_path, global_lat, global_lon):
    """
    Fallback: build the 80-km buffer mask from the coastline shapefile.
    Requires shapely, pyproj, and pyshp.
    """
    try:
        import shapefile
        import shapely.geometry as sgeom
        import shapely.ops as sops
        import pyproj
        from scipy.ndimage import binary_dilation
    except ImportError as e:
        log.error(f"Fallback requires shapely, pyproj, pyshp, scipy: {e}")
        return None

    log.info("Building mask from coastline shapefile (fallback) ...")
    log.info("  Reading shapefile ...")
    sf          = shapefile.Reader(shp_path)
    WGS84       = "EPSG:4326"
    ALBERS      = "EPSG:5070"
    BUFFER_M    = BUFFER_KM * 1_000

    project_to   = pyproj.Transformer.from_crs(WGS84, ALBERS, always_xy=True).transform
    project_back = pyproj.Transformer.from_crs(ALBERS, WGS84, always_xy=True).transform

    # Crop grid to US extent
    lat_idx = np.where((global_lat >= 17) & (global_lat <= 76))[0]
    lon_idx = np.where((global_lon >= -180) & (global_lon <= -60))[0]
    lat_c   = global_lat[lat_idx]
    lon_c   = global_lon[lon_idx]
    nlat, nlon = len(lat_c), len(lon_c)
    lat_step = abs(float(global_lat[1] - global_lat[0]))
    lon_step = abs(float(global_lon[1] - global_lon[0]))

    coast_grid = np.zeros((nlat, nlon), dtype=bool)
    log.info("  Rasterising coastline segments ...")
    n = 0
    for shape in sf.iterShapes():
        for (x, y) in shape.points:
            if 17 <= y <= 76 and -180 <= x <= -60:
                r = int(round((lat_c.max() - y) / lat_step))
                c = int(round((x - lon_c.min()) / lon_step))
                if 0 <= r < nlat and 0 <= c < nlon:
                    coast_grid[r, c] = True
        n += 1
    log.info(f"  {n} shapes, {coast_grid.sum()} coastline cells")

    # Dilate: at 0.25° resolution, 80 km ≈ 3-4 cells
    radius = 4
    Y, X   = np.ogrid[-radius:radius+1, -radius:radius+1]
    disk   = (X**2 + Y**2) <= radius**2
    mask2d = binary_dilation(coast_grid, structure=disk)

    # Remove Great Lakes
    gl_lat = np.where((lat_c >= 41) & (lat_c <= 49.5))[0]
    gl_lon = np.where((lon_c >= -93) & (lon_c <= -75))[0]
    mask2d[np.ix_(gl_lat, gl_lon)] = False
    log.info(f"  Buffer applied: {mask2d.sum():,} coastal cells")
    return mask2d, lat_c, lon_c, lat_idx, lon_idx


def process_file(fp, mask2d, lat_c, lon_c, lat_idx, lon_idx, out_dir):
    stem    = Path(fp).stem
    outname = f"{stem}_US_OCEAN_COAST_80km.nc"
    outpath = os.path.join(out_dir, outname)
    if os.path.exists(outpath):
        return "skip"
    try:
        src      = nc.Dataset(fp)
        data     = np.array(src.variables["SPEI3"][0, :, :])[np.ix_(lat_idx, lon_idx)]
        tval     = np.array(src.variables["time"][:])
        tunit    = getattr(src.variables["time"], "units",    "days since 1900-01-01")
        tcal     = getattr(src.variables["time"], "calendar", "standard")
        src.close()

        data_out = np.where(mask2d, data.astype(np.float32), FILL)

        dst = nc.Dataset(outpath, "w", format="NETCDF4")
        dst.createDimension("time", 1)
        dst.createDimension("lat",  len(lat_c))
        dst.createDimension("lon",  len(lon_c))
        tv = dst.createVariable("time", "f8", ("time",));  tv.units = tunit; tv.calendar = tcal; tv[:] = tval
        la = dst.createVariable("lat",  "f4", ("lat",));   la.units = "degrees_north"; la[:] = lat_c
        lo = dst.createVariable("lon",  "f4", ("lon",));   lo.units = "degrees_east";  lo[:] = lon_c
        sp = dst.createVariable("SPEI3","f4", ("time","lat","lon"),
                                fill_value=FILL, zlib=True, complevel=4)
        sp.long_name  = "SPEI-3 US ocean coast 80 km buffer"
        sp.source     = os.path.basename(fp)
        sp.created    = datetime.now().isoformat()
        sp[0]         = data_out
        dst.close()
        return "ok"
    except Exception as e:
        return f"ERROR: {e}"


# ---------------------------------------------------------------------------
def main():
    log.info("=" * 65)
    log.info("23-mask_spei3_coastal.py — SPEI-3 US coastal masking")
    log.info(f"Input : {INPUT_DIR}")
    log.info(f"Output: {OUTPUT_DIR}")
    log.info("=" * 65)

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    # ---- Find input files --------------------------------------------------
    files = sorted(glob.glob(os.path.join(INPUT_DIR, "SPEI3_*.nc")))
    if not files:
        log.error(f"No SPEI3 files found in {INPUT_DIR}"); sys.exit(1)
    log.info(f"Found {len(files)} SPEI3 files")

    # ---- Read global grid from first file ----------------------------------
    ref       = nc.Dataset(files[0])
    global_lat = np.array(ref.variables["lat"][:])
    global_lon = np.array(ref.variables["lon"][:])
    ref.close()

    # ---- Build coastal mask ------------------------------------------------
    result = build_mask_from_spei48(SPEI48_DIR, global_lat, global_lon)
    if result is None:
        log.warning("SPEI-48 reference not found — falling back to shapefile")
        result = build_mask_from_shapefile(COAST_SHP, global_lat, global_lon)
    if result is None:
        log.error("Could not build coastal mask"); sys.exit(1)

    mask2d, lat_c, lon_c, lat_idx, lon_idx = result

    # ---- Process files -----------------------------------------------------
    t0 = time.time()
    ok = skip = err = 0
    for i, fp in enumerate(files, 1):
        status = process_file(fp, mask2d, lat_c, lon_c, lat_idx, lon_idx, OUTPUT_DIR)
        bn     = os.path.basename(fp)
        if status == "ok":
            ok += 1
            if ok % 25 == 0 or i == len(files):
                elapsed = time.time() - t0
                rate    = ok / elapsed
                remain  = (len(files) - i) / rate if rate > 0 else 0
                log.info(f"  {i:3d}/{len(files)}  {ok} written  "
                         f"({rate:.1f} files/s, ~{remain/60:.0f} min remaining)")
        elif status == "skip":
            skip += 1
        else:
            log.error(f"  {bn}: {status}")
            err += 1

    log.info("=" * 65)
    log.info(f"Done: {ok} written, {skip} skipped, {err} errors")
    log.info(f"Output: {OUTPUT_DIR}")
    log.info(f"Total time: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()

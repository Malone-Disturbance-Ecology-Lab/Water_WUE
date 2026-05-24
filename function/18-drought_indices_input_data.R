# ============================================================
# Daily met drivers for ALL sites in site_lat_long_climate.csv
# Source: NASA POWER (AG) — with HTTP fallback if package absent
# Writes one CSV per site (named as site_name) to the output dir
# ============================================================

# ---- Setup & dependencies ----
.use_power_pkg <- requireNamespace("nasapower", quietly = TRUE)
if (.use_power_pkg) library(nasapower)

need_pkgs <- c("jsonlite","httr","dplyr","readr","lubridate","purrr","stringr","tidyr")
for (p in need_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
library(jsonlite); library(httr); library(dplyr); library(readr)
library(lubridate); library(purrr); library(stringr); library(tidyr)

# ---- Constants & helpers ----
POWER_FIRST_YEAR <- 1981L
POWER_LAST_YEAR  <- year(Sys.Date())

.es_hPa <- function(Tc) { 6.112 * exp(17.67 * Tc / (Tc + 243.5)) }  # Tetens (hPa)
sigma_SB <- 5.670374419e-8
MJday_to_Wm2 <- function(x) x * 1e6 / 86400

.compute_Rn <- function(SW_MJ, LW_MJ, Tair_C, alpha = 0.23, emiss = 0.98){
  Rsd_Wm2 <- MJday_to_Wm2(SW_MJ)           # SW down (W m-2)
  Rld_Wm2 <- MJday_to_Wm2(LW_MJ)           # LW down (W m-2)
  Ts_K    <- Tair_C + 273.15
  Rlu_Wm2 <- emiss * sigma_SB * (Ts_K^4)   # LW up (W m-2)
  (Rsd_Wm2 * (1 - alpha)) + (Rld_Wm2 - Rlu_Wm2)
}

# Simple retry wrapper for HTTP calls (handles transient network hiccups)
.with_retry <- function(expr, attempts = 3, wait_sec = 2){
  for (i in seq_len(attempts)){
    ok <- try(eval.parent(substitute(expr)), silent = TRUE)
    if (!inherits(ok, "try-error")) return(ok)
    if (i < attempts) Sys.sleep(wait_sec * i)
  }
  stop(ok)
}

# ---------- Option A: nasapower route ----------
.fetch_power_pkg <- function(lon, lat, start_date, end_date){
  nasapower::get_power(
    community    = "AG",
    lonlat       = c(lon, lat),
    pars         = c("T2M","RH2M","PS","WS2M",
                     "ALLSKY_SFC_SW_DWN","ALLSKY_SFC_LW_DWN","PRECTOTCORR"),
    dates        = c(start_date, end_date),
    temporal_api = "daily"
  )
}

# ---------- Option B: HTTP fallback (no nasapower needed) ----------
.fetch_power_http <- function(lon, lat, start_date, end_date){
  start_ymd <- gsub("-", "", as.character(as.Date(start_date)))
  end_ymd   <- gsub("-", "", as.character(as.Date(end_date)))
  base <- "https://power.larc.nasa.gov/api/temporal/daily/point"
  params <- paste(
    "parameters=T2M,RH2M,PS,WS2M,ALLSKY_SFC_SW_DWN,ALLSKY_SFC_LW_DWN,PRECTOTCORR",
    "community=AG",
    paste0("longitude=", lon),
    paste0("latitude=", lat),
    paste0("start=", start_ymd),
    paste0("end=", end_ymd),
    "format=JSON",
    sep="&"
  )
  url <- paste0(base, "?", params)
  resp <- .with_retry(httr::GET(url, httr::user_agent("R-nasapower-fallback")))
  httr::stop_for_status(resp)
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  js  <- jsonlite::fromJSON(txt)
  pars <- js$properties$parameter
  dates <- names(pars$T2M)
  tibble(
    DateTime = as.Date(dates, format="%Y%m%d"),
    T2M      = as.numeric(unlist(pars$T2M)[dates]),
    RH2M     = as.numeric(unlist(pars$RH2M)[dates]),
    PS       = as.numeric(unlist(pars$PS)[dates]),
    WS2M     = as.numeric(unlist(pars$WS2M)[dates]),
    SW_MJ    = as.numeric(unlist(pars$ALLSKY_SFC_SW_DWN)[dates]),
    LW_MJ    = as.numeric(unlist(pars$ALLSKY_SFC_LW_DWN)[dates]),
    PRECTOTCORR = as.numeric(unlist(pars$PRECTOTCORR)[dates])
  )
}

# ---------- One-site worker (selects pkg or HTTP) ----------
.fetch_one_site_POWER <- function(site, lat, lon, start_year, end_year, out_dir,
                                  alpha = 0.23, emiss = 0.98){
  # Validate & clamp inputs
  if (!is.finite(lat) || !is.finite(lon)) stop(sprintf("[%s] invalid lat/lon.", site))
  sy <- as.integer(max(min(start_year, POWER_LAST_YEAR), POWER_FIRST_YEAR))
  ey <- as.integer(max(min(end_year,   POWER_LAST_YEAR), POWER_FIRST_YEAR))
  if (ey < sy) stop(sprintf("[%s] end_year < start_year after clamping to POWER range.", site))
  
  message(sprintf("\n[%-7s] POWER lat=%.4f lon=%.4f years=%d–%d", site, lat, lon, sy, ey))
  start_date <- sprintf("%d-01-01", sy)
  end_date   <- sprintf("%d-12-31", ey)
  
  if (.use_power_pkg) {
    pw <- .with_retry(.fetch_power_pkg(lon, lat, start_date, end_date))
    if (!is.data.frame(pw) || nrow(pw) == 0) stop(sprintf("[%s] no POWER rows.", site))
    df <- pw %>%
      transmute(
        DateTime = as.Date(sprintf("%04d-%02d-%02d", YEAR, MM, DD)),
        Tair_C   = as.numeric(T2M),
        RH       = as.numeric(RH2M),
        PA_kPa   = as.numeric(PS),
        WS_ms    = as.numeric(WS2M),
        SW_MJ    = as.numeric(ALLSKY_SFC_SW_DWN),
        LW_MJ    = as.numeric(ALLSKY_SFC_LW_DWN),
        precip_mm= as.numeric(PRECTOTCORR)
      ) %>% arrange(DateTime)
  } else {
    message("  (nasapower not installed → using HTTP fallback)")
    df_raw <- .with_retry(.fetch_power_http(lon, lat, start_date, end_date))
    df <- df_raw %>%
      transmute(
        DateTime  = as.Date(DateTime),
        Tair_C    = as.numeric(T2M),
        RH        = as.numeric(RH2M),
        PA_kPa    = as.numeric(PS),
        WS_ms     = as.numeric(WS2M),
        SW_MJ     = as.numeric(SW_MJ),
        LW_MJ     = as.numeric(LW_MJ),
        precip_mm = as.numeric(PRECTOTCORR)
      ) %>% arrange(DateTime)
  }
  
  # Compute VPD (hPa)
  es <- .es_hPa(df$Tair_C)
  ea <- es * pmax(pmin(df$RH, 100), 0) / 100
  VPD_hPa <- pmax(es - ea, 0)
  
  # Net radiation (W m-2), daily mean
  Rn_Wm2 <- .compute_Rn(df$SW_MJ, df$LW_MJ, df$Tair_C, alpha = alpha, emiss = emiss)
  
  out <- tibble(
    DateTime  = df$DateTime,
    Tair_f    = df$Tair_C,   # °C
    VPD_f     = VPD_hPa,     # hPa
    PA_f      = df$PA_kPa,   # kPa
    NETRAD_f  = Rn_Wm2,      # W m^-2
    WS        = df$WS_ms,    # m s^-1
    precip_mm = df$precip_mm # mm/day
  )
  
  out_path <- file.path(out_dir, paste0(site, ".csv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(out, out_path)
  message(sprintf("[%s] Wrote %s (%d rows).", site, out_path, nrow(out)))
  tibble(site = site, rows = nrow(out), file = out_path)
}

# ---------- MAIN: process ALL sites in CSV ----------
fetch_all_sites_power_to_csv <- function(
    info_csv = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\data_products\\info\\site_lat_long_climate.csv",
    out_dir  = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\drivers\\ameri_drivers\\PET_drought\\input_data"
){
  message("[SETUP] Output dir: ", out_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- readr::read_csv(info_csv, show_col_types = FALSE)
  need <- c("site_name","lat","long","start_date","end_date")
  miss <- setdiff(need, names(meta))
  if (length(miss)) stop("Missing columns in info CSV: ", paste(miss, collapse=", "))
  
  # Optional per-site albedo/emissivity (columns may be absent; not required)
  if (!"albedo" %in% names(meta))     meta$albedo <- NA_real_
  if (!"emissivity" %in% names(meta)) meta$emissivity <- NA_real_
  
  meta <- meta %>%
    mutate(
      start_date = as.integer(start_date),
      end_date   = as.integer(end_date),
      lat  = as.numeric(lat),
      long = as.numeric(long)
    ) %>%
    filter(is.finite(lat), is.finite(long),
           is.finite(start_date), is.finite(end_date)) %>%
    distinct(site_name, .keep_all = TRUE) %>%
    arrange(site_name)
  
  message(sprintf("[SETUP] %d site(s) to process.", nrow(meta)))
  
  # Iterate all sites; continue if one fails (log status)
  results <- meta %>%
    mutate(
      status = NA_character_,
      file   = NA_character_,
      rows   = NA_integer_
    )
  
  for (i in seq_len(nrow(meta))){
    s <- meta$site_name[i]
    la <- meta$lat[i]
    lo <- meta$long[i]
    sy <- meta$start_date[i]
    ey <- meta$end_date[i]
    al <- suppressWarnings(as.numeric(meta$albedo[i]))
    em <- suppressWarnings(as.numeric(meta$emissivity[i]))
    
    # defaults if missing
    if (!is.finite(al)) al <- 0.23
    if (!is.finite(em)) em <- 0.98
    
    message(sprintf("\n---- [%d/%d] %s ----", i, nrow(meta), s))
    res <- try(
      .fetch_one_site_POWER(
        site = s, lat = la, lon = lo,
        start_year = sy, end_year = ey,
        out_dir = out_dir,
        alpha = al, emiss = em
      ),
      silent = TRUE
    )
    
    if (inherits(res, "try-error")){
      message(sprintf("[WARN] %s failed: %s", s, as.character(res)))
      results$status[i] <- "failed"
    } else {
      results$status[i] <- "ok"
      results$file[i]   <- res$file[1]
      results$rows[i]   <- res$rows[1]
    }
  }
  
  message("\n[SUMMARY]")
  print(results %>% select(site_name, status, rows, file), n = nrow(results))
  invisible(results)
}

# ---------- CALL: process ALL sites in the CSV ----------
results <- fetch_all_sites_power_to_csv(
  info_csv = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\data_products\\info\\site_lat_long_climate.csv",
  out_dir  = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\drivers\\ameri_drivers\\PET_drought\\input_data"
)

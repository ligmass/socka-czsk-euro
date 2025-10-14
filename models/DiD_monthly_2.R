#!/usr/bin/env Rscript
# ===== DiD monthly (CSV-only, fixest-or-lm, with log imports/exports) =====

options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"),
  Ncpus = max(1, parallel::detectCores() %/% 2),
  warn  = 1
)
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
need <- c("data.table", "lubridate", "sandwich", "lmtest", "broom", "fixest")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) try(install.packages(miss, dependencies = FALSE), silent = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(sandwich)
  library(lmtest);     library(broom)
  has_fixest <- require(fixest, quietly = TRUE)
})

dir.create("tables2", showWarnings = FALSE, recursive = TRUE)

# ---- load data (CSV) ---------------------------------------------------------
paths <- c(
  Sys.getenv("DID_DATA", ""),
  "data/monthly_panel_clean.csv",
  "data/panel.csv",
  "monthly_panel_clean.csv"
)
paths <- paths[nzchar(paths)]
path  <- paths[file.exists(paths)][1]
if (is.na(path) || !nzchar(path)) stop("no data csv found")

dat <- fread(path)
setnames(
  dat, names(dat),
  tolower(gsub("[^A-Za-z0-9]+", "_", names(dat)))
)
stopifnot("country" %in% names(dat), "date" %in% names(dat))

# ---- date -> base Date -------------------------------------------------------
if (!inherits(dat$date, "Date")) {
  if (is.numeric(dat$date)) {
    y <- dat$date %/% 100
    mo <- sprintf("%02d", dat$date %% 100)
    dat[, date := as.Date(paste0(y, "-", mo, "-01"))]
  } else {
    s <- sub("^([0-9]{4}-[0-9]{2}).*$", "\\1-01", as.character(dat$date))
    dat[, date := as.Date(s)]
  }
}
if (!inherits(dat$date, "Date")) stop("failed to parse date")
dat[, country := as.character(country)]

# ---- add log imports / exports (auto-detect columns) -------------------------
# We try common names; if multiple, pick the first match of each type.
imp_candidates <- grep("(^|_)imp(|ort|orts)?$", names(dat), value = TRUE)
imp_candidates <- c(imp_candidates,
                    grep("imports?$|trade_?imports?$", names(dat), value = TRUE))
exp_candidates <- grep("(^|_)exp(|ort|orts)?$", names(dat), value = TRUE)
exp_candidates <- c(exp_candidates,
                    grep("exports?$|trade_?exports?$", names(dat), value = TRUE))

imp_candidates <- unique(imp_candidates)
exp_candidates <- unique(exp_candidates)

if (length(imp_candidates)) {
  imp_col <- imp_candidates[1]
  dat[, ln_imp := log1p(pmax(get(imp_col), 0))]
}
if (length(exp_candidates)) {
  exp_col <- exp_candidates[1]
  dat[, ln_exp := log1p(pmax(get(exp_col), 0))]
}

# ---- treatment + rel ---------------------------------------------------------
adopt_date    <- as.Date("2009-01-01")
placebo_dates <- as.Date(c("2007-01-01", "2011-01-01"))  # set to character(0) to skip
dat[, treated := as.integer(country %in% c("sk", "SK"))]
dat[, post    := as.integer(date >= adopt_date)]
dat[, D       := treated * post]
dat[, rel     := (year(date) - year(adopt_date)) * 12 +
                 (month(date) - month(adopt_date))]

# ---- outcomes & controls -----------------------------------------------------
controls <- intersect(
  c("energy_yoy", "vat_change", "gfs_index", "ip_gap", "covid", "war22"),
  names(dat)
)

# prefer these outcomes if present; then add ln_imp/ln_exp if we created them
cands <- c(
  "hicp_yoy", "infl_yoy", "cpi_yoy", "unemp_rate", "ip_yoy",
  "minw_real", "credit_yoy", "ca_gdp", "nx_gdp"
)
extra_out <- intersect(c("ln_imp", "ln_exp"), names(dat))
outcomes <- unique(c(intersect(cands, names(dat)), extra_out))
if (!length(outcomes)) {
  bad <- c("country", "date", "treated", "post", "d", "rel", controls)
  outcomes <- setdiff(names(dat)[vapply(dat, is.numeric, logical(1))], bad)
}
cat("outcomes:", paste(outcomes, collapse = ", "), "\n")

# ---- helpers -----------------------------------------------------------------
rhs_join <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"

fit_twfe <- function(y) {
  if (has_fixest) {
    rhs <- rhs_join(c("D", controls))
    fml <- as.formula(paste(y, "~", rhs, "| country + date"))
    feols(fml, data = dat, cluster = ~ country + date)
  } else {
    rhs <- c("D", controls, "factor(country)", "factor(date)")
    fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
    m   <- lm(fml, data = dat)
    attr(m, "vcovCL") <- vcovCL(m, cluster = dat$date)
    m
  }
}

tidy_any <- function(m) {
  vc <- attr(m, "vcovCL")
  if (!has_fixest && !is.null(vc)) tidy(m, conf.int = TRUE, vcov = vc)
  else tidy(m, conf.int = TRUE)
}

gof_row <- function(m, y) {
  if (has_fixest) {
    r2  <- tryCatch(as.numeric(fixest::fitstat(m, "r2")),  error = function(e) NA_real_)
    wr2 <- tryCatch(as.numeric(fixest::fitstat(m, "wr2")), error = function(e) NA_real_)
    data.table(outcome = y, n = nobs(m), r2 = r2, r2_within = wr2,
               se_type = "cluster: country+date")
  } else {
    sm <- summary(m)
    data.table(outcome = y, n = nobs(m), r2 = sm$r.squared,
               adj_r2 = sm$adj.r.squared, se_type = "vcovCL by month")
  }
}





# ---- MAIN DiD (CSV) ----------------------------------------------------------
main_rows <- vector("list", length(outcomes)); names(main_rows) <- outcomes
gof_rows  <- vector("list", length(outcomes)); names(gof_rows)  <- outcomes
for (y in outcomes) {
  m  <- fit_twfe(y)
  tt <- tidy_any(m); tt$outcome <- y
  main_rows[[y]] <- tt
  gof_rows[[y]]  <- gof_row(m, y)
}
fwrite(rbindlist(main_rows, fill = TRUE), "tables2/did_main.csv")
fwrite(rbindlist(gof_rows,  fill = TRUE), "tables2/did_main_gof.csv")

# ---- EVENT STUDY (CSV) -------------------------------------------------------
es_rows <- vector("list", length(outcomes)); names(es_rows) <- outcomes
for (y in outcomes) {
  if (has_fixest) {
    rhs <- rhs_join(c("i(rel, treated, ref = -1)", controls))
    fml <- as.formula(paste(y, "~", rhs, "| country + date"))
    es  <- feols(fml, data = dat, cluster = ~ country + date)
    tt  <- tidy(es, conf.int = TRUE)
  } else {
    rel_fac <- factor(dat$rel)
    if (!("-1" %in% levels(rel_fac))) {
      rel_fac <- factor(rel_fac, levels = c("-1", levels(rel_fac)))
    }
    rel_fac <- stats::relevel(rel_fac, ref = "-1")
    df  <- data.frame(dat, rel_fac = rel_fac)
    rhs <- c("rel_fac:treated", controls, "factor(country)", "factor(date)")
    fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
    es  <- lm(fml, data = df)
    vc  <- vcovCL(es, cluster = df$date)
    tt  <- tidy(es, conf.int = TRUE, vcov = vc)
  }
  tt$outcome <- y
  es_rows[[y]] <- tt
}
fwrite(rbindlist(es_rows, fill = TRUE), "tables2/eventstudy_all.csv")

# ---- PLACEBOS (safe) ---------------------------------------------------------
placebo_dates <- as.Date(
  placebo_dates,
  tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y%m%d")
)
if (length(placebo_dates)) {
  pb_all <- list()
  for (pd in placebo_dates) {
    post_pb <- as.integer(dat$date >= pd)
    D_pb    <- dat$treated * post_pb
    dat_pb  <- data.frame(dat, D_pb = D_pb)
    for (y in outcomes) {
      if (has_fixest) {
        rhs <- rhs_join(c("D_pb", controls))
        fml <- as.formula(paste(y, "~", rhs, "| country + date"))
        fit <- feols(fml, data = dat_pb, cluster = ~ country + date)
        tt  <- tidy(fit, conf.int = TRUE)
      } else {
        rhs <- c("D_pb", controls, "factor(country)", "factor(date)")
        fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
        fit <- lm(fml, data = dat_pb)
        vc  <- vcovCL(fit, cluster = dat_pb$date)
        tt  <- tidy(fit, conf.int = TRUE, vcov = vc)
      }
      pd_date <- suppressWarnings(as.Date(pd))
      pd_str  <- if (!is.na(pd_date)) strftime(pd_date, "%Y-%m-%d")
                 else as.character(pd)
      pd_key  <- if (!is.na(pd_date)) strftime(pd_date, "%Y%m")
                 else gsub("\\D", "", as.character(pd))
      tt$outcome <- y
      tt$placebo_date <- pd_str
      pb_all[[paste0(y, "_", pd_key)]] <- tt
    }
  }
  fwrite(rbindlist(pb_all, fill = TRUE), "tables2/did_placebos.csv")
}

# ---- TREAT CZ (CSV) ----------------------------------------------------------
dat[, treated_cz := as.integer(country %in% c("cz", "CZ"))]
dat[, D_cz := treated_cz * post]
cz_rows <- vector("list", length(outcomes)); names(cz_rows) <- outcomes
for (y in outcomes) {
  if (has_fixest) {
    rhs <- rhs_join(c("D_cz", controls))
    fml <- as.formula(paste(y, "~", rhs, "| country + date"))
    fit <- feols(fml, data = dat, cluster = ~ country + date)
    tt  <- tidy(fit, conf.int = TRUE)
  } else {
    rhs <- c("D_cz", controls, "factor(country)", "factor(date)")
    fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
    fit <- lm(fml, data = dat)
    vc  <- vcovCL(fit, cluster = dat$date)
    tt  <- tidy(fit, conf.int = TRUE, vcov = vc)
  }
  tt$outcome <- y
  cz_rows[[y]] <- tt
}
fwrite(rbindlist(cz_rows, fill = TRUE), "tables2/did_treat_cz.csv")

library(fixest)

print_full_summary <- function(m, outcome, cluster_spec = ~ id + time) {
  cat("\n========================\nOutcome:", outcome, "\n")

  # 1) coefficients + three SE choices
  cat("\n-- Coefs w/ classical SEs --\n")
  print(summary(m, vcov = "iid"))

  cat("\n-- Coefs w/ HC (heteroskedasticity-robust) SEs --\n")
  print(summary(m, vcov = "hetero"))

  cat("\n-- Coefs w/ two-way cluster SEs --\n")
  print(summary(m, cluster = cluster_spec))

  # 2) GOF: AIC/BIC, R2 variants, LL, df, n, RMSE, etc.
  fs <- fixest::fitstat(m,
    ~ aic + bic + ll + df + n + rmse + r2 + r2_within + r2_adj
  )
  cat("\n-- Fit stats --\n")
  print(fs)
}

for (y in outcomes) {
  m <- fit_twfe(y)
  print_full_summary(m, y)
}


# ---- Optional OLS (set env vars to enable) -----------------------------------
ols_y <- Sys.getenv("OLS_Y", unset = "")
ols_x <- Sys.getenv("OLS_X", unset = "")
if (nzchar(ols_y) && nzchar(ols_x) && all(c(ols_y, ols_x) %in% names(dat))) {
  ss <- na.omit(dat[, .(date, y = get(ols_y), x = get(ols_x))])
  last_d   <- max(ss$date)
  cs       <- ss[date == last_d]
  ols_data <- if (nrow(cs) >= 20) cs else ss
  ols <- lm(y ~ x, data = ols_data)

  dw <- tryCatch(lmtest::dwtest(ols)$statistic[[1]],
                 error = function(e) NA_real_)
  jb <- tryCatch(tseries::jarque.bera.test(residuals(ols)),
                 error = function(e) NULL)
  sk <- tryCatch(moments::skewness(residuals(ols)),
                 error = function(e) NA_real_)
  ku <- tryCatch(moments::kurtosis(residuals(ols)),
                 error = function(e) NA_real_)

  jb_stat <- if (!is.null(jb)) unname(jb$statistic) else NA_real_
  jb_p    <- if (!is.null(jb)) unname(jb$p.value)   else NA_real_

  fwrite(tidy(ols, conf.int = TRUE), "tables2/ols_coeffs.csv")
  fwrite(
    data.table(
      metric = c("Durbin-Watson", "Jarque-Bera",
                 "Prob(JB)", "Skew", "Kurtosis"),
      value  = c(dw, jb_stat, jb_p, sk, ku)
    ),
    "tables2/ols_diag.csv"
  )
}

cat(
  "done. see /tables2: did_main.csv, did_main_gof.csv, eventstudy_all.csv, ",
  "did_placebos.csv (if any), did_treat_cz.csv, OLS files if enabled.\n",
  sep = ""
)

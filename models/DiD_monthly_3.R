#!/usr/bin/env Rscript
# ===== DiD monthly (CSV-only; prints+files; full metrics; spec log) =====

options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"),
  Ncpus = max(1, parallel::detectCores() %/% 2),
  warn  = 1
)
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
need <- c("data.table","lubridate","sandwich","lmtest","broom","fixest")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) try(install.packages(miss, dependencies = FALSE), silent = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(sandwich)
  library(lmtest);     library(broom)
  has_fixest <- require(fixest, quietly = TRUE)
  if (has_fixest) fixest::setFixest_notes(FALSE)
})

outdir <- "tables3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------ data ------------------------------------------------------------
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
setnames(dat, names(dat), tolower(gsub("[^A-Za-z0-9]+", "_", names(dat))))
stopifnot("country" %in% names(dat), "date" %in% names(dat))

# date to base Date
if (!inherits(dat$date, "Date")) {
  if (is.numeric(dat$date)) {
    y  <- dat$date %/% 100
    mo <- sprintf("%02d", dat$date %% 100)
    dat[, date := as.Date(paste0(y, "-", mo, "-01"))]
  } else {
    s <- sub("^([0-9]{4}-[0-9]{2}).*$", "\\1-01", as.character(dat$date))
    dat[, date := as.Date(s)]
  }
}
if (!inherits(dat$date, "Date")) stop("failed to parse date")
dat[, country := as.character(country)]

### CHANGED: always create a numeric time index for FE + clustering
dat[, time_id := as.integer(date)]

# ------------ alias logs if only log_* exist ---------------------------------
if (!"ln_imp" %in% names(dat) && "log_imp" %in% names(dat)) dat[, ln_imp := log_imp]
if (!"ln_exp" %in% names(dat) && "log_exp" %in% names(dat)) dat[, ln_exp := log_exp]

# ------------ treatment + rel -------------------------------------------------
adopt_date    <- as.Date("2009-01-01")
placebo_dates <- as.Date(c("2007-01-01", "2011-01-01"))
dat[, treated := as.integer(country %in% c("sk", "SK"))]
dat[, post    := as.integer(date >= adopt_date)]
dat[, D       := treated * post]
dat[, rel     := (year(date) - year(adopt_date)) * 12 + (month(date) - month(adopt_date))]

# ------------ controls & outcomes --------------------------------------------
controls <- intersect(c("energy_yoy","vat_change","gfs_index","ip_gap","covid","war22"), names(dat))

preferred <- c(
  "hicp_yoy","unemp_rate","ln_imp","ln_exp","log_imp","log_exp",
  "hicp","imports_world_meur","exports_world_meur",
  "repo","mro","fx_eur","fx_vol",
  "ip_yoy","infl_yoy","cpi_yoy","minw_real","credit_yoy","ca_gdp","nx_gdp"
)
if (all(c("ln_imp","log_imp") %in% names(dat))) preferred <- setdiff(preferred, "log_imp")
if (all(c("ln_exp","log_exp") %in% names(dat))) preferred <- setdiff(preferred, "log_exp")

exclude <- c("country","date","time_id","treated","post","rel","d","D","df","treated_cz","D_cz")
num_vars <- names(dat)[vapply(dat, is.numeric, logical(1))]
num_vars_clean <- setdiff(num_vars, exclude)

outcomes <- unique(c(
  intersect(preferred, names(dat)),
  setdiff(num_vars_clean, controls)
))
outcomes <- setdiff(outcomes, "date")

cat("outcomes used:", paste(outcomes, collapse = ", "), "\n")
cat("controls used:", ifelse(length(controls), paste(controls, collapse=", "), "(none)"), "\n")

# log the spec used
spec_txt <- file.path(outdir, "spec_used.txt")
writeLines(c(
  paste0("DATA_PATH: ", path),
  paste0("ADOPT_DATE: ", as.character(adopt_date)),
  paste0("PLACEBOS: ", paste(placebo_dates, collapse = ", ")),
  paste0("TREATED_DEF: country in {SK}"),
  paste0("FIXED_EFFECTS: country + time_id"),
  paste0("CLUSTERS: country + time_id (or time_id only)"),
  paste0("CONTROLS: ", ifelse(length(controls), paste(controls, collapse = ", "), "(none)")),
  paste0("OUTCOMES: ", paste(outcomes, collapse = ", "))
), con = spec_txt)

# ------------ helpers ---------------------------------------------------------
rhs_join <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"

### CHANGED: pre-filter to complete cases to avoid NA regressions
prepare_dfm <- function(y) {
  need <- unique(c(y, "D", "country", "time_id", controls))
  dfm  <- dat[, ..need]
  bad  <- rowSums(is.na(dfm)) > 0
  if (any(bad)) dfm <- dfm[!bad]
  if (nrow(dfm) < 10) stop("Too few complete observations for outcome: ", y)
  dfm
}

### CHANGED: use time_id for FE + clustering; robust fixest fallback ladder
fit_twfe <- function(y) {
  dfm <- prepare_dfm(y)
  rhs <- rhs_join(c("D", controls))
  fml_no_fe <- as.formula(paste(y, "~", rhs))

  if (has_fixest) {
    try1 <- try(fixest::feols(as.formula(paste(y, "~", rhs, "| country + time_id")),
                               data = dfm, cluster = ~ country + time_id), silent = TRUE)
    if (!inherits(try1, "try-error")) return(try1)

    try2 <- try(fixest::feols(fml_no_fe, data = dfm,
                               fixef = ~ country + time_id, cluster = ~ country + time_id), silent = TRUE)
    if (!inherits(try2, "try-error")) return(try2)

    try3 <- try(fixest::feols(fml_no_fe, data = dfm,
                               fixef = c("country","time_id"), cluster = c("country","time_id")), silent = TRUE)
    if (!inherits(try3, "try-error")) return(try3)

    try4 <- try(fixest::feols(fml_no_fe, data = dfm, fixef = c("country","time_id")), silent = TRUE)
    if (!inherits(try4, "try-error")) return(try4)

    message("fixest failed on outcome ", y, "; falling back to lm with clustered SE.")
  }

  # lm fallback with FE via factors; cluster by time_id
  fml_lm <- as.formula(paste(y, "~", paste(c("D", controls, "factor(country)", "factor(time_id)"), collapse = " + ")))
  m <- lm(fml_lm, data = dfm)
  attr(m, "dfm") <- dfm  # keep for clustering
  m
}

tidy_any <- function(m) {
  vc <- attr(m, "vcovCL")
  if (!has_fixest && !is.null(vc)) tidy(m, conf.int = TRUE, vcov = vc) else tidy(m, conf.int = TRUE)
}

# Only show treatment + controls in printouts (hide FE)
keep_mask <- function(nms) {
  base <- grepl("^D$", nms)
  if (length(controls)) base <- base | nms %in% controls
  base
}

# fitstat helper
fs_safe <- function(m) {
  if (!inherits(m, "fixest")) return(NULL)
  avail <- tryCatch(fixest::fitstat(m, show_types = TRUE), error = function(e) NULL)
  want  <- c("aic","bic","ll","rmse","r2","wr2","ar2","n")  # 'n' instead of 'nobs'
  use   <- intersect(want, if (is.null(avail)) want else rownames(avail))
  if (!length(use)) return(NULL)
  fixest::fitstat(m, as.formula(paste("~", paste(use, collapse = " + "))))
}

# GOF row
gof_row <- function(m, y) {
  if (inherits(m, "fixest")) {
    get <- function(type) tryCatch(as.numeric(fixest::fitstat(m, type)), error = function(e) NA_real_)
    data.table(
      outcome   = y,
      n         = get("n"),
      r2        = get("r2"),
      r2_within = { w <- get("wr2"); if (is.na(w)) NA_real_ else w },
      rmse      = get("rmse"),
      aic       = get("aic"),
      bic       = get("bic"),
      se_type   = "cluster: country+time_id"
    )
  } else {
    sm <- summary(m)
    data.table(
      outcome   = y,
      n         = nobs(m),
      r2        = sm$r.squared,
      r2_within = NA_real_,
      rmse      = sqrt(mean(sm$residuals^2)),
      aic       = AIC(m),
      bic       = BIC(m),
      se_type   = "HC/IID shown in logs"
    )
  }
}

### CHANGED: prints that filter out FE + robust clustering for lm fallback
print_full_summary <- function(m, outcome) {
  cat("\n========================\nOutcome:", outcome, "\n")

  if (inherits(m, "fixest")) {
    cat("\n-- Coefs (IID) [D + controls] --\n")
    s1 <- summary(m, vcov = "iid")
    print(s1$coeftable[keep_mask(rownames(s1$coeftable)), , drop = FALSE])

    cat("\n-- Coefs (HC) [D + controls] --\n")
    s2 <- summary(m, vcov = "hetero")
    print(s2$coeftable[keep_mask(rownames(s2$coeftable)), , drop = FALSE])

    cat("\n-- Coefs (cluster: country + time_id) [D + controls] --\n")
    s3 <- summary(m, cluster = ~ country + time_id)
    print(s3$coeftable[keep_mask(rownames(s3$coeftable)), , drop = FALSE])

    cat("\n-- Coefs (cluster: time_id) [D + controls] --\n")
    s4 <- summary(m, cluster = ~ time_id)
    print(s4$coeftable[keep_mask(rownames(s4$coeftable)), , drop = FALSE])

    cat("\n-- Fit stats --\n"); print(fs_safe(m))

  } else {
    dfm <- attr(m, "dfm")
    stopifnot(!is.null(dfm))
    # IID
    cat("\n-- Coefs (IID) [D + controls] --\n")
    ct_iid <- lmtest::coeftest(m)
    print(ct_iid[keep_mask(rownames(ct_iid)), , drop = FALSE])

    # HC
    cat("\n-- Coefs (HC) [D + controls] --\n")
    vc_hc  <- sandwich::vcovHC(m, type = "HC1")
    ct_hc  <- lmtest::coeftest(m, vcov. = vc_hc)
    print(ct_hc[keep_mask(rownames(ct_hc)), , drop = FALSE])

    # Two-way cluster (country + time_id) – lengths always match since we use model frame vars
    cat("\n-- Coefs (cluster: country + time_id) [D + controls] --\n")
    mf <- model.frame(m)
    vc2 <- sandwich::vcovCL(m, cluster = data.frame(country = mf$`factor(country)`, time = mf$`factor(time_id)`))
    ct2 <- lmtest::coeftest(m, vcov. = vc2)
    print(ct2[keep_mask(rownames(ct2)), , drop = FALSE])

    # Time-only cluster
    cat("\n-- Coefs (cluster: time_id) [D + controls] --\n")
    vc_t <- sandwich::vcovCL(m, cluster = mf$`factor(time_id)`)
    ct_t <- lmtest::coeftest(m, vcov. = vc_t)
    print(ct_t[keep_mask(rownames(ct_t)), , drop = FALSE])

    # Fit stats (lm)
    sm <- summary(m)
    print(list(
      nobs  = stats::nobs(m),
      r2    = sm$r.squared,
      r2_adj= sm$adj.r.squared,
      aic   = AIC(m),
      bic   = BIC(m)
    ))
  }
}

# ------------ head(dat) before DiD --------------------------------------------
cat("\n--- head(dat) (first 10 rows) ---\n")
print(utils::head(dat, 10))
fwrite(utils::head(dat, 50), file.path(outdir, "data_head.csv"))

# ------------ MAIN DiD (CSV + printable logs) --------------------------------
main_rows <- vector("list", length(outcomes)); names(main_rows) <- outcomes
gof_rows  <- vector("list", length(outcomes)); names(gof_rows)  <- outcomes
diag_rows <- list()

for (y in outcomes) {
  m <- fit_twfe(y)

  kept    <- names(coef(m))
  intended<- unique(c("D", controls))
  dropped <- setdiff(intended, kept)

  cat("\n[diagnostics:", y, "] kept:   ", ifelse(length(kept), paste(kept, collapse=", "), "(none)"),
      "\n[diagnostics:", y, "] dropped:", ifelse(length(dropped), paste(dropped, collapse=", "), "(none)"), "\n")

  diag_rows[[length(diag_rows) + 1]] <- data.table(outcome=y,
                                                   kept=if (length(kept)) paste(kept, collapse=";") else "(none)",
                                                   dropped=if (length(dropped)) paste(dropped, collapse=";") else "(none)")

  # per-outcome log (only D + controls printed)
  ftxt <- file.path(outdir, paste0("did_full_", y, ".txt"))
  con  <- file(ftxt, open = "wt")
  sink(con, split = TRUE); on.exit({sink(); close(con)}, add = TRUE)
  print_full_summary(m, y)
  sink(); close(con)

  # CSV outputs
  tt <- tidy_any(m); tt$outcome <- y
  main_rows[[y]] <- tt
  gof_rows[[y]]  <- gof_row(m, y)
}

fwrite(rbindlist(main_rows, fill=TRUE), file.path(outdir, "did_main.csv"))
fwrite(rbindlist(gof_rows,  fill=TRUE), file.path(outdir, "did_main_gof.csv"))
fwrite(rbindlist(diag_rows, fill=TRUE), file.path(outdir, "kept_dropped_controls.csv"))

# ------------ EVENT STUDY (CSV) ----------------------------------------------
es_rows <- vector("list", length(outcomes)); names(es_rows) <- outcomes
for (y in outcomes) {
  dfm <- prepare_dfm(y)
  if (has_fixest) {
    rhs <- rhs_join(c("i(rel, treated, ref = -1)", controls))
    fml <- as.formula(paste(y, "~", rhs, "| country + time_id"))
    es  <- feols(fml, data = dfm, cluster = ~ country + time_id)
    tt  <- tidy(es, conf.int = TRUE)
  } else {
    rel_fac <- factor(dat$rel)
    if (!("-1" %in% levels(rel_fac))) rel_fac <- factor(rel_fac, levels = c("-1", levels(rel_fac)))
    rel_fac <- stats::relevel(rel_fac, ref = "-1")
    df2 <- data.frame(dfm, rel_fac = rel_fac)
    rhs <- c("rel_fac:treated", controls, "factor(country)", "factor(time_id)")
    fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
    es  <- lm(fml, data = df2)
    vc  <- vcovCL(es, cluster = data.frame(country=df2$country, time=df2$time_id))
    tt  <- tidy(es, conf.int = TRUE, vcov = vc)
  }
  tt$outcome <- y
  es_rows[[y]] <- tt
}
fwrite(rbindlist(es_rows, fill = TRUE), file.path(outdir, "eventstudy_all.csv"))

# ------------ PLACEBOS (CSV) -------------------------------------------------
placebo_dates <- as.Date(placebo_dates, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y%m%d"))
if (length(placebo_dates)) {
  pb_all <- list()
  for (pd in placebo_dates) {
    post_pb <- as.integer(dat$date >= pd)
    D_pb    <- dat$treated * post_pb
    dat_pb  <- copy(dat)[, D_pb := D_pb]
    for (y in outcomes) {
      dfm <- dat_pb[, ..unique(c(y,"D_pb","country","time_id",controls))]
      dfm <- dfm[complete.cases(dfm)]
      if (nrow(dfm) < 10) next
      if (has_fixest) {
        rhs <- rhs_join(c("D_pb", controls))
        fml <- as.formula(paste(y,"~",rhs,"| country + time_id"))
        fit <- feols(fml, data = dfm, cluster = ~ country + time_id)
        tt  <- tidy(fit, conf.int=TRUE)
      } else {
        rhs <- c("D_pb",controls,"factor(country)","factor(time_id)")
        fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
        fit <- lm(fml, data=dfm)
        vc  <- vcovCL(fit, cluster = data.frame(country=dfm$country, time=dfm$time_id))
        tt  <- tidy(fit, conf.int=TRUE, vcov=vc)
      }
      pd_str <- strftime(pd, "%Y-%m-%d")
      pd_key <- strftime(pd, "%Y%m")
      tt$outcome <- y
      tt$placebo_date <- pd_str
      pb_all[[paste0(y,"_",pd_key)]] <- tt
    }
  }
  if (length(pb_all)) fwrite(rbindlist(pb_all, fill=TRUE), file.path(outdir, "did_placebos.csv"))
}

# ------------ TREAT CZ (CSV) -------------------------------------------------
dat[, treated_cz := as.integer(country %in% c("cz", "CZ"))]
dat[, D_cz := treated_cz * post]
cz_rows <- vector("list", length(outcomes)); names(cz_rows) <- outcomes
for (y in outcomes) {
  dfm <- dat[, ..unique(c(y,"D_cz","country","time_id",controls))]
  dfm <- dfm[complete.cases(dfm)]
  if (nrow(dfm) < 10) next
  if (has_fixest) {
    rhs <- rhs_join(c("D_cz",controls))
    fml <- as.formula(paste(y,"~",rhs,"| country + time_id"))
    fit <- feols(fml, data = dfm, cluster = ~ country + time_id)
    tt  <- tidy(fit, conf.int=TRUE)
  } else {
    rhs <- c("D_cz",controls,"factor(country)","factor(time_id)")
    fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
    fit <- lm(fml, data=dfm)
    vc  <- vcovCL(fit, cluster = data.frame(country=dfm$country, time=dfm$time_id))
    tt  <- tidy(fit, conf.int=TRUE, vcov=vc)
  }
  tt$outcome <- y; cz_rows[[y]] <- tt
}
fwrite(rbindlist(cz_rows, fill=TRUE), file.path(outdir, "did_treat_cz.csv"))

# ------------ Optional OLS via env vars --------------------------------------
ols_y <- Sys.getenv("OLS_Y", unset = "")
ols_x <- Sys.getenv("OLS_X", unset = "")
if (nzchar(ols_y) && nzchar(ols_x) && all(c(ols_y, ols_x) %in% names(dat))) {
  ss <- na.omit(dat[, .(date, y = get(ols_y), x = get(ols_x))])
  last_d   <- max(ss$date)
  cs       <- ss[date == last_d]
  ols_data <- if (nrow(cs) >= 20) cs else ss
  ols <- lm(y ~ x, data = ols_data)

  dw <- tryCatch(lmtest::dwtest(ols)$statistic[[1]], error = function(e) NA_real_)
  jb <- tryCatch(tseries::jarque.bera.test(residuals(ols)), error = function(e) NULL)
  sk <- tryCatch(moments::skewness(residuals(ols)), error = function(e) NA_real_)
  ku <- tryCatch(moments::kurtosis(residuals(ols)), error = function(e) NA_real_)
  jb_stat <- if (!is.null(jb)) unname(jb$statistic) else NA_real_
  jb_p    <- if (!is.null(jb)) unname(jb$p.value)   else NA_real_

  fwrite(tidy(ols, conf.int = TRUE), file.path(outdir, "ols_coeffs.csv"))
  fwrite(data.table(metric=c("Durbin-Watson","Jarque-Bera","Prob(JB)","Skew","Kurtosis"),
                    value=c(dw,jb_stat,jb_p,sk,ku)),
         file.path(outdir, "ols_diag.csv"))
}

cat(
  "\nDone. Outputs in /", outdir, ":\n",
  "- data_head.csv\n",
  "- spec_used.txt\n",
  "- did_main.csv, did_main_gof.csv\n",
  "- kept_dropped_controls.csv\n",
  "- eventstudy_all.csv\n",
  "- did_placebos.csv (if any)\n",
  "- did_treat_cz.csv\n",
  "- did_full_<outcome>.txt (per-outcome logs: only D + controls shown)\n",
  "- ols_coeffs.csv + ols_diag.csv (if OLS_* env vars set)\n",
  sep = ""
)

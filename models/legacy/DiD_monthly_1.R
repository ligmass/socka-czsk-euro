
#!/usr/bin/env Rscript
options(
  repos = c(CRAN="https://packagemanager.posit.co/cran/__linux__/jammy/latest"),
  Ncpus = max(1, parallel::detectCores() %/% 2), warn = 1
)
dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
need <- c("data.table","lubridate","sandwich","lmtest","broom","fixest")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) try(install.packages(miss, dependencies=FALSE), silent=TRUE)

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(sandwich)
  library(lmtest); library(broom)
  has_fixest <- require(fixest, quietly=TRUE)
})

dir.create("tables", showWarnings=FALSE, recursive=TRUE)

# ---- load data (CSV) ----
paths <- c(Sys.getenv("DID_DATA",""),
           "data/monthly_panel_clean.csv","data/panel.csv",
           "monthly_panel_clean.csv")
paths <- paths[nzchar(paths)]
path  <- paths[file.exists(paths)][1]
if (is.na(path) || !nzchar(path)) stop("no data csv found")
dat <- fread(path)
setnames(dat, names(dat), tolower(gsub("[^A-Za-z0-9]+","_",names(dat))))
stopifnot("country" %in% names(dat), "date" %in% names(dat))

# ---- date -> base Date ----
if (!inherits(dat$date,"Date")) {
  if (is.numeric(dat$date)) {
    y <- dat$date %/% 100; mo <- sprintf("%02d", dat$date %% 100)
    dat[, date := as.Date(paste0(y,"-",mo,"-01"))]
  } else {
    s <- sub("^([0-9]{4}-[0-9]{2}).*$","\\1-01", as.character(dat$date))
    dat[, date := as.Date(s)]
  }
}
if (!inherits(dat$date,"Date")) stop("failed to parse date")
dat[, country := as.character(country)]

# ---- treatment + rel ----
adopt_date    <- as.Date("2009-01-01")
placebo_dates <- as.Date(c("2007-01-01","2011-01-01"))  # set to character(0) to skip
dat[, treated := as.integer(country %in% c("sk","SK"))]
dat[, post    := as.integer(date >= adopt_date)]
dat[, D       := treated * post]
dat[, rel     := (year(date)-year(adopt_date))*12 + (month(date)-month(adopt_date))]

# ---- outcomes & controls ----
controls <- intersect(c("energy_yoy","vat_change","gfs_index","ip_gap","covid","war22"),
                      names(dat))
cands <- c("hicp_yoy","infl_yoy","cpi_yoy","unemp_rate","ip_yoy",
           "minw_real","credit_yoy","ca_gdp","nx_gdp")
outcomes <- intersect(cands, names(dat))
if (!length(outcomes)) {
  bad <- c("country","date","treated","post","d","rel",controls)
  outcomes <- setdiff(names(dat)[vapply(dat,is.numeric,logical(1))], bad)
}
cat("outcomes:", paste(outcomes, collapse=", "), "\n")

# ---- helpers ----
rhs_join <- function(x) if (length(x)) paste(x, collapse=" + ") else "1"
fit_twfe <- function(y){
  if (has_fixest){
    rhs <- rhs_join(c("D",controls))
    fml <- as.formula(paste(y,"~",rhs,"| country + date"))
    feols(fml, data=dat, cluster=~country+date)
  } else {
    rhs <- c("D",controls,"factor(country)","factor(date)")
    fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
    m <- lm(fml, data=dat); attr(m,"vcovCL") <- vcovCL(m, cluster=dat$date); m
  }
}
tidy_any <- function(m){
  vc <- attr(m,"vcovCL")
  if (!has_fixest && !is.null(vc)) tidy(m, conf.int=TRUE, vcov=vc)
  else tidy(m, conf.int=TRUE)
}

# ---- MAIN DiD (CSV) ----
main_rows <- vector("list", length(outcomes)); names(main_rows) <- outcomes
for (y in outcomes){ m <- fit_twfe(y); tt <- tidy_any(m); tt$outcome <- y; main_rows[[y]] <- tt }
fwrite(rbindlist(main_rows, fill=TRUE), "tables/did_main.csv")

# ---- EVENT STUDY (CSV) ----
es_rows <- vector("list", length(outcomes)); names(es_rows) <- outcomes
for (y in outcomes){
  if (has_fixest){
    rhs <- rhs_join(c("i(rel, treated, ref = -1)",controls))
    fml <- as.formula(paste(y,"~",rhs,"| country + date"))
    es <- feols(fml, data=dat, cluster=~country+date)
    tt <- tidy(es, conf.int=TRUE)
  } else {
    rel_fac <- factor(dat$rel)
    if (!("-1" %in% levels(rel_fac))) rel_fac <- factor(rel_fac, levels=c("-1",levels(rel_fac)))
    rel_fac <- stats::relevel(rel_fac, ref="-1")
    df <- data.frame(dat, rel_fac=rel_fac)
    rhs <- c("rel_fac:treated",controls,"factor(country)","factor(date)")
    fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
    es <- lm(fml, data=df); vc <- vcovCL(es, cluster=df$date)
    tt <- tidy(es, conf.int=TRUE, vcov=vc)
  }
  tt$outcome <- y; es_rows[[y]] <- tt
}
fwrite(rbindlist(es_rows, fill=TRUE), "tables/eventstudy_all.csv")

# ---- PLACEBOS (CSV) ----
placebo_dates <- as.Date(
  placebo_dates,
  tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y%m%d")
)

if (length(placebo_dates)){
  pb_all <- list()
  for (pd in placebo_dates){
    post_pb <- as.integer(dat$date >= pd)
    D_pb    <- dat$treated * post_pb
    dat_pb  <- data.frame(dat, D_pb=D_pb)
    for (y in outcomes){
      if (has_fixest){
        rhs <- rhs_join(c("D_pb",controls))
        fml <- as.formula(paste(y,"~",rhs,"| country + date"))
        fit <- feols(fml, data=dat_pb, cluster=~country+date)
        tt  <- tidy(fit, conf.int=TRUE)
      } else {
        rhs <- c("D_pb",controls,"factor(country)","factor(date)")
        fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
        fit <- lm(fml, data=dat_pb); vc <- vcovCL(fit, cluster=dat_pb$date)
        tt  <- tidy(fit, conf.int=TRUE, vcov=vc)
      }
      tt$outcome <- y; pd_date <- suppressWarnings(as.Date(pd))
      pd_str  <- if (!is.na(pd_date)) strftime(pd_date, "%Y-%m-%d") else as.character(pd)
      tt$placebo_date <- pd_str
      key <- paste0(y, "_", gsub("-", "", pd_str))  # e.g., "200701"
      pb_all[[key]] <- tt

    }
  }
  fwrite(rbindlist(pb_all, fill=TRUE), "tables/did_placebos.csv")
}

# ---- TREAT CZ (CSV) ----
dat[, treated_cz := as.integer(country %in% c("cz","CZ"))]
dat[, D_cz := treated_cz * post]
cz_rows <- vector("list", length(outcomes)); names(cz_rows) <- outcomes
for (y in outcomes){
  if (has_fixest){
    rhs <- rhs_join(c("D_cz",controls))
    fml <- as.formula(paste(y,"~",rhs,"| country + date"))
    fit <- feols(fml, data=dat, cluster=~country+date); tt <- tidy(fit, conf.int=TRUE)
  } else {
    rhs <- c("D_cz",controls,"factor(country)","factor(date)")
    fml <- as.formula(paste(y,"~",paste(rhs,collapse=" + ")))
    fit <- lm(fml, data=dat); vc <- vcovCL(fit, cluster=dat$date)
    tt  <- tidy(fit, conf.int=TRUE, vcov=vc)
  }
  tt$outcome <- y; cz_rows[[y]] <- tt
}
fwrite(rbindlist(cz_rows, fill=TRUE), "tables/did_treat_cz.csv")

# ---- OLS snapshot (optional; CSV) ----
if (all(c("mthly_hh_expense","mthly_hh_income") %in% names(dat))){
  ss <- na.omit(dat[, .(date, y=mthly_hh_expense, x=mthly_hh_income)])
  last_d <- max(ss$date); cs <- ss[date==last_d]
  ols_data <- if (nrow(cs) >= 20) cs else ss
  ols <- lm(y ~ x, data=ols_data)
  dw <- tryCatch(lmtest::dwtest(ols)$statistic[[1]], error=function(e) NA_real_)
  jb <- tryCatch(tseries::jarque.bera.test(residuals(ols)), error=function(e) NULL)
  sk <- tryCatch(moments::skewness(residuals(ols)), error=function(e) NA_real_)
  ku <- tryCatch(moments::kurtosis(residuals(ols)), error=function(e) NA_real_)
  jb_stat <- if (!is.null(jb)) unname(jb$statistic) else NA_real_
  jb_p    <- if (!is.null(jb)) unname(jb$p.value)   else NA_real_
  fwrite(tidy(ols, conf.int=TRUE), "tables/ols_coeffs.csv")
  fwrite(data.table(metric=c("Durbin-Watson","Jarque-Bera","Prob(JB)","Skew","Kurtosis"),
                    value=c(dw,jb_stat,jb_p,sk,ku)),
         "tables/ols_diag.csv")
}
cat("done. outputs in /tables: did_main.csv, eventstudy_all.csv, ",
    "did_placebos.csv (if any), did_treat_cz.csv, ",
    "ols_coeffs/ols_diag if present\n", sep="")



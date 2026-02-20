#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(lmtest)
  library(sandwich)
})

# ---------- CONFIG ----------
in_path  <- "/workspaces/socka-czsk-euro/data/monthly_panel_clean.csv"

# keep only outcomes that actually exist in your file
outcomes <- c("hicp_yoy","unemp_rate","hicp",
              "imports_world_meur","exports_world_meur",
              "mro","log_imp","log_exp")
controls_candidate <- c("fx_eur","fx_vol","repo","mro")

min_n_obs <- 60L
outdir <- "tables3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- CLUSTERING DEFAULT (PRIMARY = TIME) ----------
VCOV_TIME <- ~ time_id

# ---------- LOAD DATA ----------
dat <- fread(in_path)

# robust time_id if missing
if (!("time_id" %in% names(dat))) {
  if ("date" %in% names(dat)) {
    if (!inherits(dat$date, "Date")) dat[, date := as.IDate(date)]
    setorder(dat, date)
    udates <- sort(unique(dat$date))
    dat[, time_id := match(date, udates)]
  } else stop("No 'time_id' or 'date' column found; one is required.")
}

# ensure D exists: treated*post for CZ w/ counterfactual date; SK untreated
if (!("D" %in% names(dat))) {
  # treatment map (counterfactual adoption date for CZ; SK never adopts)
  treat_map <- data.table(country = c("CZ","SK"),
                          treat_date = as.IDate(c("2009-01-01", NA_character_)))
  dat[, date := as.IDate(date)]
  dat <- merge(dat, treat_map, by = "country", all.x = TRUE, sort = FALSE)
  dat[, treated := as.integer(!is.na(treat_date))]
  dat[, post    := as.integer(!is.na(treat_date) & date >= treat_date)]
  dat[, D       := as.integer(treated == 1L & post == 1L)]
  # event time (optional)
  dat[, rel := ifelse(!is.na(treat_date), as.integer(time_id - time_id[date == treat_date][1]), NA_integer_),
      by = country]
}

# types
dat[, country := as.character(country)]
dat[, time_id := as.integer(time_id)]
dat[, D := as.integer(D)]

# sanity
cat("\n-- treat/post sanity --\n")
print(dat[, .(n=.N,
              post_ones=sum(post, na.rm=TRUE),
              D_ones=sum(D, na.rm=TRUE)), by=country])

# ---------- NA POLICY FOR CONTROLS ----------
na_audit <- function(cols) {
  rbindlist(lapply(cols, function(v) {
    data.table(var=v,
               pNA_all = mean(is.na(dat[[v]])),
               pNA_treatpost = mean(is.na(dat[D==1, get(v)])))
  }))
}
if (length(controls_candidate)) {
  cat("\n-- missingness audit (pre-fill) --\n"); print(na_audit(controls_candidate))
  # fill only inside (treated & post) to avoid dropping those rows
  for (v in intersect(controls_candidate, names(dat))) {
    set(dat, i = which(dat$D==1 & is.na(dat[[v]])), j = v,
        value = if (is.numeric(dat[[v]])) 0 else dat[[v]])
  }
  cat("\n-- missingness audit (post-fill) --\n"); print(na_audit(controls_candidate))
  drop_due_to_na <- controls_candidate[vapply(controls_candidate, function(v) mean(is.na(dat[[v]])) > 0.25, TRUE)]
  if (length(drop_due_to_na)) {
    message("\n-- controls dropped due to NA after fill: ", paste(drop_due_to_na, collapse=", "))
    controls_candidate <- setdiff(controls_candidate, drop_due_to_na)
  } else {
    message("\n-- controls dropped due to NA after fill: (none)")
  }
  message("-- controls kept after NA policy: ", ifelse(length(controls_candidate)==0,"(none)", paste(controls_candidate, collapse=", ")))
}

# log heads
kept_outcomes <- outcomes[outcomes %in% names(dat)]
cat("outcomes used:", paste(kept_outcomes, collapse=", "), "\n")
cat("controls used (candidate):", ifelse(length(controls_candidate)==0,"(none)", paste(controls_candidate, collapse=", ")), "\n\n")

# --- add this helper above select_controls() ---
balance_two_country <- function(dt, y, must_have = c("country","time_id","D")) {
  stopifnot(all(must_have %in% names(dt)))
  # keep rows where y isn’t NA
  dt <- dt[!is.na(get(y))]
  # for each (time_id), keep only months where both CZ & SK are present
  good_t <- dt[, .N, by = .(time_id)][N >= 2L, time_id]
  dt <- dt[time_id %in% good_t]
  # also ensure within CZ we have both pre & post
  cz_has <- dt[country == "CZ", .(pre = any(D == 0L), post = any(D == 1L))]
  if (nrow(cz_has) == 0L || !cz_has$pre || !cz_has$post) {
    return(dt[0])  # empty => will trigger skip
  }
  # ensure SK actually exists in the trimmed sample
  if (!("SK" %in% dt$country)) return(dt[0])
  dt
}

# --- replace your select_controls() with this ---
select_controls <- function(dat, y, candidates, min_n = 100L) {
  # order candidates by fewest NAs globally
  cand_ord <- if (length(candidates)) {
    na_counts <- sapply(candidates, function(v) sum(is.na(dat[[v]])))
    names(sort(na_counts, decreasing = FALSE))
  } else character(0)

  base_cols <- c("country","time_id","D", y)
  # start from balanced months (both CZ & SK present for outcome y)
  df_base <- balance_two_country(dat[, ..base_cols], y)
  df_base <- na.omit(df_base)  # still remove any stragglers
  if (nrow(df_base) < min_n ||
      data.table::uniqueN(df_base$D) < 2L ||
      df_base[, .N, by=time_id][, min(N)] < 2L) {
    return(list(kept = character(0), dfm = df_base))
  }

  kept <- character(0); best_dfm <- df_base

  for (cx in cand_ord) {
    try_cols <- c(base_cols, kept, cx)
    df_try  <- dat[, ..try_cols]
    # enforce balance again after adding a control (because NA in control can unbalance)
    df_try  <- balance_two_country(df_try, y)
    df_try  <- na.omit(df_try)

    ok <- (nrow(df_try) >= min_n) &&
          (data.table::uniqueN(df_try$D) >= 2L) &&
          (df_try[, .N, by=time_id][, min(N)] >= 2L) &&
          all(c("CZ","SK") %in% df_try$country)

    if (ok) { kept <- c(kept, cx); best_dfm <- df_try }
  }

  list(kept = kept, dfm = best_dfm)
}


print_compact <- function(m, dfm, label_outcome) {
  cat("\n========================\nOutcome:", label_outcome, "\n")

  show_row <- function(label, ct) {
    cat(sprintf("\n-- Coefs (%s) [D + controls] --\n", label))
    if (!is.null(ct) && "D" %in% rownames(ct)) {
      out <- data.frame(Estimate = ct["D",1],
                        `Std. Error` = ct["D",2],
                        `t value` = ct["D",3],
                        `Pr(>|t|)` = ct["D",4])
      print(out, row.names = FALSE)
    } else cat("D not estimated.\n")
  }

  if (inherits(m,"fixest")) {
    # PRIMARY: cluster by time_id
    ct_time <- tryCatch(coeftable(summary(m, cluster = VCOV_TIME)), error=function(e) NULL)
    show_row("cluster: time_id", ct_time)

    # Secondary views
    ct_iid  <- tryCatch(coeftable(summary(m, vcov="iid")), error=function(e) NULL)
    show_row("IID", ct_iid)
    ct_het  <- tryCatch(coeftable(summary(m, vcov="hetero")), error=function(e) NULL)
    show_row("HC", ct_het)
    ct_tw   <- tryCatch(coeftable(summary(m, cluster = ~ country + time_id)), error=function(e) NULL)
    show_row("cluster: country + time_id", ct_tw)

    cat("\n-- Fit stats --\n")
    print(try(fixest::fitstat(m, show_types = TRUE), silent = TRUE))
  } else {
    # OLS fallback paths (rare here)
    ct_time <- tryCatch(lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = dfm$time_id)), error=function(e) NULL)
    show_row("cluster: time_id", ct_time)
  }
}

fit_one <- function(dfm, y, kept_ctrls) {
  rhs <- c("D", kept_ctrls)
  rhs_txt <- if (length(rhs)) paste(rhs, collapse = " + ") else "D"
  fml_fixest <- as.formula(sprintf("%s ~ %s | country + time_id", y, rhs_txt))
  m <- tryCatch(feols(fml_fixest, data = dfm), error=function(e) NULL)
  if (!is.null(m)) return(list(model=m, rhs_txt=rhs_txt))
  fml_lm <- as.formula(sprintf("%s ~ %s + factor(country) + factor(time_id)", y, rhs_txt))
  list(model = tryCatch(stats::lm(fml_lm, data=dfm), error=function(e) NULL),
       rhs_txt = rhs_txt)
}

get_beta_time <- function(m, dfm) {
  if (inherits(m,"fixest")) {
    sm <- tryCatch(summary(m, cluster = VCOV_TIME), error=function(e) NULL)
    if (is.null(sm)) return(NA_real_)
    ct <- tryCatch(coeftable(sm), error=function(e) NULL)
    if (is.null(ct) || !("D" %in% rownames(ct))) return(NA_real_)
    return(unname(ct["D","Estimate"]))
  } else {
    vct <- sandwich::vcovCL(m, cluster = dfm$time_id)
    ct  <- tryCatch(lmtest::coeftest(m, vcov. = vct), error=function(e) NULL)
    if (is.null(ct)) return(NA_real_)
    rn <- rownames(coef(summary(m))); rownames(ct) <- rn[seq_len(nrow(ct))]
    if (!("D" %in% rownames(ct))) return(NA_real_)
    return(unname(ct["D",1]))
  }
}

get_beta <- function(m, dfm, vcov_arg = NULL) {
  if (inherits(m,"fixest")) {
    sm <- tryCatch(summary(m, vcov = vcov_arg), error=function(e) NULL)
    if (is.null(sm)) return(NA_real_)
    ct <- tryCatch(coeftable(sm), error=function(e) NULL)
    if (is.null(ct) || !("D" %in% rownames(ct))) return(NA_real_)
    return(unname(ct["D","Estimate"]))
  } else {
    vv <- switch(as.character(vcov_arg),
                 iid = sandwich::vcovHC(m, type="const"),
                 hetero = sandwich::vcovHC(m, type="HC1"),
                 twoway = sandwich::vcovCL(m, cluster = interaction(dfm$country, dfm$time_id)),
                 sandwich::vcovHC(m, type="const"))
    ct <- tryCatch(lmtest::coeftest(m, vcov.=vv), error=function(e) NULL)
    if (is.null(ct)) return(NA_real_)
    rn <- rownames(coef(summary(m))); rownames(ct) <- rn[seq_len(nrow(ct))]
    if (!("D" %in% rownames(ct))) return(NA_real_)
    return(unname(ct["D",1]))
  }
}

# ---------- MAIN LOOP ----------
main_rows <- list(); gof_rows <- list(); diag_rows <- list()

cat("\n")
for (y in kept_outcomes) {

  sel <- select_controls(dat, y, controls_candidate, min_n = min_n_obs)
  kept <- sel$kept; dfm <- sel$dfm

  cat(sprintf("\n[diagnostics: %s ] kept:  %s\n", y,
              ifelse(length(kept)==0,"(none)", paste(kept, collapse=", "))))
  cat(sprintf("[diagnostics: %s ] dropped: %s\n", y,
              ifelse(length(controls_candidate)==0,"(none)",
                     ifelse(length(setdiff(controls_candidate, kept))==0,"(none)",
                            paste(setdiff(controls_candidate, kept), collapse=", ")))))

  if (nrow(dfm) < min_n_obs || data.table::uniqueN(dfm$country) < 2L || data.table::uniqueN(dfm$time_id) < 2L) {
    cat(sprintf("\n[skip:%s] Too few complete observations after control selection.\n", y)); next
  }
  if (data.table::uniqueN(dfm$D) < 2L) {
    cat(sprintf("\n[skip:%s] D has no variation in the complete sample.\n", y)); next
  }

  fit <- fit_one(dfm, y, kept); m <- fit$model
  if (is.null(m)) { cat(sprintf("\n[skip:%s] Estimation failed.\n", y)); next }

  print_compact(m, dfm, y)

  # collect estimates
  est_time <- get_beta_time(m, dfm)
  est_iid  <- get_beta(m, dfm, "iid")
  est_hc   <- get_beta(m, dfm, "hetero")
  est_cl2  <- if (inherits(m,"fixest")) get_beta(m, dfm, "twoway") else get_beta(m, dfm, "twoway")

  main_rows[[length(main_rows)+1]] <- data.table(
    outcome = y,
    kept_controls = ifelse(length(kept)==0,"(none)", paste(kept, collapse=",")),
    beta_D_time = est_time,   # PRIMARY
    beta_D_iid  = est_iid,
    beta_D_hc   = est_hc,
    beta_D_cl2  = est_cl2
  )

  if (inherits(m,"fixest")) {
    gof_rows[[length(gof_rows)+1]] <- data.table(
      outcome = y,
      n = stats::nobs(m),
      aic = tryCatch(AIC(m), error=function(e) NA_real_),
      bic = tryCatch(BIC(m), error=function(e) NA_real_),
      r2  = tryCatch(fixest::fitstat(m)["r2"], error=function(e) NA_real_),
      ar2 = tryCatch(fixest::fitstat(m)["ar2"], error=function(e) NA_real_)
    )
  }

  diag_rows[[length(diag_rows)+1]] <- data.table(
    outcome = y,
    kept = ifelse(length(kept)==0,"(none)", paste(kept, collapse=",")),
    dropped = ifelse(length(setdiff(controls_candidate, kept))==0,"(none)",
                     paste(setdiff(controls_candidate, kept), collapse=","))
  )
}

# ---------- WRITE OUTPUTS ----------
if (length(main_rows)) fwrite(rbindlist(main_rows), file.path(outdir, "did_main.csv"))
if (length(gof_rows))  fwrite(rbindlist(gof_rows),  file.path(outdir, "did_main_gof.csv"))
if (length(diag_rows)) fwrite(rbindlist(diag_rows), file.path(outdir, "kept_dropped_controls.csv"))
fwrite(head(dat, 10), file.path(outdir, "data_head.csv"))

cat("\nDone. Outputs in", outdir, ":\n- data_head.csv\n- did_main.csv, did_main_gof.csv\n- kept_dropped_controls.csv\n")

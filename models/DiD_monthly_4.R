#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(lmtest)
  library(sandwich)
  library(ggplot2)
})

# ----------------- CONFIG -----------------
in_path  <- "/workspaces/socka-czsk-euro/data/monthly_panel_clean.csv"
outdir   <- "tables4"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outcomes <- c("hicp_yoy","unemp_rate","hicp",
              "imports_world_meur","exports_world_meur",
              "mro","log_imp","log_exp")

# keep the list broad; selection will drop NA-heavy ones automatically
controls_candidate <- c("fx_eur","fx_vol","repo","mro")
min_n_obs <- 60L

# ----------------- LOAD + NORMALIZE -----------------
dat <- fread(in_path)

# lower-case names
setnames(dat, names(dat), tolower(names(dat)))

# try to normalize common alt names
if (!"country" %in% names(dat)) {
  for (alt in c("geo","cn","iso2","country_code")) {
    if (alt %in% names(dat)) { setnames(dat, alt, "country"); break }
  }
}
if (!"date" %in% names(dat)) {
  for (alt in c("period","month","time","date_month")) {
    if (alt %in% names(dat)) { setnames(dat, alt, "date"); break }
  }
}

# types
dat[, country := as.character(country)]
if (!inherits(dat$date, "Date")) suppressWarnings(dat[, date := as.IDate(date)])

if (!"time_id" %in% names(dat)) {
  setorder(dat, date)
  udates <- sort(unique(dat$date))
  dat[, time_id := match(date, udates)]
} else {
  dat[, time_id := as.integer(time_id)]
}

cat("\n-- columns present --\n"); print(names(dat))
if (!all(c("country","date","time_id") %in% names(dat))) {
  stop("Required columns missing after normalization. Have: ",
       paste(names(dat), collapse=", "))
}

# ----------------- TREATMENT CODING -----------------
# CZ 'treated' (counterfactual adoption effect) from 2009-01-01; SK untreated
treat_date_map <- data.table(
  country = c("CZ","SK"),
  treat_date = as.IDate(c("2009-01-01", NA))
)

dat <- treat_date_map[dat, on="country"]
dat[, treated := as.integer(!is.na(treat_date))]
dat[, post    := as.integer(!is.na(treat_date) & date >= treat_date)]
dat[, D       := as.integer(treated==1 & post==1)]
dat[, rel     := ifelse(treated==1, as.integer(time_id - time_id[date==treat_date][1]), NA_integer_),
    by = country]

# sanity
san <- dat[, .(n=.N, post_ones = sum(post, na.rm=TRUE), D_ones=sum(D, na.rm=TRUE)), by=country]
cat("\n-- treat/post sanity --\n"); print(san)

# ----------------- MISSINGNESS POLICY FOR CONTROLS -----------------
miss_audit <- function(d, vars) {
  rbindlist(lapply(vars, \(v) {
    data.table(
      var = v,
      pNA_all        = mean(is.na(d[[v]])),
      pNA_treatpost  = mean(is.na(d[treated==1 & post==1][[v]])),
      pNA_treatpre   = mean(is.na(d[treated==1 & post==0][[v]]))  # NEW: this catches fx_eur problem
    )
  }))
}

aud_pre  <- miss_audit(dat, controls_candidate)
cat("\n-- missingness audit (pre-fill) --\n"); print(aud_pre)

# (no fills here)

aud_post <- miss_audit(dat, controls_candidate)
cat("\n-- missingness audit (post-fill) --\n"); print(aud_post)

# drop controls with heavy NA anywhere OR any NA in treated-pre or treated-post
bad <- aud_post[pNA_all > .40 | pNA_treatpost > .05 | pNA_treatpre > .05, var]
kept_controls <- setdiff(controls_candidate, bad)

cat("\n-- controls dropped due to NA after filter: ",
    ifelse(length(bad)==0,"(none)", paste(bad, collapse=", ")), "\n")
cat("-- controls kept after NA policy: ",
    ifelse(length(kept_controls)==0,"(none)", paste(kept_controls, collapse=", ")), "\n")


# ----------------- HELPERS -----------------
`%+%` <- function(a,b) c(a,b)

select_controls <- function(dat, y, candidates, min_n = 100L) {
  if (length(candidates)) {
    na_counts <- sapply(candidates, \(v) sum(is.na(dat[[v]])))
    cand_ord  <- names(sort(na_counts, decreasing = FALSE))
  } else cand_ord <- character(0)

  base_cols <- c("country","time_id","D", y)
  df_base <- na.omit(dat[, ..base_cols])
  if (nrow(df_base) < min_n || uniqueN(df_base$D) < 2L) {
    return(list(kept = character(0), dfm = df_base))
  }
  kept <- character(0); best_dfm <- df_base
  for (cx in cand_ord) {
    try_cols <- c(base_cols, kept, cx)
    df_try  <- na.omit(dat[, ..try_cols])
    if (nrow(df_try) >= min_n && uniqueN(df_try$D) >= 2L) {
      kept <- c(kept, cx); best_dfm <- df_try
    }
  }
  list(kept = kept, dfm = best_dfm)
}

print_compact <- function(m, dfm, label_outcome) {
  cat("\n========================\nOutcome:", label_outcome, "\n")

  show_row <- function(label, vcov = NULL) {
    cat(sprintf("\n-- Coefs (%s) [D + controls] --\n", label))
    if (inherits(m,"fixest")) {
      sm <- tryCatch(summary(m, vcov = vcov), error=function(e) NULL)
      if (is.null(sm) || !("D" %in% names(coef(sm)))) {
        cat("D not estimated.\n"); return(invisible())
      }
      ct <- as.data.frame(coeftable(sm))
      if ("D" %in% rownames(ct)) print(ct["D", , drop=FALSE]) else cat("D not estimated.\n")
    } else {
      vv <- if (is.null(vcov)) sandwich::vcovHC(m, type = "const") else vcov
      ct <- tryCatch(lmtest::coeftest(m, vcov.=vv), error=function(e) NULL)
      if (is.null(ct)) { cat("coeftest failed.\n"); return(invisible()) }
      rn <- rownames(coef(summary(m))); rownames(ct) <- rn[seq_len(nrow(ct))]
      if ("D" %in% rownames(ct)) print(as.data.frame(ct["D", , drop=FALSE])) else cat("D not estimated.\n")
    }
  }

  # time clustering
  show_row("cluster: time_id", ~ time_id)
  # iid + HC
  show_row("IID", "iid")
  show_row("HC", "hetero")
  # two-way clustering
  cat("\n-- Coefs (cluster: country + time_id) [D + controls] --\n")
  if (inherits(m,"fixest")) {
    sm2 <- tryCatch(summary(m, cluster = ~ country + time_id), error=function(e) NULL)
    if (!is.null(sm2) && "D" %in% names(coef(sm2))) {
      print(as.data.frame(coeftable(sm2)["D", , drop=FALSE]))
    } else cat("D not estimated.\n")
  } else {
    vc2 <- sandwich::vcovCL(m, cluster = interaction(dfm$country, dfm$time_id))
    ct  <- tryCatch(lmtest::coeftest(m, vcov.=vc2), error=function(e) NULL)
    if (!is.null(ct) && "D" %in% rownames(ct)) {
      print(as.data.frame(ct["D", , drop=FALSE]))
    } else cat("D not estimated.\n")
  }

  # fit stats
  cat("\n-- Fit stats --\n")
  if (inherits(m,"fixest")) {
    avail <- tryCatch(fixest::fitstat(m, show_types = TRUE), error=function(e) NULL)
    if (is.null(avail)) { print(avail); return(invisible()) }
    want <- c("n","ll","aic","bic","rmse","r2","ar2","wr2")
    use  <- intersect(want, rownames(avail))
    print(fixest::fitstat(m, as.formula(paste("~", paste(use, collapse=" + ")))))
  } else {
    sm <- summary(m)
    print(list(n = nobs(m), r2 = sm$r.squared, r2_adj = sm$adj.r.squared,
               aic = tryCatch(AIC(m), error=function(e) NA_real_),
               bic = tryCatch(BIC(m), error=function(e) NA_real_)))
  }
}

fit_one <- function(dfm, y, kept_ctrls) {
  stopifnot(all(c("country","time_id") %in% names(dfm)))
  rhs <- c("D", kept_ctrls); rhs_txt <- paste(rhs, collapse = " + ")
  fml <- as.formula(sprintf("%s ~ %s | country + time_id", y, rhs_txt))
  m <- tryCatch(feols(fml, data = dfm), error=function(e) NULL)
  if (!is.null(m)) return(list(model=m, rhs_txt=rhs_txt))
  fml_lm <- as.formula(sprintf("%s ~ %s + factor(country) + factor(time_id)", y, rhs_txt))
  m2 <- tryCatch(stats::lm(fml_lm, data = dfm), error=function(e) NULL)
  list(model=m2, rhs_txt=rhs_txt)
}

# ----------------- MAIN LOOP -----------------
main_rows <- list()
gof_rows  <- list()

cat("outcomes used:", paste(outcomes, collapse=", "), "\n")
cat("controls used (candidate):", paste(controls_candidate, collapse=", "), "\n\n")

for (y in outcomes) {
  if (!(y %in% names(dat))) {
    cat(sprintf("[skip:%s] outcome not found.\n\n", y)); next
  }
  sel  <- select_controls(dat, y, kept_controls, min_n = min_n_obs)
  # ensure treated-pre exists in the model frame; if not, iteratively drop NA-heavy controls
if (nrow(sel$dfm[treated==1 & post==0]) == 0L && length(sel$kept) > 0L) {
  # re-rank kept controls by NA share within treated-pre, drop worst offenders until CZ-pre reappears
  na_rank <- sort(sapply(sel$kept, \(v) mean(is.na(dat[treated==1 & post==0][[v]]))), decreasing = TRUE)
  kept_fix <- sel$kept
  dfm_fix  <- sel$dfm
  for (cx in names(na_rank)) {
    try_keep <- setdiff(kept_fix, cx)
    try_cols <- c("country","time_id","D", y, try_keep)
    df_try <- na.omit(dat[, ..try_cols])
    if (nrow(df_try[treated==1 & post==0]) > 0L && uniqueN(df_try$D) > 1L) {
      kept_fix <- try_keep; dfm_fix <- df_try
    }
  }
  sel$kept <- kept_fix; sel$dfm <- dfm_fix
}

  kept <- sel$kept; dfm <- sel$dfm

  cat(sprintf("\n[diagnostics: %s ] kept: %s\n", y, ifelse(length(kept)==0," (none)", paste(kept, collapse=", "))))
  cat(sprintf("[diagnostics: %s ] dropped: %s\n", y,
              ifelse(length(setdiff(kept_controls, kept))==0,"(none)", paste(setdiff(kept_controls, kept), collapse=", "))))

  # guards
  if (nrow(dfm) < min_n_obs || uniqueN(dfm$country) < 2L || uniqueN(dfm$time_id) < 2L) {
    cat(sprintf("[skip:%s] Too few complete observations.\n\n", y)); next
  }
  if (uniqueN(dfm$D) < 2L) {
    cat(sprintf("[skip:%s] D has no variation after filtering.\n\n", y)); next
  }

  fit <- fit_one(dfm, y, kept)
  m <- fit$model
  if (is.null(m)) { cat(sprintf("[skip:%s] estimation failed.\n\n", y)); next }

  # report
  print_compact(m, dfm, y)

  # collect rows
  coef_safe <- function(model, vc) {
    if (inherits(model,"fixest")) {
      sm <- tryCatch(summary(model, vcov=vc), error=function(e) NULL)
      if (is.null(sm) || !("D" %in% names(coef(sm)))) return(NA_real_)
      as.numeric(coeftable(sm)["D","Estimate"])
    } else NA_real_
  }
  main_rows[[length(main_rows)+1]] <- data.table(
    outcome = y,
    kept_controls = ifelse(length(kept)==0,"(none)", paste(kept, collapse=",")),
    beta_D_iid = coef_safe(m, "iid"),
    beta_D_hc  = coef_safe(m, "hetero"),
    beta_D_clt = coef_safe(m, ~ time_id),
    beta_D_cl2 = {
      sm2 <- tryCatch(summary(m, cluster = ~ country + time_id), error=function(e) NULL)
      if (is.null(sm2) || !("D" %in% names(coef(sm2)))) NA_real_ else as.numeric(coeftable(sm2)["D","Estimate"])
    }
  )

  if (inherits(m,"fixest")) {
    gof_rows[[length(gof_rows)+1]] <- data.table(
      outcome = y,
      n = nobs(m),
      aic = tryCatch(AIC(m), error=function(e) NA_real_),
      bic = tryCatch(BIC(m), error=function(e) NA_real_),
      r2  = tryCatch(as.numeric(fitstat(m)["r2"]), error=function(e) NA_real_),
      ar2 = tryCatch(as.numeric(fitstat(m)["ar2"]), error=function(e) NA_real_)
    )
  }
}

# ----------------- WRITE -----------------
if (length(main_rows)) fwrite(rbindlist(main_rows), file.path(outdir,"did_main.csv"))
if (length(gof_rows))  fwrite(rbindlist(gof_rows),  file.path(outdir,"did_main_gof.csv"))
fwrite(head(dat,10), file.path(outdir,"data_head.csv"))

cat("\nDone. Outputs in", outdir, ":\n- data_head.csv\n- did_main.csv, did_main_gof.csv\n")

# ----------------- OPTIONAL: quick ES plot (safe) -----------------
# writes only if hicp_yoy present and rel defined for treated country
try({
  if ("hicp_yoy" %in% names(dat)) {
    es_df <- dat[treated==1 & !is.na(rel) & !is.na(hicp_yoy),
                 .(y = mean(hicp_yoy, na.rm=TRUE)), by=rel][order(rel)]
    gg <- ggplot(es_df, aes(x=rel, y=y)) +
      geom_hline(yintercept = es_df[rel<0, mean(y, na.rm=TRUE)], linetype="dashed") +
      geom_vline(xintercept = 0, linetype="dotted") +
      geom_line() + geom_point() +
      labs(x="months relative to 2009-01 (CZ)", y="avg HICP yoy (treated)", title="Event-style mean path (treated)")
    ggsave(filename = file.path(outdir,"es_hicp_yoy.png"), plot=gg, width=7, height=4, dpi=150)
  }
}, silent = TRUE)

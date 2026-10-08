suppressPackageStartupMessages(library(jsonlite))
root <- normalizePath(getwd())
stopifnot(file.exists(file.path(root, "R/closure_rcpp.cpp")))
dir.create("results/R", recursive = TRUE, showWarnings = FALSE)
dir.create("build/Rcpp", recursive = TRUE, showWarnings = FALSE)
native_sources <- c("R/closure_rcpp.cpp", "analysis/closure_core.cpp",
                    "analysis/closure_core.h", "analysis/closure_numerics.h", "analysis/closure_impl.hpp")
native_cache <- file.path("build/Rcpp", digest::digest(vapply(native_sources,
  function(p) digest::digest(file=p,algo="sha256"),character(1)),algo="xxhash64"))
local({
  target <- parent.env(environment())
  old <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_CXXFLAGS"), unset = NA_character_)
  on.exit({
    for (key in names(old)) {
      if (is.na(old[[key]])) Sys.unsetenv(key)
      else do.call(Sys.setenv, setNames(list(old[[key]]), key))
    }
  })
  bh <- system.file("include", package = "BH")
  if (!nzchar(bh)) stop("Install the BH package before loading the simulation functions.")
  Sys.setenv(PKG_CPPFLAGS = paste(Sys.getenv("PKG_CPPFLAGS"), paste0("-I", shQuote(bh))),
             PKG_CXXFLAGS = paste(Sys.getenv("PKG_CXXFLAGS"), "-std=c++17"))
  Rcpp::sourceCpp("R/closure_rcpp.cpp", cacheDir = native_cache,
                  rebuild = FALSE, env = target)
})
closure_rcpp <- function(x, mandatory, alpha = .05, mode = 0L)
  closure_native(x, mandatory, alpha, mode)
e_cap_rcpp <- function(m, alpha = .05) e_cap_native(m, alpha)
by_cutoffs_rcpp <- function(m, alpha = .05) by_cutoffs_native(m, alpha)
by_bins_rcpp <- function(p, alpha = .05) by_bins_native(p, alpha)
reject_e_rcpp <- function(e, alpha = .05) reject_e_native(e, alpha)
step_up_rcpp <- function(x, thresholds, decreasing = TRUE)
  step_up_native(x, thresholds, decreasing)

as_matrix <- function(x) if (is.matrix(x)) x else matrix(x, nrow = 1L)
log_add <- function(x, y) {
  a <- pmax(x, y)
  ans <- a + log(exp(x - a) + exp(y - a))
  ans[is.infinite(a) & a < 0] <- -Inf
  ans
}
reject_p <- function(lp, alpha = .05, by = FALSE, log.p = TRUE) {
  lp <- as_matrix(lp); m <- ncol(lp)
  cuts <- if (by) by_cutoffs_rcpp(m, alpha) else alpha * seq_len(m) / m
  if (by && log.p) return(step_up_rcpp(exp(lp), cuts, FALSE))
  step_up_rcpp(lp, if (log.p) log(cuts) else cuts, FALSE)
}
reject_loge <- function(le, alpha = .05) {
  reject_e(bounded_e(le, alpha), alpha)
}
reject_e <- function(e, alpha = .05) {
  e <- as_matrix(e); bins <- .calibrated_bins(e, alpha)
  if (!is.null(bins)) {
    scores <- bins; scores[bins == 0L] <- Inf
    return(step_up_rcpp(scores, seq_len(ncol(e)), FALSE))
  }
  reject_e_rcpp(e, alpha)
}
bounded_e <- function(le, alpha = .05) {
  le <- as_matrix(le); cap <- e_cap_rcpp(ncol(le), alpha)
  ans <- exp(pmin(le, log(cap)))
  ans[le >= log(cap)] <- cap
  ans
}
log_abf <- function(x, variance, prior_sd = .2) {
  r <- prior_sd^2 / variance
  -.5 * log1p(r) + .5 * x^2 / variance * r / (1 + r)
}
log_adjuster <- function(y) {
  ans <- y * 0 - log(2)
  small <- y > 0 & y < 1e-3; regular <- y >= 1e-3
  ans[small] <- log(.5 + y[small]/6 + y[small]^2/24 + y[small]^3/120)
  z <- y[regular]
  ans[regular] <- z + log1p(-(1 + z) * exp(-z)) - 2 * log(z)
  ans
}
by_calibrator <- function(p, alpha = .05) {
  p <- as_matrix(p); m <- ncol(p); bins <- by_bins_rcpp(p, alpha)
  q <- matrix(m / (alpha * pmax(1, bins)), nrow(p), m); q[bins == 0L] <- 0

  attr(q, "BY_calibration") <- list(alpha = alpha, bins = bins)
  q
}
.calibrated_bins <- function(e, alpha) {
  info <- attr(e, "BY_calibration")
  if (is.null(info) || !identical(info$alpha, alpha) ||
      !identical(dim(e), dim(info$bins))) return(NULL)
  ref <- ncol(e) / (alpha * pmax(1, info$bins)); ref[info$bins == 0L] <- 0
  if (!identical(as.vector(e), as.vector(ref))) return(NULL)
  info$bins
}
closed <- function(values, mandatory = NULL, alpha = .05, mode = "e") {
  values <- as_matrix(values)
  if (is.null(mandatory)) mandatory <- matrix(FALSE, nrow(values), ncol(values))
  bins <- if (mode == "e") .calibrated_bins(values, alpha) else NULL
  closure_rcpp(if (is.null(bins)) values else bins, as_matrix(mandatory), alpha,
               if (!is.null(bins)) 2L else as.integer(mode == "by"))
}
persistent_e <- function(e, previous = NULL, alpha = .05) {
  mandatory <- reject_e(e, alpha)
  if (!is.null(previous)) mandatory <- mandatory | previous
  ans <- closed(e, mandatory, alpha)
  stopifnot(all(!mandatory | ans))
  ans
}

mean_consistent_R <- function(e, set, alpha = .05) {
  k <- sum(set); if (!k) return(TRUE)
  bins <- .calibrated_bins(as_matrix(e), alpha)
  if (!is.null(bins)) {
    m <- length(e); strength <- ifelse(bins > 0, m/pmax(1,bins), 0)
    bs <- bins[set][order(strength[set])]; bu <- bins[!set][order(strength[!set])]
    cs <- c(0,cumsum(sort(strength[set]))); cu <- c(0,cumsum(sort(strength[!set])))
    for (a in seq_len(k)) for (b in 0:sum(!set)) {
      lhs <- k*(cs[a+1]+cu[b+1]); rhs <- a*(a+b)
      nu <- (a+b+19)*.Machine$double.eps/2; err <- nu/(1-nu)*max(1,abs(lhs)+rhs)
      if (lhs-rhs < -err) return(FALSE)
      if (abs(lhs-rhs) <= err) {
        z <- c(bs[seq_len(a)],bu[seq_len(b)]); z <- z[z>0]
        total <- if(length(z))m*sum(gmp::as.bigq(1,z)) else gmp::as.bigq(0)
        if(k*total < rhs)return(FALSE)
      }
    }
    return(TRUE)
  }
  s <- sort(e[set]); u <- sort(e[!set])
  cs <- c(0, cumsum(s)); cu <- c(0, cumsum(u))
  for (a in seq_len(k)) for (b in 0:sum(!set)) {
    lhs <- alpha * k * (cs[a+1] + cu[b+1]); rhs <- a * (a+b)
    nu <- (a+b+19) * .Machine$double.eps/2
    err <- nu/(1-nu) * max(1,abs(lhs)+rhs)
    if (lhs-rhs < -err) return(FALSE)
    if (abs(lhs-rhs) <= err &&
        gmp::as.bigq(alpha)*k*sum(gmp::as.bigq(c(s[seq_len(a)],u[seq_len(b)]))) < rhs)
      return(FALSE)
  }
  TRUE
}
estimate <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  mu <- mean(x); se <- sd(x) / sqrt(length(x))
  c(mean = mu, mcse = se, lower = mu - 1.96*se, upper = mu + 1.96*se)
}
method_labels <- c(
  bh = "BH snapshot", by = "BY snapshot", ebh = "e-BH snapshot", bh_union = "BH union",
  by_spend = "BY spending", bonferroni = "Bonferroni spending", anytime_by = "Anytime-p BY",
  gaussian = "Gaussian mixture", p_mixture = "P-calibrator mixture", freeze = "Freeze ledger",
  carefree = "Carefree A1", closed_gaussian = "Retained Gaussian closure",
  closed_p = "Retained calibrated-p closure", closed_by = "Retained closed BY",
  closed_carefree = "Retained carefree closure")
write_csv <- function(x, file) write.csv(x, file, row.names = FALSE, na = "")

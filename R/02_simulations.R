source("R/01_methods.R")

simulation_configs <- list(
  list(scenario = "Nested all null", n_alt = 0L), list(scenario = "Nested sparse"),
  list(scenario = "Nested correlated", rho = .4), list(scenario = "Refreshed sparse", nested = FALSE),
  list(scenario = "Inflated uncorrected", true_lambda = 1.5),
  list(scenario = "Inflated corrected", true_lambda = 1.5, used_lambda = 1.5),
  list(scenario = "Weak sparse", effect = .15), list(scenario = "Very sparse", n_alt = 5L),
  list(scenario = "Student t3 exact p all null", n_alt = 0L, nested = FALSE, distribution = "t3"),
  list(scenario = "Mean bias all causal null", n_alt = 0L, bias = .1))

run_simulation <- function(config, reps = 5000L, seed = 20260919L, closure = TRUE) {
  started <- proc.time()[3]
  p <- modifyList(list(scenario = "", n_alt = 20L, effect = .25, rho = 0, nested = TRUE,
    distribution = "normal", true_lambda = 1, used_lambda = 1, bias = 0, prior_sd = .2), config)
  m <- 200L; times <- 12L; alpha <- .05
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260920L)
  order <- t(replicate(reps, sample.int(m)))
  index <- cbind(rep(seq_len(reps), m), as.vector(order))
  truth <- order <= p$n_alt
  set.seed(seed)
  theta <- numeric(m)
  if (p$n_alt) theta[seq_len(p$n_alt)] <- p$effect * rep(c(1, -1), length.out = p$n_alt)
  theta <- matrix(rep(theta, each = reps), reps, m)
  u <- matrix(0, reps, m); empty <- matrix(FALSE, reps, m)
  bh_union <- ebh_union <- by_spend <- bonferroni <- frozen <- empty
  frozen_values <- u; minp <- u; logmix <- logpmix <- u - Inf; maxe <- u
  prev <- list(); maxfdp <- list(); power_auc <- list(); trajectories <- list()
  metric_arrays <- list()
  for (tt in seq_len(times)) {
    n <- 25 * tt
    noise <- if (p$distribution == "t3") matrix(rt(reps*m, 3), reps, m)/sqrt(3) else
      sqrt(p$rho) * rnorm(reps) + sqrt(1-p$rho) * matrix(rnorm(reps*m), reps, m)
    if (p$nested) {
      u <- u + 25*(theta + p$bias) + sqrt(25*p$true_lambda)*noise; x <- u/n
    } else x <- theta + p$bias + sqrt(p$true_lambda/n)*noise
    z <- abs(x)*sqrt(n/p$used_lambda)
    lp <- log(2) + if (p$distribution == "t3") pt(-z*sqrt(3), df = 3, log.p = TRUE) else pnorm(-z, log.p = TRUE)
    le <- log_abf(x, p$used_lambda/n, p$prior_sd)
    lp <- matrix(lp[index], reps, m); le <- matrix(le[index], reps, m)
    bh <- reject_p(lp); by <- reject_p(lp, by = TRUE); ebh <- reject_loge(le)
    bh_union <- bh_union | bh
    ebh_union <- ebh_union | ebh
    by_spend <- by_spend | reject_p(lp, alpha/times, by = TRUE)
    bonferroni <- bonferroni | (lp <= log(alpha/(m*times)))
    minp <- pmin(minp, lp + log(times))
    logmix <- log_add(logmix, le-log(times))
    logpmix <- log_add(logpmix, -log(2)-.5*lp-log(times))
    curr <- le; curr[frozen] <- frozen_values[frozen]
    declared <- reject_loge(curr); fresh <- declared & !frozen
    stopifnot(all(!frozen | declared))
    frozen_values[fresh] <- le[fresh]; frozen <- frozen | declared
    maxe <- pmax(maxe, le)
    e <- bounded_e(logmix); q <- by_calibrator(exp(minp)); adj <- bounded_e(log_adjuster(maxe))
    methods <- list(bh = bh, by = by, ebh = ebh, bh_union = bh_union, ebh_union = ebh_union,
      by_spend = by_spend, bonferroni = bonferroni, anytime_by = reject_e(q), gaussian = reject_e(e),
      p_mixture = reject_loge(logpmix), freeze = frozen, carefree = reject_e(adj))
    if (closure) methods <- c(methods, list(
      closed_gaussian = persistent_e(e, prev$closed_gaussian),
      closed_p = persistent_e(q, prev$closed_p), closed_by = closed(exp(minp), prev$closed_by, mode = "by"),
      closed_carefree = persistent_e(adj, prev$closed_carefree)))
    stopifnot(identical(methods$gaussian, reject_loge(logmix)),
      identical(methods$anytime_by, reject_p(minp, by = TRUE)),
      identical(methods$carefree, reject_loge(log_adjuster(maxe))))
    for (name in names(methods)) {
      r <- methods[[name]]
      if (!name %in% c("bh", "by", "ebh") && !is.null(prev[[name]])) stopifnot(all(!prev[[name]] | r))
      fdp <- rowSums(r & !truth)/pmax(1, rowSums(r))
      power <- if (p$n_alt) rowSums(r & truth)/p$n_alt else rep(NA_real_, reps)
      maxfdp[[name]] <- pmax(if (is.null(maxfdp[[name]])) 0 else maxfdp[[name]], fdp)
      power_auc[[name]] <- (if (is.null(power_auc[[name]])) 0 else power_auc[[name]]) + power/times
      metrics <- list(fdr = fdp, power = power, max_fdp = maxfdp[[name]], mean_discoveries = rowSums(r))
      trajectories[[length(trajectories)+1L]] <- do.call(rbind, lapply(names(metrics), function(metric)
        data.frame(scenario = p$scenario, release = tt, method = name, metric = metric,
          as.list(estimate(metrics[[metric]])), row.names = NULL)))
      metric_arrays[[name]] <- c(metrics, list(mean_time_power = power_auc[[name]]))
    }
    if (closure) stopifnot(all(!methods$gaussian | methods$closed_gaussian),
      all(!methods$anytime_by | methods$closed_p), all(!methods$carefree | methods$closed_carefree))
    prev <- methods
  }
  summary <- do.call(rbind, lapply(names(metric_arrays), function(name)
    do.call(rbind, lapply(names(metric_arrays[[name]]), function(metric)
      data.frame(scenario = p$scenario, method = name, metric = metric,
        as.list(estimate(metric_arrays[[name]][[metric]])), row.names = NULL)))))
  paired <- list()
  comparisons <- list(c("gaussian", "anytime_by"), c("freeze", "gaussian"))
  if (closure) comparisons <- c(comparisons, list(c("closed_gaussian", "gaussian"), c("closed_p", "anytime_by"),
    c("closed_by", "anytime_by"), c("closed_gaussian", "anytime_by")))
  for (pair in comparisons) for (metric in c("power", "fdr")) {
    delta <- metric_arrays[[pair[1]]][[metric]] - metric_arrays[[pair[2]]][[metric]]
    paired[[length(paired)+1]] <- data.frame(scenario = p$scenario, first = pair[1], second = pair[2],
      metric = metric, as.list(estimate(delta)), row.names = NULL)
  }
  result <- list(parameters = p, replicates = reps, seed = seed, seconds = unname(proc.time()[3]-started),
    summary = summary, trajectories = do.call(rbind, trajectories), paired = do.call(rbind, paired),
    final_replicate_metrics = metric_arrays)
  cat("Completed R simulation:", p$scenario, "in", round(result$seconds, 1), "s\n")
  result
}

run_simulations <- function(reps = 5000L, workers = 2L) {
  started_at <- format(Sys.time(), tz="UTC", usetz=TRUE)
  message("Full R simulation run started: ", started_at)
  results <- parallel::mclapply(simulation_configs, run_simulation, reps = reps, mc.cores = workers, mc.set.seed = FALSE)
  stopifnot(!any(vapply(results, inherits, logical(1), "try-error")))
  sensitivity_configs <- lapply(seq_len(9), function(i) {
    grid <- expand.grid(prior_sd = c(.1,.2,.4), used_lambda = c(1,1.5,2))
    list(scenario = paste0("lambda=",grid$used_lambda[i],", prior_sd=",grid$prior_sd[i]),
      prior_sd = grid$prior_sd[i], used_lambda = grid$used_lambda[i])
  })
  sensitivity <- parallel::mclapply(sensitivity_configs, run_simulation, reps = reps, closure = TRUE,
    mc.cores = workers, mc.set.seed = FALSE)
  stopifnot(!any(vapply(sensitivity, inherits, logical(1), "try-error")))
  metadata <- list(language = "R", R = R.version.string, alpha = .05, hypotheses = 200,
    releases = 12, replicates = reps, seed = 20260919L, tie_order_seed = 20260920L,
    RNG = RNGkind(), generated_in_R = TRUE, independent_of_Python_draws = TRUE,
    benchmark_configurations=length(results), sensitivity_configurations=length(sensitivity),
    methods_per_configuration=length(unique(results[[1]]$summary$method)),
    started_at=started_at, completed_at=format(Sys.time(), tz="UTC", usetz=TRUE),
    source_sha256=setNames(lapply(c("R/01_methods.R","R/02_simulations.R",native_sources),
      function(p)digest::digest(file=p,algo="sha256")),c("R/01_methods.R","R/02_simulations.R",native_sources)),
    note = "Full rerun with the recorded design and seeds. Every benchmark and sensitivity configuration evaluates all 16 procedures. Rcpp accelerates deterministic closure decisions.")
  saveRDS(list(metadata = metadata, experiments = results, sensitivity = sensitivity), "results/R/simulations.rds")
  write_csv(do.call(rbind, lapply(results, `[[`, "summary")), "results/R/simulation_summary.csv")
  write_csv(do.call(rbind, lapply(results, `[[`, "trajectories")), "results/R/simulation_trajectories.csv")
  write_csv(do.call(rbind, lapply(results, `[[`, "paired")), "results/R/paired_comparisons.csv")
  write_csv(do.call(rbind, lapply(sensitivity, `[[`, "summary")), "results/R/sensitivity_summary.csv")
  write_json(metadata, "results/R/simulation_metadata.json", pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), "results/R/sessionInfo.txt")
  invisible(results)
}
if (sys.nframe() == 0L) run_simulations()

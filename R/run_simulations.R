args <- commandArgs(trailingOnly = TRUE)
study <- if (length(args)) args[1] else 'all'
choices <- c('benchmark', 'direct', 'factorial', 'extension', 'cumulative', 'all')
if (!study %in% choices) stop('Study must be one of: ', paste(choices, collapse = ', '))
workers <- if (length(args) >= 2L) suppressWarnings(as.integer(args[2])) else
  if (.Platform$OS.type == 'unix') 2L else 1L
stopifnot(length(workers) == 1L, !is.na(workers), workers > 0L)
if (.Platform$OS.type != 'unix') workers <- 1L
options(simulation_workers = workers)
if (study %in% c('all','direct')) {
  needed <- c('results/R/real_data_ledger.csv','results/R/endpoint_metadata.csv')
  if (!all(file.exists(needed))) stop('The direct study requires the two prepared simulation-design inputs (real_data_ledger: release, available, protein_id, endpoint, variant, wald_se_first_order; endpoint_metadata: release, endpoint, cases, controls): ', paste(needed, collapse = ', '))
}
if (study %in% c('all','benchmark')) {
  source('R/02_simulations.R')
  run_simulations(workers = workers)
}
if (study %in% c('all','direct')) {
  source('R/11_retention_study.R')
  run_retention_study(workers = workers)
}
if (study %in% c('all','factorial')) source('R/23_retention_boundary_study.R')
if (study %in% c('all','extension')) {
  source('R/38_cumulative_extension.R')
  run_cumulative_extension(workers = workers)
}
if (study %in% c('all','cumulative')) {
  needed <- c('results/R/boundary_study.rds','results/R/retention_study.rds')
  if (!all(file.exists(needed))) stop('Run the direct and factorial studies before the cumulative replay.')
  source('R/36_cumulative_reporting.R')
  run_cumulative_reporting(workers = workers)
}

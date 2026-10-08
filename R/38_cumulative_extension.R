source('R/01_methods.R')
source('R/retention_diagnostics.R')

extension_grid <- function() {
  g<-expand.grid(releases=c(6L,12L),final_snr=c(5,7),signals=c(60L,100L,140L),
    KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  g<-g[,c('signals','final_snr','releases')]
  g$overlap<-0;g$growth<-'Balanced';g$id<-seq_len(nrow(g))
  g
}

simulate_extension <- function(config,reps=1000L,m=200L,prior_sd=3,alpha=.05) {
  T<-config$releases;set.seed(2026092300L+config$id)
  info<-outer(rep(1,m),(1:T)/T,function(a,t)t^a)
  truth<-matrix(FALSE,reps,m);theta<-matrix(0,reps,m)
  for(i in seq_len(reps)) {
    alt<-sample.int(m,config$signals);truth[i,alt]<-TRUE
    theta[i,alt]<-config$final_snr*sample(c(-1,1),config$signals,TRUE)
  }
  noise<-function()sqrt(.2)*matrix(rnorm(reps),reps,m)+sqrt(.8)*matrix(rnorm(reps*m),reps,m)
  signals<-config$signals
  walk<-matrix(0,reps,m);mix<-matrix(-Inf,reps,m)
  names_methods<-c('gaussian','fresh_gaussian','retained_gaussian','union_gaussian')
  prev<-first<-mx<-list()
  for(n in names_methods){prev[[n]]<-matrix(FALSE,reps,m);first[[n]]<-matrix(T+1L,reps,m);mx[[n]]<-numeric(reps)}
  trajectories<-list()
  uncertified<-union_changed<-withdrawal<-difference<-loss_event<-rep(FALSE,reps);maxloss<-numeric(reps)
  for(t in seq_len(T)) {
    q<-info[,t];dq<-q-if(t==1)0 else info[,t-1]
    walk<-walk+sweep(noise(),2,sqrt(dq),'*')
    z<-sqrt(config$overlap)*sweep(walk,2,sqrt(q),'/')+sqrt(1-config$overlap)*noise()
    x<-theta+sweep(z,2,sqrt(q),'/');variance<-matrix(rep(1/q,each=reps),reps,m)
    le<-log_abf(x,variance,prior_sd=prior_sd)
    mix<-log_add(mix,le-log(T));e<-bounded_e(mix)
    ss<-list(gaussian=reject_e(e))
    ss$fresh_gaussian<-closed(e,ss$gaussian)
    ss$retained_gaussian<-persistent_e(e,prev$retained_gaussian)
    ss$union_gaussian<-prev$union_gaussian|ss$fresh_gaussian
    stopifnot(all(!prev$gaussian|ss$gaussian),all(!prev$retained_gaussian|ss$retained_gaussian))
    differs<-rowSums(ss$union_gaussian!=ss$fresh_gaussian)>0
    union_changed<-union_changed|differs
    for(i in which(differs))uncertified[i]<-uncertified[i]|!mean_consistent_R(e[i,],ss$union_gaussian[i,])
    a<-ss$fresh_gaussian;b<-ss$retained_gaussian;delta<-rowSums(a)-rowSums(b)
    withdrawal<-withdrawal|rowSums(prev$fresh_gaussian&!a)>0
    difference<-difference|rowSums(a!=b)>0;loss_event<-loss_event|delta>0;maxloss<-pmax(maxloss,delta)
    for(n in names_methods) {
      r<-ss[[n]];fdp<-rowSums(r&!truth)/pmax(rowSums(r),1);mx[[n]]<-pmax(mx[[n]],fdp)
      first[[n]][r&first[[n]]==T+1L]<-t
      metrics<-list(fdr=fdp,discoveries=rowSums(r),true_discoveries=rowSums(r&truth),power=rowSums(r&truth)/signals)
      for(metric in names(metrics))trajectories[[length(trajectories)+1L]]<-
        data.frame(release=t,method=n,metric=metric,t(estimate(metrics[[metric]])))
    }
    prev<-ss
  }
  stopifnot(identical(first$fresh_gaussian,first$union_gaussian))
  raw<-lapply(names_methods,function(n) {
    r<-prev[[n]];num<-rowSums(r);tp<-rowSums(r&truth)
    list(discoveries=num,true_discoveries=tp,fdr=rowSums(r&!truth)/pmax(num,1),supfdr=mx[[n]],
      power=tp/signals,first_declaration=rowSums(first[[n]]*truth)/signals)
  });names(raw)<-names_methods
  summary<-do.call(rbind,lapply(names_methods,function(n)do.call(rbind,
    lapply(names(raw[[n]]),function(metric)data.frame(method=n,metric=metric,t(estimate(raw[[n]][[metric]])))))))
  differences<-lapply(names(raw$union_gaussian),function(metric)raw$union_gaussian[[metric]]-raw$retained_gaussian[[metric]])
  names(differences)<-names(raw$union_gaussian)
  paired<-do.call(rbind,lapply(names(differences),function(metric)data.frame(metric=metric,t(estimate(differences[[metric]])))))
  list(config=config,summary=summary,paired=paired,trajectories=do.call(rbind,trajectories),raw=raw,
    uncertified=uncertified,union_changed=union_changed,withdrawal=withdrawal,difference=difference,
    size_loss=loss_event,max_loss=maxloss,
    retained_fdp_above_nominal=mx$retained_gaussian>alpha,union_fdp_above_nominal=mx$union_gaussian>alpha,
    histories=reps)
}

run_cumulative_extension <- function(workers=2L) {
  grid<-extension_grid()
  protocol<-list(
    purpose='Prespecified extension of the cumulative-union comparison along the axis where unions failed certification most often',
    question='As signal density, signal strength, and the release horizon increase with independently refreshed summaries, how often does the cumulative union of fresh mean-closure selections lose certification, and does its SupFDR approach or exceed the nominal level?',
    grid=grid,replicates_per_configuration=1000L,hypotheses=200L,alpha=.05,prior_sd=3,
    evidence='Gaussian mixture score with known marginal variance 1/q, equal release weights',
    cross_hypothesis_correlation=.2,cumulative_noise_fraction=0,information_growth='Balanced',
    seed_rule='2026092300 + configuration id',
    comparators=c('Gaussian e-BH','fresh mean closure','retained mean closure','cumulative union of fresh selections'),
    primary_outcomes=c('uncertified union frequency','SupFDR of e-BH, retained closure, and the fresh union',
      'fresh-retained withdrawal and size loss','paired union minus retained differences'),
    recorded_before_simulation=format(Sys.time(),tz='UTC',usetz=TRUE))
  write_json(protocol,'results/R/cumulative_extension_protocol.json',pretty=TRUE,auto_unbox=TRUE)
  write_csv(grid,'results/R/cumulative_extension_design.csv')
  jobs<-lapply(seq_len(nrow(grid)),function(i)grid[i,])
  runs<-if(workers>1L && .Platform$OS.type=='unix')
    parallel::mclapply(jobs,simulate_extension,mc.cores=workers,mc.set.seed=FALSE) else lapply(jobs,simulate_extension)
  stopifnot(!any(vapply(runs,inherits,logical(1),'try-error')))
  names(runs)<-paste0('Extension ',grid$id)
  saveRDS(runs,'results/R/cumulative_extension.rds')
  for(k in c('summary','paired','trajectories'))write_csv(do.call(rbind,lapply(names(runs),function(n)
    cbind(study=n,runs[[n]]$config[,c('signals','final_snr','releases')],runs[[n]][[k]]))),
    paste0('results/R/cumulative_extension_',k,'.csv'))
  events<-do.call(rbind,lapply(names(runs),function(n) {
    z<-runs[[n]];data.frame(study=n,z$config[,c('signals','final_snr','releases')],replicates=z$histories,
      withdrawals=sum(z$withdrawal),membership_differences=sum(z$difference),size_losses=sum(z$size_loss),
      union_changes=sum(z$union_changed),uncertified_unions=sum(z$uncertified),
      retained_fdp_above_nominal=sum(z$retained_fdp_above_nominal),union_fdp_above_nominal=sum(z$union_fdp_above_nominal))
  }))
  write_csv(events,'results/R/cumulative_extension_events.csv')
  files<-c('R/38_cumulative_extension.R','R/retention_diagnostics.R','R/01_methods.R',native_sources)
  write_json(list(configurations=length(runs),replicates_per_configuration=1000L,
    protocol='results/R/cumulative_extension_protocol.json',
    analysis_type='New prespecified simulation configurations with their own seeds; original studies unchanged',
    union_certification_is_not_an_observed_FDR_test=TRUE,
    source_sha256=setNames(lapply(files,function(p)digest::digest(file=p,algo='sha256')),files),
    completed_at=format(Sys.time(),tz='UTC',usetz=TRUE)),
    'results/R/cumulative_extension_validation.json',pretty=TRUE,auto_unbox=TRUE)
  invisible(runs)
}
if(sys.nframe()==0L)run_cumulative_extension()

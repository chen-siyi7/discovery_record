source('R/01_methods.R')
source('R/retention_diagnostics.R')
source('R/11_retention_study.R')

replay_cumulative <- function(kind,config=NULL,saved=NULL) {
  if(kind=='factorial') {
    reps<-1000L;m<-200L;T<-6L
    set.seed(2026092000L+config$id)
    exponent<-if(config$growth=='Balanced')rep(1,m) else rep(c(.5,2),each=m/2)
    info<-outer(exponent,(1:T)/T,function(a,t)t^a)
    truth<-matrix(FALSE,reps,m);theta<-matrix(0,reps,m)
    for(i in seq_len(reps)) {
      alt<-sample.int(m,config$signals);truth[i,alt]<-TRUE
      theta[i,alt]<-config$final_snr*sample(c(-1,1),config$signals,TRUE)
    }
    noise<-function()sqrt(.2)*matrix(rnorm(reps),reps,m)+sqrt(.8)*matrix(rnorm(reps*m),reps,m)
    signals<-config$signals
  } else {
    design<-saved$design;m<-design$m;T<-design$T;reps<-saved$reps
    RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(saved$seed)
    all_null<-grepl('All null',kind);df<-saved$variance_df;versions<-saved$versions
    truth<-matrix(FALSE,reps,m);old_only<-new_only<-truth;theta<-matrix(0,reps,m)
    for(i in seq_len(reps))if(!all_null) {
      alt<-sample.int(m,design$n_alt)
      theta[i,alt]<-design$amplitude*sample(c(-1,1),length(alt),replace=TRUE)
      truth[i,alt]<-TRUE
      if(versions){old_only[i,alt[1:20]]<-TRUE;new_only[i,alt[21:40]]<-TRUE}
    }
    perm<-t(replicate(reps,sample.int(m)))
    ix<-cbind(rep(seq_len(reps),m),as.vector(perm))
    truthp<-matrix(truth[ix],reps,m);info<-design$info;signals<-design$n_alt
  }
  walk<-matrix(0,reps,m);mix<-matrix(-Inf,reps,m)
  names_methods<-c('gaussian','fresh_gaussian','retained_gaussian','union_gaussian')
  prev<-first<-mx<-list()
  for(n in names_methods){prev[[n]]<-matrix(FALSE,reps,m);first[[n]]<-matrix(T+1L,reps,m);mx[[n]]<-numeric(reps)}
  trajectories<-list();uncertified<-rep(FALSE,reps);union_changed<-rep(FALSE,reps)
  summaries<-events<-list();withdrawal<-difference<-loss_event<-rep(FALSE,reps)
  maxloss<-numeric(reps)
  for(t in seq_len(T)) {
    q<-info[,t];dq<-q-if(t==1)0 else info[,t-1]
    if(kind=='factorial') {
      walk<-walk+sweep(noise(),2,sqrt(dq),'*')
      z<-sqrt(config$overlap)*sweep(walk,2,sqrt(q),'/')+sqrt(1-config$overlap)*noise()
      x<-theta+sweep(z,2,sqrt(q),'/');variance<-matrix(rep(1/q,each=reps),reps,m)
      le<-log_abf(x,variance,prior_sd=3)
      active_truth<-truth
    } else {
      if(design$kind=='benchmark') z<-matrix(rnorm(reps*m),reps,m) else {
        inst<-matrix(rnorm(reps*max(design$instrument)),reps)[,design$instrument,drop=FALSE]
        outcome<-matrix(rnorm(reps*4),reps)[,design$endpoint,drop=FALSE]
        z<-sqrt(.2)*inst+sqrt(.2)*outcome+sqrt(.6)*matrix(rnorm(reps*m),reps,m)
      }
      walk<-walk+sweep(z,2,sqrt(dq),'*');mu<-theta
      if(versions){if(t<=2)mu[new_only]<-0 else mu[old_only]<-0}
      se<-design$se/sqrt(q);x<-mu+sweep(walk,2,design$se/q,'*')
      vv<-matrix(rep(se^2,each=reps),reps)
      if(is.finite(df))vv<-vv*matrix(rchisq(reps*m,df)/df,reps,m)
      le<-log_abf(x,vv,.2);le<-matrix(le[ix],reps,m)
      active_truth<-truthp
    }
    mix<-log_add(mix,le-log(T));e<-bounded_e(mix)
    ss<-list(gaussian=reject_e(e))
    ss$fresh_gaussian<-closed(e,ss$gaussian)
    ss$retained_gaussian<-persistent_e(e,prev$retained_gaussian)
    ss$union_gaussian<-prev$union_gaussian|ss$fresh_gaussian
    stopifnot(all(!prev$gaussian|ss$gaussian), all(!prev$retained_gaussian|ss$retained_gaussian))
    differs<-rowSums(ss$union_gaussian!=ss$fresh_gaussian)>0
    union_changed<-union_changed|differs
    for(i in which(differs))uncertified[i]<-uncertified[i]|!mean_consistent_R(e[i,],ss$union_gaussian[i,])
    if(t>1) {
      z<-diagnose_release(e,ss$fresh_gaussian,ss$retained_gaussian,prev$retained_gaussian,t)
      summaries[[length(summaries)+1L]]<-z$summary
      if(!is.null(z$events))events[[length(events)+1L]]<-z$events
    }
    a<-ss$fresh_gaussian;b<-ss$retained_gaussian;delta<-rowSums(a)-rowSums(b)
    withdrawal<-withdrawal|rowSums(prev$fresh_gaussian&!a)>0
    difference<-difference|rowSums(a!=b)>0;loss_event<-loss_event|delta>0;maxloss<-pmax(maxloss,delta)
    for(n in names_methods) {
      r<-ss[[n]];mx[[n]]<-pmax(mx[[n]],rowSums(r&!active_truth)/pmax(rowSums(r),1))
      first[[n]][r&first[[n]]==T+1L]<-t
      metrics<-list(fdr=rowSums(r&!active_truth)/pmax(rowSums(r),1),
        discoveries=rowSums(r),true_discoveries=rowSums(r&active_truth),
        power=if(signals>0 && any(active_truth))rowSums(r&active_truth)/signals else rep(NA_real_,reps))
      for(metric in names(metrics))trajectories[[length(trajectories)+1L]]<-
        data.frame(release=t,method=n,metric=metric,t(estimate(metrics[[metric]])))
    }
    prev<-ss
  }

  for(n in setdiff(names_methods,'union_gaussian')) {
    r<-prev[[n]];num<-rowSums(r);tp<-rowSums(r&active_truth);ref<-saved$raw[[n]]
    stopifnot(identical(num,ref$discoveries),identical(tp,ref$true_discoveries),
      identical(rowSums(r&!active_truth)/pmax(num,1),ref$fdr))
    if(kind=='factorial')stopifnot(identical(mx[[n]],ref$supfdr),
      identical(rowSums(first[[n]]*active_truth)/signals,ref$first_declaration)) else {
      stopifnot(identical(mx[[n]],ref$max_fdp))
      if(!all_null)stopifnot(identical(rowSums(first[[n]]*active_truth)/signals,ref$first_discovery))
    }
  }
  pair<-saved$pair_raw$gaussian
  if(kind=='factorial')stopifnot(identical(withdrawal,pair$withdrawal),identical(difference,pair$difference),
     identical(loss_event,pair$size_loss),identical(maxloss,pair$max_loss)) else
    stopifnot(identical(withdrawal,pair$any_withdrawal),identical(difference,pair$different),
     identical(loss_event,pair$cardinality_cost),identical(maxloss,pair$max_cardinality_cost))
  message('Cumulative replay matched: ',if(kind=='factorial')paste('factorial',config$id) else kind)
  raw<-lapply(names_methods,function(n) {
    r<-prev[[n]];num<-rowSums(r);tp<-rowSums(r&active_truth)
    list(discoveries=num,true_discoveries=tp,fdr=rowSums(r&!active_truth)/pmax(num,1),
      supfdr=mx[[n]],power=if(any(active_truth))tp/signals else rep(NA_real_,reps),
      first_declaration=if(any(active_truth))rowSums(first[[n]]*active_truth)/signals else rep(NA_real_,reps))
  });names(raw)<-names_methods

  stopifnot(identical(first$fresh_gaussian,first$union_gaussian))
  summary<-do.call(rbind,lapply(names_methods,function(n)do.call(rbind,
    lapply(names(raw[[n]]),function(metric)data.frame(method=n,metric=metric,t(estimate(raw[[n]][[metric]])))))))
  differences<-lapply(names(raw$union_gaussian),function(metric)raw$union_gaussian[[metric]]-raw$retained_gaussian[[metric]])
  names(differences)<-names(raw$union_gaussian)
  paired<-do.call(rbind,lapply(names(differences),function(metric)data.frame(metric=metric,t(estimate(differences[[metric]])))))
  list(summary=summary,paired=paired,trajectories=do.call(rbind,trajectories),raw=raw,
    uncertified=uncertified,union_changed=union_changed,
    histories=reps,matched=TRUE,config=config,
    calibrated=kind=='factorial'||is.infinite(saved$variance_df))
}

run_cumulative_reporting <- function(workers=2L) {
  factorial_saved<-readRDS('results/R/boundary_study.rds')
  direct_saved<-readRDS('results/R/retention_study.rds')
  jobs<-c(lapply(seq_along(factorial_saved),function(i)list(kind='factorial',config=factorial_saved[[i]]$config,saved=factorial_saved[[i]])),
    lapply(names(direct_saved),function(n)list(kind=n,saved=direct_saved[[n]])))
  runs<-if(workers>1L && .Platform$OS.type=='unix')parallel::mclapply(jobs,function(j)do.call(replay_cumulative,j),mc.cores=workers,mc.set.seed=FALSE) else lapply(jobs,function(j)do.call(replay_cumulative,j))
  stopifnot(!any(vapply(runs,inherits,logical(1),'try-error')))
  names(runs)<-c(paste0('Factorial ',seq_along(factorial_saved)),names(direct_saved))
  saveRDS(runs,'results/R/cumulative_reporting.rds')
  for(k in c('summary','paired','trajectories'))write_csv(do.call(rbind,lapply(names(runs),function(n)cbind(study=n,calibrated=runs[[n]]$calibrated,runs[[n]][[k]]))),paste0('results/R/cumulative_reporting_',k,'.csv'))
  events<-do.call(rbind,lapply(names(runs),function(n) {
    z<-runs[[n]];data.frame(study=n,calibrated=z$calibrated,replicates=z$histories,
      union_changes=sum(z$union_changed),uncertified_unions=sum(z$uncertified))
  }))
  write_csv(events,'results/R/cumulative_reporting_events.csv')
  files<-c('R/36_cumulative_reporting.R','R/retention_diagnostics.R','R/11_retention_study.R','R/01_methods.R',native_sources)
  write_json(list(configurations=length(runs),replicates_per_configuration=1000L,
    all_original_replicate_outcomes_match=TRUE,all_original_pair_events_match=TRUE,
    first_declarations_match_fresh_union=TRUE,
    analysis_type='Post-analysis diagnostic replay with unchanged designs and seeds',
    union_certification_is_not_an_observed_FDR_test=TRUE,
    source_sha256=setNames(lapply(files,function(p)digest::digest(file=p,algo='sha256')),files),
    completed_at=format(Sys.time(),tz='UTC',usetz=TRUE)),
    'results/R/cumulative_reporting_validation.json',pretty=TRUE,auto_unbox=TRUE)
  invisible(runs)
}
if(sys.nframe()==0L)run_cumulative_reporting()

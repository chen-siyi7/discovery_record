source("R/01_methods.R")

retention_design <- function(kind = "application") {
  if(kind == "benchmark") return(list(m=200L,T=12L,se=rep(1/sqrt(300),200),
    info=matrix(rep((1:12)/12,each=200),200),endpoint=rep(1L,200),
    instrument=seq_len(200),n_alt=20L,amplitude=.25,kind=kind))
  d <- read.csv("results/R/real_data_ledger.csv")
  ends <- c("I9_CHD","I9_AF","I9_HEARTFAIL","M13_RHEUMA")
  last <- subset(d,release==13 & available)
  set.seed(202609151)
  proteins <- sample(unique(last$protein_id),300,replace=FALSE)
  ids <- unlist(lapply(ends,function(e) match(paste(proteins,e),paste(last$protein_id,last$endpoint))))
  z <- last[ids,];stopifnot(!anyNA(ids),nrow(z)==1200)
  meta <- read.csv("results/R/endpoint_metadata.csv")
  frac <- sapply(10:13,function(t) {
    a <- meta[match(paste(ends,t),paste(meta$endpoint,meta$release)),]
    b <- meta[match(paste(ends,13),paste(meta$endpoint,meta$release)),]
    (1/a$cases+1/a$controls)^(-1)/(1/b$cases+1/b$controls)^(-1)
  })
  stopifnot(all(apply(frac,1,diff)>0))
  list(m=1200L,T=4L,se=z$wald_se_first_order,info=frac[rep(1:4,each=300),],
       endpoint=rep(1:4,each=300),instrument=match(z$variant,unique(z$variant)),
       n_alt=60L,amplitude=.2,kind=kind,protein_id=z$protein_id,variant=z$variant)
}

retention_run <- function(name,kind="application",versions=FALSE,df=Inf,all_null=FALSE,reps=1000L,seed=202609152L) {
  design <- retention_design(kind);m<-design$m;T<-design$T
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)

  truth <- matrix(FALSE,reps,m); old_only<-new_only<-truth
  theta <- matrix(0,reps,m)
  for(i in seq_len(reps)) if(!all_null) {
    alt<-sample.int(m,design$n_alt)
    theta[i,alt]<-design$amplitude*sample(c(-1,1),length(alt),replace=TRUE)
    truth[i,alt]<-TRUE
    if(versions) {old_only[i,alt[1:20]]<-TRUE;new_only[i,alt[21:40]]<-TRUE}
  }

  perm <- t(replicate(reps,sample.int(m)))
  ix <- cbind(rep(seq_len(reps),m),as.vector(perm))
  truthp<-matrix(truth[ix],reps,m); currentp<-matrix((truth & !old_only)[ix],reps,m)
  methods<-c("gaussian","fresh_gaussian","retained_gaussian","anytime_by","fresh_p","retained_p","fresh_by","retained_by")
  prev<-ever<-first<-maxfdp<-maxcurrent<-sum_power<-withdrawals<-list()
  for(method in methods) {
    prev[[method]]<-ever[[method]]<-matrix(FALSE,reps,m)
    first[[method]]<-matrix(T+1L,reps,m)
    maxfdp[[method]]<-maxcurrent[[method]]<-sum_power[[method]]<-withdrawals[[method]]<-numeric(reps)
  }
  families<-c("gaussian","p","by");pair<-list()
  for(f in families)pair[[f]]<-list(any_withdrawal=rep(FALSE,reps),different=rep(FALSE,reps),cardinality_cost=rep(FALSE,reps),max_cardinality_cost=numeric(reps))
  walk<-matrix(0,reps,m);mix<-matrix(-Inf,reps,m);minp<-matrix(0,reps,m)
  trajectories<-list();clock_start<-proc.time()[[3]]
  for(t in seq_len(T)) {
    q<-design$info[,t];dq<-q-if(t>1)design$info[,t-1] else 0
    if(kind=="benchmark") noise<-matrix(rnorm(reps*m),reps,m) else {
      inst<-matrix(rnorm(reps*max(design$instrument)),reps)[,design$instrument,drop=FALSE]
      outcome<-matrix(rnorm(reps*4),reps)[,design$endpoint,drop=FALSE]
      noise<-sqrt(.2)*inst+sqrt(.2)*outcome+sqrt(.6)*matrix(rnorm(reps*m),reps,m)
    }
    walk<-walk+sweep(noise,2,sqrt(dq),"*")
    mu<-theta
    if(versions){if(t<=2)mu[new_only]<-0 else mu[old_only]<-0}
    se<-design$se/sqrt(q)
    x<-mu+sweep(walk,2,design$se/q,"*")

    vv<-matrix(rep(se^2,each=reps),reps)
    if(is.finite(df)) vv<-vv*matrix(rchisq(reps*m,df)/df,reps,m)
    zz<-abs(x)/sqrt(vv)
    lp<-log(2)+if(is.finite(df))pt(zz,df,lower.tail=FALSE,log.p=TRUE) else pnorm(zz,lower.tail=FALSE,log.p=TRUE)
    le<-log_abf(x,vv,.2)
    lp<-matrix(lp[ix],reps,m);le<-matrix(le[ix],reps,m)
    mix<-log_add(mix,le-log(T));minp<-pmin(minp,lp+log(T))
    e<-bounded_e(mix);p<-exp(minp);cal<-by_calibrator(p)
    ss<-list(gaussian=reject_e(e),anytime_by=reject_e(cal))
    ss$fresh_gaussian<-closed(e,ss$gaussian)
    ss$retained_gaussian<-persistent_e(e,prev$retained_gaussian)
    ss$fresh_p<-closed(cal,ss$anytime_by)
    ss$retained_p<-persistent_e(cal,prev$retained_p)
    ss$fresh_by<-closed(p,mode="by")
    ss$retained_by<-closed(p,prev$retained_by,mode="by")
    for(f in families) {
      fresh<-ss[[paste0("fresh_",f)]];ret<-ss[[paste0("retained_",f)]]
      delta<-rowSums(fresh)-rowSums(ret)
      stopifnot(all(delta>=0))
      pair[[f]]$any_withdrawal<-pair[[f]]$any_withdrawal | rowSums(prev[[paste0("fresh_",f)]] & !fresh)>0
      pair[[f]]$different<-pair[[f]]$different | rowSums(fresh!=ret)>0
      pair[[f]]$cardinality_cost<-pair[[f]]$cardinality_cost | delta>0
      pair[[f]]$max_cardinality_cost<-pmax(pair[[f]]$max_cardinality_cost,delta)
    }
    for(method in methods) {
      r<-ss[[method]];n<-rowSums(r);tp<-rowSums(r & truthp);fdp<-rowSums(r & !truthp)/pmax(n,1)
      cf<-rowSums(r & !currentp)/pmax(n,1)
      if(grepl("retained",method))stopifnot(all(!prev[[method]] | r))
      withdrawals[[method]]<-withdrawals[[method]]+rowSums(prev[[method]] & !r)
      first[[method]][r & !ever[[method]]]<-t;ever[[method]]<-ever[[method]] | r
      maxfdp[[method]]<-pmax(maxfdp[[method]],fdp);maxcurrent[[method]]<-pmax(maxcurrent[[method]],cf)
      pow<-if(all_null)rep(NA_real_,reps) else tp/design$n_alt
      sum_power[[method]]<-sum_power[[method]]+pow
      for(metric in c("discoveries","true_discoveries","power","fdr")) {
        xmetric<-switch(metric,discoveries=n,true_discoveries=tp,power=pow,fdr=fdp)
        trajectories[[length(trajectories)+1L]]<-data.frame(scenario=name,method=method,release=t,metric=metric,t(estimate(xmetric)))
      }
    }
    prev<-ss
    message(name,": release ",t,"/",T)
  }
  raw<-list();summary<-list()
  for(method in methods) {
    r<-prev[[method]];n<-rowSums(r)
    delay<-if(all_null)rep(NA_real_,reps) else rowSums(first[[method]]*truthp)/design$n_alt

    metrics<-list(discoveries=n,true_discoveries=rowSums(r & truthp),
      power=if(all_null)rep(NA_real_,reps) else rowSums(r & truthp)/design$n_alt,
      fdr=rowSums(r & !truthp)/pmax(n,1),max_fdp=maxfdp[[method]],
      current_label_fdp=rowSums(r & !currentp)/pmax(n,1),withdrawals=withdrawals[[method]],
      first_discovery=delay,power_auc=sum_power[[method]]/T)
    raw[[method]]<-metrics
    for(metric in names(metrics))summary[[length(summary)+1L]]<-data.frame(scenario=name,method=method,metric=metric,t(estimate(metrics[[metric]])))
  }
  paired<-list()
  for(f in families) {
    for(metric in names(pair[[f]]))paired[[length(paired)+1L]]<-data.frame(scenario=name,family=f,metric=metric,t(estimate(as.numeric(pair[[f]][[metric]]))))
    for(metric in c("discoveries","true_discoveries","first_discovery","power_auc","fdr")) {
      delta<-raw[[paste0("retained_",f)]][[metric]]-raw[[paste0("fresh_",f)]][[metric]]
      paired[[length(paired)+1L]]<-data.frame(scenario=name,family=f,metric=paste0("retained_minus_fresh_",metric),t(estimate(delta)))
    }
  }
  list(summary=do.call(rbind,summary),paired=do.call(rbind,paired),trajectories=do.call(rbind,trajectories),
       raw=raw,pair_raw=pair,design=design,versions=versions,variance_df=df,reps=reps,seed=seed,elapsed=proc.time()[[3]]-clock_start)
}

run_retention_study <- function(reps=1000L,workers=2L) {
  configs<-list(list(name="Homogeneous cumulative",kind="benchmark"),
    list(name="Heterogeneous cumulative"),list(name="Endpoint versions",versions=TRUE),
    list(name="Estimated variance: df 30",df=30),list(name="Estimated variance: df 200",df=200),
    list(name="All null: df 30",df=30,all_null=TRUE),list(name="All null: df 200",df=200,all_null=TRUE))
  jobs<-function(x)do.call(retention_run,c(x,list(reps=reps)))
  runs<-if(workers>1L && .Platform$OS.type=="unix")parallel::mclapply(configs,jobs,mc.cores=workers,mc.set.seed=FALSE) else lapply(configs,jobs)
  stopifnot(!any(vapply(runs,inherits,logical(1),"try-error")))
  names(runs)<-vapply(configs,`[[`,character(1),"name")
  saveRDS(runs,"results/R/retention_study.rds")
  for(k in c("summary","paired","trajectories"))write_csv(do.call(rbind,lapply(runs,`[[`,k)),paste0("results/R/retention_",k,".csv"))
  write_json(list(replicates_per_scenario=reps,workers=workers,design_seed=202609151,data_seed=202609152,
    scenarios=lapply(runs,function(z)list(m=z$design$m,releases=z$design$T,versions=z$versions,
      variance_df=if(is.finite(z$variance_df))z$variance_df else "fixed",elapsed_seconds=z$elapsed)),
    source_sha256=digest::digest(file="R/11_retention_study.R",algo="sha256"),
    method_source_sha256=setNames(lapply(c("R/01_methods.R",native_sources),
      function(p)digest::digest(file=p,algo="sha256")),c("R/01_methods.R",native_sources)),
    completed_at=format(Sys.time(),tz="UTC",usetz=TRUE)),"results/R/retention_metadata.json",pretty=TRUE,auto_unbox=TRUE)
  invisible(runs)
}
if(sys.nframe()==0L)run_retention_study()

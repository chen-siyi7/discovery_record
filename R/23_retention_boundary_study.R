source('R/01_methods.R')
grid <- expand.grid(overlap=c(0,.5,1), growth=c('Balanced','Staggered'),
                    signals=c(20L,60L), final_snr=c(3,5), stringsAsFactors=FALSE)
grid$id <- seq_len(nrow(grid))
reps <- 1000L; m <- 200L; T <- 6L
protocol <- list(grid=grid,replicates=reps,m=m,releases=T,alpha=.05,
  prior_sd=3,release_weights=rep(1/T,T),cross_hypothesis_correlation=.2,
  information_exponents=list(Balanced=1,Staggered=c(.5,2)),
  signal_locations='Sampled independently of observations',
  overlap='Fraction of marginal noise variance from cumulative observations; the remainder is independently refreshed',
  seed=2026092000L,estimation='Known marginal variances; exact Gaussian p-values',
  outcomes=c('any rank change','withdrawal','composition difference','size loss',
             'discoveries','power','FDR','SupFDR','first declaration'),
  created_at=format(Sys.time(),tz='UTC',usetz=TRUE))
write_json(protocol,'results/R/boundary_protocol.json',pretty=TRUE,auto_unbox=TRUE)
write_csv(grid,'results/R/boundary_design.csv')

run_one <- function(config) {
 set.seed(2026092000L+config$id)
 alpha <- .05
 exponent <- if(config$growth=='Balanced')rep(1,m) else rep(c(.5,2),each=m/2)
 info <- outer(exponent,(1:T)/T,function(a,t)t^a)
 truth <- matrix(FALSE,reps,m); theta <- matrix(0,reps,m)
 for(i in seq_len(reps)) {
  alt <- sample.int(m,config$signals);truth[i,alt]<-TRUE
  theta[i,alt]<-config$final_snr*sample(c(-1,1),config$signals,TRUE)
 }
 methods <- c('gaussian','fresh_gaussian','retained_gaussian','anytime_by',
              'fresh_p','retained_p','fresh_by','retained_by')
 families <- c('gaussian','p','by')
 prev <- first <- mx <- list()
 for(n in methods) {prev[[n]]<-matrix(FALSE,reps,m);first[[n]]<-matrix(T+1L,reps,m);mx[[n]]<-numeric(reps)}
 pair <- lapply(families,function(f)list(withdrawal=rep(FALSE,reps),difference=rep(FALSE,reps),
             size_loss=rep(FALSE,reps),max_loss=numeric(reps)))
 names(pair)<-families
 rank_change <- rep(FALSE,reps); previous_rank <- NULL
 walk <- matrix(0,reps,m); mix <- matrix(-Inf,reps,m); minp <- matrix(0,reps,m)
 noise <- function() sqrt(.2)*matrix(rnorm(reps),reps,m)+sqrt(.8)*matrix(rnorm(reps*m),reps,m)
 for(t in seq_len(T)) {
  q<-info[,t];dq<-q-if(t==1)0 else info[,t-1]
  walk<-walk+sweep(noise(),2,sqrt(dq),'*')
  z<-sqrt(config$overlap)*sweep(walk,2,sqrt(q),'/')+sqrt(1-config$overlap)*noise()
  x<-theta+sweep(z,2,sqrt(q),'/');variance<-matrix(rep(1/q,each=reps),reps,m)
  lp<-log(2)+pnorm(abs(x)/sqrt(variance),lower.tail=FALSE,log.p=TRUE)
  le<-log_abf(x,variance,prior_sd=3)
  mix<-log_add(mix,le-log(T));minp<-pmin(minp,lp+log(T))
  e<-bounded_e(mix);p<-exp(minp);cal<-by_calibrator(p)
  ranks<-t(apply(mix,1,order))
  if(t>1)rank_change<-rank_change|rowSums(ranks!=previous_rank)>0
  previous_rank<-ranks
  ss<-list(gaussian=reject_e(e),anytime_by=reject_e(cal))
  ss$fresh_gaussian<-closed(e,ss$gaussian)
  ss$retained_gaussian<-persistent_e(e,prev$retained_gaussian)
  ss$fresh_p<-closed(cal,ss$anytime_by)
  ss$retained_p<-persistent_e(cal,prev$retained_p)
  ss$fresh_by<-closed(p,mode='by')
  ss$retained_by<-closed(p,prev$retained_by,mode='by')
  for(f in families) {
   a<-ss[[paste0('fresh_',f)]];b<-ss[[paste0('retained_',f)]]
   delta<-rowSums(a)-rowSums(b)
   stopifnot(all(delta>=0))
   pair[[f]]$withdrawal<-pair[[f]]$withdrawal|rowSums(prev[[paste0('fresh_',f)]]&!a)>0
   pair[[f]]$difference<-pair[[f]]$difference|rowSums(a!=b)>0
   pair[[f]]$size_loss<-pair[[f]]$size_loss|delta>0
   pair[[f]]$max_loss<-pmax(pair[[f]]$max_loss,delta)
  }
  for(n in methods) {
   r<-ss[[n]];fdp<-rowSums(r&!truth)/pmax(rowSums(r),1)
   mx[[n]]<-pmax(mx[[n]],fdp)
   first[[n]][r&first[[n]]==T+1L]<-t
   if(grepl('retained',n))stopifnot(all(!prev[[n]]|r))
  }
  stopifnot(all(!ss$gaussian|ss$retained_gaussian),all(!ss$anytime_by|ss$retained_p))
  prev<-ss
 }
 raw<-list();summary<-paired<-list()
 for(n in methods) {
  r<-prev[[n]];num<-rowSums(r);tp<-rowSums(r&truth)
  raw[[n]]<-list(discoveries=num,true_discoveries=tp,power=tp/config$signals,
                 fdr=rowSums(r&!truth)/pmax(num,1),supfdr=mx[[n]],
                 first_declaration=rowSums(first[[n]]*truth)/config$signals)
  for(metric in names(raw[[n]]))summary[[length(summary)+1L]]<-data.frame(config,method=n,metric=metric,t(estimate(raw[[n]][[metric]])))
 }
 for(f in families) {
  stopifnot(identical(pair[[f]]$withdrawal,pair[[f]]$difference))
  for(metric in names(pair[[f]]))paired[[length(paired)+1L]]<-data.frame(config,family=f,metric=metric,t(estimate(as.numeric(pair[[f]][[metric]]))))
  for(metric in names(raw[[1]]))paired[[length(paired)+1L]]<-data.frame(config,family=f,metric=paste0('retained_minus_fresh_',metric),t(estimate(raw[[paste0('retained_',f)]][[metric]]-raw[[paste0('fresh_',f)]][[metric]])))
 }
 out<-list(config=config,summary=do.call(rbind,summary),paired=do.call(rbind,paired),
           pair_raw=pair,raw=raw,any_rank_change=rank_change,reps=reps)
 message('Completed boundary configuration ',config$id,'/',nrow(grid))
 out
}
jobs<-split(grid,seq_len(nrow(grid)))
runs<-parallel::mclapply(jobs,run_one,mc.cores=getOption("simulation_workers", if (.Platform$OS.type == "unix") 2L else 1L),mc.set.seed=FALSE)
stopifnot(!any(vapply(runs,inherits,logical(1),'try-error')))
saveRDS(runs,'results/R/boundary_study.rds')
write_csv(do.call(rbind,lapply(runs,`[[`,'summary')),'results/R/boundary_summary.csv')
write_csv(do.call(rbind,lapply(runs,`[[`,'paired')),'results/R/boundary_pairs.csv')
write_json(list(completed_at=format(Sys.time(),tz='UTC',usetz=TRUE),configs=nrow(grid),
  reps=reps,seed_base=2026092000L,source_sha256=digest::digest(file='R/23_retention_boundary_study.R',algo='sha256'),
  method_source_sha256=setNames(lapply(c('R/01_methods.R',native_sources),
    function(p)digest::digest(file=p,algo='sha256')),c('R/01_methods.R',native_sources))),
  'results/R/boundary_metadata.json',pretty=TRUE,auto_unbox=TRUE)

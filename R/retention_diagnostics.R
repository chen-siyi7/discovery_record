boundary_completion <- function(e, previous, k) {
  stopifnot(length(e)==length(previous),sum(previous)<=k,k<=length(e))
  order_e <- order(e,seq_along(e),decreasing=TRUE)
  ans <- previous
  if(k>sum(ans))ans[head(order_e[!ans[order_e]],k-sum(ans))]<-TRUE
  ans
}

boundary_profile <- function(e, fresh, previous, alpha=.05) {
  k<-sum(fresh);m<-length(e)
  stopifnot(k>0,sum(previous)<=k,all(is.finite(e)),all(e>=0))
  candidate<-boundary_completion(e,previous,k)
  a<-seq_len(k);threshold<-a/(alpha*k)
  sp<-cumsum(sort(e[fresh]));sc<-cumsum(sort(e[candidate]))

  unselected_term<-function(set) vapply(threshold,function(c)sum(pmin(e[!set]-c,0)),numeric(1))
  up<-unselected_term(fresh);uc<-unselected_term(candidate)
  reserve<-sp+up-a*threshold
  displacement<-sp-sc;gain<-uc-up
  residual<-sc+uc-a*threshold
  sufficient<-reserve-displacement

  scale<-abs(sp)+abs(sc)+abs(up)+abs(uc)+a*threshold
  nu<-(4*m+64)*.Machine$double.eps/2
  err<-nu/(1-nu)*pmax(1,scale)
  exact_ok<-residual>=0;sufficient_ok<-sufficient>=0
  ambiguous<-which(abs(residual)<=err|abs(sufficient)<=err)
  if(length(ambiguous)) {
    q<-gmp::as.bigq(e);aq<-gmp::as.bigq(alpha)
    ps<-q[which(fresh)[order(e[fresh])]]
    cs<-q[which(candidate)[order(e[candidate])]]
    for(i in ambiguous) {
      c<-gmp::as.bigq(i)/(aq*k)
      negsum<-function(z) {d<-z-c;sum(d[d<0])}
      r<-sum(ps[seq_len(i)])+negsum(q[!fresh])-i*c
      delta<-sum(ps[seq_len(i)])-sum(cs[seq_len(i)])
      z<-sum(cs[seq_len(i)])+negsum(q[!candidate])-i*c
      exact_ok[i]<-z>=0;sufficient_ok[i]<-r-delta>=0
      reserve[i]<-as.numeric(r);displacement[i]<-as.numeric(delta)
      residual[i]<-as.numeric(z);sufficient[i]<-as.numeric(r-delta)
    }
  }
  list(candidate=candidate,no_loss=all(exact_ok),
       sufficient=all(sufficient_ok),
       profile=data.frame(a=a,reserve=reserve,displacement=displacement,
         unselected_gain=gain,residual=residual,sufficient_margin=sufficient))
}

diagnose_release <- function(e,fresh,retained,previous,release,alpha=.05) {
  n<-nrow(e);omitted<-rowSums(previous&!fresh);loss<-rowSums(fresh)-rowSums(retained)
  stopifnot(all(loss>=0),all(!previous|retained))
  events<-list()
  for(i in which(omitted>0)) {
    z<-boundary_profile(e[i,],fresh[i,],previous[i,],alpha)

    stopifnot(z$no_loss==mean_consistent_R(e[i,],z$candidate,alpha),
      z$no_loss==(loss[i]==0),!z$sufficient||z$no_loss)
    p<-z$profile; j<-which.min(p$residual)
    events[[length(events)+1L]]<-data.frame(replicate=i,release=release,
      fresh_size=sum(fresh[i,]),retained_size=sum(retained[i,]),
      omitted=omitted[i],size_loss=loss[i],sufficient=z$sufficient,
      limiting_a=p$a[j],reserve=p$reserve[j],displacement=p$displacement[j],
      unselected_gain=p$unselected_gain[j],residual=p$residual[j],
      minimum_sufficient_margin=min(p$sufficient_margin))
  }
  stopifnot(all(loss[omitted==0]==0),all(fresh[omitted==0,,drop=FALSE]==retained[omitted==0,,drop=FALSE]))
  list(summary=data.frame(release=release,comparisons=n,
      boundary_crossings=sum(omitted>0),size_losses=sum(loss>0),
      sufficient_crossings=if(length(events))sum(vapply(events,function(z)z$sufficient,logical(1))) else 0L),
    events=if(length(events))do.call(rbind,events) else NULL)
}

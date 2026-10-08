.next_double <- function(x, up = TRUE) {
  stopifnot(length(x)==1L,is.finite(x),x>0)
  b <- as.integer(writeBin(as.double(x),raw(),size=8,endian='little'))
  for(i in seq_along(b)) {
    if(up) {if(b[i]<255L){b[i]<-b[i]+1L;break};b[i]<-0L}
    else {if(b[i]>0L){b[i]<-b[i]-1L;break};b[i]<-255L}
  }
  readBin(as.raw(b),double(),n=1,size=8,endian='little')
}
.by_cutoffs_exact <- local({
  cache <- new.env(parent=emptyenv())
  function(m,alpha=.05) {
    key <- paste(m,sprintf('%a',alpha))
    if(exists(key,cache,inherits=FALSE))return(get(key,cache,inherits=FALSE))
    h <- sum(gmp::as.bigq(1,seq_len(m)))
    exact <- gmp::as.bigq(alpha)*seq_len(m)/(m*h)
    cuts <- as.numeric(exact)
    for(i in seq_len(m)) {
      while(gmp::as.bigq(cuts[i])>exact[i])cuts[i]<-.next_double(cuts[i],FALSE)
      while(gmp::as.bigq(.next_double(cuts[i],TRUE))<=exact[i])cuts[i]<-.next_double(cuts[i],TRUE)
    }
    assign(key,cuts,cache);cuts
  }
})
.by_bins_exact <- function(p,m,alpha=.05) {
  b <- findInterval(p,.by_cutoffs_exact(m,alpha),left.open=TRUE)+1L
  b[b>m] <- 0L;b
}
.roundoff <- function(terms,magnitude) {
  nu <- (terms+16)*.Machine$double.eps/2
  nu/(1-nu)*pmax(1,magnitude)
}
.e_pass <- function(e,k,m,alpha) {
  lhs <- (alpha*e)*k;rhs <- rep(m,length(lhs))
  err <- .roundoff(3,abs(lhs)+m);out <- lhs>=rhs
  amb <- which(is.finite(lhs) & abs(lhs-rhs)<=err)
  for(i in amb)out[i]<-gmp::as.bigq(alpha)*gmp::as.bigq(e[i])*rep(k,length.out=length(e))[i]>=m
  out
}
.e_cap <- function(m,alpha) {
  cap<-m/alpha
  if(gmp::as.bigq(alpha)*gmp::as.bigq(cap)<m)cap<-.next_double(cap,TRUE)
  cap
}

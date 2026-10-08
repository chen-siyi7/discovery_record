source('R/numerical_helpers.R')
certified_sets_exact <- function(values,sets,alpha=.05,mode='e') {
  m<-length(values);sizes<-rowSums(sets);A<-sets[-1,,drop=FALSE];d<-rowSums(A)
  calibrated<-NULL;info<-attr(values,'BY_calibration')
  if(mode=='e' && !is.null(info) && identical(info$alpha,alpha) &&
     identical(dim(values),dim(info$bins))) {
    expected<-m/(alpha*pmax(1,info$bins));expected[info$bins==0L]<-0
    if(identical(as.vector(values),as.vector(expected)))calibrated<-info$bins
  }
  bounds<-lapply(seq_len(nrow(A)),function(i) {
    if(!is.null(calibrated)) {
      bins<-calibrated[A[i,]];bins<-bins[bins>0L]
      return(if(length(bins))m*sum(gmp::as.bigq(1,bins))/d[i] else gmp::as.bigq(0))
    }
    if(mode=='e')return(gmp::as.bigq(alpha)*sum(gmp::as.bigq(values[A[i,]]))/d[i])
    bins<-.by_bins_exact(values[A[i,]],d[i],alpha);bins<-bins[bins>0L]
    if(length(bins))sum(gmp::as.bigq(1,bins)) else gmp::as.bigq(0)
  })
  b<-vapply(bounds,as.numeric,numeric(1));counts<-tcrossprod(A,sets)
  good<-rep(FALSE,nrow(sets));good[sizes==0L]<-TRUE
  for(r in which(sizes>0L)) {
    f<-counts[,r]/sizes[r];err<-.roundoff(3,abs(f)+abs(b))
    if(any(f-b>err))next
    amb<-which(abs(f-b)<=err)
    good[r]<-all(vapply(amb,function(i)gmp::as.bigq(counts[i,r],sizes[r])<=bounds[[i]],logical(1)))
  }
  good
}

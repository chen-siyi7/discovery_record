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

.calibrated_bins <- function(e, alpha) {
  info <- attr(e, "BY_calibration")
  if (is.null(info) || !identical(info$alpha, alpha) ||
      !identical(dim(e), dim(info$bins))) return(NULL)
  ref <- ncol(e) / (alpha * pmax(1, info$bins)); ref[info$bins == 0L] <- 0
  if (!identical(as.vector(e), as.vector(ref))) return(NULL)
  info$bins
}

as_matrix <- function(x) {
  if (is.matrix(x)) x else matrix(x, nrow = 1L)
}

.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between zero and one.", call. = FALSE)
  }
  invisible(alpha)
}

.numeric_matrix <- function(x, name) {
  x <- as_matrix(x)
  if (!is.numeric(x) || anyNA(x) || nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must be a nonempty numeric vector or matrix without NA/NaN.",
         call. = FALSE)
  }
  x
}

.logical_matrix <- function(x, template, name) {
  x <- as_matrix(x)
  if (!is.logical(x) || anyNA(x) || !identical(dim(x), dim(template))) {
    stop(name, " must be logical, without NA, and match the input dimensions.",
         call. = FALSE)
  }
  x
}

log_add <- function(x, y) {
  a <- pmax(x, y)
  ans <- a + log(exp(x - a) + exp(y - a))
  ans[is.infinite(a) & a < 0] <- -Inf
  ans[is.infinite(a) & a > 0] <- Inf
  ans
}

.step_up <- function(x, thresholds, decreasing) {
  out <- matrix(FALSE, nrow(x), ncol(x))
  for (r in seq_len(nrow(x))) {
    sorted <- sort(x[r, ], decreasing = decreasing)
    passed <- if (decreasing) sorted >= thresholds else sorted <= thresholds
    k <- if (any(passed)) max(which(passed)) else 0L
    if (k > 0L) {
      out[r, ] <- if (decreasing) x[r, ] >= thresholds[k] else
        x[r, ] <= thresholds[k]
    }
  }
  out
}

reject_p <- function(lp, alpha=.05, by=FALSE, log.p=TRUE) {
  .check_alpha(alpha);lp<-.numeric_matrix(lp,'lp')
  if(any(lp>if(log.p)0 else 1) || (!log.p && any(lp<0)))stop('Invalid p-values.')
  m<-ncol(lp);cuts<-if(by).by_cutoffs_exact(m,alpha) else alpha*seq_len(m)/m
  if(by && log.p)return(.step_up(exp(lp),cuts,FALSE))
  .step_up(lp,if(log.p)log(cuts) else cuts,FALSE)
}

reject_loge <- function(le,alpha=.05) {reject_e(bounded_e(le,alpha),alpha)}

reject_e <- function(e,alpha=.05) {
  .check_alpha(alpha);e<-.numeric_matrix(e,'e')
  if(any(e<0))stop('e must be nonnegative.')
  bins<-.calibrated_bins(e,alpha);m<-ncol(e)
  if(!is.null(bins)){scores<-bins;scores[bins==0L]<-Inf;return(.step_up(scores,seq_len(m),FALSE))}
  out<-matrix(FALSE,nrow(e),m)
  for(r in seq_len(nrow(e))) {
    passed<-.e_pass(sort(e[r,],decreasing=TRUE),seq_len(m),m,alpha)
    k<-if(any(passed))max(which(passed)) else 0L
    if(k)out[r,]<-.e_pass(e[r,],k,m,alpha)
  }
  out
}

bounded_e <- function(le,alpha=.05) {
  .check_alpha(alpha);le<-.numeric_matrix(le,'le');cap<-.e_cap(ncol(le),alpha)
  ans<-exp(pmin(le,log(cap)));ans[le>=log(cap)]<-cap;ans
}

log_abf <- function(x, variance, prior_sd = .2) {
  if (!is.numeric(x) || any(!is.finite(x)) || !is.numeric(variance) ||
      !length(variance) || any(!is.finite(variance) | variance <= 0) ||
      !is.numeric(prior_sd) || length(prior_sd) != 1L ||
      !is.finite(prior_sd) || prior_sd < 0) {
    stop("Use finite x, positive finite variance, and a nonnegative prior_sd.",
         call. = FALSE)
  }
  r <- prior_sd^2 / variance
  -.5 * log1p(r) + .5 * x^2 / variance * r / (1 + r)
}

log_adjuster <- function(y) {
  if (!is.numeric(y) || anyNA(y) || any(y < 0)) {
    stop("y must be nonnegative and contain no NA/NaN.", call. = FALSE)
  }
  ans <- y * 0 - log(2)
  small <- y > 0 & y < 1e-3
  regular <- is.finite(y) & y >= 1e-3
  ans[small] <- log(.5 + y[small]/6 + y[small]^2/24 + y[small]^3/120)
  z <- y[regular]
  ans[regular] <- z + log1p(-(1 + z) * exp(-z)) - 2 * log(z)
  ans[is.infinite(y)] <- Inf
  ans
}

by_calibrator <- function(p,alpha=.05) {
  .check_alpha(alpha);p<-.numeric_matrix(p,'p')
  if(any(!is.finite(p)|p<0|p>1))stop('p must lie in [0,1].')
  m<-ncol(p);bins<-matrix(.by_bins_exact(as.vector(p),m,alpha),nrow(p),m)
  q<-matrix(m/(alpha*pmax(1,bins)),nrow(p),m);q[bins==0L]<-0
  attr(q,'BY_calibration')<-list(alpha=alpha,bins=bins);q
}

.closure_feasible <- function(values,selected,alpha,mode,harmonic=NULL,bins=NULL) {
  k<-sum(selected);m<-length(values);if(!k)return(TRUE)
  if(mode=='e') {
    exact<-function(ix) {
      if(!is.null(bins)) {
        b<-bins[ix];z<-gmp::as.bigq(rep(0,length(b)));nz<-b>0L
        z[nz]<-gmp::as.bigq(m)/b[nz];return(z)
      }
      z<-gmp::as.bigq(alpha)*gmp::as.bigq(values[ix]);z[z>m]<-m;z
    }
    x<-if(is.null(bins))pmin(alpha*values,m) else ifelse(bins>0,m/pmax(1,bins),0)
    si<-which(selected);ui<-which(!selected)
    s<-x[si];u<-x[ui];cs<-c(0,cumsum(s));cu<-c(0,cumsum(u))
    for(a in seq_len(k)) {
      d<-k*u-a;below<-d<0;amb<-which(abs(d)<=.roundoff(3,abs(k*u)+a))
      if(length(amb))below[amb]<-as.logical(k*exact(ui[amb])<a)
      b<-sum(below);lhs<-k*(cs[a+1]+cu[b+1]);rhs<-a*(a+b)
      err<-.roundoff(a+b+3,abs(lhs)+rhs)
      if(lhs-rhs < -err)return(FALSE)
      if(abs(lhs-rhs)<=err && k*sum(exact(c(si[seq_len(a)],ui[seq_len(b)])))<rhs)return(FALSE)
    }
  } else {
    s<-values[selected];u<-values[!selected];q<-length(u)
    for(n in seq_len(m)) {
      bs<-.by_bins_exact(s,n,alpha);bu<-.by_bins_exact(u,n,alpha)
      cs<-c(0,cumsum(ifelse(bs>0,1/pmax(1,bs),0)))
      cu<-c(0,cumsum(ifelse(bu>0,1/pmax(1,bu),0)))
      for(a in seq.int(max(1L,n-q),min(k,n))) {
        b<-n-a;lhs<-k*(cs[a+1]+cu[b+1]);err<-.roundoff(n+3,abs(lhs)+a)
        if(lhs-a < -err)return(FALSE)
        if(abs(lhs-a)<=err) {
          z<-c(bs[seq_len(a)],bu[seq_len(b)]);z<-z[z>0]
          if(k*sum(gmp::as.bigq(1,z))<a)return(FALSE)
        }
      }
    }
  }
  TRUE
}

closed <- function(values,mandatory=NULL,alpha=.05,mode='e') {
  .check_alpha(alpha);mode<-match.arg(mode,c('e','by'))
  values<-.numeric_matrix(values,'values')
  if(any(!is.finite(values)|values<0)||(mode=='by'&&any(values>1)))stop('Invalid closure input.')
  nr<-nrow(values);m<-ncol(values)
  bins<-if(mode=='e').calibrated_bins(values,alpha) else NULL
  if(is.null(mandatory))mandatory<-matrix(FALSE,nr,m)
  mandatory<-.logical_matrix(mandatory,values,'mandatory');out<-matrix(FALSE,nr,m)
  for(r in seq_len(nr)) {
    ord<-order(values[r,],decreasing=mode=='by',method='radix')
    v<-values[r,ord];forced<-mandatory[r,ord];min_k<-sum(forced)
    available<-rev(which(!forced));found<-FALSE
    for(k in seq.int(m,min_k)) {
      selected<-forced;extra<-k-min_k
      if(extra>0L)selected[available[seq_len(extra)]]<-TRUE
      if(.closure_feasible(v,selected,alpha,mode,bins=if(is.null(bins))NULL else bins[r,ord])){found<-TRUE;break}
    }
    if(!found)stop('No certified superset at row ',r,'.')
    out[r,ord]<-selected
  }
  out
}

persistent_e <- function(e, previous = NULL, alpha = .05) {
  mandatory <- reject_e(e, alpha)
  if (!is.null(previous)) {
    previous <- .logical_matrix(previous, mandatory, "previous")
    mandatory <- mandatory | previous
  }
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
  if (!length(x)) {
    return(c(mean = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  mu <- mean(x)
  se <- stats::sd(x) / sqrt(length(x))
  c(mean = mu, mcse = se, lower = mu - 1.96 * se, upper = mu + 1.96 * se)
}

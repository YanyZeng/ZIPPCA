#' @title ZIPPCA
#' @description Zero-Inflated Probabilistic PCA (ZIPPCA),for dimension reduction and data ordination of multivariate abundance data,
#' and propose an efficient variational approximation method for estimation,inference, and prediction.
#'
#' @param X matrix of observations.
#' @param V vector of the sample covariate.
#' @param family distribution of models. Two options are "poisson" and "negative.binomial". Defaults to "negative.binomial".
#' @param n.factors the rank or number of factors, after dimensional reduction. Defaults to 2.
#' @param rank logical, if TRUE, the rank or number of factors, is chosen from 1 to 5 by HIC (hybrid information criterion). Defaults to FALSE.
#' @param trace logical, defaults to \code{FALSE}. if \code{TRUE} each current iteration step information will be printed.
#' @param maxit maximum number of iterations within \code{optim} and \code{constrOptim} function, defaults to 100.
#' @param parallel logical, if TRUE, use parallel toolbox to accelerate.

#' @return
#'
#'
#'  \item{VLB }{ variational lower bound of log likelihood}
#'  \item{lvs}{list of latent variables
#'  \itemize{
#'    \item{pi }{ the probabilities of excess zeros}
#'    \item{factor_scores }{ coordinates or factor scores in low-dimensional subspace}
#'    }}
#'  \item{params}{list of model parameters
#'  \itemize{
#'    \item{factor_coefs_j }{ coefficients of latent variables fator scores or factor loadings}
#'    \item{factor_coefs_0 }{ taxon-specific intercepts}
#'    \item{alpha }{ sample-specifc intercepts that adjusts for the sequencing depth}
#'    \item{gamma }{ coeffcients of sample covariate}
#'    \item{dispersion }{ taxon-specific over-dispersion parameter for negative binomial distribution}
#'    }}
#'  \item{hic}{ the number of the rank selection, chosen by HIC type information criterion}

#'  
#' @examples
#'  n.n = 100
#'  n.m = 50
#'  n.factors = 2
#'  set.seed(37)
#'  f <- matrix(0,nrow = n.n, ncol = n.factors)
#'  for(i in 1:n.n){
#'    f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
#'  }
#'  betaj <- matrix(0,nrow = n.m, ncol = n.factors)
#'  for(j in 1:n.m){
#'    betaj[j,] <- rnorm(n.factors, mean =0, sd = 1)
#'  }
#'  beta0 <- rep(1,n.m)
#'  alpha <- rep(0,n.n)
#'  l <- matrix(alpha,n.n,n.m) +matrix(beta0,n.n,n.m,byrow=TRUE) + f %*% t(betaj)
#'  lambda <- exp(l)
#'  eta_j <- runif(n.m, 0,1)
#'  z <- matrix(0,n.n,n.m)
#'  for(i in 1:n.n){
#'    z[i,] <- rbinom(n.m, size=1, prob=eta_j)
#'  }
#'  X <- matrix(0,n.n,n.m,byrow = TRUE)
#'  for(i in 1:n.n){
#'    for(j in 1:n.m){
#'      X[i,j] <- rnbinom(n=1,size=10,mu=lambda[i,j])
#'    }
#'  }
#'  X[z==1]=0
#'  zerorow <- which(rowSums(X)==0)
#'  if(length(zerorow) >0 ){
#'    X <- X[-zerorow,];f <- f[-zerorow,]
#'  }
#'  zerocol <- which(colSums(X)==0)
#'  if(length(zerocol) >0 ){
#'    X <- X[,-zerocol];betaj <- t(t(betaj)[,-zerocol])
#'  }
#'  result <- ZIPPCA::ZIPPCApn(X)
#'  f_coordinates <- result$lvs$factor_scores
#' @export
#'

ZIPPCApn <- function(X, V=NULL, family = "negative.binomial", n.factors=2, rank=FALSE,
                     trace = FALSE, maxit = 100, parallel=TRUE){

  ZIPNVA <- function(X,V,family = "negative.binomial",n.factors,trace,maxit) {

    n.s<-dim(X)[1]; n.f<-dim(X)[2];
    if(is.null(V)){Y <- 0
    }else if(is.numeric(V)){Y <- V
    }else{
      Y <- as.numeric(as.factor(V))-1}

    out.list <- list()

    ### Initialization
    X_sc <- scale(log2(1+X),center = T,scale = T)
    re <- svd(X_sc,n.factors,n.factors)
    if(n.factors==1){
      factor_scores <- new.factor_scores <- re$u * (re$d[1])
    }else{factor_scores <- new.factor_scores <- re$u %*% diag(re$d[1:n.factors])}

    factor_coefs_j <- new.factor_coefs_j <- re$v
    factor_coefs_0 <-  new.factor_coefs_0 <-  rep(1,n.f)
    factor_coefs <- cbind(factor_coefs_0,factor_coefs_j)
    pzero.col <- apply(X, 2, function(x) {sum(x==0)/n.s})
    z.hat <- pi <- new.pi <- t(ifelse(t(X)==0, pzero.col, 0))
    eta <- new.eta <- round(apply(pi, 2, mean), 6)
    fit <- mvabund::manyglm(X ~ 1, family = family)
    phi <- fit$phi + 1e-5
    dispersion <- new.dispersion <- NULL
    if(any(phi>100))phi[phi>100]=100; if(any(phi<0.10))phi[phi<0.10]=0.10;
    new.dispersion <- dispersion <- 1/phi
    sigma <- new.sigma <- matrix(0.01,n.s,n.factors)
    alpha <- new.alpha <- rep(0,n.s)
    gamma <-  new.gamma <-  rep(1,n.f)

    ### VA iteration about variational lower bound
    cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff=1e5;eps = 1e-4;max.iter = 100;
    m.cur.logfunc <- -1e6; b.cur.logfunc <- -1e6; d.cur.logfunc <- -1e6; s.cur.logfunc <- -1e6;

    while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
      if(trace) cat("Iteration:", iter, "\n")
      ## VLB
      VLB <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
        x2 <- x
        new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
        new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
        new.alpha <- v[1:n.s]
        new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
        new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
        new.gamma <- x2[1:n.f]
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- pi
        new.eta <- eta
        new.dispersion <- d

        pi.new.mat <- matrix(0,n.s,n.f,byrow = TRUE)
        dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y +matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
        if(family == "negative.binomial"){
          e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          a <- (1-new.pi)*(lgamma(X+dispersion.mat)-lfactorial(X)-lgamma(dispersion.mat)-(X+dispersion.mat)*log(1+dispersion.mat/exp(e.mat))+dispersion.mat*log(dispersion.mat)-dispersion.mat*ll)
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          b <- -(1-new.pi)*dispersion.mat*log(1+exp(e.mat)/dispersion.mat)
        }
        if(family == "poisson"){
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          a <- (1-new.pi)*(X * ll - lfactorial(X) - exp(e.mat))
          b <- (1-new.pi)*(- exp(e.mat))
        }
        y1 <- ifelse(X>0,a,b)
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2)) }
        fun3 <- function(i) {
          new.eta <- new.eta+1e-8
          new.pi <- new.pi+1e-8
          pi.mat <-  (1-new.pi[i,])*log((1-new.eta)/(1-new.pi[i,]))+ new.pi[i,]*log(new.eta/new.pi[i,])
          sum(pi.mat)
        }
        y <- sum(y1)+ sum(sapply(1:n.s,fun2))+sum(sapply(1:n.s,fun3))
        return(y)
      }

      log_base <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
        x2 <- x
        new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
        new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
        new.alpha <- v[1:n.s]
        new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
        new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
        new.gamma <- x2[1:n.f]
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- pi
        new.eta <- eta
        new.dispersion <- d

        pi.new.mat <- matrix(0,n.s,n.f,byrow = TRUE)
        dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y +matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
        if(family == "negative.binomial"){
          c <- (1-new.pi)*(lgamma(X+dispersion.mat)-lfactorial(X)-lgamma(dispersion.mat)+X*log(exp(ll)/(dispersion.mat+exp(ll)))+dispersion.mat*log(dispersion.mat/(dispersion.mat+exp(ll))))

        }
        if(family == "poisson"){
          c <- (1-new.pi)*(X * ll - lfactorial(X) - exp(ll))
        }

        y <- sum(c)
        return(y)
      }

      ## update eta
      new.eta <- round(apply(new.pi, 2, mean), 6)

      ## update beta,beta0,gamma
      beta_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
        x2 <- x
        new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
        new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
        new.alpha <- v[1:n.s]
        new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
        new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
        new.gamma <- x2[1:n.f]
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- pi
        new.eta <- eta
        new.dispersion <- d

        dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
        bl<- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) +matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + new.factor_scores %*% t(new.factor_coefs_j)
        if(family=="negative.binomial"){
          e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          a <- (1-new.pi)*(-(X+dispersion.mat)*log(1+dispersion.mat/exp(e.mat))-dispersion.mat*bl)
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          b <- -(1-new.pi)*dispersion.mat*log(1+exp(e.mat)/dispersion.mat)
        }
        if(family=="poisson"){
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          a <- (1-new.pi)*(X*bl-exp(e.mat))
          b <- (1-new.pi)*(-exp(e.mat))
        }
        y1 <- ifelse(X>0,a,b)
        y <- sum(y1)
        return(y)
      }
      beta_grad_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
        x2 <- x
        new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
        new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
        new.alpha <- v[1:n.s]
        new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
        new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
        new.gamma <- x2[1:n.f]
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- pi
        new.eta <- eta
        new.dispersion <- d

        b3 <- b2 <-  b30 <- b20 <- beta_grad <- beta0_grad <- ga <- ga0 <- gamma_grad <- NULL
        dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y +matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
        if(family=="negative.binomial"){
          e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          e.mat2 <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          for(w in 1:n.factors){
            b3 <- c(b3,(1-new.pi)*(sweep(-dispersion.mat,1,new.factor_scores[,w],"*") + sweep((new.sigma[,w]%*%t(new.factor_coefs_j[,w])),1,-new.factor_scores[,w],"+") * (-(X + dispersion.mat) * (dispersion.mat / (exp(e.mat) + dispersion.mat)))))
            b30 <- c(b30,-(1-new.pi)*(sweep((new.sigma[,w]%*%t(new.factor_coefs_j[,w])),1,new.factor_scores[,w],"+") * (exp(e.mat2)/(1+exp(e.mat2)/dispersion.mat))))
          }
          b2 <- (1-new.pi)*((X+dispersion.mat)*(dispersion.mat/(exp(e.mat)+dispersion.mat))-dispersion.mat)
          b20 <- -(1-new.pi)*(exp(e.mat2)/(1+exp(e.mat2)/dispersion.mat))
          ga <- (1-new.pi)*((X+dispersion.mat)*(dispersion.mat/(exp(e.mat)+dispersion.mat))-dispersion.mat)*Y
          ga0 <- -(1-new.pi)*(exp(e.mat2)/(1+exp(e.mat2)/dispersion.mat))*Y
        }
        if(family=="poisson"){
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          for(w in 1:n.factors){
            b3 <- c(b3,(1-new.pi)*(sweep(X,1,new.factor_scores[,w],"*") - sweep((new.sigma[,w]%*%t(new.factor_coefs_j[,w])),1,new.factor_scores[,w],"+") * exp(e.mat)))
            b30 <- c(b30,(1-new.pi)*(- sweep((new.sigma[,w]%*%t(new.factor_coefs_j[,w])),1,new.factor_scores[,w],"+") * exp(e.mat)))
          }
          b2 <- (1-new.pi)*(X-exp(e.mat))
          b20 <- (1-new.pi)*(-exp(e.mat))
          ga <- (1-new.pi)*(X-exp(e.mat))*Y
          ga0 <- (1-new.pi)*(-exp(e.mat))*Y
        }
        X2 <- matrix(X,n.s,n.f*n.factors)
        b3 <- matrix(b3,n.s,n.f*n.factors)
        b30 <- matrix(b30,n.s,n.f*n.factors)
        beta_grad <- ifelse(X2>0,b3,b30)
        beta0_grad <- ifelse(X>0,b2,b20)
        gamma_grad <- ifelse(X>0,ga,ga0)
        return(c(colSums(beta_grad), colSums(beta0_grad), colSums(gamma_grad)))
      }
      q <- try(optim(c(factor_coefs_j,factor_coefs_0,gamma), v = c(new.factor_scores,new.alpha), s = c(new.sigma), pi = new.pi,eta = new.eta,d=new.dispersion, method = "BFGS", fn = beta_f, gr = beta_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j; new.factor_coefs_0 <- factor_coefs_0;new.gamma <- gamma
      }else{
        if(iter > 1 && b.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_j <- factor_coefs_j; new.factor_coefs_0 <- factor_coefs_0;new.gamma <- gamma
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_j <- matrix(q$par[1:(n.f * n.factors)],n.f,n.factors); x2 <- q$par[-(1:(n.f * n.factors))];
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)];
          new.gamma <- x2[1:n.f]
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ## update dispersion
      if(family=="negative.binomial") {
        dispersion_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          a <- (1-new.pi)*(lgamma(X+dispersion.mat)-lgamma(dispersion.mat)-(X+dispersion.mat)*log(1+dispersion.mat/exp(e.mat))+dispersion.mat*log(dispersion.mat)-dispersion.mat*ll)
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          b <- -(1-new.pi)*dispersion.mat*log(1+exp(e.mat)/dispersion.mat)
          y1 <- ifelse(X>0,a,b)
          y <- sum(y1)
          return(y)
        }

        dispersion_grad_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          grad.dispersion <-grad.d <- grad.d0 <- NULL
          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          grad.d <- (1-new.pi)*(digamma(X+dispersion.mat)-digamma(dispersion.mat)-log(1+dispersion.mat/exp(e.mat))-(X+dispersion.mat)/
                                  (exp(e.mat)+dispersion.mat)+log(dispersion.mat)+1-ll)
          e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
          grad.d0 <- -(1-new.pi)*(log(1+exp(e.mat)/dispersion.mat)-exp(e.mat)/(1+exp(e.mat)/dispersion.mat)/dispersion.mat)
          grad.dispersion <- ifelse(X>0,grad.d,grad.d0)
          return(colSums(grad.dispersion))
        }
        q <- try(optim(c(dispersion),x=c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), s = c(new.sigma), pi = new.pi,eta = new.eta, method = "BFGS", fn = dispersion_f, gr = dispersion_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
        if("try-error" %in% class(q)){ new.dispersion <- dispersion
        }else{
          if(iter > 1 && d.cur.logfunc > q$value){if(trace)
            cat("Optimization of dispersion did not improve on iteration step ",iter,"\n");
            new.dispersion <- dispersion;
          }else{
            if(trace) cat("Model parameters dispersion updated","\n")
            new.dispersion <- q$par[1:n.f];
            if(q$convergence != 0) { if(trace) cat("Optimization of dispersion did not converge on iteration step ", iter,"\n") }
          }
        }
      }

      delta.pi.required <- 1e-3; p.iter <- 1; p.max.iter <- 10; delta.pi <- abs(pi);
      while(!all(delta.pi < delta.pi.required) & (p.iter < p.max.iter)){
        ## update pi
        new.pi <- matrix(0,n.s,n.f,byrow = TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
        e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
        if(family=="negative.binomial"){
          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          sum1 <- exp(-dispersion.mat*log(1+exp(e.mat)/dispersion.mat))
        }
        if(family=="poisson"){
          sum1 <- exp( - exp(e.mat))
        }
        for(i in 1:n.s){
          new.pi[i,] <- new.eta/(new.eta+(1-new.eta)*sum1[i,]+1e-8)
        }
        new.pi[X!=0]=0

        ## update sigma
        s_f <-function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL){
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          y1 <- rep(0,n.s)
          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          # for(i in 1:n.s){
          #   y1[i] <- -1/2*(sum(new.sigma[i,])-log(det(diag(new.sigma[i,]))))
          # }

          for(i in 1:n.s){
            if(n.factors==1){
              y1[i] <- -1/2*((new.sigma[i])-log(new.sigma[i]))
            }else{
              y1[i] <- -1/2*(sum(new.sigma[i,])-log(det(diag(new.sigma[i,]))))
            }}

          if(family=="negative.binomial"){
            e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            a <- (1-new.pi)*(-(X+dispersion.mat)*log(1+dispersion.mat/exp(e.mat)))
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            b <- -(1-new.pi)*dispersion.mat*log(1+exp(e.mat)/dispersion.mat)
            y2 <- ifelse(X>0,a,b)
          }
          if(family=="poisson"){
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            y2 <- -(1-new.pi)*exp(e.mat)
          }
          y <- sum(y1)+sum(y2)
          return(y)
        }
        s_grad_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL){
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          grad.sigma <- matrix(0,n.s,n.factors,byrow = TRUE)
          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y +matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          if(family=="negative.binomial"){
            e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            mu <- -0.5*(1-new.pi)*((X+dispersion.mat)*(dispersion.mat/(exp(e.mat)+dispersion.mat)))
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            mu2 <- -0.5*(1-new.pi)*(exp(e.mat)/(1+exp(e.mat)/dispersion.mat))
            mu.mat <- ifelse(X>0,mu,mu2)
          }
          if(family=="poisson"){
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            mu.mat <- -0.5*(1-new.pi)*exp(e.mat)
          }
          beta2 <- new.factor_coefs_j^2
          # for(i in 1:n.s){
          #   grad.sigma[i,] <- diag(-0.5*(diag(rep(1,n.factors))-(diag(new.sigma[i,]))^-1) + diag(apply(mu.mat[i,]*beta2,2,sum)))
          # }

          for(i in 1:n.s){
            if(n.factors==1){
              grad.sigma[i] <- diag(-0.5*(diag(rep(1,n.factors))-(new.sigma[i])^-1) + sum(mu.mat[i,]*beta2))
            }else{
              grad.sigma[i,] <- diag(-0.5*(diag(rep(1,n.factors))-(diag(new.sigma[i,]))^-1) + diag(apply(mu.mat[i,]*beta2,2,sum)))
            }}

          return(c(grad.sigma))
        }
        q <- try(constrOptim(c(sigma), method = "BFGS", f = s_f, grad = s_grad_f, x = c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), pi = new.pi,eta = new.eta,d=new.dispersion,ui=diag(1,n.s*n.factors,n.s*n.factors),ci=rep(1e-8,n.s*n.factors),control = list(trace = 0, fnscale = -1, maxit = maxit,reltol=1e-6)), silent = TRUE)
        if("try-error" %in% class(q)){ new.sigma <- sigma
        }else{
          if(iter > 1 && s.cur.logfunc > q$value){if(trace)
            cat("Optimization of sigma did not improve on iteration step ",iter,"\n");
            new.sigma <- sigma;
          }else{
            if(trace) cat("Variational parameters sigma updated","\n")
            new.sigma <- q$par[1:(n.factors*n.s)]; new.sigma <- matrix(new.sigma,n.s,n.factors);
            if(q$convergence != 0) { if(trace) cat("Optimization of sigma did not converge on iteration step ", iter,"\n") }
          }
        }

        ## update factor_scores and alpha
        m_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL){
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          bl<- new.factor_scores %*% t(new.factor_coefs_j)+ matrix(new.alpha,n.s,n.f)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          if(family=="negative.binomial"){
            e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            a <- (1-new.pi)*(-(X+dispersion.mat)*log(1+dispersion.mat/exp(e.mat))-dispersion.mat*bl)
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            b <- -(1-new.pi)*dispersion.mat*log(1+exp(e.mat)/dispersion.mat)
          }
          if(family=="poisson"){
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            a <- (1-new.pi)*(X * bl - exp(e.mat))
            b <- (1-new.pi)*(- exp(e.mat))
          }
          y3 <- ifelse(X>0,a,b)
          foo2 <- function(i) { -0.5 *  sum(factor_scores[i,]^2) }
          y <- sum(y3) + sum(sapply(1:n.s,foo2))
          return(y)
        }
        m_grad_f <- function(x,v=NULL,s=NULL,pi=NULL,eta=NULL,d=NULL) {
          x2 <- x
          new.factor_scores <- new.factor_coefs_j <- new.factor_coefs_0 <- new.alpha <- new.sigma <- NULL
          new.factor_scores <- matrix(c(v[1:(n.s*n.factors)]),n.s,n.factors); v <- v[-(1:(n.s*n.factors))]
          new.alpha <- v[1:n.s]
          new.factor_coefs_j <- matrix(c(x2[1:(n.f * n.factors)]),n.f,n.factors); x2 <- x2[-(1:(n.f * n.factors))]
          new.factor_coefs_0 <- x2[1:n.f]; x2 <- x2[-(1:n.f)]
          new.gamma <- x2[1:n.f]
          new.sigma <- matrix(s,n.s,n.factors)
          new.pi <- pi
          new.eta <- eta
          new.dispersion <- d

          grad.m <- grad.alpha <- NULL
          dispersion.mat <- matrix(new.dispersion,n.s,n.f,byrow=TRUE)
          ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y +matrix(new.alpha,n.s,n.f) + new.factor_scores %*% t(new.factor_coefs_j)
          if(family=="negative.binomial"){
            e.mat <- ll - 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            sum1 <- (1-new.pi)*((X+dispersion.mat)*(dispersion.mat/(exp(e.mat)+dispersion.mat))-dispersion.mat)
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            sum2 <- -(1-new.pi)*(exp(e.mat)/(1+exp(e.mat)/dispersion.mat))
          }
          if(family=="poisson"){
            e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)
            sum1 <- (1-new.pi)*(X-exp(e.mat))
            sum2 <- (1-new.pi)*(-exp(e.mat))
          }
          sum <- ifelse(X>0,sum1,sum2)
          for(l in 1:n.factors) {
            grad.m <- c(grad.m,rowSums(sweep(sum,2,new.factor_coefs_j[,l],"*"))-new.factor_scores[,l])
          }
          grad.alpha <-rowSums(sum)

          return(c(grad.m,grad.alpha))
        }
        q <- try(optim(c(factor_scores,alpha), method = "BFGS", fn = m_f, gr = m_grad_f, x = c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), s = c(new.sigma), pi = new.pi,eta = new.eta,d=new.dispersion,control = list(trace = 0, fnscale = -1, maxit = maxit,reltol=1e-6)), silent = TRUE)
        if("try-error" %in% class(q)){ new.factor_scores <- factor_scores; new.alpha <- alpha
        }else{
          if(iter > 1 && m.cur.logfunc > q$value){if(trace)
            cat("Optimization of m and alpha did not improve on iteration step",iter,"\n");
            new.factor_scores <- factor_scores; new.alpha <- alpha;
          }else{
            if(trace) cat("Variational parameters m and alpha updated","\n")
            new.factor_scores <- matrix(q$par[1:(n.factors*n.s)],n.s,n.factors); x2 <- q$par[-(1:(n.s * n.factors))];
            new.alpha <- x2[1:n.s]; x2 <- x2[-(1:n.s)];
            if(q$convergence != 0) { if(trace) cat("Optimization of m and alpha did not converge on iteration step", iter,"\n") }
          }
        }

        q_m <- list(value = m_f(c(new.factor_scores,new.alpha), x = c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), s = c(new.sigma), pi = new.pi,eta = new.eta,d=new.dispersion))
        m.new.logfunc <- q_m$value
        m.cur.logfunc <- m.new.logfunc

        q_s <- list(value = s_f(c(new.sigma), x = c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), pi = new.pi,eta = new.eta,d=new.dispersion))
        s.new.logfunc <- q_s$value
        s.cur.logfunc <- s.new.logfunc

        delta.pi <- abs(new.pi-pi)
        sigma <- new.sigma
        factor_scores <- new.factor_scores
        alpha <- new.alpha
        pi <- new.pi
        p.iter <- p.iter+1
      }

      q_b <- list(value = beta_f(c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), s = c(new.sigma), pi = new.pi,eta = new.eta,d=new.dispersion))
      b.new.logfunc <- q_b$value
      b.cur.logfunc <- b.new.logfunc

      if(family=="negative.binomial"){
        q_d <- list(value = dispersion_f(c(new.dispersion), x=c(new.factor_coefs_j,new.factor_coefs_0,new.gamma),v = c(new.factor_scores,new.alpha), s = c(new.sigma), pi = new.pi,eta = new.eta))
        d.new.logfunc <- q_d$value
        d.cur.logfunc <- d.new.logfunc
      }

      ## Take values of VLB to define stopping rule
      q <- list(value = VLB(c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), s = c(new.sigma),pi = new.pi,eta = new.eta,d=new.dispersion))
      new.VLB <- q$value
      diff=abs(new.VLB-cur.VLB)
      ratio <- abs(new.VLB/cur.VLB);
      if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, ". Difference in VLB:",diff,"\n")
      cur.VLB <- new.VLB

      q <- list(value = log_base(c(new.factor_coefs_j,new.factor_coefs_0,new.gamma), v = c(new.factor_scores,new.alpha), s = c(new.sigma),pi = new.pi,eta = new.eta,d=new.dispersion))
      cur.log<- q$value

      eta <- new.eta;
      alpha <- new.alpha;
      factor_coefs_0 <-new.factor_coefs_0;
      factor_coefs_j <- new.factor_coefs_j;
      gamma <- new.gamma
      factor_scores <- new.factor_scores;
      sigma <- new.sigma;
      pi <- new.pi;
      if(family=="negative.binomial"){dispersion <- new.dispersion}

      iter <- iter + 1
    }

    if(iter > 99){
      print(paste(family,"ZIPPCA Not converging!"))
    }

    if(family=="negative.binomial"){
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(alpha,n.s,n.f) + factor_scores %*% t(factor_coefs_j)
      dispersion.mat <- matrix(dispersion,n.s,n.f,byrow=TRUE)
      e.mat <- ll + 0.5 * (sigma) %*% t(factor_coefs_j^2)
      e.mat2 <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y  + matrix(alpha,n.s,n.f)+ factor_scores %*% t(factor_coefs_j) - 0.5 * (sigma) %*% t(factor_coefs_j^2)
      e.mat3 <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y  + factor_scores %*% t(factor_coefs_j) + 0.5 * (sigma) %*% t(factor_coefs_j^2)
      
      mu_p <- (exp(e.mat)+dispersion.mat)/(1+dispersion.mat/exp(e.mat2))
      Q_nb <- mu_p/rowSums(mu_p)
      Q_nb_z <- (1-new.pi)*mu_p/(rowSums((1-new.pi)*mu_p))
      mu_z <- (1-new.pi)*exp(e.mat3)
    }
    if(family=="poisson"){
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y + matrix(alpha,n.s,n.f) + factor_scores %*% t(factor_coefs_j)
      e.mat <- ll + 0.5 * (sigma) %*% t(factor_coefs_j^2)
      Q_poi <- exp(e.mat)/rowSums(exp(e.mat))
      Q_poi_z <- (1-new.pi)*exp(e.mat)/(rowSums((1-new.pi)*exp(e.mat)))
      e.mat2 <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y  + factor_scores %*% t(factor_coefs_j) + 0.5 * (sigma) %*% t(factor_coefs_j^2)
      mu_z <- (1-new.pi)*exp(e.mat2)
      
    }
    
    ## print the output
    out.list$VLB <- cur.VLB
    out.list$CLL <- cur.log
    out.list$iter=iter-1
    out.list$lvs$pi <- pi
    out.list$lvs$factor_scores <- factor_scores
    out.list$lvs$sigma <- sigma
    out.list$params$eta <- eta
    out.list$params$factor_coefs_j <- factor_coefs_j
    out.list$params$factor_coefs_0 <- factor_coefs_0
    out.list$params$gamma <- gamma
    out.list$params$alpha <- alpha
    if(family=="negative.binomial"){
      out.list$params$dispersion <- dispersion
      out.list$mu <- exp(e.mat3)
      out.list$Q <- Q_nb
      out.list$muz <- mu_z
      out.list$Qz <- Q_nb_z
    }
    if(family=="poisson"){
      out.list$Q <- Q_poi
      out.list$mu <- exp(e.mat2)
      out.list$muz <- mu_z
      out.list$Qz <- Q_poi_z
    }
    return(out.list)
  }

  if(rank==FALSE & family == "negative.binomial"){
    re <- ZIPNVA(X,V,family = "negative.binomial",n.factors,trace,maxit)
  }else if(rank==FALSE & family == "poisson"){
    re <- ZIPNVA(X,V,family = "poisson",n.factors,trace,maxit)
  }else if(rank==TRUE){
    if (parallel){
      #cl <- parallel::makeCluster(detectCores(logical = FALSE))
      cl <- parallel::makeCluster(getOption("cl.cores", 4))
      doParallel::registerDoParallel(cl)
    }
    out.list <- list()
    p <- ncol(X)
    n <- nrow(X)
    r <- 5
    fold <- 5
    beta <- list()
    beta0 <- list()
    f <- list()
    pi <- list()
    gamma <- list()
    alpha <- list()
    sigma <- list()
    eta <- list()
    Q <- list()
    mu <- list()
    Qz <- list()
    muz <- list()
    L <- rep(0,r); lob <- rep(0,r)
    iter <- rep(0,r)
    hic <- rep(0,r);
    
    if(family == "negative.binomial"){
      dispersion <- list()
      if (parallel){
        Mres <- foreach::foreach(w=1:r) %dopar% {
          re <- ZIPNVA(X,V,family = "negative.binomial",n.factors=w,trace,maxit)
          re
        }
        for(w in 1:r){
          L[w] <- Mres[[w]]$VLB
          lob[w] <-  Mres[[w]]$CLL
          iter[w] <- Mres[[w]]$iter
          beta[[w]] <- Mres[[w]]$params$factor_coefs_j
          beta0[[w]] <- Mres[[w]]$params$factor_coefs_0
          gamma[[w]] <- Mres[[w]]$params$gamma
          alpha[[w]] <- Mres[[w]]$params$alpha
          dispersion[[w]] <- Mres[[w]]$params$dispersion
          f[[w]] <- Mres[[w]]$lvs$factor_scores
          sigma[[w]] <- Mres[[w]]$lvs$sigma
          pi[[w]] <- Mres[[w]]$lvs$pi
          eta[[w]] <- Mres[[w]]$params$eta
          mu[[w]] <- Mres[[w]]$mu
          Qz[[w]] <- Mres[[w]]$Qz
          Q[[w]] <- Mres[[w]]$Q
          muz[[w]] <- Mres[[w]]$muz
          hic[w] <- -2*lob[w]+log(n)*(w*p-w^2)+2*w*n
        }
      }else{
        for(w in 1:r){
          re <- ZIPNVA(X,V,family = "negative.binomial",n.factors=w,trace,maxit)
          L[w] <- re$VLB
          lob[w] <- re$CLL
          iter[w] <- re$iter
          beta[[w]] <- re$params$factor_coefs_j
          beta0[[w]] <- re$params$factor_coefs_0
          gamma[[w]] <- re$params$gamma
          alpha[[w]] <- re$params$alpha
          dispersion[[w]] <- re$params$dispersion
          sigma[[w]] <- re$lvs$sigma
          f[[w]] <- re$lvs$factor_scores
          eta[[w]] <- re$params$eta
          pi[[w]] <- re$lvs$pi
          Q[[w]] <- re$Q
          mu[[w]] <- re$mu
          Qz[[w]] <- re$Qz
          muz[[w]] <- re$muz
          hic[w] <- -2*lob[w]+log(n)*(w*p-w^2)+2*w*n
        }
      }}else if (family == "poisson"){
        if (parallel){
          Mres <- foreach::foreach(w=1:r) %dopar% {
            re <- ZIPNVA(X,V,family = "poisson",n.factors=w,trace,maxit)
            re
          }
          for(w in 1:r){
            L[w] <- Mres[[w]]$VLB
            lob[w] <-  Mres[[w]]$CLL
            iter[w] <- Mres[[w]]$iter
            beta[[w]] <- Mres[[w]]$params$factor_coefs_j
            beta0[[w]] <- Mres[[w]]$params$factor_coefs_0
            gamma[[w]] <- Mres[[w]]$params$gamma
            alpha[[w]] <- Mres[[w]]$params$alpha
            f[[w]] <- Mres[[w]]$lvs$factor_scores
            sigma[[w]] <- Mres[[w]]$lvs$sigma
            pi[[w]] <- Mres[[w]]$lvs$pi
            eta[[w]] <- re$params$eta
            mu[[w]] <- Mres[[w]]$mu
            Qz[[w]] <- Mres[[w]]$Qz
            Q[[w]] <- Mres[[w]]$Q
            muz[[w]] <- Mres[[w]]$muz
            hic[w] <- -2*lob[w]+log(n)*(w*p-w^2)+2*w*n
          }
        }else{
          for(w in 1:r){
            re <- ZIPNVA(X,V,family = "poisson",n.factors=w,trace,maxit)
            L[w] <- re$VLB
            lob[w] <- re$CLL
            iter[w] <- re$iter
            beta[[w]] <- re$params$factor_coefs_j
            beta0[[w]] <- re$params$factor_coefs_0
            gamma[[w]] <- re$params$gamma
            alpha[[w]] <- re$params$alpha
            sigma[[w]] <- re$lvs$sigma
            f[[w]] <- re$lvs$factor_scores
            eta[[w]] <- re$params$eta
            pi[[w]] <- re$lvs$pi
            Q[[w]] <- re$Q
            mu[[w]] <- re$mu
            Qz[[w]] <- re$Qz
            muz[[w]] <- re$muz
            hic[w] <- -2*lob[w]+log(n)*(w*p-w^2)+2*w*n
          }
        }}
    
    out.list$hic <- which.min(hic)
    out.list$VLB <- L[[which.min(hic)]]
    out.list$CLL <- lob[[which.min(hic)]]
    out.list$iter <- iter[[which.min(hic)]]
    out.list$lvs$pi <- pi[[which.min(hic)]]
    out.list$lvs$factor_scores <- f[[which.min(hic)]]
    out.list$lvs$sigma <- sigma[[which.min(hic)]]
    out.list$params$factor_coefs_j <- beta[[which.min(hic)]]
    out.list$params$factor_coefs_0 <- beta0[[which.min(hic)]]
    out.list$params$alpha <- alpha[[which.min(hic)]]
    out.list$params$gamma <- gamma[[which.min(hic)]]
    out.list$params$eta <- eta[[which.min(hic)]]
    
    if(family=="negative.binomial"){
      out.list$params$dispersion <- dispersion[[which.min(hic)]]
    }
    
    if (parallel){
      parallel::stopCluster(cl = cl)
    }
    return(out.list)
    
  }
}


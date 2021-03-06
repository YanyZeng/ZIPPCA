# ZIPPCA
A general framework, Zero-Inflated Probabilistic PCA (ZIPPCA), for dimension reduction and data ordination of multivariate abundance data, and propose an efficient variational approximation method for estimation,inference, and prediction.

# Installation
```r
install.packages("devtools")  
devtools::install_github("YanyZeng/ZIPPCA")  
library(ZIPPCA) 
```

# Description
Dimension reduction and ordination analysis are often applied to multivariate abundance data to visualize broad trends of how similar or different microbial communities are. However, this process is complicated by several statistical challenges. In particular, microbiome data produced by high-throughput sequencing are count-valued, correlated, high-dimensional, and over-dispersed with excess zeros. To address these challenges, we introduce a probabilistic framework by extending the standard factor analysis model, and propose a variational approximation algorithm for estimation, inference, and data ordination. We demonstrate the superior performance of the proposed method over existing ones on a number of simulated examples and a gut microbiome dataset.

# Usage
```r
ZIPPCApn <- function(X, V=NULL, family = "negative.binomial", n.factors=2, rank=FALSE, trace = FALSE, maxit = 100, parallel=TRUE)
```
* X: matrix of observations.
* V: vector of the sample covariate.
* family: distribution of models. Two options are "poisson" and "negative.binomial". Defaults to "negative.binomial".
* n.factors: the rank or number of factors, after dimensional reduction. Defaults to 2.
* rank: logical, if TRUE, the rank or number of factors, is chosen from 1 to 5 by HIC (hybrid information criterion). Defaults to FALSE.
* trace: logical, defaults to FALSE. if TRUE each current iteration step information will be printed.
* maxit: maximum number of iterations within optim and constrOptim function, defaults to 100.
* parallel: logical, if TRUE, use parallel toolbox to accelerate.

# examples
```
 n.n = 100
 n.m = 50
 n.factors = 2
 set.seed(37)
 f <- matrix(0,nrow = n.n, ncol = n.factors)
 for(i in 1:n.n){
   f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
 }
 betaj <- matrix(0,nrow = n.m, ncol = n.factors)
 for(j in 1:n.m){
   betaj[j,] <- rnorm(n.factors, mean =0, sd = 1)
 }
 beta0 <- rep(1,n.m)
 alpha <- rep(0,n.n)
 l <- matrix(alpha,n.n,n.m) +matrix(beta0,n.n,n.m,byrow=TRUE) + f %*% t(betaj)
 lambda <- exp(l)
 eta_j <- runif(n.m, 0,1)
 z <- matrix(0,n.n,n.m)
 for(i in 1:n.n){
   z[i,] <- rbinom(n.m, size=1, prob=eta_j)
 }
 X <- matrix(0,n.n,n.m,byrow = TRUE)
 for(i in 1:n.n){
   for(j in 1:n.m){
     X[i,j] <- rnbinom(n=1,size=10,mu=lambda[i,j])
   }
 }
 X[z==1]=0
 zerorow <- which(rowSums(X)==0)
 if(length(zerorow) >0 ){
   X <- X[-zerorow,];f <- f[-zerorow,]
 }
 zerocol <- which(colSums(X)==0)
 if(length(zerocol) >0 ){
   X <- X[,-zerocol];betaj <- t(t(betaj)[,-zerocol])
 }
 result <- ZIPPCA::ZIPPCApn(X)
 f_coordinates <- result$lvs$factor_scores
 }
 ```

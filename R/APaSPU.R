
#' The function for adaptive sum of powered score (aSPU) test
#'
#' It gives the simulation-based p-values of the sum of powered score (SPU) with different power and aSPU test.
#' @param betas the vector (length = p) of maximum likelihood estimates; p is the number of paramters to be tested.
#' @param vcov_mat the p by p variance-covariance matrix of MLE.
#' @param pow the power vector for aSPU tests.
#' @param B number of simulations
#' @return  simulation-based p-values of the sum of powered score (SPU) with different power and aSPU test
#' @details the suggested power term is c(1:8, Inf)
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#' @references
#' Zhiyuan Xu, Wei Pan (2014) Approximate score-based testing with application to multivariate
#' trait association analysis. Genetic Epidemiology. 39(6): 469-479.
#'
#' Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014) Testing for association with multiple
#' traits in generalized estimation equations, with application to neuroimaging data. Neuroimage.
#' 96:309-25
#' @examples
#' library(lme4)
#' data("exdat_GLMM")
#' fit1<-glmer(Y~.-ID + (1|ID),family="binomial",
#' data=exdat_GLMM,control=glmerControl(optimizer="bobyqa"))
#' betas<-fixef(fit1)[-1]
#' Cov<-vcov(fit1)[-1,-1]
#' APaSPU(betas,Cov,pow=c(1:8,Inf),B=1000)
#' @export
#' @importFrom MASS mvrnorm ginv
#' @importFrom lme4 glmer
APaSPU <- function(betas, vcov_mat, pow = c(1:8,Inf), B){

  V = ginv(as.matrix(vcov_mat))
  U = V %*% betas
  K = length(betas)
  gamma2 = pow
  gamma1 = 1
  weight = F
  spu=spuval(U=U,V=V, gamma1=gamma1,gamma2=gamma2,K,weight=weight)
  p=rep(NA,length(gamma1)*length(gamma2))
  ###### NULL Distribution
  T0s=matrix(NA,nrow=B,ncol=length(gamma1)*length(gamma2))
  u.null<-mvrnorm(B,rep(0,length(U)),V)
  for (b in 1:B) T0s[b,]=spuval(U=u.null[b,],V=V,gamma1=gamma1,gamma2=gamma2,K=K,weight=weight)
  for (g in 1:ncol(T0s))  p[g]<-(1+sum(abs(spu[g])<abs(T0s[,g])))/(1+B)

  ##### p-values
  P0s<-PermPvs(T0s,B=B)
  minp0<-apply(P0s,1,min)
  Paspu<-(sum(minp0<=min(p))+1)/(B+1)
  pvs<-data.frame(matrix(c(p,Paspu),nrow=1))

  names(pvs)=c(paste0("SPU.",pow),"aSPU")
  return(pvs)
}









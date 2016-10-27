
#' The function gee.approxscore to construct the approximate score vector from GEE model
#' @param dat, the data table for analysis
#' @param y, response name
#' @param par, variables of interest for testing
#' @param cov, covariates needed to be adjusted
#' @param groupID, group ID name
#' @param corstr, working correlation structure for GEE
#' @param family type, common ones are "gaussian" and "binomial".
#' @return Score vector, variance-covariance matrix of Score vector and p-value by Score test
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#'
#' @importFrom gee gee
#' @importFrom stats as.formula pchisq
gee.approxscore = function(dat,y,par,cov, groupID, corstr="independence",family){

  if(is.null(cov)) form<-as.formula(paste(y,"~0+",paste(par,collapse="+"),sep="")) else
  form<-as.formula(paste(y,"~0+",paste(par,collapse="+"),
                         "+",paste(cov,collapse="+"),sep=""))
  id = dat[,groupID]
  fit1 = gee(form, id=id,data=dat,family=family,corstr=corstr)
  p.num = length(par)     # number of par
  betas<-fit1$coef[1:p.num]
  Rcov<-fit1$robust.variance[1:p.num,1:p.num]
  Mcov<-fit1$naive.variance[1:p.num,1:p.num]
  VmInv = solve(Mcov)
  U = VmInv %*% betas
  CovS=VmInv %*% Rcov %*% VmInv


  ss<-t(U)%*%ginv(CovS)%*%U
  pval<-1-pchisq(ss,p.num)

  return(list(U=U, CovS = CovS,pval=pval))
}

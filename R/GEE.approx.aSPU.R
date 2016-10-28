#' The function GEE.approx.aSPU compute the p-values of aSPU test
#' based on approximated Score vectors from GEE
#' @param dat, the data table for analysis
#' @param y, response name
#' @param par, variables of interest for testing
#' @param cov, covariates needed to be adjusted
#' @param groupID, group ID name
#' @param corstr, working correlation structure for GEE
#' @param family type, common ones are "gaussian" and "binomial".
#' @param pow, power integer candidates
#' @param B, simulation number to calculate P-values
#' @return P-values for SPU and aSPU test based on approximated Score vector from GEE
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#' @references
#' Zhiyuan Xu, Wei Pan (2014) Approximate score-based testing with application to multivariate
#' trait association analysis. Genetic Epidemiology. 39(6): 469-479.
#'
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014) A powerful and adaptive
#' association test for rare variants, Genetics, 197(4), 1081-95
#'
#' Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014) Testing for association with multiple
#' traits in generalized estimation equations, with application to neuroimaging data. Neuroimage.
#' 96:309-25
#' @examples
#'
#' data("exdat_GLMM")
#' gee.aSPU <- GEE.approx.aSPU(exdat_GLMM, "Y", paste0("X.",1:10), cov = NULL,
#' groupID = "ID",corstr="independence",family="binomial", B = 1000 )
#' gee.aSPU
#' @export
#
GEE.approx.aSPU = function(dat,y,par,cov, groupID, corstr="independence",family="binomial", pow=c(1:8, Inf), B){

  gee1 = gee.approxscore(dat = dat,y=y, par=par,cov=cov, groupID = groupID, corstr=corstr, family=family)
  U = as.vector(gee1$U)
  CovU = as.matrix(gee1$CovS)
  pout <- SPU.mulT(U=U,V=CovU,pow=pow,B=B)
  return(pout)
}

#' PermPvs is the function to calcualte p-values for SPU test by permutation
#'
#' It gives the simulation-based p-values of the sum of powered score (SPU) with different power and aSPU test.
#' @param T0s, Null statistics matrix, the SPU statistic is matrix, rows for gamma2 and columns for gamma1, where gamma2 is for traits and gamma1 is for SNPs
#' @param B number of simulations
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#' @return  simulation-based p-values of the sum of powered score (SPU) with different power
#' @details the suggested power term is c(1:8, Inf)
#' @export


PermPvs<-function(T0s,B){
  P0s=apply(T0s,2,function(z)( B-rank(abs(z))+1)/(B) )
  return(P0s)
}

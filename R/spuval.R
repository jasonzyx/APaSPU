#' The function is to calculate the SPU statistics
#'
#' @param U, the Score vector
#' @param V, the variance-covariance matrix of U
#' @param gamma1, power candidates
#' @param gamma2, power candidates
#' @param K, number of traits
#' @param weight, TRUE or FALSE, default is FALSE
#' @return  the SPU statistics
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014) A powerful and adaptive
#' association test for rare variants, Genetics, 197(4), 1081-95
#'
#' Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014) Testing for association with multiple
#' traits in generalized estimation equations, with application to neuroimaging data. Neuroimage.
#' 96:309-25
#'


spuval<-function(U,V,gamma1,gamma2, K,weight=F){
  if (weight==F) U2=matrix(U,nrow=K,byrow=T)  ## rows for traits and columns for snps
  if (weight==T) U2=matrix(U/sqrt(diag(V)),nrow=K,byrow=T)
  spumat<-matrix(NA,ncol=length(gamma1),nrow=length(gamma2))
  for (g1 in 1:length(gamma1)){
    if (gamma1[g1]<Inf) spu.g1<-apply(U2,1,function(z)sum(z^gamma1[g1])) else spu.g1<-apply(U2,1,function(z)max(abs(z)))
    for (g2 in 1:length(gamma2)){
      if (gamma2[g2]<Inf) spumat[g2,g1]=sum(spu.g1^(gamma2[g2]/gamma1[g1])) else spumat[g2,g1]=max(abs(spu.g1))
    }
  }
  spu0=c(t(spumat))
  spu0
}

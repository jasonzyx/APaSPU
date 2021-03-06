% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spuval.R
\name{spuval}
\alias{spuval}
\title{The function is to calculate the SPU statistics}
\usage{
spuval(U, V, gamma1, gamma2, K, weight = F)
}
\arguments{
\item{U, }{the Score vector}

\item{V, }{the variance-covariance matrix of U}

\item{gamma1, }{power candidates}

\item{gamma2, }{power candidates}

\item{K, }{number of traits}

\item{weight, }{TRUE or FALSE, default is FALSE}
}
\value{
the SPU statistics
}
\description{
The function is to calculate the SPU statistics
}
\author{
Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
}
\references{
Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014) A powerful and adaptive
association test for rare variants, Genetics, 197(4), 1081-95

Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014) Testing for association with multiple
traits in generalized estimation equations, with application to neuroimaging data. Neuroimage.
96:309-25
}


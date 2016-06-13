\name{ComLikMulCom{mvtnorm, multcomp,MASS}}
\alias{ComLikMulCom}

\title{Multiple Comparisons using Composite Likelihood
}
\description{Multiple hypothesis testing \eqn{H_0 :C\beta = 0}, based on composite likelihood inference for clustered data, with common covariance structure across the clusters. The critical value is based on multivariate normal quantile, using the R-package \code{multcomp} by Hothorn, et.al.
}
\usage{
ComLikMulCom(data, yno, idno, y.type = c("normal", "probit", "quadratic exponential")[1], type = c("many_one", "all_pair")[1], f = 1)
}

\arguments{
  \item{data}{data set in long format containing response, covariates and one variable that
  identifies the cluster}
  \item{yno}{the response column number}
  \item{idno}{column number of the variable that identify the clusters}
  \item{y.type}{the model of the response variable, it could be normal, probit, or quadratic exponential}
  \item{type}{type of the test, it could be many-to-one or all-pair comparisons}
  \item{f}{the baseline variable for many to one type}
}


\details{The current functions handle the clustered normal and probit data with equal size. For quadratic exponential data, clusters can have different sizes.

}
\value{
\item{m}{number of individuals in each cluster, if equal}
\item{n}{number of clusters}
\item{p}{number of covariates}
\item{Cov_beta}{covariance matrix of the covariates}
\item{T_stat}{the test statistic}
\item{normal_quantile}{the simultaneous normal quantile }
\item{MNQ}{results of individual null hypotheses based on composite likelihood}
\item{Bon}{results of individual null hypotheses based on Bonferroni approach}
\item{Sidak}{results of individual null hypotheses based on Dunn-Sidak approach}
\item{Scheffe}{results of individual null hypotheses based on Scheffe method}
}

\references{
Azadbakhsh,M.,  Gao,X., Jankowski,H. (2016)
{Multiple Comparisons Using Composite Likelihood in Clustered Data}. Submitted.

Hothorn, T., Bretz, F., Westfall, P., Heiberger, R. M., Schutzenmeister, A.(2010)
{ multcomp: Simultaneous Inference for General Linear Hypotheses}, R package version 1.1-7, http://CRAN.R-project.org/package=multcomp
}


\author{
Mahdis Azadbakhsh,
Xin Gao,
Hanna Jankowski.
}
\examples{
### artificial data set (normal)
m <- 3
n <- 20
p <- 5

v <- sample(seq(-4,4,0.05),360,replace= TRUE)
cx <- matrix(v, nrow=(n*m))

ind <- rep(0,n*m)
for(i in 1: n) {ind[(m*(i-1)+1):(m*i)]<-   i}

cxx <- cbind(cx,ind)
colnames(cxx) <- c("res", "c1","c2","c3","c4", "c5", "ind")

idno <- which( colnames(cxx)=="ind" )
yno <- which( colnames(cxx)=="res" )

ComLikMulCom(cxx, yno , idno , y.type ="normal" ,  type = "many_one", f=3)
ComLikMulCom(cxx, yno , idno , y.type ="normal" ,  type = "all_pair")



### probit


v <- sample(seq(-4,4,0.05),300,replace= TRUE)
y <- sample(c(0,1),60,replace= TRUE)
cx <- matrix(v, nrow=60)

colnames(cx) <- c("c1","c2","c3","c4", "c5")


id <- rep(0,n*m)
for(i in 1: n) {id[(m*(i-1)+1):(m*i)]<-   i}


cxx <- as.data.frame(cbind(y,cx,id))
idno <- which( colnames(cxx)=="id" )
yno <- which( colnames(cxx)=="y" )


ComLikMulCom(cxx, yno , idno , y.type ="probit" ,  type = "all_pair")



### Normal
m <- 3
n <- 40
p <- 5
require(MASS)
library(MASS)
library(mvtnorm)
beta <-   c(0,2,0,0,1)
sigma  <- matrix(c(0.8, 0.4, 0.4, 0.4, 0.8, 0.4, 0.4, 0.4, 0.8), nrow=3)

desx      <-    array(dim=c(m,p,n))
y         <-    matrix(0,nrow=m,ncol=n)
mu        <-    matrix(0,nrow=m,ncol=n)

for (k in 1:n){
  for(j in 1:p){
    desx[,j,k] <-   rnorm(m,1,3.5)
  }
  mu[,k]      <-   tcrossprod(desx[,,k], t(beta))
  y[,k]        <-   mvrnorm(1,mu[,k], sigma, tol = 1e-6, empirical = FALSE)
}
sss    <-    array(dim=c(m,(p+1),n))


for(i in 1:n){
  ind <- rep(i,m)
  sss[,,i]  <- cbind(ind,desx[,,i])
}

yyy         <-   rep(0,120)

for( i in 1:40){
  yyy[((i-1)*3+1):(i*3)]   <-   y[,i]
}

xxx       <-   matrix(0,nrow=120, ncol=6)

for( i in 1:40){
  xxx[((i-1)*3+1):(i*3),] 	<- 	sss[,,i]
}
dat <- cbind(yyy,xxx)
colnames(dat) <- c("response","id", "c1","c2","c3","c4", "c5")

idno <- which( colnames(dat)=="id" )
yno <- which( colnames(dat)=="response" )
ComLikMulCom(dat, yno , idno , y.type ="normal" ,  type = "many_one", f=3)
ComLikMulCom(dat, yno , idno , y.type ="normal" ,  type = "all_pair")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multiple comparisons }
\keyword{ composite likelihood }% __ONLY ONE__ keyword per line



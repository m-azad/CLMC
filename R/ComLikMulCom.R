


ComLikMulCom <- function(data, yno , idno , y.type =c("normal", "probit", "quadratic exponential")[1] ,  type = c("many_one", "all_pair")[1], f=1){
  
  
  require(mvtnorm)
  
  N   <-  length((data[,idno]))               #total number of measurements
  n   <-  length(unique(data[,idno]))         #number of individuals
  p   <-  ncol(data)-2                        # number of factors
  m   <-  N/n
  
  if("many_one" %in% type){
    
    mt1     <- matrix(0, nrow=(p-1), ncol=p)
    mt1[,f] <- -1
    
    for (i in 1:(p-1)){
      if(mt1[i,i]==0) {        
        mt1[i,i]<- 1}
      else if (mt1[i,i]==-1){
        k <- i
{break}}
    }

for (j in k:(p-1)){
  mt1[j,(j+1)]<- 1}

C     <-    mt1

  }


if("all_pair" %in% type){
  
  C     <-   matrix(0,nrow=p,ncol=p*(p-1)/2)
  comb   <- 	combn(p,2)
  C[cbind(comb[1,],1:(p*(p-1)/2))] <- 1
  C[cbind(comb[2,],1:(p*(p-1)/2))] <- -1
  C     <-    t(C)
}


CLcov <- function(data, yno , idno, N, n, p , y.type = c("normal", "probit", "quadratic exponential")[1] ,  type = c("many_one", "all_pair")[1]){
  
   
  data <- data.frame(data)
  
  yy <- as.matrix(data[,yno])
  xx <- as.matrix(data[,-c(yno,idno)])
  
  
  
  if("normal" %in% y.type){
    ###    Estimating beta and sigma from marginal normal log likelihood
    
    sum1     <-   t(xx)%*%xx
    sum2     <-   t(xx)%*%yy
    beta_hat <-   solve(sum1)%*%sum2  
    
    
    yhat      <-   matrix(0,nrow=m,ncol=n)
    e         <-   matrix(0,nrow=m,ncol=n)
    
    for(i in 1: n){
      index <- which(data[, idno]==i)
      d <- as.matrix(data[index,])
      e[,i] <- d[,yno] - d[,-c(idno,yno)]%*%beta_hat
    }
    
    
    sigma_hat <-   matrix(0, nrow=m, ncol=m)
    
    
    sii     <-   sum(apply(e,1,function(x)sum(x^2)))/(n*m-p)  
    
    E         <-   rep(0,n)
    for (i in 1:n){
      for( j in 1:(m-1)){
        E[i] 	<- 	E[i] + sum(e[j,i]*e[(j+1):m,i])
      }  
    }
    r 		<- 	sum(E)/((n-p)*choose(m,2)*sii)
    
    Sigma_hat  <- 	array(rep(r,m^2),dim=c(m,m))
    Sigma_hat  <- 	Sigma_hat+ (1-r)*diag(m)
    Sigma_hat  <- 	sii*Sigma_hat 		                 
    
    ###  	computing covariance of coefficients
    
    c1		<-	matrix(0,nrow=p,ncol=p)
    c2		<- 	matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      
      index <- which(data[, idno]==i)
      d <- as.matrix(data[index,])
      
      #d <- as.matrix(data[data$id==i,])
      c1	<-	c1 +(t(d[,-c(idno,yno)])%*%d[,-c(idno,yno)])
      c2	<-	c2 +t(d[,-c(idno,yno)])%*%Sigma_hat%*%d[,-c(idno,yno)]
    }
    
    Cov_beta  <- 	solve(c1)%*%c2 %*%t(solve(c1)) 
    
    
  } 
  
  
  
  if("probit" %in% y.type){
    
    fun1_b <- function(s,x,y){  
      
      ps   <- as.numeric(pnorm(s,0,1,lower.tail=TRUE))
      pms  <- as.numeric(pnorm(-s,0,1,lower.tail=TRUE))
      ds   <- as.numeric(dnorm(s,0,1))  
      if  (abs(s)<= 5) {z1 <- ((y-ps)*(ds*x))/(ps*pms)
      } else if (s>5 & y==1) {z1 <- ds*x/ps
      } else if (s>5 & y==0) {z1 <- -s*x
      } else if (s< -5 & y==1) {z1 <- -s*x
      } else if (s< -5 & y==0) {z1 <- -ds*x/pms}
      return(z1)
    }
    
    
    fun2_b <- function(s,x,y){  
      
      ps   <- as.numeric(pnorm(s,0,1,lower.tail=TRUE))
      pms  <- as.numeric(pnorm(-s,0,1,lower.tail=TRUE))
      ds   <- as.numeric(dnorm(s,0,1))
      if  (y==1 & s>-5 ) {z2 <-as.numeric((-ds^2/ps^2) - (ds*s/ps))*x%*%t(x)
      }else if (y==1 & s< -5) {z2 <- 0
      }else if (y==0 & s< 5 ) {z2 <-as.numeric((-ds^2/pms^2) + (ds*s/pms))*x%*%t(x)
      }else if (y==0 & s> 5) {z2 <- 0}
      return(z2)
    }
    ###Initial value
    
    fit       <-    glm(yy~ xx ,family=binomial(link = "probit"), data , 
                        control = list(epsilon = 1e-6, maxit = 30, trace = FALSE))
    
    beta_hat  <-    as.vector(fit$coefficients[-1])
    
    
    ###Newton Raphson for beta
    
    sig     <-   1
    
    for (k in 1:1000){ 
      
      beta_hat_n   <-  beta_hat 
      sc_b         <-  rep(0,p)
      hs_b         <-  matrix(0,p,p)
      
      for(i in 1:nrow(xx)){
        
        X    <-  xx[i,]%*%beta_hat_n/sig
        
        sc_b <- sc_b + fun1_b(X,xx[i,],yy[i,])
        hs_b <- hs_b + fun2_b(X,xx[i,],yy[i,])
      }
      
      beta_hat <- ifelse(sig*abs(solve(hs_b)%*%sc_b) < 10^(-5),beta_hat_n,beta_hat_n - sig*solve(hs_b)%*%sc_b)
      
      if (all(beta_hat== beta_hat_n)){break}
    } 
    
    v_s  <- matrix(0,p,p)
    H    <- matrix(0,p,p)
    
    
    for ( i in 1:n){
      
      index <- which(data[, idno]==i)
      d <- as.matrix(data[index,])
      
      dmu   <-   as.vector(dnorm(d[,-c(idno,yno)]%*%beta_hat))*d[,-c(idno,yno)]
      mu    <-   pnorm(d[,-c(idno,yno)]%*%beta_hat)
      V     <-   diag(as.vector(pnorm(d[,-c(idno,yno)]%*%beta_hat))*as.vector(pnorm(-d[,-c(idno,yno)]%*%beta_hat)))
      sij   <-   (d[,yno]-mu)%*%t(d[,yno]-mu)
      
      v_s   <-   v_s + t(dmu)%*%solve(V)%*%sij%*%solve(V)%*%dmu
      H     <-   H + t(dmu)%*%solve(V)%*%dmu
      
    }
    
    Cov_beta  <-  solve(H)%*%v_s%*%solve(H)
    
  }
  
  if("quadratic exponential" %in% y.type){
    ## creating the column of dependence factor (w)
    
    pat      <-  array(0,dim=c(n,2))  
    pat[,2]  <-  table(data[,idno])
    ct       <-  1
    breakdown<-  vector("list", n)
    
    
    for (i in 1:(N-1)){
      breakdown[[ct]]  <- c(breakdown[[ct]],data[i,yno])
      if (data[,idno][i+1]>data[,idno][i]){ct<-ct+1}
    }
    for (i in N:N){ 
      breakdown[[ct]]  <- c(breakdown[[ct]],data[i,yno])
    }
    for ( i in 1: n){
      brk <- breakdown[[i]]
      brk <- unlist(brk)
      pat[i,1]         <- sum(brk)
    }
    
    
    ss     <-  pat[,1]
    ff     <-  pat[,2]
    data$z <-  rep(ss,ff)
    data$m <-  rep(ff,ff)
    data$w <-  ifelse(data[,yno]==1,  -(data$m-2*data$z+1),-(data$m-2*data$z-1))
    
    
    r     <-  dim(data)[2]
    xdat  <-  data[,-c((r-1),(r-2))] 
    
    factnam  <-  names(xdat[,-c(idno,yno)])
    
    resp     <-  xdat[,yno]  
    
    fmla     <-  as.formula(paste("resp ~ ", paste(factnam, collapse= "+")))
    fit      <-  glm(fmla ,data=xdat ,family=binomial(link=logit))
    beta <-  fit$coefficients[-1]
    
    beta_hat     <-   as.matrix(fit$coefficients[2:(p+1)])
    what     <-   fit$coefficients[(p+2)]
    
    
    X  <- split(data, data[,idno])
    X  <- lapply(X, function(x) { x[c("id")] <- NULL; x })
    X  <- lapply(X, function(x) { x[c("w")] <- NULL; x })
    
    
    scorevec  <-    array(0,dim=c(n,p))
    
    for ( i in 1: n){
      
      xxs      <-    as.data.frame(X[[i]])
      xxs$indc <-    ifelse(xxs[,yno]==1, 1, -1)
      e        <-    rep(0,nrow(xxs))
      no3      <-    which( colnames(xxs)=="indc" )
      
      for(q in 1: nrow(xxs)){
        med <- xxs[,no3]
        e[q] <- ifelse(med[q]==1,exp(as.matrix(xxs[q,-c(yno,no3,(p+2),(p+3))])%*%beta_hat-what*(xxs$m[q]-2*xxs$z[q]+1))
                       ,exp(-as.matrix(xxs[q,-c(yno,no3,(p+2),(p+3))])%*%beta_hat +what*(xxs$m[q]-2*xxs$z[q]-1)) ) 
      }
      
      
      xxj           <- cbind(xxs,e)
      xxj$prob      <- 1/(1+e)
      
      no4 <- which( colnames(xxj)=="z" )
      no5 <- which( colnames(xxj)=="m" )
      no6 <- which( colnames(xxj)=="indc" )
      no7 <- which( colnames(xxj)=="e" )
      no8 <- which( colnames(xxj)=="prob" )
      
      multmat       <- xxj[,-c(yno,no4,no5,no6,no7,no8)]*xxj$prob*xxj$indc
      scorevec[i,]  <- apply(multmat,2,sum)
    }
    
    no1 <- which( colnames(xdat)=="id" )
    no2 <- which( colnames(xdat)=="w" )
    
    xdes     <-  xdat[,-c(yno,no1,no2)]
    xdes     <-  as.matrix(xdes)
    
    
    vmat     <-  var(scorevec)
    hmat     <-  t(xdes*(fit$weights))%*%xdes
    Cov_beta <-  n*solve(hmat)%*%vmat%*%solve(hmat)
    
    std_bet  <- sqrt(diag(Cov_beta))
    
    Z        <- beta_hat/std_bet
    pValue   <- 2*pnorm(-abs(Z))  
    
    
    
  }
  
  res <- list(beta_hat = beta_hat ,Cov_beta = Cov_beta,  T_stat= T,normal_quantile =quantile )
  return(res)
  
}


covcl <- CLcov(data = data, yno = yno, idno = idno, N = N, n = n, p= p , y.type = y.type, type=type)



bstr      <- C%*%covcl$beta_hat
cov_str   <- C%*%covcl$Cov_beta%*%t(C)
se_str    <- sqrt(diag(cov_str))^(-1)
v_mat     <- diag(se_str)
corr_beta <- v_mat %*% cov_str %*% v_mat

T         <-	 bstr*se_str
cov_T 	  <- 	 v_mat%*%cov_str%*%t(v_mat)  


###	computing the m-variate normal quantile to consider it as threshold and performing hypothesis test	

qu <- qmvnorm(0.95, tail ="both.tails", mean = 0,  corr=cov_T)
quantile <- qu$quantile

rm <- abs(T)<quantile
sM     <-    ifelse(rm ==TRUE,"Accept","Reject") 

###Bonferroni

Bon <- 0.05/nrow(C)
qBon <- qnorm(1-Bon/2, lower.tail = TRUE)  

rb <- abs(T)<qBon
sB    <-  ifelse(rb==TRUE,"Accept","Reject")                      

###Dunn-Sidak

dunn <- 1-(1-0.05)^(1/nrow(C))
qDunn <- qnorm(1-dunn/2)  

rd <- abs(T)<qDunn
sD   <-   ifelse(rd==TRUE,"Accept","Reject") 

### Scheffe

qscheffe  <- sqrt(qchisq(0.95,(p-1)))  

rs <- abs(T)< qscheffe
sS     <-  ifelse( rs ==TRUE,"Accept","Reject") 


res <- list(m = m, n = n, p= p, Cov_beta= covcl$Cov_beta , T_stat= T,normal_quantile =quantile , MNQ = sM  , Bon = sB ,  Sidak = sD , Scheffe = sS)
return(res)

}


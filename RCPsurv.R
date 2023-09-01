RCPsurvSIM<-function(seed=NA,n,gamma,beta,alpha1,alpha2,mu,sigma){
  if(is.numeric(seed)){set.seed(seed)}
  
  u1 <-runif(n,0,1)
  X <- rnorm(n,0,1)
  Z <- runif(n,min = 0,max = 4)
  eta<-rnorm(n,mean = mu,sd = sigma)
  
  Ti.star    <-sqrt(-log(u1)/exp(gamma*X + beta*Z + alpha1*I(Z>=eta) + alpha2*I(Z>=eta)*(Z-eta)))
  Ci         <-runif(n=n, min = 0,max = 5)
  cen        <-ifelse(Ti.star<Ci,1,0)
  Yi         <-pmin(Ti.star,Ci)
  
  dat<-cbind(Yi,cen,X,Z)
  dat<-cbind(c(1:n),dat[order(dat[,1]),])
  dat<-as.data.frame(dat)
  names(dat)<-c("id","Yi","cen","X","Z")
  return(dat)
}

truncated_normal_b_moment<-function(order, mu, sigma, b){
  b_h = (b-mu)/sigma
  b_pdf = dnorm(b_h)
  b_cdf = pnorm(b_h)
  if(b_cdf<10^{-100}){f = NA}else{f = b_pdf/b_cdf}
  
  value<-irm2<-irm1<-0
  
  for(r in c(0:order)){
    if(r==0){ir=1}else if(r==1){ir=-f}else{ir=-b_h^{r-1}*f+(r-1)*irm2}
    value = value + choose(order, r)*mu^{order-r}*sigma^r*ir
    irm2 = irm1
    irm1 = ir
  }
  return(value)
}

truncated_normal_a_moments<-function(order_max, mu, sigma, a){
  XX<-numeric(order_max+1)
  for(o in c(0:order_max)){XX[o+1]<-ifelse(o%%2==0,1,-1)*truncated_normal_b_moment(o, -mu, sigma, -a)}
  return(XX)
}

truncated_normal_b_moments<-function(order_max, mu, sigma, b){
  XX<-numeric(order_max+1)
  for(o in c(0:order_max)){XX[o+1]<-truncated_normal_b_moment(o, mu, sigma, b)}
  return(XX)
}

moment_method<-function(m, moments, boundary, LB){
  if(sum(is.finite(moments))==1){x<-w<-rep(0,m)}else{
    N<-min(m, (sum(is.finite(moments))-1)/2)
    repeat{
      if(N==1){print("ERROR: insufficient number of node"); return(NULL);break}
      #Hankel matrix H
      h<-matrix(0,nrow=N+1,ncol=N+1)
      for(i in c(1:(N+1))){
        for(j in c(1:(N+1))){
          h[i,j] = moments[i+j-1]
        }}
      
      #compute R, the lower triangular Cholesky factor of H.
      XER<-try(chol(h),silent = TRUE)
      if(inherits(XER, "try-error")){N<-(N-1)}else{
        r<-XER
        alpha<-numeric(N)
        beta <-numeric(N-1)
        alpha[1]<-r[1,2]/r[1,1]
        for(o in c(2:N)){alpha[o]<-r[o,o+1]/r[o,o]-r[o-1,o]/r[o-1,o-1]}
        for(o in c(1:(N-1))){beta[o]<-r[o+1,o+1]/r[o,o]}
        
        if(N==2){jacobi <- diag(alpha)+rbind(c(0,beta),c(beta,0))}else{
          jacobi <- diag(alpha)+rbind(cbind(0,diag(beta)),0)+rbind(0,cbind(diag(beta),0))}
        jacobi[which(is.na(jacobi))]<-0
        
        #get the eigen-decomposition of the Jacobi matrix.
        eigen_jacobi.test<-eigen(jacobi)
        x.test           <-eigen_jacobi.test$values
        if((sum(x.test<boundary)>=1)&(LB==TRUE)){N<-(N-1)}else 
          if((sum(x.test>boundary)>=1)&(LB==FALSE)){N<-(N-1)}else{
            eigen_jacobi<-eigen_jacobi.test; x<-x.test ; break
          }}}
    x<-c(x,rep(0,m-N))
    eigvec <-eigen_jacobi$vectors
    w<-numeric(m)
    for(o in c(1:N)){w[o]<-moments[1]*eigvec[1,o]^2}
  }
  return(list(node=x[order(x)], weight=w[order(x)]))
}

fi_given_theta<-function(eta, xi.old, lambda, cumhaz, X, Z, cen, Ind){
  H  <-cbind(X,Z,Ind,Ind*(Z-eta))%*%xi.old
  return(ifelse(cen==1,lambda*exp(H),1)*exp(-cumhaz*exp(H)))
}

#----------------EM algorithm -----------------------------------------------------------#
EM_step <- function(mu, sigma, estimate, cumhaz, rs, X, Z, cen, estimate_update, m, n, P){
  
  node_a<-node_b<-wt_a<-wt_b<-matrix(0,nrow=n,ncol=m)
  for(i in c(1:n)){
    moments_a<-truncated_normal_a_moments(order_max = 2*m+1,mu = mu,sigma = sigma,a = Z[i])
    moments_b<-truncated_normal_b_moments(order_max = 2*m+1,mu = mu,sigma = sigma,b = Z[i])
    rule_a<-moment_method(m=m, moments=moments_a, boundary=Z[i], LB=TRUE)
    rule_b<-moment_method(m=m, moments=moments_b, boundary=Z[i], LB=FALSE)
    node_a[i,]<-rule_a$node   
    node_b[i,]<-rule_b$node
    wt_a[i,]  <-rule_a$weight*(1-pnorm(q = Z[i],mean = mu ,sd = sigma))
    wt_b[i,]  <-rule_b$weight*(pnorm(q = Z[i],mean = mu ,sd = sigma))
  }
  
  lambda <-c(cumhaz[1],diff(cumhaz))
  fi_numerB<-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_b[,o1], xi.old=estimate, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z, Ind=1)})
  fi_numerA<-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_a[,o1], xi.old=estimate, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z, Ind=0)})
  fi_denom<-apply(fi_numerA*wt_a,1,sum)+apply(fi_numerB*wt_b,1,sum)
  fi_B    <-fi_numerB/replicate(fi_denom,n = m)
  fi_A    <-fi_numerA/replicate(fi_denom,n = m)
  
  IZ_A      <-matrix(0,nrow=n,ncol=m)
  IZ_B      <-matrix(1,nrow=n,ncol=m)
  HJ_A      <-sapply(c(1:m),function(o){cbind(X,Z,0,0)%*%estimate})
  HJ_B      <-sapply(c(1:m),function(o){cbind(X,Z,1,Z-node_b[,o])%*%estimate})
  IZ        <-apply(IZ_A*fi_A*wt_a,1,sum)+apply(IZ_B*fi_B*wt_b,1,sum)
  exp_eta   <-apply(node_a*fi_A*wt_a,1,sum)+apply(node_b*fi_B*wt_b,1,sum)
  exp_eta_IZ<-apply(node_a*IZ_A*fi_A*wt_a,1,sum)+apply(node_b*IZ_B*fi_B*wt_b,1,sum)
  exp_eta_sq<-apply(node_a^2*fi_A*wt_a,1,sum)+apply(node_b^2*fi_B*wt_b,1,sum)
  exp_hj    <-apply(exp(HJ_A)*fi_A*wt_a,1,sum)+apply(exp(HJ_B)*fi_B*wt_b,1,sum)
  exp_hj_IZ <-apply(IZ_A*exp(HJ_A)*fi_A*wt_a,1,sum)+apply(IZ_B*exp(HJ_B)*fi_B*wt_b,1,sum)
  exp_hj_IZe<-apply(node_a*IZ_A*exp(HJ_A)*fi_A*wt_a,1,sum)+apply(node_b*IZ_B*exp(HJ_B)*fi_B*wt_b,1,sum)
  exp_hj_IZesq<-apply(node_a^2*IZ_A*exp(HJ_A)*fi_A*wt_a,1,sum)+apply(node_b^2*IZ_B*exp(HJ_B)*fi_B*wt_b,1,sum)
  
  if(estimate_update==TRUE){
    XZmatrix    <-cbind(X,Z,1,Z)
    XZmatrix_str<-(XZmatrix*cbind(replicate(P+1,exp_hj),replicate(2,exp_hj_IZ))-cbind(matrix(0,nrow = nrow(XZmatrix),ncol = P+2),exp_hj_IZe))
    den_1       <-sapply(c(1:n),function(i){(sum(exp_hj[rs[[i]]]))^{-1}})
    den_2       <-sapply(c(1:n),function(i){(sum(exp_hj[rs[[i]]]))^{-2}})
    
    grad_xi     <-apply(replicate(P+3,cen)*(cbind(X,Z,IZ,(IZ*Z-exp_eta_IZ))-
                                              replicate(P+3,den_1)*rbind(t(sapply(c(1:(n-1)),function(i){apply(XZmatrix_str[rs[[i]],],2,sum)})), XZmatrix_str[n,])),2,sum)
    
    hess_xi<-sapply(1:(P+3), function(o1){sapply(1:(P+3), function(o2){
      num_1<-(XZmatrix[,o1]*XZmatrix[,o2]*((o1<=(P+1))*(o2<=(P+1))*exp_hj+(1-(o1<=(P+1))*(o2<=(P+1)))*exp_hj_IZ)-
                (o1==(P+3))*XZmatrix[,o2]*exp_hj_IZe-(o2==(P+3))*XZmatrix[,o1]*exp_hj_IZe+(o1==(P+3))*(o2==(P+3))*exp_hj_IZesq)
      sum(cen*sapply(c(1:n),function(i){
        -den_1[i]*{sum(num_1[rs[[i]]])}+
          den_2[i]*{sum(XZmatrix_str[rs[[i]],o1])}*{sum(XZmatrix_str[rs[[i]],o2])} })
      )})})
    
    estimate_M<-estimate-solve(hess_xi)%*%grad_xi
  }else{estimate_M<-estimate}
  
  fi_numerB <-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_b[,o1], xi.old=estimate_M, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z, Ind=1)})
  fi_numerA <-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_a[,o1], xi.old=estimate_M, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z, Ind=0)})
  fi_denom  <-apply(fi_numerA*wt_a,1,sum)+apply(fi_numerB*wt_b,1,sum)
  fi_B      <-fi_numerB/replicate(fi_denom,n = m)
  fi_A      <-fi_numerA/replicate(fi_denom,n = m)
  HJ_A      <-sapply(c(1:m),function(o){cbind(X,Z,0,0)%*%estimate_M})
  HJ_B      <-sapply(c(1:m),function(o){cbind(X,Z,1,Z-node_b[,o])%*%estimate_M})
  exp_eta   <-apply(node_a*fi_A*wt_a,1,sum)+apply(node_b*fi_B*wt_b,1,sum)
  exp_eta_sq<-apply(node_a^2*fi_A*wt_a,1,sum)+apply(node_b^2*fi_B*wt_b,1,sum)
  exp_hj    <-apply(exp(HJ_A)*fi_A*wt_a,1,sum)+apply(exp(HJ_B)*fi_B*wt_b,1,sum)
  cumhaz_M  <-cumsum(sapply(1:n, function(o2){cen[o2]/sum(exp_hj[rs[[o2]]])})) 
  
  if(estimate_update==TRUE){
    mu_M      <-sum(exp_eta)/n                                                    
    sigma_M   <-sqrt(sum(exp_eta_sq)/n-mu_M^2)
    return(list(mu_M=mu_M,sigma_M=sigma_M,estimate_M=estimate_M,cumhaz_M=cumhaz_M))
  }else{
    return(list(cumhaz_M=cumhaz_M))
  }
}

pl <- function(mu, sigma, estimate, cumhaz, X, Z, cen, n, m){
  
  node_a<-node_b<-wt_a<-wt_b<-matrix(0,nrow=n,ncol=m)
  for(i in c(1:n)){
    moments_a<-truncated_normal_a_moments(order_max = 2*m+1,mu = mu,sigma = sigma,a = Z[i])
    moments_b<-truncated_normal_b_moments(order_max = 2*m+1,mu = mu,sigma = sigma,b = Z[i])
    rule_a<-moment_method(m=m, moments=moments_a, boundary=Z[i], LB=TRUE)
    rule_b<-moment_method(m=m, moments=moments_b, boundary=Z[i], LB=FALSE)
    node_a[i,]<-rule_a$node   
    node_b[i,]<-rule_b$node
    wt_a[i,]  <-rule_a$weight*(1-pnorm(q = Z[i],mean = mu ,sd = sigma))
    wt_b[i,]  <-rule_b$weight*(pnorm(q = Z[i],mean = mu ,sd = sigma))
  }
  
  lambda <-c(cumhaz[1],diff(cumhaz))
  fi_numerB<-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_b[,o1], xi.old=estimate, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z,Ind=1)})
  fi_numerA<-sapply(c(1:m),function(o1){
    fi_given_theta(eta=node_a[,o1], xi.old=estimate, lambda=lambda, cumhaz=cumhaz, cen=cen, X=X, Z=Z,Ind=0)})
  fi_denom<-apply(fi_numerA*wt_a,1,sum)+apply(fi_numerB*wt_b,1,sum)
  
  return(log(fi_denom))
}


RCPsurvEST<-function(data, P, m=10, tolerance=10^{-3}, gamma0=NA, beta0=NA, alpha10=NA, alpha20=NA, mu0=NA, sigma0=NA, TRACE=FALSE){
  n      <-nrow(data)
  data    <- data[order(data[,2]),]
  data[,1]<-c(1:nrow(data))
  Ti_d   <- data[,2]; cen_d<-data[,3]; Z_d  <-data[,(4+P)]
  if(P==1){X_d  <-as.numeric(data[,c(4:(4+P-1))])}else{X_d  <-as.matrix(data[,c(4:(4+P-1))])}
  rs     <- apply(as.matrix(Ti_d),1,function(t){which(Ti_d >= t)}) #riskset indicators
  
  if(is.numeric(gamma0)==FALSE){gamma0 <- rep(0,P)}else
    if(length(gamma0)!=P){print("ERROR: The length of gamma0 does not equal to P.");break}
  if(is.numeric(beta0)==FALSE){beta0 <- 0}else
    if(length(beta0)!=1){print("ERROR: beta0 is not a constant.");break}
  if(is.numeric(alpha10)==FALSE){alpha10 <- 0}else
    if(length(alpha10)!=1){print("ERROR: alpha10 is not a constant.");break}
  if(is.numeric(alpha20)==FALSE){alpha20 <- 0}else
    if(length(alpha20)!=1){print("ERROR: alpha20 is not a constant.");break}
  if(is.numeric(mu0)==FALSE){mu0 <- median(Z_d)}else
    if(length(mu0)!=1){print("ERROR: mu0 is not a constant.");break}else if((mu0<min(Z_d))|(mu0>max(Z_d))){
      mu0<-max(min(mu0,max(Z_d)),min(Z_d)); print("ERROR: The initial value for mu0 is out of the range of Z.")}
  if(is.numeric(sigma0)==FALSE){sigma0 <- 2}else
    if((length(sigma0)!=1)|(sigma0<=0)){print("ERROR: sigma0 is not a positive constant.");break}
  
  if(length(gamma0)!=P){print("ERROR: The length of gamma0 is not equal to P."); break}
  cumhaz   <-seq(10^{-4},1,length.out=n)
  estimate <-c(gamma0, beta0, alpha10, alpha20)
  mu       <-mu0
  sigma    <-max(sigma0,0.05)
  
  conv<-1
  I<-R<-tick<-0
  repeat{
    I<-I+1
    
    if(tick==0){theta0<-c(cumhaz,estimate,mu,sigma)}
    if(tick==2){
      r<-(theta1-theta0)
      v<-(theta2-theta1)-r
      a<- -sqrt(sum(r^2))/sqrt(sum(v^2))
      theta.dash<-theta0-2*a*r+a^2*v
      
      if((max(abs(2*a*r+a^2*v))<5)&prod(theta.dash[1:n]>=0)&(theta.dash[n+4+P+1]>=0.05)){
        cumhaz <-theta.dash[1:n]; estimate<-theta.dash[(n+1):(n+4+P-1)];
        mu     <-theta.dash[n+4+P]; sigma   <-theta.dash[n+4+P+1]
      }}
    
    #X = X_d; Z = Z_d; cen = cen_d
    EM_step_result<- EM_step(mu = mu,sigma = sigma,estimate = estimate,cumhaz = cumhaz,rs = rs, m=m, n=n, P=P,
                             X = X_d, Z = Z_d, cen = cen_d, estimate_update = TRUE)
    
    d0<-cumhaz-EM_step_result$cumhaz_M
    d1<-estimate-EM_step_result$estimate_M
    d2<-c(mu,sigma)-c(EM_step_result$mu_M,EM_step_result$sigma_M)
    
    mu       <- EM_step_result$mu_M
    sigma    <- max(EM_step_result$sigma_M,0.05)
    estimate <- EM_step_result$estimate_M
    cumhaz   <- EM_step_result$cumhaz_M
    
    obs.like<-sum(pl(estimate = estimate,mu =mu, sigma=sigma,cumhaz = cumhaz, X = X_d, Z = Z_d,cen = cen_d, n=n, m=m))
    
    dist_MLE<-max(abs(c(d0,d1,d2)))
    if((R==1)&(tick==2)&(dist_MLE>=tolerance)){R<-0}
    if((R==1)&(tick==2)&(dist_MLE<tolerance)){R<-2}
    if((R==0)&(tick==2)&(dist_MLE<tolerance)){R<-1}
    
    if(TRACE==TRUE){print(round(c(I, estimate, mu, sigma, dist_MLE,obs.like,R),4))}
    
    if(tick==0){theta1<-c(cumhaz,estimate,mu,sigma)}
    if(tick==1){theta2<-c(cumhaz,estimate,mu,sigma)}
    if(tick==2){tick<- -1}
    tick<-tick+1
    
    if(R==2){
      conv        <-0
      xi.hat      <-c(estimate, mu, sigma)
      Lambda.hat  <-cumhaz
      pl.hat      <-pl(estimate = estimate,mu =mu, sigma=sigma,cumhaz = cumhaz, X = X_d, Z = Z_d,cen = cen_d, n=n, m=m)
      like.hat    <-sum(pl.hat)
      break}
  }
  
  hn<-5*n^{-1}
  storeinf_ph<-storeinf_mh<-c()
  for (j in 1:(P+5)){
    Iinf        <-0
    hn.vec      <-numeric(P+5)
    hn.vec[j]   <-hn
    xi.p        <-xi.hat+hn.vec
    Lambda.star <-Lambda.hat
    repeat{
      Iinf <-Iinf+1
      EM_step_result<-EM_step(mu = xi.p[4+P],sigma = xi.p[4+P+1],estimate = xi.p[1:(4+P-1)],cumhaz = Lambda.star,rs = rs, m=m, n=n, P=P,
                              X = X_d, Z = Z_d,cen = cen_d, estimate_update = FALSE)
      
      d0_inf     <-Lambda.star-EM_step_result$cumhaz_M
      Lambda.star<-EM_step_result$cumhaz_M
      
      if((max(abs(c(d0_inf)))<tolerance)|Iinf>50){
        storeinf_ph<-rbind(storeinf_ph, c(xi.p, Lambda.star))
        break}
    }
  }
  
  for (j in 1:(P+5)){
    Iinf        <- 0
    hn.vec      <-numeric(P+5)
    hn.vec[j]   <-hn
    xi.p        <-xi.hat-hn.vec
    Lambda.star <-Lambda.hat
    repeat{
      Iinf <-Iinf+1
      EM_step_result<-EM_step(mu = xi.p[4+P],sigma = xi.p[4+P+1],estimate = xi.p[1:(4+P-1)],cumhaz = Lambda.star,rs = rs, m=m, n=n, P=P,
                              X = X_d, Z = Z_d, cen = cen_d, estimate_update = FALSE)
      
      d0_inf     <-Lambda.star-EM_step_result$cumhaz_M
      Lambda.star<-EM_step_result$cumhaz_M
      
      if((max(abs(c(d0_inf)))<tolerance)|Iinf>50){
        storeinf_mh<-rbind(storeinf_mh, c(xi.p, Lambda.star))
        break}
    }
  }
  
  dli<-matrix(0,nrow=n,ncol=(P+5))
  for(j in c(1:(P+5))){dli[,j]<-(pl(estimate = storeinf_ph[j,c(1:(4+P-1))],mu = storeinf_ph[j,(4+P)], sigma = storeinf_ph[j,(4+P+1)],
                                    cumhaz = storeinf_ph[j,c((4+P+2):ncol(storeinf_ph))], X = X_d, Z = Z_d, cen = cen_d, m=m, n=n)-
                                   pl(estimate = storeinf_mh[j,c(1:(4+P-1))],mu = storeinf_mh[j,(4+P)], sigma = storeinf_mh[j,(4+P+1)],
                                      cumhaz = storeinf_mh[j,c((4+P+2):ncol(storeinf_mh))], X = X_d, Z = Z_d, cen = cen_d, m=m, n=n))/2/hn}
  dli2<-matrix(0,nrow=(P+5),ncol=(P+5))
  for(i in c(1:n)){dli2<-dli2+as.vector(dli[i,])%o%as.vector(dli[i,])}
  
  est.se_f<-sqrt(diag(solve(dli2)))
  
  if((abs(xi.hat[P+2])<0.5)&(est.se_f[P+2]>5)){
    est.se_f_sub<-sqrt(diag(solve(dli2[-(P+2),-(P+2)])))
    se.vec      <-c(est.se_f_sub[1:(P+1)],NA,est.se_f_sub[(P+2):(P+4)])
  }else{se.vec  <-est.se_f}
  
  return(list(loglik=obs.like,gamma.hat=xi.hat[1:P], beta.hat=xi.hat[P+1], alpha1.hat=xi.hat[P+2], alpha2.hat=xi.hat[P+3], 
              mu.hat=xi.hat[P+4], sigma.hat=xi.hat[P+5], gamma.hat.se=se.vec[1:P], beta.hat.se=se.vec[P+1], 
              alpha1.hat.se=se.vec[P+2], alpha2.hat.se=se.vec[P+3], mu.hat.se=se.vec[P+4], sigma.hat.se=se.vec[P+5],
              conv=conv))
}




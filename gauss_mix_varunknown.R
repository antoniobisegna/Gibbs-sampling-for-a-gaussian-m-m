data= mix.gauss(N=100,
                std = c(3,3),
                K=2,
                pi.c = c(.5,.5),
                mu=c(-10,10),
                plot = T)

library('mvtnorm')
library('bpp')
library('gtools')


sample_z = function(x,pi,mu,sig2){
  z = rep(0, length(x))
  for(i in 1:length(z)){
    std=sqrt(sig2)
    p.z.given.x = as.vector(pi) * dnorm(x[i],mu,std)
    z[i] = sample(1:length(pi), size=1,prob=p.z.given.x)
  }
  return(z)
}


sample_pi = function(z,k,lambda){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi = gtools::rdirichlet(1,counts+lambda)
  return(pi)
}

Gibbs_gaus_mix = function(G,burnin,thin,K,data,lambda=1){
  
  #Initializing storage
  
  Z         = matrix(NA,nrow = G+burnin, ncol = length(data))
  sig2.post = matrix(NA,nrow = G+burnin, ncol = K)
  mu.post   = matrix(NA,nrow = G+burnin, ncol = K)
  pi.post   = matrix(NA,nrow = G+burnin, ncol = K)
  
  Iterations = burnin+(thin*G)
  
  #Inizializing parameters
  sig2.post[1,] = 1/(rgamma(K,1,1))
  mu.post[1,] = rnorm(K,0,sqrt(sig2.post[1,])) 
  pi.post[1,] = rep(1/K,K)
  Z[1,]       = sample_z(x=data,mu = mu.post[1,],pi=pi.post[1,],sig2 = sig2.post[1,])
  
  
  #Inizzializing tmp
  mu.tmp   = mu.post[1,]
  sig2.tmp = sig2.post[1,]
  pi.tmp   = pi.post[1,]
  Z.tmp    = Z[1,]
  
  for (iter in 2:Iterations) {
    
    #############################
    ### STEP 1: SAMPLE MU ###
    #############################
    for (i in 1:K) {
      # prior info
      mu0 = mu.tmp[i]
      sd0 = sqrt(sig2.tmp[i])
      k0  = 1
      nu0 = 1
      
      # First set n, ybar and s^2
      ybar  = mean(data[Z.tmp==i])
      s = sd(data[Z.tmp==i])
      n = length(which(Z.tmp==i))
      
      # posterior inference
      kn  <- k0+n
      nun <- nu0+n
      mun <- (k0*mu0+n*ybar)/kn
      s2n <- (nu0*(sd0^2)+(n-1)*(s^2)+k0*n*(ybar-mu0)^2/(kn))/(nun)

      sig2.tmp[i] <- 1/rgamma(1,nun/2,s2n*nun/2)
      mu.tmp[i]   <- rnorm(1,mun,sqrt(sig2.tmp[i]/kn) )
    }
    if(iter<= burnin | iter%%thin == 0){
      sig2.post[iter,] = sig2.tmp
      mu.post[iter,]   = mu.tmp
    }
    
    #############################
    ### STEP 2: SAMPLE pi ###
    #############################
    pi.tmp= sample_pi(z=Z.tmp,
                      k=max(unique(Z.tmp)),
                      lambda = lambda)
    if(iter<= burnin | iter%%thin == 0){
      pi.post[iter,]=pi.tmp
    }
    
    #############################
    ### STEP 3: SAMPLE Z  ###
    #############################
    Z.tmp = sample_z(x=data,
                     pi = pi.tmp,
                     mu = mu.tmp,
                     sig2 = sig2.tmp)
    if(iter<= burnin | iter%%thin == 0){
      Z[iter,]=Z.tmp
    }
  }
  
  results      = list()
  results$sig2 = sig2.post
  results$mu   = mu.post
  results$Z    = Z
  results$pi   = pi.post
  return(results)
}

proto  = Gibbs_gaus_mix(G=5000,
                        burnin = 1000,
                        thin = 1,
                        K=2,
                        data=data$samples)

plot(proto$mu[,1],type='l')
plot(proto$mu[,2],type='l')



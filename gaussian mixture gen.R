####################################################
##### Code to generate a mixture of Gaussians  #####
####################################################

mix.gauss = function(N, K, pi.c, mu, std, plot=F){
  
  #Sample N components based on their pi (mixing proportion)
  components = sample(1:K, prob=pi.c, size=N, replace=TRUE)
  
  #Generate N random samples from a 1D-Gaussian, with the corresponding weights
  samples    = rnorm(n=N,mean=mu[components],sd=std[components])
  
  if(plot == T){
    par(mfrow=c(1,2))
    
    #Density plot
    plot(density(samples),main=paste0('Gaussian mixture of ',K,' components'))
    y=rep(0,length(samples))
    points(samples,y,col=components,cex=0.75)
    
    #Barplot
    barplot(table(components)/length(components),
            xlab = 'Components',
            ylab = 'Frequency')
    par(mfrow=c(1,1))
  }
  results            = list()
  results$components = components
  results$samples    = samples
  
  return (results)
}




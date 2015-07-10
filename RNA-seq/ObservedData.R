#generate observed RNA-seq
ObservedData<- function(){
  
  p <- 0.05 #state noise
  q=.1 #observation noise
  n <- 100 # number of observations
  ch <-50 # chaange ponit
  seq_depth = 10.75 #
  mu = .01
  #state noise
  s_noise1 = rbinom(n, 1, p)
  s_noise2 = rbinom(n, 1, p)
  s_noise3 = rbinom(n, 1, p)
  s_noise4 = rbinom(n, 1, p)
  ATM<- vector("numeric",n)#initial states
  P53 <- vector("numeric",n)
  WIP1 <- vector("numeric",n)
  MDM2 <- vector("numeric",n)
  ATM[1] = 0
  P53[1] = 1
  WIP1[1] = 0
  MDM2[1] = 1
  for (k in 2:n){
    
    ATM[k] = !WIP1[k-1] & (ATM[k-1] | dna_dsb)+ s_noise1[k]
    if (ATM[k]==2)
      ATM[k] =0
    
    if (k>ch) ATM[k] =0
    
    if (ATM[k]==0)
      l1 = seq_depth* exp(mu +theta1[k])
    if(ATM[k]==1)
      l1 = seq_depth* exp(mu + delta1+theta1[k])
    
    
    P53[k] = !MDM2[k-1] & (ATM[k-1] | WIP1[k-1]) + s_noise2[k]
    if (P53[k]==2)
      P53[k]=0
  #  if(k>ch)  P53[k] =1
    
    if (P53[k]==0)
      l2 = seq_depth* exp(mu +theta2[k])
    if(P53[k]==1)
      l2 = seq_depth* exp(mu + delta2+theta2[k])
    
    
    WIP1[k] = P53[k-1] + s_noise3[k]
    if (WIP1[k]==2)
      WIP1[k]=0
    
    # if (k>ch) WIP1[k] =0
    if (WIP1[k]==0)
      l3 = seq_depth* exp(mu +theta3[k])
    if(WIP1[k]==1)
      l3 = seq_depth* exp(mu + delta3+theta3[k])
    
    
    MDM2[k] = !ATM[k-1] & (P53[k-1] | WIP1[k-1]) + s_noise4[k]
    if (MDM2[k]==2)
      MDM2[k] =0
    
   
    if (MDM2[k]==0)
      l4 = seq_depth* exp(mu +theta4[k])
    if(MDM2[k]==1)
      l4 = seq_depth* exp(mu + delta4+theta4[k])
    
    #observations
    Y1[k] = rpois(1, l1)
    Y2[k] = rpois(1, l2)
    Y3[k] = rpois(1, l3)
    Y4[k] = rpois(1, l4)
    
  }
  return(list(Y1,Y2,Y3,Y4))
  
}
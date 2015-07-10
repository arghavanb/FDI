whitenesstest<- function(residual,lag){
    nr <- nrow(residual)
    ccc <- vector("numeric",)
   
    window <- 150
   Q <- vector("numeric",)
    k=1
      for (i in 1:(nr-window)){
          # for (i in 51:nr){
          # if(i <= (nr-window)){
        LBtest<- mq(residual[i:(i+window),],lag)
        Q[i] <- LBtest[lag,2]
        #  LBtest <- mq(residual[1:i,],lag)
        if(LBtest[lag,3] < .05){
             ccc[k] <- i +window
            # ccc[k] <- i
            k<- k+1
            
        }
        
            }
    
    #   plot(Q,type="l")
    return(ccc)
  
}
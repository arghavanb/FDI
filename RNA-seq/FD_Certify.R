FD_Certify <- function( xx){
    FCT <- vector("numeric",)
    if(length(which(xx<200))>6){
        
    chk =10 # interval around the faulty point that needed to be tested.
    k = -chk/2
    FDT<- xx[xx < 200]
	l<- length(FDT)
	t=1
	
    
	for (i in (chk/2+1):l){
		cnt <-0
		# checking 10 points befor and after point QQ[a] if 5 point out of this 20 points are faulty so we decide that this point is valid as a faulty point.
      
		for (j in 1:chk){
			a <- FDT[i]
			if(any(FDT,a)){
			cnt<- cnt+1
			}
            
			k <- k+1
		}
		if (cnt>5){
            
            FCT[t] = FDT[i]
        t<- t+1
		}
		
	}
    }
   
    
	return(FCT)
}
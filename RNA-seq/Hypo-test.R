Hypothesis <- function(fnrm,fp1,fa1,fw1,fm1,fa0,fp0,fw0,fm0,fault){
    m <- length(fp1)
hyponrm <- vector("numeric",m)
hypoa1 <- vector("numeric",m)
hypop1 <- vector("numeric",m)
hypow1 <- vector("numeric",m)
hypom1 <- vector("numeric",m)
hypoa0 <- vector("numeric",m)
hypop0 <- vector("numeric",m)
hypow0 <- vector("numeric",m)
hypom0 <- vector("numeric",m)
K<- 8
hyponrm[fault-1] <- 1/K
hypoa1[fault-1] = 1/K
hypop1[fault-1] = 1/K
hypow1[fault-1] = 1/K
hypom1[fault-1] = 1/K
hypoa0[fault-1] = 1/K
hypop0[fault-1] = 1/K
hypow0[fault-1] = 1/K
hypom0[fault-1] = 1/K
cnt <-0
aa<- vector("numeric",)

for (l in fault:m){
    
    #	hyponrm[l] <- (fnrm[l] * hyponrm[l-1])/(fnrm[l] * hyponrm[l-1]+fp1[l] *hypop1[l-1]+fa1[l]*hypoa1[l-1]+fw1[l] * hypow1[l-1]+fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
    
    hypop1[l] <- (fp1[l] * hypop1[l-1])/(fnrm[l] * hyponrm[l-1]+ fp1[l] *hypop1[l-1] +fa1[l]*hypoa1[l-1]+fw1[l] * hypow1[l-1]+fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
    
    hypoa1[l] <- (fa1[l] * hypoa1[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1]+fw1[l] * hypow1[l-1]+fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
    
     hypow1[l] <- (fw1[l] * hypow1[l-1])/(fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1]+fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
    
    
    hypom1[l] <- (fm1[l] * hypom1[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1] +fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
    
    hypoa0[l] <- (fa0[l] * hypoa0[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1] +fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1]+fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])

 hypop0[l] <- (fp0[l] * hypop0[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1] +fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1] +fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
 
     hypow0[l] <- (fw0[l] * hypow0[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1] +fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1] +fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
     
     hypom0[l] <- (fm0[l] * hypom0[l-1])/( fp1[l] *hypop1[l-1]+ fa1[l]*hypoa1[l-1] + fw1[l] * hypow1[l-1] +fm1[l] * hypom1[l-1]+fa0[l] * hypoa0[l-1] +fp0[l] * hypop0[l-1]+fw0[l] * hypow0[l-1]+fm0[l] * hypom0[l-1])
     
     temp <- c(hypoa1[l],hypop1[l],hypow1[l],hypom1[l],hypoa0[l],hypop0[l],hypow0[l],hypom0[l])
     #aa[l]<- which.max(temp)
      if( length(which(temp>.8))!=0){
         cntrl<- which(temp>.8)
         if (cntrl ==8){#cntrl shows which gene is stuck
     flag=1
    diag = c(flag,l)
    return(diag)
      }
     }
     # if(5[l]<.5)
     # if(aa[l] !=1)
     # cnt=cnt+1
}

#par(mfrow=c(4,2))
#plot(hyponrm,type="l")
# plot(hypoa1,xlim=c(199,400),ylim=c(0:1),ylab='ATM s-a-1',type="l")


#plot(hypop1,xlim=c(199,400),ylim=c(0:1),ylab='P53 s-a-1',type="l")

#plot(hypow1,xlim=c(199,400),ylim=c(0:1),ylab='WIP1 s-a-1',type="l")

#plot(hypom1,xlim=c(199,400),ylim=c(0:1),ylab='MDM2 s-a-1',type="l")

#plot(hypoa0,xlim=c(199,400),ylim=c(0:1),ylab='ATM s-a-0',type="l")

#plot(hypop0,xlim=c(199,400),ylim=c(0:1),ylab='P53 s-a-0',type="l")

#plot(hypow0,xlim=c(199,400),ylim=c(0:1),ylab='WIP1 s-a-0',type="l")

#plot(hypom0,xlim=c(199,400),ylim=c(0:1),ylab='MDM2 s-a-0',type="l")

#return(cnt)
}


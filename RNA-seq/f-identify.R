#source("Prediction_M.R")
#source("Update.R")

#source("Hypo-test.R")
fault_identify<- function(fault,Y1,Y2,Y3,Y4,seq_depth,theta,delta,PVD){

#Staes
gene <- 4 # [ATM,P53,WIP1,MDM2]
num_states <- 2^gene #withot dna_dsb
A <- matrix(nrow=num_states,ncol =gene)
for (i in 1:gene){
    temp = rbind(array(0,c(num_states/2^i,1)),array(1,c(num_states/2^i,1)))
    A[,i] = kronecker(matrix(1,2^(i-1)),temp)
}
Aa1 = A
Aa1[,1]=1
Aa1 =t(Aa1)
Ap1 = A
Ap1[,2]=1
Ap1 =t(Ap1)
Aw1 = A
Aw1[,3]=1
Aw1 =t(Aw1)
Am1 = A
Am1[,4]=1
Am1 =t(Am1)
Aa0 = A
Aa0[,1]=0
Aa0 =t(Aa0)

Ap0 = A
Ap0[,2]=0
Ap0 =t(Ap0)
Aw0 = A
Aw0[,3]=0
Aw0 =t(Aw0)
Am0 = A
Am0[,4]=0
Am0 =t(Am0)
A= t(A)
p <- 0.05
#
#seq_depth = runif(1, min=2, max=3.75)# sequencing depth
 # for [250K-300K]
# seq_depth <- (0.035+2)/2 #1K-50K: [0.035, 2.00]
mu = .01
delta1 =  rnorm(1,1,.05)
delta2 = rnorm(1,1,.05)
delta3 =  rnorm(1,1,.05)
delta4 =  rnorm(1,1,.05)
delta = c(delta1,delta2,delta3,delta4)
##
#P <- #matrix(1/16,num_states,1)#initial PVD I assume it is uniform
Pp1 <- PVD[,fault]#matrix(1/16,16,1)
Pa1 <- PVD[,fault]
Pw1 <- PVD[,fault]
Pm1 <- PVD[,fault]
Pm1 <- PVD[,fault]
Pa0 <- PVD[,fault]
Pp0 <- PVD[,fault]
Pw0 <- PVD[,fault]
Pm0 <- PVD[,fault]
ATM<- vector("numeric",n)#initial statexxs
P53 <- vector("numeric",n)
WIP1 <- vector("numeric",n)
MDM2 <- vector("numeric",n)
proby <- vector("numeric",n)
res_erra1 <- matrix(nrow=4,ncol=n)
res_errp1 <- matrix(nrow=4,ncol=n)
res_errw1 <- matrix(nrow=4,ncol=n)
res_errm1 <- matrix(nrow=4,ncol=n)
res_erra0 <- matrix(nrow=4,ncol=n)
res_errp0 <- matrix(nrow=4,ncol=n)
res_errw0 <- matrix(nrow=4,ncol=n)
res_errm0 <- matrix(nrow=4,ncol=n)
n<- 400
#dna_dsb <- cbind(array(0,c(1,1320)),array(1,c(1,n-1320)))
#dna_dsb <- array(1,c(1,n))

fnrm <- vector("numeric",n)
fp1 <- vector("numeric",n)
fa1 <- vector("numeric",n)
fw1 <- vector("numeric",n)
fm1<-vector("numeric",n)
fa0 <-vector("numeric",n)
fp0 <-vector("numeric",n)
fw0 <-vector("numeric",n)
fm0 <-vector("numeric",n)
B <- vector("numeric",num_states)
#D <- D_matrix(A,mu,delta,theta,seq_depth)
Mk = prediction_M(A,p)
Mp1= Stuck_P531(Ap1,p)
Ma1 =Stuck_ATM1(Aa1,p)
Mw1 =Stuck_wip1(Aw1,p)
Mm1 =Stuck_mdm1(Am1,p)
Ma0 =Stuck_ATM0 (Aa0,p)
Mp0 =Stuck_P530(Ap0,p)
Mw0 =Stuck_wip0(Aw0,p)
Mm0 =Stuck_mdm0(Am0,p)

for (k in fault:n){
    

   
    theta = c(theta1[k],theta2[k],theta3[k],theta4[k])
    y <- rbind(Y1[k],Y2[k],Y3[k],Y4[k])
    #  t <-Update_M(y,A,delta,theta,seq_depth)
    
    # Tk <- diag(t)
    # PP <- Mk %*% P
    #  B<- Tk %*% PP
    # proby[k] <- t %*% PP
    # P <- B/norm(as.matrix(B))
    # fnrm[k] = t%*%PP
    
    ###
    tp1 <-Update_M(y,Ap1,delta,theta,seq_depth)
    
    Tp1 <- diag(tp1)
    PPp1 <- Mp1 %*% Pp1
    Bp1<- Tp1 %*% PPp1
    Pp1 <- Bp1/norm(as.matrix(Bp1))
    fp1[k] = tp1%*%PPp1
  
    ####
    ta1 <-Update_M(y,Aa1,delta,theta,seq_depth)
    
    Ta1 <- diag(ta1)
    PPa1 <- Ma1 %*% Pa1
    Ba1<- Ta1 %*% PPa1
    Pa1 <- Ba1/norm(as.matrix(Ba1))
    fa1[k] = ta1%*%PPa1
   
  

    ####
    tw1 <-Update_M(y,Aw1,delta,theta,seq_depth)
    
    Tw1 <- diag(tw1)
    PPw1 <- Mw1 %*% Pw1
    Bw1<- Tw1 %*% PPw1
    Pw1 <- Bw1/norm(as.matrix(Bw1))
    fw1[k] = tw1%*%PPw1
   
    ####
    tm1 <-Update_M(y,Am1,delta,theta,seq_depth)
    
    Tm1 <- diag(tm1)
    PPm1 <- Mm1 %*% Pm1
    Bm1<- Tm1 %*% PPm1
    Pm1 <- Bm1/norm(as.matrix(Bm1))
    fm1[k] = tm1%*%PPm1
 
    
    ###
    ta0 <-Update_M(y,Aa0,delta,theta,seq_depth)
    
    Ta0 <- diag(ta0)
    PPa0 <- Ma0 %*% Pa0
    Ba0<- Ta0 %*% PPa0
    Pa0 <- Ba0/norm(as.matrix(Ba0))
    fa0[k] = ta0%*%PPa0
   
    
    ###
    tp0 <-Update_M(y,Ap0,delta,theta,seq_depth)
    Tp0 <- diag(tp0)
    PPp0 <- Mp0 %*% Pp0
    Bp0<- Tp0 %*% PPp0
    Pp0 <- Bp0/norm(as.matrix(Bp0))
    fp0[k] = tp0%*%PPp0
  
    
    ###
    tw0 <-Update_M(y,Aw0,delta,theta,seq_depth)
    Tw0 <- diag(tw0)
    PPw0 <- Mw0 %*% Pw0
    Bw0<- Tw0 %*% PPw0
    Pw0 <- Bw0/norm(as.matrix(Bw0))
    fw0[k] = tw0%*%PPw0

    
    ###
    tm0 <-Update_M(y,Am0,delta,theta,seq_depth)
    Tm0 <- diag(tm0)
    PPm0 <- Mm0 %*% Pm0
    Bm0<- Tm0 %*% PPm0
    Pm0 <- Bm0/norm(as.matrix(Bm0))
    fm0[k] = tm0%*%PPm0
 
    
    
}

ss <- Hypothesis(fnrm,fp1,fa1,fw1,fm1,fa0,fp0,fw0,fm0,fault)
return(ss)

}
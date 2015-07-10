source("Prediction_M.R")
source("Update.R")

source("Hypo-test.R")
fmr <- vector("numeric",100)
#for (l in 1:100){
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
    p <- 0.08 #state noise
    q=.1
    n <- 400 # number of observations
    ch <-200 # chaange ponit
    s_noise1 = rbinom(n, 1, p)
    s_noise2 = rbinom(n, 1, p)
    s_noise3 = rbinom(n, 1, p)
    s_noise4 = rbinom(n, 1, p)
    #
    #seq_depth = runif(1, min=2, max=3.75)# sequencing depth
     seq_depth =  (9.75+11.75)/2 # for [250K-300K]
    # seq_depth <- (0.035+2)/2 #1K-50K: [0.035, 2.00]
    mu = .01
    delta1 =  rnorm(1,1,.05)
    delta2 = rnorm(1,1,.05)
    delta3 =  rnorm(1,1,.05)
    delta4 =  rnorm(1,1,.05)
    delta = c(delta1,delta2,delta3,delta4)
    theta1 = rnorm(n,0,q)
    theta2 =  rnorm(n,0,q)
    theta3= rnorm(n,0,q)
    theta4 = rnorm(n,0,q)
    #theta  = c(theta1[1],theta2[1],theta3[1],theta4[1])
    #
    P <- matrix(1/16,16,1)#initial PVD I assume it is uniform
    Pp1 <- matrix(1/16,16,1)
    Pa1 <- matrix(1/16,16,1)
    Pw1 <- matrix(1/16,16,1)
    Pm1 <- matrix(1/16,16,1)
    Pm1 <- matrix(1/16,16,1)
    Pa0 <- matrix(1/16,16,1)
    Pp0 <- matrix(1/16,16,1)
    Pw0 <- matrix(1/16,16,1)
    Pm0 <- matrix(1/16,16,1)
    ATM<- vector("numeric",n)#initial statexxs
    P53 <- vector("numeric",n)
    WIP1 <- vector("numeric",n)
    MDM2 <- vector("numeric",n)
    proby <- vector("numeric",n)
    ATM[1] = 1
    P53[1] = 1
    WIP1[1] = 0
    MDM2[1] = 0
    #context
    
    #dna_dsb <- cbind(array(0,c(1,1320)),array(1,c(1,n-1320)))
    dna_dsb <- array(1,c(1,n))
    Y1 <- vector("numeric",n)
    Y2 <- vector("numeric",n)
    Y3 <- vector("numeric",n)
    Y4 <- vector("numeric",n)
    res_err <- matrix(,nrow=gene,ncol=n)
    res_erra1 <- matrix(,nrow=gene,ncol=n)
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
    Mp1=prediction_M(Ap1,p)
    Ma1 =prediction_M(Aa1,p)
    Mw1 =prediction_M(Aw1,p)
    Mm1 =prediction_M(Am1,p)
    Ma0 =prediction_M(Aa0,p)
    Mp0 =prediction_M(Ap0,p)
    Mw0 =prediction_M(Aw0,p)
    Mm0 =prediction_M(Am0,p)


for (k in 2:n){
        
        ATM[k] = !WIP1[k-1] & (ATM[k-1] | dna_dsb[k-1])+ s_noise1[k]
        if (ATM[k]==2)
        ATM[k] =0
        
        # if (k>ch) { ATM[k] =1}
        
        if (ATM[k]==0)
        l1 = seq_depth* exp(mu +theta1[k])
        if(ATM[k]==1)
        l1 = seq_depth* exp(mu + delta1+theta1[k])
        
        
        P53[k] = !MDM2[k-1] & (ATM[k-1] | WIP1[k-1]) + s_noise2[k]
        if (P53[k]==2)
        P53[k]=0
        
        #  if(k>ch) P53[k] =0
        
        if (P53[k]==0)
        l2 = seq_depth* exp(mu +theta2[k])
        if(P53[k]==1)
        l2 = seq_depth* exp(mu + delta2+theta2[k])
        
        
        
        WIP1[k] = P53[k-1] + s_noise3[k]
        if (WIP1[k]==2)
        WIP1[k]=0
        
        if (k>ch) WIP1[k] =0
         
        if (WIP1[k]==0)
        l3 = seq_depth* exp(mu +theta3[k])
        if(WIP1[k]==1)
        l3 = seq_depth* exp(mu + delta3+theta3[k])
        
        
        MDM2[k] = !ATM[k-1] & (P53[k-1] | WIP1[k-1]) + s_noise4[k]
        if (MDM2[k]==2)
        MDM2[k] =0
        
        
        #  if (k>ch) MDM2[k] = 0
         
        if (MDM2[k]==0)
        l4 = seq_depth* exp(mu +theta4[k])
        if(MDM2[k]==1)
        l4 = seq_depth* exp(mu + delta4+theta4[k])
        
        #observations
        Y1[k] = rpois(1, l1)
        Y2[k] = rpois(1, l2)
        Y3[k] = rpois(1, l3)
        Y4[k] = rpois(1, l4)
        theta = c(theta1[k],theta2[k],theta3[k],theta4[k])
		y <- rbind(Y1[k],Y2[k],Y3[k],Y4[k])
        t <-Update_M(y,A,delta,theta,seq_depth)
      
        Tk <- diag(t)
        PP <- Mk %*% P
        B<- Tk %*% PP
        P <- B/norm(as.matrix(B))
        fnrm[k] = t%*%PP
        D <- D_matrix(A,mu,delta,theta,seq_depth)
        Apriori = D %*% PP#E[y[k]|y[k-1],...,y[1]]
    
        res_err[1,k] = Y1[k] - Apriori[1]
        res_err[2,k] = Y2[k] - Apriori[2]
        res_err[3,k] = Y3[k] - Apriori[3]
        res_err[4,k] = Y4[k] - Apriori[4]
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
        Da1 <- D_matrix(Aa1,mu,delta,theta,seq_depth)
        Aprioria1 = Da1 %*% PPa1#E[y[k]|y[k-1],...,y[1]]
        
        res_erra1[1,k] = Y1[k] - Aprioria1[1]
        res_erra1[2,k] = Y2[k] - Aprioria1[2]
        res_erra1[3,k] = Y3[k] - Aprioria1[3]
        res_erra1[4,k] = Y4[k] - Aprioria1[4]
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
#fa1<-fa1/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fp1<-fp1/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fw1<-fw1/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fm1<-fm1/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fa0<-fa0/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fp0<-fp0/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fw0<-fw0/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
#fm0<-fm0/(fa1+fp1+fw1+fm1+fa0+fp0+fw0+fm0)
fault=200

ss <- Hypothesis(fnrm,fp1,fa1,fw1,fm1,fa0,fp0,fw0,fm0,fault)
#fmr[l] <- ss/200
#}

#mfmr <- mean(fmr)
#mfmr
#sd<- sqrt(var(fmr))
#sd

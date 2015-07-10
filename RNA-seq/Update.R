Update_M <- function(y,Z,delta,theta,seq_depth){
    mu = .01
   
    num_states =ncol(Z)
	t <- vector("numeric",num_states)
	
	for (i in 1:num_states){
        
        if (Z[1,i] ==0){
              l1 = seq_depth*exp(mu+theta[1])
              # l1 = exp(mu+theta[1])
        }
        if (Z[1,i]==1){
            l1 = seq_depth*exp(mu+delta[1]+theta[1])
            # l1 = exp(mu+delta[1]+theta[1])
        }
        if (Z[2,i] ==0){
              l2 = seq_depth*exp(mu+theta[2])
              #l2 = exp(mu+theta[2])
        }
        if (Z[2,i] ==1){
            l2 = seq_depth*exp(mu+delta[2]+theta[2])
            # l2 = exp(mu+delta[2]+theta[2])
        }
        if (Z[3,i] ==0){
            l3 = seq_depth*exp(mu+theta[3])
            # l3 = exp(mu+theta[3])
        }
        if (Z[3,i] ==1){
             l3 = seq_depth*exp(mu+delta[3]+theta[3])
             #  l3 = exp(mu+delta[3]+theta[3])

        }
        if (Z[4,i] ==0){
            l4 = seq_depth*exp(mu+theta[4])
            # l4 = exp(mu+theta[4])
        }
        if (Z[4,i] ==1){
             l4 = seq_depth*exp(mu+delta[4]+theta[4])
             # l4 = exp(mu+delta[4]+theta[4])
        }
        
        # t[i] = exp(-l1-l2-l3-l4+ log(l1/seq_depth)*y[1]+ log(l2/seq_depth)*y[2]+log(l3/seq_depth)*y[3]+log(l4/seq_depth)*y[4])
       t[i] = dpois(y[1],l1)*dpois(y[2],l2)*dpois(y[3],l3)*dpois(y[4],l4)
    }
	return(t)
}
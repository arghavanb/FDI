#this fuction estimates expectation of y given x 
D_matrix <- function(A,mu,delta,theta,seq_depth){
	num_states <- ncol(A)
	d<- nrow(A)
	D <-matrix(nrow=d,ncol=num_states)
	for (i in 1:num_states){
        
		for (j in 1:d){
			if(A[j,i]==0){	D[j,i] = seq_depth* exp(mu +theta[j])	}
			if(A[j,i] ==1){ D[j,i] =  seq_depth* exp(mu + delta[j]+theta[j]) }
		}
	}
	return(D)
}
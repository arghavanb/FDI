prediction_M <- function(C,p){
	#p <- .08
	num_states <- ncol(C)
	M <- matrix(,nrow =num_states,ncol = num_states)
    
	for (i in 1:num_states){
		for(j in 1:num_states){
            
			X1 = (C[1,j] | dna_dsb) & !C[3,j]
			X2 = !C[4,j] & (C[1,j] | C[3,j])
            X3 = C[2,j]
            X4 = !C[1,j] & (C[2,j] | C[3,j])
            n1 = abs(C[1,i] - X1)
			n2 = abs(C[2,i] - X2)
            n3 = abs(C[3,i] - X3)
            n4 = abs(C[4,i] - X4)
			M[i,j] = (n1*p+(1-n1)*(1-p))*(n2*p+(1-n2)*(1-p))*(n3*p+(1-n3)*(1-p))*(n4*p+(1-n4)*(1-p))
		}
	}
	return(M)
}
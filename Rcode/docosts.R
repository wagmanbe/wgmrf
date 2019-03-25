###########################################################
### FUNCTIONS FOR CALCULATING TRAD,FIELDS,AND F-S COST
### DATE: OCT 1ST, 2016
### Benjamin Moore Wagman
### Adapted from new-costs.R (Alvaro Nosedal-Sanchez)
###########################################################

library('spam')

###  COST.nfields = Function that computes cost for 11 fields, 
###               ...1 region at a time, 1 season at a time. 

###  INPUTS: 1) List containing fields of climate data for a region and season
###          2) List containing fields of model data ordered same as obs.  
###          3) Q1 = precision matrix corresponding to the given region
###          4) S1 = estimated covariance matrix for the given region
### 	     26) Alpha = parameter for scaling GMRF cost. 

###  OUTPUTS: 11x11 matrix of cost. 
###       1) traditional cost:	  For S1=diagonal matrix & alpha=1
###	      2) fields cost: 	    	For S1 with off-diagonal values & alpha=1
###	      3) fields-space costs. 	For S1 with off-diagonals & alpha<1
###  	      
###  CALLS TO: nobody
#########################################################################################################

cost.nfields<-function(obs,mod,Q1,S1,alpha){

  rows<-dim(Q1)[1]

  Q1.spam<-(1-alpha)*as.spam(Q1)+(alpha)*diag.spam(rows)  
  
  NEW.Q1<-kronecker(S1, Q1.spam) # defines the NEW.Q1 operator. 
  # When alpha = 1 Q1 is the ID matrix. 
  # When alpha is not 1, actual Q is used. 
  # Kroenicker: basically it's each element of the S matrix * Q, e.g. S11*Q S12*Q...etc. 
    
  cost_detail_mat<-matrix(data = NA,nrow=length(obs),ncol=length(obs)) #initialize empty nfields by nfields matrix. 
  
  for (i in 1:length(obs)){
    for (j in 1:length(obs)){
      #NEW.Q1<-kronecker(S1[i,j],Q1.spam) # each S(ij) in the entire S becomes S(ij)*Q
      FIELDS.Q<-S1[i,j]*Q1.spam # S^-1(ij) becomes S^-1(i,j)*Q. Here we're looking at one S^-1(i,j)*Q at a time. Dim = nobs x nobs
       
      #ans<-FIELDS.Q-NEW.Q1[start[i]:fin[i],start[j]:fin[j]]
      
      #define field 1
      obs_in1<-matrix(t(obs[[i]])) #for i=1 this is 1:19008, which is the first obs vector
      mod_in1<-matrix(t(mod[[i]]))
      obs_in1[is.na(obs_in1)] <- 0
      mod_in1[is.na(mod_in1)] <- 0
      
      #define field 2
      obs_in2<-matrix(t(obs[[j]]))
      mod_in2<-matrix(t(mod[[j]])) 
      obs_in2[is.na(obs_in2)] <- 0
      mod_in2[is.na(mod_in2)] <- 0
      
      cost_detail_mat[i,j]<-t(obs_in1-mod_in1)%*%FIELDS.Q%*%(obs_in2-mod_in2) # Cost = V1t * S[1,1]*Q * V2...
    
    }
  }
  return(cost_detail_mat)
  
}

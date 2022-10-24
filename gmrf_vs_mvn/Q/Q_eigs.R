# Calculate the eigenvalues of the Q matrix. 
# Then use the solver to find the ideal alpha. 


load('Q_5x19.Rdata')

result<-eigen(Q_5x19,only.values=TRUE)
eigs_5x19<-result$values

save(eigs_5x19,file="Q_eigs.Rdata")
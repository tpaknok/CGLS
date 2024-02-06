likelihood.lambda2<-function(lambda,formula,y,x,C){
  require(spaMM)
  lambda <- round(lambda,3)
  V<-diag(diag(C))
  C<-C-V
  n<-nrow(C)
  C.lambda<-(V+lambda*C)
  
  sim_data <- data.frame(y=y,x=x,comm=paste0("comm",1:nrow(comm)))
  rownames(C.lambda) <- colnames(C.lambda) <- paste0("comm",1:nrow(comm))

  m2_spamm <- fitme(formula,corrMatrix=C.lambda,
                    data=sim_data,method = "REML",
                    #init=list(lambda=NaN,phi=NaN),
                    control.HLfit=list(max.iter=10000))
  #logL <- logLik(m2_spamm)
  #cat(lambda,logL,"\n")
  phi <- m2_spamm$phi
  
  #return(-logL) 
  return(phi)
}
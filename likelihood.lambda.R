likelihood.lambda<-function(lambda,y,X,V,comm){
  message(lambda)
  V_lambda <- V*lambda
  diag(V_lambda) <- diag(V)
  
  C.lambda<- get_comm_pair_r_3(comm,V_lambda)
  n <- ncol(C.lambda)
  
  beta<-solve(t(X)%*%solve(C.lambda)%*%X)%*%(t(X)%*%solve(C.lambda)%*%y)
  sig2e<-as.double((1/n)*(t(y-X%*%beta)%*%solve(C.lambda)%*%(y-X%*%beta)))
  logL<--(1/2)*t(y-X%*%beta)%*%solve(sig2e*C.lambda)%*%(y-X%*%beta)-(1/2)*
    determinant(sig2e*C.lambda,logarithm=TRUE)$modulus-(n/2)*log(2*pi)
  
  AIC <- 2*(ncol(as.matrix(X))+2) - 2*logL
  return(AIC) 
}
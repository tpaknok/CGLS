likelihood.lambda.INLA<-function(lambda,formula,data,V,comm,prior="pc.prior.auto",hyper = NULL,fixed=T){
  require(INLA)
  message(lambda, "-INLA")
  V_lambda <- V*lambda
  diag(V_lambda) <- diag(V)
  C.lambda<- get_comm_pair_r_3(comm,V_lambda)
  
  P.lambda <- solve(C.lambda)
  res_f <- gsub("\\bf\\([^)]*\\)", "",Reduce(paste, deparse(formula)))
  res_f <- gsub("\\+\\s*", "", res_f)

  if (prior == "pc.prior.auto") {
  m <- lm(res_f,data=data)
  sdres <- sd(residuals(m))
  hyper_param <- c(3*sdres,0.01)
  } else {
    hyper_param <- hyper
  }
  
  #argus <- c(list(formula = as.formula(formula),
                  #data = data))  
  
  #m_INLA_lambda <- do.call(inla,argus)
  
  formula <- as.formula(gsub("hyper_param","hyper_param",Reduce(paste, deparse(formula)))) #so replace it with a new hyper_param works??

  m_INLA_lambda <- inla(formula,
                          data=data,
                         control.compute = list(dic = TRUE,waic=T,return.marginals.predictor=TRUE),
                         control.inla = list(tolerance = 1e-10))
                         #control.family=list(hyper=list(theta1=list(param=c(1e8),fixed=T))))
  
  m_INLA_lambda <- inla.rerun(m_INLA_lambda)
  #summary(m_INLA_lambda)

  wAIC <- m_INLA_lambda$waic$waic
  return(wAIC)
}

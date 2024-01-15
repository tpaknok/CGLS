sim_CGLS <- function(comm,C,ef_mean,b1=0,signals.X="sr",signals.error=T) {
  #sp <- rnorm(ncol(comm),0,1)
  #ef_mean <- as.matrix(comm)%*%sp
  
  if (signals.X == "sr") {
    x <- rowSums(comm) #species richness
  } 
  
  if (signals.X == "phy_cor") {
    x <- c(t(chol(C)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=1)) #signals in x
  }
  
  if (signals.X == "no_phy_cor") {
    x <- rnorm(nrow(comm),0,1) #no signals in X
  }
  
  if (signals.error == T) {
    ef <- b1*x+c(t(chol(C)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=1)) #signals in error
  } else {
    ef <- b1*x+rnorm(n=nrow(comm),0,sd=1) #no signals in error
  }
  
  sim_data <- data.frame(ef=ef,x=x,comm=1:nrow(comm))
  
  m <- lm(ef~x,data=sim_data)
  p <- summary(m)$coefficients[2,4]
  effect_lm <- summary(m)$coefficients[2,1]
  #type1 <- ifelse(p < 0.05, type1+1,type1)
  
  m2 <- gls(ef~x,data=sim_data,correlation=corSymm(C[lower.tri(C)], fixed = T))
  p_gls <- summary(m2)$tTable[2,4]
  effect_gls <- summary(m2)$tTable[2,1]
  #type1_gls <- ifelse(p_gls < 0.05, type1_gls+1,type1_gls)
  #summary(m2)
  
  ML.opt<-optim(runif(1),likelihood.lambda,y=ef,X=x,C=C,method="L-BFGS-B",
                lower=0.0,upper=1.0)
  V<-diag(diag(C))
  C_temp<-C-V
  n<-nrow(C_temp)
  lambda<-ML.opt$par
  logL<--ML.opt$value
  C.lambda<-(V+lambda*C_temp) 
  
  m2_optim <- gls(ef~x,data=sim_data,correlation=corSymm(C.lambda[lower.tri(C.lambda)], fixed = T))
  p_gls_optim <- 2*pt(abs(summary(m2_optim)$tTable[2,3]),m2_optim$dims$N-m2_optim$dims$p-1,lower.tail=F) #estimated lambda, thus df = n-3. 
  effect_gls_optim <- summary(m2_optim)$tTable[2,1]
  #type1_gls_optim <- ifelse(p_gls_optim < 0.05, type1_gls_optim+1,type1_gls_optim)
  
  return(data.frame(p,p_gls,p_gls_optim,effect_lm,effect_gls,effect_gls_optim)) 
  # m2_spamm <- fitme(ef~x+corrMatrix(1|comm),corrMatrix=C,data=sim_data)
  # m2_spamm <- fitme(ef~x+corrMatrix(1|comm),corrMatrix=C,data=sim_data,control.HLfit=list(sparse_precision=FALSE))
  # p_mm <- summary(m2_spamm,details=T)$beta_table[2,4]
  # effect_mm[i] <- summary(m2_spamm,details=T)$beta_table[2,1]
  # type1_mm <- ifelse(p_mm < 0.05, type1_mm+1,type1_mm)
}

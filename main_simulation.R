comm <- read.csv("gnfurey-e120_soils-0beb647/Data/comm.csv")
data <- read.csv("gnfurey-e120_soils-0beb647/SubmissionData/Dataset_S1.csv")
sp_list <- read.csv("gnfurey-e120_soils-0beb647/Data/species_list.csv")

comm <- comm[match(data$Plot,comm$Plot),12:29]

names <- sp_list[match(colnames(comm),sp_list$Species5),"Species"] 
colnames(comm)[[3]] <- "Amorpha canescens"
colnames(comm)[[11]] <- "Monarda fistulosa"
names[[3]] <- "Amorpha canescens"
names[[11]] <- "Monarda fistulosa"
names[[1]] <- "Achillea millefolium"
names[[2]] <- "Agropyron smithii"
names[[4]] <- "Andropogon gerardii"
names[[7]] <- "Koeleria pyramidata"
names[[11]] <- "Monarda fistulosa"
names[[13]] <- "Dalea purpurea"

colnames(comm) <- names       
comm <- comm[,order(colnames(comm))]
library(rtrees)
library(ape)
#library(taxize)
#taxonomy <- tax_name(names,c("genus"))
#taxonomy <- taxonomy[,c(2,4,3)]
#colnames(taxonomy)[[1]] <- "species" 

plant_names <- sp_list_df(data.frame(species=names),"plant")
tree <- get_tree(sp_list=plant_names,taxon="plant",scenario = "at_basal_node")
plot(tree)
V_sp <- vcv(tree)
cor_V <- vcv(tree,corr=T)
V_sp <- V_sp[order(rownames(V_sp)),order(colnames(V_sp))]
colnames(V_sp) <- rownames(V_sp) <- gsub("_"," ",rownames(V_sp))

###
library(phytools)
library(Matrix)
source("get_comm_pair_r.R") #best version. reduced computational time from 16s to 0.3s
source("likelihood.lambda.R")
source("likelihood.lambda.INLA.R")
source("sim_CGLS.R")

set.seed(1000)

b1 <- 0
n <- 500
optim_n <- 100
optim_seq <- sample(c(rep(TRUE, optim_n), rep(FALSE, n-optim_n)), n ,replace = F)

#prior1 <- list(prior = "pc.prec")
f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,hyper=list(prec = prior1))

sim_results_sr_corE0.5 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=0.5,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_corE0.5<- data.frame(do.call(rbind,sim_results_sr_corE0.5),data="sr_corE0.5")
sim_results_sr_corE <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=1,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_corE <- data.frame(do.call(rbind,sim_results_sr_corE),data="sr_corError")
sim_results_sr_nocorE <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=0,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_nocorE <- data.frame(do.call(rbind,sim_results_sr_nocorE),data="sr_nocorError")

lapply(sim_results_sr_corE, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_corE0.5, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_nocorE, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_corE, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_sr_corE0.5, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_sr_nocorE, function(x) range(x,na.rm=T))[8:14]

##################
b1 <- 0.5
set.seed(1000)
n <- 500
optim_n <- 100
optim_seq <- sample(c(rep(TRUE, optim_n), rep(FALSE, n-optim_n)), n ,replace = F)

sim_results_sr_corE0.5_2 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=0.5,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_corE0.5_2<- data.frame(do.call(rbind,sim_results_sr_corE0.5_2),data="sr_corE0.5")
sim_results_sr_corE_2 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=1,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_corE_2 <- data.frame(do.call(rbind,sim_results_sr_corE_2),data="sr_corError")
sim_results_sr_nocorE_2 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="sr",signals_intercept=T,signals_slope=F,lambda_true=0,true_model=T,optim_model=x,INLA_formula = f))
sim_results_sr_nocorE_2 <- data.frame(do.call(rbind,sim_results_sr_nocorE_2),data="sr_nocorError")

lapply(sim_results_sr_corE_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_corE0.5_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_nocorE_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_sr_corE_2, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_sr_corE0.5_2, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_sr_nocorE_2, function(x) range(x,na.rm=T))[8:14]
###################
b1 <- 0
n <- 500
optim_n <- 100
optim_seq <- sample(c(rep(TRUE, optim_n), rep(FALSE, n-optim_n)), n ,replace = F)

set.seed(1000)

#f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)+f(comm2,x,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)
f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)+f(comm2,x,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)
sim_results_slope_sr_corE0.5 <-lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=0.5,true_model=T,optim_model=x,
                                                                     INLA_formula = f))
sim_results_slope_sr_corE0.5<- data.frame(do.call(rbind,sim_results_slope_sr_corE0.5),data="phy_cor_corE0.5")
sim_results_slope_sr_corE <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=1,true_model=T,optim_model=x,
                                                                   INLA_formula = f))
sim_results_slope_sr_corE <- data.frame(do.call(rbind,sim_results_slope_sr_corE),data="phy_cor_corError")
sim_results_slope_sr_nocorE <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=0,true_model=T,optim_model=x,
                                                                     INLA_formula = f))
sim_results_slope_sr_nocorE <- data.frame(do.call(rbind,sim_results_slope_sr_nocorE),data="phy_cor_nocorError")

lapply(sim_results_slope_sr_corE, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_corE0.5, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_nocorE, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_corE, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_slope_sr_corE0.5, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_slope_sr_nocorE, function(x) range(x,na.rm=T))[8:14]

###################
b1 <- 0.5
n <- 500
optim_n <- 100
optim_seq <- sample(c(rep(TRUE, optim_n), rep(FALSE, n-optim_n)), n ,replace = F)

set.seed(1000)

#f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)+f(comm2,x,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)
f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)+f(comm2,x,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)
sim_results_slope_sr_corE0.5_2 <-lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=0.5,true_model=T,optim_model=x,
                                                                     INLA_formula = f))
sim_results_slope_sr_corE0.5_2<- data.frame(do.call(rbind,sim_results_slope_sr_corE0.5_2),data="phy_cor_corE0.5")
sim_results_slope_sr_corE_2 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=1,true_model=T,optim_model=x,
                                                                   INLA_formula = f))
sim_results_slope_sr_corE_2 <- data.frame(do.call(rbind,sim_results_slope_sr_corE_2),data="phy_cor_corError")
sim_results_slope_sr_nocorE_2 <- lapply(optim_seq,function(x) sim_CGLS(comm,V_sp,0,1,b1=b1,signals_X="phy_cor",signals_intercept=T,signals_slope=T,lambda_true=0,true_model=T,optim_model=x,
                                                                     INLA_formula = f))
sim_results_slope_sr_nocorE_2 <- data.frame(do.call(rbind,sim_results_slope_sr_nocorE_2),data="phy_cor_nocorError")

lapply(sim_results_slope_sr_corE_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_corE0.5_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_nocorE_2, function(x) sum(x<0.05,na.rm=T)/length(x[!is.na(x)]))[1:7]
lapply(sim_results_slope_sr_corE_2, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_slope_sr_corE0.5_2, function(x) range(x,na.rm=T))[8:14]
lapply(sim_results_slope_sr_nocorE_2, function(x) range(x,na.rm=T))[8:14]

#####
####################
lapply(sim_results_cor[1:5],function(x) sum(x < 0.05))
lapply(sim_results_no_cor[1:5],function(x) sum(x < 0.05))
lapply(sim_results_corX_nocorE[1:5],function(x) sum(x < 0.05))
lapply(sim_results_nocorX_corE[1:5],function(x) sum(x < 0.05))

###

library(ggplot2)

p_sim <- ggplot(data=combined_sim_results,aes(y=effect_gls_optim,x=effect_lm))+
  geom_point()+
  geom_abline()+
  xlim(-2,2)+
  ylim(-2,2)+
  labs(y="Optimized GLS estimates",x="OLS estimates")+
  facet_wrap(~data)
plot(p_sim)


plot(effect_lm,effect_gls)
plot(effect_lm,effect_gls_optim)
t.test(abs(effect_lm),abs(effect_gls_optim),paired=T)

### empirical
library(tidyverse)
library(nlme)
library(phylosignal)
cor(data[,c(22:29)])
soil <- data[,c(22:29)]
colnames(soil) <- word(colnames(soil),4,sep="_")
C <- get_comm_pair_r_3(comm,V_sp)

library(AICcmodavg)
f <- y~x+f(comm,model="generic0",Cmatrix=P.lambda,prior = "pc.prec",hyper=hyper_param)

soil_empiricial <- function(soil_data,x,V,names,soil_properties) {
  print(soil_properties)
  y <- soil_data
  df <- data.frame(y=y,x=x)
  df$comm <- df$comm2 <- 1:nrow(df)
  m <- gls(y~scale(x),data=df)
  summary(m)
  m_sig <- summary(m)$tTable[2,4]
  #m2 <- gls(y~log(NumSp),data=data,correlation=corSymm(C[lower.tri(C)], fixed = T))
  #summary(m2)
  init_lambda <- lambdaTest(df$y,C)$Lambda
  
  ML.opt<-optim(runif(1),likelihood.lambda.INLA,
                formula= f,
                data=df,V=V_sp,comm=comm,method="L-BFGS-B",
                lower=0.0,upper=1.0,control=list(factr=1e7,pgtol=1e-9))

  lambda_INLA<-ML.opt$par
  wAIC <- ML.opt$value
  V_INLA <- V*lambda_INLA
  diag(V_INLA) <- diag(V)
  C.lambda.INLA<- get_comm_pair_r_3(comm,V_INLA) 
  
  #m2_spamm_optim <- fitme(#y~scale(X)+corrMatrix(1|comm)+corrMatrix(0+scale(X)|comm),
                         # y~scale(X)+corrMatrix(1|comm),
                         # corrMatrix=C.lambda,
                         # data=df, method="REML",
                         # init=list(lambda=NaN,phi=NaN),
                         # control.HLfit=list(max.iter=100000))
  #m2_optim <- gls(y~log(X),correlation=corSymm(C.lambda[lower.tri(C.lambda)], fixed = T),data=df)
  #summary(m2_optim)
  #m2_optim_sig <- 2*pt(abs(fixef(m2_spamm_optim)/sqrt(diag(vcov( m2_spamm_optim)))),df.residual(m2_spamm_optim),lower.tail=F)[2]
  newobs <- data.frame(y=NA,x=unique(x),comm=NA,comm2=NA)
  INLA_df <- rbind(df,newobs)
  
  m2_INLA_LM <- inla(y~x,data=INLA_df,control.compute = list(dic=T,waic=T))
  m2_INLA_LM <- inla.rerun(m2_INLA_LM)
  wAIC_LM <- m2_INLA_LM$waic$waic
  
  prec.mat.INLA <- solve(C.lambda.INLA)
  m2_INLA_optim <- inla(y~x+f(comm,model="generic0",Cmatrix=prec.mat.INLA,prior = "pc.prec",hyper= c(3*sd(residuals(m)),0.01)),data=INLA_df,control.predictor=list(compute=TRUE),control.compute = list(dic=T,waic=T),safe=T,control.inla = list(tolerance = 1e-10))
  m2_INLA_optim <- inla.rerun(m2_INLA_optim) # https://groups.google.com/g/r-inla-discussion-group/c/qkoV9ZtA1Wo
  m2_INLA_optim <- inla.rerun(m2_INLA_optim) # https://groups.google.com/g/r-inla-discussion-group/c/qkoV9ZtA1Wo
  
  summary(m2_INLA_optim)
  m2_INLA_optim_predict <- m2_INLA_optim$summary.fitted.values[-1:-(nrow(df)),]
  m2_INLA_optim_predict <- data.frame(newobs$x,m2_INLA_optim_predict[,c("mean","0.025quant","0.975quant")])
  m2_INLA_optim_predict$sig <- ifelse(sign(m2_INLA_optim$summary.fixed)[2,3] == sign(m2_INLA_optim$summary.fixed)[2,5],-1,1)
  
  library(ggeffects)
  m_predict <- data.frame(ggeffect(m,terms="x"),sig=m_sig)

  m_predict <- m_predict[,c(1,2,4,5,7)]
  #m2_optim_predict <- data.frame(ggeffect(m2_spamm_optim,terms="X"),sig=m2_optim_sig)
  #m2_optim_predict <- as.data.frame(predict(m2_spamm_optim,newdata=data.frame(X=seq(min(df$X),max(df$X),length.out=10),comm=1),re.form=NA))
  #m2_optim_predict <- cbind(X=seq(min(df$X),max(df$X),length.out=10),m2_optim_predict,get_intervals(m2_spamm_optim,data.frame(X=seq(min(df$X),max(df$X),length.out=10),comm=1),re.form=NA,intervals="fixefVar"),sig= m2_optim_sig)
  #colnames(m2_optim_predict) <- colnames(m_predict)
  
  colnames(m2_INLA_optim_predict) <- colnames(m_predict)
  #m_predict <- data.frame(predictSE(m,newdata=data.frame(X=1:16)))
  #m_predict$conf.low <- m_predict$fit-m_predict$se.fit*1.96
  #m_predict$conf.high <- m_predict$fit+m_predict$se.fit*1.96
  #m_predict$sig <- m_sig
  
 #m2_optim_predict <- data.frame(predictSE(m2_optim,newdata=data.frame(X=1:16)))
 #m2_optim_predict$conf.low <- m2_optim_predict$fit-m2_optim_predict$se.fit*1.96
  #m2_optim_predict$conf.high <- m2_optim_predict$fit+m2_optim_predict$se.fit*1.96
  #m2_optim_predict$sig <- m2_optim_sig
  
  predict_df <- rbind(cbind(m_predict,model="OLS",soil=soil_properties,lambda=NA,hyper_Gaussian = NA, hyper_Comm = NA,wAIC = NA,wAIC_LM=NA),
                      cbind(m2_INLA_optim_predict,model="Optimized INLA",soil=soil_properties,lambda=lambda_INLA,hyper_Gaussian=m2_INLA_optim$summary.hyperpar[1,1],hyper_Comm =m2_INLA_optim$summary.hyperpar[2,1],wAIC=wAIC,wAIC_LM = wAIC_LM))
  
  return(predict_df)
}

df <- data.frame(x=log(data$NumSp))
predict_df <- lapply(1:ncol(soil),function(x) soil_empiricial(soil_data=soil[,x],x=log(data$NumSp),V=V_sp,soil_properties=colnames(soil)[x]))

predict_df <- do.call(rbind,predict_df)
predict_df$x <- exp(predict_df$x)
predict_df$sig_binary <- ifelse(predict_df$sig < 0.05,"Significant","Insignificant")
lapply(1:ncol(comm),function(x) table(rowSums(comm[comm[,x]>0,])))
### plot
library(ggplot2)
library(tidyverse)
plot_data <- data %>% 
  select(starts_with("gm2_2017") | contains("NumSp")) %>%
  #select(contains("2015.2017") | contains("NumSp")) %>%
  pivot_longer(!NumSp,names_to="soil")

plot_data$soil <- word(plot_data$soil,4,sep="_")

plot_data$soil <- fct_relevel(plot_data$soil,"Carbon","Nitrogen","Potassium","Calcium","Magnesium","CEC","pH","Phosphorus")
predict_df$soil <- fct_relevel(predict_df$soil,"Carbon","Nitrogen","Potassium","Calcium","Magnesium","CEC","pH","Phosphorus")
predict_df$ratio <- predict_df$hyper_Gaussian/predict_df$hyper_Comm
predict_df$ratio <- (1/predict_df$hyper_Gaussian)/(1/predict_df$hyper_Comm)

p <- ggplot(data=predict_df)+
  geom_jitter(data=plot_data,aes(y=value,x=NumSp),width=0.25)+
  geom_line(aes(y=predicted,x=x,colour=model,linetype=sig_binary))+
  geom_ribbon(aes(y=predicted,x=x,fill=model,ymin=conf.low,ymax=conf.high),alpha=0.3,colour="transparent")+
  facet_wrap(~soil,scales="free",nrow=2,ncol=4)+
  theme_classic()+
  labs(x="Species richness",y="Soil properties",colour="Model",fill="Model",linetype="Significance")+
  scale_linetype_manual(values=c(2,1))+
  scale_colour_manual(values=c("#000000","#E69F00"))+
  scale_fill_manual(values=c("#000000","#E69F00"))+
  theme(legend.position="bottom")

plot(p)

ggsave("Figure/soil.tiff",height=8.4*1.5,width=16.8,units="cm",compression="lzw",dpi=600)

vcv2phylo(C_orig)
######### useless = assume no evolutinoary history
C <- matrix(0,nrow=nrow(comm),ncol=nrow(comm))

for (j in 1:nrow(comm)){
  for (k in 1:nrow(comm)) {
    comm1 <- comm[j,]
    comm2 <- comm[k,]
    
    species <- specnumber(comm1)*specnumber(comm2)
    overlap <- sum(comm1+comm2 ==2)
    
    cor <- overlap/sqrt(species)
    
    C[j,k] <- cor
  }
}

C <- as.matrix(nearPD(C,corr=T,keepDiag=T)$mat)

colnames(C) <- rownames(C) <- 1:nrow(comm)
#C <- corSymm(C[lower.tri(C)], fixed = T)

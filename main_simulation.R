comm <- read.csv("C:/Users/pakno/OneDrive/Desktop/BEF/gnfurey-e120_soils-0beb647/Data/comm.csv")
data <- read.csv("C:/Users/pakno/OneDrive/Desktop/BEF/gnfurey-e120_soils-0beb647/SubmissionData/Dataset_S1.csv")
sp_list <- read.csv("C:/Users/pakno/OneDrive/Desktop/BEF/gnfurey-e120_soils-0beb647/Data/species_list.csv")

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
source("get_comm_pair_r.R")
source("likelihood.lambda.R")

b1 <- 0
set.seed(999)
C_orig <- get_comm_pair_r(comm,V_sp)
C_orig[upper.tri(C_orig)] <- t(C_orig)[upper.tri(C_orig)]
C <- C_orig
C <- as.matrix(nearPD(C_orig,corr=T,keepDiag=T,maxit=100000)$mat)

#type1 <- type1_mm <- type1_gls <- type1_gls_optim <- 0
#effect_lm <- effect_gls <- effect_mm <- effect_gls_optim <- 0

b1 <- 0.25
n <- 1000
sim_results_no_cor <- lapply(1:n,function(x) sim_CGLS(comm,C,0,b1=b1,signals.X="no_phy_cor",signals.error=F))
sim_results_no_cor <- data.frame(do.call(rbind,sim_results_no_cor),data="nocorX_nocorError")
sim_results_cor <- lapply(1:n,function(x) sim_CGLS(comm,C,0,b1=b1,signals.X="phy_cor",signals.error=T))
sim_results_cor <- data.frame(do.call(rbind,sim_results_cor),data="corX_corError")
sim_results_corX_nocorE <- lapply(1:n,function(x) sim_CGLS(comm,C,0,b1=b1,signals.X="phy_cor",signals.error=F))
sim_results_corX_nocorE <- data.frame(do.call(rbind,sim_results_corX_nocorE),data="corX_nocorError")
sim_results_nocorX_corE <- lapply(1:n,function(x) sim_CGLS(comm,C,0,b1=b1,signals.X="no_phy_cor",signals.error=T))
sim_results_nocorX_corE <- data.frame(do.call(rbind,sim_results_nocorX_corE),data="nocorX_corError")
combined_sim_results <- rbind(sim_results_no_cor,sim_results_cor,sim_results_corX_nocorE,sim_results_nocorX_corE)

lapply(sim_results_cor[1:3],function(x) sum(x < 0.05))
lapply(sim_results_no_cor[1:3],function(x) sum(x < 0.05))
lapply(sim_results_corX_nocorE[1:3],function(x) sum(x < 0.05))
lapply(sim_results_nocorX_corE[1:3],function(x) sum(x < 0.05))

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

cor(data[,c(22:29)])

y <- data$AbvBioAnnProd.2015.2017
m <- gls(y~log(NumSp),data=data)
summary(m)

m2 <- gls(y~log(NumSp),data=data,correlation=corSymm(C[lower.tri(C)], fixed = T))
summary(m2)

ML.opt<-optim(runif(1),likelihood.lambda,y=y,X=log(data$NumSp),C=C,method="L-BFGS-B",
              lower=0.0,upper=1.0)
V<-diag(diag(C))
C_temp<-C-V
n<-nrow(C_temp)
lambda<-ML.opt$par
logL<--ML.opt$value
C.lambda<-(V+lambda*C_temp) 

m2_optim <- gls(y~log(NumSp),data=data,correlation=corSymm(C.lambda[lower.tri(C.lambda)], fixed = T))
summary(m2_optim)
p_gls_optim <- 2*pt(summary(m2_optim)$tTable[2,3],m2_optim$dims$N-m2_optim$dims$p-1,lower.tail=F) #estimated lambda, thus df = n-3. 
effect_gls_optim[i] <- summary(m2_optim)$tTable[2,1]
type1_gls_optim <- ifelse(p_gls_optim < 0.05, type1_gls_optim+1,type1_gls_optim)

library(ggeffects)
m2_predict <- ggeffect(m2,terms=c("NumSp[1:16,by=1]"))
m2_optim_predict <- ggeffect(m2_optim,terms=c("NumSp[1:16,by=1]"))
m_predict <- ggeffect(m,terms=c("NumSp[1:16,by=1]"))

predict_df <- rbind(cbind(m_predict,model="lm"),cbind(m2_predict,model="gls"),cbind(m2_optim_predict,model="gls_optim"))

library(ggplot2)

p <- ggplot(data=predict_df)+
  geom_point(data=data,aes(y=y,x=NumSp))+
  geom_line(aes(y=predicted,x=x,colour=model))+
  geom_ribbon(aes(y=predicted,x=x,fill=model,ymin=conf.low,ymax=conf.high),alpha=0.1,colour="transparent")+
  theme_classic()

plot(p)

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

get_comm_pair_r <- function (comm,V) {
  
  pairwise_comm_r <- function (pairwise_comm,V) {
    pairwise_comm <- pairwise_comm[colSums(pairwise_comm) > 0]
    V <- V[rownames(V) %in% colnames(pairwise_comm),colnames(V) %in% colnames(pairwise_comm)]
    
    if (length(V) == 1) {
      r <- 1
    } else {
      pairwise_comm <- as(pairwise_comm,"sparseMatrix") #much much much quicker!
      list_df <- list(cov_component=expand.grid(C1=1:ncol(pairwise_comm),C2=1:ncol(pairwise_comm)),
                      C1=expand.grid(C1=1:ncol(pairwise_comm),C2=1:ncol(pairwise_comm)),
                      C2=expand.grid(C1=1:ncol(pairwise_comm),C2=1:ncol(pairwise_comm)))
      
      cov <- lapply(list_df, function(x) x$cov <- V[as.matrix(x)])
      numerator <- sum(cov[[1]]*t(pairwise_comm[1,list_df[[1]]$C1]*pairwise_comm[2,list_df[[1]]$C2]))
      denom1 <- sum(cov[[2]]*t(pairwise_comm[1,list_df[[2]]$C1]*pairwise_comm[1,list_df[[2]]$C2]))
      denom2 <- sum(cov[[3]]*t(pairwise_comm[2,list_df[[3]]$C1]*pairwise_comm[2,list_df[[3]]$C2]))
      r <- numerator/sqrt(denom1*denom2)
    }
    return(r)
  }
  
  require(nlme)
  pairs <- expand.grid(C1 = 1:nrow(comm), C2= 1:nrow(comm))
  pairs <- subset(pairs,C1>=C2)
  pairwise_comm <- split(pairs,seq(nrow(pairs)))
  
  pairs$r <- sapply(pairwise_comm, function(x) pairwise_comm_r(comm[unlist(x),],V))
  
  cor_M <- matrix(NA,nrow(comm),nrow(comm))
  cor_M[as.matrix(pairs[,1:2])] <- pairs$r
  
  return(cor_M)
}

#rm(list=ls())


#EXECUTER MRM_CRITERE_IS_MCMC_STAGE.R AU PREALABLE



#TESTS
options(warn=-1)




N.calls <- 210
REP=20


PROBAS = list(list(p= 5*1e-2,ordre.p=2), list(p= 5*1e-3,ordre.p=3), list(p= 5*1e-4,ordre.p=4))
DIMENSIONS = c(2, 4, 6, 8, 10, 12, 15)


EXAMPLE = "EXAMPLE_1"

df <- data.frame(EXAMPLE=character(), p=double(), ordre.p=integer(), DIM=double(), MRM_time=double(), MRM_um=double(), MRM_uM=double(), MRM_mcmc_time=double(), MRM_mcmc_um=double(), MRM_mcmc_uM=double())


for (k in 1:length(PROBAS)) {
  
  p=PROBAS[[k]]$p
  ordre.p=PROBAS[[k]]$ordre.p
  
  for (i in DIMENSIONS) {
    nDim=i
    
    tmp<-rep(0.0, REP)
    tmp2<- rep(0.0, REP)
    um<-rep(0.0, REP)
    uM<-rep(0.0, REP)
    um2<-rep(0.0, REP)
    uM2<-rep(0.0, REP)
    
    alpha_start= (1/nDim)*(2.38)^2         #PEUT ETRE MODIFIE 
    
    
    choix.loi = list()
    InputDist <- list()
    dir.monotony = c(1,rep(-1, nDim-1))
    for(i in 1:nDim){
      choix.loi[[i]] = list("gamma",c(i+1,1))
    }
    InputDist <- choix.loi
    for(i in 1:nDim){
      nparam <- length(choix.loi[[i]][[2]])
      for(j in 1:nparam){
        InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
      }
    }
    f  <- function(X){
      result =  (X[1]/sum(X)) - qbeta((p),2,sum(c(2:nDim)+1))
      return(result)
    }
    
    j <- 1
    while(j <= REP) {
      
      an.error.occured <- FALSE
      time1 <- Sys.time()
      
      
      res.MRM <- MRM(f, nDim, choix.loi, dir.monotony, N.calls, Method = "MRM", code ="C", ordre.p = ordre.p, MAXIMIN="MLE")
      browser()
      
      tryCatch( { res.MRM <- MRM(f, nDim, choix.loi, dir.monotony, N.calls, Method = "MRM", code ="C", ordre.p = ordre.p, MAXIMIN="MLE") }
                , error = function(e) {an.error.occured <<- TRUE ; print("ERROR")})
      if (an.error.occured) {break}
  
      time2 <- Sys.time()
      tmp[j] <- as.double(difftime(time2, time1, tz, units = "secs"))
      um[j] <- res.MRM$tab[dim(res.MRM$tab)[1], "Um" ]/p
      uM[j] <- res.MRM$tab[dim(res.MRM$tab)[1], "UM" ]/p
      
      
      
      time3 <- Sys.time()
      tryCatch( { res.MRM2 <- MRM(f, nDim, choix.loi = choix.loi, dir.monot = dir.monotony, N.calls, Method = "MRM_MCMC2", code ="C", ordre.p = ordre.p,  alpha_start = alpha_start) }
                , error = function(e) {an.error.occured <<- TRUE ; print("ERROR")})
      if (an.error.occured) {break}
      time4 <- Sys.time()
      tmp2[j] <- as.double(difftime(time4, time3, tz, units = "secs"))
      tmp[j] <- as.double(difftime(time2, time1, tz, units = "secs"))
      um2[j] <- res.MRM2$tab[dim(res.MRM2$tab)[1], "Um" ]/p
      uM2[j] <- res.MRM2$tab[dim(res.MRM2$tab)[1], "UM" ]/p
      
      j <-j+1
      
    }
    
    de <- list(EXAMPLE=EXAMPLE, p=p, ordre.p=ordre.p, DIM=nDim, MRM_time=mean(tmp), MRM_um=mean(um), MRM_uM=mean(uM), MRM_mcmc_time=mean(tmp2), MRM_mcmc_um=mean(um2), MRM_mcmc_uM=mean(uM2) )
    df = rbind(df,de, stringsAsFactors=FALSE)
    
  }
  
  
}





EXAMPLE = "EXAMPLE_2"


df2 <- data.frame(EXAMPLE=character(), p=double(), ordre.p=integer(), DIM=double(), MRM_time=double(), MRM_um=double(), MRM_uM=double(), MRM_mcmc_time=double(), MRM_mcmc_um=double(), MRM_mcmc_uM=double())


for (k in 1:length(PROBAS)) {
  
  p=PROBAS[[k]]$p
  ordre.p=PROBAS[[k]]$ordre.p
  
  for (i in DIMENSIONS) {
    nDim=i
    tmp<-rep(0.0, REP)
    tmp2<- rep(0.0, REP)
    um<-rep(0.0, REP)
    uM<-rep(0.0, REP)
    um2<-rep(0.0, REP)
    uM2<-rep(0.0, REP)
    
    alpha_start= (1/nDim)*(2.38)^2         #PEUT ETRE MODIFIE 
    
    choix.loi = list()
    InputDist <- list()
    dir.monotony = rep(c(-1, 1), nDim%/%2)
    if (nDim%%2==1) {dir.monotony=c(dir.monotony, -1)}
    for(i in 1:nDim){
      choix.loi[[i]] = list("norm", c(0, sqrt(i)) ) 
    }
    InputDist <- choix.loi
    for(i in 1:nDim){
      nparam <- length(choix.loi[[i]][[2]])
      for(j in 1:nparam){
        InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
      }
    }
    f  <- function(X){
      result =  sum(dir.monotony*X)-qnorm(p, mean = 0 , sd = sqrt((nDim+1)*nDim*0.5))
      return(result)
    }
    
    
    j <- 1
    while(j <= REP) {
      
      an.error.occured <- FALSE
      time1 <- Sys.time()
      tryCatch( { res.MRM <- MRM(f, nDim, choix.loi, dir.monotony, N.calls, Method = "MRM", code ="C", ordre.p = ordre.p, MAXIMIN="MLE") }
                , error = function(e) {an.error.occured <<- TRUE ; print("ERROR")})
      if (an.error.occured) {break}
      
      time2 <- Sys.time()
      tmp[j] <- as.double(difftime(time2, time1, tz, units = "secs"))
      um[j] <- res.MRM$tab[dim(res.MRM$tab)[1], "Um" ]/p
      uM[j] <- res.MRM$tab[dim(res.MRM$tab)[1], "UM" ]/p
      
      
      time3 <- Sys.time()
      tryCatch( { res.MRM2 <- MRM(f, nDim, choix.loi = choix.loi, dir.monot = dir.monotony, N.calls, Method = "MRM_MCMC2", code ="C", ordre.p = ordre.p,  alpha_start = alpha_start) }
                , error = function(e) {an.error.occured <<- TRUE ; print("ERROR")})
      if (an.error.occured) {break}
      time4 <- Sys.time()
      tmp2[j] <- as.double(difftime(time4, time3, tz, units = "secs"))
      tmp[j] <- as.double(difftime(time2, time1, tz, units = "secs"))
      um2[j] <- res.MRM2$tab[dim(res.MRM2$tab)[1], "Um" ]/p
      uM2[j] <- res.MRM2$tab[dim(res.MRM2$tab)[1], "UM" ]/p
      
      j <-j+1
      
    }
    
    de <- list(EXAMPLE=EXAMPLE, p=p, ordre.p=ordre.p, DIM=nDim, MRM_time=mean(tmp), MRM_um=mean(um), MRM_uM=mean(uM), MRM_mcmc_time=mean(tmp2), MRM_mcmc_um=mean(um2), MRM_mcmc_uM=mean(uM2) )
    df2 = rbind(df2,de, stringsAsFactors=FALSE)
    
  }
  
  
}







# SAUVEGARDE LES 2 DATAFRAMES CREES DANS LE WORKING DIRECTORY
save(df, df2, file = "My_Two_Objects.RData") 




#load("My_Two_Objects.RData")





#Il faut installer le package numDeriv!!

#-------------------------------------------------------#
#--------------------Fonctions test :-------------------# 
#-------------------------------------------------------#

#-----------#
#------I----#
#-----------#

#true.value = 1e-14
#nDim = 2

#N = 1e4
#set.seed(1986)
#X=cbind(runif(N),runif(N))

Method = "HLRF"
IS = TRUE
N.appels = 300
q = 100/N.appels

#u.dep=c(rep(0,nDim-1),0.01)
u.dep = rep(0, nDim)
choix.loi = list()
dir.monotony = c(1,rep(-1, nDim-1))

for(i in 1:nDim){
    choix.loi[[i]] = list("gamma",c(i+1,1))
}

f  <- function(X){
 result =  (X[1]/sum(X)) - qbeta((true.value),2,sum(c(2:nDim)+1))
 return(result)
}

#res=FORM(f,u.dep,choix.loi,N.appels,Method,IS,q);res



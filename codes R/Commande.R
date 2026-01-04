nDim = 5
ordre.p = 2
true.value = 10^(-ordre.p)
BUDGET = 150
source("/home/C23072/Documents/MRM/MRM_Critere/fonction_test_1.R")
source("/home/C23072/Documents/MRM/MRM_Critere/MRM_Critere.R")

res.MLE=MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "MRM", code = "C", ordre.p, MAXIMIN ="MLE")


#res.subset=MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "M-Subset", code = "C", ordre.p, MAXIMIN ="essai")


#res.sort=MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "MRM", code = "C", ordre.p, MAXIMIN ="essai")
#res.sort=MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "M-Subset", code = "C", ordre.p, MAXIMIN ="essai")

B <- dim(res.MLE)[1]


#res.MLE[B, ]
#res.subset[dim(res.subset)[1], ]



NN <-50

M <- res.MLE[NN, 2]

x.M <- B

plot(1:B, res.MLE[,1], type='l',col="black",ylim=c(0, nDim*true.value),lwd=2, xlim=c(0, x.M))
lines(1:B, res.MLE[,2],col="black",lwd=2)
lines(1:B, res.MLE[,3],col="blue",lty=2,lwd=2)
lines(1:B, res.MLE[,4],col="blue",lty=1,lwd=2)
lines(1:B, res.MLE[,5],col="blue",lty=1,lwd=2)
lines(1:B, res.MLE[,7],col="red",lty=2,lwd=2)
lines(1:B, res.MLE[,8],col="red",lty=1,lwd=2)
lines(1:B, res.MLE[,9],col="red",lty=1,lwd=2)
lines(1:B, rep(true.value, B), col="orange",lwd=2)



##########################################################################################################
##########################################################################################################
##########################################################################################################

nDim = 6
ordre.p = 3
true.value = 10^(-ordre.p)
BUDGET = 150
source("/home/C23072/Documents/MRM/MRM_Critere/fonction_test_1.R")
source("/home/C23072/Documents/MRM/MRM_Critere/MRM_Critere.R")


N <- 1780


RES = MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "MRM", code = "C", ordre.p, MAXIMIN ="MLE")

B <- dim(RES)[1]

RESULT <- matrix(0, nrow = N, ncol = 8)

RESULT[1,] <- RES[B, c(1,2,3,4,5,7,8,9)]


for(i in (N+1):(N+180)){
  print(i);flush.console();
  res.MLE=MRM(f, nDim, choix.loi, dir.monotony,BUDGET, Method = "MRM", code = "C", ordre.p, MAXIMIN ="MLE")
  RESULT[i, ] <- res.MLE[B, c(1,2,3,4,5,7,8,9)]
  RES <- RES + res.MLE
}

Estim.old <- RESULT[, 3]
Estim.new <- RESULT[, 6]

DENS.old <- density(Estim.old)
DENS.new <- density(Estim.new)

quantile.old <- quantile(Estim.old, c(0.025,0.975))
quantile.new <- quantile(Estim.new, c(0.025,0.975))

par(mfrow = c(2,2))

plot(1:B, RES[,1]/N, type='l',col="black",lwd=2,ylim=c(0, 0.0015))
lines(1:B, RES[,2]/N,col="black",lwd=2)
lines(1:B, RES[,3]/N,col="blue",lty=2,lwd=2)
lines(1:B, RES[,4]/N,col="blue",lty=1,lwd=2)
lines(1:B, RES[,5]/N,col="blue",lty=1,lwd=2)
lines(1:B, RES[,7]/N,col="red",lty=2,lwd=2)
#lines(1:B, RES[,8]/N,col="red",lty=1,lwd=2)
#lines(1:B, RES[,9]/N,col="red",lty=1,lwd=2)
lines(1:B, rep(true.value, B), col="orange",lwd=2)

plot(1:N , abs(cumsum(Estim.new)/(1:N) - true.value) ,type='l',col="red" , log='y')
lines(1:N, abs(cumsum(Estim.old)/(1:N) - true.value),type='l',col="blue"  )

xm <- min(c(Estim.old ,Estim.new) )
xM <- max(c(Estim.old ,Estim.new) )

xx <- seq(from = 0, to = 0.003,length=1e3)

yy.old <- dnorm(xx, mean(Estim.old), sd = sd(Estim.old) )
yy.new <- dnorm(xx, mean(Estim.new), sd = sd(Estim.new) )

Hist.old <- hist(Estim.old,xlim=c(xm,xM),freq = FALSE,ylim=c(0,max(c(DENS.old[[2]] ,yy.old, 900))))
yM.old <- max(Hist.old$density)
lines( rep(true.value, 1e3), seq(from = 0, to = yM.old,length = 1e3),col="orange",lwd=2)

lines( rep(quantile.old[1], 1e3), seq(from = 0, to = yM.old,length = 1e3),col="blue",lwd=2,lty=1)
lines( rep(quantile.old[2], 1e3), seq(from = 0, to = yM.old,length = 1e3),col="blue",lwd=2,lty=1)
lines( rep(mean(RESULT[,3]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="blue",lwd=2,lty=2)
lines( rep(mean(RESULT[,4]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="green",lwd=2,lty=2)
lines( rep(mean(RESULT[,5]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="green",lwd=2,lty=2)
lines(xx,yy.old,col="blue")
lines(DENS.old,col="blue",lty=2,lwd=2)



Hist.new <- hist(Estim.new,xlim=c(xm,xM),freq = FALSE,ylim=c(0,max(c(DENS.new[[2]] ,yy.old, 900))))
yM.new <- max(Hist.new$density)
lines( rep(true.value, 1e3), seq(from = 0, to = yM.new,length = 1e3),col="orange",lwd=2)
lines( rep(quantile.new[1], 1e3), seq(from = 0, to = yM.new,length = 1e3),col="red",lwd=2,lty=1)
lines( rep(quantile.new[2], 1e3), seq(from = 0, to = yM.new,length = 1e3),col="red",lwd=2,lty=1)
lines( rep(mean(RESULT[,6]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="red",lwd=2,lty=2)
lines( rep(mean(RESULT[,7]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="green",lwd=2,lty=2)
lines( rep(mean(RESULT[,8]), 1e3), seq(from = 0, to = yM.old,length = 1e3),col="green",lwd=2,lty=2)
lines(DENS.new,col="red",lty=2,lwd=2)
lines(xx, yy.new,col="red")

var.old <- mean( (Estim.old - true.value )^2 )
var.new <- mean( (Estim.new - true.value )^2 )

var.old; var.new



#-----------------------------------------------------------------------------------------------------------------#
#
#
#
#                                                  Algorithme MRM (MLE/IS/IS-MCMC)  :
#
#
#
#-----------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------------------#
# Input :											    	    
#-----------------------------------------------------------------------------------------------------------------#
# f : Fonction de d?faillance.								    
# ndim : Dimension de l'espace
# choix.loi : liste contenant le nom des lois d'entr?es et leur(s) param?tre(s).  
# dir.monot : vecteur appartenant ? {-1,1}^ndim. dirmonot[i] = -1 : code d?croissant ; 1 : code croissant en X[i]
# N.appels : Budget d'appel au code donn? par l'utilisateur
# Code = {"R", "C"}, choix du langage de programmation pour certaine fonction
# ordre.p : Ordre de la proba cherch?e : si la proba est de l'ordre de 1e-k, ordre.p = k.
#           On r?servera alors floor((ordre.p+2)*log(10)/log(2)) appels pour la recherche par dichotomie
# BUDGET : Nombre d'appels au code autoris? par l'utilisateur
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Output :	Renvoi une matrice de taille BUDGETx6 : 
#-----------------------------------------------------------------------------------------------------------------#
#  um = bornes inf?rieures 
#  uM = bornes sup?rieures 
#  MLE = Estimateur du maximum de vraissemblance
#  u.m.Fisher = Bornes inf?rieures de l'intervalle de confiance li? au MLE
#  u.M.Fisher = Bornes sup?rieures de l'intervalle de confiance li? au MLE
#  CV = coefficient de variation de l'estimateur MLE
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

#dyn.load("C:/Users/loicb/Desktop/MRM/Code_C/Domination/IS_DOMINANT.so")
dyn.load("IS_DOMINANT.so")
library(e1071)
library(monmlp)
library(MASS)
library(coda)
library(numDeriv)



MRM <- function(f, ndim, choix.loi, dir.monot, Method, N.appels, code, ordre.p, MAXIMIN="MLE", choix.classif="MonMLP", epsilon=0.1, alpha_start = (1/ndim)*(2.38)^2, transfMCMC=TRUE, Burn_in="GelmanRubin"){
  
  if(missing(MAXIMIN)){print("Missing MAXIMIN : default MAXIMIN MLE chosen")}
  if(missing(choix.classif)){print("Missing classifier choice : default MonMLP classifier chosen")}
  if(missing(epsilon) && (Method=="MRM_IS")){print("Missing epsilon choice (MRM_IS): default espilon=0.1 chosen")}
  if(missing(alpha_start)){print("Missing alpha_start : default alpha_start=(1/ndim)*(2.38)^2 chosen")}
  if(missing(transfMCMC)){print("Missing transfMCMC : default transfMCMC=TRUE chosen")}
  if(missing(Burn_in)){print("Missing Burn_in : default Burn_in=GelmanRubin chosen")}

  
  #-----------------------------------------------------------------------------------------------------------------#
  #                         Cr?ation de la liste contenant le nom des lois d'entr?es et			
  #                                    des param?ters de chacune des lois
  #-----------------------------------------------------------------------------------------------------------------#
  InputDistribution <- function(choix.loi){
    
    InputDist <- list()
    InputDist <- choix.loi
    
    for(i in 1:ndim){
      nparam <- length(choix.loi[[i]][[2]])
      for(j in 1:nparam){
        InputDist[[i]]$q <- paste("q",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$p <- paste("p",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$d <- paste("d",InputDist[[i]][[1]], sep = "");
        InputDist[[i]]$r <- paste("r",InputDist[[i]][[1]], sep = "");
      }
    }
    InputDist
  }
  
  InputDist <- InputDistribution(choix.loi)
  #print(InputDist)
  #-----------------------------------------------------------------------------------------------------------------#
  #	                           Transformation dans l'espace uniforme	  (Y_i suppos?s ind?pendants)
  #-----------------------------------------------------------------------------------------------------------------#
  G <- function(X){
    XU <- numeric()
    for(i in 1:ndim){
      if(dir.monot[i] == -1){X[i] <- 1 - X[i]}
      XU[i] <- do.call(InputDist[[i]]$q,c(list(X[i,drop = FALSE]), InputDist[[i]][[2]]))
    }
    return(f(XU))
  }
  #-----------------------------------------------------------------------------------------------------------------#
  #                Calcul de l'intersection entre la diagonale et l'etat limite de defaillance (non v?rifi? par Bastien)
  #-----------------------------------------------------------------------------------------------------------------#
  Intersect <- function(ndim, FUNC, N.dicho){
    
    #eps   <- 1
    a     <- 2
    k     <- 2
    res   <- list()
    u.new <- 0
    temp  <- 0
    u.dep <- list()
    out   <- list()
    comp  <- 2
    
    u.dep[[1]]  <- rep(1/2, ndim)
    temp        <- FUNC(u.dep[[1]])
    u.dep[[2]]  <- sign(temp)
    cp 	    <- 1 						# compteur d'appel ? G
    u.other     <- u.dep
    u.new       <- u.dep[[1]]
    LIST        <- list()
    LIST[[1]]   <- u.dep
    list.set    <- 0
    list.set[1] <- LIST[[1]][2]
    
    if(temp > 0){
      u.new <- u.dep[[1]] - 1/(a^k)
    }else{
      u.new <- u.dep[[1]] + 1/(a^k)
    }
    
    eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )
    
    if( ( u.other[[2]] != sign(temp)) & (eps > 1e-7) ){
      exist   <- exist + 1
      u.other <- list( u.dep[[1]], sign(temp) )
    }
    k <- k + 1
    
    while(cp < N.dicho){         #N.dicho d?finie en d?but de programme
      
      u.dep[[1]] <- u.new 
      temp       <- FUNC(u.dep[[1]])
      u.dep[[2]] <- sign(temp)
      cp 	     <- cp + 1
      
      if(temp > 0){
        u.new <- u.dep[[1]] - 1/(a^k)
      }else{
        u.new <- u.dep[[1]] + 1/(a^k)
      }
      
      eps <- ( u.new - u.dep[[1]] )%*%( u.new - u.dep[[1]] )
      k              <- k + 1
      LIST[[comp]]   <- u.dep
      list.set[comp] <- LIST[[comp]][[2]]
      comp           <- comp + 1
    }
    return(LIST)
    list.set <- as.numeric(list.set)
    
    # Si il n'y a qu'un seul point : 
    if( abs(sum(list.set)) == length(LIST)){
      res[[1]] <- LIST[[length(LIST)]][[1]]
      res[[2]] <- LIST[[length(LIST)]][[2]]
      out <- list(res, cp)
      return(out)
    }else{
      u.dep[[1]] <- LIST[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- LIST[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- LIST[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- LIST[[max(which(list.set == 1))]][[2]]
      res[[1]] <- rbind(u.dep[[1]], u.other[[1]])
      res[[2]] <- c(u.dep[[2]], u.other[[2]])
      out <- list(res, cp)
      return(out)
    }
    
  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #                         		 Fonction qui test les points de x qui sont 
  #				                  "plus petit" ou "plus grand"  que y.
  #					             Si set == 1 : Renvoi TRUE si x[i]<= y				   
  #					             Si set == 2 : Renvoi TRUE si x[i]=> y		
  #		   
  #-----------------------------------------------------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------------------------------------------------#
  #				                          Code en C
  #-----------------------------------------------------------------------------------------------------------------#
  # COMM_BASTIEN : Renvoie une liste de longueur celle de l'echantillon X, avec pour chaque point x_i (dans U) la valeur 1(TRUE) si x domine tous les 
  # points de Y,  0 sinon (si set == 2). Sinon (set==1) renvoie pour chaque x_i 1 si x_i est domin? par  tous les points de Y, 0 sinon.
  is.dominant.C <- function(x, y, ndim, set){
    L.x <- 0
    
    if(is.null(dim(x)) ){
      L.x <- 1;
    }else{
      L.x <- dim(x)[1]
    }
    
    #   set <- ifelse(set == 1, 2, -1)
    
    dens <- .C("IS_DOMINANT", 
               as.double(x), 
               as.double(y), 
               as.integer(ndim), 
               as.integer(set),
               as.integer(L.x), 
               result = integer(L.x));
    
    result <- dens$result
    
    rm(dens);rm(L.x);
    
    return(result);
  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #                         				Code en R
  #-----------------------------------------------------------------------------------------------------------------#
  is.dominant.R <- function(x, y, ndim, set){
    
    if((set != 1)&(set != 2)){
      print("ERROR : set must to be equal to 1 or 2.")
      break()
    } 
    
    dominant <- NULL;
    
    # Cas o? x est un vecteur
    if( is.null(dim(x)) ){
      if(set == 2){
        if ( sum(x >= y) == ndim ){
          return(TRUE)
        }else{
          return(FALSE)
        } 
      }else{
        if( sum(x <= y) == ndim ){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }
    
    y.1 	 <- NULL
    Y.2  	 <- NULL
    y.1 	 <- rep(y,dim(x)[1])
    y.2 	 <- matrix(y.1, ncol = ndim, byrow = TRUE)
    
    if(set == 2){
      dominant <- apply(x >= y.2, 1, sum) == ndim
    }else{
      dominant <- apply(x <= y.2, 1, sum) == ndim
    }
    
    return(dominant)
  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #		                Fonction qui retourne la valeur des volume
  #		                     Un- ou  Un+ (d?termin?s par les points S)
  #			                par la methode de Monte Carlo
  #-----------------------------------------------------------------------------------------------------------------#
  # COMM_BASTIEN : NON ca renvoie une valeur approch?e par MC (en limitant le nb de tests de domination) des volumes Vol(Un-) (set = 2)
  # ou 1-Vol(Un+) (set=1), qu'on suppose d?finis ? partir de l'ensemble de points S. Si Un- : On tire des unifs sur U puis on compte d'abord
  # celles domin?es par le point uu de S le + eloign? de l'origine. On les enl?ve de l'echantillon d'unifs de meme que celle qui dominent
  # uu (car inutiles pour la suite), puis on recommence avec le 2eme point le + eloign? de 0, et ainsi de suite. De m?me pour Un+.
  #-----------------------------------------------------------------------------------------------------------------#
  VOLVOL <- function(X.MC, S, set){ 
    
    if(set == 1){S <- 1 - S}
    
    if(is.null(dim(S))){
      if(set == 1){
        return(1 - prod(S))
      }
      if(set == 2){
        return(prod(S))
      }
    }
    
    DS <- dim(S)[1] 
    if(ndim == 2){
      S    <- S[order(S[,1]),]
      res  <- diag(outer(S[,1], c(0,S[1:(DS - 1), 1]), "-"))
      res1 <- res%*%S[,2]
      res1 <- ifelse(set == 1, 1 - res1, res1)
      return(res1) 
    }
    
    
    MC.VOL <- 0
    RES.VOL <- 0
    MM <- dim(X.MC)[1]
    
    for(i in 1:(DS-1)){
      u  <- apply(S, MARGIN = 1, prod)
      u.max <- which.max(u)
      uu    <- S[u.max, ]
      S     <- S[-u.max, ]
      ss.1 <- is.dominant.C(X.MC, uu, ndim, 1) 
      ss.2 <- is.dominant.C(X.MC, uu, ndim, 2) 
      RES.VOL <- RES.VOL + sum(ss.1)
      X.MC <- X.MC[which((ss.1 == 0)&(ss.2 == 0) ) , ]
    }
    
    uu <- S
    if(set == 1){
      tt <- is.dominant.C(X.MC, uu, ndim, 2)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- 1 - RES.VOL/MM
    }else{
      tt <- is.dominant.C(X.MC, uu, ndim, 1)
      RES.VOL <- RES.VOL + sum(tt)
      RES.VOL <- RES.VOL/MM
    }
    return(RES.VOL) 
    
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #				                     Choix du langage
  #-----------------------------------------------------------------------------------------------------------------
  if(code == "R"){
    is.dominant <- is.dominant.R
    #Volume.MC   <- Volume.MC.R
  }
  if(code == "C"){
    is.dominant <- is.dominant.C
    #Volume.MC   <- Volume.MC.C
  }
  
  #-----------------------------------------------------------------------------------------------------------------
  #				                     Extrait la fronti?re d'un ensemble
  #-----------------------------------------------------------------------------------------------------------------
  # COMM_BASTIEN : Utilis?e dans Monte-Carlo (MC.monotone)
  Frontier <- function(S, set){
    if(is.null(dim(S)) |( dim(S)[1] == 1)){
      return(S)
    }
    R <- NULL
    if(set == 1){
      S <- 1 - S
    }
    while(!is.null(dim(S))){
      aa  <- apply(S, MARGIN = 1, prod)
      temp <- S[which.max(aa), ]
      R <- rbind(R, temp)
      S <- S[-which.max(aa), ]
      ss <- is.dominant(S, temp, ndim, -1)
      if(!is.null(dim(S))){
        S <- S[which(ss == FALSE),]
      }else{
        S <- matrix(S, ncol= ndim)
        S <- S[which(ss == FALSE), ]
        R = rbind(R, S)
        if(set == 1){
          return(1 - R)
        }else{
          return(R)
        }
      }
    }
    
    if(set == 1){
      return(1 - R)
    }else{
      return(R)
    }
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  # 		         				methode de Monte Carlo (Monotone)
  #
  #-----------------------------------------------------------------------------------------------------------------
  Method.MC.monotone <- function(N.appels){ 
    
    if(ndim > 2){
      UU <- runif(ndim*10^(ordre.p + 2))
      UU <- matrix(UU, ncol=ndim, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }
    
    NN <- 0
    N.tot <- 0
    res <- NULL
    Z.safe <- NULL
    Z.fail <- NULL
    TADAM <- 0
    while(NN < N.appels){
      if(N.tot%%100 == 0){print(NN);flush.console();}
      U <- runif(ndim)
      N.tot <- N.tot + 1
      
      if( is.null(Z.safe)& is.null(Z.fail) ){
        t.u <- G(U)
        NN <- NN + 1
        TADAM[NN] <- N.tot
      }
      
      if(is.null(Z.safe)&( !is.null(Z.fail)) ){
        ttf <- is.dominant(Z.fail, U, ndim, set = 2)
        if( sum(ttf) == 0 ){
          t.u <- G(U)
          NN <- NN + 1
          TADAM[NN] <- N.tot          
        }else{
          t.u <- -1
        }
      }
      
      if(!is.null(Z.safe) & is.null(Z.fail) ){       
        tts <- is.dominant(Z.safe, U, ndim, set = 1)
        if(sum(tts) == 0){
          t.u <- G(U)
          NN <- NN + 1
          TADAM[NN] <- N.tot          
        }else{
          t.u <- 1
        }
      }
      
      if((!is.null(Z.safe)) &( !is.null(Z.fail)) ){      
        ttf <- is.dominant(Z.fail, U, ndim, set = 2)       
        tts <- is.dominant(Z.safe, U, ndim, set = 1)
        if( (sum(tts) == 0) & (sum(ttf) == 0) ){
          t.u <- G(U)
          NN <- NN + 1
          TADAM[NN] <- N.tot
        }
        if( (sum(tts) == 0)& (sum(ttf) != 0) ){
          t.u <- -1
        }
        if( (sum(tts) != 0)& (sum(ttf) == 0) ){
          t.u <- 1
        }
        
      }
      
      if(t.u <= 0){
        Z.fail <- rbind(Z.fail, U)
        res <- c(res, 1)
      }else{        
        Z.safe <- rbind(Z.safe, U)
        res <- c(res, 0)
      }
    }
    I <- 1:N.tot
    alpha <- 0.05
    cum.res <- cumsum(res)
    M.res <- cum.res/I
    Var   <- (M.res)*(1 - M.res)/I
    u.m   <- M.res - qnorm(1 - alpha/2)*sqrt(Var)
    u.M   <- M.res + qnorm(1 - alpha/2)*sqrt(Var)
    CV <- 100*sqrt(Var)/M.res
    
    #browser();
    if(is.null(Z.fail)){
      um <- 0
    }else{
      ZF <- Frontier(Z.fail, -1)
      um <- VOLVOL(UU, ZF, -1)
    }
    
    if(is.null(Z.safe)){
      uM <-1
    }else{ 
      ZS <- Frontier(Z.safe, set = 1)
      uM <- VOLVOL(UU, ZS, 1)
    }
    
    return(list(cbind(u.m, u.M, M.res, CV, Var)[TADAM, ], um, uM, N.tot))
  }
  
  Method.MC <- function(N.appels){
    
    X <- matrix(runif(N.appels*ndim), ncol = ndim)
    Y <- apply(X, MARGIN = 1, G)
    I <- 1:N.appels
    
    alpha <- 0.05
    
    res     <- 1*(Y <= 0)
    cum.res <- cumsum(res)
    
    M.res <- cum.res/I
    Var   <- (M.res)*(1 - M.res)/I
    u.m   <- M.res - qnorm(1 - alpha/2)*sqrt(Var)
    u.M   <- M.res + qnorm(1 - alpha/2)*sqrt(Var)
    
    return(cbind(u.m,u.M,M.res))
  }
  
  ################################################################
  #             SUBSET MONOTONE
  ################################################################
  
  M.Subset <- function(N.appels, H, N.dicho){
    
    
    if(ndim > 2){
      UU <- runif(ndim*10^(ordre.p + 2))
      UU <- matrix(UU, ncol=ndim, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }
    
    V <- list()
    V <- Intersect(ndim, H, N.dicho)  
    
    u.dep   <- list()
    u.other <- list()
    
    
    list.set <- 0
    for(i in 1:length(V)){
      list.set[i] <- V[[i]][[2]]
    }
    
    u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
    u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
    
    u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
    u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
    
    ZF <- Z.fail <- t(as.matrix(u.dep[[1]]))
    ZS <- Z.safe <- t(as.matrix(u.other[[1]]))
    
    cp <- length(V)   
    um <-prod(V[[cp]][[1]])
    uM <- 1 - prod(1 - V[[cp]][[1]])
    p.hat <- 0  
    SIGN <- 0
    j <- 1
    Um <- 0
    UM <- 1
    
    W <- 1:(2^(ndim) - 1)
    
    while(cp < N.appels){
      U <- Sim.non.dominated.space(1, ZS, ZF, W)
      tts1 <- is.dominant(Z.safe, U, ndim, 1)
      ttf1 <- is.dominant(Z.fail, U, ndim, 2)
      TEST <-  (sum(tts1) == 0 ) & ( sum(ttf1) == 0) 
      
      if( sum(tts1) > 0 ){
        SIGN[j] <- 0
      }
      if( sum(ttf1) > 0 ){
        SIGN[j] <- 1
      }
      
      if(!TEST){
        Um[j] <- ifelse(j==1,um , Um[j-1])
        UM[j] <- ifelse(j==1,uM , UM[j-1])
      }
      
      if(TEST){
        SIGN[j] <- (1 - sign(H(U)))/2
        
        if(SIGN[j] == 0){
          
          ss     <- is.dominant(Z.safe, U, ndim, 2)
          Z.safe <- Z.safe[which(ss == FALSE), ]
          Z.safe <- rbind(U, Z.safe)
          
          UM[j] <- VOLVOL(UU, Z.safe , 1)
          Um[j] <- ifelse(j==1,um , Um[j-1])
        }
        if(SIGN[j] == 1){
          
          ff     <- is.dominant(Z.fail, U, ndim, 1)
          Z.fail <- Z.fail[which(ff == FALSE), ]
          Z.fail <- rbind(U, Z.fail)
          
          Um[j] <- VOLVOL(UU, Z.fail, 2)
          UM[j] <- ifelse(j==1,uM , UM[j-1])
        }
        
        cp <- cp + 1
      }
      
      p.hat[j] <- um + (uM - um)*mean(SIGN)
      
      j <- j + 1
      print(cp);flush.console();
    }
    alpha <- 0.025
    Var   <- ((uM - um)^2) *mean(SIGN)*(1 - mean(SIGN))/length(SIGN)
    ICinf   <- p.hat - qnorm(1 - alpha/2)*sqrt(Var)
    ICsup   <- p.hat + qnorm(1 - alpha/2)*sqrt(Var)
    return(cbind(Um, UM, p.hat,ICinf,ICsup ,SIGN))
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #	               	              Construction du maximum de vraissemblance
  #-----------------------------------------------------------------------------------------------------------------x
  
  log.likehood <- function(p , p.k, signature){
    gamma <- (p - p.k[,1])/(p.k[,2] - p.k[,1])
    u     <- (gamma^signature)*((1 - gamma)^(1 - signature))
    return(prod(u))
  }
  
  
  #max de vraisemblance avec ancienne simulation
  log.likehood.with.old <- function(p , p.k, signature){
    pp <- 0
    dd <- dim(signature)[1]
    for(i in 1:dd){
      gamma <- (p - p.k[i,1])/(p.k[1:i,2] - p.k[1:i,1])
      pp <- pp + sum(    signature[1:i,i]*log( gamma) + (1 - signature[1:i,i])*log(1 - gamma) ) 
    }
    return(pp)
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #		                            Convertit un entier en binaire
  #-----------------------------------------------------------------------------------------------------------------
  #COMM_BASTIEN: Reste dans la division par 2 (parit?) pour chaque composante de x. Utilis?e dans SIM.
  as.binary <- function (x) { 
    base <- 2;
    r <- numeric(ndim)
    for (i in ndim:1){ 
      r[i] <- x%%base 
      x <- x%/%base
    } 
    return(r) 
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #		       On tire 1 point uniform?ment dans U\U^+(x)  (cf.Klein MOUOUSSAMY)
  #-----------------------------------------------------------------------------------------------------------------
  #COMM_BASTIEN :Fonction SIM. Utilis?e pour simuler dans l'espace non domin?. MOdifi?e pour g?nerer un ?chantillon. Marche bien.
  
  SIM <- function(CP, x, W){
    AA <- NULL
    #      W <- 1:(2^(ndim) - 1) (on partitionne l'espace autour de x en W zones)
    B <- 0
    # On s?pare l'espace autour de x en 2^d - 1 zones, et l'on calcul le volume de chaque zone.
    B <- apply( matrix(W, ncol = 1), 
                MARGIN = 1, 
                function(v){
                  Z <- as.binary(v)
                  v <- 0
                  u <- 0
                  for(j in 1:ndim){
                    u[j] <- ifelse(Z[j] == 0, 1 - x[j], x[j])
                  }     
                  return(prod(u))  
                }
    )
    
    B <- cumsum(B)
    
    for (i in 1:CP){
      U <- runif(1, 0, max(B))
      #On choisit al?atoirement dans quelle "zone" autour de x on va tirer des points
      pos <- ifelse(U < B[1], 1, which.max(B[B <= U]) + 1)
      Z <- as.binary(W[pos])
      A <- 0
      #On tire uniform?ment dans la zone choisit al?atoirement
      for(i in 1:ndim){
        A[i] = ifelse(Z[i] == 0, runif(1, x[i], 1), runif(1, 0, x[i]))  
      }
      AA= rbind(A, AA)
      
    }
    return(AA)
  }
  
  #-----------------------------------------------------------------------------------------------------------------
  #                   Simulation de CP points uniform?ment dans l'espace non domin?
  #-----------------------------------------------------------------------------------------------------------------
  #Comm_BASTIEN : On utilise la methode du rejet. Sauf que plut?t que de simuler des unif sur tout U puis de garder que celles dans U_n (non domin?)
  #                on fait mieux en simulant avec SIM sur U\U^+(x), o? x est le point de Z.safe le plus "proche" de l'origine.
  Sim.non.dominated.space <- function(CP, Z.safe, Z.fail, W){
    #print(c("Z.safe size :", dim(Z.safe)[1], "Z.fail size :", dim(Z.fail)[1] )) #BASTIEN !!!
    CP1 <- 0;
    Y   <- NULL
    Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
    Y.temp1 <- Z.safe[which.max(Y.temp),]
    while(CP1 < CP){
      Y.temp2 <- SIM(1, Y.temp1, W)[1,]
      tts1 <- is.dominant(Z.safe, Y.temp2, ndim, 1)
      ttf1 <- is.dominant(Z.fail, Y.temp2, ndim, 2)
      if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) ){
        Y <- rbind(Y,Y.temp2)
        CP1 <- CP1 + 1
      } 
    }
    return(Y)
  }
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  #
  # 		     				   PARTIE 1 : ESTIMATEUR MRM MLE
  #
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  Choix.method.1 <- function(N.appels, H, N.dicho){
   
    start_time <- Sys.time()
    # V contient le premier tri de points avec la methode de la dichotomie.
    # Z.fail : ensemble des points (definissant le) du domaine de d?faillance
    # Z.safe : ensemble des points (definissant le) du domaine de s?ret?
    # cp : compteur d'appel ? H (=G)
    # um : vecteur contenant les bornes exactes inf?rieures
    # uM : vecteur contenant les bornes exactes sup?rieures
    # MLE :Estimateur du maximum de vraissemblance
    # I.Fisher : information de Fisher associ? au MLE
    # CV : vecteur contenant les coefficients de variations
    # SIGN : Signature des points test?s
    
    Um <- 0
    UM <- 1
    eps <- 1e-7
    alpha <- 0.025   # CHOIX DE ALPHA ?
    SIGN  <- 0
    ICinf <- 0
    ICsup <- 0
    VAR   <- 0
    CV.MLE    <- 0
    CV.max    <- 0
    MLE <- 0
    p.bar <- 0
    ZS <- list()
    ZF <- list()
    um <- 0
    uM   <- 1
    # Cr?ation du tirage uniforme si ndim > 2
    if(ndim > 2){
      UU <- runif(ndim*10^(ordre.p + 2))
      UU <- matrix(UU, ncol=ndim, byrow = TRUE)
      D.UU <- dim(UU)[1]
    }
    
    if(ndim == 2){UU <- 0; D.UU <- 0}
    
    #DICHOTOMIE
    if(TRUE){
      V <- list()
      V <- Intersect(ndim, H, N.dicho)    
      list.set <- 0
      for(i in 1:length(V)){
        list.set[i] <- V[[i]][[2]] }
      #initialisation de Z.fail et Z.safe a partir des points calcul?s (V) pdt la dichotomie initiale.
      u.dep   <- list()
      u.other <- list()
      u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
      Z.fail <- t(as.matrix(u.dep[[1]]))
      Z.safe <- t(as.matrix(u.other[[1]]))
      cp     <- length(V)  
      #initialisation des bornes a partir du 1er elt de U- (resp U+) calcul? par dichotomie
      um <- prod(V[[cp]][[1]]) 
      uM   <- 1 - prod(1 - V[[cp]][[1]])
    }
    
    j = 1
    
    # si on choisit de calculer le MLE et les IC avec likehood.with.old
    SIGN.2.old <- 0
    v = 0
    VAR.2 <- 0   
    ICinf.2 <- 0
    ICsup.2 <- 0
    
    W <- 1:(2^(ndim) - 1)
    #library("stats")
    
    SIMU <- NULL
    
    time100 <-0 ;  Time100 <-0
    while(cp < N.appels){
      
      if (j==110) {Time100 <- Sys.time() ; time100 <-as.double(difftime(Time100, start_time, tz,units = "secs"))}
      
      # CP : Nombre de points que l'on veut avoir dans l'espace non domin?
      # Y : ensembles des points dans l'espace non domin?
      
      cat(paste("  compteur =",cp, "--"));flush.console()
      # if (cp%%10 == 0){cat(paste("  compteur =",cp, "--"));flush.console() }
      
      CCPP <- 1
      u <- Sim.non.dominated.space(CCPP, Z.safe, Z.fail, W)
      
      # On trie les points de Y : 
      SIMU <- rbind(SIMU, u)
      Z.safe.old <- 0
      Z.fail.old <- 0
      # On regarde dans quel ensemble appartient u
      # On rajoute u ? l'ensemble auquel il appartient. On garde sa signature
      # et on calcule le volume du nouvel ensemble construit
      H.u = H(u)
      if(H.u > 0){
        Z.safe.old <- rbind(u, Z.safe)
        ss     <- is.dominant(Z.safe, u, ndim, 2)
        # on enleve les points definissant Z.safe qui sont domin?s par le nouveau point u 
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(u, Z.safe)
        vol <- VOLVOL(UU,Z.safe, 1)
        #vol <- Volume.MC.R(Z.safe, 1)  #on utilise le m?me nb de tirage d'unifs : ndim*10^(ordre.p + 2)
        
        if(j == 1){
          Um[j] <- um
          if(vol >= uM){
            UM[j] <- uM
          }else{
            UM[j] <- vol
          }
        }else{
          Um[j] <- Um[j - 1]
          if(vol >= UM[j-1]){
            UM[j] <- UM[j-1]
          }else{
            UM[j] <- vol
          }
        }
        
        #CC <- ifelse(cp == N.appels, 1, 0)
        SIGN[j] <- 0
        
      }else{
        Z.fail.old <- rbind(u, Z.fail)
        
        ff     <- is.dominant(Z.fail, u, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(u, Z.fail)
        vol = VOLVOL(UU, Z.fail, 2)
        #vol <- Volume.MC.R(Z.fail, 2)
        
        if(j == 1){
          UM[j] <- uM
          if(vol <= um){
            Um[j] <- um
          }else{
            Um[j] <- vol
          }
        }else{
          UM[j] <- UM[j - 1]
          if(vol <= Um[j-1]){
            Um[j] = Um[j-1]
          }else{
            Um[j] <- vol
          }
        }
        
        
        CC <- ifelse(cp == N.appels, 1, 0)
        SIGN[j] <- 1
      }
      cp <- cp + 1
      
      ZS[[j]] <- Z.safe
      ZF[[j]] <- Z.fail
      
      if(MAXIMIN == "MLE"){
        # MLE.test <- optimize(f = log.likehood,
        #                      interval = c(Um[j],UM[j]),
        #                      maximum = TRUE,
        #                      tol = 10^(-3*ordre.p),
        #                      signature = SIGN,
        #                      p.k = cbind(Um, UM)
        # )   
        # MLE[j] <- as.numeric(MLE.test[1])
        # 
        # # VAR Corresponds a \hat{J}_n(\hat{p_n}) cf. prop 3.2 article Bousquet! (o?  \hat{p_n} est le MLE)
        # VAR <- sum( 1/((MLE - Um)*(UM- MLE)))
        # bn  <- 1/VAR
        # 
        # # cf. Proposition 3.4. Approximation de (\hat{J}_n(\hat{p_n}))**(5/2) \ [ |\hat{J}_n(\hat{p_n})|]
        # an  <- eps*VAR^(5/2)/abs( sum( 1/((MLE + eps - Um)*(UM - MLE - eps))) - sum( 1/((MLE - Um)*(UM- MLE)))  )
        # 
        # # Comment on calcule cet I.C ??
        # ICinf[j]  <- MLE[j] - qnorm(1 - alpha)/sqrt(VAR + alpha/an)
        # ICsup[j]  <- MLE[j] + qnorm(1 - alpha)/sqrt(VAR - alpha/an) 
        # 
        # #Pourquoi on calcule le  CV ainsi ??
        # CV.MLE[j] <- 100/(sqrt(VAR)*MLE[j])   
        
        # if(FALSE){
        #   SIGN.2 <- matrix(0, ncol = j, nrow = j)
        #   diag(SIGN.2) <- SIGN
        #   
        #   if(j > 1){
        #     SIGN.2[1:(j-1),1:(j-1)] <- SIGN.2.old 
        #     for(jj in 1:(j-1)){
        #       if( SIGN.2[jj, j-1] == 1){
        #         if( SIGN[j] == 1){
        #           if( sum( SIMU[jj,  ] <= SIMU[j, ] ) != ndim ){
        #             SIGN.2[jj, j] <- 1
        #           }
        #         }else{
        #           SIGN.2[jj, j] <- 1
        #         }
        #       }
        #     }#Fin de la boucle for
        #   }
        #   
        #   #browser();
        #   SIGN.2.old <- SIGN.2
        #   
        #   #browser();
        #   #MLE.test.2 <- optimize(f = log.likehood.with.old, interval = c(Um[j],UM[j]),maximum = TRUE,signature = SIGN.2,  p.k = cbind(Um, UM) )
        #   MLE.test.2 <- optimize(f = log.likehood.with.old, interval = c(Um[j],UM[j]),maximum = TRUE,tol = 10^(-3*ordre.p) ,signature = SIGN.2,  p.k = cbind(Um, UM) )      
        #   
        #   
        #   p.bar[j] <- as.numeric(MLE.test.2[1])
        #   
        #   VAR.2 <- VAR.2 + sum(1/( (UM - Um - p.bar[j] + UM[j])*(p.bar[j] - Um[j] )  )  )
        #   ICinf.2[j]  <- p.bar[j] - qnorm(1 - alpha)/sqrt(VAR.2)
        #   ICsup.2[j]  <- p.bar[j] + qnorm(1 - alpha)/sqrt(VAR.2)
        # } # si on choisit de calculer aussi le MLE et les IC avec likehood.with.old (je n'ai pas encore lu cette partie)
        
        
      }
      
      
      j <- j + 1 
      
    } #fin du while cp < N.appels
    # browser();
    if(MAXIMIN  == "MLE"){
      print("On est dans le MLE")
      #RR1 <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE)
      #RR1 <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE, p.bar.1, p.bar.2, p.bar.3 ) # si on choisit de calculer le MLE aveclikehood.with.old
      #RR1 <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE, p.bar, ICinf.2, ICsup.2)
      #RR1 <- cbind(Um, UM,  MLE)
      RR1 <- cbind(Um, UM)
      #browser()
      RR <- list(tab=RR1, ZF=ZF, ZS=ZS, time100=time100)
      #return(RR) 
      return(RR) 
      
    }
    return(RR)
  }
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  # 					               FIN	METHODE MRM MLE
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #                  Entra?nement du classifieur
  #
  #-----------------------------------------------------------------------------------------------------------------
  classifieur.train<- function(choix.classif, ZF, ZS){
    #IL FAUDRAIT demander en entree des specifications optionnelles sur les classifieurs
    #VERIFIE POUR SVM ET MonMLP
    len_ZF=dim(ZF)[1]
    len_ZS=dim(ZS)[1]
    ZZ=rbind(ZS, ZF)
    y_val <- NULL
    for (i  in 1:(len_ZS+len_ZF)){
      if(i<=len_ZS){y_val<-c(y_val, 1)}
      else{y_val<-c(y_val, 0)}
    }
    y_val=as.matrix(y_val, nrow=(len_ZS+len_ZF))
    
    
    if (choix.classif=="SVM"){
      y_val <- factor(y_val, labels = c(0, 1))
      #mod_lin <- svm(x=UU, y=y_val, scale = FALSE, type = "eps-regression", kernel = "linear" )
      mod_poly <- svm(x=ZZ, y=y_val, scale = TRUE, type = "C-classification", coef0 = 1,
                      kernel = "polynomial", degree = 3, cost = 100, epsilon = 0.1 )
      return(mod_poly)
    }
    
    if (choix.classif=="SVM-reg"){
      
      mod_poly <- svm(x=ZZ, y=y_val, scale = TRUE, type = "eps-regression", coef0 = 1,
                      kernel = "polynomial", degree = 3, cost = 1, epsilon = 0.01 )
      return(mod_poly)
    }
    
    if (choix.classif=="MonMLP"){
      
      w.mon <- monmlp.fit(x = ZZ, y = y_val, hidden1 = 4, hidden2 = 4, monotone = 1,
                          bag = FALSE, iter.max = 500, iter.stopped = 10, silent = TRUE)
      return(w.mon)
    }
    
    if (choix.classif=="MonSVM"){
      
      
    }
    
    
  }
  #-----------------------------------------------------------------------------------------------------------------
  #
  #                  Predicition du classifieur
  #
  #-----------------------------------------------------------------------------------------------------------------
  #X est un ensemble de points qu'on cherche ? evaluer avec le classifieur
  #Renvoie 1 si on est dans hat{U}_k^+ (0 sinon)
  #weights=classifieur.train(choix.classif, ZF, ZS)
  #VERIFIE POUR SVM ET MonMLP
  classifieur.pred<- function(X, choix.classif, weights){
    
    if (choix.classif=="SVM"){
      ypoly.test <- as.numeric(predict(weights,X))-1
      return(ypoly.test)
    }
    
    if (choix.classif=="SVM-reg"){
      ypoly.test <- ifelse( (predict(weights,X)>0.5), 1, 0 )
      return(ypoly.test)
      
      
    }
    
    if (choix.classif=="MonMLP"){
      ypoly.test <- ifelse( (monmlp.predict(X, weights=weights)>0.5), 1, 0 )    
      return(ypoly.test)
    }
    
    if (choix.classif=="MonSVM"){
      
      
      
    }
  } 
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #                   Simulation de CP points uniform?ment dans l'espace d?fini par l'intersection de
  #                   la zone lab?lis?e >0 (dans hat\U^+) par le classifieur et l'espace non domin?.
  #                  
  #-----------------------------------------------------------------------------------------------------------------
  #Comm_BASTIEN : On utilise la methode du rejet. On peut  utiliser la m?me approche acc?l?r?e que dans le cas uniforme.
  Sim.class <- function(CP, choix.classif, weights, Z.fail, Z.safe, W){
    CP1 <- 0;
    Y   <- NULL
    Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
    Y.temp1 <- Z.safe[which.max(Y.temp),]
    while(CP1 < CP){
      Y.temp2 <- SIM(1, Y.temp1, W)[1,]
      Y.temp2=matrix(Y.temp2, ncol=ndim)
      tts <- is.dominant(Z.safe, Y.temp2, ndim, 1)
      ttf <- is.dominant(Z.fail, Y.temp2, ndim, 2)
      if( (classifieur.pred(Y.temp2, choix.classif, weights)==1) & ( sum(tts) == 0) & ( sum(ttf) == 0)){
        Y <- rbind(Y,Y.temp2)
        CP1 <- CP1 + 1
      }
    }
    return(Y)
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonction qui retourne la valeur approch?e par MC du volume d?fini par le classifieur  hat{U}_n^+ (set = 1), 
  #         ou hat{U}_n^- (set = 2).
  #
  #-----------------------------------------------------------------------------------------------------------------#
  #COMM_BASTIEN : On utilise la methode de Monte Carlo avec rejet pour estimer le volume de hat{U}_n^+ inter Un^+.
  VOLCLASS <- function(X.MC, choix.classif, weights, Z.safe, Z.fail){
    CP1 <- 0 
    for (i in 1:dim(X.MC)[1]) {
      Y.temp2 <- X.MC[i,]
      Y.temp2 <- matrix(Y.temp2, ncol=ndim)
      tts1 <- is.dominant(Z.safe, Y.temp2, ndim, 1)
      ttf1 <- is.dominant(Z.fail, Y.temp2, ndim, 2)
      if( (classifieur.pred(Y.temp2, choix.classif, weights)==1) & ( sum(tts1) == 0) & ( sum(ttf1) == 0)){
        CP1 <- CP1 + 1
      }
    }
    volume <- CP1/dim(X.MC)[1]
    return(volume)
  }
  
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  #
  # 		     				   PARTIE 2 : ESTIMATEUR MRM IS
  #
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  Choix.method.2 <- function(N.appels, H, choix.classif, N.dicho){
    
    # V contient le premier tri de points avec la methode de la dichotomie.
    # Z.fail : ensemble des points (definissant le) du domaine de d?faillance
    # Z.safe : ensemble des points (definissant le) du domaine de s?ret?
    # cp : compteur d'appel ? H (=G)
    # um : vecteur contenant les bornes exactes inf?rieures
    # uM : vecteur contenant les bornes exactes sup?rieures
    # MLE :Estimateur du maximum de vraissemblance
    # I.Fisher : information de Fisher associ? au MLE
    # CV : vecteur contenant les coefficients de variations
    # SIGN : Signature des points test?s
    # choix.classif : choix du classifieur
    
    Um <- 0
    UM <- 1
    eps <- 1e-7
    alpha <- 0.025   # CHOIX DE ALPHA ?
    SIGN  <- 0
    Classif  <- list()  #liste de listes contenant les donn?es du classifieur ? chaque ?tape
    val_f <- 0     #valeur de f_k-1(X_k:=u)
    ICinf <- 0
    ICsup <- 0
    VAR   <- 0
    val_f <- 0
    p.bar <- 0  # ICI p.bar correspond au vecteur des p.bar tels que definis dans Bousquet 2021.
    um <- 0
    uM   <- 1
    old_weights <- 0
    classif_init <- 0
    Z.safe<- NULL #Fronti?re de U_k+
    Z.fail<- NULL #Fronti?re de U_k-
    List.Z.safe <- list() # Liste des ensembles Z.safe selon l'it?ration
    List.Z.fail <- list() # Liste des ensembles Z.fail selon l'it?ration
    ZS <- NULL  #Contient  TOUT les points ?valu?s jusqu'a alors qui conduisent a la d?faillance
    ZF <- NULL  # Symm?trique
    
    
    #Cr?ation du tirage uniforme si ndim > 2 OU ndim=2 ! (pour estimer les bornes ou le volume lab?lis? par le classifieur si VOLCASS)
    UU <- runif(ndim*10^(ordre.p + 2))  
    # writeBin(UU, zz, size = 8)
    UU <- matrix(UU, ncol=ndim, byrow = TRUE)
    
    #if(ndim == 2){UU <- 0; D.UU <- 0}
    
    #DICHOTOMIE ET INITIALISATION DES BORNES/DU CLASSIFIEUR
    if(TRUE){
      V <- list()
      V <- Intersect(ndim, H, N.dicho)    
      list.set <- 0
      for(i in 1:length(V)){
        list.set[i] <- V[[i]][[2]] }
      #initialisation de Z.fail et Z.safe a partir des points calcul?s (V) pdt la dichotomie initiale.
      u.dep   <- list()
      u.other <- list()
      u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
      Z.fail <- t(as.matrix(u.dep[[1]]))
      Z.safe <- t(as.matrix(u.other[[1]]))
      ZS <- matrix(data = u.other[[1]] , ncol=ndim, byrow=TRUE)  
      ZF <- matrix(data = u.dep[[1]] , ncol=ndim, byrow=TRUE)  
      cp     <- length(V)  
      #initialisation des bornes a partir du 1er elt de U- (resp U+) caclul? par dichotomie
      um <- prod(V[[cp]][[1]]) 
      uM   <- 1 - prod(1 - V[[cp]][[1]])
      classif_init <- classifieur.train(choix.classif, ZF, ZS)
      #cat(paste("  classif_init =",classif_init));flush.console()
      #cat(paste("  type of classif_init =",typeof(classif_init)));flush.console()
    }
    
    j = 1
    
    # si on choisit de calculer le MLE et les IC avec likehood.with.old
    SIGN.2.old <- 0
    v = 0
    VAR.2 <- 0   
    ICinf.2 <- 0
    ICsup.2 <- 0
    
    W <- 1:(2^(ndim) - 1)
    #library("stats")
  
    while(cp < N.appels){
      # cp : Nombre de points que l'on veut avoir dans l'espace non domin?
      cat(paste("  j =",j));flush.console()
      cat(paste("  compteur =",cp, "--"));flush.console()
      #if (cp%%10 == 0){ cat(paste("  compteur =",cp, "--"));flush.console() }
      
      ##u <- Sim.non.dominated.space(1, Z.safe, Z.fail, W) (MRM classique)
      
      #Simulation de X_k:=u avec proba eps_{k-1} dans l'espace \hat{U}_{k-1}^+
      #   Pour simuler avec proba 1-eps_{k-1} dans U_{k-1}\hat{U}_{k-1}^+ on simule dans l'espace
      #   non domin? avec Sim.non.dominated.space puis on utilise la methode du rejet. PROBLEME si l'espace non
      #   domin? tout entier est inclus dans \hat{U}_{k-1}^+ (impossible avec un classif monotone)
      ru <- runif(n = 1)
      if(j>1){
        old_weights=Classif[[j-1]]  
      }else{
        old_weights=classif_init }
      #cat(paste("  old_weights =",old_weights));flush.console()
      if (ru<1-epsilon){
        cat(" simu.unif. U^+_j-1 ");flush.console()
        u <- Sim.class(CP=1 ,choix.classif, weights = old_weights, Z.fail, Z.safe, W)
        u <- matrix(u, ncol=ndim)
      }else{
        cat(" simu. unif ds U_j-1 prive de U^+_j-1 ")
        cpt <- 0
        cpt2 <- 0
        while(cpt<1  & cpt2<1000 ){
          u <- Sim.non.dominated.space(1, Z.safe, Z.fail, W)
          u <- matrix(u, ncol=ndim)
          #cat(paste("u", u))
          #cat(paste("dim_u", dim(u)))
          #cat(paste("  classifieur.pred_u =",classifieur.pred(u, choix.classif, weights=old_weights) )) 
          #cat(paste("  cpt2", cpt2 ));flush.console()
          #classifieur.pred(u, choix.classif, weights=ifelse(j==1,classif_init,classif[j-1]))==0
          if( classifieur.pred(u, choix.classif, weights=old_weights)==0 ){
            cpt <- 1
          }else{cpt2 <-cpt2+1 }
        }
        if (cpt2==1000){cat("\n TROP D'ITERATIONS (>1000) POUR SIMULER DANS U_k prive de Uchapeau_k+ ");flush.console()}
      }
      #cat("A ");flush.console()
      
      #u <- Sim.non.dominated.space(1, Z.safe, Z.fail, W)
      # On trie les points de l'espace non domin?: 
      Z.safe.old <- 0
      Z.fail.old <- 0
      # On regarde dans quel ensemble appartient u
      # On rajoute u ? l'ensemble auquel il appartient. On garde sa signature
      # et on calcule le volume du nouvel ensemble construit
      H.u = H(u)
      if(H.u > 0){
        #Z.safe.old <- rbind(u, Z.safe)
        ss     <- is.dominant(Z.safe, u, ndim, 2)
        # on enleve les points definissant Z.safe qui sont domin?s par le nouveau point u 
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(u, Z.safe)
        vol <- VOLVOL(UU,Z.safe, 1)
        # On mets a jour les bornes. Attention, ici la taille de l'echantillon UU n'a d'impact que sur la pr?cision du caclul
        # des bornes (calcul?es par VOLVOL) et non sur le caclul des estimateurs p.bar et hat{p}, 
        # CONTRAIREMENT au cas MRM MLE. (a part si on utilise VOLCLASS)
        if(j == 1){
          Um[j] <- um
          if(vol >= uM){
            UM[j] <- uM
          }else{
            UM[j] <- vol
          }
        }else{
          Um[j] <- Um[j - 1]
          if(vol >= UM[j-1]){
            UM[j] <- UM[j-1]
          }else{
            UM[j] <- vol
          }
        }
        
        
        SIGN[j] <- 0
        ZS<- rbind(ZS, u)
        
        
      }else{
        #Z.fail.old <- rbind(u, Z.fail)
        ff     <- is.dominant(Z.fail, u, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(u, Z.fail)
        vol = VOLVOL(UU, Z.fail, 2)
        
        if(j == 1){
          UM[j] <- uM
          if(vol <= um){
            Um[j] <- um
          }else{
            Um[j] <- vol
          }
        }else{
          UM[j] <- UM[j - 1]
          if(vol <= Um[j-1]){
            Um[j] = Um[j-1]
          }else{
            Um[j] <- vol
          }
        }
        
        
        SIGN[j] <- 1
        ZF<- rbind(ZF, u)
      }
      
      List.Z.safe[[j]] <- Z.safe
      List.Z.fail[[j]] <- Z.fail
      
      
      #cat("B ");flush.console()
      #la nouvelle real. X_k dont on a calcul? la signature permet de cacluler un nouveau classifieur :
      if(j==1){
        Classif <- list(classifieur.train(choix.classif, ZF, ZS))
      }else{
        Classif <-append(Classif, list(classifieur.train(choix.classif, ZF, ZS)) ) } #classif[j]
      
      #cat("C ");flush.console()
      #CALCUL DE p.bar. Volumes estim?s avec des uniformes sur tout l'espace U! Utilise VOLCLASS.
      if(j==1){
        volclass <- VOLCLASS(UU, choix.classif, classif_init, Z.safe = Z.safe, Z.fail = Z.fail )
      }else{
        volclass <- VOLCLASS(UU, choix.classif, Classif[[j-1]],  Z.safe = Z.safe, Z.fail = Z.fail ) # si on utilise des unifs sur U pour estimer le volume
      }
      #cat(paste("  volclass =",round(volclass, digits = 3)));flush.console()
      
      #CALCUL DE p.bar alternatif.
      #BCP PLUS PRECIS QU'AVEC VOLCLASS MAIS PRENDS BCP PLUS DE TEMPS CAR ON DOIT (?) (PERIODIQUEMENT) RECALCULER UN
      #  ECHANTILLON UNIFORME SUR UNE REGION ADAPTEE.
      # CCP <- 10^(ordre.p + 2)
      # if(j%%100 == 1){  #Avec quelle fr?quence on recalcule le point de U_n+ le plus proche de l'origine ?
      #   Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
      #   Y.temp1 <- Z.safe[which.max(Y.temp),]
      #   YY <- SIM(CCP, Y.temp1, W) #MULTISIM
      # }
      # weights.temp <- 0
      # if(j>1){ weights.temp = classif[[j-1]]}else{weights.temp = classif_init }
      # volclass2 <- 0
      # cp1 <- 1
      # while(cp1 < CCP+1){
      #   Y.temp3 <- matrix(YY[cp1, ], ncol=ndim)
      #   tts <- is.dominant(Z.safe, Y.temp3, ndim, 1)
      #   ttf <- is.dominant(Z.fail, Y.temp3, ndim, 2)
      #   if( (classifieur.pred(Y.temp3, choix.classif, weights.temp)==1) & ( sum(tts) == 0) & ( sum(ttf) == 0)){
      #     volclass2 <- volclass2 + 1
      #   }
      #   cp1 <- cp1+1
      # }
      # cat(paste("  prod(1-Y.temp1) =",round(prod(1-Y.temp1), digits = 3)));flush.console()
      # volclass2 <- volclass2*(1/CCP)*(1-prod(1-Y.temp1))
      # cat(paste("  volclass2 =",round(volclass2, digits = 3)));flush.console()
      
      
      
      if(j==1){
        val_f[j] <- epsilon*(1/((uM-um)-volclass))*ifelse(classifieur.pred(u, choix.classif , classif_init)==0, 1, 0) +
          (1-epsilon)*(1/(volclass))*ifelse(classifieur.pred(u, choix.classif , classif_init)==1, 1, 0)
      }else{
        val_f[j] <- epsilon*(1/((UM[j-1]-Um[j-1])-volclass))*ifelse(classifieur.pred(u, choix.classif , Classif[[j-1]])==0, 1, 0) +
          (1-epsilon)*(1/(volclass))*ifelse(classifieur.pred(u, choix.classif , Classif[[j-1]])==1, 1, 0)
      }
      cat(paste("  volnondom =",round(UM[j-1]-Um[j-1], digits = 3)));flush.console()
      cat(paste("  volclassif =",round(volclass, digits = 3)));flush.console()
      cat(paste("  val_f[j] =",round(val_f[j], digits = 3)));flush.console()
      
      
      if(j==1){
        p.bar[j] <- uM - (1/(val_f[j]))*ifelse(SIGN[j]==0, 1, 0)
      }else{
        p.bar[j] <- UM[j-1] - (1/(val_f[j]))*ifelse(SIGN[j]==0, 1, 0) #attention SIGN==1 si on est dans U-
      }
      #cat(" D \n");flush.console()
      cat(" \n");flush.console()
      
      
      cp<- cp + 1
      j <- j + 1 
      
    } #fin du while cp < N.appels
    
    
    
    print("On est dans la methode MRM IS ")
    print(c("Classifieur choisi :", choix.classif))
    #      RR <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE)
    #      RR <- list(Um, UM, MLE, ICinf, ICsup, CV.MLE, Z.fail, Z.safe)
    #      RR <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE, p.bar.1, p.bar.2, p.bar.3 ) # si on choisit de calculer le MLE aveclikehood.with.old
    RR1 <- cbind(Um, UM,  p.bar)
    #browser()
    RR <- list(tab=RR1, Classif=Classif, List.Z.fail=List.Z.fail, List.Z.safe=List.Z.safe, ZF=ZF, ZS=ZS )
    return(RR) 
    #return(list(ZS,ZF))
  }
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  # 					               FIN	METHODE MRM IS
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Algo de Metropolis-Hastings (avec marche al?atoire)
  #
  #-----------------------------------------------------------------------------------------------------------------#
  #COMM_BASTIEN :Sigma doit etre du type matrix (et positive def symm).
  RW_Metropolis<- function(startvalue, nb_iter, Z.safe, Z.fail, alpha=(2.38)^2, Sigma){
    chain <- startvalue
    accept_rate <- 0
    normal_shift <- mvrnorm(n=nb_iter, mu=rep(0, ndim), Sigma=(alpha/ndim)*Sigma)
    for (i in 1:nb_iter){
      #browser()
      if (i==1){previous_proposal = chain}else{previous_proposal = chain[i,]}
      proposal = previous_proposal + normal_shift[i,]
      tts1 <- is.dominant(Z.safe, proposal, ndim, 1)
      ttf1 <- is.dominant(Z.fail, proposal, ndim, 2)
      if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) & all(proposal>=0) & all(1-proposal>=0) ){
        chain <- rbind(chain, proposal)
        accept_rate <- accept_rate + 1
      }else{
        chain <- rbind(chain, previous_proposal)
      }
    }
    accept_rate <- accept_rate/nb_iter
    return(list(chain=mcmc(chain), accept_rate=accept_rate ))
  }
  
  
  #Utilis?e dans Metropolis-Hastings avec transformation
  Inverse_normal_transf <- function(x, mean=0, sd=1){
    return(qnorm(x , mean, sd))
  }
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Algo de Metropolis-Hastings (avec marche al?atoire dans espace image)
  #
  #-----------------------------------------------------------------------------------------------------------------#
  Inverse_normal_transf <- function(x, mean=0, sd=1){
    return(qnorm(x , mean, sd))
  }
  
  #matrice Sigma a la place de sigma ?
  Metropolis_Hastings<- function(startvalue, nb_iter,  Z.safe, Z.fail , sigma=(1/ndim)*(2.38)^2, Sigma=diag(ndim), transf_params=c("normal", 0, 1) ){ 
    chain <- startvalue
    mean= transf_params[2]
    sd = transf_params[3]
    accept_rate <- 0
    normal_shift <- mvrnorm(n=nb_iter, mu=rep(0, ndim), sigma*Sigma)
    for (i in 1:nb_iter){
      if (i==1){previous_proposal = chain}else{previous_proposal = chain[i,]}
      if (transf_params[1]=="normal") {
        proposal = sapply(previous_proposal, FUN=qnorm) + normal_shift[i,]
        proposal=sapply(proposal, FUN=pnorm) # on effectue la transform?e inverse
      }
      tts1 <- is.dominant(Z.safe, proposal, ndim, 1)
      ttf1 <- is.dominant(Z.fail, proposal, ndim, 2)
      if( (sum(tts1) == 0 ) & ( sum(ttf1) == 0) & all(proposal>= 5e-04) & all(1-proposal>= 5e-04)  ){
        #browser()
        # alpha =  log(abs(det(jacobian(Inverse_normal_transf, previous_proposal, method="simple") )))  -
        #   log(abs(det(jacobian(Inverse_normal_transf, proposal, method="simple") )))     
        det1=1;  det2=1 
        for (i in 1:ndim) {
          det1=det1*grad(Inverse_normal_transf, previous_proposal[i], method="simple")
          det2=det2*grad(Inverse_normal_transf, proposal[i], method="simple")
        }
        alpha = log(abs(det1))-log(abs(det2))
        alpha= min(alpha, 0)
      }else{
        alpha = -Inf}
      #print(paste("alpha=", alpha))
      if (is.na(alpha)) {alpha = -Inf}
      
      U=runif(n=1)
      #print(paste( "logU=", log(U)))
      if(log(U) <= alpha){
        chain <- rbind(chain, proposal)
        accept_rate <- accept_rate + 1
      }else{
        chain <- rbind(chain, previous_proposal)
      }
    }
    accept_rate <- accept_rate/nb_iter
    return(list(chain=mcmc(chain), accept_rate=accept_rate ))
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonction qui renvoie le param?tre d'?chelle pour Metropolis
  #
  #-----------------------------------------------------------------------------------------------------------------#
  #COMM_BASTIEN :On choisit le param?tre d'?chelle de la loi de la marche al?atoire dans Metropolis pour obtenir un taux d'acceptation de l'ordre de 1/4. 
  alpha_adapt <- function(startvalue,  Z.safe, Z.fail, transfMCMC=TRUE, alpha_start=(1/ndim)*(2.38)^2, Sigma=diag(ndim), transf_params = c("normal", 0, 1)  ){
    if (transfMCMC==TRUE) {
      alpha=alpha_start
      #Metropolis_Hastings<- function(startvalue, nb_iter,  Z.safe, Z.fail , sigma=(1/ndim)*(2.38)^2,  transf_params=c("normal", 0, 1) )
      chain_data <- Metropolis_Hastings(startvalue, nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail ,  sigma=alpha , Sigma=Sigma)
      accept_rate <- chain_data$accept_rate
      while ( (accept_rate<0.2) || (accept_rate>0.3) ) {
        chainstart <- chain_data$chain[dim(chain_data$chain)[1], ]
        alpha <- ifelse((accept_rate<0.2), (1/3)*alpha, 3* alpha )
        chain_data <- Metropolis_Hastings(chainstart, nb_iter=1000, Z.safe= Z.safe, Z.fail=Z.fail , sigma=alpha, Sigma=Sigma)
        accept_rate <- chain_data$accept_rate
      }
    }else{
      chainstart <- startvalue
      alpha=(1/ndim)*(2.38)^2
      chain_data <- RW_Metropolis(startvalue, 1000,  Z.safe, Z.fail , alpha=alpha)
      accept_rate <- chain_data$accept_rate
      while ( (accept_rate<0.2) || (accept_rate>0.3) ) {
        chainstart <- chain_data$chain[dim(chain_data$chain)[1], ]
        alpha <- ifelse((accept_rate<0.2), 3*alpha, (1/3)* alpha )
        chain_data <- RW_Metropolis(chainstart, 1000, transf_params,  Z.safe, Z.fail , alpha=alpha)
        accept_rate <- chain_data$accept_rate
      }

    }
    return(list(chain=chain_data$chain, alpha=alpha))

  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonction qui renvoie la sample covariance matrix d'un echantillon
  #
  #-----------------------------------------------------------------------------------------------------------------#
  #COMM_BASTIEN : data  correspond ? l'?chantillon des points, array (matrix) de dimension (nb_points , ndim)
  Sigma_estimate<- function(data){
    nb_points= dim(data)[1]
    M= matrix(data=data , nrow = nb_points , ncol = ndim)
    M= M - matrix(data = colMeans(M), nrow = nb_points, ncol = ndim , byrow = TRUE)
    return( (1/(nb_points-1))*(t(M)%*% M) )
  }
  
  
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonction qui renvoie la sample covariance matrix d'un echantillon
  #
  #-----------------------------------------------------------------------------------------------------------------#
  #COMM_BASTIEN : data  correspond ? l'?chantillon des points, array (matrix) de dimension (nb_points , ndim)
  Sigma_estimate<- function(data, regularize=FALSE){
    nb_points= dim(data)[1]
    M= matrix(data=data , nrow = nb_points , ncol = ndim)
    M= M - matrix(data = colMeans(M), nrow = nb_points, ncol = ndim , byrow = TRUE)
    if (regularize==FALSE) {
      return( (1/(nb_points-1))*(t(M)%*% M) )
    }else{
      return( 0.9*(1/(nb_points-1))*(t(M)%*% M) + 0.1*diag(ndim))
    }
  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonctions qui 
  #
  #-----------------------------------------------------------------------------------------------------------------#
  
  #Fonction qui renvoie un des points de la chaine de chain appartenant au nouvel espace non-domin? (chain doit ?tre du type ...)
  rej_sampl_nondom <- function(chain, start_point, Z.safe, Z.fail){
    #browser()
    cpt <- start_point
    M <- dim(chain)[1]
    if(cpt>M){stop("startpoint not in chain")}
    #if(chain){warning("")}
    
    cpt <- start_point
    chainstart <-NULL
    tts1 <- 1 ; ttf1 <- 1
    while( ((sum(tts1) != 0 ) || ( sum(ttf1) != 0)) & (cpt<M) ){
      tts1 <-is.dominant(Z.safe, chain[cpt,], ndim, 1)
      ttf1 <-is.dominant(Z.fail,chain[cpt,], ndim, 2)
      cpt=cpt+1
    }
    if(cpt==M){stop("No value found by rejection sampling")
    }else{val <- matrix(chain[cpt-1,], ncol=ndim)}
    return(list(val=val, cpt=cpt-1))
  }
  
  
  
  #Fonction qui renvoie un des points de la chaine de chain appartenant ? la zone lab?lis?e par le classifieur 
  # (set = 1 si on est dans hat{U}_k^+, 0 sinon) au sein de l'espace non-domin? (chain doit ?tre du type ...)
  rej_sampl_class <- function(chain, start_point, Z.safe, Z.fail, choix.classif, weights, set){
    cpt <- start_point
    M <- dim(chain)[1]
    if(cpt>M){stop("startpoint not in chain")}
    #if(chain){warning("")}
    
    cpt <- start_point
    chainstart <-NULL
    tts1 <- 1 ; ttf1 <- 1
    class_pred <- 2
    while( ( (sum(tts1) != 0)  || ( sum(ttf1) != 0) || (class_pred != set) ) & (cpt<M) ){
      cpt=cpt+1
      tts1 <-is.dominant(Z.safe, chain[cpt,], ndim, 1)
      ttf1 <-is.dominant(Z.fail,chain[cpt,], ndim, 2)
      class_pred <- classifieur.pred( matrix(chain[cpt,], ncol=ndim), choix.classif, weights)
    }
    if(cpt==M){stop("No value found by rejection sampling")
    }else{val <- matrix(chain[cpt,], ncol=ndim)}
    return(list(val=val, cpt=cpt))
  }
  
  #-----------------------------------------------------------------------------------------------------------------#
  #
  #		      Fonction qui
  #
  #-----------------------------------------------------------------------------------------------------------------#
  
  #Fonction qui renvoie un des points de la chaine de chain appartenant au nouvel espace non-domin? (chain doit ?tre du type ...)
  Gelman_Rubins <- function(chainstart, alpha_start,  Z.safe, Z.fail){
    #browser()
    chain1 <- NULL ; chain2 <- NULL ; chain3 <- NULL ; chain1supp <- NULL ; chain2supp <- NULL ; chain3supp <- NULL
    chain1 = Metropolis_Hastings(startvalue=chainstart[[1]]$val, nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain  
    chain2 = Metropolis_Hastings(startvalue=chainstart[[2]]$val, nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain
    chain3 = Metropolis_Hastings(startvalue=chainstart[[3]]$val, nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain
    combinedchains = mcmc.list(chain1, chain2, chain3)
    mpsrf <- gelman.diag(combinedchains)$mpsrf
    
    cpt <- 1
    while ((mpsrf>1.1) & (cpt < 10)) {
      chain1supp = Metropolis_Hastings(startvalue=chain1[dim(chain1)[1], ], nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain  
      chain2supp = Metropolis_Hastings(startvalue=chain1[dim(chain2)[1], ], nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain
      chain3supp = Metropolis_Hastings(startvalue=chain1[dim(chain3)[1], ], nb_iter=1000, Z.safe=Z.safe, Z.fail=Z.fail, sigma=alpha_start)$chain
      chain1 = mcmc(rbind( as.matrix(chain1), as.matrix(chain1supp))) 
      chain2 = mcmc(rbind( as.matrix(chain2), as.matrix(chain2supp)))
      chain3 = mcmc(rbind( as.matrix(chain3), as.matrix(chain3supp)))
      
      combinedchains = mcmc.list(chain1, chain2, chain3)
      mpsrf <- gelman.diag(combinedchains)$mpsrf
      cpt <- cpt +1
    }
    
    if(cpt==10){print("Stationnarity not reached after 30000 iterations")}
    return(list(chain1 = chain1, chain2= chain2, chain3 = chain3 ))
    #return(combinedchains)
    
  }  

  
  

  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  #
  # 		     				   PARTIE 3 : ESTIMATEUR MRM IS-MCMC
  #
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  Choix.method.3 <- function(N.appels, H, choix.classif, N.dicho, mcmc_period){
    
    # cp : compteur d'appel ? H (=G)
    # CV : vecteur contenant les coefficients de variations
    # choix.classif : choix du classifieur
    start_time <- Sys.time()
    um <- 0               
    uM   <- 1 
    Um <- 0               #Um : vecteur contenant les bornes exactes inf?rieures
    UM <- 1
    eps <- 1e-7
    alpha <- 0.025        # CHOIX DE ALPHA ?
    SIGN  <- 0            # SIGN : Signature des points test?s
    val_f <- 0
    p.bar <- 0            # uM : vecteur contenant les bornes exactes sup?rieures
    Z.safe<- NULL         #Fronti?re de U_k+
    Z.fail<- NULL         #Fronti?re de U_k-
    List.Z.safe <- list() # Liste des ensembles Z.safe selon l'it?ration
    List.Z.fail <- list() # Liste des ensembles Z.fail selon l'it?ration
    Classif  <- list()  #liste de listes contenant les donn?es du classifieur ? chaque ?tape
    Volnondom <- list()
    Volclassif <- list()
    ZS <- NULL            #Contient  TOUT les points ?valu?s jusqu'a alors qui conduisent a la d?faillance
    ZF <- NULL            # Pareil pour la non-d?faillance
    W <- 1:(2^(ndim) - 1)
    
    
    
    #DICHOTOMIE ET INITIALISATION DES BORNES/DU CLASSIFIEUR
    if(TRUE){
      V <- list() # V contient le premier tri de points avec la methode de la dichotomie.
      V <- Intersect(ndim, H, N.dicho)    
      list.set <- 0
      for(i in 1:length(V)){
        list.set[i] <- V[[i]][[2]] }
      #initialisation de Z.fail et Z.safe a partir des points calcul?s (V) pdt la dichotomie initiale.
      u.dep   <- list()
      u.other <- list()
      u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
      Z.fail <- t(as.matrix(u.dep[[1]]))
      Z.safe <- t(as.matrix(u.other[[1]]))
      Z.safe_init <- Z.safe  
      Z.fail_init <- Z.fail
      ZS <- matrix(data = u.other[[1]] , ncol=ndim, byrow=TRUE)  
      ZF <- matrix(data = u.dep[[1]] , ncol=ndim, byrow=TRUE)  
      cp     <- length(V)  
      #initialisation des bornes a partir du 1er elt de U- (resp U+) caclul? par dichotomie
      um <- prod(V[[cp]][[1]]) 
      uM   <- 1 - prod(1 - V[[cp]][[1]])
      classif_init <- classifieur.train(choix.classif, ZF, ZS)
    }
    #Cr?ation du tirage uniforme (pour initaliser les realisations dans l'espace non-domin? par la methode du rejet)
    UU <- runif(ndim*10^(ordre.p + 2))  
    UU <- matrix(UU, ncol=ndim, byrow = TRUE)
    indep_chain <- NULL # Correspond ? un ?chantilonnage i.i.d ? l'?tape 1
    for(i in 1:dim(UU)[1]) {
      tts1 <- is.dominant(Z.safe, UU[i,], ndim, 1) #ou l'inverse ?
      ttf1 <- is.dominant(Z.fail,  UU[i,], ndim, 2)
      if ((sum(tts1)==0) & (sum(ttf1)==0) ) { 
        indep_chain <- rbind(indep_chain, UU[i,] )
      }
    }
    ttsf <- rep(1,dim(indep_chain)[1])
    CPT_CHAIN <- 1
    #initialisation des volumes U_0 et ?_0
    volnondom_chain <-  1 -( prod(1-Z.safe[1, ]) + prod(Z.fail[1, ]) )
    volclass_init <-  VOLCLASS(UU, choix.classif, classif_init, Z.safe=Z.safe, Z.fail= Z.fail ) # A VERIFIER ! 
    if (volclass_init>=volnondom_chain) {warning("volclass_init (MC) >= >=volnondom_init (exact_calc) ")}
    time_0 <- Sys.time() ; cat(paste(" initialisation time : ", as.double(difftime(time_0, start_time, tz,units = "secs")), "s", "\n"  ));flush.console()
    
    
    M <- 2*1e4 
    j = 1
    
    
    while(cp < N.appels){
      #browser()
      
      #Simulation d'une nouvelle observation X_n+k-1 par la methode du rejet a partir de indep_chain
      ru <- runif(n = 1)
      if(j==1){
        old_weights=classif_init
      }else{
        old_weights=Classif[[j-1]]}
      if (ru<1-epsilon){
        tmp <- rej_sampl_class(chain=indep_chain, start_point = CPT_CHAIN , Z.safe, Z.fail, choix.classif, weights=old_weights, set=1)
      }else{
        tmp <- rej_sampl_class(chain=indep_chain, start_point = CPT_CHAIN , Z.safe, Z.fail, choix.classif, weights=old_weights, set=0)
      }
      u <- tmp$val
      CPT_CHAIN <- tmp$cpt+1
      
      
      # On rajoute u ? l'ensemble auquel il appartient. On garde sa signature et on calcule le volume du nouvel ensemble construit
      # On trie les points de l'espace non domin?:
      H.u = H(u)
      if(H.u > 0){
        ss     <- is.dominant(Z.safe, u, ndim, 2)
        # on enleve les points definissant Z.safe qui sont domin?s par le nouveau point u 
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(Z.safe, u )
        vol <- VOLVOL(UU,Z.safe, 1)
        if(j == 1){
          Um[j] <- um
          if(vol >= uM){
            UM[j] <- uM
          }else{
            UM[j] <- vol
          }
        }else{
          Um[j] <- Um[j - 1]
          if(vol >= UM[j-1]){
            UM[j] <- UM[j-1]
          }else{
            UM[j] <- vol
          }
        }
        
        SIGN[j] <- 0
        ZS<- rbind(ZS, u)
        
      }else{
        ff     <- is.dominant(Z.fail, u, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(Z.fail, u)
        vol = VOLVOL(UU, Z.fail, 2)
        
        if(j == 1){
          UM[j] <- uM
          if(vol <= um){
            Um[j] <- um
          }else{
            Um[j] <- vol
          }
        }else{
          UM[j] <- UM[j - 1]
          if(vol <= Um[j-1]){
            Um[j] = Um[j-1]
          }else{
            Um[j] <- vol
          }
        }
        
        
        SIGN[j] <- 1
        ZF<- rbind(ZF, u)
      }
      
      List.Z.safe[[j]] <- Z.safe
      List.Z.fail[[j]] <- Z.fail
      
      Classif[[j]] <- classifieur.train(choix.classif, ZF, ZS)            #Creation du nouveau classifieur
      
      #CALCUL DES BORNES ET DES ESTIMATEURS
      if (j<=mcmc_period) {
        for (i in 1:dim(indep_chain)[1]){ttsf[i] <- ttsf[i]*ifelse(is.dominant(u ,indep_chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 )}
        Volnondom[[j]] <- (sum(ttsf==1)/dim(indep_chain)[1])*volnondom_chain
        ttsfc <- ttsf
        for (i in 1:dim(indep_chain)[1]){ ttsfc[i] <- ttsfc[i]*ifelse(classifieur.pred( matrix(indep_chain[i,], ncol=ndim), choix.classif , Classif[[j]] )==1, 1, 0)  }
        Volclassif[[j]] <- (sum(ttsfc==1)/M)*volnondom_chain
      }else{
        for (i in 1:M){ ttsf[i] <- ttsf[i]*ifelse(is.dominant(u ,chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 ) }
        Volnondom[[j]] <- (sum(ttsf==1)/M)*volnondom_chain
        ttsfc <- ttsf
        for (i in 1:M){ ttsfc[i] <- ttsfc[i]*ifelse(classifieur.pred( matrix(chain[i,], ncol=ndim), choix.classif , Classif[[j]] )==1, 1, 0)  }
        Volclassif[[j]] <- (sum(ttsfc==1)/M)*volnondom_chain
      }
      if (Volclassif[[j]]>=Volnondom[[j]]) {print(" ERROR : Volclassif[[j]]  >= Volnondom[[j]] ")}
      cat(paste("Volnondom",j, "=", round(Volnondom[[j]], digits=5), "  Volclassif",j, "=", round(Volclassif[[j]], digits=5), "\n"));flush.console()
      
      
      if (j==1) {
        UM[j] <- uM-(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)
        Um[j] <- um+(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)
      }else{
        UM[j] <- UM[j-1]-(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)    #(1-UM[j]) = (1-UM[j-1])+(volnondom[j-1]-volnondom[j])*ifelse(SIGN[j]==0, 1, 0)
        Um[j] <- Um[j-1]+(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)    #equiv. ? Um[j] <- UM[j]-Volnondom[[j]] ?
      }
      
      #CALCUL DES ESTIMATEURS P.BAR
      if (j==1) {
        val_f[j] <- epsilon*(1/(volnondom_chain-volclass_init))*ifelse(classifieur.pred(u, choix.classif , classif_init)==0, 1, 0) +
          (1-epsilon)*(1/(volclass_init))*ifelse(classifieur.pred(u, choix.classif , classif_init)==1, 1, 0)
        p.bar[j] <-  uM - (1/(val_f[j]))*ifelse( SIGN[j]==0, 1, 0)
      }else{
        val_f[j] <- epsilon*(1/(Volnondom[[j-1]]-Volclassif[[j-1]]))*ifelse(classifieur.pred(u, choix.classif , Classif[[j-1]])==0, 1, 0) +
          (1-epsilon)*(1/(Volclassif[[j-1]] ))*ifelse(classifieur.pred(u, choix.classif , Classif[[j-1]])==1, 1, 0)
        p.bar[j] <-  UM[j-1] - (1/(val_f[j]))*ifelse( SIGN[j]==0, 1, 0)
        }
      
      
      
      #CONSTRUCTION DE LA CHAINE DE MARKOV 
      # TANT QUE rho= Vol(U_n+k)/Vol(U_n)>0.? et CPT_CHAIN<dim(chain), on a pas besoin de cr?er une nouvelle cha?ne de Markov !
      # CPDT le probl?me c'est que l'erreur de pk+ et celle de pk+1 est corr?l?e car on calcule a partir de la m?me trajectoire de la cha?ne!
      # C'est POURQUOI il faut mettre a jour la trajectoire malgr? tout p?riodiquement  (mcmc_period = 40 ?)
      
      if ((j%%mcmc_period == 0) ){
        
        #browser()
        time1<-Sys.time()
        volnondom_chain <- Volnondom[[j]]
        Sigma_choice <- diag(ndim) #Sigma_choice=Sigma_estimate(indep_chain)  # or Y[n-p])=combinedchains
        chainstart <- rej_sampl_nondom(indep_chain, start_point = 1, Z.safe, Z.fail)$val
        
        M <- 2*M                        #car on divise par 2 apr?s
        if(transfMCMC==FALSE){
          #chain = RW_Metropolis(startvalue=chainstart, nb_iter=M, Z.safe=Z.safe, Z.fail=Z.fail, Sigma=Sigma_choice)$chain
          
          # M <- 1.5*M
          # chain1 = RW_Metropolis(startvalue=chainstart[1,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # chain2 = RW_Metropolis(startvalue=chainstart[2,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # chain3 = RW_Metropolis(startvalue=chainstart[3,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # combinedchains = mcmc.list(chain1, chain2, chain3) 
          # #plot(combinedchains)
          # #gelman.diag(combinedchains)
          # #gelman.plot(combinedchains)
          # 
          # #ON FAIT LE BURN-IN POUR LES 3 CHAINES DE LA COMBINED CHAIN !!
          # #PAS BESOIN SI ON DEBUTE DANS UNE ZONE DE HAUTE DENSITE
          # m=as.integer(M/2)
          # cat(paste("  Burn-In ="));flush.console()
          # chain1=window(chain1, start=m, thin=1)
          # chain2=window(chain2, start=m, thin=1)
          # chain3=window(chain3, start=m, thin=1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # chain = mcmc(rbind( as.matrix(chain1), as.matrix(chain2), as.matrix(chain2) )) 
          # 
          # #Deduction d'un echantillon quasi-indep. par echantill. par paquets (de taille t=10 ou autre ?)
          # indep_chain1=window(chain1, thin=10)
          # indep_chain2=window(chain2, thin=10)
          # indep_chain3=window(chain3, thin=10)
          # indep_chains=mcmc(rbind( as.matrix(indep_chain1), as.matrix(indep_chain2), as.matrix(indep_chain2) ) )
          # 
        }else{
          #browser() #On choisit un param?tre d'?chelle proche de 1/4 par dichotomie avec alpha_adapt
          test_chain_data <- alpha_adapt(startvalue=chainstart,  Z.safe=Z.safe, Z.fail=Z.fail, transfMCMC=TRUE)
          time2<-Sys.time()
          chainstart <- test_chain_data$chain[dim(test_chain_data$chain)[1], ]
          chain_data <- Metropolis_Hastings(startvalue=chainstart, nb_iter=M, Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha )
          chain <- mcmc(rbind( as.matrix(test_chain_data$chain), as.matrix(chain_data$chain)))
          accept_rate <- chain_data$accept_rate
        }
        #ON FAIT LE BURN-IN
        chain=window(chain, start=as.integer(M/2), thin=1)  #chainbis=matrix(data = chain, ncol=ndim)
        ACT=as.integer(dim(chain)[1]/effectiveSize(chain))+1
        if (TRUE) { # La variance de l'estimateur MC classique pour estimer un volume rho=Vol(V)/Vol(U_n) est rho(1-rho)/M. 
          # Si on a un echantillon (quasi-)indep unif. sur U_n, l'estimateur de MC classique aura ainsi un Coeff de Var = sqrt((1-rho)/rho)*sqrt(1/M)
          # Dans le cas o? rho>=1/5, alors CV <= 2*sqrt(1/M) (ce qui conduit ? M=1600)
          # Dans le cas o? rho>=1/2, alors CV <= sqrt(1/M)
          # Dans le cas o? rho>=3/4, alors CV <= 0.58*sqrt(1/M)  (ce qui conduit ? M=133)
          # Dans le cas o? rho>=9/10, alors CV <= 0.33*sqrt(1/M)  (ce qui conduit ? M=45)
          # On veut donc prendre sqrt(1/M)=0.05, soit M=400.  
          # Il faut ensuite multiplier M par l'ACT si on utilise le thm ergodique. (on utilise la moyenne des 3 estimatuers caclul?s pour chque chaine, ce qui divise la variance par 3 et l'?cart-type par sqrt(3)). 
        }
        M_MC <- 100
        cat(paste("  ACT =", max(ACT), "\n", "sigma =", test_chain_data$alpha,  "accept_rate=", accept_rate, "\n",  " size of chain  (old_M)", dim(chain)[1], " size needed  (new_M)", M_MC *max(ACT),  "\n"  ));flush.console()
        if (M_MC *max(ACT)<=dim(chain)[1]) {
          chain <- chain[1:(M_MC *max(ACT)), ]
        }else{
          chain_supp = Metropolis_Hastings(startvalue=chain[dim(chain)[1], ], nb_iter=M_MC*max(ACT)-dim(chain)[1], Z.safe=Z.safe, Z.fail=Z.fail)$chain
          chain=mcmc(rbind( as.matrix(chain), as.matrix(chain_supp)))
        }
        indep_chain=window(chain, thin=40)
        M <- dim(chain)[1]
        ttsf <- rep(1,M)
        CPT_CHAIN <- 1
        time3 <- Sys.time()
        cat(paste(" choosing chainstart and scale param for MC time : ", as.double(difftime(time2, time1, tz,units = "secs")), "s", "\n"  ));flush.console()
        cat(paste(" building MC time : ", as.double(difftime(time3, time2, tz,units = "secs")), "s", "\n"  ));flush.console()
        
        
        
        
      }
      
      j <- j+1
      cp <- cp+1
      
      
    } #fin du while cp < N.appels
    
    
    
    print("On est dans la methode MRM_IS_MCMC ")
    RR1 <- cbind(Um, UM, Volnondom, Volclassif,   p.bar)
    RR <- list(tab=RR1, Classif=Classif, List.Z.fail=List.Z.fail, List.Z.safe=List.Z.safe, ZF=ZF, ZS=ZS )
    return(RR) 

    
  }
  
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  # 					               FIN	METHODE MRM IS-MCMC 
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  #
  # 		     				   PARTIE 4 : ESTIMATEUR MRM MLE-MCMC 
  #
  #
  #
  #-----------------------------------------------------------------------------------------------------------------
  
  Choix.method.4 <- function(N.appels, H, choix.classif, N.dicho, mcmc_period){
    
    # cp : compteur d'appel ? H (=G)
    # MLE :Estimateur du maximum de vraissemblance
    # CV : vecteur contenant les coefficients de variations
    # SIGN : Signature des points test?s
    # choix.classif : choix du classifieur
    start_time <- Sys.time()
    Um <- 0
    UM <- 1
    eps <- 1e-7
    MLE <- 0
    SIGN  <- 0
    beta <- 0             # vecteur contenant la suite des rapports des volumes non-domin?s (estim?s)
    um <- 0               # um : vecteur contenant les bornes exactes inf?rieures
    uM   <- 1             # uM : vecteur contenant les bornes exactes sup?rieures
    Z.safe<- NULL         #Fronti?re de U_k+
    Z.fail<- NULL         #Fronti?re de U_k-
    List.Z.safe <- list() # Liste des ensembles Z.safe selon l'it?ration
    List.Z.fail <- list() # Liste des ensembles Z.fail selon l'it?ration
    Volnondom <- list()
    ZS <- NULL            #Contient  TOUT les points ?valu?s jusqu'a alors qui conduisent a la d?faillance
    ZF <- NULL            # Symm?trique
    W <- 1:(2^(ndim) - 1)
    
    
    
    #DICHOTOMIE ET INITIALISATION DES BORNES/DU CLASSIFIEUR
    if(TRUE){
      V <- list() # V contient le premier tri de points avec la methode de la dichotomie.
      V <- Intersect(ndim, H, N.dicho)    
      list.set <- 0
      for(i in 1:length(V)){
        list.set[i] <- V[[i]][[2]] }
      #initialisation de Z.fail et Z.safe a partir des points calcul?s (V) pdt la dichotomie initiale.
      u.dep   <- list()
      u.other <- list()
      u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
      Z.fail <- t(as.matrix(u.dep[[1]]))
      Z.safe <- t(as.matrix(u.other[[1]]))
      Z.safe_init <- Z.safe  
      Z.fail_init <- Z.fail
      ZS <- matrix(data = u.other[[1]] , ncol=ndim, byrow=TRUE)  
      ZF <- matrix(data = u.dep[[1]] , ncol=ndim, byrow=TRUE)  
      cp     <- length(V)  
      #initialisation des bornes a partir du 1er elt de U- (resp U+) caclul? par dichotomie
      um <- prod(V[[cp]][[1]]) 
      uM   <- 1 - prod(1 - V[[cp]][[1]])
    }
    #Cr?ation du tirage uniforme (pour initaliser les realisations dans l'espace non-domin? par la methode du rejet)
    UU <- runif(ndim*10^(ordre.p + 2))  
    UU <- matrix(UU, ncol=ndim, byrow = TRUE)
    indep_chain <- NULL # Correspond ? un ?chantilonnage i.i.d ? l'?tape 1
    for(i in 1:dim(UU)[1]) {
      tts1 <- is.dominant(Z.safe, UU[i,], ndim, 1) #ou l'inverse ?
      ttf1 <- is.dominant(Z.fail,  UU[i,], ndim, 2)
      if ((sum(tts1)==0) & (sum(ttf1)==0) ) { 
        indep_chain <- rbind(indep_chain, UU[i,] )
      }
    }
    ttsf <- rep(1,dim(indep_chain)[1])
    CPT_CHAIN <- 1
    #initialisation des volumes U_0 et ?_0
    #browser()
    volnondom_chain <-  1 -( prod(1-Z.safe[1, ]) + prod(Z.fail[1, ]) )
    
    time_1 <- Sys.time() ; cat(paste(" initialisation time : ", time_1 - start_time   ));flush.console()
    
    M <- 2*1e4 
    j = 1
    
    
    
    time100 <-0 ;  Time100 <-0
    while(cp < N.appels){
      if (j==110) {Time100 <- Sys.time() ; time100 <-as.double(difftime(Time100, start_time, tz,units = "secs"))}
      
      #cat(paste("  j =",j,  "\n"));flush.console()
      #Simulation d'une nouvelle observation X_n+k-1 par la methode du rejet a partir de indep_chain
      an.error.occured <- FALSE
      tryCatch( { tmp<-rej_sampl_nondom(chain = indep_chain, start_point = CPT_CHAIN, Z.safe, Z.fail) }
                , error = function(e) {an.error.occured <<- TRUE ; print("ERROR")})
      if (an.error.occured) {return(NULL)}
      
      
      
      
      u <- tmp$val
      CPT_CHAIN <- tmp$cpt+1
      
      
      # On rajoute u ? l'ensemble auquel il appartient. On garde sa signature et on calcule le volume du nouvel ensemble construit
      # On trie les points de l'espace non domin?:
      #cat(paste(" u =", round(u[1, ], digits=6); "\n"));flush.console()
      H.u = H(u)
      if(H.u > 0){
        ss     <- is.dominant(Z.safe, u, ndim, 2)
        # on enleve les points definissant Z.safe qui sont domin?s par le nouveau point u 
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(Z.safe, u )
        vol <- VOLVOL(UU,Z.safe, 1)
        if(j == 1){
          Um[j] <- um
          if(vol >= uM){
            UM[j] <- uM
          }else{
            UM[j] <- vol
          }
        }else{
          Um[j] <- Um[j - 1]
          if(vol >= UM[j-1]){
            UM[j] <- UM[j-1]
          }else{
            UM[j] <- vol
          }
        }
        
        SIGN[j] <- 0
        ZS<- rbind(ZS, u)
        
      }else{
        ff     <- is.dominant(Z.fail, u, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(Z.fail, u)
        vol = VOLVOL(UU, Z.fail, 2)
        
        if(j == 1){
          UM[j] <- uM
          if(vol <= um){
            Um[j] <- um
          }else{
            Um[j] <- vol
          }
        }else{
          UM[j] <- UM[j - 1]
          if(vol <= Um[j-1]){
            Um[j] = Um[j-1]
          }else{
            Um[j] <- vol
          }
        }
        
        
        SIGN[j] <- 1
        ZF<- rbind(ZF, u)
      }
      
      List.Z.safe[[j]] <- Z.safe
      List.Z.fail[[j]] <- Z.fail
      
 
      #CALCUL DES ESTIMATEURS 
      if (j<=mcmc_period) {
        for (i in 1:dim(indep_chain)[1]){ttsf[i] <- ttsf[i]*ifelse(is.dominant(u ,indep_chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 )}
        Volnondom[[j]] <- (sum(ttsf==1)/dim(indep_chain)[1])*volnondom_chain
      }else{
        for (i in 1:M){
          ttsf[i] <- ttsf[i]*ifelse(is.dominant(u ,chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 )
          }
        Volnondom[[j]] <- (sum(ttsf==1)/M)*volnondom_chain
        }
      #Comment calculer UM[j] ? (1-UM[j]) = (1-UM[j-1])+(volnondom[j-1]-volnondom[j])*ifelse(SIGN[j]==0, 1, 0)
      if (j==1) {
        UM[j] <- uM-(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)
        Um[j] <- um+(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)
      }else{
        UM[j] <- UM[j-1]-(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)
        #Um[j] <- UM[j]-Volnondom[[j]]
        Um[j] <- Um[j-1]+(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)
      }
      
      MLE.test <- optimize(f = log.likehood,
                           interval = c(Um[j],UM[j]),
                           maximum = TRUE,
                           tol = 10^(-3*ordre.p),
                           signature = SIGN,
                           p.k = cbind(Um, UM)    )
      MLE[j] <- as.numeric(MLE.test[1])
      cat(paste("  volnondom",j, "=", round(Volnondom[[j]], digits=6), "  MLE",j, "=", round(MLE[j], digits=6),  "\n"));flush.console()
      
      
      #CONSTRUCTION DE LA CHAINE DE MARKOV
      # TANT QUE rho= Vol(U_n+k)/Vol(U_n)>0.5 et CPT_CHAIN<dim(chain), on a pas besoin de cr?er une nouvelle cha?ne de Markov !
      # CPDT le probl?me c'est que l'erreur de pk+ et celle de pk+1 est corr?l?e car on calcule a partir de la m?me trajectoire de la cha?ne!
      # C'est POURQUOI il faut mettre a jour la trajectoire malgr? tout p?riodiquement  (mcmc_period = 40 ?)
      if ((j%%mcmc_period == 0) ){
        
        time1<-Sys.time()
        old_volnondom_chain <-  volnondom_chain 
        volnondom_chain <- Volnondom[[j]]
        Sigma_choice <- diag(ndim) #Sigma_choice=Sigma_estimate(indep_chain)  # or Y[n-p])=combinedchains
        chainstart <- rej_sampl_nondom(indep_chain, start_point = 1, Z.safe, Z.fail)$val
        
        M <- 2*M                        #car on divise par 2 apr?s
        if(transfMCMC==FALSE){
          #chain = RW_Metropolis(startvalue=chainstart, nb_iter=M, Z.safe=Z.safe, Z.fail=Z.fail, Sigma=Sigma_choice)$chain
          
          # M <- 1.5*M
          # chain1 = RW_Metropolis(startvalue=chainstart[1,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # chain2 = RW_Metropolis(startvalue=chainstart[2,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # chain3 = RW_Metropolis(startvalue=chainstart[3,], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, alpha=(2.38)^2, Sigma=Sigma_choice)
          # combinedchains = mcmc.list(chain1, chain2, chain3) 
          # #plot(combinedchains)
          # #browser()
          # #gelman.diag(combinedchains)
          # #gelman.plot(combinedchains)
          # 
          # #ON FAIT LE BURN-IN POUR LES 3 CHAINES DE LA COMBINED CHAIN !!
          # #PAS BESOIN SI ON DEBUTE DANS UNE ZONE DE HAUTE DENSITE
          # m=as.integer(M/2)
          # cat(paste("  Burn-In ="));flush.console()
          # chain1=window(chain1, start=m, thin=1)
          # chain2=window(chain2, start=m, thin=1)
          # chain3=window(chain3, start=m, thin=1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # ACT1=dim(chain1)[1]/effectiveSize(chain1)
          # chain = mcmc(rbind( as.matrix(chain1), as.matrix(chain2), as.matrix(chain2) )) 
          # 
          # #Deduction d'un echantillon quasi-indep. par echantill. par paquets (de taille t=10 ou autre ?)
          # indep_chain1=window(chain1, thin=10)
          # indep_chain2=window(chain2, thin=10)
          # indep_chain3=window(chain3, thin=10)
          # indep_chains=mcmc(rbind( as.matrix(indep_chain1), as.matrix(indep_chain2), as.matrix(indep_chain2) ) )
          # 
        }else{
          #browser() #On choisit un param?tre d'?chelle proche de 1/4 par dichotomie avec alpha_adapt
          
          test_chain_data <- alpha_adapt(startvalue=chainstart,  Z.safe=Z.safe, Z.fail=Z.fail,  transfMCMC=TRUE)
          time2<-Sys.time()
          chainstart <- test_chain_data$chain[dim(test_chain_data$chain)[1], ]
          chain_data <- Metropolis_Hastings(startvalue=chainstart, nb_iter=M, Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha )
                                    chain <- mcmc(rbind( as.matrix(test_chain_data$chain), as.matrix(chain_data$chain)))
          accept_rate <- chain_data$accept_rate
        }
        #ON FAIT LE BURN-IN
        chain=window(chain, start=as.integer(M/2), thin=1)  #chainbis=matrix(data = chain, ncol=ndim)
        ACT=as.integer(dim(chain)[1]/effectiveSize(chain))+1
        #browser()
        #ACT <- max (ACT)   #ou quantile(ACT)[4] si on veut le 3e quartile slmnt (en grande dim); prendre[3] pour la mediane
        ACT <- mean(ACT)
        # La variance de l'estimateur MC classique pour estimer un volume rho=Vol(V)/Vol(U_n) est rho(1-rho)/M. 
        # Si on a un echantillon (quasi-)indep unif. sur U_n, l'estimateur de MC classique aura ainsi un Coeff de Var = sqrt((1-rho)/rho)*sqrt(1/M)
        # Dans le cas o? rho>=1/5, alors CV <= 2*sqrt(1/M) (ce qui conduit ? M=1600)
        # Dans le cas o? rho>=1/2, alors CV <= sqrt(1/M)
        # Dans le cas o? rho>=3/4, alors CV <= 0.58*sqrt(1/M)  (ce qui conduit ? M=133)
        # Dans le cas o? rho>=9/10, alors CV <= 0.33*sqrt(1/M)  (ce qui conduit ? M=45)
        # On veut donc prendre sqrt(1/M)=0.05, soit M=400.  
        # Il faut ensuite multiplier M par l'ACT si on utiliqe le thm ergodique. (on utilise la moyenne des 3 estimatuers caclul?s pour chque chaine, ce qui divise la variance par 3 et l'?cart-type par sqrt(3)).
        M_MC <- 60
        cat(paste("  ACT =", ACT, "\n", "accept_rate=", accept_rate, "\n",  " size of chain  (old_M)", dim(chain)[1], " size needed  (new_M)", M_MC *ACT, "sigma=", test_chain_data$alpha, "\n"  ));flush.console()
        if (ACT>=1000) { ACT <- 1000 ;  cat(paste(" size fixed  (new_M)", M_MC *ACT, "\n"  ));flush.console()  }
        if (M_MC *ACT<=dim(chain)[1]) {
          chain <- chain[1:(M_MC *ACT), ]
        }else{
          chain_supp = Metropolis_Hastings(startvalue=chain[dim(chain)[1], ], nb_iter=M_MC*ACT-dim(chain)[1], Z.safe=Z.safe, Z.fail=Z.fail)$chain
          chain=mcmc(rbind( as.matrix(chain), as.matrix(chain_supp)))
        }
        
        #indep_chain=window(chain, thin=ACT)
        indep_chain=window(chain, thin = ACT) 
        M <- dim(chain)[1]
        ttsf <- rep(1,M)
        CPT_CHAIN <- 1
        beta <- c(beta, (old_volnondom_chain - volnondom_chain)/(M_MC*volnondom_chain) )
        time3 <- Sys.time()
        cat(paste(" choosing chainstart and scale param for MC time : ", as.double(difftime(time2, time1, tz,units = "secs")), "s", "\n"  ));flush.console()
        cat(paste(" building MC time : ", as.double(difftime(time3, time2, tz,units = "secs")), "s", "\n"  ));flush.console()
      }
      
      j <- j+1
      cp <- cp+1
      

    } #fin du while cp < N.appels
    
    
    
    print("On est dans la methode MRM_MCMC ")
    RR1 <- cbind(Um, UM,  MLE)
    #browser()
    RR <- list(tab=RR1, List.Z.fail=List.Z.fail, List.Z.safe=List.Z.safe, ZF=ZF, ZS=ZS, beta=beta ) 
    #return(RR) 
    return(RR) 
  
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------
  #
  #
  #
  # 					               FIN	METHODE MRM MLE-MCMC 
  #
  #
  #-----------------------------------------------------------------------------------------------------------------

  
  Choix.method.5 <- function(N.appels, H, choix.classif, N.dicho, alpha_start){
    
    # cp : compteur d'appel ? H (=G)
    # MLE :Estimateur du maximum de vraissemblance
    # CV : vecteur contenant les coefficients de variations
    # SIGN : Signature des points test?s
    # choix.classif : choix du classifieur
    start_time <- Sys.time()
    Um <- 0
    UM <- 1
    MLE <- 0
    SIGN  <- 0
    beta <- 0             # vecteur contenant la suite des rapports des volumes non-domin?s (estim?s)
    um <- 0               # um : vecteur contenant les bornes exactes inf?rieures
    uM   <- 1             # uM : vecteur contenant les bornes exactes sup?rieures
    Z.safe<- NULL         #Fronti?re de U_k+
    Z.fail<- NULL         #Fronti?re de U_k-
    List.Z.safe <- list() # Liste des ensembles Z.safe selon l'it?ration
    List.Z.fail <- list() # Liste des ensembles Z.fail selon l'it?ration
    Volnondom <- list()
    ZS <- NULL            #Contient  TOUT les points ?valu?s jusqu'a alors qui conduisent a la d?faillance
    ZF <- NULL            # Symm?trique
    W <- 1:(2^(ndim) - 1)
    first_loop <- TRUE
    
    
    #DICHOTOMIE ET INITIALISATION DES BORNES/DU CLASSIFIEUR
    if(TRUE){
      V <- list() # V contient le premier tri de points avec la methode de la dichotomie.
      V <- Intersect(ndim, H, N.dicho)    
      list.set <- 0
      for(i in 1:length(V)){
        list.set[i] <- V[[i]][[2]] }
      #initialisation de Z.fail et Z.safe a partir des points calcul?s (V) pdt la dichotomie initiale.
      u.dep   <- list()
      u.other <- list()
      u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
      u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]
      u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
      u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]
      Z.fail <- t(as.matrix(u.dep[[1]]))
      Z.safe <- t(as.matrix(u.other[[1]]))
      Z.safe_init <- Z.safe  
      Z.fail_init <- Z.fail
      ZS <- matrix(data = u.other[[1]] , ncol=ndim, byrow=TRUE)  
      ZF <- matrix(data = u.dep[[1]] , ncol=ndim, byrow=TRUE)  
      cp     <- length(V)  
      #initialisation des bornes a partir du 1er elt de U- (resp U+) caclul? par dichotomie
      um <- prod(V[[cp]][[1]]) 
      uM   <- 1 - prod(1 - V[[cp]][[1]])
    }
    #Cr?ation du tirage uniforme (pour initaliser les realisations dans l'espace non-domin? par la methode du rejet)
    UU <- runif(ndim*10^(ordre.p + 2))  
    UU <- matrix(UU, ncol=ndim, byrow = TRUE)

    tts1 <- is.dominant.C(UU, Z.safe, ndim, 1)
    ttf1 <- is.dominant.C(UU, Z.fail, ndim, 2)
    indep_chain <- UU[which((tts1 == 0)&(ttf1 == 0) ) , ]
    
    ttsf <- rep(1,dim(indep_chain)[1])
    CPT_CHAIN <- 1
    #initialisation des volumes U_0 et ?_0
    #browser()
    volnondom_chain <-  1 -( prod(1-Z.safe[1, ]) + prod(Z.fail[1, ]) )
    
    time_1 <- Sys.time() ; cat(paste(" initialisation time : ", as.double(difftime(time_1, start_time, tz,units = "secs")), "s", "\n"  ));flush.console()
    
    
    M <- 2*1e4 
    j = 1
    j_chain = 1
    
    
    #time100 <-0 ;  Time100 <-0
    
    while(cp < N.appels){
      
      #if (j==100+ length(V)) {Time100 <- Sys.time() ; time100 <-as.double(difftime(Time100, start_time, tz,units = "secs"))}
      
      #Simulation d'une nouvelle observation X_n+k-1 par la methode du rejet a partir de indep_chain
      #On recalcule la chaine si on peut plus faire du rejection sampling ou si  Volnondom/Volnondom_chain<=1/2 (pour la pr?cision) ou si on a deja calcul? trop de sous-volumes (par ex 25) avec la meme trajectoire
      next_chain <- FALSE          #<
      if (j>1){ if(Volnondom[[j-1]]<0.5 *volnondom_chain) {next_chain <- TRUE ;  print("time to recompute chain 2") }      }
      #if (next_chain==FALSE & (j-j_chain==25) ) {  next_chain <- TRUE ;  print("time to recompute chain 1")     }
      if (next_chain==FALSE){if (first_loop & (j==30)) {  next_chain <- TRUE ;  print("time to recompute chain 1")     }}
      if (next_chain==FALSE) {                 
        tryCatch( { tmp<-rej_sampl_nondom(chain = indep_chain, start_point = CPT_CHAIN, Z.safe, Z.fail) }
                  , error = function(e) {next_chain <<- TRUE ;  print("time to recompute chain 3") })
      }
      

      if (next_chain==FALSE) {
        u <- tmp$val            #cat(paste(" u =", round(u[1, ], digits=6); "\n"));flush.console()
        CPT_CHAIN <- tmp$cpt+1
        
        # On rajoute u ? l'ensemble auquel il appartient. On garde sa signature et on calcule le volume du nouvel ensemble construit
        # On trie les points de l'espace non domin?:
        H.u = H(u)
        if(H.u > 0){
          ss     <- is.dominant(Z.safe, u, ndim, 2)
          # on enleve les points definissant Z.safe qui sont domin?s par le nouveau point u 
          Z.safe <- Z.safe[which(ss == FALSE), ]
          Z.safe <- rbind(Z.safe, u )
          SIGN[j] <- 0
          ZS<- rbind(ZS, u)
          
        }else{
          ff     <- is.dominant(Z.fail, u, ndim, 1)
          Z.fail <- Z.fail[which(ff == FALSE), ]
          Z.fail <- rbind(Z.fail, u)
          SIGN[j] <- 1
          ZF<- rbind(ZF, u)
        }
        
        List.Z.safe[[j]] <- Z.safe
        List.Z.fail[[j]] <- Z.fail
        #CALCUL DES ESTIMATEURS 
        if (first_loop) {
          # for (i in 1:dim(indep_chain)[1]){
          #   ttsf[i] <- ttsf[i]*ifelse(is.dominant(u ,indep_chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 )
          #   }
          ttsf <- ttsf*(1-is.dominant.C(indep_chain, u, ndim,(1-SIGN[j])+1))
          Volnondom[[j]] <- (sum(ttsf==1)/dim(indep_chain)[1])*volnondom_chain
          
          
          
        }else{
          # browser()
          # for (i in 1:M){
          #   ttsf[i] <- ttsf[i]*ifelse(is.dominant(u , chain[i,] , ndim, set=SIGN[j]+1)==TRUE , 0 , 1 )
          # }
          ttsf <- ttsf*(1-is.dominant.C(chain, u, ndim,(1-SIGN[j])+1))
          Volnondom[[j]] <- (sum(ttsf==1)/M)*volnondom_chain
          
        }
        #Comment calculer UM[j] ? (1-UM[j]) = (1-UM[j-1])+(volnondom[j-1]-volnondom[j])*ifelse(SIGN[j]==0, 1, 0)
        if (j==1) {
          UM[j] <- uM-(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)
          Um[j] <- um+(volnondom_chain-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)
        }else{
          UM[j] <- UM[j-1]-(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==0, 1, 0)
          #Um[j] <- UM[j]-Volnondom[[j]]
          Um[j] <- Um[j-1]+(Volnondom[[j-1]]-Volnondom[[j]])*ifelse(SIGN[j]==1, 1, 0)
        }
        
        # MLE.test <- optimize(f = log.likehood,
        #                      interval = c(Um[j],UM[j]),
        #                      maximum = TRUE,
        #                      tol = 10^(-3*ordre.p),
        #                      signature = SIGN,
        #                      p.k = cbind(Um, UM)    )
        # MLE[j] <- as.numeric(MLE.test[1])
        # cat(paste("  volnondom",j, "=", round(Volnondom[[j]], digits=6), "  MLE",j, "=", round(MLE[j], digits=6),  "\n"));flush.console()
        cat(paste("  compteur = ",cp, "--"));flush.console()

        j <- j+1
        cp <- cp+1
        
        
        }else{ 
          j_chain <-j          
          #(if next_chain == TRUE)
          #CONSTRUCTION DE LA NOUVELLE CHAINE DE MARKOV
          # TANT QUE rho= Vol(U_n+k)/Vol(U_n)>0.5 et CPT_CHAIN<dim(chain), on a pas besoin de cr?er une nouvelle cha?ne de Markov !
          # CPDT le probl?me c'est que l'erreur de pk+ et celle de pk+1 est corr?l?e car on calcule a partir de la m?me trajectoire de la cha?ne!
          # C'est POURQUOI il faut mettre a jour la trajectoire malgr? tout p?riodiquement  (M_MC = 50)
          
          
          time1<-Sys.time()
          old_volnondom_chain <-  volnondom_chain 
          volnondom_chain <- Volnondom[[j-1]]
          # tryCatch( { Sigma_choice <- Sigma_estimate(indep_chain, regularize = FALSE)}
          #           , error = function(e) {Sigma_choice <<- Sigma_estimate(indep_chain, regularize = TRUE) ;  print(paste("Sigma_estimate not able to compute : regularization needed", "\n") ) })
          Sigma_choice <- Sigma_estimate(indep_chain, regularize = TRUE)
          
          if(Burn_in=="GelmanRubin"){
            M <- 2*M # car on effectue un thinning de taille 2 ensuite
            chainstart <- list()
            #debug(rej_sampl_nondom)
            if (first_loop) {chain <- indep_chain}
            chainstart[[1]]<- rej_sampl_nondom(chain, start_point = 1, Z.safe, Z.fail)
            chainstart[[2]]<- rej_sampl_nondom(chain, start_point = chainstart[[1]]$cpt+1, Z.safe, Z.fail)
            chainstart[[3]]<- rej_sampl_nondom(chain, start_point = chainstart[[2]]$cpt+1, Z.safe, Z.fail)
            #unbug(rej_sampl_nondom)
            tryCatch({alpha_start <-  test_chain_data$alpha}
                     , error = function(e) {alpha_start <<- (1/ndim)*(2.38)^2 })
        
            cat(paste(" Gelman-Rubins "));flush.console()
            chains <- Gelman_Rubins(chainstart, alpha_start,  Z.safe, Z.fail)
            
            test_chain_data <- alpha_adapt(startvalue=chains$chain1[dim(chains$chain1)[1],],  Z.safe=Z.safe, Z.fail=Z.fail,  alpha_start=alpha_start, transfMCMC=TRUE)
            time2<-Sys.time()
            cat(paste(" Building chain "));flush.console()
            chain1 = Metropolis_Hastings(startvalue=chains$chain1[dim(chains$chain1)[1],], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha)
            chain2 = Metropolis_Hastings(startvalue=chains$chain2[dim(chains$chain2)[1],], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha)
            chain3 = Metropolis_Hastings(startvalue=chains$chain3[dim(chains$chain3)[1],], nb_iter=as.integer(M/3), Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha)
            #gelman.plot(combinedchains)
            chain = mcmc(rbind( as.matrix(chain1$chain), as.matrix(chain2$chain), as.matrix(chain3$chain) ))
            accept_rate <- (chain1$accept_rate + chain2$accept_rate + chain3$accept_rate)/3
            chain=window(chain, thin=2) 
            
          }else{
            
            M <- 2*M                        #car on divise par 2 apr?s
            chainstart <- rej_sampl_nondom(indep_chain, start_point = 1, Z.safe, Z.fail)$val
            #browser() #On choisit un param?tre d'?chelle proche de 1/4 par dichotomie avec alpha_adapt
            try(alpha_start <-  test_chain_data$alpha)
            test_chain_data <- alpha_adapt(startvalue=chainstart,  Z.safe=Z.safe, Z.fail=Z.fail,  alpha_start=alpha_start, transfMCMC=TRUE)
            chainstart <- test_chain_data$chain[dim(test_chain_data$chain)[1], ]
            chain_data <- Metropolis_Hastings(startvalue=chainstart, nb_iter=M, Z.safe=Z.safe, Z.fail=Z.fail, sigma=test_chain_data$alpha )
            chain <- mcmc(rbind( as.matrix(test_chain_data$chain), as.matrix(chain_data$chain)))
            accept_rate <- chain_data$accept_rate
            time2<-Sys.time()
            chain=window(chain, start=as.integer(M/2), thin=1) #ON FAIT LE BURN-IN
          }
          
          ACT=as.integer(dim(chain)[1]/effectiveSize(chain))+1
          ACT <- as.integer(mean(ACT))   #ACT <- max (ACT)   #ou quantile(ACT)[4] si on veut le 3e quartile slmnt (en grande dim); prendre[3] pour la mediane
          # Il faut ensuite multiplier M par l'ACT si on utiliqe le thm ergodique. (on utilise la moyenne des 3 estimatuers caclul?s pour chque chaine, ce qui divise la variance par 3 et l'?cart-type par sqrt(3)).
          M_MC <- 50 #100
          cat(paste("  ACT =", ACT, "\n", "accept_rate=", accept_rate, "\n",  " size of chain  (old_M)", dim(chain)[1], " size needed  (new_M)", M_MC *ACT, "alpha=", test_chain_data$alpha, "\n"  ));flush.console()
          if (ACT>=1000) { ACT <- 1000 ;  cat(paste(" size fixed  (new_M)", M_MC *ACT, "\n"  ));flush.console()  }
          if (M_MC *ACT<=dim(chain)[1]) {
            chain <- mcmc(chain[1:(M_MC *ACT), ])
          }else{
            chain_supp = Metropolis_Hastings(startvalue=chain[dim(chain)[1], ], nb_iter=M_MC*ACT-dim(chain)[1], Z.safe=Z.safe, Z.fail=Z.fail)$chain
            chain=mcmc(rbind( as.matrix(chain), as.matrix(chain_supp)))
          }
          
          indep_chain=window(chain, thin=ACT) 
          M <- dim(chain)[1]
          ttsf <- rep(1,M)
          CPT_CHAIN <- 1
          beta <- c(beta, (old_volnondom_chain - volnondom_chain)/(M_MC*volnondom_chain) )
          first_loop <- FALSE
          time3 <- Sys.time()
          cat(paste(" (Burn-in) and  chainstart + scale param for MC total time  : ", as.double(difftime(time2, time1, tz,units = "secs")), "s", "\n"  ));flush.console()
          cat(paste(" building MC time : ", as.double(difftime(time3, time2, tz,units = "secs")), "s", "\n"  ));flush.console()
        
      }
   

      
    } #fin du while cp < N.appels
    
    
    
    print("On est dans la methode MRM_MCMC ")
    #RR1 <- cbind(Um, UM,  MLE)
    RR1 <- cbind(Um, UM)
    #browser()
    RR <- list(tab=RR1, List.Z.fail=List.Z.fail, List.Z.safe=List.Z.safe, ZF=ZF, ZS=ZS, beta=beta)
    #return(RR) 
    return(RR) 
    
  }  
  
  
  
  
  
  #Budget de point pour l'initialisation dichotomique
  N.DICHO   <- floor(1 + (ordre.p+2)*log(10)/log(2))
  
  if(Method == "MRM"){
    RESULT <- Choix.method.1(N.appels, G, N.dicho = N.DICHO)
  }
  
  if(Method == "MRM_IS"){
    RESULT <- Choix.method.2(N.appels, G, choix.classif, N.dicho = N.DICHO)
  }
  
  if(Method == "MRM_IS_MCMC"){
    RESULT <- Choix.method.3(N.appels, G, choix.classif, N.dicho = N.DICHO, mcmc_period = mcmc_period )
  }
  
  if(Method == "MRM_MCMC"){
    RESULT <- Choix.method.4(N.appels, G, choix.classif, N.dicho = N.DICHO, mcmc_period = mcmc_period )
    #RESULT <- Choix.method.5(N.appels, G, choix.classif, N.dicho = N.DICHO, mcmc_period = mcmc_period )
  }

  if(Method == "MRM_MCMC2"){
    RESULT <- Choix.method.5(N.appels, G, choix.classif, N.dicho = N.DICHO, alpha_start = alpha_start )
  }
  
  if(Method == "MC"){
    RESULT <- Method.MC(N.appels)
  }
  
  if(Method == "MC_monotone"){
    RESULT <- Method.MC.monotone(N.appels)
  }
  
  if(Method == "M-Subset"){
    RESULT <-   M.Subset(N.appels, G, N.dicho = N.DICHO)
  }
  
  return(RESULT)
  
}


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#						Fin de MRM	
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------




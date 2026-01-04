#-----------------------------------------------------------------------------------------------------------------#
#
#
#
#                                                  Algorithme MRM :
#
#
#
#-----------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------------------#
# Input :											    	    
#-----------------------------------------------------------------------------------------------------------------#
# f : Fonction de défaillance.								    
# ndim : Dimension de l'espace
# choix.loi : liste contenant le nom des lois d'entrées et leur(s) paramètre(s).  
# dir.monot : vecteur appartenant à {-1,1}^ndim. dirmonot[i] = -1 : code décroissant ; 1 : code croissant en X[i]
# N.appels : Budget d'appel au code donné par l'utilisateur
# Method = {1, 2, 3,...} choix de la méthode
# Code = {"R", "C"}, choix du langage de programmation pour certaine fonction
# ordre.p : Ordre de la proba cherchée : si la proba est de l'ordre de 1e-k, ordre.p = k.
#           On réservera alors floor((ordre.p+2)*log(10)/log(2)) appels pour la recherche par dichotomie
# BUDGET : Nombre d'appels au code autorisé par l'utilisateur
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#
# Output :	Renvoi une matrice de taille BUDGETx6 : 
#-----------------------------------------------------------------------------------------------------------------#
#  um = bornes inférieures 
#  uM = bornes supérieures 
#  MLE = Estimateur du maximum de vraissemblance
#  u.m.Fisher = Bornes inférieures de l'intervalle de confiance lié au MLE
#  u.M.Fisher = Bornes supérieures de l'intervalle de confiance lié au MLE
#  CV = coefficient de variantion de l'estimateur MLE
#-----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------#

dyn.load("/home/C23072/Documents/MRM/Code_C/Domination/IS_DOMINANT.so")
#dyn.load("D:\\MRM\\MRM_Icossar\\IS_DOMINANT.dll")


MRM <- function(f, ndim, choix.loi, dir.monot, N.appels, Method = 1, code, ordre.p, MAXIMIN ){

#-----------------------------------------------------------------------------------------------------------------#
#                         Création de la liste contenant le nom des lois d'entrées et			
#                                    des paramèters de chacune des lois
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

#-----------------------------------------------------------------------------------------------------------------#
#	                           Transformation dans l'espace uniforme	  
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
#                Calcul de l'intersection entre la diagonale et l'etat limite de defaillance
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
    cp 	    <- 1 						# compteur d'appel à G
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

    while(cp < N.dicho){         #N.dicho définie en début de programme
 
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

  is.dominant.C <- function(x, y, ndim, set){
    L.x <- 0

    if(is.null(dim(x)) ){
	L.x <- 1;
    }else{
      L.x <- dim(x)[1]
    }

 #   set <- ifelse(set == 1, 2, -1)

    dens <- .C("is_dominant", 
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

# Cas où x est un vecteur
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
#		                Fonction qui crée le jeu de donnée
#-----------------------------------------------------------------------------------------------------------------#
  Create.data.C <- function(ORDRE.P){
    K <- 10^(ORDRE.P + 2)

    .C("Sauvegarde", 
       as.integer(K),
       as.integer(ndim))
  }
#-----------------------------------------------------------------------------------------------------------------#
#		              Fonction qui retourne les points
#                   qui participe au volume d'un ensemble S de points
#-----------------------------------------------------------------------------------------------------------------#
  VOLVOL <- function(X.MC, S, set){ 

   if(set == 1){S <- 1 - S}
 
    if(is.null(dim(S))){
      if(set == 1){
        #return(prod(1 - S))
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

#-----------------------------------------------------------------------------------------------------------------#
#		                Fonction qui retourne la valeur du volume
#		                     déterminés par les points de S
#			                par la méthode de Monte Carlo
#-----------------------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------------------------------#
#				                    Code en R
#-----------------------------------------------------------------------------------------------------------------#
    Volume.MC.R <- function(S, set){
	
	if(ordre.p <= 3){
        K <- 10^(ordre.p + 2)
      }else{
        K <- 1e6
      }

      Y 	 <- NULL
	Y.res  <- NULL

      if(set == 1){S <- 1-S}
	if(is.null(dim(S))){       
	  return(prod(S))
	}     

      # On construit Y, un échantillon uniforme sur [0,1]^DIM

      for(d in 1:ndim){
        Y.temp <- runif(K)
        Y 	   <- cbind(Y, Y.temp)
      }

      while(!is.null(dim(Y))){
	  R <- is.dominant(S, Y[1,], ndim, 2)
	  if( sum(R) != 0 ){
		R.bis  <- is.dominant(Y, Y[1, ], ndim, 1)
		Y.temp <- Y[which(R.bis == TRUE), ]
		Y.res  <- rbind(Y.res, Y.temp)
		Y 	 <- Y[which(R.bis == FALSE), ]
	  }else{
	    R.bis   <- is.dominant(Y, Y[1,], ndim, 2)
	    Y 	<- Y[which(R.bis == FALSE), ]
	  }
	  if( (is.null(dim(Y)))| (length(Y) == 0)){
		test <- is.dominant(Y, S, ndim, 1)
		if(sum(test) != 0){Y.res <- rbind(Y.res, Y)}
		P    <- dim(Y.res)[1]/K
	      return(P)
        }
	}
      P <- dim(Y.res)[1]/K
	return(P)
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
#				                     Extrait la frontière d'un ensemble
#-----------------------------------------------------------------------------------------------------------------

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
#
#
# 		         				Méthode de Monte Carlo
#
#
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
################################################################
#             SUSET MONOTONE
################################################################
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
#-----------------------------------------------------------------------------------------------------------------

    # p.k = (p.k^-, p.k^+)
    # p : variables
    # signature[1] = 1 si H(Y[k]) <= 0, 0 sinon
    # On calcul ici  : (- l'estimateur du maximum de vraissemblance), car on utilise
    #                    la fonction "optimize" de R, qui
    #                    renvoi l'argmin de la fonction

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

#    log.likehood.with.old <- function(p , p.k, signature){
#      pp <- 1
#      dd <- dim(signature)[1]
#      for(i in 1:dd){
#        gamma <- (p - p.k[i,1])/(p.k[1:i,2] - p.k[1:i,1])
#        pp <- pp*prod( ((gamma)^signature[1:i,i])*((1 - gamma)^(1 - signature[1:i,i])))
#      }
#      return(pp)
#    }

#maximum de vraisemblance pour la classification
            log.likehood.C <- function(signature, p){
            u     <- signature*log(p/(1 - p)) + log(1 - p)
            return(sum(u))
            }
#-----------------------------------------------------------------------------------------------------------------
#		                         Fin de la construction du maximum de vraissemblance
#-----------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------
#		                            Convertit un entier en binaire
#-----------------------------------------------------------------------------------------------------------------
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
#		                            Fin de la conversion en binaire
#-----------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------
#		                             On tire 1 point uniformément autour de x
#-----------------------------------------------------------------------------------------------------------------
    
    SIM <- function(x, W){

#      W <- 1:(2^(ndim) - 1)
      B <- 0
      # On sépare l'espace autour de x en 2^d - 1 zones, et l'on calcul le volume de chaque zone.
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
      U <- runif(1, 0, max(B))

      #On choisit aléatoirement dans quelle "zone" autour de x on va tirer des points
      pos <- ifelse(U < B[1], 1, which.max(B[B <= U]) + 1)

      Z <- as.binary(W[pos])
      A <- 0
      
      #On tire uniformément dans la zone choisit aléatoirement
      for(i in 1:ndim){
        A[i] = ifelse(Z[i] == 0, runif(1, x[i], 1), runif(1, 0, x[i]))  
      }
      return(A)      
    }
#-----------------------------------------------------------------------------------------------------------------
#		Fin du tirage de 1 point uniformément autour de x
#-----------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------
#                   Simulation de CP points dans l'espace non dominé
#-----------------------------------------------------------------------------------------------------------------

    Sim.non.dominated.space <- function(CP, Z.safe, Z.fail, W){
      CP1 <- 0;
      Y   <- NULL
      Y.temp  <- apply(1 - Z.safe, MARGIN = 1, prod)
      Y.temp1 <- Z.safe[which.max(Y.temp),]
      while(CP1 < CP){
        Y.temp2 <- SIM(Y.temp1, W)
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
# 		     				   PARTIE 1 : Sans procédure de classification
#
#
#
#-----------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------
#
#
#
# 						            METHODE 1.1 : Monte Carlo progressif
#
#
#
#-----------------------------------------------------------------------------------------------------------------

  Choix.method.1.1 <- function(N.appels, H, N.dicho){

    # V contient le premier tri de points avec la méthode de la dichotomie.
    # Z.fail : ensemble des points du domaines de défaillance
    # Z.safe : ensemble des points du domaines de sûreté
    # cp : compteur d'appel à H
    # um : vecteur contenant les bornes exactes inférieures
    # uM : vecteur contenant les bornes exactes supérieures
    # MLE :Estimateur du maximum de vraissemblance
    # I.Fisher : information de Fisher associé au MLE
    # CV : vecteur contenant les coefficients de variations
    # SIGN : Signature des points testés

    #Création du tirage si ndim > 2

    if(ndim > 2){
     # zz <- file("C:\\Users\\C23072\\Documents\\MRM\\Code_C\\Tirage\\echantillon.bin", "wb") 
      UU <- runif(ndim*10^(ordre.p + 2))
     # writeBin(UU, zz, size = 8)
      UU <- matrix(UU, ncol=ndim, byrow = TRUE)
      #Create.data.C(ordre.p)
      D.UU <- dim(UU)[1]
    }

    if(ndim == 2){UU <- 0; D.UU <- 0}
    
    V <- list()
    V <- Intersect(ndim, H, N.dicho)    

    list.set <- 0
    for(i in 1:length(V)){
      list.set[i] <- V[[i]][[2]]
    }

    u.dep   <- list()
    u.other <- list()

    u.dep[[1]] <- V[[max(which(list.set == -1))]][[1]]
    u.dep[[2]] <- V[[max(which(list.set == -1))]][[2]]

    u.other[[1]] <- V[[max(which(list.set == 1))]][[1]]
    u.other[[2]] <- V[[max(which(list.set == 1))]][[2]]

    Z.fail <- t(as.matrix(u.dep[[1]]))
    Z.safe <- t(as.matrix(u.other[[1]]))

    cp     <- length(V)  
    
    um <- 0
    uM <- 1

    Um <- 0
    UM <- 1
    eps <- 1e-7
    alpha <- 0.025
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

#    if(V[[1]][[2]] == 1){
#      uM[1]   <- 1 - prod(1 - V[[1]][[1]])
#      UM[1]   <- 1 - prod(1 - V[[1]][[1]])
#      SIGN[1] <- 0
#      ZS <- V[[1]][[1]]
#      ZF <- NULL
#    }else{
#      um[1]   <- prod(V[[1]][[1]])
#      Um[1]   <- prod(V[[1]][[1]])
#      SIGN[1] <- 1
#      ZF <- V[[1]][[1]]
#      ZS <- NULL
#    }


    # N.set : Nombre de points théoriques qui tomberais dans 
    # l'ensemble "set" pour trouver la proba um[1] ou uM[1]




#    for(i in 1:cp){
#      if(V[[i]][[2]] == 1){
#        UM[i]   <- 1 - prod(1 - V[[i]][[1]])

#        if(i == 1){
#          Um[i] = 0
#        }else{
#          Um[i]   <- Um[i - 1]
#        }

#        SIGN[i] <- 0

#        ZS <- rbind(ZS, V[[i]][[1]])
#      }else{
#        Um[i]   <- prod(V[[i]][[1]])
#        if(i == 1){
#          UM[i] = 1
#        }else{
#          UM[i]   <- UM[i - 1]
#        }
        
#        SIGN[i] <- 1
       
#        ZF <- rbind(ZF, V[[i]][[1]])
#      }



#      if(MAXIMIN == "MLE"){
#        MLE.test <- optimize(f = log.likehood,
#                             interval = c(Um[i],UM[i]),
#                             signature = SIGN,
#                             p.k = cbind(Um, UM)
#                             )   
#        MLE[i] <- as.numeric(MLE.test[1])
#        VAR <- sum( 1/((MLE[i] - um)*(uM- MLE[i])))
#        bn <- 1/VAR
#        an <- eps*VAR^(5/2)/abs( sum( 1/((MLE[i] + eps - Um[i])*(UM[i] - MLE[i] - eps))) - sum( 1/((MLE[i] - Um[i])*(UM[i]- MLE[i])))  )
#        ICinf[i] <- MLE[i] - qnorm(1 - alpha)/sqrt(VAR + alpha/an)
#        ICsup[i] <- MLE[i] + qnorm(1 - alpha)/sqrt(VAR - alpha/an) 
#        CV[i] <- 100/(sqrt(VAR)*MLE[i])
#      }
#else{
#        p.bar[i] = (Um[i] + UM[i])/2
#      }
      #print("On test le nouvel estimateur!");browser()
#       p.bar[i] = (Um[i] + UM[i])/2
#    }



#    j <- cp + 1
um <- prod(V[[cp]][[1]])
uM   <- 1 - prod(1 - V[[cp]][[1]])

j = 1

SIGN.2.old <- 0

v = 0
M.f = 0
var.p.bar.theo = 0
var.p.hat.theo = 0
omega.theo = 0
p.hat.theo = 0
p.hat.bis = 0

VV <- 0
VAR.2 <- 0
ICinf.2 <- 0
ICsup.2 <- 0
p.bar.1 = 0
p.bar.2 = 0
p.bar.3 = 0

M.f = 0
var.p.bar = 0
omega = 0
p.hat = 0
p.hat.1 = 0
var.p.hat = 0
IC.inf.max <- 0
IC.sup.max <- 0
var <- 0
 W <- 1:(2^(ndim) - 1)
    if(MAXIMIN == "classification"){
#      library("monmlp")
source("package_R/monmlp/R/gam.style.R")
source("package_R/monmlp/R/linear.prime.R")
source("package_R/monmlp/R/linear.R")
source("package_R/monmlp/R/logistic.prime.R")
source("package_R/monmlp/R/logistic.R")
source("package_R/monmlp/R/monmlp.cost.R")
source("package_R/monmlp/R/monmlp.fit.R")
source("package_R/monmlp/R/monmlp.initialize.R")
source("package_R/monmlp/R/monmlp.nlm.R")
source("package_R/monmlp/R/monmlp.predict.R")
source("package_R/monmlp/R/monmlp.reshape.R")
source("package_R/monmlp/R/tansig.prime.R")
source("package_R/monmlp/R/tansig.R")
    }
    if(MAXIMIN != "MLE"){
      library("stats")
    }

    SIMU <- NULL

    while(cp < N.appels){
      # CP : Nombre de points que l'on veut avoir dans l'espace non dominé
      # Y : ensembles des points dans l'espace non dominé

#      print(paste("compteur =",cp));flush.console()
      cat(paste("  compteur =",cp, "--"));flush.console()


      # temp.vol : on récupère pour chaque point de Y, sa contribution au volume de
      #            de Z.safe et Z.fail, en récupérant le min des deux.


      if(MAXIMIN != "MLE"){

        temp.vol <- 0;

        CCPP <- 40
        Y <- Sim.non.dominated.space(CCPP, Z.safe, Z.fail, W)

        if(MAXIMIN == "classification"){

          X = rbind(Z.safe, Z.fail)
          indic.s <- ifelse(is.null(dim(Z.safe)), 1, dim(Z.safe)[1])
          indic.f <- ifelse(is.null(dim(Z.fail)), 1, dim(Z.fail)[1])
          T <- cbind(  c(rep(1, indic.s), rep(0, indic.f)), c(rep(0, indic.s), rep(1, indic.f)))

#         W.MON <- monmlp.fit(x = X, y = T, hidden1 = 3, monotone = c(1,2),
#                               n.ensemble = 2, bag = TRUE, silent = TRUE)

          #W.MON <- monmlp.fit(x = X, y = T, hidden1 = 3, monotone = 2, n.ensemble = 2, bag = TRUE, silent = TRUE)

          W.MON <- monmlp.fit(x = X, y = T, hidden1 = 3, monotone = 2, n.ensemble = 2, bag = FALSE, silent = TRUE)

        }

        temp.vol <- apply(Y, MARGIN = 1,
			      function(u){
                                if(MAXIMIN == "sort"){
                                  v.s <- is.dominant(Y, u, ndim, 2)
                                  v.f <- is.dominant(Y, u, ndim, 1)
                                  vol.s <- sum(v.s)
                                  vol.f <- sum(v.f)
                                  vol <- min( c(vol.s, vol.f) )
                                  return(vol)
                                }
                                if(MAXIMIN == "vol"){
                                  ttt     <- 0
                                  test.s  <- rbind(Z.safe, u)
                                  test.f  <- rbind(Z.fail, u)
                                  VS      <- VOLVOL(UU, test.s, 1)
                                  VF      <- VOLVOL(UU, test.f, 2)
                                  vol.test.s <- UM[j-1] - VS
                                  vol.test.f <- VF - Um[j-1]
                                  vol.test <- c(vol.test.s, vol.test.f)
                                  ttt      <- min(vol.test)
                                  return(ttt) 
                               }
                               if(MAXIMIN == "classification"){
                                 ttt     <- 0
                                 test.s  <- rbind(Z.safe, u)
                                 VS      <- VOLVOL(UU, test.s, 1)
                                 #p.mon   <- monmlp.predict(x = matrix(u, nrow = 1), weights = W.MON)
                                 bound <- ifelse(j == 1, uM, UM[j-1])
                                 #vol.test.s <- (bound - VS)*p.mon[1]
                                 vol.test.s <- (bound - VS)
                                 return(vol.test.s) 
                               }
                               if(MAXIMIN == "essai"){
                                 ttt     <- 0
                                 test.s  <- rbind(Z.safe, u)
                                 test.f  <- rbind(Z.fail, u)
                                 VS      <- VOLVOL(UU, test.s, 1)
                                 VF      <- VOLVOL(UU, test.f, 2)
                                 bound <- ifelse(j == 1, uM - um, UM[j-1] - Um[j-1])
                                 vol.test.s <- bound/(VS - VF)
                                 return(vol.test.s) 
                               }
			      })
          if(MAXIMIN == "classification"){
            temp.vol <- temp.vol*apply(Y, MARGIN = 1, function(z){monmlp.predict(x = matrix(z, nrow = 1), weights = W.MON)[1]})
          }

          u <- Y[which.max(temp.vol), ]
#          uu <- Y[which.max(temp.vol), ]
#          b <-  ifelse(j == 1, um, Um[j-1])
#          B <-  ifelse(j == 1, uM, UM[j-1])
##          v[j] <- (B - b)/ndim
#          v[j] <- (B - b)

 #         CCPP <- 50

#          SIM.NORM <- function(Compt, Z.safe, Z.fail, uu, var){
#            res <- matrix(ncol= ndim, nrow = Compt)
#            for(i in 1:Compt){
#              vv <- -1
#              f <- -1
#              s <- -1
#              while( !(  (sum(s) == 0) & (sum(f) == 0) & (sum(vv <= 1) == ndim) & (sum(vv>= 0) == ndim )) ){
#                vv <- uu + rnorm(ndim, 0 , sqrt(var) )
#                s <- is.dominant(Z.safe, vv, ndim, 1)
#                f <- is.dominant(Z.fail, vv, ndim, 2)
#              }
#               res[i, ] <- vv
#            }
#            return(res)
#          }

#          u <- SIM.NORM(1, Z.safe, Z.fail, uu, v[j])

      }else{
        CCPP <- 1
        u <- Sim.non.dominated.space(CCPP, Z.safe, Z.fail, W)
      }#fin du else
      # On tri les points de Y : 
      SIMU <- rbind(SIMU, u)
      Z.safe.old <- 0
      Z.fail.old <- 0
      # On récupère la contribution maximale parmi toute les contributions

      # On regarde dans quel ensemble appartient u
      # On rajoute u à l'ensemble auquel il appartient. On garde sa signature
      # et on calcule le volume du nouvel ensemble construit
      H.u = H(u)
      if(H.u > 0){
        Z.safe.old <- rbind(u, Z.safe)
        ss     <- is.dominant(Z.safe, u, ndim, 2)
        Z.safe <- Z.safe[which(ss == FALSE), ]
        Z.safe <- rbind(u, Z.safe)
      # tt     <- is.dominant(Y, u, ndim, 2)

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

        CC <- ifelse(cp == N.appels, 1, 0)
        SIGN[j] <- 0

      }else{
        Z.fail.old <- rbind(u, Z.fail)

        ff     <- is.dominant(Z.fail, u, ndim, 1)
        Z.fail <- Z.fail[which(ff == FALSE), ]
        Z.fail <- rbind(u, Z.fail)
       # tt     <- is.dominant(Y, u, ndim, 1)

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


        CC <- ifelse(cp == N.appels, 1, 0)
        SIGN[j] <- 1
      }
      cp <- cp + 1

      ZS[[j]] <- Z.safe
      ZF[[j]] <- Z.fail

      if(MAXIMIN == "MLE"){
        MLE.test <- optimize(f = log.likehood,
                               interval = c(Um[j],UM[j]),
                               maximum = TRUE,
                               tol = 10^(-3*ordre.p),
                               signature = SIGN,
                               p.k = cbind(Um, UM)
                               )   
        MLE[j] <- as.numeric(MLE.test[1])

        VAR <- sum( 1/((MLE - Um)*(UM- MLE)))
        bn  <- 1/VAR
        an  <- eps*VAR^(5/2)/abs( sum( 1/((MLE + eps - Um)*(UM - MLE - eps))) - sum( 1/((MLE - Um)*(UM- MLE)))  )

        ICinf[j]  <- MLE[j] - qnorm(1 - alpha)/sqrt(VAR + alpha/an)
        ICsup[j]  <- MLE[j] + qnorm(1 - alpha)/sqrt(VAR - alpha/an) 
        CV.MLE[j] <- 100/(sqrt(VAR)*MLE[j])   


        SIGN.2 <- matrix(0, ncol = j, nrow = j)
        diag(SIGN.2) <- SIGN

        if(j > 1){
          SIGN.2[1:(j-1),1:(j-1)] <- SIGN.2.old 
          for(jj in 1:(j-1)){
            if( SIGN.2[jj, j-1] == 1){
              if( SIGN[j] == 1){
                if( sum( SIMU[jj,  ] <= SIMU[j, ] ) != ndim ){
                  SIGN.2[jj, j] <- 1
                }
              }else{
                  SIGN.2[jj, j] <- 1
              }
            }
          }#Fin de la boucle for
        }

#browser();
        SIGN.2.old <- SIGN.2

#browser();
#        MLE.test.2 <- optimize(f = log.likehood.with.old, interval = c(Um[j],UM[j]),maximum = TRUE,signature = SIGN.2,  p.k = cbind(Um, UM) )
        MLE.test.2 <- optimize(f = log.likehood.with.old, interval = c(Um[j],UM[j]),maximum = TRUE,tol = 10^(-3*ordre.p) ,signature = SIGN.2,  p.k = cbind(Um, UM) )      


        p.bar[j] <- as.numeric(MLE.test.2[1])

        VAR.2 <- VAR.2 + sum(1/( (UM - Um - p.bar[j] + UM[j])*(p.bar[j] - Um[j] )  )  )
        ICinf.2[j]  <- p.bar[j] - qnorm(1 - alpha)/sqrt(VAR.2)
        ICsup.2[j]  <- p.bar[j] + qnorm(1 - alpha)/sqrt(VAR.2) 


#yp <- 0
#yp1 <- 0
#MM <- 1000
#p <- seq(from = Um[j], to = UM[j],length=MM)
#for(iii in 1:MM){
#  yp[iii] <- log.likehood.with.old(p[iii], cbind(Um, UM),SIGN.2 )
#  yp1[iii] <- log.likehood(p[iii], cbind(Um, UM),SIGN )
#}
#par(mfrow=c(2,2))
#plot(p, yp, type='l')
#lines(rep(p.bar[j], 1e3), seq(from = min(yp[-c(1,MM) ]),to=max(yp[-c(1,MM)]),length=1e3 ),col=1 )
#lines(rep(  p[ which.max( yp[-c(1,MM)] ) ] , 1e3), seq(from = min(yp[-c(1,MM) ]),to=max(yp[-c(1,MM)]),length=1e3 ),col=1,lwd=2,lty=2 )
#plot(p, log(yp1),col="red", type='l')
#lines(rep(MLE[j], 1e3), seq(from = min(log(yp1[-c(1,MM) ])),to=max(log(yp1[-c(1,MM)])),length=1e3 ),col=2 )
#lines(rep(  p[ which.max(yp[-c(1,MM)] ) ], 1e3) , seq(from = min(log(yp1[-c(1,MM) ])),to=max(log(yp1[-c(1,MM)])),length=1e3 ),col=2 ,lwd=2,lty=2 )
#plot(1:j,  cumsum(( true.value -  p.bar)^2),col=1,type='l' ,lwd=2,log="y")
#lines(1:j, cumsum(( true.value - MLE)^2), col=2 ,lwd=2)

#print(c( p[which.max(yp[-c(1,500)] ),max(yp1[-c(1,500)] ))  )

#print(p.bar[j]);print(SIGN[j]);
#browser();
#        temp.estim <- 0
#        for(kk in 1:j){
#          temp.estim.1 <- 0
#          cp.temp <- 0
#          Bjk <- ((MLE[j] - Um[kk])*(UM[1:kk] - Um[1:kk] + Um[kk] - MLE[j] ))^(-1)/sum( ((MLE[j] - Um[kk])*(UM[1:kk] - Um[1:kk] + Um[kk] - MLE[j] ))^(-1))
#          for(jj in 1:kk){
#            ts <- is.dominant(ZS[[jj]],SIMU[jj,] , ndim, 1)
#            tf <- is.dominant(ZF[[jj]], SIMU[jj,], ndim, 2)
#            if( (sum(ts) == 0 ) & ( sum(tf) == 0) ){
##              temp.estim.1 <- temp.estim.1 + Um[jj] + (UM[jj] - Um[jj])*SIGN[jj]
##              temp.estim.1 <- temp.estim.1 + Um[jj] + ((UM[jj] - Um[jj])/(UM[kk] - Um[kk]))*( Um[kk] + (UM[kk] - Um[kk])*SIGN[jj])
#              temp.estim.1[jj] <- Um[kk] + (UM[jj] - Um[jj])*SIGN[jj]
##              cp.temp <- cp.temp + 1
#            }else{
#              temp.estim.1[jj] <- Um[kk]
#            }
#          }
##          temp.estim[kk] <- ifelse(cp.temp == 0, 0 ,temp.estim.1/cp.temp)
#          temp.estim[kk] <- sum(Bjk*temp.estim.1)
#        }
#        p.bar.1[j] <- mean(temp.estim)
#        p.bar.2[j] <- mean(Um + (UM - Um)*SIGN)
#        OMEGA <- ((MLE - Um)*(UM - MLE[j]))^(-1)/sum( ( (MLE - Um)*(UM - MLE[j]))^(-1))
#        p.bar.3[j] <- sum(OMEGA* (Um + (UM - Um)*SIGN))

      }else{

#        M  <- NULL
#        M1 <- NULL
#        RR <- Sim.non.dominated.space(300, Z.safe, Z.fail, W)
#        M1 <- c(M1, apply(RR, MARGIN = 1, function(rr){exp(- sum(abs(rr - uu)^2)/(2*v[j]))}))
#        MMM1 <- mean(M1)*(B-b)

##        NN <- 100
#        #S <- SIM.NORM(NN, Z.safe, Z.fail, uu, v[j])
#        #S1 <- apply(S, MARGIN = 1, function(ss){exp(sum(abs(ss - uu)^2)/(2*v[j]))})
#        p.bar[j] <- (Um[j] + UM[j])/2  + (2*SIGN[j] - 1)*exp( sum(abs(u - uu)^2)/(2*v[j]))*MMM1/2
#        if(  !( (p.bar[j] >= Um[j])&(p.bar[j])<=UM[j]) ){
#          if(p.bar[j] < Um[j]){
#            p.bar[j] <- Um[j]
#          }
#          if(p.bar[j] > UM[j]){
#            p.bar[j] <- UM[j]
#          }
#        }
#        if(SIGN[j] == 0){
#          p.bar[j] <- Um[j]
#        }else{
#          M  <- NULL
#          M1 <- NULL
#          RR <- Sim.non.dominated.space(300, Z.safe, Z.fail)
#          M1 <- c(M1, apply(RR, MARGIN = 1, function(rr){exp(- sum(abs(rr - uu)^2)/(2*v[j]))}))
#          MMM1 <- mean(M1)*(B-b)
#          p.bar[j] <- UM[j] - exp( sum(abs(u - uu)^2)/(2*v[j]))*MMM1/2
#        } 


 #       #var.p.bar[j] <- 0.25*mean(S1^2)*MMM1^2
 #       var.p.bar[j] <-  (UM[j] - Um[j])^2

#        omega <- (j/var.p.bar)/sum(1/var.p.bar)

#        p.hat[j] <- mean(omega*p.bar)

#        var.p.hat[j] <- sum(1/var.p.bar)^(-1)
#        CV.max[j] <- 100*sqrt(var.p.hat[j])/p.hat[j]

#        IC.inf.max[j] <- p.hat[j] - qnorm(1-alpha/2)*sqrt(var.p.hat[j])
#        IC.sup.max[j] <- p.hat[j] + qnorm(1-alpha/2)*sqrt(var.p.hat[j])
##print((Um[j] + UM[j])/2);flush.console();
##print(p.bar[j]);flush.console();
#print(p.hat[j]);flush.console();
##print(var.p.bar[j]);flush.console();
#print(var.p.hat[j]);flush.console();
#if(j>=10){browser();}
#browser();
        }#fin du if maximin == "sort"

        j <- j + 1 

    } #fin du while cp < N.appels
# browser();
    if(MAXIMIN  == "MLE"){
print("On est dans le MLE")
#      RR <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE)
      #RR <- list(Um, UM, MLE, ICinf, ICsup, CV.MLE, Z.fail, Z.safe)
#      RR <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE, p.bar.1, p.bar.2, p.bar.3 )
      RR <- cbind(Um, UM, MLE, ICinf, ICsup, CV.MLE, p.bar, ICinf.2, ICsup.2)
#lll = length(Um);plot(1:lll,Um,type='l',ylim=c(0,0.2));lines(1:lll,UM);lines(1:lll,p.bar,col="orange");lines(1:lll,p.hat,col="blue");
#lines(1:lll,rep(true.value,lll),col="purple");lines(1:lll,MLE,col="red")
#browser();
       #RR <- list(Um, UM, MLE, p.hat, ICinf, ICsup, CV.MLE, CV.max, Z.fail, Z.safe)
    }else{
print("On est dans le sort")
     # RR <- cbind(Um, UM, p.hat, IC.inf.max, IC.sup.max, CV.max)
 RR <- cbind(Um, UM)
#plot(1:(j-1),omega, col="red",type='l');lines(1:(j-1),omega.theo, col="blue")
#lll = length(Um);lines(1:lll,Um,type='l',ylim=c(0,0.2),col="orange");lines(1:lll,UM,col="orange");lines(1:lll,p.bar,col="green");
#lines(1:lll,p.hat,col="blue");
#lines(1:lll,rep(true.value,lll),col="purple")

#browser();
#      RR <- list(Um, UM, p.hat, IC.inf.max, IC.sup.max, CV.max)
    } 

    return(RR)
  }

#MLE1 = 0
#for(i in 1:lll){
#a <- optimize(f = log.likehood,
#               interval = c(Um[i],UM[i]),
#                maximum = TRUE,
#                signature = SIGN[1:i],
#                p.k = cbind(Um[1:i], UM[1:i]))   
#    MLE1[i] <- as.numeric(a[1])}
#lines(1:lll,MLE1,col="brown")

#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					               FIN	METHODE 1.1
#
#
#
#-----------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
# 					            	FIN DE LA PARTIE 2
#
#
#
#-----------------------------------------------------------------------------------------------------------------  
  #Budget de point pour l'initialisation dichotomique
  N.DICHO   <- floor(1 + (ordre.p+2)*log(10)/log(2))

  if(Method == "MRM"){
    RESULT <- Choix.method.1.1(N.appels, G, N.dicho = N.DICHO)
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



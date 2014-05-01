cat("\014")
source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/extraction data.r")
#source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/tests stats.r")
#source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/opti pf.r")
source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/previsions.r")
library(nleqslv)
##On dispose des donnees suivantes:
# X vecteur des rendements
# rf rendement moyen de l'actif sans risque
##Les anticipations nous donnent les elements suivants:
# la matrice P informe sur quels actifs sont concernes par quelle prediction
# le vecteur Q (derive de P) renvoie les valeurs de sur/sous-performance absolue pour chaque actif
# le coefficient "weight on views" thau
# la matrice omega donne les niveaux de confiance sur chaque vue

######################## Etape 0 : fonctions de base ( servent pour les 2 modeles) ##############################################################################

meanreturn <- function(X,p){
  # Renvoie le rendement historique moyen du portefeuille
  #X est le vecteur des rendements des actifs du portefeuille, p le vecteur des poids attribués
  somme = 0
  for (i in 1:length(p)){
    somme = somme + mean(X[,i])*p[i]
  }
  return(somme)
}

volatility <- function(X,p){
  #Renvoie la volatilite historique du portefeuille, calculee comme l'ecart type de la var "rendement du portefeuille"
  # X est le vecteur des rendements et p le vecteur des poids
  Y <- X
  for (i in 1:length(p)){
    Y[,i] <- p[i]*Y[,i]
  }
  return(sqrt(var(rowMeans(Y))))
}

r_av_coeff <- function(X,p,rf){
  #Renvoie le coefficient d'aversion au risque determine par le choix du portefeuille
  # X et p determinent le portefeuille, rf est le rendement de l'actif sans risque sur la periode
  return((meanreturn(X,p)-rf)/(volatility(X,p)))
}

f_Er <- function(X,tau,impexc,P,Q,omega){
  # Renvoie E(R) en fonction du implied excess returns et des previsions
  sigma <- sigma(X)
  return(impexc + tau*sigma%*%t(P)%*%solve(omega+tau*P%*%sigma%*%t(P))%*%(Q - P%*%impexc))
}

sigma <- function(X){
  #Renvoie la matrice de variance historique
  return(var(X))
}
sigma_ant <- function(X,v_ant){
  # Renvoie la matrice de variance anticipee
  n <- length(v_ant)
  cor_hist <-cor(X)
  y <- matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      y[i,j] <- cor_hist[i,j]*v_ant[i]*v_ant[j]
    }
  }
  return(y)
}

coskewness <- function(X){
  # Renvoie la matrice de coskewness pour le portefeuille dont on a les rendements historiques X
  n <- ncol(X)
  T <- nrow(X)
  skew <- array(data = 0, dim = c(n,n*n))
  sous_produit<-function(t){
    #Renvoie la matrice de coskewness pour une date t donnee
    x <- t(as.numeric(X[t,] - t(colMeans(X))))
    #return(m_kronecker(list(t(x),x,x)))
    return(kronecker(kronecker(t(x),x),x))
  }
  for (t in 1:T){
    skew <- skew + sous_produit(t)
  }
  return(skew/T)
}
coskewness_ant <- function(X,v_ant){
  #Renvoie la matrice de coskewness anticipee pour les volatilites anticipees (memes hypoth?ses de stationnarite que pour la matrice de variance anticipee)
  #X : rendements historiques
  # v_ant : vecteur de taille celle du portefeuille renvoyant les volatilites anticipees pour chaque actif
  n <- length(v_ant)
  skew_hist=skewness(X)
  std <- apply(X,2,sd)
  y <- matrix(0,n,n*n)
  for (i in 1:n){
    for (j in 1:(n*n)){
      y[i,j] <- skew_hist[i,j]*v_ant[j%/%n]*v_ant[i]*v_ant(j%%n)/(std[j%/%n]*std[i]*std(j%%n))
    }
  }
  return(y)
}

cokurtosis <- function(X){
  # Renvoie la matrice de cokurtosis pour le portefeuille dont on a les rendements historiques X
  n <- ncol(X)
  T <- nrow(X)
  kur <- array(data = 0, dim = c(n,n^3))
  sous_produit<-function(t){
    #Renvoie la matrice de coskewness pour une date t donnee
    x <- t(as.numeric(X[t,] - t(colMeans(X))))
    return(kronecker(kronecker(kronecker(t(x),x),x),x))
  }
  for (t in 1:T){
    kur <- kur + sous_produit(t)
  }
  return(kur/T)
}
cokurtosis_ant <- function(X,v_ant){
  #Renvoie la matrice de coskewness anticipee pour les volatilites anticipees (memes hypoth?ses de stationnarite que pour la matrice de variance anticipee)
  #X : rendements historiques
  # v_ant : vecteur de taille celle du portefeuille renvoyant les volatilites anticipees pour chaque actif
  decomposition <- function(a,b){
    n <- max(4,floor(log(a,b)) +1) #Pour info, donne le nombre de chiffres dans la decomposition
    y <- seq(1,1,length=n)
    print(y)
    q <- a
    i <- 1
    while(q >0){
      y[i] <- y[i]+q%%b
      q <- q %/% b
      i <- i+1
    }
    print(y)
    return(y)
  }
  n <- length(v_ant)
  kur_hist=kurtosis(X)
  std <- apply(X,2,sd)
  y <- matrix(0,n,n^3)
  for (i in 1:n){
    for (J in 1:((n^3)-1)){
      decomp <- decomposition(J,n)
      l <- decomp[1]
      j <- decomp[2]
      k <- decomp[3]
      y[i,j+1] <- (kur_hist[i,j+1]*v_ant[i]*v_ant[j]*v_ant[k]*v_ant[l])/(kur_hist[i,j+1]*std[i]*std[j]*std[k]*std[l])  
    }
    y[i,1] <- kur_hist[i,1]*std[i]*std[1]*std[1]*std[1]
  }
  return(y)
}

######################## Etape 1 : modele de B-L simple et amelioration ##############################################################################

expected_simple <- function(X,p,rf){
  #Renvoie le vecteur "implied excess returns" qui represente les rendements optimaux 
  #etant donnes le choix d'allocation et la var historique
  #NB: si p est historiquement optimal, la fonction renvoie meanreturn(X,p)
  #X, p et rf comme d'habitude, sigma matrice de covariance
  return(r-av_coeff(X,p,rf)%*%sigma(X)%*%p)
}

f_wBL<-function(X,lambda,Er){
  # Renvoie l'allocation optimale (Black-Litterman basique)
  sigma <- sigma(X)
  w=(1/as.numeric(lambda))*solve(sigma)%*%Er
  return(w)
}

f_wBL_somme1<-function(X,lambda,Er){
  # Renvoie l'allocation optimale (B-L avec somme des poids valant 1)
  m=length(Er)
  sigma <- sigma(X)
  v=numeric(m)+1
  A=t(v)%*%solve(sigma)%*%Er
  B=t(v)%*%solve(sigma)%*%v
  gamma=(A-as.numeric(lambda))/B
  w=(1/as.numeric(lambda))*solve(sigma)%*%(Er-as.numeric(gamma)*v)
  return(w)
}

######################## Etape 2 : modele de B-L etendu aux moments d'ordre 2 et 3 ##############################################################################

expected <- function(X,p,lambda){
  # Renvoie le rendement espere etant donne X et p 
  # en integrant les moments d'ordre superieur
  # lambda est le coeff d'aversion au risque
  w <- t(p)
  Rbarre <- colMeans(X) # Vecteur des rendements historiques moyens
  mu <- p%*%Rbarre
  print(paste('mu: ', mu))
  sigmaw <- sigma(X)%*%w
  mu2 <- p%*%sigmaw
  print(paste('mu2: ',mu2))
  omegaw <- coskewness(X)%*%kronecker(w,w)
  mu3 <- p%*%omegaw
  print(paste('mu3: ',mu3))
  psiw <- cokurtosis(X)%*%kronecker(kronecker(w,w),w)
  mu4 <- p%*%psiw
  print(paste('mu4:', mu4))
  mus <- list(mu,mu2,mu3,mu4)
  A <- 1 + ((lambda^2)*mu2/2) - ((lambda^3)*mu3/6) + ((lambda^4)*mu4/24)
  print(paste('A:',A))
  delta <-function(k){return(as.numeric((lambda^k)/(A*factorial(k))))}
  print(delta(1))
  print(dim(sigmaw))
  return(delta(1)*sigmaw - delta(2)*omegaw + delta(3)*psiw)
}

f_wBL_ext <- function(X,Er,lambda,w0){
  # Renvoie le portefeuille p de BL apres la premiere etape
  # w0 sert de root pour la resolution du probleme d'inversion, on aimerait prendre le wbl original par exemple
  v_ant <- volatiliteGARCH(X)
  expected_ant <- function(w){
    # Meme fonction que expected mais avec les valeurs anticipees
    p <- t(w)
    Rbarre <- colMeans(X) # Vecteur des rendements historiques moyens
    mu <- p%*%Rbarre
    print(paste('mu: ', mu))
    sigmaw <- sigma_ant(X,v_ant)%*%w
    mu2 <- p%*%sigmaw
    print(paste('mu2: ',mu2))
    omegaw <- coskewness_ant(X,v_ant)%*%kronecker(w,w)
    mu3 <- p%*%omegaw
    print(paste('mu3: ',mu3))
    psiw <- cokurtosis_ant(X,v_ant)%*%kronecker(kronecker(w,w),w)
    mu4 <- p%*%psiw
    print(paste('mu4:', mu4))
    mus <- list(mu,mu2,mu3,mu4)
    A <- 1 + ((lambda^2)*mu2/2) - ((lambda^3)*mu3/6) + ((lambda^4)*mu4/24)
    print(paste('A:',A))
    delta <-function(k){return(as.numeric((lambda^k)/(A*factorial(k))))}
    print(delta(1))
    print(dim(sigmaw))
    return(delta(1)*sigmaw - delta(2)*omegaw + delta(3)*psiw - Er)
  }
  # On calcule ensuite les poids w qui verifient expected_ant(w) = Er grace au package nleqslv
  return(nleqslv(w0,expected_ant)$x)
}

######################## Etape 3 : Pour tester la performance d'un portefeuille ##############################################################################
titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
rdt_j = rendements(datap_j)
rdt_j = global_return(datap_j)
rdt_j=rdt_j[,titres_selec_oat]
######################## Etape 4 : main.cpp ##############################################################################

main_bl <- function(){
  titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
  titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")
  
  #datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
  #rdt_j = rendements(datap_j)
  
  sigma<- sigma(rdt_j) #matrice des covariances
  fopt2<-function(Q,X,r) {
    
    n=length(X)
    v1=numeric(n)+1
    M=matrix(c(0,0,0,0),2,2)
    
    M[1,1]=t(X)%*%solve(Q)%*%X             #étape de calcul intermédiaire pour calculer les poids
    M[1,2]=t(X)%*%solve(Q)%*%v1  
    M[2,1]=t(v1)%*%solve(Q)%*%v1  
    M[2,2]=M[1,1]*M[2,1]-M[1,2]^2  
    
    poids=((M[2,1]*r-M[1,2])/M[2,2])*solve(Q)%*%X +((M[1,1]-M[1,2]*r)/M[2,2])*solve(Q)%*%v1  #vecteur de poids
    return(poids) 
    
  }
  
  poids_strat0 <- fopt2(sigma,colMeans(rdt_j),0.0012) #poids du portefeuille stratégique calculé via la fonction fopt2
  moy_rdt_hist=t(poids_strat0)%*%colMeans(rdt_j) #rendement historique du portefeuille 
  lambda <- (moy_rdt_hist-0.0012)/(t(poids_strat0)%*%sigma%*%poids_strat0) #coefficient d'aversion au risque
  mu <- expected(rdt_j,t(poids_strat0),lambda)
  P <- f_P(rdt_j)
  Q <- f_Q_garch(rdt_j)
  omega <- f_omega_garch(rdt_j)
  Er <- f_Er(rdt_j,1,mu,P,Q,omega)
  w0 <- f_wBL_somme1(rdt_j,lambda,Er)
  wBL_final <- f_wBL_ext(rdt_j,Er,lambda,w0)
  print(wBL_final)
}

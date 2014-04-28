source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/extraction data.r")
#source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/tests stats.r")
source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/opti pf.r")

#On dispose des donnees suivantes:
# X vecteur des rendements
# la matrice des corrélations historiques
# sigma la matrice de covariance
# rf rendement moyen de l'actif sans risque

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

r-av_coeff <- function(X,p,rf){
  #Renvoie le coefficient d'aversion au risque determine par le choix du portefeuille
  # X et p determinent le portefeuille, rf est le rendement de l'actif sans risque sur la periode
  return((meanreturn(X,p)-rf)/(volatility(X,p)))
}

pi <- function(p,rf,sigma){
  #Renvoie le vecteur "implied excess returns" qui represente les rendements optimaux 
  #etant donnes le choix d'allocation et la var historique
  #NB: si p est historiquement optimal, la fonction renvoie meanreturn(X,p)
  #X, p et rf comme d'habitude, sigma matrice de covariance
  return(r-av_coeff(X,p,rf)%*%sigma%*%p)
}

#Les anticipations nous donnent les elements suivants:
# la matrice P informe sur quels actifs sont concernes par quelle prediction
# le vecteur Q (derive de P) renvoie les valeurs de sur/sous-performance absolue pour chaque actif
# le coefficient "weight on views" thau
# la matrice omega donne les niveaux de confiance sur chaque vue

#on a alors (MODELE STANDARD):
p_bl <- function(X,p,rf,sigma,P,Q,thau,omega){
  # Renvoie l'allocation optimale etant donnees les anticipations
  impexc <- pi(p,rf,sigma)
  expect <- impexc + thau*sigma%*%t(P)%*%solve(omega+thau%*%P%*%sigma%*%t(P))%*%(Q - P%*%impexc)
  return(solve(r-av_coeff(X,p,rf)*sigma)%*%expect)
  #Attention, la somme des poids n'est pas generalement egale a 1, il peut y avoir des poids negatifs
  # et on suppose que le coeff d'aversion au risque est constant au cours du temps
}

#(MODELE ETENDU)

# On utilise:
# la matrice sigmma_prime des anticipations de covariances

m_kronecker <- function(v){
  souskro <- function(A,B){
    # Calcule le produit de 2 matrices
    m <- nrow(A)
    n <- ncol(A)
    p <- nrow(B)
    q <- ncol(B)
    print('PUSSY/')
    print(m)
    print(n)
    print(p)
    print(q)
    print('/Pussy')
    I <- m*p
    print(I)
    J <- n*q
    print(J)
    y <- array(dim=c(I,J))
    for (i in 1:I){
      for (j in 1:J){
        y[i,j] <- A[1+floor((i-1)/p),1+floor((j-1)/q)]*B[1+(i-1)%%p,1+(j-1)%%q]
      }
    }
    return(y)
  }
  # Le produit de Kronecker est associatif, cette fonction sert donc a calculer
  # le produit de 3 ou plus elements
  n <- length(v) # si n < 2 la fonction plante
  print(v[1])
  print(v[2])
  y <- kronecker(v[1],v[2])
  if (n > 2){
    for (i in 3:n){
      y <- kronecker(y,v[i])
    }
  }
  return(y)
}

coskewness <- function(X){
  # Renvoie la matrice de coskewness pour le portefeuille dont on a les rendements historiques X
  n <- ncol(X)
  T <- nrow(X)
  skew <- array(dim = c(n,n*n))
  sous_produit<-function(t){
    #Renvoie la matrice de coskewness pour une date t donnee
    x <- X[t,] - t(colMeans(X))
    return(m_kronecker(c(t(x),x,x)))
  }
  for (t in 1:T){
    skew <- skew + sous_produit(t)
  }
  return(skew/T)
}

cokurtosis <- function(X){
  # Renvoie la matrice de cokurtosis pour le portefeuille dont on a les rendements historiques X
  n <- ncol(X)
  T <- nrow(X)
  kur <- array(dim = c(n,n^3))
  sous_produit<-function(t){
    #Renvoie la matrice de coskewness pour une date t donnee
    x <- X[t,] - t(colMeans(X))
    return(m_kronecker(c(t(x),x,x,x)))
  }
  for (t in 1:T){
    skew <- skew + sous_produit(i)
  }
  return(skew/T)
}

expected <- function(X,p,lambda, sigma){
  # Renvoie le rendement espere etant donne X et p 
  # en integrant les moments d'ordre superieur
  # lambda est le coeff d'aversion au risque
  w <- t(p)
  Rbarre <- colMeans(X) # Vecteur des rendements historiques moyens
  mu <- p%*%Rbarre
  sigmaw <- sigma%*%w
  mu2 <- p%*%sigmaw
  omegaw <- coskewness(X)%*%m_kronecker(c(w,w))
  mu3 <- p%*%omegaw
  psiw <- cokurtosis(X)%*%m_kronecker(c(w,w,w))
  mu4 <- p%*%psiw
  mus <- c(mu,mu2,mu3,mu4)
  A <- 1 + ((lambda^2)*mu2/2) - ((lambda^3)*mu3/6) + ((lambda^4)*mu4/24)
  delta <-function(k){return((lambda^k)/(A*factorial(k)))}
  return(delta(1)*sigmaw - delta(2)*omegaw + delta(3)*psiw)
}

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
datap_m = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "M")

rdt_j = rendements(datap_j)
rdt_m = rendements(datap_m) 

SIGMA=cov(rdt_j) #matrice des covariances
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

poids_strat0=fopt2(SIGMA,colMeans(rdt_j),0.05) #poids du portefeuille stratégique calculé via la fonction fopt2
moy_rdt_hist=t(poids_strat0)%*%colMeans(rdt_j) #rendement historique du portefeuille 
lambda=(moy_rdt_hist-0.05)/(t(poids_strat0)%*%SIGMA%*%poids_strat0) #coefficient d'aversion au risque

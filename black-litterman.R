source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/extraction data.r")
source("/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/code/tests stats.r")
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

#(MODELE AMELIORE)

# On utilise:
# la matrice sigmma_prime des anticipations de covariances

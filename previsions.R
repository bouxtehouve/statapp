library(forecast)
library(fGarch)
library(fUnitRoots)

f_P <- function(X){
  # Fonction inutile mais c'est plus lisible. On ne donne que des vues absolues donc P est la matrice identite
  return(diag(ncol(X)))
}

######################## Methode 1: prevision avec ARMA ##############################################################################

# fonction de prevision basee sur des modeles ARIMA : vecteur de donnees x, horizon de prevision: n ( par défaut 60 ds la fonction)
previsionARMA<-function(X, n = 60) {
  u=ts(X)
  modele=auto.arima(u)    # retourne le meilleur ARIMA selon les critÃ¨res AIC, BIC 
  prev=forecast(modele,n) # retourne les prÃ©visions jusqu'Ã  t=n
  return(mean(prev$mean)) # on retourne la moyenne des prÃ©visions afin de construire le vecteur Q
}

f_Q <- function(X){
  # Renvoie le vecteur Q des niveaux d' "excess performance"
  return(apply(X,2,"previsionARMA"))
}

f_omega<-function(X,P){
  # construction de la matrice de confiance omega selon la methode automatique d'Idzorek:pk%*%sigma%*%pk'
  # pk vecteur ligne de la matrice P, sigma matrice des correlations
  sigma <- var(X)
  N=dim(P)[1]
  omega=matrix(0,N,N)
  for (i in 1:N) { 
    p=P[i,]
    omega[i,i]=t(p)%*%sigma%*%p
  }
  return(omega)
}

####################### Methode 2: prevision avec modele  GARCH #####################

perron_phillips <- function(X){
  # Verification de la stationnarite de nos series de rendements
  # Check par le test de Perron Phillips
  W=c()
  for (i in 1:dim(X)[2]) { 
    if (PP.test(ts(X[,i]),lshort=TRUE)$p.value>0.05) {W=c(W,i)}
  }
  return(W)
}

# -->tous les rendements sont stationnaires
# On peut donc calculer des GARCH sans differencier la serie

# modele GARCH permettant d'avoir des previsions sur la volatilite
# on postule des GARCH(1,1) pour chacune de nos actions (courant)

previsionGARCH <-function(X, m = 60){
  fit=garchFit(~garch(1, 1),X,trace = FALSE )
  a=predict(fit,n.ahead=m)$meanForecast
  return(mean(a))
}

volatiliteGARCH <-function(X, m = 60){
  fit=garchFit(~garch(1, 1),X,trace = FALSE )
  a=predict(fit,n.ahead=m)$standardDeviation
  return(mean(a))
}

f_Q_garch <- function(X){
  # Renvoie le vecteur Q des niveaux d' "excess performance" avec prevision garch
  return(apply(X,2,"previsionGARCH"))
}

f_omega_garch<-function(X){
  # Matrice de niveau de confiance dans les vues Omega
  # on a la variance pour m dates via le modele GARCH--> on prend la moyenne des variances comme estimateurs
  # Methode Mankert --> matrice incorporant les variances des vues
  N=dim(X)[2]
  omega=matrix(0,N,N)
  for (i in 1:N) {omega[i,i]= mean(volatiliteGARCH(X[,i])^2)}
  return(omega)
}
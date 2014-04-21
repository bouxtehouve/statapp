#A coller dans application

######################## Black litterman ########################################################################

library(forecast)
library(rugarch)
library(fGarch)

###################### FONCTIONS OPTIMISATIONS #################################################################
fopt<-function(Q,X,r,rf) {
  n=length(X)
  e=numeric(n)+1
  v=c()
  m=c()
  
  u=X-rf*e
  v[1]=t(e)%*%solve(Q)%*%u
  v[2]=t(u)%*%solve(Q)%*%u
  
  poidsrisque=((r-rf)*solve(Q)%*%u)/v[2]    #vecteur de poids des actifs risqu?s
  poidssansrisque=1-t(e)%*%poidsrisque       #poids de l'actif sans risque
  m=c(poidssansrisque,poidsrisque)
  
  return(m)
}

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
#méthode d'optimisation via la fonction d'utilité utilisée ds Black litterman
#fopt_utilite=function(sigma,mu){
#  m=length(mu)
#  v=numeric(m)+1
#  A=t(v)%*%solve(sigma)%*%mu
#  B=t(v)%*%solve(sigma)%*%v
#  gamma=(A-as.numeric(lambda))/B
#  w=(1/as.numeric(lambda))*solve(sigma)%*%(mu-as.numeric(gamma)*v)
# return(w)
#}
#opti=fopt_utilite(SIGMA,colMeans(rdt_j[,titres_selec_oat]))

####################### Première étape ##########################################################################

# Vector of implied excess return (rendement portefeuille stratégique)

rf=mean(rdt_j[,"OAT"]) #rendement de l'actif sans risque
SIGMA=var(rdt_j[,titres_selec_oat]) #matrice des covariances

# je pense qu'on optimise sans prendre en compte l'actif sans risque
# m=length(fopt(SIGMA,colMeans(rdt_j[,titres_selec_oat]),0.003,mean(rdt_j[,"OAT"])))
# poids_strat=fopt(SIGMA,colMeans(rdt_j[,titres_selec_oat]),0.003,mean(rdt_j[,"OAT"]))[1:(m-1)] #on retire le dernier poids correspondant à celui de l'acitf sans risque

poids_strat=fopt2(SIGMA,colMeans(rdt_j[,titres_selec_oat]),0.003)
moy_rdt_hist=t(poids_strat)%*%colMeans(rdt_j[,titres_selec_oat]) #rendement historique du portefeuille 

lambda=(moy_rdt_hist-rf)/(t(poids_strat)%*%SIGMA%*%poids_strat)
as.numeric(lambda)
rdt_strat=as.numeric(lambda)*(SIGMA%*%poids_strat) #vector of implied excess return de Black Litterman


######################### Deuxième étape ########################################################################
## on utilise la matrice des rendements net de l'actif sans risque pour pouvoir faire nos prévisions
N=dim(rdt_j)
rdt_net=(rdt_j-rdt_j[,N[2]]) #on soustrait aux rendements celui de l'actif sans risque
rdt_net_trunc=rdt_net[1:N[1],1:(N[2]-1)]

# matrice P: on ne donne que des vues absolues --> P=Identité(16,16)

######################## Méthode 1: prévision avec ARMA ##############################################################################

# fonction de prévision basée sur des modèles ARIMA : vecteur de données x, horizon de prévision: n (ds la fonction)
previsionARMA<-function(x) {
  n=20
  u=ts(x)
  modele=auto.arima(u)    # retourne le meilleur ARIMA selon les critères AIC, BIC 
  prev=forecast(modele,n) # retourne les prévisions jusqu'à t=n
  return(mean(prev$mean)) # on retourne la moyenne des prévisions afin de construire le vecteur Q
}

Q=apply(rdt_net_trunc,2,"previsionARMA")  #construction du vecteur Q contenant nos prev pr chaque actif
summary(Q)


# construction de la matrice de confiance omega selon la méthode d'Idzorek:pk%*%sigma%*%pk'
# pk vecteur ligne de la matrice P, sigma matrice des corrélations

f_omega<-function(P,Sigma){
  N=dim(P)[1]
  omega=matrix(0,N,N)
  for (i in 1:N) { 
    p=P[i,]
    omega[i,i]=t(p)%*%Sigma%*%p
  }
  return(omega)
}

OMEGA=f_omega(diag(15),cov(datap_j[,titres_selec_oat]))

# CF littérature , on fixe tau à 1 (Satchell et Scowcroft)
tau=1

########################### Troisième étape ####################################################################
# rendement espéré ds le modele de Black Litterman

rdt_prev=tau*SIGMA%*%solve(OMEGA+tau*SIGMA)%*%(Q-rdt_strat)
rdt_litterman=rdt_strat+rdt_prev

# portefeuille Black Litterman

WBL=solve(as.numeric(lambda)*SIGMA)%*%rdt_litterman
t(WBL)%*%colMeans(rdt_j[,titres_selec_oat]) 
sum(WBL) #la somme des poids est différente de 1 ds le modèle de BL

# portefeuille BL amélioré pour que les poids soient égaux à 1
WBL2=fopt2(SIGMA,rdt_litterman,0.003) 
t(WBL2)%*%colMeans(rdt_j[,titres_selec_oat])


####################### Modele 2:prévision avec modèle  GARCH et matrice de covariance anticipée ################

# Vérification de la stationnarité de nos séries de rendements
# Check par le test de Perron Phillips

library(fUnitRoots)

W=c()
for (i in 1:dim(rdt_net_trunc)[2]) { 
  
  if (PP.test(ts(rdt_net_trunc[,i]),lshort=TRUE)$p.value>0.05) {W=c(W,i)}
  
}
W

# -->tous les rendements sont stationnaires
# On peut donc calculer des GARCH sans différencier la série

# modèle GARCH permettant d'avoir des previsions sur la volatilité
# on postule des GARCH(1,1) pour chacune de nos actions (courant)

# prévision à l'horizon m=100
previsionGARCH<-function(x){
  m=20
  spec=ugarchspec()
  modele=ugarchfit(spec,x)
  prev=ugarchforecast(modele,n.ahead=m)
  return(mean(fitted(prev)))
}

volatiliteGARCH<-function(x){
  m=20
  spec=ugarchspec()
  modele=ugarchfit(spec,x)
  prev=ugarchforecast(modele,n.ahead=m)
  return(sigma(prev))
}

# Nouveau vecteur Q avec prev garch
Q_garch=apply(rdt_net_trunc,2,"previsionGARCH")  #construction du vecteur Q contenant nos prev pr chaque actif
summary(Q_garch)

# Matrice de niveau de confiance dans les vues Omega
# on a la variance pour m dates via le modèle GARCH--> on prend la moyenne des variances comme estimateurs
# Méthode Mankert --> matrice incorporant les variances des vues
f_omega_garch<-function(x){
  N2=dim(x)[2]
  omega=matrix(0,N2,N2)
  for (i in 1:N2) { omega[i,i]= mean(volatiliteGARCH(x[,i])^2)}
  return(omega)
  
}

omega_garch = f_omega_garch(rdt_net_trunc)


# Construction de la matrice de variance anticipée
# hypothèse: stabilité des corrélations dans le temps

# Corrélations historiques
rho_hist=cor(rdt_net_trunc)
# cf p13 WERLé 
SIGMA_forecast=matrix(0,dim(rdt_net_trunc)[2],dim(rdt_net_trunc)[2])
for (i in 1:dim(rdt_net_trunc)[2]) {
  for (j in 1:dim(rdt_net_trunc)[2]) {
    SIGMA_forecast[i,j]=rho_hist[i,j]*sqrt(omega_garch[i,i])*sqrt(omega_garch[j,j])
  }
}

# rendement espéré ds le modele de Black Litterman

rdt_prev_garch=tau*SIGMA_forecast%*%solve(omega_garch+tau*SIGMA_forecast)%*%(Q_garch-rdt_strat)
rdt_litterman_garch=rdt_strat+rdt_prev_garch

WBL_garch=solve(as.numeric(lambda)*SIGMA_forecast)%*%rdt_litterman_garch # la somme des poids est différente de 1

WBL2_garch=fopt2(SIGMA_forecast,rdt_litterman_garch,0.003) # tq la somme des poids vaut 1


# test sur les données de 2007

datap_j_2007=load_data("01/01/2006", "20/01/2006", titres = titres_selec_oat, type = "J")
rdt_j_2007 = rendements(datap_j_2007)

# Résultats Sur données réelles impacte 

t(poids_strat)%*%colMeans(rdt_j_2007)
t(WBL2)%*%colMeans(rdt_j_2007)
t(WBL2_garch)%*%colMeans(rdt_j_2007)

# ratio de Sharpe

S_garch=(t(WBL2_garch)%*%colMeans(rdt_j_2007))/sqrt(t(WBL2_garch)%*%cov(datap_j_2007)%*%WBL2_garch)
S_WBL2=(t(WBL2)%*%colMeans(rdt_j_2007))/sqrt(t(WBL2)%*%cov(datap_j_2007)%*%WBL2)
S_strat=(t(poids_strat)%*%colMeans(rdt_j_2007))/sqrt(t(poids_strat)%*%cov(datap_j_2007)%*%poids_strat)

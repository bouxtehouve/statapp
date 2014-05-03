source("/Users/dano/Desktop/essai victor/extraction data.r")
source("/Users/dano/Desktop/essai victor/tests stats.r")
source("/Users/dano/Desktop/essai victor/opti pf.r")

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","oat","cac40")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")

rdt_j = global_return(datap_j)
rdt_j_a=rdt_j[,titres_selec_oat]
#rdt_m = global_return(datap_m) 

######################## Black litterman ########################################################################

library(forecast)
library(fGarch)

###################### FONCTIONS OPTIMISATIONS #######################################################################
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
####################### Première étape ###############################################################################
#on a vu precedemment que l'on ne prend pas en compte l'actif sans risque dans notre portefeuille
# Vector of implied excess return (rendement portefeuille stratégique)
rf=mean(rdt_j[,"OAT"])
SIGMA=cov(rdt_j_a) #matrice des covariances

poids_strat0=fopt2(SIGMA,colMeans(rdt_j_a),0.001169648) #poids du portefeuille stratégique calculé via la fonction fopt2/ on fixe le rendement à r=0.01169648 (optimal selon la frontière efficiente)
moy_rdt_hist=t(poids_strat0)%*%colMeans(rdt_j_a) #rendement historique du portefeuille 
lambda=(moy_rdt_hist-rf)/(t(poids_strat0)%*%SIGMA%*%poids_strat0) #coefficient d'aversion au risque

############################## FONCTIONS OPTIMISATIONS BLACK LITTERMAN ################################################
#méthode d'optimisation via la fonction d'utilité utilisée ds Black litterman
# -> mzximiser en w , U = t(w)%*%mu -(lambda)/2*t(w)%*%sigma%*%w
fopt_u<-function(sigma,mu) {
  w=(1/as.numeric(lambda))*solve(sigma)%*%mu
  return(w)
}

#2 ème fonction d'optimisation avec contrainte que la somme des poids soit égale à 1
fopt_uc<-function(sigma,mu){
  m=length(mu)
  v=numeric(m)+1
  A=t(v)%*%solve(sigma)%*%mu
  B=t(v)%*%solve(sigma)%*%v
  gamma=(A-as.numeric(lambda))/B
  w=(1/as.numeric(lambda))*solve(sigma)%*%(mu-as.numeric(gamma)*v)
  return(w)
}
######################################################################################################################

poids_strat1=fopt_uc(SIGMA,colMeans(rdt_j_a))
rdt_strat1=as.numeric(lambda)*(SIGMA%*%poids_strat1)
mean(rdt_strat1)
######################### Deuxième étape #############################################################################

# matrice P: on ne donne que des vues absolues --> P = Identité(16,16)

######################## Méthode 1: prévision avec ARMA ##############################################################################

# fonction de prévision basée sur des modèles ARIMA : vecteur de données x, horizon de prévision: n (ds la fonction)
previsionARMA<-function(x) {
  n=60
  u=ts(x)
  modele=auto.arima(u)    # retourne le meilleur ARIMA selon les critères AIC, BIC 
  prev=forecast(modele,n) # retourne les prévisions jusqu'à t=n
  return(mean(prev$mean)) # on retourne la moyenne des prévisions afin de construire le vecteur Q
}

Q=apply(rdt_j_a,2,"previsionARMA")  #construction du vecteur Q contenant nos previsions pr chaque actif
summary(Q)


# construction de la matrice de confiance omega selon la méthode automatique d'Idzorek:pk%*%sigma%*%pk'
# pk vecteur ligne de la matrice P, sigma matrice des corrélations

f_omega<-function(P,Sigma){
  N=dim(P)[1]
  omega=matrix(0,N,N)
  for (i in 1:N) { 
    p=P[i,]
    omega[i,i]=t(p)%*%SIGMA%*%p
  }
  return(omega)
}

OMEGA=f_omega(diag(15),cov(rdt_j_a))

# CF littérature , on fixe tau à 1 (Satchell et Scowcroft)
tau=1

########################### Troisième étape ##########################################################################
# rendement espéré ds le modele de Black Litterman

rdt_prev=tau*SIGMA%*%solve(OMEGA+tau*SIGMA)%*%(Q-rdt_strat1)
rdt_litterman=rdt_strat1+rdt_prev

# portefeuille Black Litterman

WBL=solve(as.numeric(lambda)*SIGMA)%*%rdt_litterman 
sum(WBL) #la somme des poids est différente de 1 ds le modèle de BL

# portefeuille BL amélioré pour que les poids soient égaux à 1
WBL2=fopt_uc(SIGMA,rdt_litterman) 

####################### Modele 2:prévision avec modèle  GARCH et matrice de covariance anticipée #####################

# Vérification de la stationnarité de nos séries de rendements
# Check par le test de Perron Phillips

library(fUnitRoots)

W=c()
for (i in 1:dim(rdt_j_a)[2]) { 
  
  if (PP.test(ts(rdt_j_a[,i]),lshort=TRUE)$p.value>0.05) {W=c(W,i)}
  
}
W

# -->tous les rendements sont stationnaires
# On peut donc calculer des GARCH sans différencier la série

# modèle GARCH permettant d'avoir des previsions sur la volatilité
# on postule des GARCH(1,1) pour chacune de nos actions (courant)

previsionGARCH<-function(x){
  m=60
  fit=garchFit(~garch(1, 1),x,trace = FALSE )
  a=predict(fit,n.ahead=m)$meanForecast
  return(mean(a))
}

volatiliteGARCH<-function(x){
  m=60
  fit=garchFit(~garch(1, 1),x,trace = FALSE )
  a=predict(fit,n.ahead=m)$standardDeviation
  return(mean(a))
}

# Nouveau vecteur Q avec prev garch
Q_garch=apply(rdt_j_a,2,"previsionGARCH") #construction du vecteur Q contenant nos previsions pr chaque actif (égale à la moyenne pour un GARCH)
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

omega_garch = f_omega_garch(rdt_j_a)


# Construction de la matrice de variance anticipée
# hypothèse: stabilité des corrélations dans le temps

# Corrélations historiques
rho_hist=cor(rdt_j_a)
# cf p13 WERLé 
g=c() #vecteur qui contient les volatilités prédites de chaque actif
for (i in 1:dim(rdt_j_a)[2]) {g=c(g,mean(volatiliteGARCH(rdt_j_a[,i])))}

SIGMA_forecast=matrix(0,dim(rdt_j_a)[2],dim(rdt_j_a)[2])
for (i in 1:dim(rdt_j_a)[2]) {
  for (j in 1:dim(rdt_j_a)[2]) {
    SIGMA_forecast[i,j]=rho_hist[i,j]*g[i]*g[j]
  }
}

# rendement espéré ds le modele de Black Litterman

rdt_prev_garch=tau*SIGMA_forecast%*%solve(omega_garch+tau*SIGMA_forecast)%*%(Q_garch-rdt_strat1)
rdt_litterman_garch=rdt_strat1+rdt_prev_garch

WBL_garch=solve(as.numeric(lambda)*SIGMA_forecast)%*%rdt_litterman_garch # la somme des poids est différente de 1

WBL2_garch=fopt_uc(SIGMA_forecast,rdt_litterman_garch) # tq la somme des poids vaut 1

# test sur les données de 2007

datap_j_2007=load_data("01/01/2007", "01/05/2007", titres = titres_selec, type = "J")
rdt_j_2007 = rendements(datap_j_2007)
rdt_j_a_2007=rdt_j_2007[,titres_selec_oat]
SIGMA_2007=cov(rdt_j_a_2007)

# Rendements sur la période 2007

rdt_strat_2007=t(poids_strat1)%*%colMeans(rdt_j_a_2007)
rdt_arma_2007=t(WBL2)%*%colMeans(rdt_j_a_2007)
rdt_garch_2007=t(WBL2_garch)%*%colMeans(rdt_j_a_2007)

# ratios de Sharpe

#rdt=matrice des rendements , r=rendement de l'actif sans risque, sigma= matrice des corrélations
f_sharpe<-function(poids,rdt,r,sigma){                      
  a=(t(poids)%*%colMeans(rdt)-r)/sqrt(t(poids)%*%sigma%*%poids)
  return(a)
}

f_sharpe(poids_strat1,rdt_j_a_2007,mean(rdt_j_2007[,"OAT"]),SIGMA_2007)
f_sharpe(WBL2,rdt_j_a_2007,mean(rdt_j_2007[,"OAT"]),SIGMA_2007)
f_sharpe(WBL2_garch,rdt_j_a_2007,mean(rdt_j_2007[,"OAT"]),SIGMA_2007)

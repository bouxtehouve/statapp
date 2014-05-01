source("/Users/dano/Desktop/essai victor/extraction data.r")
source("/Users/dano/Desktop/essai victor/tests stats.r")
source("/Users/dano/Desktop/essai victor/opti pf.r")

library(tseries)
library(fUnitRoots)
library(forecast)

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
datap_m = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "M")

rdt_j = global_return(datap_j)
rdt_j_a=rdt_j[,titres_selec_oat]

#Test de racine unité de Perron Phillips
W=c()
for (i in 1:dim(rdt_j_a)[2]) { 
  
  if (PP.test(ts(rdt_j_a[,i]),lshort=TRUE)$p.value>0.05) {W=c(W,i)}
  
}
W
# conclusion:test de racine unité rejeté --> stationnarité

#test de racine unité Augmented Dickey Fuller
W1=c()

for (i in 1:dim(rdt_j_a)[2]) {
  if (adfTest(ts(rdt_j_a[,i]),lags=10,"nc")@test$p.value>0.05) {W1=c(W1,i)}
}
W1
#conclusion: test de racine unité rejeté--> stationnarité

#Test de Blancheur des résidus et de normalité

#résidus ARMA
residusARMA<-function(x) {
  u=ts(x)
  modele=auto.arima(u)    # retourne le meilleur ARIMA selon le critère AIC
  return(modele$res) 
}

#calcul de p+q=nombre de coeff du modèle hors constante
nbcoeffARMA<-function(x){
  u=ts(x)
  modele=auto.arima(u)
  return(length(modele$coef))
}

M1=apply(rdt_j_a,2,"residusARMA") #on récupère tous les résidus des modèles ARMA
M2=apply(rdt_j_a,2,"nbcoeffARMA") #nombre de coeff hors constante

#Test de Ljung Box
W2=c()
for (i in 1:dim(M1)[2]){
  if (Box.test(M1[,i], lag = 10, type = c( "Ljung-Box"), fitdf = M2[i])$p.value>0.05) {W2=c(W2,i)}
}
W2 #length(W2)>0 on rejette la blancheur des résidus globale


#Test de normalité des résidus : Jarque Bera
W3=c()
for (i in 1:dim(M1)[2]){
  if (jarque.bera.test(M1[,i])$p.value>0.05) {W3=c(W3,i)}
}
W3
#conclusion: test de normalité rejeté

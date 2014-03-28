source("D:/cours/2A/Statap/Scripts R/extraction data.r")
source("D:/cours/2A/Statap/Scripts R/tests stats.r")
source("D:/cours/2A/Statap/Scripts R/opti pf.r")

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
datap_m = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "M")

rdt_j = rendements(datap_j)
rdt_m = rendements(datap_m)

# Produit une ACP sur les cours des actions pour obtenir le cercle des corrélations grâce au package suivant
library("FactoMineR")
#PCA(rdt_j,graph=T)

plot(datap_j[,"vinci"],type='l')

P=fopt(var(datap_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),0.003,mean(rdt_j[,"OAT"]))
f_esp(datap_j,P)

fronteffi(var(datap_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),mean(rdt_j[,"OAT"]))

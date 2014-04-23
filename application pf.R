source("D:/cours/2A/Statap/Scripts R/extraction data.r")
source("D:/cours/2A/Statap/Scripts R/tests stats.r")
source("D:/cours/2A/Statap/Scripts R/opti pf.r")

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","OAT")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")
datap_m = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "M")

rdt_j = rendements(datap_j)
rdt_m = rendements(datap_m)

# Produit une ACP sur les cours des actions pour obtenir le cercle des corr?lations gr?ce au package suivant
library("FactoMineR")
#PCA(rdt_j,graph=T)

plot(datap_j[,"vinci"],type='l')

P=fopt(var(datap_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),0.003,mean(rdt_j[,"OAT"]))
f_esp(datap_j,P)

fronteffi(var(datap_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),mean(rdt_j[,"OAT"]))

#Frontière efficiente
library(fPortfolio)
rdt_journalier=as.timeSeries(rdt_j[,titres_selec_oat]) #necessaire de transformer rdt_j en séries temporelles

spec=portfolioSpec() #description des spécificités de notre portefeuille
setTargetReturn(spec)=mean(colMeans(X)) # spécifier un rendement objectif
spec

constraints=c("minW[1:n]=-1") #spécification des contraintes, ici on autorise des poids négatifs
portfolioConstraints(X, spec,constraints) # description des contraintes

frontier <- portfolioFrontier(X, spec,constraints) # calcul l'ensemble des points de la frontière efficiente
print(frontier)

tailoredFrontierPlot(object = frontier) # tracé de la frontière efficiente

weightsPlot(frontier,col= palette(rainbow(20))) # graphique permettant de voir les poids de chaque actifs pour chaque portefeuille de la frontière efficiente

efficientPortfolio(X, spec,constraints) #  calcul le portefeuille optimal compte tenu du rendement objectif et des contraintes imposées
tangencyPortfolio(X,spec,constraints) # calcul le portefeuille ayant le meilleur rapport rendement/risque de la frontière efficiente


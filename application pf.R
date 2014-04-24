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

rdt_journalier=as.timeSeries(rdt_j[,titres_selec_oat]) #necessaire de transformer rdt_j en séries temporelles avec la library fPortfolio

spec=portfolioSpec() #description des spécificités de notre portefeuille
setTargetReturn(spec)=mean(colMeans(rdt_journalier)) # spécifier un rendement objectif
spec

constraints=c("minW[1:n]=-1") #spécification des contraintes, ici on autorise des poids négatifs
portfolioConstraints(rdt_journalier, spec,constraints) # description des contraintes

frontier <- portfolioFrontier(rdt_journalier, spec,constraints) # calcul l'ensemble des points de la frontière efficiente
print(frontier)

tailoredFrontierPlot(object = frontier) # tracé de la frontière efficiente

setTargetReturn(spec)=0.0012 #on peut augmenter le rendement objectif , cf résultat frontier et frontière efficiente

weightsPlot(frontier,col= palette(rainbow(20))) # graphique permettant de voir les poids de chaque actifs pour chaque portefeuille de la frontière efficiente

efficientPortfolio(rdt_journalier, spec,constraints) #  calcul le portefeuille optimal compte tenu du rendement objectif et des contraintes imposées
tangencyPortfolio(rdt_journalier,spec,constraints) # calcul le portefeuille ayant le meilleur rapport rendement/risque de la frontière efficiente 

# on peu également utiliser nos propres fonctions pour trouver les poids optimaux
# d'après le graphe de la frontière efficiente, le meilleur rendement est 0.12%

fopt2(cov(rdt_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),0.0012) #on obtient des poids identiques à ceux obtenus avec efficientPortfolio

# Frontière efficiente Var-mean

# Il peut aussi être interessant de comparer les rendements en fonction d' une autre mesure de risque: VAR
# l'option frontier nous retourne 50 (nb modifiable) portefeuilles sur la frontière efficiente. on peut utiliser les poids de ces portefeuilles pour calculer les VAR95%

poids_effi=attributes(frontier@portfolio)$portfolio$weights #retourne les poids de chacun des 50 portefeuilles sur la frontière efficiente

#fonctions utiles

somme_col<-function(M){ #somme le colonnes d'une matrice
  n=dim(M)[2]
  v=numeric(n)
  for (i in 1:n) {
    v=v+M[,i]
  }
  return(v)
}

combine<-function(w,M){ #multiplie chaque colonne d'une matrice par les poids d'un vecteur 
  n=dim(M)[1]
  m=dim(M)[2]
  A=matrix(0,n,m)
  for (i in 1:m) {A[,i]=as.numeric(w[i])*M[,i]}
  return(A)
}

#Matrice des rendements de tous les portefeuilles de la frontière efficiente

rdt_fronteffi<-function(W,R){ 
  n=dim(W)[1]
  s=dim(R)[1]
  M=matrix(0,s,n)
  for (i in 1:n) {M[,i]=somme_col(combine(W[i,],R))}
  return(M)
}

rdt_effi=rdt_fronteffi(poids_effi,rdt_j[,titres_selec_oat]) 

#calcul de la VAR95% pour chaque portefeuille de la frontière efficiente

frontierRETURN_VAR<-function(M){
  
  Var=c()
  rendement=c()
  
  for (i in 1:dim(M)[2]) {Var=c(Var,f_var(M,i,0.95))
                          rendement=c(rendement,mean(M[,i]))
  }
  M=cbind(Var,rendement)
  return(M)
}

#plot frontière efficiente Return-Var95%
plot(frontierRETURN_VAR(rdt_effi)[,"Var"],frontierRETURN_VAR(rdt_effi)[,"rendement"],type="l",col="red",main=expression("Frontière efficiente rendement/Var95%"),ylab="rendement",xlab="Var95%"))




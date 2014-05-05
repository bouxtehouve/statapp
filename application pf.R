path="C:/Users/Yann/Documents/GitHub/statapp"
source(paste(path,"extraction data.r",sep="/"))
source(paste(path,"tests stats.r",sep="/"))
source(paste(path,"opti pf.r",sep="/"))
source(paste(path,"black_litterman_code_kevin.r",sep="/"))
source(paste(path,"black-litterman.r",sep="/"))
source(paste(path,"frontiere_efficiente.r",sep="/"))

titres_selec=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini","oat","cac40")
titres_selec_oat=c("orange","tf1","airliquide","alcatel-lucent","carrefour","kering","loreal","peugeot","thales","bnp","bouygues","soge","total","vinci","capgemini")

datap_j = load_data("01/01/2003", "31/12/2006", titres = titres_selec, type = "J")

rdt_j = global_return(datap_j)

# Produit une ACP sur les cours des actions pour obtenir le cercle des corr?lations gr?ce au package suivant
library(FactoMineR)
PCA(rdt_j,graph=T)

#Frontière efficiente
library(fPortfolio)

rdt_journalier=as.timeSeries(global_return(datap_j)) # necessaire de transformer rdt_j en séries temporelles avec la library fPortfolio
rdt_journalier_oat=rdt_journalier[,titres_selec_oat] # sans l'oAT

spec=portfolioSpec() #description des spécificités de notre portefeuille
setTargetReturn(spec)=mean(colMeans(rdt_journalier)) # spécifier un rendement objectif
setRiskFreeRate(spec)=mean(rdt_journalier[,"oat"]) # rendement moyen de l'actif sans risque
spec

n=ncol(rdt_journaliers_oat)
constraints=c("minW[1:n]=-1") #spécification des contraintes, ici on autorise des poids négatifs
portfolioConstraints(rdt_journalier_oat, spec,constraints) # description des contraintes

frontier <- portfolioFrontier(rdt_journalier_oat, spec,constraints) # calcul l'ensemble des points de la frontière efficiente
print(frontier)

#tailoredFrontierPlot(object = frontier) # tracé de la frontière efficiente
#ne pas oublier de mettre le fichier frontiere_efficiente en source
frontiere_efficiente1(object = frontier) # tracé de la frontière efficiente

setTargetReturn(spec)=0.0012 #on peut augmenter le rendement objectif , cf résultat frontier et frontière efficiente

weightsPlot(frontier,col= palette(rainbow(20))) # graphique permettant de voir les poids de chaque actifs pour chaque portefeuille de la frontière efficiente

efficientPortfolio(rdt_journalier_oat, spec,constraints) #  calcul le portefeuille optimal compte tenu du rendement objectif et des contraintes imposées
tangencyPortfolio(rdt_journalier_oat,spec,constraints) # calcul le portefeuille ayant le meilleur rapport rendement/risque de la frontière efficiente (c-a-d le meilleur ratio de Sharpe)
                                                       # on selectionne ce portefeuille pariculier
# on peu également utiliser nos propres fonctions pour trouver les poids optimaux
# d'après  la frontière efficiente, le meilleur rendement est 0.001169648

fopt2(cov(rdt_j[,titres_selec_oat]),colMeans(rdt_j[,titres_selec_oat]),0.001169648) #on obtient des poids identiques à ceux obtenus avec efficientPortfolio

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
  
  for (i in 1:dim(M)[2]) {Var=c(Var,f_var(M[,i],0.05))
                          rendement=c(rendement,mean(M[,i]))
  }
  M=cbind(Var,rendement)
  return(M)
}

#plot frontière efficiente Return-Var95%
plot(frontierRETURN_VAR(rdt_effi)[,"Var"],frontierRETURN_VAR(rdt_effi)[,"rendement"],type="l",col="red",main=expression("Frontière efficiente rendement/Var95%"),ylab="rendement",xlab="Var95%"))

#Frontière efficiente CVar-Mean

spec1=portfolioSpec() # spécification du portefeuille
setType(spec1)="CVaR" # on choisi un critère de minimisation de Cvar
setAlpha(spec1)=0.05 # niveau de confiance de la Cvar
setSolver(spec1)="solveRglpk" # l'optimisation de la Cvar est un pb linéaire --> on change de solveur
constraints1=c("minW[1:6]=-999", "maxW[1:6]=+999")

frontier_Cvar= portfolioFrontier(rdt_journalier_oat,spec1,constraints1) # calcul l'ensemble des points de la frontière efficiente
print(frontier_Cvar)

tailoredFrontierPlot(frontier_Cvar,risk = "CVaR") # tracé de la frontière efficiente Mean-Cvar



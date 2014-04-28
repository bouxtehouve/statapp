###### Ce fichier produit toutes les statistiques et tests utilisés dans notre mémoire

library("moments")

## TESTS
# 1) Test de Kolmogoroff-Smirnoff
f_ks <- function(Y,i){
  return(ks.test(Y[,i],pnorm(mean(Y[,i]),sqrt(var(Y[,i]))))[2])
}

# 2) Test de Shapiro
f_shapiro <- function(Y,i){
  return(shapiro.test(Y[,i])[2])
}

# 3) Test de Jarque-Bera
f_jb<- function(Y,i){
  return(jarque.test(Y[,i])[2])
}

# Tracé de l'histogramme des rendements, visuel par rapport au suivi d'une loi normale
f_graph <- function(Y,i,n){
  plot.new()
  hist(Y[,i],prob=T,nclass=n,main=paste("Histogramme et densités du rendement de ",toupper(titres_presents[i]),sep=""))
  curve(dnorm(x,mean(Y[,i]),sqrt(var(Y[,i]))),from=-0.5,to=0.5,col='red',add=T)
  lines(density(Y[,i]),col='blue')
  legend("topleft",c("Densité empirique","Densité théorique"),col=c('blue','red'),bty="n",lty=c("solid","solid"))
}


## STATISTIQUES
# 1) Ratio de Sharpe
f_sharpe <- function(Y,i,r){
  return((mean(Y[,i])-r)/sqrt(var(V[i,i])))
}

# 2) VaR
f_var <- function(Y,i,alpha){
  return(quantile(Y[,i],alpha))
}

# 3) CVaR
f_cvar <- function(Y,i,alpha){
  var=quantile(Y[,i],1-alpha)
  return(sum(Y[,i]*(Y[,i]<var))/sum((Y[,i]<var)))
}
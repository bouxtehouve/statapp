###### Ce fichier produit toutes les statistiques et tests utilis�s dans notre m�moire

library(moments)

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

# Trac� de l'histogramme des rendements, visuel par rapport au suivi d'une loi normale
f_graph <- function(Y,i,n){
  plot.new()
  hist(Y[,i],prob=T,nclass=n,main=paste("Histogramme et densit�s du rendement de ",toupper(titres_selec[i]),sep=""))
  curve(dnorm(x,mean(Y[,i]),sqrt(var(Y[,i]))),from=-0.6,to=0.6,n=10000,col='red',add=T)
  lines(density(Y[,i]),col='blue')
  legend("topleft",c("Densit� empirique","Densit� th�orique"),col=c('blue','red'),bty="n",lty=c("solid","solid"))
}


## STATISTIQUES
# 1) Ratio de Sharpe
f_sharpe <- function(Y,i,r){
  return((mean(Y[,i])-r)/sqrt(var(V[i,i])))
}

# 2) VaR
f_var <- function(Y,alpha){
  return(quantile(Y,alpha))
}

# 2bis) VaR approxim�e de Cornish Fisher
f_var_cf <- function(Y,alpha,T=1){
  q=qnorm(alpha)
  qcf=q+skewness(Y)*(q^-1)/6+kurtosis(Y)*(q^3-3*q)/24-skewness(Y)^2*(2*q^3-5*q)/36
  return(T*mean(Y)+qcf*sqrt(T*var(Y)))
}

# 3) CVaR
f_cvar <- function(Y,alpha){
  var=quantile(Y,alpha)
  return(sum(Y*(Y<var))/sum((Y<var)))
}

# 3bis) CVaR par calcul int�gral
f_cvar_int<-function(Y,alpha){
  f<-function(a){quantile(rdt_pf,a)}
  return(integrate(f,lower=0,upper=alpha)$value/alpha)
}

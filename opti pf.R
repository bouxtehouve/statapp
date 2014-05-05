###### Ce fichier contienrt tous les calculs nécessaires à l'optimisation de notre portefeuille

## Les arguments utilisés sont les suivants :
# Q matrice de variance
# X vecteur rendement
# r rendement attendu
# rf rendement de l'actif sans risque


# 1) Calculs des volatilités

# Expression de la volatilité en fonction du rendement avec actifs risqués
volat<-function(Q,X,r) {
  n=length(X)
  e=numeric(n)+1
  A=t(e)%*%solve(Q)%*%X
  B=t(X)%*%solve(Q)%*%X
  C=t(e)%*%solve(Q)%*%e
  D=C*B-A^2
  sigma=sqrt((C*r^2-2*A*r+B)/D)
  
  return(sigma)
}

# Expression de la volatilité en fonction du rendement attendu avec actif sans risque
volat2<-function(Q,X,r,rf) {
  n=length(X)
  e=numeric(n)+1
  u=X-rf*e
  B=t(u)%*%solve(Q)%*%u
  
  sigma=sqrt(1/B)*abs(r-rf)
  return(sigma)
}


# 2) Optimisation avec actif sans risque

# Renvoie les poids optima pour un rendement attendu
fopt<-function(Q,X,r,rf) {
  n=length(X)
  e=numeric(n)+1
  v=c()
  m=c()
  
  u=X-rf*e
  v[1]=t(e)%*%solve(Q)%*%u
  v[2]=t(u)%*%solve(Q)%*%u
  
  poidsrisque=((r-rf)*solve(Q)%*%u)/v[2]    #vecteur de poids des actifs risqués
  poidssansrisque=1-t(e)%*%poidsrisque       #poids de l'actif sans risque
  m=c(poidssansrisque,poidsrisque)
  
  return(m)
}


# 3) Tracé de la frontière efficiente (risque=f(rendement)) avec actifs risqués et actif sans risque

fronteffi<-function(Q,X,rf) {
  x=seq(0,1,0.01)
  y1=volat(Q,X,x)
  y2=volat2(Q,X,x,rf)
  plot(y1,x,col="red",type="l",xlab="rendement",ylab="volatilité",main="Frontière efficiente")
  lines(y2,x,type="l",col="blue")
}

#Poids optimaux avec uniquement des actifs risqués
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

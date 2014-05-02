## Etude du coefficient beta issu de la régression linéaire du rendement du portefeuille sur celui du marché
source("C:/Users/Yann/Documents/GitHub/statapp/black_litterman_code_kevin.r")

p=t(fopt2(var(rdt_j_a),colMeans(rdt_j_a),0.012))
rdt_pf=t(p%*%t(rdt_j_a))
rdt_reg=rdt_j[,"cac40"]-mean(rdt_j[,"oat"])
lm.beta=lm(rdt_pf ~ rdt_reg )

summary(lm.beta)
library(xtable)
print(xtable(summary(lm.beta)), type="latex", file="D:/cours/2A/Statap/Mémoire/summary_beta.tex")

plot(rdt_j[,"cac40"],rdt_pf,xlab="R(marché)-r(sans risque)",ylab="R(portefeuille)",main="Régression du MEDAF")
abline(lm.beta,col='red')

# Tracés sur les résidus
par(mfrow=c(2,2))
plot(lm.beta)

# Tests de Breusch-Pagan pour l'hétéroscédasticité
library(lmtest)
bptest(lm.beta)

library(car)
ncvTest(lm.beta)

# on a de l'hétéroscédasticité, on corrige alors la matrice des covariances
library(sandwich)
coeftest(lm.beta, vcov = vcovHC(lm.beta, type = "HC3"))

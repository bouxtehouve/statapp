## Etude du coefficient beta issu de la régression linéaire du rendement du portefeuille sur celui du marché
source("C:/Users/Yann/Documents/GitHub/statapp/black_litterman_code_kevin.r")

p=t(fopt2(sigma(rdt_j_a),colMeans(rdt_j_a),0.012))
rdt_pf=t(p%*%t(rdt_j[,titres_cac40]))

lm.beta=lm(rdt_pf ~ rdt_j[,"cac40"])

summary(lm.beta)

plot(rdt_j[,"cac40"],rdt_pf,xlab="Rendement du marché",ylab="Rendement du portefeuille",main="Régression linéaire du rendement du portefeuille sur celui du marché")
abline(lm.beta,col='red')

# Tracés sur les résidus
par(mfrow=c(2,2))
plot(lm.beta)

# Tests de Breusch-Pagan pour l'hétéroscédasticité
library(lmtest)
bptest(lm.beta)

library(car)
ncvTest(lm.beta)

print(xtable(summary(lm.beta)), type="latex", file="D:/cours/2A/Statap/Mémoire/summary_beta.tex")
## Etude du coefficient beta issu de la r�gression lin�aire du rendement du portefeuille sur celui du march�
source(paste(path,"black_litterman_code_kevin.r",sep="/"))

p=t(fopt2(sigma(rdt_j_a),colMeans(rdt_j_a),0.001169648))
rdt_pf=t(p%*%t(rdt_j_a))

lm.beta=lm(rdt_pf ~ rdt_j[,"cac40"])

summary(lm.beta)

plot(rdt_j[,"cac40"],rdt_pf,xlab="Rendement du march�",ylab="Rendement du portefeuille",main="R�gression lin�aire du rendement du portefeuille sur celui du march�")
abline(lm.beta,col='red')

# Trac�s sur les r�sidus
par(mfrow=c(2,2))
plot(lm.beta)

# Tests de Breusch-Pagan pour l'h�t�rosc�dasticit�
library(lmtest)
bptest(lm.beta)

library(car)
ncvTest(lm.beta)

print(xtable(summary(lm.beta)), type="latex", file="D:/cours/2A/Statap/M�moire/summary_beta.tex")
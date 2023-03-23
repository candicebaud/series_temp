require(zoo)
require(tseries)
library(dplyr)

#base de données
setwd("C:/Users/candi/Desktop/ETUDES/ENSAE2A/semestre 2/séries temporelles/series temp/series_temp")
data <- read.csv('valeurs_mensuelles_pesticides.csv', sep=";")
data <- data[2]

data <- as.data.frame(data[-(1:3),])#on enlève les premières lignes qui ne sont pas des données
indice <- as.data.frame(as.numeric(unlist(data)))

xm.source <- zoo(indice[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm.source)
test <- tail(xm.source, n=2) #pour comparer nos prévisions avec les vraies données
xm <- xm.source[(250):(T-2)] #pour le modèle

mean(xm.source)
plot(xm, xaxt="n") #plot des données
axis(side=1,at=seq(0,400,12)) #pour mettre l'axe x

par(mfrow=c(1,2))
acf(xm)
pacf(xm) #saisonnalité apparente : saisonnalité de 12 donc annuelle, avec été et hiver différenciés (corr pos et neg)
dev.off()

pp.test(xm) #test de philippe perron, on rejette à 1% l'hypothèse que la série n'est pas stationnaire

#on retire la moyenne de xm
xm <- xm - mean(xm)

#on enlève la saisonnalité apparente
xm <- diff(xm, lag = 12)
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #la saisonnalité a bien disparu

#on identifie avec l'acf et la pacf les ordres maximums à tester
pmax = 12
qmax = 11

#modèle
#tous les modèles
selection_print <- function(pmax, qmax){
  res_aic <- matrix(nrow=pmax, ncol=qmax)
  for (p in 1:pmax){
    for (q in 1:qmax){
      model <- arima(xm, c(p,0, q))
      res_aic[p,q] <- model$aic
      print(c(p,q))
    }}
  return(res_aic)
}

selec <- selection_print(pmax,qmax) #matrice qui renvoie les AIC de tous les modèles
min(selec)#l'aic minimal est de 1088.943
which(selec == min(selec),  arr.ind=TRUE) #on choisit l'AIC le plus petit

p = 5
q = 11

#au dessus on n'a pas exploré les modèles avec p = 0 ou q = 0 donc on les explore ci-dessous
selection_bis <- function(pmax, qmax){
  p_0 <- numeric(length = qmax)
  q_0 <- numeric(length = pmax)
  for (q in 1:qmax){
    p_0[q] <- arima(xm, c(0,0,q))$aic
  }
  for (p in 1:pmax){
    q_0[p] <- arima(xm, c(p,0,0))$aic
  }
  return(list(p_0, q_0))
}

res_0 <- selection_bis(12, 11) #résultat de notre exploration rendu sous forme de liste
min(res_0[[1]])#minimum pour p=0
min(res_0[[2]])#minimum pour q=0
#les deux aic minimums obtenus pour lorsque p est fixé à 0 et q est fixé à 0 respectivement sont de 1097.329 et 1097.295


#On choisit donc un modèle ARMA(5,11)
arma_fit <- arima(xm, c(5,0,11))
arma_fit

#résidus
plot(arma_fit$residuals)
acf(arma_fit$residuals)
pacf(arma_fit$residuals)

hist(arma_fit$residuals)
library(forecast)
checkresiduals(arma_fit)

#Q test
#test
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arma_fit$residuals, 24, 5) #tests de LB pour les ordres 1 a 24
#on rejette le fait que les résidus soient corrélés 

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signific <- as.data.frame(signif(arma_fit))

#prévision
model_pred <- predict(arma_fit, n.ahead=2)
serie_pred <- zoo(c(xm, model_pred$pred))

#graphiques
xm_all <- xm.source[250:T] - mean(xm.source[250:(T-2)])
xm_all <- diff(xm_all, lag = 12)

plot(xm_all, col = 'black', ylab = 'Série', main = 'Prévision des 2 prochaines valeurs de la série')
#lines(xm_all, col = 'black', type = 'p') pour avoir des ronds à chaque valeur de la série temporelle
U = model_pred$pred + 1.96*model_pred$se
L = model_pred$pred - 1.96*model_pred$se
xx = c(time (U), rev (time (U)))
yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray (0.6, alpha=0.2))
lines(model_pred$pred, type = "p", col = "red")
lines(model_pred$pred, type = 'l', col = 'red')
legend("topleft", legend=c("Données réelles", "Prédiction"), col=c("red", "black"), lty=1:2, cex=0.4)


#export de la table de significativité des modèles pour le document latex
library(xtable)
xtable(signific)

xtable(signific %>% select(ar1, ar2, ar3, ar4, ar5, ma1, ma2, ma3, ma4, ma5, ma6))

xtable(signific %>% select(ma7, ma8, ma9, ma10, ma11, intercept))

xtable((as.data.frame(selec))%>%select(V1, V2, V3, V4, V5, V6, V7, V8, V9))

xtable((as.data.frame(selec))%>%select(V8, V9, V10, V11))

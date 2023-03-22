require(zoo)
require(tseries)
library(dplyr)

#base de données
data <- read.csv('valeurs_mensuelles_pesticides.csv', sep=";")
data <- data[2]

data <- as.data.frame(data[-(1:3),])#on enlève les premières lignes qui ne sont pas des données
indice <- as.data.frame(as.numeric(unlist(data)))

xm.source <- zoo(indice[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm.source)
test <- tail(xm.source, n=4) #pour comparer nos prévisions avec les vraies données
xm <- xm.source[(250):(T-4)] #pour le modèle

plot(xm, xaxt="n") #plot des données
axis(side=1,at=seq(0,400,12)) #pour mettre l'axe x

par(mfrow=c(1,2))
acf(xm)
pacf(xm) #saisonnalité apparente : saisonnalité de 12 donc annuelle, avec été et hiver différenciés (corr pos et neg)
dev.off()

pp.test(xm) #test de philippe perron, on rejette à 1% l'hypothèse que la série n'est pas stationnaire


#on enlève la saisonnalité apparente
xm <- diff(xm, lag = 12)
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #la saisonnalité a bien disparu sauf sur la pacf


#on retire la moyenne de xm
xm <- xm - mean(xm)
plot(xm, xaxt="n")
axis(side=1,at=seq(0,400,12))
acf(xm, lag = 30) #q peut aller jusqu'à 12??? ou bien on dit qu'on s'arrête en gros à 3?
pacf(xm, lag = 50) #p= 26 ?  ou pareil on s'arrête vers 4-5 : pb à 12 aussi gros pic

pmax = 12
qmax = 11

#modèle
#tous les modèles
selection_print <- function(pmax, qmax){
  res_aic <- matrix(nrow=pmax, ncol=qmax)
  for (p in 0:pmax){
    for (q in 0:qmax){
      model <- arima(xm, c(p,0, q))
      res_aic[p,q] <- model$aic
      print(c(p,q))
    }}
  return(res_aic)
}

selec <- selection_print(pmax,qmax)
which(selec == min(selec),  arr.ind=TRUE)

p = 5
q = 11

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
model_pred <- predict(arma_fit, n.ahead=4)

serie_pred <- zoo(c(xm, model_pred$pred))

xm_all <- xm.source[250:T] - mean(xm.source[250:(T-4)])
xm_all <- diff(xm_all, lag = 12)
plot(xm_all, col = 'red')
lines(serie_pred, col = 'black')

#intervalle de confiance juste pour voir à modif 
IC95_sup <- zoo(model_pred$pred + 1.96*model_pred$se/2)
IC95_low <- zoo(model_pred$pred - 1.96*model_pred$se/2)

lines(IC95_sup, col = 'blue')
lines(IC95_low, col = 'blue')

#calcul rmse  
rmse <- sqrt(sum((model_pred$pred - tail(xm_all, n=4))**2)/4)









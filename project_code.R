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
xm <- xm.source[1:(T-100)] #pour faire de la prévision

plot(xm, xaxt="n") #plot des données
axis(side=1,at=seq(0,300,12)) #pour mettre l'axe x

#xm <- xm - mean(xm) #pour enlever la moyenne
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #saisonnalité apparente : saisonnalité de 12 donc annuelle
dev.off()

#on enlève la saisonnalité apparente
xm <- diff(xm, lag = 12)
plot(xm, xaxt="n")
axis(side=1,at=seq(0,300,12))
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #la saisonnalité a bien disparu sauf sur la pacf

pp.test(xm) #test de philippe perron, on rejette à 1% l'hypothèse que la série n'est pas stationnaire

#au niveau acf : on voit qu'on rentre environ dans l'itv de confiance pour l'acf au lag 5 environ (donc on prend q = 4)
#au niveau pacf : on rentre complètement à p = 2 mais on ressort à p = 12
#dans les deux cas on ressort à p=12 même après différenciation ... 

#deuxième saisonnalité?
#xm_bis <- diff(xm, lag = 12)
#plot(xm_bis)
#par(mfrow=c(1,2))
#acf(xm_bis)
#pacf(xm_bis)

#on retire la moyenne de xm
xm <- xm - mean(xm)
plot(xm)
acf(xm)
pacf(xm)

arima_1 <- arima(xm, c(2,0,5))
arima_1

Box.test(arima_1$residuals, lag=6, type="Ljung-Box", fitdf=5)

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima_1$residuals, 24, 5) #tests de LB pour les ordres 1 a 24
round(Qtests(arima_1$residuals,24,fitdf=5),3)

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arima_1)






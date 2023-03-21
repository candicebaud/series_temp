setwd("C:/Users/candi/Desktop/ETUDES/ENSAE2A/semestre 2/séries temporelles/series temp")

require(zoo)
require(tseries)
library(dplyr)

data <- read.csv('valeurs_mensuelles.csv', sep=";")
data <- data[2]


data <- as.data.frame(data[-(1:3),])#on enlève les premières lignes qui ne sont pas des données
indice <- as.data.frame(as.numeric(unlist(data)))

xm.source <- zoo(indice[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm.source)
xm <- xm.source[1:(T-4)]

#plot(xm, ylim = 100)
plot(xm, xaxt="n", ylim= 100) #plot des données
axis(side=1,at=seq(0,244,12)) #mettre les axes

xm <- xm - mean(xm) #pour enlever la moyenne
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #pas de saisonnalité apparente
dev.off()

pp.test(xm) #test de philippe perron ok, on rejette à 1% l'hypothèse que la série n'est pas stationnaire


#au niveau acf : on voit qu'on rentre complètement dans l'itv de confiance pour l'acf au lag 17 (donc on prend q =16)
#au niveau pacf : on rentre complètement à p = 4

#donc il faudrait prendre un ARIMA(4, 0, 17)

arima4017 <- arima(xm, c(4,0,17))

Box.test(arima4017$residuals, lag=6, type="Ljung-Box", fitdf=5) #

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima4017$residuals, 24, 5) #tests de LB pour les ordres 1 a 24
round(Qtests(arima4017$residuals,24,fitdf=5),3)

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arima4017)









#avec une autre base ? 
hydro <- read.csv('valeurs_mensuelles_hydrocarbures.csv', sep = ';')
hydro <- hydro[2]

hydro <- as.data.frame(data[-(1:3),])#on enlève les premières lignes qui ne sont pas des données
indice_hydro <- as.data.frame(as.numeric(unlist(hydro)))

serie_hydro <- zoo(indice_hydro[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(serie_hydro)
serie_hydro <- serie_hydro[1:(T-4)]

plot(serie_hydro)
acf(serie_hydro)
pacf(serie_hydro)

pp.test(serie_hydro)



#encore une autre base ? 
elec <- read.csv('valeurs_mensuelles_elec.csv', sep = ';')
elec <- elec[2]

elec <- as.data.frame(data[-(1:3),])#on enlève les premières lignes qui ne sont pas des données
indice_elec <- as.data.frame(as.numeric(unlist(hydro)))

serie_elec <- zoo(indice_elec[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(serie_elec)
serie_elec <- serie_elec[1:(T-4)]

plot(serie_elec)
acf(serie_elec)
pacf(serie_elec)

pp.test(serie_elec)


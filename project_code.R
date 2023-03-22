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
xm <- xm.source[1:(T-4)] #pour le modèle

plot(xm, xaxt="n") #plot des données
axis(side=1,at=seq(0,300,12)) #pour mettre l'axe x

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

#test
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


#sélection du modèle ? 
library(weakARMA)
pmax=12
qmax=5
selec <- ARMA.selec(xm, 1, 1)

#tester tous les modèles à la main
selection <- function(pmax, qmax){
  res <- c(NA, NA)
  AIC <- -Inf
  for (p in 0:pmax){
    for (q in 0:qmax){
      model <- arima(xm, c(p,0, q))
      if(model$aic>=AIC){
        AIC <- model$aic
        res[1] = p
        res[2] = q
      }
    }
  }
  return(res)
}

#ne marche pas lol
selection(pmax, qmax)


#tous les modèles
selection_print <- function(pmax, qmax){
  res_aic <- matrix(nrow=pmax, ncol=qmax)
  for (p in 0:pmax){
    for (q in 0:qmax){
      model <- arima(xm, c(p,0, q))
      res_aic[p,q] <- model$aic
      }}
  return(res)
}

selec <- selection_print(pmax,qmax)

#le modèle sélectionné est celui avec p=12 et q=5... 

arima_000 <- arima(xm, c(0,0,0))
arima_205 <- arima(xm, c(2,0,5))
arima_105 <- arima(xm, c(1,0,5))

arima_selec <- arima(xm, c(12,0,5))
Box.test(arima_selec$residuals, lag=6, type="Ljung-Box", fitdf=5)
#pas rejetté donc ok 

#forecast 
model_pred <- predict(arima_selec, n.ahead=4)
ts.plot(xm, model_pred$pred, ylim = c(0, 200), 
        ylab = "CPI", col = "blue")

serie_pred <- zoo(c(xm, model_pred$pred))

xm_all <- xm.source - mean(xm.source[1:(T-4)])
xm_all <- diff(xm_all, lag = 12)
plot(xm_all)
lines(serie_pred, col = 'red')



















#sélection du modèle
## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=xm, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}


## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}



#on fait tous les modèles
armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima (patienter...)

#on sélectionne ceux qui sont bien
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec

#on évalue avec le AIC, BIC et le RMSE sur la prédiction
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cree une liste des ordres p et q des modeles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les elements de la liste
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]]))) #cree une liste des modeles candidats estimes
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et BIC des modeles candidats
### L'ARMA(2,1) minimise les criteres d'information.

trend <- 1:length(xm)
lt <- lm(xm ~ trend) 

rps <- lapply(models, function(m) as.zoo(predict(m,4)$pred)) #previsions de r
#xmps <- lapply(rps, function(rp) rp+cbind(1,c((T-3):T))%*%lt$coefficients) #previsions de xm
#rmse <- vapply(xmps, FUN.VALUE=numeric(1), function(xmp) sqrt(sum((as.zoo(xmp)-tail(xm.source,4))^2))) #calcule les rmse out-of-sample











#prévision
library(forecast)
f <- forecast(arima_1,h=100)
plot(xm, xlim=c(290, 350), ylim=c(min(xm), max(xm, f$upper)), xlab="Années", ylab="Valeur", main="Prévisions avec intervalle de confiance à 95%")
lines(f$mean, col="red")
lines(f$upper, col="blue", lty=2)
lines(f$lower, col="blue", lty=2)
points(xm, col="green")
legend("topright", legend=c("Données réelles", "Prévisions", "Intervalles de confiance"), col=c("green", "red", "blue"), lty=c(1, 1, 2), cex=0.8)

#on fait les mêmes opérations à notre série test de comparaison
test <- test - mean(xm)
test <- diff(test, lag=12)

#pb on ne voit pas les itvl et la prediction est naze
plot(xm, xlim=c(250,350))
lines(test, col = 'blue')
lines(f$mean, col='red')


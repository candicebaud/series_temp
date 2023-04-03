require(zoo)
require(tseries)
library(dplyr)

#fonctions utilisées dans le code
signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}

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

armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}


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

#mean(xm.source)
plot(xm, xaxt="n") #plot des données
axis(side=1,at=seq(0,400,12)) #pour mettre l'axe x

par(mfrow=c(1,2))
acf(xm)
pacf(xm) #saisonnalité apparente : saisonnalité de 12 donc annuelle, avec été et hiver différenciés (corr pos et neg)
dev.off()

#on retire la moyenne de xm
xm <- xm - mean(xm)

#on enlève la saisonnalité apparente
xm <- diff(xm, lag = 12)
par(mfrow=c(1,2))
acf(xm)
pacf(xm) #la saisonnalité a bien disparu

pp.test(xm) #test de philippe perron, on rejette à 1% l'hypothèse que la série n'est pas stationnaire

#on identifie avec l'acf et la pacf les ordres maximums à tester
pmax = 12
qmax = 11

#modèle
#estimation de tous les modèles et sélection des modèles valides
arma_valid <- armamodelchoice(12,11)
selec <- arma_valid[arma_valid[,"ok"]==1&!is.na(arma_valid[,"ok"]),]

#les modèles possibles sont donnés par
selec
#on peut donc choisir p=12, q=2 ou p=12, q=9, ou p=5, q=11

#on fit les trois modèles et on calcule les aic
arma_12_2 <- arima(xm, c(12,0,2))
arma_12_9 <- arima(xm, c(12,0,9))
arma_5_11 <- arima(xm, c(5,0,11))

#aic
arma_12_2$aic
arma_12_9$aic
arma_5_11$aic #modèle ayant le plus petit aic

#bic
BIC(arma_12_2)
BIC(arma_12_9)
BIC(arma_5_11) #modèle ayant le plus petit bic


#valeurs sélectionnées pour notre modèle
p = 5
q = 11

#fit du modèle
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
Qtests(arma_fit$residuals, 24, 5) #tests de LB pour les ordres 1 a 24
#on rejette le fait que les résidus soient corrélés 

signific <- as.data.frame(signif(arma_fit))

#causalité
roots <- polyroot(sort(arma_fit$coef[c('ar1', 'ar2', 'ar3', 'ar4', 'ar5')]))
modulus_roots <- Mod(roots)
modulus_roots #les coefficients sont bien plus grands que 1 donc le modèle est causal

#prévision
model_pred <- predict(arma_fit, n.ahead=2)
serie_pred <- zoo(c(xm, model_pred$pred))

#graphiques
xm_all <- xm.source[250:T] - mean(xm.source[250:(T-2)])
xm_all <- diff(xm_all, lag = 12)

dev.off()

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


#bonus: table avec les aic et bic de tous les modèles 
selection_print_aic <- function(pmax, qmax){
  res_aic <- matrix(nrow=pmax, ncol=qmax)
  for (p in 1:pmax){
    for (q in 1:qmax){
      model <- arima(xm, c(p,0, q))
      res_aic[p,q] <- model$aic
      print(c(p,q))
    }}
  return(res_aic)
}
selec_aic <- selection_print_aic(pmax,qmax)

xtable((as.data.frame(selec_aic))%>%select(V1, V2, V3, V4, V5, V6, V7, V8, V9))

xtable((as.data.frame(selec_aic))%>%select(V8, V9, V10, V11))

selection_print_bic <- function(pmax, qmax){
  res_bic <- matrix(nrow=pmax, ncol=qmax)
  for (p in 1:pmax){
    for (q in 1:qmax){
      model <- arima(xm, c(p,0, q))
      res_bic[p,q] <- BIC(model)
      print(c(p,q))
    }}
  return(res_bic)
}
selec_bic <- selection_print_bic(pmax,qmax)

xtable((as.data.frame(selec_bic))%>%select(V1, V2, V3, V4, V5, V6, V7, V8, V9))

xtable((as.data.frame(selec_bic))%>%select(V8, V9, V10, V11))



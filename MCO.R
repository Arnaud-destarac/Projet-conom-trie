
#                 REGRESSION LINEAIRE MULTIPLE

rm(list=ls())  # Effacer les donnees en memoire
graphics.off() # Fermer les graphiques
library(readxl)

# Lecture des données

données <- read_excel('données.xlsx')
attach(données)

Y = Celec_menages/Pop1
X1 = (PIB2020/`IPC(base100=2015)`)/Pop1
X2 = Pelec/`IPC(base100=2015)`
X3 = DJU

n=length(Y)

# Choix des variables 

vec <- c(X1,X2,X3)
X <- matrix( vec, ncol=3) 
Y=matrix(Y,n,1)
k=ncol(X);
K=k+1

# Transformation des variables en logarithme

y=log(Y) 
x=log(X)

# Estimation MCO

OLS=lm(formula = y ~ x)

# Affichage des résultats

summary(OLS)

xc = cbind(1,x) 
bhat = OLS$coefficients 
yf = xc %*% bhat
res = y - yf
scr = t(res) %*% res

# Graphique de la serie Y(i) en fonction de i
nobs=cbind(1:n)

dev.new()
xbnd=c(1*range(nobs)[1], 1*range(nobs)[2])
ybnd=c(0.9*range(y)[1], 1.1*range(y)[2])
plot(nobs,y,xlab=" ",ylab="serie Y",col="green",xlim=xbnd,ylim=ybnd,type="l")
title(main="Observed values")


## NORMALITE DES ERREURS ?

#QQ-Plot
qqnorm(y=res, xlab='Quantiles loi normales', ylab='Quantiles échantillon')

# Test de Kolomorov-Smirnov
ks.test(x=res,y="pnorm",mean=0, sd=5.905)

# Histogramme des résidus 

dev.new()
hist(residuals(OLS))


## AUTOCORRELATION ?

# Statistique de Durbin-Watson

d1 = scr
d2 =  t(res[2:n]-res[1:n-1]) %*% (res[2:n]-res[1:n-1])
dw = d2/scr
print (dw)

# Test de Ljung-Box

#Pour savoir le nombre de retards (lags) indiqués dans le test Ljung-Box
pacf(res)
dev.print(device= jpeg, file="pacf.jpeg", width=600)
acf(res)
dev.print(device= jpeg, file="acf.jpeg", width=600)
#
# lag : déterminé à partir des graphes de PACF et ACF, ici 2
# fitdf : K, ici 4

Box.test(res, lag = 4, type = c( "Ljung-Box"), fitdf = 4)


## HOMOSCEDASTICITE ?

# Regression auxiliaire de White 

e2=res*res 
xcarre=x*x 
x1x2=X1*X2 
x1x3=X1*X3 
x2x3=X2*X3 
WAUX=lm(formula = e2 ~ x+xcarre+x1x2+x1x3+x2x3) 
summary(WAUX) 
resw = WAUX$residuals
SCRw = t(resw) %*% resw
Rw = 1 - SCRw/(var(e2)*(n-1))
Rw_sq = Rw^2
print(n*Rw_sq)
# la valeur obtenue est inférieure à la valeur seuil chi-deux(9) pour alpha = 0.05 (=16.919), donc on ne rejette pas l'hyp. d'homoscédasticité

# Graphique des résidus
dev.new()
ybnd=c(0.9*range(res)[1], 1.1*range(res)[2])
plot(nobs,res,xlab=" ",ylab="residuals",col="blue",xlim=xbnd,ylim=ybnd,type="l")
title(main="Residuals")


## MULTICOLINEARITE ?

# Coefficients de correlation entre les variables x
cc<-cor(x)
print(cc)

# calcul du VIF

cat("Critere du VIF","\n")
VIF <-matrix(0,k,1)
xc=cbind(x,1)
xtx1=solve(t(xc)%*%xc)
for (j in 1:k) {
  VIF[j]= var(x[1:n,j])*(n-1)*xtx1[j,j]
  cat(j,VIF[j],"\n")
}

## STABILITE TEMPORELLE

# Test de Chow

library(strucchange)

# Chow (en faisant l'hypothèse d'une rupture à l'observation "point")
sctest(y ~ x, type = "Chow", point = )

# Chow pas à pas car pas de rupture observée sur le graphe des résidus 
# Rq: le test commence à 5 car il faut plus d'obersvations que de variables dans chaque sous-période
for(i in 5:20) {
  print(sctest(y ~ x, type = "Chow", point = i) )
}
# pour i = 11, la stats de Fischer est la plus élevée

# Test Cusum

Wr <- efp(y ~ x, type = "Rec-CUSUM")
plot(Wr)
# résultat inintéressant

# Test Cusum Square

rr <- (recresid(y ~ x))
rr <- rr^2
cumrr <- cumsum(rr)/scr

# Valeurs seuil de la distribution Cusum

c0 = 0.197 # cf table avec ici n-K = 28
Kp1=K+1

t2 <- ts(Kp1:n)
t2 = c(5:32)

smin <-((t2-K)/(n-K))-c0
smax <- ((t2-K)/(n-K))+c0
#

vec2 <- c(smin, cumrr, smax)
cusum2 <- matrix(vec2, ncol = 3); 
matplot(t2, cusum2, type ="l")
# même interprétation que pour le test de chow, rupture pour i = 11
# on refait le test pour la régression de 2000 (i=11) à 2021

#régression
y_2 = y[11:n]
x_2 = x[11:n]
n_2 = length(y_2)
OLS_2=lm(formula = y_2 ~ x_2)
summary(OLS_2)

xc_2 = cbind(1,x_2) 
bhat_2 = OLS_2$coefficients 
yf_2 = xc_2 %*% bhat_2
res_2 = y_2 - yf_2
scr_2 = t(res_2) %*% res_2

# Test Cusum Square

rr_2 <- (recresid(y_2 ~ x_2))
rr_2 <- rr_2^2
cumrr_2 <- cumsum(rr_2)/scr_2

# Valeurs seuil de la distribution Cusum

c0_2 = 0.23298 # cf table avec ici n-K = 18
Kp1=K+1

t2_2 = c(11:32)

smin_2 <-((t2_2-K)/(n_2-K))-c0_2
smax_2 <- ((t2_2-K)/(n_2-K))+c0_2
#

vec2_2 <- c(smin_2, cumrr_2, smax_2)
cusum2_2 <- matrix(vec2_2, ncol = 3); 
matplot(t2_2, cusum2_2, type ="l")
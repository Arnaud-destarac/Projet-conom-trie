
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
x = x[1:n,1:3]

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

##implémenter méthode de Cochrane-Orcutt pour obtenir les coeffs sans autocorrélation

library(orcutt)
coch = cochrane.orcutt(OLS)
coch
summary(coch)

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

# Chow pas à pas car pas de rupture observée sur le graphe des résidus 
# Rq: le test commence à 5 et finit à n-K-1=27 car il faut plus d'obersvations que de variables dans chaque sous-période
for(i in 5:(n-K-1)) {
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

t2p <- ts(Kp1:n)
t2 = c(5:32)

smin <-((t2-K)/(n-K))-c0
smax <- ((t2-K)/(n-K))+c0
#

vec2 <- c(smin, cumrr, smax)
cusum2 <- matrix(vec2, ncol = 3); 
matplot(t2, cusum2, type ="l")
# rupture pour i = 10
# on refait le test pour la régression de 1999 (i=10) à 2021


## Test de rupture sur le 2è sous-échantillon

rupture = 10
#régression
y_2 = y[rupture:n]
x_2 = x[rupture:n,1:3]
n_new = length(y_2)

OLS_2=lm(formula = y_2 ~ x_2)
summary(OLS_2)

xc_2 = cbind(1,x_2) 
xt_2 = t(xc_2)
bmco_2 = OLS_2$coefficients 
ycalc_2 = xc_2 %*% bmco_2
xtx_2= xt_2 %*% xc_2
xtx1_2 = solve(xtx_2) 
u_2=y_2-xc_2%*%bmco_2
scr_2 = t(u_2) %*% u_2
#Test de Chow
for(i in 5:(n_new-K-1)) {
  print(sctest(y_2 ~ x_2, type = "Chow", point = i) )
}
# Pas d'autres ruptures observées

#Test de Cusum Square
rr_2 <- (recresid(y_2 ~ x_2))
rr_2 <- rr_2^2
cumrr_2 <- cumsum(rr_2)/scr_2

c0 = 0.23298 # cf table avec ici n_new-K = 19
Kp1=K+1

t2_2 <- ts(Kp1:n)
t2_2 = c(Kp1:n_new)

smin_2 <-((t2_2-K)/(n_new-K))-c0
smax_2 <- ((t2_2-K)/(n_new-K))+c0
#

vec2_2 <- c(smin_2, cumrr_2, smax_2)
cusum2_2 <- matrix(vec2_2, ncol = 3); 
matplot(t2_2+rupture, cusum2_2, type ="l")



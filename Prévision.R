
rm(list=ls())  # Effacer les donnees en memoire
graphics.off() # Fermer les graphiques
library(readxl)

# Lecture des données

données <- read_excel('données.xlsx')
attach(données)

n_y = 32

Y = Celec_menages[1:n_y]/Pop1[1:n_y]
X1 = (PIB2020/`IPC(base100=2015)`)/Pop1
X2 = Pelec/`IPC(base100=2015)`
X3 = DJU

n_x = length(X1)

# Choix des variables 

vec <- c(X1,X2,X3)
X <- matrix( vec, ncol=3) 
Y=matrix(Y,n_y,1)
k=ncol(X);
K=k+1

# Transformation des variables en logarithme

y=log(Y) 
x=log(X)

#on ne garde que les données post-rupture
rupture = 10
y_2 = y[rupture:n_y]
x_2 = x[rupture:n_y,1:3]
n = n_y - rupture + 1

# Estimation MCO

OLS=lm(formula = y_2 ~ x_2)

# Affichage des résultats

summary(OLS)

# On enlève l'AR
#library(orcutt)
#coch = cochrane.orcutt(OLS)
#coch
#summary(coch)

xc = cbind(1,x_2) 
bmco = OLS$coefficients
# bmco = coch$coefficients
xt = t(xc)
xtx= xt %*% xc
xtx1 = solve(xtx)
yf = xc %*% bmco
u = y_2 - yf
scr = t(u) %*% u

# histogramme des résidus MCO
dev.new()
hist(residuals(OLS))
d1 <- density(u)
plot(d1, main="Histogramme : noyau des résidus MCO")


## Prévision

nprev = n_x - n_y # =9, de 2022 à 2030

for( j in 1:nprev){
  xf<-x[n_y+j,1:3]
  nvx <-matrix(c(1,xf))
nvx1 = t(nvx);
nvy= nvx1 %*% bmco 
#
# Ecart-type de prevision
#
sprev= sqrt( (scr/(n-K)) * (1+ nvx1 %*% xtx1 %*% nvx ))
print(sprev)
print(xf)

# Prévision Bootstrap méthode percentile-t et Bootstrap des résidus centrés réduits
#
utild=u
#
# mettre des # devant les lignes suivantes jusqu'? "d?but du boostrap"
# pour ne pas normaliser les erreurs
#
h=xc %*% xtx1 %*% xt
for (i in 1:n) {
  utild[i] = u[i]/sqrt(1-h[i,i])
}
utildmoy=0.
for (i in 1:n) {
  utildmoy = utildmoy+utild[i]
}
utildmoy=utildmoy/n
utild=utild-utildmoy
#
# debut du bootstrap
#
# nombre de replications
#
m=1000;
nvybb<-matrix(0,m,1)
zb<-matrix(0,m,1)
#
# Debut de la boucle de bootstrap
#
for (b in 1:m) {
  #
  indboot=sample(1:n,n,replace = T)
  ub=utild
  for (i in 1:n) {
    ub[i] = utild[indboot[i]]
  }
  yb = xc %*% bmco+ub
  xtyb = xt %*% yb
  bb = xtx1 %*% xtyb
  yhatb = xc%*%bb 
  scrb = t(yb - yhatb) %*% (yb - yhatb) 
  s2b = scrb / (n-K) 
  sb = sqrt (s2b)
  #
  # Prévision bootstrap de Y
  #
  nvyb=nvx1 %*% bb
  pdt = 1 + nvx1 %*% xtx1 %*% nvx
  sprevb = sb*sqrt(pdt)
  # construction de Z*(i)
  indboot2=sample(1:n,1,replace = T)
  uprevb=utild[indboot2]
  z = (nvyb - nvy-uprevb)/ sprevb 
  #
  nvybb[b]=nvyb
  zb[b]=z
}
#
# fin de la boucle sur b
#
zitri = sort (zb)
#
# histogramme zitri
#

dev.new()
d2 <- density(zitri)
plot(d2, main= "Histogramme de z*")

i025=floor(0.025*m)
i05=floor(0.05*m)
i95=floor(0.95*m)
i975=floor(0.975*m)
#construction de l'IC bootstrapt
inft025 = nvy - zitri[i975]*sprev
inft05 = nvy - zitri[i95]*sprev
supt95 = nvy - zitri[i05]*sprev
supt975 = nvy - zitri[i025]*sprev

#
# Intervalle de prédiction MCO - Borne inf et sup à 95%

# T de Student 
#
tstud=qt(0.975,n-K)
#
prevymin=nvy-tstud*sprev
prevymax=nvy+tstud*sprev

# Affichage 

cat("Intervalles de confiance",j, ":",exp(nvy), exp(prevymin), exp(prevymax), exp(inft025), exp(supt975),sprev, "\n")

}


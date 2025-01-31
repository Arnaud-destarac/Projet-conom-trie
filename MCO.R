
#                 REGRESSION LINEAIRE MULTIPLE

rm(list=ls())  # Effacer les donnees en memoire
graphics.off() # Fermer les graphiques

library(ggplot2)
library(car)          # Pour tests de multicolinéarité
library(lmtest)       # Pour tests d'hypothèses
library(sandwich)     # Pour erreurs robustes
library(forecast)     # Pour les projections
library(readxl)

######################################################
####"#---------- PRESENTATION DU MODELE ---------#####
######################################################

# Lecture des données

donnees <- read_excel('données-conso_elec.xlsx')

# Exploration des données
str(donnees)
summary(donnees)

# Visualisation des tendances historiques
ggplot(donnees, aes(x = Date, y = Celec_menages/Pop1*10**6)) +
  geom_line() +
  labs(title = "Évolution de la consommation électrique agrégée des ménages français", x = "Année", y = "Consommation (en GWh)")

# Operations sur les variables
Y <- (donnees$Celec_menages * 10^3) / donnees$Pop1  # Consommation en MWh par habitant
X1 <- (donnees$PIB2020 * 10^9) / donnees$Pop1      # PIB réel en euros par habitantBASE 2020
X2 <- donnees$Pelec / (donnees$`IPC(base100=2015)` / 100)  # Prix électricité corrigé en euro / MWh 
X3 <- donnees$DJU                               # Indice climatique

n <- length(Y)  # Nombre d'observations
            
# Choix des variables 
vec <- c(X1, X2, X3)
X <- matrix(vec, ncol = 3, byrow = FALSE)  # Matrice des variables explicatives
Y <- matrix(Y, nrow = n, ncol = 1)   # Vecteur de la variable dépendante

# Transformation en logarithme
y <- log(Y)
x <- log(X)
x <- x[1:n, 1:3]  
k=ncol(X)
K=k+1

# Estimation du modèle MCO
OLS <- lm(formula = y ~ x)
OLS1 <- lm(log(Y) ~ log(X1) + log(X2) + log(X3))

# Affichage des résultats du modèle de regression

summary(OLS)

xc = cbind(1,x) 
bhat = OLS$coefficients 
yf = xc %*% bhat
res = y - yf                        #=residuals(OLS)   plus simple
scr <- sum(residuals(OLS1)^2)

# Graphique de la serie Y(i) et des X(i) en fonction de i(années)
nobs=cbind(1:n)
années = 1990:2021

dev.new()
xbnd=c(1*range(nobs)[1], 1*range(nobs)[2])
ybnd=c(0.9*range(y)[1], 1.1*range(y)[2])
plot(années,y,xlab="années",ylab="log(C_elec/Pop1)",col="green",xlim=range(années),ylim=ybnd,type="p", lwd=2)
lines(années, yf, col="blue", lty=2, lwd=2, type='o')  # Ajout des valeurs prédites
legend("topright", legend=c("valeurs observées", "prédiction"), col=c("green", "blue"), lty=c(1, 2), lwd=2)
title(main="Valeurs observées et prédites")

ybnd=c(0.9*range(log(X1))[1], 1.1*range(log(X1))[2])
plot(années,log(X1),xlab="années",ylab="log(PIB/Pop1)",col="red",xlim=range(années),ylim=ybnd,type="p", lwd=2)

ybnd=c(0.9*range(X2)[1], 1.1*range(X2)[2])
plot(années,X2,xlab="années",ylab="Pelec (€/MWh - base 2015)",col="red",xlim=range(années),ylim=ybnd,type="p", lwd=2)
# On voit clairement une rupture --> c'est sur cette variable qu'on va créer une variable muette

ybnd=c(0.9*range(X3)[1], 1.1*range(X3)[2])
plot(années,X3,xlab="années",ylab="DJU",col="red",xlim=range(années),ylim=ybnd,type="p", lwd=2)

######################################################
#####---------- TESTS SUR LE MODELE ---------#####
######################################################

##### MULTICOLINEARITE ? (H2) -> CHAPITRE III #####
library(corrplot)

# Calcul de la matrice de corrélation (log des variables explicatives)
correlation_matrix <- cor(x)
correlation_matrix

library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 90)
## coefficients de correlation entre les variables x peu eleves en valeur absolue, c'est bon signe 

# Calcul du VIF (Variance Inflation Factor) : 1ere methode
vif_results <- vif(OLS1)
print(vif_results)
## VIFs peu éleves donc c'est bon 

# Calcul du VIF : 2eme méhode

cat("Critere du VIF","\n")
VIF <-matrix(0,k,1)
xc=cbind(x,1)
xtx1=solve(t(xc)%*%xc)
for (j in 1:k) {
  VIF[j]= var(x[1:n,j])*(n-1)*xtx1[j,j]
  cat(j,VIF[j],"\n")
}


##### SUR LES ERREURS (NORMALITE avec espe nulle) ? -> CHAPITRE II #####


#ols_plot_resid_fit(OLS1) #residuals vs fitted values plot -> indique ok 

# Graphique des résidus
dev.new()
ybnd=c(0.9*range(residuals(OLS1))[1], 1.1*range(residuals(OLS1))[2])
plot(nobs,residuals(OLS1),xlab=" ",ylab="residuals",col="blue",xlim=xbnd,ylim=ybnd,type="l")
title(main="Residuals")

library(olsrr)
ols_test_normality(OLS1)

#QQ-Plot
qqnorm(y=residuals(OLS1), xlab='Quantiles loi normales', ylab='Quantiles échantillon', main="QQ plot des résidus")
qqline(residuals(OLS1), col = "red", lwd = 2)  # Ajout d'une ligne de référence
## La majorité des points suivent bien la ligne de référence (ligne rouge). Cela suggère que les résidus se rapprochent d'une distribution normale.

# Histogramme des résidus 
dev.new()
hist(residuals(OLS1))


##### SUR L' AUTOCORRELATION ->  #####

# Statistique de Durbin-Watson

d1 = scr
d2 =  t(res[2:n]-res[1:n-1]) %*% (res[2:n]-res[1:n-1])
dw = d2/scr
print (dw)


#Pour savoir le nombre de retards (lags) indiqués dans le test Ljung-Box
pacf(res, main="pacf des résidus")
dev.print(device= jpeg, file="pacf.jpeg", width=600)
acf(res, main="acf des résidus")
dev.print(device= jpeg, file="acf.jpeg", width=600)

# Test de Ljung-Box
# lag : déterminé à partir des graphes de PACF et ACF
# fitdf : K, ici 4 --> Ljung-Box pas possible
#Box.test(res, lag = 4, type = c( "Ljung-Box"), fitdf = 3)

# ou tester avec Breusch-Godfrey avec un lag de 1, 2, .... 
bgtest(OLS1, order = 1) 
## p value élevée suggère absence d'autocorrelation 


##### SUR HOMOSCEDASTICITE ? #####

# Graphique des résidus
dev.new()
ybnd=c(0.9*range(residuals(OLS1))[1], 1.1*range(residuals(OLS1))[2])
plot(nobs,residuals(OLS1),xlab=" ",ylab="residuals",col="blue",xlim=xbnd,ylim=ybnd,type="l")
title(main="Residuals")
## pas de dispersion croissante : bon signe 

#test de Breusch Pagan
?bptest
bptest(OLS1)

# Regression auxiliaire de White 
e2=res*res 
xcarre=x*x 
x1x2=X1*X2 
x1x3=X1*X3 
x2x3=X2*X3 
WAUX=lm(formula = e2 ~ x+xcarre+x1x2+x1x3+x2x3) 
summary(WAUX) 

 # methode 2 de white avec bibliotheque
library("skedastic")
skedastic_package_white  <- white(mainlm =  OLS1, interactions = TRUE)
skedastic_package_white
# on obtient p=0,342, on peut donc pas rejeter H0 : absence d'heteroscedasticité

#observation visuelle de l'homoscédasticité
plot(fitted(OLS1), rstandard(OLS1), 
     xlab = "Valeurs ajustées", ylab = "Résidus standardisés")
abline(h = 0, col = "red")


## STABILITE TEMPORELLE

# Test de Chow

library(strucchange)

# Chow pas à pas car pas de rupture observée sur le graphe des résidus 
# Rq: le test commence à 5 et finit à n-K-1=27 car il faut plus d'obersvations que de variables dans chaque sous-période
for(i in 5:(n-K-1)) {
  print(sctest(y ~ x, type = "Chow", point = i) )
}

# pour i = 19, la stats de Fischer est la plus élevée


## Test Cusum Square

rr <- (recresid(y ~ x))
rr <- rr^2
cumrr <- cumsum(rr)/scr

# Valeurs seuil de la distribution Cusum

c0 = 0.197 # cf table avec ici n-K = 28
Kp1=K+1

#t2 <- ts(Kp1:n)
t2 = c(5:32)

smin <-((t2-K)/(n-K))-c0
smax <- ((t2-K)/(n-K))+c0
#

vec2 <- c(smin, cumrr, smax)
cusum2 <- matrix(vec2, ncol = 3); 
matplot(t2+1989, cusum2, xlab ="années",  type ="l")
# rupture pour i = 19
# on refait le test pour la régression de 2008 (i=19) à 2021


## Test de rupture sur le 2è sous-échantillon

rupture = 19
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

#Test de Cusum Square
rr_2 <- recresid(y_2 ~ x_2)
rr_2 <- rr_2^2
cumrr_2 <- cumsum(rr_2)/scr_2

c0 = 0.30221 # cf table avec ici n_new-K = 10
Kp1=K+1

t2_2 <- ts(Kp1:n)
t2_2 = c(Kp1:n_new)

smin_2 <-((t2_2-K)/(n_new-K))-c0
smax_2 <- ((t2_2-K)/(n_new-K))+c0
#

vec2_2 <- c(smin_2, cumrr_2, smax_2)
cusum2_2 <- matrix(vec2_2, ncol = 3); 
matplot(t2_2+2007, cusum2_2, xlab ="années",  type ="l")

# nouvelle rupture pour i = 25, ie 2014

rupture = 25
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


#Test de Cusum Square
rr_2 <- (recresid(y_2 ~ x_2))
rr_2 <- rr_2^2
cumrr_2 <- cumsum(rr_2)/scr_2

c0 = 0.37359 # cf table avec ici n_new-K = 4
Kp1=K+1

t2_2 <- ts(Kp1:n)
t2_2 = c(Kp1:n_new)

smin_2 <-((t2_2-K)/(n_new-K))-c0
smax_2 <- ((t2_2-K)/(n_new-K))+c0

vec2_2 <- c(smin_2, cumrr_2, smax_2)
cusum2_3 <- matrix(vec2_2, ncol = 3); 
matplot(t2_2+2013, cusum2_3, xlab ="années",  type ="l")
# pas d'autre rupture


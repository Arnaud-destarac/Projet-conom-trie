
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

# Graphique de la serie Y(i) en fonction de i(années)
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
## coefficients de correlation entre les variables x peu eleves en valeur absolue ie inférieurs à 1, c'est bon signe 

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
## ouais ok 

#test de shapiro 
#shapiro.test(residuals(OLS1))
# p value de ~0.32 => on rejette pas H0 : les residus suivent une loi normale 

# Test de Kolomorov-Smirnov
#?ks.test
#ks.test(x=res,y="pnorm",mean=0, sd=5.905)
# blabla


##### SUR L' AUTOCORRELATION ->  #####

# Statistique de Durbin-Watson

d1 = scr
d2 =  t(res[2:n]-res[1:n-1]) %*% (res[2:n]-res[1:n-1])
dw = d2/scr
print (dw)
# on obtient d=0.987 <= dU  ie il y'a autocorrelation positive ...
# ou tester avec Ljung-box avec un lag de 1, 2, .... 
#?Box.test
#Box.test(residuals(OLS), lag = 1, type = "Ljung-Box")
## p value élevée suggère absence d'autocorrelation 

# ou tester avec Breusch-Godfrey avec un lag de 1, 2, .... 
bgtest(OLS1, order = 1) 
## p value élevée suggère absence d'autocorrelation 

# Test de Ljung-Box

#Pour savoir le nombre de retards (lags) indiqués dans le test Ljung-Box
pacf(res, main="pacf des résidus")
dev.print(device= jpeg, file="pacf.jpeg", width=600)
acf(res, main="acf des résidus")
dev.print(device= jpeg, file="acf.jpeg", width=600)
# lag : déterminé à partir des graphes de PACF et ACF, ici 2! => modele AR(2) ? 
# fitdf : K, ici 4
Box.test(res, lag = 4, type = c( "Ljung-Box"), fitdf = 3)


##implémenter méthode de Cochrane-Orcutt pour obtenir les coeffs sans autocorrélation

    # fonction intégrée
#library(orcutt)
#OLS_corrige = cochrane.orcutt(OLS)
#OLS_corrige
#summary(OLS_corrige)
#avec ce test, on a bien une valeur de d=2,9 comrise entre dU et 4-dU => absence d'autocorrelation

    # fonction manuelle, suivant les étapes dans le poly 
cochrane_orcutt_moi <- function(OLS, seuil = 1e-4, max_iter = 100) {
   Initialisation
  iter <- 0
  diff <- seuil + 1  # Pour entrer dans la boucle
  residus <- residuals(OLS)  # Résidus initiaux du modèle OLS
  y_current <- y  # Variable dépendante
  x_current <- x  # Variables explicatives
  while (diff > seuil && iter < max_iter) {
    # Étape 2 : Calcul de rho_chapeau
    rho_chap <- sum(residus[-1] * residus[-length(residus)]) / sum(residus^2)
    
    # Étape 3 : Transformation des variables (données différenciées)
    y_trans <- y_current[-1] - rho_chap * y_current[-length(y_current)]
    X_trans <- x_current[-1, ] - rho_chap * x_current[-nrow(x_current), ]
    
    # Étape 4 : Régression des variables transformées
    model_trans <- lm(y_trans ~ X_trans - 1)  # Régression sans intercept (X_trans contient déjà l'intercept ajusté)
    
    # Mise à jour des résidus
    new_residus <- residuals(model_trans)
    
    # Calcul de la différence entre rho courant et précédent (critère de convergence)
    diff <- abs(rho_chap - sum(new_residus[-1] * new_residus[-length(new_residus)]) / sum(new_residus^2))
    
    # Mise à jour pour la prochaine itération
    residus <- new_residus
    iter <- iter + 1
  }
  
  # Résultats finaux
  list(
    model = model_trans,
    rho = rho_chap,
    iterations = iter,
    converged = iter < max_iter
  )
}

# application de la méthode manuelle au modèle initial
#OLS_corrige_moi <- cochrane_orcutt_moi(OLS)

# Résumé du modèle corrigé
#OLS_corrige_moi
#summary(OLS_corrige_moi$model)
# re calcul de la statistique de Durbin Watson
#res_corr <- residuals(OLS_corrige_moi)  # Extraire les résidus corrigés
#n_corr <- length(res_corr)
#d1_1 = sum(res_corr^2)
#d2_1 =  t(res_corr[2:n_corr]-res_corr[1:n_corr-1]) %*% (res_corr[2:n_corr]-res_corr[1:n_corr-1])
#dw_1 = d2_1/d1_1
#print (dw_1)

# Affichage de rho et nombre d'itérations
#cat("Rho:", OLS_corrige_moi$rho, "\n")
#cat("Nombre d'itérations:", OLS_corrige_moi$iterations, "\n")



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
#resw = WAUX$residuals
#SCRw = t(resw) %*% resw
#Rw = 1 - SCRw/(var(e2)*(n-1))
#Rw_sq = Rw^2
#print(n*Rw_sq)
# la valeur obtenue est inférieure à la valeur seuil chi-deux(9) pour alpha = 0.05 (=16.919), donc on ne rejette pas l'hyp. d'homoscédasticité

 # methode 2 de white avec bibliotheque
library("skedastic")
skedastic_package_white  <- white(mainlm =  OLS1, interactions = TRUE)
skedastic_package_white
# on obtient p=0,342, on peut donc pas rejeter H0 : absence d'heteroscedasticité
 # methode 3 
bptest(OLS1, ~ fitted(OLS1) + I(fitted(OLS1)^2))

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

# Test Cusum

#Wr <- efp(y ~ x, type = "Rec-CUSUM")
#plot(Wr)
# résultat inintéressant

# Test Cusum Square

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


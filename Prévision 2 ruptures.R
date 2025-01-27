# PREVISION en prenant en compte 2 ruptures

rm(list=ls())  # Effacer les donnees en memoire
graphics.off() # Fermer les graphiques
library(readxl)

# Lecture des données

donnees <- read_excel('données.xlsx')
attach(donnees)

n_y = 32 # nombre d'années dans les données de base

# Operations sur les variables
Y <- (donnees$Celec_menages * 10^3)[1:n_y] / donnees$Pop1[1:n_y]  # Consommation en MWh par habitant
X1 <- (donnees$PIB2020 * 10^9) / donnees$Pop1      # PIB réel en euros par habitantBASE 2020
X2 <- donnees$Pelec / (donnees$`IPC(base100=2015)` / 100)  # Prix électricité corrigé en euro / MWh / hab
X3 <- donnees$DJU                               # Indice climatique

n_x = length(X1) # nombre d'années totales dans les variables explicatives(en comptant celles qu'on aura rajouté jusqu'à 2030)
n_test = 4 # nombre d'années utilisées pour calculer la qualité de la prédiction
n = n_y - n_test # nombre d'années utilisées pour la régression


# Introduction de variables muettes
rupture_1 = 19 # 1re rupture due à la libéralisation du marché après le traité de Lisbonne
rupture_2 = 25 # 2è rupture due à l'essor des renouvelables ?

# passage au log pour la régression
vec <- c(X1,X2,X3)
X <- matrix(vec, ncol=3) 
Y=matrix(Y,n_y,1)
k=ncol(X);
K=k+1
y=log(Y) 
x=log(X)

# création de la variable muette 
muet_fin <- c(rep(0, rupture_2-1), rep(1, n_x - rupture_2 + 1))
muet_mid <- c(rep(0, rupture_1-1), rep(1, rupture_2 - rupture_1), rep(0, n_x - rupture_2 + 1))
x2_1 <- muet_mid * log(X2)
x2_2 <- muet_fin * log(X2)

# création de la matrice de données utilisé pour la régression (on a ajouté la variable muette aux données de base)
x_2 <- cbind(x,x2_1, x2_2)
x_2_reg <- x_2[1:n,1:(k+2)]
y_reg <- y[1:n]

# Estimation MCO
OLS=lm(formula = y_reg ~ x_2_reg)

# Affichage des résultats
summary(OLS)
bmco = OLS$coefficients

# on en déduit les coefficients estimés pour les vrais variables (sans la muette)
bmco_f <- c(bmco[1:2],bmco[3]+bmco[5]+bmco[6], bmco[4])
print(bmco_f)

# on calcule les y estimés et les résidus, en utilisant les données de base et les coefficients qu'on vient de calculer
x_reg <- x[1:n,1:k]
xc = cbind(1,x_reg)
xt = t(xc)
xtx= xt %*% xc
xtx1 = solve(xtx)
yf = xc %*% bmco_f
u = y_reg - yf
scr = t(u) %*% u


## Calcul d'indcateurs de qualité de la prédiction qu'on veut faire avec cette régression

rmse=0
mae=0
mape=0
for(j in 1:n_test){
  xf<-x[n+j,1:3]
  yfobs<-y[n+j]
  nvx <-matrix(c(1,xf))
  nvx1 = t(nvx);
  nvy= nvx1 %*% bmco_f
  #
  # Ecart-type de prevision
  #
  sprev= sqrt( (scr/(n-K)) * (1+ nvx1 %*% xtx1 %*% nvx ))
  print(sprev)
  print(xf)
  #
  # Intervalle de prédiction MCO - Borne inf et sup à 95%
  
  # T de Student 
  #
  tstud=qt(0.975,n-K)
  #
  prevymin=nvy-tstud*sprev
  prevymax=nvy+tstud*sprev
  
  # Affichage 
  
  cat("Intervalles de confiance",j, ":",nvy, prevymin, prevymax,sprev, "\n")
  rmse=rmse+(yfobs-nvy)*(yfobs-nvy)
  mae = mae + abs(yfobs-nvy)
  mape = mape + abs((yfobs-nvy)/yfobs)
}
rmse=sqrt(rmse/n_test)
mae = mae/n_test
mape = (mape/n_test) * 100
cat("RMSE", rmse, "\n")
cat("MAE", mae, "\n")
cat("MAPE", mape, "%", "\n")


##prédiciton de 2022 à 2030


nprev = n_x - n_y # =9, de 2022 à 2030
for(j in 1:nprev){
  xf<-x[n_y+j,1:3]
  yfobs<-y[n_y+j]
  nvx <-matrix(c(1,xf))
  nvx1 = t(nvx);
  nvy= nvx1 %*% bmco_f
  #
  # Ecart-type de prevision
  #
  sprev= sqrt( (scr/(n-K)) * (1+ nvx1 %*% xtx1 %*% nvx ))
  print(sprev)
  print(xf)
  #
  # Intervalle de prédiction MCO - Borne inf et sup à 95%
  
  # T de Student 
  #
  tstud=qt(0.975,n-K)
  #
  prevymin=nvy-tstud*sprev
  prevymax=nvy+tstud*sprev
  
  # Affichage 
  
  cat("Intervalles de confiance",2021 +j, ":",nvy, prevymin, prevymax,sprev, "\n")
}


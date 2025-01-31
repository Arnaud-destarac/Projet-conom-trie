# PREVISION en prenant en compte 1 rupture

####SCENARIO "Pessimiste" : même fichier que "Prévision scén opt.R" sauf en ligne 12, ici sheet 1

rm(list=ls())  # Effacer les donnees en memoire
graphics.off() # Fermer les graphiques
library(readxl)

# Lecture des données

donnees <- read_excel('données-conso_elec.xlsx')
donnees_prev <- read_excel('donnees_prev.xlsx', sheet=1)
attach(donnees)
attach(donnees_prev)

n_y = 32 # nombre d'années dans les données de base

# Operations sur les variables
Y <- (donnees$Celec_menages * 10^3)[1:n_y] / donnees$Pop1[1:n_y]  # Consommation en MWh par habitant
X1 <- (donnees$PIB2020 * 10^9) / donnees$Pop1      # PIB réel en euros par habitantBASE 2020
X2 <- donnees$Pelec / (donnees$`IPC(base100=2015)` / 100)  # Prix électricité corrigé en euro / MWh
X3 <- donnees$DJU                               # Indice climatique

n_x = length(X1) # nombre d'années totales dans les variables explicatives(en comptant celles qu'on aura rajouté jusqu'à 2030)
n_test = 4 # nombre d'années utilisées pour calculer la qualité de la prédiction
n = n_y - n_test # nombre d'années utilisées pour la régression


rupture = 19 # rupture en 2008 cf fichier MCO

# Passage au log pour la régression
vec <- c(X1,X2,X3)
X <- matrix(vec, ncol=3) 
Y=matrix(Y,n_y,1)
k=ncol(X);
K=k+1
y=log(Y) 
x=log(X)

# Récupération des données prévisionnelles
X1_prev <- donnees_prev$PIB_par_hab 
X2_prev <- donnees_prev$Pelec_base_2015
X3_prev <- donnees_prev$DJU 
vec_p <- c(X1_prev,X2_prev,X3_prev)
X_prev <- matrix(vec_p, ncol=3) 
x_prev = log(X_prev)

# Création de la variable muette 
muet_fin <- c(rep(0, rupture-1), rep(1, n_x - rupture + 1))
x2_b <- muet_fin * log(X2)

# Création de la matrice de données utilisé pour la régression (avec variable muette)
x_2 <- cbind(x,x2_b)
x_2_reg <- x_2[1:n,1:(k+1)]
y_reg <- y[1:n]

# Estimation MCO
OLS=lm(formula = y_reg ~ x_2_reg)

# Affichage des résultats
summary(OLS)
bmco = OLS$coefficients

# Calcul des coefficients pour les vrais variables (sans la muette)
bmco_f <- c(bmco[1:2],bmco[3]+bmco[5], bmco[4])
print(bmco_f)

# Estimation de la variable endogène à partir des coeffs qu'on vient d'estimer, et calcul des résidus
x_reg <- x[1:n,1:k]
xc = cbind(1,x_reg)
xt = t(xc)
xtx= xt %*% xc
xtx1 = solve(xtx)
yf = xc %*% bmco_f
u = y_reg - yf
scr = t(u) %*% u


## Calcul d'indicateurs de qualité de la prédiction

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

Pop <- Population
nprev = length(X1_prev) # =9, de 2022 à 2030
y_prev <- rep(0,nprev)
for(j in 1:nprev){
  xf<-x_prev[j,1:3]
  nvx <-matrix(c(1,xf))
  nvx1 = t(nvx);
  nvy= nvx1 %*% bmco_f
  pop = Pop[j]*1e6
  y_prev[j] = nvy
  #
  # Ecart-type de prevision
  #
  sprev= sqrt( (scr/(n-K)) * (1+ nvx1 %*% xtx1 %*% nvx ))
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
  
  cat("Intervalles de confiance",2021 +j, ":", exp(prevymin)*1e-3*pop, exp(nvy)*1e-3*pop, exp(prevymax)*1e-3*pop,sprev, "\n")
}

# Plot des résultats de prédiction

années <- c(Date,Année)
Y_prev = exp(y_prev)
Conso_prev_2 <- (Y_prev*1e-3) * (Pop*1e6) # car Y_prev est en MWh/hab et Pop en Millions d'habitants
Conso <- Celec_menages
Conso_prol <- c(Conso,Conso_prev_2)
Conso_print<- (Conso_prev_1)

ybnd=c(0.9*range(Conso_prol)[1], 1.1*range(Conso_prol)[2])
plot(années,Conso_prol,xlab="années",ylab="Consommation (GWh)",col="red",xlim=range(Années),ylim=ybnd,type="p", lwd=2)

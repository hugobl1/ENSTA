###########################################################
#Exercice 2

#Question 1
set.seed(999)
help("set.seed")

#Initialisation des donn�es
n=50
lambda=0.3
a=5
b=1
var=2
x_i=rnorm(n)
eps_i=rnorm(n,0,sqrt(var))
hist(x_i,main=paste("R�partition de l'�chantillon X"))
hist(eps_i,main=paste("R�partition des erreurs"))
#On calcule l'�chantillon Z, d'apr�s les informations de l'�nonc�
z_i=a+b*x_i+eps_i
hist(z_i,main=paste("R�partition de l'�chantillon Z"))
#Puis l'�chantillon Y, gr�ce � la relation h_lambda, formule de transformation de
#Bickel et Docksum:
y_i=sign(lambda*z_i+1)*abs(lambda*z_i+1)^(1/lambda)
hist(y_i,main=paste("R�partition de l'�chantillon Y"))
#V�rification:
z_i-sign(y_i)*((abs(y_i)^lambda)-1)/lambda
#[1]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#[9]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.665335e-16  0.000000e+00
#[17] -8.881784e-16 -8.881784e-16  0.000000e+00  0.000000e+00 -8.881784e-16  0.000000e+00 -8.881784e-16  0.000000e+00
#[25]  0.000000e+00  0.000000e+00 -4.440892e-16 -8.881784e-16  4.440892e-16  0.000000e+00  0.000000e+00  0.000000e+00
#[33] -8.881784e-16  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  4.440892e-16
#[41]  0.000000e+00  8.881784e-16  0.000000e+00  0.000000e+00  0.000000e+00 -4.440892e-16  0.000000e+00  4.440892e-16
#[49]  0.000000e+00  8.881784e-16
#Ok, c'est bon!

#R�gression lin�aire simple de z en fonction de x:
res_z=lm(z_i~x_i)
summary(res_z)


#R�gression lin�aire simple de y en fonction de x:
res_y=lm(y_i~x_i)
summary(res_y)

#On repr�sente les r�sultats:
#Pour z en fonction de x:

plot(x_i,z_i, main="R�gression lin�aire de Z en fonction de X", xlab = "x",ylab = "z")
abline(res_z$coefficients[1],res_z$coefficients[2],col="red")

#Pour y en fonction de x:

plot(x_i,y_i,main="R�gression lin�aire de Y en fonction de X", xlab = "x",ylab = "y")
abline(res_y$coefficients[1],res_y$coefficients[2],col="red")

#Etude des r�sidus

library(MASS) # pour la fonction stdres

#Cas de la r�gression lin�aire simple de z fonction de x

plot(res_z$fitted,z_i,main="Valeurs ajust�s/observ�s (R�gression z~x)")

abline(0,1,col="red")
# les points alignent assez bien suivant la premi�re bissectrice,

plot(res_y$fitted,y_i,main="Valeurs ajust�s/observ�s (R�gression y~x)")

abline(0,1,col="red")
# les points alignent moyennement bien suivant la premi�re bissectrice,
# ce qui est un gage visuel que l'ajustement assez bon mais l�g�rement moins bon qu'avant


plot(res_z$fitted,studres(res_z),col=3,pch=3,main="Diff�rents r�sidus en fonction des valeurs ajust�es(z~x)") #studentis�s par validation crois�e
points(res_z$fitted,stdres(res_z), col=2,pch=2 )# fitted to variance 1

points(res_z$fitted,res_z$residuals,main="diff�rents r�sidus")# non norm�s
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)

# les r�sidus studentis�s ne montrent pas de d�pendance ou de tendance,
# ils sont raisonnablement compris entre -2 et 2, sauf deux points
plot(res_z$fitted,res_z$residuals,main="diff�rents r�sidus")# non norm�s

plot(res_y$fitted,studres(res_y),col=3,pch=3,main="Diff�rents r�sidus en fonction des valeurs ajust�es(y~x)") #studentis�s par validation crois�e
points(res_y$fitted,stdres(res_y), col=2,pch=2 )# fitted to variance 1
points(res_y$fitted,res_y$residuals,main="diff�rents r�sidus")# non norm�s
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)
# quelques r�sidus ne sont pas compris entre -2 et 2
##
plot(res_y$fitted,res_y$residuals,main="diff�rents r�sidus")# non norm�s

qqnorm(studres(res_z),main="graphe quantile-quantile")
qqline(stdres(res_z)) # passe par le 1er et le 3�me quartile
# les points s'alignent bien,
# il n'y a pas de raison de mettre en doute la gaussiannit� des r�sidus.
shapiro.test(studres(res_z))
##
qqnorm(studres(res_y),main="graphe quantile-quantile")
qqline(stdres(res_y)) # passe par le 1er et le 3�me quartile
# les points s'alignent moins bien que sur l'exemple pr�c�dent
shapiro.test(studres(res_y))
##

#> shapiro.test(studres(res_y))
#
#Shapiro-Wilk normality test

#data:  studres(res_y)
#W = 0.93397, p-value = 0.007836
# Le test de Shapiro rejette l'hypoth�se de gaussianit�
# (avec une risque de premi�re esp�ce 0.05)
# outer pour titre g�n�ral hors du cadre du dernier plot
title(main="Validation en r�gression lin�aire multiple", outer=TRUE)

#Partie 2

#Question 1
library(MASS)
library(car)
X = as.matrix(cbind(1,x_i))
(t(X)%*%X)/n

#Question 2

Q = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X)
Lmle = function(Z){
  n = length(Z)
  sig2 = ( t(Z)%*%Q%*%Z )/n
  -n/2*log(sig2)
}

h_lambda=function(lambda,y){
  (sign(y)*abs(y)^(lambda)-1)/lambda
}



lmin= function(lambda,Y){
  n=length(Y)
  lmax=(lambda-1)*sum(log(abs(Y)))+Lmle(h_lambda(lambda,Y))-(n/2)*(1+log(2*pi))
  -lmax
}



Vlmin = Vectorize(lmin,"lambda")

curve(expr=(Vlmin(x,y_i)),from=0, to= 2)

#Question 3

resopt = nlm(Vlmin,p=0.4,Y=y_i,hessian=TRUE) 

lambda=resopt$estimate
hessian=resopt$hessian
V=(hessian)^{-1}

#Question 4

qn=qnorm(0.975)
Ic_b=lambda-sqrt(V)*qn
Ic_h=lambda+sqrt(V)*qn

lambda_0=c(1 ,1/2, 0.3 ,0)
W_vect=V^{-1}*(lambda*c(1,1,1,1)-lambda_0)^2
W_vect>qchisq(0.95,1)
1-pchisq(W_vect,1)
#Question 5
#On modifie les fonctions pour le test faisant intervenir h_tilde_lambda

h_tilde_lambda=function(lambda,y){
  if(lambda==0){
    rep=log(y)
  }
  else {
    rep=(y^(lambda)-1)/lambda
  }
  rep
}
lmin_tilde= function(lambda,Y){
  n=length(Y)
  lmax=(lambda-1)*sum(log(abs(Y)))+Lmle(h_tilde_lambda(lambda,Y))-(n/2)*(1+log(2*pi))
  -lmax
}
lambda_0=c(1 ,1/2, 0.3 )
LRV_vect=2*(-Vlmin(lambda,y_i)+Vlmin(lambda_0,y_i))
LRV_vect>qchisq(0.95,1)
1-pchisq(LRV_vect,1)

LRV_0=2*(-lmin_tilde(lambda,y_i)+lmin_tilde(0,y_i))
LRV_0>qchisq(0.95,1)
1-pchisq(LRV_0,1)
#Question 6

rep=powerTransform(y_i~x_i)
summary(rep)

#Question 7

n=50
lambda=0.3
a=5
b=1
var=2
compteur=0
for (i in 1:100){
  x_test=rnorm(n)
  eps_test=rnorm(n,0,sqrt(var))
  z_test=a+b*x_test+eps_test
  y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
  reptest=powerTransform(y_test~x_test)
  if (testTransform(reptest,0.3)$pval<0.05)
  {
    compteur=compteur+1
  }
}
compteur
Nb_moyen_rejets=compteur/100

#Pour la puissance
n=50
lambda=0.1
a=5
b=1
var=2
compteur=0
for (i in 1:100){
  x_test=rnorm(n)
  eps_test=rnorm(n,0,sqrt(var))
  z_test=a+b*x_test+eps_test
  y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
  reptest=powerTransform(y_test~x_test)
  if (testTransform(reptest,0.3)$pval<0.05)
  {
    compteur=compteur+1
  }
}
compteur
Nb_moyen_rejets=compteur/100

#On transforme cela en fonction:
puissance=function(lambda)
{
  n=50
  a=5
  b=1
  var=2
  compteur=0
  for (i in 1:100){
    x_test=rnorm(n)
    eps_test=rnorm(n,0,sqrt(var))
    z_test=a+b*x_test+eps_test
    y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
    reptest=powerTransform(y_test~x_test,family = "yjPower")
    if (testTransform(reptest,0.3)$pval<0.05)
    {
      compteur=compteur+1
    }
  }
  compteur
  Nb_moyen_rejets=compteur/100
  Nb_moyen_rejets
}
puissvectoriz = Vectorize(puissance,"lambda")

curve(expr=(puissvectoriz(x)),from=0.05, to= 1.5)
lambda=0.9
x_test=rnorm(n)
eps_test=rnorm(n,0,sqrt(var))
z_test=a+b*x_test+eps_test
y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
reptest=powerTransform(y_test~x_test)

#Bonus

n=50
lambda=0
a=5
b=1
var=2
compteur=0
for (i in 1:100){
  x_test=rnorm(n)
  eps_test=rnorm(n,0,sqrt(var))
  z_test=a+b*x_test+eps_test
  y_test=exp(z_test)
  reptest=powerTransform(y_test~x_test)
  if (testTransform(reptest,0)$pval<0.05)
  {
    compteur=compteur+1
  }
}
compteur
Nb_moyen_rejets=compteur/100

#Pour la puissance
n=50
lambda=0.1
a=5
b=1
var=2
compteur=0
for (i in 1:100){
  x_test=rnorm(n)
  eps_test=rnorm(n,0,sqrt(var))
  z_test=a+b*x_test+eps_test
  y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
  reptest=powerTransform(y_test~x_test)
  if (testTransform(reptest,0)$pval<0.05)
  {
    compteur=compteur+1
  }
}
compteur
Nb_moyen_rejets=compteur/100

#On transforme cela en fonction:
puissancebis=function(lambda)
{
  n=50
  a=5
  b=1
  var=2
  compteur=0
  for (i in 1:100){
    x_test=rnorm(n)
    eps_test=rnorm(n,0,sqrt(var))
    z_test=a+b*x_test+eps_test
    y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
    reptest=powerTransform(y_test~x_test)
    if (testTransform(reptest,0)$pval<0.05)
    {
      compteur=compteur+1
    }
  }
  compteur
  Nb_moyen_rejets=compteur/100
  Nb_moyen_rejets
}
puissvectorizbis = Vectorize(puissancebis,"lambda")

curve(expr=(puissvectorizbis(x)),from=0.1, to= 0.3)
lambda=0.9
x_test=rnorm(n)
eps_test=rnorm(n,0,sqrt(var))
z_test=a+b*x_test+eps_test
y_test=sign(lambda*z_test+1)*abs(lambda*z_test+1)^(1/lambda)
reptest=powerTransform(y_test~x_test)

###########################################################
#Exercice 3

library(MASS)
library(car)
##---Cas pratique

#Question 1 : Lecture des donnees 

cycle <- read.csv('NbCycleRupture.csv',header = TRUE, sep=";")

#verification des donnees importees, on a bien 27 observations et 4 variables 
cycle

#Def du modele : Y =alpha+ beta1 a1 +  beta_2 a2 + beta_3 a3 + epsilon
a1 = (50*cycle$x1 + 300)
a2 = cycle$x2 + 9
a3 = (5*cycle$x3 + 45)
df = data.frame(cycle$y,a1,a2,a3)

reg_multi <-lm(cycle$y~a1 + a2 + a3, data = df) #Ajustement du modele


summary(reg_multi)



reg_multi_bis <- lm(y~x1+x2+x3,data=cycle)
summary(reg_multi_bis)

reg_multi2<-lm(y~x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3),data=cycle)
summary(reg_multi2)
#On obtient un R tres grand... On a trop ajuste ! 

#V�rification des hypoth�ses du mod�le
#Test de normalit� des r�sidus:
#Pour le premier mod�le:
qqnorm(studres(reg_multi),main="graphe quantile-quantile")
qqline(stdres(reg_multi)) # passe par le 1er et le 3�me quartile
shapiro.test(stdres(reg_multi))
#Shapiro-Wilk normality test

#data:  stdres(reg_multi)
#W = 0.86097, p-value = 0.001914
#La p-value est inf�rieure � 0.05, on peut donc rejeter l'hypoth�se des r�sidus gaussiens

#Pour le second mod�le:
qqnorm(studres(reg_multi_bis),main="graphe quantile-quantile")
qqline(stdres(reg_multi_bis)) # passe par le 1er et le 3�me quartile
shapiro.test(stdres(reg_multi_bis))
#Shapiro-Wilk normality test
#data:  stdres(reg_multi_bis)
#W = 0.86097, p-value = 0.001914
#La p-value est inf�rieure � 0.05, on peut donc rejeter l'hypoth�se des r�sidus gaussiens

#Affiches des figures pour la validation:
#Mod�le 1
par(mfrow=c(2,2))
plot(reg_multi)

par(mfrow=c(1,1))
plot(reg_multi_bis$fitted,studres(reg_multi_bis),col=3,pch=3) #studentis�s par validation crois�e
points(reg_multi_bis$fitted,stdres(reg_multi_bis), col=2,pch=2 )# fitted to variance 1
points(reg_multi_bis$fitted,reg_multi_bis$residuals,main="diff�rents r�sidus")# non norm�s
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)


par(mfrow=c(2,2))
plot(reg_multi_bis$fitted,reg_multi_bis$residuals,main="diff�rents r�sidus")# non norm�s
plot(reg_multi$fitted,reg_multi$residuals,main="diff�rents r�sidus")# non norm�s

#Mod�le 2
par(mfrow=c(2,2))
plot(reg_multi2)
points(reg_multi2$fitted,stdres(reg_multi2), col=2,pch=2 )# fitted to variance 1
points(reg_multi2$fitted,reg_multi2$residuals,main="diff�rents r�sidus")# non norm�s

par(mfrow=c(1,1))
plot(reg_multi2$fitted,studres(reg_multi2),col=3,pch=3) #studentis�s par validation crois�e
points(reg_multi2$fitted,stdres(reg_multi2), col=2,pch=2 )# fitted to variance 1
points(reg_multi2$fitted,reg_multi2$residuals,main="diff�rents r�sidus")# non norm�s
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)


qqnorm(studres(reg_multi2),main="graphe quantile-quantile")
qqline(stdres(reg_multi2)) # passe par le 1er et le 3�me quartile
shapiro.test(stdres(reg_multi2))
#Shapiro-Wilk normality test

#data:  stdres(reg_multi)
#W = 0.86097, p-value = 0.2904
#La p-value est sup�rieure � 0.05, on ne rejete donc pas l'hypoth�se des r�sidus gaussiens



#Question 2
p=10
n=27
r=6
A=matrix(c(0,0,0,0,1,0,0,0,0,0,
           0,0,0,0,0,1,0,0,0,0,
           0,0,0,0,0,0,1,0,0,0,
           0,0,0,0,0,0,0,1,0,0,
           0,0,0,0,0,0,0,0,1,0,
           0,0,0,0,0,0,0,0,0,1
),nrow=6,ncol=10,byrow=TRUE)
theta=reg_multi2$coefficients
V=vcov(reg_multi2)

sigma=sqrt(sum((cycle$y-reg_multi2$fitted.values)^2)/(n-p))

W=(t(A%*%theta)%*%solve(A%*%V%*%t(A))%*% A%*%theta)/r

qf(0.95,r,n-p)

abs(W)>qf(0.95,r,n-p)

#Wbis=(t(A%*%theta)%*%solve(A%*%%*%t(A))%*% A%*%theta/r)/((sigma)^2/(n-p))
plot(reg_multi2)




#Question 3

#Pour le mod�le (M1bis)
repex3_1bis=powerTransform(reg_multi_bis)
summary(repex3_1bis) 
trans_y=log(cycle$y)
trans_reg_1bis=lm(trans_y~cycle$x1+cycle$x2+cycle$x3)
summary(trans_reg_1bis)
plot(trans_reg_1bis$fitted,stdres(trans_reg_1bis), col=2,pch=2 ,main="R�sidus Studentis�s/Valeurs ajust�es (Mod�le M1bis)")
shapiro.test(studres(trans_reg_1bis))
#Pour le mod�le (M2bis)
repex3_2bis=powerTransform(reg_multi2)
summary(repex3_2bis)
trans_reg_2bis=lm(log(y)~x1+x2+x3+I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3),data = cycle)
summary(trans_reg_2bis)
plot(trans_reg_2bis$fitted,stdres(trans_reg_2bis), col=2,pch=2, main="R�sidus Studentis�s/Valeurs ajust�es (Mod�le M2bis)" )
shapiro.test(studres(trans_reg_2bis))
#Question 4

p=10
n=27
r=6
A=matrix(c(0,0,0,0,1,0,0,0,0,0,
           0,0,0,0,0,1,0,0,0,0,
           0,0,0,0,0,0,1,0,0,0,
           0,0,0,0,0,0,0,1,0,0,
           0,0,0,0,0,0,0,0,1,0,
           0,0,0,0,0,0,0,0,0,1
),nrow=6,ncol=10,byrow=TRUE)
theta=trans_reg_2bis$coefficients
V=vcov(trans_reg_2bis)

sigma=sqrt(sum((log(cycle$y)-trans_reg_2bis$fitted.values)^2)/(n-p))

W=(t(A%*%theta)%*%solve(A%*%V%*%t(A))%*% A%*%theta)/r

qf(0.95,r,n-p)

abs(W)>qf(0.95,r,n-p)

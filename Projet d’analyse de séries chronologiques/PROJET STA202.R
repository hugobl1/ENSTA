setwd("/Users/33628/Desktop/STA202")
rm(list=objects())
library(zoo)
library(timeDate)
library(forecast)
library(xts)
library(mgcv)

####MISE EN FORME DES DONNEES

###IMPORT DES DONNEES

#lecture du fichier conso_quotidienne
data<-read.csv("conso_quotidienne_brute_triee.csv", sep = ";")

#creation des dates pour construire une serie temporelle
date1<-strptime(c("01/01/2013 00:00"),"%d/%m/%Y %H:%M")
date2<-strptime(c("31/12/2020 23:30"),"%d/%m/%Y %H:%M")
date<-seq(date1,date2, by = "30 min")

#recuperation des donnees de conso_electrique pour la region PACA
data1<-data[data$region == "Provence-Alpes-CÃ´te d'Azur",]$consommation_brute_electricite_rte

head(data1)

#construction de la serie temporelle
X<-as.vector(na.locf(data1))
X<-c(X[1],X)
conso<-c(X)
conso.xts<-xts(conso, order.by = date)



head(conso.xts)



conso.ts<-ts(conso,start = 1, frequency = 48)
plot(conso.ts, type ="l",col = "red")

#Creation des variables de conso journalieres
#on utilise la methode du TP1
Date = format(date,"%Y/%m/%d")
mean.day<-tapply(conso,as.factor(Date),mean)
head(mean.day)
Y<-as.vector(mean.day)
Y.train=Y[1:2000]
Y.predict=Y[2001:2800]


new_date1<-strptime(c("01/01/2013"),"%d/%m/%Y")
new_date2<-strptime(c("31/12/2020"),"%d/%m/%Y")
new_date<-seq(new_date1,new_date2, by = "1 day")

week_date<-seq(new_date1)

conso_frame<-data.frame(Y,new_date)
conso_jour.xts<-xts(Y, order.by = new_date)

head(conso_frame)

### ANALYSE DESCRIPTIVE SUR LA SERIE TEMPORELLE BRUTE

##profils
#mensuel
month<-format(new_date,"%m")
mean.month<-tapply(conso_jour.xts,as.factor(month),mean)
plot(mean.month, type ="b")
#annuel
year<-format(new_date,"%Y")
mean.year<-tapply(conso_jour.xts,as.factor(year),mean)
plot(mean.year, type ="b")
#journalier
day<-format(new_date,"%d")
mean.day<-tapply(conso_jour.xts,as.factor(day),mean)
plot(mean.day,type='b')
#semaine
dow<-as.factor(.indexwday(conso.xts))
conso_day<-tapply(conso.xts,dow,mean)
par(mfrow=c(1,2))
plot(conso_day,type='b',xlab="Jour de la semaine",ylab="Consommation électrique sur 30 minutes")
boxplot(X~ dow, col = "plum4",xlab="Jour de la semaine",ylab="Consommation électrique sur 30 minutes ")

mean.week=conso_day


boxplot(conso_frame$Y~ month, col = "lightblue")
title("Répartition des moyennes par tranche de 30 minutes par mois")


#statistique de base, boxplot et graphe des conso par mois puis par an

par(mfrow=c(3,2))
plot(conso_day,type='b',xlab="Jour de la semaine",ylab="Consommation électrique sur moyenne 30 minutes")
title("Répartition des moyennes hebdommadaires par tranche de 30 minutes")
boxplot(X~ dow, col = "plum4",xlab="Jour de la semaine",ylab="Consommation électrique moyenne sur 30 minutes ")
title("Répartition des moyennes hebdommadaires par tranche de 30 minutes")
plot(mean.month,type='b', xlab="Mois",ylab="Consommation électrique moyenne sur 30 minutes")
title("Répartition des moyennes mensuelles pour des tranche de 30 minutes")
boxplot(conso_frame$Y~ month, col = "lightblue",xlab="Mois",ylab="Consommation électrique moyenne sur 30 minutes")
title("Répartition des moyennes mensuelles par tranche de 30 minutes")
plot(mean.year, type ='b', xlab="Mois",ylab="Consommation électrique moyenne sur 30 minutes" )
title("Répartition des moyennes annuelles par tranche de 30 minutes")
boxplot(conso_frame$Y ~ year, col = "plum4", xlab="Mois",ylab="Consommation électrique moyenne sur 30 minutes")
title("Répartition des moyennes annuelles par tranche de 30 minutes")

#resume des statistiques de base 
summary(mean.month)
summary(mean.year)

par(mfrow=c(2,2))
hist(mean.month,col = "lightblue")
title("Histogramme des moyennes mensuelles par tranche de 30 min")
hist(mean.year,col = "plum4")
boxplot(conso_frame$Y~ day, col = "red")

dev.off()



#correlations
acf(conso_frame$Y, lag.max = 21)
pacf(conso_frame$Y, lag.max = 21)

autocorr<-function(x,h)
{
  x.lag<-lag.xts(x, k = h, na.pad = TRUE)
  return(cor(x.lag,x,use="pairwise.complete.obs"))
}


resultat_corr<-sapply(c(1:15),autocorr,x = conso_jour.xts)

plot(resultat_corr)




####MODELISATION ET PREVISIONS
conso_jour.xts.train=conso_jour.xts[1:2000]
new_date.train=new_date[1:2000]
new_date.predict=new_date[2001:2800]
###Analyse de la tendance 
##Moyenne mobile

l<-365
MA<-filter(conso_jour.xts.train, filter = array(1/l, dim = l), method = c("convolution"),sides = 2, circular = F)
MA.xts<-xts(MA, order.by = new_date.train)

plot(conso_jour.xts.train, type = "l", col ="green")
lines(MA.xts, col = "red")
title("moyenne mobile pour l = 365")


##Estimateur Ã noyau  
n<-length(Y.train)
t<-seq(1,n, by = 1)
x<-seq(1,n, by = 1)

#noyau_gaussien
noyau_gaussien<-function(x1)
{
  K<-dnorm((x1-t),mean = 0, sd = sqrt(h/2))/sum(dnorm((x1-t),mean = 0, sd = sqrt(h/2)))
  f<-sum(Y.train*K)
  return(f)
}

plot(conso_jour.xts.train, type = "l", col ="green")

par(mfrow=c(2,2))
dev.off()

h<-1000
F1<-lapply(x,noyau_gaussien)
F1<-unlist(F1)
F1.xts<-xts(F1,order.by = new_date.train)
lines(F1.xts, type ="l",col = "plum4")
title("noyau gaussien h = 1000")

h<-5000
F2<-lapply(x,noyau_gaussien)
F2<-unlist(F2)
F2.xts<-xts(F2,order.by = new_date.train)
lines(F2.xts, type ="l",col = "plum4")
title("noyau gaussien h = 5000")


h<-20000
F3<-lapply(x,noyau_gaussien)
F3<-unlist(F3)
F3.xts<-xts(F3,order.by = new_date.train)
lines(F3.xts, type ="l",col = "plum4")
title("noyau gaussien h = 20000")


h<-50000
F4<-lapply(x,noyau_gaussien)
F4<-unlist(F4)
F4.xts<-xts(F4,order.by = new_date.train)
lines(F4.xts, type ="l",col = "plum4")
title("noyau gaussien h = 50000")

#noyau



##Polynome locaux 
conso.lo<-loess(Y.train~t,degree = 2, span = 0.95)
conso.lo.xts<-xts(conso.lo$fitted, order.by = new_date.train)
plot(conso_jour.xts.train,type ="l", col ="green")
lines(conso.lo.xts, col ="red")
title("polynome locaux span = 0.95")


#regression lineaire pour tendance polynomiale
reg<-lm(Y.train~1+t)
reg.xts<-xts(reg$fitted.values,order.by = new_date.train)
plot(conso_jour.xts.train,type ="l", col ="green")
lines(reg.xts, col = "blue")
title("RÃ©gression linÃ©aire")
summary(reg)


###Analyse de la saisonnalite
##Modele moyenne mobile
l<-15
MA<-filter(conso_jour.xts.train, filter = array(1/l, dim = l), method = c("convolution"),sides = 2, circular = F)
MA.xts<-xts(MA, order.by = new_date.train)

plot(conso_jour.xts.train, type = "l", col ="green")
lines(MA.xts, col = "red")
title("moyenne mobile pour l = 15")

##Modele noyau gaussien
Y.trendless.gaussien<-Y.train-F4
Y.trendless.gaussien.xts<-xts(Y.trendless.gaussien, order.by = new_date.train)


#noyau_gaussien_bis
noyau_gaussien_bis<-function(x1)
{
  K<-dnorm((x1-t),mean = 0, sd = sqrt(h/2))/sum(dnorm((x1-t),mean = 0, sd = sqrt(h/2)))
  f<-sum(Y.trendless.gaussien*K)
  return(f)
}

h<-100
F1.trendless<-lapply(x,noyau_gaussien_bis)
F1.trendless<-unlist(F1.trendless)
F1.trendless.xts<-xts(F1.trendless,order.by = new_date.train)
plot(Y.trendless.gaussien.xts, col = "green")
lines(F1.trendless.xts, type ="l",col = "plum4")





##Modele regression lineaire 
Y.trendless.lineaire<-Y.train - reg$fitted.values
Y.trendless.lineaire.xts<-xts(Y.trendless.lineaire, order.by = new_date.train)
plot(Y.trendless.gaussien.xts, col = "green")


##Modele polynome locaux 

Y.trendless<-Y.train-conso.lo$fitted
conso.trendless<-data.frame(t,Y.trendless,new_date.train)
Y.trendless.xts<-xts(Y.trendless,order.by = new_date.train)
l.trendless<-7
month<-as.numeric(format(new_date.train,"%m"))







###Predictions 

##estimateur parametrique - serie de Fourier
w=2*pi/365
fourier<-cbind(cos(w*t[1:2000]), sin(w*t[1:2000]))
N<-20
for(i in c(2:N))
{
  fourier<-cbind(fourier,cos(i*w*t), sin(i*w*t))
}
matplot(fourier[,1:10],type='l')
dim(fourier)

head(fourier[,1:10])

reg.fourier<-lm(Y.trendless.gaussien[1:2000] ~ fourier[,1:16]-1)

summary(reg.fourier)

plot(Y, type ="l", col ="green")
lines(reg.fourier$fitted.values+reg$fitted.values, col = "red")

coefficients=reg.fourier$coefficients


predict=rep(0,800)
for(i in c(1:8))
{

  predict<-predict+coefficients[2*i-1]*cos(i*w*c(2001:2800))+ coefficients[2*i]*sin(i*w*c(2001:2800))
  head(predict)
}

##On prolonge la régression linéaire sur l'intervalle de test
regpredict=reg$coefficients[2]*c(2001:2800)+reg$coefficients[1]
predict=predict+regpredict
lines(c(2001:2800),predict, col = "blue")
title("Modèle paramétrique des saisonnalités par les séries de Fourier ")



## Prédiction
dev.off()
conso_ts=ts(conso_frame[1:2000,1], start=2013, frequency=365)

conso_components=decompose(conso_ts)
plot(conso_components)


##Lissage exponentiel simple
exposimple <- HoltWinters(conso_ts, beta=FALSE, gamma = FALSE)

plot(exposimple)

View(m)

predictsimple <- forecast(exposimple, h=800)
plot(predictsimple)
lines(2013+c(2001:2800)/365,conso_frame[2001:2800,1], col = "blue")

expodouble<- HoltWinters(conso_ts, gamma = FALSE)

predictdouble <- forecast(expodouble,h=800)


##Lissage exponentiel double
plot(predictdouble)
lines(2013+c(2001:2800)/365,conso_frame[2001:2800,1], col = "blue")

plot(2013+c(2001:2800)/365,predictdouble$mean,type='l')
lines(2013+c(2001:2800)/365,conso_frame[2001:2800,1], col = "blue")


##Prédiction grâce à un lissage exponentiel de HolWinters
expoHolt<- HoltWinters(conso_ts)
plot(expoHolt,xlab="Année",ylab="Observées/Ajustées")
predictHolt <- forecast(expoHolt,h=800)
predictHoltbis<-predict(expoHolt,n.ahead=800)
##Affichage données entrainement/données test
plot(2013+c(1:2800)/365,conso_frame[1:2800,1],type='l',xlab="Année",ylab="Consomation électrique sur un créneau de 30 minutes en MW")
lines(2013+c(2001:2800)/365,conso_frame[2001:2800,1], col = "blue")


#Affichage prédiction
plot(2013+c(2001:2800)/365,predictHoltbis$mean,type='l',xlab="Année",ylab="Consomation électrique sur un créneau de 30 minutes en MW")
plot(2013+c(2001:2800)/365,predictHolt$mean,type='l',xlab="Année",ylab="Consomation électrique sur un créneau de 30 minutes en MW", main="Holt-Winters Forecasting")
lines(2013+c(2001:2800)/365,conso_frame[2001:2800,1], col = "blue")

#ERREUR PREDICTION
plot(2013+c(2001:2800)/365,(conso_frame[2001:2800,1]-predictHolt$mean)/conso_frame[2001:2800,1],type='l',xlab="Année",ylab="Erreur relative",main="Erreur relative de la prédiction par lissage exponentiel de Holt-Winters")

plot(2013+c(2001:2800)/365,(conso_frame[2001:2800,1]-predict)/conso_frame[2001:2800,1],type='l',xlab="Année",ylab="Erreur relative",main="Erreur relative de la prédiction par série de Fourier") ##Fourier

#Prédiction précise pour du très court-terme (type 1 semaine)
##La prédiction va désormais être fait en modifiant les données d'entrainement et de test
short=data1[1:7056]
short<-c(short[1],short)
short=na.locf(short)
short.train=short[1:6720]
short.predict=short[6721:7056]

short_ts=ts(as.numeric(short.train), start=1, frequency=48)
plot(short_ts)
short_ts=na.locf(short_ts)

##Affichage données entrainement/données test
plot(+c(1:7056)/(48*30),data1[1:7056] ,type='l',xlab="Mois",ylab="Consomation électrique sur un créneau de 30 minutes en MW")
lines(+c(6721:7056)/(48*30),short.predict, col = "blue")

shortexpoHolt<- HoltWinters(short_ts)
shortpredictHolt <- forecast(shortexpoHolt,h=336)
plot(shortexpoHolt, xlab="Jour",ylab="Observées/Ajustées")

plot(shortpredictHolt)
plot(c(1:336)/48, short.predict,type='l',main="Holt-Winters Forecasting",xlab="Jour",ylab="Consomation électrique sur un créneau de 30 minutes en MW")
lines(c(1:336)/48,shortpredictHolt$mean, col = "blue")

plot((short.predict-shortpredictHolt$mean)/short.predict)

rm(list = ls()) #clear list

#automatic installation of required packages
packages <- c("xlsx","calibrate","stargazer","sandwich","lmtest","getopt","CausalGAM","ggplot2","reshape2","xts",
              "lattice","gridExtra","gtable","plm","lfe","lmtest","car","tis","foreign","MASS","quantreg","ggrepel",
              "dplyr","stringr","datasets","rio","psych","systemfit","MatchIt","CRTgeeDR","eurostat","plyr","zoo","ggthemes",
              "robumeta","metafor","dplyr","clubSandwich","Hmisc","metafor","pracma","pkgs","broom","sjPlot", "here", "data.table")

#load packages
library(xlsx)
library(calibrate) 
library (stargazer) 
library(sandwich)
library(lmtest)
library(getopt)
library(CausalGAM)
library(ggplot2)
library(reshape2)
library(xts)
library(lattice)
library(gridExtra)
library(gtable)
library(plm)
library(lfe)
library(lmtest)
library(car)
library(tis)
library(foreign)
library(MASS)
library(quantreg)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(datasets)
library(rio)
library(psych)
library(systemfit)
library(foreign)
library(MatchIt)
library(CRTgeeDR)
library(eurostat)
library(plyr)
library(zoo)
library(ggthemes)
library("robumeta")
library("metafor")
library("dplyr")
library(clubSandwich)
library(Hmisc)
library(metafor)
library(pracma)
library(pkgs)
library(broom)
library(sjPlot)
library(here)
library(data.table)

#load data
dat <- fread(here("data/Coding-ineq-World-Economy-final.csv"))

#calculate the partial correlation coefficient
dat$PartialCorrelationCoefficient <- dat$Tstatistic / (sqrt((dat$Tstatistic^2)+dat$DegreesofFreedom))

dat$PartialCorrelationCoefficient <- dat$PartialCorrelationCoefficient*dat$Transform
dat$Tstatistic <- dat$Tstatistic*dat$Transform

#calculate the standard error of the partial correlation coefficient
dat$StandardErrorPartialCorrelation <- sqrt((1-(dat$PartialCorrelationCoefficient)^2)/dat$DegreesofFreedom)

#Precision
dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation

#Variance 
dat$Variance <- dat$StandardErrorPartialCorrelation^2
#PrecVariance
dat$PrecVariance <- 1 / dat$Variance

dat <- escalc(measure="ZCOR", ri=PartialCorrelationCoefficient, ni=Observations, data=dat) 

#average year
dat$MeanYearData<- (dat$StartYear+dat$EndYear)/2

dat_long <- melt(dat, id=1:52)
dat_DevelopingCountriesOnly <- subset(dat_long, DevelopingCountriesOnly %in% c('1'))
dat_AdvancedCountriesOnly <- subset(dat_long, AdvancedCountriesOnly %in% c('1'))
dat_MixofCountries <- subset(dat_long, MixofCountries %in% c('1'))
dat_TradeGlobalizationOnly <- subset(dat_long, TradeGlobalization %in% c('1'))
dat_FinancialGlobalizationOnly <- subset(dat_long, FinancialGlobalization %in% c('1'))
dat_OverallEconomicGlobalizationOnly <- subset(dat_long, OverallEconomicGlobalization %in% c('1'))

min(dat_long$MeanYearData)
max(dat_long$MeanYearData)
dat_long$MeanYearData <- as.numeric(dat_long$MeanYearData)

#Table 1: average published globalization-spending partial correlaton

#all-set
#unweighted average
uwa <- sum(dat$PartialCorrelationCoefficient, na.rm=TRUE) / 1254
uwa

#(precision-weighted) average
dat$PreSE <- 1/dat$Variance
numerator <- dat$PartialCorrelationCoefficient*dat$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(dat$PreSE, na.rm=TRUE)
wa

#median
median(dat$PartialCorrelationCoefficient, na.rm=TRUE)

res <- rma(yi, vi, data=dat, method="REML") #Random Effecs
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  

fes <- rma(yi, vi, data=dat, method="FE") #Fixed Effecs
fes 

hes <- rma(yi, vi, data=dat, method="HS") #Hunter-Schmidt
hes 
predict(hes, digits=3, transf=transf.ztor)
confint(hes)  

#exclude top and bottom 10%
topbottom <- group_by(dat_long, id) %>%
  mutate(rank = rank(desc(PartialCorrelationCoefficient))) %>%
  filter(rank >=125) %>%
  filter(rank <=1128) %>%
  arrange(rank)

#unweighted average
uwa <- sum(topbottom$PartialCorrelationCoefficient, na.rm=TRUE) / 1004
uwa

#(precision-weighted) average
topbottom$PreSE <- 1/topbottom$Variance
numerator <- topbottom$PartialCorrelationCoefficient*topbottom$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(topbottom$PreSE, na.rm=TRUE)
wa

#median
median(topbottom$PartialCorrelationCoefficient, na.rm=TRUE)

res <- rma(yi, vi, data=topbottom, method="REML") #Random Effecs
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  

fes <- rma(yi, vi, data=topbottom, method="FE") #Fixed Effecs
fes 

hes <- rma(yi, vi, data=topbottom, method="HS") #Hunter-Schmidt
hes 
predict(hes, digits=3, transf=transf.ztor)
confint(hes)

#Trade globalization only
#unweighted average
uwa <- sum(dat_TradeGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE) / 732
uwa

#(precision-weighted) average
dat_TradeGlobalizationOnly$PreSE <- 1/dat_TradeGlobalizationOnly$Variance
numerator <- dat_TradeGlobalizationOnly$PartialCorrelationCoefficient*dat_TradeGlobalizationOnly$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(dat_TradeGlobalizationOnly$PreSE, na.rm=TRUE)
wa

#median
median(dat_TradeGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE)

res <- rma(yi, vi, data=dat_TradeGlobalizationOnly, method="REML") #Random Effecs
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  

fes <- rma(yi, vi, data=dat_TradeGlobalizationOnly, method="FE") #Fixed Effecs
fes 

hes <- rma(yi, vi, data=dat_TradeGlobalizationOnly, method="HS") #Hunter-Schmidt
hes 
predict(hes, digits=3, transf=transf.ztor)
confint(hes)

#Financial globalization only
#unweighted average
uwa <- sum(dat_FinancialGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE) / 473
uwa

#(precision-weighted) average
dat_FinancialGlobalizationOnly$PreSE <- 1/dat_FinancialGlobalizationOnly$Variance
numerator <- dat_FinancialGlobalizationOnly$PartialCorrelationCoefficient*dat_FinancialGlobalizationOnly$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(dat_FinancialGlobalizationOnly$PreSE, na.rm=TRUE)
wa

#median
median(dat_FinancialGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE)

res <- rma(yi, vi, data=dat_FinancialGlobalizationOnly, method="REML") #Random Effecs
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  

fes <- rma(yi, vi, data=dat_FinancialGlobalizationOnly, method="FE") #Fixed Effecs
fes 

hes <- rma(yi, vi, data=dat_FinancialGlobalizationOnly, method="HS") #Hunter-Schmidt
hes 
predict(hes, digits=3, transf=transf.ztor)
confint(hes)

#Overall economic globalization only
#unweighted average
uwa <- sum(dat_OverallEconomicGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE) / 49
uwa

#(precision-weighted) average
dat_OverallEconomicGlobalizationOnly$PreSE <- 1/dat_OverallEconomicGlobalizationOnly$Variance
numerator <- dat_OverallEconomicGlobalizationOnly$PartialCorrelationCoefficient*dat_OverallEconomicGlobalizationOnly$PreSE
wa <- sum(numerator, na.rm=TRUE)/sum(dat_OverallEconomicGlobalizationOnly$PreSE, na.rm=TRUE)
wa

#median
median(dat_OverallEconomicGlobalizationOnly$PartialCorrelationCoefficient, na.rm=TRUE)

res <- rma(yi, vi, data=dat_OverallEconomicGlobalizationOnly, method="REML") #Random Effecs
res 
predict(res, digits=3, transf=transf.ztor)
confint(res)  

fes <- rma(yi, vi, data=dat_OverallEconomicGlobalizationOnly, method="FE") #Fixed Effecs
fes 

hes <- rma(yi, vi, data=dat_OverallEconomicGlobalizationOnly, method="HS") #Hunter-Schmidt
hes 
predict(hes, digits=3, transf=transf.ztor)
confint(hes)

#Figure 1

partialvector <- dat_long$PartialCorrelationCoefficient
h<-hist(partialvector, breaks=10, col="red", xlab="partial correlation coefficient", 
        main="Distribution of globalization-spending\n partial correlation coefficients") 
xfit<-seq(min(partialvector),max(partialvector),length=40) 
yfit<-dnorm(xfit,mean=mean(partialvector),sd=sd(partialvector)) 
yfit <- yfit*diff(h$mids[1:2])*length(partialvector) 
lines(xfit, yfit, col="blue", lwd=2)

#kernel density plot
d <- density(dat_long$PartialCorrelationCoefficient) # returns the density data 
plot(d) # plots the results

hist(partialvector,breaks = 10, freq=F,main="Distribution of partial correlations:\n globalization and income inequality\n (n=1254)",xlab="partial correlation coefficient\n economic globalization-income inequality",ylab="density", ylim=c(0,3), xlim=c(-1,1))
lines(density(partialvector), col="red", lwd=2) 
density(partialvector)

dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation

#distributional statistics
max(dat_long$PartialCorrelationCoefficient)
min(dat_long$PartialCorrelationCoefficient)
sd(dat_long$PartialCorrelationCoefficient)

#Table 3
#Multi-variate meta-regression results

#column (1)
WLS <- lm(PartialCorrelationCoefficient ~ TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, weights=PrecVariance, data=dat_long)
summary(WLS)
coef_test(WLS, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#column (2)
#Random effects panel
RE <- rma(PartialCorrelationCoefficient, Variance, mods = ~ TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, method="REML", weights=PrecSE, data=dat_long) 
summary(RE)

coef_test(RE, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#column (3)
#Robust regression
Robust_app <- rlm(PartialCorrelationCoefficient ~  TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, data=dat_long)
summary(Robust_app)

#column (4)
#Fisher's z
pubbias_2_var_gts_Fisher <- lm(yi ~ TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, weights=PrecVariance, data=dat_long)
summary(pubbias_2_var_gts_Fisher)

coef_test(pubbias_2_var_gts_Fisher, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#stargazer table
#WLS with standard errors clustered at the study level
ses.WLS_app <- list(coef_test(WLS, vcov = "CR1", 
                               cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.WLS_app <- list(coef_test(WLS, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,3]) 


#WLS (Fisher's z) with standard errors clustered at the study level
ses.WLS.Fisher_app <- list(coef_test(pubbias_2_var_gts_Fisher, vcov = "CR1", 
                                     cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.WLS.Fisher_app <- list(coef_test(pubbias_2_var_gts_Fisher, vcov = "CR1", 
                                       cluster = dat_long$id, test = "naive-t")[,3]) 

#Robust regression with standard errors clustered at the study level
ses.Robust_app <- list(coeftest(Robust_app,vcov=NeweyWest(Robust_app, verbose=T))[,2]) 
tvals.Robust_app <- list(coeftest(Robust_app,vcov=NeweyWest(Robust_app, verbose=T))[,3]) 
tvals.Robust_app

#note that column (2) is not integrated in this stargazer output because of compatibility problems; i.e. the stargazer code only shows reesults from Table 3 for the columns (1), (3) and (4). Results for column (2) are available from the RE object
stargazer(WLS, Robust_app, pubbias_2_var_gts_Fisher, t=list(unlist(tvals.WLS_app), unlist(tvals.Robust_app), unlist(tvals.WLS.Fisher_app)), se=list(unlist(ses.WLS_app), unlist(ses.Robust_app), unlist(ses.WLS.Fisher_app)))

####
#Publication bias: funnel plot

#Figure 2

plot_funnel <- ggplot(data=dat,
                      aes(x=PartialCorrelationCoefficient, y=PrecSE)) +
  geom_point() +
  xlab("Partial correlation coefficient") +
  ylab("Precision (1 / Standard error)") +
  ggtitle("Funnel plot of globalization-inequality\n partial correlations (n=1254)")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  geom_vline(xintercept=0, colour="black", linetype=2)+
  theme(legend.text = element_text(colour="black", size = 6))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.title.x=element_text(size=14)) +
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))
plot_funnel

#normalize impact factor
dat_long$JournalImpactFactor <- as.numeric(dat_long$JournalImpactFactor)
dat_long$MaxImpactFactor <- max(dat_long$JournalImpactFactor)
dat_long$MaxImpactFactor <- as.numeric(dat_long$MaxImpactFactor)
dat_long$NormalizedImpactFactor <- dat_long$JournalImpactFactor / max(dat_long$JournalImpactFactor)

#Table 4
#column (1)
pubbias_1 <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_long)
summary(pubbias_1)

coef_test(pubbias_1, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#column(2)
pubbias_2 <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, weights=PrecVariance, data=dat_long)
summary(pubbias_2)

coef_test(pubbias_2, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#column (2)
pubbias_3 <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + NormalizedImpactFactor + TopIncomeShare + BottomIncomeShare + Theil + IncomeShareRatio + HighLowSkilled + OtherIneqVar + FinancialGlobalization + OverallEconomicGlobalization + DevelopingCountriesOnly + MixofCountries  + Primary + Crossauthor + Prior + Technology + Education, weights=PrecVariance, data=dat_long)
summary(pubbias_3)

coef_test(pubbias_3, vcov = "CR1", 
          cluster = dat_long$id, test = "naive-t")

#values for stargazer table
ses.pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", 
                          cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", 
                            cluster = dat_long$id, test = "naive-t")[,3]) 

ses.pubbias_2 <- list(coef_test(pubbias_2, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_2 <- list(coef_test(pubbias_2, vcov = "CR1", 
                                  cluster = dat_long$id, test = "naive-t")[,3]) 

ses.pubbias_3 <- list(coef_test(pubbias_3, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_3 <- list(coef_test(pubbias_3, vcov = "CR1", 
                                  cluster = dat_long$id, test = "naive-t")[,3]) 

ses.pubbias_4 <- list(coef_test(pubbias_4, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_4 <- list(coef_test(pubbias_4, vcov = "CR1", 
                                  cluster = dat_long$id, test = "naive-t")[,3]) 

ses.pubbias_5 <- list(coef_test(pubbias_5, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_5 <- list(coef_test(pubbias_5, vcov = "CR1", 
                                  cluster = dat_long$id, test = "naive-t")[,3]) 

ses.pubbias_6 <- list(coef_test(pubbias_6, vcov = "CR1", 
                                cluster = dat_long$id, test = "naive-t")[,2]) 
tvals.pubbias_6 <- list(coef_test(pubbias_6, vcov = "CR1", 
                                  cluster = dat_long$id, test = "naive-t")[,3]) 


stargazer(pubbias_1, pubbias_2, pubbias_3, pubbias_4, pubbias_5, pubbias_6, t=list(unlist(tvals.pubbias_1), unlist(tvals.pubbias_2), unlist(tvals.pubbias_3), unlist(tvals.pubbias_4), unlist(tvals.pubbias_5), unlist(tvals.pubbias_6)), se=list(unlist(ses.pubbias_1), unlist(ses.pubbias_2),unlist(ses.pubbias_3),unlist(ses.pubbias_4),unlist(ses.pubbias_5),unlist(ses.pubbias_6)))

#summary statistics for Table 2
mean(dat_long$PartialCorrelationCoefficient)
sd(dat_long$PartialCorrelationCoefficient)

mean(dat_long$Gini)
sd(dat_long$Gini)

mean(dat_long$TopIncomeShare)
sd(dat_long$TopIncomeShare)

mean(dat_long$BottomIncomeShare)
sd(dat_long$BottomIncomeShare)

mean(dat_long$IncomeShareRatio)
sd(dat_long$IncomeShareRatio)

mean(dat_long$Theil)
sd(dat_long$Theil)

mean(dat_long$HighLowSkilled)
sd(dat_long$HighLowSkilled)

mean(dat_long$OtherIneqVar)
sd(dat_long$OtherIneqVar)

mean(dat_long$TradeGlobalization)
sd(dat_long$TradeGlobalization)

mean(dat_long$FinancialGlobalization)
sd(dat_long$FinancialGlobalization)

mean(dat_long$OverallEconomicGlobalization)
sd(dat_long$OverallEconomicGlobalization)

mean(dat_long$CrossSection)
sd(dat_long$CrossSection)

mean(dat_long$DataHyperGlobalizationOnly)
sd(dat_long$DataHyperGlobalizationOnly)

mean(dat_long$AdvancedCountriesOnly)
sd(dat_long$AdvancedCountriesOnly)

mean(dat_long$DevelopingCountriesOnly)
sd(dat_long$DevelopingCountriesOnly)

mean(dat_long$MixofCountries)
sd(dat_long$MixofCountries)

mean(dat_long$CountryFixedEffects)
sd(dat_long$CountryFixedEffects)

mean(dat_long$NonOLS)
sd(dat_long$NonOLS)

mean(dat_long$StandardErrorPartialCorrelation)
sd(dat_long$StandardErrorPartialCorrelation)

mean(dat_long$NormalizedImpactFactor)
sd(dat_long$NormalizedImpactFactor)

mean(dat_long$DevelopmentJournal)
sd(dat_long$DevelopmentJournal)

mean(dat_long$Primary)
sd(dat_long$Primary)

mean(dat_long$Crossauthor)
sd(dat_long$Crossauthor)

mean(dat_long$Prior)
sd(dat_long$Prior)

mean(dat_long$Technology)
sd(dat_long$Technology)

mean(dat_long$GDPgrowth)
sd(dat_long$GDPgrowth)

mean(dat_long$Unemployment)
sd(dat_long$Unemployment)

mean(dat_long$IncomeLevel)
sd(dat_long$IncomeLevel)

mean(dat_long$Population)
sd(dat_long$Population)

mean(dat_long$Education)
sd(dat_long$Education)

mean(dat_long$TradeUnion)
sd(dat_long$TradeUnion)

mean(dat_long$SocialSpending)
sd(dat_long$SocialSpending)

mean(dat_long$Democracy)
sd(dat_long$Democracy)

mean(dat_long$PartisanPolitics)
sd(dat_long$PartisanPolitics)

#distributional statistics
max(dat_long$PartialCorrelationCoefficient)
min(dat_long$PartialCorrelationCoefficient)
sd(dat_long$PartialCorrelationCoefficient)

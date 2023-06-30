#This script conducts the driver analysis using 2013-2019 survey data


rm(list=ls())

#LOAD LIBRARY FUNCTIONS ... 
source("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/scripts/Functions_Juveniles.R")
#source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/GIS_functions.R")

#Calculate Difference in range of values
rg<-function(m) {
  max(m, na.rm=TRUE) - min(m, na.rm=TRUE)
}


library(forcats)
library(geosphere)
library(rgdal)
library(stringr)
library(mgcv)
detach(package:dplyr) #dplyr has issues if plyr is loaded first
library(dplyr)
library(MASS)
library(performance)
library(see)
library(patchwork)
library(emmeans)
library(rcompanion)
library(survey)
library(gridExtra)
library(ggpubr)
library(lemon)
library(arm)
library(purrr)
library(tibble)
library(car)
library(extrafont)
library(remotes)
library(corrplot)
library(jtools)


# remotes::install_version("Rttf2pt1", version = "1.3.8")
# extrafont::font_import()
#loadfonts(device="win")


setwd("T:/Benthic/Projects/Juvenile Project")


#LOAD DATA
df<-read.csv("T:/Benthic/Projects/Juvenile Project/Data/JuvDen_Pred_SITE_AllYears.csv")#Combined juvenile delta density and all predictors
cover_sec<-read.csv("T:/Benthic/Projects/Juvenile Project/BenthicCover_JuvenileProject_Tier1_SECTOR.csv")

#remove columns
df<-subset(df,select=c(DATE_,OBS_YEAR,REGION,ISLAND,SEC_NAME,DEPTH_BIN,REEF_ZONE,STRATANAME,SITE,n,NH,sw,TRANSECTAREA_j,JuvColCount,JuvColDen,
                       LATITUDE,LONGITUDE,Depth_Median,CORAL,CORALst,CCA,SAND_RUB,TURF,EMA_MA, YearSinceDHW4,DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR01,DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR03,
                       DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR05,DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR10,DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR10YR01,
                       DHW.Np10y_Major_Degree_Heating_Weeks_CRW_Daily_YR10,mean_SST_CRW_CoralTemp_Daily_YR10,sd_SST_CRW_CoralTemp_Daily_YR10,mean_weekly_range_SST_CRW_CoralTemp_Daily_YR10,mean_biweekly_range_SST_CRW_CoralTemp_Daily_YR10,mean_Chlorophyll_A_ESA_OC_CCI_8Day_YR10,sd_Chlorophyll_A_ESA_OC_CCI_8Day_YR10,
                       mean_annual_range_Chlorophyll_A_ESA_OC_CCI_8Day_YR10, WavePower,HumanDen,HerbivoreBio))

#Combine site level data with sector cover data
cover_sec$OBS_YEAR<-cover_sec$ANALYSIS_YEAR
df<-left_join(df,cover_sec)

#Convert NH (number of 50x50m grids to area)
df$Area<-df$NH*250


df$CORAL_sec<-ifelse(df$SEC_NAME=="Baker",23.62746339,df$CORAL_sec)#no 2018 benthic cover sites so using fish sector data
df$CoralSec_A<-df$Area*df$CORAL_sec

df<-filter(df,ISLAND !="Guguan") #only 1 year

#reformat dates
df$DATE_<-lubridate::ymd(df$DATE_)
#Generate table of date ranges and n for each island
df$OBS_MONTH<-lubridate::month(df$DATE_,label=TRUE)
df$OBS_DAY<-day(df$DATE_)

meta<-df %>%
  group_by(REGION,ISLAND,OBS_YEAR) %>%
  summarize(MinDate=min(DATE_),MaxDate=max(DATE_),n=length(n))

meta<-df %>%
  group_by(REGION,ISLAND,OBS_YEAR) %>%
  summarize(MinMonth=min(OBS_MONTH),MaxMonth=max(OBS_MONTH),
            MinDay=min(OBS_DAY),MaxDay=max(OBS_DAY),
            n=length(n))

meta$Month<-ifelse(meta$MinMonth==meta$MaxMonth,as.character(meta$MaxMonth),paste(meta$MinMonth,meta$MaxMonth,sep=" - "))
meta$DateRange<-paste(meta$Month,meta$OBS_YEAR,sep=" ")


meta<-meta[,c("REGION","ISLAND","OBS_YEAR","DateRange","n")]
head(meta)
meta<-meta %>% mutate(T1_T2= dplyr::recode(OBS_YEAR,
                                           `2013`="T1",     
                                           `2014`="T1",
                                                `2015`="T1",
                                                `2016`="T2",
                                                `2017`="T2",
                                                `2018`="T2",
                                                `2019`="T3"))
meta$R_Y<-paste(meta$REGION,meta$OBS_YEAR,sep="-")
meta$T1_T2<-ifelse(meta$R_Y=="NWHI-2017","T3",meta$T1_T2)
View(meta)

coord<-df %>%
  group_by(REGION,ISLAND) %>%
  summarize(Latitude=median(LATITUDE),Longitude=median(LONGITUDE))

wide.date<-meta %>%
  dplyr::select(REGION,ISLAND,T1_T2,DateRange) %>%
  pivot_wider(names_from = T1_T2,values_from = DateRange)

wide.n<-meta %>%
  dplyr::select(REGION,ISLAND,T1_T2,n) %>%
  pivot_wider(names_from = T1_T2,values_from = n)
View(wide.n)

wide<-left_join(wide.date,wide.n,by=c("REGION","ISLAND"))

wide$T1<-paste(wide$T1.x,"(",wide$T1.y,")")
wide$T2<-paste(wide$T2.x,"(",wide$T2.y,")")
wide$T3<-paste(wide$T3.x,"(",wide$T3.y,")")

wide<-left_join(wide,coord,by=c("REGION","ISLAND"))

head(wide)
wide<-wide[,c("REGION","ISLAND","Latitude","Longitude","T1","T2","T3")]

wide$ISLAND <- factor(wide$ISLAND, levels = c("Kure","Pearl_&_Hermes","Lisianski","French_Frigate","Kauai",
                                              "Oahu","Molokai","Maui","Lanai","Kahoolawe","Hawaii","Wake",
                                              "Howland","Baker","Kingman","Palmyra","Jarvis","Saipan","Tinian",
                                              "Rota","Aguijan","Guam","Farallon_de_Pajaros","Maug","Pagan","Asuncion",
                                              "Sarigan","Swains","Tutuila","Ofu_&_Olosega","Tau","Rose"))
wide<- wide[order(wide$ISLAND),];View(wide)


write.csv(wide,file="T:/Benthic/Projects/Juvenile Project/Tables/TableS1.csv")


#Dealing with missing Year Since DHW4 data:
#If a site never experienced a >=4 DHW event then set YearSinceDHW4 to the most recent survey date - 1st recorded DHW data (1/1/1985)
df$DATE_<-lubridate::ymd(df$DATE_)
dhw_start<-lubridate::ymd("1985-01-01")

#how many NA values for year since = 198 (14.3%)
summary(df$YearSinceDHW4) 

#Assign year since DHW 4 event that has NA value the max year for the climatology = 32 years
df$YearSinceDHW4<-ifelse(is.na(df$YearSinceDHW4),32,df$YearSinceDHW4)

head(subset(df,SEC_NAME=="OAH_NORTH")) #check that column was changed correctly

View(df)

#Remove strata that have no corallivore or herbivore biomass data
df.new<-df %>% filter(!is.na(HerbivoreBio)) #Now 1294-dropping 93 sites with NA values

table(df.new$REGION,df.new$OBS_YEAR)

nrow(df.new)

#Convert latitude to absolute value
df.new$LATITUDE<-abs(df.new$LATITUDE)


#Rename Predictors
colnames(df.new)[colnames(df.new)=="DHW.MeanMax_Degree_Heating_Weeks_CRW_Daily_YR10"]<-"MeanDHW10"
colnames(df.new)[colnames(df.new)=="mean_SST_CRW_CoralTemp_Daily_YR10"]<-"MeanSST"
colnames(df.new)[colnames(df.new)=="mean_Chlorophyll_A_ESA_OC_CCI_8Day_YR10"]<-"Meanchla"
colnames(df.new)[colnames(df.new)=="mean_biweekly_range_SST_CRW_CoralTemp_Daily_YR10"]<-"SST_Range"

hist(log10(df.new$HumanDen+0.5))
df.new$logHumanDen<-log10(df.new$HumanDen+0.5)


#Extract predictors and merge with new survey weights dataset
pcols<-c("SITE","CORAL","CoralSec_A","CORALst","CCA","TURF","EMA_MA","SAND_RUB","Depth_Median","LATITUDE",
         "MeanDHW10","Meanchla","SST_Range",
         "MeanSST","WavePower","YearSinceDHW4","logHumanDen","HerbivoreBio")

p<-df.new[,pcols]

#Combine survey weighted juvenile data and predictors
rcols<-c("OBS_YEAR","REGION","SITE","TRANSECTAREA_j","JuvColCount","n","NH","sw")

r<-df.new[,rcols]

nrow(r) #should be 1294
r<-left_join(r,p)
nrow(r);View(r) #should be 1294


#Testing for Multicolinarity
which(colnames(r)=="CORAL")
preds<-r[,9:ncol(r)]
# library(GGally)
# ggpairs(preds)


par(mfrow=c(1,1))
M = cor(preds)
png(width = 750, height = 750, filename = "T:/Benthic/Projects/Juvenile Project/Figures/Drivers/JuvenilePredictorsCorPlot_AllYears.png")
corrplot(M, method = 'number')
dev.off()

#We used a correlation coefficent of 0.65 which is fairly conservative. 

#Confirmed with VIF - a priori cut off 3, but all less than 2.
fit1 <- lm(JuvColDen ~ CORAL + CoralSec_A +  CCA +  EMA_MA + SAND_RUB + Depth_Median +  
             MeanDHW10 + Meanchla + MeanSST +
             WavePower + YearSinceDHW4 + logHumanDen +HerbivoreBio, data = df.new)

car::vif(fit1)

#Turf and CCA correlated, latitude and Mean SST correlated, MeanMaxDHW and SST Range correlated
#Dropping turf, latitutde, and SST range

preds <- scale(preds, center = T, scale = T);colnames(preds)<-paste("scaled",colnames(preds),sep="_")

final.df<-cbind(df.new,preds)



# Quick plots of top predictors -------------------------------------------

par(mfrow=c(2,2))
plot(final.df$JuvColDen~final.df$MeanDHW10)
plot(final.df$JuvColDen~final.df$Depth_Median)
#plot(final.df$JuvColDen~final.df$LATITUDE)
plot(final.df$JuvColDen~final.df$Meanchla) 

plot(final.df$JuvColDen~final.df$scaled_CORAL) 
plot(final.df$JuvColDen~final.df$scaled_CCA) 
plot(final.df$JuvColDen~final.df$scaled_EMA_MA) 
plot(final.df$JuvColDen~final.df$scaled_SAND_RUB) 

plot(final.df$JuvColDen~final.df$YearSinceDHW4)
plot(final.df$JuvColDen~final.df$scaled_SST_Range)
plot(final.df$JuvColDen~final.df$scaled_HerbivoreBio)

par(mfrow=c(1,1))
plot(final.df$JuvColDen~final.df$YearSinceDHW4)
plot(final.df$JuvColDen~final.df$scaled_CoralSec_A) 
plot(final.df$JuvColDen~final.df$scaled_logHumanDen) 

head(final.df)

#There are sites that never experienced a DHW 4 heating event = NA. Set these NAs to max of climatolgoical time series = 32 years
final.df$YearSinceDHW4<-ifelse(is.na(final.df$YearSinceDHW4),32,final.df$YearSinceDHW4)


# #Backwards Model selection with Wald Tests (similar to LRTs) ------------
data.cols<-c("OBS_YEAR","REGION","ISLAND","SEC_NAME","STRATANAME","SITE","TRANSECTAREA_j","JuvColCount","n","NH","sw","SITE","CORAL","CoralSec_A","CCA","EMA_MA","SAND_RUB","Depth_Median",
                    "MeanDHW10","Meanchla","MeanSST","WavePower","YearSinceDHW4","scaled_CORAL","scaled_CoralSec_A","scaled_CCA","scaled_EMA_MA","scaled_SAND_RUB","scaled_Depth_Median",
"scaled_MeanDHW10","scaled_Meanchla","scaled_MeanSST","scaled_WavePower","scaled_YearSinceDHW4","scaled_logHumanDen","logHumanDen","scaled_HerbivoreBio","HerbivoreBio")

final.df<-final.df[,data.cols]


#StRS design
#concantate all nested variables into 1 column otherwise you will overestimate the variance with this method only accounting for the first variable
final.df$Strat_conc<-paste(final.df$OBS_YEAR, final.df$REGION,final.df$ISLAND,final.df$STRATANAME,sep = "_")

des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=final.df)


# Testing for polynomial relationships -----------------------------------------------------

#Testing polynomial relationships with DEPTH
d_poly3<-svyglm(JuvColCount ~  
                  poly(scaled_Depth_Median,3),
                design=des, family="poisson",offset=log(TRANSECTAREA_j))

d_poly2<-svyglm(JuvColCount ~  
                    poly(scaled_Depth_Median,2),
                    design=des, family="poisson",offset=log(TRANSECTAREA_j))
d<-svyglm(JuvColCount ~  
            scaled_Depth_Median,
          design=des, family="poisson",offset=log(TRANSECTAREA_j))

anova(d,d_poly2) 
anova(d_poly3,d_poly2) 

AIC(d_poly3)
AIC(d_poly2)
AIC(d)
#2nd order polynomial is best fit

#Testing polynomial relationships with site-level CORAL COVER
d_poly3<-svyglm(JuvColCount ~  
                  poly(scaled_CORAL,3),
                design=des, family="poisson",offset=log(TRANSECTAREA_j))

d_poly2<-svyglm(JuvColCount ~  
                  poly(scaled_CORAL,2),
                design=des, family="poisson",offset=log(TRANSECTAREA_j))
d<-svyglm(JuvColCount ~  
            scaled_CORAL,
          design=des, family="poisson",offset=log(TRANSECTAREA_j))

anova(d,d_poly2) 
anova(d_poly3,d_poly2) 

AIC(d_poly3)
AIC(d_poly2)
AIC(d)


#Tested other variables- no polynomial relationships

#testing influence of high herbivore biomass stratum
drop.hi<-subset(final.df,HerbivoreBio <90)
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=drop.hi)


d<-svyglm(JuvColCount ~  
            scaled_HerbivoreBio,
          design=des, family="poisson",offset=log(TRANSECTAREA_j))


#Tested dropping sites with outlier herbivore biomass over 90 g/m2. - herbivore still retained as signficant in best fit model. 
df.d<-drop.hi
df.d$scaled_HerbivoreBio<- seq(min(drop.hi$scaled_HerbivoreBio),max(drop.hi$scaled_HerbivoreBio),
                             by=round(rg(drop.hi$scaled_HerbivoreBio),5)/nrow(drop.hi))


p <- predict(d, newdata = df.d, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata<-cbind(df.d,p)
newdata$Predict.lwr <- newdata$Predicted_Juv - 1.96 * newdata$SE_Juv # confidence interval upper bound
newdata$Predict.upr <- newdata$Predicted_Juv + 1.96 * newdata$SE_Juv # confidence interval lower bound
head(newdata)


att <- attributes(scale(final.df$HerbivoreBio))
mylabels <- seq(0,95,10)
mybreaks <- scale(mylabels, att$`scaled:center`, att$`scaled:scale`)[,1]

#Plot
ggplot(newdata, aes(x = scaled_HerbivoreBio, y = Predicted_Juv)) +
  geom_line() +
  geom_ribbon(data = newdata,
              aes(ymin = Predict.lwr, ymax = Predict.upr),
              alpha = 0.1)+
  geom_rug(data=newdata,mapping=aes(x=scaled_HerbivoreBio,y=0,color=REGION))+
  ylab("Predicted Juvenile Abudance") +
  xlab("Herbivore Biomass (g/m2)")+ 
  scale_x_continuous(labels=mylabels,breaks=mybreaks)


##### Model Selection ####
#Remove French Frigate Shoals to test whether sector coral cover area is retained in best fit model.
#Note- analyses reveal that sector-level coral cover x area is still signficant even after removing FFS which has anolomously high coral cover area.
#Backwards model selection order is consistent with model including FFS.
# no.ffs<-subset(final.df,ISLAND != "French Frigate")
# no.ffs$Strat_conc<-paste(no.ffs$OBS_YEAR, no.ffs$REGION,no.ffs$ISLAND,no.ffs$STRATANAME,sep = "_")
# des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=no.ffs)

final.df$Strat_conc<-paste(final.df$OBS_YEAR, final.df$REGION,final.df$ISLAND,final.df$STRATANAME,sep = "_")

des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=final.df)

#Global model - INCLUDING HERBIVORES
global.mod1<-svyglm(JuvColCount ~
                      scaled_CORAL*scaled_MeanDHW10+ 
                      scaled_CCA*poly(scaled_Depth_Median,2,raw=TRUE)+
                      scaled_CoralSec_A*scaled_MeanDHW10 +
                      scaled_EMA_MA*scaled_MeanDHW10 +
                      scaled_SAND_RUB*scaled_MeanDHW10 +
                      scaled_HerbivoreBio*scaled_MeanDHW10 +
                      poly(scaled_Depth_Median,2,raw=TRUE)*scaled_MeanDHW10 +
                      scaled_Meanchla +
                      scaled_MeanSST +
                      scaled_WavePower*scaled_MeanDHW10+
                      scaled_YearSinceDHW4*scaled_MeanDHW10+
                      scaled_logHumanDen*scaled_MeanDHW10,
                    design=des, family="poisson",offset=log(TRANSECTAREA_j)) #also tried quasipoisson -no change in model??

summary(global.mod1)

#Only option to generate a R2 like metric for these kinds of models
cor(global.mod1$y, fitted(global.mod1))^2
AIC(global.mod1)

summ(global.mod1) #jtools


#Backwards model selection
RED.MOD1 <- update(global.mod1, .~. -scaled_MeanDHW10:scaled_EMA_MA) #drop term
anova(global.mod1, RED.MOD1,method="Wald") #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD1)


RED.MOD2 <- update(RED.MOD1, .~. -scaled_MeanDHW10:scaled_SAND_RUB) #drop  term
anova(RED.MOD1, RED.MOD2) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD2)

RED.MOD3 <- update(RED.MOD2, .~. -scaled_CCA:poly(scaled_Depth_Median, 2, raw = TRUE)) #drop term
anova(RED.MOD2, RED.MOD3,test = "Chisq") #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD3)

RED.MOD4 <- update(RED.MOD3, .~. -scaled_MeanDHW10:scaled_HerbivoreBio) #drop term
anova(RED.MOD3, RED.MOD4) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD4)

RED.MOD5 <- update(RED.MOD4, .~. -scaled_Meanchla) #drop term
anova(RED.MOD4, RED.MOD5) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD5)

RED.MOD6 <- update(RED.MOD5, .~. -scaled_MeanDHW10:poly(scaled_Depth_Median, 2, raw = TRUE)) #drop term
anova(RED.MOD5, RED.MOD6) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD6)

RED.MOD7 <- update(RED.MOD6, .~. -scaled_CCA) #drop term
anova(RED.MOD6, RED.MOD7) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD7)

RED.MOD8 <- update(RED.MOD7, .~. -) #drop term
anova(RED.MOD7, RED.MOD8) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD8)

RED.MOD9 <- update(RED.MOD8, .~. -scaled_CORAL:scaled_MeanDHW10) #drop term
anova(RED.MOD8, RED.MOD9) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD9)

RED.MOD10 <- update(RED.MOD9, .~. -caled_HerbivoreBio) #drop term
anova(RED.MOD9, RED.MOD10) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD10)

RED.MOD11 <- update(RED.MOD10, .~. -scaled_HerbivoreBio) #drop term
anova(RED.MOD10, RED.MOD11) #LRT --> move forward w/ whichever model keeps/removes term
summary(RED.MOD11)

AIC(RED.MOD9)
AIC(RED.MOD10)
AIC(RED.MOD11)


best.mod<-RED.MOD10
summary(best.mod)

#Caculate McFadden's R2
summ(best.mod) #jtools


psrsq<-function(object, method=c("Cox-Snell","Nagelkerke"),...){
  UseMethod("psrsq",object)
}

psrsq.glm<-function(object, method=c("Cox-Snell","Nagelkerke"),...){
  nullmodel<-update(object,.~1)
  method<-match.arg(method)
  ell0<-as.vector(logLik(nullmodel))
  ell1<-as.vector(logLik(object))
  n<-object$df.null+1
  
  mutualinf<-  -2*(ell1-ell0)/n
  r2cs<-1-exp(mutualinf)
  if (method == "Cox-Snell") 
    return(r2cs)
  scaling<-1-exp(2*ell0/n)
  r2cs/scaling
}

psrsq.svyglm<-function(object, method=c("Cox-Snell", "Nagelkerke"),...){
  method<-match.arg(method)
  if (!(object$family$family %in% c("binomial","quasibinomial","poisson","quasipoisson")))
    stop("Only implemented for discrete data")
  w<-weights(object$survey.design,"sampling")
  N<-sum(w)
  n<-sum(object$prior.weights)
  minus2ell0<-object$null.deviance*(N/n)
  minus2ell1<-object$deviance*(N/n)
  mutualinf<-(minus2ell1-minus2ell0)/N
  r2cs<-1-exp(mutualinf)
  if (method =="Cox-Snell") 
    return(r2cs)
  if (any(w<1)) warning("Weights appear to be scaled: rsquared may be wrong")
  scaling<-1-exp(-minus2ell0/N)
  r2cs/scaling
}



psrsq.svyglm(best.mod,method= "Nagelkerke")
psrsq(best.mod,method= "Nagelkerke")







#Model diagnostics
svystdres(best.mod,doplot=TRUE)
#We have overdispersion
res <- residuals(best.mod, type="deviance")
plot(log(predict(best.mod)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)





# PARAMETER ESTIMATES +/- SE ----------------------------------------------

sum <- summary(best.mod)
sum.co <- data.frame(sum$coefficients)
sum.co$Variable <- rownames(sum.co)
sum.co <- data.frame(sum.co, row.names = NULL)
sum.co <- sum.co[ order(abs(sum.co$Estimate), decreasing = T),]
var_ord <- sum.co$Variable


sum.co$lwr.CI <- sum.co$Estimate - 1.96 * sum.co$Std..Error # confidence interval lower bound
sum.co$upr.CI <- sum.co$Estimate + 1.96 * sum.co$Std..Error # confidence interval upper bound
head(sum.co)

sum.co <- sum.co[ which(sum.co$Variable != "(Intercept)"),]


# sum.co$Variable <- factor(sum.co$Variable, levels = var_ord)
# sum.co <- sum.co[order(factor(sum.co$Variable, levels = var_ord)),]
sum.co$Variable_plot <- factor(c("HS_ts x HS_sev",
                                 "Depth",
                                 "Coral Cover^2",
                                 "Coral Cover",
                                 "HS_ts",
                                 "Depth^2",
                                 "Sector-level Coral Cover",
                                 "Human Density",
                                 "Unconsolidated Cover",
                                 "Macroalgae Cover",
                                 "Herbivore Biomass",
                                 "Sector-level Coral Cover x HSsev",
                                 "Wave Power x HSsev",
                                 "Coral Cover^3",
                                 "Heat Stress",
                                 "Wave Power"),
                               levels = c("HS_ts x HS_sev",
                                          "Depth",
                                          "Coral Cover^2",
                                          "Coral Cover",
                                          "HS_ts",
                                          "Depth^2",
                                          "Sector-level Coral Cover",
                                          "Human Density",
                                          "Unconsolidated Cover",
                                          "Macroalgae Cover",
                                          "Herbivore Biomass",
                                          "Sector-level Coral Cover x HSsev",
                                          "Wave Power x HSsev",
                                          "Coral Cover^3",
                                          "Heat Stress",
                                          "Wave Power"))

write.csv(sum.co,file="T:/Benthic/Projects/Juvenile Project/Manuscript/Tables/Table S3_unformatted.csv")


sum.co$Sig <- NA
sum.co <- transform(sum.co,
                    Sig=ifelse(Pr...t..<0.05,"p<0.05","p>0.05"))

var_plot <-
  ggplot(sum.co,aes(x = reorder(Variable_plot,Estimate), y = Estimate)) +
  geom_point(aes(color = Sig),size = 3) +
  geom_hline(yintercept = 0, color = "lightgray") +
  geom_errorbar(aes(ymin = lwr.CI, ymax = upr.CI,color = Sig), width=.2,
                position=position_dodge(.9))  +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.background = element_rect(colour = "black", fill = "white"),
    text = element_text(size = 18)) +
  xlab("") +
  ylab("\nParameter Estimate") +
  scale_y_continuous(limits = c(-0.6,0.7)) +
  scale_x_discrete(limits = rev(levels(sum.co$Variable_plot))) +
  scale_color_manual(values = c("#009E73","black"))

var_plot


setwd("T:/Benthic/Projects/Juvenile Project/Manuscript/Figures")
pdf(width = 12, height = 12, file = "Figure 3.pdf")
var_plot
dev.off()


# function to predict bleaching and create a plot based on variable of interest
Predictplot <- function(mod=best.mod,dat=newdata, us_pred="CORAL",predictor="scaled_CORAL", predictor_name,sigcol="black",bks=2){
  dat$s_X<-dat[,predictor] #scaled predictor
  dat$X<-dat[,us_pred] #unscaled predictor
  
  p <- predict(mod, newdata = dat, type = "response",se.fit=TRUE)
  p<-as.data.frame(p)
  colnames(p)<-c("Predicted_Juv","SE_Juv")
  dat<-cbind(dat,p)
  dat$Predict.lwr <- dat$Predicted_Juv - 1.96 * dat$SE_Juv # confidence interval upper bound
  dat$Predict.upr <- dat$Predicted_Juv + 1.96 * dat$SE_Juv # confidence interval lower bound
  head(dat)
  
  mx_val<-max(dat$X, na.rm = TRUE)
  
  #Unscaling predictor to plot on x axis
  att <- attributes(scale(dat$X))
  mylabels <- seq(0,mx_val,bks)
  mybreaks <- scale(mylabels, att$`scaled:center`, att$`scaled:scale`)[,1]
  
  
  #Try mapping geom_rug(dat2,aes(x=s_X)) - dat2= final.df
  
  #Plot
  plot1<-ggplot(data=dat, aes(x = s_X, y = Predicted_Juv)) + 
    geom_line(color=sigcol,size=1) +
    geom_ribbon(data = dat,
                aes(ymin = Predict.lwr, ymax = Predict.upr),
                alpha = 0.1)+
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title = element_text(face = "bold"),
      text = element_text(size = 18),
      panel.grid = element_blank()
    ) +
    ylab("Predicted Juvenile Abudance") +
    xlab(predictor_name)  +
    scale_x_continuous(labels = mylabels,breaks=mybreaks) 
  
  return(plot1)
  
}


#10yr Meam Max DHW
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-seq(min(final.df$scaled_MeanDHW10),max(final.df$scaled_MeanDHW10),
                              by=round(rg(final.df$scaled_MeanDHW10),5)/nrow(final.df))
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4,na.rm=T)

hs.plot<-Predictplot(best.mod,dat=newdata,"MeanDHW10","scaled_MeanDHW10",expression(bold('Mean Max '^o*'C-weeks')),"black",2)+
  geom_rug(data=final.df,mapping=aes(x=scaled_MeanDHW10,y=0))


#Depth
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- seq(min(final.df$scaled_Depth_Median),max(final.df$scaled_Depth_Median),
                                  by=round(rg(final.df$scaled_Depth_Median),4)/nrow(final.df))
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)


depth.plot<-Predictplot(best.mod,newdata,"Depth_Median","scaled_Depth_Median","Median Depth (m)","#009E73",2)+
  geom_rug(data=final.df,mapping=aes(x=scaled_Depth_Median,y=0))


#Coral Cover
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- seq(min(final.df$scaled_CORAL),max(final.df$scaled_CORAL),
                            by=round(rg(final.df$scaled_CORAL),5)/nrow(final.df))
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

coral.plot<-Predictplot(best.mod,newdata,"CORAL","scaled_CORAL","% Coral Cover","#009E73",10)+
  geom_rug(data=final.df,mapping=aes(x=scaled_CORAL,y=0))

#Sector-level Coral Cover
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- seq(min(final.df$scaled_CoralSec_A),max(final.df$scaled_CoralSec_A),
                                 by=round(rg(final.df$scaled_CoralSec_A),4)/nrow(final.df))
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

coralsec.plot<-Predictplot(best.mod,newdata,"CoralSec_A","scaled_CoralSec_A","Sector % Coral Cover x Area (ha)","#009E73",50000000)+
  geom_rug(data=final.df,mapping=aes(x=scaled_CoralSec_A,y=0))


#Sand/Rubble Cover
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- seq(min(final.df$scaled_SAND_RUB),max(final.df$scaled_SAND_RUB),
                               by=round(rg(final.df$scaled_SAND_RUB),3)/nrow(final.df))
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

sandrub.plot<-Predictplot(best.mod,newdata,"SAND_RUB","scaled_SAND_RUB","% Unconsolidated Cover","#009E73",10)+
  geom_rug(data=final.df,mapping=aes(x=scaled_SAND_RUB,y=0))

#Macroalgae Cover
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- seq(min(final.df$scaled_EMA_MA),max(final.df$scaled_EMA_MA),
                             by=round(rg(final.df$scaled_EMA_MA),7)/nrow(final.df))
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

ma.plot<-Predictplot(best.mod,newdata,"EMA_MA","scaled_EMA_MA","% Macroalgae Cover","#009E73",10)+
  geom_rug(data=final.df,mapping=aes(x=scaled_EMA_MA,y=0))

#Herbivore Biomass
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- seq(min(final.df$scaled_HerbivoreBio),max(final.df$scaled_HerbivoreBio),
                             by=round(rg(final.df$scaled_HerbivoreBio),7)/nrow(final.df))
newdata$scaled_WavePower <-mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

herb.plot<-Predictplot(best.mod,newdata,"HerbivoreBio","scaled_HerbivoreBio",expression(bold(Herbivore~Biomass~(g~m^{-2}))),"#009E73",10)+
  geom_rug(data=final.df,mapping=aes(x=scaled_HerbivoreBio,y=0))


#Wave Power
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- seq(min(final.df$scaled_WavePower),max(final.df$scaled_WavePower),
                                by=round(rg(final.df$scaled_WavePower),6)/nrow(final.df))
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

wave.plot<-Predictplot(best.mod,newdata,"WavePower","scaled_WavePower",expression(bold(Wave~Power~(kWh~m^{-1}))),"black", 60000)+
  geom_rug(data=final.df,mapping=aes(x=scaled_WavePower,y=0))


#LogHumanDensity
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- seq(min(final.df$scaled_logHumanDen),max(final.df$scaled_logHumanDen),
                                  by=round(rg(final.df$scaled_logHumanDen),3)/nrow(final.df))
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-mean(final.df$scaled_YearSinceDHW4)

human.plot<-Predictplot(best.mod,newdata,"logHumanDen","scaled_logHumanDen",expression(bold(Log~Human~Density~(km^{-2}))),"#009E73",0.5)+
  geom_rug(data=final.df,mapping=aes(x=scaled_logHumanDen,y=0))

#Year Since DHW4 
newdata <- final.df
newdata$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata$scaled_CORAL <- mean(final.df$scaled_CORAL)
newdata$scaled_CoralSec_A <- mean(final.df$scaled_CoralSec_A)
newdata$scaled_SAND_RUB <- mean(final.df$scaled_SAND_RUB)
newdata$scaled_EMA_MA <- mean(final.df$scaled_EMA_MA)
newdata$scaled_HerbivoreBio <- mean(final.df$scaled_HerbivoreBio)
newdata$scaled_WavePower <- mean(final.df$scaled_WavePower)
newdata$scaled_Depth_Median<- mean(final.df$scaled_Depth_Median)
newdata$scaled_logHumanDen <- mean(final.df$scaled_logHumanDen)
newdata$scaled_MeanDHW10<-mean(final.df$scaled_MeanDHW10)
newdata$scaled_YearSinceDHW4<-seq(min(final.df$scaled_YearSinceDHW4),max(final.df$scaled_YearSinceDHW4),
                                  by=round(rg(final.df$scaled_YearSinceDHW4),4)/nrow(final.df))

tsdhw.plot<-Predictplot(best.mod,newdata,"YearSinceDHW4","scaled_YearSinceDHW4","Years Since Heat Stress Event","#009E73",2)+
  geom_rug(data=final.df,mapping=aes(x=scaled_YearSinceDHW4,y=0))


# save full plot
setwd("T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/")
# ytitle <- text_grob("Predicted Juvenile Colonies/m2", size = 18, face = "bold", rot = 90)

ytitle <- text_grob(expression(bold(paste("Predicted Juvenile Colonies",m^-2))), size = 18, face = "bold", rot = 90)
pdf(width = 14.5, height = 16,file = "Figure 4.pdf") #Then convert to jpeg in Adobe for best quality

grid.arrange(arrangeGrob(depth.plot + ggtitle("A)"),
                         coral.plot + ggtitle("B)"),
                         tsdhw.plot + ggtitle("C)"), 
                         coralsec.plot + ggtitle("D)"),
                         human.plot + ggtitle("E)"), 
                         sandrub.plot + ggtitle("F)"), 
                         ma.plot + ggtitle("G)"),
                         herb.plot + ggtitle("H)"),
                         hs.plot + ggtitle("I)"),
                         wave.plot + ggtitle("J)"),
                         nrow = 4), 
             nrow = 2, heights = c(10,1),
             left = ytitle)
dev.off()


# Visualizing Interactions----------------------------------------

#Heat Stress Severity x year since heat stress event
r <- subset(final.df,YearSinceDHW4<=3)
newdata1<-r
newdata1$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata1$scaled_CORAL <- mean(r$scaled_CORAL)
newdata1$scaled_CoralSec_A <- mean(r$scaled_CoralSec_A)
newdata1$scaled_SAND_RUB <- mean(r$scaled_SAND_RUB)
newdata1$scaled_EMA_MA <- mean(r$scaled_EMA_MA)
newdata1$scaled_HerbivoreBio <- mean(r$scaled_HerbivoreBio)
newdata1$scaled_WavePower <- mean(r$scaled_WavePower)
newdata1$scaled_Depth_Median<- mean(r$scaled_Depth_Median)
newdata1$scaled_logHumanDen <- mean(r$scaled_logHumanDen)
newdata1$scaled_MeanDHW10<-seq(min(r$scaled_MeanDHW10),max(r$scaled_MeanDHW10),
                              by=round(rg(r$scaled_MeanDHW10),3)/nrow(r))
newdata1$scaled_YearSinceDHW4<-mean(r$scaled_YearSinceDHW4)


m <- subset(final.df,YearSinceDHW4>3 & YearSinceDHW4 <=10)
newdata2<-m
newdata2$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata2$scaled_CORAL <- mean(m$scaled_CORAL)
newdata2$scaled_CoralSec_A <- mean(m$scaled_CoralSec_A)
newdata2$scaled_SAND_RUB <- mean(m$scaled_SAND_RUB)
newdata2$scaled_EMA_MA <- mean(m$scaled_EMA_MA)
newdata2$scaled_HerbivoreBio <- mean(m$scaled_HerbivoreBio)
newdata2$scaled_WavePower <- mean(m$scaled_WavePower)
newdata2$scaled_Depth_Median<- mean(m$scaled_Depth_Median)
newdata2$scaled_logHumanDen <- mean(m$scaled_logHumanDen)
newdata2$scaled_MeanDHW10<-seq(min(m$scaled_MeanDHW10),max(m$scaled_MeanDHW10),
                              by=round(rg(m$scaled_MeanDHW10),3)/nrow(m))
newdata2$scaled_YearSinceDHW4<-mean(m$scaled_YearSinceDHW4)


o <- subset(final.df,YearSinceDHW4>10)
newdata3<-o
newdata3$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata3$scaled_CORAL <- mean(o$scaled_CORAL)
newdata3$scaled_CoralSec_A <- mean(o$scaled_CoralSec_A)
newdata3$scaled_SAND_RUB <- mean(o$scaled_SAND_RUB)
newdata3$scaled_EMA_MA <- mean(o$scaled_EMA_MA)
newdata3$scaled_HerbivoreBio <- mean(o$scaled_HerbivoreBio)
newdata3$scaled_WavePower <- mean(o$scaled_WavePower)
newdata3$scaled_Depth_Median<- mean(o$scaled_Depth_Median)
newdata3$scaled_logHumanDen <- mean(o$scaled_logHumanDen)
newdata3$scaled_MeanDHW10<-seq(min(o$scaled_MeanDHW10),max(o$scaled_MeanDHW10),
                               by=round(rg(o$scaled_MeanDHW10),5)/nrow(o))
newdata3$scaled_YearSinceDHW4<-mean(o$scaled_YearSinceDHW4)



p <- predict(best.mod, newdata = newdata1, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata1<-cbind(newdata1,p)
newdata1$Predict.lwr <- newdata1$Predicted_Juv - 1.96 * newdata1$SE_Juv # confidence interval upper bound
newdata1$Predict.upr <- newdata1$Predicted_Juv + 1.96 * newdata1$SE_Juv # confidence interval lower bound
newdata1$HSts_cat<-"0-3 years"

p <- predict(best.mod, newdata = newdata2, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata2<-cbind(newdata2,p)
newdata2$Predict.lwr <- newdata2$Predicted_Juv - 1.96 * newdata2$SE_Juv # confidence interval upper bound
newdata2$Predict.upr <- newdata2$Predicted_Juv + 1.96 * newdata2$SE_Juv # confidence interval lower bound
newdata2$HSts_cat<-"3-10 years"

p <- predict(best.mod, newdata = newdata3, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata3<-cbind(newdata3,p)
newdata3$Predict.lwr <- newdata3$Predicted_Juv - 1.96 * newdata3$SE_Juv # confidence interval upper bound
newdata3$Predict.upr <- newdata3$Predicted_Juv + 1.96 * newdata3$SE_Juv # confidence interval lower bound
newdata3$HSts_cat<-"> 10 years"

#Merge into 1 dataframe and add column for each HSts category. 
all.newdata<-rbind(newdata1,newdata2,newdata3)


#Unscaling predictor to plot on x axis
att <- attributes(scale(final.df$MeanDHW10))
mylabels <- seq(0,14,2)
mybreaks <- scale(mylabels, att$`scaled:center`, att$`scaled:scale`)[,1]

#colors<-c("cyan4","purple3","goldenrod2")
colors<-c("cyan4","purple3","gray67")

#Reorder HSts variables
all.newdata$HSts_cat <- factor(all.newdata$HSts_cat, levels = c("0-3 years","3-10 years","> 10 years"))
all.newdata<- all.newdata[order(all.newdata$HSts_cat),];head(all.newdata)



#Plot
plot1<-ggplot() +
  geom_line(data=all.newdata,aes(x = scaled_MeanDHW10, y = Predicted_Juv,color=HSts_cat),size=1) +
  geom_ribbon(data = all.newdata,aes(x = scaled_MeanDHW10,ymin = Predict.lwr, ymax = Predict.upr,fill=HSts_cat),alpha = 0.1)+
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title=element_blank(),
    legend.key.size = unit(1.5, 'cm'),
    legend.text = element_text(size=16),
    text = element_text(size = 18),
    panel.grid = element_blank()
  ) +
  ylab(expression(bold(paste("Predicted Juvenile Colonies",m^-2)))) +
  xlab(expression(bold('Mean Max '^o*'C-weeks')))  +
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  scale_x_continuous(labels = comma(mylabels),breaks=mybreaks)+
  scale_y_continuous(limits=c(0,20))+
  geom_rug(data=final.df,mapping=aes(x=scaled_MeanDHW10,y=0))


plot1

#### Wave Power x Heat stress- 
l <- subset(final.df,MeanDHW10<4)

newdata1<-l
newdata1$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata1$scaled_CORAL <- mean(l$scaled_CORAL)
newdata1$scaled_CoralSec_A <- mean(l$scaled_CoralSec_A)
newdata1$scaled_SAND_RUB <- mean(l$scaled_SAND_RUB)
newdata1$scaled_EMA_MA <- mean(l$scaled_EMA_MA)
newdata1$scaled_HerbivoreBio <- mean(l$scaled_HerbivoreBio)
newdata1$scaled_Depth_Median<- mean(l$scaled_Depth_Median)
newdata1$scaled_logHumanDen <- mean(l$scaled_logHumanDen)
newdata1$scaled_WavePower <- seq(min(l$scaled_WavePower),max(l$scaled_WavePower),
                            by=round(rg(l$scaled_WavePower),3)/nrow(l))
newdata1$scaled_MeanDHW10<-mean(l$scaled_MeanDHW10)
newdata1$scaled_YearSinceDHW4<-mean(l$scaled_YearSinceDHW4)


h <- subset(final.df,MeanDHW10>=4)

newdata2<-h
newdata2$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata2$scaled_CORAL <- mean(h$scaled_CORAL)
newdata2$scaled_CoralSec_A <- mean(h$scaled_CoralSec_A)
newdata2$scaled_SAND_RUB <- mean(h$scaled_SAND_RUB)
newdata2$scaled_EMA_MA <- mean(h$scaled_EMA_MA)
newdata2$scaled_HerbivoreBio <- mean(h$scaled_HerbivoreBio)
newdata2$scaled_Depth_Median<- mean(h$scaled_Depth_Median)
newdata2$scaled_logHumanDen <- mean(h$scaled_logHumanDen)
newdata2$scaled_WavePower <- seq(min(h$scaled_WavePower),max(h$scaled_WavePower),
                                by=round(rg(h$scaled_WavePower),6)/nrow(h))
newdata2$scaled_MeanDHW10<-mean(h$scaled_MeanDHW10)
newdata2$scaled_YearSinceDHW4<-mean(h$scaled_YearSinceDHW4)



p <- predict(best.mod, newdata = newdata1, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata1<-cbind(newdata1,p)
newdata1$Predict.lwr <- newdata1$Predicted_Juv - 1.96 * newdata1$SE_Juv # confidence interval upper bound
newdata1$Predict.upr <- newdata1$Predicted_Juv + 1.96 * newdata1$SE_Juv # confidence interval lower bound
newdata1$HSsev_cat<-"low"

p <- predict(best.mod, newdata = newdata2, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata2<-cbind(newdata2,p)
newdata2$Predict.lwr <- newdata2$Predicted_Juv - 1.96 * newdata2$SE_Juv # confidence interval upper bound
newdata2$Predict.upr <- newdata2$Predicted_Juv + 1.96 * newdata2$SE_Juv # confidence interval lower bound
newdata2$HSsev_cat<-"high"


#Merge into 1 dataframe and add column for each HStsev category. 
all.newdata<-rbind(newdata1,newdata2)


#Unscaling predictor to plot on x axis
att <- attributes(scale(final.df$WavePower))
mylabels <- seq(0,max(final.df$WavePower),75000)
mybreaks <- scale(mylabels, att$`scaled:center`, att$`scaled:scale`)[,1]

colors<-c("orange1","red3")

#Reorder HSts variables
all.newdata$HSsev_cat <- factor(all.newdata$HSsev_cat, levels = c("low","high"))
all.newdata<- all.newdata[order(all.newdata$HSsev_cat),];head(all.newdata)

#Plot
plot3<-ggplot() +
  geom_line(data=all.newdata,aes(x = scaled_WavePower, y = Predicted_Juv,color=HSsev_cat),size=1) +
  geom_ribbon(data = all.newdata,aes(x = scaled_WavePower,ymin = Predict.lwr, ymax = Predict.upr,fill=HSsev_cat),alpha = 0.1)+
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title=element_blank(),
    legend.key.size = unit(1.5, 'cm'),
    legend.text = element_text(size=16),
    text = element_text(size = 18),
    panel.grid = element_blank()
  ) +
  ylab(expression(bold(paste("Predicted Juvenile Colonies",m^-2)))) +
  xlab(expression(bold(Wave~Power~(kWh~m^{-1}))))  +
  scale_color_manual(labels=c(expression('< 4 '^o*'C-wk'),expression(""> '4 '^o*'C-wk')),values = colors)+
  scale_fill_manual(labels=c(expression('< 4 '^o*'C-wk'),expression(""> '4 '^o*'C-wk')),values = colors)+
  scale_x_continuous(labels = comma(mylabels),breaks=mybreaks)+
  scale_y_continuous(limits=c(0,20))+
  geom_rug(data=final.df,mapping=aes(x=scaled_WavePower,y=0))

plot3


#### Sector-level Cover x Heat stress- 
l <- subset(final.df,MeanDHW10<4)

newdata1<-l
newdata1$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata1$scaled_CORAL <- mean(l$scaled_CORAL)
newdata1$scaled_CoralSec_A <-seq(min(l$scaled_CoralSec_A),max(l$scaled_CoralSec_A),
                                 by=round(rg(l$scaled_CoralSec_A),4)/nrow(l))
newdata1$scaled_SAND_RUB <- mean(l$scaled_SAND_RUB)
newdata1$scaled_EMA_MA <- mean(l$scaled_EMA_MA)
newdata1$scaled_HerbivoreBio <- mean(l$scaled_HerbivoreBio)
newdata1$scaled_Depth_Median<- mean(l$scaled_Depth_Median)
newdata1$scaled_logHumanDen <- mean(l$scaled_logHumanDen)
newdata1$scaled_WavePower <- mean(l$scaled_WavePower)
newdata1$scaled_MeanDHW10<-mean(l$scaled_MeanDHW10)
newdata1$scaled_YearSinceDHW4<-mean(l$scaled_YearSinceDHW4)


h <- subset(final.df,MeanDHW10>=4)

newdata2<-h
newdata2$TRANSECTAREA_j <- 1 #Need to keep survey area constant
newdata2$scaled_CORAL <- mean(h$scaled_CORAL)
newdata2$scaled_CoralSec_A <- seq(min(h$scaled_CoralSec_A),max(h$scaled_CoralSec_A),
                                  by=round(rg(h$scaled_CoralSec_A),5)/nrow(h))
newdata2$scaled_SAND_RUB <- mean(h$scaled_SAND_RUB)
newdata2$scaled_EMA_MA <- mean(h$scaled_EMA_MA)
newdata2$scaled_HerbivoreBio <- mean(h$scaled_HerbivoreBio)
newdata2$scaled_Depth_Median<- mean(h$scaled_Depth_Median)
newdata2$scaled_logHumanDen <- mean(h$scaled_logHumanDen)
newdata2$scaled_WavePower <- mean(h$scaled_WavePower)
newdata2$scaled_MeanDHW10<-mean(h$scaled_MeanDHW10)
newdata2$scaled_YearSinceDHW4<-mean(h$scaled_YearSinceDHW4)



p <- predict(best.mod, newdata = newdata1, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata1<-cbind(newdata1,p)
newdata1$Predict.lwr <- newdata1$Predicted_Juv - 1.96 * newdata1$SE_Juv # confidence interval upper bound
newdata1$Predict.upr <- newdata1$Predicted_Juv + 1.96 * newdata1$SE_Juv # confidence interval lower bound
newdata1$HSsev_cat<-"low"

p <- predict(best.mod, newdata = newdata2, type = "response",se.fit=TRUE)
p<-as.data.frame(p)
colnames(p)<-c("Predicted_Juv","SE_Juv")
newdata2<-cbind(newdata2,p)
newdata2$Predict.lwr <- newdata2$Predicted_Juv - 1.96 * newdata2$SE_Juv # confidence interval upper bound
newdata2$Predict.upr <- newdata2$Predicted_Juv + 1.96 * newdata2$SE_Juv # confidence interval lower bound
newdata2$HSsev_cat<-"high"


#Merge into 1 dataframe and add column for each HStsev category. 
all.newdata<-rbind(newdata1,newdata2)


#Unscaling predictor to plot on x axis
att <- attributes(scale(final.df$scaled_CoralSec_A))
mylabels <- seq(0,max(final.df$CoralSec_A),50000000);mylabels
mybreaks <- scale(mylabels, att$`scaled:center`, att$`scaled:scale`)[,1]

colors<-c("orange1","red3")

#Reorder HSts variables
all.newdata$HSsev_cat <- factor(all.newdata$HSsev_cat, levels = c("low","high"))
all.newdata<- all.newdata[order(all.newdata$HSsev_cat),];head(all.newdata)

#Plot
plot2<-ggplot() +
  geom_line(data=all.newdata,aes(x = scaled_CoralSec_A, y = Predicted_Juv,color=HSsev_cat),size=1) +
  geom_ribbon(data = all.newdata,aes(x = scaled_CoralSec_A,ymin = Predict.lwr, ymax = Predict.upr,fill=HSsev_cat),alpha = 0.1)+
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title=element_blank(),
    legend.key.size = unit(1.5, 'cm'),
    legend.text = element_text(size=16),
    text = element_text(size = 18),
    panel.grid = element_blank()
  ) +
  ylab(expression(bold(paste("Predicted Juvenile Colonies",m^-2)))) +
  xlab(expression(bold("Sector % Coral Cover x Area (ha)")))  +
  scale_color_manual(labels=c(expression('< 4 '^o*'C-wk'),expression(""> '4 '^o*'C-wk')),values = colors)+
  scale_fill_manual(labels=c(expression('< 4 '^o*'C-wk'),expression(""> '4 '^o*'C-wk')),values = colors)+
  scale_x_continuous(labels = comma(mylabels),breaks=mybreaks)+
  scale_y_continuous(limits=c(-1,85))+
  geom_rug(data=final.df,mapping=aes(x=scaled_CoralSec_A,y=0))

plot2


# save full plot
setwd("T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/")
ytitle <- text_grob(expression(bold(paste("Predicted Juvenile Colonies ",m^-2))), size = 18,
                    face = "bold", rot = 90,hjust=0.2)
pdf(width = 16, height = 9.25, file = "Fig.5.pdf")
grid.arrange(arrangeGrob(plot1 + ggtitle("A)"),
                         plot2 + ggtitle("B)"),
                         plot3 + ggtitle("C)"),
                         nrow = 1), 
             nrow = 2, heights = c(10,1),
             left = ytitle)
dev.off()




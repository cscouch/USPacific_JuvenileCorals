# This script will clean the raw benthic REA data using method E that comes directly from the new data base application.
#Note- these data represent the revised data structure insituted in November 2018. Several recent dead and condition columns were added
#These data only include surveys conducted between 2013-2019

#Uses all forereef site-level data between 2013-2019 


###Add 2015 and 2017 Laysan cover data before finishing script
# 2017 sites "LAY-05112" "LAY-05120" "LAY-05131" "LAY-05142" "LAY-05157" "LAY-05177" "LAY-05178"

rm(list=ls())

#LOAD LIBRARY FUNCTIONS ... 
source("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/scripts/Functions_Juveniles.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/fish_team_functions.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/Islandwide Mean&Variance Functions.R")

library(forcats)
library(geosphere)
library(rgdal)
library(stringr)
library(mgcv)
detach(package:dplyr) #dplyr has issues if plyr is loaded first
library(dplyr)
library(tidyr)
library(RCurl)

#LOAD DATA
jwd_site<-read.csv("T:/Benthic/Projects/Juvenile Project/JuvProject_SITE_weights_AllYears_wHD.csv") #This dataframe has both juvenile and human density data

tsdhw<-read.csv("T:/Benthic/Projects/Juvenile Project/Juvenile_TimeSinceDHW4_8_v2.csv"); tsdhw<-subset(tsdhw,select= -c(X,ISLAND))
wave<-read.csv("T:/Benthic/Projects/Juvenile Project/Pacific_WaveActionData_v5.csv")

githubURL <- "https://github.com/cscouch/USPacific_JuvenileCorals/blob/main/Data/Survey_Master_Timeseries_2022-04-04.Rdata?raw=true"
load(url(githubURL))

cover1<-read.csv("T:/Benthic/Projects/Juvenile Project/Data/BenthicCover_2010-2019_Tier1_SITE.csv")#Cover from all sites
cover3<-read.csv("T:/Benthic/Projects/Juvenile Project/Data/BenthicCover_2010-2019_Tier3_SITE.csv")#Cover from all sites
sectors<-read.csv("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)


#Change Region Names to correspond to juvenile data
Convert_Region<-function(data){
  data$REGION<-ifelse(data$ISLAND %in% c("Farallon de Pajaros", "Maug", "Asuncion", "Alamagan", "Pagan", "Agrihan", "Guguan", "Sarigan","Farallon_de_Pajaros")
                      ,"NMI", as.character(data$REGION))
  data$REGION<-ifelse(data$ISLAND %in% c("Saipan", "Tinian", "Aguijan", "Rota", "Guam")
                    ,"SMI", as.character(data$REGION))
  data$REGION<-ifelse(data$ISLAND %in% c("Howland","Baker")
                    ,"PHOENIX", as.character(data$REGION))
  data$REGION<-ifelse(data$ISLAND =="Wake"
                    ,"WAKE", as.character(data$REGION))
  data$REGION<-ifelse(data$ISLAND %in% c("Kingman","Palmyra","Jarvis")
                    ,"LINE", as.character(data$REGION))
  return(data$REGION)
}

SM$REGION<-Convert_Region(SM)
cover1$REGION<-Convert_Region(cover1)
cover3$REGION<-Convert_Region(cover3)
sectors$REGION<-Convert_Region(sectors)
jwd_site$REGION<-Convert_Region(jwd_site)
SM$STRATANAME<-paste(SM$SEC_NAME,SM$REEF_ZONE,SM$DEPTH_BIN,sep="_")


#Remove spaces in island and sec names
jwd_site <- mutate_if(jwd_site, 
                is.character, 
                str_replace_all, pattern = " ", replacement = "_")

wave<-mutate_if(wave, 
                is.character, 
                str_replace_all, pattern = " ", replacement = "_")
cover1<-mutate_if(cover1, 
              is.character, 
              str_replace_all, pattern = " ", replacement = "_")

cover3<-mutate_if(cover3, 
                  is.character, 
                  str_replace_all, pattern = " ", replacement = "_")

sectors<-mutate_if(sectors, 
                  is.character, 
                  str_replace_all, pattern = " ", replacement = "_")

SM<-mutate_if(SM, 
                is.character, 
              str_replace_all, pattern = " ", replacement = "_")


#View(jwd_site)
nrow(jwd_site) #Should be 1405

#Columns to keep
jcols<-c("DATE_","SITEVISITID", "OBS_YEAR", "REGION", "ISLAND","SEC_NAME", "SITE","HABITAT_CODE","REEF_ZONE","MIN_DEPTH_M","MAX_DEPTH_M",
         "DEPTH_BIN","TRANSECTAREA_j","JuvColCount","JuvColDen","STRATANAME","LATITUDE", "LONGITUDE",
         "n","NH","sw","HumanDen")
jwd_site<-jwd_site[,jcols]
head(jwd_site)


#Calculate median depth
jwd_site<-jwd_site %>%
  rowwise() %>%
  mutate(Depth_Median=median(c(MIN_DEPTH_M,MAX_DEPTH_M)))


#Subset survey master and env columns of interest -not prior to 10/14/22, I was using the 750m VIIRS Chla data, but the data only goes back to 2012. Swaping in 4km data

cols<-c("MISSIONID","DATE_","SITEVISITID", "OBS_YEAR", "REGION", "ISLAND","SEC_NAME", "SITE","REEF_ZONE",
               "DEPTH_BIN","STRATANAME", "DHW.MeanMax_Degree_Heating_Weeks_YR01","DHW.MeanMax_Degree_Heating_Weeks_YR03",
        "DHW.MeanMax_Degree_Heating_Weeks_YR05","DHW.MeanMax_Degree_Heating_Weeks_YR10","DHW.MeanMax_Degree_Heating_Weeks_YR10YR01",
        "DHW.MaxMax_Degree_Heating_Weeks_YR03","DHW.MaxMax_Degree_Heating_Weeks_YR05","DHW.MaxMax_Degree_Heating_Weeks_YR10","DHW.Np10y_Major_Degree_Heating_Weeks_YR10",
        "mean_SST_CRW_Daily_YR10","sd_SST_CRW_Daily_YR10","mean_weekly_range_SST_CRW_Daily_YR10","mean_biweekly_range_SST_CRW_Daily_YR10","mean_Chlorophyll_A_ESAOCCCI_8Day_YR10","sd_Chlorophyll_A_ESAOCCCI_8Day_YR10",
        "mean_annual_range_Chlorophyll_A_ESAOCCCI_8Day_YR10")
sm_env<-SM[,cols]


#Not enough sampling in each sector- pool them together
sm_env$SEC_NAME<-ifelse(sm_env$SEC_NAME %in% c("TAU_OPEN","TAU_SANCTUARY"),"TAU",as.character(sm_env$SEC_NAME))
sm_env$SEC_NAME<-ifelse(sm_env$SEC_NAME %in% c("SWA_OPEN","SWA_SANCTUARY"),"SWA",as.character(sm_env$SEC_NAME))
View(sm_env)

#Calculate CVchla and CVsst
# sm_env$CVsst<-sm_env$sd_SST_CRW_Daily_YR10/sm_env$mean_SST_CRW_Daily_YR10
# sm_env$CVchla<-sm_env$sd_Chlorophyll_A_ESAOCCCI_8Day_YR05/sm_env$mean_Chlorophyll_A_ESAOCCCI_8Day_YR05

# #Change NaN to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

sm_env[is.nan(sm_env)] <- NA

#Combine Tier 1 and 3 cover
cover<-left_join(cover3,cover1[,c("SITEVISITID","CORAL","MA","TURF","SED")])  
nrow(cover);View(cover)

#Create summarized benthic columns
cover$RUBBLE<-cover$CCAR+cover$RUB+cover$TURFR
cover$TURF_BARE<-cover$TURFH+cover$HARD
cover$CCA<-cover$CCAH+cover$CCAR
cover$SAND<-cover$SED
cover$TURF<-cover$TURFH
cover$EMA_MA<-cover$EMA+cover$PESP+cover$MA

#Consolidate columns and combine rubble and sand
cols<-c("SITEVISITID", "OBS_YEAR", "REGION", "ISLAND","SEC_NAME", "SITE","REEF_ZONE",
        "DEPTH_BIN", "LATITUDE", "LONGITUDE","CORAL","CCA","RUBBLE","SAND","TURF","MA","EMA_MA")
cover<-cover[,cols]
cover$SAND_RUB<-cover$RUBBLE+cover$SAND
head(cover)

#Make tweaks to pooling sector pooling
#Not enough sampling in each sector- pool them together
cover$SEC_NAME<-ifelse(cover$SEC_NAME %in% c("TAU_OPEN","TAU_SANCTUARY"),"TAU",as.character(cover$SEC_NAME))
cover$SEC_NAME<-ifelse(cover$SEC_NAME %in% c("SWA_OPEN","SWA_SANCTUARY"),"SWA",as.character(cover$SEC_NAME))
View(cover)

#Convert Protected Reef Slope to Forereef and Subset just Forereef sites
cover$REEF_ZONE<-ifelse(cover$REEF_ZONE=="Protected_Slope","Forereef",as.character(cover$REEF_ZONE))
cover<-subset(cover,REEF_ZONE=="Forereef") #only include forereef


cover$STRATANAME<-paste(cover$SEC_NAME,cover$REEF_ZONE,cover$DEPTH_BIN,sep="_") #Create stratum

#Calculate strata mean cover for each stratum

cover_st<-cover %>%
  group_by(OBS_YEAR,STRATANAME) %>%
  summarize(CORALst=mean(CORAL,na.rm=T),
            CCAst=mean(CCA,na.rm=T),
            SAND_RUBst=mean(SAND_RUB,na.rm=T),
            TURFst=mean(TURF,na.rm=T),
            EMA_MAst=mean(EMA_MA,na.rm=T),)


#remove NAs from tsdhw
tsdhw<-tsdhw %>% filter(!is.na(YearSinceDHW4))

# Combine juvenile and predictor data at site-level -----------------------

cover_forsite<-cover %>% select(SITE,CORAL,CCA,SAND_RUB,TURF,EMA_MA)
sm_env<-sm_env %>% select(SITEVISITID, SITE,DHW.MeanMax_Degree_Heating_Weeks_YR01,DHW.MeanMax_Degree_Heating_Weeks_YR03,
                          DHW.MeanMax_Degree_Heating_Weeks_YR05,DHW.MeanMax_Degree_Heating_Weeks_YR10,DHW.MeanMax_Degree_Heating_Weeks_YR10YR01,
                          DHW.MaxMax_Degree_Heating_Weeks_YR03,DHW.MaxMax_Degree_Heating_Weeks_YR05,DHW.MaxMax_Degree_Heating_Weeks_YR10,DHW.Np10y_Major_Degree_Heating_Weeks_YR10,
                          mean_SST_CRW_Daily_YR10,sd_SST_CRW_Daily_YR10,mean_weekly_range_SST_CRW_Daily_YR10,mean_biweekly_range_SST_CRW_Daily_YR10,mean_Chlorophyll_A_ESAOCCCI_8Day_YR10,sd_Chlorophyll_A_ESAOCCCI_8Day_YR10,
                          mean_annual_range_Chlorophyll_A_ESAOCCCI_8Day_YR10)
wave<-wave %>% select(SITE,WavePower)

all_pred_site<- jwd_site   %>% 
  left_join(cover_forsite) %>%
  left_join(tsdhw) %>%
  left_join(sm_env) %>%
  left_join(wave)


nrow(jwd_site)
nrow(all_pred_site)

all_pred_site<-left_join(all_pred_site,cover_st)
nrow(all_pred_site)

#More manual tweaks
#there are some juvenile sites that don't have cover- sub with strata means for a given year- the fish team conducted surveys in the same strata
all_pred_site$CORAL<-ifelse(is.na(all_pred_site$CORAL),all_pred_site$CORALst,all_pred_site$CORAL)
all_pred_site$CCA<-ifelse(is.na(all_pred_site$CCA),all_pred_site$CCAst,all_pred_site$CCA)
all_pred_site$SAND_RUB<-ifelse(is.na(all_pred_site$SAND_RUB),all_pred_site$SAND_RUBst,all_pred_site$SAND_RUB)
all_pred_site$TURF<-ifelse(is.na(all_pred_site$TURF),all_pred_site$TURFst,all_pred_site$TURF)
all_pred_site$EMA_MA<-ifelse(is.na(all_pred_site$EMA_MA),all_pred_site$EMA_MAst,all_pred_site$EMA_MA)


#Missing depths- went back to datasheets or use max depth
# all_pred_site$Depth_Median<-ifelse(all_pred_site$SITE=="TIN-00646",21.3,all_pred_site$Depth_Median)
# all_pred_site$Depth_Median<-ifelse(all_pred_site$SITE=="TIN-00581",13.7,all_pred_site$Depth_Median)
# all_pred_site$Depth_Median<-ifelse(all_pred_site$SITE=="WAK-00430",9.44,all_pred_site$Depth_Median)
# all_pred_site$Depth_Median<-ifelse(all_pred_site$SITE=="TUT-01987",22.1,all_pred_site$Depth_Median)
all_pred_site$Depth_Median<-ifelse(is.na(all_pred_site$Depth_Median),all_pred_site$MAX_DEPTH_M,all_pred_site$Depth_Median)#missing min depth from several 2013 MHI sites


View(all_pred_site)

#Check for duplicate sites
all_pred_site %>% 
  group_by(SITE) %>% 
  filter(n()>1)

#Last Clean-up
all_pred_site[all_pred_site==-9991] <- NA #change -9991 in environemtnal data to NA

all_pred_site<-dplyr::filter(all_pred_site,SITE !="LAY-05016")

write.csv(all_pred_site, file="T:/Benthic/Projects/Juvenile Project/Data/JuvDen_Pred_SITE_AllYears.csv",row.names = F)



# Calculating sector-level BENTHIC COVER -----------------------------------------------------------


#We are missing some cover data from some benthic sites so use all benthic and fish sites
#Merge together wsd and sectors

sectors<-read.csv("C:/Users/courtney.s.couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SectorArea_Juveniles.csv", stringsAsFactors=FALSE)

strataKEEP<-unique(all_pred_site$STRATANAME)
cover<-dplyr::filter(cover,STRATANAME %in% strataKEEP)

#Add NH for forereef and protected reef sites
#Not enough sampling in each sector- pool them together
sectors$SEC_NAME<-ifelse(sectors$SEC_NAME %in% c("TAU_OPEN","TAU_SANCTUARY"),"TAU",as.character(sectors$SEC_NAME))
sectors$SEC_NAME<-ifelse(sectors$SEC_NAME %in% c("SWA_OPEN","SWA_SANCTUARY"),"SWA",as.character(sectors$SEC_NAME))

tmp<-ddply(sectors,.(SEC_NAME,REEF_ZONE,DEPTH_BIN,AREA_HA),
           summarize, temp=length(ISLAND))
new.NH<-ddply(tmp,.(SEC_NAME,REEF_ZONE,DEPTH_BIN),
              summarize, new.Area=sum(AREA_HA))
nrow(sectors)
sectors<-left_join(sectors,new.NH)
sectors<-subset(sectors,select=-c(AREA_HA))
sectors<-sectors %>% rename(AREA_HA=new.Area)
nrow(sectors)

sectors<-mutate_if(sectors, 
                  is.character, 
                  str_replace_all, pattern = " ", replacement = "_")



wsd<-left_join(cover,sectors[,c("SEC_NAME","REEF_ZONE","DEPTH_BIN","AREA_HA")]);nrow(wsd);head(wsd)

#Subset just Forereef sites
# wsd$REEF_ZONE<-ifelse(wsd$REEF_ZONE=="Protected Slope","Forereef",as.character(wsd$REEF_ZONE))
# wsd<-subset(wsd,REEF_ZONE=="Forereef")

wsd$SEC_NAME<-wsd$SEC_NAME
wsd$ANALYSIS_YEAR<-wsd$OBS_YEAR

data.cols<-c("CORAL","CCA","RUBBLE","SAND","TURF","MA","SAND_RUB")

### CALCULATE MEAN AND VARIANCE WITHIN STRATA ###
SPATIAL_POOLING_BASE<-c("REGION","ISLAND", "SEC_NAME", "REEF_ZONE", "STRATANAME","DEPTH_BIN")    
ADDITIONAL_POOLING_BY<-c("ANALYSIS_YEAR")                                    # additional fields that we want to break data at, but which do not relate to physical areas (eg survey year or method)

#generate within strata means and vars
POOLING_LEVEL<-c(SPATIAL_POOLING_BASE, ADDITIONAL_POOLING_BY)
dps<-Calc_PerStrata(wsd, data.cols, c(POOLING_LEVEL, "AREA_HA"))
head(dps$Mean)

###### REMOVE STRATA with N=1 (cannot pool those up)
dps$Mean<-dps$Mean[dps$Mean$N>1,]
dps$SampleVar<-dps$SampleVar[dps$SampleVar$N>1,]
dps$SampleSE<-dps$SampleSE[dps$SampleSE$N>1,]

# e.g. SAVE BY ISLAND AND REEF_ZONE PER YEAR
OUTPUT_LEVEL<-c("REGION","ISLAND","SEC_NAME","STRATANAME","DEPTH_BIN","ANALYSIS_YEAR") 
dpst<-Calc_Pooled_Simple(dps$Mean, dps$SampleVar, data.cols, OUTPUT_LEVEL, "AREA_HA");dpst<-as.data.frame(dpst)

#Clean up- remove SE columns and remove "Mean" from column names
dpst<-dpst %>% dplyr::select(Mean.REGION:Mean.SAND_RUB,-c(Mean.N))

dpst<-dpst %>%
  dplyr::rename_all(funs(stringr::str_replace_all(., "Mean.", "")))
head(dpst)

cover_sum<-dpst
head(cover_sum)

cover_sum3<-cover_sum

cover_sum2<-cover_sum%>% dplyr::filter(ANALYSIS_YEAR>2013)

#SAVE BY SECTOR PER YEAR
OUTPUT_LEVEL<-c("REGION","ISLAND","SEC_NAME","ANALYSIS_YEAR") 
dpsec<-Calc_Pooled_Simple(dps$Mean, dps$SampleVar, data.cols, OUTPUT_LEVEL, "AREA_HA")

dpsec<-as.data.frame(dpsec)
dpsec<-dpsec %>% dplyr::select(Mean.REGION:Mean.SAND_RUB,-c(Mean.N))

dpsec<-dpsec %>%
  dplyr::rename_all(funs(stringr::str_replace_all(., "Mean.", "")))
head(dpsec)

colnames(dpsec)[c(5:11)]<-paste(colnames(dpsec)[c(5:11)],"sec",sep="_")

write.csv(dpsec, file="T:/Benthic/Projects/Juvenile Project/BenthicCover_JuvenileProject_Tier1_SECTOR.csv",row.names = F)







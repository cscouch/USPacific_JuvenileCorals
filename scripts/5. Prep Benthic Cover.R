rm(list=ls())


library(gdata)             # needed for drop_levels()
library(reshape)           # reshape library inclues the cast() function used below
library(RODBC)            # to connect to oracle

#LOAD LIBRARY FUNCTIONS ... 
source("C:/Users/Courtney.S.Couch/Documents/GitHub/Benthic-Scripts/Functions/Benthic_Functions_newApp.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/core_functions.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/fish_team_functions.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/Islandwide Mean&Variance Functions.R")

#BIA data - this is from CPCE
load("T:/Benthic/Data/REA Coral Demography & Cover/Raw from Oracle/ALL_BIA_STR_RAW_NEW.rdata")   #bia

bia$SITE<-SiteNumLeadingZeros(bia$SITE)

#CNET data - from CoralNet
#These data contain human annotated data. There may be a small subset of robot annotated data. 
#The robot annotations are included because the confidence threshold in CoralNet was set to 70-90% allowing the robot to annotate points when it was 70-90% certain.
load("T:/Benthic/Data/REA Coral Demography & Cover/Raw from Oracle/ALL_BIA_STR_CNET.rdata") #load data

cnet$SITE<-SiteNumLeadingZeros(cnet$SITE)

#Temporary work around for merging in 2015 and 2017 NWHI data that hasn't been uploaded to Oracle yet- remove this once Michael has incorporated data
new.nw<-read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Raw Data from CoralNet/2015_2017_NWHI_CnetAnnotations_formatted.csv")
new.nw<-new.nw %>% drop_na(ROUNDID) #remove blank rows
new.laysan<-read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Raw Data from CoralNet/2015_2017_NWHI_CnetAnnotations_Laysan_formatted.csv")
new.laysan<-new.laysan %>% drop_na(ROUNDID) #remove blank rows

new.cnet<-rbind(new.nw,new.laysan)

class(new.cnet$DATE_)
class(new.cnet$DATE_TAKEN)

#Date conversations still not working
new.cnet$DATE_<-lubridate::mdy(new.cnet$DATE_)
new.cnet$DATE_TAKEN<-lubridate::ymd(new.cnet$DATE_TAKEN);head(new.cnet$DATE_TAKEN)
new.cnet$DATE_ANNOTATED<-lubridate::ymd_hms(new.cnet$DATE_ANNOTATED);head(new.cnet$DATE_ANNOTATED)

#combine old cnet and 2015 & 2017 nwhi cnet data
cnet<-rbind(cnet,new.cnet) 
table(cnet$REGION,cnet$OBS_YEAR)
table(new.cnet$ISLAND,new.cnet$OBS_YEAR)


##Generate Table of all the bia categories to review
head(bia)
bia_tab<-ddply(bia,.(TIER_1, CATEGORY_NAME, TIER_2, SUBCATEGORY_NAME, TIER_3, GENERA_NAME),summarize,count=sum(POINTS))
#write.csv(bia_tab, file="BIA categories.csv")
table(bia$TIER_1)
table(bia$TIER_2)

##Generate Table of all the bia categories to review
head(cnet)
cnet_tab<-ddply(cnet,.(CATEGORY_CODE,CATEGORY_NAME,SUBCATEGORY_CODE,SUBCATEGORY_NAME,GENERA_CODE,GENERA_NAME,FUNCTIONAL_GROUP),summarize,count=length(ROUNDID))

sm$SITE<-SiteNumLeadingZeros(sm$SITE)
sectors<-read.csv("C:/Users/courtney.s.couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SectorArea_Juveniles.csv")

test<-subset(sm,OBS_YEAR=="2019",TRANSECT_PHOTOS=="-1");nrow(test)

# Merge together all Photoquad Datasets & make sure columns match ---------------------------------------
bia$METHOD<-"CPCE"

cnet$POINTS<-1
cnet$METHOD<-"CNET"
cnet$REP<-cnet$REPLICATE
cnet$IMAGE_NAME<-cnet$ORIGINAL_FILE_NAME
cnet$PHOTOID<-cnet$IMAGE_NUMBER
cnet$TIER_1<-cnet$CATEGORY_CODE
cnet$TIER_2<-cnet$SUBCATEGORY_CODE
cnet$TIER_3<-cnet$GENERA_CODE


#Combine cpc and coralnet
FIELDS_TO_RETAIN<-c("MISSIONID","METHOD", "REGION", "OBS_YEAR","ISLAND", "SITEVISITID","SITE", "LATITUDE", "LONGITUDE", "REEF_ZONE", "DEPTH_BIN", "PERM_SITE", 
                    "CLIMATE_STATION_YN", "MIN_DEPTH", "MAX_DEPTH", "HABITAT_CODE", "REP", "IMAGE_NAME", "PHOTOID", "ANALYST", "TIER_1", "CATEGORY_NAME", 
                    "TIER_2", "SUBCATEGORY_NAME", "TIER_3", "GENERA_NAME", "POINTS")
x<-bia[,FIELDS_TO_RETAIN]; head(x)
y<-cnet[,FIELDS_TO_RETAIN]; head(y)

ab<-rbind(x,y)


#Flag sites that have more than 33 and less than 15 images
#With the exception of OCC 2012 sites, there should be 30 images/site/10 points/image
test<-ddply(ab,.(OBS_YEAR,SITEVISITID,SITE),summarize,count=sum(POINTS))
test2<-test[test$count<150 |test$count>330,]
View(test2)

#Remove sites with less than 150 points
test3<-test[test$count<150,];test3
ab<-ab[!(ab$SITE %in% test3$SITE),];head(ab)
subset(ab,SITE %in% c("TUT-00210","TUT-00275","OAH-00558")) #double check that sites were dropped properly

#Generate a table of # of sites/region and year from original datasets before data cleaning takes place
#use this later in the script to make sure sites haven't been dropped after data clean up.
oracle.site<-ddply(ab,.(REGION,OBS_YEAR),summarize,nSite=length(unique(SITE)))
oracle.site

#Check this against site master list
table(sm$REGION,sm$OBS)
ab.site<-ddply(subset(cnet,OBS_YEAR=="2019"),.(REGION,OBS_YEAR),summarize,nSite=length(unique(SITE)));ab.site

SURVEY_INFO<-c("OBS_YEAR", "REGION",  "ISLAND")
survey_island<-Aggregate_InputTable(cnet, SURVEY_INFO)

#There are some missing Tier3 information. If these data are missing then fill it with tier2 code
CATEGORY_FIELDS<-c("METHOD", "TIER_1", "CATEGORY_NAME", "TIER_2", "SUBCATEGORY_NAME", "TIER_3", "GENERA_NAME")
summary(ab[,CATEGORY_FIELDS])
levels(ab$TIER_3)<-c(levels(ab$TIER_3), levels(ab$TIER_2))
levels(ab$GENERA_NAME)<-c(levels(ab$GENERA_NAME), levels(ab$SUBCATEGORY_NAME))
ab[is.na(ab$TIER_3), ]$TIER_3<-ab[is.na(ab$TIER_3), ]$TIER_2
ab[is.na(ab$GENERA_NAME), ]$GENERA_NAME<-ab[is.na(ab$GENERA_NAME), ]$SUBCATEGORY_NAME
ab<-droplevels(ab)

# Reclassify EMA and Halimeda --------------------------------------------

#CREATING CLASS EMA "Encrusting Macroalgae
levels(ab$TIER_1)<-c(levels(ab$TIER_1), "EMA")
levels(ab$CATEGORY_NAME)<-c(levels(ab$CATEGORY_NAME), "Encrusting macroalga")
ab[ab$GENERA_NAME %in% c("Lobophora sp","Peyssonnelia sp", "Encrusting macroalga"),]$TIER_1<-"EMA"
ab[ab$GENERA_NAME %in% c("Lobophora sp","Peyssonnelia sp", "Encrusting macroalga"),]$TIER_2<-"EMA"
ab[ab$GENERA_NAME %in% c("Lobophora sp","Peyssonnelia sp", "Encrusting macroalga"),]$SUBCATEGORY_NAME<-"Encrusting macroalga"
ab[ab$GENERA_NAME %in% c("Lobophora sp","Peyssonnelia sp", "Encrusting macroalga"),]$CATEGORY_NAME<-"Encrusting macroalga"


###Create a Halimeda class
ab$GENERA_NAME<-as.character(ab$GENERA_NAME)
ab$TIER_1<-as.character(ab$TIER_1)

ab$TIER_3<-ifelse(ab$TIER_3=="HALI","HAL",as.character(ab$TIER_3))
ab$TIER_1<-ifelse(ab$TIER_3=="HAL","HAL",as.character(ab$TIER_1))
ab$CATEGORY_NAME<-ifelse(ab$TIER_3=="HAL","Halimeda sp",ab$CATEGORY_NAME)

hal<-subset(ab,TIER_1=="HAL")
head(hal)


test<-ddply(ab,.(REGION,OBS_YEAR),summarize,nSite=length(unique(SITE)))
test

#### WORKING WITH CLEAN DATA FILE AT THIS POINT  
ab<-droplevels(ab)
table(ab$ISLAND, ab$OBS_YEAR)

summary(ab)

#We are missing depth bin, reef zone and habitat_code information from some sites.
#This information is also missing from the SURVEY MASTER file

levels(ab$DEPTH_BIN)<-c(levels(ab$DEPTH_BIN), "UNKNOWN")
levels(ab$REEF_ZONE)<-c(levels(ab$REEF_ZONE), "UNKNOWN")
levels(ab$HABITAT_CODE)<-c(levels(ab$HABITAT_CODE), "UNKNOWN")

ab[is.na(ab$DEPTH_BIN),]$DEPTH_BIN<-"UNKNOWN"
ab[is.na(ab$REEF_ZONE),]$REEF_ZONE<-"UNKNOWN"
ab[is.na(ab$HABITAT_CODE),]$HABITAT_CODE<-"UNKNOWN"


#Generate a SITE table
sites<-unique(ab[,c("METHOD","REGION","OBS_YEAR","ISLAND","PERM_SITE","CLIMATE_STATION_YN","SITEVISITID","LATITUDE","LONGITUDE","REEF_ZONE","DEPTH_BIN")])

dim(sites)


# Generate Site-level Data at TIER 1 level--------------

photo<-dcast(ab, formula=METHOD + OBS_YEAR + SITEVISITID + SITE  ~ TIER_1, value.var="POINTS", sum, fill=0)
head(photo)

r_levels<-c(unique(as.character(ab$TIER_1)))
photo$N<-rowSums(photo[,r_levels])

#Subtract mobile inverts and tape wand shallow and unclassified
photo$new.N<-photo$N-(photo$MF+photo$UC+photo$TW)

r_levels<-c(unique(as.character(ab$TIER_1)))
data.cols<-c(r_levels)
data.cols

#Calculate proportion
photo[,data.cols]<-(photo[,data.cols]/photo$new.N)*100
head(photo)

r_levels<-c(unique(as.character(ab$TIER_1)))
T1data.cols<-c(r_levels)
T1data.cols<-T1data.cols[!T1data.cols %in% c("TW","UC","MF")]

wsd<-merge(sites, photo, by=c("METHOD", "OBS_YEAR", "SITEVISITID"), all.y=T)

#Make sure that you have the correct # of sites/region and year
test1<-ddply(wsd,.(REGION,OBS_YEAR),summarize,nSite_wsd=length(unique(SITE)))

#check against original number of sites pulled from oracle
full_join(test1,oracle.site)

#Merge Tier 1 data with SURVEY MASTER FILE
sm<-read.csv("C:/Users/courtney.s.couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SURVEY MASTER_Juveniles.csv")

#Convert date formats
sm$DATE_<-lubridate::mdy(sm$DATE_)
class(sm$DATE_)

head(sm)

sm<-sm[,c("DATE_","MISSIONID","SITEVISITID","SITE","ANALYSIS_YEAR","ANALYSIS_SCHEME","OBS_YEAR","SEC_NAME","EXCLUDE_FLAG","new_MIN_DEPTH_M","new_MAX_DEPTH_M")]
wsd_t1<-merge(sm,wsd,by=c("SITEVISITID","SITE","OBS_YEAR"),all.y=TRUE)
head(wsd_t1)

#remove 
wsd_t1

#Remove the unknowns and TWS columns
wsd_t1<-subset(wsd_t1,select= -c(MF,UC,TW))

#Remove sites that have less than 150 points after removing MF,UC, TW (a lot of our early imagery was poor quality and will be removed)
wsd_t1<-subset(wsd_t1,new.N >=150)

#Remove special mission data with exclude flag =-1
wsd_t1<-subset(wsd_t1,EXCLUDE_FLAG!="-1")

#Save Tier 1 site data to t drive. This file has all sites (fish, benthic and OCC) that were annoated between 2010 and 2018
write.csv(wsd_t1, file="T:/Benthic/Projects/Juvenile Project/Data/BenthicCover_2010-2019_Tier1_SITE.csv",row.names=F)


# Generate Site-level Data at TIER 3 level--------------
photo<-dcast(ab, formula=METHOD + OBS_YEAR + SITEVISITID + SITE  ~ TIER_3, value.var="POINTS", sum, fill=0)
head(photo)

r_levels<-c(unique(as.character(ab$TIER_3)))
photo$N<-rowSums(photo[,r_levels])
data.cols<-c(r_levels)

#Substract mobile inverts and tape wand shallow and uclassified
photo$new.N<-photo$N-(photo$WAND+photo$UNK+photo$TAPE+photo$MOBF+photo$SHAD)


#Calculate proportion
photo[,data.cols]<-(photo[,data.cols]/photo$new.N)*100
head(photo)

r_levels<-c(unique(as.character(ab$TIER_3)))
T3data.cols<-c(r_levels)
T3data.cols<-T3data.cols[!T3data.cols %in% c("WAND","UNK","TAPE","MOBF","SHAD")]

wsd<-merge(sites, photo, by=c("METHOD", "OBS_YEAR", "SITEVISITID"), all.y=T)

#Merge Tier 3 data with SURVEY MASTER data
wsd_t3<-merge(sm,wsd,by=c("SITEVISITID","SITE","OBS_YEAR"),all.y=TRUE)
head(wsd_t3)

test<-wsd_t3[is.na(wsd_t3$TRANSECT_PHOTOS),]
View(test) # none of the 2010 imagery has TRANSECT_PHOTOS assigned - ASK MICHAEL TO FIX
wsd_t3$TRANSECT_PHOTOS<-"-1" #make sure that all rows = -1


#Remove the unknowns and TWS columns
wsd_t3<-subset(wsd_t3,select= -c(WAND,UNK,TAPE,MOBF,SHAD))

#Remove sites that have less than 150 points after removing MF,UC, TW (a lot of our early imagery was poor quality and will be removed)
wsd_t3<-subset(wsd_t3,new.N >=150)

#Remove special mission data with exclude flag =-1
wsd_t3<-subset(wsd_t3,EXCLUDE_FLAG!="-1")


#Save Tier 1 site data to t drive. This file has all sites (fish, benthic and OCC) that were annoated between 2010 and 2018
write.csv(wsd_t3, file="T:/Benthic/Projects/Juvenile Project/Data/BenthicCover_2010-2019_Tier3_SITE.csv",row.names = F)


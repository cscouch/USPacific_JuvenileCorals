#This script imports raw point data from CPCe and CoralNet and calculates site-level cover data at the functional and tier 3 level.


rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/USPacific_JuvenileCorals/"))


#Load Functions and helper files 
source("scripts/Functions_Juveniles.R")

sm<-read.csv("SupportFiles/SURVEY MASTER_Juveniles.csv") #The survey master file is a full list of all of the sites with metadata
sectors<-read.csv("SupportFiles/SectorArea_Juveniles.csv") 

sm$SITE<-as.factor(sm$SITE)
sm$SITE<-SiteNumLeadingZeros(sm$SITE)

#Convert date formats
sm$DATE_<-lubridate::mdy(sm$DATE_)
class(sm$DATE_)

head(sm)

# Load and format Point Data ----------------------------------------------

#CPCe data (2010-2014)
load("Data/PacificNCRMP_CPCE_StRS.rdata")   #bia

# Change site number such as MAR-22 to MAR-0022 to avoid changing to "March 22"
bia$SITE<- as.factor(bia$SITE)
bia$SITE<-SiteNumLeadingZeros(bia$SITE)

#CoralNet data (2014-2019)
#These data contain human annotated data. There may be a small subset of robot annotated data. 
#The robot annotations are included because the confidence threshold in CoralNet was set to 70-90% allowing the robot to annotate points when it was 70-90% certain.
load("Data/PacificNCRMP_CNET_StRS.rdata") #load data

cnet$SITE<- as.factor(cnet$SITE)
cnet$SITE<-SiteNumLeadingZeros(cnet$SITE)

#CoralNet data from NWHI 2015 and 2017 that hasn't been uploaded NCEI yet
new.cnet<-read.csv("Data/2015_2017_NWHI_CnetAnnotations.csv")
new.cnet<-new.cnet %>% drop_na(ROUNDID) #remove blank rows

class(new.cnet$DATE_)
class(new.cnet$DATE_ANNOTATED)

#Convert Date formats
new.cnet$DATE_<-lubridate::mdy(new.cnet$DATE_)
new.cnet$DATE_TAKEN<-lubridate::mdy(new.cnet$DATE_TAKEN);head(new.cnet$DATE_TAKEN)
new.cnet$DATE_ANNOTATED<-as.Date(mdy_hms(new.cnet$DATE_ANNOTATED));head(new.cnet$DATE_ANNOTATED)

#Combine all CoralNet dataframes
cnet<-rbind(cnet,new.cnet) 
table(cnet$REGION,cnet$OBS_YEAR)


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


#Combine cpce and coralnet
FIELDS_TO_RETAIN<-c("MISSIONID","METHOD", "REGION", "OBS_YEAR","ISLAND", "SITEVISITID","SITE", "LATITUDE", "LONGITUDE", "REEF_ZONE", "DEPTH_BIN", "PERM_SITE", 
                    "CLIMATE_STATION_YN", "MIN_DEPTH", "MAX_DEPTH", "HABITAT_CODE", "REP", "IMAGE_NAME", "PHOTOID", "ANALYST", "TIER_1", "CATEGORY_NAME", 
                    "TIER_2", "SUBCATEGORY_NAME", "TIER_3", "GENERA_NAME", "POINTS")
x<-bia[,FIELDS_TO_RETAIN]; head(x)
y<-cnet[,FIELDS_TO_RETAIN]; head(y)

ab<-rbind(x,y) #Full raw point data


#Generate a table of # of sites/region and year from original datasets before data cleaning takes place
#use this later in the script to make sure sites haven't been dropped after data clean up.
oracle.site<-ddply(ab,.(REGION,OBS_YEAR),summarize,nSite=length(unique(SITE)))
oracle.site


# Reclassify Encrusting Macroalgae and Halimeda --------------------------------------------

#There are some missing Tier3 information. If these data are missing then fill it with tier2 code
CATEGORY_FIELDS<-c("METHOD", "TIER_1", "CATEGORY_NAME", "TIER_2", "SUBCATEGORY_NAME", "TIER_3", "GENERA_NAME")
summary(ab[,CATEGORY_FIELDS])
levels(ab$TIER_3)<-c(levels(ab$TIER_3), levels(ab$TIER_2))
levels(ab$GENERA_NAME)<-c(levels(ab$GENERA_NAME), levels(ab$SUBCATEGORY_NAME))
ab[is.na(ab$TIER_3), ]$TIER_3<-ab[is.na(ab$TIER_3), ]$TIER_2
ab[is.na(ab$GENERA_NAME), ]$GENERA_NAME<-ab[is.na(ab$GENERA_NAME), ]$SUBCATEGORY_NAME
ab<-droplevels(ab)


#Create EMA class "Encrusting Macroalgae
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


test<-ddply(ab,.(REGION,OBS_YEAR),summarize,nSite=length(unique(SITE)))
test

#### Clean Data 
ab<-droplevels(ab)
table(ab$ISLAND, ab$OBS_YEAR)

summary(ab)
# 
# #We are missing depth bin, reef zone and habitat_code information from some sites.
# #This information is also missing from the SURVEY MASTER file
# 
# levels(ab$DEPTH_BIN)<-c(levels(ab$DEPTH_BIN), "UNKNOWN")
# levels(ab$REEF_ZONE)<-c(levels(ab$REEF_ZONE), "UNKNOWN")
# levels(ab$HABITAT_CODE)<-c(levels(ab$HABITAT_CODE), "UNKNOWN")
# 
# ab[is.na(ab$DEPTH_BIN),]$DEPTH_BIN<-"UNKNOWN"
# ab[is.na(ab$REEF_ZONE),]$REEF_ZONE<-"UNKNOWN"
# ab[is.na(ab$HABITAT_CODE),]$HABITAT_CODE<-"UNKNOWN"


#Generate a SITE table
sites<-unique(ab[,c("METHOD","REGION","OBS_YEAR","ISLAND","PERM_SITE","CLIMATE_STATION_YN","SITEVISITID","LATITUDE","LONGITUDE","REEF_ZONE","DEPTH_BIN")])

dim(sites)


# Generate Site-level Data at TIER 1 (functional) level--------------

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
write.csv(wsd_t1, file="Data/outputs/BenthicCover_2010-2019_Tier1_SITE.csv",row.names=F)


# Generate Site-level Data at TIER 3 (genus-morphology for corals, dominant genus level for algae) level--------------
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

wsd_t3$TRANSECT_PHOTOS<-"-1" #make sure that all rows = -1


#Remove the unknowns and TWS columns
wsd_t3<-subset(wsd_t3,select= -c(WAND,UNK,TAPE,MOBF,SHAD))

#Remove sites that have less than 150 points after removing MF,UC, TW (a lot of our early imagery was poor quality and will be removed)
wsd_t3<-subset(wsd_t3,new.N >=150)

#Remove special mission data with exclude flag =-1
wsd_t3<-subset(wsd_t3,EXCLUDE_FLAG!="-1")


#Save Tier 1 site data to t drive. This file has all sites (fish, benthic and OCC) that were annoated between 2010 and 2018
write.csv(wsd_t3, file="Data/outputs/BenthicCover_2010-2019_Tier3_SITE.csv",row.names = F)


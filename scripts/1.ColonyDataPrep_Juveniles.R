# Using R version 4.1.0 (2021-05-18)

# This script will clean the raw colony-level NOAA NCRMP data
rm(list=ls())

#LOAD LIBRARY FUNCTION ...
source("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/scripts/Functions_Juveniles.R")

## LOAD raw data downloaded from NCEI:
# https://www.fisheries.noaa.gov/inport/item/36165
# https://www.fisheries.noaa.gov/inport/item/22790
# https://www.fisheries.noaa.gov/inport/item/36164
# https://www.fisheries.noaa.gov/inport/item/36166

#The 4 files were combined then read in
df<-read.csv("T:/Benthic/Projects/Juvenile Project/Data/ALL_REA_JUVCORAL_RAW_2013-2020_JuvenileProject.csv")

x<-df

#I was having issues with NA values in the 2017 NWHI data- had to remerge data. Delete this commented text before going public
#load("T:/Benthic/Data/REA Coral Demography & Cover/Raw from Oracle/ALL_REA_JUVCORAL_RAW_2013-2020.rdata") 
# x<-df #leave this as df
# 
# x$RY<-paste(x$REGION,x$OBS_YEAR,sep="_")
# x<-subset(x,RY!="NWHI_2017")
# 
# #Convert date formats
# class(x$DATE_)
# x$DATE_ <- as.Date(x$DATE_, format = "%Y-%m-%d")
# 
# #Create vector of column names to include then exclude unwanted columns from dataframe
# DATA_COLS<-c("MISSIONID","REGION","REGION_NAME","ISLAND","ISLANDCODE","SITE","REEF_ZONE","DEPTH_BIN","OBS_YEAR",
#              "DATE_","NO_SURVEY_YN","EXCLUDE_FLAG","SITEVISITID","HABITAT_CODE","DIVER","TRANSECTNUM","SEGMENT","SEGWIDTH","SEGLENGTH",
#              "COLONYID","TAXONCODE","MORPH_CODE","COLONYLENGTH","GENUS_CODE","S_ORDER","TAXONNAME")
# 
# 
# #remove extraneous columns
# head(x[,DATA_COLS])
# x<-x[,DATA_COLS]
# 
# #Cleanup 2017 NWHI data to merge with the rest of the juvenile data -temporary workaround until data team migrates data to Oracle
# nw<-read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Raw from Oracle/PMNM2017_JUVENILECOLONY_QCd.csv")
# 
# head(nw[,DATA_COLS])
# nw<-nw[,DATA_COLS]
# 
# nw$COLONYID<-nrow(x)+1:length(nw$SITEVISITID)
# nw$COLONYID<-ifelse(nw$TAXONCODE=="AAAA",NA,nw$COLONYID)
# 
# #Create vector of column names to include then exclude unwanted columns from dataframe
# DATA_COLS<-c("MISSIONID","REGION","REGION_NAME","ISLAND","ISLANDCODE","SITE","REEF_ZONE","DEPTH_BIN","OBS_YEAR",
#              "DATE_","NO_SURVEY_YN","EXCLUDE_FLAG","SITEVISITID","HABITAT_CODE","DIVER","TRANSECTNUM","SEGMENT","SEGWIDTH","SEGLENGTH",
#              "COLONYID","TAXONCODE","MORPH_CODE","COLONYLENGTH","GENUS_CODE","S_ORDER","TAXONNAME")
# 
# 
# #remove extraneous columns
# head(nw[,DATA_COLS])
# nw<-nw[,DATA_COLS]
# 
# #Convert date formats
# class(nw$DATE_)
# nw$DATE_ <- dmy(nw$DATE_)
# nw$DATE_ <- as.Date(nw$DATE_, format = "%Y-%m-%d")
# 
# x<-rbind(x,nw)
# 
# write.csv(x,file="T:/Benthic/Data/REA Coral Demography & Cover/Raw from Oracle/ALL_REA_JUVCORAL_RAW_2013-2020_JuvenileProject.csv")

# Change site number such as MAR-22 to MAR-0022 to avoid changing to "March 22"
x$SITE<- as.factor(x$SITE)
x$SITE<-SiteNumLeadingZeros(x$SITE) 

### Use these functions to look at data
head(x)
tail(x)
table(x$REGION, x$OBS_YEAR) #review years and regions in dataframe


#Double check level and class of variables to make sure there aren't any errors
sapply(x,levels)
sapply(x,class)##Change column names to make code easier to code

colnames(x)[colnames(x)=="TRANSECTNUM"]<-"TRANSECT" #Change column name

# Merge Juvenile data and SITE MASTER -------------------------------------
#The survey master file is a full list of all of the sites with 
# load site master to merge with colony data
survey_master<-read.csv("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SURVEY MASTER_Juveniles.csv")

#Use SM coordinates-some coordinates are wrong in data and need to be updated
colnames(survey_master)[colnames(survey_master)=="LATITUDE_LOV"]<-"LATITUDE" #Change column name- we will eventually change this column back to "taxoncode" after we modify the spcode names to match the taxalist we all feel comfortable identifying
colnames(survey_master)[colnames(survey_master)=="LONGITUDE_LOV"]<-"LONGITUDE" #Change column name- we will eventually change this column back to "taxoncode" after we modify the spcode names to match the taxalist we all feel comfortable identifying

#Check that OBS_YEAR, SITEVISITID, and SITE are all the same in both x and survey master
OYerror<-which(x$OBS_YEAR!=survey_master$OBS_YEAR[match(x$SITEVISITID,survey_master$SITEVISITID)])
SIerror<-which(as.vector(x$SITE)!=survey_master$SITE[match(x$SITEVISITID,survey_master$SITEVISITID)])
SIOYerrors<-unique(c(OYerror,SIerror))
if(length(SIOYerrors)>0){print(paste0("Warning: Raw Data disagree with Survey Master for sitevisitids: ",x$SITEVISITID[SIOYerrors]))}

#merge colony data and survey master
x<-left_join(x, survey_master[,c("OBS_YEAR","SITEVISITID","SITE","LATITUDE","LONGITUDE","SEC_NAME","ANALYSIS_YEAR","new_MIN_DEPTH_M","new_MAX_DEPTH_M")])

colnames(x)[colnames(x)=="new_MIN_DEPTH_M"]<-"MIN_DEPTH_M" #Change column name
colnames(x)[colnames(x)=="new_MAX_DEPTH_M"]<-"MAX_DEPTH_M" #Change column name

#CHECK THAT all SEC_NAME are present in the survey_master file
test<-x[is.na(x$SEC_NAME), c("MISSIONID","REGION", "SITE","OBS_YEAR"),]
test<-droplevels(test);table(test$SITE,test$MISSIONID) #create a table of missing sites by missionid
if(dim(test)[1]>0) {cat("sites with MISSING SECTORS present")}   # should be 0
#These are sites from a special mission that will be deleted

# CLEAN UP ----------------------------------------------------------------

##Remove sites that were only surveyed for photoquads but not demographics
#Note-photoquad only sites are not included in data prior to 2018
#Test whether there are missing values in the NO_SURVEY_YN column. The value should be 0 or -1
x.na<-x[is.na(x$NO_SURVEY_YN)&x$OBS_YEAR>2013,]
View(x.na)

x$NO_SURVEY_YN<-is.na(x$NO_SURVEY_YN)<-0 #Change NAs (blank cells) to 0
x<-subset(x,NO_SURVEY_YN==0)
x<-subset(x,SEGLENGTH!="NA") #Remove segments that were not surveyed for coral demography

##Clean-up taxa information
#AAAA denotes segments where no colonies were observed
#Check to see whether S_ORDER is NA and not AAAA (the code for no colonies observed on the segment)
x[x$TAXONCODE!="AAAA"& is.na(x$S_ORDER),] #this dataframe should be empty

#Change columns to character
x$GENUS_CODE<-as.character(x$GENUS_CODE)
x$TAXONCODE<-as.character(x$TAXONCODE)
x$S_ORDER<-as.character(x$S_ORDER)

#Make sure there are no NA values in genus code or taxoncode if it's supposed to be a scleractinian
subset(x,S_ORDER=="Scleractinia" & GENUS_CODE=="NA") #this dataframe should be empty
subset(x,S_ORDER=="Scleractinia" & TAXONCODE=="NA") #this dataframe should be empty

#There are some old SPCODES that were a combination of taxa and weren't included in the complete taxa list
#Change these unknown genera or taxoncodes to the spcode and the remaining NAs in the Taxon and genus code to AAAA
x$GENUS_CODE<-ifelse(x$TAXONCODE=="UNKN","UNKN",x$GENUS_CODE)
x$GENUS_CODE<-ifelse(x$TAXONCODE=="AAAA","AAAA",x$GENUS_CODE)
x$GENUS_CODE<-ifelse(x$TAXONCODE %in% c("MOAS","LEPA"),"UNKN",x$GENUS_CODE)

View(x) #view data in separate window

#Check that Unknown scl were changed correctly
head(subset(x,TAXONCODE=="UNKN"&S_ORDER=="Scleractinia"),40)
head(subset(x,GENUS_CODE=="UNKN"&S_ORDER=="Scleractinia"))
head(subset(x,GENUS_CODE=="AAAA"))


##Calcuating segment and transect area and add column for transect area
x$TRANSECTAREA<-Transectarea(x)
summary(x$TRANSECTAREA)
head(x)
nrow(x)


## CLEAN UP NAs ##
NegNineCheckCols=c("S_ORDER","TAXONNAME","COLONYLENGTH")
x[,NegNineCheckCols][x[,NegNineCheckCols] ==-9] <- NA #Convert missing numeric values to NA (they are entered as -9 in Oracle)


jwd<-droplevels(x)
write.csv(jwd,file="T:/Benthic/Projects/Juvenile Project/Data/CoralBelt_Juveniles_raw_CLEANED.csv",row.names = FALSE)


# Using R version 4.1.0 (2021-05-18)

# This script will clean the raw colony-level NOAA NCRMP data
rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/USPacific_JuvenileCorals/"))
      
#LOAD LIBRARY FUNCTION ...
source("scripts/Functions_Juveniles.R")

#Load Raw Juvenile Data
## These data can also be downloaded from the National Center for Environmental Information (see links in manuscript)
load("Data/ALL_REA_JUVCORAL_RAW_2013-2022.rdata")

x<-subset(df,OBS_YEAR<2020)


# Data clean-up -----------------------------------------------------------

# Change site number such as MAR-22 to MAR-0022 to avoid changing to "March 22"
x$SITE<- as.factor(x$SITE)
x$SITE<-SiteNumLeadingZeros(x$SITE) 


#Convert date formats
class(x$DATE_)
x$DATE_ <- as.Date(x$DATE_, format = "%Y-%m-%d")

### Use these functions to look at data
head(x)
tail(x)
table(x$REGION, x$OBS_YEAR) #review years and regions in dataframe


#Double check level and class of variables to make sure there aren't any errors
sapply(x,levels)
sapply(x,class)##Change column names to make code easier to code

colnames(x)[colnames(x)=="TRANSECTNUM"]<-"TRANSECT" #Change column name

#remove Depth and coordinates columns from juvenile dataset- some of them values are missing. Use depths and coords from the SURVEY MASTER file
x<-subset(x,select=-c(SITE_MAX_DEPTH,SITE_MIN_DEPTH,LATITUDE,LONGITUDE))

# Merge Juvenile data and SITE MASTER -------------------------------------
#The survey master file is a full list of all of the sites with metadata
survey_master<-read.csv("SupportFiles/SURVEY MASTER_Juveniles.csv") 

#Use SM coordinates
colnames(survey_master)[colnames(survey_master)=="LATITUDE_LOV"]<-"LATITUDE" #Change column name
colnames(survey_master)[colnames(survey_master)=="LONGITUDE_LOV"]<-"LONGITUDE" #Change column name

#Check that OBS_YEAR, SITEVISITID, and SITE are all the same in both x and survey master
OYerror<-which(x$OBS_YEAR!=survey_master$OBS_YEAR[match(x$SITEVISITID,survey_master$SITEVISITID)])
SIerror<-which(as.vector(x$SITE)!=survey_master$SITE[match(x$SITEVISITID,survey_master$SITEVISITID)])
SIOYerrors<-unique(c(OYerror,SIerror))
if(length(SIOYerrors)>0){print(paste0("Warning: Raw Data disagree with Survey Master for sitevisitids: ",x$SITEVISITID[SIOYerrors]))}

#merge colony data and survey master
x<-left_join(x, survey_master[,c("OBS_YEAR","SITEVISITID","SITE","LATITUDE","LONGITUDE","SEC_NAME","new_MIN_DEPTH_M","new_MAX_DEPTH_M")])

colnames(x)[colnames(x)=="new_MIN_DEPTH_M"]<-"MIN_DEPTH_M" #Change column name
colnames(x)[colnames(x)=="new_MAX_DEPTH_M"]<-"MAX_DEPTH_M" #Change column name

#CHECK THAT all SEC_NAME are present in the survey_master file
test<-x[is.na(x$SEC_NAME), c("MISSIONID","REGION", "SITE","OBS_YEAR"),]
test<-droplevels(test);table(test$SITE,test$MISSIONID) #create a table of missing sites by missionid
if(dim(test)[1]>0) {cat("sites with MISSING SECTORS present")}   # should be 0
#These are sites from a special mission that will be deleted


# Final Clean-up ----------------------------------------------------------------

##Remove sites that were only surveyed for benthic cover but juveniles
#Test whether there are missing values in the NO_SURVEY_YN column. The value should be 0 or -1
x.na<-x[is.na(x$NO_SURVEY_YN)&x$OBS_YEAR>2013,]
head(x.na)

x$NO_SURVEY_YN<-is.na(x$NO_SURVEY_YN)<-0 #Change NAs (blank cells) to 0
x<-subset(x,NO_SURVEY_YN==0)
x<-subset(x,SEGLENGTH!="NA") #Remove segments that were not surveyed for juveniles

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



## CLEAN UP NAs ##
NegNineCheckCols=c("S_ORDER","TAXONNAME","COLONYLENGTH")
x[,NegNineCheckCols][x[,NegNineCheckCols] ==-9] <- NA #Convert missing numeric values to NA (they are entered as -9 in Oracle)


jwd<-droplevels(x)
write.csv(jwd,file="Data/outputs/CoralBelt_Juveniles_raw_CLEANED.csv",row.names = FALSE)


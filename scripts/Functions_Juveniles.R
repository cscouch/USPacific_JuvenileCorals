# Using R version 4.1.0 (2021-05-18)



library(reshape2)
library(ggplot2) 
library(data.table)
library(plyr)
library(dplyr)
library(gdata)
library(tidyr)
library(plotrix)
library(scales)  # for pretty_breaks() function
library(splitstackshape)
library(lsmeans)
library(multcomp)
library(lubridate)
library(forcats)
library(stringr)
library(rcompanion)
library(car)
library(gridExtra)
library(sp)
library(sf)
library(ggsn)
library(ggspatial)
library(ggrepel)
library(rnaturalearth)
library(rgeos)
library(cowplot)
library(survey)
library(svydiags)
library(multcompView)



######################################################
# function SiteNumLeadingZeros, adds leasing zeros to SITE code numeric parts to make eg OAH-1 become OAH-001
# this function therefore makes it easier to sort site names meaningfully and also removes the problem of eg MAR-22 site being treated in csv output as if it means March 22nd
# some site names have letter in the second 3 portion eg GAR-R3 .. those sites are not changed, because there re very few of those and those are generally well sorted anyway
# (and, it seems harder to work out which situations those are and how to deal with all possible variants .. therefore code just runs for situations where there are only digits in the part of the site name after the hyphen) 
#####################################################
SiteNumLeadingZeros <- function(site_names)
{	
  tmp<-levels(as.factor(site_names))
  for (i in 1:length(tmp)) {
    s<-tmp[i]
    if (nchar(s)<9) {   # only change values where name length is too short ()
      ss<-strsplit(as.character(tmp[i]),"-")
      s1<-ss[[1]][1]
      s2<-ss[[1]][2]
      if (length(x=grep("[A-Z]",unlist(strsplit(toupper(s2),""))))==0) 
      {
        tmp[i]<-paste(s1, formatC(as.numeric(s2), width=5, flag="0"), sep="-")
      }
    }
  }
  levels(site_names)<-tmp
  
  return(site_names)
} #SiteNumLeadingZeros



# GENERAL FUNCTIONS -------------------------------------------------------
#We have changed the way we number our segments over the course of NCRMP monitoring.
#Historically segments 0,5,10 and 15 were converted to 1,3,5,7 when migrated to oracle, We found this confusing and are converting them back to 0-15
ConvertSegNumber<-function(data){
  data$SEGMENT<-as.factor(data$SEGMENT)
  data<-data %>% mutate(SEGMENT.new=recode(SEGMENT,
                                           `1`="0",
                                           `3`="5",
                                           `5`="10",
                                           `7`="15",
                                           `NA`="NA"))
  return(data$SEGMENT.new)
}

#Generate depth bin column from DB_RZ column (depth bin/reef zone)
Generate_DB<-function(data){
  data$DB_RZ<-as.factor(data$DB_RZ)
  data<-data %>% mutate(DEPTH_BIN=recode(DB_RZ,
                                         `BA`="All",
                                         `BM`="Mid",
                                         `FM`="Mid",
                                         `LM`="Mid",
                                         `BS`="Shallow",
                                         `FS`="Shallow",
                                         `LS`="Shallow",
                                         `FD`="Deep",
                                         `LD`="Deep"))
  return(data$DEPTH_BIN)
}



##Calcuate segment and transect area and add column for transect area for methods c,e,f
Transectarea<-function(data){
  data$SEGAREA<-data$SEGLENGTH*data$SEGWIDTH # Calculate segment area
  
  #Calculate total transect area then merge back to a dataframe
  s.df<-ddply(data, .(MISSIONID,REGION,ISLAND,OBS_YEAR,SITE,TRANSECT,SEGMENT,SITEVISITID),
              summarise,
              SEGAREA=median(SEGAREA))
  tr.df<-ddply(s.df, .(MISSIONID,REGION,ISLAND,OBS_YEAR,SITE,TRANSECT,SITEVISITID),
               summarise,
               TRANSECTAREA=sum(SEGAREA))
  
  data<-left_join(data,tr.df, by=c("MISSIONID","REGION","ISLAND","OBS_YEAR","SITE","SITEVISITID","TRANSECT"))
  
  
  return(data$TRANSECTAREA)
}


####Functions for benthic summary metrics

#This function calculates total area surveyed per transect
Calc_SurveyArea_By_Transect<-function(data){
  
  tr.df<-ddply(data, .(SITE,SITEVISITID,TRANSECT),
               summarise,
               TRANSECTAREA=median(TRANSECTAREA))
  
  return(tr.df)
}


#This function calculates total area surveyed per site
Calc_SurveyArea_By_Site<-function(data){
  
  tr.df<-ddply(data, .(SITE,TRANSECT,SITEVISITID),
               summarise,
               TRANSECTAREA=median(TRANSECTAREA))
  
  tr.df2<-ddply(tr.df, .(SITE,SITEVISITID),
                summarise,
                TRANSECTAREA=sum(TRANSECTAREA))
  return(tr.df2)
}

## TRANSECT LEVEL SUMMARY FUNCTIONS #######

#This function calculates colony density at the transect scale by first calculating the total survey area (using Calc_SurveyArea_By_Transect) then calcuating colony density
Calc_ColDen_Transect<-function(data, grouping_field="GENUS_CODE"){
  #require(tidyr)
  
  data$GROUP<-data[,grouping_field] #assign a grouping field for taxa
  
  #Remove Tubastrea
  data$S_ORDER<-ifelse(data$GROUP=="TUSP",NA,as.character(data$S_ORDER))
  data$GROUP<-ifelse(data$GROUP=="TUSP","AAAA",as.character(data$GROUP))
  
  
  #Calculate # of colonies for each variable. You need to have S_ORDER and Fragment here so you can incorporate zeros properly later in the code
  a<-ddply(data, .(METHOD,SITE,SITEVISITID,TRANSECT, S_ORDER,GROUP,Fragment),
           summarise,
           ColCount=length(COLONYID)) #change to count
  
  #Convert from long to wide and insert 0s for taxa that weren't found at each site.
  #ca<-dcast(a, formula=SITE + SITEVISITID +TRANSECT +Fragment+S_ORDER~ GROUP, value.var="ColCount",fill=0)
  ca0=spread(a,key=GROUP,value=ColCount,fill=0) #Keepin' it TIDYR
  data.cols<-names(ca0[7:dim(ca0)[2]]) #define your data coloumns- note you need to index by column number not by names-you may have situations where there are no AAAA
  field.cols<-c("METHOD","SITE", "SITEVISITID", "TRANSECT","Fragment") #define field columns
  
  ### Drop all fragments, but don't drop a fragment-only transect... ###
  #change colony counts for fragments to 0 so that we account for the transects that only had fragments
  ca0[which(ca0$Fragment <0), data.cols]<-0
  
  #At this point you will have multiple rows for each site/transect so sum data by site and transect. This will help you properly insert 0s
  field.cols<-c("METHOD","SITE", "SITEVISITID", "TRANSECT")
  ca1<-aggregate(ca0[,data.cols], by=ca0[,field.cols], sum)
  rm(list='ca0')
  
  #Create a list of scleractinian taxa that are in the dataframe as columns then sum across just those taxa to get total scl
  b<-subset(data,S_ORDER=="Scleractinia");taxalist<-as.character(unique(b$GROUP))
  ca1$SSSS<-rowSums(ca1[,taxalist,drop=FALSE]) #calculate total colony density
  ca2 <- tidyr::gather(ca1, GROUP, ColCount, names(ca1[5:dim(ca1)[2]]), factor_key=TRUE) #convert wide to long format
  rm(list='ca1')
  
  #Remove everything that isn't a scleractinian
  taxalist2<-c(taxalist,"SSSS")
  ca3<-ca2[ca2$GROUP %in% taxalist2,]
  rm(list='ca2')
  
  #Calculate transect area surveyed and colony density***
  #trarea<-Calc_SurveyArea_By_Transect(data)
  uTA=unique(data[,c("METHOD","SITE","SITEVISITID","TRANSECT","TRANSECTAREA")])
  out<-join(uTA,ca3, by=c("METHOD","SITE","SITEVISITID","TRANSECT"))
  out$ColDen<-out$ColCount/out$TRANSECTAREA
  colnames(out)[which(colnames(out) == 'GROUP')] <- grouping_field #change group to whatever your grouping field is.
  
  return(out)
}



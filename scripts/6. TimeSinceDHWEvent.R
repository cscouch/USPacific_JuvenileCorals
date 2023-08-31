#This script reads in the raw DHW data that were merged with environmental summary data set
#The script calculates the time between the last survey date and the previously most recent 8 DHW event 

#BEFORE RUNNING SCRIPT- 
#The raw Coral Reef Watch DHW data files are very large and can not be stored in github. To run this script, you will need to
#navigate to https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_dhw_v1_0.html to download the raw data and store them locally.

#This script uses a for loop to read in the raw data file for each island individually then remove the df from the workspace before reading in the next- you can run this script on your laptop


rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/USPacific_JuvenileCorals/"))

library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(lubridate)
library(ggplot2)


#Define the list of island files you will be loading -1 folder with 1 DHW file for each island stored as .Rdata
file_list <- list.files("M:/Environmental Data Summary/DataDownload/Degree_Heating_Weeks/RAW")
# file.miss<-c("Swains_raw_Degree_Heating_Weeks.RData","Kauai_raw_Degree_Heating_Weeks.RData","Kingman_raw_Degree_Heating_Weeks.RData",
#              "Molokai_raw_Degree_Heating_Weeks.RData","Niihau_raw_Degree_Heating_Weeks.RData","Oahu_raw_Degree_Heating_Weeks.RData",
#              "Palmyra_raw_Degree_Heating_Weeks.RData")
#file_list<-("Lisianski_raw_Degree_Heating_Weeks.RData")
#load("M:/Environmental Data Summary/DataDownload/Degree_Heating_Weeks/RAW/Kure_raw_Degree_Heating_Weeks.RData")

#Read in raw juvenile data (this file is used to identify which sites you want to extract DHW for- the file needs to have a date column)
juvdata<-read.csv("Data/outputs/CoralBelt_Juveniles_raw_CLEANED.csv")
#juvdata_recent<-dplyr::filter(juvdata,OBS_YEAR>=2017) #Only use most recent survey years

juvdata <- mutate_if(juvdata,
                is.character,
                stringr::str_replace_all, pattern = " ", replacement = "_")
levels(as.factor(juvdata$ISLAND))

#Calculate minimum positive value
minpositive = function(x) min(x[x > 0])


# Calculate time since last DHW  4 and 8 event ----------------------------

TimeSinceDHW8<-function(data,jwd){

  #Remove columns with column names = NA and rows with NaN
  data <- data[!is.na(names(data))]
  data <- na.omit(data)
  
  #Subset data to make it easier to work with
  dhw<-pivot_longer(data, cols = 6:ncol(data), names_to = "DATE_", values_to = "DHW")
  
  dhw$DATE_<-as.Date(dhw$DATE_, format = "%Y-%m-%d")
  head(dhw)
  
  dhw4<-as.data.frame(dplyr::filter(dhw,DHW >=4))
  class(dhw4$DATE_)
  class(dhw4$DHW4_DATE)
  head(dhw4)

  dhw8<-as.data.frame(dplyr::filter(dhw,DHW >=8));dhw8<-subset(dhw8,select= -c(LATITUDE_LOV,LONGITUDE_LOV))
  class(dhw8$DATE_)
  head(dhw8)
  

  #Only include juvenile sites
  jwd.meta<-ddply(jwd,.(SITE,SITEVISITID,DATE_),summarize,x=mean(SITEVISITID,na.rm=T))
  jwd.meta<-subset(jwd.meta,select=-c(x));colnames(jwd.meta)<-c("SITE","SITEVISITID","J.DATE")
  jwd.meta$J.DATE<-as.Date(jwd.meta$J.DATE, format = "%Y-%m-%d")
  class(jwd.meta$J.DATE)
  head(jwd.meta)

  #Merge 4 DHW data and Juvenile data
  dhw4<-merge(jwd.meta,dhw4,by=c("SITE","SITEVISITID"),all.x=T) #tried using left_join- it doesn't like the date columns
  head(dhw4)
  
  #Calculate difference between dates
  dhw4$Tdiff4<-difftime(as.Date(dhw4$J.DATE) ,as.Date(dhw4$DATE_) , units = c("weeks")) 
  dhw4$YearSinceDHW4<-as.numeric(dhw4$Tdiff/52) #transform Tdiff into years

  #Merge 8 DHW data and Juvenile data
  dhw8<-merge(jwd.meta,dhw8,by=c("SITE","SITEVISITID"),all.x=T) #tried using left_join- it doesn't like the date columns

  #Calculate difference between dates
  dhw8$Tdiff8<-difftime(as.Date(dhw8$J.DATE) ,as.Date(dhw8$DATE_) , units = c("weeks")) 
  dhw8$YearSinceDHW8<-as.numeric(dhw8$Tdiff/52) #transform Tdiff into years
  head(dhw8)
  
  
  #Select most recent date for 4DHW
  df4<-ddply(dhw4,.(SITE,SITEVISITID,ISLAND),summarize,
            YearSinceDHW4=minpositive(YearSinceDHW4)) #Identify minimum positive value (aka most recent date prior to survey)
  df8<-ddply(dhw8,.(SITE,SITEVISITID,ISLAND),summarize,
             YearSinceDHW8=minpositive(YearSinceDHW8))
  head(df4);head(df8)
  
  df<-merge(df4,df8,by=c("SITE","SITEVISITID","ISLAND"),all=T) #merge DHW4 and DHW8
  
  
  #Note- NA values are sites where no DHW events have ever been observed, Inf values indicate sites that saw DHW events only after the most recent surveys
  #Change Inf to NA
  df<-df%>% mutate_if(is.numeric, ~ifelse(abs(.) == Inf,NA,.))
    
  return(df)

}

df.all<-NULL
for(i in 1:length(file_list)){
  
  #i = 1
  filename = file_list[i]
  load(filename)
  df<-TimeSinceDHW8(df_i,juvdata)
  rm(list='df_i') #remove the dhw df for a given island from work space- files are too large to load all at once
  df.all<-rbind(df.all,df)
}
View(df.all)

#df.all3<-rbind(df.all,df.all2)

write.csv(df.all,file="Data/outputs/Juvenile_TimeSinceDHW4_8_v2.csv")



# Calculate Date of Peak DHW for each year --------------------------------


DHWPeak<-function(data,jwd){
  
  #Remove columns with column names = NA and rows with NaN
  data <- data[!is.na(names(data))]

  #Subset data to make it easier to work with
  dhw<-pivot_longer(data, cols = 6:ncol(data), names_to = "DATE_", values_to = "DHW")
  
  dhwR<-filter(dhw,DATE_>="2013-01-01")

  # dhw$DATE_<-lubridate::ymd(dhw$DATE_)
  # head(dhw)
  
  dhwR$OBS_YEAR<-year(dhwR$DATE_)
  
  #Identify Date of peak DHW
  peak <- dhwR %>%
    group_by(SITE,OBS_YEAR) %>%
    filter(DHW == max(DHW)) %>%
    select(dhw=DHW, PeakDHW=DHW, DATE_) %>%
    filter(DATE_==min(DATE_))
  peak$PeakDate<-peak$DATE_
  peak <- peak[,c("SITE","OBS_YEAR","PeakDHW","PeakDate")]
  
  peak<-subset(peak,SITE %in% c(jwd$SITE))
  
  return(peak)
  
}

peak.all<-NULL
for(i in 1:length(file_list)){
  
  #i = 1
  filename = file_list[i]
  load(filename)
  df<-DHWPeak(df_i,juvdata)
  rm(list='df_i') #remove the dhw df for a given island from work space- files are too large to load all at once
  peak.all<-rbind(peak.all,df)
}
View(peak.all)

write.csv(peak.all,file="Data/outputs/Juvenile_PeakDHW.csv")

  
data=df_i
jwd=juvdata


# Create a timeline of heat stress for Pacific NCRMP regions  -------------

#Calculate average DHW for each date  
HS_Timeline<-function(data,jwd){
  
  #Remove columns with column names = NA and rows with NaN
  data <- data[!is.na(names(data))]
 
  #Only include juvenile sites
  jwd.meta<-ddply(jwd,.(SITE,SITEVISITID,DATE_),summarize,x=mean(SITEVISITID,na.rm=T))
  jwd.meta<-subset(jwd.meta,select=-c(x));colnames(jwd.meta)<-c("SITE","SITEVISITID","J.DATE")
  jwd.meta$J.DATE<-as.Date(jwd.meta$J.DATE, format = "%Y-%m-%d")
  class(jwd.meta$J.DATE)
  head(jwd.meta) 
  
  dhw<-pivot_longer(data, cols = 6:ncol(data), names_to = "DATE_", values_to = "DHW")
  
  dhwR<-filter(dhw,DATE_>="2008-01-01")
  
  #Merge 4 DHW data and Juvenile data
  dhwR<-merge(jwd.meta,dhwR,by=c("SITE","SITEVISITID"),all.x=T) #tried using left_join- it doesn't like the date columns
  head(dhwR)
  
  dhwM <- dhwR %>%
    group_by(ISLAND,DATE_) %>%
    summarise(MeanDHW=mean(DHW,na.rm=TRUE))
  
  return(dhwM)
  
}

hs.all<-NULL
for(i in 1:length(file_list)){
  
  #i = 1
  filename = file_list[i]
  load(filename)
  df<-HS_Timeline(df_i,juvdata)
  rm(list='df_i') #remove the dhw df for a given island from work space- files are too large to load all at once
  hs.all<-rbind(hs.all,df)
}
View(hs.all)


write.csv(hs.all,file="Data/outputs/Juvenile_HSTimeline.csv",row.names=FALSE)


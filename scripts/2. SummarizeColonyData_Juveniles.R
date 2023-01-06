#This script summarizes juvenile data from NCRMP 2013-2019 at the site, stratum, island and sector level
#It also identifies which sectors and strata have been surveyed in all years
#It calculate delta density

# Using R version 4.1.0 (2021-05-18)

rm(list=ls())

#LOAD LIBRARY FUNCTION ...
source("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/scripts/Functions_Juveniles.R")

#LOAD DATA
jwd<-read.csv("T:/Benthic/Projects/Juvenile Project/Data/CoralBelt_Juveniles_raw_CLEANED.csv")

#Tweaks before calculating Site-level data-------------------------------------------------
#Colony fragments and scleractinans are subsetted in the functions 
#Add a column for adult fragments so we can remove them from the dataset later (-1 indicates fragment)
jwd$Fragment <- 0 # you need to add this column so that you can use the site level functions correctly
jwd$DATE_ <- ymd(jwd$DATE_)
jwd$METHOD<-"DIVER"
jwd$ANALYST<-jwd$DIVER
jwd$SEGAREA<-jwd$SEGLENGTH*jwd$SEGWIDTH

# #Round depth to nearest 4 decimal places to avoid mergings issues
# jwd$MIN_DEPTH_M<-round(jwd$MIN_DEPTH_M, 4)
# jwd$MAX_DEPTH_M<-round(jwd$MAX_DEPTH_M, 4)
# jwd$LATITUDE<-round(jwd$LATITUDE,7)
# jwd$LONGITUDE<-round(jwd$LONGITUDE,5)


#Fix Errors in Data:
#2 sites with incorrect Reef zone and depth bin
jwd$REEF_ZONE<-ifelse(jwd$SITE=="HAW-04285","Forereef",as.character(jwd$REEF_ZONE))
jwd$DEPTH_BIN<-ifelse(jwd$SITE=="FFS-04155","Shallow",as.character(jwd$DEPTH_BIN))

#3 Maug sites miscoded as Lagoon, should be FRF
mau.sites<-c("MAU-00603","MAU-00539","MAU-00551")
jwd$REEF_ZONE<-ifelse(jwd$SITE %in% mau.sites,"Forereef",as.character(jwd$REEF_ZONE))

#Remove 2 sites that weren't surveyed for juveniles
jwd<-jwd[!(jwd$SITE %in% c("OFU-01012","PAG-00596")),]

#Remove special missions (not NCRMP surveys)
#Change all special missions to exclude flag =-1, right now they are 0. Then exclude these sites
jwd$MISSIONID<-as.factor(jwd$MISSIONID)
levels(jwd$MISSIONID)
jwd<-jwd[!jwd$MISSIONID %in% c("MP1410","MP1512","MP1602","SE1602","MP2006"),] 

#Exclude PRIA 2017 sites because we want similar time intervals following bleaching events for all regions
jwd$Year_Island<-paste(jwd$OBS_YEAR,jwd$ISLAND,sep="_")
jwd<-jwd[!jwd$Year_Island %in% c("2017_Baker","2017_Jarvis","2017_Howland"),] 

jwd<-droplevels(jwd);levels(jwd$MISSIONID) #special missions should be gone
View(jwd)

#Dealing with colonies <1cm- divers haven't recorded sub cm colonies consistently through time.
#Change colonies that are <1cm to NA. I'm not subsetting these data because I need to keep the placeholder in the dataframe in case a site only had colonies <1cm or >5cm
View(subset(jwd,COLONYLENGTH<1))
nrow(subset(jwd,COLONYLENGTH<1))
nrow(jwd)
jwd$S_ORDER<-ifelse(jwd$COLONYLENGTH<1|jwd$COLONYLENGTH==5,NA,as.character(jwd$S_ORDER))
jwd$GENUS_CODE<-ifelse(jwd$COLONYLENGTH<1|jwd$COLONYLENGTH==5,"AAAA",as.character(jwd$GENUS_CODE))
jwd$TAXONCODE<-ifelse(jwd$COLONYLENGTH<1|jwd$COLONYLENGTH==5,"AAAA",as.character(jwd$TAXONCODE))
jwd$COLONYLENGTH<-ifelse(jwd$COLONYLENGTH<1|jwd$COLONYLENGTH==5,NA,jwd$COLONYLENGTH)
jwd$TAXONNAME<-ifelse(jwd$COLONYLENGTH<1|jwd$COLONYLENGTH==5,NA,as.character(jwd$TAXONNAME))

nrow(subset(jwd,COLONYLENGTH>1))
nrow(subset(jwd,COLONYLENGTH<1)) #should be 0



#Create list of sites with metadata - doesn't include depth, lat and long because there are issues with slightly different decimal places that cause issues with merging. merge in later
SURVEY_SITE<-c("METHOD","MISSIONID","DATE_","SITEVISITID", "OBS_YEAR", "REGION", "REGION_NAME", "ISLAND","ISLANDCODE","SEC_NAME", "SITE","HABITAT_CODE","REEF_ZONE",
"DEPTH_BIN")
survey_site<-unique(jwd[,SURVEY_SITE])

n_occur <- data.frame(table(survey_site$SITE)) 
n_occur[n_occur$Freq > 1,]

# Generate Juvenile Density at the TRANSECT & SITE-LEVEL BY GENUS--------------------------------------------------
jcd.gen<-Calc_ColDen_Transect(jwd,"GENUS_CODE"); colnames(jcd.gen)[colnames(jcd.gen)=="ColCount"]<-"JuvColCount";colnames(jcd.gen)[colnames(jcd.gen)=="ColDen"]<-"JuvColDen";colnames(jcd.gen)[colnames(jcd.gen)=="TRANSECTAREA"]<-"TRANSECTAREA_j"

#Drop 2nd transect since we only surveyed 1 transect after 2017 and the 2nd transect wasn't surveyed consistently prior to 2018
jcd.gen$TRANSECT<-as.factor(jcd.gen$TRANSECT)

site.data.gen<-dplyr::filter(jcd.gen,TRANSECT %in% c("1","3")) #subseting first transect (different transect numbering was used over the years 1 & 3 refer to the 1st transects- I know it doesn't make sense)
summary(site.data.gen$TRANSECT)

site.data.gen2<-dplyr::select(site.data.gen, -c(TRANSECT)) #Create a copy (it takes a while to run the transect function above)

# Merge Site level data with sectors file and export site data ------------
sectors<-read.csv("C:/Users/Courtney.S.Couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SectorArea_Juveniles.csv", stringsAsFactors=FALSE)

#Merge together survey meta data and sector area files and check for missmatches 
meta<-left_join(survey_site,sectors)
nrow(survey_site)
nrow(meta)
nrow(subset(site.data.gen2,GENUS_CODE=="SSSS")) #1611 sites

head(survey_site)


#Merge site level data and metadata
site.data.gen2<-left_join(site.data.gen2,meta)
site.data.gen2$Juvpres.abs<-ifelse(site.data.gen2$JuvColDen>0,1,0) #add presence/absense column 
nrow(subset(site.data.gen2,GENUS_CODE=="SSSS")) #should be 1611
site.data.gen2[which(is.na(site.data.gen2$AREA_HA)),] #The NA values are from special missions and maug Lagoon that will be dropped later in the script - ok
nrow(meta)

#Make tweaks to pooling sector pooling and drop specific islands because they were not surveyed more than once or didn't have enough sampling across years 
isl.drop<-c("Alamagan","Midway","Maro","Johnston","Guguan","Agrihan") 
site.data.gen2<-site.data.gen2[!site.data.gen2$ISLAND %in% c(isl.drop),] #Remove islands that don't have enough strata sampled across the years 

#Not enough sampling in each sector- pool them together
site.data.gen2$SEC_NAME<-ifelse(site.data.gen2$SEC_NAME %in% c("TAU_OPEN","TAU_SANCTUARY"),"TAU",as.character(site.data.gen2$SEC_NAME))
site.data.gen2$SEC_NAME<-ifelse(site.data.gen2$SEC_NAME %in% c("SWA_OPEN","SWA_SANCTUARY"),"SWA",as.character(site.data.gen2$SEC_NAME))
View(site.data.gen2)


#Convert Protected Reef Slope to Forereef and Subset just Forereef sites
site.data.gen2$REEF_ZONE<-ifelse(site.data.gen2$REEF_ZONE=="Protected Slope","Forereef",as.character(site.data.gen2$REEF_ZONE))
site.data.gen2<-subset(site.data.gen2,REEF_ZONE=="Forereef") #only include forereef

#Add up NH (number of 250m2 grid areas in stratum) for forereef and protected reef sites
tmp<-unique(site.data.gen2[,c("OBS_YEAR","SEC_NAME","REEF_ZONE","DEPTH_BIN","NH")])

new.NH<-ddply(tmp,.(OBS_YEAR,SEC_NAME,REEF_ZONE,DEPTH_BIN),
              summarize, new.NH=sum(NH))
nrow(site.data.gen2) #double check that the same number of rows are present before and after joining
site.data.gen2<-left_join(site.data.gen2,new.NH)
site.data.gen2<-subset(site.data.gen2,select=-c(NH)) #remove old NH column
site.data.gen2<-site.data.gen2 %>% dplyr::rename(NH=new.NH) #rename NH
nrow(site.data.gen2)#double check that the same number of rows are present before and after joining

table(site.data.gen2$ISLAND,site.data.gen2$OBS_YEAR)


#Change Region Names -breaking up marianas and PRIAs into subregions because of broad geography and thermal history
site.data.gen2$REGION<-ifelse(site.data.gen2$ISLAND %in% c("Maug", "Asuncion", "Alamagan", "Pagan", "Agrihan", "Guguan", "Sarigan","Farallon de Pajaros")
                              ,"NMI", as.character(site.data.gen2$REGION))
site.data.gen2$REGION<-ifelse(site.data.gen2$ISLAND %in% c("Saipan", "Tinian", "Aguijan", "Rota", "Guam")
                              ,"SMI", as.character(site.data.gen2$REGION))
site.data.gen2$REGION<-ifelse(site.data.gen2$ISLAND %in% c("Howland","Baker")
                              ,"PHOENIX", as.character(site.data.gen2$REGION))
site.data.gen2$REGION<-ifelse(site.data.gen2$ISLAND =="Wake"
                              ,"WAKE", as.character(site.data.gen2$REGION))
site.data.gen2$REGION<-ifelse(site.data.gen2$ISLAND %in% c("Kingman","Palmyra","Jarvis")
                              ,"LINE", as.character(site.data.gen2$REGION))

site.data.gen2$STRATANAME<- paste(site.data.gen2$SEC_NAME,site.data.gen2$REEF_ZONE,site.data.gen2$DEPTH_BIN,sep="_")
site.data.gen2$REGION_YEAR<-paste(site.data.gen2$REGION,site.data.gen2$OBS_YEAR,sep="_")
site.data.gen2$DB_RZ<- paste(site.data.gen2$REEF_ZONE,site.data.gen2$DEPTH_BIN,sep="_")

#Remove the 2014 NWHI data - we have very low sampling in 2014 & we do not have benthic cover data for 2014 
site.data.gen2<-site.data.gen2[site.data.gen2$REGION_YEAR !="NWHI_2014",]

length(unique(site.data.gen2$SITE))


#Calculate survey weights (inverse proportion weighting)
w.df<-ddply(site.data.gen2,.(OBS_YEAR,SEC_NAME,REEF_ZONE,DEPTH_BIN,NH),
            summarize,
            n=length(unique(SITE))) #quantify # of sites/stratum

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(site.data.gen2,w.df) #merge weights with site-level data
head(site.sw)
site.swS<-dplyr::filter(site.sw,GENUS_CODE=="SSSS")#only include total scleractinans

length(unique(site.swS$SITE))

#remove strata that have less than 2 sites
site.swS<-subset(site.swS,n>1)
summary(site.swS$n)
length(unique(site.swS$SITE))


#Check for duplicate sites
site.swS %>% 
  group_by(SITE) %>% 
  filter(n()>1)

#We have 1405 sites for spatial and correlative analysis

#Merge Lat and Long and depths back in
survey_master<-read.csv("C:/Users/courtney.s.couch/Documents/GitHub/USPacific_JuvenileCorals/SupportFiles/SURVEY MASTER_Juveniles.csv")

#Use SM coordinates-some coordinates are wrong in data and need to be updated
colnames(survey_master)[colnames(survey_master)=="LATITUDE_LOV"]<-"LATITUDE" #Change column name- we will eventually change this column back to "taxoncode" after we modify the spcode names to match the taxalist we all feel comfortable identifying
colnames(survey_master)[colnames(survey_master)=="LONGITUDE_LOV"]<-"LONGITUDE" #Change column name- we will eventually change this column back to "taxoncode" after we modify the spcode names to match the taxalist we all feel comfortable identifying

colnames(survey_master)[colnames(survey_master)=="new_MIN_DEPTH_M"]<-"MIN_DEPTH_M" #Change column name
colnames(survey_master)[colnames(survey_master)=="new_MAX_DEPTH_M"]<-"MAX_DEPTH_M" #Change column name

#merge colony data and survey master
site.swS<-left_join(site.swS, survey_master[,c("SITEVISITID","SITE","LATITUDE","LONGITUDE","MIN_DEPTH_M","MAX_DEPTH_M")])
View(site.swS)

#Save site-level data
write.csv(site.swS,file="T:/Benthic/Projects/Juvenile Project/Data/JuvProject_SITE_weights_AllYears.csv",row.names = F)


# GENERATE DATA FOR TEMPORAL ANALYSIS ---------------------------------------------------
site.swS$STRATANAME<- paste(site.swS$SEC_NAME,site.swS$REEF_ZONE,site.swS$DEPTH_BIN,sep="_") #create strata column
st.list<-ddply(site.swS,.(OBS_YEAR,REGION,ISLAND,SEC_NAME,STRATANAME),summarize,n=length(unique(SITE)))#number of sites/stratum

#Generate list of strata that were surveyed in all years for a given region and had at least 2 sites/stratum
st.list_w<-reshape2::dcast(st.list, formula=REGION+ISLAND+SEC_NAME+STRATANAME~ OBS_YEAR, value.var="n",fill=0)
dCOLS<-c("2013","2014","2015","2016","2017","2018","2019")
st.list_w$year_n<-rowSums(st.list_w[,dCOLS] > 0, na.rm=T) #count # of years of data
st.list_w2<-subset(st.list_w,REGION %in% c("NMI","SMI","LINE","PHOENIX","WAKE","SAMOA") & year_n>=2) #subset strata that were surveyed twice for these regions
st.list_w3<-subset(st.list_w,REGION %in% c("NWHI","MHI") & year_n>=3) #subset strata that were surveyed three times for these regions
st.list_w4<-rbind(st.list_w2,st.list_w3) #combine strata list to use for temporal analysis

head(st.list_w4);st.list_w4<-droplevels(st.list_w4) #generate the list

data.gen_temp<-site.swS[site.swS$STRATANAME %in% c(st.list_w4$STRATANAME),] #Subset juv data to only include strata of interest 
isl.drop<-c("Maui","Niihau","Kure","French Frigate","Midway","Pearl & Hermes","Lisianski") #excluding all NWHI islands- poor sampling across years
data.temporal<-data.gen_temp[!data.gen_temp$ISLAND %in% c(isl.drop),] #Remove islands that don't have enough strata sampled across the years 
length(unique(data.temporal$SITE))

#We have 1010 sites for the temporal analysis- note,only strata that were surveyed in years were include, BUT not all possible strata for a given sector are included. E.g.
#we may be missing the shallow stratum because it was only sampled in one year. Care should be taken when comparing temporal patterns between islands

write.csv(data.temporal,file="T:/Benthic/Projects/Juvenile Project/JuvDen_Temporal.csv",row.names = FALSE)



# TEMPORAL Patterns in Juveniles- Analysis and Plots --------------------------------------------

#Use survey package to calculate mean SE and conduct statistical analyses
data.temporal$OBS_YEAR<-as.factor(data.temporal$OBS_YEAR)

#Create contactenated Strata variable and establish survey design with survey weights
data.temporal$Strat_conc<-paste(data.temporal$OBS_YEAR, data.temporal$REGION,data.temporal$ISLAND,data.temporal$SEC_NAME,data.temporal$DB_RZ,sep = "_")
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=data.temporal)

#Calculate regional mean and SE
temp_Rmean<-svyby(~JuvColDen,~OBS_YEAR+REGION,des,svymean)

#Test fixed effects of region and year
modR<-svyglm(JuvColCount ~ REGION*OBS_YEAR, design=des,offset= TRANSECTAREA_j, family="poisson")
anova(modR)

#Run separate post hoc tests for each region to test for differences between years- I don't care about comparing all possible combinations of year and region
library(multcomp)

mhi.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="MHI"))
mhi<-svyglm(JuvColCount ~ OBS_YEAR, design=mhi.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(mhi, mcp(OBS_YEAR="Tukey"))) 

wake.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="WAKE"))
wake<-svyglm(JuvColCount ~ OBS_YEAR, design=wake.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(wake, mcp(OBS_YEAR="Tukey"))) 

ph.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="PHOENIX"))
ph<-svyglm(JuvColCount ~ OBS_YEAR, design=ph.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(ph, mcp(OBS_YEAR="Tukey"))) 

line.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="LINE"))
l<-svyglm(JuvColCount ~ OBS_YEAR, design=line.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(l, mcp(OBS_YEAR="Tukey"))) 

sm.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="SMI"))
sm<-svyglm(JuvColCount ~ OBS_YEAR, design=sm.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(sm, mcp(OBS_YEAR="Tukey"))) 

nm.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="NMI"))
nm<-svyglm(JuvColCount ~ OBS_YEAR, design=nm.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(nm, mcp(OBS_YEAR="Tukey")))

sam.des<-svydesign(id=~1, strata=~Strat_conc, weights=~sw,data=subset(data.temporal,REGION=="SAMOA"))
sam<-svyglm(JuvColCount ~ OBS_YEAR, design=sam.des,offset= TRANSECTAREA_j, family="poisson")
summary(glht(sam, mcp(OBS_YEAR="Tukey"))) 


#Calculate adjusted pvalues for multiple test corrections
pvals<-c(0.54508,0.01609,0.00325,0.366,0.312,0.000001,0.176,0.000041,0.305)
round(p.adjust(pvals, "BH"), 3) #0.568 0.082 0.010 0.447 0.447 0.000 0.447 0.001 0.568


# PLOTTING Figure 2 ----------------------------------------------------------------

#bar plot of juv by region by year with post hoc tests 
temp_Rmean$REGION <- factor(temp_Rmean$REGION, levels = c("MHI","NMI","WAKE","SMI","SAMOA","LINE","PHOENIX"))
temp_Rmean$OBS_YEAR<-as.factor(temp_Rmean$OBS_YEAR)
#Add Posthoc groupings from glms
temp_Rmean<- temp_Rmean[order(temp_Rmean$REGION),];temp_Rmean
temp_Rmean$sig<-c("ab","a","b","a","b","","","","","","","a","b","","")

#scale_fill_manual(values = c("#CC79A7","#D55E00","#E69F00","#F0E442","#009E73","#56B4E9","#0072B2","#999999")) +
  
p8 <- ggplot(temp_Rmean, aes(x=OBS_YEAR, y=JuvColDen,fill=REGION)) +
  #geom_bar(stat = "identity", position = position_dodge2(preserve='single'), width = 1, color="black") +
  geom_errorbar(aes(y=JuvColDen, x=OBS_YEAR,ymin=JuvColDen-se, ymax=JuvColDen+se), width=.2)+
  geom_point(color="black",pch=21,size=4)+
  facet_grid(~REGION, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 12),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(angle = 90,size=14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("#D55E00","#E69F00","#F0E442","#009E73","#56B4E9","#0072B2","#999999")) +
  labs(x="",y=expression(paste("Mean Juvenile Colonies ",m^-2)))+
  scale_y_continuous(expand = c(0,0), limits = c(0,16)) +
  geom_text(aes(x=OBS_YEAR,y=JuvColDen+se,label=sig, group = REGION),
            position = position_dodge(),
            vjust = -0.5) 
p8

ggsave(plot=p8,file="T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/Figure2.jpg",width=10,height=5)



# Spatial Trends in Juveniles averaged across years --------------------------------------------

table(site.swS$ISLAND,site.swS$DB_RZ)
nrow(site.swS) #total number of sites used in the spatial analysis

#Use survey package to calculate mean SE and conduct statistical analyses
site.swS$OBS_YEAR<-as.factor(site.swS$OBS_YEAR)
site.swS$REGION<-as.factor(site.swS$REGION)

#Establish survey design
site.swS$Strat_conc<-paste(site.swS$OBS_YEAR, site.swS$REGION,site.swS$ISLAND,site.swS$SEC_NAME,site.swS$DB_RZ,sep = "_")
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.swS)

#Calculate regional mean and SE
spatial_Rmean<-svyby(~JuvColDen,~REGION,des,svymean)

#Test fixed effects of region and year
modR<-svyglm(JuvColCount ~ REGION, design=des,offset= TRANSECTAREA_j, family="poisson")
summary(modR)
tuk2<-glht(modR, mcp(REGION="Tukey")) 
tuk.cld2 <- cld(tuk2)
sig<-c("bd","b","a","ab","c","d","a","c")

spatial_Rmean<-cbind(spatial_Rmean,sig)


#svystdres(modR,stvar="DB_RZ",doplot=TRUE)
null.mod<-svyglm(JuvColCount ~ 1, design=des,offset= TRANSECTAREA_j, family="poisson")
anova(null.mod,modR)


#bar plot of juv by region by year with post hoc tests 
spatial_Rmean$REGION <- factor(spatial_Rmean$REGION, levels = c("NWHI","MHI","NMI","WAKE","SMI","SAMOA","LINE","PHOENIX"))
#Add Posthoc groupings from glms
spatial_Rmean<- spatial_Rmean[order(spatial_Rmean$REGION),];spatial_Rmean


spatialR <- ggplot(spatial_Rmean, aes(x=REGION, y=JuvColDen,fill=REGION)) +
  #geom_bar(stat = "identity", position = position_dodge2(preserve='single'), width = 1, color="black") +
  geom_errorbar(aes(y=JuvColDen, x=REGION,ymin=JuvColDen-se, ymax=JuvColDen+se), width=.2)+
  geom_point(color="black",pch=21,size=5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 14),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(angle = 90,size=14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  scale_fill_manual(name=NULL, values = c("#CC79A7","#D55E00","#E69F00","#F0E442","#009E73","#56B4E9","#0072B2","#999999")) +
  labs(x="",y=expression(bold(paste("Mean Juvenile Colonies ",m^-2))),title = "B")+
  scale_y_continuous(expand = c(0,0), limits = c(0,16)) +
  geom_text(aes(x=REGION,y=JuvColDen+se,label=sig, group = REGION),
            position = position_dodge(),
            vjust = -0.5,size = 5) 
spatialR

#ggsave(plot=spatialR,file="T:/Benthic/Projects/Juvenile Project/Figures/DensityRegionalSpatial.jpg",width=10,height=5)


spatial_Imean<-svyby(~JuvColDen,~REGION + ISLAND,des,svymean)
spatial_Imean$REGION <- factor(spatial_Imean$REGION, levels = c("NWHI","MHI","NMI","WAKE","SMI","SAMOA","LINE","PHOENIX"))


p10 <- ggplot(spatial_Imean, aes(x=reorder(ISLAND,-JuvColDen), y=JuvColDen,color=REGION)) +
  #geom_bar(stat = "identity", position = position_dodge2(preserve='single'), width = 1, color="black") +
  geom_point(size=2)+
  geom_errorbar(aes(y=JuvColDen, x=ISLAND,ymin=JuvColDen-se, ymax=JuvColDen+se), width=.2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        #legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 12),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(vjust=0,angle = 90)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,25)) +
  scale_color_manual(name=NULL,values = c("#CC79A7","#D55E00","#E69F00","#F0E442","#009E73","#56B4E9","#0072B2","#999999")) +
  labs(x="Island",y=expression(paste("Mean Juvenile Colonies  ",m^-2)))
p10


ggsave(plot=p10,file="T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/FigureS2.jpg",width=10,height=5)


# Plot Pacific-wide Map most recent Juvenile Density -------------------------

##Helpful website for plotting maps with ggplot https://r-spatial.org/r/2018/10/25/ggplot2-sf-2.html

#https://rpubs.com/valentin/pacific-centered-map-voronoi-tessellation


#Identify median lat and long that surveys were conducted for each island and year
lat.sum<-ddply(site.swS,.(REGION,ISLAND),
               summarize,
               Y=median(LATITUDE,na.rm=T),
               X=median(LONGITUDE,na.rm = T))

j.sum<-ddply(spatial_Imean,.(REGION,ISLAND),
             summarize,
             JuvColDen=mean(JuvColDen,na.rm=T))
juv_coords<-left_join(j.sum,lat.sum)
head(juv_coords)


#Theme for maps
theme_set(
  theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 8)))

#Get world spatial polygons from the rnaturalearth package
#Cut out area of world that doesn't include the Pacific and bind E and W Pacific and shift the geographical coordinates for a Pacific view (see website for illustration)
#Note- you will get several warnings about "world is invalid" and an issues with rgeos- ignore these
world <- rnaturalearth::ne_countries(scale = 'medium', returnclass = "sp")
box_cut <- bbox2SP(n = 90, s = -90, w = -150, e = 140, proj4string = world@proj4string) #you can tweak the W and E coords to zoom in and out, but adjust the N and S coords in the function below using ymin and ymax
world_crop <- gDifference(world, box_cut)

pacific_crop <- world_crop %>% 
  st_as_sf() %>% # change from sp to sf object/class
  st_shift_longitude() %>% 
  st_crop(c(xmin = st_bbox(.)[["xmin"]],
            xmax = st_bbox(.)[["xmax"]],
            ymin = -15,
            ymax = 30))

#convert delta juvenile density data to spatial points df
xy <- juv_coords[,c(5,4)]
juv_coords_sp <- SpatialPointsDataFrame(coords = xy, data = juv_coords,
                                        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#Crop the juv spdf
delta_shift <- juv_coords_sp %>% 
  st_as_sf() %>%
  st_shift_longitude() %>% 
  st_crop(c(xmin = 120, xmax = 250, ymin = -50, ymax = 30)) %>% 
  # Also adds the coordinates to be used for labeling with geom_text_repel
  bind_cols(st_coordinates(.) %>% as.data.frame())

#Create an Inset map using a similar process described above
box_cut2 <- bbox2SP(n = 90, s = -90, w = -120, e = 110, proj4string = world@proj4string)
world_crop2 <- gDifference(world, box_cut2)

###https://github.com/r-spatial/sf/issues/1902
#Having issues with the 

pacific_crop2 <- world_crop2 %>% 
  st_as_sf() %>% # change from sp to sf object/class
  st_shift_longitude() %>% 
  st_crop(c(xmin = st_bbox(.)[["xmin"]],
            xmax = st_bbox(.)[["xmax"]],
            ymin = -40,
            ymax = 50))
pacific_crop_bb = st_as_sfc(st_bbox(pacific_crop)) #draw box for the area of main map
pacific_box = st_as_sfc(st_bbox(pacific_crop2)) #draw box for the area of inset map

#plot inset map
insetmap<-ggplot() +
  geom_sf(data = pacific_crop2)+
  geom_sf(data = pacific_crop_bb, fill = NA, color = "black", size = 1.2)+
  geom_sf(data = pacific_box, fill = NA, color = "black", size = 0.4)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank())


#plot main data map with island colors = gradient of juvenile density
deltamap<-ggplot() +
  geom_sf(data = pacific_crop)+ #basemap
  geom_sf(data = delta_shift,aes(fill = JuvColDen), size = 3, shape = 21)+ #data
  geom_text_repel(data = delta_shift, #add island labels 
                  aes(x = X...7, y = Y...8, label = ISLAND),
                  size = 3,
                  fontface = "bold",
                  segment.size = 0.25,
                  box.padding = 0.4,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  seed = 2020-5-16)+
  annotation_scale(location = "bl", width_hint = 0.4)+ #add scale bar
  scale_fill_gradient2(midpoint = 8, #Color scheme
                        high = 'forestgreen',
                        mid = 'yellow2',
                        low = 'red2',
                        na.value = 'gray95',
                        name=expression(paste("Mean Juvenile Colonies  ",m^-2)))+  #you can add a legend title here if you want
  #theme(legend.position = c(0.9,0.15))
  theme(legend.position = "bottom")
  

#Combine main and inset maps
finalmap = ggdraw() +
  draw_plot(deltamap) +
  draw_plot(spatialR, x = 0.05, y = 0.14, width = 0.3, height = 0.3)+
  #draw_plot(insetmap, x = 0.02, y = 0.07, width = 0.3, height = 0.3)
  draw_plot(insetmap, x = 0.77, y = 0.2, width = 0.23, height = 0.23)
  
finalmap

ggsave(plot=finalmap,file="T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/new.Figure1.jpg",width=13,height=9)

#OLD version
#plot main data map with island colors = REGION
delta_shift$REGION <- factor(delta_shift$REGION, levels = c("NWHI","MHI","NMI","WAKE","SMI","SAMOA","LINE","PHOENIX"))
delta_shift<- delta_shift[order(delta_shift$REGION),];delta_shift


deltamap<-ggplot() +
  geom_sf(data = pacific_crop)+ #basemap
  geom_sf(data = delta_shift, aes(fill = REGION),color="black",size = 3, pch=21)+ #data
  geom_text_repel(data = delta_shift, #add island labels 
                  aes(x = X...7, y = Y...8, label = ISLAND),
                  size = 5,
                  fontface = "bold",
                  segment.size = 0.25,
                  box.padding = 0.4,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  seed = 2020-5-16)+
  labs(title="A")+
  scale_fill_manual(name=NULL,values = c("#CC79A7","#D55E00","#E69F00","#F0E442","#009E73","#56B4E9","#0072B2","#999999")) +
  annotation_scale(location = "bl", width_hint = 0.4)+ #add scale bar
  #theme(legend.position = c(0.9,0.15))
  theme(legend.position = "none")


#Combine main and inset maps
finalmap = ggdraw() +
  draw_plot(deltamap) +
  draw_plot(spatialR, x = 0.05, y = 0.075, width = 0.45, height = 0.45)+
  #draw_plot(insetmap, x = 0.02, y = 0.07, width = 0.3, height = 0.3)
  draw_plot(insetmap, x = 0.77, y = 0.1, width = 0.23, height = 0.23)

finalmap

ggsave(plot=finalmap,file="T:/Benthic/Projects/Juvenile Project/Manuscript/Figures/Figure 1.jpg",width=13,height=9)


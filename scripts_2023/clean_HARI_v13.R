# check problem with same value for a series of countri in the barplot of totla infections per country. 

library(gdata)
library(raster)
library(FME)
library(RColorBrewer)
library(usdm)
library(mice)
setwd("/Users/ruchita/Documents/HARI/")
source("scripts_2023/00_Library.r")

# Imputation for countries that did not report a single survey? 
###############################################################
ImpTAG = T

# Import data extracted by RAs
##############################
myT = read.csv("/Users/ruchita/Documents/HARI/May_9_2022_Consolidated_Sheet_tvb.csv")

# Remove columns with empty spaces  
myT = myT[which(names(myT) %in%  c("X","X.1","X.2","X.3","X.4","X.5","X.6","X.7","X.8","X.9","X.10") == F)]

#Conver variables to numeric if neded
myT$Bed = as.numeric(myT$Bed)
myT$StudyPeriod = as.numeric(myT$StudyPeriod)

# read shapefile 
mySHP = shapefile("/Users/ruchita/Documents/HARI/data/World.shp")
mySHP@data = mySHP@data[c("ISO3")]
myFAO = mySHP@data

# aggregate all hospitalization rates 
myFAO$HRNTYR = myT$HRNTYR[match(mySHP@data$ISO3,myT$ISO3)]

# Read latest population data from world bank
myPop2019 = read.csv("/Users/ruchita/Documents/HARI/data/pop2019_world_bank.csv")
myFAO$Pop2019 = myPop2019$Pop2019[match(myFAO$ISO3,myPop2019$ISO3)]

# Import Income Group Information
############################################
myGNI2019 = read.csv("data/GNI_Capita_current_USD_2019.csv")
myGNI2019$GNI_2019 = ifelse(is.na(myGNI2019$GNI_2019),myGNI2019$GNI_2018,myGNI2019$GNI_2019)
myGNI2019$GNI_2019 = ifelse(is.na(myGNI2019$GNI_2019),myGNI2019$GNI_2017,myGNI2019$GNI_2019)

# Classification: # low income = 1; middle income = 2; high-income = 3; 
myT$GNI_2019 = myGNI2019$GNI_2019[match(myT$ISO3,myGNI2019$ISO3)]
myT$IncomeGroup = ifelse(myT$GNI_2019 > 1036,2,1) # low income = 1; low-middle income = 2; high-middle-income = 3; high-income = 4. 
myT$IncomeGroup = ifelse(myT$GNI_2019 > 12535,myT$IncomeGroup+1,myT$IncomeGroup)

myFAO$GNI_2019 = myGNI2019$GNI_2019[match(myFAO$ISO3,myGNI2019$ISO3)]
myFAO$IncomeGroup = ifelse(myFAO$GNI_2019 > 1036,2,1) 
myFAO$IncomeGroup = ifelse(myFAO$GNI_2019 > 12535,myFAO$IncomeGroup+1,myFAO$IncomeGroup)

# clean names of bacteria and define pathogens groups
######################################################
source("scripts_2023/clean_pathogens.R")

# Remove undefined pathogen group
myT = myT[which(is.na(myT$PathogensG)==F),]

# Remove undefined Drug group 
myT = myT[which(is.na(myT$Drug)==F),]

# Exclude Salmonella 
myT = subset(myT, PathogensG != "Salmonella spp")

# clean drug names 
########################################################
myT$Drug = gsub(" ","",myT$Drug)
myT$Drug = gsub("AMI/PEN","PEN/AMI",myT$Drug)
myT$Drug = gsub("PEM","PEN",myT$Drug)
myT$Drug = gsub("QUI/CAR","CAR/QUI",myT$Drug)
myT$Drug = gsub("AMP","PEN",myT$Drug)

# Remove undefined Drug group 
myT = myT[which(is.na(myT$Drug)==F),]

# clean resistance and Infections 
myT$HRes = as.numeric(myT$HRes)
myT$HNinf = as.numeric(myT$HNinf)

# Remove undefined Drug group 
myT = myT[which(is.na(myT$HRes)==F),]

# clean hospitalization rates 
myT$HRNTYR = as.numeric(myT$HRNTYR)

# Impute missing rates from covariates of the world Bank (physician per capita, Urban poulation, etc)
source("scripts_2023/Hopitalization_Rates_Imputation.r")
print(paste0("Countries with hospitalization Rates : ",length(na.exclude(myWB$HRNTYR))))
myT$HRNTYR_p = myWB$HRNTYR_p[match(myT$ISO3,myWB$ISO3)]
myT$HRNTYR = ifelse(is.na(myT$HRNTYR)==T,myT$HRNTYR_p,myT$HRNTYR)

# Add population  
myT$Pop2019 = myFAO$Pop2019[match(myT$ISO3,myFAO$ISO3)]

# number of publications 
print("number of publications")
print(length(unique(myT$DOI)))

# number of countries with publications
print("number of countries with one publication")
length(unique(myT$ISO3))

# rates by pathogens 
myTG = aggregate(count~myT$PathogensG, data = myT,FUN = sum)
myTG = myTG[rev(order(myTG$count)),]
print(myTG)

# total rates 
print(sum(myTG$count))

# Import Income Group Information
############################################
myGNI2019 = read.csv("data/GNI_Capita_current_USD_2019.csv")
myGNI2019$GNI_2019 = ifelse(is.na(myGNI2019$GNI_2019),myGNI2019$GNI_2018,myGNI2019$GNI_2019)
myGNI2019$GNI_2019 = ifelse(is.na(myGNI2019$GNI_2019),myGNI2019$GNI_2017,myGNI2019$GNI_2019)

# Classification: # low income = 1; middle income = 2; high-income = 3; 
myT$GNI_2019 = myGNI2019$GNI_2019[match(myT$ISO3,myGNI2019$ISO3)]
myT$IncomeGroup = ifelse(myT$GNI_2019 > 1036,2,1) # low income = 1; low-middle income = 2; high-income = 3; 
myT$IncomeGroup = ifelse(myT$GNI_2019 > 12535,myT$IncomeGroup+1,myT$IncomeGroup)

myFAO$GNI_2019 = myGNI2019$GNI_2019[match(myFAO$ISO3,myGNI2019$ISO3)]
myFAO$IncomeGroup = ifelse(myFAO$GNI_2019 > 1036,2,1) 
myFAO$IncomeGroup = ifelse(myFAO$GNI_2019 > 12535,myFAO$IncomeGroup+1,myFAO$IncomeGroup)

# Estimate hospitalization rates per country 
######################################################################
# First aggregate hospitalization rates by country 
myTHR = aggregate(cbind(Pop2019,HRNTYR)~ISO3, FUN = mean, na.rm = T,data = myT)
myTHR = myTHR[rev(order(myTHR$HRNTYR)),]
myTHR$IncomeGroup = myT$IncomeGroup[match(myTHR$ISO3,myT$ISO3)]

# Only select 'big enough' countries (5,000,000 people)
myTHR = subset(myTHR, Pop2019 > 5e6) 

# Colors for income categories 
ColsCat = c("khaki","orange1","navyblue")
par(mfrow = c(1,1))
par(mar = c(5,5,5,5))
barplot(myTHR$HRNTYR, names = myTHR$ISO3, space = 0, col = ColsCat[myTHR$IncomeGroup], border = F,
        ylab = "Hospitalized Pop / Year (%)", width = 1, las = 2)

# Impute hospitalization rates by income group 
##############################################
myTHR_IC = aggregate(HRNTYR~IncomeGroup, data = myT, FUN = median)
myFAO$HRNTYR_IC = myTHR_IC$HRNTYR[match(myFAO$IncomeGroup,myTHR_IC$IncomeGroup)]
myFAO$HRNTYR = ifelse(is.na(myFAO$HRNTYR),myFAO$HRNTYR_IC,myFAO$HRNTYR)

# Estimate the number of drug resistance infections in MC analysis
#######################################################################
# Parameters distributions
# Define a dataframe with parameters distribution used for the sensitivity analysis (Latin Hyper Cube) 
# Order of parameters: 
# LOS; length of stay in hospital
# Occupancy Low Income Countries, 
# Occupancy Middle Income countries 
# Occupancy High Income countries 
# Occupancy percentile of high values to be excluded 
ParRange <- data.frame(
 min = c(2, 1, 0.8, 0.6,0.98), 
 max = c(7, 1.5, 1, 1,1))
rownames(ParRange) <- c("LOS", "OLIC","OMIC","OHIC","H_PCT")

# Create a drug-bug combination Index 
myT$DBComb = paste(myT$PathogensG,myT$Drug, sep = "_")
myT$DBCombCTY = paste(myT$PathogensG,myT$Drug,myT$ISO3, sep = "_")

#Identify countries without minimum 2 surveys 
DFSbyCountry = as.data.frame(cbind(unique(myT$ISO3),rep(NA,length(unique(myT$ISO3)))))
names(DFSbyCountry) = c("ISO3","SurveybyCountry")

# This is to identify the country that have less than 2 surveys for a given drug-bug combination. 
for (ii in 1:nrow(DFSbyCountry)){
 DFSbyCountry$SurveybyCountry[ii] = length(unique(subset(myT, ISO3 == DFSbyCountry$ISO3[ii])$DOI))
} 
DFSbyCountry2Surveys = subset(DFSbyCountry, SurveybyCountry >=2)

# Number of countries with at least two surveys 
print(paste0("Number countries with 2 surveys = ",nrow(DFSbyCountry2Surveys)))

NoSurveyID = which(is.na(match(myFAO$ISO3,DFSbyCountry2Surveys$ISO3))==T)

# Number of drug-bug comninations
nDBC = length(unique(myT$DBComb))

set.seed(999)
# Define number of BS
nBS = 1000

# define an array to store output 
RinF_Y =  array(0, dim=c(nrow(myFAO),nDBC,nBS))



for(bs in 1:nBS){
 
 # Sampling of parameters for each bootstrap
 pars = as.data.frame(Latinhyper(ParRange, 1))
 
 # Make assumptions about occupancy rates
 myT$OccRate[myT$IncomeGroup == 1] = pars$OLIC
 myT$OccRate[myT$IncomeGroup == 2] = pars$OMIC
 myT$OccRate[myT$IncomeGroup == 2] = pars$OMIC
 myT$OccRate[myT$IncomeGroup == 3] = pars$OHIC
 myT$OccRate[is.na(myT$OccRate)] = 1
 
 # Estimate length of stay in hospital 
 myT$LOS = rep(pars$LOS,nrow(myT)) 
 
 # Calculate number of hospital visits over the study period 
 myT$Visits =  ((myT$Bed * myT$OccRate) / myT$LOS) * myT$StudyPeriod
 
 # number of drug-resistant infections per drug/bug combination per visit at hospital 
 myT$RInfHospV = myT$HNinf * myT$HRes / myT$Visits 
 
 # Adjust this for the length of the study period 
 myT$StudyPeriodYear = as.numeric(myT$StudyPeriod) / 365
 myT$RInfperHosp_Y = myT$RInfHospV / myT$StudyPeriodYear
 
 for (i in 1:length(unique(myT$DBComb))){
  
  # Subset on a specific drug-bug combination 
  myST = subset(myT, DBComb == unique(myT$DBComb[i]))
  
  # Remove the potential outliers (extreme high values) according to percentile subject to sensitivity analysis (H_PCT)
  myST = myST[myST$RInfperHosp_Y < quantile(na.exclude(myST$RInfperHosp_Y), prob = pars$H_PCT),]
  
  # Count the number of surveys for THIS drug-bug combination per country
  DFSbyCountry = as.data.frame(cbind(unique(myST$ISO3),rep(NA,length(unique(myST$ISO3)))))
  names(DFSbyCountry) = c("ISO3","SurveybyCountry")
  
  # Iidentify countries that have less than 2 surveys for a given drug-bug combination. 
  for (ii in 1:nrow(DFSbyCountry)){
   DFSbyCountry$SurveybyCountry[ii] = length(unique(subset(myST, ISO3 == DFSbyCountry$ISO3[ii])$DOI))
  } 
  DFSbyCountryEnoughSurveys = subset(DFSbyCountry,SurveybyCountry >= 2) # This is arbitrary
  NotEnoughSurveys = which(is.na(match(myFAO$ISO3,DFSbyCountryEnoughSurveys$ISO3))==T)
  
  # Build Table of prevalence of drug-resistant infection per hospitalization by country per year
  ####################################################################################################
  myAT = aggregate(RInfperHosp_Y ~ ISO3, data = myST, FUN = mean, na.rm = T)
  myFAO$RInfperHosp_Y_Run = myAT$RInfperHosp_Y[match(myFAO$ISO3,myAT$ISO3)]
  
  # Remove national average for countries that have less than 2 surveys 
  myFAO$RInfperHosp_Y_Run[NotEnoughSurveys] = NA
  
  #  Calculate average rate per income group
  myAT$IncomeGroup = myFAO$IncomeGroup[match(myAT$ISO3,myFAO$ISO3)]
  myAAT = aggregate(RInfperHosp_Y ~ IncomeGroup, data = myAT, FUN = mean, na.rm = T)
  
  # Calculate resistance rates by income groups 
  myFAO$RInfperHosp_Y_IC_Run = myAAT$RInfperHosp_Y[match(myFAO$IncomeGroup,myAAT$IncomeGroup)]
  
  # Impute missing value by income groups resistant infection rates 
  myFAO$RInfperHosp_Y_Run = ifelse(is.na(myFAO$RInfperHosp_Y_Run)==T,myFAO$RInfperHosp_Y_IC_Run,myFAO$RInfperHosp_Y_Run)
  
  # Calculate  number of infections per country. Each given drug-bug combination is a column of the matrix RInfDF
  if(ImpTAG == F){myFAO$RInfperHosp_Y_Run[NoSurveyID]=NA}
  
  # Fill the output array with result     
  RinF_Y[,i,bs] = myFAO$Pop2019 * myFAO$HRNTYR * myFAO$RInfperHosp_Y_Run
  
 } # enf of loop on DBC
 
} # end of loop on bootstraps

# Calculate mean number of infection, and sd by country 
#######################################################
RinF_Y_M = apply(RinF_Y, MARGIN = c(1,3), FUN = sum, na.rm = T) 

DRI_World = apply(RinF_Y_M,MARGIN = 1,FUN = mean, na.rm = T) / 1e6 # million of drug-resistant infections
DRI_World_CIup = DRI_World + 1.96 * apply(RinF_Y_M,MARGIN = 1,FUN = sd) /  1e6
DRI_World_CIlw = DRI_World - 1.96 * apply(RinF_Y_M,MARGIN = 1,FUN = sd) /  1e6
myFAO$RinF_Y_all = DRI_World
myFAO$RinF_Y_all_up = DRI_World_CIup
myFAO$RinF_Y_all_lw = DRI_World_CIlw

# Calculate number of infections by country per DBC
###################################################
# Mean 
Inf_by_ISO_DBC = apply(RinF_Y, MARGIN = c(1,2), FUN = mean, na.rm = T)
Inf_by_ISO_DBC_df = as.data.frame(Inf_by_ISO_DBC)
names(Inf_by_ISO_DBC_df) = unique(myT$DBComb)

# Add country names
ISO3 = myFAO$ISO3
Inf_by_ISO_DBC_df = cbind(ISO3,Inf_by_ISO_DBC_df)

#write.csv(Inf_by_ISO_DBC_df, "/Users/thomavan/Dropbox/HARI/ms/Supp_Table_2.csv")

#Total number of infections in the world (millions)
##################################################
MeanEst = format(sum(DRI_World, na.rm = T), digits = 4)
LwEst = format(sum(DRI_World_CIlw, na.rm = T), digits = 4)
UpEst = format(sum(DRI_World_CIup, na.rm = T), digits = 4)
print(paste0("Global DR Infections = ",MeanEst," [",LwEst,"-",UpEst,"]"))

# Calculate mean number of infection, and sd per capita by country per capita
#############################################################################
DRI_World_cap = (DRI_World * 1e6) / (1e-3 * myFAO$Pop2019)# million of drug-resistant infections per 1000 poeple 
DRI_World_CIup_cap = (DRI_World_CIup * 1e6) / (1e-3 * myFAO$Pop2019)
DRI_World_CIlw_cap = (DRI_World_CIlw * 1e6) / (1e-3 * myFAO$Pop2019)

# output HARI per country with CI
Table_HARI_CNY_CI = as.data.frame(cbind(myFAO$ISO3,DRI_World * 1e6,DRI_World_CIlw * 1e6,DRI_World_CIup * 1e6))
names(Table_HARI_CNY_CI) = c("ISO3","HARI","HARI_95_lw","HARI_95_up")
Table_HARI_CNY_CI$HARI = Table_HARI_CNY_CI$HARI = signif(as.numeric(Table_HARI_CNY_CI$HARI), digits = 3)
Table_HARI_CNY_CI$HARI_95_lw = Table_HARI_CNY_CI$HARI_95_lw = signif(as.numeric(Table_HARI_CNY_CI$HARI_95_lw), digits = 3)
Table_HARI_CNY_CI$HARI_95_up = Table_HARI_CNY_CI$HARI_95_up = signif(as.numeric(Table_HARI_CNY_CI$HARI_95_up), digits = 3)
Table_HARI_CNY_CI$HARI = as.character(Table_HARI_CNY_CI$HARI)
Table_HARI_CNY_CI$HARI_95_lw = as.character(Table_HARI_CNY_CI$HARI_95_lw)
Table_HARI_CNY_CI$HARI_95_up = as.character(Table_HARI_CNY_CI$HARI_95_up)
Table_HARI_CNY_CI$HARI_95_lw[Table_HARI_CNY_CI$HARI_95_lw<0] = 0 
#write.table(Table_HARI_CNY_CI, "/Users/thomavan/Dropbox/HARI/ms/Supp_Table_3.csv")
#Table_HARI_CNY_CI_Final = Table_HARI_CNY_CI
# Add value to my FAO dataframe 
myFAO$RinF_Y_all_Cap = DRI_World_cap
myFAO$RinF_Y_all_Cap_lw = DRI_World_CIlw_cap
myFAO$RinF_Y_all_Cap_up = DRI_World_CIup_cap


################################################################################
# Crude map of total infections
################################################################################
# This now need to be multiplied by the hospitalization rate
mySHP@data$res = myFAO$RinF_Y_all[match(mySHP$ISO3,myFAO$ISO3)] * 1e6
mySHP@data$resPlot = sqrt(mySHP@data$res) / 400

X = coordinates(mySHP)[,1]
Y = coordinates(mySHP)[,2]

png(width = 1600, height = 1200, filename = "current_plot.png")
par(mfrow = c(1,1))
plot((mySHP), lwd = 0.5, border = "black")

points(X,Y, col = "red3", pch = 16, cex = mySHP@data$resPlot)

dev.off()

# Show the top 10 Countries 
top10 = myFAO[rev(order(myFAO$RinF_Y_all)),][1:20,]
top10$DRI = top10$RinF_Y_all / 1e6
top10[c("ISO3","RinF_Y_all","RinF_Y_all_lw","RinF_Y_all_up")]

# Show total DRI by income group
myFAOIncome = aggregate(myFAO[c("RinF_Y_all","RinF_Y_all_lw","RinF_Y_all_up")], by = list(myFAO$IncomeGroup), FUN = sum)
names(myFAOIncome) = c("IncomeGroup","RinF_Y_all","RinF_Y_all_lw","RinF_Y_all_up")
myFAOIncome[myFAOIncome<0] = 0

# output as table 
CRTable = mySHP@data[c("ISO3","res")] 
CRTable = CRTable[rev(order(CRTable$res)),]
CRTable$res = round(CRTable$res)
CRTable$res = as.character(format(CRTable$res, big.mark = ","))
names(CRTable) = c("Country ISO3","Infections")
#write.table(CRTable, "/Users/thomavan/Dropbox/HARI/ms/Supp_Table_1.txt", sep = "\t", row.names = F, quote = F)

# Present similar information as barplot pooled by region
#########################################################
RinF_Y_DBC_M = apply(RinF_Y, MARGIN = c(1,2), FUN = mean, na.rm = T) 
par(mar = c(11,5,5,5))
par(mfrow = c(1,3))
RInfAVG_df = as.data.frame(RinF_Y_DBC_M)
names(RInfAVG_df) = unique(myT$DBComb)
RInfAVG_dfA = aggregate(RInfAVG_df, by = list(myFAO$IncomeGroup), FUN = sum, na.rm = T)
RInfAVG_dfAM = as.matrix(RInfAVG_dfA)

# Clean resulting matrix 
NamesCatclean = substr(dimnames(RInfAVG_dfAM)[[2]],1,3)

# Remove first column with income group 
RInfAVG_dfAM = RInfAVG_dfAM[,-1]

# Define colors for each DBC valid for each income group by pathogens
colDF = as.data.frame(cbind(dimnames(RInfAVG_dfAM)[[2]],substr(dimnames(RInfAVG_dfAM)[[2]],1,3)))
names(colDF) = c("DBC","ColsCat")

colDF = colDF[order(colDF$ColsCat),]

# Define color vector
ColsCat = brewer.pal(length(unique(colDF$ColsCat)),"Set1")

Cols = NULL

for (i in 1:length(colDF$ColsCat)){
 coltmp = rep(ColsCat[i],length(which(colDF$ColsCat == unique(colDF$ColsCat)[i])))
 if(length(coltmp) > 1){
  coltmp = addalpha(coltmp,seq(0.5,1,0.5/(length(coltmp)-1)))
 } 
 Cols = c(Cols,coltmp) 
} # categories 

colDF = cbind(colDF,Cols)
names(colDF) = c("DBC","ColsCat","Cols")

otherDF = as.data.frame(cbind("others","others","grey"))
names(otherDF) =  c("DBC","ColsCat","Cols")

colDF = rbind(colDF,otherDF)

colDF$PathogenG = gsub( "_.*$", "", colDF$DBC)

# Sanity check on colors for each DBC
par(mfrow = c(1,1))
barplot(rep(10,nrow(colDF)), names.arg = colDF$DBC, las = 2, col = colDF$Cols, space = 0, border = F)

# Define figure meta parameters: Ymax and percentage for pooling in 'other' class
#########################################
Pct = 0.005
Ymax = 5e6

# Low-Income Countries
######################
# Barplot by drug-bug combination
par(mar = c(4,11,2,1))
par(mfrow = c(3,2))

DCBs_ICv = RInfAVG_dfAM[1,]
Othersv = DCBs_ICv[DCBs_ICv < (Pct * sum(DCBs_ICv))]
DCBs_ICvBig = DCBs_ICv[DCBs_ICv >= (Pct * sum(DCBs_ICv))]
NamesBig = dimnames(RInfAVG_dfAM)[[2]][DCBs_ICv >= Pct * sum(DCBs_ICv)]

DCBs_ICvMerge = c(DCBs_ICvBig,sum(Othersv))
NamesMerge = c(NamesBig,"others")
ColsMerge = colDF$Cols[match(NamesMerge,colDF$DBC)]

DCB_IC_df = as.data.frame(cbind(NamesMerge,DCBs_ICvMerge,ColsMerge))
DCB_IC_df$DCBs_ICvMerge = as.numeric(DCB_IC_df$DCBs_ICvMerge)
DCB_IC_df = DCB_IC_df[order(DCB_IC_df$NamesMerge),]

barplot(DCB_IC_df$DCBs_ICvMerge, names.arg = DCB_IC_df$NamesMerge, las = 2, space = 0, cex.names = 1,
        xlim = c(0,Ymax), col = DCB_IC_df$ColsMerge, main = "Low-Income Countries", cex.main = 2, horiz=T, border = F)
abline(v = seq(1e6,5e6,1e6), col = "grey", lty = 3)

# Pie by Pathogens
par(mar = c(1,1,1,1))
DCB_IC_df$PathogenG = gsub( "_.*$", "", DCB_IC_df$NamesMerge)
DCB_IC_dfG = aggregate(DCBs_ICvMerge~PathogenG,data = DCB_IC_df, FUN = sum,na.rm = T)
DCB_IC_dfG$ColsGMerge = colDF$Cols[match(DCB_IC_dfG$PathogenG,colDF$PathogenG)]

TotalInf = format(round(1000*round(sum(DCBs_ICv)/1000)), big.mark = ",")
pie(DCB_IC_dfG$DCBs_ICvMerge, col = DCB_IC_dfG$ColsGMerge, border = F, labels = DCB_IC_dfG$PathogenG, radius = 1.4, main = paste0("\n",TotalInf), cex.main =2)

# Middle-Income Countries
######################
# Barplot by drug-bug combination
par(mar = c(4,11,2,1))

DCBs_ICv = RInfAVG_dfAM[2,]
Othersv = DCBs_ICv[DCBs_ICv < (Pct * sum(DCBs_ICv))]
DCBs_ICvBig = DCBs_ICv[DCBs_ICv >= (Pct * sum(DCBs_ICv))]
NamesBig = dimnames(RInfAVG_dfAM)[[2]][DCBs_ICv >= (Pct * sum(DCBs_ICv))]

DCBs_ICvMerge = c(DCBs_ICvBig,sum(Othersv))
NamesMerge = c(NamesBig,"others")
ColsMerge = colDF$Cols[match(NamesMerge,colDF$DBC)]

DCB_IC_df = as.data.frame(cbind(NamesMerge,DCBs_ICvMerge,ColsMerge))
DCB_IC_df$DCBs_ICvMerge = as.numeric(DCB_IC_df$DCBs_ICvMerge)
DCB_IC_df = DCB_IC_df[order(DCB_IC_df$NamesMerge),]

barplot(DCB_IC_df$DCBs_ICvMerge, names.arg = DCB_IC_df$NamesMerge, las = 2, space = 0, cex.names = 0.6,
        xlim = c(0,Ymax), col = DCB_IC_df$ColsMerge, axes = T, main = "Middle-Income Countries", cex.main = 2, horiz = T, border = F)
abline(v = seq(1e6,5e6,1e6), col = "grey", lty = 3)

# Pie by Pathogens
par(mar = c(1,1,1,1))
DCB_IC_df$PathogenG = gsub( "_.*$", "", DCB_IC_df$NamesMerge)
DCB_IC_dfG = aggregate(DCBs_ICvMerge~PathogenG,data = DCB_IC_df, FUN = sum,na.rm = T)
#DCB_IC_dfG$ColsGMerge = DCB_IC_df$ColsMerge[match(DCB_IC_dfG$PathogenG,DCB_IC_df$PathogenG)]
DCB_IC_dfG$ColsGMerge = colDF$Cols[match(DCB_IC_dfG$PathogenG,colDF$PathogenG)]

TotalInf = format(round(1000*round(sum(DCBs_ICv)/1000)), big.mark = ",")
pie(DCB_IC_dfG$DCBs_ICvMerge, col = DCB_IC_dfG$ColsGMerge, border = F, labels = DCB_IC_dfG$PathogenG, radius = 1.4, main = paste0("\n",TotalInf), cex.main =2)


# High-Income Countries
######################
# Barplot by drug-bug combination
par(mar = c(4,11,2,1))

DCBs_ICv = RInfAVG_dfAM[3,]
Othersv = DCBs_ICv[DCBs_ICv < (Pct * sum(DCBs_ICv))]
DCBs_ICvBig = DCBs_ICv[DCBs_ICv >= (Pct * sum(DCBs_ICv))]
NamesBig = dimnames(RInfAVG_dfAM)[[2]][DCBs_ICv >= (Pct * sum(DCBs_ICv))]

DCBs_ICvMerge = c(DCBs_ICvBig,sum(Othersv))
NamesMerge = c(NamesBig,"others")
ColsMerge = colDF$Cols[match(NamesMerge,colDF$DBC)]

DCB_IC_df = as.data.frame(cbind(NamesMerge,DCBs_ICvMerge,ColsMerge))
DCB_IC_df$DCBs_ICvMerge = as.numeric(DCB_IC_df$DCBs_ICvMerge)
DCB_IC_df = DCB_IC_df[order(DCB_IC_df$NamesMerge),]

barplot(DCB_IC_df$DCBs_ICvMerge, names.arg = DCB_IC_df$NamesMerge, las = 2, space = 0, cex.names = 0.6,
        xlim = c(0,Ymax), col = DCB_IC_df$ColsMerge, axes = T,main = "High-Income Countries", cex.main = 2, horiz = T, border = F)
abline(v = seq(1e6,5e6,1e6), col = "grey", lty = 3)

# Pie by Pathogens
par(mar = c(1,1,1,1))
DCB_IC_df$PathogenG = gsub( "_.*$", "", DCB_IC_df$NamesMerge)
DCB_IC_dfG = aggregate(DCBs_ICvMerge~PathogenG,data = DCB_IC_df, FUN = sum,na.rm = T)
#DCB_IC_dfG$ColsGMerge = DCB_IC_df$ColsMerge[match(DCB_IC_dfG$PathogenG,DCB_IC_df$PathogenG)]
DCB_IC_dfG$ColsGMerge = colDF$Cols[match(DCB_IC_dfG$PathogenG,colDF$PathogenG)]
TotalInf = format(round(1000*round(sum(DCBs_ICv)/1000)), big.mark = ",")
pie(DCB_IC_dfG$DCBs_ICvMerge, col = DCB_IC_dfG$ColsGMerge, border = F, labels = DCB_IC_dfG$PathogenG, radius = 1.4, main = paste0("\n",TotalInf), cex.main =2)


# Calculate total by income group with CI
DRI_Income_agg = aggregate(myFAO[c("RinF_Y_all","RinF_Y_all_lw","RinF_Y_all_up")], by = list(myFAO$IncomeGroup), FUN = sum, na.rm = T)


# Number of low-income countries with data 
myTL = subset(myT, IncomeGroup == 1)

print(paste0("Low-Income ountries with data n = ",length(unique(myTL$ISO3))))


setEPS()
#postscript("/Users/ruchita/Documents/HARI/figures/fig3_eps.eps",width = 12, height = 5)
par(mfrow = c(1,1))
par(mar = c(4,4,4,4))
plot((mySHP), lwd = 0.5, border = "black")
#rbPal <- colorRampPalette(c('red','blue'))
#mySHP@data$Col[which(mySHP@data$resPlot>0)] <- rbPal(10)[as.numeric(cut(mySHP@data$resPlot[which(mySHP@data$resPlot>0)],breaks = 10))]
points(X,Y, col = "red3", pch = 16, cex = 0.25 * mySHP@data$resPlot)
#mySHP@data$ResCol<- ifelse(mySHP@data$resPlot>0,1,0)
#points(X,Y, col = mySHP@data$Col, pch = 16, cex = 1.2*mySHP@data$ResCol)

XXX = (15*400)^2

#image(1, mySHP@data$resPlot, t(seq_along(mySHP@data$resPlot)), col = mySHP@data$Col, axes=FALSE)


text(-160,30,"Infections/Year", font = 2)
points(-160,0,cex = 4.9, lwd = 3)
text(-110,0,"30,000,000", cex = 1)
points(-160,-38,cex = 2.9, lwd = 3)
text(-110,-38,"10,000,000", cex = 1)
points(-160,-58,cex = .9, lwd = 3)
text(-110,-58,"1,000,000", cex = 1)

dev.off()





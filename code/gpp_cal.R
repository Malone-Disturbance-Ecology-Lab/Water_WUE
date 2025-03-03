rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)

#packageVersion("REddyProc") 1.3.3

#https://rdrr.io/cran/flux/man/gpp.html
#https://rdrr.io/rforge/REddyProc/man/sEddyProc.example.html

# florida sites
data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/florida/gaps_US-Skr.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/florida/gaps_US-EvM.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/florida/gaps_US-TaS.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/florida/gaps_US-KS3.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/florida/gaps_US-KS4.csv",header=TRUE)

# california sites
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-Dmg.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-EDN.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-EKH.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-EKP.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-Myb.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-Srr.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/gaps_US-Tw1.csv",header=TRUE)

# massachusetts

#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/massachusetts/gaps_US-PHM.csv",header=TRUE)

# south carolina

#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/south carolina/gaps_US-HB1.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/south carolina/gaps_US-HB2.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/south carolina/gaps_US-HB3.csv",header=TRUE)

# New Jersey

#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/new jersey/gaps_US-Hpy.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/new jersey/gaps_US-Mrm.csv",header=TRUE)

## lousiana 

#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/lousiana/gaps_US-LA1.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/lousiana/gaps_US-LA2.csv",header=TRUE)
#data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/lousiana/gaps_US-LA3.csv",header=TRUE)


oct<-(data)

oct=oct[,-1]

EddyDataWithPosix.F <- fConvertTimeToPosix(oct, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')

#EddyProc.C<- sEddyProc$new('oct', EddyDataWithPosix.F, c('NEE','LE','H', 'Rg','Tair','VPD','Ustar','PA','NETRAD'),DTS=48)
EddyProc.C<- sEddyProc$new('oct', EddyDataWithPosix.F, c('NEE','LE','H', 'Rg','Tair','VPD','Ustar','PA','NETRAD'),DTS=48)
#EddyProc.C<- sEddyProc$new('oct', EddyDataWithPosix.F, c('NEE','LE','H', 'Rg','Tair','VPD','Ustar','PA'),DTS=48)

#EddyProc.C<- sEddyProc$new('oct', EddyDataWithPosix.F, c('NEE','LE','H', 'Rg','Tair','VPD','Ustar'),DTS=48)

uStarTh <- EddyProc.C$sEstUstarThreshold()$uStarTh

EddyProc.C$sMDSGapFillAfterUstar('NEE') # it uses threshold 
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFillAfterUstar('LE')
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFill('Rg')
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE) # Only fills gaps where sufficient meteorological variables are available, dont fill all gaps
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFill('VPD', FillAll.b=FALSE)
EddyProc.C$sExportResults()

filled1<- EddyProc.C$sExportResults()
##########################################################################################

### GPP part
################
##Florida

EddyProc.C$sSetLocationInfo(Lat_deg.n=25.3629, Long_deg.n=-81.0776, TimeZone_h.n=1)#US-Skr 
#EddyProc.C$sSetLocationInfo(Lat_deg.n=25.3539, Long_deg.n=-80.3810, TimeZone_h.n=1)#US-EvM
#EddyProc.C$sSetLocationInfo(Lat_deg.n=25.1908, Long_deg.n=-80.6391, TimeZone_h.n=1)#US-tas
#EddyProc.C$sSetLocationInfo(Lat_deg.n=28.7084, Long_deg.n=-80.7427, TimeZone_h.n=1)#US-ks3
#EddyProc.C$sSetLocationInfo(Lat_deg.n=28.6000, Long_deg.n=-80.7207, TimeZone_h.n=1)#US-ks4

##################

## california
#EddyProc.C$sSetLocationInfo(Lat_deg.n=38.0015, Long_deg.n=-121.6691, TimeZone_h.n=1)#US-dmg
#EddyProc.C$sSetLocationInfo(Lat_deg.n=37.6156, Long_deg.n=-122.1140, TimeZone_h.n=1)#US-EDN
#EddyProc.C$sSetLocationInfo(Lat_deg.n=36.8094, Long_deg.n=-121.7523, TimeZone_h.n=1)#US-EKH
#EddyProc.C$sSetLocationInfo(Lat_deg.n=36.8558, Long_deg.n=-121.7488, TimeZone_h.n=1)#US-EKP
#EddyProc.C$sSetLocationInfo(Lat_deg.n=38.2006, Long_deg.n=-122.0264, TimeZone_h.n=1)#US-Srr
#EddyProc.C$sSetLocationInfo(Lat_deg.n=38.1074, Long_deg.n=-121.6469, TimeZone_h.n=1)#US-Tw1

#########################################################################################

##Massachusetts 
#EddyProc.C$sSetLocationInfo(Lat_deg.n=42.7423, Long_deg.n=-70.8301, TimeZone_h.n=1)#US-PhM
###############################################################################################
# south carolina

#EddyProc.C$sSetLocationInfo(Lat_deg.n=33.3455, Long_deg.n=-79.1957, TimeZone_h.n=1)#US-HB1
#EddyProc.C$sSetLocationInfo(Lat_deg.n=33.3242, Long_deg.n=--79.2440, TimeZone_h.n=1)#US-HB2
#EddyProc.C$sSetLocationInfo(Lat_deg.n=33.3482, Long_deg.n=--79.2322, TimeZone_h.n=1)#US-HB3

##############################################################################################

# New Jersey 

#EddyProc.C$sSetLocationInfo(Lat_deg.n=40.7692, Long_deg.n=-74.0853, TimeZone_h.n=1)#US-Hpy
#EddyProc.C$sSetLocationInfo(Lat_deg.n=40.8164, Long_deg.n=-74.0435, TimeZone_h.n=1)#US-Mrm

###############################################################################################
# lousiana 

#EddyProc.C$sSetLocationInfo(Lat_deg.n=29.5013, Long_deg.n=-90.4449, TimeZone_h.n=1)#US-LA1
#EddyProc.C$sSetLocationInfo(Lat_deg.n=29.8587, Long_deg.n=-90.2869, TimeZone_h.n=1)#US-LA2
#EddyProc.C$sSetLocationInfo(Lat_deg.n=29.4936, Long_deg.n= -89.9153, TimeZone_h.n=1)#US-LA3

################################################################################################
#EddyProc.C$sMDSGapFillAfterUstar('NEE')
#EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE)
#EddyProc.C$sMDSGapFill('VPD', FillAll.b=FALSE)

EddyProc.C$sMRFluxPartition(Suffix.s='uStar')  # night time #Reichstein 2005
EddyProc.C$sGLFluxPartition(Suffix.s='uStar') # day time  Reco_DT, GPP_DT #Lasslop, light response curve
#EddyProc.C$sTKFluxPartition(Suffix.s='uStar')# modified day 


filled1<- EddyProc.C$sExportResults()

EddyDataWithPosix.F$DateTime <- strftime(EddyDataWithPosix.F$DateTime, "%Y-%m-%d %H:%M:%S")

filled2=data.frame(cbind(DateTime=EddyDataWithPosix.F$DateTime, Year=EddyDataWithPosix.F$Year,
DoY=EddyDataWithPosix.F$DoY,Hour=EddyDataWithPosix.F$Hour,NEE_f=filled1$NEE_uStar_f,LE_f=filled1$LE_uStar_f,
Tair_f=filled1$Tair_f,VPD_f=filled1$VPD_f,Rg_f=filled1$Rg_f,GPP_DT=filled1$GPP_DT_uStar,
GPP_nt=filled1$GPP_uStar_f,Reco_DT=filled1$Reco_DT_uStar,Reco_nt=filled1$Reco_uStar))


# florida
write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/florida/fill_US-Skr.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/fill_US-EvM.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/fill_US-TaS.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/fill_US-KS3.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/florida/fill_US-KS4.csv")

# california
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-Dmg.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-EDN.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-EKH.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-EKP.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-Srr.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/fill_US-Tw1.csv")

## massachusetts
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/massachusetts/fill_US-PHM.csv")

## south carolina
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/south carolina/fill_US-HB1.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/south carolina/fill_US-HB2.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/south carolina/fill_US-HB3.csv")

# New Jersey
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/new jersey/fill_US-Hpy.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/new jersey/fill_US-Mrm.csv")

# luoisina 

#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/lousiana/fill_US-LA1.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/lousiana/fill_US-LA2.csv")
#write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/lousiana/fill_US-LA3.csv")
#########################################################################################################

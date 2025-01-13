rm(list=ls(all=TRUE)) 
library(Rcpp)
library(REddyProc)

#packageVersion("REddyProc") 1.3.3

#https://rdrr.io/cran/flux/man/gpp.html
#https://rdrr.io/rforge/REddyProc/man/sEddyProc.example.html

data <-read.csv(file="C:/ammara_MD/a_yale/Research Project/data/california/fill_US-EKH_.csv",header=TRUE)

df<-(data)

df2=df[,-1]

## get proper datetime format
EddyDataWithPosix.F <- fConvertTimeToPosix(df2, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')

# include variables that you want to gap fill along with meteorological variables
EddyProc.C<- sEddyProc$new('df2', EddyDataWithPosix.F, c('NEE','LE','H', 'Rg','Tair','VPD','Ustar'),DTS=48)


#############################################################################################333
#ustar filtering step
uStarTh <- EddyProc.C$sEstUstarThreshold()$uStarTh

## gap fill individual variables
EddyProc.C$sMDSGapFillAfterUstar('NEE') # 
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFillAfterUstar('LE')
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFill('Rg')
EddyProc.C$sExportResults()

#### need to gap fill Tair and VPD to use for GPP calculation

EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE)
EddyProc.C$sExportResults()

EddyProc.C$sMDSGapFill('VPD', FillAll.b=FALSE)
EddyProc.C$sExportResults()

filled1<- EddyProc.C$sExportResults()
##########################################################################################

### GPP part

EddyProc.C$sSetLocationInfo(Lat_deg.n=36.8094, Long_deg.n=-121.7523, TimeZone_h.n=1)#US-EKH

EddyProc.C$sMDSGapFillAfterUstar('NEE')
EddyProc.C$sMDSGapFill('Tair', FillAll.b=FALSE)
EddyProc.C$sMDSGapFill('VPD', FillAll.b=FALSE)

EddyProc.C$sMRFluxPartition(Suffix.s='uStar')  # night time #Reichstein 2005
EddyProc.C$sGLFluxPartition(Suffix.s='uStar') # day time  Reco_DT, GPP_DT #Lasslop, light response curve (preferred method where more sunlight)
################################################################################################################

## save gap filled variables in one dataframe
filled1<- EddyProc.C$sExportResults()

filled2=data.frame(cbind(DateTime=EddyDataWithPosix.F$DateTime, Year=EddyDataWithPosix.F$Year,
                         DoY=EddyDataWithPosix.F$DoY,Hour=EddyDataWithPosix.F$Hour,NEE=filled1$NEE_uStar_f,GPP_DT=filled1$GPP_DT_uStar,
                         GPP_nt=filled1$GPP_uStar_f,Reco_DT=filled1$Reco_DT_uStar,Reco_nt=filled1$Reco_uStar))



write.csv(filled2,"C:/ammara_MD/a_yale/Research Project/data/california/check.csv")

########################################################################################################3

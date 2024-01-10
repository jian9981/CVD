
setwd("P:/ORD_Bress_202112038D/Jian")

.libPaths(c(.libPaths(),"D:/R_Temp/Packages/3.4.0"))



library(RODBC)
library(ggplot2)
library(stringr)
library(dplyr)
library(gtsummary)
library(mice)
library(lubridate)
odbcDataSources()
conn_odbc <- odbcConnect("VHACDWRB02_17")
sqlTables(conn_odbc,catalog="ORD_Bress_202112038D",schema = "dflt") 
data0=sqlQuery(conn_odbc, "select  * from ORD_Bress_202112038D.dflt.ArbData_Selected_aa_3 where inclusion=1 and BL_CVD=0 and BL_Endstage_kidney_disease=0 and BL_Kidney_transplant=0 ")  #exclusion criteria on page 5 of the analysis plan
table(data0$BL_CVD )
table(data0$BL_Endstage_kidney_disease)
table(data0$BL_Kidney_transplant)
mean(is.na(data0$BL_PatientFIPS))
CAINC1=read.csv("./data/CAINC1__ALL_AREAS_1969_2018.csv")
CAINC1$GeoFIPS1=as.integer(CAINC1$GeoFIPS)
CAINC1=CAINC1%>%filter(LineCode==3)%>%dplyr::select(-one_of(c("GeoName","Region","TableName", "LineCode", "IndustryClassification","Description" , "Unit")))
CAINC1=CAINC1%>%tidyr::pivot_longer(-c(GeoFIPS,GeoFIPS1),values_to =  "area.median.income", names_to = "year")
CAINC1=CAINC1%>%mutate(year=as.numeric(gsub("X","",year)),
                       area.median.income=as.numeric(area.median.income),
                       GeoFIPS=str_trim(GeoFIPS,"both"))



data0 <- data0 %>%
  mutate(region=as.factor(region),
         BL_Priority = case_when(BL_Priority =="1" ~ "1(highest need)",
                                 is.na(BL_Priority)~NA,
                                 TRUE ~ "2 through 9"),
         BL_FrailtyIndexMaxValue_grp = case_when(BL_FrailtyIndexMaxValue<7 ~ "Non-frail(number of conditions<7",
                                                 BL_FrailtyIndexMaxValue>=7 ~ "Frail(number of conditions>=7"),
         race_ethnicity = case_when(ethnicity=="NOT HISPANIC OR LATINO" & race=="BLACK OR AFRICAN AMERICAN" ~ "Non-Hispanic Black",
                                    ethnicity=="NOT HISPANIC OR LATINO" & race=="WHITE" ~ "Non-Hispanic White",
                                    ethnicity=="HISPANIC OR LATINO"   ~ "Hispanic",
                                    TRUE~"Other"),
         OtherAntiHyper=BL_M_ARNI+BL_M_Central+BL_M_DirectVas+
           BL_M_DirectReInh+BL_M_AldoRecep+BL_M_PostDiu, 
         year=year(drugdateTime),
         GeoFIPS=BL_patientFIPS
  )

data0 <- data0 %>%left_join(CAINC1%>%dplyr::select(year, GeoFIPS, area.median.income))

bsl.covariates=c("ageAtIndex", "gender","race_ethnicity",
                 "area.median.income", "region","BL_PatientSupplementalIns","BL_Priority", "BL_Smoker","BL_Homeless",
                 "BL_BMI","BL_Systolic" , "BL_Diastolic","BL_HeartRate",
                 "BL_total_cholesterol_value", "BL_HDL_C_value" ,"BL_LDL_C_value","BL_triglyceride_value",
                 "BL_hemoglobin_A1c_value","BL_eGFR_value","BL_urineACR_value", "BL_serum_sodium_value","BL_serum_potassium_value" ,
                 "BL_LVEF_value","BL_Alcohol",
                 "BL_Arrhythmia","BL_Cancer","BL_Chronic_Liver_Disease","BL_Chronic_Lung_Disease",
                 "BL_Cirrhosis","BL_Diabetes","BL_Depression","BL_Drug_substance_abuse",
                 "BL_Obstructive_sleep_apnea",
                 "BL_FrailtyIndexMaxValue",
                 "BL_FrailtyIndexMaxValue_grp","BL_Syncope",
                 "BL_M_AspirinUsed","BL_M_StainUsed","BL_M_SGLT2Used","BL_M_CalcCB",
                 "BL_M_ThiaDiu","BL_M_LoopDiu","BL_M_BetaB", "BL_M_AlphaB",
                 "OtherAntiHyper",  "BL_M_NumOfAntiHyperMed",
                 "BL_PrimaryCareVisitCount","BL_HospitalizationsCount",
                 "BL_ERCareVisitCount")
data2 <- data0[,c("patientICN" ,"group_ACEI_ARB", bsl.covariates)]
data2=data2%>%mutate(BL_M_AlphaB=BL_M_AlphaB>0,
                     BL_M_BetaB=BL_M_BetaB>0,
                     BL_M_CalcCB=BL_M_CalcCB>0,
                     BL_M_LoopDiu=BL_M_LoopDiu>0,
                     BL_M_ThiaDiu=BL_M_ThiaDiu>0)
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###       MI, Weighting, and all analysis of effect - point estimate using original data (no resampling..)                 ##############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2=data2 %>% mutate_at(c("race_ethnicity", "BL_Priority"), as.factor) ## do we impute outcomes? - no
miss=findmiss(data2)
miss/nuni(data2$patientICN)
mvars=names(miss);length(mvars) 
glimpse(data2[,mvars])

outcomemiss=intersect(mvars, c("CVD.flp.year", "CVD")) 
outcomemiss

## MI now
ini=mice(data2, maxit=0)
(meth1=ini$meth)
meth1[mvars] 
meth1[c(outcomemiss, "patientICN" ,"group_ACEI_ARB") ]="" # do not impute these

puse=c( bsl.covariates) ## only use these
pdrop=setdiff(names(data2), puse)
pred1=quickpred(data2, mincor=0.01, exclude=pdrop)

dim(pred1)
pred1[mvars,] %>% rowSums ## total predictors per missing variable
colSums(pred1)[colSums(pred1)>0] ## all variables ever used as predictor

mvars1=setdiff(mvars, outcomemiss)
pred1[mvars1, ] %>% rowSums %>% sort
meth1[meth1!=""] 

t1=Sys.time()
imputed.data=mice(data2, m=1, predictorMatrix=pred1, meth=meth1, seed=1977, maxit=10)
t2=Sys.time()
t2-t1
findmiss(complete(imputed.data,2))
findmiss(data0)
table(data2$group_ACEI_ARB)
data2=data2%>%mutate(group_ACEI_ARB=(group_ACEI_ARB==2)*1) #Encode 0=ACEI, 1=ARB. According to Katie, ARB is less frequent
save(data2, pred1, meth1, bsl.covariates,file="./data/data for imputation.RData")

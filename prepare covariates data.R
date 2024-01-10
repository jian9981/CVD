
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
save(data0, file="./data/source covariate data.RData")
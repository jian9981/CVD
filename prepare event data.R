
setwd("P:/ORD_Bress_202112038D/Jian")
set.seed(123)

.libPaths(c(.libPaths(),"D:/R_Temp/Packages/3.4.0"))



library(RODBC)
library(ggplot2)
library(stringr)
library(dplyr)

odbcDataSources()
conn_odbc <- odbcConnect("VHACDWRB02_17")
sqlTables(conn_odbc,catalog="ORD_Bress_202112038D",schema = "dflt") 
flp=sqlQuery(conn_odbc, "select  * from ORD_Bress_202112038D.dflt.FOLLOWUP_CVD_LIFETEST_AA")  
#check one row per subject
flp%>%group_by(patientICN)%>%summarise(n=n())%>%ungroup()%>%summarise(max.n=max(n))
table(flp$CVD_Date>flp$Dateofdeath) #I found 248 subjects have CVD date after deathe date, chnage that to be the same as death date
flp=flp%>%mutate(death=(!is.na(Dateofdeath))*1, death.days=as.numeric(Dateofdeath-drugdateTime), last.follow.up=as.numeric(as.Date("2022-12-31") -drugdateTime), death.or.cvd=(death|CVDFlag)*1 )
flp=flp%>%mutate(CVD=if_else(!is.na(FU_CVD),1,0), CVD.flp.days=pmin(FU_CVD_Days, death.days,last.follow.up, na.rm = T ), death.or.cvd.days=CVD.flp.days) 
flp=flp%>%select(patientICN,Dateofdeath,drugdateTime,death, CVD,FU_CVD_Days, death.days,last.follow.up,death.or.cvd,  CVD.flp.days,death.or.cvd.days)
flp$multi.status.CVD.death=flp$death
flp$multi.status.CVD.death[flp$CVD==1] =2 #~40000 patients had CVD and death on the same day. Here I set CVD statue for them on this variable. 
#flp$multi.status.CVD.death=factor(flp$multi.status.CVD.death, 0:2, labels = c("censor","death","CVD"))
#checking
table(flp$multi.status.CVD.death,flp$CVD)
table(flp$multi.status.CVD.death,flp$death)

save(flp,file="./data/event data.RData")

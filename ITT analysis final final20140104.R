setwd("P:/ORD_Bress_202112038D/Jian")
source("./code/support_functions_12152023.R")
source("./code/from Chong/functions2.R")

.libPaths(c(.libPaths(),"D:/R_Temp/Packages/3.4.0"))
library(RODBC)
library(ggplot2)
library(stringr)
library(dplyr)
library(gtsummary)
library(mice)
library(PSweight)
library(glmnet)
library(cobalt)
library(WeightIt)
library(haven)
library(survminer)
library(survival)
library(cmprsk)
library(cmprskcoxmsm)
require(doParallel);
require(doRNG);
library(twang)
library(forcats)
require(survey);require(MatchIt); require(mice);library(boot);require(tidyverse);require(plyr);library(viridis); require(car);require(rlist);require(twang)
require(broom.mixed);require(ggpubr);require(survminer);require(glmnet)

load(file="./data/event data.RData")
load("./data/data for imputation.RData")

covariates=names(data2)[3:51]
fml=as.formula(paste("group_ACEI_ARB~", paste(c(covariates), collapse="+")))

#imputation
# imputed.data=mice(data2, m=10, predictorMatrix=pred1, meth=meth1, seed=1977, maxit=10)
# complete.data=complete(imputed.data, "long", include = F)
# save(complete.data,file="./data/imputated data.RData")
load("./data/imputated data.RData")


mycif <- function( data, fml ) {
  weight1=WeightIt::weightit(fml, data=data,estimand = "ATE")
  weight.trimed=WeightIt::trim(weight1,at=.99 )
  data$weight=weight.trimed$weights
  
  data=data%>%mutate(trt=group_ACEI_ARB)
  
  cif0=getcif_faster(dat=data,tx=0, evtvar="multi.status.CVD.death", dayvar="death.or.cvd.days", weightvar="weight" )
  cif1=getcif_faster(dat=data,tx=1, evtvar="multi.status.CVD.death", dayvar="death.or.cvd.days", weightvar="weight" )
  cif0=data.frame(cif0)
  cif1=data.frame(cif1)
  cif.cvd=rbind(cif0 %>%mutate(Event="CVD", Treatment="ACEI"), cif1%>%mutate(Event="CVD", Treatment="ARB"))
  
  cif.death=wtkm(dat=data, evtvar="death" , dayvar="death.days",treatvar="trt", wt= "weight")
  names(cif.death)=c( "tm1"  ,     "cif"     ,   "Treatment")
  cif.death=cif.death%>%mutate(Event="Death",Treatment=if_else(Treatment==1, "ARB","ACEI" ))
  
  cif.all=rbind(cif.cvd,  cif.death)
  temp=expand.grid(Event=c("CVD","Death"), Treatment=c("ACEI","ARB"), tm1=1:8000)
  cif.all=temp%>%left_join(cif.all)%>%arrange(Event, Treatment, tm1)%>%tidyr::fill(cif,.direction="down")
  cif.all
}




cl <- makeCluster(10) # don't use all your cores 
registerDoParallel(cl)
t1=Sys.time()
cif.ori<- foreach(i=1:10) %dopar%{
  library(dplyr)
  library(survey)
  library(questionr)
  library(survival)
  data1=complete.data%>%filter(.imp==i)
  data1=data1%>%left_join(flp)
  
  mycif(data = data1,fml=fml)
  }
stopCluster(cl)
t2=Sys.time()
t2-t1



#bootstrapping parallel parallel
cl <- makeCluster(20) # don't use all your cores 
registerDoParallel(cl)
t1=Sys.time()
cif.btr<- foreach(i=1:10) %dopar%{
  library(dplyr)
  library(survey)
  library(questionr)
  library(survival)
  require(doParallel);
  require(doRNG);
  data1=complete.data%>%filter(.imp==i)
  data1=data1%>%left_join(flp)
  all.cif<- foreach(j=1:100) %dopar%{
    library(dplyr)
    library(survey)
    library(questionr)
    library(survival)
    data.btr=data1[sample.int(dim(data1)[1],replace=T),]
    mycif(data = data.btr,fml=fml)
  }
  all.cif
}
stopCluster(cl)
t2=Sys.time()
t2-t1

save(cif.ori, cif.btr,file="./data/cif bootstrapping.RData")

cif.btr.sum=NULL
for (i in 1:10){
  cif.btr1=do.call(rbind, cif.btr[[i]])
  cif.btr1=cif.btr1%>%dplyr::group_by(Event,Treatment, tm1 )%>%dplyr::summarise(n=n(), mean.cif=mean(cif), se.cif=sd(cif))%>%dplyr::rename(cif=mean.cif)
  cif.btr1$.imp=i
  cif.btr.sum=rbind(cif.btr.sum, cif.btr1)
}
cif.btr.com=cif.btr.sum%>%dplyr::group_by(Event,Treatment, tm1 )%>%dplyr::summarise(mean.cif=mean(cif), V_w=mean(se.cif^2), V_b=(sd(cif))^2)%>%dplyr::rename(cif=mean.cif)
cif.btr.com=cif.btr.com%>%dplyr::mutate(V=V_w+(1+1/10)*V_b, se.cif=sqrt(V), lower=cif-1.96*se.cif,upper=cif+1.96*se.cif )


cif.btr.com%>%ggplot()+
  geom_step(aes(x=tm1/365, y=cif, group=interaction(Treatment, Event), color=Event,linetype=Treatment))+
  geom_step(aes(x=tm1/365, y=lower, group=interaction(Treatment, Event), color=Event,linetype=Treatment), alpha=0.5)+
  geom_step(aes(x=tm1/365, y=upper, group=interaction(Treatment, Event), color=Event,linetype=Treatment), alpha=0.2)+
  xlab("Year")+ylab("CIF")

ggsave("./results/cif for CVD and death.png")


complete.data$weight=c(rep(NA, nrow(data0)),unlist(weight.list))

cl <- makeCluster(5) # don't use all your cores 
registerDoParallel(cl)

fit.weight<- foreach(i=1:10) %dopar%{
  library(dplyr)
  WeightIt::weightit(fml, data=complete.data%>%filter(.imp==i),estimand = "ATE")
}
stopCluster(cl)





weight.list=lapply(fit.weight, "[[", "weights") 


data1=complete.data%>%filter(.imp==1)%>%dplyr::select(-one_of(c(".imp", ".id")))
# data1$Y=rnorm(n=dim(data1)[1]) #Create a random variable Y, it is needed by the PSweight function but it's values don't affact the PS weight
# PS=PSweight(fml, yname="Y", data=data1)
# propensity=PS$propensity
# data1$p0=propensity[,1]
# data1$p1=propensity[,2]
# data1=data1%>%mutate(propensity.score=if_else(group_ACEI_ARB==T, p1,p0), weight=1/propensity.score)
# data1%>%ggplot()+geom_histogram(aes(x=propensity.score))+facet_grid(group_ACEI_ARB)
# fit=glm(fml, family = binomial(link = "logit"),data=data1)
# pred=predict(fit, type="response")
# dt=data1
# dt$pred=pred
# dt=dt%>%mutate(ps=(pred^group_ACEI_ARB)*((1-pred)^(1-group_ACEI_ARB)), weight=1/ps)

weight1=WeightIt::weightit(fml, data=data1,estimand = "ATE")
#summary(data1$weight-weight1$weights) #Yeah, the two ways get same weights
#summary(dt$weight-weight1$weight)  #yeah, the two ways get same weights too
summary(weight.list[[1]])
summary(weight1$weights) 
table(weight.list[[1]]==weight1$weights) #yeah, they are the same;

complete.data$ATE=c(rep(NA, nrow(data0)),unlist(weight.list))
#check allignment
table(complete.data$ATE[complete.data$.imp==1]==weight.list[[1]]) #yeah they align

complete.data%>%filter(.imp!=0) %>%group_by(.imp)%>%summarise(quantile(ATE, 0.99))

complete.data$ATE1=unlist(tapply(complete.data$ATE, complete.data$.imp, function(x){y=x; if(sum(!is.na(x))>0){cut=quantile(x, 0.99); y[x>cut]=cut};return(y)}))
check=complete.data%>%filter(ATE>=19)%>%select(.imp, ATE, ATE1)
complete.data=complete.data%>%mutate(ATE=ATE1)%>%select(-ATE1)



cl <- makeCluster(5) # don't use all your cores 
registerDoParallel(cl)
t1=Sys.time()
bal.table<- foreach(i=1:10) %dopar%{
  library(dplyr)
  library(cobalt)
  library(WeightIt)
  weight.trimed=trim(fit.weight[[i]],at=.99 )
  t01=bal.tab(weight.trimed,abs=T,un=T, disp=c("means","sds"))
  t01$Balance
}
stopCluster(cl)
t2=Sys.time()
t2-t1 # 6min

bal.table.mean=Reduce("+",lapply(bal.table, "[", -1) )/10

write.csv(round(bal.table.mean,3), "./results/covariates balance table mean of 10 imputations after trimming.csv", row.names = T, na="")

bal.table.mean$var=row.names(bal.table.mean)
bal.table.mean=bal.table.mean%>%arrange(Diff.Un)
bal.table.mean$var=fct_inorder(bal.table.mean$var)
bal.table.mean$order=1:nrow(bal.table.mean)

m1=bal.table.mean%>%select(var, Diff.Un)%>%rename(Diff=Diff.Un)%>%mutate(weighting="Before weighting",Diff=abs(Diff))
m2=bal.table.mean%>%select(var, Diff.Adj)%>%rename(Diff=Diff.Adj)%>%mutate(weighting="After weighting",Diff=abs(Diff))
m=rbind(m1,m2)

bal1=ggplot(data=m, aes(var, Diff, group=weighting))+geom_point(aes(color=weighting))+geom_hline(yintercept = c(0.1, 0.2), color="red", lty=2)+xlab(" ")+coord_flip()
bal.plot1=bal1+scale_color_manual(values=c("black","red"),name="")+theme(legend.position = "bottom")+ylab("Absolute Standardized differences")
(fig2=bal.plot1)

png("./results/covariate balance mean of 10 imputations after trimming.png", res=300, width=2000, height=2200)
fig2
dev.off()


dt1=complete.data%>%filter(.imp==1)






t01=bal.tab(weight1,abs=T,un=T, disp=c("means","sds"))
t01a=t01$Balance
write.csv(t01$Balance, "./results/covariates balance table before trimming.csv", row.names = T)

bal.plot(weight1, var.name = "prop.score", which="both", mirror=T)
png(filename = ".\\results\\propensity score distribution before triming.png",width=4,height = 4,units="in",res = 200)
bal.plot(weight1, var.name = "prop.score", which="both", mirror=T)+xlab("Propensity Score")+ggtitle("")
dev.off()

png(filename = ".\\results\\unadjusted propensity score distribution before triming.png",width=4,height = 4,units="in",res = 200)
bal.plot(weight1, var.name = "prop.score", which="unadjusted", mirror=T)+xlab("Propensity Score")+ggtitle("")
dev.off()

varnames=read.csv(".\\data\\varnames.csv")
varnames=varnames[,c('old','new')]
png(filename = ".\\results\\covariates balance plot before trimming.png",width=8,height = 12,units="in",res = 200)
love.plot(weight1,abs=T, var.names = varnames, var.order = "unadjusted",thresholds = c(m=0.1))
dev.off()


q=quantile(weight1$weights, probs = (1:100)/100)
q=as.data.frame(q)
names(q)="percentile.before.truncation"

weight.trimed=trim(weight1,at=.99 )

q1=quantile(weight.trimed$weights, probs = (1:100)/100)
q1=as.data.frame(q1)
names(q1)="percentile.after.truncation"
q2=cbind(q,q1)
write.csv(q2, "./results/percentils of weight.csv" )


summary(weight.trimed$weights) 
summary(weight.trimed)
t1=bal.tab(weight.trimed,abs=T,un=T, disp=c("means","sds"))
t1a=t1$Balance
write.csv(t1$Balance, "./results/covariates balance table after trimming top 1 percent weight.csv", row.names = T)

png(filename = ".\\results\\covariates balance plot after top 1 percent weight trimed.png",width=8,height = 12,units="in",res = 200)
love.plot(weight.trimed,abs=T, var.names = varnames, var.order = "unadjusted",binary="std",thresholds = c(m=0.1))
dev.off()

png(filename = ".\\results\\propensity score distribution after top 1 percent weight trimed.png",width=4,height = 4,units="in",res = 200)
bal.plot(weight.trimed, var.name = "prop.score", which="both", mirror=T)+xlab("Propensity Score")+ggtitle("")
dev.off()

png(filename = ".\\results\\unadjusted propensity score distribution after top 1 percent weight trimed.png",width=4,height = 4,units="in",res = 200)
bal.plot(weight.trimed, var.name = "prop.score", which="unadjusted", mirror=T)+xlab("Propensity Score")+ggtitle("")
dev.off()





data1$weight=weight.trimed$weights

data1$ps=weight.trimed$ps


data=data1%>%left_join(flp)
data=data%>%mutate(treatment=if_else(group_ACEI_ARB==T, "ACEI","ARB"), CVD.flp.year=CVD.flp.days/365, death.or.cvd.year=death.or.cvd.days/365)
table(data$treatment)

fit=survfit(Surv(death.or.cvd.year, multi.status.CVD.death)~treatment,weights = weight, data=data)
surv.prob=tbl_survfit(fit, time=(1:20))
surv.prob1=summary(fit, time=(1:20))
surv.prob1a=as.data.frame(surv.prob1[c("strata", "time", "n.risk",  "cumhaz", "lower","upper")])
write.csv(surv.prob1a,".\\results\\cumulative incidence CVD competing death  with weight after top 1 percent weight trimed.csv", row.names = F )

plot(fit, col=c(1,2,1,2), lty=c(2,2,1,1),
     mark.time=F, lwd=2,
     xlab="Years post-index",ylab="Cumulative Incidence Function")
legend(1,0.45,c("Death:ACEI", "Death:ARB", "CVD:ACEI", "CVD:ARB"), col=c(1,2,1,2), lty=c(2,2,1,1),
       lwd=2)
png(".\\results\\cumulative incidence CVD competing death with weight after top 1 percent weight trimed.png", width = 10, height = 8, units = "in", res=200)
plot(fit, col=c(1,2,1,2), lty=c(2,2,1,1),
     mark.time=F, lwd=2,
     xlab="Years post-index",ylab="Cumulative Incidence Function")
legend(0,0.45,c("CVD:ACEI", "CVD:ARB", "Death:ACEI", "Death:ARB"), col=c(1,2,1,2), lty=c(2,2,1,1),
       lwd=2)
dev.off()


fit.cox=coxph(Surv(CVD.flp.year, CVD)~treatment,weights = weight, data=data)
ss=summary(fit.cox)
ss$conf.int
write.csv(ss$conf.int,".\\results\\CVD HR with weight after top 1 percent weight trimed.csv", row.names = F )





fit1=survfit(Surv(death.or.cvd.year, death.or.cvd)~treatment,weights = weight, data=data)
surv.prob=tbl_survfit(fit1, time=(1:20)*12)
surv.prob1=summary(fit1, time=(1:20)*12)
surv.prob1=as.data.frame(surv.prob1[c("strata", "time", "n.risk", "surv", "lower","upper")])
write.csv(surv.prob1,".\\results\\CVD or Death survival probability  with weight after top 1 percent weight trimed.csv", row.names = F )

ggsurvplot(fit1, data=data, risk.table = F, break.time.by=12, ggtheme=theme_minimal(), xlab="Months")

png(".\\results\\CVD or Death KM curve with weight after top 1 percent weight trimed.png", width = 10, height = 8, units = "in", res=200)
ggsurvplot(fit1, data=data, risk.table = F, break.time.by=12, ggtheme=theme_minimal(), xlab="Months")
dev.off()


fit1.cox=coxph(Surv(death.or.cvd.year, death.or.cvd)~treatment,weights = weight, data=data)
summary(fit1.cox)
ss=summary(fit1.cox)
ss$conf.int
write.csv(ss$conf.int,".\\results\\CVD or death HR with weight after top 1 percent weight trimed.csv", row.names = F )

data%>%group_by(treatment)%>%summarise(n=n(),total.CVD.flp.days=sum(CVD.flp.days), n.CVD=sum(CVD))%>%mutate(n.cvd.per100.yr=n.CVD*365*100/total.CVD.flp.days)
median(data$CVD.flp.days)/365.2

#try cmprskcoxmsm package
#convert logical variables to numeric variables
cols=sapply(data,is.logical)
data[,cols]=lapply(data[,cols],as.numeric)


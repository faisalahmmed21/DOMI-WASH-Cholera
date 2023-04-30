# OCV WASH DOMI Trial (Faisal Ahmmed, 2023)

library(dplyr) # For data management
library(foreign)
library(gtsummary)
library (survival) # AG model, Merginal mean model, PWP-TT, PWP-GT model

setwd("E:/E/EADII2/Analysis#9/Datasets")
analytic_0 <- read.csv("DOMI_WASH_journal.csv") 



# Utilities functions
hrci<- function(model){
  sumcof<- summary(model)$coef
  cicof_2tail<- exp(confint(model))
  cof<- c(paste(round(exp(sumcof[1,1]),2),"(",
                round(cicof_2tail[1,1],2),",",
                round(cicof_2tail[1,2],2),")", sep=""),
          round(sumcof[1,6],3))
  return(cof)
}


PE<- function(dat,pt,out,conf){
  dat[["survobj1"]]<- Surv(dat[[pt]],dat[[out]])
  model<- coxph(survobj1~dat[[conf]], cluster=cluster,data=dat)
  pe<-hrci(model)
  return(pe)
}

# Function to calculate incidence and the PE using cox-propotional hazard model

CE<- function(dat,pt,out,conf){
  participant<- table(dat[[conf]]) 
  case<- tapply(dat[[out]],dat[[conf]],sum,na.rm=T)
  persontime<- tapply(dat[[pt]],dat[[conf]],sum,na.rm=T)/365.25
  incidence<- case*1000/persontime
  pe<- PE(dat,pt,out,conf)
  row1<- c(participant[1],case[1],round(persontime[1],0),round(incidence[1],1))
  row2<-c(participant[2],case[2],round(persontime[2],0),round(incidence[2],1))
  #all_print<- rbind(cbind(row2,row1),pe)
  all_print<- c(row2,row1,pe)
  return(all_print)
}






# Bivariate relationship of available wash factors with Cholera
cl<-c(" ","Yes"," Group"," "," ","No","Group"," ","Hazard","Ratio (HR)" )
cl1<-c("n","Case","PY","IR/1000PY","n","Case","PY","IR/1000PY","Crude HR","p-val")

## Training Set
tab<-as.matrix(rbind(cl1,
                     CE(filter(analytic_0,Type=="Train"),"time","case","Flush_Toilet"),
                     CE(filter(analytic_0,Type=="Train"),"time","case","Drinking_Water"),
                     CE(filter(analytic_0,Type=="Train"),"time","case","Daily_Water_Use"),
                     CE(filter(analytic_0,Type=="Train"),"time","case","Hand_Wash"),
                     CE(filter(analytic_0,Type=="Train"),"time","case","Place_disposal")))

colnames(tab)<-cl
rownames(tab)<-c(" ","Have Flush Toilet","Safe Drinking Water","Safe Daily Water Use","Hand Wash using Soap","Have Waste Disposal Place")
kable(tab,caption = "Incidence rates and HR  of Culture Positive for Cholera")





## Incidence rates and Hazard Ratio (HR) efficacy  
analytic_0<- analytic_0 %>% mutate(wash=ifelse(wash_c=="Better",1,0))

### FUll Set
tab<-as.matrix(rbind(cl1,
                     CE(analytic_0,"time","case","wash"),
                     CE(filter(analytic_0, age_c=="0-4,yrs"),"time","case","wash"),
                     CE(filter(analytic_0, age_c=="5-14,yrs"),"time","case","wash"),
                     CE(filter(analytic_0, age_c=="15+,yrs"),"time","case","wash")))

colnames(tab)<-c("Better","WasH"," "," ","Not Better","WasH"," "," ","Hazards Ratio","(HR)")
rownames(tab)<-c(" ", "All","0-4,yrs","5-14,yrs","15+, yrs")
kable(tab,caption = "")


### Adjusted PE
coxph(Surv(time,case) ~ wash + age_c + occup + expend_m + clinic_dis_m, cluster=cluster ,data=analytic_0)
coxph(Surv(time,case) ~ wash + expend_m, cluster=cluster ,data=filter(analytic_0, age_c=="0-4,yrs"))
coxph(Surv(time,case) ~ wash + occup + expend_m + clinic_dis_m, cluster=cluster ,data=filter(analytic_0, age_c=="5-14,yrs"))
coxph(Surv(time,case) ~ wash + hindu +  expend_m + clinic_dis_m, cluster=cluster ,data=filter(analytic_0, age_c=="15+,yrs"))


## Supplement Table 1.Baseline characteristics of the training, and validation subpopulations in the control population

analytic_0 %>% select(Type, age_c, sex, hindu, own, occup, expend_m, clinic_dis_m) %>%
  tbl_summary(by=Type,
              type = all_dichotomous() ~ "categorical") %>%
  add_overall() %>% 
  bold_labels() %>%
  add_p() 



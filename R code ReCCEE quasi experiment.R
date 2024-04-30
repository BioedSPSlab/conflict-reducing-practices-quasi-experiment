#Data preparation
library(nlme)
library(lme4)
library(dplyr)
library(MASS)
library(gtsummary)
library(magrittr)
library(tidyverse)
library(jtools)
library(huxtable)
library(broom.mixed) 
library(readr)
library(WebPower) 
library(interactions)
library(ltm)
library(sjPlot)
library(TAM)
library(sjmisc)
library(emmeans)
library(data.table)
library(readr)
library(dplyr)
#DATA
#Administrated to >11,000 undergraduate biology students from 72  group (courses).
#Analyzed and consented to both the pre- and post-course survey for four semester.
#In other words, participants that were missing any of the variables were removed from the analysis (n=6896).
#Students are nested within courses (random effect group). We used minimum 20 students each random effect groups as recommended by Simmons et al (2011), resulting 55 groups (n=6719) as data1
#Note: This data only include courses that minimun of 20 students 
X23_05_19_ReCCEE_quasi_data <- read_csv("C:/Users/rqa2a/OneDrive - Middle Tennessee State University/Research - ReCCEE Projects/ReCCEE quasi experiment/ReCCEE Manuscript/Next Journal/23 05 19 ReCCEE quasi data.csv")
data1<-X23_05_19_ReCCEE_quasi_data
rm(X23_05_19_ReCCEE_quasi_data)
################ CFA ################ 
cceeps.model <- ' rsrm =~ rsrm1 + rsrm2 + rsrm3 + rsrm4
aut =~ aut1 + aut2 + aut3'
cceeps.fit <- cfa(cceeps.model, data=data1, std.lv=TRUE, estimator = "WLS")
summary(cceeps.fit, standardized=TRUE, fit.measures=TRUE)
fitMeasures(cceeps.fit, c("chisq", "df", "pvalue", "agfi", "nfi", "tli", "cfi", "rmsea", "srmr"))
lavaanPlot(model = cceeps.fit, node_options = list(shape = "box", fontname = "Helvetica"), edge_options = list(color = "grey"), coefs = TRUE, sig = .05)

######################################
pre.hum<-data1[c("human1", "human2", "human3", "human4", "human5", "human6", "human7", "human8")]
post.hum<-data1[c("phuman1", "phuman2", "phuman3", "phuman4", "phuman5", "phuman6", "phuman7", "phuman8")]
cronbach.alpha(pre.hum)
cronbach.alpha(post.hum)
#rasch
mod.pre<- tam.mml(pre.hum, irtmodel= "PCM2", control= list(snodes=1000, maxiter=50))
mod.post<- tam.mml(post.hum,irtmodel= "PCM2", control=list(snodes=1000,maxiter=50))

mod.pre$EAP.rel
ability.pre<- tam.wle(mod.pre)
fit.pre.accept<- TAM::msq.itemfit( mod.pre )
summary(fit.pre.accept)

mod.post$EAP.rel
ability.post<- tam.wle(mod.post)
fit.post.accept<- TAM::msq.itemfit( mod.post )
summary(fit.post.accept)

data1$humtheta2<-ability.pre$theta
data1$phumtheta2<-ability.post$theta

# Validity evolution understanding
evound_pre<-data1[c("evound1", "evound2", "evound3", "evound4", "evound5", "evound6", "evound7",
                                            "evound8", "evound9", "evound10", "evound11", "evound12", "evound13", "evound14")]
cronbach.alpha(evound_pre)

mod.evound<-tam.mml(evound_pre, irtmodel= "PCM", control= list(snodes=1000,  maxiter=50))
mod.evound$EAP.rel
fit.evound<- TAM::msq.itemfit( mod.evound )
summary(fit.evound)
ability.evound<- tam.wle(mod.evound)
data1$evoundtheta2<-ability.evound$theta

#validity religiosity
religiosity<-data1[c("rel1", "rel2", "rel3", "rel4")]
cronbach.alpha(religiosity)

mod.rel<-tam.mml(religiosity, irtmodel= "PCM", control= list(snodes=1000,  maxiter=50))
mod.rel$EAP.rel
fit.rel<- TAM::msq.itemfit( mod.rel )
summary(fit.rel)
wle.rel<- tam.wle(mod.rel)
data1$reltheta2<-wle.rel$theta

#autonomy
aut_data<-data1[c("aut1", "aut2", "aut3")]
cronbach.alpha(aut_data)

mod.aut<-tam.mml(aut_data, irtmodel= "PCM", control= list(snodes=1000, maxiter=50))
plot(mod.aut, items=1:3, export=FALSE)
mod.aut$EAP.rel
fit.aut <- TAM::msq.itemfit( mod.aut )
summary(fit.aut)
rel.aut<- tam.wle(mod.aut)
data1$auttheta2<-rel.aut$theta
#compute thresholds
thresh3 <- TAM::IRT.threshold(mod.aut)
print(thresh3)
IRT.WrightMap(thresh3)

#Religious role model
rsrm_all<-data1[c("rsrm1", "rsrm2", "rsrm3", "rsrm4")]
cronbach.alpha(rsrm_all)

mod.rsrm<-tam.mml(rsrm_all, irtmodel= "PCM", control= list(snodes=1000, maxiter=50))
mod.rsrm$EAP.rel
plot(mod.rsrm, items=1:4, export=FALSE)

fit.rsrm <- TAM::msq.itemfit( mod.rsrm )
summary (fit.rsrm)
rs.rm<- tam.wle(mod.rsrm)
data1$rsrmtheta2<-rs.rm$theta
# compute thresholds
thresh2 <- TAM::IRT.threshold(mod.rsrm)
print(thresh2)
IRT.WrightMap(thresh2)

dev.off()
#ANALYSES
#Recode
data1$religion_recode<-data1$religion
data1$religion_recode<-dplyr::recode(data1$religion_recode, '1'="Other",'2'="Christian", '3'= "Other",'4'="Other",'5'= "Other", '6'="No religion", '7'="Other", '8'="Other")
data1$gender<- dplyr::recode(data1$gender, '1'="Female",'2'="Male", '3'= "non binary",'4'="other",'5'= "other")
data1$race<- dplyr::recode(data1$race, '1'="Black/Hispanic/Native",'2'="Asian",'3'= "Multiracial", '4'="White",)
data1$biomajor<- dplyr::recode(data1$biomajor, '1'="Biology major", '0'="Non Biology Major")
data1$semester<- dplyr::recode(data1$semester, '2'="Spring 2020",'3'="Fall 2020",'4'= "Spring 2021", '5'="Fall 2018")

#make new religiosity category 
x<- data1$reltheta2
data1$religiosity2<- case_when(x > mean(x)+sd(x) ~ "+ 1 SD",
                               x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "Mean",
                               x < mean(x)-sd(x) ~ "- 1 SD") 
data1$religiosity2 <- factor(data1$religiosity2, levels= c("+ 1 SD", "Mean", "- 1 SD"))

data1$relcat<-paste(data1$religion_recode, data1$religiosity2) #Merge together

#Create new category student according to 
data1$religion_recode2<-dplyr::recode(data1$relcat,
                                      'Christian - 1 SD'="No religion",
                                      'Other - 1 SD'="No religion",
                                      'Christian Mean'="Christian",
                                      'Christian + 1 SD'="Christian",
                                      'Other Mean'="Other",
                                      'Other + 1 SD'="Other",
                                      'No religion Mean'="No religion",
                                      'No religion + 1 SD'="No religion", 
                                      'No religion - 1 SD'="No religion")

dem.acc <- data1 %>% select(semester, gender, race, biomajor, religion_recode, relcat, religion_recode2)


dem.acc %>% tbl_summary(label = list(gender ~ "Gender", semester ~ "Semester",
                                     race ~ "Race/ethnicity", biomajor ~ "Bio Major", religion_recode2 ~ "Religion final"))
#see course
as.factor(data1$course)
frequency.course<-table(data1$course)
frequency.course <- as.data.frame(frequency.course)

#change reference group
data1$gender <- as.factor(data1$gender)
data1$gender<- relevel (data1$gender, ref =  "Male")
data1$race <- as.factor(data1$race) 
data1$race <- relevel(data1$race, ref = "White")
data1$semester <- as.factor(data1$semester)
data1$semester <- relevel(data1$semester, ref = "Spring 2020")
data1$biomajor <- as.factor(data1$biomajor) 
data1$biomajor <- relevel(data1$biomajor, ref = "Biology major")
data1$religion_recode2 <- as.factor(data1$religion_recode2) 
data1$religion_recode2<- relevel(data1$religion_recode2, ref =  "No religion")
data1$course<- as.factor (data1$course)

#Check Power
#To calculate the statistical power given sample size and effect size:
wp.crt2arm(f = 0.25, n = 20 , J =55 , icc = 0.02, alpha = 0.05, power = )

#We performed a power analysis within a multilevel model structure using the WebPower package in R and found we have 93% Power to detect effects of d = 0.25 assuming a moderate amount of the overall variance will occur at the classroom level (between-cluster variance = 0.02)
library(Hmisc)
deff(data1$phumtheta2, data1$course)


###############


#Human evolution acceptance
#Model 0:
hum.semester<- lmer(phumtheta2 ~ gender + race + semester + biomajor + humtheta2 + evoundtheta2 + (1|course), data=data1, REML = F)
hum<- lmer(phumtheta2 ~ gender + race  + biomajor + humtheta2 + evoundtheta2 + (1|course), data=data1, REML = F)
summ(hum.semester)
#check if semester significant, if no, next model exclude semester 
#Model 1:
hum<- lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2+ (1|course), data=data1, REML = F)
#Model - rsrm
hum.rsrm1<- lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2 + rsrmtheta2 + (1|course), data=data1, REML = F)
hum.rsrm2<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   +  rsrmtheta2*religion_recode2  + (1|course), data=data1, REML = F)
hum.rsrm3<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   + reltheta2*rsrmtheta2+ (1|course), data=data1, REML = F)
hum.rsrm4<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   + reltheta2*rsrmtheta2*religion_recode2 + (1|course), data=data1, REML = F)

#Model - aut
hum.aut1<- lmer(phumtheta2 ~ gender + race  + biomajor + humtheta2 + evoundtheta2 + auttheta2 + (1|course), data=data1, REML = F)
hum.aut2<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   +  auttheta2*religion_recode2  + (1|course), data=data1, REML = F)
hum.aut3<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2  + reltheta2*auttheta2+ (1|course), data=data1, REML = F)
hum.aut4<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   + reltheta2*auttheta2*religion_recode2 + (1|course), data=data1, REML = F)

##estimate marginal result
hum.rsrm.emm4<-emmeans(hum.rsrm4, "religion_recode2")
hum.rsrm.emm4
hum.aut.emm4<-emmeans(hum.aut4, "religion_recode2")
hum.aut.emm4


#Export result
tab_model(hum.rsrm1, hum.rsrm2, hum.rsrm3, hum.rsrm4, show.se = TRUE, show.stat = TRUE, pred.labels = c("intercept", "gender"))
tab_model(hum, hum.aut1, hum.aut2, hum.aut3, hum.aut4, show.se = TRUE, show.stat = TRUE, pred.labels = c("intercept", "gender"))

#plot result
###RSRSM
jpeg("hum rsrsm 1.jpeg", 
     width = 24, height = 12, units = "cm", res = 200, 
     quality = 100)
interactions::interact_plot(hum.rsrm4, pred =rsrmtheta2, 
                            modx = reltheta2, mod2=religion_recode2, 
                            data = data1,
                            interval = TRUE,
                            modx.values = c(0,2.16),
                            mod2.labels =c("Non-Religious Students", "Christian Students", "Other-Religion Students"),
                            x.label = "Students' Perception of Religious Role Model",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-10,10) + theme(axis.text = element_text(size = 12))    
# Close device
dev.off()
#no rel
jpeg("hum rsrsm 2.jpeg", 
     width = 24, height = 12, units = "cm", res = 200, 
     quality = 100)
interactions::interact_plot(hum.rsrm4, pred =rsrmtheta2, 
                            modx = reltheta2, mod2=religion_recode2, 
                            data = data1,
                            interval = TRUE,
                            modx.values = c(0,0),
                            colors = "CUD Bright",
                            mod2.labels =c("Non-Religious Students", "Christian Students", "Other-Religion Students"),
                            x.label = "Students' Perception of Religious Role Model",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-10,10)+ theme(axis.text = element_text(size = 12))    
dev.off()
#all
jpeg("hum rsrsm 3.jpeg", 
     width = 10, height = 8, units = "cm", res = 200, 
     quality = 100)
interactions::interact_plot(hum.rsrm4, pred =rsrmtheta2, 
                            modx = reltheta2,
                            data = data1,
                            pred.labels = c("Religiosity"),
                            interval = TRUE,
                            colors = "Greys",
                            modx.values = c(0,0),
                            x.label = "Students' Perception of Religious Role Model",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-10,10)
dev.off()


###AUTONOMY
#christian and other
jpeg("hum aut 1.jpeg", 
     width = 24, height = 12, units = "cm", res = 200, 
     quality = 100)

interactions::interact_plot(hum.aut4, pred =autheta2, 
                            modx = reltheta2, mod2=religion_recode2, 
                            data = data1,
                            interval = TRUE,
                            mod2.labels =c("Non-Religious Students", "Christian Students", "Other-Religion Students"),
                            modx.values = c(0,2.16),
                            x.label = "Students' Perception of Autonomy",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-6,2) + theme(axis.text = element_text(size = 12))    
dev.off()

jpeg("hum aut 2.jpeg", 
     width = 24, height = 12, units = "cm", res = 200, 
     quality = 100)
interactions::interact_plot(hum.aut4, pred =autheta2, 
                            modx = reltheta2, mod2=religion_recode2, 
                            data = data1,
                            interval = TRUE,
                            modx.values = c(0,0),
                            colors = "CUD Bright",
                            mod2.labels =c("Non-Religious Students", "Christian Students", "Other-Religion Students"),
                            x.label = "Students' Perception of Autonomy",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-6,2)+ theme(axis.text = element_text(size = 12))    
dev.off()

jpeg("hum aut 3.jpeg", 
     width = 10, height = 8, units = "cm", res = 200, 
     quality = 100)
interactions::interact_plot(hum.aut4, pred =autheta2, 
                            modx = reltheta2,
                            data = data1,
                            pred.labels = c("Religiosity"),
                            interval = TRUE,
                            colors = "Greys",
                            modx.values = c(0,0),
                            x.label = "Students' Perception of Autonomy",
                            y.label = "Predicted Human Evolution Acceptance")  +ylim(-1,1) + xlim(-6,2)
dev.off()
#hum relneg
#christian and other

#non rel
#all


#slope
#plot result

library(psych)
library(ggplot2)
library(nlme) #for mixed effects models
library(reghelper)  #for probing interactions

##SLOPE ANALYSIS

round(sd(x), 2)

sim.hum.rsrm<-sim_slopes(hum.rsrm4, pred=rsrmtheta2, modx =reltheta2 , mod2 = religion_recode2,  modx.values = c(0,2.16))

sim.hum.aut<-sim_slopes(hum.aut4, pred=autheta2, modx =reltheta2 , mod2 = religion_recode2,  modx.values = c(0,2.16))

t1<-tidy(sim.hum.rsrm)
t5<-tidy(sim.hum.aut)

list.slope<-list(t1,t5)
slope.an<-rbindlist(list.slope)
slope.an$p.value<-round(slope.an$p.value, 3)
write.csv(slope.an, "slope analysis 23 05 11.csv")


#####
aut_data<-data1[c("aut1", "aut2", "aut3")]
aut_data$course<-data1$course
aut_data$course<-as.factor(aut_data$course)
aut_data$aut<-rowMeans(aut_data[1:3], na.rm = FALSE)
aut.melt<-aut_data %>% select(aut, course) %>%
  melt(id='course')
str(aut.melt)
# Get the unique levels of the course variable
unique_courses <- unique(aut.melt$course)
# Create a named vector for the new levels
new_levels <- setNames(as.character(1:length(unique_courses)), levels(unique_courses))
# Relevel the course variable
aut.melt$course <- factor(aut.melt$course, levels = unique_courses, labels = new_levels)
# Calculate mean and standard deviation for each course
course_stats_aut<- aut.melt %>%
  group_by(course) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))
overall_mean.aut <- mean(aut.melt$value, na.rm = TRUE)
# Create the bar plot
ggplot(course_stats_aut, aes(x = course, y = mean)) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Course", y = "Mean Autonomy") +
  theme_minimal() +geom_hline(yintercept = overall_mean.aut, linetype = "dashed", color = "black")

ggplot(course_stats_aut, aes(x = course, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "black", size=0.5) +
  geom_point(color = "black", size = 2) +
  labs(x = "Course", y = "Inclusive teaching through\n decision autonomy") +
  theme_minimal() + 
  geom_hline(yintercept = mean(course_stats_aut$mean), linetype = "dashed", color = "black") +
  scale_y_continuous(breaks=c(1,2,3,4,5,6), limits=c(1,6.5))  + theme(axis.title.y = element_text(size = 13))
#boxplot intead
ggplot(aut.melt, aes(x = course, y = value)) +
  geom_boxplot(fatten=NULL, outlier.shape=NA) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.75, size = 1, linetype = "solid") +
  labs(x = "Course", y = "Inclusive teaching through\n decision autonomy") + 
  theme(axis.title.x = element_text(size = 14, face = "bold", family = "Arial", color = "gray"),
        axis.title.y = element_text(size = 12, face = "bold", family = "Arial", color = "gray"),
        axis.text.x = element_text(size = 12, face = "bold", family = "Arial", color = "gray"),
        axis.text.y = element_text(size = 12, face = "bold", family = "Arial", color = "gray")) +
  theme_minimal() +
  geom_hline(yintercept = mean(course_stats_rsrm$mean), linetype = "dashed", color = "black")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6), limits=c(1,6)) + theme(axis.title.y = element_text(size = 16))


rsrm_all$course<-data1$course
rsrm_all$course<-as.factor(rsrm_all$course)
rsrm_all$rsrm<-rowMeans(rsrm_all[1:4], na.rm = FALSE)
rsrm.melt <- melt(rsrm_all, id='course')
rsrm.melt<-rsrm.melt%>%filter(variable=="rsrm")
# Get the unique levels of the course variable
unique_courses_rsrm <- unique(rsrm.melt$course)
# Create a named vector for the new levels
new_levels_rsrm <- setNames(as.character(1:length(unique_courses_rsrm)), levels(unique_courses_rsrm))
# Relevel the course variable
rsrm.melt$course <- factor(rsrm.melt$course, levels = unique_courses_rsrm, labels = new_levels_rsrm)

# Calculate mean and standard deviation for each course
course_stats_rsrm <- rsrm.melt %>%
  group_by(course) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))
# Calculate the overall mean
overall_mean.rsrm <- mean(rsrm.melt$value, na.rm = TRUE)

# Create the bar plot
ggplot(course_stats_rsrm, aes(x = course, y = mean)) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Course", y = "Mean Compatibility") +
  theme_minimal() + geom_hline(yintercept = overall_mean.rsrm, linetype = "dashed", color = "black")
#mean wishker plot
ggplot(course_stats_rsrm, aes(x = course, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "black", size=0.5) +
  geom_point(color = "black", size = 2) +
  labs(x = "Course", y = "Compatibility with \n role models") + 
  theme_minimal() +
  geom_hline(yintercept = mean(course_stats_rsrm$mean), linetype = "dashed", color = "black")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6), limits=c(1,6.5)) + theme(axis.title.y = element_text(size = 16))
#boxplot
ggplot(rsrm.melt, aes(x = course, y = value)) +
  geom_boxplot(fatten=NULL, outlier.shape=NA) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.75, size = 1, linetype = "solid") +
  labs(x = "Course", y = "Compatibility with \n role models") + 
  theme(axis.title.x = element_text(size = 14, face = "bold", family = "Arial", color = "gray"),
        axis.title.y = element_text(size = 16, face = "bold", family = "Arial", color = "gray"),
        axis.text.x = element_text(size = 12, face = "bold", family = "Arial", color = "gray"),
        axis.text.y = element_text(size = 12, face = "bold", family = "Arial", color = "gray")) +
  theme_minimal() +
  geom_hline(yintercept = mean(course_stats_rsrm$mean), linetype = "dashed", color = "black")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6), limits=c(1,6)) + theme(axis.title.y = element_text(size = 16))

#median and mean
ggplot(rsrm.melt, aes(x = course, y = value)) +
  geom_boxplot(outlier.shape=NA, color="gray")+
  stat_summary(fun.y = mean, geom = "point", shape=16, size = 2, color="red") +
  labs(x = "Course", y = "Compatibility with \n role models") + 
  theme_minimal() +
  geom_hline(yintercept = mean(course_stats_rsrm$mean), linetype = "dashed", color = "red")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6), limits=c(1,6.5)) + theme(axis.title.y = element_text(size = 16))


#
hum.rsrm4<-lmer(phumtheta2 ~ gender + race + biomajor + humtheta2 + evoundtheta2   + reltheta2*rsrmtheta2*religion_recode2 + (1|course), data=data1, REML = F)
qqnorm(residuals(hum.rsrm4))
data1$predict.rsrm<-predict(hum.rsrm4)

ggplot(data1, aes(rsrmtheta2, predict.rsrm,col = course)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = F, alpha=0.5) +
  geom_smooth(method = 'lm', se = F, aes(group = 1)) + theme_minimal()+
  facet_wrap(~religion_recode2) # one facet per class

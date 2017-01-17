#####################################################################################
################## Code for analysing data on vegetation plots 2015 #################
#####################################################################################

#Reading data 
data_veg<-read.table("./data/clean/datafile_vegplots.txt",header=T,sep="\t",dec=",")
head(data_veg)
tail(data_veg)
str(data_veg)

attach(data_veg)
data_veg$plot<-as.factor(data_veg$plot)

# #Construct new variable
# #Is plant higher than vegetation?
# data_veg$higher_veg_max<-as.factor(with(data_veg,ifelse(diff_veg_h_max_shoot_h<0,"1","0")))
# data_veg$higher_veg_mean<-as.factor(with(data_veg,ifelse(diff_veg_h_mean_shoot_h<0,"1","0")))
# 
#Construct new variable: attack
data_veg$attack<-as.factor(with(data_veg,ifelse(n_eggs_max>0,"1","0")))
# 
# #Construct new variable: redants_pres
# data_veg$redants_pres<-as.factor(with(data_veg,ifelse(n_redants>0,"1","0")))
# 
# #Construct new variable: prop_pred
# data_veg$prop_pred<-with(data_veg,n_predated/n_fl_corrected)
# 
#Construct new variable: observation-level random effect
data_veg$id<-1:100

# summary(data_veg)
# names(data_veg)
# 
# plot(attack,most_adv)
# summary(lm(fruit_set~attack))
# plot(attack,n_fruits_total)

#Correlations
library(PerformanceAnalytics)
data_vegcorr<-data_veg[c(4:10,17:34)]
cor(data_vegcorr,use="pairwise.complete.obs" )
write.table(cor(data_vegcorr,use="pairwise.complete.obs" ),"clipboard",sep="\t")
chart.Correlation(data_vegcorr,pch=20)

#PCA
PCA1<-prcomp(~n_fl_corrected+most_adv+shoot_h,center=T,scale=T)
plot(PCA1)
biplot(PCA1)
PCA1
summary(PCA1)

###########################################################################

#Boxplots
for (i in 4:ncol(data_veg)){boxplot(as.numeric(data_veg[,i]),main=paste(colnames(data_veg)[i]))}

#Histograms
for (i in 4:ncol(data_veg)){hist(as.numeric(data_veg[,i]),main=paste(colnames(data_veg)[i]))}

############################################################################

#Standardize X variables
z <- scale(data_veg[,c("bud_h","shoot_h","veg_h_mean","diff_veg_h_max_shoot_h","diff_veg_h_mean_shoot_h",
                    "phen_index","most_adv","n_intact_fruits","n_predated","fruit_set",
                    "n_fl_corrected","n_eggs_max","n_redants","dist_closest_redants")])
data_veg$z.bud_h <- z[,1]
data_veg$z.shoot_h <- z[,2]
data_veg$z.veg_h_mean <- z[,3]
data_veg$z.diff_veg_h_max_shoot_h <- z[,4]
data_veg$z.diff_veg_h_mean_shoot_h <- z[,5]
data_veg$z.phen_index <- z[,6]
data_veg$z.most_adv <- z[,7]
data_veg$z.n_intact_fruits <- z[,8]
data_veg$z.n_predated <- z[,9]
data_veg$z.fruit_set <- z[,10]
data_veg$z.n_fl_corrected <- z[,11]
data_veg$z.n_eggs_max <- z[,12]
data_veg$z.n_redants <- z[,13]
data_veg$z.dist_closest_redants <- z[,14]

head(data_veg)
summary(data_veg)

###############################################################################
#Differences in variables between treatments?

summary(lm(phen_index~treat))
summary(lm(most_adv~treat))
summary(lm(shoot_h~treat))
summary(lm(n_intact_fruits~treat))
summary(lm(fr_int_h_max~treat))
summary(lm(n_predated~treat))
summary(lm(n_redants~treat))
plot(treat,n_redants,xlab="Treatment",ylab="Number of Myrmica ants")  #Do ants prefere to forage where there is higher vegetation?
summary(lm(n_eggs_max~treat))
summary(lm(n_fl_corrected~treat))
summary(lm(fruit_set~treat))

###############################################################################
#Construct models for interaction (intensity)
#data_veg_comp<-data_veg[complete.cases(data_veg[37:50]),]

#attach(data_veg_comp)
library(lme4)
library(MuMIn)
library(car)
model1v<-glm(attack~(z.shoot_h+z.most_adv+z.n_fl_corrected+z.n_redants)
             *treat,family="binomial",na.action="na.fail")
model1v<-glm(attack~(z.shoot_h+z.most_adv+z.n_fl_corrected)
             *treat,family="binomial",na.action="na.fail")
summary(model1v)
r.squaredLR(model1v) 
vif(model1v) #High vif for ants

#Model selection
models1v<-dredge(model1v)
subset(models1v, delta <2) 
importance(models1v) 

model1vavg<-model.avg(models1v, subset=delta <2,revised.var = TRUE)
summary(model1vavg)

#Separate
model1CUT<-glm(attack~z.shoot_h+z.most_adv+z.n_fl_corrected+z.n_redants,
               family="binomial",na.action="na.fail",data=subset(data_veg,treat=="Cut"))
summary(model1CUT)

model1UNCUT<-glm(attack~z.shoot_h+z.most_adv+z.n_fl_corrected+z.n_redants,
                 family="binomial",na.action="na.fail",data=subset(data_veg,treat=="Uncut"))
summary(model1UNCUT)

#Response=n_eggs_max

model2v<-glm(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected+z.n_redants)
             *treat,family="poisson",na.action="na.fail")
summary(model2v)
vif(model2v) #High vif for ants
model2v<-glm(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
             *treat,family="poisson",na.action="na.fail")
model2v<-glm(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
             *treat,family="quasipoisson",na.action="na.fail")
model2v<-glm.nb(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
             *treat,na.action="na.fail")
model2v<-hurdle(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
                *treat,dist="negbin",na.action="na.fail")
summary(model2v)
vif(model2v) 
1 - pchisq(summary(model2v)$deviance, 
           summary(model2v)$df.residual
) #Model does not fit the data! Both with poisson and quasipoisson

with(subset(data_veg,treat=="Cut"),plot(most_adv,n_eggs_max))
with(subset(data_veg,treat=="Cut"),abline(lm(n_eggs_max~most_adv)))
with(subset(data_veg,treat=="Uncut"),plot(most_adv,n_eggs_max))
with(subset(data_veg,treat=="Uncut"),abline(lm(n_eggs_max~most_adv)))


model2vnb<-glm.nb(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
             *treat,na.action="na.fail")
summary(model2vnb)
vif(model2vnb)
1 - pchisq(summary(model2vnb)$deviance, 
           summary(model2vnb)$df.residual
) #Model fits the data! But no significant effects :(

#Include plot nested in treatment
#Poisson model
model2v_n<-glmer(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
           *treat+(1|plot),na.action="na.fail",data=data_veg,family="poisson")
#model2v_n<-glmer(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
#*treat+(1|treat:plot),na.action="na.fail",data=data_veg,family="poisson")
#Gives exactly the same output: Function is smart enought to see that since 
#plot 1 is only found in conjunction with treatment Cut, plot is nested in treatment
summary(model2v_n)
r.squaredGLMM(model2v_n) 
overdisp.glmer(model2v_n) #Overdispersion
library(blmeco) 
dispersion_glmer(model2v_n) #it shouldn't be over 1.4
overdisp_fun(model2v_n)
library(aods3)
gof(model2v_n)
(sum(residuals(model2v_n,"pearson")^2))/df.residual(model2v_n) #Overdisp. param.

#Lognormal with observation-level random effect
model2v_n<-glmer(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
                 *treat+(1|plot)+(1|id),na.action="na.fail",data=data_veg,family="poisson")
summary(model2v_n) #NS
r.squaredGLMM(model2v_n) 
overdisp.glmer(model2v_n) #No overdispersion

#Negative binomial
model2vnb_n<-glmer.nb(n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
                      *treat+(1|plot),na.action="na.fail",data=data_veg)
summary(model2vnb_n) #NS
#Same smartness here

#Quasipoisson
model2v_n<-glmmPQL(fixed=n_eggs_max~(z.shoot_h+z.most_adv+z.n_fl_corrected)
                   *treat,random=~1|plot,na.action="na.fail",data=data_veg,
                   family="quasipoisson")
#And here
summary(model2v_n) #Only most_adv is significant

#Model selection in negative binomial
models2vnb_n<-dredge(model2vnb_n)
subset(models2vnb_n, delta <2) 
importance(models2vnb_n) 

model2vnb_n_avg<-model.avg(models2vnb_n, subset=delta <2,revised.var = TRUE)
summary(model2vnb_n_avg)

#With attack
model2v_n<-glmer(attack~(z.shoot_h+z.most_adv+z.n_fl_corrected)
                 *n_redants+(1|plot),na.action="na.fail",data=data_veg,family="binomial",
                 control=glmerControl(optimizer="bobyqa"))


###################################################################################33

xyplot(n_eggs_max~z.most_adv | treat, 
       panel=function(x, y){ 
         panel.xyplot(x, y) 
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       } 
)

xyplot(n_eggs_max~z.most_adv | n_redants, 
       panel=function(x, y){ 
         panel.xyplot(x, y) 
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       } 
)

xyplot(n_eggs_max~z.shoot_h | treat, 
       panel=function(x, y){ 
         panel.xyplot(x, y) 
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       } 
)

xyplot(n_eggs_max~z.n_redants | treat, 
       panel=function(x, y){ 
         panel.xyplot(x, y) 
         panel.lmline(x, y, lty = 2)  # Least squares broken line
       } 
)

##########################################################################
#Try piecewiseSEM
library(piecewiseSEM)

modList<-list(
  glmer.nb(n_eggs_max~(most_adv+n_fl_corrected+shoot_h)*n_redants+
             (1|plot),na.action="na.fail",data=data_veg,
           control=glmerControl(optimizer="bobyqa")),
  glmer.nb(n_redants~treat+(1|plot),na.action="na.fail",data=data_veg)
)

sem.fit(modList, data_veg,
        corr.errors = c("most_adv~~n_fl_corrected","shoot_h~~n_fl_corrected","most_adv~~shoot_h"))
sem.model.fits(modList)
sem.coefs(modList, data_veg,
          corr.errors = c("most_adv~~n_fl_corrected","shoot_h~~n_fl_corrected","most_adv~~shoot_h"))

model1<-glmer(n_eggs_max~n_redants+(1|id),family="poisson") #Nearly * effect
(sum(residuals(model1,"pearson")^2))/df.residual(model1) #Overdisp. param.

model1<-glmer.nb(n_redants~treat+(1|plot),na.action="na.fail",data=data_veg)
model2<-glmer(n_redants~treat+(1|plot),na.action="na.fail",data=data_veg,family="poisson")
summary(model1) #Significant effect of treatment on ants
summary(model2)
AIC(model1,model2)

model1<-glm.nb(n_eggs_max~(most_adv+n_fl_corrected+shoot_h)*n_redants,na.action="na.fail",data=data_veg)
model2<-glmer.nb(n_eggs_max~(z.most_adv+z.n_fl_corrected+z.shoot_h)*z.n_redants+
                   (1|plot),na.action="na.fail",data=data_veg,
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
model2<-glmer.nb(n_eggs_max~(most_adv+n_fl_corrected+shoot_h)*n_redants+
                   (1|plot),na.action="na.fail",data=data_veg,
                 control=glmerControl(optimizer="bobyqa"))
model3<-glmer.nb(n_eggs_max~most_adv*n_redants+(1|plot),na.action="na.fail",data=data_veg)

summary(model1)
summary(model2)
summary(model3)
AIC(model1,model2,model3)

(sum(residuals(model2,"pearson")^2))/df.residual(model2) #Overdisp. param.

#Convergence problems!

#Try with PCA for traits

PCA1<-prcomp(~most_adv+n_fl_corrected+shoot_h,center=T,scale=T)
plot(PCA1)
biplot(PCA1)
PCA1
summary(PCA1)
data_veg$pca1<-predict(PCA1)[,1]

#Models using pca1 instead of indidividual
model2v_n<-glmer(n_eggs_max~pca1*n_redants+(1|plot),na.action="na.fail",data=data_veg,family="poisson")
summary(model2v_n)
(sum(residuals(model2v_n,"pearson")^2))/df.residual(model2v_n) #Overdisp. param.

model2v_n<-glmer(n_eggs_max~pca1*n_redants+(1|plot)+(1|id),na.action="na.fail",data=data_veg,family="poisson")
summary(model2v_n)
(sum(residuals(model2v_n,"pearson")^2))/df.residual(model2v_n) #Overdisp. param.

model2vnb_n<-glmer.nb(n_eggs_max~pca1*n_redants+(1|plot),na.action="na.fail",data=data_veg)
summary(model2vnb_n)
(sum(residuals(model2vnb_n,"pearson")^2))/df.residual(model2vnb_n) #Overdisp. param.


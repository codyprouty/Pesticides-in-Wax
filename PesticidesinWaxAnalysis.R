##

#Analysis for lithium chloride manuscript

#Contact: cprouty@ufl.edu    https://scholar.google.com/citations?user=PpeDx78AAAAJ&hl=en

##

#load packages
library(lme4)
library(lsmeans)
library(multcomp)
library(survival)
library(brglm2)
##

#Load datasheet
Wax <- read.csv("PesticidesinWax.csv")
##

#Data organization
names(Wax)[1] <- "Treatment"
Wax$SurvObj <- with(Wax, Surv(DeathDay, Survived == 0))
Survival <- survfit(SurvObj ~ Treatment, data = Wax)

#Subsetted for the analysis of day of death
Dead <- subset(Wax, Survived == 0)
##

#Effect of treatment on death day
DD <- lm(DeathDay ~ Treatment, data= Dead)
anova(DD)
summary(DD)

lsm<-lsmeans (DD, list( ~ Treatment))
cld(lsm)
###

#Effect of treatment on survival
Sur <- glm(Survived ~ Treatment, 
          family = binomial, Wax, method=brglmFit)
anova(Sur, test = "Chisq")
summary(Sur)
###

#Graph of survival across treatments
Survival <- survfit(SurvObj ~ Treatment, data = Wax)
plot(Survival, col=c("#BBBBBB","#66CCEE", "#4477AA", "#228833","#999933", "#CCBB44", "#EE6677", "#AA3377", "#000000"), 
     lty=1, lwd = 4, xlim=c(-5,30), xlab="Time (days)", ylab="Proportion Survived")
legend("bottomleft", c("Neg Control", "Acetone", "Amitraz", "Chlorothalonil", "Clothianidin", "Coumaphos", "Imidacloprid", "TauFluvalinate", "Dimethoate"), 
       col=c("#BBBBBB","#66CCEE", "#4477AA", "#228833","#999933", "#CCBB44", "#EE6677", "#AA3377", "#000000"), lty=1, lwd = 3, title="Pesticide")
###

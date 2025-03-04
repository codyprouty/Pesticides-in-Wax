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
library(survival)
library(survminer)
library(coxme)
library(car)
library(visreg)

##

#Load datasheet
Wax <- read.csv("PesticidesinWax.csv")
Wax$Count = rep(1, length(Wax$Replicate))
Wax$Replicate = as.factor(Wax$Replicate)
Aggregated = aggregate(Count ~ Replicate, data = Wax, sum)
Wax = merge(Wax, Aggregated, by = "Replicate")
names(Wax)[names(Wax) == "Count.x"] <- "Count"
names(Wax)[names(Wax) == "Count.y"] <- "Total"

#Increase count to be equal to total, count then becomes bee number
#Add a row for every day less than death day for a given bee while i < length of datasheet.
i = 1
while(i < 2717){
  if(i > 1 && Wax$Count[i] <= Wax$Count[i-1] && Wax$Count[i-1] < Wax$Total[i-1]){
    Wax$Count[i] = Wax$Count[i-1] + Wax$Count[i]
    i=i+1
    }
  else{
    i=i+1
  }
}
names(Wax)[names(Wax) == "Count"] <- "BeeNumber"

#For every individual bee, we need to show all the times they were alive. 
#Add a row for each death day where the bee was alive.

i = 1
while(i < 2717){
  if(Wax$DeathDay[i] > 1){
    j = 1
    while(j < Wax$DeathDay[i]){
      Row <- c(Wax$Replicate[i], Wax$Treatment[i], j, 1, Wax$BeeNumber[i], Wax$Total[i])
      Wax <- rbind(Wax, Row)
      j=j+1
    }
    i=i+1
  }
  else{
    i=i+1
  }
}

Wax$DeathDay = as.numeric(Wax$DeathDay)
Wax$Survived = as.numeric(Wax$Survived)

Wax10 = Wax

##Only check survival of first 10 days
Wax10$Survived[which(Wax10$DeathDay > 10)] <- 1
##

#Data organization
Wax10$SurvObj <- with(Wax10, Surv(DeathDay, Survived == 0))
Wax$SurvObj <- with(Wax, Surv(DeathDay, Survived == 0))
Survival <- survfit(SurvObj ~ Treatment, data = Wax10)

#Subsetted for the analysis of day of death
Dead <- subset(Wax, Survived == 0)
##

#Effect of treatment on death day
DD <- lmer(DeathDay ~ Treatment + (1|Replicate), data= Dead)
anova(DD)

resids <- resid(DD)
hist(resids, breaks = 10, xlab="residuals")
shapiro.test(resids)

Dead$DeathDay.t <- Dead$DeathDay + 1
qqp(Dead$DeathDay, "norm")

poisson <- fitdistr(Dead$DeathDay, "Poisson")
qqp(Dead$DeathDay, "pois", lambda=poisson$estimate)

nbinom <- fitdistr(Dead$DeathDay, "Negative Binomial")
qqp(Dead$DeathDay, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

gamma <- fitdistr(Dead$DeathDay, "gamma")
qqp(Dead$DeathDay, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#Not normally distributed, transformations ineffective

m1 <- coxme(Surv(DeathDay, Survived) ~ Treatment + (1|Replicate), data = Wax3)


DD <- kruskal.test(DeathDay ~ Treatment, data= Dead)
DD

pairwise.wilcox.test(Dead$DeathDay, Dead$Treatment,
                     p.adjust.method = "BH")
###

#Effect of treatment on survival
Sur <- glm(Survived ~ Treatment, 
          family = binomial, Wax, method=brglmFit)
anova(Sur, test = "Chisq")
summary(Sur)

Wax1 <- subset(Wax, Treatment != "zDimethoate")
Wax2 <- subset(Wax1, Treatment != "1NegControl")
Wax3 <- subset(Wax2, Treatment != "Acetone")

Wax11 = subset(Wax10, Treatment != "zDimethoate")

m1 <- coxme(Surv(DeathDay, Survived) ~ Treatment + (1|Replicate), data = Wax1)
anova(m1)
summary(m1)

res <- pairwise_survdiff(Surv(DeathDay, Survived) ~ Treatment, data = Wax3, p.adjust.method = "bonferroni")
res

Surv <- glmer(Survived ~ Treatment + (1|Replicate), family = binomial(link="logit"), Wax)
Anova(Surv, test="Chisq")

lsm<-lsmeans (m1, list( ~ Treatment))
cld(lsm)


m1 <- coxme(Surv(DeathDay, Survived) ~ Treatment + (1|Replicate), data = Wax11)
anova(m1)
summary(m1)
lsm<-lsmeans (m1, list( ~ Treatment))
cld(lsm)
plot(lsm)

m1 <- coxph(Surv(DeathDay, Survived) ~ Treatment, data = Wax)
m1
anova(m1)
summary(m1)
lsm<-lsmeans (m1, list( ~ Treatment))
cld(lsm)

res <- pairwise_survdiff(Surv(DeathDay, Survived) ~ Treatment, data = Wax10, p.adjust.method = "bonferroni")
res

new_df <- with(Wax,
               data.frame(Treatment = unique(Wax$Treatment))
)

fit <- survfit(m1, newdata = new_df)

#png("SmallSurvival.png", width=8, height=6, units="in", res=400)
ggsurvplot(fit, conf.int = TRUE, palette = "Dark3", 
           censor = FALSE, surv.median.line = "hv", data=new_df$Treatment)


library(ecotox)

results <- LT_probit(survfit(Surv(DeathDay, Survived == 0) ~ Treatment, data = Wax),
                     p = c(50,99),
                     data = Wax)
results




Wax10Imi = subset(Wax10, Treatment == "Imidacloprid")
km_fit <- survfit(Surv(DeathDay, Survived == 0) ~ Treatment, data = Wax)
plot(km_fit, col=c("#BBBBBB","#66CCEE", "#4477AA", "#228833","#999933", "#CCBB44", "#EE6677", "#AA3377", "#000000"), 
     lty=1, lwd = 4, xlim=c(-5,30), xlab="Time (days)", ylab="Proportion Survived")
summary(km_fit)
###

#Graph of survival across treatments
Survival <- survfit(SurvObj ~ Treatment, data = Wax)
plot(Survival, col=c("#BBBBBB","#66CCEE", "#4477AA", "#228833","#999933", "#CCBB44", "#EE6677", "#AA3377", "#000000"), 
     lty=1, lwd = 4, xlim=c(-5,30), xlab="Time (days)", ylab="Proportion Survived")
legend("bottomleft", c("Neg Control", "Acetone", "Amitraz", "Chlorothalonil", "Clothianidin", "Coumaphos", "Imidacloprid", "TauFluvalinate", "Dimethoate"), 
       col=c("#BBBBBB","#66CCEE", "#4477AA", "#228833","#999933", "#CCBB44", "#EE6677", "#AA3377", "#000000"), lty=1, lwd = 3, title="Pesticide")
###

tgcWax <- summarySE(Wax, measurevar="Survived", groupvars=c("Treatment"))
tgcWax10 <- summarySE(Wax10, measurevar="Survived", groupvars=c("Treatment"))
write.csv(tgcWax10, "waxTable.csv")

###
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


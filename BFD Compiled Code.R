# Full BFD Analyses 
# Updates 8/21/2020

library("lme4")
library("car")
library("multcomp")
library("ggplot2")
library("emmeans")
library("piecewiseSEM")
library("glmmTMB")
library("nlme")
library("lmerTest")
library("tidyverse")
library("rcompanion")
library("blmeco")
library("ggpubr")
library("plotrix")
library("epiDisplay")


#1. SEM #### 
# Load and organize Data
setwd("C:/Users/bwl42/Desktop/sem") 
quaddat= read.csv("BFDFinal2.csv", header=T, sep=",")
quaddat <- quaddat %>% mutate(uninf = (25-totinf),  # Just adding some new variables
                              inf.ratio = (totinf/25),
                              w.six = (25),addist = (adddist),
                              inf.ratio.off = ((totinf-1)/24),
                              w.off = (24), 
                              avgdistinf.off = ((avgdistinf *25)/24), totpop = (totad + totnym), 
                              logpop = (log(totpop)),
                              block = if_else(Block == 1 , "one", "two"), lognym = (log(totnym)),
                              rand=as.factor(paste(Block,Rep)),propsource = (sourcetot/totpop))
sem_subset = quaddat %>% 
  dplyr::select(lethal,risk, Block, block, totad,totnym,addist,nymdist,inf.ratio,w.six)
head(sem_subset)

quaddat %>% 
  group_by(Treat) %>% 
  dplyr::summarize(mean_inf = mean(inf.ratio), SE_inf = std.error(inf.ratio))

# a. Trend Verification and Visualization #### 
# Titer
titercheck <- glm(titer ~ inf.ratio,  data = quaddat)
cor.test(quaddat$titer,quaddat$inf.ratio)
ggplot(quaddat, mapping = aes(x = totinf, y = titer)) + geom_point(aes(color = as.factor(Block)))+
  geom_smooth(method = "glm") + labs(title = "Visual Identification correlated with Quadrant Titer") 
          # Trend validates our technique, though block effect complicates interpretation
# Univariate Plots
codes <- list("C" = "Control", "R" = "Risk", "L" = "Lethal")
plotdat <- quaddat %>% 
  mutate(Treat = dplyr::recode(Treat, !!!codes)) 
plotdat$Treat <- factor(plotdat$Treat, levels=c("Control", "Risk", "Lethal"))
a= ggplot(plotdat, mapping = aes(x = Treat, y = inf.ratio)) + geom_boxplot() + xlab("")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
b= ggplot(plotdat, mapping = aes(x = Treat, y = totnym)) + geom_boxplot() + xlab("")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
c= ggplot(plotdat, mapping = aes(x = Treat, y = totad)) + geom_boxplot() + xlab("")
d= ggplot(plotdat, mapping = aes(x = Treat, y = addist)) + geom_boxplot() + xlab("")
e= ggplot(plotdat, mapping = aes(x = totad, y = inf.ratio, weight = 24)) + geom_smooth(method = "glm",method.args = list(family = "binomial")) + geom_point(aes(color = Treat))
f= ggplot(plotdat, mapping = aes(x = totnym, y = inf.ratio, weight = 24)) + geom_smooth(method = "glm",method.args = list(family = "binomial"))+ geom_point(aes(color = Treat))
g= ggplot(plotdat, mapping = aes(x = addist, y = inf.ratio, weight = 24)) + geom_smooth(method = "glm",method.args = list(family = "binomial"))+ geom_point(aes(color = Treat))
h= ggplot(plotdat, mapping = aes(x = nymdist, y = inf.ratio, weight = 24)) + geom_smooth(method = "glm",method.args = list(family = "binomial"))+ geom_point(aes(color = Treat))

ggarrange(a,b,c,d) # Some validation, some big questions about validity of these trends
ggarrange(e,f,g,h) # These are wildly oversimplified vs the actual component models, but not great

# Univariate Model for OVerall Infection
head(quaddat)
first.mod <- glm(inf.ratio ~ Treat + Block, family = "binomial", weights = w.six, data = quaddat)
plotNormalHistogram(quaddat$inf.ratio)
summary(first.mod)
Anova(first.mod)
cld(emmeans(first.mod, ~ Treat))
residualPlot(first.mod)
with(summary(first.mod), 1 - deviance/null.deviance)

# Few extras
plotNormalHistogram(quaddat$totnym)
p2nym.mod <- glm(totnym ~ Treat + Block, data = quaddat)
summary(p2nym.mod)
nym2inf.mod <- glm(inf.ratio ~ totnym + Block, data = quaddat)
summary(nym2inf.mod)

# Likelihood ratio calculation 
quaddat %>% 
  dplyr::summarize(sum_inf = sum(totinf), sum_uninf = sum(uninf))
# INF = 520, UNINF = 680
520/680
# Total prevalence is 0.764
quaddat %>% 
  dplyr::group_by(Treat) %>% 
  dplyr::summarize(sum_inf = sum(totinf), sum_uninf = sum(uninf),
                   prev = sum_inf/(sum_inf+sum_uninf))

# Prevalence ratios
# L to C
0.445/0.39
# R to C 
0.465/0.39
# Overall 
((178+186)/(178+186+222+214))/0.39

# These essentially just confirm the individual pathway effects.

# b. Create SEM ####
# run 1 #
ad.pop.1 <- lm(totad ~ lethal + risk + block, data = quaddat)
nym.pop.1 <- lm(totnym ~ lethal + risk + block, data = quaddat)
ad.dist.1 <- lm(addist ~ lethal + risk + block, data = quaddat)
nym.dist.1 <- glm(nymdist ~ lethal + risk + block, data = quaddat)
inf.rat.1 <- glm(inf.ratio ~ totad + totnym + addist + nymdist+ block, family="binomial", weights=w.six, data=quaddat)
sem.1 <- psem(ad.pop.1, nym.pop.1, ad.dist.1,nym.dist.1,inf.rat.1, data = quaddat)
summary(sem.1, standardize = "none", conserve = TRUE)
summary(update(sem.1, totad %~~% totnym, addist %~~% nymdist), standardize="scale")
# run 2 final #
ad.pop.1.1 <- lm(totad ~ lethal , data = quaddat)
nym.pop.1.1 <- lm(totnym ~ lethal + risk+ block, data = quaddat)
ad.dist.1.1 <- lm(addist ~ lethal + risk + totad + totnym, data = quaddat)
nym.dist.1.1 <- lm(nymdist ~ lethal , data = quaddat)
inf.rat.1.1 <- glm(inf.ratio ~ totnym  + nymdist+ block + lethal, family="binomial", weights=w.six, data=quaddat)
sem.1.1 <- psem(ad.pop.1.1, nym.pop.1.1, ad.dist.1.1, nym.dist.1.1, inf.rat.1.1, data = quaddat)
summary(sem.1.1, standardize = "scale", conserve = TRUE)
summary(update(sem.1.1, totnym %~~% totad, nymdist %~~% addist), standardize="scale") # Final Model

# For exporting into table
#coeftable1 <- as.data.frame(coefs(sem.1.1))
#View(coeftable1)
#write.csv(coeftable, "coeftable.csv")
#write.csv(dSep(sem.1.1), "dSeptable.csv")

# 2. 5P Analyses #### 
# Load data 
setwd("C:/Users/bwl42/Desktop/Projects/Big Field Dispersal/5p")
fivetime <- read.csv("5ptimenew.csv", sep=",", header=TRUE) # Additional non-time series data available too
fivetime <- fivetime %>% 
  filter(Treat != "LL", Treat != "LR") %>% 
  mutate(pop = (nymtot + adtot)) %>% mutate(logpop = log(pop))

fivesummary <- fivetime %>% 
  group_by(Treat) %>% 
  dplyr::summarize(mean_adtop = mean(adtop), SE_adtop = std.error(adtop),
                   mean_nymtop = mean(nymtop), SE_nymtop = std.error(nymtop),
                   mean_adtot = mean(adtot), SE_adtot = std.error(adtot),
                   mean_nymtot = mean(nymtot), SE_nymtot = std.error(nymtot))
fivesummary

# a. GLMs ####  

# Adult Position
plotNormalDensity(fivetime$adtop)
ad.top <- glmer(adtop ~ Treat*Time + (1|Cage), family = "binomial", data= fivetime, weight = adtot)   
dispersion_glmer(ad.top)
summary(ad.top)
Anova(ad.top) 
cld(emmeans(ad.top, ~ Treat*Time, type= "response", adjust = "none")) 
View(as.data.frame(Anova(ad.top)))

# Nymph Position
nym.top <- glmer(nymtop ~ Treat*Time +  (1|Cage) , family= "binomial", data= fivetime, weight = nymtot)   
dispersion_glmer(nym.top)
summary(nym.top)
Anova(nym.top)
cld(emmeans(nym.top, ~ Treat*Time, type= "response", adjust = "none")) 

# Adult Population
adultpop <- lmer(adtot ~ Treat*Time + (1|Cage), data= fivetime)   
overdispersion.test(adultpop)
summary(adultpop)
predict(adultpop,data.frame(Treat="HL",Time=2, Cage = "HL1" ))
Anova(adultpop, type = 3) 
cld(emmeans(adultpop, ~ Treat*Time))


# Nymph Population
nymphpop <- lmer(nymtot ~ Treat*Time + (1|Cage), data= fivetime)  
summary(nymphpop)
Anova(nymphpop, type = 3, test.statistic = "F") 
cld(emmeans(nymphpop, ~Treat*Time))
View(fivetime %>% dplyr::select(Treat,Time,nymtot))


# Adult Movement
adultdist <- lmer(addist ~ Treat*Time + (1|Cage),  data= fivetime)   
summary(adultdist)
Anova(adultdist) 
cld(emmeans(adultdist, ~Treat))

# Nymph Movement
nymphdist <- lmer(nymdist ~ Treat*Time + (1|Cage),  data= fivetime)   
summary(nymphdist)
Anova(nymphdist) 
cld(emmeans(nymphdist, ~Treat))


# b. Position Boxplots ####
codes <- list("C" = "Control", "HR" = "Risk", "HL" = "Lethal")
BPdat <- fivetime %>% 
  mutate(Treat = dplyr::recode(Treat, !!!codes)) 
BPdat$Treat <- factor(BPdat$Treat, levels=c("Control", "Risk", "Lethal"))

adultBP <- ggplot(BPdat, aes(x = Treat, y= prop)) + 
  geom_boxplot(mapping = aes(x = Treat, y= adtop*100), size = 1,  position=position_dodge(), stat = "boxplot") +         
  theme_classic() + ylab("% Aphid Adults on Top Half of Host") +
  geom_hline(aes(yintercept= 50), linetype = "dotted") + xlab("") + 
  theme(axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12))

nymphBP <- ggplot(BPdat, aes(x = Treat, y= prop)) + 
  geom_boxplot(mapping = aes(x = Treat, y= nymtop*100), size = 1,  position=position_dodge(), stat = "boxplot") +         
  theme_classic()  + ylab("% Aphid Nymphs on Top Half of Host") + 
  geom_hline(aes(yintercept= 50), linetype = "dotted") + xlab("")+ 
  theme(axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12))

ggarrange(adultBP,nymphBP, nrow=2) #figure size is 450x600

# c. Population and Movement Figures ####
head(fivetime)
codes <- list("C" = "Control", "HR" = "Risk", "HL" = "Lethal")
popmovesummary <- BPdat %>% 
  group_by(Treat, Time) %>% 
  dplyr::summarize(meanadpop = mean(adtot), SEadpop = std.error(adtot),
                   meannympop = mean(nymtot), SEnympop = std.error(nymtot),
                   meanaddist = mean(addist), SEaddist = std.error(addist),
                   meannymdist = mean(nymdist), SEnymdist = std.error(nymdist),
                   meanadtop = mean(adtop), SEadtop = std.error(adtop),
                   meannymtop = mean(nymtop), SEnymtop = std.error(nymtop))
popmovesummary
pd = position_dodge(.1)

# adult abundance
adultPopfig <- ggplot(popmovesummary, aes(x=Time, y=meanadpop, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meanadpop-SEadpop, ymax=meanadpop+SEadpop), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Aphid Adult Abundance")
adultPopfig

# nymph abundance
nymphPopfig <- ggplot(popmovesummary, aes(x=Time, y=meannympop, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meannympop-SEnympop, ymax=meannympop+SEnympop), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Aphid Nymph Abundance")
nymphPopfig


# adult Distance From Center
adultdistfig <- ggplot(popmovesummary, aes(x=Time, y=meanaddist, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meanaddist-SEaddist, ymax=meanaddist+SEaddist), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Aphid Adult Dist. From Center")
adultdistfig

# nymph Distance From Center
nymphdistfig <- ggplot(popmovesummary, aes(x=Time, y=meannymdist, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meannymdist-SEnymdist, ymax=meannymdist+SEnymdist), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Aphid Nymph Dist. From Center")
nymphdistfig

ggarrange(adultPopfig,nymphPopfig,
          adultdistfig,nymphdistfig, ncol = 2,nrow=2)

# adult abundance
adulttopfig <- ggplot(popmovesummary, aes(x=Time, y=meanadtop, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meanadtop-SEadtop, ymax=meanadtop+SEadtop), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Prop. Adults Feeding Top Half")
adulttopfig

# nymph abundance
nymphtopfig <- ggplot(popmovesummary, aes(x=Time, y=meannymtop, colour=Treat)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 4, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=meannymtop-SEnymtop, ymax=meannymtop+SEnymtop), width=.1, position = pd) +
  theme_classic() + 
  theme(axis.text= element_text(size = 12), axis.title = element_text(size = 12),
        legend.position = "bottom")+
  labs(x = "Day", y = "Prop. Nymphs Feeding Top Half")
nymphtopfig

ggplot(data= fivetime %>% filter(Treat=="HR"), mapping = aes(x = Time, y = nymtot, color = Treat)) + 
  geom_point() + geom_smooth( method = "glm")+
  facet_wrap(~Cage)

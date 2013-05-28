library(glmmLDTS)
#library(lme4)
setwd("/media/Hitachi2GB/00NMML/glmp.ts/glmmLDTS")
load("Test5Seals.RData")


#-------------------------------------------------------------------------------
# FIT THE MODEL
#-------------------------------------------------------------------------------

Vhour <- c("sixAMness", "sixAMness2", "noonness", "noonness2")
Vdate <- c("sprness", "sprness2", "sumness", "sumness2")

Vform <- c("HOStatus ~ HrFLo + HrFLo2 ")
for(i in 1:length(Vhour)) Vform <- paste(as.character(Vform), "+", Vhour[i])
for(i in 1:length(Vdate)) Vform <- paste(as.character(Vform), "+", Vdate[i])
for(i in 1:length(Vhour)) {
	for(j in 1:length(Vdate)) {
		Vform <- paste(as.character(Vform), " + ", Vhour[i], ":", 
			Vdate[j], sep = "")
	}
}
Vform <- as.formula(Vform)

#undebug(glmmLDTS1)
all.fit.1 <- glmmLDTS1(Vform,
	random = HOStatus ~ SpeNo,
	data = Test5SealsData,
	timecol = "time.vec", group.vec = "SpeNo")
all.fit.1$formula
all.fit.1$random.terms
all.fit.1$timecol
all.fit.1$link.func
all.fit.1$trialscol
all.fit.1$group.vec
all.fit.1$ridge.reg
all.fit.1$lambda
all.fit.1$start.time
all.fit.1$end.time
all.fit.1$R.cov.parameters
all.fit.1$G.cov.parameters
all.fit.1$typeIII.hypoth
all.fit.1$fixed.effects
all.fit.1$random.effects
all.fit.1$outer.iterations
all.fit.1$inner.iterations2

attach(Test5SealsData)
glmfit <- glm(Vform,
    data = Test5SealsData,
    family = binomial)



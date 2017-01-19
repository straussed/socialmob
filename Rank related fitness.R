##################################################################
##                        Eli Strauss                           ##
##                  Rank related fitness                        ##
##                    January 5th, 2017                         ##
##################################################################
library(lme4)
library(ggplot2)
library(MASS)
library(glmmADMB)
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")
ranks <- ranks[ranks$Year != 1988,]
###testing

ranks$fitness2yo <- NA
for(row in 1:length(ranks[,1])){
  if(ranks[row, 'Clan'] %in% c('north', 'south', 'hz') & ranks[row,'Year'] < 2010){next}
  year <- ranks[row,'Year']
  ranks[row,'fitness2yo'] <- length(filter(tblHyenas, Mom == ranks[row,'ID'] & (is.na(Disappeared) | Disappeared - Birthdate > 730) & format(Birthdate + 730, '%Y') == year)[,1])
}


ranks$stan.rank.model <- ranks$stan.rank+1
ranks$is.alpha <- ifelse(ranks$Rank == 1, 1,0)

- 

m1 <- lm(data = ranks, fitness2yo ~ Rank)
m2 <- lm(data = ranks, fitness2yo ~ Rank + I(Rank^2))
m3 <- lmer(data = ranks, fitness2yo ~ stan.rank.model + I(stan.rank.model^2) + is.alpha + (1|ID))
m4 <- lmer(data = ranks, fitness2yo ~ stan.rank.model + is.alpha + (1|ID))
m5 <- lmer(data = ranks, fitness2yo ~ stan.rank.model + (1|ID))
m6 <- lm(data = ranks, fitness2yo ~ stan.rank.model + is.alpha)

p3 <- glmer(data = ranks, fitness2yo ~ stan.rank.model + I(stan.rank.model^2) + is.alpha + (1|ID), family = 'poisson')
p4 <- glmer(data = ranks, fitness2yo ~ stan.rank.model + is.alpha + (1|ID), family = 'poisson')
p5 <- glmer(data = ranks, fitness2yo ~ stan.rank.model + (1|ID), family = 'poisson')
p6 <- glm(data = ranks, fitness2yo ~ stan.rank.model + is.alpha, family = 'poisson')

#####bin by stan rank for plotting
numbins <- 20
binedge <- c(seq(0,1.999999, length.out = numbins-1), 2)
ranks$rank.bin <- cut(ranks$stan.rank.model, binedge)

fitRanksBin <- aggregate(x = ranks, by = list(ranks$rank.bin), FUN = mean, na.rm = T)[,c('Group.1', 'fitness2yo', 'stan.rank.model')]
ranks.ordered <- ranks[order(ranks$stan.rank.model),]  
 
####Model with quadratic 
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(ranks.ordered$stan.rank.model, exp(fixef(p3)[1] + fixef(p3)[2]*ranks.ordered$stan.rank.model + fixef(p3)[3]*(ranks.ordered$stan.rank.model^2) + fixef(p3)[4]*ranks.ordered$is.alpha), lwd = 3)

#####Model without quadratic
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(ranks.ordered$stan.rank.model, exp(fixef(p4)[1] + fixef(p4)[2]*ranks.ordered$stan.rank.model + fixef(p4)[3]*ranks.ordered$is.alpha), lwd = 3)


##################################Testosterone#########################
ttt <- read.csv("~/Documents/Fisibase/testosterone.csv")
ttt$hyenaID <- tolower(ttt$hyenaID)
ttt$poop_date <- as.Date(ttt$poop_date, format = '%d-%b-%y')
ttt$poop_year <- as.numeric(format(ttt$poop_date, '%Y'))
ttt$am.pm <- as.factor(ttt$am.pm)
ttt$ng.g <- as.numeric(ttt$ng.g)
ttt[ttt$poop_time == '','am.pm'] <- NA

ttt.am <- filter(ttt, am.pm == 'AM')
ttt.am.mean <- aggregate(ng.g ~ poop_year + hyenaID, data = ttt.am, mean)

ranks$fecalT <- left_join(x = ranks, y = ttt.am.mean, by = c("Year" = "poop_year", "ID" = "hyenaID"))$ng.g
boxplot(ranks$fecalT ~ ranks$DiffBinary)
#######################################################################

#######################################################################

confints <- confint.merMod(m3)[3:5,]

fitdif <- function(start, end){
  est <- exp(fixef(p4)[1] + fixef(p4)[2]*end + fixef(p4)[3]*(ifelse(end == 2, 1,0))) - exp(fixef(p4)[1] + fixef(p4)[2]*start + fixef(p4)[3]*(ifelse(start == 2, 1,0)))
  #est <- (fixef(p5)[1] + fixef(p5)[2]*end) - (fixef(p5)[1] + fixef(p5)[2]*start)
  #upbound <- (confints[1,1] + confints[2,1]*end + confints[3,1]*(end^2)) - (confints[1,1] + confints[2,1]*start + confints[3,1]*(start^2))
  #lowbound <- (confints[1,2] + confints[2,2]*end + confints[3,2]*(end^2)) - (confints[1,2] + confints[2,2]*start + confints[3,2]*(start^2))
  #return(c(est, upbound, lowbound))
  return(est)
}

ranks$Diff <- NA
ranks$Move <- "None"
for(row in 1:length(ranks[,1])){
  if(ranks[row,'ID']==ranks[row,'IDold']){
    start = end
    next
  }
  end <- ranks[row,]$stan.rank.model
  start <- ranks[ranks$IDold == ranks[row,'ID'] & ranks$Year == ranks[row,'Year'], 'stan.rank.model']
  ranks[row,'Diff'] <- fitdif(start, end)
  if(start < end){
    ranks[row,'Move'] <- "Up"
  }else if(start > end){ranks[row,'Move'] <- "Down"}
}
ranks$WhichHalf <- ifelse(ranks$stan.rank < 0, 0, 1)
ranks$DiffBinary  <- ifelse(ranks$Diff > 0, 1, 0)
ggplot(data = ranks, aes(x = Diff, fill  = as.factor(WhichHalf))) +
  geom_bar(position = 'dodge', binwidth = .1)

########################Add coalition info##################################
for(row in 1:length(ranks[,1])){
  id <- ranks[row,'ID']
  year <- ranks[row,'Year']
  clan <- ranks[row,'Clan']
  ranks[row,'coals'] <- length(filter(aggsFull,
                         Agg == id,
                         Clan == clan,
                         Year == year,
                         Seq > 0)[,1])
  ranks[row,'obsTime'] <- sum(sessions[grep(id, filter(sessions, format(Date, '%Y') == year)$Hyenas),'Time'], na.rm = T)
  ranks[row,'coalPart'] <- length(unique(unlist(strsplit(filter(aggsFull,
                                            Agg == id,
                                            Clan == clan,
                                            Year == year,
                                            Seq > 0)[,'Group'], ','))))-2
  if(ranks[row,'coals'] == 0){ranks[row, 'coalPart'] <- 0}
  ranks[row, 'coalRate'] <- ranks[row,'coals']/ranks[row,'obsTime']
  ranks[row,'aggs'] <- length(filter(aggsFull,
                           Agg == id,
                           Clan == clan,
                           Year == year)[,1])
  ranks[row,'aggRate'] <- ranks[row,'aggs']/ranks[row,'obsTime']
  
}
ranks$Move <- factor(ranks$Move, levels = c('None', 'Up', 'Down'))
############################################################################

boxplot(ranks$coalRate ~ ranks$Move)
boxplot(ranks$coalPart ~ ranks$Move)
plot(data = filter(ranks, Move == 'None'), coalPart ~ Rank)
ggplot(data = filter(ranks, obsTime != 0), aes(y = coalPart, x = Move, col = Rank)) + 
  geom_point()


################################Modeling####################################
####coalition partners
ranks$Year <- as.factor(ranks$Year)
ranks$Clan <- as.factor(ranks$Clan)
cpm0 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank)
cpm1 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank, offset = log(obsTime))
cpm2 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank + Move, offset = log(obsTime))
cpm3 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank + Move + (1|Year), offset = log(obsTime))
cpm4 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank + Move + (1|Year) + (1|Clan), offset = log(obsTime))
cpm5 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coalPart ~ Rank + Move + Move*Rank + (1|Year) + (1|Clan), offset = log(obsTime))

cpmnb0 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coalPart ~ Rank)
cpmnb1 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coalPart ~ Rank + offset(log(obsTime)))
cpmnb2 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coalPart ~ Rank + Move + offset(log(obsTime)))
cpmnb3 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coalPart ~ Rank + (1|Clan/Year) + offset(log(obsTime)))
cpmnb4 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coalPart ~ Rank + Move + (1|Clan) + offset(log(obsTime)))
cpmnb5 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coalPart ~ Rank + Move + (1|Clan/Year) + offset(log(obsTime)))
cpmnb6 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coalPart ~ Rank * Move + (1|Clan/Year) + offset(log(obsTime)))

AIC(cpm0, cpm1, cpm2, cpm3, cpm4, cpm5, cpmnb0, cpmnb1,cpmnb2,cpmnb3, cpmnb4, cpmnb5, cpmnb6)

####coalition rate
cm0 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank)
cm1 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank, offset = log(obsTime))
cm2 <- glm(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank + Move, offset = log(obsTime))
cm3 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank + Move + (1|Year), offset = log(obsTime))
cm4 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank + Move + (1|Year) + (1|Clan), offset = log(obsTime))
cm5 <- glmer(data = filter(ranks, obsTime != 0), family = poisson, formula = coals ~ Rank + Move + Move*Rank + (1|Year) + (1|Clan), offset = log(obsTime))

cmnb0 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coals ~ Rank)
cmnb1 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coals ~ Rank + offset(log(obsTime)))
cmnb2 <- glm.nb(data = filter(ranks, obsTime != 0), formula = coals ~ Rank + Move + offset(log(obsTime)))
cmnb3 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coals ~ Rank + (1|Clan/Year) + offset(log(obsTime)))
cmnb4 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coals ~ Rank + Move + (1|Clan) + offset(log(obsTime)))
cmnb5 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coals ~ Rank + Move + (1|Clan/Year) + offset(log(obsTime)))
cmnb6 <- glmmadmb(data = filter(ranks, obsTime != 0), family = "nbinom", formula = coals ~ Rank * Move + (1|Clan/Year) + offset(log(obsTime)))

AIC(cm0, cm1, cm2, cm3, cm4, cm5, cmnb0, cmnb1,cmnb2,cmnb3, cmnb4, cmnb5, cmnb6)


############################################################################

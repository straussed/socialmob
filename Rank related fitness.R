##################################################################
##                        Eli Strauss                           ##
##                  Rank related fitness                        ##
##                    January 5th, 2017                         ##
##################################################################
library(lme4)
library(ggplot2)
library(MASS)
library(glmmADMB)
library(rethinking)
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")
source("~/Documents/Research/RCIISI/minicrank/IdentifyRankChanges.R")
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


ranksDecay <- data.frame()

for(year in unique(ranks$Year)){
  yrranks <- ranks[ranks$Year == year,]
  for(id in unique(yrranks$ID)){
    if(id %in% ranks[ranks$Year == (year+1), 'ID']){
      ranksDecay <- rbind(ranksDecay, cbind(Year = year, ID = id, Rank = round(yrranks[yrranks$ID == id, 'stan.rank'], 4),
                                            Delta = round(ranks[ranks$Year == (year + 1) & ranks$ID == id, 'stan.rank'] - yrranks[yrranks$ID == id, 'stan.rank'], 4),
                                            absRank = yrranks[yrranks$ID == id, 'Rank'],
                                            absDelta = yrranks[yrranks$ID == id, 'Rank'] - ranks[ranks$Year == (year + 1) & ranks$ID == id, 'Rank'],
                                            Clan = unique(ranks[ranks$ID == id,'Clan'])))
    }
  }
}
ranksDecay$Rank <- as.numeric(ranksDecay$Rank)
ranksDecay$Delta <- as.numeric(ranksDecay$Delta)
ranksDecay$absRank <- as.numeric(ranksDecay$absRank)
ranksDecay$absDelta <- as.numeric(ranksDecay$absDelta)

###Integrate rank.changes object in order to remove individuals who change ranks
for(row in 1:length(rank.changes[,1])){
  ranksDecay <- filter(ranksDecay, Year != (rank.changes[row,'Year']-1) | (ID != rank.changes[row,'Upmover'] & ID !=rank.changes[row,'Downmover']))
}

hist(ranksDecay$absDelta)


# absM <- lm(data = ranksDecay, absDelta ~ I(poly(absRank, 2)))
# 
# plot(ranksDecay$absDelta ~ ranksDecay$absRank)
# rseq <- seq(0,60)
# rseq.s <- seq(0:60)^2
# lines(coef(absM)[1] + coef(absM)[2]*rseq + coef(absM)[3]*rseq.s)

#ranksDecay$absRank.s <- ranksDecay$absRank ^ 2
#polynomial seems inappropriate
rnkfit <- map(alist(
            Delta ~ dnorm(u, s),
            u <- a + b1*Rank,
            a ~ dnorm(0,10),
            b1 ~ dnorm(0,10),
            s ~ dunif(0,20)
            ), data = ranksDecay)
plot(precis(rnkfit))

plot(ranksDecay$absDelta ~ jitter(ranksDecay$absRank))
with(filter(ranksDecay, Clan == 'talek'), plot(Delta ~ jitter(Rank), main = 'Talek'))
with(filter(ranksDecay, Clan != 'talek'), plot(Delta ~ jitter(Rank), main = 'Not Talek'))
abline(coef(rnkfit)[1], coef(rnkfit)[2])

alphaDescent <- function(id){
  if(length(id) == 0 ){
    return(F)
  }else if(is.na(id) | id == ''){
    return(F)
  }else if(id == 'kb'){
    return(T)
  }else(alphaDescent(tblHyenas[tblHyenas$ID == id,'Mom']))
}


ranks$Move <- "None"
for(row in 1:length(ranks[,1])){
  ranks[row,'clanSize'] = length(filter(ranks, Year == ranks[row,'Year'], Clan == ranks[row,'Clan'])[,1])
  if(ranks[row,'ID']==ranks[row,'IDold']){
    start = end
    next
  }
  end <- ranks[row,]$stan.rank.model
  start <- ranks[ranks$IDold == ranks[row,'ID'] & ranks$Year == ranks[row,'Year'], 'stan.rank.model']
  if(start < end){
    ranks[row,'Move'] <- "Up"
  }else if(start > end){ranks[row,'Move'] <- "Down"}
}


#######plot lines for each individual
#ranks$ID <- ranks$NewOrder
par(fig = c(0,1,0,1))
plot(data = ranks, Rank ~ Year, type = 'n', ylim = c(50,0), main = 'Talek')
for(id in unique(ranks[ranks$Clan == 'talek',]$ID)){
  #clr <- ifelse(alphaDescent(id), 'blue', 'black')
  #if(id == 'nav'){clr <- 'red'}
  clr <- 'black'
  if('Up' %in% ranks[ranks$ID == id,'Move'] | 'Down' %in% ranks[ranks$ID == id, 'Move']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = clr))}
}

par(fig = c(0.07,(.07+(6/25)+.02),0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'south',], Rank ~ Year, type = 'n', ylim = c(20,0), main = 'South')
for(id in unique(ranks[ranks$Clan == 'south',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'Move'] | 'Down' %in% ranks[ranks$ID == id, 'Move']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}

par(fig = c(0,1,0,1), new = T)
par(fig = c((.07+(6/25)+.02),(.33+.26), 0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'north',], Rank ~ Year, type = 'n', ylim = c(20,0), main = 'North')
for(id in unique(ranks[ranks$Clan == 'north',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'Move'] | 'Down' %in% ranks[ranks$ID == id, 'Move']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}

par(fig = c(0,1,0,1), new = T)
par(fig = c(0.59,.59+.26, 0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'hz',], Rank ~ Year, type = 'n', ylim = c(20,0), main = 'Happy Zebra', axes = F, frame.plot = T)
axis(1, at = seq(2009, 2012))
axis(2, at = c(0, 5, 10, 15, 20))
for(id in unique(ranks[ranks$Clan == 'hz',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'Move'] | 'Down' %in% ranks[ranks$ID == id, 'Move']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}


##################Compare fitness######################
latIds <- unique(filter(ranks, Year >= 2008, Year <= 2012, stan.rank < 0, Clan == 'talek')$ID)
earlyPer <- filter(ranks, ID %in% latIds & Clan == 'talek' & Year < 2008 & Year >= 2003)
latePer <- filter(ranks, ID %in% earlyPer$ID & Clan == 'talek' & Year >= 2008 & Year <= 2012)
dev.off()
par(mfrow = c(1,2))
hist(earlyPer$fitness2yo)
hist(latePer$fitness2yo)

#######Plot ranks by age#######
#ranks$Age <- NA
par(fig = c(0,1,0,1))
plot(data = ranks, Rank ~ Age, type = 'n', main = 'Talek', ylim = c(50,0))
for(id in unique(ranks[ranks$Clan == 'talek',]$ID)){
  irank <- ranks[ranks$ID == id & ranks$Year == min(ranks[ranks$ID == id,'Year']),'Rank']
  ranks[ranks$ID == id,'dRank'] <- irank - ranks[ranks$ID == id, 'Rank']
  if(!is.na(tblHyenas[tblHyenas$ID == id,'Birthdate'])){
    ranks[ranks$ID == id,'Age'] <- as.numeric(ranks[ranks$ID == id,'Year']) - as.numeric(format(tblHyenas[tblHyenas$ID == id,'Birthdate'], '%Y')) 
    }
  clr <- 'black'
  if(id %in% c(rank.changes$Upmover, rank.changes$Downmover)){
    with(filter(ranks, ID == id), lines(Rank ~ Age, lwd = 1.5, col = 'black'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Age, lwd = 1.5, col = 'black'))}
}

####delta rank 
dev.off()
par(fig = c(0,1,0,1))
plot(data = ranks, dRank ~ Age, type = 'n', axes = T, ylab= 'Net change in rank since onset of adulthood')
for(id in unique(ranks[ranks$Clan == 'talek',]$ID)){
  irank <- ranks[ranks$ID == id & ranks$Year == min(ranks[ranks$ID == id,'Year']),'Rank']
  ranks[ranks$ID == id,'dRank'] <- irank - ranks[ranks$ID == id, 'Rank']
  if(!is.na(tblHyenas[tblHyenas$ID == id,'Birthdate'])){
    ranks[ranks$ID == id,'Age'] <- as.numeric(ranks[ranks$ID == id,'Year']) - as.numeric(format(tblHyenas[tblHyenas$ID == id,'Birthdate'], '%Y')) 
  }
  clr <- 'black'
  if(id %in% c(rank.changes$Upmover, rank.changes$Downmover)){
    with(filter(ranks, ID == id), lines(dRank ~ Age, lwd = 1.5, col = 'black'))
  }else{with(filter(ranks, ID == id), lines(dRank ~ Age, lwd = 1.5, col = 'black'))}
}
abline(h = 0, col = 'red', lty = 2)
mRankDecay <- glm(data = ranks, exp(dRank) ~ Age)
plot(mRankDecay)

##################################################################
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
binedge <- c(seq(0,2.00001, length.out = numbins-1), 2)
ranks$rank.bin <- cut(ranks$stan.rank.model, binedge)

fitRanksBin <- aggregate(x = ranks, by = list(ranks$rank.bin), FUN = mean, na.rm = T)[,c('Group.1', 'fitness2yo', 'stan.rank.model')]
ranks.ordered <- ranks[order(ranks$stan.rank.model),]  

####Model with quadratic 
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(ranks.ordered$stan.rank.model, exp(fixef(p3)[1] + fixef(p3)[2]*ranks.ordered$stan.rank.model + fixef(p3)[3]*(ranks.ordered$stan.rank.model^2) + fixef(p3)[4]*ranks.ordered$is.alpha), lwd = 3)

#####Model without quadratic
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(ranks.ordered$stan.rank.model, exp(fixef(p4)[1] + fixef(p4)[2]*ranks.ordered$stan.rank.model + fixef(p4)[3]*ranks.ordered$is.alpha), lwd = 3)


###Variance
fitRanksBinVar <- aggregate(x = ranks, by = list(ranks$rank.bin), FUN = sd, na.rm = T)[,c('Group.1', 'fitness2yo', 'stan.rank.model')]
plot(fitRanksBinVar$stan.rank.model, fitRanksBinVar$fitness2yo, ylim = c(0,1), lwd = 3)


#####Model with rethinking
mreth <- glimmer(fitness2yo ~ (1|ID) + stan.rank.model + is.alpha, ranks[!is.na.data.frame(ranks$fitness2yo),], family = poisson)
mrs <- map2stan(mreth$f, mreth$d)
srm.seq <- seq(from = 0, to = 2, length.out = 1000)
mu <- link(mrs, data.frame(stan_rank_model = srm.seq, is_alpha = c(rep(0, 900), rep(0,100)), ID = as.factor(sample(ranks$ID, 1000, replace = T))))
mu.PI <- apply(mu, 2, PI)
mu.M <- apply(mu, 2, mean)
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(srm.seq, mu.M, col = col.alpha(rangi2, .4))
shade(mu.PI, srm.seq)##


#####Model with rethinking without alpha
mreth <- glimmer(fitness2yo ~ (1|ID) + stan.rank.model, ranks[!is.na.data.frame(ranks$fitness2yo),], family = poisson)
mrs <- map2stan(mreth$f, mreth$d)
srm.seq <- seq(from = 0, to = 2, length.out = 1000)
mu <- link(mrs, data.frame(stan_rank_model = srm.seq, rep(0,100)), ID = as.factor(sample(ranks$ID, 1000, replace = T)))
mu.PI <- apply(mu, 2, PI)
mu.M <- apply(mu, 2, mean)
plot(fitRanksBin$stan.rank.model, fitRanksBin$fitness2yo, ylim = c(0,1), lwd = 3)
lines(srm.seq, mu.M, col = col.alpha(rangi2, .4))
shade(mu.PI, srm.seq)##



##############################Observed Fitness########################
obsFit <- data.frame()
fitCutoff <- 3
rank.move <- filter(ranks, Move != 'None')
for(row in 1:length(rank.move[,1])){
  if(!length(ranks[ranks$Clan == rank.move[row,'Clan'] & ranks$Year == (rank.move[row,'Year']+fitCutoff),1])){next}
  #if(!rank.move[row,'ID'] %in% ranks[ranks$Year == (rank.move[row,'Year']+fitCutoff),'ID']){next}
  fit <- mean(filter(ranks, ID == rank.move[row,'ID'],
                     Year > rank.move[row, 'Year'],
                     Year <= rank.move[row,'Year']+fitCutoff
                     )$fitness2yo)
  longevity <- length(unique(ranks[rank.move[row,'ID'] %in% ranks$ID]$Year))
  obsFit <- rbind(obsFit, c(rank.move[row,c('Year', 'ID', 'Rank', 'stan.rank', 'dRank','Age', 'Move')], fit = fit))
}

boxplot(data = obsFit, fit ~ Move)
ggplot(data = obsFit, aes(y = fit, x = stan.rank, col = Move))+
  geom_point()+
  geom_smooth(method = lm)
######################################################################






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
par(mfrow = c(1,2))
boxplot(ranks$coalRate ~ ranks$Move, main = 'Rate of coalitions')
boxplot(ranks$coalPart ~ ranks$Move, main = 'Number of coalition partners')
plot(data = filter(ranks, Move == 'None'), coalPart ~ Rank)
ggplot(data = filter(ranks, obsTime != 0), aes(y = coalPart, x = Move, col = Rank)) + 
  geom_point()


summary(glm.nb(data = ranks, fitness2yo ~ stan.rank * clanSize))
summary(glm.nb(data = ranks, fitness2yo ~ stan.rank))
















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

r <- residuals(cmnb3)
summary(glm(exp(r) ~ filter(ranks, obsTime != 0)$Move, family = 'Gamma'))
############################################################################

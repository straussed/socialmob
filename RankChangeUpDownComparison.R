##################################################################
##                        Eli Strauss                           ##
##                Rank change modeling non-SNA                  ##
##                    January 2nd, 2017                         ##
##################################################################

library(dplyr)
library(sna)
library(network)
options(stringsAsFactors = FALSE)

source("~/Documents/Research/RCIISI/minicrank/IdentifyRankChanges.R")
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")
rank.changes.full <- rank.changes
ranks <- ranks[ranks$Year != 1988,]


##########################Add Kin################################
rank.changes$upKin <- NA
rank.changes$downKin <- NA
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    curRanks <- clanRanks[clanRanks$Year == year,]
    for(id in curRanks$ID){
      mom <- filter(tblHyenas, ID == id)$Mom
      sibs <- filter(tblHyenas, Mom == mom)$ID
      offspring <- filter(tblHyenas, Mom == id)$ID
      l <- curRanks[curRanks$ID %in% c(mom,sibs,offspring), 'ID']
      if(length(rank.changes[rank.changes$Year == year & rank.changes$Upmover == id,]$upKin)){
        rank.changes[rank.changes$Year == year & rank.changes$Upmover == id,]$upKin <- length(l)-1
      }
      if(length(rank.changes[rank.changes$Year == year & rank.changes$Downmover == id,]$downKin)){
        rank.changes[rank.changes$Year == year & rank.changes$Downmover == id,]$downKin <- length(l)-1
      }
    }
  }
}
################################################################



#######################Add coals###############################  
rank.changes$upCoal <- NA
rank.changes$downCoal <- NA
rank.changes$upNumCoal <- NA
rank.changes$downNumCoal <- NA
rank.changes$finalRank <- NA
rank.changes$initRank <- NA
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    curRanks <- clanRanks[clanRanks$Year == year,]
    fullCoals <- filter(aggsFull, Agg %in% curRanks$ID & Seq >= 1, as.numeric(Year) == year)
    upmovers <- unique(rank.changes[rank.changes$Year == year & rank.changes$Clan == clan,'Upmover'])
    downmovers <- unique(rank.changes[rank.changes$Year == year & rank.changes$Clan == clan,'Downmover'])
    if(length(c(upmovers,downmovers)) == 0){next}
    for(focal in upmovers){
      curCoals <- filter(fullCoals, Agg == focal)
      rank.changes[rank.changes$Upmover == focal & rank.changes$Year == year,'finalRank'] <- curRanks[curRanks$ID == focal,'Rank']
      rank.changes[rank.changes$Upmover == focal & rank.changes$Year == year,'upNumCoal'] <- length(curCoals[,1])
      if(length(curCoals[,1]) == 0){
        rank.changes[rank.changes$Upmover == focal & rank.changes$Year == year,'upCoal'] <- 0
        next
      }
      mates <- c()
      for(i in 1:length(curCoals[,1])){
        mates <- c(mates, strsplit(curCoals[i,'Group'], split = ',')[[1]][-1])
      }
      mates <- unique(mates)
      rank.changes[rank.changes$Upmover == focal & rank.changes$Year == year,'upCoal'] <- length(mates)-1
    }
    
    for(focal in downmovers){
      curCoals <- filter(fullCoals, Agg == focal)
      rank.changes[rank.changes$Downmover == focal & rank.changes$Year == year,'initRank'] <- curRanks[curRanks$IDold == focal,'Rank']
      rank.changes[rank.changes$Downmover == focal & rank.changes$Year == year,'downNumCoal'] <- length(curCoals[,1])
      if(length(curCoals[,1]) == 0){
        rank.changes[rank.changes$Downmover == focal & rank.changes$Year == year,'downCoal'] <- 0
        next
      }
      mates <- c()
      for(i in 1:length(curCoals[,1])){
        mates <- c(mates, strsplit(curCoals[i,'Group'], split = ',')[[1]][-1])
      }
      mates <- unique(mates)
      rank.changes[rank.changes$Downmover == focal & rank.changes$Year == year,'downCoal'] <- length(mates)-1
    }
  }
}
################################################################





######################Remove littermates#########################
remove <- c()
for(row in 1:length(rank.changes[,1])){
  up <- rank.changes[row,'Upmover']
  down <- rank.changes[row,'Downmover']
  ub <- filter(tblHyenas, ID == up)$Birthdate
  db <- filter(tblHyenas, ID == down)$Birthdate
  if(length(ub) & length(db) && ub == db && ub != '1/1/1900'){
    remove <- c(remove, row)
  }
}
rank.changes <- rank.changes[-remove,]
################################################################



#######################Presence of young########################
rank.changes$upNumKids <- NA
rank.changes$upNumMales <- NA
rank.changes$upNumFemales <- NA
rank.changes$downNumKids <- NA
rank.changes$downNumMales <- NA
rank.changes$downNumFemales <- NA
rank.changes$upDaughters <- NA
rank.changes$upSexRatio <- NA
rank.changes$downSexRatio <- NA
for(row in 1:length(rank.changes[,1])){
  kids <- filter(tblHyenas, Mom == rank.changes[row,'Upmover'], format(Birthdate, '%Y') == rank.changes[row,'Year'] | format(Birthdate, '%Y') == rank.changes[row,'Year']-1)
  kids <- filter(kids, is.na(Disappeared) | format(Disappeared, '%Y') < rank.changes[row,'Year'])
  rank.changes[row,'upDaughters'] <- length(filter(tblHyenas, Mom == rank.changes[row,'Upmover'], 
                                            format(Birthdate, '%Y') <= rank.changes[row,'Year'] & Sex == 'f',
                                            (is.na(Disappeared) | format(Disappeared, '%Y') < rank.changes[row,'Year']))[,1])
  rank.changes[row,'upSexRatio'] <- length(filter(tblHyenas, Mom == rank.changes[row,'Upmover'], 
                                          format(Birthdate, '%Y') <= rank.changes[row,'Year'] & Sex == 'f')[,1])/length(filter(tblHyenas, 
                                                                                                                        Mom == rank.changes[row,'Upmover'], 
                                                                                                                        format(Birthdate, '%Y') <= rank.changes[row,'Year'],
                                                                                                                        Sex == 'f' | Sex == 'm')[,1])
  rank.changes[row,'upNumKids'] <- length(kids[,1])
  rank.changes[row,'upNumMales'] <- length(filter(kids, Sex == 'm')[,1])
  rank.changes[row,'upNumFemales'] <- length(filter(kids, Sex == 'f')[,1])
  
  kids <- filter(tblHyenas, Mom == rank.changes[row,'Downmover'], format(Birthdate, '%Y') == rank.changes[row,'Year'] | format(Birthdate, '%Y') == rank.changes[row,'Year']-1)
  kids <- filter(kids, is.na(Disappeared) | format(Disappeared, '%Y') < rank.changes[row,'Year'])
  rank.changes[row,'downDaughters'] <- length(filter(tblHyenas, Mom == rank.changes[row,'Downmover'], 
                                              format(Birthdate, '%Y') <= rank.changes[row,'Year'] & Sex == 'f',
                                              (is.na(Disappeared) | format(Disappeared, '%Y') < rank.changes[row,'Year']))[,1])
  rank.changes[row,'downSexRatio'] <- length(filter(tblHyenas, Mom == rank.changes[row,'Downmover'], 
                                                    format(Birthdate, '%Y') <= rank.changes[row,'Year'] & Sex == 'f')[,1])/length(filter(tblHyenas, 
                                                                                                                                        Mom == rank.changes[row,'Downmover'], 
                                                                                                                                        format(Birthdate, '%Y') <= rank.changes[row,'Year'],
                                                                                                                                        Sex == 'f' | Sex == 'm')[,1])
  rank.changes[row,'downNumKids'] <- length(kids[,1])
  rank.changes[row,'downNumMales'] <- length(filter(kids, Sex == 'm')[,1])
  rank.changes[row,'downNumFemales'] <- length(filter(kids, Sex == 'f')[,1])
}
################################################################


########################Aggression intensity####################
rank.changes$upIntensity <- NA
rank.changes$downIntensity <- NA
for(row in 1:length(rank.changes[,1])){
  rank.changes[row,'upIntensity'] <- mean(as.numeric(filter(aggsFull, Agg == rank.changes[row,'Upmover'])$ThreatLevel), na.rm = T)
  rank.changes[row,'downIntensity'] <- mean(as.numeric(filter(aggsFull, Agg == rank.changes[row,'Downmover'])$ThreatLevel), na.rm = T)
  
}
################################################################



###################Emit/Receive Up hierarchy agg################
rank.changes$upEmit <- NA
rank.changes$downEmit <- NA
rank.changes$upRec <- NA
rank.changes$downRec <- NA
rank.changes$upUnprov <- NA
rank.changes$downUnprov <- NA
rank.changes$revCoal <- NA
for(row in 1:length(rank.changes[,1])){
  upmover <- rank.changes[row,'Upmover']
  downmover <- rank.changes[row,'Downmover']
  
  upmoverAggsRec <- filter(aggsFull, Recip == upmover, as.numeric(Year) >= min(ranks[ranks$ID == upmover,'Year']))
  upmoverUpAggsRec <- filter(upmoverAggsRec, AggRank > RecipRank)
  if(length(upmoverUpAggsRec[,1] == 0)){
    rank.changes[row,'upRec'] <- 0
  }else{
    rank.changes[row,'upRec'] <- 1
  }
  
  upmoverAggsEmit <- filter(aggsFull, Agg == upmover, as.numeric(Year) >= min(ranks[ranks$ID == upmover,'Year']))
  rank.changes[row,'upUnprov'] <- length(filter(upmoverAggsEmit, Context %in% c('unprov', 'unprovoked'))[,1])/length(upmoverAggsEmit[,1])
  upmoverUpAggsEmit <- filter(upmoverAggsEmit, AggRank > RecipRank)
  if(length(upmoverUpAggsEmit[,1] == 0)){
    rank.changes[row,'upEmit'] <- 0
  }else{
    rank.changes[row,'upEmit'] <- 1
  }
  
  downmoverAggsRec <- filter(aggsFull, Recip == downmover, as.numeric(Year) >= min(ranks[ranks$ID == downmover,'Year']))
  downmoverUpAggsRec <- filter(downmoverAggsRec, AggRank > RecipRank)
  if(length(downmoverUpAggsRec[,1] == 0)){
    rank.changes[row,'downRec'] <- 0
  }else{
    rank.changes[row,'downRec'] <- 1
  }
  rank.changes[row,'revCoal'] <- ifelse(length(filter(aggsFull, Recip == downmover, Agg == upmover, Seq > 0, Year == rank.changes[row,'Year'])[,1]),
                                        1,0)
    
  downmoverAggsEmit <- filter(aggsFull, Agg == downmover, as.numeric(Year) >= min(ranks[ranks$ID == downmover,'Year']))
  rank.changes[row,'downUnprov'] <- length(filter(downmoverAggsEmit, Context %in% c('unprov', 'unprovoked'))[,1])/length(downmoverAggsEmit[,1])
  downmoverUpAggsEmit <- filter(downmoverAggsEmit, AggRank > RecipRank)
  if(length(downmoverUpAggsEmit[,1] == 0)){
    rank.changes[row,'downEmit'] <- 0
  }else{
    rank.changes[row,'downEmit'] <- 1
  }
}
################################################################
boxplot(rank.changes$finalRank ~ rank.changes$revCoal)



###############Collapese data for comparison####################
rc.collapse <- data.frame()
for(year in unique(rank.changes$Year)){
  rcYear <- filter(rank.changes, Year == year)
  for(clan in unique(rcYear$Clan)){
    curRC <- filter(rcYear, Clan == clan)
    for(upmover in unique(curRC$Upmover)){
      temp <- curRC[curRC$Upmover == upmover,]
      attach(temp, warn.conflicts = FALSE)
      downKinAverage <- round(mean(downKin, na.rm = TRUE), digits = 2)
      downCoalAverage <- round(mean(downCoal, na.rm = TRUE), digits = 2)
      downNumCoalAverage <- round(mean(downNumCoal, na.rm = TRUE), digits = 2)
      downNumKidsAverage <- round(mean(downNumKids, na.rm = TRUE), digits = 2)
      downNumMalesAverage <- round(mean(downNumMales, na.rm = TRUE), digits = 2)
      downNumFemalesAverage <- round(mean(downNumFemales, na.rm = TRUE), digits = 2)
      downIntensityAverage <- round(mean(downIntensity, na.rm = TRUE), digits = 2)
     
      rand <- sample(length(Downmover),1)
      downKinRandom <- downKin[rand]
      downCoalRandom <- downCoal[rand]
      downNumCoalRandom <- downNumCoal[rand]
      
      rc.collapse <- rbind(rc.collapse, cbind(Year, Upmover, finalRank, Clan, upKin, upCoal, upNumCoal,
                                              upNumKids, upNumMales, upNumFemales, upIntensity,downIntensityAverage,
                                              downKinAverage, downCoalAverage, downNumCoalAverage,upUnprov, downUnprov,
                                              upDaughters, upSexRatio, downDaughters, downSexRatio,
                                              downKinRandom, downCoalRandom, downNumCoalRandom, upEmit, upRec, downEmit, downRec,
                                              downNumKidsAverage, downNumMalesAverage, downNumFemalesAverage))
      detach(temp)
    }
  }
}
rc.collapse[,c(-2,-4)] <- sapply(rc.collapse[,c(-2,-4)], as.numeric)
rc.collapse <- unique(rc.collapse)

################################################################




#######################Explore data#############################

hist(rc.collapse$upCoal - rc.collapse$downCoalAverage)
with(filter(rc.collapse, finalRank >= 5), hist(upCoal-downCoalAverage))

hist(rc.collapse$upNumCoal - rc.collapse$downNumCoalAverage)
with(filter(rc.collapse, finalRank <= 5), hist(upNumFemales - downNumFemalesAverage))

hist(rc.collapse$upNumKids - rc.collapse$downNumKidsAverage)
hist(rc.collapse$upNumFemales - rc.collapse$downNumFemalesAverage)

hist(rc.collapse$upIntensity - rc.collapse$downIntensityAverage)
################################################################



#########################Model things###########################
rc.dat.up <- data.frame()
rc.dat.down <- data.frame()
for(row in 1:length(rank.changes[,1])){
  rc.dat.up <- rbind(rc.dat.up, cbind(rank.changes[row,c("Year", 'Upmover', 'Clan', 'upKin', 'upIntensity', 'upUnprov', 'upCoal', 'upNumCoal', 'finalRank', 'upNumKids', 'upNumMales', 'upNumFemales', 'upEmit', 'upRec', 'upDaughters', 'upSexRatio')],1))
  rc.dat.down <- rbind(rc.dat.down, cbind(rank.changes[row,c('Year', 'Downmover', 'Clan', 'downKin', 'downIntensity', 'downUnprov', 'downCoal', 'downNumCoal', 'initRank', 'downNumKids', 'downNumMales', 'downNumFemales', 'downEmit', 'downRec', 'downDaughters', 'downSexRatio')],0))
}

names(rc.dat.up) <- c('Year', 'ID', 'Clan', 'Kin', 'Intensity', 'Unprov', 'CoalPartners', 'NumCoal', 'Rank', 'NumKids', 'NumMales', 'NumFemales', 'Emit', 'Rec', 'Daughters', 'SexRatio', 'upDown')
names(rc.dat.down) <- c('Year', 'ID', 'Clan', 'Kin', 'Intensity', 'Unprov', 'CoalPartners', 'NumCoal', 'Rank', 'NumKids', 'NumMales', 'NumFemales', 'Emit','Rec', 'Daughters', 'SexRatio', 'upDown')
rc.dat <- rbind(rc.dat.up, rc.dat.down)
rc.dat <- unique(rc.dat)

library(lme4)
summary(glmer(data = filter(rc.dat, Rank <= 5), formula = upDown ~ CoalPartners + (1|Year) ,family = binomial))
summary(glm(dat = rc.dat, formula = upDown ~ Unprov +Unprov:Rank, family = binomial))
################################################################


#########################behavior rates#########################

################################################################


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

rc.dat$fecalT <- left_join(x = rc.dat, y = ttt.am.mean, by = c("Year" = "poop_year", "ID" = "hyenaID"))$ng.g
boxplot(rc.dat$fecalT ~ rc.dat$upDown)
#######################################################################



########################By year#################################
yearly <- data.frame()
for(clan in unique(ranks[,'Clan'])){
  for(year in unique(ranks[ranks$Clan == clan,'Year'])){
    curRanks <- ranks[ranks$Clan == clan & ranks$Year == year,]
    alpha <- curRanks[curRanks$Rank == 1, 'IDold']
    attackOnAlpha <- length(filter(aggsWinner, Recip == alpha, Year == year)[,1])
    alphaDaughters <- length(filter(tblHyenas, Mom == alpha, ID %in% curRanks$ID)[,1])
    alphaIntensity <- mean(as.numeric(filter(aggsFull, Agg == alpha, Year == year)$ThreatLevel), na.rm = T)/mean(as.numeric(filter(aggsFull, Agg == alpha, Year != year)$ThreatLevel), na.rm = T)
    clanSize <- length(curRanks[,1])
    rankChanges <- length(rc.dat[rc.dat$Year == year & rc.dat$Clan == clan,1])
    movers <- min(length(rc.dat[rc.dat$Year == year & rc.dat$Clan == clan & rc.dat$upDown == 1, 1]), length(rc.dat[rc.dat$Year == year & rc.dat$Clan == clan & rc.dat$upDown == 0, 1]))
    yearly <- rbind(yearly, cbind(clan, year, alpha, attackOnAlpha, alphaDaughters, alphaIntensity, clanSize, rankChanges, movers))
  }
}
yearly[,c(-1,-3)] <- sapply(yearly[,c(-1,-3)], as.numeric)

talekDeath <- as.numeric(c(format(filter(tblHyenas, ID == 'bsh')$Disappeared, '%Y'),
                          format(filter(tblHyenas, ID == 'mrph')$Disappeared, '%Y'),
                          format(filter(tblHyenas, ID == 'dion')$Disappeared, '%Y')))

yearly$alphaDeath <- NA
for(y in yearly[yearly$clan == 'talek',]$year){
  ans <- y - talekDeath
  if(length(ans[ans >= 0])){yearly[yearly$year == y & yearly$clan == 'talek', 'alphaDeath'] <- min(ans[ans >= 0])}
}

yearly$FitChange <- 0
for(clan in unique(ranks$Clan)){
  for(year in unique(ranks[ranks$Clan == clan,'Year'])){
    curRanks <- filter(ranks, Clan == clan, Year == year)
    curRanks[is.na(curRanks$Diff), 'Diff'] <- 0
    yearly[yearly$clan == clan & yearly$year == year,'FitChange'] <- mean(abs(curRanks$Diff))
  }
}


 
################################################################
summary(glm(data = yearly[yearly$clan == 'talek',], as.numeric(rankChanges) ~ as.numeric(alphaDaughters), family = poisson))
summary(glm(data = yearly, FitChange ~ alphaIntensity, family= beta))

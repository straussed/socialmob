####
setwd("~/Documents/Research/RCIISI")
options(stringsAsFactors = FALSE)
library(network)
library(dplyr)
library(pscl)
library(reshape2)

ranks <- read.csv("~/Documents/Research/RCIISI/Working_19_Oct_2016_12/Results.csv")
####Aggressive interactions####
aggs <- read.csv("~/Documents/Research/RCIISI/AggsISIAllClans.csv")
names(aggs) <- c("Agg", "Recip", "Date", "R1", "R2", "R3", "Seq", "ThreatLevel", "Context", "Clan")
aggs$Date <- fix.dates(aggs$Date)
aggs$Agg <- tolower(aggs$Agg)
aggs$Recip <- tolower(aggs$Recip)

###tblSessions###
tsess <- read.csv('~/Documents/Fisibase/tblSessions.csv')
tsess$clan = rep('talek')
ssess <- read.csv('~/Documents/Fisibase/tblSessions_serena.csv')
columns <- c('SESSION', 'LOCATION', 'DATE', 'START', 'STOP', 'HYENAS', 'TIME', 'clan')
sessions <- rbind(tsess[,columns], ssess[,columns])
names(sessions) <- tolower(names(sessions))
sessions$date <- fix.dates.longy(sessions$date)
sessions$start <- as.POSIXct(fix.times(sessions$start))
sessions$stop <- as.POSIXct(fix.times(sessions$stop))
sessions$time <- as.numeric(sessions$stop - sessions$start)/60
sessions$numHyenas <- apply(sessions[,'hyenas', drop = F], 1, function(x){length(unlist(strsplit(x, ',')))})-1

hps <- sapply(sessions$hyenas, strsplit, ',')
names(hps) <- sessions$session

checkSes <- function(ses, hyena){
  return(hyena %in% strsplit(ses$hyenas, ','))
}



####Rank changes from IdentifyRankChanges
rank.changes
l <- length(rank.changes[,1])
rc.indiv <- data.frame()
for(row in 1:length(rank.changes[,1])){
  upmover <- rank.changes[row, 'Upmover']
  numCoals <- length(filter(aggs, Agg == upmover, format(Date, '%Y') == rank.changes[row,'Year'], Seq >= 1)[,1])
  numCoalsPrev <- length(filter(aggs, Agg == upmover, format(Date, '%Y') == rank.changes[row,'Year']-1, Seq >= 1)[,1])
  numCoalsNext <- length(filter(aggs, Agg == upmover, format(Date, '%Y') == rank.changes[row,'Year']+1, Seq >= 1)[,1])
  rank <- as.numeric(ranks[ranks$NewOrder == upmover & ranks$Year == rank.changes[row,'Year'], 'Rank'])
  ##rates for coalitions
  sessSeen <- sessions[grep(paste(',', upmover, ',', sep = ''), sessions$hyenas),]
  sessSeenYear <- sessSeen[format(sessSeen$date, '%Y') == rank.changes[row,'Year'],]
  timeSeen <- sum(filter(sessSeenYear, numHyenas >= 3)[,'time'])
  rateCoals <- numCoals/timeSeen
  rc.indiv <- rbind(rc.indiv, c(ID = upmover, timeSeen, rateCoals, numCoals, numCoalsPrev, numCoalsNext, Year = rank.changes[row,'Year'], Rank = rank, Direction = 'Upmover'))
  ####Downmover
  downmover <- rank.changes[row, 'Downmover']
  numCoals <- length(filter(aggs, Agg == downmover, format(Date, '%Y') == rank.changes[row,'Year'], Seq >= 1)[,1])
  numCoalsPrev <- length(filter(aggs, Agg == downmover, format(Date, '%Y') == rank.changes[row,'Year']-1, Seq >= 1)[,1])
  numCoalsNext <- length(filter(aggs, Agg == downmover, format(Date, '%Y') == rank.changes[row,'Year']+1, Seq >= 1)[,1])
  rank <- as.numeric(ranks[ranks$NewOrder == downmover & ranks$Year == rank.changes[row,'Year'], 'Rank'])
  sessSeen <- sessions[grep(paste(',', downmover, ',', sep = ''), sessions$hyenas),]
  sessSeenYear <- sessSeen[format(sessSeen$date, '%Y') == rank.changes[row,'Year'],]
  timeSeen <- sum(filter(sessSeenYear, numHyenas >= 3)[,'time'])
  rateCoals <- numCoals/timeSeen
  rc.indiv <- rbind(rc.indiv, c(ID = downmover, timeSeen, rateCoals, numCoals, numCoalsPrev, numCoalsNext, Year = rank.changes[row,'Year'], Rank = rank, Direction = 'Downmover'))
}

names(rc.indiv) <- c('ID', 'timeSeen', 'rateCoals', 'numCoals', 'numCoalsPrev', 'numCoalsNext', 'Year', 'Rank', 'Direction')
rc.indiv <- unique(rc.indiv)
ggplot(data = unique(rc.indiv), aes(y = scale(as.numeric(rateCoals)), x = factor(Direction)))+
  geom_boxplot()
  #geom_abline(slope = 1, intercept = 0)
  #geom_jitter(aes(x = as.numeric(numCoalsNext), y = as.numeric(numCoals), col= 'blue'))
  
summary(glm(data = filter(rc.indiv), as.numeric(factor(Direction))-1 ~ as.numeric(rateCoals) * as.numeric(Rank), family = binomial))

ggplot(data = rc.indiv, aes(y = as.numeric(rateCoals), x = as.numeric(Rank), col = Direction))+
  geom_point()
hist(log(as.numeric(unique(rc.indiv$rateCoals))))

lm(data = rc.indiv, numCoals ~ factor(Direction))
rc.indiv$Direction <- as.factor(rc.indiv$Direction)


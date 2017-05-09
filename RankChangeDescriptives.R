##################################################################
##                        Eli Strauss                           ##
##                 Rank Change Descriptives                     ##
##                     March 27th, 2017                         ##
##################################################################

library(dplyr)
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")

############################Add rank change info to ranks##############
for(row in 1:length(ranks[,1])){
  ranks.temp <- filter(ranks, Clan == ranks[row,'Clan'], Year == ranks[row,'Year'])
  if(ranks[row,'Rank'] == ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'None'
  }else if(ranks[row,'Rank'] < ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'Up'
  }else if(ranks[row,'Rank'] > ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'Down'
  }
}
ranks$RankChange <- as.factor(ranks$RankChange)


#######plot rank changes with lines for each individual######
#ranks$ID <- ranks$NewOrder
par(fig = c(0,1,0,1))
plot(data = ranks, Rank ~ Year, type = 'n', ylim = c(55,0), main = 'Talek', xlim = c(1987,2015))
for(id in unique(ranks[ranks$Clan == 'talek',]$ID)){
  #clr <- ifelse(alphaDescent(id), 'blue', 'black')
  #if(id == 'nav'){clr <- 'red'}
  clr <- 'black'
  if('Up' %in% ranks[ranks$ID == id,'RankChange'] | 'Down' %in% ranks[ranks$ID == id, 'RankChange']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = clr))}
}

par(fig = c(0.07,(.07+(6/25)+.02),0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'south',], Rank ~ Year, type = 'n', ylim = c(25,0), main = 'South')
for(id in unique(ranks[ranks$Clan == 'south',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'RankChange'] | 'Down' %in% ranks[ranks$ID == id, 'RankChange']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}

par(fig = c(0,1,0,1), new = T)
par(fig = c((.07+(6/25)+.02),(.33+.26), 0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'north',], Rank ~ Year, type = 'n', ylim = c(20,0), main = 'North')
for(id in unique(ranks[ranks$Clan == 'north',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'RankChange'] | 'Down' %in% ranks[ranks$ID == id, 'RankChange']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}

par(fig = c(0,1,0,1), new = T)
par(fig = c(0.59,.59+.26, 0.10,.40), new = T)
plot(data = ranks[ranks$Clan == 'hz',], Rank ~ Year, type = 'n', ylim = c(20,0), main = 'Happy Zebra', axes = F, frame.plot = T)
axis(1, at = seq(2009, 2015))
axis(2, at = c(0, 5, 10, 15, 20))
for(id in unique(ranks[ranks$Clan == 'hz',]$ID)){
  if('Up' %in% ranks[ranks$ID == id,'RankChange'] | 'Down' %in% ranks[ranks$ID == id, 'RankChange']){
    with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5, col = 'red'))
  }else{with(filter(ranks, ID == id), lines(Rank ~ Year, lwd = 1.5))}
}


##################Identify rank change period#######################
rank.changers <- ranks[ranks$RankChange != 'None',]
for(row in 1:length(rank.changers[,1])){
  up <- ifelse(rank.changers[row,'RankChange'] == 'Up', 1, 0)
  year <- rank.changers[row,'Year']
  clan <- rank.changers[row,'Clan']
  yr.ranks <- ranks[ranks$Year == year & ranks$Clan == clan,]
  id <- rank.changers[row,'ID']
  ##identify the individuals ID swapped places with
  if(length(yr.ranks[yr.ranks$RankChange == 'Up',1]) == 1 | length(yr.ranks[yr.ranks$RankChange == 'Down',1]) == 1){
    if(up){swappers = yr.ranks[yr.ranks$RankChange == 'Down','ID']}else{swappers <- yr.ranks[yr.ranks$RankChange == 'Up','ID']}
  }else{
    if(up){
      swappers <- filter(yr.ranks, ID %in% yr.ranks$ID[which(yr.ranks$ID == id)+1:length(yr.ranks$ID)], ID %in% yr.ranks$IDold[1:which(yr.ranks$IDold == id)])$ID
    }else{
      swappers <- filter(yr.ranks, ID %in% yr.ranks$ID[1:which(yr.ranks$ID == id)], ID %in% yr.ranks$IDold[which(yr.ranks$IDold == id)+1:length(yr.ranks$ID)])$ID
    }
  }
  
  ##calculate end of swap period as first aggression indicating reversal with highest/lowest ranking of the swappers (depending on direction).
  ##If rank change is implied by interaction with other individuals, return NA
  if(up){
    end <- min(filter(aggsWinner, Year == year, Agg == id, Recip == swappers[which.min(yr.ranks[yr.ranks$ID %in% swappers,]$Rank)])$Date)
  }else{
    end <- min(filter(aggsWinner, Year == year, Agg %in% swappers[which.max(yr.ranks[yr.ranks$ID %in% swappers,]$Rank)], Recip == id)$Date)
  }
  
  ##calculate start of swap period as first aggression indicating reversal with any of the swappers.
  if(up){
    start <- max(filter(aggsWinner, Year <= year, Agg %in% swappers, Recip == id, Date < end)$Date)
  }else{
    start <- max(filter(aggsWinner, Year <= year, Agg == id, Recip %in% swappers, Date < end)$Date)
  }
  
  ranks[ranks$ID == id & ranks$Year == year, 'rc.start']<- start
  ranks[ranks$ID == id & ranks$Year == year, 'rc.end']<- end
  ranks[ranks$ID == id & ranks$Year == year, 'rc.duration'] <- end-start
}

rank.changers.all <- ranks[!is.na(ranks$rc.duration),]
rank.changers.full.info <- ranks[is.finite(ranks$rc.duration),]
rank.changers.implied <- ranks[is.infinite(ranks$rc.end),]
rank.changers.juvenile <- ranks[is.finite(ranks$rc.end) & is.infinite(ranks$rc.start),]
rank.changers.end <- ranks[is.finite(ranks$rc.end),]

rank.changers.juvenile

rank.changers.all$Birthdate <- left_join(rank.changers.all, tblHyenas, by = c('ID'))$Birthdate
rank.changers.known.birth <- rank.changers.all[complete.cases(rank.changers.all$Birthdate),]
rank.changers.known.birth$AgeAtStart <- as.Date(rank.changers.known.birth$rc.start, origin = '1970-01-01') - rank.changers.known.birth$Birthdate
rank.changers.known.birth$AgeAtEnd <- as.Date(rank.changers.known.birth$rc.end, origin = '1970-01-01') - rank.changers.known.birth$Birthdate
rank.changers.age <- rank.changers.known.birth
rank.changers.age[is.infinite(rank.changers.age$AgeAtStart),'AgeAtStart'] <- 0
rank.changers.age <- rank.changers.age[is.finite(rank.changers.age$AgeAtEnd),]
hist(as.numeric(rank.changers.age$AgeAtStart))
rank.changers.age$AgeAtStart <- as.numeric(rank.changers.age$AgeAtStart)/365
rank.changers.age$AgeAtEnd <- as.numeric(rank.changers.age$AgeAtEnd)/365
rank.changers.age$yaxis <- NA
rank.changers.age <- rank.changers.age[order(rank.changers.age$AgeAtStart, decreasing = T),]
rank.changers.age[rank.changers.age$RankChange == 'Up','yaxis'] <- seq(from = 2, to = 1.05, length.out = length(rank.changers.age[rank.changers.age$RankChange == 'Up','yaxis']))
rank.changers.age[rank.changers.age$RankChange == 'Down','yaxis'] <- seq(from = 0, to = .95, length.out = length(rank.changers.age[rank.changers.age$RankChange == 'Down','yaxis']))



######Plot lines for age####
change.colors <- c('red', 'blue')
plot(data = rank.changers.age, yaxis ~ AgeAtStart, type = 'n', xlim = c(0, 15), yaxt = 'n', ylab = 'Direction of Rank Change', xlab = 'Age range of period of change')
abline(h = 1, lty = 3, lwd = 2)
rect(xleft = -1, xright = 1.5, ybottom = -1, ytop = 3, angle = 45, col = rgb(.3,.3,.3, alpha = .5))
mtext(text = c('Up'), side = 2, at = 1.5, cex = 1.2)
mtext(text = c('Down'), side = 2, at = 0.5, cex = 1.2)
for(row in 1:length(rank.changers.age[,1])){
  color = change.colors[ifelse(rank.changers.age[row,'RankChange'] == 'Up',1,2)]
  lines(x = c(rank.changers.age[row,'AgeAtStart'], rank.changers.age[row,'AgeAtEnd']), 
        y = c(rank.changers.age[row,'yaxis'], rank.changers.age[row,'yaxis']), 
        col = color,
        lwd = 2)
  text(x = rank.changers.age[row,'AgeAtEnd'] + 1, y = rank.changers.age[row,'yaxis'], labels = rank.changers.age[row,'ID'])
  lines(x = c(rank.changers.age[row,'AgeAtStart'], rank.changers.age[row,'AgeAtStart']),
        y = c(rank.changers.age[row,'yaxis']-.01, rank.changers.age[row,'yaxis']+.01),
        col = color,
        lwd = 2)
  lines(x = c(rank.changers.age[row,'AgeAtEnd'], rank.changers.age[row,'AgeAtEnd']),
        y = c(rank.changers.age[row,'yaxis']-.01, rank.changers.age[row,'yaxis']+.01),
        col = color,
        lwd = 2)
}

hist(as.numeric(rank.changers.age$AgeAtEnd), breaks = 10, main = 'Age at rank change', xlab = 'Years', col = 'grey')

########Calculate  proportions of ranks assigned where there was a reversal####
rank.assignments <- 0
tot.changes <- 0
for(clan in unique(ranks$Clan)){
  r <- filter(ranks, Clan == clan)
  r <- r[r$Year != min(r$Year) & r$Year != max(r$Year),]
  rank.assignments <- rank.assignments + length(r[,1])
  tot.changes <- tot.changes+ length(r[r$RankChange != 'None',1])
}
tot.changes/rank.assignments

#####Calculate proportion of rank inheritance corresponding to Kawamura's rules#####
newf.count <- 0
krule.count <- 0
for(clan in unique(ranks$Clan)){
  r <- filter(ranks, Clan == clan)
  r <- r[r$Year != min(r$Year) & r$Year != max(r$Year),]
  for(row in 1:length(r[,1])){
    id <- r[row,'ID']
    if(r[row,'Year'] == min(r[r$ID == id,'Year'])){
      newf.count <- newf.count + 1
      if(r[row,'ID'] == r[row,'IDold']){
        krule.count <- krule.count + 1
      }
    }
  }
}
krule.count/newf.count


result <- c()
alphas <- c(unique(filter(ranks, ranks$Rank == 1)$ID), 'kb')
alpha_data <- tibble(ID = NA, Year = NA, Rate = NA)
alphas_that_recuse <- tibble(ID = c('kb', 'clov', 'hel'),
                             Year = c(1988, 2012, 2014))
for(row in 1:tibdim(alphas_that_recuse)){
  alpha = alphas_that_recuse[row,]$ID
  birth_year <- as.numeric(format(tblHyenas[tblHyenas$ID == alpha,]$FirstSeen, '%Y'))
  last_year <- as.numeric(format(tblHyenas[tblHyenas$ID == alpha,]$LastSeen, '%Y'))
  change_year <- alphas_that_recuse[row,]$Year
  
  for(year in birth_year:change_year){
    repro_rate <- tibdim(filter(tblHyenas, Mom == alpha, format(Birthdate, '%Y')  < year, format(Birthdate, '%Y')  >= year-3))/3
    alpha_data <- add_row(alpha_data, ID = alpha, Year = year, Rate = repro_rate)
  }
}


##########
ggplot(data = filter(ranks, RankChange != 'None'), aes(stan.rank))+
  geom_histogram(bins = 22)+
  xlab('Rank') + 
  main('Rank changes are concen')

ggplot(data = filter(ranks, RankCategory == 'High'), aes(y=abs(RankDiffAbs),x=ai_top3_deg))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_classic()


####Rank changes after death of alpha###
tibble(ID = c('koi', 'bsh', 'mrph', 'rbc'),
        ChangeAfterDeath = c('yes', 'yes', 'yes', 'yes'))



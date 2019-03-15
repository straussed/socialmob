################################################################################
#           Miscellaneous descriptives and plots for publication               #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)

load('8.ranks_with_rank_change.RData')
source('0.multiplot.R')

###Read in data####
aggressions <- read.csv('0.rawdata/tblAggression.csv')
aggressions$date <- as.Date(aggressions$date)
aggressions$Year <- as.numeric(format(aggressions$date, '%Y'))

aggsFull <- aggressions
aggsFull$enteredby <- NULL
aggsFull$doublecheck <- NULL

excludeResponse <- c("ignores", "ignore", "ct", "counterattack", "counter", "counters", "counterattacks")
aggsWinner <- filter(aggsFull, !response1 %in% excludeResponse, !response2 %in% excludeResponse, !response3 %in% excludeResponse) 
aggsWinner <- filter(aggsWinner, !context %in% c('ct', 'counter', 'counterattack'))

names(aggsFull) <- c('Clan', 'AggID', 'Session', 'Date', 'Time', 'Agg', 'Recip',
                     'Group','GroupComp', 'Context', 'AltContext', 'B1', 'B2',
                     'R1', 'R2', 'R3', 'Seq','Year')

names(aggsWinner) <- c('Clan', 'AggID', 'Session', 'Date', 'Time', 'Agg', 'Recip',
                     'Group','GroupComp', 'Context', 'AltContext', 'B1', 'B2',
                     'R1', 'R2', 'R3', 'Seq','Year')

aggsFull$Agg <- tolower(aggsFull$Agg)
aggsFull$Recip <- tolower(aggsFull$Recip)

aggsWinner$Agg <- tolower(aggsWinner$Agg)
aggsWinner$Recip <- tolower(aggsWinner$Recip)



#demographic information
hyenas <- read.csv('0.rawdata/tblHyenas.csv')
hyenas[hyenas$clan == 'serena s',]$clan <- 'south'
hyenas[hyenas$clan == 'serena n',]$clan <- 'north'
hyenas[hyenas$clan == 'happy zebra',]$clan <- 'hz'
hyenas$birthdate <- as.Date(hyenas$birthdate, format = '%d-%b-%y')
hyenas$disappeared <- as.Date(hyenas$disappeared, format = '%d-%b-%y')

tblHyenas <- hyenas
names(tblHyenas) <- c('ID', 'Last.Updated', 'SampleID', 'Eartage', 'Name', 'PrevID',
                      'Sex', 'AgeClass', 'Status', 'FirstSeen', 'DenGrad', 'Disappeared',
                      'Mom', 'Birthdate','NumLittermates',
                      'LitRank', 'ArrivedDen', 'LeaveDen', 'Fate', 'MortalitySource',
                      'DeathDate', 'Weaned', 'Clan', 'Park', 'Notes')

###Descriptives###
####Number of aggressions
all.aggs <- inner_join(aggsWinner, ranks, by = c('Agg' = 'ID', 'Year')) %>% 
  inner_join(ranks, by = c('Recip' = 'ID', 'Year'))

## Coalitions used to measure bond strength - From 5.create_coalition_networks.R
## sum(coal.summary$num.coals) )  ###1913 coalitions used to measure bond strength

yearly.aggs <- all.aggs %>% 
  group_by(Clan, Year) %>% 
  summarize(num.aggs = length(Agg))


ids.per.year <- ranks %>% 
  group_by(Clan, Year) %>% 
  summarize(num.ids = length(ID))

ids.per.year %>% group_by(Clan) %>% summarise(clan.group.size = mean(num.ids))

yearly.aggs <- left_join(yearly.aggs, ids.per.year)
yearly.aggs$aggs.per.id <- yearly.aggs$num.aggs/yearly.aggs$num.ids

###Rank changes per year###
rc.per.year <- ranks %>% 
  group_by(Clan, Year) %>%
  summarize(prop.changes = sum(RankChange != 'None')/length(RankChange),
            num.changes = sum(RankChange != 'None'),
            no.changes = sum(RankChange == 'None'))



yearly.aggs <- left_join(yearly.aggs, rc.per.year)

pdf(file = 'plots/13.FigS1.pdf', 4, 4)
ggplot(data = yearly.aggs, aes(x = aggs.per.id, y = prop.changes)) + 
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_classic()+
  xlab('Aggressive acts observed per individual')+
  ylab('Proportion of individuals experiencing rank reversal')+
  annotate('text', x = 15, y = 0.75, label = 'Pearson correlation = 0.184')

dev.off()

glm(data = yearly.aggs, cbind(num.changes, no.changes) ~ aggs.per.id, family = 'binomial') %>% summary()
cor(yearly.aggs$prop.changes, yearly.aggs$aggs.per.id)


##Average group size
ranks %>% group_by(Clan, Year) %>% summarize(num_inds = length(ID)) %>% 
  ungroup() %>% group_by(Clan) %>% summarize(mean_size = mean(num_inds))

###Proportion of individuals who inherit rank according to MRI and YA
mri <- 0
notmri <- 0
for(id in unique(ranks$ID)){
  clan <- ranks[ranks$ID == id,]$Clan[1]
  first.year <- min(ranks[ranks$ID == id,]$Year)
  if(first.year == min(ranks[ranks$Clan == clan,]$Year)){next}
  mom <- filter(tblHyenas, ID == id)$Mom
  if(is.na(mom) | !mom %in% ranks[ranks$Year == first.year,]$ID){next}
  ranks.year <- filter(ranks, Year == first.year, Clan == clan)
  if(which(ranks.year$ID == id) == 1){
    notmri <- notmri + 1
    next
  }
  up.hyena <- ranks.year[which(ranks.year$ID == id)-1,]$ID
  if(up.hyena == mom){
    mri <- mri + 1
  }else{
    if(which(ranks.year$ID == id) == 2){
      notmri <- notmri + 1
      next
    }
    up.hyena2 <- ranks[which(ranks.year$ID == id)-2,]$ID
    if(up.hyena2 != mom){
      notmri <- notmri + 1
    }else{
      if(!is.na(filter(tblHyenas, ID == up.hyena)$Birthdate) && 
         !is.na(filter(tblHyenas, ID == id)$Birthdate)){
        if(filter(tblHyenas, ID == up.hyena)$Birthdate == filter(tblHyenas, ID == id)$Birthdate){
          mri <- mri+1
        }else{
          notmri <- notmri+1
        }
      }
    }
  }
}
mri/(mri+notmri)

####Proportion of individuals involved in a rank change
length(unique(filter(ranks, RankChange != 'None')$ID))/length(unique(ranks$ID))

###Number of rank changes
table(ranks$RankChange)
sum(table(ranks$RankChange)[c(1,3,4)])/sum(table(ranks$RankChange))


##Proportion of rank reversals where daughter passes mother##

up.movers <- filter(ranks, RankDiffAbs > 0)
not.pass.mom <- 0
pass.mom <- 0


mom.daughter <- data.frame()
for(row in 1:nrow(up.movers)){
  mom <- filter(tblHyenas, ID == up.movers[row,]$ID)$Mom
  cranks <- filter(ranks, Clan == up.movers[row,]$Clan, Year == up.movers[row,]$Year)
  if(length(mom) > 1) print(mom)
  if(is.na(mom) | mom == ''){
    next
  }
  if(mom %in% cranks$ID){
    if(filter(cranks, ID == up.movers[row,]$ID)$Rank < filter(cranks, ID == mom)$Rank & 
       filter(cranks, OldOrder == up.movers[row,]$ID)$Rank > filter(cranks, OldOrder == mom)$Rank){
      pass.mom <- pass.mom + 1
      print(paste0(mom, 'pass'))
      mom.daughter <- rbind(mom.daughter, data.frame(up.movers[row,], Pass = 'pass'))
    }else{
      not.pass.mom <- not.pass.mom + 1
      mom.daughter <- rbind(mom.daughter, data.frame(up.movers[row,], Pass = 'no pass'))
    }
  }else{
    not.pass.mom <- not.pass.mom + 1
    mom.daughter <- rbind(mom.daughter, data.frame(up.movers[row,], Pass = 'no pass'))
  }
}
not.pass.mom
pass.mom
pass.mom/(pass.mom + not.pass.mom)




##Figure 1
talek <- ggplot(data = filter(ranks, Clan == 'talek'), aes(y = Rank, x = Year)) + 
  ylim(52,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black', size = 0.3)+
  theme(plot.margin = unit(x = c(0,8.5,5.5,5.5), 'pt'))+
  ggtitle(label = NULL, subtitle = 'Talek')+
  theme(axis.text = element_text(size = 9), plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 10))

north <- ggplot(data = filter(ranks, Clan == 'north'), aes(y = Rank, x = Year)) + 
  ylim(22,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black', size = 0.3) + 
  theme(plot.margin = unit(x = c(0,0,0,5.5), 'pt'))+
  ggtitle(label = NULL, subtitle = 'Serena North')+
  theme(axis.text = element_text(size = 9), plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 10))

south <- ggplot(data = filter(ranks, Clan == 'south'), aes(y = Rank, x = Year)) + 
  ylim(25,0)+
  theme_classic() + 
  ylab(label = NULL)+
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black', size = 0.3)+
  theme(plot.margin = unit(x = c(0,0,0,5.5), 'pt'))+
  ggtitle(label = NULL, subtitle = 'Serena South')+
  theme(axis.text = element_text(size = 9), plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 10))

hz <- ggplot(data = filter(ranks, Clan == 'hz'), aes(y = Rank, x = Year)) + 
  ylim(20,0)+
  theme_classic() + 
  ylab(label= NULL)+
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black', size= 0.3)+
  theme(plot.margin = unit(x = c(0,5.5,0,5.5), 'pt'))+
  ggtitle(label = NULL, subtitle = 'Happy Zebra')+
  scale_x_continuous(breaks = c(2010, 2012, 2014))+
  theme(axis.text = element_text(size = 9), plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 10))

pdf(file = "plots/13.Fig1.pdf",
    width = 4.5,
    height = 3.5)
multiplot(talek, north, south, hz, layout = matrix(c(1,1,1,1,1,1,2,3,4), nrow = 3, byrow = TRUE))
dev.off()


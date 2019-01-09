################################################################################
#                    Prep raw data for use with DynaRankR                      #
#                                                                              #
#                                                                              #
#                            By Eli Strauss                                    #
#                                                                              #
#                            December 2018                                     #
################################################################################

library(dplyr)
options(stringsAsFactors = FALSE)

setwd('~/Documents/Research/socialmob/ranksDynaRank/0.rawdata/')

#####Set seed####
set.seed(1989)

####Read in data####
aggressions <- read.csv('tblAggression.csv')
aggressions$date <- as.Date(aggressions$date)
aggressions$Year <- as.numeric(format(aggressions$date, '%Y'))

#demographic information
#hyenas <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblHyenas.csv')
hyenas <- read.csv('tblHyenas.csv')
hyenas[hyenas$clan == 'serena s',]$clan <- 'south'
hyenas[hyenas$clan == 'serena n',]$clan <- 'north'
hyenas[hyenas$clan == 'happy zebra',]$clan <- 'hz'


#Fix id loda -> luda
hyenas[hyenas$id == 'loda',]$id <- 'luda'

#Fix bd
hyenas[hyenas$id == 'bd',]$birthdate <- hyenas[hyenas$id == 'bd',]$first.seen

#Fix tru
hyenas[hyenas$id == 'tru',]$mom <- 'bor'
  
  
#remove hyena 44 who disappeared before our data begin
hyenas <- hyenas[-which(hyenas$id == '44'),]

#cash appears twice
hyenas <- filter(hyenas, !(id == 'cash' & clan == 'north'))

hyenas$birthdate <- as.Date(hyenas$birthdate)
hyenas$death.date <- as.Date(hyenas$death.date)
hyenas$disappeared <- as.Date(hyenas$disappeared)

#sessions and hyenas per session to fix some missing disappeared dates
hps <- read.csv('tblHyenasPerSession.csv', colClasses = 'character')
sessions <- read.csv('tblSessions.csv', colClasses = 'character')
hps$date <- left_join(hps, sessions, by = 'session')$date
hps$date <- as.Date(hps$date)

last.seen <- hps %>% group_by(hyena) %>% summarize(last.seen = max(date, na.rm = T))

hyenas[is.na(hyenas$disappeared),]$disappeared <- 
  left_join(hyenas[is.na(hyenas$disappeared),],
            last.seen, by = c('id' = 'hyena'))$last.seen

#Clan membership after fission
talekMembership <- read.csv('ClanMembership.csv')

excludeResponse <- c("ignores", "ignore", "ct", "counterattack", "counter", "counters", "counterattacks")
aggsWinner <- filter(aggressions, !response1 %in% excludeResponse, !response2 %in% excludeResponse, !response3 %in% excludeResponse) 
aggsWinner <- filter(aggsWinner, !context %in% c('ct', 'counter', 'counterattack'))

####List of clans####
clans <- c('talek', 'north', 'south', 'hz')
ranks <- data.frame()
all.intx <- data.frame()
all.conts <- data.frame()

for(currentClan in clans){
  ####Prep data for current clan
  #starting ranks
  
  initial.ranks <- read.csv(paste0('iranks_', currentClan, '.csv'))
  #aggressions
  aggs <- filter(aggsWinner, clan == currentClan)
  #contestants
  females <- filter(hyenas, clan == currentClan, sex == 'f', !is.na(birthdate) | id %in% initial.ranks$ID)
  
  ##hyenas of current clan
  currentHyenas <- filter(hyenas, clan == currentClan)
  
  females$EndYear <- format(do.call(pmin, c(females[,c('death.date', 'disappeared')], na.rm = T)), '%Y')
  
  #Females are added the first *complete* year that they are at least 1.5 years old
  females$StartYear <- format(females$birthdate + 365*2.5, '%Y')
  females <-females[,c('id', 'StartYear', 'EndYear')]
  names(females) <- c('ID', 'StartYear', 'EndYear')
  
  ##If the clan is talek, remove females that leave in 2000
  if(currentClan == 'talek'){
    ##Make list of talek east
    easties <- filter(talekMembership, Membership == 'e')
    for(eh in unique(females$ID)){
      ehmom <- filter(currentHyenas, id == eh)$mom
      if(ehmom %in% easties$ID){easties <- rbind(easties, c(eh, 'e', 'kid', 'EDS'))}
    }
    ####remove talek east
    females[females$ID %in% easties$ID,]$EndYear <- 1999
  }
  
  first.year <- min(initial.ranks$Year)
  if(currentClan %in% c('talek', 'hz')){
    last.year <- 2015
  }else{
    last.year <- 2016
  }
  
  ###Assign final tables to be used by DynaRank package
  #Initial ranks
  assign(paste0('initial.ranks.', currentClan), initial.ranks$ID)
  females[females$ID %in% initial.ranks$ID,'StartYear'] <- first.year
  
  #Contestants
  females <- filter(females, StartYear <= last.year,
                    StartYear <= EndYear)
  contestants <- data.frame()
  for(id in females$ID){
    contestants <- rbind(contestants, data.frame(id, period = seq(from =filter(females, ID == id)$StartYear, to = filter(females, ID == id)$EndYear, by = 1)))
  }
  contestants <- filter(contestants, period >= first.year, period <= last.year)
  contestants$convention1 <- left_join(contestants, hyenas, by = 'id')$mom
  contestants$convention2 <- left_join(contestants, hyenas, by = 'id')$litrank
  contestants <- filter(contestants, convention1 != '' | id %in% initial.ranks$ID)
  contestants <- arrange(contestants, period)
  assign(paste0('contestants.', currentClan), contestants)
  
  all.conts <- rbind(all.conts, data.frame(contestants, clan = currentClan))
  
  #Interactions
  interactions <- filter(aggs,
                         aggressor %in% c(initial.ranks$ID, contestants$id),
                         recip %in% c(initial.ranks$ID, contestants$id),
                         Year <= last.year,
                         Year >= first.year) %>%
    semi_join(contestants, by = c('aggressor' = 'id', 'Year' = 'period')) %>%
    semi_join(contestants, by = c('recip' = 'id', 'Year' = 'period')) %>%
    rename(winner = aggressor, loser = recip, period = Year) %>%
    dplyr::select(winner, loser, period) %>%
    arrange(period)
  assign(paste0('interactions.', currentClan), interactions)
  
  all.intx <- rbind(all.intx, filter(aggs,
                                     aggressor %in% c(initial.ranks$ID, contestants$id),
                                     recip %in% c(initial.ranks$ID, contestants$id),
                                     Year <= last.year,
                                     Year >= first.year) %>%
                      semi_join(contestants, by = c('aggressor' = 'id', 'Year' = 'period')) %>%
                      semi_join(contestants, by = c('recip' = 'id', 'Year' = 'period')))
}

###Interactions used to infer ranks:
length(unique(all.intx$aggid)) ##12505
length(unique(all.intx[all.intx$groupcomp != '','aggid'])) ##2966 coalitions
length(unique(all.intx[all.intx$groupcomp != '','aggid']))/length(unique(all.intx$aggid)) ##23.72% coalitions



##Yearly aggression summary
yearly.aggs <- all.intx %>% 
  group_by(clan, Year) %>% 
  summarize(num.aggs = length(aggressor))


ids.per.year <- all.conts %>% 
  group_by(clan, period) %>% 
  summarize(num.ids = length(id))

yearly.aggs <- left_join(yearly.aggs, ids.per.year, by = c('clan', 'Year' = 'period'))
yearly.aggs$aggs.per.id <- yearly.aggs$num.aggs/yearly.aggs$num.ids

##Summary of aggs/id 
mean(yearly.aggs$aggs.per.id)
sd(yearly.aggs$aggs.per.id)
min(yearly.aggs$aggs.per.id)
max(yearly.aggs$aggs.per.id)



save(file = '../2.hyena_data.RData',
     list = c('contestants.hz', 'contestants.north', 'contestants.south', 'contestants.talek',
              'interactions.hz', 'interactions.north', 'interactions.south', 'interactions.talek',
              'initial.ranks.hz', 'initial.ranks.north', 'initial.ranks.south', 'initial.ranks.talek'))


################################################################################
#                       Create coalition networks                               #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################

rm(list = ls())
setwd('~/Documents/Research/socialmob/ranksDynaRank/')
options(stringsAsFactors = FALSE)
library(dplyr)

###Read in data####
aggressions <- read.csv('0.rawdata/tblAggression.csv')
aggressions$date <- as.Date(aggressions$date)
aggressions$Year <- as.numeric(format(aggressions$date, '%Y'))

aggsFull <- aggressions
aggsFull$enteredby <- NULL
aggsFull$doublecheck <- NULL
names(aggsFull) <- c('Clan', 'AggID', 'Session', 'Date', 'Time', 'Agg', 'Recip',
                     'Group','GroupComp', 'Context', 'AltContext', 'B1', 'B2',
                     'R1', 'R2', 'R3', 'Seq','Year')

aggsFull$Agg <- tolower(aggsFull$Agg)
aggsFull$Recip <- tolower(aggsFull$Recip)

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

hps <- read.csv('0.rawdata/tblHyenasPerSession.csv', colClasses = 'character')
sessions <- read.csv('0.rawdata/tblSessions.csv', colClasses = 'character')
hps$date <- left_join(hps, sessions, by = 'session')$date
hps$date <- as.Date(hps$date)
hps$Year <- format(hps$date, '%Y')
names(hps) <- c('Session', 'Hyena', 'FAS', 'FASStart', 'FASStop' ,'FeedingFAS',
                'Follow', 'Tracked', 'Date', 'Year')


##ranks
load('4.ranks.RData')
names(ranks) <- c('Year', 'ID', 'Rank', 'Stan.Rank', 'OldOrder', 'Clan')
ranks[ranks$Clan == 'happy.zebra',]$Clan <- 'hz'
ranks[ranks$Clan == 'serena.s',]$Clan <- 'south'
ranks[ranks$Clan == 'serena.n',]$Clan <- 'north'

coal.summary <- data.frame()

####Calculate coalition strength networks####
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    ids <- yearRanks$ID
    
    ####Coalition network
    cnet <- matrix(nrow = nrow(yearRanks), ncol = nrow(yearRanks), dimnames = list(yearRanks$ID, yearRanks$ID), data = 0)
    coalsTemp <- filter(aggsFull, Year == year, Clan == clan, Group != 'n')
    coalPartners <- cbind(expand.grid(yearRanks$ID, yearRanks$ID), rep(0))
    names(coalPartners) <- c('Focal', 'Alter', 'NumCoals')
    coalPartners$Focal <- as.character(coalPartners$Focal)
    coalPartners$Alter <- as.character(coalPartners$Alter)
    if(!nrow(coalsTemp)){next}
    coal.ids <- c()
    for(row in 1:nrow(coalsTemp)){
      focal <- coalsTemp[row,]$Agg
      for(alter in unique(strsplit(coalsTemp[row,]$GroupComp, ',')[[1]][-1])){
        if(alter %in% yearRanks$ID & alter != focal & focal %in% yearRanks$ID){
          coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,]$NumCoals <- coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,]$NumCoals + 1
          coal.ids <- c(coal.ids, coalsTemp[row,'AggID'])
        }
      }
    }
    coalPartnersPaired <- coalPartners[coalPartners$NumCoals != 0,]
    cnet[as.matrix(coalPartnersPaired[,1:2])] <- coalPartnersPaired[,3]
    
    assign(paste('coalNet', clan, year, sep = '_'), cnet)
    coal.summary <- rbind(coal.summary, data.frame(num.coals = length(unique(coal.ids)),
                                                   clan, year))
  }
}

##
sum(coal.summary$num.coals)  ###1913 coalitions used to measure bond strength

##add metrics to ranks table
ranks$coal_deg <- NA
ranks$coal_top3_deg <- NA

for(row in 1:nrow(ranks)){
  year <- ranks[row,]$Year
  clan <- ranks[row,]$Clan
  focal <- ranks[row,]$ID
  mom <- filter(tblHyenas, ID == focal)$Mom
  
  if(exists(paste('coalNet', clan, year, sep = '_'))){
    coal_net <- get(paste('coalNet', clan, year, sep = '_'))
    obs <- length(unique(filter(hps, Hyena == focal, Year == year)$Session))
    if(nrow(coal_net)){
      ranks[row,]$coal_deg <- sum(coal_net[focal,])
      ranks[row,]$coal_top3_deg <- coal_net[focal,order(coal_net[focal,], decreasing = T)[1:3]] %>% 
        sum()
    } 
  }
}

###Get observations for each individual
hps$Year <- as.numeric(format(hps$Date, '%Y'))
obs_counts <- hps %>% group_by(Hyena, Year) %>% summarize(obs_counts = length(Session))
ranks$obs_counts <- left_join(ranks, obs_counts, by = c('ID' = 'Hyena', 'Year'))$obs_counts


save(ranks, file = '6.ranks_with_coalition_data.RData')
save(list = ls()[starts_with(vars = ls(), 'coalNet')], file = '6.12.coalition_networks.RData')

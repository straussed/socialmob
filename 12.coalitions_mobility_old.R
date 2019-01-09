################################################################################
#                    Coalitions and upward mobility                            #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################

rm(list =ls())
library(dplyr)
options(stringsAsFactors = FALSE)
load('8.ranks_with_rank_change.RData')
load('6.12.coalition_networks.RData')

###Read in data####
aggressions <- read.csv('0.rawdata/tblAggression.csv')
aggressions$date <- as.Date(aggressions$date)
aggressions$Year <- as.numeric(format(aggressions$date, '%Y'))

aggsFull <- aggressions
names(aggsFull) <- c('Clan', 'Session', 'Date', 'Time', 'Agg', 'Recip',
                     'Group','GroupComp', 'Context', 'AltContext', 'B1', 'B2',
                     'R1', 'R2', 'R3', 'Seq', 'Notes', 'Year')


aggsFull$Agg <- tolower(aggsFull$Agg)
aggsFull$Recip <- tolower(aggsFull$Recip)


######Do coalitions increase likelihood of aggressing together up the hierarchy?####

### identify up-hierarchy **triadic** coalitions

full.coals <- filter(aggsFull, GroupComp != '')
full.coals$Year <- as.numeric(format(full.coals$Date, '%Y'))
full.coals$AggRank <- left_join(full.coals, ranks, by = c('Agg' = 'OldOrder', 'Year'))$Stan.Rank
full.coals$RecipRank <- left_join(full.coals, ranks, by = c('Recip' = 'OldOrder', 'Year'))$Stan.Rank

##How many coalitions have more than 2 aggressors?
triadic <- filter(full.coals, lengths(gregexpr(',', full.coals$GroupComp)) == 3,
                  !is.na(AggRank), !is.na(RecipRank))

non_adults <- c()
for(row in 1:nrow(triadic)){
  ids <- strsplit(triadic[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == triadic[row,]$Year,
                          Clan == triadic[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
triadic <- triadic[-non_adults,]
triadic <- unique(triadic[,c('Date', 'Year', 'GroupComp', 'Recip', 'Time', 'Clan')])


more.than.two <- filter(full.coals, lengths(gregexpr(',', full.coals$GroupComp)) > 3,
                        !is.na(AggRank), !is.na(RecipRank))

non_adults <- c()
for(row in 1:nrow(more.than.two)){
  ids <- strsplit(more.than.two[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == more.than.two[row,]$Year,
                          Clan == more.than.two[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
more.than.two <- more.than.two[-non_adults,]
more.than.two <- unique(more.than.two[,c('Date', 'Year', 'GroupComp', 'Recip', 'Time', 'Clan')])




###What proportion of coalitions are triadic?
nrow(triadic)/(nrow(triadic) + nrow(more.than.two))

up.coals <- filter(full.coals, AggRank < RecipRank)
non.up.coals <- filter(full.coals, !is.na(AggRank), !is.na(RecipRank), AggRank > RecipRank)
non.up.coals <- filter(non.up.coals, lengths(gregexpr(',', non.up.coals$GroupComp)) == 3)

#Make sure all coalitions are between adults
non_adults <- c()
for(row in 1:nrow(non.up.coals)){
  ids <- strsplit(non.up.coals[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == non.up.coals[row,]$Year,
                          Clan == non.up.coals[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
non.up.coals <- non.up.coals[-non_adults,]



rev.coals <- semi_join(up.coals, up.coals[which(duplicated(up.coals[,c('GroupComp', 'Date','Time')])),
                                          c('GroupComp', 'Date','Time')], by = c('Date', 'Time', 'GroupComp'))
bridge.coals <- semi_join(full.coals, up.coals, by = c('GroupComp', 'Date','Time')) %>%
  anti_join(rev.coals, by = c('Date', 'Time', 'GroupComp'))

rev.coals <- filter(rev.coals, lengths(gregexpr(',', rev.coals$GroupComp)) == 3)
bridge.coals <- filter(bridge.coals, lengths(gregexpr(',', bridge.coals$GroupComp)) == 3)

#Make sure all coalitions are between adults
non_adults <- c()
for(row in 1:nrow(rev.coals)){
  ids <- strsplit(rev.coals[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == rev.coals[row,]$Year,
                          Clan == rev.coals[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
rev.coals <- rev.coals[-non_adults,]


#Make sure all bridging coalitions are between adults
non_adults <- c()
for(row in 1:nrow(bridge.coals)){
  ids <- strsplit(bridge.coals[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == bridge.coals[row,]$Year,
                          Clan == bridge.coals[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
bridge.coals <- bridge.coals[-non_adults,]


#####
rev.coal.count <- unique(rev.coals[,c('Date', 'Year', 'GroupComp', 'Recip', 'Time', 'Clan')])
non.up.count <- unique(non.up.coals[,c('Date', 'Year', 'GroupComp', 'Recip', 'Time', 'Clan')])

both <- 0
none <- 0
one <- 0

for(row in 1:nrow(rev.coal.count)){
  ids <- strsplit(rev.coal.count[row,]$GroupComp, ',')[[1]][-1]
  
  changers <- sum(filter(ranks, ID %in% ids, Year == rev.coal.count[row,]$Year)$Stan.Rank > 
                    filter(ranks, ID == rev.coal.count[row,]$Recip, Year == rev.coal.count[row,]$Year)$Stan.Rank)
  
  changers.ny <- sum(filter(ranks, ID %in% ids, Year == rev.coal.count[row,]$Year+1)$Stan.Rank > 
                       filter(ranks, ID == rev.coal.count[row,]$Recip, Year == rev.coal.count[row,]$Year+1)$Stan.Rank)
  
  
  if(changers.ny == 0){
    none <- none+1
  }else if(changers.ny == 1){
    one <- one+1
  }else if(changers.ny == 2){
    both <- both+1
  }
}



both/(both+one+none)
one/(both+one+none)

### 68.75% of revolutionary coalitions result in at least 1 individual moving up the hierarchy
(both+one)/(both+one+none)

### 31.25% of revolutionary coalitions result 0 individuals moving up the hierarchy
none/(none+one+both)

##Are revolutionary coalitoins associated with rank reversals? chi-squared test
matrix(nrow =2, ncol = 2,
       dimnames = list(c('reversal', 'no reversal'),
                       c('revolutionary', 'down-hierarchy')),
       data = 0)



##Do 'revolutionary' coalitions involve top allies?###
coal_cont <- matrix(nrow = 2, ncol = 2, dimnames = list(c('top', 'not top'),
                                                        c('up-hierarchy', 'down-hierarchy')),
                    data = 0)

rev.coal.count$CoalCount <- NA
rev.coal.count$Direction <- 'Up'
non.up.count$SumCoals <- NA
for(row in 1:nrow(rev.coal.count)){
  coal_net<- get(paste('coalNet', rev.coal.count[row,]$Clan, rev.coal.count[row,]$Year, sep = '_'))
  ids <- strsplit(rev.coal.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3]){
    coal_cont['top','up-hierarchy'] <- coal_cont['top','up-hierarchy']+1
  }else{
    coal_cont['not top','up-hierarchy'] <- coal_cont['not top','up-hierarchy']+1
  }
  rev.coal.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
  rev.coal.count[row,]$SumCoals <- sum(coal_net[ids,])
}

non.up.count$CoalCount <- NA
non.up.count$Direction <- 'Down'
non.up.count$SumCoals <- NA
for(row in 1:nrow(non.up.count)){
  coal_net<- get(paste('coalNet', non.up.count[row,]$Clan, non.up.count[row,]$Year, sep = '_'))
  ids <- strsplit(non.up.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3]){
    coal_cont['top','down-hierarchy'] <- coal_cont['top','down-hierarchy']+1
  }else{
    coal_cont['not top','down-hierarchy'] <- coal_cont['not top','down-hierarchy']+1
  }
  non.up.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
  non.up.count[row,]$SumCoals <- sum(coal_net[ids,])
}

####top allies engage in up-hierarchy coalitions
#######p = 0.0039
coal_cont
chisq.test(coal_cont)

coal.type <- rbind(non.up.count, rev.coal.count)
coal.type$Direction <- as.factor(coal.type$Direction)


##Using logistf package - strength of coalition predict up- vs down-directed coalition
bond.strength.mod <- logistf::logistf((as.numeric(coal.type$Direction)-1) ~ coal.type$CoalCount)
tot.coal.mod <- logistf::logistf((as.numeric(coal.type$Direction)-1) ~ coal.type$SumCoals)
betas.dyad <- vector(length = 1000)
betas.dyad[1] <- bond.strength.mod$coefficients[2]
betas.total <- vector(length = 1000)
betas.total[1] <- tot.coal.mod$coefficients[2]

coal.type.perm <- coal.type
for(i in 2:1000){
  coal.type.perm <- coal.type.perm %>%
    group_by(Clan, Year)%>%
    mutate(Direction = sample(Direction, replace = FALSE))
  
  bond.strength.perm <- logistf::logistf((as.numeric(coal.type.perm$Direction)-1) ~ coal.type.perm$CoalCount)
  tot.coal.perm <- logistf::logistf((as.numeric(coal.type.perm$Direction)-1) ~ coal.type.perm$SumCoals)
  betas.dyad[i] <- bond.strength.perm$coefficients[2]
  betas.total[i] <- tot.coal.perm$coefficients[2]
}

sum(betas.dyad >= betas.dyad[1])/1000
sum(betas.total >= betas.total[1])/1000



###Plot


##Figure 2

coef(bond.strength.mod)
plogit <- coef(bond.strength.mod)[1] + coef(bond.strength.mod)[2]*(1:25)
pred.up <- data.frame(coal.count = 1:25,
                      predictions = exp(plogit)/(1+exp(plogit)))

# Summarise data to create histogram counts
h = coal.type %>% group_by(Direction) %>%
  mutate(breaks = cut(CoalCount, breaks=seq(0,20,1), labels=seq(0.5,20,1), 
                      include.lowest=TRUE),
         breaks = as.numeric(as.character(breaks))) %>%
  group_by(Direction, breaks) %>% 
  summarise(n = n()) %>%
  mutate(pct = ifelse((as.numeric(Direction)-1)==0, n/sum(n), 1 - n/sum(n))) 

h$Direction <- as.numeric(h$Direction)-1

downs <- coal.type[coal.type$Direction=='Down',]
ups <- coal.type[coal.type$Direction=='Up',]
ups$CoalCount <- ups$CoalCount + runif(nrow(ups), -0.9, -0.1)
downs$CoalCount <- downs$CoalCount + runif(nrow(downs), -0.9, -0.1)

pdf("plots/12.Fig2.pdf",
    height = 4,
    width = 4)
ggplot(data = pred.up, aes(y = predictions, x = coal.count))+
  geom_line(size = 1)+
  geom_segment(data=h, size=4, show.legend=FALSE,
               aes(x=breaks, xend=breaks, y=Direction, yend=pct, colour=factor(Direction)))+
  geom_segment(data=downs, aes(x=CoalCount, xend=CoalCount, y=0, yend=-0.02), size=0.2, colour="grey30") +
  geom_segment(data=ups, aes(x=CoalCount, xend=CoalCount, y=1, yend=1.02), size=0.2, colour="grey30") +
  theme_classic()+
  ylab('Probability coalition is directed up the hierarchy')+
  xlab('Strength of bond between coalition partners')+
  scale_color_manual(values = c('black', 'black'))

dev.off()


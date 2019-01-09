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
library(ggplot2)
options(stringsAsFactors = FALSE)
load('8.ranks_with_rank_change.RData')
load('6.12.coalition_networks.RData')
set.seed(1989)

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

#Exclude final years for each clan because ranks can't change in final year because
#of validation rule
aggsFull <- filter(aggsFull, !(Clan == 'talek' & Year == 2015),
                   !(Clan == 'south' & Year %in% c(2015, 2016)),
                   !(Clan == 'north' & Year == 2015),
                   !(Clan == 'hz' & Year == 2015))


######Do coalitions increase likelihood of aggressing together up the hierarchy?####

### identify up-hierarchy **triadic** coalitions

full.coals <- filter(aggsFull, GroupComp != '')
full.coals$Year <- as.numeric(format(full.coals$Date, '%Y'))
full.coals$AggRank <- left_join(full.coals, ranks, by = c('Agg' = 'OldOrder', 'Year'))$Stan.Rank
full.coals$RecipRank <- left_join(full.coals, ranks, by = c('Recip' = 'OldOrder', 'Year'))$Stan.Rank

##How many coalitions have more than 2 aggressors?
triadic <- filter(full.coals, lengths(gregexpr(',', full.coals$GroupComp)) == 3,
                  !is.na(AggRank), !is.na(RecipRank))

##Only include triadic coalitions where both members are adult females
non_adults <- c()
for(row in 1:nrow(triadic)){
  ids <- strsplit(triadic[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == triadic[row,]$Year,
                          Clan == triadic[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
triadic <- triadic[-non_adults,]

##Polyadic coalitions with more than two allies
more.than.two <- filter(full.coals, lengths(gregexpr(',', full.coals$GroupComp)) > 3,
                        !is.na(AggRank), !is.na(RecipRank))
##Only include polyadic coalitions where both members are adult females
non_adults <- c()
for(row in 1:nrow(more.than.two)){
  ids <- strsplit(more.than.two[row,]$GroupComp, split = ',')[[1]][-1]
  if(!all(ids %in% filter(ranks, Year == more.than.two[row,]$Year,
                          Clan == more.than.two[row,]$Clan)$ID))
    non_adults <- c(non_adults, row)
}
more.than.two <- more.than.two[-non_adults,]

###What proportion of coalitions are triadic?
length(unique(triadic$AggID))/(length(unique(triadic$AggID)) + length(unique(more.than.two$AggID)))

###How many coalitions have more than 2 aggressors?
dyadic <- filter(aggsFull, GroupComp == '')
dyadic$Year <- as.numeric(format(dyadic$Date, '%Y'))
dyadic$AggRank <- left_join(dyadic, ranks, by = c('Agg' = 'OldOrder', 'Year'))$Stan.Rank
dyadic$RecipRank <- left_join(dyadic, ranks, by = c('Recip' = 'OldOrder', 'Year'))$Stan.Rank

up.dyadic <- filter(dyadic, AggRank < RecipRank)
down.dyadic <- filter(dyadic, RecipRank < AggRank)



up.coals <- filter(triadic, AggRank < RecipRank)
down.coals <- filter(triadic, AggRank > RecipRank)


##table 
up.tab <- table(up.coals$AggID)
down.tab <- table(down.coals$AggID)
#both aggressors are lower ranked
rev.coals <- filter(up.coals, AggID %in% names(up.tab[up.tab == 2]))
bridge.coals <- filter(up.coals, AggID %in% names(up.tab[up.tab == 1]))
down.coals <- filter(down.coals, AggID %in% names(down.tab[down.tab == 2]))


###Are coalitions more likely in up- or down-hierarchy aggressions?
dyadic.triadic <- matrix(nrow =2, ncol = 2,
                         dimnames = list(c('down', 'up'),
                                         c('dyadic', 'triadic')),
                         data = c(length(unique(down.dyadic$AggID)), length(unique(up.dyadic$AggID)),
                                  length(unique(down.coals$AggID)), length(unique(rev.coals$AggID))))

chisq.test(dyadic.triadic)
###Proportion of down hierarchy interactions that are triadic/dyadic + triadic
dyadic.triadic[1,2]/sum(dyadic.triadic[1,]) ##5.8%

###Proportion of up hierarchy interactions that are triadic/dyadic + triadic
dyadic.triadic[2,2]/sum(dyadic.triadic[2,]) ##4.2 % 

sum(dyadic.triadic[2,])/sum(dyadic.triadic)


## Confirm all triadic coalitions are accounted for
length(unique(rev.coals$AggID)) + 
  length(unique(bridge.coals$AggID)) + 
  length(unique(down.coals$AggID)) == length(unique(triadic$AggID))

rev.coal.count <- unique(rev.coals[,c('Date', 'AggID', 'Year', 'GroupComp', 'Recip', 'Clan')])
bridge.coal.count <- unique(bridge.coals[,c('Date', 'AggID', 'Year', 'GroupComp', 'Recip', 'Clan')])
down.coal.count <- unique(down.coals[,c('Date', 'AggID', 'Year', 'GroupComp', 'Recip', 'Clan')])

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

### 65.63% of revolutionary coalitions result in at least 1 individual moving up the hierarchy
(both+one)/(both+one+none)

### 34.37% of revolutionary coalitions result 0 individuals moving up the hierarchy
none/(none+one+both)

none.down <- 0
one.down <- 0
both.down <- 0
for(row in 1:nrow(down.coal.count)){
  ids <- strsplit(down.coal.count[row,]$GroupComp, ',')[[1]][-1]
  
  changers <- sum(filter(ranks, ID %in% ids, Year == down.coal.count[row,]$Year)$Stan.Rank < 
                    filter(ranks, ID == down.coal.count[row,]$Recip, Year == down.coal.count[row,]$Year)$Stan.Rank)
  
  changers.ny <- sum(filter(ranks, ID %in% ids, Year == down.coal.count[row,]$Year+1)$Stan.Rank < 
                       filter(ranks, ID == down.coal.count[row,]$Recip, Year == down.coal.count[row,]$Year+1)$Stan.Rank)
  
  
  if(changers.ny == 0){
    none.down <- none.down+1
  }else if(changers.ny == 1){
    one.down <- one.down+1
  }else if(changers.ny == 2){
    both.down <- both.down+1
  }
}

### 0.43% of revolutionary coalitions result in at least 1 individual moving up the hierarchy
(both.down+one.down)/(both.down+one.down+none.down)

### 99.57% of revolutionary coalitions result 0 individuals moving up the hierarchy
none.down/(none.down+one.down+both.down)


##Are revolutionary coalitions associated with rank reversals? chi-squared test
rev.chisq <- matrix(nrow =2, ncol = 2,
       dimnames = list(c('reversal', 'no reversal'),
                       c('revolutionary', 'down-hierarchy')),
       data = c(both+one, none, both.down+one.down, none.down))

chisq.test(rev.chisq)


##Are revolutionary coalitions associated with rank reversals? chi-squared test
rev.chisq <- matrix(nrow =2, ncol = 2,
                    dimnames = list(c('reversal', 'no reversal'),
                                    c('revolutionary', 'down-hierarchy')),
                    data = c(both+one, none, both.down+one.down, none.down))

chisq.test(rev.chisq)



##Do 'revolutionary' coalitions involve top allies?###
coal_cont <- matrix(nrow = 2, ncol = 2, dimnames = list(c('top', 'not top'),
                                                        c('up-hierarchy', 'down-hierarchy')),
                    data = 0)

rev.coal.count$CoalCount <- NA
rev.coal.count$Direction <- 'Up'
rev.coal.count$SumCoals <- NA
for(row in 1:nrow(rev.coal.count)){
  coal_net<- get(paste('coalNet', rev.coal.count[row,]$Clan, rev.coal.count[row,]$Year, sep = '_'))
  ids <- strsplit(rev.coal.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3] &
     ids[1] %in% names(sort(coal_net[ids[2],], decreasing = TRUE))[1:3]){
    coal_cont['top','up-hierarchy'] <- coal_cont['top','up-hierarchy']+1
  }else{
    coal_cont['not top','up-hierarchy'] <- coal_cont['not top','up-hierarchy']+1
  }
  rev.coal.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
  rev.coal.count[row,]$SumCoals <- sum(coal_net[ids,])
}

down.coal.count$CoalCount <- NA
down.coal.count$Direction <- 'Down'
down.coal.count$SumCoals <- NA
for(row in 1:nrow(down.coal.count)){
  coal_net<- get(paste('coalNet', down.coal.count[row,]$Clan, down.coal.count[row,]$Year, sep = '_'))
  ids <- strsplit(down.coal.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3] &
     ids[1] %in% names(sort(coal_net[ids[2],], decreasing = TRUE))[1:3]){
    coal_cont['top','down-hierarchy'] <- coal_cont['top','down-hierarchy']+1
  }else{
    coal_cont['not top','down-hierarchy'] <- coal_cont['not top','down-hierarchy']+1
  }
  down.coal.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
  down.coal.count[row,]$SumCoals <- sum(coal_net[ids,])
}

####top allies engage in up-hierarchy coalitions
#######p = 0.0005
coal_cont
chisq.test(coal_cont)

coal.type <- rbind(down.coal.count, rev.coal.count)
coal.type$Direction <- as.factor(coal.type$Direction)

###Ensure dyads appear in a consisent order (i.e., A,B always rather than some A,B and others B,A)
alphabetize <- function(x){
  ids <- strsplit(x, ',')[[1]][-1]
  return(paste0(',',paste(ids[order(ids)], collapse = ','), ','))
}

coal.type$GroupComp <- sapply(coal.type$GroupComp, alphabetize)

bond.strength.mod <- lme4::glmer(data = coal.type, Direction ~ CoalCount + (1|Clan), family = 'binomial')
summary(bond.strength.mod)

tot.coal.mod <- lme4::glmer(data = coal.type, Direction ~ SumCoals + (1|Clan), family = 'binomial')
summary(tot.coal.mod)
aic <- MuMIn::AICc(bond.strength.mod, tot.coal.mod)
aic$AICc[1] - aic$AICc[2]


betas.dyad <- vector(length = 1000)
betas.dyad[1] <- bond.strength.mod@beta[2]
betas.total <- vector(length = 1000)
betas.total[1] <- tot.coal.mod@beta[2]

coal.type.perm <- coal.type
for(i in 2:1000){
  coal.type.perm <- coal.type.perm %>%
    group_by(Clan, Year)%>%
    mutate(Direction = sample(Direction, replace = FALSE))
  
  bond.strength.perm <- lme4::glmer(data = coal.type.perm, Direction ~ CoalCount + (1|Clan), family = 'binomial')
  tot.coal.perm <- lme4::glmer(data = coal.type.perm, Direction ~ SumCoals + (1|Clan), family = 'binomial')
  betas.dyad[i] <- bond.strength.perm@beta[2]
  betas.total[i] <- tot.coal.perm@beta[2]
}

sum(betas.dyad >= betas.dyad[1])/1000
sum(betas.total >= betas.total[1])/1000


###Plot


##Figure 2 - excluded from final paper for space reasons

coef(bond.strength.mod)
plogit <- bond.strength.mod@beta[1] + bond.strength.mod@beta[2]*(1:25)
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
  scale_color_manual(values = c('black', 'black'))+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9))

dev.off()


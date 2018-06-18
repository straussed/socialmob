##################################################################
##                        Eli Strauss                           ##
##      Predict magnitude of rank change by allies and kin      ##
##                       Nov 28th, 2017                         ##
##################################################################
options(stringsAsFactors = FALSE)
library(magrittr)
library(tidyverse)
library(igraph)
library(grid)
library(lme4)
library(logistf)

###Read in data####
aggressions <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblAggression.csv')
aggressions$date <- as.Date(aggressions$date)
aggressions$Year <- as.numeric(format(aggressions$date, '%Y'))

aggsFull <- aggressions
names(aggsFull) <- c('Clan', 'Session', 'Date', 'Time', 'Agg', 'Recip',
                     'Group','GroupComp', 'Context', 'AltContext', 'B1', 'B2',
                     'R1', 'R2', 'R3', 'Seq', 'Notes','Year')

#demographic information
hyenas <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblHyenas.csv')
hyenas[hyenas$clan == 'serena s',]$clan <- 'south'
hyenas[hyenas$clan == 'serena n',]$clan <- 'north'
hyenas[hyenas$clan == 'happy zebra',]$clan <- 'hz'
hyenas$birthdate <- as.Date(hyenas$birthdate)
hyenas$disappeared <- as.Date(hyenas$disappeared)

tblHyenas <- hyenas
names(tblHyenas) <- c('ID', 'Last.Updated', 'SampleID', 'Eartage', 'Name', 'PrevID',
                   'Sex', 'AgeClass', 'Status', 'FirstSeen', 'DenGrad', 'Disappeared',
                   'Mom', 'Birthdate','NumLittermates',
                   'LitRank', 'ArrivedDen', 'LeaveDen', 'Fate', 'MortalitySource',
                   'DeathDate', 'Weaned', 'Clan', 'Park', 'Notes')

hps <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblHyenasPerSession.csv', colClasses = 'character')
sessions <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblSessions.csv', colClasses = 'character')
hps$date <- left_join(hps, sessions, by = 'session')$date
hps$date <- as.Date(hps$date)
hps$Year <- format(hps$date, '%Y')
names(hps) <- c('Session', 'Hyena', 'FAS', 'FASStart', 'FASStop' ,'FeedingFAS',
                'Follow', 'Tracked', 'Date', 'Year')




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
    for(row in 1:nrow(coalsTemp)){
      focal <- coalsTemp[row,]$Agg
      for(alter in unique(strsplit(coalsTemp[row,]$GroupComp, ',')[[1]][-1])){
        if(alter %in% yearRanks$ID & alter != focal & focal %in% yearRanks$ID){
          coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,]$NumCoals <- coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,]$NumCoals + 1
        }
      }
    }
    coalPartnersPaired <- coalPartners[coalPartners$NumCoals != 0,]
    cnet[as.matrix(coalPartnersPaired[,1:2])] <- coalPartnersPaired[,3]
    
    assign(paste('coalNet', clan, year, sep = '_'), cnet)
  }
}


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


####Examine coalition degree in top3 partners wtd by obs####
ranks.perm <- na.omit(dplyr::select(ranks, ID, Clan, Year, RankDiffAbs, coal_top3_deg, obs_counts))
top3_mod <- ranks.perm %>% lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + (1|ID) + (1|Clan/Year) + log(offset(obs_counts)))
top3_poly <- ranks.perm %>% lme4::lmer(formula = RankDiffAbs ~ poly(coal_top3_deg,2) + (1|ID) + (1|Clan/Year) + log(offset(obs_counts)))
summary(top3_poly)
AIC(top3_mod, top3_poly)

#removing 3 largest rank changes doesn't change the effect
filter(ranks, abs(RankDiffAbs) <10) %>% lme4::lmer(formula = RankDiffAbs ~ poly(coal_top3_deg,2) + (1|ID) + (1|Clan/Year) + log(offset(obs_counts))) %>% summary()

####Analysis
predictions.list <- list()
betas.list <- list()

for(i in 1:999){
  ranks.perm$RankDiffAbs <- ranks.perm %>% group_by(Clan, Year) %>%
    sample_frac(replace = FALSE) %>% 
    ungroup() %>% dplyr::select(RankDiffAbs) %>% .[[1]]
  
  mod.perm <- ranks.perm %>% lme4::lmer(formula = RankDiffAbs ~ poly(coal_top3_deg, 2) + (1|ID) + (1|Clan/Year) + log(offset(obs_counts)))
  
  betas.list[[i]] <- data.frame(Intercept = mod.perm@beta[1], coal_top3_deg = mod.perm@beta[2], coal_top3_poly = mod.perm@beta[3], offset = mod.perm@beta[4])
  
  if(i <= 100){
    predictions.list[[i]] <- data.frame(rc = predict(mod.perm, type = 'response'), 
                                        coal_top3_deg = ranks.perm$coal_top3_deg,
                                        obs_counts = ranks.perm$obs_counts)
  }
}
predictions <- do.call(rbind, predictions.list)
betas <- do.call(rbind, betas.list)




tblHyenas$BirthYear <- as.numeric(format(tblHyenas$Birthdate, '%Y'))
tblHyenas$Survive2 <- ifelse(is.na(tblHyenas$Disappeared) | 
                               tblHyenas$Disappeared - tblHyenas$Birthdate > (365*2),
                             TRUE,
                             FALSE)

left_join(ranks, tblHyenas, by = c('ID' = 'Mom', 'Year' = 'BirthYear')) %>%
  group_by(ID, Year) %>% summarize(ARS = sum(!is.na(Name)), 
                                   Surv2 = sum(as.numeric(Survive2))) %>%
  left_join(ranks, ., by = c('ID', 'Year')) -> ranks

ranks[is.na(ranks$Surv2),]$Surv2 <- 0

ranks %>% group_by(ID) %>% summarize(mean.rank = mean(stan.rank), 
                                     total.cubs = sum(Surv2),
                                     mean.cubs = mean(Surv2),
                                     mean.ars = mean(ARS)) -> lrs

lrs <- filter(lrs, ID %in% filter(tblHyenas, !is.na(Disappeared))$ID)
lrs$repro.lifespan <- as.numeric(left_join(lrs, tblHyenas, by = 'ID')$Disappeared - 
                                   left_join(lrs, tblHyenas, by = 'ID')$Birthdate - 2*365)/365

survive_to_4 <- filter(tblHyenas, !is.na(Disappeared), Disappeared - Birthdate >= 4*365)

nrow(filter(lrs, ID %in% survive_to_4$ID))
pois.quad.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ poly(mean.rank, 2), family = 'poisson')
summary(pois.quad.lrs)

pois.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ mean.rank, family = 'poisson')
summary(pois.lrs)
pois.exp.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ I(exp(1)^mean.rank), family = 'poisson')
summary(pois.exp.lrs)

##Mean lifetime reproductive success
mean(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)
sd(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)/sqrt(nrow(filter(lrs, ID %in% survive_to_4$ID)))
min(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)
max(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)

AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)
AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[1,2] - AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[3,2]
AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[2,2] - AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[3,2]



quad.ars <- glm(data = filter(lrs, ID %in% survive_to_4$ID), 
                formula = total.cubs ~ poly(mean.rank, 2) + offset(log(repro.lifespan)),
                family = poisson(link = "log"))
summary(quad.ars)

ars <- glm(data = filter(lrs, ID %in% survive_to_4$ID), 
           formula = total.cubs ~ mean.rank + offset(log(repro.lifespan)),
           family = poisson(link = "log"))
summary(ars)
exp.ars <- glm(data = filter(lrs, ID %in% survive_to_4$ID), 
               formula = total.cubs ~ I(exp(1)^mean.rank) + offset(log(repro.lifespan)),
               family = poisson(link = "log"))

AIC(quad.ars, ars, exp.ars)

##predict consequence of rank change on lrs based on models

predicted.lrs.final.rank <- 
  exp(predict(pois.exp.lrs, newdata = data.frame(mean.rank = ranks$stan.rank)))
predicted.lrs.initial.rank <- 
  exp(predict(pois.exp.lrs, 
              newdata = data.frame(mean.rank = ranks$stan.rank - ranks$RankDiff)))

ranks$lrs.diff <- predicted.lrs.final.rank - predicted.lrs.initial.rank


####Amplifications of small rank differences###
matriarchs <- c('kb', 'dj', '03', 'coch', '40',
                'rbc', 'wafl', 'digs', 'shrm',
                'clov', 'java', 'coel', 'pike')
ranks$matriline <- NA
for(matriarch in matriarchs){
  assign(paste0('matriarch_', matriarch), matriarch)
  for(i in 1:20){
    assign(paste0('matriarch_', matriarch), 
           unique(c(get(paste0('matriarch_', matriarch)),
                    filter(tblHyenas, Mom %in% get(paste0('matriarch_', matriarch)))$ID)))
  }
  ranks[ranks$ID %in% get(paste0('matriarch_', matriarch)),]$matriline <- matriarch
}
ranks[is.na(ranks$matriline),]$matriline <- 'other'
ranks[ranks$matriline == 'coch',]$matriline <- '03'

first_year <- filter(ranks, Clan == 'talek') %>% 
  group_by(ID) %>% summarize(first_year = min(Year), 
                             first_rank = Rank[which.min(Year)],
                             matriline = matriline[1]) 
ribbon <- filter(ranks, Clan == 'talek') %>%
  group_by(matriline, Year) %>%
  summarize(ymin = min(Rank), ymax = max(Rank), max.stan = max(stan.rank),
            mean.stan = mean(stan.rank))

ribbon[ribbon$matriline == '40' & ribbon$Year > 2007,c('ymin', 'ymax')] <- NA


matriline_labels <- data.frame(x.start = rep(1987), 
                               matriline = c('kb', 'dj', '03', '40', 'other'),
                               y.start = c(0.7, 3.3, 5.5, 8, 15))
matriline_labels <- arrange(matriline_labels, matriline)

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <-  c("#4477AA", viridis::viridis(5)[5], viridis::viridis(5)[4], 
           lighten(desat(viridis::viridis(1), 0.5), 1.6), "#CC6677")

matriline.colors <- ggplot(filter(ranks, Clan == 'talek'), 
                           aes(x = Year, y = Rank, group = ID, col = matriline)) +
  geom_ribbon(data = ribbon, aes(x = Year, ymin = ymin - 0.5, ymax = ymax +0.5,
                                 group = matriline, fill = matriline),
              inherit.aes = FALSE)+
  geom_ribbon(data = filter(ranks, Clan == 'talek'), aes(x = Year, ymin = Rank - 0.5, ymax = Rank+0.5,
                                                         group = ID,  fill = matriline),
              inherit.aes = FALSE) +
  geom_line(size = 0.5) +
  theme_classic() +
  ylim(55, 0) + 
  theme(legend.position = 'none') +
  geom_label(data = matriline_labels,
             aes(x = x.start, y = y.start, label = matriline, group = matriline, fill = matriline),
             label.size = NA,
             label.padding = unit(0.2, 'lines'),
             fontface = 2,
             color = c('black', 'black', 'black', 'black', 'black'))+
  geom_point(data = data.frame(hex = c('hex', 'hex'), x = c(2005, 2014), y = c(19, 40)),
             aes(x = x, y = y), inherit.aes = FALSE,
             fill = desat(cols[2], sat = 0.5), shape = 23, size =2.5, stroke = 1.4,
             col = darken(cols[2]))+
  scale_color_manual(values = darken(cols)) + 
  scale_fill_manual(values = desat(cols, sat = 0.5))+
  ggtitle('a)')


##Track difference between females over time##
matriline.pairs <- data.frame(id1 = c('kb', 'dj', '03'),
                              id2 = c('dj', '03', '40'))


get.descendants <- function(mom){
  if(nrow(filter(tblHyenas, Mom == mom))){
    her.descendants <- c(mom)
    for(kid in filter(tblHyenas, Mom == mom)$ID){
      her.descendants <- c(her.descendants, get.descendants(kid))
    }
    return(her.descendants)
  }else{
    return(mom)
  }
}

# #set coch's mom to be 03 in tblHyenas
# tblHyenas[tblHyenas$ID == 'coch',]$Mom <- '03'
#set coch's mom to be 03 in tblHyenas
tblHyenas[tblHyenas$ID == 'coch',]$Mom <- '03'

descendants <- list('03' = get.descendants('03'),
                    '40' = get.descendants('40'),
                    dj = get.descendants('dj'),
                    kb = get.descendants('kb'))

yearly.data.no.change <- list()
ids.that.change <- unique(filter(ranks, RankChange != 'None')$ID)
counter <- 1
for(row in 1:nrow(matriline.pairs)){
  for(year in sort(unique(ranks$Year))){
    id1.desc <- descendants[[matriline.pairs[row,1]]]
    id2.desc <- descendants[[matriline.pairs[row,2]]]
    
    mean.stan.rank.diff <- filter(ranks, Year == year, ID %in% id1.desc,
                                  !ID %in% ids.that.change)$stan.rank %>%
      mean() - 
      mean(filter(ranks, Year == year, ID %in% id2.desc,
                  !ID %in% ids.that.change)$stan.rank)
    
    mean.abs.rank.diff <- filter(ranks, Year == year, ID %in% id2.desc,
                                 !ID %in% ids.that.change)$Rank %>%
      mean() - 
      mean(filter(ranks, Year == year, ID %in% id1.desc,
                  !ID %in% ids.that.change)$Rank)
    
    
    yearly.data.no.change[[counter]] <- data.frame(id1 = matriline.pairs[[row,1]],
                                                   id2 = matriline.pairs[[row,2]],
                                                   mean.stan.rank.diff,
                                                   mean.abs.rank.diff,
                                                   year)
    counter <- counter+1
  }
}
rank.diff.over.time.no.change <- do.call(rbind, yearly.data.no.change)

rank.diff.over.time.no.change$id1 <- factor(rank.diff.over.time.no.change$id1, 
                                            levels = c('kb', 'dj', '03'))

rank.diff.over.time.no.change$id2 <- factor(rank.diff.over.time.no.change$id2, 
                                            levels = c('dj', '03', '40'))

mean.rank.diff <- 
  ggplot(rank.diff.over.time.no.change, 
         aes(x = year, y = mean.abs.rank.diff, 
             col = id1, fill = id2))+
  geom_smooth(method = 'lm', se = FALSE, col = 'black', 
              aes(x = year, y = mean.abs.rank.diff), inherit.aes = F)+
  geom_jitter(shape = 21, stroke = 1.5, size = 2, width = 1, height = 1)+
  theme_classic()+
  ylab('Difference in mean rank of adjacent matrilines')+
  xlab('Year')+
  scale_fill_manual(name = 'Difference between\nadjacent matrilines',
                    labels = c('kb-dj', 'dj-03', '03-40'), 
                    values = cols[c(3,1,2)])+
  scale_color_manual(name = 'Difference between\nadjacent matrilines',
                     labels = c('kb-dj', 'dj-03', '03-40'),
                     values = cols[c(4,3,1)])+
  theme(legend.position = c(0.35, 0.65),
        legend.title.align = 0.5)+
  ggtitle('b)')


######Do coalitions increase likelihood of aggressing together up the hierarchy?####

### identify up-hierarchy **dyadic** coalitions


full.coals <- filter(aggsFull, GroupComp != '')
full.coals$Year <- as.numeric(format(full.coals$Date, '%Y'))
full.coals$AggRank <- left_join(full.coals, ranks, by = c('Agg' = 'OldOrder', 'Year'))$stan.rank
full.coals$RecipRank <- left_join(full.coals, ranks, by = c('Recip' = 'OldOrder', 'Year'))$stan.rank

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
  
  changers <- sum(filter(ranks, ID %in% ids, Year == rev.coal.count[row,]$Year)$stan.rank > 
                    filter(ranks, ID == rev.coal.count[row,]$Recip, Year == rev.coal.count[row,]$Year)$stan.rank)
  
  changers.ny <- sum(filter(ranks, ID %in% ids, Year == rev.coal.count[row,]$Year+1)$stan.rank > 
                       filter(ranks, ID == rev.coal.count[row,]$Recip, Year == rev.coal.count[row,]$Year+1)$stan.rank)
  
  
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
(both+one)/(both+one+none)
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
for(row in 1:nrow(rev.coal.count)){
  coal_net<- get(paste('coalNet', rev.coal.count[row,]$Clan, rev.coal.count[row,]$Year, sep = '_'))
  ids <- strsplit(rev.coal.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3]){
    coal_cont['top','up-hierarchy'] <- coal_cont['top','up-hierarchy']+1
  }else{
    coal_cont['not top','up-hierarchy'] <- coal_cont['not top','up-hierarchy']+1
  }
  rev.coal.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
}

non.up.count$CoalCount <- NA
non.up.count$Direction <- 'Down'
for(row in 1:nrow(non.up.count)){
  coal_net<- get(paste('coalNet', non.up.count[row,]$Clan, non.up.count[row,]$Year, sep = '_'))
  ids <- strsplit(non.up.count[row,]$GroupComp, ',')[[1]][-1]
  if(ids[2] %in% names(sort(coal_net[ids[1],], decreasing = TRUE))[1:3]){
    coal_cont['top','down-hierarchy'] <- coal_cont['top','down-hierarchy']+1
  }else{
    coal_cont['not top','down-hierarchy'] <- coal_cont['not top','down-hierarchy']+1
  }
  non.up.count[row,]$CoalCount <- coal_net[ids[1], ids[2]]
}

chisq.test(coal_cont)

coal.type <- rbind(non.up.count, rev.coal.count)
coal.type$Direction <- as.factor(coal.type$Direction)

##Using logistf package
bond.strength.mod <- logistf::logistf((as.numeric(coal.type$Direction)-1) ~ coal.type$CoalCount)



#####Plots!######

##Figure 1
talek <- ggplot(data = filter(ranks, Clan == 'talek'), aes(y = Rank, x = Year)) + 
  ylim(52,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black')+
  ggtitle('Talek Clan')

north <- ggplot(data = filter(ranks, Clan == 'north'), aes(y = Rank, x = Year)) + 
  ylim(22,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black') + 
  ggtitle('Serena North Clan')

south <- ggplot(data = filter(ranks, Clan == 'south'), aes(y = Rank, x = Year)) + 
  ylim(25,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black')+
  ggtitle('Serena South Clan')

hz <- ggplot(data = filter(ranks, Clan == 'hz'), aes(y = Rank, x = Year)) + 
  ylim(20,0)+
  theme_classic() + 
  geom_line(aes(y = Rank, x = Year, group = ID), col = 'black')+
  ggtitle('Happy Zebra Clan')

pdf(file = "~/Documents/Figures/Rank changes/Figure 1 - rank reversals.pdf",
    width = 6,
    height = 4)
multiplot(talek, north, south, hz, layout = matrix(c(1,1,1,1,1,1,2,3,4), nrow = 3, byrow = TRUE))
dev.off()
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

pdf("~/Documents/Figures/Rank changes/Figure 2 - bond strength predicts coalition direction.pdf",
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



##Figure 3

betas.summary = rbind(data.frame(estimate = names(betas), 
                                 ci = apply(X = betas, MARGIN = 2, FUN = quantile, probs = 0.025)),
                      data.frame(estimate = names(betas), 
                                 ci = apply(X = betas, MARGIN = 2, FUN = quantile, probs = 0.975)))

observed <- data.frame(estimate = names(betas),
                       mean = top3_poly@beta)

observed$estimate <- factor(observed$estimate, levels = c('coal_top3_deg', 'coal_top3_poly', 'Intercept', 
                                                          'offset'),
                            labels = c('Coalitions with top allies',
                                       'Coalitions with top allies (squared)',
                                       'Intercept',
                                       'Offset'))

betas.summary$estimate <- factor(betas.summary$estimate, levels = c('coal_top3_deg', 'coal_top3_poly', 'Intercept', 
                                                                    'offset'),
                                 labels = c('Coalitions with top allies',
                                            'Coalitions with top allies (squared)',
                                            'Intercept',
                                            'Offset'))

inset <- ggplot(data = betas.summary, aes(x = ci, y = estimate))+
  geom_line(size = 1) + 
  geom_point(data = observed, 
             aes(x = mean, y = estimate), size = 4, col = 'black',
             shape = 21) + 
  theme_bw() + 
  geom_vline(aes(xintercept = 0), col = 'firebrick', lty = 2) + 
  ylab('')+
  xlab('')
# theme(axis.text.y = element_text(size = 16))+
# theme(axis.text.x = element_text(size = 16))


top3_ranks <- data.frame(rc = predict(top3_poly, type = 'response'),
                         coal_top3_deg = na.omit(ranks[,c('coal_top3_deg', 'obs_counts')])$coal_top3_deg,
                         obs_counts = na.omit(ranks[,c('coal_top3_deg', 'obs_counts')])$obs_counts,
                         category = 'Model Predictions')
obs <- data.frame(rc = ranks$RankDiffAbs,
                  coal_top3_deg = ranks$coal_top3_deg,
                  obs_counts = ranks$obs_counts,
                  category = 'Observed')

predictions$category <- 'Permuted'

fig3_ranks <- rbind(predictions,
                    obs,
                    top3_ranks)

fig3_ranks$category <- factor(fig3_ranks$category, 
                              levels = c('Model Predictions',
                                         'Observed',
                                         'Permuted'))


main <- ggplot(data=fig3_ranks, aes(x = coal_top3_deg, y = rc, col = category))+
  geom_jitter(size = 2)+
  theme_classic() + 
  ylab('Rank change due to rank reversal') +
  xlab('Strength of coalitionary ties with top partners') + 
  theme(legend.position = c(0.6,0.2))+
  ylim(c(-10, 25))+
  scale_color_manual(values = c('firebrick', 'dodgerblue', 'grey'),
                     guide = guide_legend(title = ''))


vp = viewport(width = 1, height = 0.2, x = 0.5, y = 0.9)

pdf("~/Documents/Figures/Rank changes/Figure 3 - rank change due to rank reversal.pdf",
    height= 6,
    width = 8)
print(main)
print(inset, vp = vp)
dev.off()



##Figure 4a
pois.exp.pred <- predict.glm(pois.exp.lrs, se.fit = TRUE, type = "response")
lrs.mean <- mean(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)

lrs.plot <- ggplot(filter(lrs, ID %in% survive_to_4$ID), aes(x = mean.rank, y = total.cubs, label = ID))+
  geom_point(size =2) +
  geom_line(aes(y = pois.exp.pred$fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), size = 1.2, col = 'firebrick')+
  geom_line(aes(y = pois.exp.pred$fit + 2*pois.exp.pred$se.fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), lty = 2, size = 1, col = 'firebrick')+
  geom_line(aes(y = pois.exp.pred$fit - 2*pois.exp.pred$se.fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), lty = 2, size = 1, col = 'firebrick')+
  theme_classic()+
  xlab('Mean lifetime standardized rank') + 
  ylab('Total cubs produced')+
  ggtitle('a)')

##Figure 4b

predicted.lrs.final.rank <- 
  exp(predict(pois.exp.lrs, newdata = data.frame(mean.rank = ranks$stan.rank)))
predicted.lrs.initial.rank <- 
  exp(predict(pois.exp.lrs, 
              newdata = data.frame(mean.rank = ranks$stan.rank - ranks$RankDiff)))

ranks$lrs.ratio <- predicted.lrs.final.rank/predicted.lrs.initial.rank
ranks$lrs.diff <- predicted.lrs.final.rank - predicted.lrs.initial.rank

predict.consequences <- ggplot(filter(ranks, lrs.diff != 0), aes(x = lrs.diff/lrs.mean))+
  geom_histogram(fill = 'firebrick')+
  theme_classic()+
  xlab('Predicted change in LRS due to rank change')+
  ylab('')+
  ggtitle('b)')

## Figure 4c
exp.change.lrs <- ggplot(filter(ranks, lrs.ratio != 1), aes(x = stan.rank, y = RankDiffAbs, fill = lrs.ratio, size = exp(abs(1-lrs.ratio))))+
  geom_jitter(shape = 21, col = 'grey80')+
  theme_classic()+
  geom_hline(yintercept = 0, lty = 2)+
  xlab('Rank after change')+
  ylab('Number of positions moved') +
  scale_fill_gradientn(colors = c("firebrick", "white","dodgerblue"),
                       values = scales::rescale(c(0.5, 1 , 3)),
                       guide = "colorbar", limits = c(0.5, 3))+
  guides(size = FALSE, 
         fill= guide_colorbar(title = 'New LRS ÷ Old LRS', direction = 'horizontal',
                              title.position = 'top', title.hjust = 0.5, barwidth = 7, nbin = 100))+
  theme(legend.position = c(0.15, 0.75))+
  ggtitle('c)')

pdf("~/Documents/Figures/Rank changes/Figure 4 - predicted effects of rank change.pdf",
    height = 8,
    width = 8)
multiplot(lrs.plot, predict.consequences,
          exp.change.lrs, layout = matrix(c(1,2,3,3), byrow = TRUE, nrow = 2))
dev.off()

#Figure 5
pdf("~/Documents/Figures/Rank changes/Figure 5 - intergenerational effects.pdf",
    height = 8.5,
    width = 7)
gridExtra::grid.arrange(matriline.colors, mean.rank.diff, nrow = 2, heights = c(5,3))
dev.off()

###Descriptives###
####Number of aggressions
all.aggs <- inner_join(aggsFull, ranks, by = c('Agg' = 'ID', 'Year')) %>% 
  inner_join(ranks, by = c('Recip' = 'ID', 'Year'))
nrow(filter(all.aggs, GroupComp != ''))
nrow(filter(all.aggs, GroupComp != ''))/nrow(all.aggs)

####Number of females
length(unique(ranks$ID))

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
sum(table(ranks$RankChange)[c(1,3)])/sum(table(ranks$RankChange))

###Number of aggressions involving females#####
aggs_between_females <-  inner_join(aggsFull, ranks, by = c('Agg' = 'ID', 'Year')) %>% 
  inner_join(ranks, by = c('Recip' = 'ID', 'Year'))
nrow(aggs_between_females)
nrow(filter(aggs_between_females, GroupComp != ''))


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

##################################################################
##                        Eli Strauss                           ##
##      Predict magnitude of rank change by allies and kin      ##
##                       Nov 28th, 2017                         ##
##################################################################
#library(igraph)
#library(mlogit)
#library(ggplot2)
#source("/mnt/home/straus46/RankChangeModels/Fisibase/fisibasetidy/ReadTidyData.R")
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")
library(magrittr)
library(tidyverse)
library(igraph)

########Calculate AI degree####################
#Create ai networks for each clan and year
for(clan in unique(ranks$Clan)){
  for(year in unique(ranks[ranks$Clan == clan,]$Year)){
    ids <- ranks[ranks$Clan == clan & ranks$Year == year,]$ID
    ###session data
    sess.year <- filter(hps, Hyena %in% ids, format(Date, '%Y') == year)
    if(tibdim(sess.year)){
      sess.year <- add_column(sess.year, Clan = rep(clan, tibdim(sess.year)), Year = rep(year, tibdim(sess.year)))
      ai.count.net <- matrix(data = 0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
      for(sess.num in unique(sess.year$Session)){
        ids.sess <- as.matrix(expand.grid(sess.year[sess.year$Session == sess.num,]$Hyena, sess.year[sess.year$Session == sess.num,]$Hyena))
        ai.count.net[ids.sess] <- ai.count.net[ids.sess]+1
      }
      ai.net <- ai.count.net
      for(a in rownames(ai.count.net)){
        for(b in rownames(ai.count.net)){
          ai.net[a,b] <- ai.count.net[a,b]/(ai.count.net[a,a] + ai.count.net[b,b] - ai.count.net[a,b])
        }
      }
      diag(ai.net) <- 0
      assign(paste('ai_net', year, clan, sep = '_'), ai.net)
      assign(paste('ai_count_net', year, clan, sep = '_'), ai.count.net)
    }else{
      assign(paste('ai_net', year, clan, sep = '_'), NULL)
    }
  }
}

####Calculate coalition degree####
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    ####Coalition network
    cnet <- matrix(nrow = tibdim(yearRanks), ncol = tibdim(yearRanks), dimnames = list(yearRanks$ID, yearRanks$ID), data = 0)
    coalsTemp <- filter(aggsFull, Year == year, Clan == clan, Group != 'n')
    coalPartners <- cbind(expand.grid(yearRanks$ID, yearRanks$ID), rep(0))
    names(coalPartners) <- c('Focal', 'Alter', 'NumCoals')
    coalPartners$Focal <- as.character(coalPartners$Focal)
    coalPartners$Alter <- as.character(coalPartners$Alter)
    if(!tibdim(coalsTemp)){next}
    for(row in 1:tibdim(coalsTemp)){
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

##add metrics to ranks tibble
ranks <- add_column(ranks, coal_deg = rep(NA, tibdim(ranks)),
                    coal_top3_deg = rep(NA, tibdim(ranks)),
                    coal_deg_wtd = rep(NA),
                    coal_top3_wtd = rep(NA),
                    ally_ranks = rep(NA),
                    mom_change = rep(NA),
                    num_kin = rep(NA),
                    ally_kin = rep(NA),
                    ally_kin_ai = NA,
                    coals_with_kin = rep(NA),
                    ai_deg = rep(NA),
                    ai_top3_deg = rep(NA),
                    coal_gini = rep(NA)
)


for(row in 1:tibdim(ranks)){
  year <- ranks[row,]$Year
  clan <- ranks[row,]$Clan
  focal <- ranks[row,]$ID
  mom <- filter(tblHyenas, ID == focal)$Mom
  if(length(mom) && mom %in% filter(ranks, Year == year, Clan == clan)$ID){
    ranks[row,]$mom_change <- filter(ranks, Year == year, Clan == clan, ID == mom)$RankDiffAbs
  }
  if(length(mom) && !is.na(mom)){
    if(!any(is.na(filter(tblHyenas, ID %in% filter(ranks, Year == year, Clan == clan)$ID)$Mom))){
      ##count number of half sibs, moms, daughters in hierarchy
      filter(tblHyenas, ID %in% filter(ranks, Year == year, Clan == clan)$ID & ID != focal) %>% 
        filter(Mom %in% c(mom, focal) | ID == mom) %>% .[['ID']] -> kin
      ranks[row,]$num_kin <- length(kin)
    }else{kin <- NA}
  }else{kin <- NA}
  if(exists(paste('coalNet', clan, year, sep = '_')) & exists(paste('ai_net', year, clan, sep = '_')) & !is.null(get(paste('ai_net', year, clan, sep = '_')))){
    coal_net <- get(paste('coalNet', clan, year, sep = '_'))
    ai_net <- get(paste('ai_net', year, clan, sep = '_'))
    coal_net_wtd <- coal_net/ai_net
    coal_net_wtd[is.nan(coal_net_wtd)] <- 0
    if(tibdim(coal_net)){
      ranks[row,]$coal_deg <- sum(coal_net[focal,])
      ranks[row,]$coal_deg_wtd <- sum(coal_net_wtd[focal,])
      ranks[row,]$coal_top3_deg <- coal_net[focal,order(coal_net[focal,], decreasing = T)[1:3]] %>% 
        sum()
      coal_net[focal,order(coal_net[focal,], decreasing = T)[1:3]] %>% names() -> allies
      coal_net[focal,order(ai_net[focal,], decreasing = T)[1:3]] %>% names() -> allies_ai
      ranks[row,]$ally_ranks <- sum(filter(ranks, ID %in% allies, Year == year)$Rank)
      ranks[row,]$ally_kin <- ifelse(all(is.na(kin)), NA,
                                     sum(allies %in% kin)/3)
      ranks[row,]$coal_gini <- ineq::Gini(coal_net[focal,])
      ranks[row,]$ally_kin_ai <- ifelse(all(is.na(kin)), NA,
                                     sum(allies_ai %in% kin)/3)
      ranks[row,]$coals_with_kin <- ifelse(all(is.na(kin)), NA,
                                           sum(coal_net[kin,focal]))
      ranks[row,]$coal_top3_wtd <- coal_net_wtd[focal,order(coal_net_wtd[focal,], decreasing = T)[1:3]] %>% 
        sum() 
    } 
  }
}


for(row in 1:tibdim(ranks)){
  year <- ranks[row,]$Year
  clan <- ranks[row,]$Clan
  focal <- ranks[row,]$ID
  ai_net <- get(paste('ai_net', year, clan, sep = '_'))
  if(tibdim(ai_net)){
    ranks[row,]$ai_deg <- sum(ai_net[focal,])
    ranks[row,]$ai_top3_deg <- ai_net[focal,order(ai_net[focal,], decreasing = T)[1:3]] %>% 
      sum() 
  }
}





################Analyze coalitions and AIs together for whole hierarchy ###############

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

coal_perm_stacked <- tibble()
#par(mfrow = c(2,2))
coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = ranks_mod$coal_deg,
                                      coal_top_deg = ranks_mod$coal_top3_deg,
                                      ai_all_deg = ranks_mod$ai_deg,
                                      ai_top_deg = ranks_mod$ai_top3_deg,
                                      rank_change = ranks_mod$RankDiffAbs,
                                      network = 'Observed'))

iterations <- 1000

rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

both_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_deg + ai_deg + (1|ID))
both_top3_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + ai_top3_deg + (1|ID))


obs_coef <- as_tibble(cbind(both_deg_intercept = both_mod %>% .@beta %>% .[1],
                            both_coal_deg = both_mod %>% .@beta %>% .[2],
                            both_ai_deg = both_mod %>% .@beta %>% .[3],
                            both_top3_intercept = both_top3_mod %>% .@beta %>% .[1],
                            both_coal_top3 = both_top3_mod %>% .@beta %>% .[2],
                            both_ai_top3 = both_top3_mod %>% .@beta %>% .[3]))

obs_coef <- as_tibble(cbind(both_deg_intercept = both_mod %>% .@beta %>% .[1],
                            both_coal_deg = both_mod %>% .@beta %>% .[2],
                            both_ai_deg = both_mod %>% .@beta %>% .[3],
                            both_top3_intercept = both_top3_mod %>% .@beta %>% .[1],
                            both_coal_top3 = both_top3_mod %>% .@beta %>% .[2],
                            both_ai_top3 = both_top3_mod %>% .@beta %>% .[3]))

perm_coef <- as_tibble(cbind(both_deg_intercept = rep(NA, iterations),
                             both_coal_deg = rep(NA, iterations),
                             both_ai_deg = rep(NA, iterations),
                             both_top3_intercept = rep(NA, iterations),
                             both_coal_top3 = rep(NA, iterations),
                             both_ai_top3 = rep(NA, iterations)))


ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks_mod$Clan, ranks_mod$Year)

for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    dplyr::select(RankDiffAbs) %>% .[[1]] ->
    ranks.perm$RankDiffAbs
  
  ranks.perm %>%
    lme4::lmer(formula = RankDiffAbs ~ coal_deg + ai_deg + (1|ID)) -> all_perm_mod
  all_perm_mod%>% 
    .@beta %>% .[1] -> perm_coef[i,1]
  all_perm_mod%>% 
    .@beta %>% .[2] -> perm_coef[i,2]
  all_perm_mod%>% 
    .@beta %>% .[3] -> perm_coef[i,3]
  
  ranks.perm %>%
    lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + ai_top3_deg + (1|ID)) -> top_perm_mod
  top_perm_mod %>%
    .@beta %>% .[1] -> perm_coef[i,4]
  top_perm_mod %>%
    .@beta %>% .[2] -> perm_coef[i,5]
  top_perm_mod %>%
    .@beta %>% .[3] -> perm_coef[i,6]
  
  rank_perm_change[,i] <- as.character(ranks.perm$RankDiffAbs)
}

coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = rep(ranks_mod$coal_deg, iterations),
                                      coal_top_deg = rep(ranks_mod$coal_top3_deg, iterations),
                                      ai_all_deg = rep(ranks_mod$ai_deg, iterations),
                                      ai_top_deg = rep(ranks_mod$ai_top3_deg, iterations),
                                      rank_change = as.numeric(as.vector(rank_perm_change)),
                                      network = 'Permuted'))



par(mfrow = c(ceiling(length(perm_coef)/2), 2))
for(param in names(perm_coef)){
  range <- c(min(perm_coef[param])-.005, max(perm_coef[param], obs_coef[param]) + .005)
  mar.default <- c(5,4,4,2) + 0.1
  par(family = 'Trebuchet MS', mar = mar.default + c(0, 1, 0, 0))
  d <- density(perm_coef[[param]])
  p = ifelse(2*tibdim(which(perm_coef[[param]] >= obs_coef[[param]]))/iterations > 1,
             round(2*tibdim(which(perm_coef[[param]] <= obs_coef[[param]]))/iterations, digits = 3),
             round(2*tibdim(which(perm_coef[[param]] >= obs_coef[[param]]))/iterations, digits = 3))
  coal_top3_plot <- plot(d, lwd = 6, col = 'grey87',
                         main = paste0('Effect size = ', round(obs_coef[param], digits = 3), ', ', 'p = ', p),
                         xlab = param,
                         cex.lab=2, cex.axis=2, cex.main=2,
                         xlim = range)
  polygon(d, col = 'grey87', border = 'grey87')
  abline(v = obs_coef[param], col = 'dodgerblue1', lty=2, lwd=3)
}



facet_labels = c(High = 'High ranking individuals', Low = 'Low ranking individuals')
colors = c('dodgerblue1', 'grey87')

coal_perm_stacked$rank_change_category <- ifelse(coal_perm_stacked$rank_change == 0,'None', 
                                                 ifelse(coal_perm_stacked$rank_change < 0,'Down','Up'))
coal_perm_stacked$rank_change_category %<>% factor(., levels = c('Down', 'None', 'Up'))


###Coalitions###
ggplot(data = filter(coal_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = coal_top_deg, y = rank_change, col = network)) + 
  geom_jitter(col = alpha('grey', .1))+
  geom_jitter(data = filter(coal_perm_stacked, network == 'Observed'), aes(x = coal_top_deg, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$both_top3_intercept, slope = perm_coef$both_coal_top3, col = 'grey', alpha = 0.02, size = 1) +
  geom_abline(intercept = obs_coef$both_top3_intercept, slope = obs_coef$both_coal_top3, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('Number of coalitions with top partners')+
  ylab('Amount of rank change')


top_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change_category, y = coal_top_deg, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with top allies')+
  coord_cartesian(ylim= c(0,35))

top_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold', size = 12,family = 'Trebuchet MS'),
        axis.title = element_text(size = 14,family = 'Trebuchet MS'),
        axis.text = element_text(size = 12, family = 'Trebuchet MS'),
        legend.text = element_text(size = 12, family = 'Trebuchet MS'),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank())

###AI#####
ggplot(data = filter(coal_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = ai_top_deg, y = rank_change, col = network)) + 
  geom_jitter(col = alpha('grey', .1))+
  geom_jitter(data = filter(coal_perm_stacked, network == 'Observed'), aes(x = ai_top_deg, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$both_top3_intercept, slope = perm_coef$both_ai_top3, col = 'grey', alpha = 0.02, size = 1) +
  geom_abline(intercept = obs_coef$both_top3_intercept, slope = obs_coef$both_ai_top3, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('Sum of association indices with top partners') + 
  ylab('Amount of rank change')
  



########Relationship between kin and rank change######

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

kin_perm_stacked <- tibble()
#par(mfrow = c(2,2))
kin_perm_stacked <- bind_rows(kin_perm_stacked,
                               tibble(mom_change = ranks_mod$mom_change,
                                      num_kin = ranks_mod$num_kin,
                                      ally_kin = ranks_mod$ally_kin,
                                      coals_with_kin = ranks_mod$coals_with_kin,
                                      rank_change = ranks_mod$RankDiffAbs,
                                      network = 'Observed'))

iterations <- 1000

rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

#kin_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ mom_change + num_kin + ally_kin + coals_with_kin + (1|ID))
kin_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ mom_change + num_kin + (1|ID))

obs_coef <- as_tibble(cbind(kin_intercept = kin_mod %>% .@beta %>% .[1],
                            mom_change = kin_mod %>% .@beta %>% .[2],
                            num_kin = kin_mod %>% .@beta %>% .[3],
                            ally_kin = kin_mod %>% .@beta %>% .[4],
                            coals_with_kin = kin_mod%>% .@beta %>% .[5]))
perm_coef <- as_tibble(cbind(kin_intercept = rep(NA, iterations),
                            mom_change = rep(NA, iterations),
                            num_kin = rep(NA, iterations),
                            ally_kin = rep(NA, iterations),
                            coals_with_kin = rep(NA, iterations)))

ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks_mod$Clan, ranks_mod$Year)

for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    dplyr::select(RankDiffAbs) %>% .[[1]] ->
    ranks.perm$RankDiffAbs
  
  ranks.perm %>%
    #lme4::lmer(formula = RankDiffAbs ~ mom_change + num_kin + ally_kin + coals_with_kin + (1|ID)) -> kin_perm_mod
    lme4::lmer(formula = RankDiffAbs ~ mom_change + num_kin + (1|ID)) -> kin_perm_mod
  for(col in 1:length(perm_coef)){
    kin_perm_mod%>% 
      .@beta %>% .[col] -> perm_coef[i,col]
  }
  
  rank_perm_change[,i] <- as.character(ranks.perm$RankDiffAbs)
}

kin_perm_stacked <- bind_rows(kin_perm_stacked,
                              tibble(mom_change = rep(ranks_mod$mom_change,iterations),
                                     num_kin = rep(ranks_mod$num_kin,iterations),
                                     ally_kin = rep(ranks_mod$ally_kin,iterations),
                                     coals_with_kin = rep(ranks_mod$coals_with_kin,iterations),
                                     rank_change = as.numeric(as.vector(rank_perm_change)),
                                     network = 'Permuted'))



par(mfrow = c(ceiling(length(perm_coef)/2), 2))
for(param in names(perm_coef)){
  range <- c(min(perm_coef[param])-.005, max(perm_coef[param], obs_coef[param]) + .005)
  mar.default <- c(5,4,4,2) + 0.1
  par(family = 'Trebuchet MS', mar = mar.default + c(0, 1, 0, 0))
  d <- density(perm_coef[[param]])
  p = ifelse(2*tibdim(which(perm_coef[[param]] >= obs_coef[[param]]))/iterations > 1,
             round(2*tibdim(which(perm_coef[[param]] <= obs_coef[[param]]))/iterations, digits = 3),
             round(2*tibdim(which(perm_coef[[param]] >= obs_coef[[param]]))/iterations, digits = 3))
  coal_top3_plot <- plot(d, lwd = 6, col = 'grey87',
                         main = paste0('Effect size = ', round(obs_coef[param], digits = 3), ', ', 'p = ', p),
                         xlab = param,
                         cex.lab=2, cex.axis=2, cex.main=2,
                         xlim = range)
  polygon(d, col = 'grey87', border = 'grey87')
  abline(v = obs_coef[param], col = 'dodgerblue1', lty=2, lwd=6)
}


###Mom Change###
ggplot(data = filter(kin_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = mom_change, y = rank_change)) + 
  geom_jitter(col = alpha('grey', .2))+
  geom_jitter(data = filter(kin_perm_stacked, network == 'Observed'), aes(x = mom_change, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$kin_intercept, slope = perm_coef$mom_change, col = 'grey', alpha = 0.05, size = 1) +
  geom_abline(intercept = obs_coef$kin_intercept, slope = obs_coef$mom_change, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('Rank change of mother')+
  ylab('Rank change of daughter')


###Number of kin###
ggplot(data = filter(kin_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = num_kin, y = rank_change)) + 
  geom_jitter(col = alpha('grey', .2))+
  geom_jitter(data = filter(kin_perm_stacked, network == 'Observed'), aes(x = num_kin, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$kin_intercept, slope = perm_coef$num_kin, col = 'grey', alpha = 0.05, size = 1) +
  geom_abline(intercept = obs_coef$kin_intercept, slope = obs_coef$num_kin, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('Number of adult close kin (mothers, daughters, and maternal sisters)')+
  ylab('Amount of rank change')

ranks$clan_size <- NA
clan_size <- left_join(ranks, 
                             ranks %>% group_by(Clan, Year) %>% summarize(clan_size = length(Year)),
                             by = c('Clan', 'Year'))$clan_size
ranks$clan_size <- clan_size

kin_coals <- na.omit(data.frame(ID = ranks$ID, Clan = ranks$Clan, Year = ranks$Year, RankDiffAbs = ranks$RankDiffAbs, kin = ranks$coals_with_kin, non_kin = (ranks$coal_deg - ranks$coals_with_kin), 
                                prop_kin= ranks$num_kin/clan_size))

chisq <- chisq.test(kin_coals$kin, p = kin_coals$prop_kin)
chisq$p.value

kin_coals <- filter(kin_coals, kin + non_kin > 0)


Xsq <- apply(X = kin_coals[,c('kin', 'non_kin', 'prop_kin')], FUN = function(x)(chisq.test(unlist(x[1:2]), p = unlist(c(x[3], (1-x[3]))))$statistic), MARGIN = 1)
kin_coals$Xsq <- Xsq

plot(kin_coals$Xsq, kin_coals$RankDiffAbs)
summary(lm(data  = kin_coals, RankDiffAbs ~ Xsq))

##########Predicting rank change yearly#####
yearly_rank_change <- ranks %>% group_by(Year, Clan) %>% summarize(total_change = sum(abs(RankDiffAbs))/2, coal_skew = moments::skewness(coal_deg, na.rm = T), clan_size = length(coal_top3_deg)) %>% ungroup()

yearly_rank_change <- yearly_rank_change[order(yearly_rank_change$total_change),]

plot(yearly_rank_change$total_change ~ yearly_rank_change$coal_skew)
par(mfrow = c(1,1))


ranks_yearly <- left_join(ranks, yearly_rank_change, by = c('Clan', 'Year'))
ranks_yearly <- ranks_yearly[order(ranks_yearly$total_change),]
ranks_yearly$Clan_Year <- paste(ranks_yearly$Clan, ranks_yearly$Year, sep = ' ')
ranks_yearly$Clan_Year <- factor(ranks_yearly$Clan_Year, levels = unique(ranks_yearly$Clan_Year))


summary(glm(data= yearly_rank_change, total_change ~ coal_skew, family = 'poisson'))


ggplot(data = ranks_yearly, aes(x = coal_deg, y = Clan_Year, fill = total_change)) + 
  geom_joy(rel_min_height = 0.01)+
  theme_classic()+
  scale_fill_continuous(guide = guide_legend(title = '# Rank changes'), low = viridis(6)[1], high = viridis(6)[6])  + 
  xlab('Number of coalitions by each individual') + 
  theme(legend.position = c(.8, .2))+
  ylab('(Fewer changes)                      Clan and Year                 (More changes)')
  

# #yearly_rank_change$skew <- NA
# #for(row in 1:nrow(yearly_rank_change)){
#   year = yearly_rank_change[row,]$Year
#   clan = yearly_rank_change[row,]$Clan
#   coal_net <- get(paste('coalNet', clan, year, sep = '_'))
#   yearly_rank_change[row,]$skew <- moments::skewness(coal_net[,])
#     #coal_net %>% 
#     #graph_from_adjacency_matrix(mode = 'undirected', weighted = T) %>% 
#     #strength() %>% moments::skewness()
# }

summary(glm(data= yearly_rank_change, total_change ~ coal_skew, family = 'poisson'))


####Standard deviation of vertex strength
ggplot(data = perm_data, aes(x = coal_skew, y = total_change))+
  geom_jitter(col = alpha('grey', .01))+
  geom_abline(data = perm_coef, aes(intercept = intercept, slope = coal_skew), col = alpha('grey', .1))+
  geom_jitter(col = 'dodgerblue', data=yearly_rank_change, aes(x = coal_skew, y = total_change))+
  geom_abline(data = obs_coef, aes(intercept = intercept, slope = coal_skew), col = 'dodgerblue', size = 1.5)+
  theme_classic()+
  ylab('Number of dyadic rank relationships changed per network')+
  xlab('Skew in # relationships per individual')

p = 2*length(which(perm_coef$coal_skew > obs_coef$coal_skew))/iterations





####Example of permutation####

obs <- coef(lm(data = ranks, coal_deg ~ Rank))
perm <- list()
coal_perm <- matrix(NA, nrow=nrow(ranks), ncol=100)

ranks.perm <- ranks
for(i in 1:100){
  ranks.perm %>% 
    group_by(Clan, Year) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    dplyr::select(coal_deg) %>% .[[1]] ->
    ranks.perm$coal_deg
  perm[[i]] <- coef(lm(data = ranks.perm, coal_deg ~ Rank))
  
  coal_perm[,i] <- ranks.perm$coal_deg
}
perm <- do.call(rbind, perm)



p <- ggplot(data = ranks, aes(x = Rank, y = coal_deg)) + 
  geom_point(col = 'dodgerblue')+
  geom_abline(intercept = obs[1], slope = obs[2], col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('')+
  ylab('')

ggsave(plot = p, filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/observed.png'),
       device = 'png', dpi = 300, width = 5, height = 5, units = 'in')

p <- ggplot(data = data.frame(Rank = rep(ranks$Rank, i), coal_deg = coal_deg), aes(x = Rank, y = coal_deg)) + 
  geom_blank()+
  theme_classic()+
  xlab('')+
  ylab('')
ggsave(p, filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/permutation_0.png'),
       device = 'png', dpi = 300, width = 5, height = 5, units = 'in')

for(i in 1:20){
  coal_deg <- c()
  for(ii in 1:i){
    coal_deg <- c(coal_deg, coal_perm[,ii])
  }
  #png(filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/permutation_', i, '.png'),res = 300, width = 4000, height = 4000)
  p <- ggplot(data = data.frame(Rank = rep(ranks$Rank, i), coal_deg = coal_deg), aes(x = Rank, y = coal_deg)) + 
    geom_point(col = alpha('grey', .1))+
    geom_abline(intercept = perm[1:i,1], slope = perm[1:i,2], col = alpha('grey', .2), size = 2) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_classic()+
    xlab('')+
    ylab('')
  #name <- ifelse(i < 10, paste0('0', i), i)
  ggsave(plot = p, filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/permutation_', i, '.png'),
        device = 'png', dpi = 300, width = 5, height = 5, units = 'in')
        
  #dev.off()
}

p <- ggplot(data = data.frame(Rank = rep(ranks$Rank, i), coal_deg = coal_deg), aes(x = Rank, y = coal_deg)) + 
  geom_point(col = alpha('grey', .1))+
  geom_abline(intercept = perm[1:i,1], slope = perm[1:i,2], col = alpha('grey', .2), size = 2) +
  geom_abline(intercept = obs[1], slope = obs[2], col = 'dodgerblue', size = 1.5) +
  geom_point(data = ranks, aes(x = Rank, y = coal_deg), col = 'dodgerblue')+
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  xlab('')+
  ylab('')
ggsave(plot = p, filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/permutation_21.png'),
       device = 'png', dpi = 300, width = 5, height = 5, units = 'in')


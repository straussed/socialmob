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
library(grid)

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

####Calculate coalition degree (instead of rate)####
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    ids <- yearRanks$ID
    ##Rate not included in analysis. Low rankers have inflated rates due to 
    ##low #s of observation
    
    # sess.year <- filter(hps, Hyena %in% ids, format(Date, '%Y') == year)
    # time_matrix <- matrix(data = 0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
    # for(sess.num in unique(sess.year$Session)){
    #   ids.sess <- as.matrix(expand.grid(sess.year[sess.year$Session == sess.num,]$Hyena, sess.year[sess.year$Session == sess.num,]$Hyena))
    #   time <- as.numeric(filter(sessions, Session == sess.num)$Stop) - as.numeric(filter(sessions, Session == sess.num)$Start)
    #   if(!is.na(time) & time >= 0) time_matrix[ids.sess] <- time_matrix[ids.sess]+time
    # }
    ##Convert time matrix to hours
    time_matrix <- time_matrix/(60*60)
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
    
    #divide by time spent to gether for rates
    #cnet[1:length(ids), 1:length(ids)] <- as.vector(cnet) / as.vector(time_matrix)
    #cnet[is.nan(cnet)] <- 0
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


yearly_rank_change$density <- NA
yearly_rank_change$num_coals <- NA
yearly_rank_change$modularity <- NA

for(clan in unique(yearly_rank_change$Clan)){
  for(year in unique(filter(yearly_rank_change, Clan == clan)$Year)){
    coal_net <- get(paste('coalNet', clan, year, sep = '_')) %>% graph_from_adjacency_matrix(mode = 'undirected', weighted = TRUE)
    if(!ecount(coal_net)){next}
    sum(edge_attr(coal_net)$weight) ->
      yearly_rank_change[yearly_rank_change$Clan == clan & yearly_rank_change$Year == year,]$num_coals
    graph.density(coal_net) ->
      yearly_rank_change[yearly_rank_change$Clan == clan & yearly_rank_change$Year == year,]$density
    
    if(ecount(coal_net))
    coal_net_connected <- delete_vertices(coal_net, which(degree(coal_net) == 0))
    cluster_fast_greedy(coal_net_connected) %>% modularity() ->
      yearly_rank_change[yearly_rank_change$Clan == clan & yearly_rank_change$Year == year,]$modularity
  }
}

plot(yearly_rank_change$total_change ~ yearly_rank_change$modularity)

summary(glm(yearly_rank_change, formula = total_change ~ modularity, family = 'poisson'))


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

for(i in 1:100){
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
ggsave(plot = p, filename = paste0('~/Documents/Presentations/BEACON Talk 2017/permutation/permutation_101.png'),
       device = 'png', dpi = 300, width = 5, height = 5, units = 'in')





#####Look at extreme cases######

old_ranks <- ranks[,c('Rank', 'Year', 'Clan', 'IDold')]
old_ranks$RankChange <- left_join(old_ranks, ranks, by = c('IDold' = 'ID', 'Clan', 'Year'))$RankChange

fall_far <- filter(ranks, RankDiffAbs <= -5)
fall_from_top <- filter(old_ranks, RankChange == 'Down', Rank < 4) %>% semi_join(ranks, ., by = c('Year', 'Clan', 'ID' = 'IDold'))
extreme_down <- unique(rbind(fall_far, fall_from_top))
                      

rise_far <- filter(ranks, RankDiffAbs >= 5)
rise_to_top <- filter(ranks, Rank < 4, RankChange == 'Up')
extreme_up <- unique(rbind(rise_far, rise_to_top))

top <- ggplot(data = extreme_up, aes(x = coal_top3_deg, y = ai_top3_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down, aes(x = coal_top3_deg, y = ai_top3_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year of Change')+
  xlab('Coalition degree')+
  ylab('AI degree')

all <- ggplot(data = extreme_up, aes(x = coal_deg, y = ai_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down, aes(x = coal_deg, y = ai_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year of Change')+
  xlab('Coalition degree')+
  ylab('AI degree')

multiplot(all, top)



extreme_up_previous <- semi_join(ranks, mutate(extreme_up, Prev_Year = Year - 1), by = c('ID', 'Clan', 'Year' = 'Prev_Year'))

extreme_down_previous <- semi_join(ranks, mutate(extreme_down, Prev_Year = Year - 1), by = c('ID', 'Clan', 'Year' = 'Prev_Year'))

previous_top <- ggplot(data = extreme_up_previous, aes(x = coal_top3_deg, y = ai_top3_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down_previous, aes(x = coal_top3_deg, y = ai_top3_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year before Change')+
  xlab('Coalition degree')+
  ylab('AI degree')

previous_all <- ggplot(data = extreme_up_previous, aes(x = coal_deg, y = ai_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down_previous, aes(x = coal_deg, y = ai_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year before Change')+
  xlab('Coalition degree')+
  ylab('AI degree')


extreme_up_next <- semi_join(ranks, mutate(extreme_up, Prev_Year = Year + 1), by = c('ID', 'Clan', 'Year' = 'Prev_Year'))

extreme_down_next <- semi_join(ranks, mutate(extreme_down, Prev_Year = Year + 1), by = c('ID', 'Clan', 'Year' = 'Prev_Year'))

next_top <- ggplot(data = extreme_up_next, aes(x = coal_top3_deg, y = ai_top3_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down_next, aes(x = coal_top3_deg, y = ai_top3_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year after Change')+
  xlab('Coalition degree')+
  ylab('AI degree')

next_all <- ggplot(data = extreme_up_next, aes(x = coal_deg, y = ai_deg))+
  geom_point(col = 'red') + 
  geom_point(data = extreme_down_next, aes(x = coal_deg, y = ai_deg), col = 'blue')+
  theme_classic()+
  xlim(0,40)+
  ggtitle('Year after Change')+
  xlab('Coalition degree')+
  ylab('AI degree')



multiplot(previous_top, top, next_top)


####Fitness####
tblHyenas$Survive2 <- FALSE 
tblHyenas$Survive2 <- ifelse(is.na(tblHyenas$Disappeared), TRUE, FALSE)
tblHyenas[tblHyenas$Survive2 == FALSE,]$Survive2 <- with(tblHyenas[tblHyenas$Survive2 == FALSE,], 
                                                         ifelse(!is.na(Birthdate) & Disappeared - Birthdate >= (365*2), TRUE, FALSE))
tblHyenas$BirthYear <- format(tblHyenas$Birthdate, '%Y')
yearly_cubs <- tblHyenas %>% group_by(BirthYear, Mom) %>% summarize(YearlyCubs = sum(Survive2))


rise_to_top$yearly.cubs.after.change <- NA
for(row in 1:nrow(rise_to_top)){
  last_year <- filter(tblHyenas, ID == rise_to_top[row,]$ID)$LastSeen %>%
    format('%Y') %>% as.numeric()
  rise_to_top[row,]$yearly.cubs.after.change <- 
    sum(filter(yearly_cubs, 
               BirthYear <= last_year, 
               Mom == rise_to_top[row,]$ID)$YearlyCubs) / 
    (last_year - rise_to_top[row,]$Year)  
    
  
}


fall_from_top$yearly.cubs.after.change <- NA
for(row in 1:nrow(fall_from_top)){
  last_year <- filter(tblHyenas, ID == fall_from_top[row,]$ID)$LastSeen %>%
    format('%Y') %>% as.numeric()
  fall_from_top[row,]$yearly.cubs.after.change <- 
    sum(filter(yearly_cubs, 
               BirthYear <= last_year, 
               Mom == fall_from_top[row,]$ID)$YearlyCubs) / 
    (last_year - fall_from_top[row,]$Year)
}


yearly_rs <- rbind(data.frame(ID = rise_to_top$ID,
                              yearly.cubs.after.change = rise_to_top$yearly.cubs.after.change,
                              RankChange = 'Up'),
                   data.frame(ID = fall_from_top$ID,
                              yearly.cubs.after.change = fall_from_top$yearly.cubs.after.change,
                              RankChange = 'Down'))

plot(yearly_rs$yearly.cubs.after.change ~ factor(yearly_rs$RankChange))

###Survival### -- Summary: RankChange has no effect on survival
library(survival)
library(survminer)
##Deal with two 'Both' conditions
ranks[ranks$RankChange == 'Both' & ranks$ID == 'peep',]$RankChange <- 'Up'
ranks[ranks$RankChange == 'Both' & ranks$ID == 'lg',]$RankChange <- 'None'
survival_df <- data.frame(ID = ranks$ID, Year = ranks$Year, Rank = ranks$Rank,
                          RankChange = ranks$RankChange,
                          RankDiff = ranks$RankDiffAbs,
                          Disappeared = left_join(ranks, tblHyenas, by = 'ID')$Disappeared,
                          Birthdate = left_join(ranks, tblHyenas, by = 'ID')$Birthdate,
                          LastSeen = left_join(ranks, tblHyenas, by = 'ID')$LastSeen,
                          FirstSeen = left_join(ranks, tblHyenas, by = 'ID')$FirstSeen,
                          type = 'right',
                          time = NA,
                          event = NA)

survival_df$RankChange <- factor(survival_df$RankChange, levels = c('Down', 'None', 'Up'))

for(row in 1:nrow(survival_df)){
  if(is.na(survival_df$Birthdate[row])){ ##Birthdate missing
    if(is.na(survival_df$Disappeared[row])){
      survival_df$time[row] <- survival_df$LastSeen[row] - survival_df$FirstSeen[row]
      survival_df$event[row] <- 0
    }else{
      survival_df$time[row] <- survival_df$Disappeared[row] - survival_df$FirstSeen[row]
      survival_df$event[row] <- 1
    }
  }else{###Birthday present
    if(is.na(survival_df$Disappeared[row])){
      survival_df$time[row] <- survival_df$LastSeen[row] - survival_df$Birthdate[row]
      survival_df$event[row] <- 0
    }else{
      survival_df$time[row] <- survival_df$Disappeared[row] - survival_df$Birthdate[row]
      survival_df$event[row] <- 1
    }
  }
}

surv.fit <- coxph(Surv(time = survival_df$time, 
                       event = survival_df$event) ~
                    survival_df$RankChange)

summary(surv.fit)

ggsurvplot(survfit(Surv(time, event) ~ RankChange, data = survival_df), conf.int = T, pval = T, data = survival_df, risk.table = F,
           palette = c('dodgerblue4', 'dodgerblue', 'darkorange'), size = 2)



####LRS for animals who have died##### - Summary: no difference between upmovers and downmovers

full_life <- filter(tblHyenas, !is.na(Birthdate) & !is.na(Disappeared), 
                    ID %in% ranks$ID)

full_life$rank.change <- NA
full_life$cubs.produced <- NA
for(row in 1:nrow(full_life)){
  full_life$cubs.produced[row] <- nrow(filter(tblHyenas, Mom == full_life$ID[row]))
  net.rank.change <- sum(filter(ranks, ID == full_life$ID[row])$RankDiffAbs)
  full_life$rank.change[row] <- ifelse(net.rank.change > 0, 'Up', ifelse(net.rank.change < 0, 'Down', 'None'))
}

plot(full_life$cubs.produced ~ as.factor(full_life$rank.change))



#######Look at fitness change within individuals####
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


changers.list <- list()
counter <- 1
for(id in unique(filter(ranks, RankChange != 'None')$ID)){
  change_year <- filter(ranks, ID == id, RankChange != 'None')$Year %>% max()
  if((change_year - 2) %in% filter(ranks, ID == id)$Year & 
     (change_year + 2) %in% filter(ranks, ID == id)$Year){
    changers.list[[counter]] <- data.frame(id,
                                    change.year = change_year,
                                    pre.ars = mean(filter(ranks, ID == id,
                                                   Year < change_year)$ARS),
                                    post.ars = mean(filter(ranks, ID == id,
                                                           Year > change_year)$ARS),
                                    net.change = sum(filter(ranks, ID == id)$RankDiff),
                                    net.direction = ifelse(sum(filter(ranks, ID == id)$RankDiff) < 0,
                                                           'Down',
                                                           'Up'),
                                    rank.group = ifelse(filter(ranks, ID == id,
                                                               Year == change_year)$stan.rank >= 0,
                                                        'High',
                                                        'Low'),
                                    pre.surv2 = mean(filter(ranks, ID == id,
                                                           Year > change_year)$Surv2),
                                    post.surv2 = mean(filter(ranks, ID == id,
                                                           Year < change_year)$Surv2))
    counter <- counter+1
  }
}
changers <- do.call(rbind, changers.list)
changers$ars.diff <- changers$post.ars - changers$pre.ars
#changers$rank.group <- as.factor(changers$rank.group)

ggplot(data = changers,
       aes( y = ars.diff, x= net.direction)) +
  geom_label(aes(label = id), position = position_jitter(width = 0.3))+
  geom_violin(alpha = 0.2)+
  facet_wrap( ~ rank.group)+
  theme_classic()+
  xlab('Direction of rank change') + 
  ylab('Consequence of change in yearly offpsring produced')

ggplot(data = changers,
       aes( y = post.surv2 - pre.surv2, x= net.direction)) +
  geom_label(aes(label = id), position = position_jitter(width = 0.3))+
  geom_violin(alpha = 0.2)+
  facet_grid(. ~ rank.group)+
  theme_classic()+
  xlab('Direction of rank change') + 
  ylab('Consequence of change in yearly offpsring that survive to 2yo')




ranks %>% group_by(ID) %>% summarize(mean.rank = mean(stan.rank), 
                                     total.cubs = sum(Surv2),
                                     mean.cubs = mean(Surv2),
                                     mean.ars = mean(ARS)) -> lrs

lrs <- filter(lrs, ID %in% filter(tblHyenas, !is.na(Disappeared))$ID)

survive_to_4 <- filter(tblHyenas, !is.na(Disappeared), Disappeared - Birthdate >= 4*365)

ggplot(filter(lrs, ID %in% survive_to_4$ID), aes(x = mean.rank, y = total.cubs, label = ID))+
  geom_label() +
  geom_smooth(method = 'lm')+
  theme_classic()


pois.quad.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ poly(mean.rank, 2), family = 'poisson')
summary(pois.quad.lrs)

pois.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ mean.rank, family = 'poisson')
summary(pois.lrs)
AIC(pois.quad.lrs, pois.lrs)


quad.ars <- lm(data = filter(lrs, ID %in% survive_to_4$ID), formula = mean.cubs ~ poly(mean.rank, 2) )
summary(quad.ars)

ars <- lm(data = filter(lrs, ID %in% survive_to_4$ID), formula = mean.cubs ~ mean.rank)
summary(ars)

AIC(quad.ars, ars)


ggplot(filter(lrs, ID %in% survive_to_4$ID), aes(x = mean.rank, y = mean.cubs, label = ID))+
  geom_label() +
  geom_smooth(method = 'loess')+
  theme_classic()



##Feeding##
scans <- read_csv('~/Documents/Fisibase/tblScans.csv', col_types = 'ccccccc')
scans$Session <- as.character(scans$Session)
scans$Date <- left_join(scans, tblSessions, by = c('Session' = 'session'))$date
format(scans$Date, '%Y') %>% as.numeric() %>% hist()

changers$pre.fd <- NA
changers$post.fd <- NA
changers$pre.fd.obs <- NA
changers$post.fd.obs <- NA
for(row in 1:nrow(changers)){
  change_year <- changers[row,]$change.year
  pre.id.scans <- filter(scans, format(Date, '%Y') < change_year, Hyena == changers[row,]$id)
  post.id.scans <- filter(scans, format(Date, '%Y') > change_year, Hyena == changers[row,]$id)
  changers$pre.fd.obs[row] <- nrow(pre.id.scans)
  changers$post.fd.obs[row] <- nrow(post.id.scans)
  changers$pre.fd[row] <- nrow(filter(pre.id.scans, BehaviorCode %in% c('fd', 'fd?'))) / 
    nrow(pre.id.scans)
  changers$post.fd[row] <- nrow(filter(post.id.scans, BehaviorCode %in% c('fd', 'fd?'))) / 
    nrow(post.id.scans)
}

changers$fd.dif <- changers$post.fd - changers$pre.fd

ggplot(data = filter(changers, pre.fd.obs >= 10, post.fd.obs >= 10),
       aes( y = fd.dif, x= net.direction)) +
  geom_label(aes(label = id), position = position_jitter(width = 0.3))+
  geom_violin(alpha = 0.2)+
  theme_classic()+
  xlab('Direction of rank change') + 
  ylab('Difference in feeding before vs after rank change')



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

# matriline_labels <- ribbon %>% group_by(matriline) %>% 
#   summarize (x.start = 1997,
#              y.start = mean(c(ymin[which(Year == 1997)], ymax[which(Year == 1997)])))

matriline_labels <- data.frame(x.start = rep(1987), 
                               matriline = c('kb', 'dj', 'coch', '03', '40', 'other'),
                               y.start = c(1.5, 3, 5, 6.5, 8, 15))
matriline_labels <- arrange(matriline_labels, matriline)

matriline_labels <- data.frame(x.start = rep(1987), 
                               matriline = c('kb', 'dj', '03', '40', 'other'),
                               y.start = c(0.7, 3.3, 5.5, 8, 15))
matriline_labels <- arrange(matriline_labels, matriline)


desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(5)
cols[2] <- 'darkgoldenrod1'

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
             color = 'gray25')+
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = desat(cols, sat= 0.5))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())



stan.rank.labels <- data.frame(matriline = sort(unique(filter(ranks, Clan == 'talek')$matriline)),
                               start = rep(1988),
                               start.stan = filter(ribbon, Year == 1988)$max.stan,
                               mean.start.stan = filter(ribbon, Year == 1988)$mean.stan)

#mean standardized rank
mean.matriline.rank <- ggplot(ribbon, aes(x = Year, y= mean.stan, group = matriline))+
  geom_line(data = ribbon, aes(col = matriline), size = 2)+
  theme_classic()+
  geom_label(data = stan.rank.labels,
             aes(x = start, y = mean.start.stan, label = matriline, group = matriline, fill = matriline),
             label.size = NA,
             label.padding = unit(0.2, 'lines'),
             color = 'gray25')+
  scale_color_manual(values = cols)  + 
  scale_fill_manual(values = desat(cols, sat= 0.5))+
  theme(legend.position = 'none')+
  ylab('Mean standardized rank')

grid.newpage()
grid.draw(rbind(ggplotGrob(matriline.colors), 
                ggplotGrob(mean.matriline.rank),
          size = 'last'))
gridExtra::grid.arrange(matriline.colors, mean.matriline.rank, nrow = 2, heights = c(2,1))

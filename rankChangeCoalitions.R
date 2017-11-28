##################################################################
##                        Eli Strauss                           ##
##      Predict magnitude of rank change using coalitions       ##
##                     April 11th, 2017                         ##
##################################################################
#library(igraph)
#library(mlogit)
#library(ggplot2)
#source("/mnt/home/straus46/RankChangeModels/Fisibase/fisibasetidy/ReadTidyData.R")
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")
library(magrittr)
library(tidyverse)

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
                    ally_ranks <- rep(NA)
)

for(row in 1:tibdim(ranks)){
  year <- ranks[row,]$Year
  clan <- ranks[row,]$Clan
  focal <- ranks[row,]$ID
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
      ranks[row,]$ally_ranks <- sum(filter(ranks, ID %in% allies, Year == year)$Rank)
      ranks[row,]$coal_top3_wtd <- coal_net_wtd[focal,order(coal_net_wtd[focal,], decreasing = T)[1:3]] %>% 
        sum() 
    } 
  }
}



##add metrics to ranks tibble
ranks <- add_column(ranks, ai_deg = rep(NA, tibdim(ranks)),
                    ai_top3_deg = rep(NA, tibdim(ranks))
)

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



################Do permutations all ranks together - unweighted coalitions###############

ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

coal_perm_stacked <- tibble()
#par(mfrow = c(2,2))
for(category in list('High', 'Low')){
  ranks_mod <- filter(ranks, RankCategory %in% category)
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                                 tibble(coal_all_deg = ranks_mod$coal_deg,
                                        coal_top_deg = ranks_mod$coal_top3_deg,
                                        rank_change = ranks_mod$RankDiffAbs,
                                        rank_category = category,
                                        network = 'Observed'))
  
  iterations <- 100
  
  rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

  coal_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_deg + (1|ID))
  coal_top3_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + (1|ID))
  
  obs_coef <- as_tibble(cbind(coal_deg = coal_mod %>% .@beta %>% .[2],
                              coal_top3_deg = coal_top3_mod %>% .@beta %>% .[2]))
  perm_coef <- tibble(coal_deg = rep(NA, iterations), 
                      coal_top3_deg = rep(NA, iterations))
  
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
      lme4::lmer(formula = RankDiffAbs ~ coal_deg + (1|ID)) %>% 
      .@beta %>% .[2] -> perm_coef[i,1]
    
    ranks.perm %>%
      lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + + (1|ID)) %>%
      .@beta %>% .[2] -> perm_coef[i,2]
    
    rank_perm_change[,i] <- as.character(ranks.perm$RankDiffAbs)
  }
  
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = rep(ranks_mod$coal_deg, iterations),
                                      coal_top_deg = rep(ranks_mod$coal_top3_deg, iterations),
                                      rank_change = as.numeric(as.vector(rank_perm_change)),
                                      rank_category = category,
                                      network = 'Permuted'))
  
  # coal_plot <- hist(perm_coef$coal_deg, breaks = 50, col = 'black',
  #                   main = paste0(category, '--All'),
  #                   xlab = 'Coefficient estimate')
  # segments(x1 = obs_coef$coal_deg, x0 = obs_coef$coal_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
  # text(x = -.018, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$coal_deg >= obs_coef$coal_deg))/iterations))
  # 
  # coal_top3_plot <- hist(perm_coef$coal_top3_deg, breaks = 50, col = 'black',
  #                        main = paste0(category, '--Top'),
  #                        xlab = 'Total Degree of Coalitionary Aggression')
  # segments(x1 = obs_coef$coal_top3_deg, x0 = obs_coef$coal_top3_deg, y0 = 0, y1 = max(coal_top3_plot$counts), col = 'red', lty=2, lwd=2)
  # text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/iterations))
  # 
  range <- c(min(perm_coef$coal_top3_deg)-.005, max(perm_coef$coal_top3_deg, obs_coef$coal_top3_deg) + .005)
  mar.default <- c(5,4,4,2) + 0.1
  par(family = 'Trebuchet MS', mar = mar.default + c(0, 1, 0, 0))
  d <- density(perm_coef$coal_top3_deg)
  coal_top3_plot <- plot(d, lwd = 6, col = 'grey87',
                         main = paste0('Effect size = ', round(obs_coef$coal_top3_deg, digits = 3), ', ', 'p = ', round(2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/iterations, digits = 3)),
                         xlab = 'Degree of coalitionary aggression',
                         cex.lab=2, cex.axis=2, cex.main=2,
                         xlim = range)
  polygon(d, col = 'grey87', border = 'grey87')
  abline(v = obs_coef$coal_top3_deg, col = 'dodgerblue1', lty=2, lwd=6)
}


facet_labels = c(High = 'High ranking individuals', Low = 'Low ranking individuals')
colors = c('dodgerblue1', 'grey87')

coal_perm_stacked$rank_change_category <- ifelse(coal_perm_stacked$rank_change == 0,'None', 
                                                 ifelse(coal_perm_stacked$rank_change < 0,'Down','Up'))
coal_perm_stacked$rank_change_category %<>% factor(., levels = c('Down', 'None', 'Up'))

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



all_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change, y = coal_all_deg, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with all allies')
all_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold'))





################Analyze whole hierarchy together - unweighted coalitions###############

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

coal_perm_stacked <- tibble()
#par(mfrow = c(2,2))
coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = ranks_mod$coal_deg,
                                      coal_top_deg = ranks_mod$coal_top3_deg,
                                      rank_change = ranks_mod$RankDiffAbs,
                                      network = 'Observed'))

iterations <- 100

rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

coal_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_deg + (1|ID))
coal_top3_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + (1|ID))
  
obs_coef <- as_tibble(cbind(coal_deg_intercept = coal_mod %>% .@beta %>% .[1],
                            coal_deg = coal_mod %>% .@beta %>% .[2],
                            coal_top3_intercept = coal_top3_mod %>% .@beta %>% .[1],
                            coal_top3_deg = coal_top3_mod %>% .@beta %>% .[2]))
perm_coef <- tibble(coal_deg_intercept = rep(NA, iterations),
                    coal_deg = rep(NA, iterations), 
                    coal_top3_intercept = rep(NA, iterations),
                    coal_top3_deg = rep(NA, iterations))

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
    lme4::lmer(formula = RankDiffAbs ~ coal_deg + (1|ID)) -> all_perm_mod
  all_perm_mod%>% 
    .@beta %>% .[1] -> perm_coef[i,1]
  all_perm_mod%>% 
    .@beta %>% .[2] -> perm_coef[i,2]
  
  ranks.perm %>%
    lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + (1|ID)) -> top_perm_mod
  top_perm_mod %>%
    .@beta %>% .[1] -> perm_coef[i,3]
  top_perm_mod %>%
    .@beta %>% .[2] -> perm_coef[i,4]
  
  rank_perm_change[,i] <- as.character(ranks.perm$RankDiffAbs)
}
  
coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = rep(ranks_mod$coal_deg, iterations),
                                      coal_top_deg = rep(ranks_mod$coal_top3_deg, iterations),
                                      rank_change = as.numeric(as.vector(rank_perm_change)),
                                      network = 'Permuted'))

# coal_plot <- hist(perm_coef$coal_deg, breaks = 50, col = 'black',
#                   main = paste0(category, '--All'),
#                   xlab = 'Coefficient estimate')
# segments(x1 = obs_coef$coal_deg, x0 = obs_coef$coal_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
# text(x = -.018, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$coal_deg >= obs_coef$coal_deg))/iterations))
# 
# coal_top3_plot <- hist(perm_coef$coal_top3_deg, breaks = 50, col = 'black',
#                        main = paste0(category, '--Top'),
#                        xlab = 'Total Degree of Coalitionary Aggression')
# segments(x1 = obs_coef$coal_top3_deg, x0 = obs_coef$coal_top3_deg, y0 = 0, y1 = max(coal_top3_plot$counts), col = 'red', lty=2, lwd=2)
# text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/iterations))
# 
range <- c(min(perm_coef$coal_top3_deg)-.005, max(perm_coef$coal_top3_deg, obs_coef$coal_top3_deg) + .005)
mar.default <- c(5,4,4,2) + 0.1
par(family = 'Trebuchet MS', mar = mar.default + c(0, 1, 0, 0))
d <- density(perm_coef$coal_top3_deg)
coal_top3_plot <- plot(d, lwd = 6, col = 'grey87',
                       main = paste0('Effect size = ', round(obs_coef$coal_top3_deg, digits = 3), ', ', 'p = ', round(2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/iterations, digits = 3)),
                       xlab = 'Degree of coalitionary aggression',
                       cex.lab=2, cex.axis=2, cex.main=2,
                       xlim = range)
polygon(d, col = 'grey87', border = 'grey87')
abline(v = obs_coef$coal_top3_deg, col = 'dodgerblue1', lty=2, lwd=6)


colors = c('dodgerblue1', 'grey87')

coal_perm_stacked$rank_change_category <- ifelse(coal_perm_stacked$rank_change == 0,'None', 
                                                 ifelse(coal_perm_stacked$rank_change < 0,'Down','Up'))
coal_perm_stacked$rank_change_category %<>% factor(., levels = c('Down', 'None', 'Up'))

ggplot(data = filter(coal_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = coal_top_deg, y = rank_change, col = network)) + 
  geom_jitter(col = alpha('grey', .2))+
  geom_jitter(data = filter(coal_perm_stacked, network == 'Observed'), aes(x = coal_top_deg, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$coal_top3_intercept, slope = perm_coef$coal_top3_deg, col = 'grey', alpha = 0.05, size = 1) +
  geom_abline(intercept = obs_coef$coal_top3_intercept, slope = obs_coef$coal_top3_deg, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()



#######old plots######

top_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change_category, y = coal_top_deg, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with top allies')+
  coord_cartesian(ylim= c(0,35))+
  coord_flip()

top_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold', size = 12,family = 'Trebuchet MS'),
        axis.title = element_text(size = 14,family = 'Trebuchet MS'),
        axis.text = element_text(size = 12, family = 'Trebuchet MS'),
        legend.text = element_text(size = 12, family = 'Trebuchet MS'),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank())



all_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change_category, y = coal_all_deg, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with all allies')
all_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold'))



################Analyze coalitions and AIs together for whole hierarchy ###############

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

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

iterations <- 100

rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

both_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_deg + ai_deg + (1|ID))
both_top3_mod <- ranks_mod %>% lme4::lmer(formula = RankDiffAbs ~ coal_top3_deg + ai_top3_deg + (1|ID))

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

# coal_plot <- hist(perm_coef$coal_deg, breaks = 50, col = 'black',
#                   main = paste0(category, '--All'),
#                   xlab = 'Coefficient estimate')
# segments(x1 = obs_coef$coal_deg, x0 = obs_coef$coal_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
# text(x = -.018, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$coal_deg >= obs_coef$coal_deg))/iterations))
# 
# coal_top3_plot <- hist(perm_coef$coal_top3_deg, breaks = 50, col = 'black',
#                        main = paste0(category, '--Top'),
#                        xlab = 'Total Degree of Coalitionary Aggression')
# segments(x1 = obs_coef$coal_top3_deg, x0 = obs_coef$coal_top3_deg, y0 = 0, y1 = max(coal_top3_plot$counts), col = 'red', lty=2, lwd=2)
# text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/iterations))
# 

par(mfrow = c(ceiling(length(perm_coef)/2), 2))
for(param in names(perm_coef)){
  range <- c(min(perm_coef[param])-.005, max(perm_coef[param], obs_coef[param]) + .005)
  mar.default <- c(5,4,4,2) + 0.1
  par(family = 'Trebuchet MS', mar = mar.default + c(0, 1, 0, 0))
  d <- density(perm_coef[[param]])
  coal_top3_plot <- plot(d, lwd = 6, col = 'grey87',
                         main = paste0('Effect size = ', round(obs_coef[param], digits = 3), ', ', 'p = ', round(2*tibdim(which(perm_coef[[param]] >= obs_coef[[param]]))/iterations, digits = 3)),
                         xlab = param,
                         cex.lab=2, cex.axis=2, cex.main=2,
                         xlim = range)
  polygon(d, col = 'grey87', border = 'grey87')
  abline(v = obs_coef[param], col = 'dodgerblue1', lty=2, lwd=6)
}



colors = c('dodgerblue1', 'grey87')

coal_perm_stacked$rank_change_category <- ifelse(coal_perm_stacked$rank_change == 0,'None', 
                                                 ifelse(coal_perm_stacked$rank_change < 0,'Down','Up'))
coal_perm_stacked$rank_change_category %<>% factor(., levels = c('Down', 'None', 'Up'))


###Coalitions###
ggplot(data = filter(coal_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = coal_top_deg, y = rank_change, col = network)) + 
  geom_jitter(col = alpha('grey', .2))+
  geom_jitter(data = filter(coal_perm_stacked, network == 'Observed'), aes(x = coal_top_deg, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$both_top3_intercept, slope = perm_coef$both_coal_top3, col = 'grey', alpha = 0.05, size = 1) +
  geom_abline(intercept = obs_coef$both_top3_intercept, slope = obs_coef$both_coal_top3, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()

###AI#####
ggplot(data = filter(coal_perm_stacked, network == 'Permuted')[1:(100*length(rank_perm_change)/iterations),], aes(x = ai_top_deg, y = rank_change, col = network)) + 
  geom_jitter(col = alpha('grey', .2))+
  geom_jitter(data = filter(coal_perm_stacked, network == 'Observed'), aes(x = ai_top_deg, y = rank_change), col = 'dodgerblue')+
  geom_abline(intercept = perm_coef$both_top3_intercept, slope = perm_coef$both_ai_top3, col = 'grey', alpha = 0.05, size = 1) +
  geom_abline(intercept = obs_coef$both_top3_intercept, slope = obs_coef$both_ai_top3, col = 'dodgerblue', size = 1.5) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()








################Coalitions weighted by AI###############

ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

coal_perm_stacked <- tibble()
par(mfrow = c(2,2))
for(category in list('High', 'Low')){
  ranks_mod <- filter(ranks, RankCategory %in% category)
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                                 tibble(coal_all_wtd = ranks_mod$coal_deg_wtd,
                                        coal_top_wtd = ranks_mod$coal_top3_wtd,
                                        rank_change = ranks_mod$RankDiffAbs,
                                        rank_category = category,
                                        network = 'Observed'))
  
  iterations <- 1000
  
  rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)
  
  coal_mod <- ranks_mod %>% glm(formula = RankDiffAbs ~ coal_deg_wtd, family = gaussian)
  coal_top3_mod <- ranks_mod %>% glm(formula = RankDiffAbs ~ coal_top3_wtd, family = gaussian)
  
  obs_coef <- as_tibble(cbind(coal_deg_wtd = coal_mod %>% coef() %>% .[2],
                              coal_top3_wtd = coal_top3_mod %>% coef() %>% .[2]))
  perm_coef <- tibble(coal_deg_wtd = rep(NA, iterations), 
                      coal_top3_wtd = rep(NA, iterations))
  
  ranks.perm <- ranks_mod
  ###Select variables to permute within##
  ranks.perm$PermCategory <- paste(ranks_mod$Clan, ranks_mod$Year)
  
  for(i in 1:iterations){
    ranks.perm %>% 
      group_by(PermCategory) %>% 
      sample_frac(replace = F) %>% 
      ungroup() %>% 
      select(RankDiffAbs) %>% .[[1]] ->
      ranks.perm$RankDiffAbs
    
    ranks.perm %>%
      glm(formula = RankDiffAbs ~ coal_deg_wtd, family = gaussian) %>%
      coef() %>% .[2] -> perm_coef[i,1]
    
    ranks.perm %>%
      glm(formula = RankDiffAbs ~ coal_top3_wtd, family = gaussian) %>%
      coef() %>% .[2] -> perm_coef[i,2]
    
    rank_perm_change[,i] <- as.character(ranks.perm$RankDiffAbs)
  }
  
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                                 tibble(coal_all_wtd = rep(ranks_mod$coal_deg_wtd, iterations),
                                        coal_top_wtd = rep(ranks_mod$coal_top3_wtd, iterations),
                                        rank_change = as.numeric(as.vector(rank_perm_change)),
                                        rank_category = category,
                                        network = 'Permuted'))
  
  coal_plot <- hist(perm_coef$coal_deg_wtd, breaks = 50, col = 'black',
                    main = paste0(category, '--All'),
                    xlab = 'Coefficient estimate')
  segments(x1 = obs_coef$coal_deg_wtd, x0 = obs_coef$coal_deg_wtd, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
  text(x = -.018, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$coal_deg_wtd >= obs_coef$coal_deg_wtd))/iterations))
  
  coal_top3_plot <- hist(perm_coef$coal_top3_wtd, breaks = 50, col = 'black',
                         main = paste0(category, '--Top'),
                         xlab = 'Coefficient estimate')
  segments(x1 = obs_coef$coal_top3_wtd, x0 = obs_coef$coal_top3_wtd, y0 = 0, y1 = max(coal_top3_plot$counts), col = 'red', lty=2, lwd=2)
  text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$coal_top3_wtd >= obs_coef$coal_top3_wtd))/iterations))
}


facet_labels = c(High = 'High ranking individuals', Low = 'Low ranking individuals')
colors = c('dodgerblue1', 'grey87')

coal_perm_stacked$rank_change_category <- ifelse(coal_perm_stacked$rank_change == 0,'None', 
                                                 ifelse(coal_perm_stacked$rank_change < 0,'Down','Up'))
coal_perm_stacked$rank_change_category %<>% factor(., levels = c('Up', 'None', 'Down'))

top_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change_category, y = coal_top_wtd, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with top allies')#+
  #coord_cartesian(ylim= c(0,35))

top_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold', size = 12,family = 'Trebuchet MS'),
        axis.title = element_text(size = 14,family = 'Trebuchet MS'),
        axis.text = element_text(size = 12, family = 'Trebuchet MS'),
        legend.text = element_text(size = 12, family = 'Trebuchet MS'))



all_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change_category, y = coal_all_wtd, fill = network)) + 
  geom_boxplot(alpha = .7, outlier.shape = NA)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with all allies')
all_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold'))



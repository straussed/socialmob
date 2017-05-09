##################################################################
##                        Eli Strauss                           ##
##      Modeling magnitude of rank change using permutation     ##
##                     April 11th, 2017                         ##
##################################################################
#library(igraph)
#library(mlogit)
#library(ggplot2)
#source("/mnt/home/straus46/RankChangeModels/Fisibase/fisibasetidy/ReadTidyData.R")
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")


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

################Do permutations all ranks together###############

##Determine start and end years for each clan
clans <- unique(ranks$Clan)
clan_years <- tibble(Clan = clans, First = NA, Last = NA)
for(row in 1:tibdim(clan_years$Clan)){
  clan_years[row,'First'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% min(na.rm = T)
  clan_years[row,'Last'] <- filter(ranks, Clan == clan_years[row,][[1]])$Year %>% max(na.rm = T)
}

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

iterations <- 10000
ai_mod <- ranks_mod %>% lm(formula = RankDiff ~ ai_deg)
ai_top3_mod <- ranks_mod %>% lm(formula = RankDiff ~ ai_top3_deg)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                   ai_top3_deg = ai_top3_mod %>% coef() %>% .[2]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    ai_top3_deg = rep(NA, iterations))

ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks_mod$Clan, ranks_mod$Year)

for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiff) %>% .[[1]] ->
    ranks.perm$RankDiff

  ranks.perm %>%
    lm(formula = RankDiff ~ ai_deg) %>%
    coef() %>% .[2] -> perm_coef[i,1]
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_top3_deg) %>%
    coef() %>% .[2] -> perm_coef[i,2]
}

par(mfrow = c(1,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black',
                main = 'Effects from model with all allies',
                xlab = 'Effect of AI with all allies')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
text(x = -.018, y = 500, paste0('p = ', tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/10000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',
                     main = 'Effects from model with all allies',
                     xlab = 'Effect of AI with top allies')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red', lty=2, lwd=2)
text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/10000))

###############################################


################Test effect of AI for high rankers only###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

#ranks_mod <- ranks_mod %>% group_by(Clan, Year) %>% mutate(ai_deg = (ai_deg-mean(ai_deg))/sd(ai_deg), 
                                               #            ai_top3_deg = (ai_top3_deg-mean(ai_top3_deg))/sd(ai_top3_deg))

                    
iterations <- 10000
ranks_mod <- filter(ranks_mod, RankCategory == 'High')
mod_length <- tibdim(ranks_mod)
###matrix store permuted values
rank_diff_perm <- matrix(NA, nrow=mod_length, ncol=iterations)
ranks_mod$RankDev <- abs(ranks_mod$RankDiffAbs)

ai_mod <-  ranks_mod %>% lm(formula = RankDiff ~ ai_deg)
ai_top3_mod <- ranks_mod %>% lm(formula = RankDiff ~ ai_top3_deg)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    ai_top3_deg = rep(NA, iterations))

ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks.perm$Clan, ranks.perm$Year)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiffAbs) %>% .[[1]] ->
    ranks.perm$RankDiffAbs
  
  ranks.perm %>%
    lm(formula = RankDiffAbs ~ ai_deg) %>%
    coef() %>% .[2] -> perm_coef[i,1]
  
  ranks.perm %>%
    lm(formula = RankDiffAbs ~ ai_top3_deg) %>%
    coef() %>% .[2] -> perm_coef[i,2]
  
  rank_diff_perm[,i] <- ranks.perm$RankDiffAbs
}

par(mfrow = c(1,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black',
                main = 'Effects from model with all allies\n>High status<',
                xlab = 'Effect of AI with all allies')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red',lty=2,lwd=2)
text(x = -.015, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/iterations))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',
                     xlab = 'Effect of AI with top allies',
                     main = '')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red',lty=2,lwd=2)
text(x = .5, y = 500, paste0('p = ', 2*tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/iterations))

ai_perm_stacked <- tibble(ai_all_deg = rep(ranks_mod$ai_deg, iterations),
                          ai_top_deg = rep(ranks_mod$ai_top3_deg, iterations),
                          rank_diff = as.vector(rank_diff_perm),
                          rank_change = ifelse(rank_diff == 0, 'None',
                                               ifelse(rank_diff < 0, 'Down',
                                                      'Up')),
                          rank_category = 'High',
                          network = 'Permuted')

ai_perm_stacked <- bind_rows(ai_perm_stacked,
                             tibble(ai_all_deg = ranks_mod$ai_deg,
                             ai_top_deg = ranks_mod$ai_top3_deg,
                             rank_diff = ranks_mod$RankDiff,
                             rank_change = ifelse(rank_diff == 0, 'None',
                                                  ifelse(rank_diff < 0, 'Down',
                                                         'Up')),
                             rank_category = 'High',
                             network = 'Observed'))



################Test effect of AI for low rankers only###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

#ranks_mod <- ranks_mod %>% group_by(Clan, Year) %>% mutate(ai_deg = (ai_deg-mean(ai_deg))/sd(ai_deg), 
                                                           #ai_top3_deg = (ai_top3_deg-mean(ai_top3_deg))/sd(ai_top3_deg))
ranks_mod <- filter(ranks_mod, RankCategory == 'Low')
mod_length <- tibdim(ranks_mod)
###matrix store permuted values
rank_diff_perm <- matrix(NA, nrow=mod_length, ncol=iterations)
iterations <- 10000
ai_mod <- ranks_mod %>% lm(formula = RankDiff ~ ai_deg)
ai_top3_mod <- ranks_mod %>% lm(formula = RankDiff ~ ai_top3_deg)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    ai_top3_deg = rep(NA, iterations))

ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks.perm$Clan, ranks.perm$Year)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiffAbs) %>% .[[1]] ->
    ranks.perm$RankDiffAbs
  
  ranks.perm %>%
    lm(formula = RankDiffAbs ~ ai_deg) %>%
    coef() %>% .[2] -> perm_coef[i,1]
  
  ranks.perm %>%
    lm(formula = RankDiffAbs ~ ai_top3_deg) %>%
    coef() %>% .[2] -> perm_coef[i,2]
  
  rank_diff_perm[,i] <- ranks.perm$RankDiffAbs
}

par(mfrow = c(1,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black',
                main = 'Effects from model with all allies\n>Low status<',
                xlab = 'Effect of AI with all allies')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty = 2, lwd = 2)
text(x = -.022, y = 600, paste0('p = ', 2*tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/iterations))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',
                     main = 'Effects from model with top allies\n>Low status<',
                     xlab = 'Effect of AI with top allies')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red', lty = 2, lwd = 2)
text(x = -.09, y = 450, paste0('p = ', 2*tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/iterations))


ai_perm_stacked <- bind_rows(ai_perm_stacked,
                             tibble(ai_all_deg = rep(ranks_mod$ai_deg, iterations),
                                    ai_top_deg = rep(ranks_mod$ai_top3_deg, iterations),
                                    rank_diff = as.vector(rank_diff_perm),
                                    rank_change = ifelse(rank_diff == 0, 'None',
                                                         ifelse(rank_diff < 0, 'Down',
                                                                'Up')),
                                    rank_category = 'Low',
                                    network = 'Permuted'))

ai_perm_stacked <- bind_rows(ai_perm_stacked,
                             tibble(ai_all_deg = ranks_mod$ai_deg,
                                    ai_top_deg = ranks_mod$ai_top3_deg,
                                    rank_diff = ranks_mod$RankDiffAbs,
                                    rank_change = ifelse(rank_diff == 0, 'None',
                                                         ifelse(rank_diff < 0, 'Down',
                                                                'Up')),
                                    rank_category = 'Low',
                                    network = 'Observed'))

ai_perm_stacked$rank_change <- factor(ai_perm_stacked$rank_change,
                                      levels = c('Up', 'None', 'Down'))


facet_labels = c(High = 'High ranking individuals', Low = 'Low ranking individuals')
colors = c('dodgerblue1', 'grey87')
  
top_plot <- ggplot(data = ai_perm_stacked, aes(x = rank_change, y = ai_top_deg, fill = network)) + 
  geom_boxplot(alpha = .7)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of association with top allies')
top_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold', size = 12,family = 'Trebuchet MS'),
        axis.title = element_text(size = 14,family = 'Trebuchet MS'),
        axis.text = element_text(size = 12, family = 'Trebuchet MS'),
        legend.text = element_text(size = 12, family = 'Trebuchet MS'))
  


all_plot <- ggplot(data = ai_perm_stacked, aes(x = rank_change, y = ai_all_deg, fill = network)) + 
  geom_boxplot(alpha = .7)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of association with all allies')
all_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold'))



################Test effect of AI on direction of rank movement among movers of all ranks###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

#ranks_mod <- ranks_mod %>% group_by(Clan, Year) %>% mutate(ai_deg = (ai_deg-mean(ai_deg))/sd(ai_deg), 
#            ai_top3_deg = (ai_top3_deg-mean(ai_top3_deg))/sd(ai_top3_deg))


iterations <- 1000
ranks_mod <- filter(ranks_mod, RankChange != 'None', RankCategory == 'High')
mod_length <- tibdim(ranks_mod)
###matrix store permuted values
rank_change_perm <- matrix(NA, nrow=mod_length, ncol=iterations)
ranks_mod$RankDev <- abs(ranks_mod$RankDiffAbs)

ai_mod <-  ranks_mod %>% glm(formula = RankChange ~ ai_deg, family = binomial)
ai_top3_mod <- ranks_mod %>% glm(formula = RankChange ~ ai_top3_deg, family = binomial)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    ai_top3_deg = rep(NA, iterations))

ranks.perm <- ranks_mod
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks.perm$Clan, ranks.perm$Year)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankChange) %>% .[[1]] ->
    ranks.perm$RankChange
  
  ranks.perm %>%
    glm(formula = RankChange ~ ai_deg, family = binomial) %>%
    coef() %>% .[2] -> perm_coef[i,1]
  
  ranks.perm %>%
    glm(formula = RankChange ~ ai_top3_deg, family =binomial) %>%
    coef() %>% .[2] -> perm_coef[i,2]
  
  rank_change_perm[,i] <- ranks.perm$RankChange
}

par(mfrow = c(1,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black',
                main = 'Effects from model with all allies\n>Low status<',
                xlab = 'Effect of AI with all allies')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty = 2, lwd = 2)
text(x = -.022, y = 600, paste0('p = ', 2*tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/10000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',
                     main = 'Effects from model with top allies\n>Low status<',
                     xlab = 'Effect of AI with top allies')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red', lty = 2, lwd = 2)
text(x = -.09, y = 450, paste0('p = ', 2*tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/1000))



####################AI by dyad - are dyads that associate less more likely to change ranks?########


###build dyad tibble with AI for each dyad in each year
get_ai_previous<- function(df){
  Year = unique(df$Year)-1
  Clan = unique(df$Clan)
  if(exists(paste('ai_net', Year, Clan, sep = '_'))){
    ai_net <- get(paste('ai_net', Year, Clan, sep = '_'))
    df[df$Focal %in% dimnames(ai_net)[[1]] & df$Alter %in% dimnames(ai_net)[[1]],]$AI_previous <- 
      df[,c('Focal', 'Alter')] %>% 
      filter(Focal %in% dimnames(ai_net)[[1]] & Alter %in% dimnames(ai_net)[[1]]) %>%
      as.matrix() %>% 
      ai_net[.]
  }
  return(df)
}

ranks_mod <- anti_join(ranks, clan_years, by = c('Clan', 'Year' = 'First')) %>%
  anti_join(clan_years, by = c('Clan', 'Year' = 'Last'))

dyads <- ranks_mod %>%
  group_by(Clan, Year) %>% 
  do(expand.grid(Focal = .$ID, Alter = .$ID)) %>%
  ungroup()

dyads$AI_previous <- NA
for(clan in unique(dyads$Clan)){
  for(year in unique(dyads[dyads$Clan == clan,]$Year)){
    dyads[dyads$Clan == clan & dyads$Year == year,] <- filter(dyads, Clan == clan, Year == year) %>% get_ai_previous
  }
}

dyads <- dyads[-which(dyads$Focal == dyads$Alter),]


####Add in rank and rank change info###############
dyads$RankChange <- factor('NoChange', levels = c('NoChange', 'Change'))
rank_changers <- filter(ranks, !RankChange %in% c('None', 'Both'))
for(row in 1:tibdim(rank_changers)){
  id <- rank_changers[row,]$ID
  year <- rank_changers[row,]$Year
  clan <- rank_changers[row,]$Clan
  rank_temp <- filter(ranks, Clan == clan, Year == year)
  pos <- which(rank_temp$ID == id)
  old_pos <- which(rank_temp$IDold == id)
  swappers <- c(rank_temp$ID[1:pos][!rank_temp$ID[1:pos] %in% rank_temp$IDold[1:old_pos]],
                rank_temp$IDold[1:old_pos][!rank_temp$IDold[1:old_pos] %in% rank_temp$ID[1:pos]])
  dyads[dyads$Focal == id & 
          dyads$Alter %in% swappers & 
          dyads$Year == year,]$RankChange <- 'Change'
  dyads[dyads$Focal %in% swappers & 
          dyads$Alter == id & 
          dyads$Year == year,]$RankChange <- 'Change'
}

dyads$Focal_rank <- left_join(dyads, ranks_mod, by = c('Year', 'Focal' = 'ID'))$RankCategory
dyads$Alter_rank <- left_join(dyads, ranks_mod, by = c('Year', 'Alter' = 'ID'))$RankCategory
#############Model effect of AI on presence of rank change############

###All ranks
dyads_mod <- dyads[complete.cases(dyads[,c('RankChange', 'AI_previous')]),]
obs_coef <- coef(glm(formula = RankChange ~ AI_previous, data = dyads_mod, family = 'binomial'))[2]

iterations = 10000
perm_coef  <- tibble(AI_previous = rep(NA, iterations))
dyads_perm <- dyads_mod
for(i in 1:iterations){
  dyads_perm %>% 
    group_by(Clan, Year) %>%
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankChange) %>%
    .[[1]] ->
    dyads_perm$RankChange
  
  dyads_perm %>% 
    glm(formula = RankChange ~ AI_previous, family = 'binomial') %>%
    coef() %>% .[2] -> perm_coef[i,1]
}

par(mfrow = c(1,1))
dyad_plot <- hist(perm_coef$AI_previous, breaks = 50, col = 'black', main = 'Changes among all individuals', 
                  xlab = 'Effect of AI from previous year')
segments(x1 = obs_coef, x0 = obs_coef, y0 = 0, y1 = max(dyad_plot$counts), col = 'red', lwd = 2, lty = 2)
text(x = -3, y = 300, paste0('p = ', tibdim(which(perm_coef$AI_previous <= obs_coef))/10000))

###Mixed ranks
dyads_mod <- dyads[complete.cases(dyads[,c('RankChange', 'AI_previous')]),]
dyads_mod <- filter(dyads_mod, Focal_rank != Alter_rank)
obs_coef <- coef(glm(formula = RankChange ~ AI_previous, data = dyads_mod, family = 'binomial'))[2]

iterations = 10000
perm_coef  <- tibble(AI_previous = rep(NA, iterations))
dyads_perm <- dyads_mod
for(i in 1:iterations){
  dyads_perm %>% 
    group_by(Clan, Year) %>%
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankChange) %>%
    .[[1]] ->
    dyads_perm$RankChange
  
  dyads_perm %>% 
    glm(formula = RankChange ~ AI_previous, family = 'binomial') %>%
    coef() %>% .[2] -> perm_coef[i,1]
}

par(mfrow = c(1,1))
dyad_plot <- hist(perm_coef$AI_previous, breaks = 50, col = 'black', main = 'Changes among individuals of different rank categories', 
                  xlab = 'Effect of AI from previous year')
segments(x1 = obs_coef, x0 = obs_coef, y0 = 0, y1 = max(dyad_plot$counts), col = 'red', lwd = 2, lty = 2)
text(x = -15, y = 300, paste0('p = ', tibdim(which(perm_coef$AI_previous <= obs_coef))/10000))


###Low rankers
dyads_mod <- dyads[complete.cases(dyads[,c('RankChange', 'AI_previous')]),]
#Scale and center
dyads_mod <- dyads_mod %>% group_by(Clan, Year) %>% mutate(AI_previous = (AI_previous-mean(AI_previous))/sd(AI_previous))
#Select low rankers
dyads_mod <- filter(dyads_mod, Focal_rank == 'Low', Alter_rank == 'Low')
#observed glm coefficients
obs_coef <- coef(glm(formula = RankChange ~ AI_previous, data = dyads_mod, family = 'binomial'))[2]


iterations = 1000

#store randomized data of RankChangers
changers_perm <- matrix(nrow = tibdim(filter(dyads_mod, RankChange == 'Change')), ncol = iterations)

perm_coef  <- tibble(AI_previous = rep(NA, iterations))
dyads_perm <- dyads_mod
for(i in 1:iterations){
  dyads_perm %>% 
    group_by(Clan, Year) %>%
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankChange) %>%
    .[[1]] ->
  dyads_perm$RankChange
  
  dyads_perm %>% 
    ungroup() %>%
    glm(formula = RankChange ~ AI_previous, family = 'binomial') %>%
    coef() %>% .[2] -> 
  perm_coef[i,1]
  
  dyads_perm %>%
    ungroup() %>%
    filter(RankChange == 'Change') %>%
    select(AI_previous) %>% .[[1]] -> 
  changers_perm[,i]
}

par(mfrow = c(1,2))
dyad_plot <- hist(perm_coef$AI_previous, breaks = 50, col = 'black', 
                  main = 'Changes among low status individuals', 
                  xlab = 'Effect of AI from previous year')
segments(x1 = obs_coef, x0 = obs_coef, y0 = 0, y1 = max(dyad_plot$counts), col = 'red', lwd = 2, lty = 2)
text(x = -8, y = 500, paste0('p = ', tibdim(which(perm_coef$AI_previous <= obs_coef))/10000))

changers_plot_data <- tibble(AI_previous = as.vector(changers_perm), Network = 'Randomized', RankCategory = 'Low')
changers_plot_data <- bind_rows(changers_plot_data, 
                                tibble(AI_previous = filter(dyads_mod, RankChange == 'Change')$AI_previous, 
                                       Network = 'Observed',
                                       RankCategory = 'Low'))
changers_plot_data

###High rankers
dyads_mod <- dyads[complete.cases(dyads[,c('RankChange', 'AI_previous')]),]
#Scale and center
dyads_mod <- dyads_mod %>% group_by(Clan, Year) %>% mutate(AI_previous = (AI_previous-mean(AI_previous))/sd(AI_previous))
#Select low rankers
dyads_mod <- filter(dyads_mod, Focal_rank == 'High', Alter_rank == 'High')
#observed glm coefficients
obs_coef <- coef(glm(formula = RankChange ~ AI_previous, data = dyads_mod, family = 'binomial'))[2]
iterations = 1000

#store randomized data of RankChangers
changers_perm <- matrix(nrow = tibdim(filter(dyads_mod, RankChange == 'Change')), ncol = iterations)

perm_coef  <- tibble(AI_previous = rep(NA, iterations))
dyads_perm <- dyads_mod
for(i in 1:iterations){
  dyads_perm %>% 
    group_by(Clan, Year) %>%
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankChange) %>%
    .[[1]] ->
  dyads_perm$RankChange
  
  dyads_perm %>% 
    glm(formula = RankChange ~ AI_previous, family = 'binomial') %>%
    coef() %>% .[2] -> 
  perm_coef[i,1]
  
  dyads_perm %>%
    ungroup() %>%
    filter(RankChange == 'Change') %>%
    select(AI_previous) %>% .[[1]] -> 
  changers_perm[,i]
}


dyad_plot <- hist(perm_coef$AI_previous, breaks = 50, col = 'black',
                  main = 'Changes among high status individuals',
                  xlab = 'Effect of AI from previous year')
segments(x1 = obs_coef, x0 = obs_coef, y0 = 0, y1 = max(dyad_plot$counts), col = 'red', lwd = 2, lty = 2)
text(x = -5, y = 300, paste0('p = ', tibdim(which(perm_coef$AI_previous <= obs_coef))/10000))

changers_plot_data <- bind_rows(changers_plot_data,
                                tibble(AI_previous = as.vector(changers_perm), 
                                       Network = 'Randomized', 
                                       RankCategory = 'High'))
changers_plot_data <- bind_rows(changers_plot_data, 
                                tibble(AI_previous = filter(dyads_mod, RankChange == 'Change')$AI_previous, 
                                       Network = 'Observed',
                                       RankCategory = 'High'))

ggplot(changers_plot_data, aes(x = RankCategory, y = AI_previous, col = Network)) + 
  geom_boxplot()










##All rankers
ranks <- mutate(ranks, RankChangeBin = ifelse(RankChange != 'None', 'Change', 'NoChange'))
ranks.tab <- ranks %>% 
  group_by(Clan, Year, RankCategory) %>%
  count(Total = RankChangeBin)

ranks %>% 
  group_by(Clan, Year) %>%
  summarize(ClanSize = tibdim(ID)) %>% 
  ungroup() -> 
ranks.tab
                        
ranks.tab$Changers <- 0    
for(clan in unique(ranks.tab$Clan)){
  for(year in unique(ranks.tab[ranks.tab$Clan==clan,]$Year)){
    ranks.tab[ranks.tab$Clan==clan & ranks.tab$Year==year,]$Changers <- tibdim(filter(ranks,
                                                                             Year == year,
                                                                             Clan == clan,
                                                                             RankChange != 'None'))
  }
}       

ggplot(data = ranks.tab, aes(x = ClanSize, y = Changers/ClanSize))+
  geom_point()+
  geom_smooth(method = lm)

##High rankers
ranks %>% 
  filter(RankCategory == 'High') %>%
  group_by(Clan, Year) %>%
  summarize(ClanSize = tibdim(ID)) %>% 
  ungroup() -> 
ranks.tab.high

ranks.tab.high$Changers <- 0    
for(clan in unique(ranks.tab$Clan)){
  for(year in unique(ranks.tab.high[ranks.tab$Clan==clan,]$Year)){
    ranks.tab.high[ranks.tab.high$Clan==clan & ranks.tab.high$Year==year,]$Changers <- tibdim(filter(ranks,
                                                                                      Year == year,
                                                                                      Clan == clan,
                                                                                      RankCategory == 'High',
                                                                                      RankChange != 'None'))
  }
}    

ggplot(data = ranks.tab.high, aes(x = ClanSize, y = Changers/ClanSize))+
  geom_point()+
  geom_smooth(method = lm)


##Low rankers
ranks %>% 
  filter(RankCategory == 'Low') %>%
  group_by(Clan, Year) %>%
  summarize(ClanSize = tibdim(ID)) %>% 
  ungroup() -> 
  ranks.tab.low

ranks.tab.low$Changers <- 0    
for(clan in unique(ranks.tab.low$Clan)){
  for(year in unique(ranks.tab.low[ranks.tab$Clan==clan,]$Year)){
    ranks.tab.low[ranks.tab.low$Clan==clan & ranks.tab.low$Year==year,]$Changers <- tibdim(filter(ranks,
                                                                                                     Year == year,
                                                                                                     Clan == clan,
                                                                                                     RankCategory == 'Low',
                                                                                                     RankChange != 'None'))
  }
}    

ggplot(data = ranks.tab.low, aes(x = ClanSize, y = Changers/ClanSize))+
  geom_point()+
  geom_smooth(method = lm)

summary(lm(data =ranks.tab.low, I(Changers/ClanSize) ~ ClanSize))

library(aniDom)
clan = 'talek'
year <- 1993
ids <- tibble(ID = filter(ranks, Year == year, Clan == clan)$ID)
aggs <- filter(aggsWinner, Agg %in% ids$ID, Recip %in% ids$ID, Year == year)
r <- left_join(ids, filter(ranks, Year == year), by = 'ID')$Rank
a <- plot_hierarchy_shape(identity = ids$ID, rank = r, winners = aggs$Agg, losers = aggs$Recip)

a <- elo_scores(aggs$Agg, aggs$Recip, randomise = T)
elo <- tibble(ID = row.names(a), elo = as.vector(a[,1]))
elo$rank <- left_join(elo, filter(ranks, Year == year, Clan == clan), by = 'ID')$Rank
plot(data = elo, elo ~ rank)
plot_ranks(a , plot.CIs = T)

winners.rank <- rank[match(winners, identity)]
losers.rank <- rank[match(losers, identity)]
xx <- winners.rank - losers.rank
x <- 1:(max(abs(xx)))
y <- rep(NA, length(x))
CI.upper <- y
CI.lower <- y
for (i in 1:length(x)) {
  y[i] <- sum(xx == -x[i])/sum(abs(xx) == x[i])
  CI.upper[i] <- y[i] + sqrt(y[i] * (1 - y[i])/sum(abs(xx) == 
                                                     x[i])) + 0.5/sum(abs(xx) == x[i])
  CI.upper[i] <- min(CI.upper[i], 1)
  CI.lower[i] <- y[i] - sqrt(y[i] * (1 - y[i])/sum(abs(xx) == 
                                                     x[i])) - 0.5/sum(abs(xx) == x[i])
  CI.lower[i] <- max(CI.lower[i], 0)
}
x <- x[!is.na(y)]
y <- y[!is.na(y)]
plot(x, y, xlab = "Difference in rank", ylab = "Probability that higher rank wins", 
     ylim = c(min(0.5, min(y)), 1), pch = 20, cex = 2)
arrows(x, CI.lower, x, CI.upper, length = 0.1, angle = 90, 
       code = 3)
if (fitted) {
  l <- loess(y ~ x)
  lines(x, l$fitted, col = "red", lwd = 2)
}
invisible(data.frame(Rank.diff = x, Prob.dom.win = y, CI.upper = CI.upper, 
                     CI.lower = CI.lower))



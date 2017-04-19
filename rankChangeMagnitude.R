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

########Determine rank before change######

########################################

#####Calculate magnatude of rank change####
ranks$RankDiff <- NA


calc_diff <- function(rank_row, ranks){
  return(rank_row[6] %>% as.numeric %>% round(digits = 4) - 
           filter(ranks, IDold %in% rank_row[3], Year %in% rank_row[2])$stan.rank %>% as.numeric() %>% round(digits = 4) 
  )
}


ranks$RankDiff <- apply(ranks, 1, calc_diff, ranks)
###############################################

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
iterations <- 1000
ai_mod <- ranks %>% lm(formula = RankDiff ~ ai_deg + stan.rank)
ai_top3_mod <- ranks %>% lm(formula = RankDiff ~ ai_top3_deg + stan.rank)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                   stan.rank_ai = ai_mod %>% coef() %>% .[3],
                   ai_top3_deg = ai_top3_mod %>% coef() %>% .[2],
                   stan.rank_ai3 = ai_top3_mod %>% coef() %>% .[3]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    stan.rank_ai = rep(NA, iterations),
                    ai_top3_deg = rep(NA, iterations),
                    stan.rank_ai3 = rep(NA, iterations))


ranks.perm <- ranks
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks$Clan, ranks$Year)

for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiff) %>% .[[1]] ->
    ranks.perm$RankDiff

  ranks.perm %>%
    lm(formula = RankDiff ~ ai_deg + stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,1:2]
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_top3_deg+ stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,3:4]
}

par(mfrow = c(2,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black', xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red')
text(x = -.015, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/1000))

stan.rank_ai_plot <- hist(perm_coef$stan.rank_ai, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$stan.rank_ai, x0 = obs_coef$stan.rank_ai, y0 = 0, y1 = max(stan.rank_ai_plot$counts), col = 'red')
text(x = -.012, y = 30, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai >= obs_coef$stan.rank_ai))/1000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red')
text(x = -.042, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/1000))

stan.rank_ai3_plot <- hist(perm_coef$stan.rank_ai3, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$stan.rank_ai3, x0 = obs_coef$stan.rank_ai3, y0 = 0, y1 = max(stan.rank_ai3_plot$counts), col = 'red')
text(x = -.012, y = 40, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai3 >= obs_coef$stan.rank_ai3))/1000))

###############################################

################Do permutations within rank category###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

iterations <- 1000
ai_mod <- ranks %>% lm(formula = RankDiff ~ ai_deg + stan.rank)
ai_top3_mod <- ranks %>% lm(formula = RankDiff ~ ai_top3_deg + stan.rank)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            stan.rank_ai = ai_mod %>% coef() %>% .[3],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2],
                            stan.rank_ai3 = ai_top3_mod %>% coef() %>% .[3]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    stan.rank_ai = rep(NA, iterations),
                    ai_top3_deg = rep(NA, iterations),
                    stan.rank_ai3 = rep(NA, iterations))

ranks.perm <- ranks
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks$Clan, ranks$Year, ranks$RankCategory)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiff) %>% .[[1]] ->
    ranks.perm$RankDiff
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_deg + stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,1:2]
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_top3_deg+ stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,3:4]
}

par(mfrow = c(2,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black', xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red')
text(x = -.015, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/1000))

stan.rank_ai_plot <- hist(perm_coef$stan.rank_ai, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$stan.rank_ai, x0 = obs_coef$stan.rank_ai, y0 = 0, y1 = max(stan.rank_ai_plot$counts), col = 'red')
text(x = -.012, y = 30, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai >= obs_coef$stan.rank_ai))/1000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red')
text(x = -.042, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/1000))

stan.rank_ai3_plot <- hist(perm_coef$stan.rank_ai3, breaks = 50, col = 'black',xlim = c(-.05, .05), ylim = c(0,100))
segments(x1 = obs_coef$stan.rank_ai3, x0 = obs_coef$stan.rank_ai3, y0 = 0, y1 = max(stan.rank_ai3_plot$counts), col = 'red')
text(x = -.012, y = 40, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai3 >= obs_coef$stan.rank_ai3))/1000))

################Test effect of AI for high rankers only###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

iterations <- 1000
ai_mod <- filter(ranks, RankCategory == 'High') %>% lm(formula = RankDiff ~ ai_deg + stan.rank)
ai_top3_mod <- filter(ranks, RankCategory == 'High') %>% lm(formula = RankDiff ~ ai_top3_deg + stan.rank)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            stan.rank_ai = ai_mod %>% coef() %>% .[3],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2],
                            stan.rank_ai3 = ai_top3_mod %>% coef() %>% .[3]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    stan.rank_ai = rep(NA, iterations),
                    ai_top3_deg = rep(NA, iterations),
                    stan.rank_ai3 = rep(NA, iterations))

ranks.perm <- filter(ranks, RankCategory == 'High')
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks.perm$Clan, ranks.perm$Year)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiff) %>% .[[1]] ->
    ranks.perm$RankDiff
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_deg + stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,1:2]
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_top3_deg+ stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,3:4]
}

par(mfrow = c(2,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red')
text(x = -.015, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/1000))

stan.rank_ai_plot <- hist(perm_coef$stan.rank_ai, breaks = 50, col = 'black')
segments(x1 = obs_coef$stan.rank_ai, x0 = obs_coef$stan.rank_ai, y0 = 0, y1 = max(stan.rank_ai_plot$counts), col = 'red')
text(x = -.012, y = 30, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai >= obs_coef$stan.rank_ai))/1000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red')
text(x = -.042, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/1000))

stan.rank_ai3_plot <- hist(perm_coef$stan.rank_ai3, breaks = 50, col = 'black')
segments(x1 = obs_coef$stan.rank_ai3, x0 = obs_coef$stan.rank_ai3, y0 = 0, y1 = max(stan.rank_ai3_plot$counts), col = 'red')
text(x = -.012, y = 40, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai3 >= obs_coef$stan.rank_ai3))/1000))

################Test effect of AI for low rankers only###############

###rank category
ranks <- ranks %>%
  mutate(RankCategory = cut(stan.rank, breaks = 2, labels = c('Low', 'High')))

iterations <- 1000
ai_mod <- filter(ranks, RankCategory == 'Low') %>% lm(formula = RankDiff ~ ai_deg + stan.rank)
ai_top3_mod <- filter(ranks, RankCategory == 'Low') %>% lm(formula = RankDiff ~ ai_top3_deg + stan.rank)

obs_coef <- as_tibble(cbind(ai_deg = ai_mod %>% coef() %>% .[2],
                            stan.rank_ai = ai_mod %>% coef() %>% .[3],
                            ai_top3_deg = ai_top3_mod %>% coef() %>% .[2],
                            stan.rank_ai3 = ai_top3_mod %>% coef() %>% .[3]))
perm_coef <- tibble(ai_deg = rep(NA, iterations), 
                    stan.rank_ai = rep(NA, iterations),
                    ai_top3_deg = rep(NA, iterations),
                    stan.rank_ai3 = rep(NA, iterations))

ranks.perm <- filter(ranks, RankCategory == 'Low')
###Select variables to permute within##
ranks.perm$PermCategory <- paste(ranks.perm$Clan, ranks.perm$Year)
for(i in 1:iterations){
  ranks.perm %>% 
    group_by(PermCategory) %>% 
    sample_frac(replace = F) %>% 
    ungroup() %>% 
    select(RankDiff) %>% .[[1]] ->
    ranks.perm$RankDiff
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_deg + stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,1:2]
  
  ranks.perm %>%
    lm(formula = RankDiff ~ ai_top3_deg+ stan.rank) %>%
    coef() %>% .[2:3] -> perm_coef[i,3:4]
}

par(mfrow = c(2,2))
ai_plot <- hist(perm_coef$ai_deg, breaks = 50, col = 'black')
segments(x1 = obs_coef$ai_deg, x0 = obs_coef$ai_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red')
text(x = -.015, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_deg >= obs_coef$ai_deg))/1000))

stan.rank_ai_plot <- hist(perm_coef$stan.rank_ai, breaks = 50, col = 'black')
segments(x1 = obs_coef$stan.rank_ai, x0 = obs_coef$stan.rank_ai, y0 = 0, y1 = max(stan.rank_ai_plot$counts), col = 'red')
text(x = -.012, y = 30, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai >= obs_coef$stan.rank_ai))/1000))

ai_top3_plot <- hist(perm_coef$ai_top3_deg, breaks = 50, col = 'black')
segments(x1 = obs_coef$ai_top3_deg, x0 = obs_coef$ai_top3_deg, y0 = 0, y1 = max(ai_top3_plot$counts), col = 'red')
text(x = -.042, y = 30, paste0('p = ', tibdim(which(perm_coef$ai_top3_deg >= obs_coef$ai_top3_deg))/1000))

stan.rank_ai3_plot <- hist(perm_coef$stan.rank_ai3, breaks = 50, col = 'black')
segments(x1 = obs_coef$stan.rank_ai3, x0 = obs_coef$stan.rank_ai3, y0 = 0, y1 = max(stan.rank_ai3_plot$counts), col = 'red')
text(x = -.012, y = 40, paste0('p = ', tibdim(which(perm_coef$stan.rank_ai3 >= obs_coef$stan.rank_ai3))/1000))



####################AI by dyad - are dyads that associate less more likely to change ranks?########
dyads <- group_by(ranks, Clan, Year) %>% {expand.grid(.$ID, .$ID)}
##############################################
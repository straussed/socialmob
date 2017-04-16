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
#sink(file = 'RankChangePermutations_Output.txt')

#####Tibble length function######
tibdim <- function(tib, margin = 1){
  if(margin == 1){
    length(tib[[1]])
  }else{
    length(tib)
  }
}
###################################


############################Add rank change info to ranks##############
for(row in 1:length(ranks[,1])){
  ranks.temp <- filter(ranks, Clan == ranks[row,'Clan'], Year == ranks[row,'Year'])
  if(ranks[row,'Rank'] == ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'None'
  }else if(ranks[row,'Rank'] < ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'Up'
  }else if(ranks[row,'Rank'] > ranks.temp[ranks.temp$IDold == ranks[row, 'ID'],'Rank']){
    ranks[row,'RankChange'] <- 'Down'
  }
}
ranks$RankChange <- factor(ranks$RankChange, levels = c('None', 'Up', 'Down'))
#############################################

#####Calculate magnatude of rank change####
ranks$RankDiff <- NA

calc_diff <- function(rank_row, ranks){
  return(as.numeric(rank_row[6]) - as.numeric(filter(ranks, IDold %in% rank_row[3], Year %in% rank_row[2])$stan.rank))
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

################Do permutations#############

############################################


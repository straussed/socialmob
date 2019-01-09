################################################################################
#                  Add rank change info to ranks dataframe                     #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)
load('6.ranks_with_coalition_data.RData')
names(ranks)


assign_rank_change <- function(a){
  rc <- data.frame(RankChange = rep(NA, nrow(a)),
               RankDiff = rep(NA, nrow(a)))
  stan_rank_unit <- 2/(nrow(a)-1)
  for(i in 1:length(a$ID)){
    id <- a$ID[i]
    pos <- which(a$ID == id)
    old_pos <- which(a$OldOrder == id)
    if((pos == 1 | pos == nrow(a)) & a$OldOrder[pos] == id){
      rc[i,1] <- 'None'
      rc[i,2] <- 0
      rc[i,3] <- 0
    }else{
      passers <- a$ID[1:pos][!a$ID[1:pos] %in% a$OldOrder[1:old_pos]]
      passees <- a$OldOrder[1:old_pos][!a$OldOrder[1:old_pos] %in% a$ID[1:pos]]
      if(length(passers) & length(passees)){
        rc[i,1] <- 'Both'
        rc[i,2] <- (length(passees) - length(passers))*stan_rank_unit
        rc[i,3] <- (length(passees) - length(passers))
      }else if(length(passers)){
        rc[i,1] <- 'Down'
        rc[i,2] <- -length(passers)*stan_rank_unit
        rc[i,3] <- -length(passers)
      }else if(length(passees)){
        rc[i,1] <- 'Up'
        rc[i,2] <- length(passees)*stan_rank_unit
        rc[i,3] <- length(passees)
      }else{
        rc[i,1] <- 'None'
        rc[i,2] <- 0
        rc[i,3] <- 0
      }
    }
  }
  return(rc)
}


ranks$RankChange <- NA
ranks$RankDiff <- NA
ranks$RankDiffAbs <- NA
for(clan in unique(ranks$Clan)){
  for(year in unique(ranks[ranks$Clan == clan,]$Year)){
    ranks[ranks$Clan == clan & ranks$Year == year,c('RankChange',
                                                    'RankDiff',
                                                    'RankDiffAbs')] <-
      ranks[ranks$Clan == clan & ranks$Year == year,] %>% 
      assign_rank_change()
    
  }
}

ranks$RankChange <- factor(ranks$RankChange, levels = c('Down', 'None', 'Up', 'Both'))

sum(abs(ranks$RankDiffAbs)) ##460 initial ranks determined by Elo
                            ##278 standard initial ranks


save(ranks, file = '8.ranks_with_rank_change.RData')

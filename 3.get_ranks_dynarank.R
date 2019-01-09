################################################################################
#                       Get ranks using DynaRank                               #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################


rm(list = ls())
options(stringsAsFactors = FALSE)

library(DynaRankR)
library(dplyr)


setwd('~/Documents/Research/socialmob/ranksDynaRank/')

load('2.hyena_data.RData')

set.seed(1989)


################################################################################
#
#    Talek
#
talek_ranks <- informed_matreorder(contestants = contestants.talek,
                                   convention = 'mri',
                                   n = 100,
                                   shuffles = 30,
                                   require.corroboration = TRUE,
                                   initial.ranks = initial.ranks.talek,
                                   interactions = interactions.talek)

plot_ranks(talek_ranks)

################################################################################
#
#    North
#

north_ranks <- informed_matreorder(contestants = contestants.north,
                        convention = 'mri',
                        n = 100,
                        shuffles = 30,
                        require.corroboration = TRUE,
                        initial.ranks = initial.ranks.north,
                        interactions = interactions.north)


plot_ranks(north_ranks)

################################################################################
#
#    South
#

south_ranks <- informed_matreorder(contestants = contestants.south,
                        convention = 'mri',
                        n = 100,
                        shuffles = 30,
                        require.corroboration = TRUE,
                        initial.ranks = initial.ranks.south,
                        interactions = interactions.south)

plot_ranks(south_ranks)

################################################################################
#
#    Happy Zebra
#

hz_ranks <- informed_matreorder(contestants = contestants.hz,
                     convention = 'mri',
                     n = 100,
                     shuffles = 30,
                     require.corroboration = TRUE,
                     initial.ranks = initial.ranks.hz,
                     interactions = interactions.hz)

plot_ranks(hz_ranks)




## Combine data
talek_ranks$clan <- 'talek'
north_ranks$clan <- 'serena.n'
south_ranks$clan <- 'serena.s'
hz_ranks$clan <- 'happy.zebra'

#remove last year
ranks <- rbind(filter(talek_ranks, period <= 2014), 
               filter(north_ranks, period <= 2015), 
               filter(south_ranks, period <= 2015), 
               filter(hz_ranks,period <= 2014))

save(ranks, file = '4.ranks.RData')

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

for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    ####Coalition network
    cnet <- matrix(nrow = tibdim(yearRanks), ncol = tibdim(yearRanks), dimnames = list(yearRanks$ID, yearRanks$ID), data = 0)
    coalsTemp <- filter(aggsFull, Year == year, Clan == clan, Group != '')
    coalPartners <- cbind(expand.grid(yearRanks$ID, yearRanks$ID), rep(0))
    names(coalPartners) <- c('Focal', 'Alter', 'NumCoals')
    coalPartners$Focal <- as.character(coalPartners$Focal)
    coalPartners$Alter <- as.character(coalPartners$Alter)
    if(!tibdim(coalsTemp)){next}
    for(row in 1:tibdim(coalsTemp)){
      focal <- coalsTemp[row,]$Agg
      for(alter in strsplit(coalsTemp[row,]$Group, ',')[[1]][-1]){
        if(alter %in% yearRanks$ID & alter != focal & focal %in% yearRanks$ID){
          coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,'NumCoals'] <- coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,'NumCoals'] + 1
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
                    coal_top3_wtd = rep(NA)
)

for(row in 1:tibdim(ranks)){
  year <- ranks[row,]$Year
  clan <- ranks[row,]$Clan
  focal <- ranks[row,]$ID
  if(exists(paste('coalNet', clan, year, sep = '_'))){
    coal_net <- get(paste('coalNet', clan, year, sep = '_'))
    if(tibdim(coal_net)){
      ranks[row,]$coal_deg <- sum(coal_net[focal,])
      ranks[row,]$coal_top3_deg <- coal_net[focal,order(coal_net[focal,], decreasing = T)[1:3]] %>% 
        sum() 
    } 
  }
}


################Do permutations all ranks together###############

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
  ranks_mod <- filter(ranks, RankChange != 'None', RankCategory %in% category)
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                                 tibble(coal_all_deg = ranks_mod$coal_deg,
                                        coal_top_deg = ranks_mod$coal_top3_deg,
                                        rank_change = as.character(ranks_mod$RankChange),
                                        rank_category = category,
                                        network = 'Observed'))
  
  iterations <- 1000
  
  rank_perm_change <- matrix(NA, nrow=tibdim(ranks_mod), ncol=iterations)

  coal_mod <- ranks_mod %>% glm(formula = RankChange ~ coal_deg, family = binomial)
  coal_top3_mod <- ranks_mod %>% glm(formula = RankChange ~ coal_top3_deg, family = binomial)
  
  obs_coef <- as_tibble(cbind(coal_deg = coal_mod %>% coef() %>% .[2],
                              coal_top3_deg = coal_top3_mod %>% coef() %>% .[2]))
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
      select(RankChange) %>% .[[1]] ->
      ranks.perm$RankChange
    
    ranks.perm %>%
      glm(formula = RankChange ~ coal_deg, family = binomial) %>%
      coef() %>% .[2] -> perm_coef[i,1]
    
    ranks.perm %>%
      glm(formula = RankChange ~ coal_top3_deg, family = binomial) %>%
      coef() %>% .[2] -> perm_coef[i,2]
    
    rank_perm_change[,i] <- as.character(ranks.perm$RankChange)
  }
  
  coal_perm_stacked <- bind_rows(coal_perm_stacked,
                               tibble(coal_all_deg = rep(ranks_mod$coal_deg, iterations),
                                      coal_top_deg = rep(ranks_mod$coal_top3_deg, iterations),
                                      rank_change = as.vector(rank_perm_change),
                                      rank_category = category,
                                      network = 'Permuted'))
  
  coal_plot <- hist(perm_coef$coal_deg, breaks = 50, col = 'black',
                    main = paste0(category, '--All'),
                    xlab = 'Coefficient estimate')
  segments(x1 = obs_coef$coal_deg, x0 = obs_coef$coal_deg, y0 = 0, y1 = max(ai_plot$counts), col = 'red', lty=2, lwd=2)
  text(x = -.018, y = 500, paste0('p = ', tibdim(which(perm_coef$coal_deg >= obs_coef$coal_deg))/1000))
  
  coal_top3_plot <- hist(perm_coef$coal_top3_deg, breaks = 50, col = 'black',
                         main = paste0(category, '--Top'),
                         xlab = 'Coefficient estimate')
  segments(x1 = obs_coef$coal_top3_deg, x0 = obs_coef$coal_top3_deg, y0 = 0, y1 = max(coal_top3_plot$counts), col = 'red', lty=2, lwd=2)
  text(x = -.042, y = 300, paste0('p = ', 2*tibdim(which(perm_coef$coal_top3_deg >= obs_coef$coal_top3_deg))/1000))
}


facet_labels = c(High = 'High ranking individuals', Low = 'Low ranking individuals')
colors = c('dodgerblue1', 'grey87')

top_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change, y = coal_top_deg, fill = network)) + 
  geom_boxplot(alpha = .7)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with top allies')
top_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold', size = 12,family = 'Trebuchet MS'),
        axis.title = element_text(size = 14,family = 'Trebuchet MS'),
        axis.text = element_text(size = 12, family = 'Trebuchet MS'),
        legend.text = element_text(size = 12, family = 'Trebuchet MS'))



all_plot <- ggplot(data = coal_perm_stacked, aes(x = rank_change, y = coal_all_deg, fill = network)) + 
  geom_boxplot(alpha = .7)+
  scale_fill_manual(values = colors, name = '') +
  theme_bw() + 
  facet_grid(. ~ rank_category, labeller = labeller(rank_category = facet_labels))+
  xlab('Direction of rank change')+
  ylab('Total degree of coalitions with all allies')
all_plot +
  theme(strip.background = element_rect(fill = alpha(colors[1], .7)),
        strip.text.x = element_text(face = 'bold'))



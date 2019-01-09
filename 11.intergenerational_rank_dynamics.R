################################################################################
#                   Intergenerational rank dynamics                            #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################

rm(list =ls())
library(dplyr)
library(grid)
library(ggplot2)
options(stringsAsFactors = FALSE)

load('8.ranks_with_rank_change.RData')

#demographic information
#hyenas <- read.csv('/Volumes/Holekamp/code_repository/R/1_output_tidy_tbls/tblHyenas.csv')
hyenas <- read.csv('0.rawdata/tblHyenas.csv')
hyenas[hyenas$clan == 'serena s',]$clan <- 'south'
hyenas[hyenas$clan == 'serena n',]$clan <- 'north'
hyenas[hyenas$clan == 'happy zebra',]$clan <- 'hz'


#Fix id loda -> luda
hyenas[hyenas$id == 'loda',]$id <- 'luda'

#Fix bd
hyenas[hyenas$id == 'bd',]$birthdate <- hyenas[hyenas$id == 'bd',]$first.seen

#Fix tru
hyenas[hyenas$id == 'tru',]$mom <- 'bor'


#remove hyena 44 who disappeared before our data begin
hyenas <- hyenas[-which(hyenas$id == '44'),]

#cash appears twice
hyenas <- filter(hyenas, !(id == 'cash' & clan == 'north'))

hyenas$birthdate <- as.Date(hyenas$birthdate)
hyenas$death.date <- as.Date(hyenas$death.date)
hyenas$disappeared <- as.Date(hyenas$disappeared)
tblHyenas <- hyenas
names(tblHyenas) <- c('ID', 'Last.Updated', 'SampleID', 'Eartage', 'Name', 'PrevID',
                      'Sex', 'AgeClass', 'Status', 'FirstSeen', 'DenGrad', 'Disappeared',
                      'Mom', 'Birthdate','NumLittermates',
                      'LitRank', 'ArrivedDen', 'LeaveDen', 'Fate', 'MortalitySource',
                      'DeathDate', 'Weaned', 'Clan', 'Park', 'Notes')





####Amplifications of small rank differences###
matriarchs <- c('kb', 'dj', '03', 'coch', '79',
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
  summarize(ymin = min(Rank), ymax = max(Rank), max.stan = max(Stan.Rank),
            mean.stan = mean(Stan.Rank))

ribbon[ribbon$matriline == '79' & ribbon$Year > 2007,c('ymin', 'ymax')] <- NA


matriline_labels <- data.frame(x.start = rep(1987), 
                               matriline = c('kb', 'dj', '03', '79', 'other'),
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
             size = 3,
             color = c('black', 'black', 'black', 'black', 'black'))+
  geom_point(data = data.frame(hex = c('hex', 'hex'), x = c(2005, 2014), y = c(19, 40)),
             aes(x = x, y = y), inherit.aes = FALSE,
             fill = desat(cols[2], sat = 0.5), shape = 23, size =2.5, stroke = 1.4,
             col = darken(cols[2]))+
  scale_color_manual(values = darken(cols)) + 
  scale_fill_manual(values = desat(cols, sat = 0.5))+
  ggtitle('a)')+
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.title = element_text(size= 11),
        plot.margin = unit(c(5.5, 10, 5.5, 5.5), units = 'pt'))


##Track difference between females over time##
matriline.pairs <- data.frame(id1 = c('kb', 'dj', '03'),
                              id2 = c('dj', '03', '79'))


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
                    '79' = get.descendants('79'),
                    dj = get.descendants('dj'),
                    kb = get.descendants('kb'))

yearly.data.no.change <- list()
ids.that.change <- unique(filter(ranks, RankChange != 'None')$ID)
counter <- 1
for(row in 1:nrow(matriline.pairs)){
  for(year in sort(unique(ranks$Year))){
    id1.desc <- descendants[[matriline.pairs[row,1]]]
    id2.desc <- descendants[[matriline.pairs[row,2]]]
    
    mean.Stan.Rank.diff <- filter(ranks, Year == year, ID %in% id1.desc,
                                  !ID %in% ids.that.change)$Stan.Rank %>%
      mean() - 
      mean(filter(ranks, Year == year, ID %in% id2.desc,
                  !ID %in% ids.that.change)$Stan.Rank)
    
    mean.abs.rank.diff <- filter(ranks, Year == year, ID %in% id2.desc,
                                 !ID %in% ids.that.change)$Rank %>%
      mean() - 
      mean(filter(ranks, Year == year, ID %in% id1.desc,
                  !ID %in% ids.that.change)$Rank)
    
    
    yearly.data.no.change[[counter]] <- data.frame(id1 = matriline.pairs[[row,1]],
                                                   id2 = matriline.pairs[[row,2]],
                                                   mean.Stan.Rank.diff,
                                                   mean.abs.rank.diff,
                                                   year)
    counter <- counter+1
  }
}
rank.diff.over.time.no.change <- do.call(rbind, yearly.data.no.change)

rank.diff.over.time.no.change$id1 <- factor(rank.diff.over.time.no.change$id1, 
                                            levels = c('kb', 'dj', '03'))

rank.diff.over.time.no.change$id2 <- factor(rank.diff.over.time.no.change$id2, 
                                            levels = c('dj', '03', '79'))

mean.rank.diff <- 
  ggplot(rank.diff.over.time.no.change, 
         aes(x = year, y = mean.abs.rank.diff, 
             col = id1, fill = id2))+
  geom_smooth(method = 'lm', se = FALSE, col = 'black', 
              aes(x = year, y = mean.abs.rank.diff), inherit.aes = F)+
  geom_jitter(shape = 21, stroke = 1.5, size = 2, width = 1, height = 1)+
  theme_classic()+
  ylab('Difference in rank of adjacent matrilines')+
  xlab('Year')+
  scale_fill_manual(name = 'Difference between\nadjacent matrilines',
                    labels = c('kb-dj', 'dj-03', '03-79'), 
                    values = cols[c(3,1,2)])+
  scale_color_manual(name = 'Difference between\nadjacent matrilines',
                     labels = c('kb-dj', 'dj-03', '03-79'),
                     values = cols[c(4,3,1)])+
  theme(legend.position = c(0.30, 0.75),
        legend.title.align = 0.5,
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.title = element_text(size= 11))+
  ggtitle('b)')


#Figure 5
pdf("plots/11.Fig5.pdf",
    height = 7,
    width = 4.5)
gridExtra::grid.arrange(matriline.colors, mean.rank.diff, nrow = 2, heights = c(5,3))
dev.off()


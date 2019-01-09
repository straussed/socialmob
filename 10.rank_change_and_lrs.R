################################################################################
#    Predict consequence of rank change on lifetime reproductive success       #
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
set.seed(1989)

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

ranks %>% group_by(ID) %>% summarize(mean.rank = mean(Stan.Rank), 
                                     total.cubs = sum(Surv2),
                                     mean.cubs = mean(Surv2),
                                     mean.ars = mean(ARS)) -> lrs

lrs <- filter(lrs, ID %in% filter(tblHyenas, !is.na(Disappeared))$ID)
lrs$repro.lifespan <- as.numeric(left_join(lrs, tblHyenas, by = 'ID')$Disappeared - 
                                   left_join(lrs, tblHyenas, by = 'ID')$Birthdate - 2*365)/365

survive_to_4 <- filter(tblHyenas, !is.na(Disappeared), Disappeared - Birthdate >= 4*365)

nrow(filter(lrs, ID %in% survive_to_4$ID))
pois.quad.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ poly(mean.rank, 2), family = 'poisson')
summary(pois.quad.lrs)

pois.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ mean.rank, family = 'poisson')
summary(pois.lrs)
pois.exp.lrs <- glm(data = filter(lrs, ID %in% survive_to_4$ID), formula = total.cubs ~ I(exp(mean.rank)), family = 'poisson')
summary(pois.exp.lrs)

##Mean lifetime reproductive success
mean(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)
sd(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)/sqrt(nrow(filter(lrs, ID %in% survive_to_4$ID)))
min(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)
max(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)

AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)
AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[1,2] - AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[3,2]
AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[2,2] - AIC(pois.quad.lrs, pois.lrs, pois.exp.lrs)[3,2]


##predict consequence of rank change on lrs based on models

predicted.lrs.final.rank <- 
  exp(predict(pois.exp.lrs, newdata = data.frame(mean.rank = ranks$Stan.Rank)))
predicted.lrs.initial.rank <- 
  exp(predict(pois.exp.lrs, 
              newdata = data.frame(mean.rank = ranks$Stan.Rank - ranks$RankDiff)))

ranks$lrs.diff <- predicted.lrs.final.rank - predicted.lrs.initial.rank


######Plots


##Figure 4a
pois.exp.pred <- predict.glm(pois.exp.lrs, se.fit = TRUE, type = "response")
lrs.mean <- mean(filter(lrs, ID %in% survive_to_4$ID)$total.cubs)

lrs.plot <- ggplot(filter(lrs, ID %in% survive_to_4$ID), aes(x = mean.rank, y = total.cubs, label = ID))+
  geom_point(size =2) +
  geom_line(aes(y = pois.exp.pred$fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), size = 1.2, col = 'firebrick')+
  geom_line(aes(y = pois.exp.pred$fit + 2*pois.exp.pred$se.fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), lty = 2, size = 1, col = 'firebrick')+
  geom_line(aes(y = pois.exp.pred$fit - 2*pois.exp.pred$se.fit, x = filter(lrs, ID %in% survive_to_4$ID)$mean.rank), lty = 2, size = 1, col = 'firebrick')+
  theme_classic()+
  theme(plot.margin = unit(x = c(5.5,16.5,5.5,5.5),units = 'pt'))+
  scale_x_continuous(labels = c('Lowest', 'Highest'), breaks = c(-1, 1))+
  xlab('Mean lifetime rank') + 
  ylab('Total cubs produced')+
  ggtitle('a)') +
  theme(axis.title = element_text(size = 11), plot.title = element_text(size = 11),
        axis.text = element_text(size = 11))

##Figure 4b

predicted.lrs.final.rank <- 
  exp(predict(pois.exp.lrs, newdata = data.frame(mean.rank = ranks$Stan.Rank)))
predicted.lrs.initial.rank <- 
  exp(predict(pois.exp.lrs, 
              newdata = data.frame(mean.rank = ranks$Stan.Rank - ranks$RankDiff)))

ranks$lrs.ratio <- predicted.lrs.final.rank/predicted.lrs.initial.rank
ranks$lrs.diff <- predicted.lrs.final.rank - predicted.lrs.initial.rank

predict.consequences <- ggplot(filter(ranks, lrs.diff != 0), aes(x = lrs.diff/lrs.mean))+
  geom_histogram(fill = 'firebrick')+
  theme_classic()+
  xlab('Predicted change in LRS')+
  ylab('')+
  ggtitle('b)')+
  theme(axis.title = element_text(size = 11), plot.title = element_text(size = 11),
        axis.text = element_text(size = 11))

## Figure 4c
exp.change.lrs <- ggplot(filter(ranks, lrs.ratio != 1), aes(x = Stan.Rank, y = RankDiffAbs, fill = lrs.ratio, size = exp(abs(1-lrs.ratio))))+
  geom_jitter(shape = 21, col = 'grey80')+
  theme_classic()+
  geom_hline(yintercept = 0, lty = 2)+
  xlab('Standardize rank after change')+
  ylab('Number of positions moved') +
  ylim(-11, 22)+
  scale_fill_gradientn(colors = c("firebrick", "white","dodgerblue"),
                       values = scales::rescale(c(0.5, 1 , 3.5)),
                       guide = "colorbar", limits = c(0.5, 3.5))+
  guides(size = FALSE, 
         fill= guide_colorbar(title = 'New LRS ÷ Old LRS', direction = 'horizontal',
                              title.position = 'top', title.hjust = 0.5, barwidth = 8, nbin = 100))+
  theme(legend.position = c(0.3, 0.8), legend.key.size = unit(c(5, 5), unit = 'pt'),
        legend.title = element_text(size = 11), legend.text = element_text(size = 11),
        axis.title = element_text(size = 11), plot.title = element_text(size = 11),
        axis.text = element_text(size = 11))+
  ggtitle('c)')

pdf("plots/10.Fig4.pdf",
    height = 5,
    width = 4.5)
source('0.multiplot.R')
multiplot(lrs.plot, predict.consequences,
          exp.change.lrs, layout = matrix(c(1,2,3,3), byrow = TRUE, nrow = 2))
dev.off()


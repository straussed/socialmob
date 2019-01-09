################################################################################
#             Model and plot effect of coalitions on rank change               #
#                                                                              #
#                                                                              #
#                           By Eli Strauss                                     #
#                                                                              #
#                           December 2018                                      #
################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(grid)
load('8.ranks_with_rank_change.RData')
set.seed(1989)


####Examine coalition degree in top3 partners wtd by obs####
ranks.perm <- na.omit(dplyr::select(ranks, ID, Clan, Year, RankDiffAbs, coal_top3_deg, obs_counts))
top3_mod <- ranks.perm %>% lm(formula = RankDiffAbs ~ coal_top3_deg + log(offset(obs_counts)))
top3_poly <- ranks.perm %>% lme4::lmer(formula = RankDiffAbs ~ poly(coal_top3_deg,2) + (1|ID) + log(offset(obs_counts)))
summary(top3_poly)
mod_no_coal <- ranks.perm %>% blme::blmer(formula = RankDiffAbs ~ (1|ID)+ log(offset(obs_counts)))
aic <- MuMIn::AICc(top3_mod, top3_poly, mod_no_coal)
aic$AICc[2] - aic$AICc[3]

####Analysis
predictions.list <- list()
betas.list <- list(data.frame(Intercept = top3_poly@beta[1], coal_top3_deg = top3_poly@beta[2], coal_top3_poly = top3_poly@beta[3], offset = top3_poly@beta[4]))

mod.perm <- top3_poly
for(i in 1:999){
  ##If fit is singular, try again
  mod.perm@theta <- 0
  while(lme4::isSingular(mod.perm)){
    ranks.perm <- ranks.perm %>% group_by(Clan, Year) %>%
      mutate(RankDiffAbs = sample(RankDiffAbs, replace = FALSE))
    mod.perm <- ranks.perm %>% lme4::lmer(formula = RankDiffAbs ~ poly(coal_top3_deg, 2) + (1|ID) + log(offset(obs_counts)))
  }
  
  betas.list[[i+1]] <- data.frame(Intercept = mod.perm@beta[1], coal_top3_deg = mod.perm@beta[2], coal_top3_poly = mod.perm@beta[3], offset = mod.perm@beta[4])
  
  if(i <= 100){
    predictions.list[[i]] <- data.frame(rc = predict(mod.perm, type = 'response'), 
                                        coal_top3_deg = ranks.perm$coal_top3_deg,
                                        obs_counts = ranks.perm$obs_counts)
  }
}
predictions <- do.call(rbind, predictions.list)
betas <- do.call(rbind, betas.list)

##p-values from permutation
## beta coalitions - 0.001
sum(betas[1,]$coal_top3_deg <= betas$coal_top3_deg)/length(betas$coal_top3_deg)

## beta coalitions squared - 0.001
sum(betas[1,]$coal_top3_poly <= betas$coal_top3_poly)/length(betas$coal_top3_poly)

####Outlier analysis####
#removing largest rank changes doesn't change the effect
ranks.perm.outlier <- na.omit(dplyr::select(filter(ranks,abs(RankDiffAbs) <=10), ID, Clan, Year, RankDiffAbs, coal_top3_deg, obs_counts))
top3.poly.outlier <- ranks.perm.outlier  %>% lm(formula = RankDiffAbs ~ poly(coal_top3_deg,2) + log(offset(obs_counts)))
summary(top3.poly.outlier)

#predictions.list.outlier <- list()
betas.list.outlier <- list(data.frame(Intercept = top3.poly.outlier$coefficients[1], coal_top3_deg = top3.poly.outlier$coefficients[2], coal_top3_poly = top3.poly.outlier$coefficients[3], offset = top3.poly.outlier$coefficients[4]))

for(i in 1:999){
  ranks.perm.outlier$RankDiffAbs <- ranks.perm.outlier %>% group_by(Clan, Year) %>%
    sample_frac(replace = FALSE) %>% 
    ungroup() %>% dplyr::select(RankDiffAbs) %>% .[[1]]
  
  mod.perm <- ranks.perm.outlier %>% lm(formula = RankDiffAbs ~ poly(coal_top3_deg, 2) + log(offset(obs_counts)))
  
  betas.list.outlier[[i+1]] <- data.frame(Intercept = mod.perm$coefficients[1], coal_top3_deg = mod.perm$coefficients[2], coal_top3_poly = mod.perm$coefficients[3], offset = mod.perm$coefficients[4])
  
  # if(i <= 100){
  #   predictions.list.outlier[[i]] <- data.frame(rc = predict(mod.perm, type = 'response'), 
  #                                       coal_top3_deg = ranks.perm$coal_top3_deg,
  #                                       obs_counts = ranks.perm$obs_counts)
  # }
}
#predictions.outlier <- do.call(rbind, predictions.list.outlier)
betas.outlier <- do.call(rbind, betas.list.outlier)

##p-values from permutation -  outlier analysis
## beta coalitions - 0.011
sum(betas.outlier[1,]$coal_top3_deg <= betas.outlier$coal_top3_deg)/length(betas.outlier$coal_top3_deg)

## beta coalitions squared - 0.028
sum(betas.outlier[1,]$coal_top3_poly <= betas.outlier$coal_top3_poly)/length(betas.outlier$coal_top3_poly)



######Plots#####

##Figure 3

betas.summary = rbind(data.frame(estimate = names(betas), 
                                 ci = apply(X = betas, MARGIN = 2, FUN = quantile, probs = 0.025)),
                      data.frame(estimate = names(betas), 
                                 ci = apply(X = betas, MARGIN = 2, FUN = quantile, probs = 0.975)))

observed <- data.frame(estimate = names(betas),
                       mean = unlist(betas[1,]))

observed$estimate <- factor(observed$estimate, levels = c('coal_top3_deg', 'coal_top3_poly', 'Intercept', 
                                                          'offset'),
                            labels = c('Coalitions with top allies',
                                       'Coalitions with top allies (squared)',
                                       'Intercept',
                                       'Offset'))

betas.summary$estimate <- factor(betas.summary$estimate, levels = c('coal_top3_deg', 'coal_top3_poly', 'Intercept', 
                                                                    'offset'),
                                 labels = c('Coalitions with top allies',
                                            'Coalitions with top allies (squared)',
                                            'Intercept',
                                            'Offset'))

inset <- ggplot(data = betas.summary, aes(x = ci, y = estimate))+
  geom_line(size = 1) + 
  geom_point(data = observed, 
             aes(x = mean, y = estimate), size = 2, col = 'black',
             shape = 21) + 
  theme_bw() + 
  geom_vline(aes(xintercept = 0), col = 'firebrick', lty = 2) + 
  ylab(NULL)+
  xlab(NULL)+
  theme(text = element_text(size = 9),
                plot.margin = unit(c(0,0,0,0), 'pt'))
# theme(axis.text.y = element_text(size = 16))+
# theme(axis.text.x = element_text(size = 16))


top3_ranks <- data.frame(rc = predict(top3_poly, type = 'response'),
                         coal_top3_deg = na.omit(ranks[,c('coal_top3_deg', 'obs_counts')])$coal_top3_deg,
                         obs_counts = na.omit(ranks[,c('coal_top3_deg', 'obs_counts')])$obs_counts,
                         category = 'Model Predictions')
obs <- data.frame(rc = ranks$RankDiffAbs,
                  coal_top3_deg = ranks$coal_top3_deg,
                  obs_counts = ranks$obs_counts,
                  category = 'Observed')

predictions$category <- 'Permuted'

fig3_ranks <- rbind(predictions,
                    obs,
                    top3_ranks)

fig3_ranks$category <- factor(fig3_ranks$category, 
                              levels = c('Model Predictions',
                                         'Observed',
                                         'Permuted'))


main <- ggplot(data=fig3_ranks, aes(x = coal_top3_deg, y = rc, col = category))+
  geom_jitter(size = 0.5)+
  theme_classic() + 
  ylab('Positions moved up (+) or down (-) the hierarchy') +
  xlab('Strength of coalitionary ties with top partners') + 
  theme(legend.position = c(0.6,0.1))+
  ylim(c(-12, 20))+
  scale_color_manual(values = c('firebrick', 'dodgerblue', 'grey'),
                     guide = guide_legend(title = NULL))+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8), 
        legend.direction = 'horizontal',
        legend.key.size = unit(8, 'pt'))
  


source('0.multiplot.R')

pdf("plots/9.Fig3.pdf",
    height= 3.5,
    width = 4.5)
multiplot(inset, main, layout = matrix(c(1,1,
                                         2,2,
                                         2,2,
                                         2,2), byrow = TRUE, nrow = 4))
dev.off()


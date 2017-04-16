##################################################################
##                        Eli Strauss                           ##
##                  Rank related fitness                        ##
##                    February 5th, 2017                         ##
##################################################################
library(dplyr)
library(network)
library(sna)
library(Perc)
library(ggplot2)
source("~/Documents/Fisibase/fisibasetidy/ReadTidyData.R")


############Make networks for each year and each clan#############
for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    ###Dyadic aggression network
    aggNet <- network.initialize(n = length(yearRanks$ID), directed = T)
    network.vertex.names(aggNet) <- yearRanks$ID
    aggNet[as.matrix(filter(aggsFull, Clan == clan, Year == year,
                          Agg %in% yearRanks$ID, Recip %in% yearRanks$ID)[,c('Agg', 'Recip')])] <- 1
    aggNet[,]
    assign(paste('aggNet', clan, year, sep = '_'), aggNet)
    ####Coalition network
    cnet <- network.initialize(n = length(yearRanks$ID), directed = F)
    network.vertex.names(cnet) <- yearRanks$ID
    coalsTemp <- filter(aggsFull, Year == year, Clan == clan, Group != '')
    coalPartners <- cbind(expand.grid(yearRanks$ID, yearRanks$ID), rep(0))
    names(coalPartners) <- c('Focal', 'Alter', 'NumCoals')
    for(row in 1:length(coalsTemp[,1])){
      focal <- coalsTemp[row,'Agg']
      for(alter in strsplit(coalsTemp[row,'Group'], ',')[[1]][-1]){
        if(alter %in% yearRanks$ID & alter != focal & focal %in% yearRanks$ID){
          coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,'NumCoals'] <- coalPartners[coalPartners$Focal == focal & coalPartners$Alter == alter,'NumCoals'] + 1
        }
      }
    }
    coalPartnersPaired <- coalPartners[coalPartners$NumCoals != 0,]
    cnet[as.matrix(coalPartnersPaired[,1:2])] <- coalPartnersPaired[,3]
    numCoals <- as.matrix(cnet)
    numCoals[as.matrix(coalPartnersPaired[,1:2])] <- coalPartnersPaired[,3]
    set.edge.value(cnet, 'numCoals', numCoals)
    assign(paste('coalNet', clan, year, sep = '_'), cnet)
  
    ####AI Network
    obsnet <- matrix(nrow = length(yearRanks$ID), ncol = length(yearRanks$ID), dimnames = list(yearRanks$ID, yearRanks$ID), data = 0)
    ainet <- obsnet
    sessTemp <- filter(sessions, format(Date, '%Y') == year)
    for(id in yearRanks$ID){
      sess <- grep(paste(',',id,',', sep = ''), sessTemp$Hyenas)
      for(s in sess){
        part <- strsplit(sessTemp[s,'Hyenas'], ',')[[1]][-1]
        for(p in part){
          if(p %in% yearRanks$ID) {obsnet[id,p] <- obsnet[id,p]+1}
        }
      }
    }
    for(a in yearRanks$ID){
      for(b in yearRanks$ID){
        ainet[a,b] <- obsnet[a,b]/(obsnet[a,a] + obsnet[b,b] - obsnet[a,b])
      }
    }
    assign(paste('ainet', clan, year, sep = '_'), ainet)
    
    ######Greeting network
    grtnet <- matrix(nrow = length(yearRanks$ID), ncol = length(yearRanks$ID), dimnames = list(yearRanks$ID, yearRanks$ID), data = 0)
    grtsTemp <- filter(greetings, format(Date, "%Y") == year, ID1 %in% yearRanks$ID, ID2 %in% yearRanks$ID)
    if(length(grtsTemp[,1])){
      for(row in 1:length(grtsTemp)){
      grtnet[grtsTemp[row,'ID1'], grtsTemp[row,'ID2']] <- grtnet[grtsTemp[row,'ID1'], grtsTemp[row,'ID2']] + 1
      grtnet[grtsTemp[row,'ID2'], grtsTemp[row,'ID1']] <- grtnet[grtsTemp[row,'ID2'], grtsTemp[row,'ID1']] + 1
      }
    }
    
    assign(paste('grtnet', clan, year, sep = '_'), network(grtnet, ignore.eval = FALSE, names.eval = 'numGrt'))
  }
}
##################################################################



#########################Binarize ainets###########################
for(ainame in ls()[grep('^ainet_.*', ls())]){
  cutoff <- .05
  ainet <- get(ainame)
  diag(ainet) <- 0
  bin <- (ainet >= cutoff)+0
  ainameb <- paste('bin',ainame, sep= '_')
  assign(ainameb, network(bin, directed = FALSE))
}
##################################################################

#########################Full weighted ainets#####################
for(ainame in ls()[grep('^ainet_.*', ls())]){
  ainet <- get(ainame)
  diag(ainet) <- 0
  weights <- ainet
  ainet <- network(ainet, directed = FALSE)
  set.edge.value(ainet, 'AIw', as.matrix(weights))
  ainamef <- paste('full',ainame, sep= '_')
  assign(ainamef, ainet)
}
##################################################################


#####################compile centrality data######################
netDat <- ranks
netDat$Move <- ifelse(netDat$ID == netDat$IDold, 'None','Move')

for(row in rownames(netDat[netDat$Move == 'Move',])){
  if(netDat[row,'Rank'] > filter(netDat, IDold == netDat[row,'ID'], Year == netDat[row, 'Year'])$Rank){
    netDat[row,'Move'] <- 'Down'
  }else if(netDat[row,'Rank'] < filter(netDat, IDold == netDat[row,'ID'], Year == netDat[row, 'Year'])$Rank){
    netDat[row,'Move'] <- 'Up'
  }else{netDat[row,'Move'] <- NA}
}

for(clan in unique(ranks$Clan)){
  clanRanks <- filter(ranks, Clan == clan)
  for(year in unique(clanRanks$Year)){
    yearRanks <- filter(clanRanks, Year == year)
    
    ###Centrality
    ainet <- get(paste('full_ainet', clan, year, sep = '_'))
    coalNet <- get(paste('coalNet', clan, year, sep = '_'))
    aggNet <- get(paste('aggNet', clan, year, sep = '_'))
    grtNet <- get(paste('grtnet', clan, year, sep = '_'))
    
    ###Strongest allies
    netDat[netDat$Year == year & netDat$Clan == clan, 'top3coals'] <- mapply(1:length(yearRanks[,1]),
                                                                             FUN = function(x) mean(sort(as.sociomatrix(coalNet, attrname = 'numCoals')[,x],decreasing = T)[1:3]))
    netDat[netDat$Year == year & netDat$Clan == clan, 'top3ai'] <- mapply(1:length(yearRanks[,1]),
                                                                             FUN = function(x) mean(sort(as.sociomatrix(ainet, attrname = 'AIw')[,x],decreasing = T)[1:3]))
    netDat[netDat$Year == year & netDat$Clan == clan, 'top3ranks'] <- mapply(1:length(yearRanks[,1]),
                                                                          FUN = function(x) mean(filter(yearRanks, ID %in% names(sort(as.sociomatrix(ainet, attrname = 'AIw')[,x],decreasing = T)[1:3]))$Rank))
    
    
        
    netDat[netDat$Year == year & netDat$Clan == clan, 'aiDeg'] <- degree(as.sociomatrix(ainet, attrname = 'AIw'), gmode = 'graph', rescale = F, ignore.eval = F)
    netDat[netDat$Year == year & netDat$Clan == clan, 'coalDeg'] <- degree(as.sociomatrix(coalNet, attrname = 'numCoals'), gmode = 'graph', rescale = F, ignore.eval = F)
    netDat[netDat$Year == year & netDat$Clan == clan, 'aggIdeg'] <- degree(as.sociomatrix(aggNet), cmode = 'indegree', gmode = 'digraph', rescale = F, ignore.eval = T)
    netDat[netDat$Year == year & netDat$Clan == clan, 'aggOdeg'] <- degree(as.sociomatrix(aggNet), cmode = 'outdegree', gmode = 'digraph', rescale = F, ignore.eval = T)
    netDat[netDat$Year == year & netDat$Clan == clan, 'grtDeg'] <- degree(as.sociomatrix(grtNet), gmode = 'graph', rescale = F, ignore.eval = F)
  }
}
##################################################################



#####################Compare within individuals###################
moveIndComp <- filter(netDat, Move != 'None')
names(moveIndComp) <- c('Rank', 'Year', 'ID', 'IDold', 'Clan', 'stan.rank', 'Direction', 'top3coalsUnstable', 'top3aiUnstable', 'aiDegUnstable', 'coalDegUnstable', 'aggIdegUnstable', 'aggOdegUnstable', 'grtDegUnstable', 'top3ranksUnstable')
moveIndComp$aiDegStable <- NA
moveIndComp$coalDegStable <- NA
moveIndComp$aggIdegStable <- NA
moveIndComp$aggOdegStable <- NA
moveIndComp$top3coalsStable <- NA
moveIndComp$top3aiStable <- NA
moveIndComp$grtDegStable <- NA
moveIndComp$top3ranksStable <- NA
for(row in 1:length(moveIndComp[,1])){
  id <- moveIndComp[row,'ID']
  nextYear <- moveIndComp[row,'Year']+1
  nextDat <- filter(netDat, ID == id, Year == nextYear, Move == 'None')
  if(length(nextDat[,1])){
    moveIndComp[row,'aiDegStable'] <- nextDat$aiDeg
    moveIndComp[row,'coalDegStable'] <- nextDat$coalDeg
    moveIndComp[row,'aggIdegStable'] <- nextDat$aggIdeg
    moveIndComp[row,'aggOdegStable'] <- nextDat$aggOdeg
    moveIndComp[row,'grtDegStable'] <- nextDat$grtDeg
    moveIndComp[row,'top3aiStable'] <- nextDat$top3ai
    moveIndComp[row,'top3coalsStable'] <- nextDat$top3coals
    moveIndComp[row,'grtDegStable'] <- nextDat$grtDeg
    moveIndComp[row,'top3ranksStable'] <- nextDat$top3ranks
  }
}
par(mfrow = c(2,2))
attach(moveIndComp[complete.cases(moveIndComp),])
     boxplot(aiDegUnstable-aiDegStable,main = 'AI', ylab = 'Difference between Stable and Unstable Period')
     boxplot(coalDegUnstable-coalDegStable, main = 'Coalition Partners', ylab = 'Difference between Stable and Unstable Period')
     boxplot(aggIdegUnstable-aggIdegStable, main = 'Aggression Received', ylab = 'Difference between Stable and Unstable Period')
     boxplot(aggOdegUnstable-aggOdegStable, main = 'Aggression Emitted', ylab = 'Difference between Stable and Unstable Period')
detach()

par(mfrow = c(1,2))     
attach(moveIndComp[complete.cases(moveIndComp),])
     boxplot(top3aiUnstable-top3aiStable, main = 'Top 3 partners AI', ylab = 'Difference between Stable and Unstable Period')
     boxplot(top3coalsUnstable-top3coalsStable, main = 'Top 3 partners Coals', ylab = 'Difference between Stable and Unstable Period')
detach()
     
par(mfrow = c(1,2))     
attach(filter(moveIndComp, Direction == 'Up')[complete.cases(filter(moveIndComp, Direction == 'Up')),])
     boxplot(top3aiUnstable-top3aiStable, main = 'Top 3 partners AI', ylab = 'Difference between Stable and Unstable Period')
     boxplot(top3coalsUnstable-top3coalsStable, main = 'Top 3 partners Coals', ylab = 'Difference between Stable and Unstable Period')
detach()
     
par(mfrow = c(1,2))     
attach(filter(moveIndComp, Direction == 'Down')[complete.cases(filter(moveIndComp, Direction == 'Up')),])
     boxplot(top3aiUnstable-top3aiStable, main = 'Top 3 partners AI', ylab = 'Difference between Stable and Unstable Period')
     boxplot(top3coalsUnstable-top3coalsStable, main = 'Top 3 partners Coals', ylab = 'Difference between Stable and Unstable Period')
detach()


##################################################################




########################Make some plots###########################
highRankMove <- filter(moveIndComp, Rank <= 5)
##'remove Juno Loki because they are sibs. need to fix this fully later
highRankMove <- highRankMove[c(-4,-5),]
with(highRankMove, plot(top3aiStable ~ top3aiUnstable, col = as.factor(Direction), lwd = 3))
abline(a = 0, b = 1)
with(highRankMove, plot(top3coalsStable ~ top3coalsUnstable, col = as.factor(Direction), lwd = 3))
abline(a = 0, b = 1)
with(highRankMove, plot(aiDegStable ~ aiDegUnstable, col = as.factor(Direction), lwd = 3))
abline(a = 0, b = 1)
with(highRankMove, plot(aiDegStable ~ aiDegUnstable, col = as.factor(Direction), lwd = 3))
abline(a = 0, b = 1)


lowRankMove <- filter(moveIndComp, Rank > 5)
with(lowRankMove, plot(top3aiStable ~ top3aiUnstable, col = as.factor(Direction), lwd = 3))
abline(a = 0, b = 1)

moveIndComp$RankBin <- cut(moveIndComp$Rank, c(0,5,10,15,20,25,30,35,40,45,50))

ggplot(data = moveIndComp, aes((y = top3aiUnstable - top3aiStable), x = Rank, col  = Direction)) +
  geom_point() + 
  geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = moveIndComp, aes((y = top3coalsUnstable - top3coalsStable), x = Rank, col  = Direction)) +
  geom_point() + 
  geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = moveIndComp, aes((y = aiDegUnstable - aiDegStable), x = Rank, col  = Direction)) +
  geom_point() + 
  geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = moveIndComp, aes((y = coalDegUnstable - coalDegStable), x = RankBin, fill = Direction)) +
  geom_boxplot(notch = F) + 
  geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = moveIndComp, aes(y = top3ranksUnstable - top3ranksStable, x = Rank, col = as.factor(Direction))) +
  geom_point() + 
  #geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = filter(moveIndComp, Clan == 'talek'), aes((y = grtDegUnstable - grtDegStable), x = Rank, col  = Direction)) +
  geom_point() + 
  geom_smooth(method = lm)+
  ylab('Difference betwen Rank change year and subsequent year')

ggplot(data = moveIndComp)+
  geom_boxplot(aes(y = coalDegUnstable, x = RankBin, fill = 'red'), alpha = .5) + 
  geom_boxplot(aes(y = coalDegStable, x = RankBin, fill = 'blue'), alpha = .5) +
  scale_fill_identity()

ggplot(data = moveIndComp)+
  geom_boxplot(aes(y = aiDegUnstable, x = RankBin, fill = 'red'), alpha = .5) + 
  geom_boxplot(aes(y = aiDegStable, x = RankBin, fill = 'blue'), alpha = .5, position = 'dodge') +
  scale_fill_identity()

ggplot(data = moveIndComp)+
  geom_boxplot(aes(y = coalDegUnstable, x = RankBin, fill = 'red'), alpha = .5) + 
  geom_boxplot(aes(y = coalDegStable, x = RankBin, fill = 'blue'), alpha = .5) +
  scale_fill_identity()

melt(moveIndComp, id.vars = c('Rank', 'Year', 'ID', 'IDold', 'Clan'))


####################################Modularity#############################
library(bipartite)
library(igraph)
?computeModules
plotModuleWeb(computeModules(aggNet_talek_2012[,]))

g <- graph_from_adjacency_matrix(as.sociomatrix(full_ainet_talek_1995, attrname = 'AIw'), mode = 'undirected', weighted = T)
#g <- graph_from_adjacency_matrix(as.sociomatrix(coalNet_talek_2007, attrname = 'numCoals'), mode = 'undirected', weighted = T)


remove.isolates <- function(g){
  g <- igraph::delete.vertices(g, which(igraph::degree(g) == 0))
  return(g)
}


gplot <- remove.isolates(g)
c <- cluster_louvain(gplot)
layweight <- ifelse(crossing(c, gplot), 2, 1)
layout <- layout_with_kk(gplot, weights=layweight)
plot.igraph(gplot, edge.width = E(gplot)$weight, 
            vertex.color = membership(c),
            layout = layout,
            edge.color = c("darkgrey","tomato2")[crossing(c, gplot) + 1])


plot.igraph(gplot, edge.width = E(gplot)$weight*20, vertex.color = colors(length(igraph::degree(g))), vertex.label = "")
colors <- colorRampPalette(colors = c('black', 'grey'))
colors(length(igraph::degree(g)))


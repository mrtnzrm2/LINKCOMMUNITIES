library(ggtree)
library(ggplot2)
library(MASS)
library(DescTools)
library(pracma)
library(rayshader)
library(RColorBrewer)
library(viridis)
library(aricode)
library(dendextend)
library(scales)
library(mgcv)
library(network)
library(gridExtra)
library(igraph)
library(ape)
library(tidytree)
library(hrbrthemes)
library(RcppHungarian)
library(fields)
library(wesanderson)
library(fitdistrplus)
library(latex2exp)
library(ggraph)
library(akima)
library(cowplot)
library(rstan)
library(Cairo)
library(GraphAlignment)
library(reticulate)
library(stringr)
library(dplyr)
library(msir)

apply_munkres <- function(x, y, kmax){
  
  COST <- matrix(0, nrow = kmax, ncol = kmax)
  
  for (i in 1:kmax){
    for (j in 1:kmax){
      nodes.x <- which(x == i)
      nodes.y <- which(y == j)
      
      COST[i,j] <- intersect(nodes.y, nodes.x) %>%
        length()
    }
  }
  
  COST <- 1 / (COST + 1)
  
  perm <- HungarianSolver(COST)$pairs[,2]
  
  perm.x <- rep(0, length(x))
  for (i in 1:kmax)
    perm.x[x == i] <- perm[i]
  
  return(perm.x)
}

table.RMAE <- function(pred.x){
  
  #wcut < w < 3
  m03 <- mean(pred.x$rmae.mean[which(pred.x$cat == '0<w<1' |
                                pred.x$cat == '1<w<2' |
                                pred.x$cat == '2<w<3' &
                                pred.x$w > 0.9)])
  #wcut < w < 5
  m05 <- mean(pred.x$rmae.mean[which(pred.x$cat == '0<w<1' |
                                pred.x$cat == '1<w<2' |
                                pred.x$cat == '2<w<3' |
                                pred.x$cat == '3<w<4' |
                                pred.x$cat == '4<w<5' &
                                pred.x$w > 0.9)])
  
  #w > 3
  m3 <- mean(pred.x$rmae.mean[which(!(pred.x$cat == '0<w<1' |
                                  pred.x$cat == '1<w<2' |
                                  pred.x$cat == '2<w<3') &
                                  pred.x$w > 0.9)])
  
  #w > 5
  m5 <- mean(pred.x$rmae.mean[which(!(pred.x$cat == '0<w<1' |
                                  pred.x$cat == '1<w<2' |
                                  pred.x$cat == '2<w<3' |
                                  pred.x$cat == '3<w<4' |
                                  pred.x$cat == '4<w<5') &
                                  pred.x$w > 0.9)])
  #w > wcut
  mucut <- mean(pred.x$rmae.mean[which(pred.x$w > 0.9)])
  
  #w <= wcut
  mdcut <- mean(pred.x$rmae.mean[which(pred.x$w <= 0.9)])
  
  # all w
  ma <- mean(pred.x$rmae.mean)
  
  summary.RMAE <- data.frame(intervals=c('Weak (0.9<w<3)',
                                      'Weak-Medium (0.9<w<5)',
                                      'Medium-Strong (w>3)',
                                      'Strong (w>5)',
                                      'All links (w>0.9)',
                                      'Non-links (w<=0.9)',
                                      'Both links and non-links'),
                             RMAE=c(m03,
                                    m05,
                                    m3,
                                    m5,
                                    mucut,
                                    mdcut,
                                    ma),
                             ranges=c('0.9<w<3',
                                         '0.9<w<5',
                                         'w>3',
                                         'w>5',
                                         'w>0.9',
                                         'w<0.9',
                                         'Both'))
  
  return(summary.RMAE)
  
}

run.pred <- function(demo, k=F, null=F, X=23){

  if (strcmp(demo,'monkey')){
    model <- sprintf('mean_fit_%ix%i_flnMonkey_dis_zzz_reduced', 91, 40)
  } else{
    model <- sprintf('mean_fit_%ix%i_flnMouse_dis_zzz_reduced', 47, 19)
  }

  foldername <- model
  if (null)
    foldername <- sprintf('%s_null', foldername)
  
  l.files <- list.files(sprintf('%s/RDS', path.plot), pattern = model)
  if (null){
    l.files <- l.files[grepl('null', l.files)] 
  } else{
    l.files <- l.files[!grepl('null', l.files)] 
  }
  if (k && strcmp(demo,'mouse')){
    l.files <- l.files[grepl('k', l.files)] 
  } else if (!k && strcmp(demo,'mouse')){
    l.files <- l.files[!grepl('k', l.files)] 
  }
  
  #### Load network ####
  netx <- load.net(demo, tag = '')
  net <- netx$net
  nodes <- netx$nodes
  leaves <- netx$leaves
  labels <- netx$labels
  regionsCSV <- netx$regions
  dis.coords <- netx$dis.coords
  dis <- netx$dis
  dis <- dis[which(dis$target <= nodes & dis$source != dis$target),]
  
  regionsCSV$REGION <- regionsCSV$REGION %>%
    tolower()
  
  if (strcmp(demo,'monkey'))
    regionsCSV$REGION[which(regionsCSV$REGION %in% 'insulate')] <- 'insular'
  
  max.sim <- ceil(-log10(min(net$weight)))
  net$w <- log10(net$weight) + max.sim
  x.label <- sprintf('log(FLN)+%i', max.sim)
  
  is.ncd <- get.known.noconnections(net)
  
  w <- rep(NA, nrow(dis))
  w[!is.ncd] <- net$w
  
  dis.var <- rep(NA, nrow(dis))
  dis.var[!is.ncd] <- dis$weight[!is.ncd]
  
  cat.w <- rep(NA, nrow(net))
  for (e in 1:nrow(net)){
    floor.w <- floor(net$w[e])
    ceil.w <- ceil(net$w[e])
    cat.w[e] <- sprintf('%i<w<%i', floor.w, ceil.w)
  }
  
  w.cat <- rep(NA, nrow(dis))
  w.cat[!is.ncd] <- cat.w
  
  dis$weight <- dis$weight/max(dis$weight)
  
  nfl <- length(l.files)
  PRED <- data.frame()
  COUNT <- data.frame()
  
  for (fl in l.files){
    fl.tr <- str_split(fl, 'tr') %>% unlist()
    fl.tr <- str_split(fl.tr[2], '_v') %>% unlist()
    tr <- fl.tr[1]
    fl.tr <- str_split(fl.tr[2], '_reduced') %>% unlist()
    v <- fl.tr[1]
    
    mdl <- readRDS(sprintf('%s/RDS/%s', path.plot, fl))
    
    asym <- mdl$asym %>% from_adjacency_to_dataframe()
    asym <- asym[which(asym$source != asym$target),]
    sig <- mdl$sigma
    
    eta <- dis$weight + asym$weight
    
    if (!null && strcmp(model, 'monkey')){
      l <- mdl$l %>% from_adjacency_to_dataframe()
      l <- l[which(l$source != l$target),]
      eta <- eta + l$weight
    }
    
    eta <- eta/sig
    
    eta <- data.frame(source=asym$source,
                      target=asym$target,
                      weight=eta)
    
    hidden.nodes <- mdl$test.index
    
    if (tr == X)
      COUNT <- rbind(COUNT, data.frame(AREA=hidden.nodes))
    
    W <- net$w[which(!(net$target %in% hidden.nodes))]
    ETA <- eta$weight[which(!(eta$target %in% hidden.nodes) & !is.ncd)]
    
    lm.mode <- loess(W ~ ETA)
    
    pred.w <- predict(lm.mode, newdata = data.frame(ETA=eta$weight))
    rmae <- rep(NA, nrow(eta))
    rmae[!is.ncd] <-abs(net$w - pred.w[!is.ncd])/net$w
    
    train <- rep(F, nrow(eta))
    train[which(!(eta$target %in% hidden.nodes))] <- T
    
    PRED <- rbind(PRED, data.frame(source=eta$source,
                                   target=eta$target,
                                   pred=pred.w,
                                   tr=tr,
                                   v=v,
                                   train=train,
                                   rmae=rmae,
                                   w=w, 
                                   dis=dis.var,
                                   cat=w.cat))
    
    
  }
  
  x <- list(pred = PRED,
            net = net,
            nodes = nodes,
            labels = labels,
            regions = regionsCSV,
            foldername = foldername,
            count=COUNT)
  
  return(x)
}

get.known.noconnections <- function(net){
  net <- from_dataframe_to_adjacency(net)
  net[net == 0] <- NA
  net <- net %>% from_adjacency_to_dataframe()
  net <- net[which(net$source != net$target),]
  
  return(is.na(net$weight))
}

fln.H.community <- function(net, e, cascade){
  
  labels <- cascade[e,]
  labels[is.na(labels)] <- 0
  net$id <- ''
  net$comm <- T
  for (e in 1:nrow(net)){
    if (labels[net$source[e]] == labels[net$target[e]]){
      net$id[e] <- as.character(labels[net$source[e]])
    } else if (labels[net$source[e]]  > labels[net$target[e]]) {
      net$id[e] <- sprintf('%i-%i', labels[net$source[e]], labels[net$target[e]])
      net$comm[e] <- F
    } else if (labels[net$source[e]]  < labels[net$target[e]]) {
      net$id[e] <- sprintf('%i-%i', labels[net$target[e]], labels[net$source[e]])
      net$comm[e] <- F
      
    }
  }
  
  p <- ggplot(net, aes(weight))+
    facet_grid(comm~id)+
    geom_histogram(bins=10)+
    scale_x_continuous(trans = 'log10')+
    theme_bw()+
    theme(axis.text.x= element_text(angle=90))
  
  print(p)
}

measure.barea <- function(net.merde, nodes){
  
  cascade <- matrix(0, ncol = nodes, nrow = nodes)
  cascade[nodes,] <- as.factor(cutree(net.merde, k=nodes))
  
  height <- matrix(1, ncol = nodes, nrow = nodes)
  
  for (e in (nodes-1):1){
    cascade[e,] <- cutree(net.merde, k=e)
    old.cls <- unique(cascade[e+1,])
    new.cls <- unique(cascade[e,])
    
    for (e.o in old.cls){
      for (e.i in new.cls){
        
        new.nodes <- which(cascade[e,] == e.i)
        old.nodes <- which(cascade[e+1,] == e.o)
        
        rep.n <- length(intersect(new.nodes, old.nodes))
        
        if ( rep.n >= 1 && length(new.nodes) > length(old.nodes)){
          height[e,new.nodes] <- height[e+1,new.nodes] + 1
        } else {
          for (n in new.nodes){
            height[e,n] <- max(height[,n])
          }
        }
      }
    }
  }
  
  cascade <- matrix(0, ncol = nodes, nrow = nodes)
  cascade[nodes,] <- 1:nodes
  
  for (e in (nodes-1):1){
    
    cascade[e,] <- cascade[e+1,]
    
    new.coms <- cutree(net.merde, k=e)
    old.coms <- cutree(net.merde, k=e+1)
    
    old.cls <- unique(old.coms)
    new.cls <- unique(new.coms)
    
    for (e.i in new.cls){
      new.nodes <- which(new.coms == e.i)
      n <- 1
      cfirst <- c()
      csecond <- c()
      dog <- unlist(lapply(old.cls, function(x) 
        ifelse(length(intersect(new.nodes, which(old.coms == x))) >= 1,
               T, F)))
      if (sum(dog) == 2){
        cfirst <-  which(old.coms == old.cls[which(dog)[1]])
        csecond <-  which(old.coms == old.cls[which(dog)[2]])
      }
      if (length(cfirst) > 0 && length(csecond) > 0){
        for (n in nodes:e){
          n1 <- sum(height[n, cfirst] > 1)
          n2 <- sum(height[n, csecond] > 1)
          if (n1 > n2 && n1 > 0){
            the.nodes <- cfirst[which(height[n,cfirst] != 1)]
            break
          } else if (n2 > n1 && n2 > 0){
            the.nodes <- csecond[which(height[n,csecond] != 1)]
            break
          } else if (n2 == n1 && n2 > 0){
            the.nodes <- c(cfirst[which(height[n,cfirst] != 1)],
                           csecond[which(height[n,csecond] != 1)])
            break
          }
        }
        cascade[e, new.nodes] <- cascade[n, max(the.nodes)]
      }
    }
  }
  
  # height <- height[,order.trd]
  height <- from_adjacency_to_dataframe(height)
  
  # cascade <- cascade[,order.trd]
  df.cascade <- from_adjacency_to_dataframe(cascade)
  df.cascade$weight <- as.factor(df.cascade$weight)
  df.cascade$weight[which(height$weight == 1)] <- NA 
  df.cascade <- df.cascade[!is.na(df.cascade$weight),]
  
  colors <- df.cascade$weight %>%
    unique()
  
  colors <- lapply(colors, function(x) sum(df.cascade$weight == x)) %>%
    unlist() %>%
    sort(decreasing = T)
  
  area <- list()
  
  area$max <- max(colors, na.rm = T)
  area$mean <- mean(colors)
  area$max.2 <- sum(colors[1:2])
  area$sd <- sd(colors)
  
  return(area)
  
}

gg.plotNetwork_k <- function(net, labels, comms, coords, k, kmax,
                           path='', subfolder='', plt=F, filename=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s', path, subfolder), showWarnings = FALSE)
  }
  
  Data <- data.frame()
  Labels <- labels
  
  for (e in 1:nrow(comms)){
    
    new.idx <- order(comms[e,])
    new.comms <- comms[e, new.idx]
    
    mk.coord <- as.matrix(coords)
    mk.coord <- mk.coord[new.idx,]
    Data <- rbind(Data, data.frame(x=mk.coord[,1],
                                   y=mk.coord[,2],
                                   COMMS=as.factor(new.comms),
                                   LABELS=Labels[new.idx],
                                   k=e+1))
  }
 
  Data$k <- Data$k %>% as.factor()
  
  pp <- ggplot(Data, aes(x,y, fill=COMMS))+
    facet_wrap(~k)+
    geom_point(shape=21, size=10, stroke = 0.8)+
    geom_text(label=Data$LABELS, color='white', size=2.5, fontface='bold')+
    theme_void()
  
  if (plt){
    cairo_ps(filename = sprintf("%s/%s/%s.eps", path, subfolder, filename),
             width =22, height = 17,
             fallback_resolution = 200)
    print(pp)
    dev.off()
  } else{
    print(pp)
  }
  
}

gg.plotNetwork <- function(net, labels, comms, coords, k, kmax,
                        path='', subfolder='', plt=F, filename=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s', path, subfolder), showWarnings = FALSE)
  }
  
  new.idx <- order(comms)
  comms <- comms[new.idx]
  
  if (5 %in% comms){
    comms[comms == 5] <- 6
  }
  
  labels <- labels[new.idx]
  net.mat <- from_dataframe_to_adjacency(net)
  net.mat <- net.mat[new.idx, new.idx]
  net <- from_adjacency_to_dataframe(net.mat)
  net <- net[which(net$weight != 0),]
  
  mk.coord <- as.matrix(coords)
  mk.coord <- mk.coord[new.idx,]
  
  # region.df <- regions[new.idx,]
  color.graph <- hue_pal()(kmax)
  
  data <- data.frame(x=mk.coord[,1],
                     y=mk.coord[,2],
                     # COLOR=region.df$COLOR,
                     # REGION=as.factor(region.df$REGION),
                     COMMS=as.factor(comms),
                     COLCOM=color.graph[comms])
  
  # hulls <- data %>%
  #   group_by(REGION) %>%
  #   do(.[chull(.[1:2]),])
  
  # filt.data <- data[!duplicated(data$COMMS),]
  
  # color_regions <- region.df[, c('REGION','COLOR')]
  # color_regions <- color_regions[!duplicated(color_regions$REGION),]
  # color_regions <- rbind(color_regions, data.frame(REGION=filt.data$COMMS,
  #                                                  COLOR=filt.data$COLCOM))
  # 
  # color_regions <- color_regions[order(color_regions$REGION),]
  # sort.regions <- regions[order(regions$REGION),]
  
  # pp <- ggplot(data, aes(x,y, fill=REGION))+
  #   geom_shape(data=hulls, alpha=0.4, expand = unit(5, 'mm'), radius = unit(5, 'mm'))+
  #   geom_point(aes(fill=COMMS), shape=21, size=7, stroke = 0.8)+
  #   scale_fill_manual(values = color_regions$COLOR)+
  #   geom_text(label=labels, color='white', size=3, fontface='bold')+
  #   theme_void()
  
  pp <- ggplot(data, aes(x,y, fill=COMMS))+
    geom_point(shape=21, size=10, stroke = 0.8)+
    geom_text(label=labels, color='white', size=2.5, fontface='bold')+
    theme_void()
  
  if (plt){
    cairo_ps(filename = sprintf("%s/%s/%s.eps", path, subfolder, filename),
             width =9, height = 9,
             fallback_resolution = 200)
    print(pp)
    dev.off()
  } else{
    print(pp)
  }
  
}

plotWSBM <- function(net, labels, area.labels, size, k, path='', subfolder='',
                          plt=F, filename=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s', path, subfolder), showWarnings = FALSE)
  }
  
  new.idx <- order(labels)
  labels <- labels[new.idx]
  area.labels <- area.labels[new.idx]
  net.mat <- from_dataframe_to_adjacency(net)
  net.mat <- net.mat[new.idx, new.idx]
  net <- from_adjacency_to_dataframe(net.mat)
  net <- net[which(net$weight != 0),]
  
  net$weight <- log10(net$weight) + 7
  
  tabs <- rep(0, k)
  tabs[1] <- size[1]
  for (e in 2:k){
    tabs[e] <- tabs[e-1] + size[e]
  }
  
  tabs <- tabs + 0.5
  
  net$source.label <- area.labels[net$source]
  net$target.label <- area.labels[net$target]
  
  net$source.label <- factor(net$source.label, levels = rev(area.labels))
  net$target.label <- factor(net$target.label, levels = area.labels)
  
  pp <- ggplot(net, aes(target, source))+
    geom_tile(aes(fill=weight))+
    scale_fill_viridis(direction=-1, option = 'A')+
    scale_y_continuous(trans = 'reverse', breaks = 1:max(net$target),
                       labels = area.labels)+
    scale_x_continuous(breaks = 1:max(net$target), labels = area.labels)+
    geom_line(data=data.frame(x=c(0,tabs[1]), y=c(0,0)), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0,0), y=c(0,tabs[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0,tabs[1]), y=c(0,0)), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0,0), y=c(0,tabs[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(tabs[1],tabs[1]), y=c(0,tabs[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0,tabs[1]), y=c(tabs[1],tabs[1])), aes(x,y), color='green', size=1)
  
  xp <- tabs[1]
  yp <- xp
  
  for (e in 2:(length(tabs))){
    pp  <- pp + geom_line(data=data.frame(x=c(xp,tabs[e]), y=c(yp,yp)), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,xp), y=c(yp,tabs[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(tabs[e],tabs[e]), y=c(yp,tabs[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,tabs[e]), y=c(tabs[e],tabs[e])), aes(x,y), color='green', size=1)
    
    xp <- tabs[e]
    yp <- xp
  }
  
  leg.pp <- get_legend(pp)
  
  pp <- pp + guides(fill=guide_legend(title = 'w'))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
  
  # region.df <- regions[new.idx,]
  # region.df$position <- 1
  # region.df$NODE <- 1:nodes
  # color.regions <- regions[order(regions$REGION),]
  # 
  # gg_leg_y <- ggplot(region.df, aes(y=NODE, x=position, fill=REGION))+
  #   geom_tile()+
  #   scale_x_continuous(expand = c(0, 0))+
  #   scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
  #   scale_fill_manual(values = unique(color.regions$COLOR))+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   theme_void()
  # 
  # leg.gg_y <- get_legend(gg_leg_y)
  # 
  # gg_leg_x <- ggplot(region.df, aes(x=NODE, y=position, fill=REGION))+
  #   geom_tile()+
  #   scale_y_continuous(expand = c(0, 0))+
  #   scale_x_continuous(expand = c(0, 0))+
  #   scale_fill_manual(values = unique(color.regions$COLOR))+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   theme_void()
  # 
  # leg.gg_x <- get_legend(gg_leg_x)
  
  # legend <- plot_grid(leg.pp, leg.gg_y, nrow = 2, rel_heights = c(1,5), axis = 't')
  # legend <- leg.pp
  
  # g <- plot_grid(gg_leg_x + theme(legend.position = 'none',
  #                                 panel.grid = element_blank(),
  #                                 panel.border = element_blank()), 
  #                NULL, NULL,
  #                pp + theme(legend.position = 'none',
  #                           panel.grid = element_blank(),
  #                           panel.border = element_blank()),
  #                gg_leg_y + theme(legend.position = 'none',
  #                                 panel.grid = element_blank(),
  #                                 panel.border = element_blank()), 
  #                legend, 
  #                ncol=3, nrow = 2, rel_widths=c(15,1,2), rel_heights = c(1,15),
  #                align = 'hv' ,
  #                scale=0.98)
  if (plt){
    cairo_ps(filename = sprintf("%s/%s/%s.eps", path, subfolder, filename),
             width = 13, height = 11,
             fallback_resolution = 200)
    print(pp)
    
    dev.off()
  } else{
    print(g)
  }
  
}

with.potatoes <- function(a){
  
  element.a <- c()
  size.a <- c()
  for(e in 1:length(a)){
    if (!(a[e] %in% element.a)){
      element.a <- c(element.a, a[e])
      size.a <- c(size.a, sum(a == a[e]))
    }
  }
  
  return(list(element.a, size.a))
}

plot.HDensity <- function(net, net.merde, k.merde, k, nodes, height, labels,
                          regions, plt=F, filename='', foldername='',
                          subfolder='', path='', animal='monkey'){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.sim <- ceil(-log10(min(net$weight))) + 1
  net$weight <- log10(net$weight) + max.sim
  GH.partition <- cutree(net.merde, k=k.merde) %>%
    unname()
  GH.unique <- unique(GH.partition)
  
  net.zeros <- from_dataframe_to_adjacency(net)
  net.zeros <- from_adjacency_to_dataframe(net.zeros)
  
  net.zeros$Id.source <- 0
  net.zeros$Id.target <- 0
  for (e in 1:nrow(net.zeros)){
    net.zeros$Id.source[e] <- GH.partition[net.zeros$source[e]]
    net.zeros$Id.target[e] <- GH.partition[net.zeros$target[e]]
  }
  
  n.GH <- length(GH.unique)
  GH.density <- matrix(0, nrow = n.GH, ncol = n.GH)
  GH.weight <- matrix(0, nrow = n.GH, ncol = n.GH)
  for (i in 1:n.GH){
    for (j in 1:n.GH){
      sub.net <- net.zeros[which(net.zeros$Id.source == GH.unique[i] & 
                                 net.zeros$Id.target == GH.unique[j]
      ),
      ]
      if (i == j){
        
        n.nodes <- length(unique(c(sub.net$target,
                                   sub.net$source
                                  )
                                )
                          )
        
        sub.net <- sub.net[which(sub.net$weight != 0),]
        GH.density[i,j] <- nrow(sub.net)/(n.nodes*(n.nodes-1))
      } else{
        s.nodes <- sub.net$source %>%
          unique() %>%
          length()
        t.nodes <- sub.net$target %>%
          unique() %>%
          length()
        
        sub.net <- sub.net[which(sub.net$weight != 0),]
        GH.density[i,j] <- nrow(sub.net)/(s.nodes*t.nodes)
      }
      GH.weight[i,j] <- mean(sub.net$weight)
    }
  }
  
  NET <- matrix(0, nrow = nodes, ncol = nodes) %>%
    from_adjacency_to_dataframe()
  NET$density <- NA
  NET$intensity <- NA
  
  for (e in 1:nrow(NET)){
    Id.source <- GH.partition[NET$source[e]]
    Id.target <- GH.partition[NET$target[e]]
    src <- which(GH.unique == Id.source)
    tgt <- which(GH.unique == Id.target)
    NET$density[e] <- GH.density[src,tgt]
    NET$intensity[e] <- GH.weight[src,tgt]
  }
  
  tr <- as.phylo(net.merde)
  trd <- fortify(tr)
  trd <- subset(trd, isTip)
  trd <- trd[order(trd$angle, decreasing = T),]
  
  historia <- extract.k.partition(net, net.cluster, nrow(net), nodes, k, height)
  historia <- historia[match(trd$label, labels[1:nodes])]
  historia <- with.potholes(historia)
  
  u.hist <- historia[[1]]
  tab.t <- historia[[2]]
  tab2.t <- tab.t
  
  for (e in 2:length(tab.t)){
    tab.t[e] <- tab.t[e] + tab.t[e-1]
  }
  tab.t <- tab.t + 0.5
  tab.op <- c(0, tab.t)
  for (e in 1:length(tab2.t)){
    tab2.t[e] <- tab.t[e] - (tab.op[e+1] - tab.op[e])/2
  }
  
  tab2.t <- as.integer(tab2.t)
  
  trd.names <- trd$label
  
  region.df <- regions[which(regions$AREA %in% trd.names),]
  region.df <- region.df[match(trd.names, region.df$AREA),]
  
  A.density <- NET[,c('source', 'target', 'density')]
  colnames(A.density) <- c('source', 'target', 'weight')
  A.density <- A.density %>%
    from_dataframe_to_adjacency()
  colnames(A.density) <- labels[1:nodes]
  rownames(A.density) <- labels[1:nodes]
  A.density <- A.density[trd.names,trd.names]
  A.density <- A.density %>%
    from_adjacency_to_dataframe()
  A.density <- A.density[which(A.density$weight != 0),]
  
  A.intesity <- NET[,c('source', 'target', 'intensity')]
  colnames(A.intesity) <- c('source', 'target', 'weight')
  A.intesity <- A.intesity %>%
    from_dataframe_to_adjacency()
  colnames(A.intesity) <- labels[1:nodes]
  rownames(A.intesity) <- labels[1:nodes]
  A.intesity <- A.intesity[trd.names,trd.names]
  A.intesity <- A.intesity %>%
    from_adjacency_to_dataframe()
  A.intesity <- A.intesity[which(A.intesity$weight != 0),]
  
  region.df$position <- 1
  region.df$NODE <- 1:nodes
  
  gg_leg <- ggplot(region.df, aes(x=NODE, y=position, fill=REGION))+
    geom_tile()+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
    guides(fill=guide_legend(title = 'REGION'))+
    geom_text(label=trd.names, color='white', fontface=2, angle=90, size=3)+
    theme_void()
  
  leg.gg <- get_legend(gg_leg)
  
  pden <- ggplot(A.density, aes(target, source))+
    geom_tile(aes(fill=weight))+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(trans = 'reverse',expand = c(0, 0),
                       breaks = 1:nodes, labels = trd.names,
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    theme(axis.text.y.left = element_text(size=6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank())+
    scale_fill_viridis(direction = -1, option = 'A')
  
  pden <- pden + geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(0.5,0.5)), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,0.5), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(tab.t[1],tab.t[1]), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(tab.t[1],tab.t[1])), aes(x,y), color='green', size=1)
  
  xp <- tab.t[1]
  yp <- xp
  
  for (e in 2:(length(tab.t))){
    pden  <- pden + geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(yp,yp)), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,xp), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(tab.t[e],tab.t[e]), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(tab.t[e],tab.t[e])), aes(x,y), color='green', size=1)
    
    xp <- tab.t[e]
    yp <- xp
  }
  
  pden <- pden + guides(fill=guide_legend(title = 'density'))
  
  leg.pden <- get_legend(pden)
  
  legend <- plot_grid(leg.gg, leg.pden, nrow = 2, rel_heights = c(1,2), axis = 't')
  
  pden <- plot_grid(gg_leg + theme(legend.position = 'none',
                                   panel.grid = element_blank(),
                                   panel.border = element_blank()),
                    pden + theme(legend.position = 'none',
                                 panel.grid = element_blank(),
                                 panel.border = element_blank()), 
                    nrow = 2, 
                    rel_heights =c(1,15), 
                    align = 'v' ,
                    scale=0.98)
  
  pden <- plot_grid(pden,
                    legend,
                    ncol = 2, 
                    rel_widths=c(15,2), 
                    align = 'h' ,
                    scale=0.98)
  
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s_den_k_%i.eps", path, animal,
                                foldername, subfolder, filename, k.merde),
             width = 15, height = 10,
             fallback_resolution = 200)
    print(pden)
    dev.off()
  }
  else {
    return(pden)
  }
  
  pint <- ggplot(A.intesity, aes(target, source))+
    geom_tile(aes(fill=weight))+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(trans = 'reverse',expand = c(0, 0),
                       breaks = 1:nodes, labels = trd.names,
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    theme(axis.text.y.left = element_text(size=6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank())+
    scale_fill_viridis(direction = -1, option = 'A')
  
  pint <- pint + geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(0.5,0.5)), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,0.5), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(tab.t[1],tab.t[1]), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(tab.t[1],tab.t[1])), aes(x,y), color='green', size=1)
  
  xp <- tab.t[1]
  yp <- xp
  
  for (e in 2:(length(tab.t))){
    pint  <- pint + geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(yp,yp)), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,xp), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(tab.t[e],tab.t[e]), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(tab.t[e],tab.t[e])), aes(x,y), color='green', size=1)
    
    xp <- tab.t[e]
    yp <- xp
  }
  
  pint <- pint + guides(fill=guide_legend(title = 'intensity'))
  
  leg.pint <- get_legend(pint)
  
  legend <- plot_grid(leg.gg, leg.pint, nrow = 2, rel_heights = c(1,2), axis = 't')
  
  pint <- plot_grid(gg_leg + theme(legend.position = 'none',
                                   panel.grid = element_blank(),
                                   panel.border = element_blank()),
                    pint + theme(legend.position = 'none',
                                 panel.grid = element_blank(),
                                 panel.border = element_blank()), 
                    nrow = 2, 
                    rel_heights =c(1,15), 
                    align = 'v' ,
                    scale=0.98)
  
  pint <- plot_grid(pint,
                    legend,
                    ncol = 2, 
                    rel_widths=c(15,2), 
                    align = 'h' ,
                    scale=0.98)
  
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s_int_k_%i.eps", path, animal,
                                foldername, subfolder, filename, k.merde),
             width = 15, height = 10,
             fallback_resolution = 200)
    print(pint)
    dev.off()
  }
  else {
    return(pint)
  }
}

plotNetwork.seeds <- function(net, net.merde, net.cluster, nodes,
                              k, k.merde, labels, mk.coord, cascade,
                              plt=F,
                              filename='', subfolder='', path='', 
                              foldername='', animal='monkey'){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername,  subfolder), showWarnings = FALSE)
  }
  
  
  mk.coord <- as.matrix(mk.coord)
  
  seeds.na <- cascade[k.merde,]
  ncols <- seeds.na[!is.na(seeds.na)] %>%
    unique() %>%
    length()
  
  colors <- gg_color_hue(ncols)
  
  # regions <- regions[match(labels, regions$AREA),]
  
  cluster <- rep(-1, nodes)
  cluster[!is.na(seeds.na)] <- seeds.na[!is.na(seeds.na)]
  # data <- data.frame(x=mk.coord[,1],
  #                    y=mk.coord[,2],
  #                    cluster=as.factor(cluster),
  #                    label=labels[1:nodes],
  #                    REGION=regions$REGION[1:nodes])
  
  data <- data.frame(x=mk.coord[,1],
                     y=mk.coord[,2],
                     cluster=as.factor(cluster),
                     label=labels[1:nodes])
  
  pp <- ggplot(data, aes(x,y))+
    # geom_point(aes(color=cluster, shape=REGION), size=10)+
    geom_point(aes(color=cluster), size=10)+
    scale_color_manual(values = c('grey', colors))+
    geom_text(label=labels[1:nodes], color='black', fontface=2, size=3.5)+
    coord_fixed()+
    theme_void()
  
  if (plt){
    cairo_ps(filename =  sprintf("%s/%s/%s/%s.eps", path, foldername, subfolder, filename),
             width =7, height = 10,
             fallback_resolution = 200)
    print(pp)
    dev.off()
  } else {
    print(pp)
  }
}

le.grand.cascade <- function(net, net.merde, nodes, plt=F, path='',
                             animal='monkey', foldername='', subfolder='',
                             filename='merde'){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.sim <- ceil(-log10(min(net$weight))) + 1
  net$weight <- log10(net$weight) + max.sim
  
  cascade <- cascade.sil.vous.plaît(net, net.merde, nodes)
  
  grand.cascade <- data.frame(den=c(),
                              # in.out=c(),
                              # r.w=c(),
                              int=c(),
                              nat=c(),
                              height=c(), 
                              id=c())
  
  for (e in (nodes-1):2){
    seeds <- cascade[e,] %>%
      unique()
    
    seeds <- seeds[!is.na(seeds)]
    
    for (sds in seeds){
      nds <- which(cascade[e,] == sds)
      nd <- nds %>%
        length()
      
      v.links <- net[which(net$target %in% nds),]
      h.links <- net[which(net$source %in% nds),]
      s.links <- net[which(net$source %in% nds & net$target %in% nds),]
      
      # in.out <- nrow(s.links)/(nrow(v.links) + nrow(h.links) - 2*nrow(s.links))
      # r.w <- sum(s.links$weight)/(sum(v.links$weight) + sum(h.links$weight) + 
      #                               2*sum(s.links$weight))
      nat <- (sum(v.links$weight) - sum(s.links$weight))/(sum(h.links$weight) -
                                                            sum(s.links$weight))
      if (nd > 1)
        den <- nrow(s.links)/(nd*(nd-1))
      else
        den <- NA
      
      int <- mean(s.links$weight)
      
      grand.cascade <-  rbind(grand.cascade, data.frame(den=den,
                                                        # in.out=in.out,
                                                        # r.w=r.w,
                                                        int=int,
                                                        nat=nat,
                                                        height=e,
                                                        id=sds))
    }
    
  }
  
  grand.cascade$id <- as.factor(grand.cascade$id)
  
  GRAND <- data.frame(VALUE=c(), VAR=c(), HEIGHT=c(), IDS=c())
  VARS <- c( 'nat', 
             # 'in.out', 
             # 'r.w',
             'int',
             'den')
  
  for (var in VARS){
    
    gad <- grand.cascade[, c(var, 'height', 'id')]
    colnames(gad) <- c('VALUE', 'HEIGHT', 'IDS')
    gad$VAR <- var
    
    GRAND <- rbind(GRAND, gad)
    
    
  }
  
  p <- ggplot(GRAND, aes(HEIGHT, VALUE, color=IDS))+
    facet_wrap(~VAR, scale='free')+
    geom_line( alpha=0.8, size=1)+
    scale_x_continuous(trans = 'reverse')+
    # scale_y_continuous(trans = 'log10')+
    theme_bw()
  
  
  
  if (plt){
    cairo_ps(filename =  sprintf("%s/plots/%s/%s/%s/%s.eps", path, animal, 
                                 foldername, subfolder, filename),
             width =12, height = 5,
             fallback_resolution = 200)
    print(p)
    dev.off()
  } else {
    print(p)
  }
  
  # p <- ggplot(grand.cascade, aes(den, nat))+
  #   # facet_wrap(~VAR, scale='free')+
  #   geom_line(aes(color=id))+
  #   scale_x_continuous(trans = 'reverse')+
  #   # scale_y_continuous(trans = 'log10')+
  #   theme_bw()
  # 
  # print(p)
  
} 

Mountain.View <- function(ham.merde, labels, net.coord, path='', nodes=40,
                          animal='monkey', foldername='', subfolder='',
                          plt=F, filename=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername,
                       subfolder), showWarnings = FALSE)
  }
  
  # net.coord <- net.coord[labels[1:nodes],]
  net.coord <- as.matrix(net.coord)

  x.coord <- net.coord[,1]
  y.coord <- net.coord[,2]
  z.coord <- in.the.heights(ham.merde, nodes)[1:nodes]
  
  s <- interp(x.coord, y.coord, z.coord, nx = 300, ny = 300)
  
  # regions <- regions[which(regions$AREA %in% labels[1:nodes]),]
  # regions <- regions[match(labels[1:nodes], regions$AREA),]
  
  # data.3d <- data.frame(x=x.coord,
  #                       y=y.coord,
  #                       z=z.coord[1:nodes],
  #                       REGION=as.factor(regions$REGION))
  
  data.3d <- data.frame(x=x.coord,
                        y=y.coord,
                        z=z.coord[1:nodes])
  
  surface.3d <- data.frame(source=rep(s$x, each=300), target=rep(s$y, 300), weight=Reshape(t(s$z), 300*300, 1))
  surface.3d <- surface.3d[!is.na(surface.3d$weight),]
  
  pp <- ggplot(surface.3d) +
    geom_tile(aes(x = source, y = target, fill = weight)) +
    geom_contour(aes(x = source, y = target, z = weight), color = "gray", 
                 binwidth=10)+
    scale_fill_viridis(option = 'A')+
    geom_point(data=data.3d, aes(x=x, y=y), size=10)+
    # scale_color_brewer(palette = 'Dark2')+
    annotate('text', data.3d$x,
             data.3d$y, label=labels[1:nodes],
             color='white', fontface=2, size=4)+
    coord_fixed()+
    theme_void()
  
  
  if (plt){
    cairo_ps(filename = sprintf("%s/%s/%s/%s.eps", path, foldername, subfolder, filename),
             width =13, height = 10,
             fallback_resolution = 200)
    plot_gg(pp, width = 7, height = 7, raytrace = FALSE, preview = TRUE)
    
    dev.off()  
  } else {
    plot_gg(pp, width = 7, height = 7, raytrace = FALSE, preview = TRUE)
  }
  
}

hierarchical_bundling.merde <- function(net, merde, nodes, k, k.m, net.cluster, 
                                        leaves, labels, regions, att, 
                                        path='', animal='monkey', foldername='', 
                                        subfolder='', filename='', plt=F){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername,
                       subfolder), showWarnings = FALSE)
  }
  
  hierarchy <- data.frame(from=c(), to=c())
  
  for (e in 1:(nodes-1)){
    
    cluster.i <- cutree(merde, k=e)
    cluster.f <- cutree(merde, k=e+1)
    cls.i <- unique(unique(cluster.i))
    cls.f <- unique(unique(cluster.f))
    
    father <- rep(0, length(cls.f))
    
    for (fi in 1:length(cls.f)){
      indx.fi <- cluster.f == cls.f[fi]
      father[fi] <- unique(cluster.i[indx.fi])
    }
    
    h.step <- data.frame(from=paste(sprintf('step_%i', e), father, sep = '_'), 
                         to=paste(sprintf('step_%i', e+1), cls.f, sep = '_'))
    
    hierarchy <- rbind(hierarchy, h.step)
  }
  
  vertices <- data.frame(name = unique(c(as.character(hierarchy$from), 
                                         as.character(hierarchy$to)))) 
  regions <- regionsCSV
  regions <- regions[match(labels[1:nodes], regions$AREA),]
  regions <- c(rep(NA, nrow(vertices)-nodes), regions$REGION)
  
  new.labels <- c(rep(NA, nrow(vertices)-nodes), labels[1:nodes])
  cluster <- extract.k.partition(net, net.cluster, nrow(net), nodes, k, 0)
  cluster <- c(rep(NA, nrow(vertices)-nodes), cluster)
  
  mygraph <- graph_from_data_frame( hierarchy, vertices=vertices)
  connect <- with(net, data.frame(from=paste(sprintf('step_%i', nodes), source, sep = '_'),
                                  to=paste(sprintf('step_%i', nodes), target, sep = '_'),
                                  weight=weight))
  net.plotting <- paint.network(net, net.cluster, nrow(net), k, 0, kactive = T)
  est.para <- est.parameters(net.plotting)
  net.plotting$commship[which(net.plotting$commship %in% 
                                est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net.plotting$commship[which(net.plotting$commship %in% 
                                est.para$commship[is.na(est.para$Dc)])] <- -1
  
  from <- match(connect$from, vertices$name)
  to <-  match(connect$to, vertices$name)
  
  connect$value <- net.plotting$commship
  if (att == -2){
    connect$weight <- log10(connect$weight) + 7
    connect$weight <- connect$weight/max(connect$weight)
    connect$weight <- (connect$weight - min(connect$weight)+0.01)/(1-min(connect$weight)+0.01)
  } else{
    connect$weight <- ifelse(connect$value %in% att, 0.6, 0.01)
  }
  p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    geom_conn_bundle(data = get_con(from = from, to = to, 
                                    col=as.factor(connect$value),
                                    w=connect$weight),
                     aes(colour=col, alpha=w), tension = 0.80) + 
    # scale_edge_colour_brewer(palette = "Dark2") +
    geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, 
                       label=new.labels, 
                       angle = ifelse(node_angle(x,y) > 90 & node_angle(x,y) < 270, node_angle(x,y) + 180, node_angle(x,y) ),
                       colour=as.factor(cluster),
                       hjust='outward'), size=5, alpha=1, fontface = "bold") +
    geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=as.factor(regions)), 
                    size=3, shape=15) +
    # scale_color_brewer(palette = 'Spectral')+
    theme_void()
  
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s.eps", path, animal, foldername, subfolder, filename),
             width =13, height = 10,
             fallback_resolution = 200)
    print(p)
    
    dev.off()  
  } else {
    print(p)
  }
}

feature_vectors <- function(net, self.loop=F, mode=''){
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  
  for (i in 1:nodes.t){
    aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
    if (self.loop){
      if (strcmp(mode, 'BETA')){
        aik[i,i] <- mean(net$weight[which(net$target == i)], na.rm = T,)
      } else if (strcmp(mode, 'ALPHA')){
        aik[i,i] <- mean(net.source$weight[which(net.source$source == i)], na.rm = T)
      }
    }
  }
  
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  
  for (i in 1:nodes.t){
    aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
    if (self.loop){
      if (strcmp(mode, 'BETA')){
        aki[i,i] <- mean(net.source$weight[which(net.source$source == i)], na.rm = T,)
      } else if (strcmp(mode, 'ALPHA')){
        aki[i,i] <- mean(net$weight[which(net$target == i)], na.rm = T)
      }
    }
  }
  
  if (strcmp(mode, 'GAMMA')){
    for (i in 1:nodes.t){
      corr.i <- jaccard.p(aki[i,], aik[i,])
      aki[i,i] <- corr.i
      aik[i,i] <- corr.i
    }
  }
  
  x <- list('aik'= aik, 
            'aki'= aki)
  
  return(x)
}

similarity.competitive.fast <- function(net, features, leaves, it_min=1,
                                        it_max=100){
  
  aik <- features$aik
  aki <- features$aki
  
  if (it_max == leaves)
    it_max <- it_max + 1
  
  net.sim.in <- data.frame(source=c(), target=c(), weight=c())
  net.sim.final <- data.frame(source=c(), target=c(), weight=c())
  
  for (e in 1:(leaves-1)){
    if (mod(e+it_min, as.integer((leaves+it_min)/10)) == 0){
      print(sprintf('%.2f%%', ((e+it_min)/(leaves+it_min))*10))
      print(Sys.time())
    }
    
    if (net$id[e] >= it_min &&  net$id[e] < it_max){
      leaf.id <- net$id[e]
      initial.node <- net$source[e]
      final.node <- net$target[e]
      final.links <- net[which(net$target == final.node & 
                               net$source != initial.node &
                               net$id > leaf.id),]
      in.links <- net[which(net$source == initial.node &
                            net$target != final.node &
                            net$id > leaf.id),]
      if (nrow(in.links) > 0){
        n.in <- nrow(in.links)
        net.sim.in <- rbind(net.sim.in, 
                            data.frame(source=rep(leaf.id, n.in),
                                       target=in.links$id,
                                       weight=unlist(lapply(1:nrow(in.links), 
                                                            function(x)
                                                              jaccard.p.fast(aki[final.node,], 
                                                                             aki[in.links$target[x],])
                                                                        )
                                                                 )
                                                   )
                            )
        
      }
      
      if (nrow(final.links) > 0){
        n.final <- nrow(final.links)
        net.sim.final <- rbind(net.sim.final, 
                               data.frame(source=rep(leaf.id, n.final),
                                           target=final.links$id,
                                           weight=unlist(lapply(1:nrow(final.links), 
                                                                function(x)
                                                                  jaccard.p.fast(aik[initial.node,],
                                                                                 aik[final.links$source[x],])
                                                                              
                                                                                
                                                                              )
                                                                       )
                                                         )
                               )
      }
    }
  }
  
  net.sim <- rbind(net.sim.in, net.sim.final)
  
  return(net.sim)
}

jaccard.p.fast <- function(x, y){
  
  n <- length(x)
  
  x.m <- matrix(rep(x, n), ncol = n, byrow = T)
  x.m <- x.m/diag(x.m)
  x.m[is.nan(x.m)] <- Inf
  
  y.m <- matrix(rep(y, n), ncol = n, byrow = T)
  y.m <- y.m/diag(y.m)
  y.m[is.nan(y.m)] <- Inf
  
  return(sum(1/apply(pmax(x.m, y.m), c(1), sum)))
}

normalize.net <- function(A){
  
  ny <- ncol(A)
  
  for (i in 1:ny)
    A[,i] <- A[,i]/sum(A[,i], na.rm = T)
  
  return(A)
}

shuffle.net <- function(net, steps=10^6){
  
  A <- from_dataframe_to_adjacency(net)
  
  nx <- max(net$source)
  ny <- max(net$target)
  
  e <- 0
  while(e < steps){
    hat <- sample(1:ny, 4, replace = T)
    if (hat[1] != hat[2] && hat[3] != hat[4]){
      cte <- A[hat[1], hat[2]]
      A[hat[1], hat[2]] <- A[hat[3], hat[4]]
      A[hat[3], hat[4]] <- cte
      e <- e + 1
    }
  }
  
  A[A == 0] <- NA
  
  A <- from_adjacency_to_dataframe(A)
  A <- A[!is.na(A$weight),]
  A$id <- 1:nrow(A)
  
  return(A)
}

extract.k.overlap <- function(net, net.cluster, leaves, nodes, k, best.height){
  
  net.reduced <- cluster.network(net, net.cluster, leaves, k, best.height, kactive = T)
  est.para <- est.parameters(net.reduced)
  net.reduced$commship[which(net.reduced$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net.reduced$commship[which(net.reduced$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  net.reduced <- net.reduced[which(net.reduced$commship != -1),]
  
  node.historia <- rep(0, nodes)
  
  coms <- sort(unique(net.reduced$commship))
  
  for (com in coms){
    net.com <- net.reduced[which(net.reduced$commship == com),]
    nds <- sort(unique(intersect(net.com$source, net.com$target)))
    if (length(nds) > 1)
      node.historia[nds] <- 1 + node.historia[nds]
  }
  return(node.historia)
}

get.order.ham <- function(k.remark, loc.max, net, net.cluster, leaves,
                          nodes, animal, labels, regions){
  
  historia <- rep(0, nodes)
  
  for (e in 1:length(k.remark)){
    if (mod(e,5)==0)
      print(e)
    k <- k.remark[e]
    best.height <- loc.max[e]
      
    net.reduced <- cluster.network(net, net.cluster, leaves, k, best.height, kactive = T)
    est.para <- est.parameters(net.reduced)
    net.reduced$commship[which(net.reduced$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
    net.reduced$commship[which(net.reduced$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
    net.reduced <- net.reduced[which(net.reduced$commship != -1),]
    
    coms <- sort(unique(net.reduced$commship))
    
    for (com in coms){
      net.com <- net.reduced[which(net.reduced$commship == com),]
      nds <- sort(intersect(unique(net.com$source), unique(net.com$target)))
      if (length(nds) > 1){
        historia[nds] <- historia[nds] + 1
      }
    }
  }
  
  node.k.com <- extract.k.partition(net, net.cluster,
                                    leaves, nodes, 41, best.height)
  
  labels <- labels[1:nodes]
  regions <- regions[which(regions$AREA %in% labels),]
  regions <- regions[match(labels, regions$AREA),]
  
  tree.merde <- with(regions, data.frame(area=AREA, regions=REGION, 
                                         historia=historia,
                                         commship=node.k.com,
                                         order=1:nodes))
  tree.merde <- tree.merde[with(tree.merde, order(commship, historia, decreasing = T)),]
  tree.merde$nodes <- 1:nodes
  
  return(tree.merde)
  
}

from_dataframe_to_adjacency <- function(net){
  
  M <- max(net$source)
  N <- max(net$target)
  
  A <- matrix(0, nrow = M, ncol = N)
  
  for (i in 1:nrow(net))
    A[net$source[i], net$target[i]] <- net$weight[i]
  
  return(A)
}

paint.edges <- function(net, st){
  
  netPaint <- net
  netPaint$tag <- ''
  
  for (e in 1:nrow(net)){
    if (strcmp(st$nature[which(st$node == net$source[e])], 'source') && 
        strcmp(st$nature[which(st$node == net$target[e])], 'target')){
      netPaint$tag[e] <- 'purple'
    } else if (strcmp(st$nature[which(st$node == net$source[e])], 'source') && 
               strcmp(st$nature[which(st$node == net$target[e])], 'transition')){
      netPaint$tag[e] <- 'red'
    } else if (strcmp(st$nature[which(st$node == net$source[e])], 'transition') && 
               strcmp(st$nature[which(st$node == net$target[e])], 'target')){
      netPaint$tag[e] <- 'blue'
    } else if (strcmp(st$nature[which(st$node == net$source[e])], 'transition') && 
               strcmp(st$nature[which(st$node == net$target[e])], 'transition')){
      netPaint$tag[e] <- 'gold'
    }
  }
  return(netPaint)
}
similarity.special.commship <- function(net, net.cls, leaves, type='jacc', type.sim = 'jaccp', self.loop = F){
  
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  for (i in 1:nodes.t){
    aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
    if (self.loop){
      aik[i,i] <- sum(net.source$weight[which(net.source$source == i)])/sum(net.source$source == i)
      # aik[i,i] <- exp(sum(log(net.source$weight[which(net.source$source == i)]))/sum(net.source$source == i))
    }
  }
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  for (i in 1:nodes.t){
    aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
    if (self.loop){
      aki[i,i] <- sum(net$weight[which(net$target == i)])/sum(net$target == i)
      # aki[i,i] <- exp(sum(log(net$weight[which(net$target == i)]))/sum(net$target == i))
    }
  }
  
  net.sim <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j && net.cls$id[i] <= leaves && net.cls$id[j] <= leaves && net.cls$tag[i] != net.cls$tag[j]){
        
        if (net.cls$source[i] == net.cls$source[j] && net.cls$target[i] != net.cls$target[j]) {
          if (type == 'jacc'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jacc(aki[net.cls$target[i],], aki[net.cls$target[j],])
          } else if (type == 'cos'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- cos.index(aki[net.cls$target[i],], aki[net.cls$target[j],])
          } else if (type == 'jaccw'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.w(aki[net.cls$target[i],], aki[net.cls$target[j],])
          } else if (type == 'jaccp'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.p(aki[net.cls$target[i],], aki[net.cls$target[j],])
          } else if (type == 'sd'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- sd.index(aki[net.cls$target[i],], aki[net.cls$target[j],], type.sim=type.sim)
          }
        } 
        else if (net.cls$source[i] != net.cls$source[j] && net.cls$target[i] == net.cls$target[j]) {
          if (type == 'jacc'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jacc(aik[net.cls$source[i],], aik[net.cls$source[j],])
          } else if (type == 'cos'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- cos.index(aik[net.cls$source[i],], aik[net.cls$source[j],])
          } else if (type == 'jaccw'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.w(aik[net.cls$source[i],], aik[net.cls$source[j],])
          } else if (type == 'jaccp'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.p(aik[net.cls$source[i],], aik[net.cls$source[j],])
          } else if (type == 'sd'){
            net.sim[net.cls$id[i],net.cls$id[j]] <- sd.index(aik[net.cls$source[i],], aik[net.cls$source[j],], type.sim=type.sim)
          }
        }
      }
    }
  }
  net.sim <- net.sim + t(net.sim)
  net.sim[which(net.sim == 0)] <- NA
  return(net.sim)
}
self.similarity.comparison <- function(net, net1, nroles1, nature1){
  
  absolute.scale <- c('gold', 'red', 'blue', 'purple')
  
  A1.1 <- matrix(NA, nrow = 4, ncol = 4)
  
  for (ix in 1:nroles1){
    for (iy in 1:nroles1){
      if (ix <= iy){
        
        rx <- which(absolute.scale == nature1[ix])
        ry <- which(absolute.scale == nature1[iy])
        
        if (strcmp(nature1[ix], 'red') && strcmp(nature1[iy], 'red')){
          netRoles <- net1[which(net1$tag == 'red'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.competitive.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
          
        } else if (strcmp(nature1[ix], 'blue') && strcmp(nature1[iy], 'blue')){
          netRoles <- net1[which(net1$tag == 'blue'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.competitive.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if (strcmp(nature1[ix], 'gold') && strcmp(nature1[iy], 'gold')){
          netRoles <- net1[which(net1$tag == 'gold'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.competitive.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if (strcmp(nature1[ix], 'purple') && strcmp(nature1[iy], 'purple')){
          netRoles <- net1[which(net1$tag == 'purple'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.competitive.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if ((strcmp(nature1[ix], 'red') && strcmp(nature1[iy], 'gold')) ||
                   (strcmp(nature1[ix], 'gold') && strcmp(nature1[iy], 'red'))){
          netRoles <- net1[which(net1$tag == 'red' | net1$tag == 'gold'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if ((strcmp(nature1[ix], 'blue') && strcmp(nature1[iy], 'gold')) ||
                   (strcmp(nature1[ix], 'gold') && strcmp(nature1[iy], 'blue'))){
          netRoles <- net1[which(net1$tag == 'blue' | net1$tag == 'gold'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if ((strcmp(nature1[ix], 'purple') && strcmp(nature1[iy], 'gold')) ||
                   (strcmp(nature1[ix], 'gold') && strcmp(nature1[iy], 'purple'))){
          netRoles <- net1[which(net1$tag == 'purple' | net1$tag == 'gold'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if ((strcmp(nature1[ix], 'purple') && strcmp(nature1[iy], 'red')) ||
                   (strcmp(nature1[ix], 'red') && strcmp(nature1[iy], 'purple'))){
          netRoles <- net1[which(net1$tag == 'purple' | net1$tag == 'red'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
          
        } else if ((strcmp(nature1[ix], 'purple') && strcmp(nature1[iy], 'blue')) ||
                   (strcmp(nature1[ix], 'blue') && strcmp(nature1[iy], 'purple'))){
          netRoles <- net1[which(net1$tag == 'purple' | net1$tag == 'blue'),]
          
          leaves <- nrow(netRoles)
          netRoles$id <- 1:leaves
          
          simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
          simDf <- from_adjacency_to_dataframe(simMat)
          simDf <- simDf[!is.na(simDf$weight),]
          simDf <- simDf[which(simDf$source > simDf$target),]
          
          if (nrow(simDf) > 0){
            meanSim <- mean(log10(simDf$weight), na.rm = T)
            A1.1[rx, ry] <- -1/meanSim
          } 
        }
      }
    }
  }
  # for (i in 1:4){
  #   for (j in 1:4){
  #     if ( i < j ){
  #       if (is.na(A1.1[i,j]) && !is.na(A1.1[j,i])){
  #         A1.1[i,j] <- A1.1[j,i]
  #       } else if (!is.na(A1.1[i,j]) && is.na(A1.1[j,i])){
  #         A1.1[j,i] <- A1.1[i,j]
  #       }
  #     }
  #   }
  # }
  return(A1.1)
}
similarity.comparison <- function(net, net1, net2, nroles1, nroles2, nature1, nature2){
  absolute.scale <- c('gold', 'red', 'blue', 'purple')
  
  A3 <- matrix(NA, nrow = 4, ncol = 4)
  
  for (ix in 1:nroles1){
    for (iy in 1:nroles2){
      
      rx <- which(absolute.scale == nature1[ix])
      ry <- which(absolute.scale == nature2[iy])
      
      if (strcmp(nature1[ix], 'red') && strcmp(nature2[iy], 'red')){
        
        netRoles1 <- net1[which(net1$tag == 'red'),]
        netRoles2 <- net2[which(net2$tag == 'red'),]
        
        netRoles1$tag <- 'red1'
        netRoles2$tag <- 'red2'
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'blue') && strcmp(nature2[iy], 'blue')){
        
        netRoles1 <- net1[which(net1$tag == 'blue'),]
        netRoles2 <- net2[which(net2$tag == 'blue'),]
        
        netRoles1$tag <- 'blue1'
        netRoles2$tag <- 'blue2'
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'gold') && strcmp(nature2[iy], 'gold')){
        
        netRoles1 <- net1[which(net1$tag == 'gold'),]
        netRoles2 <- net2[which(net2$tag == 'gold'),]
        
        netRoles1$tag <- 'gold1'
        netRoles2$tag <- 'gold2'
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'purple') && strcmp(nature2[iy], 'purple')){
        
        netRoles1 <- net1[which(net1$tag == 'purple'),]
        netRoles2 <- net2[which(net2$tag == 'purple'),]
        
        netRoles1$tag <- 'purple1'
        netRoles2$tag <- 'purple2'
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'red') && strcmp(nature2[iy], 'gold')){
        
        netRoles1 <- net1[which(net1$tag == 'red'),]
        netRoles2 <- net2[which(net2$tag == 'gold'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'gold') && strcmp(nature2[iy], 'red')){
        
        netRoles1 <- net1[which(net1$tag == 'gold'),]
        netRoles2 <- net2[which(net2$tag == 'red'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'blue') && strcmp(nature2[iy], 'gold')) {
        
        netRoles1 <- net1[which(net1$tag == 'blue'),]
        netRoles2 <- net2[which(net2$tag == 'gold'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'gold') && strcmp(nature2[iy], 'blue')) {
        
        netRoles1 <- net1[which(net1$tag == 'gold'),]
        netRoles2 <- net2[which(net2$tag == 'blue'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'purple') && strcmp(nature2[iy], 'gold')){
        
        netRoles1 <- net1[which(net1$tag == 'purple'),]
        netRoles2 <- net2[which(net2$tag == 'gold'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'gold') && strcmp(nature2[iy], 'purple')){
        
        netRoles1 <- net1[which(net1$tag == 'gold'),]
        netRoles2 <- net2[which(net2$tag == 'purple'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'purple') && strcmp(nature2[iy], 'red')) {
        
        netRoles1 <- net1[which(net1$tag == 'purple'),]
        netRoles2 <- net2[which(net2$tag == 'red'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'red') && strcmp(nature2[iy], 'purple')) {
        
        netRoles1 <- net1[which(net1$tag == 'red'),]
        netRoles2 <- net2[which(net2$tag == 'purple'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'purple') && strcmp(nature2[iy], 'blue')){
        
        netRoles1 <- net1[which(net1$tag == 'purple'),]
        netRoles2 <- net2[which(net2$tag == 'blue'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'blue') && strcmp(nature2[iy], 'purple')){
        
        netRoles1 <- net1[which(net1$tag == 'blue'),]
        netRoles2 <- net2[which(net2$tag == 'purple'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'red') && strcmp(nature2[iy], 'blue')){
        
        netRoles1 <- net1[which(net1$tag == 'red'),]
        netRoles2 <- net2[which(net2$tag == 'blue'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      } else if (strcmp(nature1[ix], 'blue') && strcmp(nature2[iy], 'red')){
        
        netRoles1 <- net1[which(net1$tag == 'blue'),]
        netRoles2 <- net2[which(net2$tag == 'red'),]
        
        netRoles <- rbind(netRoles1, netRoles2)
        
        leaves <- nrow(netRoles)
        netRoles$id <- 1:leaves
        
        simMat <- similarity.special.commship(net, netRoles, leaves, type = 'jaccp', self.loop = T)
        simDf <- from_adjacency_to_dataframe(simMat)
        simDf <- simDf[!is.na(simDf$weight),]
        simDf <- simDf[which(simDf$source > simDf$target),]
        
        if (nrow(simDf) > 0){
          meanSim <- mean(log10(simDf$weight), na.rm = T)
          A3[rx, ry] <- -1/meanSim
        } 
        
      }
    }
  }
  return(A3)
}

from_adjacency_to_dataframe <- function(A){
  M <- dim(A)[1]
  N <- dim(A)[2]
  
  weights <- Reshape(t(A), M*N, 1)
  target <- rep(1:N, N=M)
  source <- rep(1:M, each=N)
  df <- data.frame(source = source, target = target, weight = weights)
  return(df)
}

sim.type.coop <- function(net, leaves, type='jacc', type.sim = 'jaccp', self.loop = F){
  
  nodes.s <- 91
  nodes.t <- 40
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  for (i in 1:nodes.t){
    if (i %in% unique(net$source)){
      aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
      if (self.loop){
        aik[i,i] <- sum(net.source$weight[which(net.source$source == i)])/sum(net.source$source == i)
      }
    }
  }
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  for (i in 1:nodes.t){
    if (i %in% unique(net$target)){
      aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
      if (self.loop){
        aki[i,i] <- sum(net$weight[which(net$target == i)])/sum(net$target == i)
      }
    }
  }
  
  net.sim.typeI <- matrix(0, nrow = leaves, ncol = leaves)
  net.sim.typeII <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j && net$id[i] <= leaves && net$id[j] <= leaves){
        if (net$source[i] == net$source[j] && net$target[i] != net$target[j]) {
          if (type == 'jacc'){
            # net.sim[net$id[i],net$id[j]] <- jacc(aki[net$source[i],], aki[net$source[j],])
            net.sim.typeI[net$id[i],net$id[j]] <- jacc(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'cos'){
            net.sim.typeI[net$id[i],net$id[j]] <- cos.index(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'jaccw'){
            net.sim.typeI[net$id[i],net$id[j]] <- jaccard.w(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'jaccp'){
            net.sim.typeI[net$id[i],net$id[j]] <- jaccard.p(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'sd'){
            net.sim.typeI[net$id[i],net$id[j]] <- sd.index(aik[net$target[i],], aik[net$target[j],], type.sim=type.sim)
          }
        } 
        else if (net$source[i] != net$source[j] && net$target[i] == net$target[j]) {
          if (type == 'jacc'){
            net.sim.typeII[net$id[i],net$id[j]] <- jacc(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'cos'){
            net.sim.typeII[net$id[i],net$id[j]] <- cos.index(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'jaccw'){
            net.sim.typeII[net$id[i],net$id[j]] <- jaccard.w(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'jaccp'){
            net.sim.typeII[net$id[i],net$id[j]] <- jaccard.p(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'sd'){
            net.sim.typeII[net$id[i],net$id[j]] <- sd.index(aki[net$source[i],], aki[net$source[j],], type.sim=type.sim)
          }
        }
      }
    }
  }
  net.sim.typeI[which(net.sim.typeI == 0)] <- NA
  net.sim.typeII[which(net.sim.typeII == 0)] <- NA
  
  types = list()
  types$typeI <- net.sim.typeI
  types$typeII <- net.sim.typeII
  
  return(types)
}

sim.type.comp <- function(net, leaves, type='jacc', type.sim = 'jaccp', self.loop = F){
  
  nodes.s <- 91
  nodes.t <- 40
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  for (i in 1:nodes.t){
    if (i %in% unique(net$source)){
      aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
      if (self.loop){
        aik[i,i] <- sum(net.source$weight[which(net.source$source == i)])/sum(net.source$source == i)
      }
    }
  }
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  for (i in 1:nodes.t){
    if (i %in% unique(net$target)){
      aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
      if (self.loop){
        aki[i,i] <- sum(net$weight[which(net$target == i)])/sum(net$target == i)
      }
    }
  }
  
  net.sim.typeI <- matrix(0, nrow = leaves, ncol = leaves)
  net.sim.typeII <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j && net$id[i] <= leaves && net$id[j] <= leaves){
        if (net$source[i] == net$source[j] && net$target[i] != net$target[j]) {
          if (type == 'jacc'){
            net.sim.typeII[net$id[i],net$id[j]] <- jacc(aki[net$target[i],], aki[net$target[j],])
          } else if (type == 'cos'){
            net.sim.typeII[net$id[i],net$id[j]] <- cos.index(aki[net$target[i],], aki[net$target[j],])
          } else if (type == 'jaccw'){
            net.sim.typeII[net$id[i],net$id[j]] <- jaccard.w(aki[net$target[i],], aki[net$target[j],])
          } else if (type == 'jaccp'){
            net.sim.typeII[net$id[i],net$id[j]] <- jaccard.p(aki[net$target[i],], aki[net$target[j],])
          } else if (type == 'sd'){
            net.sim.typeII[net$id[i],net$id[j]] <- sd.index(aki[net$target[i],], aki[net$target[j],], type.sim=type.sim)
          }
        } 
        else if (net$source[i] != net$source[j] && net$target[i] == net$target[j]) {
          if (type == 'jacc'){
            # net.sim[net$id[i],net$id[j]] <- jacc(aki[net$source[i],], aki[net$source[j],])
            net.sim.typeI[net$id[i],net$id[j]] <- jacc(aik[net$source[i],], aik[net$source[j],])
          } else if (type == 'cos'){
            net.sim.typeI[net$id[i],net$id[j]] <- cos.index(aik[net$source[i],], aik[net$source[j],])
          } else if (type == 'jaccw'){
            net.sim.typeI[net$id[i],net$id[j]] <- jaccard.w(aik[net$source[i],], aik[net$source[j],])
          } else if (type == 'jaccp'){
            net.sim.typeI[net$id[i],net$id[j]] <- jaccard.p(aik[net$source[i],], aik[net$source[j],])
          } else if (type == 'sd'){
            net.sim.typeI[net$id[i],net$id[j]] <- sd.index(aik[net$source[i],], aik[net$source[j],], type.sim=type.sim)
          }
        }
      }
    }
  }
  net.sim.typeI[which(net.sim.typeI == 0)] <- NA
  net.sim.typeII[which(net.sim.typeII == 0)] <- NA
  
  types = list()
  types$typeI <- net.sim.typeI
  types$typeII <- net.sim.typeII
  
  return(types)
}

link.to.region.map <- function(clusters, total.net, animal='monkey', filename='',
                               plt=F){
  all.id.net <- data.frame()
  for (id in sort(unique(clusters$commship))){
    id.net <- linkComm.region(clusters, id, regions, labels)
    id.net$edge <- id.net$edge/total.net$edge
    id.net$llk <- id.net$llk/total.net$llk
    id.net$commship <- id
    all.id.net <- rbind(all.id.net, id.net)
    
  }
  p <- ggplot(all.id.net, aes(target, source, fill=edge))+
    facet_wrap(~commship)+
    geom_tile()+
    scale_fill_viridis(option='C')+
    ggtitle('Link community - Regions, edge ratio')+
    theme(axis.text.x = element_text(angle = 90))
  q <- ggplot(all.id.net, aes(target, source, fill=llk))+
    facet_wrap(~commship)+
    geom_tile()+
    scale_fill_viridis(option='C')+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle('LLK ratio')
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/EdgeDensity_Regions_%s.eps", animal, filename),
             width =15, height = 7.5, pointsize = 20,
             fallback_resolution = 500)
    print(p)
    dev.off()
    cairo_ps(filename = sprintf("plots/%s/LLKRatio_Regions_%s.eps", animal, filename),
             width =15, height = 7.5, pointsize = 20,
             fallback_resolution = 500)
    print(q)
    dev.off()
  } else{
    print(p)
    print(q)
  }
  
}

linkComm.region <- function(net, id, regions, labels){
  labels <- labels[1:max(c(net$target))]
  net <- net[which(net$commship %in% id),]
  # vertices <- sort(unique(c(net$source, net$target)))
  regions <- regions[which(regions$AREA %in% labels),]
  nr <- length(unique(regions$REGION))
  region <- sort(unique(regions$REGION))
  
  linkInfo <- data.frame(source = rep(region, nr), target=rep(region, each=nr), edge = rep(0, nr*nr), llk = rep(0, nr*nr))
  
  for (i in 1:nrow(net)){
    e.source <- net$source[i]
    e.target <- net$target[i]
    
    e.source <- labels[e.source]
    e.target <- labels[e.target]
    
    e.source.r <- regions$REGION[which(regions$AREA == e.source)]
    e.target.r <- regions$REGION[which(regions$AREA == e.target)]
    
    linkInfo$edge[which(linkInfo$source == e.source.r & linkInfo$target == e.target.r)] <- linkInfo$edge[which(linkInfo$source == e.source.r & linkInfo$target == e.target.r)] + 1
    
    w <- net$weight[i]
    linkInfo$llk[which(linkInfo$source == e.source.r & linkInfo$target == e.target.r)] <- linkInfo$llk[which(linkInfo$source == e.source.r & linkInfo$target == e.target.r)] + log(w)
  }
  
  return(linkInfo)
}
effective.infoFlux.regions.id <- function(net, id, region, labels){

  labels.nodes <- data.frame(nodes=1:length(labels), area=labels)
  labels.nodes <- labels.nodes[which(labels.nodes$nodes <= max(net$target)),]
  net.cluster <- net[which(net$commship %in% id),]
  
  nodes.cluster <- sort(unique(c(net.cluster$source, net.cluster$target)))
  areas.cluster <- labels.nodes$area[nodes.cluster]
  regions.cluster <- unique(region$REGION[which(region$AREA %in% areas.cluster)])
  
  info.local.flux.t <- rep(0, length(regions.cluster))
  info.local.flux.s <- rep(0, length(regions.cluster))
  info.local.flux <- rep(0, length(regions.cluster))
  
  info.region.flux.t <- rep(0, length(regions.cluster))
  info.region.flux.s <- rep(0, length(regions.cluster))
  info.region.flux <- rep(0, length(regions.cluster))
  # print('All areas')
  # print(areas.cluster)
  # print(nodes.cluster)
  for (i in 1:length(regions.cluster)){
    print(regions.cluster[i])
    areas.region <-  region$AREA[which(region$REGION==regions.cluster[i])]
    areas <- areas.region
    
    areas.region <- areas.region[which(areas.region %in% labels.nodes$area)]
    areas.region <- labels.nodes$nodes[match(areas.region, labels.nodes$area)]
    info.region.flux.t[i] <- info.region.flux.t[i] - sum(log2(net$weight[which(net$target %in% areas.region)]))
    info.region.flux.s[i] <- info.region.flux.s[i] - sum(log2(net$weight[which(net$source %in% areas.region)]))
    
    # print(areas)
    areas <- areas[which(areas %in% areas.cluster)]
    # print(areas)
    areas <- labels.nodes$nodes[match(areas, labels.nodes$area)]
    # print(areas)
    info.local.flux.t[i] <- info.local.flux.t[i] - sum(log2(net.cluster$weight[which(net.cluster$target %in% areas)]))
    info.local.flux.s[i] <- info.local.flux.s[i] - sum(log2(net.cluster$weight[which(net.cluster$source %in% areas)]))
    
  }
  
  
  info.local.flux <- info.local.flux.t + info.local.flux.s
  info.region.flux <- info.region.flux.t + info.region.flux.s
  effective.flux <- data.frame(region=regions.cluster, 
                               info.loc.s=info.local.flux.s/info.region.flux.s, 
                               info.loc.t=info.local.flux.t/info.region.flux.t,
                               total.loc=info.local.flux/info.region.flux)
  
  return(effective.flux)
}

plotSymmetries <- function(net, labels, id=F, plt=F, filename='', anima='monkey'){
  
  
  net.igraph <- graph_from_data_frame(net[,c('source', 'target')], directed = F)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  labels <-  labels[cver]
  V(net.igraph)$label <- labels

  if (id){
    colors <- rep(rgb(.5,.5,.5,.5), length(unique(labels)))
    colors[which(labels == id)] <- rgb(0.9,0.2,0.1,0.8)
    V(net.igraph)$color <- colors
  } else {
    pal <- gg_color_hue(length(unique(labels)))
    colors <- pal[match(labels, sort(unique(labels)))]
    V(net.igraph)$color <- colors
  }
  
  net$weight <- -log(net$weight)
  E(net.igraph)$weight <- net$weight 
  net.coord <- layout_with_kk(net.igraph, maxiter = 20000)
  # E(net.igraph)$color <- rgb(.5,.5,.5,0.1)
  
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/symPlot_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 500)
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         vertex.label.family = 'Helvetica',
         vertex.size = 1.5,
         edge.curved=.1,
         edge.color= NA,
         vertex.label.cex = 0.09)
    legend('right', legend = sort(unique(labels)), col=pal, pch = rep(15, length(pal)),
           bty = "n",
           pt.cex = 2,
           cex = 1)
    
    dev.off()
  } else {
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         vertex.label.family = 'Helvetica',
         vertex.size = 1.5,
         edge.curved=.1,
         edge.color=NA,
         vertex.label.cex = 0.09)
    legend('right', legend = sort(unique(labels)), col=pal, pch = rep(15, length(pal)),
           bty = "n",
           pt.cex = 2,
           cex = 1)
  }
  
}
symmestriesNet <- function(symmetries){
  nnodes <- nrow(symmetries)
  nleafs <- length(symmetries[!is.na(symmetries)])/2
  source <- rep(0, nleafs)
  target <- rep(0, nleafs)
  weight <- rep(0, nleafs)
  cnt <- 1
  for (i in 1:nnodes){
    for (j in 1:nnodes){
      if (i < j && !is.na(symmetries[i,j])){
        source[cnt] <- i
        target[cnt] <- j
        weight[cnt] <- symmetries[i,j]
        cnt <- cnt + 1
      }
    }
  }
  return( data.frame(source=source, target=target, weight=weight))
}

infoFlux.regions.id <- function(net, id, region, labels){
  
  # meanInfo <- - sum((net$weight %*% log2(net$weight)))/max(net$target)
  meanInfo <- 1
  labels.nodes <- data.frame(nodes=1:length(labels), area=labels)
  net.cluster <- net[which(net$commship %in% id),]
  
  nodes.cluster <- sort(unique(c(net.cluster$source, net.cluster$target)))
  areas.cluster <- labels.nodes$area[nodes.cluster]
  regions.cluster <- unique(region$REGION[which(region$AREA %in% areas.cluster)])
  
  info.region.flux.t <- rep(0, length(regions.cluster))
  info.region.flux.s <- rep(0, length(regions.cluster))
  info.region.flux <- rep(0, length(regions.cluster))

  for (i in 1:length(regions.cluster)){
    areas <-  region$AREA[which(region$REGION==regions.cluster[i])]
    areas <- areas[which(areas %in% areas.cluster)]
    areas <- labels.nodes$nodes[match(areas,labels.nodes$area)]
    
    for (j in 1:length(areas)){
      info.region.flux.t[i] <- info.region.flux.t[i] - sum(log2(net.cluster$weight[which(net.cluster$target == areas[j])]))
      info.region.flux.s[i] <- info.region.flux.s[i] - sum(log2(net.cluster$weight[which(net.cluster$source == areas[j])]))
    }
  }
  info.region.flux <- info.region.flux.t + info.region.flux.s
  return(data.frame(region=regions.cluster, info.s=info.region.flux.s/meanInfo, info.t=info.region.flux.t/meanInfo,
                    total=info.region.flux/meanInfo))
}
plot.flux.region.id <- function(info, id=0, plt=F, filename='', animal='monkey'){
  
  info <- info[order(info$total),]
  pal <- wes_palette("Zissou1", length(unique(info$region)), type = "continuous")
  
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/FluxRegion_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 300)
    
    p <- ggplot(info, aes(x=region, y=total, fill=as.factor(round(total,digits = 3))))+
      geom_bar(stat="identity", width=0.5)+
      scale_fill_manual(values=pal)+
      ylab(TeX('$\\frac{H_{t} + H_{s}}{H_{0}}$'))+
      ggtitle(id)+
      theme(axis.text.x = element_text(angle = 90))
    print(p)
    dev.off()
    cairo_ps(filename = sprintf("plots/%s/FluxFlux_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 300)
    p <- ggplot(info, aes(x=info.s, y=info.t, color=region))+
      geom_point(size=5)+
      scale_fill_manual(values=pal)+
      geom_abline(intercept = 0, slope = 1)+
      xlab(TeX(sprintf('\\frac{\\sum_{A\\in n_{c},A\\in R} h(A,%i,s)}{H_{0}}',id)))+
      ylab(TeX(sprintf('$\\frac{\\sum_{A\\in n^{c},A\\in R} h(A,%i,t)}{H_{0}}$',id)))+
      ggtitle(id)
    print(p)
    dev.off()
  } else {
    p <- ggplot(info, aes(x=region, y=total, fill=as.factor(round(total,digits = 3))))+
      geom_bar(stat="identity", width=0.5)+
      scale_fill_manual(values=pal)+
      ylab(TeX('$\\frac{H_{t} + H_{s}}{H_{0}}$'))+
      ggtitle(id)+
      theme(axis.text.x = element_text(angle = 90))
    print(p)
    
    p <- ggplot(info, aes(x=info.s, y=info.t, color=region))+
      geom_point(size=5)+
      geom_abline(intercept = 0, slope = 1)+
      xlab(TeX(sprintf('\\frac{\\sum_{A\\in n_{c},A\\in R} h(A,%i,s)}{H_{0}}',id)))+
      ylab(TeX(sprintf('$\\frac{\\sum_{A\\in n^{c},A\\in R} h(A,%i,t)}{H_{0}}$',id)))+
      ggtitle(id)
    print(p)
  }
  
}
plot.effflux.region.id <- function(info, id=0, plt=F, filename='', animal='monkey'){
  
  info <- info[order(info$total),]
  pal <- wes_palette("Zissou1", length(unique(info$region)), type = "continuous")
  
  if (length(id)>1){
    id <- -2
  }
  
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/effFluxRegion_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 300)
    
    p <- ggplot(info, aes(x=region, y=total.loc, fill=as.factor(round(total.loc,digits = 3))))+
      geom_bar(stat="identity", width=0.5)+
      scale_fill_manual(values=pal)+
      ylab(TeX(sprintf('$H(R,c)$',id,id)))+
      ggtitle(id)+
      theme(axis.text.x = element_text(angle = 90))
    print(p)
    dev.off()
    cairo_ps(filename = sprintf("plots/%s/effFluxFlux_%s.eps", animal, filename),
             width =15, height = 10, pointsize = 20,
             fallback_resolution = 300)
    p <- ggplot(info, aes(x=info.loc.s, y=info.loc.t, color=region))+
      geom_point(size=5)+
      scale_fill_manual(values=pal)+
      geom_abline(intercept = 0, slope = 1)+
      xlab(TeX(sprintf('$\\frac{\\sum_{A\\in n_{c}, A\\in R} h(A,%i,s)}{\\sum_{c\\in C} \\sum_{A\\in n_{c}, A\\in R} h(A,c,s)}$',id)))+
      ylab(TeX(sprintf('$\\frac{\\sum_{A\\in n^{c}, A\\in R} h(A,%i,t)}{\\sum_{c\\in C} \\sum_{A\\in n^{c}, A\\in R} h(A,c,t)}$',id)))+
      # coord_fixed(ratio = 1)+
      ggtitle(id)
    print(p)
    dev.off()
  } else {
   
    p <- ggplot(info, aes(x=region, y=total.loc, fill=as.factor(round(total.loc,digits = 3))))+
      geom_bar(stat="identity", width=0.5)+
      scale_fill_manual(values=pal)+
      ylab(TeX(sprintf('$H(R,c)$',id,id)))+
      ggtitle(id)+
      theme(axis.text.x = element_text(angle = 90))
    print(p)
    p <- ggplot(info, aes(x=info.loc.s, y=info.loc.t, color=region))+
      geom_point(size=5)+
      scale_fill_manual(values=pal)+
      geom_abline(intercept = 0, slope = 1)+
      xlab(TeX(sprintf('$\\frac{\\sum_{A\\in n_{c}, A\\in R} h(A,%i,s)}{\\sum_{c\\in C} \\sum_{A\\in n_{c}, A\\in R} h(A,c,s)}$',id)))+
      ylab(TeX(sprintf('$\\frac{\\sum_{A\\in n^{c}, A\\in R} h(A,%i,t)}{\\sum_{c\\in C} \\sum_{A\\in n^{c}, A\\in R} h(A,c,t)}$',id)))+
      ggtitle(id)
    print(p)
  }
  
}
information.flux.id <- function(net, id, labels){
  # meanInfo <- -sum(net$weight %*% log2(net$weight))/max(net$target)
  meanInfo <- 1
  net <- net[which(net$commship %in% id),]
  nodes.cluster <- sort(unique(c(net$source, net$target)))
  
  node.info <- rep(0, length(nodes.cluster))
  for (i in 1:length(nodes.cluster)){
    info.t <- 0
    info.s <- 0
    if (length(net$weight[net$target == nodes.cluster[i]]) > 0){
      info.t <- -sum(log2(net$weight[net$target == nodes.cluster[i]]))
    }
    if (length(net$weight[net$source == nodes.cluster[i]]) > 0){
      info.s <- -sum(log2(net$weight[net$source == nodes.cluster[i]]))
    }
    node.info[i] <- info.t - info.s
  }
  return(data.frame(area=labels[nodes.cluster], nodes=nodes.cluster, info=node.info/meanInfo))
}

plotNetwork.infoFlux.id <- function(net, dataInfo, id, labels, plt=F, coords=T, animal='monkey',
                                    fr=T, filename=''){
  
  
  pal <- wes_palette("Zissou1", length(unique(dataInfo$info)), type = "continuous")
  
  color <- rep(rgb(.5,.5,.5,0.2), length(unique(c(net$source, net$target))))
  
  
  net.igraph <- graph_from_data_frame(net[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  V(net.igraph)$label <- labels[cver]
  color[dataInfo$nodes] <- pal[match(dataInfo$info,sort(unique(dataInfo$info)))]
  V(net.igraph)$color <- color[cver]
  
  edge.color <- rep(rgb(.5,.5,.5,0), nrow(net))
  nodes.cluster <- dataInfo$nodes
  
  net.cluster <- net[which(net$commship %in% id),]
  
  for (i in 1:nrow(net.cluster)){
    edge.color[net.cluster$id[i]] <- rgb(186/255,85/255,211/255, 0.6)
  }
  E(net.igraph)$color <- edge.color
  
  if (fr){
    net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
  } else {
    net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
  }
  if (coords){
    if (animal=='monkey'){
      net.coord <- readRDS('RDS/flnMonkeyCoords.rds')
    } else if (animal=='mouse'){
      net.coord <- readRDS('RDS/flnMouseCoords.rds')
    }
    net.coord <- as.matrix(net.coord[labels[cver],])
  }
  if (length(id) > 1){
    id<- -2
  }
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/infoFlux_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 300)
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.label.cex = 0.65,
         vertex.label.font = 2,
         main=id)
    legend('left', legend=round(sort(unique(dataInfo$info)),2), pch = rep(15, length(pal)),
           ncol=3, col = pal, bty = "n", cex = 0.5, pt.cex=1)
    # image.plot(legend.only=T, zlim=range(sort(unique(dataInfo$info))), col=pal,
    #            horizontal = T, 
    #            legend.lab= TeX('$H_{t}-H_{s}/H_{0}$'))
    dev.off()
  }else{
   
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.label.cex = 0.65,
         vertex.label.font = 2,
         main=id)
    legend('left', legend=round(sort(unique(dataInfo$info)),2), pch = rep(15, length(pal)),
           ncol=3, col = pal, bty = "n", cex = 0.5, pt.cex=1)
    # image.plot(legend.only=T, zlim=r, col=pal,
    #            horizontal = T,
    #            legend.lab= TeX('$H_{t}-H_{s}/H_{0}$'))
  }
}
plotNetwork.st.edge.fraction <- function(net, labels, coords=T, plt=F, filename='',
                                         fr=T){
  
  nodes <- sort(unique(c(net$source, net$target)))
  
  st.fraction <- rep(0, length(nodes))
  
  for (i in 1:length(nodes)){
    st.fraction[nodes[i]] <- nrow(net[which(net$target == nodes[i]),])/nrow(net[which(net$source == nodes[i]),])
  }
  st.fraction <- round(st.fraction, digits = 2)
  num.fractions <- length(unique(st.fraction))
  
  # colors <- viridis(num.fractions, option = 'plasma')
  colors <- wes_palette("Zissou1", num.fractions, type = "continuous")
  
  colors <- colors[match(st.fraction, sort(unique(st.fraction)))]
  
  net.igraph <- graph_from_data_frame(net[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  V(net.igraph)$label <- labels[cver]
  V(net.igraph)$color <- colors[cver]
  E(net.igraph)$color <- rep(rgb(0.5,0.5,.5,.5), nrow(net))
  
  if (fr){
    E(net.igraph)$weight <- -1/log10(net$weight)
  }else {
    E(net.igraph)$weight <- -log10(net$weight)
  }
  
  if (fr){
    net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
  } else {
    net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
  }
  if (coords){
    if (animal=='monkey'){
      net.coord <- readRDS('RDS/flnMonkeyCoords.rds')
    } else if (animal=='mouse'){
      net.coord <- readRDS('RDS/flnMouseCoords.rds')
    }
    net.coord <- as.matrix(net.coord[labels[cver],])
  }
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/stFraction_%s.eps", animal, filename),
             width =10, height = 10, pointsize = 20,
             fallback_resolution = 300)
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.label.cex = 0.75,
         vertex.label.font = 2)
    image.plot(legend.only=T, zlim=range(sort(unique(st.fraction))), 
               col=wes_palette("Zissou1", num.fractions, type = "continuous"),
               horizontal=T)
    dev.off()
  } else{
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.label.cex = 0.75,
         vertex.label.font = 2)
    
    image.plot(legend.only=T, zlim=range(sort(unique(st.fraction))), col=wes_palette("Zissou1", num.fractions, type = "continuous"))
  }
}
plotNetwork.overlap <- function(net, overlap, labels, coords=T, animal='monkey', fr=T,
                                filename='', plt=F, versize=15){
  
  num.overlap <- length(unique(overlap$communities))
  colors <- wes_palette("Zissou1", num.overlap, type = "continuous")
  sort.overlap <- sort(unique(overlap$communities))
  colors <- colors[match(overlap$communities, sort.overlap)]
  
  
  net.igraph <- graph_from_data_frame(net[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  V(net.igraph)$label <- labels[cver]
  V(net.igraph)$color <- colors[cver]
  E(net.igraph)$color <- rep(rgb(0.5,0.5,.5,0), nrow(net))
  E(net.igraph)$weight <- net$weight
  if (fr){
    net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
  } else {
    net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
  }
  if (coords){
    if (animal=='monkey'){
      net.coord <- readRDS('RDS/flnMonkeyCoords.rds')
    } else if (animal=='mouse'){
      net.coord <- readRDS('RDS/flnMouseCoords.rds')
    }
    net.coord <- as.matrix(net.coord[labels[cver],])
  }
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/stOverlap_%s.eps", animal, filename),
             width =12, height = 10, pointsize = 20,
             fallback_resolution = 300)
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         vertex.size = versize,
         edge.curved=.1,
         vertex.label.cex = 0.6,
         vertex.label.font = 2,
         main='Total number of link communities for each vertex')
    legend('topleft',
           legend=round(sort.overlap, digits = 2), 
           col = wes_palette("Zissou1", num.overlap, type = "continuous"),
           pch = rep(15, num.overlap),
           bty = "n",
           pt.cex = 2,
           cex = 1,
           ncol = 1,
           text.col = "black",
           horiz = F)
    dev.off()
  } else{
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         vertex.size = versize,
         edge.curved=.1,
         vertex.label.cex = 0.6,
         vertex.label.font = 2,
         main='Total number of link communities for each vertex')
    legend('topleft',
           legend=round(sort.overlap, digits = 2), 
           col = wes_palette("Zissou1", num.overlap, type = "continuous"),
           pch = rep(15, num.overlap),
           bty = "n",
           pt.cex = 2,
           cex = 1,
           ncol = 1,
           text.col = "black",
           horiz = F)
  }
}
overlapping <- function(net, labels){
  
  nodes <- sort(unique(c(net$source, net$target)))
  net <- net[which(net$commship != -1),]
  num.communities.nodes <- rep(0, length(nodes))
  for (i in 1:length(nodes)){
    num.communities.nodes[nodes[i]] <- length(unique(net$commship[which(net$target == nodes[i] | net$source == nodes[i])])) #+ length(unique(net$commship[which(net$source == nodes[i])]))
  }
  drivenData <- data.frame(area=labels[nodes], node=nodes, communities=num.communities.nodes)
  return(drivenData)
}

plot.cluster.effectivity <- function(data, id=0, plt=F, filename='', animal='monkey'){
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/effPlot_%s.eps", animal,filename),
             width =7, height = 5, pointsize = 20,
             fallback_resolution = 300)
    
    p <- ggplot(data, aes(region, role, fill=participation)) + 
      scale_fill_gradientn(colours = wes_palette("Zissou1", length(unique(data$participation)), type = "continuous")) + 
      ggtitle(id)+
      geom_tile()
    print(p)
    dev.off()
  } else {
    ggplot(data, aes(region, role, fill=participation)) + 
      scale_fill_gradientn(colours = wes_palette("Zissou1", length(unique(data$participation)), type = "continuous")) + 
      ggtitle(id)+
      geom_tile()
  }
}

effective.participation <- function(metaData, net, st.cluster, labels){
  drivenData <- data.frame(area=c(), region=c(), participation=c(), role=c())
  
  nodes.cluster <- st.cluster$nodes[which(st.cluster$nature %in% c('source','target','transition'))]
  
  for (i in nodes.cluster){
    area <- labels[i]
    region <- metaData$REGION[which(metaData$AREA == area)]
    role <- st.cluster$nature[which(st.cluster$nodes == i)]
    if (strcmp(role,'source')){
      total.edges.area <- sum(net$source == i)
      area.eff.edges <- sum(net$source == i & net$target %in% nodes.cluster)
    } else if (strcmp(role,'target')){
      total.edges.area <- sum(net$target == i)
      area.eff.edges <- sum(net$target == i & net$source %in% nodes.cluster)
    } else if (strcmp(role,'transition')){
      total.edges.area <- sum(net$source == i)
      total.edges.area <- total.edges.area + sum(net$target == i)
      area.eff.edges.out <- sum(net$source == i & net$target %in% nodes.cluster)
      area.eff.edges.in <- sum(net$target == i & net$source %in% nodes.cluster)
      area.eff.edges <- area.eff.edges.in + area.eff.edges.out
    }
    eff <- area.eff.edges/total.edges.area
    drivenData <- rbind(drivenData, data.frame(area=area, region=region, 
                                               participation=eff, role=role))
    
  }
  superData <- data.frame(region=c(), role=c(), participation=c())
  regions.drivenData <- unique(drivenData$region)
  for (i in 1:length(regions.drivenData)){
    region.roles <- unique(drivenData$role[which(drivenData$region == regions.drivenData[i])])
    for (rl in region.roles){
      region.part <- sum(drivenData$participation[which(drivenData$region == regions.drivenData[i] &
                                                          drivenData$role == rl)])/sum(metaData$REGION == regions.drivenData[i])
      superData <- rbind(superData, data.frame(region=regions.drivenData[i], role=rl,
                                               participation=region.part))
    }
  }
  return(superData)
}

plot.area.region <- function(data, id=0, plt =F, filename='', animal='monkey'){
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/raPlot_%s.eps", animal,filename),
             width =10, height = 7.5, pointsize = 20,
             fallback_resolution = 300)
    if (length(unique(data$role)) == 2){
      if ('source' %in% unique(data$role) && 'transition' %in% unique(data$role)){
        p <- ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8), rgb(0.8,0.7,0.2,1)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      } else if ('target' %in% unique(data$role) && 'transition' %in% unique(data$role)){
        p <- ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.9,0.2,0.1,0.8), rgb(0.8,0.7,0.2,1)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      } else if ('source' %in% unique(data$role) && 'target' %in% unique(data$role)){
        p <- ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8),  rgb(0.9,0.2,0.1,0.8)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      }
     
    } else if (length(unique(data$role)) == 3){
      p <- ggplot(data, aes(area, region, fill=role))+
        scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8), rgb(0.9,0.2,0.1,0.8), rgb(0.8,0.7,0.2,1)))+
        theme(axis.text.x = element_text(angle = 90))+
        ggtitle(id)+
        geom_tile()
    } else if (length(unique(data$role)) == 1){
      p <- ggplot(data, aes(area, region, fill=role))+
        scale_fill_manual(values=rgb(0.8,0.7,0.2,1))+
        theme(axis.text.x = element_text(angle = 90))+
        ggtitle(id)+
        geom_tile()
    }
    print(p)
    dev.off()
  } else {
    if (length(unique(data$role)) == 2){
      if ('source' %in% unique(data$role) && 'transition' %in% unique(data$role)){
        ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8), rgb(0.8,0.7,0.2,1)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      } else if ('target' %in% unique(data$role) && 'transition' %in% unique(data$role)){
        ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.9,0.2,0.1,0.8), rgb(0.8,0.7,0.2,1)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      } else if ('source' %in% unique(data$role) && 'target' %in% unique(data$role)){
        ggplot(data, aes(area, region, fill=role))+
          scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8),  rgb(0.9,0.2,0.1,0.8)))+
          theme(axis.text.x = element_text(angle = 90))+
          ggtitle(id)+
          geom_tile()
      }
      
    } else if (length(unique(data$role)) == 3){
      ggplot(data, aes(area, region, fill=role))+
        scale_fill_manual(values=c(rgb(0.4,0.7,0.9,0.8), rgb(0.9,0.2,0.1,0.8), rgb(0.8,0.7,0.2,1)))+
        theme(axis.text.x = element_text(angle = 90))+
        ggtitle(id)+
        geom_tile()
    }
    else if (length(unique(data$role)) == 1){
      ggplot(data, aes(area, region, fill=role))+
        scale_fill_manual(values=rgb(0.8,0.7,0.2,1))+
        theme(axis.text.x = element_text(angle = 90))+
        ggtitle(id)+
        geom_tile()
    }
  }
}

meta.data.analysis <- function(metaData, st.cluster, labels){
  
  source.areas <- labels[st.cluster$nodes[which(st.cluster$nature == 'source')]]
  source.regions <- metaData$REGION[metaData$AREA %in% source.areas]
  target.areas <- labels[st.cluster$nodes[which(st.cluster$nature == 'target')]]
  target.regions <- metaData$REGION[metaData$AREA %in% target.areas]
  
  drivenData <- data.frame()
  
  drivenData <- rbind(drivenData, data.frame('area'=source.areas, 'region'=source.regions, 'role'=rep('source', length(source.areas))))
  drivenData <- rbind(drivenData, data.frame('area'=target.areas, 'region'=target.regions, 'role'=rep('target', length(target.areas))))
  
  if(sum(which(st.cluster$nature == 'transition')) > 0){
    st.areas <- labels[st.cluster$nodes[which(st.cluster$nature == 'transition')]]
    st.regions <- metaData$REGION[regions$AREA %in% st.areas]
    drivenData <- rbind(drivenData, data.frame('area'=st.areas, 'region'=st.regions, 'role'=rep('transition', length(st.areas))))
  }
  
  return(drivenData)
}

plotNetwork.st.id <- function(net, st.cluster, id,labels, coords, fr=F, animal='monkey',
                              plt=F, filename = '', versize=15,
                              foldername='', subfolder='', path=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername, subfolder), showWarnings = FALSE)
  }
  
  net.igraph <- graph_from_data_frame(net[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  V(net.igraph)$label <- labels[cver]
  V(net.igraph)$color <- st.cluster$color[cver]
  
  edge.color <- rep( rgb(0.5,0.5,0.5,0), nrow(net))
  
  st.nodes <- st.cluster$nodes[which(st.cluster$nature != 'UNK')]
  
  if (fr){
    E(net.igraph)$weight <- -1/log10(net$weight)
  }else {
    E(net.igraph)$weight <- -log10(net$weight)
  }
  
  net <- net[which(net$commship == id),]
  
  for (i in 1:nrow(net)){
    if (net$source[i] %in%  st.nodes && net$target[i] %in%  st.nodes){
      edge.color[net$id[i]] <- rgb(186/255,85/255,211/255, 0.6)
    } 
    if ((strcmp(st.cluster$nature[which(st.cluster$nodes == net$source[i])], 'transition')) &&
               (net$target[i] %in% st.nodes)){
      edge.color[net$id[i]] <- rgb(0.4,0.7,0.9,0.8)
    }  
    if ((strcmp(st.cluster$nature[which(st.cluster$nodes == net$target[i])], 'transition')) &&
                (net$source[i] %in% st.nodes)){
      edge.color[net$id[i]] <- rgb(0.9,0.2,0.1,0.8)
    }
    if ((strcmp(st.cluster$nature[which(st.cluster$nodes == net$target[i])], 'transition') ) &&
          (strcmp(st.cluster$nature[which(st.cluster$nodes == net$source[i])], 'transition') )){
      edge.color[net$id[i]] <- rgb(0.8,0.7,0.2,1)
    }
  }
  E(net.igraph)$color <- edge.color
  
  
  if (fr){
    net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
  } else {
    net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
  }
  
  net.coord <- coords
  net.coord <- as.matrix(net.coord[labels[cver],])
  
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s.eps", path, animal, foldername, subfolder, filename),
             width =10, height = 10, pointsize = 12,
             fallback_resolution = 300)
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.size=versize,
         vertex.label.cex = 0.75,
         vertex.label.font = 2,
         main=id)
    legend("topleft",
           legend = c('Source', 'Target', 'Transition', 'Target-Source', 'Source-Transition', 'Target-Transition', 'No active'),
           col = c(
             rgb(0.9,0.2,0.1,0.8),
             rgb(0.4,0.7,0.9,0.8),
             rgb(0.8,0.7,0.2,1),
             rgb(186/255,85/255,211/255, 0.6),
             rgb(0.9,0.2,0.1,0.8),
             rgb(0.4,0.7,0.9,0.8),
             rgb(0.5,0.5,0.5,0.5)),
           pch = c(19,19,19,NA,NA,NA,19),
           bty = "n",
           pt.cex = 2,
           cex = 1,
           ncol = 1,
           text.col = "black",
           horiz = F,
           lty = c(NA, NA, 1, 1, 1, 1, NA),
           lwd = c(NA, NA, 3, 3, 3, 3, NA))
    dev.off()
  } else{
    plot(net.igraph,
         layout = net.coord,
         edge.arrow.size=0.2,
         edge.curved=.1,
         vertex.label.family = 'Helvetica',
         vertex.label.color = rgb(1,1,1),
         edge.curved=.1,
         vertex.size=versize,
         vertex.label.cex = 0.75,
         vertex.label.font = 2,
         main=id)
    legend("topleft",
           legend = c('Source', 'Target', 'Transition', 'Target-Source','No active'),
           col = c(
             rgb(0.9,0.2,0.1,0.8),
             rgb(0.4,0.7,0.9,0.8),
             rgb(0.8,0.7,0.2,1),
             rgb(186/255,85/255,211/255, 0.6),
             rgb(0.5,0.5,0.5,0.5)),
           pch = c(19,19,19,NA,19),
           bty = "n",
           pt.cex = 2,
           cex = 1,
           ncol = 1,
           text.col = "black",
           horiz = F,
           lty = c(NA, NA, 1, 1, NA),
           lwd = c(NA, NA, 3, 3, NA))
  }
}

st.clustering <- function(net, ids){
  source <- rgb(0.9,0.2,0.1,0.8)
  target <- rgb(0.4,0.7,0.9,0.8)
  st <- rgb(0.8,0.7,0.2,0.8)
  gray <- rgb(0.5,0.5,0.5,0.5)
  
  st.net <- net
  st.net <- st.net[which(st.net$commship == ids),]
  
  nodes <- sort(unique(c(net$source, net$target)))
  st.nodes <- sort(unique(c(st.net$source, st.net$target)))
  
  st.cluster <- data.frame(nodes = nodes, color = rep(gray, length(nodes)),
                           nature = rep('UNK', length(nodes)))
  
  for (i in st.nodes){
    s.edges <- sum(which(st.net$source == i))
    t.edges <- sum(which(st.net$target == i))
    total.edges <- s.edges + t.edges
    fs.edges <- s.edges/total.edges
    if (fs.edges == 0.0){
      st.cluster$color[which(st.cluster$nodes == i)] <- target
      st.cluster$nature[which(st.cluster$nodes == i)] <- 'target'
    } else if (fs.edges > 0 && fs.edges < 1){
      st.cluster$color[which(st.cluster$nodes == i)] <- st
      st.cluster$nature[which(st.cluster$nodes == i)] <- 'transition'
    } else if (fs.edges == 1.0){
      st.cluster$color[which(st.cluster$nodes == i)] <- source
      st.cluster$nature[which(st.cluster$nodes == i)] <- 'source'
    }
  }
  return(st.cluster)
}

plotNetwork_by_tree <- function(net, net.merde, net.cluster, 
                                leaves, nodes, k, best.height, labels, 
                                animal='monkey', k.merde=-1, h.merde=-1, plt=F,
                                filename='', foldername='', subfolder='',
                                path=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername, subfolder), showWarnings = FALSE)
  }
  
  if (h.merde>0){
    cut.tree <- cutree(net.merde, h=h.merde)
  }
  if (k.merde>0){
    cut.tree <- cutree(net.merde, k=k.merde)
  }
  
  tr <- as.phylo(net.merde)
  tr <- fortify(tr)
  tr <- subset(tr, isTip)
  tr <- tr[order(tr$angle, decreasing = T),]
  
  A <- from_dataframe_to_adjacency(with(net[,1:3],
                                        data.frame(source=source,
                                                   target=target,
                                                   weight=log10(weight) +
                                                     ceil(-log10(min(weight))) +
                                                     1)))
  A <- A[tr$node, tr$node]
  A <- from_adjacency_to_dataframe(A)
  A <- A[which(A$weight != 0),]

  historia <- extract.k.partition(net, net.cluster, leaves, nodes, k, best.height)
  historia <- historia[match(tr$label, labels[1:nodes])]
  historia <- with.potholes(historia)

  u.hist <- historia[[1]]
  tab.t <- historia[[2]]
  tab2.t <- tab.t
  
  for (e in 2:length(tab.t)){
    tab.t[e] <- tab.t[e] + tab.t[e-1]
  }
  tab.t <- tab.t + 0.5
  tab.op <- c(0, tab.t)
  for (e in 1:length(tab2.t)){
    tab2.t[e] <- tab.t[e] - (tab.op[e+1] - tab.op[e])/2
  }
  
  tab2.t <- as.integer(tab2.t)
  
  p <- ggplot(A, aes(target, source, fill=weight))+
    geom_tile(color='white')+
    scale_x_continuous(breaks = 1:nodes, labels = tr$label,
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    scale_y_continuous(breaks = 1:nodes, labels = tr$label, trans = 'reverse',
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    scale_fill_viridis_c(direction = -1)+
    geom_hline(yintercept = tab.t[1:(length(tab.t)-1)], color='red')+
    geom_vline(xintercept = tab.t[1:(length(tab.t)-1)], color='red')+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s_venom.eps", path, animal, foldername, subfolder, filename),
             width =12, height = 10,
             fallback_resolution = 200)
    print(p)
    
    dev.off()
  }
  else{
    print(p)
  }
  
  net.reduced <- cluster.network(net, net.cluster, leaves, k, best.height, kactive = T)
  est.para <- est.parameters(net.reduced)
  net.reduced$commship[which(net.reduced$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net.reduced$commship[which(net.reduced$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  # net.reduced <- net.reduced[which(net.reduced$commship != -1),]
  
  net.reduced <- with(net.reduced, data.frame(source=source, target=target, weight=commship))
  A <- from_dataframe_to_adjacency(net.reduced)
  A <- A[tr$node, tr$node]
  A <- from_adjacency_to_dataframe(A)
  A <- A[which(A$weight != 0),]
  
  q <- ggplot(A, aes(target, source, fill=as.factor(weight)))+
    geom_tile(color='white')+
    geom_text(label=A$weight, size=3)+
    scale_x_continuous(breaks = 1:nodes, labels = tr$label,
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    scale_y_continuous(breaks = 1:nodes, labels = tr$label, trans = 'reverse',
                       sec.axis = dup_axis(labels = u.hist,
                                           breaks = tab2.t))+
    geom_hline(yintercept = tab.t[1:(length(tab.t)-1)], color='red')+
    geom_vline(xintercept = tab.t[1:(length(tab.t)-1)], color='red')+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s_commship.eps", path, animal, foldername, subfolder, filename),
             width =12, height = 10,
             fallback_resolution = 200)
    print(q)
    
    dev.off()
  }
  else{
    print(q)
  }
}

plotDendrogram.circular <- function(net.cluster, h = -1, k = -1, filename='', plt = T, scale = '', alternative = T,
                                    tip.labels.size = 3, legend.plot.size = 50,
                                    strip.text.size = 2, offset.strip.text = 0.6,
                                    offset.strip = 0.1, branch.none = F, line.size=3,
                                    polar.text.size = 10, bar.line.size=0.8,
                                    tree.regions=list(F,c()),
                                    foldername='', subfolder='',
                                    animal='monkey', path=''){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/plots/%s/%s/%s', path, animal, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.height <- max(net.cluster$height)
  if (plt){
    cairo_ps(filename = sprintf("%s/plots/%s/%s/%s/%s.eps", path, animal, foldername, subfolder, filename),
             width =20, height = 10,
             fallback_resolution = 200)
    if (h>0){
      cut.tree <- cutree(net.cluster, h=h)
    }
    if (k>0){
      cut.tree <- cutree(net.cluster, k=k)
    }
    tr <- as.phylo(net.cluster)
    tr <- tibble::as_tibble(tr)
    n <- length(unique(cut.tree))
    anc <- matrix(0, nrow = n, ncol = 1)
    anc.tips <- matrix(0, nrow = n, ncol = 2)
    ss <- cut.tree[net.cluster$order]
    sort.ss <- sort(unique(cut.tree))
    tg <- tibble('node' = c(0), 'cluster'= c(0))
    
    if (alternative){
      for (i in 1:n){
        sp <- tr$node[which(cut.tree == sort.ss[i])]
        anc.tips[i,1] <- sp[1]
        anc.tips[i,2] <- sp[length(sp)]
        anc[i] <- MRCA(tr, sp)$node
      }
      tg <- tibble( 'node' = c(0), 'cluster' = c(0))
      for (i in 1:n){
        offsp <- offspring(tr, anc[i])$node
        tg <- full_join(tg, tibble('node' = offsp, 'cluster' = rep(sort.ss[i],length(offsp))), by = c('node', 'cluster'))
      }
    } else {
      for (i in 1:n){
        sp <- tr$node[which(cut.tree == sort.ss[i])]
        tg <- full_join(tg, tibble('node' = sp, 'cluster' = sort.ss[i]), by = c('node', 'cluster'))
        anc.tips[i,1] <- sp[1]
        anc.tips[i,2] <- sp[length(sp)]
      }
    }
    tg <- tg[2:nrow(tg),]
    tr <- full_join(tr, tg, by = 'node')
    tr$cluster[is.na(tr$cluster)] <- 'NA'
    n <- length(unique(tr$cluster))
  if (tree.regions[[1]]){
    
    regions.name <- tr$label
    regions.name <- regions.name[!is.na(regions.name)]
    region.df <-tree.regions[[2]] 
    region.df <- region.df[which(region.df$AREA %in% regions.name),]
    region.df <- region.df[match(regions.name, region.df$AREA),]
    n.regions <- length(unique(region.df$REGION))
    tr$regions <- NA
    tr$regions[which(tr$label %in% region.df$AREA)] <- region.df$REGION[which(tr$label %in% region.df$AREA)]
    
    
  }
    
    tr <- as.treedata(tr)
    if (alternative){
      if (branch.none){
        tree.plot <- ggtree(tr,layout='circular', branch.length = 'none', aes(color = as.factor(cluster)), size=line.size)
      } else{
        tree.plot <- ggtree(tr,layout='circular', aes(color = as.factor(cluster)), size=line.size)
      } 
      tree.plot <- tree.plot + geom_tiplab(offset = offset.strip.text, cex=tip.labels.size)
    } else {
      if (branch.none){
        tree.plot <- ggtree(tr,layout='circular', branch.length = 'none', color='blue', size=line.size)
      } else{
        tree.plot <- ggtree(tr,layout='circular', color='blue', size=line.size)
      }
      
      tree.plot <- tree.plot + geom_tiplab(offset = offset.strip.text, cex=tip.labels.size, aes(color = as.factor(cluster)))
    }
    tree.plot <- tree.plot + 
      scale_color_manual(values = c(gg_color_hue(n-1),'gray'))
    if (tree.regions[[1]])
      tree.plot <- tree.plot + 
      geom_tippoint(aes(fill=regions), size=7, shape=22)
      # scale_fill_brewer(length(n.regions), palette = 'Dark2')

    print(tree.plot)
      
    dev.off()
  }
  else{
    if (h>0){
      cut.tree <- cutree(net.cluster, h=h)
    }
    if (k>0){
      cut.tree <- cutree(net.cluster, k=k)
    }
    tr <- as.phylo(net.cluster)
    tr <- tibble::as_tibble(tr)
    n <- length(unique(cut.tree))
    anc <- matrix(0, nrow = n, ncol = 1)
    anc.tips <- matrix(0, nrow = n, ncol = 2)
    ss <- cut.tree[net.cluster$order]
    sort.ss <- sort(unique(cut.tree))
    tg <- tibble('node' = c(0), 'cluster'= c(0))
    
    if (alternative){
      for (i in 1:n){
        sp <- tr$node[which(cut.tree == sort.ss[i])]
        # anc.lng[i] <- length(sp)
        anc.tips[i,1] <- sp[1]
        anc.tips[i,2] <- sp[length(sp)]
        anc[i] <- MRCA(tr, sp)$node
      }
      tg <- tibble( 'node' = c(0), 'cluster' = c(0))
      for (i in 1:n){
        offsp <- offspring(tr, anc[i])$node
        # if (anc.lng == 1){
        #   offsp <- c(parent(tr, anc[i])$node,anc[i])
        # }
        tg <- full_join(tg, tibble('node' = offsp, 'cluster' = rep(sort.ss[i],length(offsp))), by = c('node', 'cluster'))
      }
    } else {
      for (i in 1:n){
        sp <- tr$node[which(cut.tree == sort.ss[i])]
        tg <- full_join(tg, tibble('node' = sp, 'cluster' = sort.ss[i]), by = c('node', 'cluster'))
        anc.tips[i,1] <- sp[1]
        anc.tips[i,2] <- sp[length(sp)]
      }
    }
    tg <- tg[2:nrow(tg),]
    tr <- full_join(tr, tg, by = 'node')
    tr$cluster[is.na(tr$cluster)] <- 'NA'
    n <- length(unique(tr$cluster))
    tr <- as.treedata(tr)
    if (alternative){
      if (branch.none){
        tree.plot <- ggtree(tr,layout='circular', branch.length = 'none', aes(color = as.factor(cluster)), size=line.size)
      } else{
        tree.plot <- ggtree(tr,layout='circular', aes(color = as.factor(cluster)), size=line.size)
      } 
      tree.plot <- tree.plot + geom_tiplab(offset = offset.strip.text, cex=tip.labels.size)
    } else {
      if (branch.none){
        tree.plot <- ggtree(tr,layout='circular', branch.length = 'none', color='blue', size=line.size)
      } else{
        tree.plot <- ggtree(tr,layout='circular', color='blue', size=line.size)
      }
        
      tree.plot <- tree.plot + geom_tiplab(offset = offset.strip.text, cex=tip.labels.size, aes(color = as.factor(cluster)))
    }
    tree.plot <- tree.plot + scale_color_manual(values = c(gg_color_hue(n-1),'gray'))
   
    return(tree.plot)
  }
}

plotDendrogram_Matrix <- function(net, net.merde, net.cluster, labels, regions,
                                  h.merde = -1, k.merde = -1, filename='', plt = T, 
                                  alternative = T, tip.labels.size = 3, 
                                  branch.none = F, line.size=3, foldername='', 
                                  subfolder='', animal='monkey', path='',
                                  k=-1, height=0){
  
  if (!strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.height <- max(net.merde$height)
  net <- net[,c('source', 'target', 'weight', 'id')]
  A <- from_dataframe_to_adjacency(net[1:3])
  A[A == 0] <- NA
  A <- log10(A) + 7
  A <- A %>% as.data.frame()
  colnames(A) <- labels
  rownames(A) <- labels
  
  if (h.merde>0){
    cut.tree <- cutree(net.merde, h=h.merde)
  }
  if (k.merde>0){
    cut.tree <- cutree(net.merde, k=k.merde)
  }
  
  tr <- as.phylo(net.merde)
  
  trd <- fortify(tr)
  trd <- subset(trd, isTip)
  trd <- trd[order(trd$angle, decreasing = T),]
  
  historia <- extract.k.partition(net, net.cluster, nrow(net), nodes, k, height)
  historia <- historia[match(trd$label, labels[1:nodes])]
  historia <- with.potholes(historia)
  
  u.hist <- historia[[1]]
  tab.t <- historia[[2]]
  tab2.t <- tab.t
  
  for (e in 2:length(tab.t)){
    tab.t[e] <- tab.t[e] + tab.t[e-1]
  }
  tab.t <- tab.t + 0.5
  tab.op <- c(0, tab.t)
  for (e in 1:length(tab2.t)){
    tab2.t[e] <- tab.t[e] - (tab.op[e+1] - tab.op[e])/2
  }
  
  tab2.t <- as.integer(tab2.t)
  
  trd.names <- trd$label
  tr <- tibble::as_tibble(tr)
  
  n <- length(unique(cut.tree))
  anc <- matrix(0, nrow = n, ncol = 1)
  anc.tips <- matrix(0, nrow = n, ncol = 2)
  ss <- cut.tree[net.cluster$order]
  sort.ss <- sort(unique(cut.tree))
  tg <- tibble('node' = c(0), 'cluster'= c(0))
  
  if (alternative){
    for (i in 1:n){
      sp <- tr$node[which(cut.tree == sort.ss[i])]
      anc.tips[i,1] <- sp[1]
      anc.tips[i,2] <- sp[length(sp)]
      anc[i] <- MRCA(tr, sp)$node
    }
    
    tg <- tibble( 'node' = c(0), 'cluster' = c(0))
    for (i in 1:n){
      offsp <- offspring(tr, anc[i])$node
      tg <- full_join(tg, tibble('node' = offsp, 'cluster' = rep(sort.ss[i],length(offsp))), by = c('node', 'cluster'))
    }
  } else {
    for (i in 1:n){
      sp <- tr$node[which(cut.tree == sort.ss[i])]
      tg <- full_join(tg, tibble('node' = sp, 'cluster' = sort.ss[i]), by = c('node', 'cluster'))
      anc.tips[i,1] <- sp[1]
      anc.tips[i,2] <- sp[length(sp)]
    }
  }
  tg <- tg[2:nrow(tg),]
  tr <- full_join(tr, tg, by = 'node')
  tr$cluster[is.na(tr$cluster)] <- 'NA'
  n <- length(unique(tr$cluster))
  
  # region.df <- regions[which(regions$AREA %in% trd.names),]
  # region.df <- region.df[match(trd.names, region.df$AREA),]
  # tr$regions <- NA
  # tr$regions[which(tr$label %in% region.df$AREA)] <- region.df$REGION[which(tr$label %in% region.df$AREA)]
  # color.regions <- regions[order(regions$REGION),]

  tr <- as.treedata(tr)
  if (alternative){
    if (branch.none){
      tree.plot <- ggtree(tr,layout='rectangular', branch.length = 'none', aes(color = as.factor(cluster)), size=line.size)
    } else{
      tree.plot <- ggtree(tr,layout='rectangular', aes(color = as.factor(cluster)), size=line.size)
    } 
    tree.plot <- tree.plot #+ theme(axis.text.x = element_text(angle = 90))
  } else {
    if (branch.none){
      tree.plot <- ggtree(tr,layout='rectangular', branch.length = 'none', color='blue', size=line.size)
    } else{
      tree.plot <- ggtree(tr,layout='rectangular', color='blue', size=line.size)
    }
    
    tree.plot <- tree.plot + geom_tiplab(offset = 0.03, cex=tip.labels.size, aes(color = as.factor(cluster)))
  }
  
  tree.plot <- tree.plot + 
    scale_color_manual(values = c(gg_color_hue(n-1), 'grey'))+
    guides(color=guide_legend(title="Link community"))
  
  leg.t <- get_legend(tree.plot)
  
  A <- A[trd.names, trd.names] %>%
    from_adjacency_to_dataframe()
  A <- A[!is.na(A$weight),]
  
  A <- ggplot(A, aes(target, source))+
    geom_tile(aes(fill=weight))+
    scale_x_continuous(expand = c(0, 0), breaks = 1:nodes,
                       labels = trd$label)+
    scale_y_continuous(trans = 'reverse',expand = c(0, 0),breaks = 1:nodes,
                       labels = trd$label)+
    scale_fill_viridis(direction = -1, option = 'A')
  
  A <- A + geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(0.5,0.5)), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,0.5), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(tab.t[1],tab.t[1]), y=c(0.5,tab.t[1])), aes(x,y), color='green', size=1)+
    geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(tab.t[1],tab.t[1])), aes(x,y), color='green', size=1)
  
  xp <- tab.t[1]
  yp <- xp
  
  for (e in 2:(length(tab.t))){
    A  <- A + geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(yp,yp)), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,xp), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(tab.t[e],tab.t[e]), y=c(yp,tab.t[e])), aes(x,y), color='green', size=1)+
      geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(tab.t[e],tab.t[e])), aes(x,y), color='green', size=1)
    
    xp <- tab.t[e]
    yp <- xp
  }
  
  
  leg.a <- get_legend(A)
  
  A <- A +guides(fill=guide_legend(title = 'w'))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))
  
  # region.df$position <- 1
  # region.df$NODE <- 1:nodes
  # 
  # gg_leg <- ggplot(region.df, aes(y=NODE, x=position, fill=REGION))+
  #   geom_tile()+
  #   scale_x_continuous(expand = c(0, 0))+
  #   scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
  #   # scale_fill_brewer(palette = 'Dark2')+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   geom_text(label=trd.names, color='white', fontface=2)+ 
  #   theme_void()
  # 
  # leg.gg <- get_legend(gg_leg)
  # 
  # legend <- plot_grid(leg.t, leg.a, leg.gg, nrow = 3, rel_heights = c(2,1,5), axis = 't')
  legend <- plot_grid(leg.t, leg.a, nrow = 2, rel_heights = c(2,2), axis = 't')
  
  g <- plot_grid(tree.plot + theme(legend.position = 'none'),
                 A + theme(legend.position = 'none',
                           panel.grid = element_blank(),
                           panel.border = element_blank()),
                 # gg_leg + theme(legend.position = 'none',
                 #                panel.grid = element_blank(),
                 #                panel.border = element_blank()),
                 legend, ncol=3, rel_widths=c(3,15,1),
                 align = 'h' ,
                 scale=0.98)
  
  if (plt){
    cairo_ps(filename = sprintf("%s/%s/%s/%s.eps", path, foldername, subfolder, filename),
             width =16, height = 10,
             fallback_resolution = 200)
    print(g)
    
    dev.off()
  }
  else{
    print(g)
  }
}

histograms.id <- function(net, comms, plt=F, filename=''){
  net.histogram <- net[which(net$commship %in% comms), c('weight', 'commship')]
  
  p <- ggplot(net.histogram, aes(x=weight, fill = as.factor(commship)))+
    facet_wrap(~as.factor(commship))+
    geom_histogram(color="#e9ecef",alpha=0.6, position = "identity", lwd=0.2, bins=50)+
    scale_x_continuous(trans = 'log10')+
    ggtitle('Weight histogram')+
    theme_ipsum()+
    xlab('w')
  if (plt){
    cairo_ps(filename = sprintf("plots/histogramProb%s.eps", filename),
             width =15, height = 10, pointsize = 12,
             fallback_resolution = 300)
    print(p)
    dev.off()
  }
  else {
    return(p)
  }
}

plotNetwork.id <- function(net, labels, comms, plt=F, filename='', coords=T,
                           animal='monkey',fr=F, ncol=1, cex=1, lgnd= T,
                           versize=10) {
  net.list <- net
  mygray <- rgb(0.5,0.5,0.5,0.5)
  net.list$color <- mygray
  k <- length(comms)
  colors <- gg_color_hue(k)
  for (i in 1:k){
    net.list$color[ which(net.list$commship == comms[i])] <- colors[i]
  }
  print(length(net.list$color[which(net.list$color != mygray)]))
  colors <- c(colors, mygray)
  
  net.igraph <- graph_from_data_frame(net.list[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  if (fr){
    E(net.igraph)$weight <- -1/log10(net.list$weight)
  }else {
    E(net.igraph)$weight <- -log10(net.list$weight)
  }
  
  V(net.igraph)$label <- labels[cver]
  nodes <- max(c(net$source, net$target))
  values.percentage <- list()
  for (i in 1:nodes){
    nodes.comm.percentage <- c()
    ss <- net.list[which(net.list$source == cver[i]),]
    sp <- net.list[which(net.list$target == cver[i]),]
    st <- rbind(ss,sp)
    
    st <- st[!duplicated(st),]
    
    for (j in 1:(k+1)){
      nodes.comm.percentage <- c(nodes.comm.percentage, length(st$color[which(st$color==colors[j])]))
    }
    
    values.percentage[i] <- list( nodes.comm.percentage)
  }
  if (fr){
    net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
  } else {
    net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
  }
  if (coords){
    if (animal=='monkey'){
      net.coord <- readRDS('RDS/flnMonkeyCoords.rds')
    } else if (animal=='mouse'){
      net.coord <- readRDS('RDS/flnMouseCoords.rds')
    }
    
    net.coord <- as.matrix(net.coord[labels[cver],])
  }

  if (plt){
    cairo_ps(filename = sprintf("plots/%s/idNetwork%s.eps", animal, filename),
             width =10, height = 10, pointsize = 12,
             fallback_resolution = 300)
    # out.file.name <- paste("plots/animationFLN_", save.name, ".png", sep="")
    # png(out.file.name, width=640*2, height=480*2)
    if (coords){
      plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.pie = values.percentage,
           vertex.size = versize,
           vertex.label.family = 'Helvetica',
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.list$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2)
      if (lgnd){
        legend('topleft',legend = c(comms, 'others'), col = colors,
               pch = rep(19, k+1),
               lty = rep(1, k+1),
               lwd = rep(2, k+1),
               ncol = ncol,
               cex = cex,
               horiz = F)
      }
    }else {
      plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.pie = values.percentage,
           vertex.size = versize,
           vertex.label.family = 'Helvetica',
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.list$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2,
           main=filename)
      if (lgnd){
        legend('topleft',legend = c(comms, 'others'), col = colors,
               pch = rep(19, k+1),
               lty = rep(1, k+1),
               lwd = rep(2, k+1),
               ncol = ncol,
               cex = cex,
               horiz = F)
      }
    }
    dev.off()
  } else{
    if (coords){
      p <- plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.pie = values.percentage,
           vertex.label.family = 'Helvetica',
           vertex.size = versize,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.list$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2)
      if (lgnd){
        legend('topleft',legend = c(comms, 'others'), col = colors,
               pch = rep(19, k+1),
               lty = rep(1, k+1),
               lwd = rep(2, k+1),
               ncol = ncol,
               cex = cex,
               horiz = F)
      }
    }else {
      p <- plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.label.family = 'Helvetica',
           vertex.pie = values.percentage,
           vertex.size = versize,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.list$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2)
      if (lgnd){
        legend('topleft',legend = c(comms, 'others'), col = colors,
               pch = rep(19, k+1),
               lty = rep(1, k+1),
               lwd = rep(2, k+1),
               ncol = ncol,
               cex = cex,
               horiz = F)
      }
    }
  }
}

plotNetwork.filter <- function(network, labels, plt=F , fr = T, coord = F, filename = '', ncol = 2){
  network <- network[order(network$source),]
  netPlot <- network
  
  nodes <- length(unique((c(network$source, network$target))))
  leaves <- length(netPlot$source)
  
  netPlot$color <- rep('',leaves)
  
  ids.nods <- sort(unique((c(network$source, network$target))))
  ids.comm <- unique(netPlot$commship)
  k <- length(unique(netPlot$commship))
  
  rain.color <- gg_color_hue(k)
  # alpha.vertex <- rep(0,k)
  
  for (i in 1:k){
    netPlot$color[which(netPlot$commship == ids.comm[i])] <-  rain.color[i]
    # alpha.vertex[i] <- length(netPlot$commship[which(netPlot$commship == ids.comm[i])])/leaves
  }
  
  net.igraph <- graph_from_data_frame(netPlot[,c('source', 'target')], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  E(net.igraph)$weight <- netPlot$weight
  V(net.igraph)$label <- labels[cver]
  # V(net.igraph)$alpha <- alpha.vertex[cver]

  values.percentage <- list()
  for (i in 1:nodes){
    nodes.comm.percentage <- c()
    ss <- netPlot[which(netPlot$source == cver[i]),]
    sp <- netPlot[which(netPlot$target == cver[i]),]
    st <- rbind(ss,sp)
    
    st <- st[!duplicated(st),]

    for (j in 1:k){
      nodes.comm.percentage <- c(nodes.comm.percentage, length(st$color[which(st$color==rain.color[j])]))
    }
    
    values.percentage[i] <- list( nodes.comm.percentage)
  }
  if (fr){
    use.mode <- layout_with_fr
  } else {
    use.mode <- layout_with_kk
  }
  if (coord){
    net.coord <- readRDS('RDS/flndCoord.rds')
    net.coord <- as.matrix(net.coord[labels[cver],])
  }
  if (plt){
    cairo_ps(filename = sprintf("plots/filteredNetwork%s.eps", filename),
             width =10, height = 10, pointsize = 12,
             fallback_resolution = 300)
    # out.file.name <- paste("plots/animationFLN_", save.name, ".png", sep="")
    # png(out.file.name, width=640*2, height=480*2)
    if (coord){
      plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(rain.color),
           vertex.pie = values.percentage,
           vertex.size = 15,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = netPlot$color,
           vertex.label.cex=1,
           vertex.label.dist=0,
           vertex.label.font = 2)
      legend('topleft',legend = ids.comm, col = rain.color,
             pch = rep(19, k),
             lty = rep(1, k),
             lwd = rep(2, k),
             ncol = ncol,
             horiz = F)
    }else {
      plot(net.igraph,
           mode = use.mode,
           vertex.shape="pie",
           vertex.pie.color = list(rain.color),
           vertex.pie = values.percentage,
           vertex.size = 15,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = netPlot$color,
           vertex.label.cex=1,
           vertex.label.dist=0,
           vertex.label.font = 2)
      legend('topleft',legend = ids.comm, col = rain.color,
             pch = rep(19, k),
             lty = rep(1, k),
             lwd = rep(2, k),
             ncol = ncol,
             horiz = F)
    }
    dev.off()
  } else{
    if (coord){
      plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(rain.color),
           vertex.pie = values.percentage,
           vertex.size = 15,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = netPlot$color,
           vertex.label.cex=1,
           vertex.label.dist=0,
           vertex.label.font = 2)
      legend('topleft',legend = ids.comm, col = rain.color,
             pch = rep(19, k),
             lty = rep(1, k),
             lwd = rep(2, k),
             horiz = F)
    }else {
      plot(net.igraph,
           mode = use.mode,
           vertex.shape="pie",
           vertex.pie.color = list(rain.color),
           vertex.pie = values.percentage,
           vertex.size = 15,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = netPlot$color,
           vertex.label.cex=1,
           vertex.label.dist=0,
           vertex.label.font = 2)
      legend('topleft',legend = ids.comm, col = rain.color,
             pch = rep(19, k),
             lty = rep(1, k),
             lwd = rep(2, k),
             horiz = F)
    }
  }
}

filter.network.estimates <- function(network, cluster, estimates, labels,
                                     k=-1, h=-1,min.nodes=10){
  estimates <- estimates[which(estimates$numNodes >= min.nodes),]
  ids.comm <- estimates$idComm
  if (k > 0) {
    cut.tree <- cutree(cluster, k=k)
  }
  if (h > 0) {
    cut.tree <- cutree(cluster, h=h)
  }
  network$commship <- cut.tree
  filter.net <- data.frame()
  for (i in unique(cut.tree)){
    for (ids in ids.comm){
      if (i == ids){
        filter.net <- rbind(filter.net, network[which(cut.tree == i),])
      }
    }
  }
  nodes <- sort(unique(c(filter.net$source, filter.net$target)))
  nw.nodes <- 1:length(nodes)
  nw.labels <- rep('', length(nodes))
  for (i in nw.nodes){
    filter.net$source[which(filter.net$source == nodes[i])] <- i
    filter.net$target[which(filter.net$target == nodes[i])] <- i
    nw.labels[i] <- labels[nodes[i]]
  }
  return(list('filtered' = filter.net, 'labels' = nw.labels))
}

weighted.clusterting <- function(net, n, k=10, directed=F){
  if (directed){
    t.m <- n*(n-1)
  } else {
    t.m <- n*(n-1)/2
  }
  m <- nrow(net)
  tpw <- sum(net$weight)
  d.m <- t.m - m
  b.weights <- matrix(0, nrow = k, ncol = d.m)
  n.weights <- matrix(0, nrow = k, ncol = n - 1)
  for (e in 1:k){
    b.weights[e,] <- sample(net$weight, d.m, replace = T)
    n.weights[e,] <- sample(net$weight, n-1, replace = F)
  }
  
  b.weights <- apply(b.weights, 1, sum)
  n.weights <- apply(n.weights, 1, sum)
  
  return((tpw - mean(n.weights))/(tpw + mean(b.weights) - mean(n.weights)))
}

# find.height <- function(Dnetwork, curve='Dview', max.wire = 1, l_continuous=T){
#     if (curve=='Dview'){
#       f <- Dnetwork$Dview
#     }  else if (curve == 'Wview'){
#       f <- Dnetwork$Wview
#     }
#     best.f <- -Inf
#     best.height <- 0
#     max.h <- f[length(f)]
#     if (l_continuous){
#       for (i in 2:length(Dnetwork$height)){
#         if (i < length(Dnetwork$height) && abs((f[i+1]-f[i])/max.h) < 0.07){
#           if ((f[i] >= best.f && i/length(Dnetwork$height) <= max.wire)){
#             best.height <- Dnetwork$height[i]
#             best.f <- f[i]
#           }
#         } else if (i < length(Dnetwork$height) && abs((f[i+1]-f[i])/max.h) >= 0.07){
#           break
#         }
#       }
#     }else{
#       for (i in 2:length(Dnetwork$height)){
#         if (i < length(Dnetwork$height)){
#           if ((f[i] >= best.f && i/length(Dnetwork$height) <= max.wire)){
#             best.height <- Dnetwork$height[i]
#             best.f <- f[i]
#           }
#         } 
#       }
#     }
#     
#   return(best.height)
# }

plotDendrogram <- function(net.cluster, h = -1, k = -1, filename='', plt = T, scale = ''){
  max.height <- max(net.cluster$height)
  if (plt){
    cairo_ps(filename = sprintf("plots/dendogram%s.eps",filename),
             width =50, height = 10,
             fallback_resolution = 300)
    avg_dend_obj <- as.dendrogram(net.cluster)
    if (h>0){
      avg_col_dend <- color_branches(avg_dend_obj, h = h) # <- check
    }
    if (k>0){
      avg_col_dend <- color_branches(avg_dend_obj, k = k) # <- check
    }
    par(cex=0.25,font=3)
    if (scale == 'log'){
      plot(avg_col_dend,  log="y", ylim=c(1,max.height))
    } else {
      plot(avg_col_dend)
    }
    abline(h = h, col = 'red')
    dev.off()
  }
  else{
    avg_dend_obj <- as.dendrogram(net.cluster)
    if (h>0){
      avg_col_dend <- color_branches(avg_dend_obj, h = h) # <- check
    }
    if (k>0){
      avg_col_dend <- color_branches(avg_dend_obj, k = k) # <- check
    }
    par(cex=0.5,font=3)
    if (scale == 'log'){
      plot(avg_col_dend,  log="y", ylim=c(1,max.height))
    } else {
      plot(avg_col_dend)
    }
    abline(h = h, col = 'red')
  }
}

paint.network <- function(net, net.cluster, leaves, k,  bh, kactive=F){

  net.plotting <- net
  net.plotting$color <- rep(rgb(0.5,0.5,0.5,0.2),leaves)
  
  if (kactive){
    net.plotting$commship <- cutree(net.cluster, k=k)
  } else{
    net.plotting$commship <- cutree(net.cluster, h=bh)
  }
  colors <- gg_color_hue(k)
  for (i in 1:k){
    if (kactive){
      net.plotting$color[which(cutree(net.cluster, k=k) == i)] <-  colors[i]
    } else{
      net.plotting$color[which(cutree(net.cluster, h=bh) == i)] <-  colors[i] #<- check
    }
  }
  return(net.plotting)
}

paint.frac.network <- function(net, net.cluster, leaves, k,  bh, id.old){
  net.plotting <- net[which(net$id <= leaves),]
  net.plotting$color <- rep(rgb(.5,.5,.5,.2),leaves)
  net.plotting$commship <- rep(-1,leaves)
  colors <- gg_color_hue(k)
  tree.id <- cutree(net.cluster, h=bh)
  for (i in 1:length(tree.id)){
    net.plotting$commship[which(net.plotting$id == id.old[i])] <- tree.id[i]
  }
  for (i in 1:k){
    net.plotting$color[which(net.plotting$commship == i)] <-  colors[i]
  }
  return(net.plotting)
}

plotNetwork <- function(net.plotting, labels, nodes, k, coords = T, plt = T, kactive = F, 
                        fr =T, filename='', ncol=2, cex = 1, animal='monkey', lgnd = T,
                        versize=15){
  
  colors <- gg_color_hue(k)
  if (-1 %in% unique(net.plotting$commship)){
    colors <- c(gg_color_hue(k), rgb(.5,.5,.5,.2))
  }
  net.igraph <- graph_from_data_frame(net.plotting[,1:2], directed = T)
  cver <- as.numeric(rownames(as.matrix( V(net.igraph))))
  if(fr){
    E(net.igraph)$weight <- -1/log10(net.plotting$weight)
  }else{
    E(net.igraph)$weight <- -log10(net.plotting$weight)
  }

  V(net.igraph)$label <- labels[cver]
  
  if (coords){
    if (animal == 'monkey'){
      net.coord <- readRDS('RDS/flnMonkeyCoords.rds')
    } else if (animal == 'mouse'){
      net.coord <- readRDS('RDS/flnMouseCoords.rds')
    }
    
    net.coord <- net.coord[labels[cver],]
    net.coord <- as.matrix(net.coord)
  } else {
    if (fr){
      net.coord <- layout_with_fr(net.igraph, maxiter = 2000)
    } else{
      net.coord <- layout_with_kk(net.igraph, maxiter = 2000)
      
    }
    net.coord <- as.matrix(net.coord)
  }
  
  values.percentage <- list()
  
  for (i in 1:nodes){
    nodes.comm.percentage <- c()
    ss <- net.plotting[which(net.plotting$source == cver[i]),]
    sp <- net.plotting[which(net.plotting$target == cver[i]),]
    st <- rbind(ss,sp)
    st <- st[!duplicated(st),]
    if (-1 %in% unique(net.plotting$commship)){
      for (j in 1:(k+1)){
        nodes.comm.percentage <- c(nodes.comm.percentage, length(st$color[which(st$color==colors[j])]))
      }
    }else{
      for (j in 1:k){
        nodes.comm.percentage <- c(nodes.comm.percentage, length(st$color[which(st$color==colors[j])]))
      }
    }
    
    values.percentage[i] <- list(nodes.comm.percentage)
  }
  if (plt){
    cairo_ps(filename = sprintf("plots/%s/network_%s.eps", animal, filename),
             width =10, height = 10, pointsize = 12,
             fallback_resolution = 300)
      plot(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.pie = values.percentage,
           vertex.label.family = 'Helvetica',
           vertex.size = versize,
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.plotting$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2,
           main=filename)
      if (lgnd){
        legend('topleft',legend = 1:k, col = colors,
               pch = rep(19, k),
               lty = rep(1, k),
               lwd = rep(2, k),
               ncol = ncol,
               horiz = F)
      }
    dev.off()
  } else{
      plot.igraph(net.igraph,
           layout = net.coord,
           vertex.shape="pie",
           vertex.pie.color = list(colors),
           vertex.pie = values.percentage,
           vertex.size = versize,
           vertex.label.family = 'Helvetica',
           vertex.label.color = rgb(1,1,1),
           edge.arrow.size=0.2,
           edge.curved=.1,
           edge.color = net.plotting$color,
           vertex.label.cex=0.75,
           vertex.label.dist=0,
           vertex.label.font = 2,
           main=filename)
    if (lgnd){
      legend('topleft',legend = 1:k, col = colors,
             pch = rep(19, k),
             lty = rep(1, k),
             lwd = rep(2, k),
             ncol = ncol,
             cex = cex,
             horiz = F)
    }
  }
}
plotLinkCommunity <- function(net.plotting, nodes, labels, regionsCSV, plt = F){
  
  regions.s <- regionsCSV[order(regionsCSV$REGION),]
  regions.s <- regions.s[which(regions.s$AREA %in% labels[1:91]),] ##
  regions.t <- regionsCSV[order(regionsCSV$REGION),]
  regions.t <- regions.t[which(regions.t$AREA %in% labels[1:40]),]
  
  tab.t<- unname(table(regions.t$REGION))
  for (e in 2:length(tab.t)){
    tab.t[e] <- tab.t[e] + tab.t[e-1]
  }
  tab.t <- tab.t + 0.5
  
  tab2.t <- unname(table(regions.t$REGION))
  for (e in 1:length(tab2.t)){
    tab2.t[e] <- tab.t[e] - tab2.t[e]/1.5
  }
  
  tab2.t <- as.integer(tab2.t)+1
  
  nodes.plotting <- data.frame()
  nodes.comm.count <- c()
  for (i in 1:nodes){
    ss <- net.plotting[which(net.plotting$source == i), c("commship", 'id')]
    sp <- net.plotting[which(net.plotting$target == i), c("commship", 'id')]
    st <- rbind(ss, sp)
    st <- st[!duplicated(st$id),]
    st$node <- which(regions.t[,1] == labels[i])
    st$area <- labels[i]
    nodes.comm.count <- c(nodes.comm.count, length(unique(st$commship)))
    nodes.plotting <- rbind(nodes.plotting,st)
  }
  
  comm.count <- data.frame(commship=unique(net.plotting$commship), numEdges=rep(0, length(unique(net.plotting$commship))))
  for (i in unique(net.plotting$commship)){
    comm.count$numEdges[which(comm.count$commship == i)] <- length(net.plotting$commship[which(net.plotting$commship == i)])
  }
  
  # nodes.plotting <- nodes.plotting[match(regions.t[,1], labels[1:nodes]),]
  
  p1 <- ggplot(nodes.plotting, aes(as.factor(commship), node, color=as.factor(commship)))+
    geom_point(shape=15, size=2)+
    theme_classic()+
    geom_hline(yintercept = tab.t[1:(length(tab.t)-1)], color='black')+
    scale_color_manual(values=gg_color_hue(length(unique(net.plotting$commship))))+
    scale_y_continuous(breaks = 1:nodes,labels=regions.t[,1], 
                       sec.axis = dup_axis(labels = unique(regions.t[order(regions.t$REGION),2],),
                                                                                 breaks = tab2.t))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab('link community membership')+
    theme(legend.position = "none")
  
  nodes.comm.count <- data.frame(Areas = labels[1:nodes], linkComCount=nodes.comm.count, areas=1:40)
  nodes.comm.count <- nodes.comm.count[match(regions.t[,1], nodes.comm.count$Areas),]
  
  p2 <- ggplot(nodes.comm.count, aes(areas, linkComCount))+
    geom_bar(stat="identity", fill = "#FF6666")+
    geom_vline(xintercept = tab.t[1:(length(tab.t)-1)], color='black')+
    scale_x_continuous(labels = regions.t[order(regions.t$REGION),1], breaks = 1:40,
                       sec.axis = dup_axis(labels = unique(regions.t[order(regions.t$REGION),2],),
                                           breaks = tab2.t))+
    coord_flip()+
    theme_classic()
  
  p3 <- ggplot(comm.count, aes(x=as.factor(commship), y= numEdges, fill=as.factor(commship)))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=gg_color_hue(length(unique(net.plotting$commship))))+
    # scale_x_continuous(limits=1:length(unique(net.plotting$commship))) +
    theme_classic()+
    xlab('link community membership')+
    theme(axis.text.x = element_text(angle = 90))+
    theme(legend.position="none")
  
  if (plt){
    cairo_ps(filename = sprintf("plots/membershipLinkComm%s.eps",filename),
             width =20, height = 10, pointsize = 12,
             fallback_resolution = 300)
    
    grid.arrange(p1,p2,p3, layout_matrix = rbind(c(3,3,NA),c(1,1,2)))
    
    dev.off()
  } else{
    grid.arrange(p1,p2,p3, layout_matrix = rbind(c(3,3,NA),c(1,1,2)))
  }
  
}

is.jaccard.p <- function(x, y, i, j){
  
  nl <- length(x)
  snl <- 1:nl
  
  a.x <- rep(0, nl-1)
  a.y <- rep(0, nl-1)
  
  a.x[2:(nl-1)] <- x[which(!(snl %in% c(i,j)))]
  a.y[2:(nl-1)] <- y[which(!(snl %in% c(i,j)))]
  
  a.x[1] <- x[j]
  a.y[1] <- y[i]
  
  jp <- 0
  for (e in 1:length(a.x)){
    max_xy <- 0
    if (x[e] > 0 && y[e] > 0){
      for (z in 1:length(a.x)){
        max_xy <- max_xy + max(a.x[z]/a.x[e],a.y[z]/a.y[e])
      }
      jp <- jp + 1/max_xy
    }
  }
  return(jp)
}

jaccard.p <- function(x, y){
  jp <- 0
  for (i in 1:length(x)){
    max_xy <- 0
    if (x[i] > 0 && y[i] > 0){
      for (j in 1:length(y)){
        max_xy <- max_xy + max(x[j]/x[i],y[j]/y[i])
      }
      jp <- jp + 1/max_xy
    }
  }
  return(jp)
}

jaccard.w <- function(x,y){
  num <- 0
  den <- 0
  nodes <- length(x)
  for (i in 1:nodes){
      num <- num + min(x[i],y[i])
      den <- den + max(x[i],y[i])
  }
  return(num/den)
}

jacc <- function(x,y){
  ai.j <- t(x)%*%y
  ai.i <- t(x)%*%x
  aj.j <- t(y)%*%y
  return(ai.j/(ai.i+aj.j-ai.j))
}

sd.index <- function(x,y, type.sim='jaccp'){
  num <- 0
  den <- 0
  nodes <- length(x)
  
  if (type.sim == 'jaccp'){
    n <- jaccard.p(x,y)
  } else if (type.sim == 'jacc2'){
    n <- jaccard.v2(x,y)
  } 
  
  x[which(x>0)]<-1
  y[which(y>0)]<-1
  
  dim <- !(x | y)
  kx <- sum(x)
  ky <- sum(y)
  ndim <- sum(dim)
  
  Delta <- (ndim/nodes  - (1-(kx/nodes))*(1- (ky/nodes)))
  n <- n*(1 + Delta)
  return(n)
}

cos.index <- function(x,y){
  ai.j <- t(x)%*%y
  ai.i <- t(x)%*%x
  aj.j <- t(y)%*%y
  return(ai.j/sqrt(ai.i*aj.j))
}

# similarity.competitive.commship <- function(net, net.cls, leaves, type='jacc', type.sim = 'jaccp', self.loop = F){
#   
#   nodes.s <- max(net$source)
#   nodes.t <- max(net$target)
#   aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
#   net.source <- net[which(net$source <= nodes.t),]
#   for (i in 1:nodes.t){
#     aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
#     if (self.loop){
#       aik[i,i] <- sum(net.source$weight[which(net.source$source == i)])/sum(net.source$source == i)
#       # aik[i,i] <- exp(sum(log(net.source$weight[which(net.source$source == i)]))/sum(net.source$source == i))
#     }
#   }
#   aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
#   for (i in 1:nodes.t){
#     aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
#     if (self.loop){
#       aki[i,i] <- sum(net$weight[which(net$target == i)])/sum(net$target == i)
#       # aki[i,i] <- exp(sum(log(net$weight[which(net$target == i)]))/sum(net$target == i))
#     }
#   }
#   
#   net.sim <- matrix(0, nrow = leaves, ncol = leaves)
#   
#   for (i in 1:leaves){
#     for (j in 1:leaves){
#       if (i < j && net.cls$id[i] <= leaves && net.cls$id[j] <= leaves){
#         if (net.cls$source[i] == net.cls$source[j] && net.cls$target[i] != net.cls$target[j]) {
#           if (type == 'jacc'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jacc(aki[net.cls$target[i],], aki[net.cls$target[j],])
#           } else if (type == 'cos'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- cos.index(aki[net.cls$target[i],], aki[net.cls$target[j],])
#           } else if (type == 'jaccw'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.w(aki[net.cls$target[i],], aki[net.cls$target[j],])
#           } else if (type == 'jaccp'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.p(aki[net.cls$target[i],], aki[net.cls$target[j],])
#           } else if (type == 'sd'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- sd.index(aki[net.cls$target[i],], aki[net.cls$target[j],], type.sim=type.sim)
#           }
#         } 
#         else if (net.cls$source[i] != net.cls$source[j] && net.cls$target[i] == net.cls$target[j]) {
#           if (type == 'jacc'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jacc(aik[net.cls$source[i],], aik[net.cls$source[j],])
#           } else if (type == 'cos'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- cos.index(aik[net.cls$source[i],], aik[net.cls$source[j],])
#           } else if (type == 'jaccw'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.w(aik[net.cls$source[i],], aik[net.cls$source[j],])
#           } else if (type == 'jaccp'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- jaccard.p(aik[net.cls$source[i],], aik[net.cls$source[j],])
#           } else if (type == 'sd'){
#             net.sim[net.cls$id[i],net.cls$id[j]] <- sd.index(aik[net.cls$source[i],], aik[net.cls$source[j],], type.sim=type.sim)
#           }
#         }
#       }
#     }
#   }
#   net.sim <- net.sim + t(net.sim)
#   net.sim[which(net.sim == 0)] <- NA
#   return(net.sim)
# }

similarity.cooperative <- function(net, leaves, type='jacc', type.sim = 'jaccp', self.loop = F){
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  for (i in 1:nodes.t){
    aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
    if (self.loop){
      aik[i,i] <- sum(net.source$weight[which(net.source$source == i)])/sum(net.source$source == i)
      # aik[i,i] <- exp(sum(log(net.source$weight[which(net.source$source == i)]))/sum(net.source$source == i))
    }
  }
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  for (i in 1:nodes.t){
    aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
    if (self.loop){
      aki[i,i] <- sum(net$weight[which(net$target == i)])/sum(net$target == i)
      # aki[i,i] <- exp(sum(log(net$weight[which(net$target == i)]))/sum(net$target == i))
    }
  }
  
  net.sim <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j && net$id[i] <= leaves && net$id[j] <= leaves){
        if (net$source[i] == net$source[j] && net$target[i] != net$target[j]) {
          if (type == 'jacc'){
            # net.sim[net$id[i],net$id[j]] <- jacc(aki[net$source[i],], aki[net$source[j],])
            net.sim[net$id[i],net$id[j]] <- jacc(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'cos'){
            net.sim[net$id[i],net$id[j]] <- cos.index(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'jaccw'){
            net.sim[net$id[i],net$id[j]] <- jaccard.w(aik[net$target[i],], aik[net$target[j],])
          } else if (type == 'jaccp'){
            net.sim[net$id[i],net$id[j]] <- jaccard.p(aik[net$target[i],], aik[net$target[j],])
            # net.sim[net$id[i],net$id[j]] <- jaccard.p(aki[net$target[i],], aki[net$target[j],])
          } else if (type == 'sd'){
            net.sim[net$id[i],net$id[j]] <- sd.index(aik[net$target[i],], aik[net$target[j],], type.sim=type.sim)
          }
        } 
        else if (net$source[i] != net$source[j] && net$target[i] == net$target[j]) {
          if (type == 'jacc'){
            net.sim[net$id[i],net$id[j]] <- jacc(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'cos'){
            net.sim[net$id[i],net$id[j]] <- cos.index(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'jaccw'){
            net.sim[net$id[i],net$id[j]] <- jaccard.w(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'jaccp'){
            net.sim[net$id[i],net$id[j]] <- jaccard.p(aki[net$source[i],], aki[net$source[j],])
          } else if (type == 'sd'){
            net.sim[net$id[i],net$id[j]] <- sd.index(aki[net$source[i],], aki[net$source[j],], type.sim=type.sim)
          }
        }
        #
        # else if (net$source[i] != net$source[j] && net$target[i] != net$target[j] &&
        #            net$target[i] == net$source[j] && net$source[i] != net$target[j]){
        #   ai.j <- t(aik[net$source[i],])%*%aki[net$target[j],]
        #   ai.i <- t(aik[net$source[i],])%*%aik[net$source[i],]
        #   aj.j <- t(aki[net$target[j],])%*%aki[net$target[j],]
        #   if (type == 'jacc'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
        #   } else if (type == 'cos'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
        #   }
        # } else if (net$source[i] != net$source[j] && net$target[i] != net$target[j] &&
        #            net$source[i] == net$target[j] && net$target[i] != net$source[j]){
        #   ai.j <- t(aki[net$target[i],])%*%aik[net$source[j],]
        #   ai.i <- t(aki[net$target[i],])%*%aki[net$target[i],]
        #   aj.j <- t(aik[net$source[j],])%*%aik[net$source[j],]
        #   if (type == 'jacc'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
        #   } else if (type == 'cos'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
        #   }
        # }
        #
      }
    }
  }
  net.sim <- net.sim + t(net.sim)
  net.sim[which(net.sim == 0)] <- NA
  return(net.sim)
}
similarity.coop.comp <- function(net ,nodes, type='jacc'){
  aik <- matrix(0, nrow = nodes, ncol = nodes)
  for (i in unique(net$source)){
    aik[i, net$target[which(net$source == i)]] <- net$weight[which(net$source == i)]
    aik[i,i] <- sum(net$weight[which(net$source == i)])/length(net$weight[which(net$source == i)])
  }
  aki <- matrix(0, nrow = nodes, ncol = nodes)
  for (i in unique(net$target)){
    aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
    aki[i,i] <- sum(net$weight[which(net$target == i)])/length(net$weight[which(net$target == i)])
  }
  
  net.sim <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j) {
        if (net$source[i] == net$source[j] && net$target[i] != net$target[j]) {
          ai.j <- t(aki[net$target[i],]+aik[net$target[i],])%*%(aki[net$target[j],]+aik[net$target[j],]) # +  t(aki[net$target[i],])%*%aik[net$target[j],] +  t(aik[net$target[i],])%*%aki[net$target[j],]
          # ai.j <- t(aki[net$target[i],])%*%aki[net$target[j],]
          # ai.j <- ai.j + t(aik[net$target[i],])%*%aik[net$target[j],]/sigma
          # ai.j <- t(aik[net$target[i],])%*%aik[net$target[j],]/sigma
          ai.i <- t(aki[net$target[i],]+aik[net$target[i],])%*%(aki[net$target[i],]+aik[net$target[i],]) #+  t(aki[net$target[i],])%*%aik[net$target[i],] +  t(aik[net$target[i],])%*%aki[net$target[i],]
          # ai.i <- ai.i + t(aik[net$target[i],])%*%aik[net$target[i],]/sigma
          # ai.i <- t(aik[net$target[i],])%*%aik[net$target[i],]/sigma
          aj.j <- t(aki[net$target[j],]+aik[net$target[j],])%*%(aki[net$target[j],]+aik[net$target[j],]) #+  t(aki[net$target[j],])%*%aik[net$target[j],] +  t(aik[net$target[j],])%*%aki[net$target[j],]
          # aj.j <- aj.j + t(aik[net$target[j],])%*%aik[net$target[j],]/sigma
          # aj.j <- t(aik[net$target[j],])%*%aik[net$target[j],]/sigma
          if (type == 'jacc'){
            net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
          } else if (type == 'cos'){
            net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
          }
        } else if (net$source[i] != net$source[j] && net$target[i] == net$target[j]) {
          ai.j <- t(aik[net$source[i],]+aki[net$source[i],])%*%(aik[net$source[j],]+aki[net$source[j],]) #+ t(aik[net$source[i],])%*%aki[net$source[j],] + t(aki[net$source[i],])%*%aik[net$source[j],]
          # ai.j <- t(aik[net$source[i],])%*%aik[net$source[j],]
          # ai.j <- ai.j + t(aki[net$source[i],])%*%aki[net$source[j],]/sigma
          # ai.j <- t(aki[net$source[i],])%*%aki[net$source[j],]/sigma
          ai.i <- t(aik[net$source[i],]+aki[net$source[i],])%*%(aik[net$source[i],]+aki[net$source[i],])  #+ t(aik[net$source[i],])%*%aki[net$source[i],] + t(aki[net$source[i],])%*%aik[net$source[i],]
          # ai.i <- t(aik[net$source[i],])%*%aik[net$source[i],]
          # ai.i <- ai.i + t(aki[net$source[i],])%*%aki[net$source[i],]/sigma
          # ai.i <-t(aki[net$source[i],])%*%aki[net$source[i],]/sigma
          aj.j <- t(aik[net$source[j],]+aki[net$source[j],])%*%(aik[net$source[j],]+aki[net$source[j],])  #+ t(aik[net$source[j],])%*%aki[net$source[j],] + t(aki[net$source[j],])%*%aik[net$source[j],]
          # aj.j <- t(aik[net$source[j],])%*%aik[net$source[j],]
          # aj.j <- aj.j + t(aki[net$source[j],])%*%aki[net$source[j],]/sigma
          # aj.j <- t(aki[net$source[j],])%*%aki[net$source[j],]/sigma
          if (type == 'jacc'){
            net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
          } else if (type == 'cos'){
            net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
          }
        } 
        #
        # else if (net$source[i] != net$source[j] && net$target[i] != net$target[j] &&
        #            net$target[i] == net$source[j] && net$source[i] != net$target[j]){
        #   ai.j <- t(aik[net$source[i],]+aki[net$source[i],])%*%(aki[net$target[j],]+aik[net$target[j],]) #+ t(aki[net$source[i],])%*%aki[net$target[j],] + t(aik[net$source[i],])%*%aik[net$target[j],]
        #   # ai.j <- ai.j + t(aki[net$source[i],])%*%aik[net$target[j],]/sigma
        #   # ai.j <- t(aki[net$source[i],])%*%aik[net$target[j],]/sigma
        #   ai.i <- t(aki[net$source[i],]+aik[net$source[i],])%*%(aki[net$source[i],]+aik[net$source[i],]) #+ t(aik[net$source[i],])%*%aki[net$source[i],] + t(aki[net$source[i],])%*%aik[net$source[i],]
        #   # ai.i <- ai.i + t(aki[net$source[i],])%*%aki[net$source[i],]/sigma
        #   # ai.i <- t(aki[net$source[i],])%*%aki[net$source[i],]/sigma
        #   aj.j <- t(aki[net$target[j],]+aik[net$target[j],])%*%(aki[net$target[j],]+aik[net$target[j],]) #+ t(aki[net$target[j],])%*%aik[net$target[j],] +  t(aik[net$target[j],])%*%aki[net$target[j],]
        #   # aj.j <- aj.j + t(aik[net$target[j],])%*%aik[net$target[j],]/sigma
        #   # aj.j <- t(aik[net$target[j],])%*%aik[net$target[j],]/sigma
        #   if (type == 'jacc'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
        #   } else if (type == 'cos'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
        #   }
        # } else if (net$source[i] != net$source[j] && net$target[i] != net$target[j] &&
        #            net$source[i] == net$target[j] && net$target[i] != net$source[j]){
        #   ai.j <- t(aki[net$target[i],]+aik[net$target[i],])%*%(aik[net$source[j],]+aki[net$source[j],]) #+ t(aik[net$target[i],])%*%aik[net$source[j],] + + t(aki[net$target[i],])%*%aki[net$source[j],]
        #   # ai.j <- ai.j + t(aik[net$target[i],])%*%aki[net$source[j],]/sigma
        #   # ai.j <- t(aik[net$target[i],])%*%aki[net$source[j],]/sigma
        #   ai.i <- t(aki[net$target[i],]+aik[net$target[i],])%*%(aki[net$target[i],]+aik[net$target[i],]) #+  t(aki[net$target[i],])%*%aik[net$target[i],] +  t(aik[net$target[i],])%*%aki[net$target[i],]
        #   # ai.i <- ai.i + t(aik[net$target[i],])%*%aik[net$target[i],]/sigma
        #   # ai.i <- t(aik[net$target[i],])%*%aik[net$target[i],]/sigma
        #   aj.j <- t(aik[net$source[j],]+aki[net$source[j],])%*%(aik[net$source[j],]+aki[net$source[j],]) #+ t(aik[net$source[j],])%*%aki[net$source[j],] + t(aki[net$source[j],])%*%aik[net$source[j],]
        #   # aj.j <- aj.j + t(aki[net$source[j],])%*%aki[net$source[j],]/sigma
        #   # aj.j <- t(aki[net$source[j],])%*%aki[net$source[j],]/sigma
        #   if (type == 'jacc'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/(ai.i+aj.j-ai.j)
        #   } else if (type == 'cos'){
        #     net.sim[net$id[i],net$id[j]] <- ai.j/sqrt(ai.i*aj.j)
        #   }
        # }
        #
      }
    }
  }
  net.sim <- net.sim + t(net.sim)
  net.sim[which(net.sim == 0)] <- NA
  return(net.sim)
}


#### Community fitering ###
# net.plt.flt <- net.plotting
# cuttree.flt <- cutree(net.cluster, h=best.height)
# filter <- 2
# for (i in unique(net.plt.flt$commship)){
#   if (length(net.plt.flt$commship[which(net.plt.flt$commship == i)]) <= filter){
#     net.plt.flt <- net.plt.flt[which(net.plt.flt$commship != i),]
#     cuttree.flt <- cuttree.flt[which(net.plt.flt$commship != i)]
#   }
# }
# 
# ori.labels <- sort(unique(net.plt.flt$commship))
# relabel <- 1:length(unique(net.plt.flt$commship))
# 
# for (i in 1:length(relabel)) {
#   net.plt.flt[which(net.plt.flt$commship == ori.labels[i]), c('commship', 'color')] <- rep(relabel[i],length())
# }
# 
# colors <- rainbow(length(unique(net.plt.flt$commship)))
# 
# for (i in unique(net.plt.flt$commship)){
#   net.plt.flt$color[which(cuttree.flt == i)] <-  colors[i] #<- check
#   # net.plotting$color[which(cutree(net.cluster, k=k) == i)] <-  colors[i]
# 
# }
# 
# # net.network <- as.network(net.plotting[,1:2], loops = F, directed = T,
# #                           matrix.type = 'edgelist')
# net.igraph <- graph_from_data_frame(net.plt.flt[,1:2], directed = T)
# E(net.igraph)$weight <- net.plt.flt$weight
# V(net.igraph)$label <- labels[1:nodes]
# 
# net.coord <- readRDS('RDS/flndCoord.rds')
# net.coord <- net.coord[labels[1:nodes],]
# values.percentage <- list()
# 
# for (i in 1:nodes){
#   nodes.comm.percentage <- c()
#   ss <- net.plt.flt[which(net.plt.flt$source == i),]
#   sp <- net.plt.flt[which(net.plt.flt$target == i),]
#   st <- rbind(ss,sp)
#   st <- st[!duplicated(st),]
#   for (j in 1:k){
#     nodes.comm.percentage <- c(nodes.comm.percentage, length(st$color[which(st$color==colors[j])]))
#   }
#   values.percentage[i] <- list(nodes.comm.percentage)
# }
# # cairo_ps(filename = sprintf("plots/network%s.eps", filename),
# #          width =10, height = 10, pointsize = 12,
# #          fallback_resolution = 300)
# plot(net.igraph,
#      layout = as.matrix(net.coord),
#      # mode = layout_with_fr,
#      vertex.shape="pie",
#      vertex.pie.color = list(colors),
#      vertex.pie = values.percentage,
#      vertex.size = 15,
#      vertex.label.color = rgb(1,1,1),
#      edge.arrow.size=0.2,
#      edge.curved=.1,
#      edge.color = net.plotting$color,
#      vertex.label.cex=1,
#      vertex.label.dist=0,
#      vertex.label.font = 2)
# dev.off()
######## Real community ##############
## Load communities
# gt <- read.table('benchmark_W_D/community.dat')
# gt <- as.data.frame(gt)
# colnames(gt) <- c('node', 'cluster')
# gt <- gt[order(gt$cluster),]
# gt$reorder <- 1:nodes
# #Sortnet
# df <- data.frame()
# for (i in 1:nodes){
#   df <- rbind(df, net[which(net$source==gt$node[i]),])
# }
# df$order <- rep('sorted',leaves)
# for (j in 1:leaves){
#   df$source[j] <- gt$reorder[which(gt$node==df$source[j])]
#   df$target[j] <- gt$reorder[which(gt$node==df$target[j])]
# }
# net.new <- rbind(net, df)
# 
# cairo_ps(filename = sprintf("plots/realNetwork%s.eps",filename),
#          width =20, height = 10, pointsize = 12,
#          fallback_resolution = 300)
# ggplot(net.new, aes(target, source, fill=weight))+
#   facet_wrap(~order)+
#   scale_fill_viridis(name="Strength",option ="C")+
#   geom_tile(color= "white",size=0.1)+
#   scale_y_continuous(trans = "reverse")+
#   theme_minimal()
# dev.off()
# 
# cairo_ps(filename = sprintf("plots/realNetwork%s.eps",filename),
#          width =20, height = 10, pointsize = 12,
#          fallback_resolution = 300)
# ggplot(net, aes(target, source, fill=weight))+
#   scale_fill_viridis(name="Strength",option ="C")+
#   geom_tile(color= "white",size=0.1)+
#   scale_y_continuous(trans = "reverse")+
#   theme_minimal()
# dev.off()
#### Using Link-comm ####

# lc <- getLinkCommunities(na.omit(net[, 1:3]), directed = TRUE,
#                          hcmethod = "single", dirweight=0.8)
# plot(lc,  type = "graph", layout =layout.fruchterman.reingold)


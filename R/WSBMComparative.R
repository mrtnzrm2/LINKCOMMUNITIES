path <- '/Users/admin/Box/DYNAMICS/RESEARCH/FLN_52_108'
  
source(sprintf('%s/R/functions/funcLinkCommunity.R', path))
source(sprintf('%s/R/functions/loadNetworks.R', path))

save.files.path <- '/Volumes/JMA_1/FLN_52_108/plots/WSBM/'
read.files.path <- '/Users/admin/Box/DYNAMICS/RESEARCH/FLN_52_108/WSBM/community_detection/CSV'

LOGEVCON.DF <- data.frame()
LOGEVMAX.DF <- data.frame()
NMI.DF <- data.frame()
PARTITIONS.DF <- data.frame()

K <- seq(2,7)
alpha <- c(0, 25, 50, 75, 100)

for (al in alpha){
  for (k in K){
    
    if (file.exists(sprintf('%s/al_%i/LOGEV/k_%i/logev_cons.csv', read.files.path, al, k))){
      
      logev.con <- read.csv(sprintf('%s/al_%i/LOGEV/k_%i/logev_cons.csv', read.files.path, al, k), header = F) %>% 
        as.numeric()
      
      LOGEVCON.DF <- rbind(LOGEVCON.DF, data.frame(al=al/100, k=k, x=logev.con, type='CON'))
      
    }
     
    if (file.exists(sprintf('%s/al_%i/LOGEV/k_%i/logev_max.csv', read.files.path, al, k))){
      
      logev.max <- read.csv(sprintf('%s/al_%i/LOGEV/k_%i/logev_max.csv', read.files.path, al, k), header = F) %>%
        as.numeric()
      
      LOGEVMAX.DF <- rbind(LOGEVMAX.DF, data.frame(al=al/100, k=k, x=logev.max, type='MAX'))
      
    }
     
    if (file.exists(sprintf('%s/al_%i/NMI/k_%i/nmi_max_cons.csv', read.files.path, al, k))){
      
      nmi.cos.max <- read.csv(sprintf('%s/al_%i/NMI/k_%i/nmi_max_cons.csv', read.files.path, al, k), header = F) %>%
        as.numeric()
      
      NMI.DF <- rbind(NMI.DF, data.frame(al=al/100, k=k, x=nmi.cos.max, type='NMI'))
    }
      
    if (file.exists(sprintf('%s/al_%i/K/CON/k_%i.csv', read.files.path, al, k))){
      partition.con <- read.csv(sprintf('%s/al_%i/K/CON/k_%i.csv', read.files.path, al, k), header = F) %>% 
        t() %>% 
        as.numeric()
      PARTITIONS.DF <- rbind(PARTITIONS.DF, data.frame(alpha=al/100, k=k, type='CON', partition=partition.con))
    }
      
    
    if (file.exists(sprintf('%s/al_%i/K/MAX/k_%i.csv', read.files.path, al, k))){
      partition.max <- read.csv(sprintf('%s/al_%i/K/MAX/k_%i.csv', read.files.path, al, k), header = F) %>%
        t() %>% 
        as.numeric()
      PARTITIONS.DF <- rbind(PARTITIONS.DF, data.frame(alpha=al/100, k=k, type='MAX', partition=partition.max))
    }
  }
}

GENALX <- rbind(LOGEVCON.DF, LOGEVMAX.DF)
GENALX$al <- GENALX$al %>% as.factor()

MAX_GENALX <- GENALX %>%
  group_by(al, type) %>%
  summarise(k = which(x == max(x)) + 1, x=max(x))

NMI.DF$al <- NMI.DF$al %>% as.factor()

MAX_NMI <- NMI.DF %>% 
  group_by(al) %>% 
  summarise(k = which(x == max(x)) + 1, x=max(x))

p <- ggplot(GENALX, aes(color=al, shape=type))+
  geom_point(aes(as.factor(k), x), alpha=0.6, size=3)+
  geom_point(data=MAX_GENALX, aes(as.factor(k), x, color=al), size=4, shape=6)+
  geom_line(aes(k-1, x), alpha=0.6, linetype='dashed')+
  scale_color_brewer(palette ='Set1')+
  xlab('k')+
  ylab('LogEvidence')+
  ggtitle('LogEvidence vs. K')+
  theme_classic()+
ggplot(NMI.DF, aes( color=al))+
  geom_point(aes(as.factor(k), x), alpha=0.6, size=3)+
  geom_point(data=MAX_NMI, aes(as.factor(k), x, color=al), size=4, alpha=0.7, shape=6)+
  geom_line(aes(k-1,x), alpha=0.6)+
  scale_color_brewer(palette ='Set1')+
  ylab('NMI')+
  xlab('k')+
  ggtitle('NMI MAX-CON')+
  theme_classic()

cairo_ps(filename = sprintf("%s/Summaries/logev_nmi.eps", save.files.path),
         width = 18, height = 7.5,
         fallback_resolution = 200)
print(p)
dev.off()

p <- ggplot(GENALX, aes(color=as.factor(k)))+
  facet_wrap(~type)+
  geom_point(aes(al, x), alpha=0.4, size=3)+
  geom_point(data=MAX_GENALX, aes(al, x, color=as.factor(k)), size=4, alpha=0.8, shape=6)+
  geom_line(aes(as.numeric(al), x), alpha=0.6, linetype='dashed')+
  scale_color_brewer(palette ='Set1')+
  xlab('alpha')+
  ylab('LogEvidence')+
  ggtitle('LogEvidence vs. al')+
  theme_classic()+
  guides(color=guide_legend(title="k"))

cairo_ps(filename = sprintf("%s/Summaries/logev_alpha.eps", save.files.path),
         width = 12, height = 7.5,
         fallback_resolution = 200)
print(p)
dev.off()


NMI.BIG.DF <- data.frame()

# for (type.1 in c('MAX', 'CON')){
#   for (type.2 in c('MAX', 'CON')){

type.1 <- 'MAX'
type.2 <- 'CON'

for (e in 2:7){
  for (al.1 in c(0, 0.25, 0.50, 0.75, 1)){
    for (al.2 in c(0, 0.25, 0.50, 0.75, 1)){
      if (al.1 <= al.2){
        if ((al.1 == 1 || al.2 == 1) && e == 7)
          next
        part.1 <- PARTITIONS.DF %>%
          filter(alpha == al.1, type == type.1, k == e) %>%
          select(partition) %>%
          unlist() %>%
          as.numeric()
        
        part.2 <- PARTITIONS.DF %>%
          filter(alpha == al.2, type == type.2, k == e) %>%
          select(partition) %>%
          unlist() %>%
          as.numeric()
        
        var.NMI <- NMI(part.1, part.2, variant = 'joint')
        
        NMI.BIG.DF <- rbind(NMI.BIG.DF, data.frame(alpha.source=al.1, alpha.target=al.2, k=e, type=sprintf('%s-%s', type.1, type.2), nmi=var.NMI)) 
      }
    }
  }
}
#   }
# }

NMI.BIG.DF$type <- factor(NMI.BIG.DF$type, levels = c('MAX-MAX', 'CON-CON', 'MAX-CON'))

NMI.BIG.DF <- rbind(NMI.BIG.DF, data.frame(alpha.source=NMI.BIG.DF$alpha.target, 
                                           alpha.target=NMI.BIG.DF$alpha.source,
                                           k=NMI.BIG.DF$k,
                                           nmi=NMI.BIG.DF$nmi,
                                           type=NMI.BIG.DF$type))

p <- ggplot(NMI.BIG.DF, aes(as.factor(alpha.target), as.factor(alpha.source), fill=nmi))+
  facet_grid(type~k)+
  geom_tile()+
  scale_fill_viridis(option = 'B')+
  geom_text(data=NMI.BIG.DF, aes(as.factor(alpha.target), as.factor(alpha.source)), label=round(NMI.BIG.DF$nmi,3), color=ifelse(NMI.BIG.DF$nmi==1,'black','white'))+
  ylab('alpha')+
  xlab('alpha')

cairo_ps(filename = sprintf("%s/Summaries/alpha_alpha_fill_nmi_facet_type-k.eps", save.files.path),
         width = 18, height = 9.5,
         fallback_resolution = 200)
print(p)
dev.off()
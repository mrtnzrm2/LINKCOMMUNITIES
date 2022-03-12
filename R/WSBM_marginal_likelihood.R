path <- '/Users/jmarti53/Library/CloudStorage/Box-Box/DYNAMICS/RESEARCH/FLN_52_108'

source(sprintf('%s/R/functions/funcLinkCommunity.R', path))
source(sprintf('%s/R/functions/loadNetworks.R', path))

save.files.path <- '/Volumes/JMA_1/FLN_52_108/plots/WSBM/'
read.files.path <- '/Users/jmarti53/Library/CloudStorage/Box-Box/DYNAMICS/RESEARCH/FLN_52_108/WSBM/community_detection/CSV'

MARGINAL.LLH <- data.frame()

type <- c('WSBM')
nature <- c('MAX', 'CON')

for (nat in nature){
  for (ty in type){
    mll <- read.csv(sprintf('%s/Summaries/%s/%s/E_MLL_nor_2.csv', read.files.path, nat, ty), header = F) %>% 
      as.matrix()
    
    
    
    MARGINAL.LLH <- rbind(MARGINAL.LLH, data.frame(nature=nat,
                                                   type=ty,
                                                   alpha=rep(c(0.25, 0.5, 0.75), 8),
                                                   k=rep(2:9, each=3),
                                                   mll=as.vector(mll)))
  }
}

MARGINAL.LLH <- MARGINAL.LLH[which(!(MARGINAL.LLH$alpha %in% c(0, 1))),]
MARGINAL.LLH$k <- factor(MARGINAL.LLH$k, levels = 2:9)
MARGINAL.LLH$nature <- factor(MARGINAL.LLH$nature, levels = nature)
MARGINAL.LLH$alpha <- factor(MARGINAL.LLH$alpha, levels = rev(c(0, 0.25, 0.5, 0.75, 1)))

p <- ggplot(MARGINAL.LLH %>% filter(type == 'WSBM'), aes(k, alpha, fill=mll))+
  facet_wrap(~nature)+
  geom_tile()+
  scale_fill_viridis(option = 'A')+
  ggtitle('WSBM')+
  geom_text(aes(k, alpha), label=MARGINAL.LLH$mll[MARGINAL.LLH$type == 'WSBM'] %>% formatC(format = "e", digits = 2),
            color=ifelse(MARGINAL.LLH$mll[MARGINAL.LLH$type == 'WSBM']  < -1950, 'white', 'black'),
            size=2.5)

cairo_ps(filename = sprintf("%s/Summaries/alpha_all_k_fill_mll_facet_-nat.eps", save.files.path),
         width = 10, height = 5,
         fallback_resolution = 200)
print(p)
dev.off()

MARGINAL.LLH <- MARGINAL.LLH[!is.nan(MARGINAL.LLH$mll),]

MAX_MARGINAL.LLH <- MARGINAL.LLH %>%
  group_by(type, nature, alpha) %>%
  summarise(k=which(mll == max(mll))+1,
            mll=max(mll))
MAX_MARGINAL.LLH <- MAX_MARGINAL.LLH %>%
  filter(type=='WSBM')

p <- ggplot(MARGINAL.LLH %>% filter(type == 'WSBM'), aes(color=alpha, shape=nature))+
  geom_point(aes(k, mll), alpha=0.6, size=2.5)+
  geom_point(data=MAX_MARGINAL.LLH, aes(as.factor(k), mll, color=alpha), shape=6, size=4)+
  geom_line(aes(as.numeric(k), mll), alpha=0.7, linetype='dashed')+
  scale_color_brewer(palette = 'Set1')+
  ylab(TeX("$E_q\\[\\log Pr(A|\\theta, Z)\\]+RENOR$"))+
  theme_classic()

cairo_ps(filename = sprintf("%s/Summaries/k_Emll_RENOR_color_all_alpha_shape_nat.eps", save.files.path),
         width = 8, height = 5,
         fallback_resolution = 200)
print(p)
dev.off()

p <- ggplot(MARGINAL.LLH %>% filter(type == 'MY'), aes(color=alpha, shape=nature))+
  geom_point(aes(k, mll), alpha=0.6, size=2.5)+
  geom_point(data=MAX_MARGINAL.LLH, aes(as.factor(k), mll, color=alpha), shape=6, size=4)+
  geom_line(aes(as.numeric(k), mll), alpha=0.7, linetype='dashed')+
  scale_color_brewer(palette = 'Set1')+
  ylab(TeX("$\\log Pr(A|\\theta, Z)$"))+
  theme_classic()

cairo_ps(filename = sprintf("%s/Summaries/k_mll_color_all_alpha_shape_nat.eps", save.files.path),
         width = 8, height = 5,
         fallback_resolution = 200)
print(p)
dev.off()

p <- ggplot(MARGINAL.LLH %>% filter(type == 'WSBM'), aes(color=alpha, shape=nature))+
  geom_point(aes(k, mll), alpha=0.6, size=2.5)+
  geom_point(data=MAX_MARGINAL.LLH, aes(as.factor(k), mll, color=alpha), shape=6, size=4)+
  geom_line(aes(as.numeric(k), mll), alpha=0.7, linetype='dashed')+
  scale_color_brewer(palette = 'Set1')+
  ylab(TeX("$E_q\\[\\log Pr(A|\\theta, Z)\\]$"))+
  theme_classic()

cairo_ps(filename = sprintf("%s/Summaries/k_mll_color_all_alpha_shape_nat.eps", save.files.path),
         width = 8, height = 5,
         fallback_resolution = 200)
print(p)
dev.off()
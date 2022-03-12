library(magrittr)
source('functions/sformat.R')

system <- '220126'
model <- 'zz_model'
variant <- 'NONULL'
subject <- 'MAC'
version <- 'merged'
dist.ver <- 'tracto2016'
threshold <- 0.528024

format.path <- list(system = system,
                    model = model,
                    variant = variant,
                    subject = subject,
                    version = version,
                    dist.ver = dist.ver)

get.path <- function(path, folder){
  
  path <- c(path, list(folder=folder))
  file.path <- '{folder}/{subject}/{system}/{version}/{dist.ver}/{model}/{variant}'
  
  return(sformat(file.path, path))
}

fln.name <- 'FLN_107x107_GB.csv'
adj.name <- 'Adjacency_91x91_GB.csv'

imputation <- read.csv(file.path('..',get.path(format.path, 'Imputations'), fln.name)) %>% as.matrix()
adjacency <- read.csv(file.path('..',get.path(format.path, 'Imputations'), adj.name)) %>% as.matrix()
adjacency <- ifelse(adjacency < threshold, 0, 1)

adj.col <- colnames(adjacency)
fln.col <- colnames(imputation)
reorder <- match(fln.col, adj.col)
adjacency <- adjacency[reorder, reorder]

imputation <- imputation * adjacency
imputation[imputation == 0 ] <- NA

for (i in 1:nrow(imputation)){
  imputation[,i] <- imputation[,i]/sum(imputation[,i], na.rm = T)
}

source('functions/plot_matrix.R')
png(file.path('..', 'plots','merged','{dist.ver}/{model}/{variant}' %>% sformat(format.path),'imputation_w.png'), width = 6, height = 5, units = 'in', res = 200)
plot.matrix(log10(imputation)+7)
dev.off()

imputation[is.na(imputation)] <- 0
write.csv(imputation, file.path('..','CSV','imputation_{model}_{variant}_{dist.ver}.csv' %>% sformat(format.path)), row.names = F)

for (e in 1:length(fln.col)){
  if(grepl('X',fln.col[e], fixed = T)){
    col_e <- stringr::str_split(fln.col[e], 'X') %>% unlist()
    fln.col[e] <- col_e[2]
  }
}
write.csv(fln.col, file.path('..','CSV','imputation_{model}_{variant}_{dist.ver}_clabel.csv' %>% sformat(format.path)), row.names = F)

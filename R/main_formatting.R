library(ggplot2)
library(viridis)
library(dendextend)
library(network)
library(gridExtra)
library(MASS)
library(igraph)
library(DescTools)
library(pracma)
library(dils)

path <- '/Users/jmarti53/Library/CloudStorage/Box-Box/DYNAMICS/RESEARCH/LINKCOMMUNITIES/MAC/220126'
source(sprintf('%s/R/functions/funcLinkCommunity.R', path))
source(sprintf('%s/R/functions/dis_to_2d.R', path))
source(sprintf('%s/R/functions/fln_to_2d.R', path))

dis.name <- 'DistanceMatrix_Map3Dnov2021_107x107_arr.csv'
dis.matrix <- read.csv(sprintf('%s/CSV/%s', path, dis.name)) %>% as.matrix()

foldername <- 'merged'
fln.name <- 'CountMatrix_Summed_Main43_Periphery8_220126_arr_fln.csv'
fln.matrix <- read.csv(sprintf('%s/CSV/%s/%s', path, foldername, fln.name))

N <- ncol(fln.matrix)

dis.to.2d(path, dis.matrix, N)
fln.to.2d(path, fln.matrix, N)

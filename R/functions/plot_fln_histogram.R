plot.fln.histogram <- function(net, path="", foldername="", subfolder=""){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  A.n <- net
  A.n$type <- "imputation"
  A.s <- net[net$target <= 50,]
  A.s$type <- "original"
  A <- A.n %>% rbind(A.s)
  p <- ggplot2::ggplot(A, ggplot2::aes(weight, fill=type))+
    ggplot2::geom_histogram(ggplot2::aes(y=..density..), color="white", alpha=0.5, position="identity")+
    ggplot2::stat_function(fun = dnorm, args = list(mean = mean(A.s$weight), sd = sd(A.s$weight)), color="blue")+
    ggplot2::stat_function(fun = dnorm, args = list(mean = mean(A.n$weight), sd = sd(A.n$weight)), color="red")+
    ggplot2::geom_text(ggplot2::aes(1.5, 0.4), label="Skewness: %.4f" %>% sprintf(moments::skewness(A.s$weight) %>% round(4)) , color="blue")+
    ggplot2::geom_text(ggplot2::aes(1.5, 0.38), label="Skewness: %.4f" %>% sprintf(moments::skewness(A.n$weight) %>% round(4)) , color="red")+
    ggplot2::xlab("-log(fln)")+
    ggplot2::theme_classic()
  
  png(sprintf("%s/%s/%s/histograms.png", path, foldername, subfolder), width = 6, height = 5, units = 'in', res = 200)
  print(p)
  dev.off()
}
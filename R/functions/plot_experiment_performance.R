parameter.average <- function(data){
  
  nolog.par <- colnames(data)[1:2]
  log.par <- colnames(data)[3:4]
  
  data.nolog <- dplyr::tibble()
  for (par in nolog.par)
    data.nolog <- data.nolog %>% dplyr::bind_rows(dplyr::tibble(variable=par, value=data %>% dplyr::pull(par), test.rmae=data %>% dplyr::pull("test.rmae")))
  
  data.log <- dplyr::tibble()
  for (par in log.par)
    data.log <- data.log %>% dplyr::bind_rows(dplyr::tibble(variable=par, value=data %>% dplyr::pull(par) %>% log10(), test.rmae=data %>% dplyr::pull("test.rmae")))
  
  mean.log <- data.log %>% dplyr::group_by(variable) %>% dplyr::summarise(mean=mean(value), sd=sd(value))
  mean.nolog <- data.nolog %>% dplyr::group_by(variable) %>% dplyr::summarise(mean=mean(value), sd=sd(value))
  mean.log[,"mean"] <- 10^(mean.log[,"mean"])
  mean.log[,"sd"] <- 10^(mean.log[,"sd"])
  
  mean.nolog <- mean.nolog %>% dplyr::bind_rows(mean.log)
  
  return(mean.nolog)
}

plot.experiment.performance <- function(serie, inst, subfolder=""){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/Regression/XGBOOST/%s', inst$plot, inst$folder, subfolder), showWarnings = FALSE)
    dir.create(sprintf('%s/%s/Regression/XGBOOST/%s/%i', inst$plot, inst$folder, subfolder, serie), showWarnings = F)
  }
  
  data <- read.csv("../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>% sprintf(inst$folder, inst$common, serie))
  # Data analysis ----
  ## Filter data ----
  data <- data %>% dplyr::filter(Feature ==  "dist")
  ## Skim data ----
  print(skimr::skim(data))
  ## Compute XGBOOST paremeters' averages ----
  print(parameter.average(data))
  # Plotting ----
  ## Histogram ----
  parameters.names <- colnames(data)[1:2]
  data.parameters <- dplyr::tibble()
  for (par in parameters.names)
    data.parameters <- data.parameters %>% dplyr::bind_rows(dplyr::tibble(variable=par, value=data %>% dplyr::pull(par), test.rmae=data %>% dplyr::pull("test.rmae")))
  p.par.1 <- ggplot2::ggplot(data.parameters, ggplot2::aes(value))+
    ggplot2::facet_wrap(~variable, nrow = 2, ncol = 2, scales = "free")+
    ggplot2::geom_histogram(bins = 20, fill="skyblue", color="white")+
    ggplot2::geom_density(alpha=0.5, fill="blue")
  parameters.names <- colnames(data)[3:4]
  data.parameters <- dplyr::tibble()
  for (par in parameters.names){
    data.parameters <- data.parameters %>% dplyr::bind_rows(dplyr::tibble(variable=par, value=data %>% dplyr::filter(Feature ==  "dist") %>% dplyr::pull(par)))
  }
  p.par.2 <- ggplot2::ggplot(data.parameters, ggplot2::aes(value))+
    ggplot2::facet_wrap(~variable, nrow = 2, ncol = 2, scales = "free")+
    ggplot2::scale_x_log10()+
    ggplot2::geom_histogram(bins = 20, fill="skyblue", color="white")+
    ggplot2::geom_density(alpha=0.5, fill="blue")
  p.par <- ggpubr::ggarrange(p.par.1, p.par.2, nrow = 2)
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_parameters.png", inst$plot, inst$folder, subfolder, serie), width = 6, height = 5, units = 'in', res = 200)
  print(p.par)
  dev.off()
  ## Scatter plot between parameters----
  p.min_tree <- ggplot2::ggplot(data, ggplot2::aes(tree_depth, min_n))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()
  p.tree_learn.rate <- ggplot2::ggplot(data, ggplot2::aes(tree_depth, learn_rate))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_y_log10()
  p.tree_learn.reduction <- ggplot2::ggplot(data, ggplot2::aes(tree_depth, loss_reduction))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_y_log10()
  p.min_learn.rate <- ggplot2::ggplot(data, ggplot2::aes(min_n, learn_rate))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_y_log10()
  p.min_learn.reduction <- ggplot2::ggplot(data, ggplot2::aes(min_n, loss_reduction))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_y_log10()
  p.rate.reduction <- ggplot2::ggplot(data, ggplot2::aes(learn_rate, loss_reduction))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_y_log10()+
    ggplot2::scale_x_log10()
  p.par <- ggpubr::ggarrange(p.min_tree, p.tree_learn.rate, p.tree_learn.reduction, p.min_learn.rate, p.min_learn.reduction, p.rate.reduction, nrow = 2, ncol = 3)
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_scatter_pars.png", inst$plot, inst$folder, subfolder, serie), width =10, height = 5, units = 'in', res = 200)
  print(p.par)
  dev.off()
  ## Ggally between parameters and rmaes ----
  p.gally <- data[,c("min_n", "tree_depth", "loss_reduction", "learn_rate", "test.rmae", "train.rmae")] %>% dplyr::mutate(learn_rate=log10(learn_rate), loss_reduction=log10(loss_reduction)) %>% GGally::ggpairs()
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_gally.png", inst$plot, inst$folder, subfolder, serie), width = 10, height = 10, units = 'in', res = 200)
  print(p.gally)
  dev.off()
  ## Scatter test.rmae vs parameters ----
  p.min <- ggplot2::ggplot(data, ggplot2::aes(min_n, test.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()
  p.tree <- ggplot2::ggplot(data, ggplot2::aes(tree_depth, test.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()
  p.learn.reduction <- ggplot2::ggplot(data, ggplot2::aes(loss_reduction, test.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_x_log10()
  p.learn.rate <- ggplot2::ggplot(data, ggplot2::aes(learn_rate, test.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_x_log10()
  p.test <- ggpubr::ggarrange(p.min, p.tree, p.learn.reduction, p.learn.rate, nrow = 2, ncol = 2)
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_scatter_test.png", inst$plot, inst$folder, subfolder, serie), width =5, height = 5, units = 'in', res = 200)
  print(p.test)
  dev.off()
  ## Scatter train.rmae vs paramteres ----
  p.min <- ggplot2::ggplot(data, ggplot2::aes(min_n, train.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()
  p.tree <- ggplot2::ggplot(data, ggplot2::aes(tree_depth, train.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()
  p.learn.reduction <- ggplot2::ggplot(data, ggplot2::aes(loss_reduction, train.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_x_log10()
  p.learn.rate <- ggplot2::ggplot(data, ggplot2::aes(learn_rate, train.rmae))+
    ggplot2::geom_point()+
    ggplot2::geom_smooth()+
    ggplot2::scale_x_log10()
  p.train <- ggpubr::ggarrange(p.min, p.tree, p.learn.reduction, p.learn.rate, nrow = 2, ncol = 2)
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_scatter_train.png", inst$plot, inst$folder, subfolder, serie), width =5, height = 5, units = 'in', res = 200)
  print(p.train)
  dev.off()
  ## test.rmae vs train.rmae ----
  p.tt <- ggpubr::ggscatter(data, x="train.rmae", y="test.rmae", add="reg.line")+
    ggplot2::geom_point()+ 
    ggpubr::stat_cor(label.x = 0.2, label.y = 0.223) +
    ggpubr::stat_regline_equation(label.x = 0.2, label.y = 0.22)
  png(sprintf("%s/%s/Regression/XGBOOST/%s/%i/xgboost_scatter_tt.png", inst$plot, inst$folder, subfolder, serie), width =5, height = 5, units = 'in', res = 200)
  print(p.tt)
  dev.off()
}
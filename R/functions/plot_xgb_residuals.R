plot.xgb.residuals <- function(serie, inst, subfolder="", suffix=""){
  source("functions/rmae.R")
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>% sprintf(inst$plot, inst$folder, subfolder, serie), showWarnings = F)
  # Load data ----
  data <- readRDS("../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie))
  # Marginalize over ids
  data <- data %>% dplyr::group_by(id) %>% dplyr::summarise(w=unique(w), stdist=unique(dist), .pred=mean(pred), .rmse=Metrics::rmse(w, pred), .sim=mean(sim))
  data <- data %>% dplyr::mutate(.resid=(w-.pred), .stdresid = (w-.pred)/sd((w-.pred)), .rmae=abs(w-.pred)/w)
  source("functions/plot_diagnosis.R")
  # Plot diagnosis ----
  p <- plot.diagnosis(data)
  purb <- ggpubr::ggarrange(p$res.fit, p$norm.qq, p$scale.location, p$den.res, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
  png("%s/%s/Regression/XGBOOST/%s/%i/residuals_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 9, height = 8, res = 200, units = "in")
  print(purb)
  dev.off()
  # Plot RMAE ----
  p.wrmae <- ggplot2::ggplot(data, ggplot2::aes(w, .rmae))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("w")+
    ggplot2::ylab("Relative prediction error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange")+
    ggplot2::theme_bw()
  p.distrmae <- ggplot2::ggplot(data, ggplot2::aes(stdist, .rmae))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("Relative prediction error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::theme_bw()
  p.simrmae <- ggplot2::ggplot(data, ggplot2::aes(.sim, .rmae))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("Relative prediction error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::theme_bw()
  p.wrmse <- ggplot2::ggplot(data, ggplot2::aes(w, .rmse))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("w")+
    ggplot2::ylab("Root mean square error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::theme_bw()
  p.distrmse <- ggplot2::ggplot(data, ggplot2::aes(stdist, .rmse))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("Root mean square error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::theme_bw()
  p.simrmse <- ggplot2::ggplot(data, ggplot2::aes(.sim, .rmse))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("Root mean square error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::theme_bw()
  purb <- ggpubr::ggarrange(p.wrmae, p.distrmae, p.simrmae, p.wrmse, p.distrmse, p.simrmse, nrow = 2, ncol = 3, labels = c("A", "B" , "C", "D"))
  png("%s/%s/Regression/XGBOOST/%s/%i/rmae_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 9, height = 6, res = 200, units = "in")
  print(purb)
  dev.off()
  # Plot cateorize w ----
  w.25 <- data %>% dplyr::filter(w < 2.5) %>% dplyr::select(w, stdist, .sim)
  w.45 <- data %>% dplyr::filter(w > 4.5) %>% dplyr::select(w, stdist, .sim)
  w.24 <- data %>% dplyr::filter(w < 4.5, w > 2.5) %>% dplyr::select(w,  stdist, .sim)
  p.w25 <- ggplot2::ggplot(w.25, ggplot2::aes(stdist, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::ylim(c(1,2.7))+
    ggpubr::stat_cor()+
    ggplot2::ggtitle("w < 2.5")+
    ggplot2::theme_bw()
  p.w45 <- ggplot2::ggplot(w.45, ggplot2::aes(stdist, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggpubr::stat_cor(label.x = -1.4)+
    ggplot2::ggtitle("w > 4.5")+
    ggplot2::theme_bw()
  p.w24 <- ggplot2::ggplot(w.24, ggplot2::aes(stdist, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::ggtitle("w > 2.5 & w < 4.5")+
    ggplot2::ylim(c(2.3, 4.8))+
    ggpubr::stat_cor()+
    ggplot2::theme_bw()
  p.w25sim <- ggplot2::ggplot(w.25, ggplot2::aes(.sim, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::ylim(c(1,2.7))+
    ggpubr::stat_cor()+
    ggplot2::ggtitle("w < 2.5")+
    ggplot2::theme_bw()
  p.w45sim <- ggplot2::ggplot(w.45, ggplot2::aes(.sim, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggpubr::stat_cor(label.x = -1.4)+
    ggplot2::ggtitle("w > 4.5")+
    ggplot2::theme_bw()
  p.w24sim <- ggplot2::ggplot(w.24, ggplot2::aes(.sim, w))+
    ggplot2::geom_point(size=0.5)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggplot2::ggtitle("w > 2.5 & w < 4.5")+
    ggplot2::ylim(c(2.3, 4.8))+
    ggpubr::stat_cor()+
    ggplot2::theme_bw()
  purb <- ggpubr::ggarrange(p.w25, p.w24, p.w45, p.w25sim, p.w24sim, p.w45sim, nrow = 2, ncol = 3, labels = c("A", "B" , "C", "D", "E", "F"))
  png("%s/%s/Regression/XGBOOST/%s/%i/catw_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 9, height = 6, res = 200, units = "in")
  print(purb)
  dev.off()
  
  data.sub <- data %>% dplyr::select(w, stdist, .sim)
  data.sub$cat <- "w<2.5"
  data.sub$cat[data.sub$w > 2.5 & data.sub$w < 4.5] <- "w>2.5 & w<4.5"
  data.sub$cat[data.sub$w > 4.5] <- "w>4.5"
  p.1 <- ggplot2::ggplot(data.sub, ggplot2::aes(stdist, w))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggpubr::stat_cor()+
    ggplot2::theme_bw()
  p.2 <- ggplot2::ggplot(data.sub, ggplot2::aes(.sim, w))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::xlab("Standardize similarity")+
    ggplot2::ylab("w")+
    ggplot2::stat_smooth(method = "lm", color="orange", fill="orange", formula = y ~ x)+
    ggpubr::stat_cor(label.x = -1.4)+
    ggplot2::theme_bw()
  p.3 <- ggplot2::ggplot(data.sub, ggplot2::aes(stdist, .sim, color=cat))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("Standardize similarirty")+
    ggpubr::stat_cor()+
    ggplot2::theme_bw()
  purb <- ggpubr::ggarrange(p.1, p.2, nrow = 1, ncol = 2, labels = c("A", "B"))
  purb <- ggpubr::ggarrange(purb, p.3, nrow = 2, ncol = 1, labels = c(NA, "C"))
  png("%s/%s/Regression/XGBOOST/%s/%i/wfeat_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 8, height =8 , res = 200, units = "in")
  print(purb)
  dev.off()
}
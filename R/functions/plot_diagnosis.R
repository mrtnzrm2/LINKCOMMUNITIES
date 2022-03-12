require(ggplot2)
plot.diagnosis <-function(model){
  p1<-ggplot(model, aes(.pred, .resid))+
    geom_point(size=0.1, alpha=0.5)
  p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1<-p1+xlab("Fitted values")+ylab("Residuals")
  p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()
  
  shap <- shapiro.test(model$.stdresid)
  
  p2<-ggplot(model, aes(sample=.stdresid))
  p2<-p2+stat_qq() + stat_qq_line()+
    xlab("Theoretical Quantiles")+
    ylab("Standardized Residuals")+
    geom_text(aes(-2,3), label="p-value: %s" %>% sprintf(shap$p.value %>% formatC(format = "e", digits = 2)))
  p2<-p2+ggtitle("Normal Q-Q")+theme_bw()
  
  p3<-ggplot(model, aes(.pred, sqrt(abs(.stdresid))))+
    geom_point(na.rm=TRUE, size=0.1)
  p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
  p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
  p3<-p3+ggtitle("Scale-Location")+theme_bw()
  
  p4 <- ggplot(model, aes(.resid))+
    geom_density(alpha=0.5, fill="skyblue")+
    stat_function(fun=dnorm, n=101, args= list(mean(model$.resid, na.rm=T), sd=sd(model$.resid, na.rm = T)))+
    ylab("Density")+
    xlab("Residuals")+
    theme_bw()
  
  # p4<-ggplot(model, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
  # p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
  # p4<-p4+ggtitle("Cook's distance")+theme_bw()
  # 
  # p5<-ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE)
  # p5<-p5+stat_smooth(method="loess", na.rm=TRUE)
  # p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
  # p5<-p5+ggtitle("Residual vs Leverage Plot")
  # p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
  # p5<-p5+theme_bw()+theme(legend.position="bottom")
  # 
  # p6<-ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm=TRUE)+stat_smooth(method="loess", na.rm=TRUE)
  # p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
  # p6<-p6+ggtitle("Cook's dist vs Leverage hii/(1-hii)")
  # p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
  # p6<-p6+theme_bw()
  
  return(list(res.fit.plot=p1, norm.qq=p2, scale.location=p3, den.res=p4))
}
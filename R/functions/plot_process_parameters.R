plot.process.parameters <- function(path.plot, foldername, hclust.features){
  source('functions/find_height.R')
  source('functions/plot_parameters.R')
  
  Curves <- colnames(hclust.features) %>% .[2:length(.)]
  
  for (curve in Curves){
    best.height <- find.height(hclust.features, curve=curve, max.wire=1,
                               l_continuous = T, min.wire = 0, jump = 0.2)
    if (pracma::strcmp(curve, 'Dc')){
      print(hclust.features$K[hclust.features$height==best.height])
      print(hclust.features$NAC[hclust.features$height==best.height])
    }
    
    plot.parameters(hclust.features, best.height, view = curve,
                    filename = sprintf('Evolution_%s', curve), plt = T,
                    subfolder = 'VariableEvolution', path=path.plot,
                    foldername = foldername)
  } 
  
  fig <- plotly::plot_ly(data = hclust.features, x = ~K) %>%
    plotly::add_trace(y=~NEC, name="NEC ~ K", mode= "lines+markers",
                      marker = list(size = 3,
                                    color = "rgba(255, 182, 193, .9)",
                                    line = list(color = 'rgba(152, 0, 0, .8)',
                                                width = 2)))
  htmlwidgets::saveWidget(plotly::as_widget(fig), sprintf("%s/%s/%s/%s.html", path.plot, foldername, 'VariableEvolution', 'nec_k'))
}
rmae <- function(truth, pred){
  
  return(mean(abs((truth-pred)/truth)))
  
}
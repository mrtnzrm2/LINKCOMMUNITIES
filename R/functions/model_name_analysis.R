model.name <- function(type.sim, typeLinkage, mode, fln.2.w,  EC, tag){
  if (!pracma::strcmp(tag,'')){
    filename <- sprintf('%s_%s_%s_%s', type.sim, typeLinkage, mode, tag)
  } else{
    filename <- sprintf('%s_%s_%s', type.sim,  typeLinkage, mode)
  }
  
  if (EC)
    filename <- sprintf('%s_EC', filename)
  
  if (fln.2.w)
    filename <- sprintf('%s_l10', filename)
  
  return(filename)
}
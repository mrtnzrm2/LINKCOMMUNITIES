file.formatting <- function(type.sim, typeLinkage, mode, fln.2.w, EC, tag){

  if (!pracma::strcmp(tag,'')){
    filename <- sprintf('%s_%s_%s_%s',type.sim, typeLinkage, mode, tag)
    foldername <- sprintf('%s_%s_%s_%s', type.sim, typeLinkage, mode, tag) 
  } else{
    filename <- sprintf('%s_%s_%s',type.sim, typeLinkage, mode)
    foldername <- sprintf('%s_%s_%s', type.sim, typeLinkage, mode)
  }
  
  if (EC){
    filename <- sprintf('%s_EC', filename)
    foldername <- sprintf('%s_EC', foldername) 
  }
  
  if (fln.2.w){
    filename <- sprintf('%s_2w', filename)
    foldername <- sprintf('%s_2w', foldername) 
  }

  return(list(filename, foldername))
}
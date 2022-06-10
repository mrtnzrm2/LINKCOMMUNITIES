fx <- inline::cxxfunction(
  signature(x_ = "matrix"),
  "
    NumericMatrix x(x_) ;
    int nr = x.nrow(), nc = x.ncol() ;
    std::vector< std::vector<double> > vec( nc ) ;
    for( int i=0; i<nc; i++) {
        NumericMatrix::Column col = x(_,i) ;
        vec[i].assign( col.begin() , col.end() ) ;
    }
    // now do whatever with it
    // for show here is how Rcpp::wrap can wrap vector<vector<> >
    // back to R as a list of numeric vectors
    return wrap( vec ) ;
  ",
  plugin = "Rcpp"
)
historical_decomposition <- function(IRF , structural_errors, series ) {
  
  nvar = dim(IRF)[3]
  HD_definitive = array(NA, dim = c(dim(IRF)[1],K,nvar))
  
  K = min(dim(IRF)[2], dim(structural_errors)[2])
  
  for (h in 1:dim(IRF)[1]) {
    
    HD_preliminary =  matrix(NA, nrow = K, ncol = nvar)
    
    IRF_simple = IRF[h,1:K,,]
    structural_errors_simple = structural_errors[h,1:K,]
    
    IRF_reshaped = array(IRF_simple, dim = c(K,nvar^2))
    IRF_selected = IRF_reshaped
    
    for (i in 1:K) {
      for (j in 1:nvar) {
        HD_preliminary[i, j] <- IRF_selected[1:i, j + series * nvar - nvar] %*% 
          structural_errors_simple[i:1, j]
      }
    }
    HD_definitive[h,,] = HD_preliminary  
  }
  return(HD_definitive)
}




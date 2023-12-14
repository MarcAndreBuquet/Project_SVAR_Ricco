## RWZaccept_modified


RWZAccept_modified <- function (a, nvar = nvar, zero = FALSE, constrained, impulses, FEVD_check, fevd0, swish = swish) # constrained prendra les matrices spécifiées + rajouter un test de FEVDs
  
{
  if(!zero) { # Sign restrictions only
 #   qr_object <- qr(matrix(rnorm(M ^ 2, 0, 1), M, M))
 #   Q <- qr.Q(qr_object)
     qr_object <- qr(matrix(rnorm(nvar^2, 0, 1), nvar, nvar))
     Q <- qr.Q(qr_object) 
     R <- qr.R(qr_object)
  } else { # if addition of Zero restrictions (Arias, 2018)
    Q <- matrix(0, nvar, nvar)
    sign_restr_selected = constrained[[1]]
    
    for(i in seq_len(nvar)) { # Build up Q
      
      slct_row <- which(sign_restr_selected[, i] == 0) ## A modifier en f° de comment on écrit les restrictions
      R <- rbind(t(swish)[slct_row, ], Q[seq_len(i - 1), ]) # A priori, remplacer sigma_chol par t(swish)
      qr_object <- qr(t(R))
      qr_rank <- qr_object[["rank"]]
      set <- if(qr_rank == 0) {seq_len(nvar)} else {-seq_len(qr_rank)}
      N_i <- qr.Q(qr_object, complete = TRUE)[, set, drop = FALSE]
      N_stdn <- crossprod(N_i, rnorm(nvar, 0, 1))
      q_i <- N_i %*% (N_stdn / norm(N_stdn, type = "2"))
      Q[i, ] <- q_i
    }
    Q <- t(Q)
  }    
  
  ###### New code for checking the sign restrictions
  
  
  # sign_restr <- function(sigma_chol,
  #                        sign_restr, M, zero = FALSE, sign_lim = 10000) {
  #   
  #   counter <- 0
  #   sign_vec <- as.vector(sign_restr)
  #   restricted <- which(!is.na(sign_vec) & sign_vec != 0)
  #   
  #   while(TRUE) {
  #     counter <- counter + 1
  #     Q <- draw_Q(sigma_chol, sign_restr, M, zero = zero)
  #     shock <- sigma_chol %*% Q
  #     shock[abs(shock) < 1e-12] <- 0
  #     
  #     shock_vec <- as.vector(shock)
  #     shock_vec[which(shock_vec < 0)] <- -1
  #     shock_vec[which(shock_vec > 0)] <- 1
  #     
  #     if(identical(shock_vec[restricted], sign_vec[restricted])) {return(shock)}
  #     if(counter > sign_lim) {
  #       stop("No matrix fitting the sign restrictions found.")
  #     }
  #   }
  # }
  ##### End of new code for sign restrictions
  
  ## Verify sign restrictions
  
  count <- 0
  
  for (j in 1:length(constrained)) {
    
    shock = t(impulses[j,,] %*% Q) # ici variables en colonne et chocs en ligne donc OK par rapport aux spécifications
    shock[abs(shock) < 1e-8] <- 0 # changer éventuellement l'arrondi à 0 et le réaugmenter si nécessaire
    
    shock_vec <- as.vector(shock)
    shock_vec[which(shock_vec < 0)] <- -1
    shock_vec[which(shock_vec > 0)] <- 1
    
    sign_vec <- as.vector(constrained[[j]])
    
    ## Rajouter un restricted
    
    #restricted <- which(!is.na(sign_vec) & sign_vec != 0)
    restricted <- which(!is.na(sign_vec))
    
    if(identical(shock_vec[restricted], sign_vec[restricted])) {count = count + 1}
  }
  
  
  ## Code for checking FEVDs
  
  if (!is.null(FEVD_check)) {
    
    count_FEVD = 0
    
    for (i in 1:length(FEVD_check)) {
      
      fevd = t(fevd0[i, , ] %*% (Q^2))
      
      fevd_vec = as.vector(fevd)
      fevd_sign_vec = as.vector(FEVD_check[[i]])
      
      restricted <- which(!is.na(fevd_sign_vec))
      
      fevd_vec_restricted = fevd_vec[restricted]
      fevd_sign_vec_restricted = fevd_sign_vec[restricted]
      
      fevd_ref = which(fevd_sign_vec_restricted == 1)
      
      if (fevd_vec_restricted[fevd_ref] == max(fevd_vec_restricted)) {count_FEVD = count_FEVD + 1}
      }
  }
 
  ## End of code for checking FEVDs 
  
  ## Return what is needed for the remaining of the function
  
  if (!is.null(FEVD_check)) {
      if (count == length(constrained) && count_FEVD == length(FEVD_check)) {acc = 1}
      else {acc = 0}}
  else {if (count == length(constrained)) {acc = 1}
    else {acc = 0}}
  
  rwz <- list(Q = Q, acc = acc)
  return(rwz)
}

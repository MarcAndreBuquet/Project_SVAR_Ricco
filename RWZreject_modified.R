RWZreject_modified <- function (Y = NULL, nlags = 4, draws = 200, subdraws = 200, nkeep = 1000, zero = FALSE,zero_period = 2, zero_list = NULL,  # rajouter argument zero = FALSE, fevd_check = NULL (fait)
          FEVD_check = NULL, KMIN = 1, KMAX = 4, constrained = NULL, constant = TRUE, 
          steps = 24) 
{
  # sanity.check.reject(Y = Y, nlags = nlags, draws = draws, 
  #                     subdraws = subdraws, nkeep = nkeep, KMIN = KMIN, KMAX = KMAX, 
  #                     constrained = constrained, constant = constant, steps = steps) # à changer pour enlever le constrained puisque sera dégagé pour d'autres contraintes (finalement, je l'ai enlevé donc faire attention à bien rentrer les bons arguments)
  varnames <- colnames(Y)
  n1 <- draws
  n2 <- subdraws
  nstep <- steps
  nlags <- nlags
  nvar <- ncol(Y)
  nobs <- nrow(Y)
  nnobs0 <- nlags + 1
  nnobs <- nobs - nlags
  nnvar0 <- nvar + 1
  ntot <- n1 * n2
  if (constant == FALSE) {
    CONS <- "F"
    ncoef <- nvar * nlags
    nncoef <- nvar * nlags
    nnvar1 <- nvar * (nlags + 1)
  }
  else {
    CONS <- "T"
    ncoef <- nvar * (nlags + 1)
    nncoef <- nvar * nlags + 1
    nnvar1 <- nvar * (nlags + 1) + 1
  }
  ###### Rajouter une étape de check s'il y a des zéros dans les matrices de sélection
  
  # Pas nécessaire si en argument, on entre bien le fait que l'on utilise des ZR
  
  ######
  model <- rfvar(Y, lags = nlags, const = CONS) # Use of Chris Sims code for reduced-form VAR
  bcoef <- model$By
  resid <- model$u
  data <- model$X
  xx <- model$xx
  uu <- crossprod(resid)
  sigma <- (1/nnobs) * uu # empirical VCOV matrix of the X for RF model
  sxx <- chol(xx)
  sv <- solve(uu)
  svt <- chol(sv)
  betaols <- t(bcoef)
  best <- betaols
  wishdof <- nnobs - nncoef
  goodresp <- array(NA, c(nkeep, nstep, nvar, nvar)) # à modifier pour ajouter une autre dimension (à priori c(nkeep, nstep, nvar, nvar)) (fait)
  BDraws <- array(NA, c(n1, nncoef, nvar)) # a l'air OK en terme de dimensionnalité
  SDraws <- array(NA, c(n1, nvar, nvar)) # OK dimensionnalité
  imp <- array(NA, c(nstep, nvar, nvar)) # rajouter une dimension (fait) et transformer en array (fait)
  fevd <- array(NA, c(nstep, nvar, nvar)) # rajouter une dimension (fait)
  goodfevd <- array(NA, c(nkeep, nstep, nvar , nvar)) # à modifier comme précédemment (fait)
  goodshock <- array(NA, c(nkeep, nnobs, nvar))
  uhatt <- matrix(NA, nnobs, nvar) # modifié en changeant la 2ème dimension de 1 à nvar 
  HDs <- matrix(NA, nnobs, nvar) # rajouté pour initialiser les HDs
  accept <- 0
  message("Starting MCMC, ", date(), ".", sep = "")
  pb0 <- txtProgressBar(min = 0, max = n1, style = 3)
  for (draws in 1:n1) {
    setTxtProgressBar(pb0, draws)
    sigmad <- solve(matrix(rWishart(1, wishdof, sv), nrow = nvar, # Create a NIW matrix
                           ncol = nvar))
    swish <- chol(sigmad)
    swsxx <- sigmad %x% xx
    bd <- rep(0, nrow(swsxx))
    betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow = nncoef, 
                    ncol = nvar)
    betadraw <- betaols + betau
    bhat <- betadraw
    imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
    impulses <- array(imfhat, dim = c(nstep, nvar, nvar))
    imp2 <- impulses^2
    imp2sum <- apply(imp2, c(2, 3), cumsum)
    stacked_preliminary_irf = NULL
    for (i in seq_len(zero_period)) {
      stacked_preliminary_irf = rbind(stacked_preliminary_irf, impulses[i,,])
    }
    mse <- apply(imp2sum, c(1, 2), sum)
    fevd0 <- array(apply(imp2sum, 3, "/", mse), dim = c(nstep, nvar, nvar))   # Jusque-là fevd0 et impulses sont de la bonne taille pour calculer tous les chocs
    for (subdraws in 1:n2) {
      # a <- matrix(rnorm(nvar^2, mean = 0, sd = 1), nvar,  ###### à modifier pour pouvoir prendre en compte les restrictions sur plus d'une seule variable (modification faite)
      #             nvar)
      RWZA <- RWZAccept_modified(nvar = nvar, zero = zero, FEVD_check = FEVD_check, fevd0 = fevd0, constrained = constrained, impulses = impulses, swish = swish, stacked_preliminary_irf = stacked_preliminary_irf, Z_cell = zero_list)
      RWZ <- RWZA$acc
      q <- RWZA$Q
      if (RWZ == 1) {
        for (j in 1:nstep) {
          imp[j, ,] <- t(impulses[j, , ] %*% q)
          fevd[j, , ] <- t(fevd0[j, , ] %*% (q^2))
        }
        accept <- accept + 1
        goodresp[accept, , , ] <- imp     # modifié en correspondance avec le chgt définitionnel précédent
        goodfevd[accept, , , ] <- fevd * 100 # pareil (rajout d'une virgule)
        BDraws[draws, , ] <- betadraw
        SDraws[draws, , ] <- sigmad
        uhat <- Y[nnobs0:nobs, ] - data %*% bhat # OK
        
        ## enlever car inutile et rajoute du temps de calcul

         for (i in 1:nnobs) {
           uhatt[i, ] <- uhat[i, ] %*% (solve(swish) %*%
                                          q)
         }
         goodshock[accept,, ] <- t(uhatt)
         
      }
      else {
        next
      }
      if (accept >= nkeep) {
        break
      }
    }
    if (accept >= nkeep) {
      break
    }
    ldraw <- draws
  }
  close(pb0)
  if (ldraw < n1) {
    BDraws <- BDraws[1:ldraw, , ]
    SDraws <- SDraws[1:ldraw, , ]
    dimnames(SDraws) <- list(1:ldraw, varnames, varnames)
  }
  if (accept < nkeep) {
    if (accept == 0) {
      stop("\n Not enough accepted draws to proceed!")
    }
    else {
      goodresp <- goodresp[1:accept, , ,] # à modifier (fait)
      goodfevd <- goodfevd[1:accept, , ,] # à modifier (fait)
      goodshock <- goodshock[1:accept, ,]
      message("\n Warning! Had only ", accept, " accepted draw(s) out of ", 
              ntot, ".", sep = "")
    }
  }
  nn1 <- accept
  #  dimnames(goodresp) <- list(1:nn1, 1:nstep, varnames)
  #  dimnames(goodfevd) <- list(1:nn1, 1:nstep, varnames)
  if (constant == FALSE) {
    dimnames(BDraws) <- list(1:ldraw, c(paste(varnames, rep(1:nlags, 
                                                            each = length(varnames)), sep = "")), varnames)
  }
  else {
    dimnames(BDraws) <- list(1:ldraw, c(paste(varnames, rep(1:nlags, 
                                                            each = length(varnames)), sep = ""), "const"), 
                             varnames)
  }
  message("\n MCMC finished, ", date(), ".", sep = "")
  return(list(IRFS = goodresp, FEVDS = goodfevd, SHOCKS = goodshock, 
              BDraws = BDraws, SDraws = SDraws))
}

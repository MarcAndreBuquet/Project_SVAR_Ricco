fevdplot <-
  function(fevddraws=NULL, type="median", labels=unlist(dimnames(fevddraws)[3]), save=FALSE, bands=c(0.16, 0.84), grid=TRUE, bw=FALSE, table=FALSE, periods=NULL){
    
    #--- SANITY CHECK ---#
    # sanity.check.fevdplot(fevddraws=fevddraws, type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw, table=table)
    
    
    graphics.off()
    par.def <- par(no.readonly = T)
    #graph parameter
    goodresp <- fevddraws
    irftype <- type #  0== median, 1== mean response
    gridgrph <- grid # grid in irf plots 0== none, 1== adds grid to plots
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      ebupp <- bands[2]# error bands for irf plots
      eblow <- bands[1]
    }else{
      ebupp <- 0.84# error bands for irf plots
      eblow <- 0.16
    }
    varlbl <- labels
    nstep <- dim(fevddraws)[2]
    nvar <- dim(fevddraws)[3]
    periodst <- is.null(periods)
    
    if(irftype=="mean"){
      imp_responses <- array(NA, dim=c(3, nstep, nvar))
      irfbands <- apply(goodresp,c(2,3),quantile,probs=c(eblow, ebupp))
      irfmean <-  array(apply(goodresp,c(2,3),mean), dim=c(1,nstep, nvar))
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
      imp_responses[1,,] <- irfmean
      imp_responses[2:3,,] <- irfbands
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
    }else{
      imp_responses <- apply(goodresp,c(2,3),quantile,probs=c(0.5, eblow, ebupp))
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
    }
    impt <- imp_responses
    impt <- aperm(impt,c(3,2,1))
    
    if(table==FALSE){
      #--- DETERMINE COLS AND ROWS OF PLOT ---#
      rowsize <-  ceiling(sqrt(nvar))
      colsize <- ceiling(nvar / rowsize)
      
      #-- GENERATE PLOTS ---#
      # dev.off()
      par(bty="o", mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
      if(bw==FALSE){
        for(i in 1:nvar){
          ulim <- max(impt[i,,1:3])
          llim <- min(impt[i,,1:3])
          plot(x=1:nstep, y=impt[i,,1], type="l", col="red", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
          if(bndtest!=TRUE){
            lines(1:nstep, impt[i,,2], col="blue", lwd=2)
            lines(1:nstep, impt[i,,3], col="blue", lwd=2)
          }
          abline(h=0, col="black")
          if(gridgrph==1){
            grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
          }
        }
        if(save==TRUE){
          dev.copy(postscript,'fevd.eps')
          dev.off()
        }
      }else{
        for(i in 1:nvar){
          ulim <- max(impt[i,,1:3])
          llim <- min(impt[i,,1:3])
          plot(x=1:nstep, y=impt[i,,1], type="l", col="black", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
          if(bndtest!=TRUE){
            lines(1:nstep, impt[i,,2], lty=2, col="black", lwd=2)
            lines(1:nstep, impt[i,,3], lty=2, col="black", lwd=2)
          }
          abline(h=0, col="black")
          if(gridgrph==1){
            grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
          }
        }
        if(save==TRUE){
          dev.copy(postscript,'fevd.eps')
          dev.off()
        }
        
      }
      par(par.def)
    }else{#if table ==TRUE
      if(periodst==1){# if period ==null
        fevdtable <- t(impt[,,1])
        fevdtable <- round(fevdtable, 2)
        colnames(fevdtable) <- varlbl
      }else{# if period !=null
        fevdtable <- t(impt[,periods,1])
        fevdtable <- round(fevdtable, 2)
        colnames(fevdtable) <- varlbl
      }
      #returns table
      return(fevdtable)
    }
  }

fn.impulse <-
  function(Bh,swish,nn){
    nvar <- nn[1]
    lags <- nn[2]
    imstep <- nn[3]
    #
    ll <- lags + 1
    n1 <- nvar + 1
    nl <- nvar*lags
    #
    Ah <- t(Bh)
    #
    imf <- matrix(nrow=imstep, ncol=nvar*nvar)
    #
    M <- matrix(nrow=nvar*imstep, ncol=nvar)
    #
    M[1:nvar,] <- t(swish)
    Mtem <- M[1:nvar,]
    #
    imf[1,] <- t(as.vector(Mtem))
    #
    ims2 <- imstep - 1
    ims1 <- min(c(ims2, lags))
    t <- 1
    while(t <=ims1){
      nt <- nvar*t
      ntt <- nvar*(t+1)
      tt <- t+1
      Mtem <- Ah[,1:nt] %*% M[1:nt,]
      M[n1:ntt,] <- M[1:nt,]
      M[1:nvar,] <- Mtem
      imf[tt,] <- t(as.vector(Mtem))
      t <- t+1
    }
    #
    for(t in ll:ims2){
      nt <- nvar*t
      ntt <- nvar*(t+1)
      tt <- t+1
      Mtem <- Ah[,1:nl] %*% M[1:nl,]
      M[n1:ntt,] <- M[1:nt,]
      M[1:nvar,] <- Mtem
      imf[tt,] <- t(as.vector(Mtem))
    }
    return(imf)
  }

fp.target <-
  function(Y=NULL, irfdraws=NULL, nlags=4, constant=TRUE, type="median", labels=colnames(Y), target= TRUE, save=FALSE, legend=TRUE, bands=c(0.16, 0.84), grid=TRUE, bw=FALSE, maxit=1000){ # what else do i need in terms of options?
    #
    #--- SANITY CHECKS ---#
    sanity.check.target(Y=Y, nlags=nlags, irfdraws=irfdraws, constant=constant, type=type, labels=labels, target= target, save=save, legend=legend, bands=bands, grid=grid, bw=bw, maxit=maxit)
    #
    graphics.off()
    par.def <- par(no.readonly = T)
    #--- PARAS---#
    ldg <- legend
    irftype <- type
    goodresp <- irfdraws
    gridgrph <- grid # grid in irf plots 0== none, 1== adds grid to plots
    tgt <- target
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      ebupp <- bands[2]# error bands for irf plots
      eblow <- bands[1]
    }else{
      ebupp <- 0.84# error bands for irf plots
      eblow <- 0.16
    }
    varlbl <- labels
    nstep <- dim(irfdraws)[2]
    nvar <- ncol(Y)
    nobs <- nrow(Y)
    n1 <- maxit
    nlags <- nlags
    nnobs0 <- nlags + 1
    nnobs <- nobs - nlags
    nnvar0 <- nvar + 1
    #
    if(constant == FALSE){
      CONS <- "F"
      ncoef <- nvar * nlags
      nncoef <- nvar * nlags
      nnvar1 <- nvar * (nlags + 1)
    }else{
      CONS <- "T"
      ncoef <- nvar * (nlags+1)
      nncoef <- nvar * nlags + 1
      nnvar1 <- nvar * (nlags + 1) + 1
    }
    #
    #--- SET UP MATRICES ---#
    targets <- matrix(0.0, nrow=nstep, ncol=nvar)
    rescale <- matrix(0.0,nrow=nstep, ncol=nvar)
    lower  <- matrix(0.0, nrow=nstep, ncol=nvar)
    upper  <- matrix(0.0, nrow=nstep, ncol=nvar)
    #
    #--- GET TARGETS ---#
    for(i in 1:nvar){# that can be done using apply
      for(k in 1:nstep){
        if(irftype=="mean"){
          targets[k,i] <- mean(goodresp[ ,k,i])
        }else{
          targets[k,i] <- quantile(goodresp[ ,k,i], probs=0.5)
        }
        rescale[k,i] <- 1.0/(quantile(goodresp[ ,k,i], probs=0.75) - quantile(goodresp[ ,k,i], probs=0.25))
        lower[k,i]  <- quantile(goodresp[ ,k,i], probs= eblow)
        upper[k,i]  <- quantile(goodresp[ ,k,i], probs= ebupp)
      }
    }
    #
    #--- RFVAR ---#
    model <- rfvar(Y,lags=nlags, const=CONS, breaks=NULL)
    bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
    resid <- model$u # same as above
    data <- model$X
    xx <- model$xx
    #
    #--- SIGMA and SXX ---#
    uu <- crossprod(resid)
    # sigma <- (1/(nnobs-nncoef))*uu
    sigma <- (1/nnobs)*uu
    betaols <- t(bcoef)
    #
    swish <- chol(sigma)
    imfhat <- fn.impulse(betaols, swish, c(nvar, nlags, nstep))
    impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
    #
    #--- PAGAN/FRY GAP ---#
    g <- matrix(1, nrow=nvar-1, ncol=1)
    gapfunc <- minqa::uobyqa(g, fn=FPgap, control = list(maxfun=n1), rescale=rescale, targets=targets, impulses=impulses, nstep= nstep)
    #
    paras <- gapfunc$par
    funcval <- gapfunc$fval
    conv <- gapfunc$ierr
    #
    #--- WARNING MESSAGE ---#
    if(conv>0){
      message('\n Warning! optimisation did not converge.', sep="")
    }
    #
    #--- CREATE MIN IMPULSES ---#
    imp <- matrix(NA, nrow=nstep, ncol=nvar)
    q <- stereo(paras)
    imfhat <- fn.impulse(betaols, swish, c(nvar, nlags, nstep))
    impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
    for(j in 1:nstep){
      imp[j,] <- t(impulses[j,,]%*%q)
    }
    #
    #--- DETERMINE COLS AND ROWS OF PLOT ---#
    rowsize <-  ceiling(sqrt(nvar))
    colsize <- ceiling(nvar / rowsize)
    #
    #--- GRAPH PARAS ---#
    if(ldg==TRUE){
      par(bty="o", oma = c(4, 1, 1, 1), mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
    }else{
      par(bty="o",  mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
    }
    #--- GENERATE PLOTS ---#
    if(bw==FALSE){
      for(i in 1:nvar){
        ulim <- max(c(imp[,i],upper[,i],lower[,i]))
        llim <- min(c(imp[,i],upper[,i],lower[,i]))
        plot(x=1:nstep, y=imp[,i], type="l", col="red", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
          lines(1:nstep, y=lower[,i], col="blue", lwd=2)
          lines(1:nstep,  y=upper[,i], col="blue", lwd=2)
        }
        abline(h=0, col="black")
        if(target==TRUE){
          lines(1:nstep, y=targets[,i], lty=2, col="blue", lwd=2)
        }
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
      }
      if(ldg==TRUE){
        if(target==TRUE){
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          if(bndtest==TRUE){
            legend("bottom",  legend=c("Fry-Pagan", "Impulse response"), col = c("red", "blue"), lty=c(1, 2), lwd=c(2, 2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }else{
            legend("bottom",  legend=c("Fry-Pagan", "Impulse response", "Error bands"), col = c("red", "blue", "blue"), lty=c(1, 2, 1), lwd=c(2, 2, 2),horiz=TRUE   , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5 , x.intersp=1, text.width=c(0.25,0.25,0.5))
          }
        }else{
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          if(bndtest==TRUE){
            legend("bottom",  legend=c("Fry-Pagan"), col = c("red"), lty=c(1), lwd=c(2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }else{
            legend("bottom",  legend=c("Fry-Pagan",  "Error bands "), col = c("red", "blue"), lty=c(1, 1), lwd=c(2, 2 ), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }
        }
      }
      if(save==TRUE){
        dev.copy(postscript,'fptarget.eps')
        dev.off()
      }
    }else{#start bw
      for(i in 1:nvar){
        ulim <- max(c(imp[,i],upper[,i],lower[,i]))
        llim <- min(c(imp[,i],upper[,i],lower[,i]))
        plot(x=1:nstep, y=imp[,i], type="l", col="black", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
          lines(1:nstep, y=lower[,i],  lty=2, col="black", lwd=2)
          lines(1:nstep,  y=upper[,i],  lty=2, col="black", lwd=2)
        }
        abline(h=0, col="black")
        if(target==TRUE){
          lines(1:nstep, y=targets[,i], col="darkgrey", lwd=2)
        }
        
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
      }
      if(ldg==TRUE){
        if(target==TRUE){
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          if(bndtest==TRUE){
            legend("bottom",  legend=c("Fry-Pagan", "Impulse response"), col = c("black", "darkgrey"), lty=c(1, 1), lwd=c(2, 2), horiz=TRUE   , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }else{
            legend("bottom",  legend=c("Fry-Pagan", "Impulse response", "Error bands"), col = c("black", "darkgrey", "black"), lty=c(1, 1, 3), lwd=c(2, 2, 2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5, x.intersp=1, text.width=c(0.25,0.25,0.5))
          }
        }else{
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          if(bndtest==TRUE){
            legend("bottom",  legend=c("Fry-Pagan"), col = c("black"), lty=c(1), lwd=c(2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }else{
            legend("bottom",  legend=c("Fry-Pagan",  "Error bands "), col = c("black", "black"), lty=c(1, 2), lwd=c(2, 2 ), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
          }
        }
      }
      if(save==TRUE){
        dev.copy(postscript,'fptarget.eps')
        dev.off()
      }
    }
    par(par.def)
  }

FPgap <-
  function(g, rescale, targets, impulses, nstep){
    gap <- 0.0
    a <- matrix(stereo(v=g))
    for(k in 1:nstep){
      scal  <-  matrix(rescale[k,] * (targets[k,] - impulses[k, , ] %*% a))
      scal <- sum(scal^2)
      gap <- gap + scal
    }
    return(gap)
  }

irfplot <-
  function(irfdraws=NULL,type="median", labels=unlist(dimnames(irfdraws)[3]), save=FALSE, bands=c(0.16, 0.84), grid=TRUE, bw=FALSE){
    #
    
    #--- SANITY CHECK ---#
    sanity.check.irfplot(irfdraws=irfdraws,type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw)
    
    graphics.off()
    par.def <- par(no.readonly = T)
    #graph parameter
    goodresp <- irfdraws
    irftype <- type #  0== median, 1== mean response
    gridgrph <- grid # grid in irf plots 0== none, 1== adds grid to plots
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      ebupp <- bands[2]# error bands for irf plots
      eblow <- bands[1]
    }else{
      ebupp <- 0.84# error bands for irf plots
      eblow <- 0.16
    }
    varlbl <- labels
    nstep <- dim(irfdraws)[2]
    nvar <- dim(irfdraws)[3]
    
    if(irftype=="mean"){
      imp_responses <- array(NA, dim=c(3, nstep, nvar))
      irfbands <- apply(goodresp,c(2,3),quantile,probs=c(eblow, ebupp))
      irfmean <-  array(apply(goodresp,c(2,3),mean), dim=c(1,nstep, nvar))
      dimnames(imp_responses) <- list(c("irf", "lower", "upper"),1:nstep, varlbl)
      imp_responses[1,,] <- irfmean
      imp_responses[2:3,,] <- irfbands
      dimnames(imp_responses) <- list(c("irf", "lower", "upper"),1:nstep, varlbl)
    }else{
      imp_responses <- apply(goodresp,c(2,3),quantile,probs=c(0.5, eblow, ebupp))
      dimnames(imp_responses) <- list(c("irf", "lower", "upper"),1:nstep, varlbl)
    }
    impt <- imp_responses
    impt <- aperm(impt,c(3,2,1))
    
    #--- DETERMINE COLS AND ROWS OF PLOT ---#
    rowsize <-  ceiling(sqrt(nvar))
    colsize <- ceiling(nvar / rowsize)
    
    #-- GENERATE PLOTS ---#
    # dev.off()
    par(bty="o", mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
    if(bw==FALSE){
      for(i in 1:nvar){
        ulim <- max(impt[i,,1:3])
        llim <- min(impt[i,,1:3])
        plot(x=1:nstep, y=impt[i,,1], type="l", col="red", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
          lines(1:nstep, impt[i,,2], col="blue", lwd=2)
          lines(1:nstep, impt[i,,3], col="blue", lwd=2)
        }
        abline(h=0, col="black")
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
      }
      if(save==TRUE){
        dev.copy(postscript,'irf.eps')
        dev.off()
      }
    }else{
      for(i in 1:nvar){
        ulim <- max(impt[i,,1:3])
        llim <- min(impt[i,,1:3])
        plot(x=1:nstep, y=impt[i,,1], type="l", col="black", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
          lines(1:nstep, impt[i,,2], lty=2, col="black", lwd=2)
          lines(1:nstep, impt[i,,3], lty=2, col="black", lwd=2)
        }
        abline(h=0, col="black")
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
        
      }
      if(save==TRUE){
        dev.copy(postscript,'irf.eps')
        dev.off()
      }
      
    }
    par(par.def)
  }

rfbvar <-
  function(Y=NULL, nlags=4, draws=1000, constant=TRUE, steps=24, shock=1){
    #
    #--- SANITY CHECK ---#
    sanity.check.rfbvar(Y=Y, nlags=nlags, draws=draws, constant=constant, steps=steps, shock=shock)
    
    #--- SET UP PARAS ---#
    varnames <- colnames(Y)
    n1 <- draws
    nstep <- steps
    nlags <- nlags
    shock <- shock
    nvar <- ncol(Y)
    nobs <- nrow(Y)
    nnobs0 <- nlags + 1
    nnobs <- nobs - nlags
    nnvar0 <- nvar + 1
    #
    if(constant == FALSE){
      CONS <- "F"
      ncoef <- nvar * nlags
      nncoef <- nvar * nlags
      nnvar1 <- nvar * (nlags + 1)
    }else{
      CONS <- "T"
      ncoef <- nvar * (nlags+1)
      nncoef <- nvar * nlags + 1
      nnvar1 <- nvar * (nlags + 1) + 1
    }
    #
    #---REDUCED FORM VAR MODEL ---#
    model <- rfvar(Y,lags=nlags, const=CONS)
    bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
    resid <- model$u # same as above
    data <- model$X
    xx <- model$xx
    
    #--- SIGMA and SXX ---#
    uu <- crossprod(resid)
    # sigma <- (1/(nnobs-nncoef))*uu
    sigma <- (1/nnobs)*uu
    
    #--- SET UP MCMC OF VAR ---#
    sxx <-  chol(xx)
    sv <- solve(uu)
    svt <-  chol(sv)
    betaols <- t(bcoef)
    best <- betaols
    wishdof <- nnobs-nncoef
    
    #--- MATRICES FOR DRAWS ---#
    goodresp <- array(NA, c(n1, nstep, nvar))
    BDraws <- array(NA, c(n1, nncoef, nvar))
    SDraws <- array(NA, c(n1, nvar, nvar))
    imp <- matrix(NA, nrow=nstep, ncol=nvar)
    fevd <- matrix(NA, nrow=nstep, ncol=nvar)
    goodfevd <- array(NA, c(n1, nstep, nvar))
    goodshock <- array(NA, c(n1, nnobs))
    
    #--- Monte CARLO INTEGRATION ---#
    message('\n Starting MCMC, ', date(),'.', sep="")
    pb <- txtProgressBar(min = 0, max = n1, style = 3)
    for(draws in 1:n1){
      setTxtProgressBar(pb, draws)
      
      #--- sigma draws ---#
      sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
      swish   <- chol(sigmad)
      
      #--- beta draws ---#
      swsxx <- sigmad  %x% xx
      bd <- rep(0, nrow(swsxx))
      #betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
      betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
      betadraw <- betaols + betau
      bhat <- betadraw
      
      #--- irfs ---#
      imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
      impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
      imp2 <- impulses^2
      
      imp2sum <- apply(imp2, c(2,3), cumsum)
      mse <-  apply(imp2sum, c(1,2), sum)
      fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))
      
      imp <- impulses[,,shock]
      fevd <- fevd0[,,shock]
      
      goodresp[draws, ,] <- imp
      goodfevd[draws, ,] <- fevd * 100
      BDraws[draws, , ] <- betadraw
      SDraws[draws, , ] <- sigmad
      uhat <-   Y[nnobs0:nobs ,] - data %*%bhat
      uhatt <- t(uhat %*% solve(swish))
      uhatt <-  uhatt[ abs(shock),]
      goodshock[draws, ] <-  uhatt
    }#end draws
    close(pb)
    #
    dimnames(goodresp) <- list(1:n1, 1:nstep, varnames)
    dimnames(goodfevd) <- list(1:n1, 1:nstep, varnames)
    dimnames(SDraws) <- list(1:n1, varnames, varnames)
    #
    if(constant == FALSE){
      dimnames(BDraws) <-  list(1:n1, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
        dimnames(BDraws) <- list(1:n1, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
      }
    #
    message('\n MCMC finished, ', date(),'.', sep="")
    return(list(IRFS=goodresp, FEVDS =goodfevd, BDraws=BDraws, SDraws=SDraws, SHOCKS = goodshock))
  }

rfvar <-
  function(ydata=NA, lags=2, xdata=NULL, const=TRUE, breaks=NULL){
    if(is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
    T <- dim(ydata)[1]
    nvar <- dim(ydata)[2]
    if(const){
      xdata <- cbind(xdata,matrix(1,T,1))
    }
    nox <- identical(xdata,NULL)
    if(!nox){
      T2 <- dim(xdata)[1]
      nx <- dim(xdata)[2]
    }
    else{
      T2 <- T
      nx <- 0
      xdata <- matrix(0,T2,0)
    }
    #
    if(!identical(T2,T)){
      print('xdata and ydata must be same length.')
      return()
    }
    if(identical(breaks,NULL))
      nbreaks <- 0
    else
      nbreaks<-length(breaks)
    breaks <- c(0,breaks,T)
    if(any(breaks[2:length(breaks)] <= breaks[1:(length(breaks)-1)]))
      stop("list of breaks must be in strictly increasing order\n")
    if(breaks[2]>lags)
      smpl <- (lags+1):breaks[2]
    else
      smpl <- NULL
    if(nbreaks>0){
      for (nb in 2:(nbreaks+1))
        smpl <- c(smpl,(breaks[nb]+lags+1):breaks[nb+1])
    }
    Tsmpl <- length(smpl)
    X <- array(0,dim=c(Tsmpl,nvar,lags))
    for(ix in seq(along=smpl))
      X[ix,,] <- t(ydata[smpl[ix]-(1:lags),,drop=FALSE])
    dim(X) <- c(Tsmpl,nvar*lags)
    X <- cbind(X, xdata[smpl,,drop=FALSE]) # rhs data for model
    y <- ydata[smpl,,drop=FALSE] # lhs data for model
    vldvr <- svd(X)
    di <- 1./vldvr$d
    dfx <- sum(vldvr$d > 100*.Machine$double.eps)
    di <- di[1:dfx]
    vldvr$u <- vldvr$u[, 1:dfx]
    vldvr$v <- vldvr$v[, 1:dfx]
    snglty <- dim(X)[2] - dfx
    B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
    u <-  y-X %*% B
    if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,lags),frequency=frequency(ydata))
    nX <- dim(X)[2]
    xx <-  di * t(vldvr$v)
    xx <-  crossprod(xx)
    By <-  t(B) # fix this only plus 1 for const
    # dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
    # By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdt.m
    if(!is.null(dimnames(ydata)[2]))
    {
      ynames <- dimnames(ydata)[[2]]
    }else
    {
      ynames <- rep("",times=nvar)
    }
    if(!nox)
    {
      if(!is.null(dimnames(xdata)[2]))
      {
        xnames <- dimnames(xdata)[[2]]
        xxnames <- c(paste(rep(ynames,each=lags),1:lags,sep=""), xnames)
        # dimnames(xx) <- list(xxnames,xxnames)
        dimnames(By) <- list(ynames, c(paste(ynames,rep(1:lags, each=length(ynames)), sep="")))
      }else
      {
        xnames <- rep("",times=nx)
        xxnames <- c(paste(rep(ynames,each=lags),1:lags,sep=""))
        # dimnames(xx) <- list(xxnames,xxnames)
        dimnames(By) <- list(ynames, c(paste(ynames,rep(1:lags, each=length(ynames)), sep=""),"const"))
      }
    }
    if (nox)
      Bx <-  NULL
    else
    {
      Bx <-  matrix(B[nvar*lags+(1:nx),],dim(B)[2],nx)
      dimnames(Bx) <- list(ynames,xnames)
    }
    return(list(By=By, Bx=Bx, u=u, xx = xx, singular=snglty, X=X))
  }

rwz.reject <-
  function(Y=NULL, nlags=4, draws=200, subdraws=200, nkeep=1000, KMIN=1, KMAX=4, constrained=NULL, constant=TRUE, steps=24){
    #
    #--- SANITY CHECK ---#
    sanity.check.reject(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps)
    #
    #--- SET UP PARAS ---#
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
    ntot <- n1*n2
    #
    if(constant == FALSE){
      CONS <- "F"
      ncoef <- nvar * nlags
      nncoef <- nvar * nlags
      nnvar1 <- nvar * (nlags + 1)
    }else{
      CONS <- "T"
      ncoef <- nvar * (nlags+1)
      nncoef <- nvar * nlags + 1
      nnvar1 <- nvar * (nlags + 1) + 1
    }
    #
    #---REDUCED FORM VAR MODEL ---#
    model <- rfvar(Y,lags=nlags, const=CONS)
    bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
    resid <- model$u # same as above
    data <- model$X
    xx <- model$xx
    #
    #--- SIGMA and SXX ---#
    uu <- crossprod(resid)
    # sigma <- (1/(nnobs-nncoef))*uu
    sigma <- (1/nnobs)*uu
    #
    #--- SET UP MCMC OF VAR ---#
    sxx <-  chol(xx)
    sv <- solve(uu)
    svt <-  chol(sv)
    betaols <- t(bcoef)
    best <- betaols
    wishdof <- nnobs-nncoef
    #
    #--- MATRICES FOR DRAWS ---#
    goodresp <- array(NA, c(nkeep, nstep, nvar))
    BDraws <- array(NA, c(n1, nncoef, nvar))
    SDraws <- array(NA, c(n1, nvar, nvar))
    imp <- matrix(NA, nrow=nstep, ncol=nvar)
    fevd <- matrix(NA, nrow=nstep, ncol=nvar)
    goodfevd <- array(NA, c(nkeep, nstep, nvar))
    goodshock <- array(NA, c(nkeep, nnobs))
    uhatt <- matrix(NA, nnobs, 1)
    #
    #--- MCMC INTEGRATION ---#
    accept <- 0
    message('Starting MCMC, ', date(),'.', sep="")
    pb0 <- txtProgressBar(min = 0, max = n1, style = 3)
    for(draws in 1:n1){
      setTxtProgressBar(pb0, draws)
      #
      #--- sigma draws ---#
      sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
      swish   <- chol(sigmad)
      #
      #--- beta draws ---#
      swsxx <-   sigmad  %x% xx
      bd <- rep(0, nrow(swsxx))
      # betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
      betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
      betadraw <- betaols + betau
      bhat <- betadraw
      #
      #--- irfs ---#
      imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
      impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
      imp2 <- impulses^2
      imp2sum <- apply(imp2, c(2,3), cumsum)
      mse <-  apply(imp2sum, c(1,2), sum)
      fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))
      #
      for(subdraws in 1:n2){
        a <- matrix(rnorm(nvar,mean=0,sd=1), nvar, 1)
        RWZA <- RWZAccept(a,KMIN,KMAX, constrained, impulses)
        RWZ <- RWZA$acc
        q <- RWZA$Q
        #
        if(RWZ==1){
          for(j in 1:nstep){
            imp[j,] <- t(impulses[j,,]%*%q)
            fevd[j,] <- t(fevd0[j,,]%*%(q^2))
          }
          accept <- accept+1
          goodresp[accept, ,] <-  imp
          goodfevd[accept, ,] <- fevd * 100
          BDraws[draws, , ] <- betadraw
          SDraws[draws, , ] <- sigmad
          uhat <-   Y[nnobs0:nobs ,] - data %*% bhat
          for(i in 1:nnobs){
            uhatt[i,] <-   uhat[i, ] %*%  (  solve(swish) %*% q)
          }
          goodshock[accept, ] <-  t(uhatt)
        }else{
          next
        }
        if(accept>=nkeep){
          break
        }
      } # end subdraws
      #
      if(accept>=nkeep){
        break
      }
      ldraw <- draws
    }#END DRAWS
    close(pb0)
    #
    #--- FIX PARA MATRICES ---#
    if(ldraw<n1){
      BDraws <- BDraws[1:ldraw, , ]
      SDraws <- SDraws[1:ldraw, , ]
      dimnames(SDraws) <- list(1:ldraw, varnames, varnames)
    }
    #
    #--- WARNING MESSAGE IN CASE OF TOO FEW DRAWS ---#
    if(accept<nkeep){#this can be shortened
      if(accept==0){
        stop("\n Not enough accepted draws to proceed!")
      }else{
        goodresp <- goodresp[1:accept, , ]
        goodfevd <- goodfevd[1:accept, , ]
        goodshock <- goodshock[1:accept, ]
        message('\n Warning! Had only ', accept,' accepted draw(s) out of ',ntot,'.', sep="")
      }
    }
    nn1 <- accept
    dimnames(goodresp) <- list(1:nn1, 1:nstep, varnames)
    dimnames(goodfevd) <- list(1:nn1, 1:nstep, varnames)
    #
    if(constant == FALSE){
      dimnames(BDraws) <-  list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
        dimnames(BDraws) <- list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
      }
    #
    message('\n MCMC finished, ', date(),'.', sep="")
    return(list(IRFS=goodresp, FEVDS = goodfevd,  SHOCKS = goodshock, BDraws=BDraws, SDraws=SDraws))
  }

RWZAccept <-
  function(a, first, last, constrained, impulses){
    QR <- qr(a)
    R <- qr.R(QR)
    Q <- qr.Q(QR)
    #
    if(R<0){Q  <- -1.0 * Q}
    #
    for(k in first:last){
      ik <- impulses[k, , ]%*%Q
      for(i in 1:length(constrained)){
        if(constrained[i]<0){
          value <- ik[-1.0*constrained[i]]
        }else{
          value <- -1.0 * ik[constrained[i]]
        }
        #
        if(value<0.0){
          if(k==first & i==1){
            Q <- 1.0 * Q
            ik <- 1.0 * ik
          } # comment this one out and uncomment bracket below.
        }else{
          acc <- 0
          rwz <- list(Q=Q, acc=acc, ika=ik)
          return(rwz)
          # } # comment out this
        }
      }
    }
    acc <- 1
    rwz <- list(Q=Q, acc=acc, ika=ik)
    return(rwz)
  }

sanity.check.fevdplot <-
  function(fevddraws=fevddraws, type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw, table=table){
    #fevddraws
    Idim <- is.null(dim(fevddraws))
    
    if(Idim==TRUE){
      
      message(" ")
      stop(" fevddraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Idiml <- length(dim(fevddraws))
    
    if(Idiml!=3){
      
      message(" ")
      stop(" fevddraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Ina <- any(is.na(fevddraws)==TRUE)
    Inan <- any(is.nan(fevddraws)==TRUE)
    Inum <-  is.numeric(fevddraws)
    
    if(Ina==TRUE | Inan==TRUE){
      
      message(" ")
      stop(" fevddraws must not contain missing values.\n", call. = FALSE)
    }
    
    #type
    if(type!="median" & type!="mean"){
      
      message(" ")
      stop("FEVD type must be mean or median.\n", call. = FALSE)
    }
    #labels
    lbll <- length(labels)
    
    if(lbll < dim(fevddraws)[3]){
      
      message(" ")
      stop("Number of labels must be equal number of variables in the model.\n", call. = FALSE)
    }
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      bndsl <- length(bands)
      
      if(bndsl !=2){
        
        message(" ")
        stop("Error bands must contain only two values c(lower, upper) or 'NULL'.\n", call. = FALSE)
      }
      
      
      if(max(bands)>=1 | min(bands)<=0){
        
        message(" ")
        stop("Error bands must be between 0 and 1.\n", call. = FALSE)
      }
    }
    return()
  }

sanity.check.irfplot <-
  function(irfdraws=irfdraws,type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw){
    #irfdraws
    Idim <- is.null(dim(irfdraws))
    
    if(Idim==TRUE){
      
      message(" ")
      stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Idiml <- length(dim(irfdraws))
    
    if(Idiml!=3){
      
      message(" ")
      stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Ina <- any(is.na(irfdraws)==TRUE)
    Inan <- any(is.nan(irfdraws)==TRUE)
    Inum <-  is.numeric(irfdraws)
    
    if(Ina==TRUE | Inan==TRUE){
      
      message(" ")
      stop(" Irfdraws must not contain missing values.\n", call. = FALSE)
    }
    
    #type
    if(type!="median" & type!="mean"){
      
      message(" ")
      stop("IRF type must be mean or median.\n", call. = FALSE)
    }
    #labels
    lbll <- length(labels)
    
    if(lbll < dim(irfdraws)[3]){
      
      message(" ")
      stop("Number of labels must be equal number of variables in the model.\n", call. = FALSE)
    }
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      bndsl <- length(bands)
      
      if(bndsl !=2){
        
        message(" ")
        stop("Error bands must contain only two values c(lower, upper) or 'NULL'.\n", call. = FALSE)
      }
      
      
      if(max(bands)>=1 | min(bands)<=0){
        
        message(" ")
        stop("Error bands must be between 0 and 1.\n", call. = FALSE)
      }
    }
    return()
  }

sanity.check.reject <-
  function(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps){
    
    Yts <-is.ts(Y)
    Ydf <- is.data.frame(Y)
    
    if(Yts==FALSE & Ydf==FALSE){
      
      message(" ")
      stop(" Data has to be a data.frame() or ts() object.\n", call. = FALSE)
    }
    
    if(ncol(Y)<2){
      
      message(" ")
      stop(" Need more than 1 variable.\n", call. = FALSE)
    }
    
    Yna <- any(is.na(Y)==TRUE)
    Ynan <- any(is.nan(Y)==TRUE)
    Ynum <-  is.numeric(Y)
    
    if(Yna==TRUE | Ynan==TRUE){
      
      message(" ")
      stop(" Data must not contain missing values.\n", call. = FALSE)
    }
    
    if(Ynum!=TRUE){
      
      message(" ")
      stop(" Data must not contain strings.\n", call. = FALSE)
    }
    
    nlagsint <- nlags%%1==0
    nlagsobs <- nrow(Y)-nlags
    
    if(nlagsint==FALSE){
      
      message(" ")
      stop("Number of lags must be integer.\n", call. = FALSE)
    }
    
    if(nlagsobs<=0){
      
      message(" ")
      stop(" Number of lags cannot be larger than the number of observations.\n", call. = FALSE)
    }
    
    if(nlags<1){
      
      message(" ")
      stop(" Need at least 1 lag.\n")
    }
    if(nlags>nrow(Y)){
      
      message(" ")
      stop(" Number of lags have to be smaller than number of observations.\n", call. = FALSE)
    }
    
    drawsint <- draws%%1==0
    
    if(drawsint!=TRUE){
      
      message(" ")
      stop(" Number of draws must be integer.\n", call. = FALSE)
    }
    
    if(draws<=0){
      
      message(" ")
      stop(" Number of draws must geater than zero.\n", call. = FALSE)
    }
    
    stepsint <- steps%%1==0
    
    if(stepsint!=TRUE){
      
      message(" ")
      stop(" Number of steps must be integer.\n", call. = FALSE)
    }
    
    if(steps<=0){
      
      message(" ")
      stop(" Number of steps must geater than zero.\n", call. = FALSE)
    }
    
    
    subdrawsint <- subdraws%%1==0
    
    if(subdrawsint!=TRUE){
      
      message(" ")
      stop(" Number of subdraws must be integer.\n", call. = FALSE)
    }
    
    if(subdraws<=0){
      
      message(" ")
      stop(" Number of subdraws must geater than zero.\n", call. = FALSE)
    }
    
    nkeepint <- nkeep%%1==0
    
    if(nkeepint!=TRUE){
      
      message(" ")
      stop(" Number of nkeep must be integer.\n", call. = FALSE)
    }
    
    if(nkeep<=0){
      
      message(" ")
      stop(" Number of nkeep must geater than zero.\n", call. = FALSE)
    }
    
    
    KMINint <- KMIN%%1==0
    
    if(KMINint!=TRUE){
      
      message(" ")
      stop("KMIN must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMIN must be greater than zero and smaller than KMAX.\n", call. = FALSE)
    }
    
    
    KMAXint <- KMAX%%1==0
    
    if(KMAXint!=TRUE){
      
      message(" ")
      stop("KMAX must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMAX must be greater than zero and greater than KMIN.\n", call. = FALSE)
    }
    
    
    cnsl <- length(constrained)
    
    if(cnsl<=0){
      
      message(" ")
      stop("Number of constraints must at least 1.\n", call. = FALSE)
    }
    
    if(cnsl>ncol(Y)){
      
      message(" ")
      stop("Number of constraints cannot be greater than number of variables.\n", call. = FALSE)
    }
    
    if(max(abs(constrained))>ncol(Y) | min(abs(constrained))==0){
      
      message(" ")
      stop("Constraints must be between 1 and the number of variables.\n", call. = FALSE)
    }
    
    cnscns <- all(constrained%%1==0)
    
    if(cnscns==FALSE){
      
      message(" ")
      stop("All constraints must be integers.\n", call. = FALSE)
    }
    
    cnsdup <- anyDuplicated(abs(constrained))
    if(cnsdup>0){
      
      message(" ")
      stop("Cannot provide multiple constraints for the same variable.\n", call. = FALSE)
    }
    
    return()
    
  }

sanity.check.reject <-
  function(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps){
    
    Yts <-is.ts(Y)
    Ydf <- is.data.frame(Y)
    
    if(Yts==FALSE & Ydf==FALSE){
      
      message(" ")
      stop(" Data has to be a data.frame() or ts() object.\n", call. = FALSE)
    }
    
    if(ncol(Y)<2){
      
      message(" ")
      stop(" Need more than 1 variable.\n", call. = FALSE)
    }
    
    Yna <- any(is.na(Y)==TRUE)
    Ynan <- any(is.nan(Y)==TRUE)
    Ynum <-  is.numeric(Y)
    
    if(Yna==TRUE | Ynan==TRUE){
      
      message(" ")
      stop(" Data must not contain missing values.\n", call. = FALSE)
    }
    
    if(Ynum!=TRUE){
      
      message(" ")
      stop(" Data must not contain strings.\n", call. = FALSE)
    }
    
    nlagsint <- nlags%%1==0
    nlagsobs <- nrow(Y)-nlags
    
    if(nlagsint==FALSE){
      
      message(" ")
      stop("Number of lags must be integer.\n", call. = FALSE)
    }
    
    if(nlagsobs<=0){
      
      message(" ")
      stop(" Number of lags cannot be larger than the number of observations.\n", call. = FALSE)
    }
    
    if(nlags<1){
      
      message(" ")
      stop(" Need at least 1 lag.\n")
    }
    if(nlags>nrow(Y)){
      
      message(" ")
      stop(" Number of lags have to be smaller than number of observations.\n", call. = FALSE)
    }
    
    drawsint <- draws%%1==0
    
    if(drawsint!=TRUE){
      
      message(" ")
      stop(" Number of draws must be integer.\n", call. = FALSE)
    }
    
    if(draws<=0){
      
      message(" ")
      stop(" Number of draws must geater than zero.\n", call. = FALSE)
    }
    
    stepsint <- steps%%1==0
    
    if(stepsint!=TRUE){
      
      message(" ")
      stop(" Number of steps must be integer.\n", call. = FALSE)
    }
    
    if(steps<=0){
      
      message(" ")
      stop(" Number of steps must geater than zero.\n", call. = FALSE)
    }
    
    
    subdrawsint <- subdraws%%1==0
    
    if(subdrawsint!=TRUE){
      
      message(" ")
      stop(" Number of subdraws must be integer.\n", call. = FALSE)
    }
    
    if(subdraws<=0){
      
      message(" ")
      stop(" Number of subdraws must geater than zero.\n", call. = FALSE)
    }
    
    nkeepint <- nkeep%%1==0
    
    if(nkeepint!=TRUE){
      
      message(" ")
      stop(" Number of nkeep must be integer.\n", call. = FALSE)
    }
    
    if(nkeep<=0){
      
      message(" ")
      stop(" Number of nkeep must geater than zero.\n", call. = FALSE)
    }
    
    
    KMINint <- KMIN%%1==0
    
    if(KMINint!=TRUE){
      
      message(" ")
      stop("KMIN must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMIN must be greater than zero and smaller than KMAX.\n", call. = FALSE)
    }
    
    
    KMAXint <- KMAX%%1==0
    
    if(KMAXint!=TRUE){
      
      message(" ")
      stop("KMAX must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMAX must be greater than zero and greater than KMIN.\n", call. = FALSE)
    }
    
    
    cnsl <- length(constrained)
    
    if(cnsl<=0){
      
      message(" ")
      stop("Number of constraints must at least 1.\n", call. = FALSE)
    }
    
    if(cnsl>ncol(Y)){
      
      message(" ")
      stop("Number of constraints cannot be greater than number of variables.\n", call. = FALSE)
    }
    
    if(max(abs(constrained))>ncol(Y) | min(abs(constrained))==0){
      
      message(" ")
      stop("Constraints must be between 1 and the number of variables.\n", call. = FALSE)
    }
    
    cnscns <- all(constrained%%1==0)
    
    if(cnscns==FALSE){
      
      message(" ")
      stop("All constraints must be integers.\n", call. = FALSE)
    }
    
    cnsdup <- anyDuplicated(abs(constrained))
    if(cnsdup>0){
      
      message(" ")
      stop("Cannot provide multiple constraints for the same variable.\n", call. = FALSE)
    }
    
    return()
    
  }

sanity.check.target <-
  function(Y=Y, nlags=nlags, irfdraws=irfdraws, constant=constant, type=type, labels=labels, target= target, save=save, legend=legend, bands=bands, grid=grid, bw=bw, maxit=maxit){
    
    
    Yts <-is.ts(Y)
    Ydf <- is.data.frame(Y)
    
    if(Yts==FALSE & Ydf==FALSE){
      
      message(" ")
      stop(" Data has to be a data.frame() or ts() object.\n", call. = FALSE)
    }
    
    if(ncol(Y)<2){
      
      message(" ")
      stop(" Need more than 1 variable.\n", call. = FALSE)
    }
    
    Yna <- any(is.na(Y)==TRUE)
    Ynan <- any(is.nan(Y)==TRUE)
    Ynum <-  is.numeric(Y)
    
    if(Yna==TRUE | Ynan==TRUE){
      
      message(" ")
      stop(" Data must not contain missing values.\n", call. = FALSE)
    }
    
    if(Ynum!=TRUE){
      
      message(" ")
      stop(" Data must not contain strings.\n", call. = FALSE)
    }
    
    nlagsint <- nlags%%1==0
    nlagsobs <- nrow(Y)-nlags
    
    if(nlagsint==FALSE){
      
      message(" ")
      stop("Number of lags must be integer.\n", call. = FALSE)
    }
    
    if(nlagsobs<=0){
      
      message(" ")
      stop(" Number of lags cannot be larger than the number of observations.\n", call. = FALSE)
    }
    
    if(nlags<1){
      
      message(" ")
      stop(" Need at least 1 lag.\n")
    }
    if(nlags>nrow(Y)){
      
      message(" ")
      stop(" Number of lags have to be smaller than number of observations.\n", call. = FALSE)
    }
    
    
    Idim <- is.null(dim(irfdraws))
    
    if(Idim==TRUE){
      
      message(" ")
      stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Idiml <- length(dim(irfdraws))
    
    if(Idiml!=3){
      
      message(" ")
      stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }
    
    Ina <- any(is.na(irfdraws)==TRUE)
    Inan <- any(is.nan(irfdraws)==TRUE)
    Inum <-  is.numeric(irfdraws)
    
    if(Ina==TRUE | Inan==TRUE){
      
      message(" ")
      stop(" Irfdraws must not contain missing values.\n", call. = FALSE)
    }
    
    #type
    if(type!="median" & type!="mean"){
      
      message(" ")
      stop("IRF type must be mean or median.\n", call. = FALSE)
    }
    #labels
    lbll <- length(labels)
    
    if(lbll < dim(irfdraws)[3]){
      
      message(" ")
      stop("Number of labels must be equal number of variables in the model.\n", call. = FALSE)
    }
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
      bndsl <- length(bands)
      
      if(bndsl !=2){
        
        message(" ")
        stop("Error bands must contain only two values c(lower, upper) or 'NULL'.\n", call. = FALSE)
      }
      
      
      if(max(bands)>=1 | min(bands)<=0){
        
        message(" ")
        stop("Error bands must be between 0 and 1.\n", call. = FALSE)
      }
    }
    
    
    
    
    mxint <- maxit%%1==0
    
    if(mxint!=TRUE){
      
      message(" ")
      stop(" Number of maxit must be integer.\n", call. = FALSE)
    }
    
    if(maxit<=0){
      
      message(" ")
      stop(" Number of maxit must geater than zero.\n", call. = FALSE)
    }
    return()
  }

sanity.check.uhlig.penalty <-
  function(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps, penalty=penalty, crit=crit){
    
    Yts <-is.ts(Y)
    Ydf <- is.data.frame(Y)
    
    if(Yts==FALSE & Ydf==FALSE){
      
      message(" ")
      stop(" Data has to be a data.frame() or ts() object.\n", call. = FALSE)
    }
    
    if(ncol(Y)<2){
      
      message(" ")
      stop(" Need more than 1 variable.\n", call. = FALSE)
    }
    
    Yna <- any(is.na(Y)==TRUE)
    Ynan <- any(is.nan(Y)==TRUE)
    Ynum <-  is.numeric(Y)
    
    if(Yna==TRUE | Ynan==TRUE){
      
      message(" ")
      stop(" Data must not contain missing values.\n", call. = FALSE)
    }
    
    if(Ynum!=TRUE){
      
      message(" ")
      stop(" Data must not contain strings.\n", call. = FALSE)
    }
    
    nlagsint <- nlags%%1==0
    nlagsobs <- nrow(Y)-nlags
    
    if(nlagsint==FALSE){
      
      message(" ")
      stop("Number of lags must be integer.\n", call. = FALSE)
    }
    
    if(nlagsobs<=0){
      
      message(" ")
      stop(" Number of lags cannot be larger than the number of observations.\n", call. = FALSE)
    }
    
    if(nlags<1){
      
      message(" ")
      stop(" Need at least 1 lag.\n")
    }
    if(nlags>nrow(Y)){
      
      message(" ")
      stop(" Number of lags have to be smaller than number of observations.\n", call. = FALSE)
    }
    
    drawsint <- draws%%1==0
    
    if(drawsint!=TRUE){
      
      message(" ")
      stop(" Number of draws must be integer.\n", call. = FALSE)
    }
    
    if(draws<=0){
      
      message(" ")
      stop(" Number of draws must geater than zero.\n", call. = FALSE)
    }
    
    stepsint <- steps%%1==0
    
    if(stepsint!=TRUE){
      
      message(" ")
      stop(" Number of steps must be integer.\n", call. = FALSE)
    }
    
    if(steps<=0){
      
      message(" ")
      stop(" Number of steps must geater than zero.\n", call. = FALSE)
    }
    
    
    subdrawsint <- subdraws%%1==0
    
    if(subdrawsint!=TRUE){
      
      message(" ")
      stop(" Number of subdraws must be integer.\n", call. = FALSE)
    }
    
    if(subdraws<=0){
      
      message(" ")
      stop(" Number of subdraws must geater than zero.\n", call. = FALSE)
    }
    
    nkeepint <- nkeep%%1==0
    
    if(nkeepint!=TRUE){
      
      message(" ")
      stop(" Number of nkeep must be integer.\n", call. = FALSE)
    }
    
    if(nkeep<=0){
      
      message(" ")
      stop(" Number of nkeep must geater than zero.\n", call. = FALSE)
    }
    
    
    KMINint <- KMIN%%1==0
    
    if(KMINint!=TRUE){
      
      message(" ")
      stop("KMIN must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMIN must be greater than zero and smaller than KMAX.\n", call. = FALSE)
    }
    
    
    KMAXint <- KMAX%%1==0
    
    if(KMAXint!=TRUE){
      
      message(" ")
      stop("KMAX must be integer.\n", call. = FALSE)
    }
    
    if(KMIN<=0 | KMIN>KMAX){
      
      message(" ")
      stop("KMAX must be greater than zero and greater than KMIN.\n", call. = FALSE)
    }
    
    
    cnsl <- length(constrained)
    
    if(cnsl<=0){
      
      message(" ")
      stop("Number of constraints must at least 1.\n", call. = FALSE)
    }
    
    if(cnsl>ncol(Y)){
      
      message(" ")
      stop("Number of constraints cannot be greater than number of variables.\n", call. = FALSE)
    }
    
    if(max(abs(constrained))>ncol(Y) | min(abs(constrained))==0){
      
      message(" ")
      stop("Constraints must be between 1 and the number of variables.\n", call. = FALSE)
    }
    
    cnscns <- all(constrained%%1==0)
    
    if(cnscns==FALSE){
      
      message(" ")
      stop("All constraints must be integers.\n", call. = FALSE)
    }
    
    cnsdup <- anyDuplicated(abs(constrained))
    if(cnsdup>0){
      
      message(" ")
      stop("Cannot provide multiple constraints for the same variable.\n", call. = FALSE)
    }
    
    if(penalty<=0){
      
      message(" ")
      stop("Penalty must be greater than 0.\n", call. = FALSE)
    }
    
    if(crit<0){
      
      message(" ")
      stop("Crit must be greater than 0.\n", call. = FALSE)
    }
    
    return()
  }

stereo <-
  function(v){
    llength <- length(v)
    ll <- llength + 1
    s <- matrix(NA, nrow=ll, ncol=1)
    for(i in 1:llength){
      s[i] <- 2*v[i] / (sqrt(sum(v^2))^2 +1)
    }
    s[ll] <- (sqrt(sum(v^2))^2 -1) / (sqrt(sum(v^2))^2 +1)
    return(s)
  }

uhlig.penalty <-
  function(Y=NULL, nlags=4, draws=1000, subdraws=1000, nkeep=1000, KMIN=1, KMAX=4, constrained=NULL, constant=TRUE, steps=24, penalty=100, crit=0.001){
    #
    #--- SANITY CHECK ---#
    sanity.check.uhlig.penalty(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps, penalty=penalty, crit=crit)
    #
    #--- SET UP PARAS ---#
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
    pen <- penalty
    ntot <- n1
    #
    if(constant == FALSE){
      CONS <- "F"
      ncoef <- nvar * nlags
      nncoef <- nvar * nlags
      nnvar1 <- nvar * (nlags + 1)
    }else{
      CONS <- "T"
      ncoef <- nvar * (nlags+1)
      nncoef <- nvar * nlags + 1
      nnvar1 <- nvar * (nlags + 1) + 1
    }
    #
    #---REDUCED FORM VAR MODEL ---#
    model <- rfvar(Y,lags=nlags, const=CONS, breaks=NULL)
    bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
    resid <- model$u # same as above
    data <- model$X
    xx <- model$xx
    #
    #--- SIGMA and SXX ---#
    uu <- crossprod(resid)
    # sigma <- (1/(nnobs-nncoef))*uu
    sigma <- (1/nnobs)*uu
    #
    #--- SET UP MCMC OF VAR ---#
    sxx <-  chol(xx)
    sv <- solve(uu)
    svt <-  chol(sv)
    betaols <- t(bcoef)
    best <- betaols
    wishdof <- nnobs-nncoef
    #
    #--- MATRICES FOR DRAWS ---#
    goodresp <- array(NA, c(nkeep, nstep, nvar))
    BDraws <- array(NA, c(n1, nncoef, nvar))
    SDraws <- array(NA, c(n1, nvar, nvar))
    imp <- matrix(NA, nrow=nstep, ncol=nvar)
    fevd <- matrix(NA, nrow=nstep, ncol=nvar)
    goodfevd <- array(NA, c(nkeep, nstep, nvar))
    goodshock <- array(NA, c(nkeep, nnobs))
    uhatt <- matrix(NA, nnobs, 1)
    #
    #--- WEIGHTING MATRIX ---#
    YD <- apply(Y, 2, diff)
    scales <-  matrix(apply(matrix(apply(YD, 2, var)),1,sqrt))
    #
    # PROJECTION to the unit sphere in R^n.
    g <- matrix(1, nrow=nvar-1, ncol=1)
    #
    #--- Monte CARLO INTEGRATION ---#
    accept <- 0
    convcnt <- 0
    # expt <- round(0.003858333 * n1)
    message('Starting MCMC, ', date(),'.', sep="")
    pb0 <- txtProgressBar(min = 0, max = n1, style = 3)
    for(draws in 1:n1){
      setTxtProgressBar(pb0, draws)
      #
      #--- sigma draws ---#
      sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
      swish   <- chol(sigmad)
      #
      #--- beta draws ---#
      swsxx <-   sigmad  %x% xx
      bd <- rep(0, nrow(swsxx))
      #betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
      betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
      betadraw <- betaols + betau
      bhat <- betadraw
      #
      #--- irfs ---#
      imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
      impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
      imp2 <- impulses^2
      imp2sum <- apply(imp2, c(2,3), cumsum)
      mse <-  apply(imp2sum, c(1,2), sum)
      fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))
      #
      #--- PENALTY FUNCTION EVAL ---#
      penaltyfunc <- minqa::uobyqa(g, fn=UhligPenalty,   control = list(maxfun=n2),  first=KMIN, last=KMAX, constrained=constrained, impulses=impulses, scales=scales, pen=pen)
      #
      betat <- penaltyfunc$par
      UAT <- penaltyfunc$fval
      convt <- penaltyfunc$ierr
      #
      penaltyfunc <- minqa::uobyqa(g, fn=UhligPenalty,   control = list(maxfun=n2),  first=KMIN, last=KMAX, constrained=constrained, impulses=impulses, scales=scales, pen=pen)
      #
      beta <- penaltyfunc$par
      UA <- penaltyfunc$fval
      conv <- penaltyfunc$ierr
      
      #--- convergence check ---#
      convcnt <- convcnt + max(c(convt,conv)) # that is not correct max can be greater than "1"
      if(convt>1 | conv>1){
        next
      }
      #
      if(abs(UAT-UA)<=crit){
        a <- stereo(beta)
        for(j in 1:nstep){
          imp[j,] <- t(impulses[j,,]%*%a)
          fevd[j,] <- t(fevd0[j,,]%*%(a^2))
        }
        accept <- accept+1
        goodresp[accept, ,] <-  imp
        goodfevd[accept, ,] <- fevd * 100
        BDraws[draws, , ] <- betadraw
        SDraws[draws, , ] <- sigmad
        uhat <-   Y[nnobs0:nobs ,] - data %*% bhat
        for(i in 1:nnobs){
          uhatt[i,] <-   uhat[i, ] %*%  (  solve(swish) %*% a)
        }
        goodshock[accept, ] <-  t(uhatt)
      }
      #
      if(accept>=nkeep){
        break
      }
      ldraw <- draws
    }# END DRAWS
    close(pb0)
    #
    #--- FIX PARA MATRICES ---#
    if(ldraw<n1){
      BDraws <- BDraws[1:ldraw, , ]
      SDraws <- SDraws[1:ldraw, , ]
      dimnames(SDraws) <- list(1:ldraw, varnames, varnames)
    }
    #
    #--- WARNING MESSAGE IN CASE OF TOO FEW DRAWS ---#
    if(accept<nkeep){
      if(accept==0){
        stop("\n Not enough accepted draws to proceed!")
      }else{
        goodresp <- goodresp[1:accept, , ]
        goodfevd <- goodfevd[1:accept, , ]
        goodshock <- goodshock[1:accept, ]
        message('\n Warning! Had only ', accept,' accepted draw(s) out of ',ntot,'. ', convcnt, ' draws did not converge.', sep="")
      }
    }
    #
    if(accept>=nkeep & convcnt>0){
      goodresp <- goodresp[1:accept, , ]
      goodfevd <- goodfevd[1:accept, , ]
      goodshock <- goodshock[1:accept, ]
      message('\n Warning! ', convcnt, ' draw(s) did not converge.', sep="")
    }
    nn1 <- accept
    dimnames(goodresp) <- list(1:nn1, 1:nstep, varnames)
    dimnames(goodfevd) <- list(1:nn1, 1:nstep, varnames)
    #
    if(constant == FALSE){
      dimnames(BDraws) <-  list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
        dimnames(BDraws) <- list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
      }
    #
    message('\n MCMC finished, ', date(),'.', sep="")
    return(list(IRFS=goodresp, FEVDS = goodfevd,  SHOCKS = goodshock, BDraws=BDraws, SDraws=SDraws))
  }

uhlig.reject <-
  function(Y=NULL,  nlags=4, draws=200, subdraws=200, nkeep=1000, KMIN=1, KMAX=4, constrained=NULL, constant=TRUE, steps=24){
    #
    #---SANITY CHECK ---#
    sanity.check.reject(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps)
    #
    #--- SET UP PARAS ---#
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
    ntot <- n1*n2
    #
    if(constant == FALSE){
      CONS <- "F"
      ncoef <- nvar * nlags
      nncoef <- nvar * nlags
      nnvar1 <- nvar * (nlags + 1)
    }else{
      CONS <- "T"
      ncoef <- nvar * (nlags+1)
      nncoef <- nvar * nlags + 1
      nnvar1 <- nvar * (nlags + 1) + 1
    }
    #
    #---REDUCED FORM VAR MODEL ---#
    model <- rfvar(ydata=Y, lags=nlags, const=CONS)
    bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
    resid <- model$u # same as above
    data <- model$X
    xx <- model$xx
    #
    #--- SIGMA and SXX ---#
    uu <- crossprod(resid)
    # sigma <- (1/(nnobs-nncoef))*uu
    sigma <- (1/nnobs)*uu
    #
    #--- SET UP MCMC OF VAR ---#
    sxx <-  chol(xx)
    sv <- solve(uu)
    svt <-  chol(sv)
    betaols <- t(bcoef)
    best <- betaols
    wishdof <- nnobs-nncoef
    #
    #--- MATRICES FOR DRAWS ---#
    goodresp <- array(NA, c(nkeep, nstep, nvar))
    BDraws <- array(NA, c(n1, nncoef, nvar))
    SDraws <- array(NA, c(n1, nvar, nvar))
    imp <- matrix(NA, nrow=nstep, ncol=nvar)
    fevd <- matrix(NA, nrow=nstep, ncol=nvar)
    goodfevd <- array(NA, c(nkeep, nstep, nvar))
    goodshock <- array(NA, c(nkeep, nnobs))
    uhatt <- matrix(NA, nnobs, 1)
    #
    #--- Monte CARLO INTEGRATION ---#
    accept <- 0
    message('Starting MCMC, ', date(),'.', sep="")
    pb0 <- txtProgressBar(min = 0, max = n1, style = 3)
    for(draws in 1:n1){
      setTxtProgressBar(pb0, draws)
      #
      #--- sigma draws ---#
      sigmad  <- solve(matrix(rWishart(1, wishdof, sv), nrow=nvar, ncol=nvar))
      swish   <- chol(sigmad)
      #
      #--- beta draws ---#
      swsxx <-   sigmad  %x% xx
      bd <- rep(0, nrow(swsxx))
      #betau <- matrix(mvrnormR(1,0,swsxx), nrow=nncoef, ncol=nvar)
      betau <- matrix(mvnfast::rmvn(1, bd, swsxx), nrow=nncoef, ncol=nvar)
      betadraw <- betaols + betau
      bhat <- betadraw
      #
      #--- irfs ---#
      imfhat <- fn.impulse(bhat, swish, c(nvar, nlags, nstep))
      impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
      imp2 <- impulses^2
      imp2sum <- apply(imp2, c(2,3), cumsum)
      mse <-  apply(imp2sum, c(1,2), sum)
      fevd0 <- array(apply(imp2sum, 3, "/",  mse), dim=c(nstep, nvar, nvar))
      #
      for(subdraws in 1:n2){
        a <- matrix(HI::rballunif(nvar,1), nvar, 1)
        UAR <- UhligAccept(a,KMIN,KMAX,constrained, impulses)
        UA <- UAR$acc
        q <- UAR$Q
        if(UA==1){
          for(j in 1:nstep){ # this can be done via apply
            imp[j,] <- t(impulses[j,,]%*%q)
            fevd[j,] <- t(fevd0[j,,]%*%(q^2))
          }
          accept <- accept+1
          goodresp[accept, ,] <-  imp
          goodfevd[accept, ,] <- fevd * 100
          BDraws[draws, , ] <- betadraw
          SDraws[draws, , ] <- sigmad
          uhat <-   Y[nnobs0:nobs ,] - data %*% bhat
          for(i in 1:nnobs){
            uhatt[i,] <-   uhat[i, ] %*%  (  solve(swish) %*% q)
          }
          goodshock[accept, ] <-  t(uhatt)
        }else{
          next
        }
        #
        if(accept>=nkeep){
          break
        }
        #
      } # end subdraws
      if(accept>=nkeep){
        break
      }
      ldraw <- draws
    }#end draws
    close(pb0)
    #
    #--- FIX PARA MATRICES ---#
    if(ldraw<n1){
      BDraws <- BDraws[1:ldraw, , ]
      SDraws <- SDraws[1:ldraw, , ]
      dimnames(SDraws) <- list(1:ldraw, varnames, varnames)
    }
    #
    #--- WARNING MESSAGE IN CASE OF TOO FEW DRAWS ---#
    if(accept<nkeep){
      if(accept==0){
        stop("\n Not enough accepted draws to proceed!")
      }else{
        goodresp <- goodresp[1:accept, , ]
        goodfevd <- goodfevd[1:accept, , ]
        goodshock <- goodshock[1:accept, ]
        message('\n Warning! Had only ', accept,' accepted draw(s) out of ',ntot,'.', sep="")
      }
    }
    nn1 <- accept
    dimnames(goodresp) <- list(1:nn1, 1:nstep, varnames)
    dimnames(goodfevd) <- list(1:nn1, 1:nstep, varnames)
    #
    if(constant == FALSE){
      dimnames(BDraws) <-  list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep="")) , varnames)}else{
        dimnames(BDraws) <- list(1:ldraw, c(paste(varnames,rep(1:nlags, each=length(varnames)), sep=""),"const"), varnames)
      }
    #
    message('\n MCMC finished, ', date(),'.', sep="")
    return(list(IRFS=goodresp, FEVDS = goodfevd,  SHOCKS = goodshock, BDraws=BDraws, SDraws=SDraws))
  }

UhligAccept <-
  function(Q, first, last, constrained, impulses){#ok
    for(k in first:last){#ok
      ik <- impulses[k, , ]%*%Q#ok
      for(i in 1:length(constrained)){#ok
        if(constrained[i]<0){#ok
          value <- ik[-1.0 * constrained[i]]#ok
        }else{#ok
          value <- -1.0 * ik[constrained[i]]#ok
        }#ok
        if(value>0.0){#ok
          if(k==first & i==1){#ok
            Q <- -1.0 * Q#ok
            ik <- -1.0 * ik#ok
          }else{#ok
            acc <- 0
            uar <- list(Q=Q, acc=acc, ika=ik)
            return(uar)
          }#ok
        }#ok
      }#end i #ok
    }#end k #ok
    acc <- 1
    uar <- list(Q=Q, acc=acc, ika=ik)
    return(uar)
  }

UhligPenalty <-
  function(g, first, last, constrained, impulses, scales, pen){
    func <- 0.0
    q <- matrix(stereo(v=g))
    for(k in first:last){
      ik <- (impulses[k, , ]%*%q) / scales
      for(i in 1:length(constrained)){
        if(constrained[i]<0){
          value <- ik[-1.0*constrained[i]]
        }else{
          value <- -1.0 * ik[constrained[i]]
        }
        if(value<0){
          func <- func +  value
        }else{
          func <- func + pen * value
        }
      }
    }
    acc <- func
    return(acc)
  }

.Random.seed <-
  c(403L, 624L, 1638542565L, 108172386L, -1884566405L, -1838154368L, 
    -250773631L, 919185230L, -1001918601L, -1002779316L, -321961507L, 
    1781331706L, 1166440499L, -117712936L, 58314745L, -938201242L, 
    1844978351L, -869947100L, 216221141L, 576699538L, 533687019L, 
    1819381040L, 675905393L, 2110048894L, 2136771815L, 1026265084L, 
    -1037166387L, -285653718L, 1263109283L, -1975616120L, 1576168937L, 
    -123742058L, 244715039L, 1526718932L, -1087137083L, 1442050242L, 
    676570459L, 800852192L, -823755935L, -566904402L, 1696695895L, 
    906100396L, 1479781309L, -321511590L, -1503089389L, 1068470072L, 
    2031323097L, 2053295558L, -749218417L, -2000861564L, 1655319477L, 
    -768462606L, 262110155L, 399143056L, -983338159L, -1865452322L, 
    -202515513L, 1164515676L, 80673453L, 1462142346L, 1143665027L, 
    -1238758168L, 255598025L, 1584402166L, 1601468671L, -748103884L, 
    1964374181L, -656573150L, 1708781115L, -1898462144L, 469722945L, 
    -888865778L, -807891657L, -53747700L, -1448147555L, -905087046L, 
    -208186893L, 289994392L, -2104807495L, -895837146L, -1276970897L, 
    -1949461532L, -133824107L, -327625390L, 1424620715L, -572587024L, 
    45700913L, -284602562L, 850959015L, -1739238724L, -1539126131L, 
    -1167198742L, -713765277L, -1419293624L, -857752151L, 695563606L, 
    -1621470241L, -2055832428L, 1623803525L, 4667778L, 276211483L, 
    -593809504L, -1185922271L, -1144033682L, 1545929751L, -1359974036L, 
    -1111928963L, -1511325670L, -867539245L, -979366408L, 1876477849L, 
    1515428486L, 777096015L, -961638076L, -2006005899L, -1471436366L, 
    1172761995L, -1584969904L, 2135108369L, 2097830302L, 24430983L, 
    -498582500L, 453086829L, 1122473546L, -329311421L, 936262568L, 
    1691700617L, -515372106L, 452959935L, 947966452L, -1681554331L, 
    1229530594L, -1939746821L, 836651776L, -2083451135L, 992809166L, 
    -1111561481L, -1899515188L, 751470941L, -1333348230L, -340137037L, 
    546100568L, 217337721L, 388352230L, 1029410351L, 1454915236L, 
    190610773L, 1220718098L, -724677013L, 832256688L, -705106191L, 
    -345336834L, -2116392857L, 2073679228L, -1618788275L, -1498716502L, 
    -1843275741L, -1791567096L, 555011433L, 1501549078L, -82028129L, 
    -538978476L, 2072162885L, 1023099458L, -640456485L, -1820780960L, 
    1417267937L, -1615469778L, 73286103L, -1968593876L, 1264235325L, 
    -1510432550L, 689823891L, 1374112952L, -1579825319L, 984164166L, 
    -1212612337L, -2059232252L, -1270406347L, 485876338L, -1881661621L, 
    1223876112L, -1647140143L, -1228800418L, 832665415L, 1755455196L, 
    608166445L, 668034826L, -336263933L, 1769548392L, -946455223L, 
    -1313552266L, 1147700351L, -1695839052L, -1854353371L, 1886752418L, 
    -1794936389L, -230452800L, 4355777L, 201450894L, -1682241353L, 
    1022248332L, 702663965L, -833046214L, 2007909747L, -175672296L, 
    -227201223L, 1249228198L, 1344398319L, -945444508L, -223954667L, 
    2147305170L, -1589878741L, -1905905296L, 1774732977L, 665360574L, 
    -360581593L, 1505302588L, 1471116301L, -1804494998L, 990945763L, 
    -965924408L, -1705927383L, 1434382038L, -677634209L, -1258556908L, 
    -1423974907L, -2066740478L, -165025125L, 722844960L, 1478693537L, 
    1856575470L, 1267548055L, -336750868L, -1797794051L, -37813862L, 
    -425518509L, 323308408L, 1053460249L, 404976646L, -1790036273L, 
    -1086757180L, 1911766773L, -655303886L, -818736885L, -1798490928L, 
    -925771119L, 1387684638L, -499915513L, -1422474852L, -1538656787L, 
    1185150922L, -592662845L, 703256872L, 1463741705L, -59357898L, 
    1903110719L, -1519843468L, -772809755L, 656586594L, -780217477L, 
    113744000L, 689151617L, -2109506994L, 931980919L, -1979738036L, 
    347394269L, -1807517190L, -1590403277L, 519623384L, 1120784121L, 
    -1052089754L, -235538001L, 961926180L, 430224597L, -1696030830L, 
    2079601131L, -620763088L, 1172790897L, 411262334L, -1435548697L, 
    1702042364L, 814180301L, 612403242L, 1241590691L, 2110404744L, 
    945171689L, -1439511658L, -1432771297L, 260754644L, 1264634309L, 
    277189570L, -1757795237L, 976298976L, 1007426145L, -848753490L, 
    -546177705L, -1250145876L, -302990659L, 2113806938L, -291892205L, 
    -126219720L, 913770201L, -1350401850L, -1395577713L, 793964932L, 
    221452981L, 1157403634L, -1514683701L, -875148400L, 1744884305L, 
    631736286L, 820777671L, 1019618396L, -555759187L, -1608562550L, 
    407246979L, 400771048L, -208708407L, -1370717706L, -137129985L, 
    -1028046284L, -1729452123L, 101752866L, 1402205499L, 1714052928L, 
    1443137089L, -1565405426L, 519341111L, -1195660532L, 553882781L, 
    856095418L, 839661811L, -341773928L, -843174215L, -1738289370L, 
    -392704145L, -964116764L, -1407815531L, 1648670802L, -424295509L, 
    -1104650512L, -1307167183L, -222633410L, -1084075609L, -1853366852L, 
    1305156493L, -1214760726L, -122456733L, -1173485752L, -1159562071L, 
    -1537513386L, -1445664033L, -1169397868L, 2013623685L, -656679806L, 
    -1362874853L, 395004576L, 938795553L, 548782446L, 768375575L, 
    -1978286996L, 1785028221L, -1217002726L, -356332077L, -1337620232L, 
    949700249L, 2105953670L, -1543380401L, 1347370052L, -1849248139L, 
    -1982264142L, 1965438091L, -187817392L, -1558214127L, -1001034594L, 
    -74841977L, 1880114972L, -675193491L, -184329910L, -1199488445L, 
    -1843235160L, 1011321993L, 1945599670L, -313150017L, 496778484L, 
    -500616347L, 1711229154L, -318580997L, -953424384L, -1630196223L, 
    839705550L, -1615732233L, -874349108L, 1116608605L, -1687995526L, 
    -1275735373L, 1782567000L, 587615865L, -1400767514L, -1178114769L, 
    1241409956L, -1782846379L, 1590792466L, 591467883L, -1633708624L, 
    -1240150543L, -1425070338L, -417652889L, -1867030404L, -1824878771L, 
    1853401514L, 1088913187L, 972592648L, -1581872023L, 1354286358L, 
    -888278881L, 1073791572L, 214819141L, -1768757950L, -293081125L, 
    -639356576L, 1134389729L, -1924188626L, 1683798231L, -864424148L, 
    -571096515L, -85548070L, 1155352467L, -1452816456L, -1258863015L, 
    -1091642810L, -526362609L, 1559119620L, -881978827L, -1874442894L, 
    1597924939L, -696993520L, 1642987985L, -2088758946L, -540168633L, 
    1473587676L, 1387180333L, -972019190L, -1759630333L, -1217895064L, 
    -1859683255L, -1170785418L, 666213247L, -1596852300L, 1768652581L, 
    1605284258L, 797669563L, -1601426240L, -616196671L, -1256933234L, 
    -947585097L, -2043408244L, 756309021L, 2115517498L, 1890659443L, 
    1771400984L, -1838797255L, -1704662874L, -1421559057L, 1584845924L, 
    -2108348395L, -649123374L, 861284139L, -1557820304L, 430122417L, 
    -163566658L, -1621512921L, -708730052L, -1533689075L, 702667370L, 
    -597866269L, 2085217480L, 747789353L, 2081087958L, -1006324129L, 
    -345514732L, -1518727931L, -1033196030L, -934973029L, 1476122656L, 
    382054817L, -134911250L, 1893906071L, -1620516372L, -597563907L, 
    1494221978L, 648642899L, 312526456L, -615840231L, 1887184646L, 
    -2006151729L, 1341133252L, 1172909557L, 17055282L, 1170233355L, 
    -141946928L, 1277966737L, 2011657758L, 997661703L, -959132516L, 
    -748174099L, 1409661642L, 1306318275L, 1718948904L, 300887049L, 
    -1379157962L, 1018380607L, -149261708L, -1435399451L, -874587550L, 
    1827527291L, 688599936L, -1558856319L, 1888046414L, 1680727415L, 
    1785750348L, 1414946781L, 1273363706L, 1912489523L, 2019675608L, 
    731762169L, -1093888666L, -1026568017L, 1588723492L, -576576555L, 
    -629308782L, -659228437L, -1300610256L, 1686191473L, 1225650302L, 
    635304679L, -1821989380L, -442714419L, -1970025686L, 1154798243L, 
    -1077808248L, 1530260457L, -1290682730L, 175717407L, -951994412L, 
    -1647707963L, -1792854334L, 1841083227L, 680672992L, 690862433L, 
    90726318L, 18773079L, -442329940L, -1184249411L, -1565383334L, 
    2045213459L, -625965768L, -1688828455L, 1324233670L, -2128182385L, 
    -468411260L, 1291323829L, 1154676466L, -898889269L, -1630656880L, 
    -912641711L, 1984666334L, 742803911L, 1338978140L, -1649633107L, 
    -1616639094L, 914177923L, 1040745192L, -1637966903L, 997520630L, 
    -2112968961L, -1959416524L, -520399195L, 1129300770L, -1226179525L, 
    1466497600L, 1308992833L, 1764401678L, 57441079L, -1151896052L, 
    -334224483L, 898399674L, -2120409101L, -446370664L, -1100141127L, 
    913900070L, -970414481L, 1701833188L, -748494955L, 723295058L, 
    -1793226069L, 1935522288L, -563145423L, -667388610L, 2020083879L, 
    -1134139204L, 2047830669L, -246514710L, -1274143645L, 52478536L, 
    -312394839L, 1116560214L, -735346209L, -1639034220L, 293446789L, 
    125599618L, -793922277L, -1570282080L, -1298824927L, 443028590L, 
    -2100301289L, 1105659756L, -1999803011L, 1754072602L, -396938029L, 
    -1336474632L, -1529231975L, -687538042L, 1888368975L, -1810110652L
  )


##############################################
##############################################
######New functions added and modified #######
##############################################
##############################################


## RWZaccept_modified


RWZaccept_modified <- function (a, nvar = nvar, zero = FALSE, constrained, impulses, FEVD_check, fevd0) # constrained prendra les matrices spécifiées + rajouter un test de FEVDs
  
{
  if(!zero) { # Sign restrictions only
    #   qr_object <- qr(matrix(rnorm(M ^ 2, 0, 1), M, M))
    #   Q <- qr.Q(qr_object)
    Q <- qr.Q(a) 
    R <- qr.R(a)
  } else { # if addition of Zero restrictions (Arias, 2018)
    Q <- matrix(0, nvar, nvar)
    for(i in seq_len(nvar)) { # Build up Q
      slct_row <- which(sign_restr[, i] == 0) ## A modifier en f° de comment on écrit les restrictions
      R <- rbind(t(swish)[slct_row, ], Q[seq_len(i - 1), ]) # A priori, remplacer sigma_chol par t(swish)
      qr_object <- qr(t(R))
      qr_rank <- qr_object[["rank"]]
      set <- if(qr_rank == 0) {seq_len(M)} else {-seq_len(qr_rank)}
      N_i <- qr.Q(qr_object, complete = TRUE)[, set, drop = FALSE]
      N_stdn <- crossprod(N_i, rnorm(M, 0, 1))
      q_i <- N_i %*% (N_stdn / norm(N_stdn, type = "2"))
      Q[i, ] <- q_i
    }
    Q <- t(Q)
  }
  
  
  ###### New code for checking the sign restrictions
  
  
  ## Verify sign restrictions
  
  count <- 0
  
  for (i in 1:length(constrained)) {
    
    shock = t(impulses[j,,] %*% Q)
    shock[abs(shock) < 1e-10] <- 0
    
    shock_vec <- as.vector(shock)
    shock_vec[which(shock_vec < 0)] <- -1
    shock_vec[which(shock_vec > 0)] <- 1
    
    sign_vec <- as.vector(sign_restr[[i]])
    
    if(identical(shock_vec[restricted], sign_vec[restricted])) {count = count + 1}
  }
  
  ## Code for checking FEVDs
  
  if (!is.null(FEVD_check)) {
    
    count_FEVD = 0
    
    for (i in 1:length(FEVD_check)) {
      
      fevd = t(fevd0[i, , ] %*% (Q^2))
      
      selected_column = which(apply(FEVD_check[[i]], 2, function(col) any(col == 1)))
      fevd_col = fevd[,selected_column]
      
      selected_row = which(apply(selected_column, 1, function(row) all(row %in% c(0, 1))))
      fevd_row = fevd_col[selected_row,]
      
      ref_row = which(fevd_row == 1)
      
      
      if (fevd_row[ref_row] == max(fevd_row)) {count_FEVD = count_FEVD + 1}
    }
  }
  
  ## Return the different results needed for the rest of the estimation
  
  if (!is.null(FEVD_check)) {
    if (count == length(constrained) & count_FEVD == length(FEVD_check)) {acc = 1}
    else {acc = 0}}
  else {if (count == length(constrained)) {acc = 1}
    else {acc = 0}}
  
  rwz <- list(Q = Q, acc = acc)
  return(rwz)
}
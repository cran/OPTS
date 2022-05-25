opts_th <- function(X, Y, m, crit = "aic", type = "binseg", prop_split = 0.5, prop_trim = 0.2, q_tail = 0.5, ...) {

  n <- nrow(X)
  p <- ncol(X)
  betahat <- SE <- rep(NA, p + 1)
  Jhat <- rep(NA, p)
  pval <- rep(NA, p)
  pval_split <- matrix(NA, nrow = m, ncol = p)
  pval_max <- pval_min <- rep(NA, m)
  
  for(j in 1 : p){
    
    glm.outj <- glm(Y ~ X[ , j], ...)
    pval[j] <- summary(glm.outj)$coef[2, 4]
  }

  if(crit == "aic"){

    aic_split <- matrix(NA, nrow = m, ncol = p + 1)
    aic_split_int <- matrix(NA, nrow = m, ncol = m * p)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pval_split[i, j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pval_min[i] <- min(pval_split[i, ])
      pval_max[i] <- max(pval_split[i, ])
      pval_sort <- sort(pval_split[i, ])
    
      glm.out0 <- glm(Y_split ~ 1, ...)
      aic_split[i, 1] <- AIC(glm.out0)
      
      for (k in 1 : p) {

        Jhatk <- pval_split[i, ] <= pval_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        aic_split[i, k + 1] <- AIC(glm.outk)
      }
    
    }
    
    cutpoints <- sort(pval_split)
    
    for(i in 1:m){
      
      cuts <- sort(pval_split[i, ])
      
      if(min(cutpoints) == cuts[1]){
        
        aic_split_int[i, cutpoints == cuts[1]] <- aic_split[i, 1]
        aic_split_int[i, cuts[1] < cutpoints & cutpoints < cuts[2]] <- aic_split[i, 2]
        
        for(k in 2 : (p - 1)){
          
          aic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            aic_split[i, k + 1]
        }
      } else{
        
        aic_split_int[i, cutpoints < cuts[1]] <- aic_split[i, 1]
        for(k in 1 : (p - 1)){
          
          aic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            aic_split[i, k + 1]
        }
        
      }
      
      aic_split_int[i, cuts[p] <= cutpoints] <- aic_split[i, p + 1]
    }
    
    lower_quant <- quantile(pval_min, probs = q_tail)
    upper_quant <- quantile(pval_max, probs = 1 - q_tail)
    subset_quant <- lower_quant <= cutpoints & cutpoints <= upper_quant
    
    if(sum(subset_quant) >= 1){
      
      cutpoints <- cutpoints[subset_quant]
      aic_split_int <- aic_split_int[, subset_quant]
    }
    
    aic_mean <- apply(aic_split_int, 2, mean, trim = prop_trim, na.rm = TRUE)
    minaic.out <- sort(aic_mean, index.return = TRUE)
      
    # min AIC rule
    if(type == "min"){
      cuthat <- cutpoints[minaic.out$ix[1]]
    }
      
    # min AIC + 0.5 SD rule
    if(type == "sd"){
      minsd <- min(aic_mean) + 0.5 * sd(aic_split_int[, minaic.out$ix[1]])
      cuthat <- cutpoints[aic_mean <= minsd][1]
    }
    
    # pelt
    if(type == "pelt"){
      pelt.out <- cpt.meanvar(aic_mean, method = "PELT")
      cpts_pelt <- cpts(pelt.out)
      cpts_neigh <- c(cpts_pelt[1], cpts_pelt[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # binseg
    if(type == "binseg"){
      bs.out <- cpt.meanvar(aic_mean, method = "BinSeg")
      cpts_bs <- cpts(bs.out)
      cpts_neigh <- c(cpts_bs[1], cpts_bs[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # amoc
    if(type == "amoc"){
      amoc.out <- cpt.meanvar(aic_mean, method = "AMOC")
      cpts_amoc <- cpts(amoc.out)
      cpts_neigh <- c(cpts_amoc, cpts_amoc + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
    
    Jhat <- pval <= cuthat
    
    if(sum(Jhat) >= 1){
      
      XB <- X[, Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    } else{
      
      glm.out <- glm(Y ~ 1, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
    output <- list(betahat = betahat, Jhat = Jhat, SE = SE, cuthat = cuthat, 
      pval = pval, cutpoints = cutpoints, aic_mean = aic_mean)
    
  }
  
  if(crit == "bic"){

    bic_split <- matrix(NA, nrow = m, ncol = p + 1)
    bic_split_int <- matrix(NA, nrow = m, ncol = m * p)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pval_split[i, j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pval_min[i] <- min(pval_split[i, ])
      pval_max[i] <- max(pval_split[i, ])
      pval_sort <- sort(pval_split[i, ])
    
      glm.out0 <- glm(Y_split ~ 1, ...)
      bic_split[i, 1] <- BIC(glm.out0)
      
      for (k in 1 : p) {

        Jhatk <- pval_split[i, ] <= pval_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        bic_split[i, k + 1] <- BIC(glm.outk)
      }
    
    }
    
    cutpoints <- sort(pval_split)
    
    for(i in 1:m){
      
      cuts <- sort(pval_split[i, ])
      
      if(min(cutpoints) == cuts[1]){
        
        bic_split_int[i, cutpoints == cuts[1]] <- bic_split[i, 1]
        bic_split_int[i, cuts[1] < cutpoints & cutpoints < cuts[2]] <- bic_split[i, 2]
        
        for(k in 2 : (p - 1)){
          
          bic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            bic_split[i, k + 1]
        }
      } else{
        
        bic_split_int[i, cutpoints < cuts[1]] <- bic_split[i, 1]
        for(k in 1 : (p - 1)){
          
          bic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            bic_split[i, k + 1]
        }
        
      }
      
      bic_split_int[i, cuts[p] <= cutpoints] <- bic_split[i, p + 1]
    }
    
    lower_quant <- quantile(pval_min, probs = q_tail)
    upper_quant <- quantile(pval_max, probs = 1 - q_tail)
    subset_quant <- lower_quant <= cutpoints & cutpoints <= upper_quant
    
    if(sum(subset_quant) >= 1){
      
      cutpoints <- cutpoints[subset_quant]
      bic_split_int <- bic_split_int[, subset_quant]
    }
    
    bic_mean <- apply(bic_split_int, 2, mean, trim = prop_trim, na.rm = TRUE)
    minbic.out <- sort(bic_mean, index.return = TRUE)
      
    # min BIC rule
    if(type == "min"){
      cuthat <- cutpoints[minbic.out$ix[1]]
    }
      
    # min BIC + 0.5 SD rule
    if(type == "sd"){
      minsd <- min(bic_mean) + 0.5 * sd(bic_split_int[, minbic.out$ix[1]])
      cuthat <- cutpoints[bic_mean <= minsd][1]
    }
    
    # pelt
    if(type == "pelt"){
      pelt.out <- cpt.meanvar(bic_mean, method = "PELT")
      cpts_pelt <- cpts(pelt.out)
      cpts_neigh <- c(cpts_pelt[1], cpts_pelt[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
    
    # binseg
    if(type == "binseg"){
      bs.out <- cpt.meanvar(bic_mean, method = "BinSeg")
      cpts_bs <- cpts(bs.out)
      cpts_neigh <- c(cpts_bs[1], cpts_bs[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # amoc
    if(type == "amoc"){
      amoc.out <- cpt.meanvar(bic_mean, method = "AMOC")
      cpts_amoc <- cpts(amoc.out)[1]
      cpts_neigh <- c(cpts_amoc, cpts_amoc + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    Jhat <- pval <= cuthat
    
    if(sum(Jhat) >= 1){
      
      XB <- X[, Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    } else{
      
      glm.out <- glm(Y ~ 1, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[,1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
    output <- list(betahat = betahat, Jhat = Jhat, SE = SE, cuthat = cuthat, 
      pval = pval, cutpoints = cutpoints, bic_mean = bic_mean)
  }
  
  output
}

opts <- function(X, Y, m, crit = "aic", prop_split = 0.5, cutoff = 0.75, ...) {

  n <- nrow(X)
  p <- ncol(X)
  betahat <- SE <- rep(NA, p + 1)
  pvals <- rep(NA, p)
  Jhat_split <- matrix(NA, nrow = m, ncol = p)

  if(crit == "aic"){

    aic_split <- rep(NA, p + 1)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pvals[j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pvals_sort <- sort(pvals)
      
      glm.out0 <- glm(Y_split ~ 1, ...)
      aic_split[1] <- AIC(glm.out0)
        
      for (k in 1 : p) {

        Jhatk <- pvals <= pvals_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        aic_split[k + 1] <- AIC(glm.outk)
      }
    
      aic_sort <- sort(aic_split[-1], index.return = TRUE)
      idx_min <- aic_sort$ix[1]
      aic_sort0 <- sort(aic_split, index.return = TRUE)
      idx_min0 <- aic_sort0$ix[1]
      
      if(idx_min0 != 1){
        
        Jhat_split[i, ] <- pvals <= pvals_sort[idx_min]
        
      } else {
        
        Jhat_split[i, ] <- rep(FALSE, p)
      }

    }

    freqs <- apply(Jhat_split, 2, mean)
    Jhat <- freqs >= cutoff

    if(sum(Jhat) >= 1){
    
      XB <- X[ , Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
    
    } else{
    
      glm.out <- glm(Y ~ 1, ...)
      Jhat <- rep(FALSE, p)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
  }
  
  if(crit == "bic"){

    bic_split <- rep(NA, p + 1)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pvals[j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pvals_sort <- sort(pvals)
      
      glm.out0 <- glm(Y_split ~ 1, ...)
      bic_split[1] <- BIC(glm.out0)
        
      for (k in 1 : p) {

        Jhatk <- pvals <= pvals_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        bic_split[k + 1] <- BIC(glm.outk)
      }
    
      bic_sort <- sort(bic_split[-1], index.return = TRUE)
      idx_min <- bic_sort$ix[1]
      bic_sort0 <- sort(bic_split, index.return = TRUE)
      idx_min0 <- bic_sort0$ix[1]
      
      if(idx_min0 != 1){
        
        Jhat_split[i, ] <- pvals <= pvals_sort[idx_min]
        
      } else {
        
        Jhat_split[i, ] <- rep(FALSE, p)
      }

    }

    freqs <- apply(Jhat_split, 2, mean)
    Jhat <- freqs >= cutoff

    if(sum(Jhat) >= 1){
    
      XB <- X[ , Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
    
    } else{
    
      glm.out <- glm(Y ~ 1, ...)
      Jhat <- rep(FALSE, p)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
  }
    
  list(betahat = betahat, Jhat = Jhat, SE = SE, freqs = freqs)
}

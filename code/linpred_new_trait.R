linpred_new_trait <- function(hM, Tr, X, C = NULL, ...) {
    
    require(mvtnorm)
    
    # some checks
    # FIX ME
    
    # pool posteriors
    draws <- poolMcmcChains(hM$postList)
    
    # collect covariates and indices
    ns <- nrow(Tr)
    nc <- hM$nc
    k <- length(draws)
    
    # phylogenetic covariance
    if (is.null(C)) {
        C <- diag(ns)
    }
    
    # mean of Beta
    mu_Beta <- 
        lapply(draws, function(x) {
            as.vector(x$Gamma %*% t(Tr))
        })
    
    # variance of Beta
    if (is.null(C)) {
        var_Beta <- 
            lapply(draws, function(x) {
                kronecker(diag(ns), x$V) 
            })
    } else {
        var_Beta <- 
            lapply(draws, function(x) {
                kronecker(C * x$rho + (1 - x$rho) * diag(ns), x$V) 
            })
    }
    
    # draw Beta
    Beta <- vector("list", k)
    for (i in 1:k) {
        Beta[[i]] <- rmvnorm(1, mu_Beta[[i]], var_Beta[[i]])
        Beta[[i]] <- matrix(Beta[[i]], nc, ns)
    }
    
    # calculate fixed effect
    LF <- lapply(Beta, function(b) X %*% b)
    
    # calculate random effect
    # not implement for now (i.e., only doing fixed-effect prediction)
    LR <- 0
    
    # calculate linear predictor
    L <- vector("list", k)
    for (i in 1:k) {
        L[[i]] <- LF[[i]] + LR
    }
    
    return(list(Beta = Beta, L = L))
    
}
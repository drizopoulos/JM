marginal_coefs <- function (object, std_errors = FALSE, time = NULL, M = 5000, L = 500) {
    
    if (!object$parameterization == "value") {
        stop("parameterization must be 'value'.\n")
    }
    
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsT <- object$termsT
    assignT <- object$assignT
    
    data <- rbind(object$data, object$data.id)
    
    if (!(all(names(assignT) %in% colnames(data) == TRUE))){
        ind <- which(!(names(assignT) %in% colnames(data)))
        data <- rbind(cbind(object$data, name = rep(object$x$W[,ind], count(object$id)$freq)),
                      cbind(object$data.id, name = object$x$W[,ind]))
        colnames(data)[which(colnames(data) == "name")] <- names(assignT)[ind]
    }
    
    timeVar <- object$timeVar
    max.t <- round(max(object$data[[timeVar]]),0)
    
    if (!is.null(time)) {
        data[[timeVar]] <- time
    }
    
    betas <- object$coefficients[['betas']]
    D <- object$coefficients[['D']]
    gammas <- object$coefficients[['gammas']]
    gammas.bs <- object$coefficients[['gammas.bs']]
    alpha <- object$coefficients[['alpha']]
    
    if (!is.null(gammas)) {
        W <- model.matrix(delete.response(TermsT), data = data)[, -1, drop = FALSE]
        X_tilde <- W
    }
    
    X <- model.matrix(TermsX, data = data)
    X_tilde <- cbind(X_tilde, X)
    const <- c(X_tilde %*% c(gammas, alpha * betas))
    Z <- model.matrix(TermsZ, data = data)
    nRE <- ncol(Z)
    n <- nrow(data)
    log_marg_hazard_ratio <- NULL
    log_marg_hazard_ratio.res <- NULL
    id.ind <- c(object$id, unique(object$id))
    
    for (i in unique(object$id)) {
        b <- mvrnorm(M, rep(0.0, nRE), D)
        rows <- which(id.ind == i)
        dat.subj <- data[rows, ]
        Z.subj <- model.matrix(TermsZ, data = dat.subj)
        Zb.subj <- Z.subj %*% t(b)
        const.subj <- matrix(rep(const[rows], each = M), nrow = length(rows), ncol = M, byrow = TRUE)
        ss_hazard_ratio <- exp(const.subj + alpha * Zb.subj)
        log_marg_hazard_ratio.res  <- log(rowMeans(ss_hazard_ratio))
        log_marg_hazard_ratio[rows] <- log_marg_hazard_ratio.res
        log_marg_hazard_ratio.res <- NULL
    }
    
    ss_coefs <- c(gammas, alpha * betas)
    marg_coefs <- c(ginv(crossprod(X_tilde)) %*% crossprod(X_tilde, log_marg_hazard_ratio))
    
    ses <- NULL
    if (std_errors) {
        D <- object$coefficients[['D']]
        diag_D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        list_thetas <- list(betas = betas, gammas= gammas, alpha = alpha, D = if (diag_D) log(diag(D)) else chol_transf(D))
        tht <- unlist(as.relistable(list_thetas))
        index <- c(1:length(betas), 
                   (length(betas)+2):(length(betas)+2+length(gammas)), 
                   (nrow(vcov(object))-length(list_thetas$D) + 1):nrow(vcov(object)))
        V <- vcov(object)[index, index]
        marg_coefs_star <- matrix(0, nrow = L, ncol = ncol(X_tilde))
        
        for (l in 1:L) {
            new_tht <- relist(MASS::mvrnorm(1, tht, V), skeleton = list_thetas)
            new_betas <- new_tht$betas
            new_gammas <- new_tht$gammas
            new_alpha <- new_tht$alpha
            new_D <- if (diag_D) diag(exp(new_tht$D), length(new_tht$D)) else chol_transf(new_tht$D)
            new_const <- c(X_tilde %*% c(new_gammas, new_alpha * new_betas))
            
            log_marg_hazard_ratio_star <- NULL
            log_marg_hazard_ratio.res.star <- NULL
            id.ind <- c(object$id, unique(object$id))
            
            for (i in unique(object$id)) {
                b <- mvrnorm(M, rep(0.0, nRE), new_D)
                rows <- which(id.ind == i)
                dat.subj <- data[rows, ]
                Z.subj <- model.matrix(TermsZ, data = dat.subj)
                Zb.subj <- Z.subj %*% t(b)
                const.subj <- matrix(rep(new_const[rows], each = M), nrow = length(rows), ncol = M, byrow = TRUE)
                ss_hazard_ratio_star <- exp(const.subj + new_alpha * Zb.subj)
                
                log_marg_hazard_ratio.res.star  <- log(rowMeans(ss_hazard_ratio_star))
                log_marg_hazard_ratio_star[rows] <- log_marg_hazard_ratio.res.star
                log_marg_hazard_ratio.res.star <- NULL
            }
            marg_coefs_star[l, ] <- c(ginv(crossprod(X_tilde)) %*% crossprod(X_tilde, log_marg_hazard_ratio_star))
        }
        ses <- sqrt(diag(var(marg_coefs_star))) 
    }
    
    if (std_errors) { list(ss_coefs = ss_coefs, marg_coefs = marg_coefs, ses = ses) }
    else {
        list(ss_coefs = ss_coefs, marg_coefs = marg_coefs)
    }
}


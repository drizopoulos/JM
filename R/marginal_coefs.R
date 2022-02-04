chol_transf <- function (x) {
    if (any(is.na(x) | !is.finite(x)))
        stop("NA or infinite values in 'x'.\n")
    if (is.matrix(x)) {
        k <- nrow(x)
        U <- chol(x)
        U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
        U[upper.tri(U, TRUE)]
    } else {
        nx <- length(x)
        k <- round((-1 + sqrt(1 + 8 * nx))/2)
        mat <- matrix(0, k, k)
        mat[upper.tri(mat, TRUE)] <- x
        mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
        res <- crossprod(mat)
        attr(res, "L") <- t(mat)[lower.tri(mat, TRUE)]
        res
    }
}

marginal_coefs <- function (object, std_errors = FALSE, M = 2000, L = 200) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!object$parameterization == "value") {
        stop("parameterization must be 'value'.\n")
    }
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsT <- object$termsT
    assignT <- object$assignT
    timeVar <- object$timeVar
    data <- rbind(object$data, object$data.id)
    
    if (!(all(names(assignT) %in% colnames(data) == TRUE))){
        ind <- which(!(names(assignT) %in% colnames(data)))
        data <- rbind(cbind(object$data, name = rep(object$x$W[,ind], table(object$id))),
                      cbind(object$data.id, name = object$x$W[,ind]))
        colnames(data)[which(colnames(data) == "name")] <- names(assignT)[ind]
    }
    
    times <- data[[timeVar]]
    nsubj <- length(unique(object$id))
    max.t <- round(max(data[[timeVar]]),0)
    betas <- object$coefficients[['betas']]
    D <- object$coefficients[['D']]
    gammas <- object$coefficients[['gammas']]
    gammas.bs <- object$coefficients[['gammas.bs']]
    alpha <- object$coefficients[['alpha']]
    ordSpline = object$control$ord
    knots = object$control$knots$`1`
    W0 <- splineDesign(knots, data[[timeVar]], ord = ordSpline, outer.ok = TRUE)
    W <- model.matrix(delete.response(TermsT), data = data)[, -1, drop = FALSE]
    X <- model.matrix(TermsX, data = data)
    X_tilde2 <- cbind(W0, W, X %*% betas)
    X_tilde_use <- cbind(W0, W, X)
    const <- c(X_tilde_use %*% c(gammas.bs, gammas, alpha * betas))
    Z <- model.matrix(TermsZ, data = data)
    nRE <- ncol(Z)
    n <- nrow(data)
    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    K <- length(sk)
    P <- times/2
    st <- outer(P, sk + 1)
    W0s <- splineDesign(knots, c(t(st)), ord = ordSpline,  outer.ok = TRUE)
    id.GK <- rep(seq_along(times), each = K)
    data2 <- data[id.GK, ]
    data2[[timeVar]] <- c(t(st))
    Xs <- model.matrix(TermsX, data = data2)
    Zs <- model.matrix(TermsZ, data = data2)
    
    surv.part      <- matrix(0, nrow = K, ncol = M)
    surv.res_num   <- matrix(0, nrow = nrow(data), ncol = M)
    surv.res_denum <- matrix(0, nrow = nrow(data), ncol = M)
    marg_hazard <- NULL
    id.ind <- c(object$id, unique(object$id))
    
    for (i in unique(object$id)) {
        b <- mvrnorm(M, rep(0.0, nRE), D)
        rows <- which(id.ind == i)
        dat.subj <- data[rows, ]
        Z.subj <- model.matrix(TermsZ, data = dat.subj)
        Zb.subj <- Z.subj %*% t(b)
        const.subj <- matrix(rep(const[rows], each = M), nrow = length(rows), ncol = M, byrow = TRUE)
        ss_hazard <- exp(const.subj + alpha * Zb.subj)
        
        for (r in rows) {
            for (k in 1:K) {
                log.h0.s <- W0s[K * (r - 1) + k, ] %*% gammas.bs
                aXB <- alpha * (Xs[K * (r - 1) + k, ] %*% betas) 
                const.subj2 <- rep((log.h0.s + aXB), M)
                Zbs.subj <- Zs[K * (r - 1) + k, ] %*% t(b) 
                surv.part[k, ] <- wk[k] * exp(const.subj2 + alpha * Zbs.subj) 
            }
            surv.res_num[r, ] <- colSums(surv.part)  
        }
        ss_surv_num <- exp(c(-exp(W[rows,] %*% gammas) * P[rows]) * surv.res_num[rows,])
        
        b <- mvrnorm(M, rep(0.0, nRE), D)
        for (r in rows) {
            for (k in 1:K) {
                log.h0.s <- W0s[K * (r - 1) + k, ] %*% gammas.bs
                aXB <- alpha * (Xs[K * (r - 1) + k, ] %*% betas) 
                const.subj2 <- rep((log.h0.s + aXB), M)
                Zbs.subj <- Zs[K * (r - 1) + k, ] %*% t(b) 
                surv.part[k, ] <- wk[k] * exp(const.subj2 + alpha * Zbs.subj) 
            }
            surv.res_denum[r, ] <- colSums(surv.part) 
        }
        ss_surv_denum <- exp(c(-exp(W[rows,] %*% gammas) * P[rows]) * surv.res_denum[rows,])
        marg_haz  <- rowMeans(ss_hazard * ss_surv_num) / rowMeans(ss_surv_denum)
        marg_hazard[rows] <- marg_haz
        marg_haz <- NULL
    }
    
    log_marg_hazard <- log(marg_hazard)
    ss_coefs <- c(gammas.bs, gammas, alpha)
    marg_coefs <- solve(t(X_tilde2)%*%X_tilde2) %*% t(X_tilde2) %*% log_marg_hazard
    marg_coefs <- c(marg_coefs)
    
    ses <- NULL
    if (std_errors) {
        D <- object$coefficients[['D']]
        diag_D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        list_thetas <- list(betas = betas, gammas= gammas, alpha = alpha, 
                            gammas.bs = gammas.bs, 
                            D = if (diag_D) log(diag(D)) else chol_transf(D))
        tht <- unlist(as.relistable(list_thetas))
        index <- c(1:length(betas), 
                   (length(betas)+2):(length(betas)+2+length(gammas)), 
                   (length(betas)+3+length(gammas)):(nrow(vcov(object))-length(list_thetas$D)),
                   (nrow(vcov(object))-length(list_thetas$D) + 1):nrow(vcov(object)))
        V <- vcov(object)[index, index]
        marg_coefs_star <- matrix(0, nrow = L, ncol = length(marg_coefs))
        timepnts <- seq(0, max.t, by = 0.01)
        
        for (l in 1:L) {
            new_tht <- relist(MASS::mvrnorm(1, tht, V), skeleton = list_thetas)
            new_betas <- new_tht$betas
            new_gammas <- new_tht$gammas
            new_alpha <- new_tht$alpha
            new_gammas.bs <- new_tht$gammas.bs
            new_D <- if (diag_D) diag(exp(new_tht$D), length(new_tht$D)) else chol_transf(new_tht$D)
            new_const <- c(X_tilde_use %*% c(new_gammas.bs, new_gammas, new_alpha * new_betas))
            X_tilde2_new <- cbind(W0, W, X %*% new_betas)
            log_marg_hazard.star <- NULL
            
            for (i in unique(object$id)) {
                b <- mvrnorm(M, rep(0.0, nRE), D)
                rows <- which(id.ind == i)
                dat.subj <- data[rows, ]
                Z.subj <- model.matrix(TermsZ, data = dat.subj)
                Zb.subj <- Z.subj %*% t(b)
                const.subj <- matrix(rep(new_const[rows], each = M), nrow = length(rows), ncol = M, byrow = TRUE)
                ss_hazard <- exp(const.subj + new_alpha * Zb.subj)
                
                for (r in rows) {
                    for (k in 1:K) {
                        log.h0.s <- W0s[K * (r - 1) + k, ] %*% new_gammas.bs
                        aXB <- new_alpha * (Xs[K * (r - 1) + k, ] %*% new_betas) 
                        const.subj2 <- rep((log.h0.s + aXB), M)
                        Zbs.subj <- Zs[K * (r - 1) + k, ] %*% t(b) 
                        surv.part[k, ] <- wk[k] * exp(const.subj2 + new_alpha * Zbs.subj) 
                    }
                    surv.res_num[r, ] <- colSums(surv.part)  
                }
                ss_surv_num <- exp(c(-exp(W[rows,] %*% new_gammas) * P[rows]) * surv.res_num[rows,])
                
                b <- mvrnorm(M, rep(0.0, nRE), D)
                for (r in rows) {
                    for (k in 1:K) {
                        log.h0.s <- W0s[K * (r - 1) + k, ] %*% new_gammas.bs
                        aXB <- new_alpha * (Xs[K * (r - 1) + k, ] %*% new_betas) 
                        const.subj2 <- rep((log.h0.s + aXB), M)
                        Zbs.subj <- Zs[K * (r - 1) + k, ] %*% t(b) 
                        surv.part[k, ] <- wk[k] * exp(const.subj2 + new_alpha * Zbs.subj) 
                    }
                    surv.res_denum[r, ] <- colSums(surv.part) 
                }
                ss_surv_denum <- exp(c(-exp(W[rows,] %*% new_gammas) * P[rows]) * surv.res_denum[rows,])
                marg_haz.star  <- rowMeans(ss_hazard * ss_surv_num) / rowMeans(ss_surv_denum)
                log_marg_hazard.star[rows] <- log(marg_haz.star)
                marg_haz.star <- NULL
            }
            marg_coefs_star[l, ] <- c(solve(t(X_tilde2_new)%*%X_tilde2_new) %*% t(X_tilde2_new) %*% log_marg_hazard.star)
        }
        ses <- sqrt(diag(var(marg_coefs_star, na.rm = TRUE))) 
    }
    
    if (std_errors) { list(ss_coefs = c(ss_coefs, betas), marg_coefs = c(marg_coefs, betas), ses = ses) }
    else {
        list(ss_coefs = c(ss_coefs, betas), marg_coefs = c(marg_coefs, betas))
    }
}

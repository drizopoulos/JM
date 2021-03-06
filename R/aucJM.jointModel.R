aucJM.jointModel <-
function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, idVar = "id", 
                           simulate = FALSE, M = 100, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$termsT
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    timeVar <- object$timeVar
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- newdata2[Time[ordTime] > Tstart, ]
    newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]
    pi.u.t <- survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                        simulate = simulate, M = M)
    pi.u.t <- sapply(pi.u.t$summaries, "[", 1, 2)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[!duplicated(id), 1]
    event <- SurvT[!duplicated(id), 2]
    names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + 1e-07
    }
    if (!all(names(pi.u.t) == names(Time)))
        stop("...")
    auc <- if (length(Time) > 1) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        pi.u.t.i <- pi.u.t[pairs[1, ]]
        pi.u.t.j <- pi.u.t[pairs[2, ]]
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & di == 0) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
        ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- paste(names(Ti), names(Tj), sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            pi2 <- survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], idVar = idVar, 
                             last.time = Time[unq_nams_i], survTimes = Thoriz, 
                             simulate = simulate, M = M)
            pi2 <- 1 - sapply(pi2$summaries, "[", 1, 2)
            ind[ind2] <- ind[ind2] * pi2[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            pi3 <- survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], idVar = idVar, 
                             last.time = Time[unq_nams_j], survTimes = Thoriz, 
                             simulate = simulate, M = M)
            pi3 <- sapply(pi3$summaries, "[", 1, 2)
            ind[ind3] <- ind[ind3] * pi3[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            pi4_i <- survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], idVar = idVar, 
                               last.time = Time[unq_nams_i], survTimes = Thoriz, 
                               simulate = simulate, M = M)
            pi4_i <- 1 - sapply(pi4_i$summaries, "[", 1, 2)
            pi4_j <- survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], idVar = idVar, 
                             last.time = Time[unq_nams_j], survTimes = Thoriz, 
                             simulate = simulate, M = M)
            pi4_j <- sapply(pi4_j$summaries, "[", 1, 2)
            ind[ind4] <- ind[ind4] * pi4_i[nams_i] * pi4_j[nams_j]
        }
        sum((pi.u.t.i < pi.u.t.j) * ind, na.rm = TRUE) / sum(ind, na.rm = TRUE)
    } else {
        NA
    }
    out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)), 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "aucJM"
    out
}

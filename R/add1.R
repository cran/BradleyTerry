"add1.BTm" <-
    function (object, scope, scale = 0,
              test = c("none", "Chisq", "F"),
              x = NULL, k = 2, ...) 
{
    Fstat <- function(table, rdf) {
        dev <- table$Deviance
        df <- table$Df
        diff <- pmax(0, (dev[1] - dev)/df)
        Fs <- (diff/df)/(dev/(rdf - df))
        Fs[df < .Machine$double.eps] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
        P[nnas] <- pf(Fs[nnas], df[nnas], rdf - df[nnas],
                      lower.tail = FALSE)
        list(Fs = Fs, P = P)
    }
    formula0 <- formula(object)
    if (formula0[[3]] == "..") 
        stop("no terms to add -- the model is already saturated")
    if (!is.character(scope)) 
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope)) 
        stop("no terms in scope for adding to object")
    oTerms <- attr(object$terms, "term.labels")
    ns <- length(scope)
    dfs <- dev <- numeric(ns + 1)
    names(dfs) <- names(dev) <- c("<none>", scope)
    dfs[1] <- object$rank
    dev[1] <- object$deviance
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(formula0, add.rhs)
    new.form[[2]] <- NULL
    Terms <- terms(new.form)
    y <- object$y0
    if (is.null(x)) {
        fc <- object$call
        fc$formula <- Terms
        data.arg <- object$call$data
        m <- if (is.null(data.arg)){
               model.frame(fc$formula, xlev = object$xlevels,
                           na.action = I)}
             else {
               model.frame(fc$formula, xlev = object$xlevels,
                           na.action = I,
                           data = eval(data.arg, parent.frame()))}
        x <- model.matrix(Terms, m, contrasts = object$contrasts)
        rownames(x) <- rownames(object$x0)
        temp <- attr(x, "assign")
    }
    n <- nrow(x)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0) > 0
    control <- object$control
    br <- object$bias.reduction
    if (is.null(br)) 
        br <- FALSE
    for (tt in scope) {
        usex <- match(asgn, match(tt, Terms), 0) > 0
        X <- x[order(rownames(x)), usex | ousex, drop = FALSE]
        z <- BTm(y ~ X,
                 offset = object$offset0,
                 control = object$control,
                 br = br,
                 order.effect = object$order.effect)
        dfs[tt] <- z$rank
        dev[tt] <- z$deviance
    }
    if (scale == 0) 
        dispersion <- summary(object,
                              dispersion = NULL)$dispersion
    else dispersion <- scale
    fam <- object$family$family
    if (fam == "gaussian") {
        if (scale > 0) 
            loglik <- dev/scale - n
        else loglik <- n * log(dev/n)
    }
    else loglik <- dev/dispersion
    aic <- loglik + k * dfs
    aic <- aic + (extractAIC(object, k = k)[2] - aic[1])
    dfs <- dfs - dfs[1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic,
                      row.names = names(dfs), 
                      check.names = FALSE)
    if (all(is.na(aic))) 
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- pmax(0, loglik[1] - loglik)
        dev[1] <- NA
        LRT <- if (dispersion == 1) 
            "LRT"
        else "scaled dev."
        aod[, LRT] <- dev
        nas <- !is.na(dev)
        dev[nas] <- pchisq(dev[nas], aod$Df[nas],
                           lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson") 
            warning(paste("F test assumes quasi", fam, " family", 
                          sep = ""))
        rdf <- object$df.residual
        aod[, c("F value", "Pr(F)")] <- Fstat(aod, rdf)
    }
    head <- c("Single term additions", "\nModel:",
              deparse(formula0), 
              if (scale > 0) paste("\nscale: ",
                                   format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

BTm <- function(formula, refcat = NULL, offset = NULL,
                contrasts = NULL, data = NULL,
                subset = NULL, br = FALSE, order.effect = NULL,
                ...){
## First define some utility functions
    expand.col <- function(col){
        col[is.na(col)] <- 0
        mat <- outer(col, col, "-")
        unravel(mat)}
    BTexpand <- function(design){
        if (ncol(design) == 0){
            matrix(, nrow(design)*(nrow(design)-1)/2, 0)
        }  else apply(design, 2, expand.col)}
    unravel <- function(matrix) matrix[row(matrix) > col(matrix)]
## Decode the function call
    if (!is.null(refcat)) refcat <- make.names(refcat)
    dots <- match.call(expand.dots = FALSE)$...
    if (!is.null(dots$weights))
        stop("weights argument is not implemented in BTm")
    if (!is.null(subset) &&
        !is.logical(subset) &&
        !is.numeric(subset))
        stop("subset argument must be a logical or numeric vector")
    m <- match.call(expand.dots = FALSE)
    m$formula <- theformula <- eval(m$formula, parent.frame())
    if (!inherits(theformula, "formula")) stop("invalid formula argument")
    m[[1]] <- as.name("model.frame")
    m$na.action <- I
    ymat <- eval(m$formula[[2]], parent.frame())
    if (!is.null(order.effect) &&
        length(order.effect) != nrow(ymat)){
        stop("order.effect has the wrong length")
    }
    m$formula[[2]] <- NULL
    if ("freq" %in% colnames(ymat) && !("Freq" %in% colnames(ymat))) {
        warning("Data column 'freq' ignored.  Use 'Freq' as the name for a column of frequencies, if that is what you intended")
    }
    freq <- {if ("Freq" %in% colnames(ymat)) ymat[, "Freq"]
                 else rep(1, nrow(ymat))}
    if (is.logical(subset) && (length(subset) != nrow(ymat)))
        stop("subset argument has wrong length")
    if (is.numeric(subset) && any(!(subset %in% 1:nrow(ymat))))
        stop("subset values not all valid")
    if (is.numeric(subset)) subset <- (1:nrow(ymat)) %in% subset
    if (is.null(subset)) subset <- rep(TRUE, nrow(ymat))
    valid.contest <- subset & !is.na(freq)
    freq <- freq[valid.contest]
    if (!is.numeric(freq)) stop("Freq column should be numeric")
    if (!is.null(order.effect))
        order.effect <- order.effect[valid.contest]
    winner <- if ("winner" %in% colnames(ymat)){
                   ymat[valid.contest, "winner"]} else
                   ymat[valid.contest, 1]
    loser <- if ("loser" %in% colnames(ymat)){
                  ymat[valid.contest, "loser"]} else
                  ymat[valid.contest, 2]
    winner <- make.names(winner)
    loser <- make.names(loser)
    playernames <- sort(union(winner, loser))
    nplayers <- length(playernames)
    formulaRHS <- theformula[[3]]
    standard.BTmodel <- formulaRHS == ".." ||
                        (length(formulaRHS) > 1 &&
                         ".." %in% as.list(formulaRHS))
    if (standard.BTmodel){
        if (formulaRHS != ".."){
            warning("RHS terms other than .. ignored")
            theformula[[3]] <- as.name("..")
        }
        if (!is.null(offset)) warning("offset ignored")
        .. <- factor(playernames)
        if (!is.null(refcat)) .. <- relevel(.., refcat)
        m <- model.frame(~ ..)
        row.names(m) <- m$..
    } else {
        m$refcat <- m$... <- m$contrasts <-
            m$subset <- m$br <- m$order.effect <- NULL
        m$drop.unused.levels <- TRUE
        if (is.null(m$data)) m$data <- data.frame(playernames)
        m <- eval(m, parent.frame())
    }
    row.names(m) <- m.rownames <- make.names(row.names(m))
    if (all(playernames %in% m.rownames)){
        m <- m[playernames, , drop = FALSE]
    } else {
        if (nrow(m) != length(playernames)) stop(
                "LHS and RHS data lengths differ")
    }
    row.names(m) <- playernames
    Terms <- attr(m, "terms")
## Get info on factors, to put in the result object
    xlevels <- NULL
    if (length(attr(Terms, "factors")) > 0){
        xlevels <- (attr(Terms, "order") == 1)
        xnames <- attr(Terms, "term.labels")[xlevels]
        xlevels <- sapply(xnames, function(label)
                          is.factor(eval(parse(text = label), m)))
        xnames <- xnames[xlevels]
        xlevels <- lapply(xnames, function(label)
                          levels(eval(parse(text = label), m)))
        names(xlevels) <- xnames}
## Set up some basic design columns for possible use in case
## covariates or offset have NAs
    player.cols <- diag(1, nplayers)
    rownames(player.cols) <- colnames(player.cols) <- playernames
    colnames(player.cols) <- paste("..", colnames(player.cols),
                                   sep = "")
## Now make the Bradley-Terry design matrix.
    cols.to.add <- logical(nrow(m))
    x <- x0 <- model.matrix(Terms, m, contrasts)
    colnames(x0) <- make.names(colnames(x0))
    xcont <- attr(x, "contrasts")
    x <-  x[ , colnames(x) != "(Intercept)", drop = FALSE]
    x.has.NA <- logical(nrow(x))
    for (i in 1:nrow(x)){
        if (any(is.na(x[i, ]))){
            is.na(x[i, ]) <- rep(TRUE, ncol(x))
            x.has.NA[i] <- TRUE
        }
    }
    cols.to.add <- x.has.NA
    x <- BTexpand(x)
## Next sort out the offset, if there is one
    offset <- offset0 <- model.extract(m, offset)
    if (!is.null(offset)){
        offset.has.NA <- is.na(offset)
        cols.to.add <- cols.to.add | offset.has.NA
        offset <- expand.col(offset)
    }
    cols.to.add <- BTexpand(player.cols[, cols.to.add,
                                        drop = FALSE])
    names.to.add <- colnames(cols.to.add)
    if (is.matrix(cols.to.add)) x <- cbind(cols.to.add, x)
    rownames(x) <-
        unravel(outer(playernames, playernames,
                         function(x, y) paste(x, y, sep = " vs ")))
    if (!is.null(offset)) names(offset) <- rownames(x)
## Now set up the matrix of binomial response data
    y <- matrix(0 , length(freq), 2)
    colnames(y) <- c("won", "lost")
    winnernum <- match(winner, playernames)
    losernum <-  match(loser, playernames)
    lower.triangle <- winnernum > losernum
    rownames(y) <- ifelse(lower.triangle,
                            paste(winner, loser, sep = " vs "),
                            paste(loser, winner, sep = " vs "))
    y[, "won"] <- ifelse(lower.triangle, freq, 0)
    y[, "lost"] <- ifelse(!lower.triangle, freq, 0)
    y1 <- matrix(, 0, 3)
    colnames(y1) <- c("won", "lost", ".order")
    .order <- if (is.null(order.effect)) rep(0, length(freq)) else
               ifelse(lower.triangle, order.effect, -(order.effect))
    for (i in unique(.order)){
        yi <- y[.order == i, , drop = FALSE]
        namesi <- list(as.factor(row.names(yi)))
        yi <- cbind(tapply(yi[, "won"], namesi, sum),
                    tapply(yi[, "lost"], namesi, sum),
                    i)
        y1 <- rbind(y1, yi)
    }
    BTyvar <- y1[, c("won", "lost"), drop = FALSE]
    x <- x[rownames(y1), , drop = FALSE]
    if (!is.null(offset)) offset <- offset[rownames(y1)]
    if (!is.null(order.effect)){
        rownames(BTyvar) <- rownames(x) <-
            paste(rownames(x), y1[, 3])
        xnames <- colnames(x)
        x <- cbind(x, y1[, ".order"])
        colnames(x) <- c(xnames, ".order")
    }
    if (ncol(x) == 0) x <- cbind(x, 0)
## Next, the main part: fit the model!
    fmla <- as.formula(paste(c("BTyvar ~ 0",
                             make.names(colnames(x))),
                             collapse = "+"))
    if (!br) {fitmodel <- glm} else {
        fitmodel <- brglm:::brglm
    }
    result <- fitmodel(fmla, offset = offset, family = binomial,
                  data = data.frame(x), ...)
## Now compute abiliy estimates and standard errors for those
    na.to.zero <- function(xmat) {
        if (dim(xmat)[2] == 0) return(xmat)
        result <- apply(xmat, 1, function(row){
            if (any(is.na(row))) rep(0, length(row)) else row
        })
        result <- if (is.matrix(result)) t(result)
                  else as.matrix(result)
        dimnames(result) <- dimnames(xmat)
        return(result)
    }
    refcat <- if (is.null(refcat)) 1 else match(refcat, playernames)
    coefs <- coef(result)
    coefs <- coefs[names(coefs) != ".order"]
    x.players <- cbind(player.cols, na.to.zero(x0))[, names(coefs),
                                                    drop = FALSE]
    x.players <- x.players[, !is.na(coefs), drop = FALSE]
    abilities <- x.players %*% na.omit(coefs)
    abilities <- abilities - abilities[refcat]
    se.ability <- rep(0, length(abilities))
    if (length(coefs) > 0){
        sqrt.vcov <- chol(vcov(result))
        coef.indices <- colnames(sqrt.vcov) != ".order"
        sqrt.vcov <- sqrt.vcov[coef.indices, coef.indices]
        v <- crossprod(sqrt.vcov %*% t(x.players))
        for (i in (1:length(abilities))[-refcat]){
            se.ability[i] <-
                sqrt(v[i, i] + v[refcat, refcat] - 2*v[i, refcat])
        }
    }
    abilities <- cbind(abilities, se.ability)
    colnames(abilities) <- c("ability", "s.e.")
    rownames(abilities) <- playernames
## Finally, tidy up the result object
    names(result$fitted.values) <-
        names(result$y) <-
        names(result$residuals) <-
        names(result$prior.weights) <-
        names(result$linear.predictors) <-
        names(result$weights) <-
        row.names(result$model) <-
            rownames(BTyvar)
    result$call <- match.call(expand.dots = TRUE)
    result$xlevels <- xlevels
    result$contrasts <- xcont
    result$terms <- Terms
    result$y0 <- ymat[valid.contest, ]
    result$order.effect <- order.effect[valid.contest]
    result$offset0 <- offset0
    xassign <- attr(x0, "assign")
    if (dim(x0)[2] > 0){
        x0 <- x0[, -1, drop = FALSE]
        attr(x0, "assign") <- xassign[-1]
    }
    result$x0 <- x0
    result$abilities <- abilities
    class(result) <- c("BTm", class(result))
    return(result)}






BTresiduals <- function(model){
    if (!inherits(model, "BTm")) stop(
          "not a Bradley-Terry model")
    res <- model$residuals  ## the "working" residuals
    denom <- model$weights
    nwords <- if (is.null(model$call$order.effect)) 3 else 4
    players <- matrix(unlist(strsplit(row.names(model$model),
                                      split = " ")),
                      nwords, length(res))[c(1,3),]
    playernames <- rownames(model$x0)
    nplayers <- nrow(model$x0)
    res.mat <- denom.mat <-
        matrix(0, nplayers, nplayers,
               dimnames = list(playernames, playernames))
    for (m in 1:ncol(players)){
        i <- players[1, m]
        j <- players[2, m]
        res.mat[i, j] <- res.mat[i, j] + res[m]*denom[m]
        res.mat[j, i] <- -res.mat[i, j]
        denom.mat[i, j] <- denom.mat[j, i] <-
            denom.mat[i, j] + denom[m]
    }
    result <- apply(res.mat, 1, sum) / apply(denom.mat, 1, sum)
    attr(result, "weights") <- apply(denom.mat, 1, sum)
    result
}

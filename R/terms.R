terms.BTm <- function(x, ...){
    model <- x
    v <- model$terms
    a <- attributes(v)
    yname <- eval(model$call$formula, parent.frame())[[2]]
    v <- as.formula(paste(yname, paste(as.character(v), 
    	collapse = " ")))
    attributes(v) <- a
    v
}


formula.BTm <- function(x, ...){
    model <- x
    formula(terms(model))
}

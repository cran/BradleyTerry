BTabilities <- function(model){
    if (!inherits(model, "BTm")) stop("model is not of class BTm")
    model$abilities
}

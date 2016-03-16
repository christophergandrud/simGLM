#' Internal function for finding highest probability density intervals from
#' a vector of simulations
#'
#' @param x vector of simulation values
#' @param probs numeric specifying the probability interval
#' @param side character either 'upper' or 'lower' specifying which side of the
#' interval to return.
#'
#' @importFrom coda as.mcmc HPDinterval
#' @importFrom dplyr %>%
#'
#' @export

HPD <- function(x, probs, side){
    if (length(x) <= 1) {
        out <- NA
    }
    else {
        both <- as.mcmc(x) %>% HPDinterval(prob = probs)
        if (side == 'lower') {
            out <- both[1]
        }
        else if (side == 'upper'){
            out <- both[2]
        }
    }
    return(out)
}
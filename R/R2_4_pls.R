R2 <- function(object, ...) UseMethod("R2")
#' R2_4_pls
#' 
#' NORMALLY, THIS FUNCTION IS PROVIDED BY THE PLS PACKAGE.
#' HAD SOME PROBLEM WITH A GIVEN (??) VERSION OF THE PLS PACKAGE.
#' THIS CODE WAS COPIED FROM THE PLS PACKAGE VERSION 2.7-1.
#' Helper function to get R^2 values for pls fits.
#' Depends on pls package.
#'
#' @param object 
#' @param estimate 
#' @param newdata 
#' @param ncomp 
#' @param comps 
#' @param intercept 
#' @param se 
#' @param ... 
#'
#' @export
#'
#' @return
#'
R2.mvr <- function(object, estimate, newdata, ncomp = 1:object$ncomp, comps,
                   intercept = cumulative, se = FALSE, ...) {
  ## Makes the code slightly simpler:  FIXME: maybe remove
  cumulative <- missing(comps) || is.null(comps)
  
  ## Figure out which estimate(s) to calculate:
  allEstimates <- c("all", "train", "CV", "test")
  if (missing(estimate)) {
    ## Select the `best' available estimate
    if (!missing(newdata)) {
      estimate = "test"
    } else {
      if (!is.null(object$validation)) {
        estimate = "CV"
      } else {
        estimate = "train"
      }
    }
  } else {
    estimate <- allEstimates[pmatch(estimate, allEstimates)]
    if (any(is.na(estimate))) stop("`estimate' should be a subset of ",
                                   paste(allEstimates, collapse = ", "))
    if (any(estimate == "all")) {
      estimate <- allEstimates[-1] # Try all estimates (except "all")
      if(missing(newdata)) estimate <- setdiff(estimate, "test")
      if(is.null(object$validation) || !cumulative)
        estimate <- setdiff(estimate, "CV")
    }
  }
  
  ## Get the needed validation statistics:
  cl <- match.call(expand.dots = FALSE)
  cl$estimate <- estimate             # update estimate argument
  cl[[1]] <- as.name("mvrValstats")
  valstats <- eval(cl, parent.frame())
  
  ## Calculate the R^2s:
  R2 <- 1 - valstats$SSE / c(valstats$SST)
  
  return(structure(list(val = R2, type = "R2", comps = valstats$comps,
                        cumulative = valstats$cumulative, call = match.call()),
                   class = "mvrVal"))
}


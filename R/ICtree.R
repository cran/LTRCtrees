#' Fit a survival tree for interval-censored survival data
#'
#' Recursive partition for interval-censored survival data in a
#' conditional inference framework.
#'
#' \code{ICtree} returns a \link[partykit]{party} object. This function extends
#' the conditional inference survival tree algorithm  in \code{\link[partykit]{ctree}}
#' to fit interval-censored survival data. This function requires the interval package,
#' which in turn requires the Icens package, which is not available on CRAN. To install
#' the Icens package, enter the following commands
#'
#' source("https://bioconductor.org/biocLite.R")
#'
#' biocLite("Icens")
#'
#' @param Formula A formula object, with the response be a \link[survival]{Surv}
#' object, with form Surv(time1, time2, type="interval2")
#' @param data A data frame contains the variables named in Formula.
#' @param Control A list of control parameters, see \link[partykit]{ctree_control}
#'
#' @return An object of class \link[partykit]{party}.
#'
#' @references Fu, W. and Simonoff, J.S. (2017). Survival trees for Interval Censored Survival data.
#' arXiv, URL: https://arxiv.org/abs/1702.07763
#'
#'
#' @examples
#' library(interval)
#' library(LTRCtrees)
#' data(bcos)
#'
#' ## Fit ICtree survival tree
#' ## make sure to attach survival package (by library(survival) ) before using Surv function
#' Ctree <- ICtree(Surv(left,right,type="interval2")~treatment, data = bcos)
#'
#' ## Plot the fitted tree
#' plot(Ctree)
#'
#'@export
ICtree <- function(Formula, data, Control = partykit::ctree_control()){
  #library(partykit) -- in order to get partykit::extree_data function
  requireNamespace("inum")

  DATA = data
  X <- DATA[,as.character(Formula[[2]][[2]])]
  Y <- DATA[,as.character(Formula[[2]][[3]])]

  if(sum(X==Y)){ ## add small noise to observed event case
    epsilon <- min(diff(sort(unique(c(X,Y)))))/20
    ID <- which(X==Y)
    DATA[ID,as.character(Formula[[2]][[3]])] <- epsilon + DATA[ID,as.character(Formula[[2]][[3]])]
  }

  if(sum(X == Inf)){
    stop("There are Infs in the left end of intervals, make sure the left end is bounded")
  }

  if(sum(Y == Inf)){
    if(length(Y[Y!=Inf])==0){
      stop("All Infs on the right end of censoring interval, unable to compute")
    }
    Right_end_impute <- max(Y[Y!=Inf])*100 # impute Inf with 100 times the max observed Y value
    DATA[which(Y==Inf),as.character(Formula[[2]][[3]])] <- Right_end_impute
  }

  ## x2 is Surv(Left,right,type="interval2") object
  .logrank_trafo <- function(x2){
    if(!(survival::is.Surv(x2) && isTRUE(attr(x2, "type") == "interval"))){stop("Response must be a 'Survival' object with Surv(time1,time2,event) format")}

    # Fit IC survival curve
    Curve <- interval::icfit(x2~1)

    ## get estimated survival
    Left <- interval::getsurv(x2[,1], Curve)[[1]]$S
    Right<- interval::getsurv(x2[,2], Curve)[[1]]$S

    Log_Left <- ifelse(Left<=0,0,Left*log(Left))
    Log_Right<- ifelse(Right<=0,0,Right*log(Right))
    result <- (Log_Left-Log_Right)/(Left-Right)

    return(as.double(result))
  }

  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- .logrank_trafo(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  partykit::ctree(formula = Formula, data=DATA, ytrafo=h2, control = Control)
}

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
#' Ctree <- ICtree(Surv(left,right,type="interval2")~treatment, data = bcos)
#'
#' ## Plot the fitted tree
#' plot(Ctree)
#'
#'@export
ICtree <- function(Formula, data, Control = partykit::ctree_control()){
  DATA = data
  X <- DATA[,as.character(Formula[[2]][[2]])]
  Y <- DATA[,as.character(Formula[[2]][[3]])]

  if(sum(X==Y)){ ## add small noise to observed event case
    epsilon <- min(diff(sort(unique(c(X,Y)))))/20
    ID <- which(X==Y)
    DATA[ID,as.character(Formula[[2]][[3]])] <- epsilon + DATA[ID,as.character(Formula[[2]][[3]])]
  }

  Response <- as.character(Formula)[[2]]

  ## x2 is Surv(Left,right,type="interval2") object
  .logrank_trafo <- function(x2){

    # Fit IC survival curve
    Curve <- interval::icfit(x2~1)

    ## get estimated survival
    Left <- interval::getsurv(x2[,1], Curve)[[1]]$S
    Right<- interval::getsurv(x2[,2], Curve)[[1]]$S

    Log_Left <- ifelse(Left<=0,0,Left*log(Left))
    Log_Right<- ifelse(Right<=0,0,Right*log(Right))
    result <- (Log_Left-Log_Right)/(Left-Right)

    return(matrix(as.double(result), ncol = 1))
  }

  h2 <- function(data, weights) {
    s <- data[, Response]
    s <- .logrank_trafo(s[weights > 0,])
    r <- rep(0, nrow(data))
    r[weights > 0] <- s
    matrix(as.double(r), ncol = 1)
  }

  partykit::ctree(formula = Formula, data=DATA, ytrafo=h2, control = Control)
}

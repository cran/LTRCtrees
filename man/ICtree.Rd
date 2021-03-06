% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICtree.R
\name{ICtree}
\alias{ICtree}
\title{Fit a survival tree for interval-censored survival data}
\usage{
ICtree(Formula, data, Control = partykit::ctree_control())
}
\arguments{
\item{Formula}{A formula object, with the response be a \link[survival]{Surv}
object, with form Surv(time1, time2, type="interval2")}

\item{data}{A data frame contains the variables named in Formula.}

\item{Control}{A list of control parameters, see \link[partykit]{ctree_control}}
}
\value{
An object of class \link[partykit]{party}.
}
\description{
Recursive partition for interval-censored survival data in a
conditional inference framework.
}
\details{
\code{ICtree} returns a \link[partykit]{party} object. This function extends
the conditional inference survival tree algorithm  in \code{\link[partykit]{ctree}}
to fit interval-censored survival data. This function itself not longer requires the interval package, but running the example below
requires the interval package (for bcos data), which in turn requires the Icens package, which is not available on CRAN. To install
the Icens package, enter the following commands

source("https://bioconductor.org/biocLite.R")

biocLite("Icens")
}
\examples{
library(Icens)
library(interval)
library(LTRCtrees)
data(bcos)

## Fit ICtree survival tree
## make sure to attach survival package (by library(survival) ) before using Surv function
Ctree <- ICtree(Surv(left,right,type="interval2")~treatment, data = bcos)

## Plot the fitted tree
plot(Ctree)

}
\references{
Fu, W. and Simonoff, J.S. (2017). Survival trees for Interval Censored Survival data.
Statistics in medicine 36 (30), 4831-4842
}

#' Copy the partykit::extree_data function from partykit to avoid dependency issue
#'
#' \code{extree_data} imports partykit::extree_data function
#'
#' @param formula Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param data Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param subset Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param na.action Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param weights Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param offset Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param cluster Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param strata Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param scores Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param yx Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param ytype Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param nmax Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @param ... Same as the one in \link[partykit]{extree_data}, check \link[partykit]{extree_data} for usage
#' @return check \link[partykit]{extree_data} for the return value
extree_data <- function(formula, data, subset, na.action = stats::na.pass, weights, offset, cluster,
                        strata, scores = NULL, yx = c("none", "matrix"), ytype = c("vector", "data.frame", "matrix"),
                        nmax = c("yx" = Inf, "z" = Inf), ...)
{
  ## call
  cl <- match.call()
  yx <- match.arg(yx, choices = c("none", "matrix"))
  ytype <- match.arg(ytype, choices = c("vector", "data.frame", "matrix"))

  ## 'formula' may either be a (multi-part) formula or a list
  noformula <- !inherits(formula, "formula")
  if(noformula) {

    ## formula needs to be a 'list' (if it is not a 'formula')
    if(!inherits(formula, "list")) stop("unsupported specification of 'formula'")

    ## specified formula elements and overall call elements
    fonam <- names(formula)
    clnam <- names(cl)[-1L]
    vanam <- c("y", "x", "z", "weights", "offset", "cluster", "strata")

    ## y and z (and optionally x) need to be in formula
    if(!all(c("y", "z") %in% fonam)) stop("'formula' needs to specify at least a response 'y' and partitioning variables 'z'")
    if(!("x" %in% fonam)) formula$x <- NULL

    ## furthermore weights/offset/cluster/strata may be in formula or call
    vars <- formula[vanam]
    names(vars) <- vanam
    if("weights" %in% clnam) {
      clvar <- try(weights, silent = TRUE)
      vars[["weights"]] <- c(vars[["weights"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$weights))
    }
    if("offset" %in% clnam) {
      clvar <- try(offset, silent = TRUE)
      vars[["offset"]] <- c(vars[["offset"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$offset))
    }
    if("cluster" %in% clnam) {
      clvar <- try(cluster, silent = TRUE)
      vars[["cluster"]] <- c(vars[["cluster"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$cluster))
    }
    if("strata" %in% clnam) {
      clvar <- try(strata, silent = TRUE)
      vars[["strata"]] <- c(vars[["strata"]], if(!inherits(clvar, "try-error")) clvar else deparse(cl$strata))
    }

    ## sanity checking
    for(v in vanam) {
      if(!is.null(vars[[v]]) && !(is.numeric(vars[[v]]) | is.character(vars[[v]]) | is.logical(vars[[v]]))) {
        warning(sprintf("unknown specification of '%s', must be character, numeric, or logical", v))
        vars[v] <- list(NULL)
      }
    }
    if(!missing(subset)) warning("'subset' argument ignored in list specification of 'formula'")
    if(!missing(na.action)) warning("'na.action' argument ignored in list specification of 'formula'")

    ## no terms (by default)
    mt <- NULL

  } else {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    mf$na.action <- na.action ### evaluate na.action
    if(missing(data)) data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "strata"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$dot <- "sequential"

    ## formula processing
    oformula <- stats::as.formula(formula)
    formula <- Formula::as.Formula(formula)
    mf$formula <- formula
    npart <- length(formula)
    if(any(npart < 1L)) stop("'formula' must specify at least one left-hand and one right-hand side")
    npart <- npart[2L]

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms in various combinations
    mt <- list(
      "all" = stats::terms(formula, data = data,                        dot = "sequential"),
      "y"   = stats::terms(formula, data = data, rhs = 0L,              dot = "sequential"),
      "z"   = stats::terms(formula, data = data, lhs = 0L, rhs = npart, dot = "sequential")
    )
    if(npart > 1L) {
      mt$yx <-stats::terms(formula, data = data, rhs = 1L,              dot = "sequential")
      for(i in 1L:(npart-1L)) {
        mt[[paste("x", if(i == 1L) "" else i, sep = "")]] <- stats::terms(
          formula, data = data, lhs = 0L, rhs = i,     dot = "sequential")
      }
    }

    ## extract variable lists
    vars <- list(
      y = attr(mt$y, "term.labels"),
      x = unique(unlist(lapply(grep("^x", names(mt)), function(i) attr(mt[[i]], "term.labels")))),
      z = attr(mt$z, "term.labels"),
      weights = if("(weights)" %in% names(mf)) "(weights)" else NULL,
      offset  = if("(offset)"  %in% names(mf)) "(offset)"  else NULL,
      cluster = if("(cluster)" %in% names(mf)) "(cluster)" else NULL,
      strata = if("(strata)" %in% names(mf)) "(strata)" else NULL
    )
    ymult <- length(vars$y) >= 1L
    if(!ymult) vars$y <- names(mf)[1L]
    ## FIXME: store information which variable(s) went into (weights), (offset), (cluster)
    ## (strata)
    ## idea: check (x and) z vs. deparse(cl$weights), deparse(cl$offset), deparse(cl$cluster)

    ## check wether offset was inside the formula
    if(!is.null(off <- attr(mt$x, "offset"))) {
      if(is.null(vars$offset)) mf[["(offset)"]] <- rep.int(0, nrow(mf))
      for(i in off) mf[["(offset)"]] <- mf[["(offset)"]] + mf[[i]]
      vars$offset <- "(offset)"
    }
  }

  ## canonicalize y/x/z term labels
  vanam <- if(noformula) names(data) else names(mf)
  ## z to numeric
  if(is.null(vars$z)) stop("at least one 'z' variable must be specified")
  if(is.integer(vars$z)) vars$z <- vanam[vars$z]
  if(is.character(vars$z)) vars$z <- vanam %in% vars$z
  if(is.logical(vars$z)) vars$z <- as.numeric(vars$z)
  if(is.null(names(vars$z))) names(vars$z) <- vanam
  vars$z <- vars$z[vanam]
  if(any(is.na(vars$z))) vars$z[is.na(vars$z)] <- 0
  vars$z <- as.numeric(vars$z)
  ## all others to integer
  for(v in c("y", "x", "weights", "offset", "cluster", "strata")) {
    if(!is.null(vars[[v]])) {
      if(is.character(vars[[v]])) vars[[v]] <- match(vars[[v]], vanam)
      if(is.logical(vars[[v]])) vars[[v]] <- which(vars[[v]])
      if(any(is.na(vars[[v]]))) {
        vars[[v]] <- vars[[v]][!is.na(vars[[v]])]
        warning(sprintf("only found the '%s' variables: %s", v, paste(vanam[vars[[v]]], collapse = ", ")))
      }
    }
    vars[[v]] <- unique(as.integer(vars[[v]]))
  }
  if(is.null(vars$y)) stop("at least one 'y' variable must be specified")


  ## FIXME: subsequently fitting, testing, splitting
  ## - fit: either pre-processed _and_ subsetted data --or-- full data object plus subset vector
  ## - test: additionally needs fit output --and-- fit function
  ## - split: additionally needs test output
  ## - tbd: control of all details

  ret <- list(
    data = if(noformula) data else mf,
    variables = vars,
    terms = mt
  )

  mf <- ret$data
  yxvars <- c(vars$y, vars$x, vars$offset, vars$cluster)
  zerozvars <- which(vars$z == 0)

  ret$scores <- vector(mode = "list", length = length(ret$variables$z))
  names(ret$scores) <- names(mf)
  if (!is.null(scores))
    ret$scores[names(scores)] <- scores

  if (length(nmax) == 1) nmax <- c("yx" = nmax, "z" = nmax)
  ### <FIXME> make meanlevels an argument and make sure intersplit is TRUE </FIXME>
  ret$zindex <- inum::inum(mf, ignore = names(mf)[zerozvars], total = FALSE,
                           nmax = nmax["z"], meanlevels = FALSE)
  if (is.finite(nmax["yx"])) {
    ret$yxindex <- inum::inum(mf[, yxvars, drop = FALSE], total = TRUE,
                              as.interval = names(mf)[vars$y], complete.cases.only = TRUE,
                              nmax = nmax["yx"], meanlevels = FALSE)
    yxmf <- attr(ret$yxindex, "levels")
    yxmf[["(weights)"]] <- NULL
    attr(ret$yxindex, "levels") <- yxmf
  } else {
    ret$yxindex <- NULL
    yxmf <- mf
  }

  ret$missings <- lapply(ret$data, function(x) which(is.na(x)))
  ret$yxmissings <- sort(unique(do.call("c", ret$missings[yxvars])))

  ## FIXME: separate object with options for: discretization, condensation, some NA handling
  ## below is just "proof-of-concept" implementation using plain model.matrix() which could
  ## be included as one option...
  if (yx == "matrix") {

    ## fake formula/terms if necessary
    formula <- Formula::as.Formula(sprintf("%s ~ %s | %s",
                                           paste(vanam[vars$y], collapse = " + "),
                                           if(length(vars$x) > 0L) paste(vanam[vars$x], collapse = " + ") else "0",
                                           paste(vanam[vars$z > 0], collapse = " + ")
    ))
    mt <- list(
      "all" = stats::terms(formula),
      "y"   = stats::terms(formula, data = data, rhs = 0L),
      "z"   = stats::terms(formula, data = data, lhs = 0L, rhs = 2L),
      "yx"  = stats::terms(formula, data = data, rhs = 1L),
      "x"   = stats::terms(formula, data = data, lhs = 0L, rhs = 1L)
    )
    ymult <- length(vars$y) > 1L
    npart <- 2L

    if (ytype == "vector" && !ymult) {
      yx <- list("y" = yxmf[, vanam[vars$y], drop = TRUE])
    } else if (ytype == "data.frame") {
      yx <- list("y" = yxmf[vanam[vars$y]])
    } else { ### ytype = "matrix"
      Ytmp <- stats::model.matrix(~ 0 + ., Formula::model.part(formula, yxmf, lhs = TRUE))
      ### <FIXME> are there cases where Ytmp already has missings? </FIXME>
      if (length(ret$yxmissings) == 0) {
        Ymat <- Ytmp
      } else {
        Ymat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Ytmp))
        Ymat[-ret$yxmissings,] <- Ytmp
      }
      yx <- list("y" = Ymat)
    }
    for(i in (1L:npart)[-npart]) {
      ni <- paste("x", if(i == 1L) "" else i, sep = "")
      ti <- if(!ymult & npart == 2L) mt$yx else mt[[ni]]
      Xtmp <- stats::model.matrix(ti, yxmf)
      if (length(ret$yxmissings) == 0) {
        Xmat <- Xtmp
      } else {
        Xmat <- matrix(0, nrow = NROW(yxmf), ncol = NCOL(Xtmp))
        Xmat[-ret$yxmissings,] <- Xtmp
      }
      yx[[ni]] <- Xmat
      if(ncol(yx[[ni]]) < 1L) {
        yx[[ni]] <- NULL
      } else {
        attr(yx[[ni]], "formula") <- formula(formula, rhs = i)
        attr(yx[[ni]], "terms") <- ti
        attr(yx[[ni]], "offset") <- yxmf[["(offset)"]]
      }
    }
    ret$yx <- yx
  }

  class(ret) <- "extree_data"
  ret
}


#' Logrank transformation function for LTRC data
#'
#' \code{.logrank_trafo} transforms Surv(time1, time2, event) objects into
#' logrank scores, which will be used later in the tree algorithm. It is
#' not designed to be used by users, not for internal used of \code{LTRCIT}
#' function.
#'
#' @param x2 A vector \link[survival]{Surv} (Surv(time1, time2, event)) objects
#' @return Logrank scores of LTRC objects
.logrank_trafo2 <- function(x2){

  unique.times <- unique(x2[,2][which(x2[,3]==1)])

  D <- rep(NA,length(unique.times))
  R <- rep(NA,length(unique.times))

  for(j in 1:length(unique.times)){
    D[j] = sum(unique.times[j]==x2[,2])
  }

  for(k in 1:length(unique.times) ){
    value <- unique.times[k]
    R[k]<- sum(apply(x2[,1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
  }

  Ratio <- D/R

  Ratio <- Ratio[order(unique.times)]
  Nelson.Aalen <- cumsum(Ratio)
  Event.time <- unique.times[order(unique.times)]

  Left <- sapply(x2[,1],function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
  Right <- sapply(x2[,2],function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

  result<- x2[,3]-(Right-Left)

  return(as.double(result))
}

#' Fit a conditional inference survival tree for LTRC data
#'
#' \code{LTRCIT} returns a \link[partykit]{party} object. This function extends
#' the conditional inference survival tree algorithm in \code{\link[partykit]{ctree}}
#' to fit left-truncated and right censored (LTRC) data.
#'
#' @param Formula A formula object, with the response be a \link[survival]{Surv}
#' object, with form Surv(time1, time2, event)
#' @param data A data frame contains the variables named in formula.
#' @param Control A list of control parameters, see \link[partykit]{ctree_control}
#'
#' @return An object of class \link[partykit]{party}.
#'
#' @references Fu, W. and Simonoff, J.S.(2017). Survival trees for left-truncated and right-censored data,
#' with application to time-varying covariate data. Biostatistics, URL: https://doi.org/10.1093/biostatistics/kxw047
#'
#'
#' @examples
#' ## The Assay of serum free light chain data in survival package
#' ## Adjust data & clean data
#' library(survival)
#' library(LTRCtrees)
#' Data <- flchain
#' Data <- Data[!is.na(Data$creatinine),]
#' Data$End <- Data$age + Data$futime/365
#' DATA <- Data[Data$End > Data$age,]
#' names(DATA)[6] <- "FLC"
#'
#' ## Setup training set and test set
#' Train = DATA[1:500,]
#' Test = DATA[1000:1020,]
#'
#' ## Fit LTRCIT survival tree
#' ## make sure to attach survival package (by library(survival) ) before using Surv function
#' LTRCIT.obj <-  LTRCIT(Surv(age, End, death) ~ sex + FLC + creatinine, Train)
#' plot(LTRCIT.obj)
#'
#' ## Putting Surv(End, death) in formula would result an error message
#' ## since LTRCIT is expecting Surv(time1, time2, event)
#'
#' ## Note that LTRCIT.obj is an object of class party
#' ## predict median survival time on test data
#' LTRCIT.pred <- predict(LTRCIT.obj, newdata = Test, type = "response")
#'
#' ## predict Kaplan Meier survival curve on test data,
#' ## return a list of survfit objects -- the predicted KM curves
#' LTRCIT.pred <- predict(LTRCIT.obj, newdata = Test, type = "prob")
#'
#'####################################################################
#'####### Survival tree with time-varying covariates ##################
#'####################################################################
#'## The pbcseq dataset of survival package
#'library(survival)
#'## Create the start-stop-event triplet needed for coxph and LTRC trees
#'first <- with(pbcseq, c(TRUE, diff(id) !=0)) #first id for each subject
#'last <- c(first[-1], TRUE) #last id
#'time1 <- with(pbcseq, ifelse(first, 0, day))
#'time2 <- with(pbcseq, ifelse(last, futime, c(day[-1], 0)))
#'event <- with(pbcseq, ifelse(last, status, 0))
#'event <- 1*(event==2)
#'
#'pbcseq$time1 <- time1
#'pbcseq$time2 <- time2
#'pbcseq$event <-  event
#'
#'pbcseq = pbcseq[1:1000,] ## fit on subset of the data to save fitting time
#' ## Fit the Cox model and LTRCIT tree with time-varying covariates
#'fit.cox <- coxph(Surv(time1, time2, event) ~ age + sex + log(bili), pbcseq)
#'LTRCIT.fit <- LTRCIT(Surv(time1, time2, event) ~ age + sex + log(bili), pbcseq)
#'plot(LTRCIT.fit)
#'
#'## transform the wide format data into long format data using tmerge function
#'## from survival function
#'## Stanford Heart Transplant data
#'jasa$subject <- 1:nrow(jasa)
#'
#'tdata <- with(jasa, data.frame(subject = subject,
#'                               futime= pmax(.5, fu.date - accept.dt),
#'                               txtime= ifelse(tx.date== fu.date,
#'                                              (tx.date -accept.dt) -.5,
#'                                              (tx.date - accept.dt)),
#'                               fustat = fustat))
#'
#'sdata <- tmerge(jasa, tdata, id=subject,death = event(futime, fustat),
#'                    trt = tdc(txtime), options= list(idname="subject"))
#'
#'sdata$age <- sdata$age - 48
#'
#'sdata$year <- as.numeric(sdata$accept.dt - as.Date("1967-10-01"))/365.25
#'
#'Cox.fit <- coxph(Surv(tstart, tstop, death) ~ age+ surgery, data= sdata)
#'LTRCIT.fit <- LTRCIT(Surv(tstart, tstop, death) ~ age + transplant, data = sdata)
#'plot(LTRCIT.fit)
#'
#' @export
LTRCIT <- function(Formula, data, Control = partykit::ctree_control()){
  #library(partykit) -- in order to get partykit::extree_data function
  requireNamespace("inum")
  if(length(as.list(Formula[[2]]))!=4){
    stop(" Response must be a 'survival' object with 'Surv(time1, time2, event)' ")
  }
  Response <- as.character(Formula)[[2]]
  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- .logrank_trafo2(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }
  partykit::ctree(formula = Formula, data=data, ytrafo=h2, control = Control)
}

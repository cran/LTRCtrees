#' Logrank transformation function for LTRC data
#'
#' \code{.logrank_trafo} transforms Surv(time1, time2, event) objects into
#' logrank scores, which will be used later in the tree algorithm. It is
#' not designed to be used by users, not for internal used of \code{LTRCIT}
#' function.
#'
#' @param x2 A vector \link[survival]{Surv} (Surv(time1, time2, event)) objects
#' @return Logrank scores of LTRC objects
.logrank_trafo <- function(x2){

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

  return(matrix(as.double(result), ncol = 1))
}

#' Fit a conditional inference survival tree for LTRC data
#'
#' \code{LTRCIT} returns an \link[partykit]{party} object. This function extends
#' the conditional inference survival tree algorithm in \code{\link[partykit]{ctree}}
#' to fit left-truncated and right censored (LTRC) data.
#'
#' @param Formula A formula object, with the response be a \link[survival]{Surv}
#' object, with form Surv(time1, time2, event)
#' @param Data A data frame contains the variables named in formula.
#' @param Control A list of control parameters, see \link[partykit]{ctree_control}
#'
#' @return An object of class \link[partykit]{party}.
#'
#' @references Fu, W. and Simonoff, J.S.(2016). Survival trees for left-truncated and right-censored data,
#' with application to time-varying covariate data. arXiv:1606.03033 [stat.ME]
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
#'LTRCART.fit <- LTRCART(Surv(tstart, tstop, death) ~ age + transplant, data = sdata)
#'plot(LTRCIT.fit)
#'
#' @export
LTRCIT <- function(Formula, Data, Control = partykit::ctree_control()){
  if(length(as.list(Formula[[2]]))!=4){
    stop(" Response must be a 'survival' object with 'Surv(time1, time2, event)' ")
  }
  Response <- as.character(Formula)[[2]]
  h2 <- function(data, weights) {
    s <- data[, Response]
    s <- .logrank_trafo(s[weights > 0,])
    r <- rep(0, nrow(data))
    r[weights > 0] <- s
    matrix(as.double(r), ncol = 1)
  }
  result <- partykit::ctree(formula = Formula, data=Data, ytrafo=h2, control = Control)
  return(result)
}

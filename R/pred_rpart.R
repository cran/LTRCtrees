#' Prediction function for \link[rpart]{rpart.object}
#'
#'
#' The output of \code{LTRCART} is an \link[rpart]{rpart} object, and as a result the
#' usual \link{predict} function on such an object returns the predicted
#' relative risk on the test set. \code{Pred.rpart} returns the predicted
#' Kaplan-Meier curves and median survival times on the test set,
#' which in some circumstances might be desirable in practice.
#' Note that this function can be applied to any \link[rpart]{rpart} survival tree
#' object, not just one produced by \code{LTRCART}
#'
#' @param formula A formula used to fit the survival tree. The response
#' is a \link[survival]{Surv} object. If it has the form Surv(time1, time2, event),
#' then \code{LTRCART} is called internally; if response has the form Surv(time, event),
#' then the \link[rpart]{rpart} is called internally.
#'
#' @param train Training set
#' @param test Test set
#'
#' @return A list of predicted KM curves and median survival times.
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
#' ## Predict median survival time and Kaplan Meier survival curve
#' ## on test data using Pred.rpart
#' LTRCART.pred <- Pred.rpart(Surv(age, End, death) ~ sex + FLC + creatinine, Train, Test)
#' LTRCART.pred$KMcurves  ## list of predicted KM curves
#' LTRCART.pred$Medians  ## vector of predicted median survival time
#'
#' @export
Pred.rpart <- function(formula, train, test){

  if( length(formula[[2]])==3){
    rtree <- rpart::rpart(formula, train)
    Formula = formula
    Formula[[3]] = 1

  }else if(length(formula[[2]])==4){
    rtree <- LTRCART(formula, train)
    Formula = formula
    Formula[[3]] = 1

  }else{
    stop("Not Surv object")
  }

  Train = train
  Train$ID <- stats::predict(rtree, type="vector")
  Keys <- unique(Train$ID)
  Keys.MM <- matrix(c(Keys,1:length(Keys)),ncol=2)

  List.KM <- list()
  List.Med <- list()

  for(p in Keys){
    subset <- Train[Train$ID == p,]
    KM <-  survival::survfit(Formula, data = subset)
    Median <- utils::read.table(textConnection(utils::capture.output(KM)),skip=2,header=TRUE)$median
    List.KM[[ Keys.MM[Keys.MM[,1]==p,2] ]] = KM
    List.Med[[ Keys.MM[Keys.MM[,1]==p,2] ]] = Median
  }

  Test.keys <- stats::predict(rtree, newdata = test, type="vector")
  Test.ID <- match(Test.keys,Keys.MM[,1])
  Test.KM <- List.KM[Test.ID]
  Test.Med <- unlist(List.Med[Test.ID])

  result <- list(KMcurves = Test.KM, Medians = Test.Med)
  return(result)
}

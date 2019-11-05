library(data.table)
folds.csv.vec <- Sys.glob("data/*/cv/*/folds.csv")
for(folds.csv in folds.csv.vec){
  folds.dt <- fread(folds.csv)
  cv.type.path <- dirname(folds.csv)
  testFolds.path.vec <- Sys.glob(file.path(cv.type.path, "testFolds", "*"))
  for(testFolds.path in testFolds.path.vec){
    test.fold <- basename(testFolds.path)
    test.ids <- folds.dt[fold==test.fold, sequenceID]
    pred.csv <- file.path(
      testFolds.path,
      "randomTrainOrderings/1/models/L1reg_linear_all/predictions.csv")
    if(file.exists(pred.csv)){
      pred.dt <- fread(pred.csv, header=TRUE)
      pred.all.train <- pred.dt[, c(1, ncol(pred.dt)), with=FALSE]
      names(pred.all.train)[2] <- "pred.log.lambda"
      stopifnot(identical(sort(test.ids), sort(pred.all.train$sequenceID)))
      cat("OK: ", pred.csv, "\n", sep="")
    }
  }
}

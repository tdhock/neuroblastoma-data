library(data.table)
folds.csv.vec <- Sys.glob("data/*/cv/*/folds.csv")
pred.err.list <- list()
for(folds.csv.i in seq_along(folds.csv.vec)){
  folds.csv <- folds.csv.vec[[folds.csv.i]]
  folds.dt <- fread(folds.csv)
  cv.type.dir <- dirname(folds.csv)
  cv.dir <- dirname(cv.type.dir)
  data.dir <- dirname(cv.dir)
  data.name <- basename(data.dir)
  data.list <- list()
  for(data.type in c("inputs", "outputs")){
    csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
    type.dt <- fread(cmd=paste("xzcat", csv.xz))
    keep.mat <- if(data.type=="inputs"){
      m <- as.matrix(type.dt[,-1, with=FALSE])
      cbind(TRUE, is.finite(m))
    }else{
      !is.na(type.dt)
    }
    keep <- apply(keep.mat, 2, all)
    data.list[[data.type]] <- type.dt[, keep, with=FALSE]
  }
  u.fold.vec <- sort(unique(folds.dt$fold))
  for(test.fold in u.fold.vec){
    folds.dt[, set := ifelse(fold==test.fold, "test", "train")]
    set.list <- list()
    for(set in unique(folds.dt$set)){
      set.ids <- folds.dt[set, sequenceID, on=list(set)]
      for(data.type in names(data.list)){
        set.dt <- data.list[[data.type]][set.ids, on=list(sequenceID)]
        set.mat <- as.matrix(set.dt[, -1, with=FALSE])
        set.list[[set]][[data.type]] <- set.mat
      }
    }
    cat(sprintf(
      "%d / %d data=%s, fold=%d\n",
      folds.csv.i, length(folds.csv.vec), data.name, test.fold))
    fun.list <- list(
      iregnet.scale1=function(inputs, outputs){
        iregnet::cv.iregnet(
          inputs, outputs,
          family="gaussian", estimate_scale=FALSE, scale_init=1)
      },
      iregnet.estscale=function(inputs, outputs){
        iregnet::cv.iregnet(
          inputs, outputs, family="gaussian")
      },
      penaltyLearning.scale1=function(inputs, outputs){
        penaltyLearning::IntervalRegressionCV(
          inputs, outputs)
      })
    fit.list <- list()
    seconds.vec <- c(constant=0)
    for(algorithm in names(fun.list)){
      fun <- fun.list[[algorithm]]
      seconds.vec[[algorithm]] <- system.time({
        fit.list[[algorithm]] <- with(set.list$train, fun(inputs, outputs))
      })[["elapsed"]]
    }
    plot(fit.list$iregnet.scale1)
    plot(fit.list$iregnet.estscale)
    plot(fit.list$penaltyLearning.scale1)
    pred.list <- list(
      constant=rep(0, nrow(set.list$test$inputs)))
    for(algorithm in names(fit.list)){
      fit <- fit.list[[algorithm]]
      pred.list[[algorithm]] <- predict(fit, set.list$test$inputs)
    }
    for(algorithm in names(pred.list)){
      pred.vec <- as.numeric(pred.list[[algorithm]])
      roc.list <- penaltyLearning::targetIntervalROC(
        set.list$test$outputs, pred.vec)
      pred.err.list[[paste(data.name, test.fold, algorithm)]] <- with(
        roc.list, data.table(
          data.name, test.fold, algorithm,
          seconds=seconds.vec[[algorithm]],
          thresholds[threshold=="predicted"], auc))
    }
  }
}
(pred.err <- do.call(rbind, pred.err.list))

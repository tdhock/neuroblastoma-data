source("packages.R")

if(FALSE){#copy gpppredicitons.csv to mdoels/GP/predictions.csv
  gppred.csv.vec <- Sys.glob(
    "data/*/cv/*/testFolds/*/*/*/gppredictions.csv")
  gp.dir.csv.vec <- file.path(dirname(gppred.csv.vec), "models", "GP", "predictions.csv")
  file.exists(gp.dir.csv.vec)
  for(d in unique(dirname(gp.dir.csv.vec)))dir.create(d, recursive=TRUE)
  file.copy(gppred.csv.vec, gp.dir.csv.vec)
}

pred.csv.vec <- Sys.glob(
  "data/*/cv/*/testFolds/*/*/*/models/*/predictions.csv")
model.dir.vec <- dirname(pred.csv.vec)
acc.csv.vec <- file.path(model.dir.vec, "accuracy.csv")
seed.dir.vec <- dirname(dirname(model.dir.vec))
n.models <- table(seed.dir.vec)
gp.names <- names(n.models)[n.models==6]
##pred.csv.vec[seed.dir.vec %in% gp.names]
(pred.csv.todo <- pred.csv.vec[!file.exists(acc.csv.vec)])

OnePred <- function(pred.csv.i){
  pred.csv <- pred.csv.todo[[pred.csv.i]]
  model.dir <- dirname(pred.csv)
  acc.csv <- file.path(model.dir, "accuracy.csv")
  cat(sprintf(
    "%4d / %4d %s\n", pred.csv.i, length(pred.csv.todo), pred.csv))
  models.dir <- dirname(model.dir)
  seed.dir <- dirname(models.dir)
  seed <- as.integer(basename(seed.dir))
  split.dir <- dirname(dirname(seed.dir))
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  data.dir <- dirname(dirname(cv.type.dir))
  folds.csv <- file.path(cv.type.dir, "folds.csv")
  folds.dt <- fread(folds.csv)
  test.seqs <- folds.dt[test.fold==fold, sequenceID]
  data.list <- list()
  for(data.type in c("evaluation", "outputs")){
    csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
    data.list[[data.type]] <- fread(cmd=paste("xzcat", csv.xz))
  }
  pred.wide.all <- fread(pred.csv, header=TRUE)
  pred.wide <- pred.wide.all[test.seqs, on=.(sequenceID), nomatch=0L]
  stopifnot(length(test.seqs) == nrow(pred.wide))
  pred.tall <- melt(
    pred.wide,
    id.vars="sequenceID",
    variable.name="train.size.fac",
    value.name="pred.log.lambda")
  pred.tall[, train.size := as.integer(sub("n", "", train.size.fac))]
  acc.dt <- pred.tall[, {
    if(any(is.na(pred.log.lambda))){
      data.table(percent.correct.labels=NA_real_, auc=NA_real_)
    }else{
      roc <- penaltyLearning::ROChange(data.list$evaluation, .SD, "sequenceID")
      ##browser(expr=train.size==4)
      with(roc, thresholds[threshold=="predicted", data.table(
        percent.correct.labels=100-error.percent, auc)])
    }
  }, by=list(train.size)]
  ## Also compute percent correct intervals.
  pred.mat <- as.matrix(pred.wide[, -1, with=FALSE])
  test.intervals <- data.list$outputs[pred.wide$sequenceID, on=list(sequenceID)]
  acc.dt[, N.test := nrow(pred.wide)]
  acc.dt[, percent.correct.intervals := test.intervals[, 100*colMeans(
    min.log.lambda < pred.mat & pred.mat < max.log.lambda)] ]
  fwrite(acc.dt, acc.csv)
}

future::plan("multiprocess")
future.apply::future_lapply(seq_along(pred.csv.todo), OnePred)

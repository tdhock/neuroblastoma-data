source("packages.R")

data.list <- list(
  evaluation=list(),
  outputs=list())

pred.csv.vec <- Sys.glob(
  "data/*/cv/*/testFolds/*/*/*/models/*/predictions.csv")
seed.dir.vec <- dirname(dirname(dirname(pred.csv.vec)))
n.models <- table(seed.dir.vec)
gp.names <- names(n.models)[n.models==6]
pred.csv.vec[seed.dir.vec %in% gp.names]

for(pred.csv.i in seq_along(pred.csv.vec)){
  pred.csv <- pred.csv.vec[[pred.csv.i]]
  model.dir <- dirname(pred.csv)
  acc.csv <- file.path(model.dir, "accuracy.csv")
  if(!file.exists(acc.csv)){
    cat(sprintf(
      "%4d / %4d %s\n", pred.csv.i, length(pred.csv.vec), pred.csv))
    models.dir <- dirname(model.dir)
    seed.dir <- dirname(models.dir)
    seed <- as.integer(basename(seed.dir))
    split.dir <- dirname(dirname(seed.dir))
    test.fold <- as.integer(basename(split.dir))
    cv.type.dir <- dirname(dirname(split.dir))
    data.dir <- dirname(dirname(cv.type.dir))
    L <- list()
    for(data.type in names(data.list)){
      csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
      if(!data.dir %in% data.list[[data.type]]){
        data.list[[data.type]][[data.dir]] <- fread(cmd=paste("xzcat", csv.xz))
      }
      L[[data.type]] <- data.list[[data.type]][[data.dir]]
    }
    pred.wide <- fread(pred.csv, header=TRUE)
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
        roc <- penaltyLearning::ROChange(L$evaluation, .SD, "sequenceID")
        ##browser(expr=train.size==4)
        with(roc, thresholds[threshold=="predicted", data.table(
          percent.correct.labels=100-error.percent, auc)])
      }
    }, by=list(train.size)]
    ## Also compute percent correct intervals.
    pred.mat <- as.matrix(pred.wide[, -1, with=FALSE])
    test.intervals <- L$outputs[pred.wide$sequenceID, on=list(sequenceID)]
    acc.dt[, percent.correct.intervals := test.intervals[, 100*colMeans(
      min.log.lambda < pred.mat & pred.mat < max.log.lambda)] ]
    fwrite(acc.dt, acc.csv)
  }
}

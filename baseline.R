source("packages.R")

future::plan("multiprocess")

data.dir <- file.path("data", "systematic")
data.list <- list()
for(data.type in c("inputs", "outputs")){
  csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
  dt <- fread(cmd=paste("xzcat", csv.xz))
  m <- as.matrix(dt[, -1, with=FALSE])
  rownames(m) <- dt$sequenceID
  data.list[[data.type]] <- m
}

order.csv.vec <- Sys.glob(file.path(
  data.dir, "cv", "*", "testFolds",
  "*", "randomTrainOrderings", "*", "order.csv"))

## 66 to 70 are testFold=2, small max train set size.
order.i <- 67

OneOrder <- function(order.i){
  order.csv <- order.csv.vec[[order.i]]
  baseline.csv <- file.path(dirname(order.csv), "baseline.csv")
  if(file.exists(baseline.csv)){
    cat(baseline.csv, " already computed.\n", sep="")
    return(fread(baseline.csv))
  }
  order.dt <- fread(order.csv)
  seed.dir <- dirname(order.csv)
  seed <- as.integer(basename(seed.dir))
  split.dir <- dirname(dirname(seed.dir))
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  all.id.vec <- rownames(data.list$inputs)
  id.list <- list(
    train=order.dt$sequenceID,
    test=all.id.vec[!all.id.vec %in% order.dt$sequenceID])
  set.list <- list()
  for(set.name in names(id.list)){
    set.id.vec <- id.list[[set.name]]
    set.list[[set.name]] <- lapply(data.list, function(m){
      m[set.id.vec,]
    })
  }
  baseline.df <- mmit::compute_optimal_costs(
    set.list$train$outputs, 0, "square")
  s <- c(
    seq(2, 20, by=2),
    seq(25, 100, by=5),
    seq(200, 1000, by=100),
    nrow(order.dt))
  train.size.vec <- s[s <= nrow(order.dt)]
  result.list <- list()
  for(size.i in seq_along(train.size.vec)){
    train.size <- train.size.vec[[size.i]]
    train.i <- 1:train.size
    cat(sprintf(
      "%4d / %4d trainSize=%d\n",
      size.i, length(train.size.vec), train.size))
    X.train <- set.list$train$inputs[train.i,]
    y.train <- set.list$train$outputs[train.i,]
    no.limits <- apply(!is.finite(y.train), 2, all)
    pred.list <- list(
      baseline=rep(
        baseline.df$pred[train.size], nrow(set.list$test$inputs)),
      L1reg_linear=if(any(no.limits)){
        rep(NA, nrow(set.list$test$inputs))
      }else{
        set.seed(1)
        fit <- penaltyLearning::IntervalRegressionCV(
          X.train, y.train, min.observations=train.size)
        as.numeric(fit$predict(set.list$test$inputs))
      })
    for(model in names(pred.list)){
      pred.vec <- pred.list[[model]]
      stopifnot(length(pred.vec)==nrow(set.list$test$outputs))
      res.vec <- if(is.na(pred.vec[1])){
        NA
      }else{
        penaltyLearning::targetIntervalResidual(
          set.list$test$outputs, pred.vec)
      }
      result.list[[paste(size.i, model)]] <- data.table(
        train.size, model,
        accuracy.percent=mean(res.vec==0)*100)
    }#for(model
    print(result.list)
  }#for(size.i
  (result <- do.call(rbind, result.list))
  fwrite(result, baseline.csv)
  result
}

res.list <- future.apply::future_lapply(seq_along(order.csv.vec), OneOrder)

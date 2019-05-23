source("packages.R")

order.csv.vec <- Sys.glob(file.path(
  "data", "*", "cv", "*", "testFolds",
  "*", "*", "*", "order.csv"))
order.i.vec <- seq_along(order.csv.vec)

## profileSize testFold=1:3 is test set size from 36 to 96, indices 61 to 75.
future::plan("multiprocess")

future.apply::future_lapply(order.i.vec, function(order.i){
  order.csv <- order.csv.vec[[order.i]]
  seed.dir <- dirname(order.csv)
  seed <- as.integer(basename(seed.dir))
  split.dir <- dirname(dirname(seed.dir))
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  data.dir <- dirname(dirname(cv.type.dir))
  data.list <- list()
  for(data.type in c("inputs", "outputs")){
    csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
    dt <- fread(cmd=paste("xzcat", csv.xz))
    m <- as.matrix(dt[, -1, with=FALSE])
    rownames(m) <- dt$sequenceID
    data.list[[data.type]] <- m
  }
  rep.val.vec <- c(
    log.log.bases="log2.n",
    log.sd="log.hall")
  for(old.name in names(rep.val.vec)){
    new.name <- rep.val.vec[[old.name]]
    colnames(data.list$inputs)[colnames(data.list$inputs) == old.name] <- new.name
  }
  keep.inputs <- apply(is.finite(data.list$inputs), 2, all)
  data.list$inputs <- data.list$inputs[, keep.inputs, drop=FALSE]
  order.dt <- fread(order.csv)
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
  pred.mat.list <- list()
  size.i.vec <- seq_along(train.size.vec)
  ##size.i.vec <- length(train.size.vec)
  for(size.i in size.i.vec){
    train.size <- train.size.vec[[size.i]]
    train.names <- sort(rownames(set.list$train$inputs)[1:train.size])
    cat(sprintf(
      "%4d / %4d files %4d / %4d trainSize=%d\n",
      order.i, length(order.csv.vec), size.i, length(train.size.vec), train.size))
    X.train <- set.list$train$inputs[train.names,]
    y.train <- set.list$train$outputs[train.names,]
    (finite.limits <- colSums(is.finite(y.train)))
    one.pred <- function(x)rep(x, nrow(set.list$test$inputs))
    na.pred <- one.pred(NA)
    X.logn <- X.train[, "log2.n", drop=FALSE]
    pred.list <- list(
      baseline_0=one.pred(baseline.df$pred[train.size]),
      unsup_BIC_1=as.numeric(set.list$test$inputs[, "log2.n"]),
      unreg_linear_1=if(all(0 < finite.limits) && 1 < length(table(X.logn))){
        fit <- penaltyLearning::IntervalRegressionUnregularized(
          X.logn, y.train)
        as.numeric(fit$predict(set.list$test$inputs))
      }else{
        na.pred
      },
      unreg_linear_2=if(all(0 < finite.limits)){
        fit <- penaltyLearning::IntervalRegressionUnregularized(
          X.train[, c("log2.n", "log.hall")], y.train,
          verbose=0)
        as.numeric(fit$predict(set.list$test$inputs))
      }else{
        na.pred
      },
      L1reg_linear_all=if(any(finite.limits < 2) || nrow(y.train) < 4){
        na.pred
      }else{
        n.folds <- min(finite.limits, 5)
        min.col <- which.min(finite.limits)
        is.finite.min <- is.finite(y.train[, min.col])
        fold.vec <- rep(NA, l=nrow(y.train))
        set.seed(1)
        fold.vec[is.finite.min] <- sample(rep(1:n.folds, l=sum(is.finite.min)))
        fold.vec[!is.finite.min] <- sample(rep(1:n.folds, l=sum(!is.finite.min)))
        fit <- penaltyLearning::IntervalRegressionCV(
          X.train, y.train, min.observations=train.size,
          verbose=0,
          fold.vec=fold.vec)
        as.numeric(fit$predict(set.list$test$inputs))
      })
    for(model in names(pred.list)){
      pred.vec <- pred.list[[model]]
      if(is.null(pred.mat.list[[model]])){
        pred.mat.list[[model]] <- matrix(
          NA, length(pred.vec), length(train.size.vec),
          dimnames=list(
            sequenceID=rownames(set.list$test$inputs),
            train.size=train.size.vec))
      }
      stopifnot(length(pred.vec)==nrow(set.list$test$outputs))
      pred.mat.list[[model]][, paste(train.size)] <- pred.vec
      res.vec <- if(is.na(pred.vec[1])){
        NA
      }else{
        penaltyLearning::targetIntervalResidual(
          set.list$test$outputs, pred.vec)
      }
      result.list[[paste(size.i, model)]] <- print(data.table(
        train.size, model,
        accuracy.percent=mean(res.vec==0)*100))
    }#for(model
  }#for(size.i
  for(model in names(pred.mat.list)){
    pred.mat <- pred.mat.list[[model]]
    pred.dt <- data.table(sequenceID=rownames(pred.mat), pred.mat)
    pred.csv <- file.path(seed.dir, "models", model, "predictions.csv")
    dir.create(dirname(pred.csv), showWarnings=FALSE, recursive=TRUE)
    fwrite(pred.dt, pred.csv)
  }
  (result <- do.call(rbind, result.list))
  ## fwrite(result, baseline.csv)
})


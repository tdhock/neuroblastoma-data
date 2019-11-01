source("packages.R")

testFold.dir.vec <- Sys.glob(file.path(
  "data", "*", "cv", "*", "testFolds",
  "*"))

n.pred.vec <- sapply(testFold.dir.vec, function(testFold.dir){
  path <- file.path(dirname(testFold.dir), "models", "*", "predictions.csv")
  length(Sys.glob(path))
})
table(n.pred.vec)

OneFold <- function(testFold.dir){
  library(data.table)
  test.fold <- as.integer(basename(testFold.dir))
  cv.type.dir <- dirname(dirname(testFold.dir))
  data.dir <- dirname(dirname(cv.type.dir))
  folds.csv <- file.path(cv.type.dir, "folds.csv")
  folds.dt <- fread(folds.csv)
  data.list <- list()
  for(data.type in c("inputs", "outputs")){
    csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
    dt <- fread(cmd=paste("xzcat", csv.xz))
    stopifnot(nrow(dt) == nrow(folds.dt))
    m <- as.matrix(dt[, -1, with=FALSE])
    rownames(m) <- dt$sequenceID
    data.list[[data.type]] <- m
  }
  rep.val.vec <- c(
    log.log.bases="log2.n",
    n.loglog="log2.n",
    "diff abs.identity.quantile.50%"="log.hall",
    log.sd="log.hall")
  for(old.name in names(rep.val.vec)){
    new.name <- rep.val.vec[[old.name]]
    colnames(data.list$inputs)[colnames(data.list$inputs) == old.name] <- new.name
  }
  keep.inputs <- apply(is.finite(data.list$inputs), 2, all)
  data.list$inputs <- data.list$inputs[, keep.inputs, drop=FALSE]
  id.list <- list(
    train=folds.dt[fold != test.fold, sequenceID],
    test=folds.dt[fold == test.fold, sequenceID])
  set.list <- list()
  for(set.name in names(id.list)){
    set.id.vec <- id.list[[set.name]]
    set.list[[set.name]] <- lapply(data.list, function(m){
      m[set.id.vec,]
    })
  }

  result.list <- list()
  pred.mat.list <- list()
  
  ####
  
  fit.list <- list()
  scale.i.list <- list( estimated = list( init= NA, estimate = TRUE),
                        fixed = list( init= 1, estimate = FALSE) )

  X.train <- matrix(set.list$train$inputs, nrow(set.list$train$inputs), ncol(set.list$train$inputs))
  Y.train <-  matrix( set.list$train$outputs , nrow(set.list$train$outputs) , ncol(set.list$train$outputs))
  
  for( model.type in c( "gaussian", "logistic", "extreme_value")){
    for( scale.type in scale.i.list){
      fit.list[length(fit.list) + 1] <-  cv.iregnet(X.train, Y.train , family = model.type, 
                     scale_init= scale.type$init ,estimate_scale= scale.type$estimate)
    }
  }
  

  ######
  
  for(size.i in size.i.vec){
    train.size <- train.size.vec[[size.i]]
    maybe.both.inf <- set.list$train$outputs[1:train.size, ]
    not.both.inf <- apply(is.finite(maybe.both.inf), 1, any)
    train.names <- sort(names(not.both.inf)[not.both.inf])
    y.train <- set.list$train$outputs[train.names, , drop=FALSE]
    X.train <- set.list$train$inputs[train.names, , drop=FALSE]
    (finite.limits <- colSums(is.finite(y.train)))
    limit.type <- ifelse(
      is.finite(y.train[,1]), ifelse(
        is.finite(y.train[,2]), "both", "lower"), "upper")
    limit.tab <- table(limit.type)
    one.pred <- function(x)rep(x, nrow(set.list$test$inputs))
    na.pred <- one.pred(NA)
    X.logn <- X.train[, "log2.n", drop=FALSE]
    pred.list <- list(
      baseline_0=one.pred(baseline.df$pred[train.size]),
      unsup_BIC_1=as.numeric(set.list$test$inputs[, "log2.n"]),
      unreg_linear_1=if(all(0 < finite.limits) && 1 < length(table(X.logn))){
        fit <- penaltyLearning::IntervalRegressionUnregularized(
          X.logn, y.train, verbose=0)
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
      L1reg_linear_all=if(
        any(finite.limits < 2) || any(limit.tab < 2) || nrow(y.train) < 4){
        na.pred
      }else{
        n.folds <- min(limit.tab, 5, floor(nrow(y.train)/2))
        fold.vec <- rep(NA, l=nrow(y.train))
        set.seed(1)
        for(l in names(limit.tab)){
          is.l <- limit.type == l
          fold.vec[is.l] <- sample(rep(1:n.folds, l=sum(is.l)))
        }
        fit <- penaltyLearning::IntervalRegressionCV(
          X.train, y.train, min.observations=nrow(y.train),
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
      result.list[[paste(size.i, model)]] <- data.table(
        train.size, model,
        accuracy.percent=with(set.list$test, {
          mean(outputs[,1] < pred.vec & pred.vec < outputs[,2])*100
        }))
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
}

pred.not.done <- order.csv.vec[n.pred.vec==0]
future::plan("multiprocess")
results <- future.apply::future_lapply(pred.not.done, OneSeed)

folds.csv.vec <- Sys.glob("data/systematic/cv/*/folds.csv")
consistent.dt.list <- list()
for(folds.i in seq_along(folds.csv.vec)){
  folds.csv <- folds.csv.vec[[folds.i]]
  folds.dt <- fread(folds.csv)
  type.dir <- dirname(folds.csv)
  cv.type <- basename(type.dir)
  consistent.dt.list[[paste(folds.i)]] <- folds.dt[, {
    test.seqIDs <- sequenceID
    fold.pred.csv <- Sys.glob(file.path(
      type.dir, "testFolds", fold,
      "*", "*", "models", "*", "predictions.csv"))
    if(length(fold.pred.csv)==0){
      data.table()
    }else{
      data.table(pred.csv=fold.pred.csv)[, {
        print(pred.csv)
        pred.dt <- fread(pred.csv, header=TRUE)
        only.int <- gsub("[^0-9]", "", names(pred.dt)[-1])
        pred.cols <- as.integer(only.int)
        pred.dt[, data.table(
          cv.type,
          test.seqIDs=length(test.seqIDs),
          pred.seqIDs=length(sequenceID),
          pred.in.test=sum(sequenceID %in% test.seqIDs),
          pred.not.in.test=sum(!sequenceID %in% test.seqIDs),
          pred.cols=length(pred.cols),
          pred.col.min=min(pred.cols),
          pred.col.max=max(pred.cols)
        )]
      }, by=.(pred.csv)]
    }
  }, by=.(fold)]
}
consistent.dt <- do.call(rbind, consistent.dt.list)

(acc.to.del <- consistent.dt[test.seqIDs != pred.seqIDs, sub("predictions.csv", "accuracy.csv", pred.csv)])
glob.to.del <- file.path(unique(dirname(dirname(acc.to.del))), "*", "accuracy.csv")
unlink(glob.to.del)
Sys.glob(glob.to.del)
consistent.dt[test.seqIDs != pred.seqIDs & grepl("unsup_BIC_1", pred.csv)]

consistent.dt[pred.in.test != test.seqIDs]




## fwrite(result, iregnet.csv)

source("packages.R")

testFold.dir.vec <- Sys.glob(file.path(
  "data", "*", "cv", "*", "testFolds",
  "*"))

n.pred.vec <- sapply(testFold.dir.vec, function(testFold.dir){
  path <- file.path(dirname(testFold.dir), "models", "*", "predictions.csv")
  length(Sys.glob(path))
})
table(n.pred.vec)

# test data dir for OneFold function
testFold.dir <- testFold.dir.vec[[1]]

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
  pred.list <- list()
  
  ####
  
  fit.list <- list()
  scale.i.list <- list( estimated = list( init= NA, estimate = TRUE),
                        fixed = list( init= 1, estimate = FALSE) )

  X.train <- matrix(set.list$train$inputs, nrow(set.list$train$inputs), ncol(set.list$train$inputs))
  Y.train <-  matrix( set.list$train$outputs , nrow(set.list$train$outputs) , ncol(set.list$train$outputs))
  
  
  for( model.type in c( "gaussian", "logistic")){
    for( scale.type in scale.i.list){
      stopifnot( nrow(X.train) == nrow(Y.train))
      fit.list[[length(fit.list) + 1]] <-  cv.iregnet(X.train, Y.train , family = model.type, 
                     scale_init= scale.type$init ,estimate_scale= scale.type$estimate)
    }
  }
  names(fit.list) <- c( "gaus_est" , "gaus_fixed" , "log_est" , "log_fixed")
  
  for( fit in fit.list){
    pred.list[[length(pred.list) + 1]] <- predict(fit , set.list$test$inputs)
  }
  names(pred.list) <- c( "iregnet_gaus_est" , "iregnet_gaus_fixed" , "iregnet_log_est" , "iregnet_log_fixed")

  
    for(model in names(pred.list)){
      pred.vec <- pred.list[[model]]
      if(is.null(pred.mat.list[[model]])){
        pred.mat.list[[model]] <- matrix(
          NA, length(pred.vec), ncol(pred.list[[model]]),
          dimnames=list(
            sequenceID=rownames(set.list$test$inputs) , colnames(pred.list[[model]])))
      }
      stopifnot(length(pred.vec)==nrow(set.list$test$outputs))
      pred.mat.list[[model]] <- pred.vec
      result.list[[paste(model)]] <- data.table(
         model,
        accuracy.percent=with(set.list$test, {
          mean(outputs[,1] < pred.vec & pred.vec < outputs[,2])*100
        }))
    }
  for(model in names(pred.mat.list)){
    pred.mat <- pred.mat.list[[model]]
    pred.dt <- data.table(sequenceID=rownames(pred.mat), pred.log.lambda=pred.mat[,1])
    for( fold in c( 1 , 2 , 3 , 4 , 5)){
    pred.csv <- file.path(testFold.dir,"randomTrainOrderings", fold , "models", model, "predictions.csv")
    dir.create(dirname(pred.csv), showWarnings=FALSE, recursive=TRUE)
    fwrite(pred.dt, pred.csv)
   }
  }
  (result <- do.call(rbind, result.list))
  n.pred.vec[[testFold.dir]] <- 1
}

pred.not.done <- testFold.dir.vec[n.pred.vec==0]
future::plan("multiprocess")
results <- future.apply::future_lapply(pred.not.done, OneFold)

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

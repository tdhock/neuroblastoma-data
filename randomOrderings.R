source("packages.R")

data.dir <- file.path("data", "systematic")
data.list <- list()
for(data.type in c("inputs", "outputs")){
  csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
  dt <- fread(cmd=paste("xzcat", csv.xz))
  m <- as.matrix(dt[, -1, with=FALSE])
  rownames(m) <- dt$sequenceID
  data.list[[data.type]] <- m
}

n.seeds <- 5
folds.csv.vec <- Sys.glob(file.path(
  data.dir, "cv", "*", "folds.csv"))
for(folds.csv in folds.csv.vec){
  fold.dt <- fread(folds.csv)
  cv.type.dir <- dirname(folds.csv)
  for(test.fold in unique(fold.dt$fold)){
    fold.dt[, set := ifelse(
      fold==test.fold, "test", "train")]
    set.list <- list()
    for(set.name in unique(fold.dt$set)){
      id.vec <- fold.dt[set.name, sequenceID, on=list(set)]
      set.list[[set.name]] <- lapply(data.list, function(m){
        m[id.vec,]
      })
    }
    for(seed in 1:n.seeds){
      set.seed(seed)
      seed.dir <- file.path(
        cv.type.dir, "testFolds",
        test.fold, "randomTrainOrderings", seed)
      ord.train.outputs <- with(set.list$train, outputs[sample(nrow(outputs)),])
      order.dt <- data.table(sequenceID=rownames(ord.train.outputs))
      dir.create(seed.dir, recursive=TRUE, showWarnings=FALSE)
      order.csv <- file.path(seed.dir, "order.csv")
      fwrite(order.dt, order.csv)
    }#for(seed
  }#for(test.fold
}#for(folds.csv

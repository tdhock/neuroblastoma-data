source("packages.R")

n.seeds <- 5
folds.csv.vec <- Sys.glob(file.path(
  "data", "*", "cv", "*", "folds.csv"))
for(folds.csv in folds.csv.vec){
  fold.dt <- fread(folds.csv)
  cv.type.dir <- dirname(folds.csv)
  for(test.fold in unique(fold.dt$fold)){
    fold.dt[, set := ifelse(
      fold==test.fold, "test", "train")]
    for(seed in 1:n.seeds){
      order.csv <- file.path(seed.dir, "order.csv")
      if(!file.exists(order.csv)){
        set.seed(seed)
        seed.dir <- file.path(
          cv.type.dir, "testFolds",
          test.fold, "randomTrainOrderings", seed)
        order.dt <- fold.dt[fold!=test.fold, .(
          sequenceID
        )][sample(.N)]
        dir.create(seed.dir, recursive=TRUE, showWarnings=FALSE)
        fwrite(order.dt, order.csv)
      }
    }#for(seed
  }#for(test.fold
}#for(folds.csv

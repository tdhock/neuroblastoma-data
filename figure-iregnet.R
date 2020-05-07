source("packages")
library(iregnet)

(accuracy.csv.vec <- Sys.glob(
  "data/*/cv/equal_labels/testFolds/*/randomTrainOrderings/1/models/*/accuracy.csv"))

csv.i <-1

result.list <- list()
for( csv.i in seq_along(auc.csv.vec)){
  accuracy.csv <- accuracy.csv.vec[[csv.i]]
  model.dir <- dirname(accuracy.csv)
  model <- basename(seed.dir)
  test.fold.dir <- dirname(dirname(dirname(dirname(model.dir))))
  test.fold <- as.integer(basename(test.fold.dir))
  cv.type.dir <- dirname(dirname(testFold.dir))
  data.dir <- dirname(dirname(cv.type.dir))
  data <- basename(data.dir)
  dt <- fread(accuracy.csv)
  result.list[[csv.i]] <- data.table(
  
  )
}


#for(csv.i in seq_along(baseline.csv.vec)){
  baseline.csv <- baseline.csv.vec[[csv.i]]
  seed.dir <- dirname(baseline.csv)
  seed <- as.integer(basename(seed.dir))
  split.dir <- dirname(dirname(seed.dir))
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  dt <- fread(baseline.csv)
  result.list[[csv.i]] <- data.table(
    seed, test.fold, cv.type=basename(cv.type.dir), dt)
#}
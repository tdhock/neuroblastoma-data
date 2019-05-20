source("packages.R")

data.list <- list()
for(data.type in c("possible_errors", "errors", "features", "targets", "folds")){
  f <- paste0(
    "../feature-learning-benchmark/labeled_problems_", data.type, ".csv")
  data.list[[data.type]] <- fread(f)
}

## TODO copy into new standard format.

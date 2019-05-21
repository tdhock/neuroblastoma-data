source("packages.R")

data.list <- list()
for(data.type in c("possible_errors", "errors", "features", "targets", "folds")){
  f <- paste0(
    "../feature-learning-benchmark/labeled_problems_", data.type, ".csv")
  data.list[[data.type]] <- fread(f)
}

## TODO copy into new standard format.
data.list$folds
names(data.list$features)
names(data.list$targets)
sys.outs <- fread("xzcat data/systematic/outputs.csv.xz")
names(data.list$features)[1] <- names(sys.outs)[1]
names(data.list$targets) <- names(sys.outs)

sys.err <- fread("xzcat data/systematic/errors.csv.xz")
names(sys.err)

err.poss <- with(data.list, errors[possible_errors, on=list(prob.dir)])
names(err.poss)

data.list$new.errors <- err.poss[, list(
  sequenceID=prob.dir,
  min.log.lambda=min.log.penalty,
  max.log.lambda=max.log.penalty,
  possible.fp,
  fp,
  possible.fn=possible.tp,
  fn,
  labels,
  errors)]

type.vec <- c(
  "features"="inputs",
  new.errors="evaluation",
  "targets"="outputs")
for(data.type in names(type.vec)){
  prefix <- type.vec[[data.type]]
  dt <- data.list[[data.type]]
  dt[, {
    out.path <- file.path("data", set, paste0(prefix, ".csv"))
    dir.create(dirname(out.path), showWarnings=FALSE, recursive=TRUE)
    d <- data.table(sequenceID, .SD)
    fwrite(d[order(sequenceID)], out.path)
    unlink(csv.xz <- paste0(out.path, ".xz"))
    system(paste("xz", out.path))
    print(head(fread(cmd=paste("xzcat", csv.xz))))
  }, by=list(set=sub("/.*", "", sequenceID))]
}

(set.targets <- data.table(data.list$targets))
set.targets[, set.name := sub("/.*", "", sequenceID)]
set.targets[, problem := sub(".*/", "", sequenceID)]
fold.targets <- set.targets[data.list$folds, on=list(set.name, problem)]
fold.targets[, {
  out.path <- file.path("data", set.name, "cv", "equal_labels", "folds.csv")
  dir.create(dirname(out.path), showWarnings=FALSE, recursive=TRUE)
  d <- data.table(sequenceID, fold)
  fwrite(d[order(sequenceID)], out.path)
  print(head(fread(out.path)))
}, by=list(set.name)]

set.name <- "ATAC_JV_adipose"
rbind(
  length(fread(cmd=paste("xzcat", file.path("data", set.name, "outputs.csv.xz")))[[1]]),
  length(fread(cmd=paste("xzcat", file.path("data", set.name, "inputs.csv.xz")))[[1]]),
  length(fread(file.path("data", set.name, "cv", "equal_labels", "folds.csv"))[[1]]))

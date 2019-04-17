source("packages.R")

data(neuroblastoma, package="neuroblastoma")
big.dt <- data.table(neuroblastoma$profiles)[, list(data=.N), by=list(profile.id, chromosome)][, list(min.data=min(data), max.data=max(data)), by=list(profile.id)][max.data==max(max.data)]

labels.xz.vec <- Sys.glob("data/*/labels.csv.xz")
for(set.i in seq_along(labels.xz.vec)){
  labels.xz <- labels.xz.vec[[set.i]]
  labels.dt <- fread(cmd=paste("xzcat", labels.xz))
  head(match.dt <- namedCapture::df_match_variable(
    labels.dt,
    sequenceID=list(
      profileID="[0-9]+",
      "_",
      chrom="chr.*")))
  table(match.dt$sequenceID.chrom)
  randcol <- function(dt, col.name, n.folds=6){
    unique.folds <- 1:n.folds
    col.vec <- dt[[col.name]]
    u.vec <- unique(col.vec)
    fold <- sample(rep(unique.folds, l=length(u.vec)))
    names(fold) <- u.vec
    fold[paste(col.vec)]
  }
  fun.list <- list(
    profileSize=function(dt){
      ifelse(
        dt$sequenceID.profileID %in% big.dt$profile.id, 1, 2)
    }, chrom=function(dt){
      as.integer(factor(dt$sequenceID.chrom))
    }, profileID=function(dt){
      randcol(dt, "sequenceID.profileID")
    }, sequenceID=function(dt){
      randcol(dt, "sequenceID")
    })
  for(split.name in names(fun.list)){
    fun <- fun.list[[split.name]]
    set.seed(1)
    fold.vec <- fun(match.dt)
    print(table(fold.vec))
    cv.dir <- file.path(dirname(labels.xz), "cv", split.name)
    label.folds <- labels.dt[, data.table(
      sequenceID, fold=fold.vec)]
    fold.counts <- label.folds[, list(
      folds=length(unique(fold))
    ), by=list(sequenceID)]
    bad <- fold.counts[folds != 1]
    if(nrow(bad)){
      print(bad)
      stop("some sequenceID numbers appear in more than one fold")
    }
    u.folds <- unique(label.folds)[order(sequenceID)]
    dir.create(cv.dir, showWarnings=FALSE, recursive=TRUE)
    print(folds.csv <- file.path(cv.dir, "folds.csv"))
    fwrite(u.folds, folds.csv)
  }
}

source("packages.R")

data(neuroblastoma, package="neuroblastoma")
big.dt <- data.table(neuroblastoma$profiles)[, list(data=.N), by=list(profile.id, chromosome)][, list(min.data=min(data), max.data=max(data)), by=list(profile.id)][max.data==max(max.data)]

labels.xz.vec <- Sys.glob("data/*/labels.csv.xz")
N.folds <- 6
for(set.i in seq_along(labels.xz.vec)){
  labels.xz <- labels.xz.vec[[set.i]]
  labels.cmd <- paste("xzcat", labels.xz)
  labels.dt <- fread(cmd=labels.cmd)
  prob.dt <- labels.dt[, list(
    labels=.N
  ), by=list(sequenceID)]
  eval.cmd <- sub("labels", "evaluation", labels.cmd)
  eval.dt <- fread(cmd=eval.cmd)
  head(match.dt <- namedCapture::df_match_variable(
    prob.dt,
    sequenceID=list(
      profileID="[0-9]+",
      "_",
      chrom="chr.*")))
  table(match.dt$sequenceID.chrom)
  randcol <- function(dt, col.name, n.folds=N.folds){
    unique.folds <- 1:n.folds
    col.vec <- dt[[col.name]]
    u.vec <- unique(col.vec)
    fold <- sample(rep(unique.folds, l=length(u.vec)))
    names(fold) <- u.vec
    fold[paste(col.vec)]
  }
  fun.list <- list(
    chrom=function(dt){
      as.integer(factor(dt$sequenceID.chrom))
    },
    profileSize=function(dt){
      randcol(dt, "sequenceID.profileID", N.folds/2)+3*ifelse(
        dt$sequenceID.profileID %in% big.dt$profile.id, 0, 1)
    },
    profileID=function(dt){
      randcol(dt, "sequenceID.profileID")
    },
    sequenceID=function(dt){
      randcol(dt, "sequenceID")
    })
  for(split.name in names(fun.list)){
    fun <- fun.list[[split.name]]
    set.seed(1)
    fold.vec <- fun(match.dt)
    print(table(fold.vec))
    cv.dir <- file.path(dirname(labels.xz), "cv", paste0("R-3.6.0-", split.name))
    prob.folds <- prob.dt[, data.table(
      sequenceID, fold=fold.vec)]
    fold.counts <- prob.folds[, list(
      folds=length(unique(fold))
    ), by=list(sequenceID)]
    bad <- fold.counts[folds != 1]
    if(nrow(bad)){
      print(bad)
      stop("some sequenceID numbers appear in more than one fold")
    }
    print(auc.dt <- prob.folds[, {
      pred.dt <- data.table(sequenceID, pred.log.lambda=0)
      L <- penaltyLearning::ROChange(eval.dt, pred.dt, "sequenceID")
      p <- L$thresholds[threshold=="predicted"]
      list(auc=L$auc, possible.fn=p$possible.fn, possible.fp=p$possible.fp)
    }, by=list(fold)])
    u.folds <- unique(prob.folds)[order(sequenceID)]
    dir.create(cv.dir, showWarnings=FALSE, recursive=TRUE)
    print(folds.csv <- file.path(cv.dir, "folds.csv"))
    fwrite(u.folds, folds.csv)
  }
}

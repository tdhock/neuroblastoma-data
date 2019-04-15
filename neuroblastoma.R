source("packages.R")

data(neuroblastomaDetailed, package="bams")
data(neuroblastoma, package="neuroblastoma")
data(neuroblastomaProcessed, package="penaltyLearning")

writeXZ <- function(set.name, ...){
  out.list <- list(...)
  for(out.type in names(out.list)){
    out.dt <- out.list[[out.type]]
    out.path <- file.path(
      "data", set.name, paste0(out.type, ".csv"))
    dir.create(dirname(out.path), showWarnings=FALSE, recursive=TRUE)
    fwrite(out.dt, out.path)
    system(paste("xz", out.path))
    print(out.path)
  }
}

new.sep <- "_chr"
seqID <- function(profile.id, chromosome){
  paste0(profile.id, new.sep, chromosome)
}
unlink("data", recursive=TRUE)
label.list <- list(
  detailed=neuroblastomaDetailed,
  systematic=neuroblastoma$annotations)
all.profiles <- data.table(neuroblastoma$profiles)
for(label.name in names(label.list)){
  label.dt <- penaltyLearning::change.labels[data.table(
    label.list[[label.name]]), on=list(annotation)]
  label.ids <- unique(label.dt[, .(profile.id, chromosome)])
  labeled.profiles <- all.profiles[label.ids, on=list(profile.id, chromosome)]
  writeXZ(
    label.name,
    profiles=labeled.profiles[, data.table(
      sequenceID=seqID(profile.id, chromosome),
      position, signal=logratio)],
    labels=label.dt[, data.table(
      sequenceID=seqID(profile.id, chromosome),
      labelStart=min,
      labelEnd=max,
      annotation,
      max.changes,
      min.changes,
      color,
      possible.fn,
      possible.fp)])
}

newID <- function(m){
  sub("[.]", new.sep, rownames(m))
}
with(neuroblastomaProcessed, writeXZ(
  "systematic",
  errors=errors[, data.table(
  sequenceID=seqID(profile.id, chromosome),
  n.segments,
  possible.fp, fp,
  possible.fn, fn,
  labels, errors,
  min.log.lambda, max.log.lambda,
  loss)],
  outputs=data.table(
    sequenceID=rownames(target.mat),
    min.log.lambda=target.mat[, "min.L"],
    max.log.lambda=target.mat[, "max.L"]),
  inputs=data.table(
    sequenceID=rownames(feature.mat),
    feature.mat)))
glob <- "data/*/*"
system(paste("du -ms", glob))
for(path in Sys.glob(glob)){
  dt <- fread(cmd=paste("xzcat", path))
  print(path)
  str(dt)
  print(dt)
}

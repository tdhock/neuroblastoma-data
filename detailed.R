source("packages.R")

data.dir <- file.path("data", "detailed")
data.list <- list()
for(data.type in c("profiles", "labels")){
  csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
  dt <- fread(cmd=paste("xzcat", csv.xz))
  data.list[[data.type]] <- dt
}

label.count <- data.list$labels[, .(
  labels=.N
), by=list(sequenceID)]

labeled.profiles <- data.list$profiles[label.count, on=list(sequenceID)]

loss.dt <- labeled.profiles[, {
  max.segments <- 20L
  fit <- jointseg::Fpsn(signal, max.segments)
  change.list <- list()
  seg.vec <- 1:max.segments
  for(model.i in seg.vec){
    end.vec <- fit$t.est[model.i, 1:model.i]
    before <- end.vec[-model.i]
    after <- before+1L
    change.list[[model.i]] <- floor((position[before]+position[after])/2)+0.5
  }
  print(.GRP)
  data.table(segments=seg.vec, change.list, loss=fit$J.est)
}, by=list(sequenceID)]

change.dt <- loss.dt[, {
  data.table(change=change.list[[1]])
}, by=list(sequenceID, segments)]
model.dt <- loss.dt[, {
  dt <- data.table(segments, loss)
  penaltyLearning::modelSelection(dt, "loss", "segments")
}, by=list(sequenceID)]

setkey(data.list$labels, sequenceID, labelStart)
label.dt <- data.list$labels[, {
  next.start <- labelStart[-1]
  end <- labelEnd[-.N]
  to.rep <- next.start < end
  middle <- floor((next.start+end)/2)
  next.start[to.rep] <- middle[to.rep]
  end[to.rep] <- middle[to.rep]
  data.table(
    labelStart=c(labelStart[1], next.start),
    labelEnd=c(end, labelEnd[.N]),
    annotation)
}, by=list(sequenceID)]
error.list <- penaltyLearning::labelError(
  model.dt, label.dt, change.dt,
  change.var="change",
  label.vars=c("labelStart", "labelEnd"),
  model.vars="segments",
  problem.vars="sequenceID")
interval.dt <- penaltyLearning::targetIntervals(
  error.list$model.errors, "sequenceID")

names(fread(cmd="xzcat data/ATAC_JV_adipose/evaluation.csv.xz"))
names(fread(cmd="xzcat data/ATAC_JV_adipose/inputs.csv.xz"))
names(fread(cmd="xzcat data/ATAC_JV_adipose/outputs.csv.xz"))

fmat <- penaltyLearning::featureMatrix(
  labeled.profiles, "sequenceID", "signal")

out.dt.list <- list(
  outputs=interval.dt[, .(sequenceID, min.log.lambda, max.log.lambda)],
  inputs=data.table(sequenceID=rownames(fmat), fmat),
  evaluation=error.list$model.errors[, .(
    sequenceID, min.log.lambda, max.log.lambda,
    possible.fp, fp,
    possible.fn, fn,
    labels, errors)])
for(out.name in names(out.dt.list)){
  out.dt <- out.dt.list[[out.name]][order(sequenceID)]
  out.path <- file.path(data.dir, paste0(out.name, ".csv"))
  fwrite(out.dt, out.path)
  unlink(csv.xz <- paste0(out.path, ".xz"))
  system(paste("xz", out.path))
  print(head(fread(cmd=paste("xzcat", csv.xz))))
}

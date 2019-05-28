### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,...){
  local.lib <- file.path(getwd(), "library")
  dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
  .libPaths(local.lib)
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
options(repos=c(
          "http://www.bioconductor.org/packages/release/bioc",
          "http://r-forge.r-project.org",
          "http://cloud.r-project.org",
          "http://cran.r-project.org"))
##library(PeakSegDisk)
works_with_R(
  "3.6.0",
  R.utils="2.8.0",
  data.table="1.12.2",
  RJSONIO="1.3.1.1",
  ##PeakSegDisk="2018.11.28",
  PeakSegJoint="2018.10.3",
  xtable="1.8.4",
  "tdhock/PeakSegDisk@671f888253b6fe1a753ebb0ba28069d2eeb30211",
  "tdhock/PeakSegPipeline@f880a978ab15756c5cacab7fcb0668bfba3b9ac5",
  "tdhock/animint2@7ac92dcc31ea35af6a6cd016e352e7c980e08198")
future::plan("multiprocess")

eval.dt <- fread(cmd="xzcat data/ATAC_JV_adipose/evaluation.csv.xz")

(min.dt <- eval.dt[, list(
  min.errors=min(errors)
  ), by=list(sequenceID)])
table(min.dt$min.errors)

min.dt <- eval.dt[, {
  dt <- unique(.SD[min(errors)==errors, .(fp, fn, errors)])
  data.table(dt, models=nrow(dt))
}, by=list(sequenceID)]
table(min.dt$models)
min.dt[models==2]
min.stats <- min.dt[, {
  list(
    min.fp=min(fp),
    max.fp=max(fp),
    min.fn=min(fn),
    max.fn=max(fn)
  )
}, by=list(sequenceID)]
##prob.dir.vec <- Sys.glob("../feature-learning-benchmark/data/*/samples/*/*/problems/*")
##file.rename(prob.dir.vec, sub("[^-./_a-zA-Z0-9]", "-", prob.dir.vec))
labels.dt <- min.dt[, {
  s <- sub(":", "-", sequenceID)
  fread(
    paste0("../feature-learning-benchmark/data/", s, "/labels.bed"),
    col.names=c("chrom", "labelStart", "labelEnd", "annotation"))
}, by=list(sequenceID)]

label.range <- labels.dt[, list(
  min.labelStart=min(labelStart),
  max.labelEnd=max(labelEnd),
  mean=mean(labelEnd-labelStart),
  labels=.N
), by=list(sequenceID)]
label.range[, range := max.labelEnd-min.labelStart]
label.range[order(range)]

some.seqs <- min.stats[1 < max.fn & 1 < max.fp]
some.seqs <- min.stats[label.range, on=list(sequenceID)][0==min.fp & 0==min.fn & 1==max.fn & 1==max.fp][1:2]
show.err <- eval.dt[some.seqs, on=list(sequenceID)]
tall.err <- melt(
  show.err,
  measure.vars=c("fp", "fn", "errors"))

profiles.dt <- some.seqs[, {
  expand <- range/20
  from <- min.labelStart-expand
  to <- max.labelEnd+expand
  s <- sub(":", "-", sequenceID)
  bg <- paste0(
    "../feature-learning-benchmark/data/",
    s, "/coverage.bedGraph")
  gz <- paste0(bg, ".gz")
  dt <- fread(
    gz,
    col.names=c("chrom", "chromStart", "chromEnd", "coverage"))
  fwrite(dt, bg, col.names=FALSE, sep="\t")
  dt[from < chromEnd & chromStart < to]
}, by=list(sequenceID)]
some.labels <- labels.dt[some.seqs, on=list(sequenceID)]
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
ggplot()+
  ggtitle(
    "Noisy coverage data and labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sequenceID ~ ., scales="free")+
  geom_tallrect(aes(
    xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
    data=some.labels,
    color="grey")+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    data=profiles.dt,
    color="grey50")+
  scale_x_continuous(breaks=seq(4e4, 5e4, by=5))

win <- function(windowStart, windowEnd){
  data.table(windowStart, windowEnd)
}
win.dt <- rbind(
  win(43447, 43457),
  win(43502, 43512))*1000
win.dt[, window := 1:.N]
setkey(profiles.dt, chromStart, chromEnd)
setkey(some.labels, labelStart, labelEnd)
setkey(win.dt, windowStart, windowEnd)
win.profiles <- foverlaps(profiles.dt, win.dt, nomatch=0L)
win.labels <- foverlaps(some.labels, win.dt, nomatch=0L)
gg <- ggplot()+
  ggtitle(
    "Noisy coverage data and labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sequenceID ~ window, scales="free", labeller=label_both)+
  geom_tallrect(aes(
    xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
    data=win.labels,
    color="grey")+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart/1e3, coverage),
    data=win.profiles,
    color="grey50")+
  scale_x_continuous(breaks=seq(4e4, 5e4, by=5))
print(gg)

win.err.list <- list()
win.segs.list <- list()
for(seq.i in 1:nrow(some.seqs)){
  s <- some.seqs[seq.i]
  pdir <- paste0(
    "../feature-learning-benchmark/data/",
    sub(":", "-", s$sequenceID))
  L <- PeakSegPipeline::problem.target(pdir, 1)
  plabels <- win.labels[sequenceID==s$sequenceID]
  plabels[, chromStart := labelStart]
  plabels[, chromEnd := labelEnd]
  selection.dt <- data.table(penaltyLearning::modelSelection(
    L$models, "total.loss", "segments"))
  for(model.i in 1:nrow(selection.dt)){
    model <- selection.dt[model.i]
    pen.str <- paste(model$penalty)
    pen.info <- PeakSegDisk::problem.PeakSegFPOP(pdir, pen.str)
    seg.dt <- data.table(sequenceID=s$sequenceID, model, pen.str, pen.info$segments)
    setkey(seg.dt, chromStart, chromEnd)
    over.dt <- foverlaps(seg.dt, win.dt, nomatch=0L)
    peak.dt <- over.dt[status=="peak"]
    e <- PeakError::PeakErrorChrom(peak.dt, plabels)
    win.err.list[[paste(seq.i, model.i)]] <-
      data.table(
        sequenceID=s$sequenceID,
        window=plabels$window,
        model, pen.str, e)
    win.segs.list[[paste(seq.i, model.i)]] <- over.dt
  }
}
win.segs <- do.call(rbind, win.segs.list)
win.err <- do.call(rbind, win.err.list)
win.segs[, segStart := ifelse(chromStart<windowStart, windowStart, chromStart)]
win.segs[, segEnd := ifelse(windowEnd<chromEnd, windowEnd, chromEnd)]

auc.dt.list <- list()
roc.dt.list <- list()
err.dt.list <- list()
roc.segs.list <- list()
roc.win.err.list <- list()
off.by <- 0.1
for(offset in seq(-5, 5, by=off.by)){
  pred.dt <- data.table(some.seqs, pred.log.lambda=c(0, offset))
  pred.eval <- eval.dt[pred.dt, on=list(sequenceID)]
  pred.eval[, min.thresh := min.log.lambda-pred.log.lambda]
  pred.eval[, max.thresh := max.log.lambda-pred.log.lambda]
  pred.eval[, piece := 1:.N]
  err.dt.list[[paste(offset)]] <- data.table(offset, pred.eval)
  print(offset)
  roc <- penaltyLearning::ROChange(eval.dt, pred.dt, "sequenceID")
  roc$roc[, thresh := (min.thresh+max.thresh)/2]
  pred.some.cols <- pred.dt[, list(id=1, sequenceID, pred.log.lambda)]
  roc.off.id <- data.table(offset, id=1, roc$roc)
  roc.off <- roc.off.id[pred.some.cols, on=list(
    id), allow.cartesian=TRUE]
  roc.off[, log.lambda := thresh + pred.log.lambda]
  roc.segs.list[[paste(offset)]] <-
    win.segs[roc.off, nomatch=0L, on=list(
      sequenceID,
      min.log.lambda<log.lambda,
      max.log.lambda>log.lambda)]
  roc.win.err.list[[paste(offset)]] <-
    win.err[roc.off, nomatch=0L, on=list(
      sequenceID,
      min.log.lambda<log.lambda,
      max.log.lambda>log.lambda)]
  off.min <- roc$roc[errors==min(errors)]
  auc.dt.list[[paste(offset)]] <- with(roc, data.table(
    auc, offset,
    min.errors=off.min$errors[1],
    n.min=nrow(off.min),
    thresholds[threshold=="min.error"]))
  roc.dt.list[[paste(offset)]] <- data.table(
    offset, roc$roc, piece=1:nrow(roc$roc), sequenceID="Total")
}
auc.dt <- do.call(rbind, auc.dt.list)
roc.dt <- do.call(rbind, roc.dt.list)
roc.segs <- do.call(rbind, roc.segs.list)
err.dt <- do.call(rbind, err.dt.list)
roc.win.err.dt <- do.call(rbind, roc.win.err.list)

ggplot()+
  geom_point(aes(
    offset, auc),
    data=auc.dt)

common.names <- intersect(names(roc.dt), names(err.dt))
both.dt <- rbind(
  err.dt[, common.names, with=FALSE],
  roc.dt[, common.names, with=FALSE])
err.dt.tall <- melt(
  both.dt,
  variable.name="error.type",
  measure.vars=c("fp", "fn", "errors"))
id2show <- function(seqID)gsub(
  "ATAC_JV_adipose/samples/AC1/|/problems/chrX:37148256-49242997", "",
  seqID)
roc.segs[, showSeq := id2show(sequenceID)]
win.labels[, showSeq := id2show(sequenceID)]
roc.win.err.dt[, showSeq := id2show(sequenceID)]
win.profiles[, showSeq := id2show(sequenceID)]
err.dt.tall[, showSeq := id2show(sequenceID)]
auc.dt[, thresh := (min.thresh+max.thresh)/2]
roc.dt[, Errors := ifelse(errors==min(errors), "Min", "More"), by=list(offset)]
min.err <- roc.dt[Errors=="Min"]
min.err[, piece := 1:.N, by=list(offset)]
roc.size <- 5
roc.peaks <- roc.segs[status=="peak"]
viz <- animint(
  title="Changepoint detection ROC curve alignment problem",
  ##first=list(offset=0.5),
  duration=list(offset=250),
  time=list(variable="offset", ms=250),
  profiles=ggplot()+
    ggtitle(
      "Noisy coverage data, labels, and predicted model")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=1200)+
    facet_grid(showSeq ~ window, scales="free", labeller=label_both)+
    geom_tallrect(aes(
      xmin=labelStart/1e3, xmax=labelEnd/1e3, fill=annotation),
      data=win.labels,
      alpha=0.5,
      color="grey")+
    scale_linetype_manual(
      "Error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3, linetype=status),
      data=roc.win.err.dt,
      showSelected=c("offset", "thresh"),
      fill=NA,
      size=2,
      color="black")+
    scale_fill_manual(values=ann.colors)+
    geom_step(aes(
      chromStart/1e3, coverage),
      data=win.profiles,
      color="grey50")+
    geom_segment(aes(
      segStart/1e3, mean,
      xend=segEnd/1e3, yend=mean),
      color="green",
      showSelected=c("offset", "thresh"),
      data=roc.segs)+
    geom_segment(aes(
      segStart/1e3, 0,
      xend=segEnd/1e3, yend=0),
      color="deepskyblue",
      showSelected=c("offset", "thresh"),
      size=5,
      alpha=0.5,
      data=roc.peaks)+
    scale_x_continuous(breaks=seq(4e4, 5e4, by=5)),
  auc=ggplot()+
    ggtitle(
      "AUC, select offset")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_point(aes(
      offset, auc),
      fill=NA,
      data=auc.dt)+
    geom_tallrect(aes(
      xmin=offset-off.by/2, xmax=offset+off.by/2),
      clickSelects="offset",
      alpha=0.5,
      data=auc.dt),
  error=ggplot()+
    ggtitle(
      "Error curves, select threshold")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(showSeq ~ ., scales="free")+
    scale_color_manual(values=c(
      fp="red",
      fn="deepskyblue",
      errors="black"))+
    scale_size_manual(values=c(
      fp=5,
      fn=3,
      errors=1))+
    xlab("Prediction threshold")+
    scale_y_continuous(
      "Number of incorrectly predicted labels",
      breaks=seq(0, 20, by=2))+
    geom_vline(aes(
      xintercept=thresh, key=piece),
      showSelected="offset",
      color="grey",
      data=min.err)+
    geom_text(aes(
      thresh+0.2, labels*0.9, key=1, label=paste0(
        "Min Errors=", errors)),
      showSelected="offset",
      hjust=0,
      color="grey",
      data=data.table(auc.dt, showSeq="Total"))+
    ## TODO geom_point, select min err and n.min.
    geom_segment(aes(
      min.thresh, value,
      key=paste(piece, error.type),
      color=error.type,
      size=error.type,
      xend=max.thresh, yend=value),
      showSelected="offset",
      data=err.dt.tall)+
    geom_tallrect(aes(
      xmin=min.thresh, xmax=max.thresh,
      tooltip=sprintf(
        "%.1f<thresh<%.1f FP=%d FN=%d",
        min.thresh, max.thresh, fp, fn),
      key=paste(offset, thresh)),
      showSelected="offset",
      clickSelects="thresh",
      alpha=0.5,
      data=roc.dt),
  roc=ggplot()+
    ggtitle(
      "ROC curves")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_path(aes(
      FPR, TPR, key=paste(offset, thresh)),
      showSelected="offset",
      data=roc.dt)+
    ## geom_point(aes(
    ##   FPR, TPR, key=paste(FPR, TPR)),
    ##   showSelected="offset",
    ##   data=min.err,
    ##   size=roc.size,
    ##   fill=NA)+
    geom_point(aes(
      FPR, TPR, fill=Errors,
      tooltip=sprintf(
        "%.1f<thresh<%.1f FP=%d FN=%d",
        min.thresh, max.thresh, fp, fn),
      key=paste(offset, thresh)),
      showSelected="offset",
      clickSelects="thresh",
      size=roc.size,
      alpha=0.7,
      data=roc.dt)+
    scale_fill_manual(values=c(
      Min="black",
      More="white"))+
    coord_equal()+
    geom_text(aes(
      0.75, 0.25, key=1, label=sprintf(
        "AUC=%.2f", auc)),
      showSelected="offset",
      data=auc.dt)+
    geom_abline(aes(
      slope=slope, intercept=intercept),
      color="grey",
      data=data.table(slope=1, intercept=0))
  )
animint2dir(viz, "figure-max-auc")

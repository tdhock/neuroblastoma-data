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
    color="grey50")

for(seq.i in 1:nrow(some.seqs)){
  s <- some.seqs[seq.i]
  pdir <- paste0(
    "../feature-learning-benchmark/data/",
    sub(":", "-", s$sequenceID))
  L <- PeakSegPipeline::problem.target(pdir, 1)
}

auc.dt.list <- list()
roc.dt.list <- list()
err.dt.list <- list()
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
  auc.dt.list[[paste(offset)]] <- with(roc, data.table(
    auc, offset,
    thresholds[threshold=="min.error"]))
  roc.dt.list[[paste(offset)]] <- data.table(
    offset, roc$roc, piece=1:nrow(roc$roc), sequenceID="Total")
}
auc.dt <- do.call(rbind, auc.dt.list)
roc.dt <- do.call(rbind, roc.dt.list)
err.dt <- do.call(rbind, err.dt.list)

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
err.dt.tall[, showSeq := gsub(
  "ATAC_JV_adipose/samples/AC1/|/problems/chrX:37148256-49242997", "",
  sequenceID)]
roc.dt[, thresh := (min.thresh+max.thresh)/2]
auc.dt[, thresh := (min.thresh+max.thresh)/2]
roc.dt[, Errors := ifelse(errors==min(errors), "Min", "More"), by=list(offset)]
min.err <- roc.dt[Errors=="Min"]
min.err[, piece := 1:.N, by=list(offset)]
roc.size <- 5
viz <- animint(
  title="Changepoint detection ROC curve alignment problem",
  ##first=list(offset=0.5),
  duration=list(offset=250),
  time=list(variable="offset", ms=250),
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

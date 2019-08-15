library(animint2)
library(data.table)

gp.dirs <- Sys.glob("data/systematic/cv/*/testFolds/*/*GP*")
random.dirs <- unique(file.path(dirname(gp.dirs), "randomTrainOrderings"))
selection.dirs <- c(gp.dirs, random.dirs)
acc.dt.list <- list()
for(sel.dir in selection.dirs){
  acc.csv.vec <- Sys.glob(file.path(
    sel.dir, "*", "models", "*", "accuracy.csv"))
  for(acc.csv in acc.csv.vec){
    model.dir <- dirname(acc.csv)
    models.dir <- dirname(model.dir)
    seed.dir <- dirname(models.dir)
    seed <- as.integer(basename(seed.dir))
    selection.dir <- dirname(seed.dir)
    selection.type <- basename(selection.dir)
    split.dir <- dirname(selection.dir)
    test.fold <- as.integer(basename(split.dir))
    cv.type.dir <- dirname(dirname(split.dir))
    dt <- fread(acc.csv)
    acc.dt.list[[acc.csv]] <- data.table(
      selection.type,
      model=basename(model.dir),
      seed, test.fold, cv.type=basename(cv.type.dir), dt)
  }
}

counts <- sapply(acc.dt.list, function(L)length(names(L)))
table(counts)
unlink(names(counts)[counts<max(counts)])

acc.dt <- do.call(rbind, acc.dt.list)
min.max <- acc.dt[, .(
  min=min(N.test),
  max=max(N.test)
), by=.(cv.type, test.fold)]
stopifnot(min.max[, min==max])

both.stats <- acc.dt[is.finite(percent.correct.intervals) & train.size<=60, list(
  mean=mean(percent.correct.intervals),
  sd=sd(percent.correct.intervals),
  median=median(percent.correct.intervals),
  q25=quantile(percent.correct.intervals, 0.25),
  q75=quantile(percent.correct.intervals, 0.75)
), by=list(selection.type, test.fold, cv.type, train.size, model)]
both.stats[, selection.disp := sub("TrainOrderings|sampleSelection", "", selection.type)]
both.stats[, Selection := ifelse(selection.disp=="random", "random", "GP")]
both.stats[, m := ifelse(model=="GP", selection.disp, model)]
model.levs <- rev(sort(unique(both.stats$m)))
model.colors <- ifelse(grepl("GP", model.levs), "red", "black")
names(model.colors) <- model.levs
both.stats[, Model := factor(m, model.levs)]
both.stats[, type.fold := paste0(cv.type, "/", test.fold)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type + test.fold ~ selection.disp, labeller=label_both)+
  geom_line(aes(
    train.size, median, color=Model, group=Model),
    size=1,
    data=both.stats)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75, fill=Model, group=Model),
    alpha=0.5,
    data=both.stats)+
  ylab("Percent correctly predicted intervals")+
  coord_cartesian(xlim=c(2, 200), ylim=c(50, 100))+
  scale_color_manual(values=model.colors)+
  scale_fill_manual(values=model.colors)+
  scale_x_log10(
    "Labeled sequences in train set")

##animint(gg)

size.counts <- both.stats[, .(
  sizes=.N,
  last.median=median[.N]
), by=.(Selection, cv.type, test.fold, type.fold, Model)]
model.counts <- size.counts[, .(
  models=.N,
  max.last.median=max(last.median)
), by=list(Selection, cv.type, test.fold, type.fold)]
model.wide <- dcast(
  model.counts,
  cv.type + test.fold + type.fold ~ Selection,
  value.var="max.last.median")
model.wide[, improvement := GP-random]

zoomout <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_line(aes(
    train.size, median, color=Model, linetype=selection.disp,
    group=paste(Model, selection.disp)),
    size=1,
    clickSelects="Model",
    showSelected="type.fold",
    data=both.stats)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75, fill=Model,
    group=paste(Model, selection.disp)),
    showSelected=c("type.fold", "selection.disp"),
    clickSelects="Model",
    alpha=0.5,
    data=both.stats)+
  ylab("Percent correctly predicted intervals")+
  scale_color_manual(values=model.colors)+
  scale_fill_manual(values=model.colors)+
  scale_x_log10(
    "Labeled sequences in train set",
    breaks=c(2, 5, 10, 20, 30, 60))
viz <- animint(
  out.dir="figure-2019-08-14-animint",
  overview=ggplot()+
    theme_bw()+
    geom_tile(aes(
      test.fold, cv.type, fill=improvement),
      data=model.wide,
      clickSelects="type.fold")+
    scale_fill_gradient2("GP-random"),
  zoomout=zoomout+
    ##guides(color="none", linetype="none", fill="none")+
    ggtitle("All data for selected test fold"),
  details=zoomout+
    ggtitle("Zoom to top left")+
    coord_cartesian(xlim=c(2, 20), ylim=c(70, 100))
)

animint2gist(viz)

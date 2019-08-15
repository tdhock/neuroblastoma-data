library(animint2)
library(data.table)

gp.paths <- Sys.glob("data/systematic/cv/*/testFolds/*/*GP*")
random.paths <- unique(file.path(dirname(gp.paths), "randomTrainOrderings"))
selection.paths <- c(gp.paths, random.paths)
acc.dt.list <- list()
for(sel.path in selection.paths){
  acc.csv.vec <- Sys.glob(file.path(
    sel.path, "*", "models", "*", "accuracy.csv"))
  for(acc.csv in acc.csv.vec){
    model.path <- dirname(acc.csv)
    models.path <- dirname(model.path)
    seed.path <- dirname(models.path)
    seed <- as.integer(basename(seed.path))
    selection.path <- dirname(seed.path)
    split.path <- dirname(selection.path)
    test.fold <- as.integer(basename(split.path))
    cv.type.path <- dirname(dirname(split.path))
    dt <- fread(acc.csv)
    auc.na <- dt[is.na(auc)]
    if(nrow(auc.na)){
      print(acc.csv)
      print(auc.na)
    }
    acc.dt.list[[acc.csv]] <- data.table(
      selection.dir=basename(selection.path),
      model.dir=basename(model.path),
      seed, test.fold, cv.type=basename(cv.type.path), dt)
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

max.train.size <- 60
both.stats <- acc.dt[
  is.finite(percent.correct.intervals) & train.size<=max.train.size, list(
  mean=mean(percent.correct.intervals),
  sd=sd(percent.correct.intervals),
  median=median(percent.correct.intervals),
  q25=quantile(percent.correct.intervals, 0.25),
  q75=quantile(percent.correct.intervals, 0.75),
  seeds=.N
), by=list(selection.dir, test.fold, cv.type, train.size, model.dir)]
table(both.stats$selection.dir)
both.stats[, selection.short := sub(
  "TrainOrderings|sampleSelection", "", selection.dir)]
table(both.stats$selection.short)
both.stats[, random.or.GP := ifelse(
  selection.short=="random", "random", "GP")]
table(both.stats$random.or.GP)
table(both.stats$model.dir)
both.stats[, model := ifelse(model.dir=="GP", selection.short, model.dir)]
table(both.stats$model)
both.stats[, baseline.or.GP := ifelse(model.dir=="GP", "GP", "baseline")]
table(both.stats$baseline.or.GP)
both.stats[, model.selection := ifelse(
  baseline.or.GP=="GP", "GP", paste0("baseline.", random.or.GP))]
table(both.stats$model.selection)
both.stats[, type.fold := paste0(cv.type, "/", test.fold)]
table(both.stats$type.fold)

GP.models.only <- both.stats[baseline.or.GP=="GP"]
animint(
  GPonly=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(type.fold ~ .)+
    geom_line(aes(
      train.size, median, group=selection.short),
      size=1,
      clickSelects="selection.short",
      data=GP.models.only)+
    geom_ribbon(aes(
      train.size, ymin=q25, ymax=q75, group=selection.short),
      alpha=0.55,
      clickSelects="selection.short",
      data=GP.models.only)+
    ylab("Percent correctly predicted intervals")+
    scale_x_log10(
      "Labeled sequences in train set")
)

## Compute improvement over random for max.train.size.
size.counts <- both.stats[, data.table(
  sizes=.N,
  .SD[train.size==max.train.size]
), by=.(cv.type, test.fold, type.fold, baseline.or.GP, random.or.GP, selection.short)]
model.counts <- size.counts[, data.table(
  .SD[which.max(median)],
  models=.N
), by=list(cv.type, test.fold, type.fold, baseline.or.GP, random.or.GP)]
some.model.counts <- model.counts[random.or.GP=="random" | baseline.or.GP=="GP"]
model.wide <- dcast(
  some.model.counts,
  cv.type + test.fold + type.fold ~ random.or.GP,
  value.var=c("median", "models"))
model.wide[, improvement := median_GP-median_random]
best.baseline <- some.model.counts[baseline.or.GP=="baseline"]
best.baseline.curves <- both.stats[best.baseline, on=.(
  cv.type, test.fold, type.fold, baseline.or.GP, random.or.GP,
  selection.short)]
model.wide[, tooltip := sprintf(
  "%d GP models for %s, best improvement over random=%.1f%%",
  models_GP, type.fold, improvement)]
animint(
  overview=ggplot()+
    theme_bw()+
    geom_tile(aes(
      test.fold, cv.type, fill=improvement, tooltip=tooltip),
      data=model.wide,
      clickSelects="type.fold")+
    geom_text(aes(
      test.fold, cv.type, label=models_GP, tooltip=tooltip),
      data=model.wide,
      clickSelects="type.fold")+
    coord_equal()+
    scale_fill_gradient2("GP-random", high=scales::muted("red"), low=scales::muted("blue")),
  GPonly=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_point(aes(
      train.size, median),
      showSelected="type.fold",
      data=best.baseline)+
    geom_line(aes(
      train.size, median, group=selection.short),
      size=1,
      clickSelects="selection.short",
      showSelected="type.fold",
      data=GP.models.only)+
    geom_ribbon(aes(
      train.size, ymin=q25, ymax=q75, group=selection.short),
      alpha=0.55,
      showSelected="type.fold",
      clickSelects="selection.short",
      data=GP.models.only)+
    ylab("Percent correctly predicted intervals")+
    scale_x_log10(
      "Labeled sequences in train set")
)

## Compute improvement over random for all train.size:
best.each.size <- both.stats[
  random.or.GP=="random" | baseline.or.GP=="GP", data.table(
  median=max(median)
), by=.(cv.type, test.fold, type.fold, train.size, random.or.GP)]
best.each.wide <- dcast(
  best.each.size,
  cv.type + test.fold + type.fold + train.size ~ random.or.GP,
  value.var=c("median"))[is.finite(GP) & is.finite(random)]
best.each.wide[, improvement := GP-random]
min.sizes <- both.stats[, .(
  min.train.size=min(train.size)
), by=.(baseline.or.GP)]
xsc <- scale_x_continuous(
  "Labeled sequences in train set",
  breaks=c(min.sizes$min.train.size, seq(10, max.train.size, by=10)))
ysc <- ylab("Percent correctly predicted intervals")
show.vline <- data.table(train.size=10)
GP.one.size <- GP.models.only[
  show.vline, on=.(train.size)]
best.each.random <- data.table(
  baseline.or.GP="baseline",
  best.each.size[random.or.GP=="random"])
rand.one.size <- data.table(
  selection.short="best random baseline",
  best.each.random[show.vline, on=.(train.size)])
both.one.size <- rbind(
  rand.one.size,
  GP.one.size[, names(rand.one.size), with=FALSE]
)[order(type.fold, median)]
both.one.size[, rank := 1:.N, by=.(type.fold)]
both.one.size[, rank.y := min(GP.models.only$q25)+(rank-1)*5]
zoomout <- ggplot()+
  theme_bw()+
  theme_animint(width=450)+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_vline(aes(
    xintercept=train.size),
    data=show.vline,
    alpha=0.5)+
  geom_text(aes(
    2, 100, key=1, label=sprintf(
      "%d GP models for test fold %s",
      models_GP, type.fold)),
    hjust=0,
    data=model.wide,
    showSelected="type.fold")+
  geom_point(aes(
    train.size, median,
    key=selection.short),
    data=both.one.size[random.or.GP=="GP"],
    clickSelects="selection.short",
    alpha=0.7,
    showSelected=c("random.or.GP", "type.fold"))+
  geom_point(aes(
    train.size, median,
    key=1),
    data=both.one.size[random.or.GP=="random"],
    showSelected=c("random.or.GP", "type.fold"))+
  geom_line(aes(
    train.size, median,
    key=1,
    linetype=random.or.GP),
    data=best.each.random,
    showSelected=c("random.or.GP", "type.fold"))+
  geom_line(aes(
    train.size, median,
    key=selection.short,
    group=selection.short, linetype=random.or.GP),
    size=1,
    clickSelects="selection.short",
    showSelected=c("random.or.GP", "type.fold"),
    data=GP.models.only)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75,
    key=selection.short,
    group=selection.short),
    alpha=0.55,
    showSelected=c("random.or.GP", "type.fold"),
    clickSelects="selection.short",
    data=GP.models.only)+
  ysc+
  xsc
both.stats[, tooltip := sprintf(
  "model=%s selection=%s",
  model.dir, random.or.GP)]
GP.selection.only <- both.stats[random.or.GP=="GP"]
random.selection.only <- both.stats[random.or.GP=="random"]
model.colors <- c(
  GP="black",
  "#A6CEE3", L1reg_linear_all="#1F78B4", #lite dark blue
  "#B2DF8A", "#33A02C", # lite dark green
  unreg_linear_2="#6A3D9A", #dark violet
  unreg_linear_1="#E31A1C", #dark red
  "#FDBF6F", "#FF7F00", #lite dark orange
  baseline_0="#CAB2D6", #lite violet
  unsup_BIC_1="#FB9A99", #lite red
  "#FFFF99", "#B15928" #yellow brown
)
zoom.rect <- data.table(xmin=2, xmax=30, ymin=75, ymax=100)
coord.zoom <- zoom.rect[, coord_cartesian(
  xlim=c(xmin, xmax), ylim=c(ymin, ymax))]
viz <- animint(
  title="Gaussian Process models for sample selection",
  overview=ggplot()+
    ggtitle("Best improvement, select test fold")+
    theme_bw()+
    xsc+
    ylab("Accuracy percent improvement of best GP - best random baseline")+
    geom_vline(aes(
      xintercept=train.size),
      data=show.vline,
      alpha=0.5)+
    geom_hline(aes(
      yintercept=yintercept),
      data=data.table(yintercept=0),
      alpha=0.5)+
    geom_line(aes(
      train.size, improvement, group=type.fold),
      data=best.each.wide,
      alpha=0.8,
      size=5,
      clickSelects="type.fold"),
  GPzoomout=zoomout+
    ggtitle("Select GP and compare with other GP, all data")+
    geom_rect(aes(
      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
      data=zoom.rect,
      alpha=0.1)+
    geom_text(aes(
      train.size, rank.y,
      key=selection.short,
      label=sprintf(
        "median=%.1f%% at size=%d for %s",
        median, train.size, selection.short)),
      hjust=0,
      data=both.one.size[baseline.or.GP=="GP"],
      clickSelects="selection.short",
      showSelected=c("random.or.GP", "type.fold"))+
    geom_text(aes(
      train.size, rank.y,
      key=1,
      label=sprintf(
        "median=%.1f%% at size=%d for %s",
        median, train.size, selection.short)),
      hjust=0,
      data=both.one.size[baseline.or.GP!="GP"],
      showSelected=c("random.or.GP", "type.fold")),
  GPzoomin=zoomout+
    ggtitle("Select GP and compare with other GP, zoomed")+
    coord.zoom,
  GPdetails=ggplot()+
    ggtitle("Compare selected GP with baselines, zoomed")+
    theme_bw()+
    theme_animint(width=450)+
    xsc+
    ysc+
    coord.zoom+
    scale_color_manual(values=model.colors, breaks=names(model.colors))+
    scale_fill_manual(values=model.colors, breaks=names(model.colors))+
    geom_vline(aes(
      xintercept=train.size),
      data=show.vline,
      alpha=0.5)+
    geom_line(aes(
      train.size, median,
      key=1,
      linetype=random.or.GP),
      data=best.each.random,
      showSelected=c("random.or.GP", "type.fold"))+
    geom_line(aes(
      train.size, median, linetype=random.or.GP,
      color=model.dir,
      tooltip=tooltip,
      key=paste(model.dir, random.or.GP),
      group=paste(model.dir, random.or.GP)),
      size=1,
      showSelected=c("random.or.GP", "selection.short", "type.fold"),
      data=GP.selection.only)+
    geom_ribbon(aes(
      train.size, ymin=q25, ymax=q75,
      fill=model.dir,
      tooltip=tooltip,
      key=paste(model.dir, random.or.GP),
      group=paste(model.dir, random.or.GP)),
      alpha=0.25,
      showSelected=c("random.or.GP", "selection.short", "type.fold"),
      data=GP.selection.only)+
    geom_line(aes(
      train.size, median, linetype=random.or.GP,
      color=model.dir,
      tooltip=tooltip,
      key=paste(model.dir, random.or.GP),
      group=paste(model.dir, random.or.GP)),
      size=1,
      showSelected=c("random.or.GP", "type.fold"),
      data=random.selection.only)+
    geom_text(aes(
      2, 100, key=1, label=sprintf(
        "%s and baselines for test fold %s", selection.short, type.fold)),
      data=both.one.size[baseline.or.GP=="GP"],
      showSelected=c("type.fold", "selection.short"),
      hjust=0)+
    geom_ribbon(aes(
      train.size, ymin=q25, ymax=q75,
      fill=model.dir,
      tooltip=tooltip,
      key=paste(model.dir, random.or.GP),
      group=paste(model.dir, random.or.GP)),
      alpha=0.25,
      showSelected=c("random.or.GP", "type.fold"),
      data=random.selection.only),
  duration=list(type.fold=2000, selection.short=2000),
  first=list(
    model.dir=c("GP", "unreg_linear_2")
  )
)
animint2dir(viz, "figure-2019-08-14-animint")

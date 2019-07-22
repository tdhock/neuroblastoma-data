source("packages.R")

(baseline.csv.vec <- Sys.glob(
  "data/systematic/cv/*/testFolds/*/*/*/baseline.csv"))
baseline.list <- list()
for(csv.i in seq_along(baseline.csv.vec)){
  baseline.csv <- baseline.csv.vec[[csv.i]]
  seed.dir <- dirname(baseline.csv)
  seed <- as.integer(basename(seed.dir))
  selection.dir <- dirname(seed.dir)
  selection.type <- basename(selection.dir)
  split.dir <- dirname(selection.dir)
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  dt <- fread(baseline.csv)
  baseline.list[[csv.i]] <- data.table(
    selection.type,
    seed, test.fold, cv.type=basename(cv.type.dir), dt)
}
baseline <- do.call(rbind, baseline.list)

(gpaccuracy.csv.vec <- Sys.glob(
  "data/systematic/cv/*/testFolds/*/sampleSelection*/*/gpaccuracy.csv"))
gp.list <- list()
for(csv.i in seq_along(gpaccuracy.csv.vec)){
  gpaccuracy.csv <- gpaccuracy.csv.vec[[csv.i]]
  seed.dir <- dirname(gpaccuracy.csv)
  seed <- as.integer(basename(seed.dir))
  selection.dir <- dirname(seed.dir)
  selection.type <- basename(selection.dir)
  split.dir <- dirname(selection.dir)
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  ## read accuracy for baseline models on selected seq thru train data.
  models.dir <- file.path(seed.dir, "models")
  model.vec <- dir(models.dir)
  meta.dt <- data.table(
    selection.type, seed,
    test.fold,
    cv.type=basename(cv.type.dir))
  for(model in model.vec){
    accuracy.csv <- file.path(models.dir, model, "accuracy.csv")
    model.dt <- fread(accuracy.csv)
    gp.list[[paste(csv.i, model)]] <- model.dt[, data.table(
      meta.dt,
      model,
      train.size,
      percent.correct.intervals)]
  }
  ## read accuracy for GP model.
  gp.dt <- fread(gpaccuracy.csv)
  gp.list[[paste(csv.i, "GP")]] <- gp.dt[, data.table(
    meta.dt,
    model=sub("sampleSelection", "", selection.type),
    train.size=n,
    percent.correct.intervals=accuracy)]
}

## read accuracy for baseline models on random orderings.
(test.fold.dirs <- unique(dirname(dirname(dirname(gpaccuracy.csv.vec)))))
rand.acc.csv.vec <- Sys.glob(
  file.path(test.fold.dirs, "randomTrainOrderings", "*", "models", "*", "accuracy.csv"))
rand.list <- list()
for(rand.acc.csv in rand.acc.csv.vec){
  model.dir <- dirname(rand.acc.csv)
  models.dir <- dirname(model.dir)
  seed.dir <- dirname(models.dir)
  seed <- as.integer(basename(seed.dir))
  selection.dir <- dirname(seed.dir)
  selection.type <- basename(selection.dir)
  split.dir <- dirname(selection.dir)
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  rand.acc.dt <- fread(rand.acc.csv)
  rand.list[[rand.acc.csv]] <- rand.acc.dt[, data.table(
    selection.type, seed,
    test.fold,
    cv.type=basename(cv.type.dir),
    model=basename(model.dir),
    train.size,
    percent.correct.intervals)]
}

gp.base <- rbind(
  do.call(rbind, gp.list),
  do.call(rbind, rand.list))

both.stats <- gp.base[is.finite(percent.correct.intervals) & train.size<=60, list(
  mean=mean(percent.correct.intervals),
  sd=sd(percent.correct.intervals),
  median=median(percent.correct.intervals),
  q25=quantile(percent.correct.intervals, 0.25),
  q75=quantile(percent.correct.intervals, 0.75)
), by=list(selection.type, test.fold, cv.type, train.size, model)]
both.stats[, selection.disp := sub("TrainOrderings|sampleSelection", "", selection.type)]
both.stats[, Selection := ifelse(selection.disp=="random", "random", "GP")]
model.levs <- rev(sort(unique(both.stats$model)))
model.colors <- ifelse(grepl("GP", model.levs), "red", "black")
names(model.colors) <- model.levs
both.stats[, Model := factor(model, model.levs)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type + test.fold ~ selection.disp, labeller=label_both)+
  geom_line(aes(
    train.size, median, color=Model),
    size=1,
    data=both.stats)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75, fill=Model),
    alpha=0.5,
    data=both.stats)+
  ylab("Percent correctly predicted intervals")+
  coord_cartesian(xlim=c(2, 200), ylim=c(50, 100))+
  scale_color_manual(values=model.colors)+
  scale_fill_manual(values=model.colors)+
  scale_x_log10(
    "Labeled sequences in train set")
dl.method <- list(cex=0.7, "last.polygons")
(dl <- directlabels::direct.label(gg, dl.method))

png("figure-2019-07-22-all.png", 15, 6, units="in", res=100)
print(dl)
dev.off()

some.stats <- both.stats[grepl("GP", model) | model=="L1reg_linear_all"]
some.stats[, label := ifelse(model==selection.disp, model, paste0(model, "\n", selection.disp))]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ cv.type + test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, median, color=Model, linetype=Selection, group=paste(Model, selection.disp)),
    size=1,
    data=some.stats)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75, fill=Model, group=paste(Model, selection.disp)),
    alpha=0.5,
    data=some.stats)+
  ylab("Percent correctly predicted intervals")+
  coord_cartesian(xlim=c(2, 200), ylim=c(50, 100))+
  scale_color_manual(values=model.colors)+
  scale_fill_manual(values=model.colors)+
  scale_x_log10(
    "Labeled sequences in train set")+
  directlabels::geom_dl(aes(
    train.size, median, color=Model, label=label),
    data=some.stats, method=dl.method)
print(gg)
png("figure-2019-07-22-L1reg.png", 15, 6, units="in", res=100)
print(gg)
dev.off()

gg2 <- gg+
  facet_grid(. ~ cv.type + test.fold + selection.disp, labeller=label_both)
png("figure-2019-07-22.png", 15, 6, units="in", res=100)
print(gg2)
dev.off()

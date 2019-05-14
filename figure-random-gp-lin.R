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
  gp.dt <- fread(gpaccuracy.csv)
  gp.list[[csv.i]] <- gp.dt[, data.table(
    selection.type, seed,
    test.fold,
    cv.type=basename(cv.type.dir),
    train.size=n,
    model=sub("sampleSelection", "", selection.type),
    accuracy.percent=accuracy)]
}
gp <- do.call(rbind, gp.list)

both.dt <- rbind(baseline, gp)[!is.na(accuracy.percent)]

(counts.dt <- dcast(
  both.dt,
  test.fold + cv.type + model ~ selection.type,
  length,
  value.var="accuracy.percent"))
gp.selection <- counts.dt[0 < sampleSelectionGP_LIN | 0 < sampleSelectionGP_SE, .(
  test.fold, cv.type, model)]
both.show <- both.dt[gp.selection, on=list(test.fold, cv.type, model)]

both.stats <- both.show[, list(
  mean=mean(accuracy.percent),
  sd=sd(accuracy.percent)
), by=list(selection.type, test.fold, cv.type, train.size, model)][4 <= train.size & train.size <= 100]
both.stats[, selection.disp := sub("TrainOrderings|sampleSelection", "", selection.type)]
dl.method <- list(cex=0.7, "last.polygons")
area.under.mean <- both.stats[, {
  L <- approx(c(0, train.size), c(0, mean), 4:100)
  list(total=sum(L$y))
}, by=list(cv.type, test.fold, model, selection.disp)]
best.random.area <- area.under.mean[selection.disp=="random", {
  .SD[which.max(total)]
}, by=list(cv.type, test.fold)]

gp.model.stats <- both.stats[model %in% c("GP_LIN", "GP_SE")]
best.random.stats <- both.stats[model=="unreg_linear_2" & selection.disp=="random"]
some <- rbind(
  gp.model.stats[, data.table(Model=model, .SD)],
  best.random.stats[, data.table(Model="unreg_linear_2
random", .SD)])
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ cv.type + test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, mean, color=Model),
    data=some)+
  geom_ribbon(aes(
    train.size, ymin=mean-sd, ymax=mean+sd, fill=Model),
    alpha=0.5,
    data=some)+
  ylab("Percent correctly predicted intervals")+
  scale_x_continuous(
    "Labeled sequences in train set",
    limits=c(4, 120))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-random-gp-lin-mean.png", 15, 6, units="in", res=100)
print(dl)
dev.off()

gp.model.diff <- both.show[model %in% c("GP_LIN", "GP_SE"), {
  select.dt <- data.table(cv.type, test.fold, seed, model="unreg_linear_2")
  rand <- both.show[select.dt, on=list(cv.type, test.fold, seed, model)]
  L <- approx(rand$train.size, rand$accuracy.percent, train.size)
  rand.acc <- L$y
  data.table(
    train.size,
    model.acc=accuracy.percent,
    random.acc=rand.acc,
    diff.acc=accuracy.percent-rand.acc)
}, by=list(cv.type, test.fold, seed, model)]
gg <- ggplot()+
  theme_bw()+
  geom_hline(yintercept=0, color="grey")+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ cv.type + test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, diff.acc, color=model, group=paste(seed, model)),
    data=gp.model.diff)+
  ylab("Accuracy difference,
GP-Linear model with random selection")+
  scale_x_continuous(
    "Labeled sequences in train set",
    limits=c(4, 120))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-random-gp-lin-diff.png", 15, 6, units="in", res=100)
print(dl)
dev.off()

gp.diff.stats <- gp.model.diff[, list(
  mean=mean(diff.acc),
  sd=sd(diff.acc)
  ), by=list(cv.type, test.fold, model, train.size)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_hline(yintercept=0, color="grey")+
  facet_grid(. ~ cv.type + test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, mean, color=model),
    data=gp.diff.stats)+
  geom_ribbon(aes(
    train.size, ymin=mean-sd, ymax=mean+sd, fill=model),
    alpha=0.5,
    data=gp.diff.stats)+
  ylab("Accuracy difference,
 GP-Linear model with random selection")+
  scale_x_continuous(
    "Labeled sequences in train set",
    limits=c(4, 120))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-random-gp-lin-diff-mean.png", 15, 6, units="in", res=100)
print(dl)
dev.off()

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type + test.fold ~ model, labeller=label_both)+
  geom_line(aes(
    train.size, mean, color=selection.disp),
    data=both.stats)+
  geom_ribbon(aes(
    train.size, ymin=mean-sd, ymax=mean+sd, fill=selection.disp),
    alpha=0.5,
    data=both.stats)+
  ylab("Percent correctly predicted intervals")+
  scale_x_continuous(
    "Labeled sequences in train set",
    limits=c(4, 120))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-random-gp-lin.png", 15, 6, units="in", res=100)
print(dl)
dev.off()


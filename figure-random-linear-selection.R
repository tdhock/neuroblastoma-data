source("packages.R")

(baseline.csv.vec <- Sys.glob(
  "data/systematic/cv/*/testFolds/*/*/*/baseline.csv"))

result.list <- list()
for(csv.i in seq_along(baseline.csv.vec)){
  baseline.csv <- baseline.csv.vec[[csv.i]]
  seed.dir <- dirname(baseline.csv)
  seed <- as.integer(basename(seed.dir))
  selection.dir <- dirname(seed.dir)
  selection.type <- basename(selection.dir)
  selection.disp <- gsub("([A-Z])", "\n\\1", selection.type)
  split.dir <- dirname(selection.dir)
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  dt <- fread(baseline.csv)
  result.list[[csv.i]] <- data.table(
    seed, test.fold,
    selection.type, selection.disp,
    cv.type=basename(cv.type.dir),
    unique(dt))
}
result <- do.call(rbind, result.list)[!is.na(accuracy.percent)]
result[, Model := factor(model, unique(model))]

counts.dt <- dcast(
  result,
  test.fold + cv.type + train.size + Model ~ selection.type,
  length,
  value.var="accuracy.percent")
counts.both <- counts.dt[0 < sampleSelectionLinear & 0 < randomTrainOrderings]
res.both <- result[counts.both, on=list(test.fold, cv.type, train.size, Model)]

dl.method <- list(cex=0.7, "last.polygons")
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type + test.fold ~ Model, labeller=label_both)+
  geom_line(aes(
    train.size, accuracy.percent, color=selection.disp, group=paste(selection.disp,seed)),
    data=res.both)+
  ylab("Percent correctly predicted intervals")+
  scale_x_log10(
    "Labeled sequences in train set",
    breaks=c(
      range(res.both$train.size),
      20))+
  coord_cartesian(xlim=c(1, 1000), ylim=c(80, 100))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-random-linear-selection.png", 15, 6, units="in", res=100)
print(dl)
dev.off()


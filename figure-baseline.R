source("packages.R")

(baseline.csv.vec <- Sys.glob(
  "data/systematic/cv/*/testFolds/*/randomTrainOrderings/*/baseline.csv"))

result.list <- list()
for(csv.i in seq_along(baseline.csv.vec)){
  baseline.csv <- baseline.csv.vec[[csv.i]]
  seed.dir <- dirname(baseline.csv)
  seed <- as.integer(basename(seed.dir))
  split.dir <- dirname(dirname(seed.dir))
  test.fold <- as.integer(basename(split.dir))
  cv.type.dir <- dirname(dirname(split.dir))
  dt <- fread(baseline.csv)
  result.list[[csv.i]] <- data.table(
    seed, test.fold, cv.type=basename(cv.type.dir), dt)
}
result <- do.call(rbind, result.list)

most.train <- result[, {
  .SD[train.size==max(train.size), list(
    min.accuracy=min(accuracy.percent),
    max.accuracy=max(accuracy.percent)
  )]
}, by=list(cv.type, test.fold, model)]
(bad <- most.train[min.accuracy != max.accuracy])
if(nrow(bad)){
  print(result[bad, on=list(cv.type, test.fold, model)][train.size==max(train.size)])
  stop("some accuracy numbers not equal")
}

if(FALSE){
  unlink(Sys.glob(
    "data/systematic/cv/profileSize/testFolds/2/randomTrainOrderings/*/baseline.csv"))
}

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type ~ test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, accuracy.percent, color=model, group=paste(model,seed)),
    data=result)+
  scale_x_log10(limits=c(NA, 30000), breaks=c(
    range(result$train.size),
    10^(1:3)))+
  coord_cartesian(ylim=c(80, 100))
(dl <- directlabels::direct.label(gg, "last.polygons"))

png("figure-baseline-lines.png", 8, 6, units="in", res=100)
print(dl)
dev.off()

stats.dt <- result[, list(
  mean=mean(accuracy.percent),
  sd=sd(accuracy.percent)
  ), by=list(cv.type, test.fold, train.size, model)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type ~ test.fold, labeller=label_both)+
  geom_line(aes(
    train.size, mean, color=model, group=model),
    data=stats.dt)+
  geom_ribbon(aes(
    train.size, ymin=mean-sd, ymax=mean+sd, fill=model, group=model),
    alpha=0.5,
    data=stats.dt)+
  ylab("Percent correctly predicted intervals")+
  scale_x_log10(limits=c(NA, 30000), breaks=c(
    range(result$train.size),
    10^(1:3)))+
  coord_cartesian(ylim=c(80, 100))
(dl <- directlabels::direct.label(gg, "last.polygons"))

png("figure-baseline.png", 8, 6, units="in", res=100)
print(dl)
dev.off()

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

ymin <- 75
png.width <- 15
xmax <- 5e4
max.dt <- result[, list(max.train=max(train.size)), by=list(cv.type, test.fold)]
max.color <- "grey"
break.vec <- c(
  min(result$train.size),
  10^(1:3))
dl.method <- list(cex=0.7, "last.polygons")
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type ~ test.fold, labeller=label_both)+
  geom_segment(aes(
    max.train, Inf,
    xend=max.train, yend=ymin),
    color=max.color,
    data=max.dt)+
  geom_text(aes(
    max.train, ymin, label=paste0(" ", max.train)),
    color=max.color,
    hjust=0,
    vjust=0,
    data=max.dt)+
  geom_line(aes(
    train.size, accuracy.percent, color=model, group=paste(model,seed)),
    data=result)+
  ylab("Percent correctly predicted intervals")+
  scale_x_log10(
    "Labeled sequences in train set",
    limits=c(NA, xmax), breaks=break.vec)+
  coord_cartesian(ylim=c(ymin, 100))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-baseline-lines.png", png.width, 6, units="in", res=100)
print(dl)
dev.off()

stats.dt <- result[is.finite(accuracy.percent), list(
  mean=mean(accuracy.percent),
  median=median(accuracy.percent),
  q25=quantile(accuracy.percent, 0.25),
  q75=quantile(accuracy.percent, 0.75),
  sd=sd(accuracy.percent)
  ), by=list(cv.type, test.fold, train.size, model)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(cv.type ~ test.fold, labeller=label_both)+
  geom_segment(aes(
    max.train, Inf,
    xend=max.train, yend=ymin),
    color=max.color,
    data=max.dt)+
  geom_text(aes(
    max.train, ymin, label=paste0(" ", max.train)),
    color=max.color,
    hjust=0,
    vjust=0,
    data=max.dt)+
  geom_line(aes(
    train.size, median, color=model, group=model),
    data=stats.dt)+
  geom_ribbon(aes(
    train.size, ymin=q25, ymax=q75, fill=model, group=model),
    alpha=0.5,
    data=stats.dt)+
  ## geom_line(aes(
  ##   train.size, mean, color=model, group=model),
  ##   data=stats.dt)+
  ## geom_ribbon(aes(
  ##   train.size, ymin=mean-sd, ymax=mean+sd, fill=model, group=model),
  ##   alpha=0.5,
  ##   data=stats.dt)+
  ylab("Percent correctly predicted intervals")+
  scale_x_log10(
    "Labeled sequences in train set",
    limits=c(NA, xmax), breaks=break.vec)+
  coord_cartesian(ylim=c(ymin, 100))
(dl <- directlabels::direct.label(gg, dl.method))
png("figure-baseline.png", png.width, 6, units="in", res=100)
print(dl)
dev.off()

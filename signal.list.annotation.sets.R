if(!file.exists("signal.list.annotation.sets.RData")){
  download.file("http://members.cbio.ensmp.fr/~thocking/neuroblastoma/signal.list.annotation.sets.RData", "signal.list.annotation.sets.RData")
}
(objs <- load("signal.list.annotation.sets.RData"))
library(data.table)
size.vec <- sapply(signal.list, nrow)
seq.labels <- do.call(rbind, lapply(names(annotation.sets), function(label.set){
  label.df <- annotation.sets[[label.set]]
  data.table(label.df)[, data.table(
    label.set,
    labels=.N
  ), by=.(pid.chr=paste0(profile.id, ".", chromosome), annotation)]
}))
set.sizes <- seq.labels[, {
  u.ids <- unique(pid.chr)
  u.sizes <- size.vec[paste(u.ids)]
  as.list(quantile(u.sizes, seq(0, 1, l=3)))
}, by=.(label.set)]
set.labels <- dcast(
  seq.labels,
  label.set ~ annotation,
  value.var="labels",
  fun.aggregate=sum)
set.sizes[set.labels, on=.(label.set)][order(`50%`)]

(labels.wide <- dcast(
  seq.labels,
  label.set + pid.chr ~ annotation,
  value.var="labels",
  fun.aggregate=sum)[order(`1breakpoint`)])
labels.wide[, n.data := size.vec[pid.chr] ]
print(labels.wide)

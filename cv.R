source("packages.R")

labels.xz.vec <- Sys.glob("data/*/labels.csv.xz")
for(set.i in seq_along(labels.xz.vec)){
  labels.xz <- labels.xz.vec[[set.i]]
  labels.dt <- fread(cmd=paste("xzcat", labels.xz))
  namedCapture::df_match_variable(
    labels.dt,
    sequenceID=list(
      profileID="[0-9]+",
      "_",
      chrom="chr.*"))
}

source("packages.R")

data.dir.vec <- Sys.glob("data/*")
file.exists(file.path(data.dir.vec, "evaluation.csv.xz"))

data.dir <- "data/systematic"


other.eval <- fread(cmd="xzcat data/detailed/evaluation.csv.xz")
names(other.eval)

errors.csv.xz <- file.path(data.dir, "errors.csv.xz")
errors.dt <- fread(cmd=paste("xzcat", errors.csv.xz))
eval.dt <- errors.dt[, names(other.eval), with=FALSE]

eval.csv <- file.path(data.dir, "evaluation.csv")
fwrite(eval.dt, eval.csv)
eval.csv.xz <- paste0(eval.csv, ".xz")
unlink(eval.csv.xz)
system(paste("xz", eval.csv))

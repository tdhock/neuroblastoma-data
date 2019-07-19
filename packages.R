### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,...){
  local.lib <- file.path(getwd(), "library")
  dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
  .libPaths(c(local.lib, .libPaths()))
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
options(repos=c(
          "http://www.bioconductor.org/packages/release/bioc",
          "http://r-forge.r-project.org",
          "http://cloud.r-project.org",
          "http://cran.r-project.org"))
works_with_R(
  "3.6.0",
  ggplot2="3.1.1",
  data.table="1.12.2",
  directlabels="2018.5.22",
  partykit="1.2.4", doParallel="1.0.14", # for mmit.
  future="1.13.0", future.apply="1.3.0",
  "anujkhare/iregnet@4d77f047c3a00a5524e1cbe140226417e3aedd92",
  "tdhock/penaltyLearning@fc2833c47d18f20e99648714be47829cf41cc084")
library(survival)
if(!requireNamespace("bams")){
  if(!file.exists("bams_1.6.tar.gz")){
    u <- "https://cran.r-project.org/src/contrib/Archive/bams/bams_1.6.tar.gz"
    download.file(u, "bams_1.6.tar.gz")
  }
  install.packages("bams_1.6.tar.gz", type="source", repos=NULL)
}
if(FALSE){
  requireGitHub::requireGitHub_package(
    "aldro61",
    "mmit/Rpackage",
    "b9d77ec604185725b51f13ec3922a6ba59bda145",
    "mmit")
}else{
  if(!require(mmit)){
    system("git clone https://github.com/aldro61/mmit.git")
    install.packages("~/R/mmit/Rpackage", type="source", repo=NULL)
  }
}

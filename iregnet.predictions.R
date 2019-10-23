source("packages.R")

folds.csv.vec <- Sys.glob("data/*/cv/*/folds.csv")

fit.list <- list()
scale.i.list <- list( estimated = list( init= "NA", estimate = TRUE), fixed = list( init= 1, estimate = FALSE) )
model.list <- list( "guassian", "logistic", "extreme_value")



table(n.pred.vec)
#take in test fold
#one fold
OneFold <- function(order.csv){
  
}


for( model.type in model.i.list){
  for( scale.type in scale.i.list){
   fit.list <- c( fit.list , cv.iregnet(X.train, Y.train , family = model.type, 
                  scale_init= scale.type$init ,estimate_scale = scale.type$estimate))  
  }
}
    
  
  ## fwrite(result, iregnet.csv)

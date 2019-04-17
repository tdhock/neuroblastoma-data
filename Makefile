data/systematic/cv/sequenceID/folds.csv: cv.R
	R --vanilla < $<
data/systematic/inputs.csv.xz: neuroblastoma.R
	R --vanilla < $<

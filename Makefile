data/systematic/cv/sequenceID/testFolds/6/randomTrainOrderings/5/baseline.csv: baseline.R
	R --vanilla < $<
figure-baseline.png: figure-baseline.R
	R --vanilla < $<
data/systematic/cv/sequenceID/testFolds/6/randomTrainOrderings/5/order.csv: randomOrderings.R
	R --vanilla < $<
data/systematic/cv/sequenceID/folds.csv: cv.R
	R --vanilla < $<
data/systematic/inputs.csv.xz: neuroblastoma.R
	R --vanilla < $<

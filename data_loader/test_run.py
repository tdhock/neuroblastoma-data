import dataproFunc as dpF
import sys, os
data_dir = sys.argv[1]
kwargs = {
"foldIDFile" : data_dir + "/cv/profileID/folds.csv",
"signalFile" : data_dir + '/profiles.csv.xz',
"outputFile" : data_dir + '/outputs.csv.xz',
"boxLen" : 100
}
for k, v in kwargs.iteritems():
    if "File" in k:
        if not os.path.exists(v):
            raise ValueError(v + " file does not exist")
mat2 = dpF.go_process(**kwargs)
mat.to_csv(data_dir + '/processedmat.csv')

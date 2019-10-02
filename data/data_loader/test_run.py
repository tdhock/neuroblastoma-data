import dataproFunc as dpF

mat2 = dpF.go_process(signalFile='./detailed/profiles.csv.xz',
                     boxLen=100, outputFile='./detailed/outputs.csv.xz',
                     foldIDFile='./detailed/folds.csv')
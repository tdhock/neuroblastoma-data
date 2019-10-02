import numpy as np
import pandas as pd


# the function takes the signalfile, the desired length,
# outputfile and foldIDfile as input 
# return: the processed dataset that contains the name, boxed signal,
# the foldID and the output
def go_process(signalFile:str, boxLen:int, outputFile:str, foldIDFile:str):
    # input the signal file
    data = pd.read_csv(signalFile)
    # get the num of sequenceID
    seq_list_name = data['sequenceID'].unique()

    # define the column range based on the data we have
    boxes_range = np.linspace(data['position'].min(),
     data['position'].max(), boxLen)
     # define the name of each column and also and index of each row
    columns = ['x' + str(iter) for iter in range(0, boxLen)]
    index = np.arange(0, len(seq_list_name), 1)

    # define a new frame
    frame = pd.DataFrame(columns=columns, index=index)

    # divide the the signal to each box
    # this is the step takes the most of time
    # be patient here, it is already faster than old __version
    for idx, row in frame.iterrows():
        name = seq_list_name[idx]
        seqID_data = data[data['sequenceID']==name]
        for indexx, ran_val in enumerate(boxes_range[:-1]):
            setset = np.array(seqID_data['position'] >= ran_val) \
                & np.array(seqID_data['position'] < boxes_range[indexx+1])
            if setset.any() == True:
                row[indexx] = seqID_data.loc[setset,'signal'].mean()
        
    frame2 = frame.copy()
    frame2.columns = {a for a in range(0, boxLen)}
    rever_u_d = frame2.iloc[:,::-1,]
    start_bin_idx = (frame2.notnull()).idxmax(axis=1, skipna=True)
    end_bin_idx = (rever_u_d.notnull()).idxmax(axis=1, skipna=True)
    bin_size = end_bin_idx - start_bin_idx + 1
    data_describe = pd.DataFrame({'start_bin_idx':pd.Series(start_bin_idx),\
        'end_bin_idx':pd.Series(end_bin_idx), 'bin_size':pd.Series(bin_size)})
    
    # This part of code can be modified to get the mean.
    for idx, rows in frame2.iterrows():
        temp = rows[data_describe.loc[idx,'start_bin_idx']:data_describe.loc[idx,'end_bin_idx']]
        temp.fillna(method='ffill', inplace=True)
        if temp.isnull().any():
            print(idx)
    

    # start process the label from here
    label = pd.read_csv(outputFile)
    label['label'] = np.nan
    op_idx = (label['min.log.lambda'] == -np.inf) & (label['max.log.lambda'] != np.inf)
    label.loc[op_idx, 'label'] = label.loc[op_idx, 'max.log.lambda'] - 1 
    op_idx = (label['min.log.lambda'] != -np.inf) & (label['max.log.lambda'] == np.inf)
    label.loc[op_idx, 'label'] = label.loc[op_idx, 'min.log.lambda'] + 1
    op_idx = (label['min.log.lambda'] != -np.inf) & (label['max.log.lambda'] != np.inf)

    frame2['sequenceID'] = pd.Series(seq_list_name)
    frame2['label'] = np.nan
    label_sq_sort = label.sort_values(by=['sequenceID'])
    frame2 = frame2.sort_values(by=['sequenceID']).reset_index()
    frame2['label'] = label_sq_sort['label'].values

    # assign foldID to the dataset
    frame2['foldID'] = np.nan
    foldID = pd.read_csv(foldIDFile)
    foldID = foldID.sort_values(by=['sequenceID']).reset_index()
    frame2['foldID'] = foldID['fold'].values
    return(frame2)



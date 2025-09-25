'''
Contains helper functions used to produce figures.
'''

import numpy as np



def ave_interval_1D(data, start_index, end_index):
    '''
    Computes the average population/payoff in each patch of a 1D model over a time interval defined by start_index and end_index
    Borrowed from piegy.figure_tools: exactly the same as the ave_interval_1D over there

    Inputs:
    - data: a 3D array with size 1 x M x max_record, expect model.U, model.V, model.Hpi, or model.Dpi
    - start_index and end_index: the time interval over which to compute the average

    Returns:
    - data_ave: 1D array, data_ave[i] is the average value of data[0][:][start_index : end_index]
    '''
    
    N = len(data)
    M = len(data[0])

    if start_index == end_index:
        start_index = end_index - 1
        
    data_ave = np.zeros(N * M)
    
    for i in range(N):
        for j in range(M):
            for k in range(start_index, end_index):
                data_ave[i * M + j] += data[i][j][k]
            data_ave[i * M + j] /= (end_index - start_index)
    
    return data_ave


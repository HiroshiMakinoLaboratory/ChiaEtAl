import FORCE
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat, savemat

if __name__ == '__main__':
    # additional batches of external coupling units
    addBatch = ['PPC', 'vS1']
    # train 35 models each
    for rep in range(35):
        for ab in addBatch:
            data_path = 'data/ALM_units.mat'
            addData_path = f'data/coupling_units_from_{ab}.mat'
            
            if os.path.exists('logs') is False:
                os.mkdir('logs')
            
            save_path = f'logs/ALMw{ab}_{rep}'
            
            if os.path.exists(save_path) is False:
                p = FORCE.create_hyperparameters(data_path, addData_path, save_path, rep)
                net = FORCE.RNN(p, save_path)
                net.train()

            if os.path.exists(save_path):
                p = FORCE.create_hyperparameters(data_path, addData_path, save_path, rep)
                net = FORCE.RNN(p, save_path)
                net.test(saveplots=True)
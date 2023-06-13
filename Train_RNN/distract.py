import FORCE
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat, savemat
import argparse
import pickle as pkl
import os

data_path = 'data/ALM_units.mat'
logs_path = 'logs'

temp = os.listdir(logs_path)
print(temp)

for drname in temp:
    try:
        rep = int(drname[-2:])
    except:
        rep = int(drname[-1])
    save_path = os.path.join(logs_path, drname)
    print(save_path)
    rnn_dr = os.path.join(save_path, 'RNN.pkl')
    
    add_index = save_path.find('ALM') + len('ALMw')
    add = save_path[add_index : add_index + 3]
    addData_path = f'data/coupling_units_from_{add}.mat'

    matfile_path = os.path.join(save_path, 'distracted')
    dist = 0.2
    dist_std = 0.02
    stim_std = 0.1
    abl_prop = 0.4
    abl_type = 'external'
    dist_type = 'early'
    dist_dur = 0.5
    
    if os.path.exists(matfile_path) is False:
        os.mkdir(matfile_path)
        
    p = FORCE.create_hyperparameters(data_path, addData_path, save_path, rep)
    net = FORCE.RNN(p, save_path)
    net.distract(distractor_amp=dist,
                 distractor_std=dist_std,
                 stim_std=stim_std,
                 ablation_proportion=abl_prop,
                 number_of_ablation_trials=1,
                 number_of_distraction_trials=100,
                 savedir=matfile_path,
                 ablation_type=abl_type,
                 distractor_type=dist_type,
                 distractor_duration=dist_dur)

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import pickle as pkl
import os
from scipy.io import loadmat, savemat
from math import ceil

def create_hyperparameters(data_path, addData_path, save_path, repetition_idx):
	"""
        Create/Modify hyperparameters settings for RNN.

    Args:
		data_path: 		Directory of mat file containing ALM activities data
		addData_path:	Directory of mat file containing coupling units data
		save_path: 		Directory to save RNN.pkl file

    Returns:
        p:    Dictionary containing the hyperparameters for RNN.

    """

	pklfile = os.path.join(save_path, 'RNN.pkl')
	if os.path.exists(pklfile):
		with open(pklfile, 'rb') as f:
			p = pkl.load(f)
		print('Hyperparameters loaded.')

	else:
		# Setup hyperparameters for RNN
		p = {
			'dt': 0.001,					# integration learning time step
			'tau': 0.01,					# time constant of neurons
			'g': 1.2,						# spectial radius, where g in range [1.2, 1.5]
			'p': 1.0,						# sparsity of the network, where p in range [0,1], where 1 is fully connected.
			'alpha': 0.01,					# coefficient of P0
			'training_trials': 500,			# number of training trials for RNN
			'learning_rate': 0.05,			# learning rate for weights update
			'repetition_index': repetition_idx,
		# Setup hyperparameters for inputs
			'stim_amplitude_R': 1.0,		# stimulus amplitude for R-licking trials
			'stim_amplitude_L': 0.0,		# stimulus amplitude for L-licking trials
			'additional_batch': 128,		# number of additional neurons
			}
		
		data = loadmat(data_path)
		L, R = data['L'], data['R']

		# Clipping and max normalization
		L, R = normalize(L), normalize(R)
		# Apply selection criteria
		concatenated_LR = np.hstack((L, R))
		normalized_LR = normalize(concatenated_LR)
		selection_criteria = normalized_LR > np.std(normalized_LR)
		passed = np.any(selection_criteria, axis=1)
		L, R = L[passed], R[passed]

		addData = loadmat(addData_path)
		addL, addR = addData['L'], addData['R']
		rng = npr.default_rng(p['repetition_index'])
		addBatch = rng.choice(range(len(addL)), p['additional_batch'], replace=True)
		addL, addR = addL[addBatch], addR[addBatch]
		offset = 0.01
		normalized_offset = offset / 5
		addL = np.clip(addL, 0.0 + normalized_offset, 1.0 - normalized_offset)
		addR = np.clip(addR, 0.0 + normalized_offset, 1.0 - normalized_offset)
		# Upsample data from sampling rate of 9.35Hz to 93.5Hz
		addL, addR = interpolating(addL, 9), interpolating(addR, 9)
		addL, addR = truncate(addL), truncate(addR)
		add_inv_sig_L, add_inv_sig_R = inv_sigmoid(addL), inv_sigmoid(addR)
		# 38/93.5 = 406ms Smoothing window
		addL, addR = moving_average(addL, 38), moving_average(addR, 38)
		add_inv_sig_L, add_inv_sig_R = moving_average(add_inv_sig_L, 38), moving_average(add_inv_sig_R, 38)

		# Upsample data from sampling rate of 9.35Hz to 93.5Hz
		L, R = interpolating(L, 9), interpolating(R, 9)
		# Only consider 0.5s of pre-sample epoch, 1.0s of sample epoch and 2.0s of delay epoch
		L, R = truncate(L), truncate(R)
		# Apply inverse sigmoid to obtain target functions
		inv_sig_L, inv_sig_R = inv_sigmoid(L), inv_sigmoid(R)
		# 38/93.5 = 406ms Smoothing window
		L, R = moving_average(L, 38), moving_average(R, 38)
		inv_sig_L, inv_sig_R = moving_average(inv_sig_L, 38), moving_average(inv_sig_R, 38)

		p.update({'L_activities': np.vstack((L, addL)), 'R_activities': np.vstack((R, addR))})
		p.update({'target_fn_L' : np.vstack((inv_sig_L, add_inv_sig_L)), 'target_fn_R': np.vstack((inv_sig_R, add_inv_sig_R))})
		p.update({'network_size': L.shape[0] + addL.shape[0]})

		if os.path.exists(save_path) is False:
			os.mkdir(save_path)
		with open(pklfile, 'wb') as f:
			pkl.dump(p, f)
		print('Number of cells: {}'.format(p['network_size']))
		print('Hyperparameters saved.')
	return p

def interpolating(activity, n_interp_points):
	'''
		Upsample the data
	'''
	interpolated = np.repeat(activity, n_interp_points+1).reshape(activity.shape[0], activity.shape[1]*(n_interp_points+1))
	interpolated = interpolated.astype('float')
	diff = (activity[:,1:] - activity[:,:-1])/(n_interp_points+1)
	for i in range(diff.shape[1]):
		for n in range(n_interp_points):
			interpolated[:,i*(n_interp_points+1) + n+1] += diff[:,i] * (n+1)
	for _ in range(n_interp_points):
		interpolated = np.delete(interpolated, -1, 1)
	return interpolated

def truncate(tseries):
	'''
		Crop the time series into 0.5s pre-sample epoch, 1.0s sample epoch and 2.0s delay epoch 
	'''
	start = round(0.5 * 93.5)
	end = round(4.0 * 93.5)
	if len(tseries.shape) > 1:
		tseries = tseries[:, start:end]
	else:
		tseries = tseries[start:end]
	return tseries

def normalize(tseries, min_value=0.0, max_value=5.0):
	'''
		Clip and max normalize the data
	'''
	# To avoid divide by zero during inverse sigmoid
	offset = 0.01
	## Threshold
	tseries = np.clip(tseries, min_value + offset, max_value - offset)
	## Max Normalization
	tseries = tseries / max_value
	return tseries

def moving_average(activity, window_size = 5):
	'''
		Data Smoothing
	'''
	left_tail = window_size//2
	right_tail = window_size-left_tail-1
	smoothed = np.zeros_like(activity)
	timelen = activity.shape[1]
	for i in range(timelen):
		if i + right_tail < window_size:
			smoothed[:,i] = np.mean(activity[:, :right_tail+i+1], axis=1)
		elif timelen - (i-left_tail) < window_size:
			smoothed[:,i] = np.mean(activity[:, -(timelen-(i-left_tail)):], axis=1)
		else:
			smoothed[:,i] = np.mean(activity[:, i-left_tail:i+right_tail+1], axis=1)
	return smoothed

def inv_sigmoid(y, beta=0.8, theta=3.0):
	'''
		Inverse sigmoidal function
	'''
	return theta + (1/beta) * np.log(y/(1-y))

# Preparing inputs for RNN
def prepare_ext_inputs(activity, target_fn, stim_amp, dist_amp=None, dist_type='early', dist_duration=0.5):
	'''
		Generate the time series for all the external inputs, including stimulus, cue and distractor
	'''
	stim = np.zeros(activity.shape[1])
	
	## Triangular Stim
	n_bins = round(1.0 * 93.5)
	half_t = ceil(n_bins/2)
	up_slope = (np.arange(half_t)+1)/half_t
	down_slope = np.flip(up_slope)
	if (n_bins - half_t) == half_t:
		triangle = np.concatenate((up_slope, down_slope)) * stim_amp
	else:
		triangle = np.concatenate((up_slope, down_slope[1:])) * stim_amp
	start = round(0.5 * 93.5)
	end = start + n_bins
	stim[start:end] = triangle
	
	# Triangular distractor
	if dist_amp is not None:
		n_bins = round(dist_duration * 93.5)
		half_t = ceil(n_bins/2)
		up_slope = (np.arange(half_t)+1)/half_t
		down_slope = np.flip(up_slope)
		if (n_bins - half_t) == half_t:
			triangle = np.concatenate((up_slope, down_slope)) * dist_amp
		else:
			triangle = np.concatenate((up_slope, down_slope[1:])) * dist_amp

		if dist_type == 'early':
			# Start at 1.75s (at 0.25s of the delay ep)
			start = round(1.75 * 93.5)
		elif dist_type == 'mid':
			# Start at 2.0s (at 0.5s of the delay ep)
			start = round(2.0 * 93.5)
		end = start + n_bins
		stim[start:end] = triangle
	
	# Delimiter Cue
	cue = np.zeros(activity.shape[1])

	n_bins = round(0.1 * 93.5)
	half_t = ceil(n_bins/2)
	up_slope = (np.arange(half_t)+1)/half_t
	down_slope = np.flip(up_slope)
	if (n_bins - half_t) == half_t:
		triangle = np.concatenate((up_slope, down_slope))
	else:
		triangle = np.concatenate((up_slope, down_slope[1:]))
	start = round(1.4 * 93.5)
	end = start + n_bins
	cue[start:end] = triangle
		
	stim = stim.repeat(activity.shape[0]).reshape(activity.shape[1], activity.shape[0]).T
	cue = cue.repeat(activity.shape[0]).reshape(activity.shape[1], activity.shape[0]).T

	inputs = {'orig': activity, 'target': target_fn, 'stim': stim, 'cue': cue}

	return inputs

def save_readouts(readouts, savedir, filename=None):
	'''
		Convert the readouts into mat file
	'''
	if filename is None:
		filename = os.path.join(savedir, 'readouts.mat')
	else:
		filename = os.path.join(savedir, filename)
	savemat(filename, readouts)
	print('Converted to Mat file.')
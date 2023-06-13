import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle as pkl
import os
from .utils import *

class RNN:
	"""
        RNN via FORCE learning.

	API:
		.reset(reset_activity=True, reset_weights=True)
		.save_weights()
		.load_weights()
		.activation(x, beta=0.8, theta=3.0)
		.run(inputs, t)
		.train(trial_types=['L', 'R'])
		.test(trial_types=['L', 'R'], saveplots=False)
		.plot_readouts(trial_types, ground_truth, readouts, save_folder_name):
		.distract(trial_types=['L', 'R', 'L_pos'],
				  distractor_amp=1.0,
				  distractor_std=0.1,
				  stim_std=0.1, 
				  ablation_proportion=0.0, 
				  number_of_ablation_trials=1,
				  number_of_distraction_trials=1, 
				  savedir=None,
				  ablation_type='external',
				  distractor_type='early'):
    """

	def __init__(self, hyperparameters, save_path):
		self.p = hyperparameters 	# assign hyperparameters of RNN as p
		self.display_idx = 0		# display the PSTH of neuron
		self.rep_idx = self.p['repetition_index']
		self.reset()				# initialize the parameters of RNN
		self.logdir = save_path
		self.pklfile = os.path.join(save_path, 'RNN.pkl')
		

	def reset(self, reset_activity=True, reset_weights=True):
		'''
			Re-initialize the paramaters of RNN.
		'''
		rnn_size = self.p['network_size']
		trg_L, trg_R = self.p['target_fn_L'], self.p['target_fn_R']
		
		if reset_activity:
			self.activity = np.mean(np.vstack((np.mean(trg_L[:,:10], 1), np.mean(trg_R[:,:10], 1))).T, 1).reshape(-1,1)
		# initialize all weights in RNN
		if reset_weights:
			rng = npr.default_rng(self.rep_idx)
			self.weights = {
							'rec': rng.normal(0, self.p['g']/np.sqrt(rnn_size), [rnn_size, rnn_size]),
							'stim': rng.normal(0, 1, [rnn_size, 1]),
							'cue': rng.normal(0, 0.1, [rnn_size, 1])
							}

	def save_weights(self):
		'''
			Save trained weights.
		'''
		if 'weights' in self.p.keys():
			self.p['weights'] = self.weights
			print('Weights Saved.')
		else:
			self.p.update({'weights': self.weights})
			print('Activity Saved.')
		
		if 'activity' in self.p.keys():
			self.p['activity'] = self.activity
		else:
			self.p.update({'activity': self.activity})

		with open(self.pklfile, 'wb') as f:
			pkl.dump(self.p, f)
	
	def load_weights(self):
		'''
			Load trained weights.
		'''
		with open(self.pklfile, 'rb') as f:
			self.p = pkl.load(f)

		if 'weights' in self.p.keys():
			self.weights = self.p['weights']
			print('Weights Loaded.')
		
		if 'activity' in self.p.keys():
			self.activity = self.p['activity']
			print('Activity Loaded.')

	def activation(self, x, beta=0.8, theta=3.0):
		'''
			Sigmoid function in numpy format
		'''
		r = 1/(1 + np.exp(-1 * beta * (x - theta)))
		return r

	def run(self, inputs, t):
		'''
			Calculate the activities and readouts of RNN at given time in a sequential manner.
		'''
		rnn_size = self.p['network_size']

		x = self.activity

		rng = npr.default_rng()
		# External inputs
		h = (self.weights['stim'] * inputs['stim'][:,t].reshape(-1,1) +
			self.weights['cue'] * inputs['cue'][:,t].reshape(-1,1) + 
			rng.normal(0, 0.15, [rnn_size, 1]))

		# Estimated total current inputs
		z = (np.dot(self.weights['rec'], self.activation(x)) + h)

		# Change in activity
		dx = (self.p['dt'] / self.p['tau']) * (-x + z)
		
		# New activity
		self.activity = (x + dx).reshape(rnn_size, 1)

		return self.activity, z

	def train(self, trial_types=['L', 'R']):
		'''
			Training RNN.
		'''
		# Reset the parameters of RNN
		self.reset()
		# Initialize save dictionary for losses
		## MSE_z: mean(sq(f - z)), where f is target function
		## MSE_r: mean(sq(n - r)), where n is normalized spike rate
		save_dict = {'MSE_z': [], 'MSE_r': []}
		
		# instantiate parameters
		rnn_size = self.p['network_size']
		P0 = self.p['alpha'] * np.eye(rnn_size)
		J = self.weights['rec']

		print('Training...')
		pbar = tqdm(range(self.p['training_trials']))
		rng = npr.default_rng(self.rep_idx)
		loss = []
		for ep in pbar:
			# At the start of every episode, reset the activity x
			self.reset(reset_activity=True, reset_weights=False)
			print(f'Epoch {ep+1}')
			loss_MSE_z = []
			P = P0

			# Alternate between L and R
			if ep % 2 == 0:
				trial_type = 'R'
			else:
				trial_type = 'L'
			print(trial_type)
			stim_amp = self.p[f'stim_amplitude_{trial_type}']
			stim_std = 0.1

			if stim_amp == 0:
				stochastic_stim_amp = 0.0
			else:
				# Sample stim amplitude from normal distribution
				stochastic_stim_amp = rng.normal(stim_amp, stim_std)
			
			act = self.p[f'{trial_type}_activities']
			tgt_fn = self.p[f'target_fn_{trial_type}']
			print(f'Stim: {stochastic_stim_amp}')
			inputs = prepare_ext_inputs(act, tgt_fn, stochastic_stim_amp)
			
			# total length of time for a trial
			timelen = inputs['target'].shape[1]
			r_t = np.zeros((rnn_size, timelen))
			z_t = np.zeros((rnn_size, timelen))

			learning_interval = 1

			for t in range(timelen):
				x, z = self.run(inputs, t)
				r = self.activation(x)
				r_t[:,t] = r.flatten()
			
				f = inputs['target'][:,t].reshape(-1,1)

				# Update the weights at intervals
				if t % learning_interval == 0:
					Pr = np.dot(P, r)
					c = 1/(1+r.T.dot(Pr))
					
					err = z - f
					LR = self.p['learning_rate']
					delta_J = c * np.outer(err.flatten(), Pr.flatten())
					J = J - LR * delta_J
					self.weights['rec'] = J

					P = P - c * (Pr.dot(Pr.T))
					
					loss_MSE_z.append(np.mean(np.square(err)))

				if t % 50 == 0 or t == timelen - 1:
					if len(loss_MSE_z) != 0:
						pbar.set_description(f'Time {t} | MSE of z {loss_MSE_z[-1]}')

			print(f'Overall MSE of z: {np.mean(loss_MSE_z)}')
			loss.append(np.mean(loss_MSE_z))

			plt.clf()
			plt.plot(r_t[250,:].T, 'r')
			plt.plot(act[250,:].T, 'k')
			plt.title(f'EP{ep+1}_{trial_type}')
			plt.draw()
			plt.pause(0.1)

		print('Done Training')
		self.p.update({'training_loss': loss})
		with open(self.pklfile, 'wb') as f:
			pkl.dump(self.p, f)
		self.save_weights()
		plt.figure(1)
		episode = np.arange(len(loss))
		plt.scatter(episode, loss)
		plt.savefig(os.path.join(self.logdir, f'MSE_of_z.png'))
		plt.close('all')
		
	def test(self, trial_types=['L', 'R'], saveplots=False):
		'''
			Testing RNN.
		'''
		readouts_r = {'L':[], 'R':[]}
		losses_r = {'L':[], 'R':[]}
		rng = npr.default_rng(11)

		for trial_type in trial_types:
			self.load_weights()			# load trained weights
			# Reset activity x
			self.reset(reset_activity=True, reset_weights=False)
			rnn_size = self.p['network_size']
			MAE_loss = []

			stim_amp = self.p[f'stim_amplitude_{trial_type}']
			stim_std = 0.1

			if stim_amp == 0:
				stochastic_stim_amp = 0.0
			else:
				# Sample stim amplitude from normal distribution
				stochastic_stim_amp = rng.normal(stim_amp, stim_std)
			
			act = self.p[f'{trial_type}_activities']
			tgt_fn = self.p[f'target_fn_{trial_type}']

			inputs = prepare_ext_inputs(act, tgt_fn, stochastic_stim_amp)
				
			timelen = inputs['target'].shape[1]
			
			r_t = np.zeros((rnn_size, timelen))
			
			print(f'Testing correct {trial_type} licking trial...')
			for t in range(timelen):
				x, z = self.run(inputs, t)
				r = self.activation(x)

				r_t[:,t] = r.flatten()
			MSE_for_r = np.mean(np.square(r - inputs['orig']), axis=1)
			losses_r[trial_type] = MSE_for_r

			print(f'Trial {trial_type} | MSE of r: {np.mean(MSE_for_r)}')
			readouts_r[trial_type] = r_t
			
		ground_truth_for_r = {'L': self.p['L_activities'], 'R': self.p['R_activities']}
		save_readouts(readouts_r, self.logdir, f'readouts_r.mat')
		save_readouts(ground_truth_for_r, self.logdir, f'ground_truth_for_r.mat')
		if saveplots:
			self.plot_readouts(trial_types, ground_truth_for_r, readouts_r, 'r_t')
		
		print('Done.')

	def plot_readouts(self, 
					  trial_types, 
					  ground_truth, 
					  readouts, 
					  save_folder_name):
		'''
			Plotting the individual fits of the readouts
		'''
		addBatch_size = self.p['additional_batch']
		ground_truth = {t: ground_truth[t][:-addBatch_size] for t in ['L', 'R']}
		readouts = {t: readouts[t][:-addBatch_size] for t in ['L', 'R']}
		if readouts['L'].shape[1] > ground_truth['L'].shape[1]:
			extended = readouts['L'].shape[1] - ground_truth['L'].shape[1]
			nulls = np.zeros((ground_truth['L'].shape[0], extended))
			ground_truth = {t: np.hstack((ground_truth[t], nulls)) for t in ['L', 'R']}
		for trial_type in trial_types:
			# Plot PSTH
			n_maps = ground_truth[trial_type].shape[0]//25 +1
			for n in range(n_maps):
				plt.figure(1)
				plt.xticks = []
				plt.rcParams['figure.figsize'] = 10,10
				fig, axes = plt.subplots(5, 5, constrained_layout=True)
				axes = axes.flatten()
				# Subplots of 25 maps
				if n < n_maps -1:
					for j in range(25):
						map_idx = (n*25)+j
						axes[j].plot(ground_truth[trial_type][map_idx], 'k')
						axes[j].plot(readouts[trial_type][map_idx], 'r')
						axes[j].set_title('{}'.format(map_idx))
						axes[j].set_xticks([])
						temp_list = [ground_truth['L'][map_idx],
									 readouts['L'][map_idx],
									 ground_truth['R'][map_idx],
									 readouts['R'][map_idx]]
						low_lim = np.min(temp_list)
						upp_lim = np.max(temp_list)
						axes[j].set_ylim([low_lim, upp_lim])
				else:
					# Subplots of maps less than 25
					m = readouts[trial_type].shape[0] - ((n_maps - 1)*25)
					for j in range(m):
						map_idx = (n*25)+j
						axes[j].plot(ground_truth[trial_type][map_idx], 'k')
						axes[j].plot(readouts[trial_type][map_idx], 'r')
						axes[j].set_title('{}'.format(map_idx))
						axes[j].set_xticks([])
						temp_list = [ground_truth['L'][map_idx],
									 readouts['L'][map_idx],
									 ground_truth['R'][map_idx],
									 readouts['R'][map_idx]]
						low_lim = np.min(temp_list)
						upp_lim = np.max(temp_list)
						axes[j].set_ylim([low_lim, upp_lim])
					# Fill the remaining subplots with blank
					for k in range(25 - m):
						axes[m+k].set_xticks([])
						axes[m+k].set_yticks([])
				print(f'{trial_type}_batch{n} done')

				folder = os.path.join(self.logdir, f'{save_folder_name}')
				if os.path.exists(folder) is False:
					os.mkdir(folder)
				plt.savefig(os.path.join(folder, f'{trial_type}_batch{n}.png'))
				plt.close(fig)

		avg_target_L = np.mean(ground_truth['L'].T, 1)
		avg_target_R = np.mean(ground_truth['R'].T, 1)
		avg_readout_L = np.mean(readouts['L'].T, 1)
		avg_readout_R = np.mean(readouts['R'].T, 1)

		title = 'Mean Spike Rate'

		fig, axes = plt.subplots(2)
		axes = axes.flatten()
		axes[0].plot(avg_target_L.flatten(), 'r')
		axes[0].plot(avg_target_R.flatten(), 'b')
		axes[0].legend(['L', 'R'])
		axes[0].set_title(f'{title} Data')
		axes[1].plot(avg_readout_L.flatten(), 'r')
		axes[1].plot(avg_readout_R.flatten(), 'b')
		axes[1].legend(['L', 'R'])
		axes[1].set_title(f'{title} RNN')
		plt.savefig(os.path.join(folder, f'{title}_Plots.svg'))
		plt.close(fig)
	
	def distract(self,
				 trial_types=['L', 'R', 'L_pos'],
				 distractor_amp=1.0,
				 distractor_std=0.1,
				 stim_std=0.1, 
				 ablation_proportion=0.0, 
				 number_of_ablation_trials=1,
				 number_of_distraction_trials=1, 
				 savedir=None,
				 ablation_type='external',
				 distractor_type='early',
				 distractor_duration=0.5):
		'''
			Generate readouts of each RNN units given trial type, in the presence of distractor.
		'''
		###################################################
		def ablation(ablation_type, ablation_proportion, ablation_rng):
			'''
				Silence the synapses with specified proportion (value between 0 to 1).
				"Internal" - randomly mask the synapses within ALM units with zeroes
				"External" - randomly mask the synapses between ALM and external coupling units with zeros
			'''
			J = self.weights['rec']
			if ablation_type == 'external':
				# Upper and Lower external synapses
				upper_J = J[:-128, -128:]
				lower_J = J[-128:, :-128]
				# Total number of external synapses available
				synap_count = np.product(upper_J.shape) + np.product(lower_J.shape)
				# number of ablated external synapses
				ablation_count = round(synap_count * ablation_proportion)
				# ablation boolean array
				ablation_bool = np.zeros(ablation_count), np.ones(synap_count - ablation_count)
				ablation_bool = np.concatenate(ablation_bool)
				ablation_rng.shuffle(ablation_bool)
				upper_ablation_bool, lower_ablation_bool = np.split(ablation_bool, 2)
				upper_ablation_bool = upper_ablation_bool.reshape(upper_J.shape)
				lower_ablation_bool = lower_ablation_bool.reshape(lower_J.shape)
				# ablation
				upper_J *= upper_ablation_bool.astype('bool')
				lower_J *= lower_ablation_bool.astype('bool')
				J[:-128, -128:] = upper_J
				J[-128:, :-128] = lower_J
			
			elif ablation_type == 'internal': 
				upper_J = J[:-128, :-128]
				# Total number of synapses available
				synap_count = np.product(upper_J.shape)
				# number of ablated synapses
				ablation_count = round(synap_count * ablation_proportion)
				# ablation boolean array
				ablation_bool = np.zeros(ablation_count), np.ones(synap_count - ablation_count)
				ablation_bool = np.concatenate(ablation_bool)
				ablation_rng.shuffle(ablation_bool)
				upper_ablation_bool = ablation_bool
				upper_ablation_bool = upper_ablation_bool.reshape(upper_J.shape)
				# ablation
				upper_J *= upper_ablation_bool.astype('bool')
				J[:-128, :-128] = upper_J
			# Update
			self.weights['rec'] = J
			print('Weights ablated.')
		#################################################

		# L_pos refers to L trial distracted with positive distractor
		rng = npr.default_rng(11)
		abl_rng = npr.default_rng(101)
		
		for m in range(number_of_ablation_trials):
			readouts_r = {'L':[], 'R':[], 'L_pos':[]}
			self.load_weights()			# load trained weights
			ablation(ablation_type, ablation_proportion, abl_rng)

			for n in range(number_of_distraction_trials):
				for trial_type in trial_types:
					
					self.reset(reset_activity=True, reset_weights=False)
					rnn_size = self.p['network_size']
					loss = []
					
					if trial_type == 'L':
						stim_amp = self.p[f'stim_amplitude_L']
						if stim_amp == 0:
							stochastic_stim_amp = 0.0
						else:
							stochastic_stim_amp = rng.normal(stim_amp, stim_std * 0.1)
						# Store the stimulus amplitude for L pos trial.
						stored_stim_L = stochastic_stim_amp					
						stochastic_dist_amp = None
						act = self.p['L_activities']

					elif trial_type == 'L_pos':
						stochastic_stim_amp = stored_stim_L
						stochastic_dist_amp = rng.normal(distractor_amp, distractor_std)
						act = self.p['L_activities']

					else:
						stim_amp = self.p[f'stim_amplitude_R']
						stochastic_stim_amp = rng.normal(stim_amp, stim_std)
						stochastic_dist_amp = None
						act = self.p['R_activities']
					
					print(f'Stim: {stochastic_stim_amp}\nDist: {stochastic_dist_amp}')
					inputs = prepare_ext_inputs(act, None, stochastic_stim_amp, stochastic_dist_amp, distractor_type, distractor_duration)

					timelen = inputs['stim'].shape[1]
					r_t = np.zeros((rnn_size, timelen))

					print(f'Generating readouts of {trial_type} licking trial...')
					for t in range(timelen):
						x, z = self.run(inputs, t)
						r = self.activation(x)

						r_t[:,t] = r.flatten()
					readouts_r[trial_type].append(r_t)

				print(f'Trial {n+1} Readouts generated.')
			if savedir is None:
				save_mat_path = os.getcwd()
			else:
				save_mat_path = savedir
			if os.path.exists(save_mat_path) is False:
				os.mkdir(save_mat_path)
			save_readouts(readouts_r, save_mat_path,
						  f'{number_of_distraction_trials}_distracted_trials_with_{ablation_proportion}_{ablation_type}_ablation_{m}.mat')


if __name__ == '__main__':
	net = RNN(create_parameters())
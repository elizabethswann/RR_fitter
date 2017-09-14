#Code to plot the corner plots and walker plots of the RR bayesian SSP fitting code

#####################################################
#Import modules
import matplotlib
matplotlib.use('agg')
import emcee_1
import corner
import matplotlib.pyplot as plt
import cPickle
import numpy as np
from math import floor
import re
import linecache
import os
import sys
import logging

###################################################################
logging.getLogger().setLevel(logging.ERROR)
###################################################################
#Load in input parameters from input file
if len(sys.argv)<3:
	print 'An input file and setup file are needed as arguments'
	print 'Input file path must be entered before setup file path'
	print 'Example: python run_RR_fitter.py /path/to/input_file.txt /path/to/setup_file.txt'
	sys.exit(0)

input_file=sys.argv[1]
setup_file=sys.argv[2]

if not os.path.isfile(input_file) or not os.path.isfile(setup_file):
	print 'Your input file or setup file does not exist. Please provide a valid file path.'
	sys.exit(0)

directory=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,10)).group(1))
datadir=directory+'data/'
savedir=directory+'saved_data/'
imagedir=directory+'images/'

num_gals=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,2)).group(1))

for index in range(num_gals):
	name_file=np.loadtxt(input_file,unpack=True,skiprows=11,dtype=str)[index]
	name_folder=os.path.splitext(name_file)[0]+'/'
	burn_in_RR=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,8)).group(1))

	if not os.path.exists(imagedir+name_folder):
		os.makedirs(imagedir+name_folder)

	###################################################################
	#Load in data from the RR bayesian SSP fitting code
	#Load in sampler chain (sampler points for walkers are in log space)
	sampler_chain=cPickle.load(
		open(savedir+name_folder+'sampler_chain_'+name_folder[:-1]+'.pkl','rb')
		)

	#Load in the log liklihoods
	log_liklihoods=cPickle.load(
		open(savedir+name_folder+'sampler_logliklihoods_'+name_folder[:-1]+'.pkl','rb')
		)

	#Load in the best fit parameters from the RR best fit (in real space)
	best_fit_params=cPickle.load(
		open(savedir+name_folder+'best_fit_params_'+name_folder[:-1]+'.pkl','rb')
		)

	ndim=len(sampler_chain[0,0,:].T)

	###################################################################
	#Plot the positions of the walkers (log space) for the RR fit
	plt.figure(1, figsize=(12, 16))
	counter=0
	for i in range(len(sampler_chain[0,0,:].T)):
		plt.subplot2grid(
			(len(sampler_chain[0,0,:].T)/2, 2),
			(int(floor(counter/2.)),int(counter%2))
			)
		plt.plot(sampler_chain[:,:,i].T,
			alpha=0.05, color='k')
		plt.ylabel(r'$t_{'+str(i+1)+'}$')
		plt.xlabel('step')
		plt.ylabel(r'ln(f$_{t_{'+str(i+1)+'}}$)')
		counter+=1
	plt.tight_layout()
	plt.savefig(imagedir+name_folder+'walker_log_'+name_folder[:-1]+'.pdf')
	
	###################################################################
	#Plot the positions of the walkers (real space) for the RR fit
	plt.figure(1, figsize=(12, 16))
	counter=0
	for i in range(len(sampler_chain[0,0,:].T)):
		plt.subplot2grid(
			(len(sampler_chain[0,0,:].T)/2, 2),
			(int(floor(counter/2.)), int(counter%2))
			)
		plt.plot(np.exp(sampler_chain[:, :, i].T),
			alpha=0.05, color='k')
		plt.ylabel(r'$t_{'+str(i+1)+'}$')
		plt.xlabel('step')
		plt.ylabel(r'f$_{t_{'+str(i+1)+'}}$')
		counter+=1
	plt.tight_layout()
	plt.savefig(imagedir+name_folder+'walker_normal_'+name_folder[:-1]+'.pdf')
	
	#################################################################
	#Plot triangle plots for parameters, in log and real space
	#Reshape the sampler chain to correct dimensions for plotting
	#Reshape the sampler into the correct shape for plotting
	samples = sampler_chain[:,burn_in_RR:,:].reshape((-1, ndim))
	labels=['$p_{1}$','$p_{2}$','$p_{3}$','$p_{4}$','$p_{5}$','$p_{6}$','$p_{7}$',
		'$p_{8}$','$p_{9}$','$p_{10}$','$p_{11}$','$p_{12}$','$p_{13}$','$p_{14}$',
		'$p_{15}$','$p_{16}$']
	
	#################################################################
	#Plot the corner plot for the RR fit in log space
	fig = corner.corner(samples,labels=labels,truths=np.log(best_fit_params),
		quantiles=[0.16,0.5,0.84],show_titles=True)
	fig.savefig(imagedir+name_folder+'triangle_log_'+name_folder[:-1]+'.pdf')
	
	#################################################################
	#Calculate the sampler chain in real space
	norm_samples=np.exp(samples)
	
	#################################################################
	#Plot the corner plot for the RR fit in real space
	fig = corner.corner(norm_samples,labels=labels,truths=best_fit_params,
		quantiles=[0.16,0.5,0.84],show_titles=True)
	fig.savefig(imagedir+name_folder+'triangle_normal'+name_folder[:-1]+'.pdf')

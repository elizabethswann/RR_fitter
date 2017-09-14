#Code to plot the corner plots and walker plots of the fullshape bayesian SSP fitting code
#that assumes a Calzetti dust model
###################################################################
#Import modules
import matplotlib
matplotlib.use('agg')
import corner
import matplotlib.pyplot as plt
import cPickle
import numpy as np
from math import floor,ceil
import re
import os
import linecache
import sys
import logging

###################################################################
logging.getLogger().setLevel(logging.ERROR)
###################################################################
#Load in input parameters from input file
if len(sys.argv)<3:
	print 'An input file and setup file are needed as arguments'
	print 'Input file path must be entered before setup file path'
	print 'Example: python run_fullshape_fitter.py /path/to/input_file.txt /path/to/setup_file.txt'
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

num_gals=boxSize=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,2)).group(1))

for index in range(num_gals):
	name_file=np.loadtxt(input_file,unpack=True,skiprows=11,dtype=str)[index]
	name_folder=os.path.splitext(name_file)[0]+'/'
	burn_in_fullshape=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,7)).group(1))

	if not os.path.exists(imagedir+name_folder):
		os.makedirs(imagedir+name_folder)	

	###################################################################
	#Load in the full chain from the bayesian fullshape fit
	sampler_chain=cPickle.load(
		open(savedir+name_folder+'sampler_chain_'+name_folder[:-1]+'_fullshape_correct_dust.pkl','rb')
		)
	#Load in the log liklihood for the bayesian fullshape fit
	log_liklihoods=cPickle.load(
		open(savedir+name_folder+'sampler_logliklihoods_'+name_folder[:-1]+'_fullshape_correct_dust.pkl','rb')
		)

	#Load in the factor from before used to renormalise the input spectrum so
	#the bayes fiter had an easier time
	factor=cPickle.load(
		open(savedir+name_folder+'logfactor_'+name_folder[:-1]+'_fullshape_correct_dust.pkl','rb')
		)

	###################################################################
	#Calculate the best fit to the fullshape solution
	index_max_liklihood=np.unravel_index(log_liklihoods.argmax(),log_liklihoods.shape)
	best_fit_params=np.exp(sampler_chain[index_max_liklihood[0],index_max_liklihood[1],1:])
	calzetti_dust=sampler_chain[index_max_liklihood[0],index_max_liklihood[1],1:]
	
	ndim=len(sampler_chain[0,0,:].T)

	###################################################################
	#Plot the walkers positions as a function fo step (log space)
	plt.figure(1, figsize=(12, 16))
	counter=0
	while counter<len(sampler_chain[0,0,:].T):
		for i in range(len(sampler_chain[0,0,:].T)):
			plt.subplot2grid((int(ceil(float(len(sampler_chain[0,0,:].T))/2.)), 2),
				(int(floor(counter/2.)), int(counter%2)))
			plt.plot(sampler_chain[:,:,i].T, alpha=0.05, color='k')
			plt.ylabel(r'$t_{'+str(i+1)+'}$')
			plt.xlabel('step')
			plt.ylabel(r'ln(f$_{t_{'+str(i+1)+'}}$)')
			counter+=1
	plt.tight_layout()
	plt.savefig(imagedir+name_folder+'walker_log_fullshape'+name_folder[:-1]+'_fullshape_calzetti_dust.pdf')
	
	###################################################################
	#Calculate the real space solutions to the fits
	exp_sampler_chain=np.copy(sampler_chain)
	#exponentiate all but the dust value which shouldn't be exponentiated
	exp_sampler_chain[:,:,1:]=np.exp(exp_sampler_chain[:,:,1:])*factor 
	
	###################################################################
	#Plot the walkers positions as a function fo step (real space)
	plt.figure(1, figsize=(12, 16))
	counter=0
	while counter<len(exp_sampler_chain[0,0,:].T):
		for i in range(len(exp_sampler_chain[0,0,:].T)):
			plt.subplot2grid((int(ceil(float(len(exp_sampler_chain[0,0,:].T))/2.)), 2),
				(int(floor(counter/2.)), int(counter%2)))
			plt.plot(exp_sampler_chain[:, :, i].T, alpha=0.05, color='k')
			plt.ylabel(r'$t_{'+str(i+1)+'}$')
			plt.xlabel('step')
			plt.ylabel(r'f$_{t_{'+str(i+1)+'}}$')
			counter+=1
	plt.tight_layout()
	plt.savefig(imagedir+name_folder+'walker_normal_fullshape'+name_folder[:-1]+'_fullshape_calzetti_dust.pdf')
	
	###################################################################
	#Plot triangle plot for the fullshape solution
	samples = sampler_chain[:,burn_in_fullshape:,:].reshape((-1, ndim)) #500
	for i in range(len(sampler_chain[:,0,0])):
		for j in range(len(sampler_chain[0,:,0])):
			sampler_chain[i,j,1:]=sampler_chain[i,j,1:]/np.sum(sampler_chain[i,j,1:])
	norm_samples=exp_sampler_chain[:,100:,:].reshape((-1,ndim))
	
	labels=['$t_{1}$','$t_{2}$','$t_{3}$','$t_{4}$','$t_{5}$','$t_{6}$','$t_{7}$',
		'$t_{8}$','$t_{9}$','$t_{10}$','$t_{11}$','$t_{12}$','$t_{13}$','$t_{14}$',
		'$t_{15}$','$t_{16}$']
	
	###################################################################
	#Plot the corner plots for the fullshape solutions
	fig = corner.corner(samples,labels=labels,quantiles=[0.16,0.5,0.84],show_titles=True)
	fig.savefig(imagedir+name_folder+'triangle_log_fullshape'+name_folder[:-1]+'_fullshape_calzetti_dust.pdf')
	
	
	fig = corner.corner(norm_samples,labels=labels,quantiles=[0.16,0.5,0.84],show_titles=True)
	fig.savefig(imagedir+name_folder+'triangle_normal_fullshape'+name_folder[:-1]+'_fullshape_calzetti_dust.pdf')

#Code to find the best fit bayesian solution to the running renormalisation
#SSP fitting method

#####################################################
#Import modules
import numpy as np
from scipy.io import readsav
from numpy import random
import emcee_1
import sys
import cPickle
from astropy.io import fits
from RR_functions import smoothMySpectrum, smoothMyError
import linecache
import re
import os
import time
###################################################################
#Define MCMC liklihood function
def lnlike(logtheta,model_arr,x,y,y_err):
	theta=np.exp(logtheta)
	test_spectra=np.dot(theta,model_arr)
	test_wave_RR,test_spectrum_RR=smoothMySpectrum(test_spectra,x,boxsize=boxSize)
	return np.sum(-(((y-test_spectrum_RR)**2.)/(2.*y_err**2.)))

#No boundaries in log space
def lnprior(logtheta):
		return 0.0

def lnprob(logtheta,model_arr,x,y,yerr):
	lp=lnprior(logtheta)
	#if not np.isfinite(lp):
	#	return -np.inf
	return lp+lnlike(logtheta,model_arr,x,y,yerr)

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

num_gals=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,2)).group(1))

#Load in SSP models used for reconstruction
models=readsav(datadir+'spectra_hist_m10_chabrier.sav') #Load in idl models
ssps=models['new_spectra']                              # Pull out SSP datacube
wave=models['wave']					#SSP wave vector

newspectra=cPickle.load(
	open(datadir+'M10_datacube_Chabrier_4ages_3metals.pkl','rb')
	)

ndim=len(newspectra)

boxSize=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,1)).group(1))
datacube=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,2)).group(1))

random_seed=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,9)).group(1))
nwalkers_RR=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,4)).group(1))
nsteps_RR=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(setup_file,6)).group(1))
survey=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,7)).group(1))

###################################################################
for index in range(num_gals):
	start_time=time.clock()
	name_file=np.loadtxt(input_file,unpack=True,skiprows=11,dtype=str)[index]
	name_folder=os.path.splitext(name_file)[0]+'/'

	np.random.seed(random_seed)

	print 'Reconstructing', name_file, 'with RR method'
###################################################################
	#Load in data or create data
	if survey in ['Fake']:
		if os.path.isfile(savedir+name_folder+'input_arr'+name_graphs+'.out')== True:
			input_arr=np.loadtxt(savedir+name_folder+'input_arr'+name_graphs+'.out')
			input_arr=input_arr/np.sum(input_arr)	#The normalisations of models to make our fake data (needs to be of length agesWanted*metalsWanted)
			spectrum_input=np.dot(input_arr,newspectra)
			wave_input=wave[:]     #slicing due to 0's at either end of models (M10 only)
		else:
			print 'For a fake galaxy, an input array is needed, stored in:', 'saved_data/'+name_folder+'input_arr'+name_graphs+'.out'
			sys.exit(0)

	elif survey in ['Lgal']:
		#Load in data from Lgal simulations
		inputdir=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,11)).group(1))
		head=fits.open(inputdir+name_file)[0]
		redshift=float(head.header['Z'])
		t=fits.open(inputdir+name_file)[1].data
	
		wave_input=t.field('wave')/(1.+redshift)
	
		wave_input=np.copy(wave) #DO NOT DO IN REAL LIFE, ONLY BECAUSE RITAS WAVELENGTH IS NOT EQUAL TO MINE AT THE 4TH DP

		#spectrum_input_dusty=t.field('flux')
		spectrum_input=t.field('flux_nodust')
	
	elif survey not in ['Lgal','Fake']:
		print 'Must choose a valid survey - Lgal or Fake'
		sys.exit(0)
	###################################################################
	#Add Calzetti Dust to input spectrum if desired
	#Define Calzetti dust model (different over different wavelength ranges)
	Es=0.44*0.3			#Arbitrarily Choseon value of Es (typical value for most gals)
	def k_lambda_small(lambda_wav):
		return (
			2.659*(-2.156 + (1.509/(lambda_wav*1e-4))
			- (0.198/((lambda_wav*1e-4)**2.))
			+ (0.011/((lambda_wav*1e-4)**3.)))
			+ 4.05
			)
			
	def k_lambda_large(lambda_wav):
		return (
			2.659*(-1.857+(1.04/(lambda_wav*1.e-4))) 
			+ 4.05
			)
	
	spectrum_input_dusty=np.zeros(len(wave_input))
	add_dust=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,4)).group(1))
	if add_dust in ['True']:
		print 'Adding dust to galaxy', name_file
		for gamma in range(len(wave_input)):
			if wave_input[gamma]*(10.**-4.) <=0.63:
				k_lambda=k_lambda_small(wave_input[gamma])
			else:
				k_lambda=k_lambda_large(wave_input[gamma])
			spectrum_input_dusty[gamma]=spectrum_input[gamma]*(10**(-0.4*Es*k_lambda))
	else:
		spectrum_input_dusty=spectrum_input
	######################################################################
	#Add white noise to our fake model with a given SNR
	add_noise=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,5)).group(1))
	if add_noise in ['True']:
		print 'Adding noise to galaxy', name_file
		SNR=float(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,6)).group(1))
		sigma = np.zeros(len(wave_input))
		e = np.zeros(len(wave_input))
		y = np.zeros(len(wave_input))
		noise=np.zeros(len(wave_input))
		for l in range(0,len(wave_input)):
			sigma[l] = (spectrum_input_dusty[l] / SNR)
			e[l] = random.normal(0.,1.)
			noise[l] = e[l]*sigma[l]
			spectrum_input_dusty[l] = spectrum_input_dusty[l] + noise[l]
	####################################################################
	#Perform a RR on our input spectrum
	input_wave_RR,input_spectrum_RR=smoothMySpectrum(spectrum_input_dusty,wave_input,boxsize=boxSize)
	input_sigma_RR=smoothMyError(spectrum_input_dusty,sigma,boxsize=boxSize)

	if not os.path.exists(savedir+name_folder):
		os.makedirs(savedir+name_folder)
	
	#Save the RR spectrum, the 'perfect' input spectrum if given, and the noisy & dusty spectrum
	np.savetxt(savedir+name_folder+'RR_spectrum'+name_folder[:-1]+'.out',
		(input_wave_RR,input_spectrum_RR,input_sigma_RR))
	
	np.savetxt(savedir+name_folder+'perfect_spectrum'+name_folder[:-1]+'.out',
		(wave_input,spectrum_input))
	
	np.savetxt(savedir+name_folder+'noisy_dusty_spectrum'+name_folder[:-1]+'.out',
		(wave_input,spectrum_input_dusty,sigma))	
	###################################################################
	#Guess initial parameters (all equal and summing to one)
	initial_guess=np.ones(ndim)
	#Initial guess is in log space
	initial_guess=np.log(initial_guess/np.sum(initial_guess))
	
	###################################################################
	#Perform the MCMC with a given number of walkers and steps to get better input parameters
	#Number of dimensions, number of walkers
	print 'Performing RR reconstruction'
	nwalkers_init=100
	#Walker initial positions start in a ball of radius 1 in log space	
	pos = [initial_guess + 1.*(2.*np.random.random(ndim)-1.) for i in range(nwalkers_init)]
	#Run the sampler for 100 steps
	sampler = emcee_1.EnsembleSampler(nwalkers_init, ndim, lnprob, 
		args=(newspectra,input_wave_RR, input_spectrum_RR, input_sigma_RR))
	sampler.run_mcmc(pos,100)
	
	#Calculate the best likelihood solution from the initial 100 steps.
	log_liklihoods=sampler.lnprobability
	index_max_liklihood=np.unravel_index(log_liklihoods.argmax(),log_liklihoods.shape)
	best_fit_params_log=sampler.chain[index_max_liklihood[0],index_max_liklihood[1],:]
	
	#Create new starting positions for the walkers in a ball about the best fit values
	#from the first 100 steps, now with a ball radius of 0.1 in log space
	new_pos=[best_fit_params_log +
		0.1*(2.*np.random.random(ndim)-1.) for i in range(nwalkers_RR)]	
	###################################################################
	#Perform the full MCMC for 1000 steps, with the new walker positions
	sampler = emcee_1.EnsembleSampler(nwalkers_RR, ndim, lnprob,
		args=(newspectra,input_wave_RR, input_spectrum_RR, input_sigma_RR))
	sampler.run_mcmc(new_pos,nsteps_RR)
	print 'RR finished on galaxy', name_file, ', saving results.'
	###################################################################
	#Save the output of the MCMC run
	
	#Save the log likelihoods for the entire chain
	f=open(savedir+name_folder+'sampler_logliklihoods_'+name_folder[:-1]+'.pkl','wb')
	cPickle.dump(sampler.lnprobability,f)
	f.close()
	
	#Save the entire chain (positions of the walkers for each step in log space) 
	f=open(savedir+name_folder+'sampler_chain_'+name_folder[:-1]+'.pkl','wb')
	cPickle.dump(sampler.chain,f)
	f.close()
	###################################################################
	#Calculate the best fit parameters (in real space)
	index_max_liklihood=np.unravel_index(log_liklihoods.argmax(),log_liklihoods.shape)
	best_fit_params=np.exp(sampler.chain[index_max_liklihood[0],index_max_liklihood[1],:])
	
	#Save the best fit parameters (in real space)
	f=open(savedir+name_folder+'best_fit_params_'+name_folder[:-1]+'.pkl','wb')
	cPickle.dump(best_fit_params,f)
	f.close()
	print 'That took',time.clock()-start_time,'seconds.'
	print '------------------------------------------------------------------------'

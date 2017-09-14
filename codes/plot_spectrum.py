#Code that plots the RR reconstruction to the spectrum , including the RR fit and the best
#fit solution to the problem in real space

###################################################################
#Import modules
import matplotlib
matplotlib.use('agg')
import emcee_1
import matplotlib.pyplot as plt
import cPickle
from scipy.io import readsav
import numpy as np
from RR_functions import smoothMySpectrum
import re
import linecache
import os
import sys
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

#Load in SSP models used for reconstruction
models=readsav(datadir+'spectra_hist_m10_chabrier.sav') #Load in idl models
ssps=models['new_spectra']                              # Pull out SSP datacube
wave=models['wave']					#SSP wave vector

newspectra=cPickle.load(
	open(datadir+'M10_datacube_Chabrier_4ages_3metals.pkl','rb')
	)

for index in range(num_gals):
	name_file=np.loadtxt(input_file,unpack=True,skiprows=11,dtype=str)[index]
	name_folder=os.path.splitext(name_file)[0]+'/'
	###################################################################
	#Load in the spectra to plot
	#Load the RR version of the input spectrum
	input_wave_RR,input_spectrum_RR,input_sigma_RR=np.loadtxt(savedir+name_folder+'RR_spectrum'+name_folder[:-1]+'.out')
	
	#Load the 'perfect' input spectrum
	wave_input,spectrum_input=np.loadtxt(savedir+name_folder+'perfect_spectrum'+name_folder[:-1]+'.out')
	
	#Load the 'noisy and dusty' input spectrum
	wave_input,spectrum_input_dusty,sigma=np.loadtxt(savedir+name_folder+'noisy_dusty_spectrum'+name_folder[:-1]+'.out')
	
	#Load in the best fit parameters from the RR best fit
	best_fit_params=cPickle.load(open(savedir+name_folder+'best_fit_params_'+name_folder[:-1]+'.pkl','rb'))
	
	###################################################################
	#Define Calzetti dust model
	Es=0.44*0.3                     #Arbitrarily Choseon value of Es (typical value for most gals)
	def k_lambda_small(lambda_wav):
		return 2.659*(-2.156 + (1.509/(lambda_wav*1e-4)) - (0.198/((lambda_wav*1e-4)**2.)) +(0.011/	((lambda_wav*1e-4)**3.))) + 4.05
	def k_lambda_large(lambda_wav):
		return 2.659*(-1.857+(1.04/(lambda_wav*1.e-4))) + 4.05
	
	frac_attenuation=np.zeros(len(wave_input))
	
	#Note: This is probably a super slow way to be doing this
	for gamma in range(len(wave_input)):
		if wave_input[gamma]*(10.**-4.) <=0.63:
			k_lambda=k_lambda_small(wave_input[gamma])
		else:
			k_lambda=k_lambda_large(wave_input[gamma])
		frac_attenuation[gamma]=10**(-0.4*Es*k_lambda)
	
	###################################################################
	#If we know the normalisation of our unattenuated spectrum (ie. fake spectra we have created)
	#then normalise it to be the same height in the range 5400AA<lambda<5600AA
	match_waverange_input=np.where([5400.<i<5600. for i in wave_input])
	normalisation_input=np.mean(spectrum_input[match_waverange_input])
	
	#Create the best fit solution spectrum and multiply it to be the same normalisation
	best_fit_spectrum=np.dot(best_fit_params,newspectra)
	match_waverange_bestfit=np.where([5400.<i<5600. for i in wave])
	normalisation_bestfit=np.mean(best_fit_spectrum[match_waverange_bestfit])
	
	normalisation=normalisation_input/normalisation_bestfit
	
	best_fit_spectrum=best_fit_spectrum*normalisation
	
	###################################################################
	#Plot spectrum and the RR reconstruction of the spectrum
	plt.figure(figsize=(8,11))
	ax1=plt.subplot(311)
	ax1.plot(wave_input,spectrum_input,'g',label='No Dust')
	ax1.plot(wave_input,spectrum_input_dusty,'r',label=r'Dust & Noise')
	ax1.plot(wave,best_fit_spectrum, color='k',label='Best Fit Spectrum')
	# plot the new best fit line
	ax1.set_ylabel('Flux (Arbitrary Scale)')
	plt.title(name_folder[:-1])
	plt.legend(loc='best')
	plt.setp(ax1.get_xticklabels(),visible=False)
	
	#Plot the ratio between the best fit and perfect input spectrum
	ax2=plt.subplot(312,sharex=ax1)
	ax2.plot(wave,best_fit_spectrum/spectrum_input)
	ax2.set_ylabel('Ratio (Best fit to Perfect Input)')
	plt.setp(ax2.get_xticklabels(),visible=False)
	
	#Plot the transmission curve of the spectrum as calculated by RR
	ax3=plt.subplot(313)
	ax3.plot(wave_input,spectrum_input_dusty/best_fit_spectrum,'r',label=r'Recovered Dust')
	ax3.plot(wave_input,frac_attenuation,'g',label='Input Dust')
	ax3.set_xlabel(r'Wavelength $\AA$')
	ax3.set_ylabel(r'Fractional Attenutation')
	plt.legend(loc='best')
	plt.savefig(imagedir+name_folder+'Spectrum_normal_recontruction_Calzetti'+name_folder[:-1]+'.pdf')
	plt.close()
	
	###################################################################
	#Recover RR input_spectrum for the bayesian fit
	fit_wave_RR,fit_spectrum_RR=smoothMySpectrum(best_fit_spectrum,wave,boxsize=41)
	
	###################################################################
	#Plot RR spectrum fit to the input RR spectrum
	ax1=plt.subplot(211)
	ax1.plot(input_wave_RR,input_spectrum_RR,'g',label='Input RR Spectrum')
	#Plot the RR best fit line
	ax1.plot(fit_wave_RR,fit_spectrum_RR, color='k',label='Best Fit RR Spectrum')
	plt.legend(loc='best')
	ax1.set_ylabel('Running Renormalised Flux')
	plt.setp(ax1.get_xticklabels(),visible=False)
	
	#Plot the residuals of the fit to the RR input spectrum
	ax2=plt.subplot(212,sharex=ax1)
	plt.plot(input_wave_RR,input_spectrum_RR-fit_spectrum_RR,color='k')
	ax2.set_xlabel(r'Wavelength $\AA$')
	ax2.set_ylabel('Residuals')
	plt.savefig(imagedir+name_folder+'Spectrum_RR_recontruction_'+name_folder[:-1]+'.pdf')
	plt.close()

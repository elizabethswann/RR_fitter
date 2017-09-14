#produces spectrum as a function
#Keywords:
#gal; Input from make_obs_galaxy.py
#setplot; can be set to 'plot' or to 'noplot'
#ps; Type of output for image. png, pdf, ps, eps and svg are supported

#imports___________________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt

#variables to change_______________________________________________________________
c=299792458 #m/s
SSP='M10'
RecycleFraction = 0.43
h=0.693
filename='testgal.sav'
wavedir='/mnt/lustre/eswann/mocks/data/' #dir where wave_M10.sav etc are stored
base_dir = '/mnt/lustre/eswann/mocks/fourierM10/' #where you want the vespa input file to be made

#Finds positions of NaNs and Infs in inarray and turns the values in those positions
#to zeros in outarray_______________________________________________________________

def rmNanInf(arr):
	arr=np.where(np.isnan(arr),0,arr)
	arr=np.where(np.isinf(arr),0,arr)
	return arr

#Crops np.ndarrays (sfh_bulge, sfh_disk etc) to have length len(ages) (number of bins?)
#and reverses order so youngest age at index 0______________________________________
	
def crop(arr,agearray):
	newarr=np.array([])
	newarr=np.append(newarr,arr[:len(agearray)])#[::-1]
	return newarr


#The produce_spectrum module________________________________________________________

def produce_spectrum(z_mass,sfh_mass,ages,Dim,Zed,models,wave_gal,w_models,ps='ps',setplot='plot'):
	#In Msun
	sfh_mass_gal=(sfh_mass)*1e10/5.2462873e+24
	z_mass=rmNanInf(z_mass)
	return_z_mass=z_mass
	
	#reconstruct the spectrum - interpolate for the metallicity in each bin
	flux_mass_rec=np.zeros(Dim)

	for i in range(0,len(ages)):
		if sfh_mass_gal[i] > 0.:
			if z_mass[i] > Zed[-1]:
				z_mass[i]=Zed[-1]
			if z_mass[i] < Zed[0]:
				z_mass[i] = Zed[0]
			#JESS LOOK AT THIS BIT BELOW!!!
			for l in range(0,Dim):
				flux_mass_rec[l] += np.interp(np.log10(z_mass[i]),np.log10(Zed),models[:,l,i])*sfh_mass_gal[i]

	flux_total_rec=flux_mass_rec

	flux=flux_total_rec
	return flux, return_z_mass

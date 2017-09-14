import subprocess
import sys
import os
import linecache
import re
import numpy as np

if len(sys.argv)<3:
	print 'An input file and setup file are needed as arguments'
	print 'Input and setup files must be stored in folder RR_fitter/input/'
	print 'Input file path name must be entered before setup file name'
	print 'Example: python run_fullshape_fitter.py /input/input_file.txt /input/setup_file.txt'
	sys.exit(0)

input_file='input/'+sys.argv[1]
setup_file='input/'+sys.argv[2]

if not os.path.isfile(input_file) or not os.path.isfile(setup_file):
	print 'Your input file or setup file does not exist. Please provide a valid file name and ensure your files are stored in the RR_fitter/input/ directory.'
	sys.exit(0)

path_base_dir=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,10)).group(1))
path_code=path_base_dir+'codes/'

num_gals=int(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,2)).group(1))
name_gals=np.loadtxt(input_file,unpack=True,skiprows=11,dtype=str)

if num_gals != len(name_gals):
	print 'Number of galaxies listed in input_file:', num_gals
	print 'Number of galaxy file names listed in input_file:', len(name_gals)
	print 'Number of galaxies given and number of listed galaxy files must be the same.'
	sys.exit(0)

input_dir=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,11)).group(1))

for i in range(len(name_gals)):
	if not os.path.isfile(input_dir+name_gals[i]):
		print 'Galaxy file:', input_dir+name_gals[i], 'does not exist.'
		print 'Please provide a valid path and filename in input_file.'
		sys.exit(0)

subprocess.check_call(['python',path_code+'bayesian_full_spectral_reconstruction.py',input_file, setup_file])

plotting=str(re.search(r"\'\s*(.*?)\s*\'",linecache.getline(input_file,9)).group(1))
if plotting in ['True']:
	print '------------------------------------------------------------------------'
	print 'Plotting best fit spectra for all galaxies'
	subprocess.check_call(['python',path_code+'plot_spectrum_fullshape.py',input_file, setup_file])
	print 'Plotting MCMC results for all galaxies'
	subprocess.check_call(['python',path_code+'plot_corner_walkers_fullshape.py',input_file, setup_file])

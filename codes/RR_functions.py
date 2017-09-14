#Code that holds the functions needed to perform RR on a spectrum

###################################################################
#Import modules
import numpy as np

###################################################################
#Define boxcar smooth functions, and error on continuum function

def doBoxSmooth(a,error_a,size):
        return (boxSmooth(a,size),boxError(a,error_a,size))

def boxSmooth(a,size):
        if size % 2 == 0.:
                print 'Choose an odd number for the boxcar smoothing width or I wont work properly'
        filter_foo=np.ones(size)/size
        return np.convolve(a, filter_foo)[size-1:-size+1]

def boxError(a,error_a,size):
        if size % 2 == 0.:
                print 'Choose an odd number for the boxcar smoothing width or I wont work properly'
        if size <2 : return a
        a=np.array(a)
        error_a=np.array(error_a)
        multiplied = (a**2.) * (error_a**2.)
        return (np.convolve(multiplied, np.ones(size)/size)[size-1:-size+1])**0.5

###################################################################
#Define smoothing function
def smoothMySpectrum(spectrum,wavelength,boxsize=41):
	spectrum_continuum=boxSmooth(spectrum,boxsize)
	spectrum_RR=np.array(spectrum[(boxsize-1)/2:-((boxsize-1)/2)])/np.array(spectrum_continuum)
	model_x=wavelength[(boxsize-1)/2:-((boxsize-1)/2)]
	return model_x,spectrum_RR

def smoothMyError(spectrum,spectrum_err,boxsize=41):
	spectrum_continuum=boxSmooth(spectrum,boxsize)
	continuum_error=np.array(boxError(spectrum,spectrum_err,boxsize))
	spectrum_RR_err=(
		(1./spectrum_continuum)
		*((spectrum_err[(boxsize-1)/2:-((boxsize-1)/2)])**2.
		+ ((spectrum[(boxsize-1)/2:-((boxsize-1)/2)] * continuum_error)**2. / spectrum_continuum**2.)
		)**0.5
		)
	return spectrum_RR_err
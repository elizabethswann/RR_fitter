# RR_fitter
Code to perform RR fitting method on spectra to recover fractional star formation histories and transmission functions

## Getting Started
Download the entire RR_fitter folder to your local machine.

The input spectra to be reconstructed using the RR_fitter should be placed in input/name_of_galaxy_spectra_folder/

The input files (also contained in input/) should be edited to contain the name of the spectra files you wish to reconstruct, their paths and a series of optional parameters.

The input and setup files should **not have their structure changed**. The input and setup parameters are read in line by line, do not change the order of parameters in the files.

The setup file should be edited if you wish to change number of walkers used in the RR MCMC fit, the burn in period, the datacube of SSP models used to reconstruct the spectrum, etc.

The RR or fullshape fitters can then be run using
```
python run_RR_fitter.py input_file.txt setup_file.txt
```

### Prerequisite python modules

* astropy
* corner
* matplotlib
* numpy
* scipy

## Testing RR_fitter

To test the RR_fitter works, run on the example Lgal galaxies ```gal_spectrum_364_MGS_M10.fits``` and ```gal_spectrum_750_MGS_M10.fits``` using the provided example setup files. You can check you get the same output as that contained in the folder example_output.

Implement this test using the RR MCMC spectral fitter with the command

```
python run_RR_fitter.py input_Lgal_example.txt setup_Lgal_example.txt
```

and the below for the fullshape MCMC spectral fitter

```
python run_fullshape_fitter.py input_Lgal_example.txt setup_Lgal_example.txt
```

## Authors

* **Elizabeth Swann**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* This code uses the emcee python package, created by Dan Foreman-Mackey and contributors, and is free software made avaliable under the MIT license (https://github.com/dfm/emcee). Note that this module is not listed as a prerequisite as it is contained as part of the RR_fitter in ```code/emcee``` and ```code/emcee_1``` as slight changes have been made for use with RR_fitter.
* Thanks to Prof Peter Thomas and Dr Rita Tojeiro for all their contributions to this work.

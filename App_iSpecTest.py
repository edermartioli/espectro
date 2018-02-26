# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Feb 23 2018
    
    Description: This routine tests the use of ipsec software
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_iSpecTest.py --input=spectrum.m.fits
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum

import numpy as np
import matplotlib.pyplot as plt

################################################################################
#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = './iSpec_v20161118/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec
################################################################################

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-p", action="store_true", dest="polar", help="polar spectrum", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_iSpecTest.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output wavelength range: ', options.wlrange
    print 'Spectrum type: ', options.spectype
    print 'Polar option: ', options.polar
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

spc = Spectrum(options.input, options.spectype, options.polar, options.telluric, options.helio)
spc.continuumOrderByOrder(binsize = 30, overlap = 10)
spc.normalizeByContinuum(constant=10.0)
spc.binningOrderByOrder(rvsampling_kps=2.5, median=False)
spc.sortByWavelength()
ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)

resolution = 80000
smoothed_star_spectrum = ispec.convolve_spectrum(ispectrum, resolution)

mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
ccf_mask = ispec.read_cross_correlation_mask(mask_file)

models, ccf = ispec.cross_correlate_with_mask(smoothed_star_spectrum, ccf_mask, \
                                              lower_velocity_limit=-50, upper_velocity_limit=50, \
                                              velocity_step=0.05, mask_depth=0.01, \
                                              fourier=False)

components = len(models)

rv = np.round(models[0].mu(), 3) # km/s
rv_err = np.round(models[0].emu(), 3) # km/s

print components, rv, rv_err

'''
suntemplate = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
smoothed_suntemplate = ispec.convolve_spectrum(suntemplate, resolution)

models2, ccf2 = ispec.cross_correlate_with_template(smoothed_star_spectrum, smoothed_suntemplate, \
                                                  lower_velocity_limit=-50, upper_velocity_limit=50, \
                                                  velocity_step=0.1, fourier=False)

components2 = len(models2)
# First component:
rv2 = np.round(models2[0].mu(), 2) # km/s
rv_err2 = np.round(models2[0].emu(), 2) # km/s

print components2, rv2, rv_err2
'''

smoothed_star_spectrum = ispec.correct_velocity(smoothed_star_spectrum, rv)

velocity = ccf['x']
xcorr = ccf['y']
xcorr_err = ccf['err']

#ispec.plot_spectra([smoothed_star_spectrum, smoothed_suntemplate])

plt.plot(velocity, xcorr)
plt.show()

# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python

    Created on Feb 23 2018
    
    Description: Radial velocity time series
    
    @author: Eder Martioli
    
    INPE / Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_RadialVelocity.py --inputdir=./espectrosflux/ --spectype=norm --object="AM Her" -tr
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from spectralclass import Spectrum
from spectralclass import SpectrumChunk
import espectrolib

import numpy as np

################################################################################
#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = './iSpec_v20161118/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec
################################################################################

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)
parser.add_option("-r", action="store_true", dest="helio", help="heliocentric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_RadialVelocity.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Object name: ', options.object
    print 'Spectrum type: ', options.spectype
    print 'Telluric correction: ', options.telluric
    print 'Heliocentric correction: ', options.helio

filelist = espectrolib.generateList(options.inputdir, options.object)

mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
ccf_mask = ispec.read_cross_correlation_mask(mask_file)
resolution = 80000

def barycentric_velocity_corrected_spectrum(spc):
    #--- Barycentric velocity correction from observation date/coordinates ---------

    day = spc.MJDTime.datetime.day
    month = spc.MJDTime.datetime.month
    year = spc.MJDTime.datetime.year
    hours = spc.MJDTime.datetime.hour
    minutes = spc.MJDTime.datetime.minute
    seconds = spc.MJDTime.datetime.second
    
    ra_hours = spc.obj_coord.ra.hms[0]
    ra_minutes = spc.obj_coord.ra.hms[1]
    ra_seconds = spc.obj_coord.ra.hms[2]
    dec_degrees = spc.obj_coord.dec.dms[0]
    dec_minutes = spc.obj_coord.dec.dms[1]
    dec_seconds = spc.obj_coord.dec.dms[2]
    
    ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)
    
    # Project velocity toward star
    barycentric_vel = ispec.calculate_barycentric_velocity_correction((year, month, day, \
                                        hours, minutes, seconds), (ra_hours, ra_minutes, \
                                        ra_seconds, dec_degrees, dec_minutes, dec_seconds))
    #--- Correcting barycentric velocity -------------------------------------------
    corrected_spectrum = ispec.correct_velocity(ispectrum, barycentric_vel)

    return corrected_spectrum

output = []

for filepath in filelist :
    head, filename = os.path.split(filepath)
    spc = Spectrum(filepath, options.spectype, False, options.telluric, options.helio)

    barycentricTime = spc.barycentricTime()

    spc.continuumOrderByOrder(binsize = 30, overlap = 10)
    spc.normalizeByContinuum(constant=10.0)
    spc.binningOrderByOrder(rvsampling_kps=2.4, median=False)
    spc.sortByWavelength()
    
    ispectrum = barycentric_velocity_corrected_spectrum(spc)
    #ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)
    
    smoothed_star_spectrum = ispec.convolve_spectrum(ispectrum, resolution)

    models, ccf = ispec.cross_correlate_with_mask(smoothed_star_spectrum, ccf_mask, \
                                              lower_velocity_limit=-100, upper_velocity_limit=100, \
                                              velocity_step=1.0, mask_depth=0.01, \
                                              fourier=False)
    '''
    rv0 = np.round(models0[0].mu()) # km/s


    models, ccf = ispec.cross_correlate_with_mask(smoothed_star_spectrum, ccf_mask, \
                                                lower_velocity_limit=rv0-2.0, upper_velocity_limit=rv0+2.0, \
                                                velocity_step=0.005, mask_depth=0.01, \
                                                fourier=False)
    '''
    
    rv = models[0].mu()
    rv_err = models[0].emu()

    output.append([filename, spc.object.replace(" ",""), barycentricTime.jd, rv, rv_err])

for result in output :
    print result[0], result[1], np.round(result[2], 6), np.round(result[3], 3), np.round(result[4], 3)



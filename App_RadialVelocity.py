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

#--- Calculate barycentric velocity correction from observation date/coordinates ---------
def barycentric_velocity(spc):
    
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
    
    # Project velocity toward star
    barycentric_vel = ispec.calculate_barycentric_velocity_correction((year, month, day, \
                                                                       hours, minutes, seconds), (ra_hours, ra_minutes, \
                                                                                                  ra_seconds, dec_degrees, dec_minutes, dec_seconds))
                                                                                                  
    return barycentric_vel
#-------

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-t", action="store_true", dest="telluric", help="telluric correction", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_RadialVelocity.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Object name: ', options.object
    print 'Spectrum type: ', options.spectype
    print 'Telluric correction: ', options.telluric

filelist = espectrolib.generateList(options.inputdir, options.object)

mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
ccf_mask = ispec.read_cross_correlation_mask(mask_file)
resolution = 80000

output = []

for filepath in filelist :
    head, filename = os.path.split(filepath)
    
    if options.verbose: print "Loading spectrum file: ", filename

    spc = Spectrum(filepath, options.spectype, False, options.telluric, False)

    barycentricTime = spc.barycentricTime()
    if options.verbose: print "Barycentric Time: ",barycentricTime

    barycentric_vel = barycentric_velocity(spc)
    if options.verbose: print "Barycentric velocity correction: ",barycentric_vel

    if options.verbose: print "Calculating continuum, binning and sorting data ..."
    spc.continuumOrderByOrder(binsize = 30, overlap = 10)
    spc.normalizeByContinuum(constant=10.0)
    spc.binningOrderByOrder(rvsampling_kps=2.4, median=False)
    spc.sortByWavelength()
    
    #ispectrum = barycentric_velocity_corrected_spectrum(spc)
    ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)
    
    if options.verbose: print "Smoothing out spectrum to resolution: ", resolution
    smoothed_star_spectrum = ispec.convolve_spectrum(ispectrum, resolution)

    #--- Correcting barycentric velocity -------------------------------------------
    if options.verbose: print "Correcting barycentric velocity ..."
    corrected_spectrum = ispec.correct_velocity(smoothed_star_spectrum, -barycentric_vel)
    #-----

    if options.verbose: print "Calculating radial velocity by cross-correlation with mask ..."
    models, ccf = ispec.cross_correlate_with_mask(smoothed_star_spectrum, ccf_mask, \
                                              lower_velocity_limit=-100, upper_velocity_limit=100, \
                                              velocity_step=0.1, mask_depth=0.01, \
                                              fourier=False)
    rv = models[0].mu()
    rv_err = models[0].emu()

    if options.verbose:
        print "Writing output:"
        print filename, spc.object.replace(" ",""), np.round(barycentricTime.jd, 6), np.round(rv, 5), np.round(rv_err, 3)
    
    output.append([filename, spc.object.replace(" ",""), barycentricTime.jd, rv, rv_err])


for result in output :
    print result[0], result[1], np.round(result[2], 6), np.round(result[3], 5), np.round(result[4], 3)



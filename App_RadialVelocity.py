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
    python $PATH/App_RadialVelocity.py --inputdir=./espectrosflux/ --wavemask=./auxdata/ESPaDOnS_OrderWavelengthRangesForRV.txt --spectype=raw --object="AM Her" -tr
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
import matplotlib.pyplot as plt
from astropy.io import fits

################################################################################
#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = './iSpec_v20161118/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec
################################################################################

#--- Function to write cross-correlation plot to output fits ---------
def writeoutputfits(ccf,filename,rv,rv_err,bjd,snr):
    
    if ".fits.gz" in filename :
        basename = os.path.splitext(os.path.splitext(filename)[0])[0]
    elif ".fits" in filename :
        basename = os.path.splitext(filename)[0]

    outputfits = basename + ".rv.fits"

    scidata = []
    scidata.append(np.array(ccf['x']))
    scidata.append(np.array(ccf['y']))
    scidata.append(np.array(ccf['err']))

    hdu = fits.PrimaryHDU(scidata)
    hdulist = fits.HDUList([hdu])
    hdulist[0].header.set('BJD',np.round(bjd,7), 'barycentric julian date (d)')

    hdulist[0].header.set('RV',rv, 'radial velocity (kps)')
    hdulist[0].header.set('RVERR',rv_err, 'radial velocity error (kps)')

    hdulist[0].header.set('SNR',snr, 'mean signal-to-noise ratio')

    hdulist[0].header.set('COL1','Radial Velocity', 'radial velocity (kps)')
    hdulist[0].header.set('COL2','Cross-correlation', 'cross-correlation')
    hdulist[0].header.set('COL3','Cross-correlation Error', 'cross-correlation error')

    hdulist[0].header['COMMENT'] = "Data processed by ESPECTRO 1.0 -> App_RadialVelocity.py"

    if os.path.exists(outputfits) :
        os.remove(outputfits)
    hdulist.writeto(outputfits)
#-------


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

#--- Load template spectrum for cross correlation ---------
def loadtemplatespectrum(filename, resolution, telluric_linelist, ccf_mask, velocity_step, wavemask):
    
    spc = Spectrum(filename)
    
    if wavemask :
        spc.applyMask(wavemask)

    barycentricTime = spc.barycentricTime()
    barycentric_vel = barycentric_velocity(spc)
    
    spc.continuumOrderByOrder(binsize = 30, overlap = 10)
    spc.normalizeByContinuum(constant=10.0)
    
    spc.binningOrderByOrder(rvsampling_kps=2.4, median=False)
    spc.sortByWavelength()

    ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)

    smoothed_spectrum = ispec.convolve_spectrum(ispectrum, resolution)

    tel_models, tel_ccf = ispec.cross_correlate_with_mask(smoothed_spectrum, telluric_linelist, \
                                                        lower_velocity_limit=-10, upper_velocity_limit=10, \
                                                        velocity_step=0.1, mask_depth=0.01, \
                                                        fourier = False, only_one_peak = True)
    tel_rv = tel_models[0].mu()
    
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(smoothed_spectrum['waveobs'], \
                                                                    telluric_linelist[dfilter], min_velocity=-tel_rv+min_vel, \
                                                                    max_velocity=-tel_rv+max_vel)
    
    
    clean_spectrum = smoothed_spectrum[~tfilter]


    corrected_spectrum = ispec.correct_velocity(clean_spectrum, tel_rv - barycentric_vel)

    models, ccf = ispec.cross_correlate_with_mask(clean_spectrum, ccf_mask, \
                                                  lower_velocity_limit=-100, upper_velocity_limit=100, \
                                                  velocity_step=velocity_step, mask_depth=0.01, \
                                                  fourier=False)
    rv = models[0].mu()

    clean_corrected_spectrum = ispec.correct_velocity(corrected_spectrum, rv)

    return clean_corrected_spectrum
#-------

parser = OptionParser()
parser.add_option("-i", "--inputdir", dest="inputdir", help="Input directory with spectral data",type='string', default="")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string', default="")
parser.add_option("-s", "--spectype", dest="spectype", help="Spectrum type: raw, norm, or fcal",type='string', default="raw")
parser.add_option("-m", "--wavemask", dest="wavemask", help="File with each order wavelength range",type='string', default="")
parser.add_option("-t", "--template", dest="template", help="Input template spectrum",type='string', default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_RadialVelocity.py -h "; sys.exit(1);

if options.verbose:
    print 'Input directory: ', options.inputdir
    print 'Object name: ', options.object
    print 'Spectrum type: ', options.spectype
    print 'Input template spectrum: ', options.template
    print 'File with wavelength ranges: ', options.wavemask

velocity_step = 0.2

filelist = espectrolib.generateList(options.inputdir, options.object)

#mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/Atlas.Arcturus.372_926nm/mask.lst""
#mask_file = ispec_dir + "input/linelists/CCF/Atlas.Sun.372_926nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.A0.350_1095nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.F0.360_698nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.G2.375_679nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K0.378_679nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K5.378_680nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.M5.400_687nm/mask.lst"
#mask_file = ispec_dir + "input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst"
mask_file = ispec_dir + "input/linelists/CCF/VALD.Sun.300_1100nm/mask.lst"

ccf_mask = ispec.read_cross_correlation_mask(mask_file)

telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

resolution = 80000

#--- Below it loads and sets up the template spectrum ---------
#templatefile = "/Users/edermartioli/Desktop/51Peg/testdata/1603681o.m.fits.gz"
templatefile = "/Users/edermartioli/Desktop/51Peg/testdata/1604052o.m.fits.gz"

if options.template :
    templatefile = options.template
    templatespectrum = loadtemplatespectrum(templatefile, resolution, telluric_linelist, ccf_mask, velocity_step, options.wavemask)
#---------------

output = []

for filepath in filelist :
    head, filename = os.path.split(filepath)
    
    if options.verbose: print "Loading spectrum file: ", filename

    spc = Spectrum(filepath, options.spectype)
    
    if options.wavemask :
        spc.applyMask(options.wavemask)

    snrmask = np.where(spc.fluxerr > 1.0)
    snr = np.mean(spc.flux[snrmask] / spc.fluxerr[snrmask])
    if options.verbose: print "Mean SNR: ",snr

    barycentricTime = spc.barycentricTime()
    if options.verbose: print "Barycentric Time: ",barycentricTime

    barycentric_vel = barycentric_velocity(spc)
    if options.verbose: print "Barycentric velocity correction: ",barycentric_vel

    if options.spectype != "norm" :
        if options.verbose: print "Calculating continuum and normalizing spectrum ..."
        spc.continuumOrderByOrder(binsize = 30, overlap = 10)
        spc.normalizeByContinuum(constant=10.0)

    if options.verbose: print "Binning and sorting the data ..."
    spc.binningOrderByOrder(rvsampling_kps=2.4, median=False)
    spc.sortByWavelength()

    ispectrum = ispec.create_spectrum_structure(spc.wl,spc.flux)

    if options.verbose: print "Smoothing out spectrum to resolution: ", resolution
    smoothed_star_spectrum = ispec.convolve_spectrum(ispectrum, resolution)


    #--- Calculating RV shift from telluric lines -------------------------------------------
    if options.verbose: print "Calculating RV shift from telluric lines ..."
    tel_models, tel_ccf = ispec.cross_correlate_with_mask(smoothed_star_spectrum, telluric_linelist, \
                                              lower_velocity_limit=-10, upper_velocity_limit=10, \
                                              velocity_step=0.1, mask_depth=0.01, \
                                              fourier = False, only_one_peak = True)
    
    tel_rv = tel_models[0].mu()
    tel_rv_err = tel_models[0].emu()
    #plt.errorbar(tel_ccf['x'],tel_ccf['y'],yerr=tel_ccf['err'], fmt='o')
    if options.verbose: print "Telluric RV shift = ",tel_rv,"+/-",tel_rv_err
    #-----

    #--- Cleaning telluric regions -------------------------------------------
    if options.verbose: print "Cleaning telluric regions ..."
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(smoothed_star_spectrum['waveobs'], \
                                                                    telluric_linelist[dfilter], min_velocity=-tel_rv+min_vel, \
                                                                    max_velocity=-tel_rv+max_vel)
    clean_spectrum = smoothed_star_spectrum[~tfilter]
    #-----

    #--- Correcting telluric shift -------------------------------------------
    if options.verbose: print "Correcting RV shift calculated from telluric lines ..."
    corrected_spectrum = ispec.correct_velocity(clean_spectrum, tel_rv)
    #-----

    #--- Correcting barycentric velocity -------------------------------------------
    if options.verbose: print "Correcting barycentric velocity ..."
    corrected_spectrum = ispec.correct_velocity(corrected_spectrum, -barycentric_vel)
    #-----

    #--- Calculating cross-correlation for RV measurements ------------------------
    if options.template :
        if options.verbose: print "Calculating radial velocity by cross-correlation with template spectrum ..."
        models, ccf = ispec.cross_correlate_with_template(corrected_spectrum, templatespectrum,  \
                                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                                            velocity_step=velocity_step, fourier=False)
    
    else :
        if options.verbose: print "Calculating radial velocity by cross-correlation with lines mask ..."
        models, ccf = ispec.cross_correlate_with_mask(corrected_spectrum, ccf_mask, \
                                                  lower_velocity_limit=-100, upper_velocity_limit=100, \
                                                  velocity_step=velocity_step, mask_depth=0.01, \
                                                  fourier=False)
    rv = models[0].mu()
    rv_err = models[0].emu()

    if options.verbose: print "RV shift = ",rv,"+/-",rv_err

    #plt.errorbar(ccf['x'],ccf['y'],yerr=ccf['err'], fmt='o')
    #plt.show()
    #-----

    writeoutputfits(ccf, filename, rv, rv_err, barycentricTime.jd, snr)

    if options.verbose:
        print "Writing output:"
        print filename, spc.object.replace(" ",""), np.round(barycentricTime.jd, 6), np.round(rv, 5), np.round(rv_err, 3), snr
    
    output.append([filename, spc.object.replace(" ",""), barycentricTime.jd, rv, rv_err, snr])

for result in output :
    print result[0], result[1], np.round(result[2], 6), np.round(result[3], 5), np.round(result[4], 3), np.round(result[5], 1)



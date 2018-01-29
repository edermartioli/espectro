# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 2017
@author: Eder Martioli
Description: ESPECTRO Library
Laboratorio Nacional de Astrofisica, Brazil
"""

import os
import numpy as np
from astropy.io.fits import getheader
from scipy import stats
from scipy.interpolate import UnivariateSpline
import glob
import gzip

import matplotlib.pyplot as plt

######################
def get_fitsfilepaths(directory):
    
    """
    Generates a list of file names in a directory tree
    by walking the tree either top-down or bottom-up.
        
    Parameters
    ----------
    directory : directory path
        
    Returns
    -------
    file_paths: a list of file paths
    """
        
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree:
    for root, directories, files in os.walk(directory):
        for filename in files:
			# Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            if filename.endswith(".fits") :
                file_paths.append(filepath)

    file_paths.sort()
    return file_paths
######################

############# get basename from FITS file path ###############
def getbasename(filepath) :
    
    basename = ''
    base = os.path.basename(filepath)
    
    dir = os.path.dirname(filepath)
    
    if ".fits.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    elif ".fits" in base :
        basename = os.path.splitext(base)[0]
    elif ".txt" in base :
        basename = os.path.splitext(base)[0]
    elif ".spc.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    else :
        print "Error: unknown extension in file: ", filepath
        exit()

    return basename, dir
###########################################

############# get Heliocentric Radial Velocity from header ####
def getEspaonsHelioRV(header) :
    comments = header['COMMENT']
    rv = 0.0
    for c in comments :
        if "Heliocentric velocity" in c :
            rv = float(c[c.find(": ")+1:c.find(" km")])
            break
    return rv
###########################################

############# FFT filtering ###############
def fft_filter(y, step=0.01, samp=20):
    """
        Perform Fast Fourier Transform filtering of data
        
        Parameters
        ----------
        y : data array
        step : step to define frequency sampling
        samp : number of sample points

        Returns
        -------
        yfftclean: filtered data array
    """

    fourier = np.fft.fft(y)
    n = len(fourier)
    freq = np.fft.fftfreq(n, d=step)
    iy = np.zeros(n, dtype=complex)
    
    
    for j in xrange(n):
        if -samp < freq[j] < samp:
            iy[j] = fourier[j]
        else:
            iy[j] = 0
    yfftclean = np.real(np.fft.ifft(iy, n))
    
    return yfftclean
#####

############# Remove Background  ###############
def removeBackground(wl,flux,wlranges):
    """
        Remove background of data
        
        Parameters
        ----------
        wl : wavelength data array
        flux : flux data array
        wlranges : wavelength ranges where to estimate background
        
        Returns
        -------
        wl,flux: reduced data arrays
    """

    mask1 = np.where(np.logical_and(wl > wlranges[0][0], wl < wlranges[0][1]))
    mask2 = np.where(np.logical_and(wl > wlranges[1][0], wl < wlranges[1][0]))
    
    flux_bkg = np.append(flux[mask1],flux[mask2])
    
    background = np.median(flux_bkg)
    
    flux -= background
    
    return wl,flux
#####

############# Multiple gaussian function  ###############
def gaussfunc(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y
#####

######################
def generateList(directory, objectname, objectkey="OBJECT"):
    
    """
    Generates a list of file names in a directory tree
    by walking the tree either top-down or bottom-up.
        
    Parameters
    ----------
    directory : directory path
        
    Returns
    -------
    file_paths: a list of file paths
    """
        
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree:
    for root, directories, files in os.walk(directory):
        for filename in files:
	    # Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            
            if filename.endswith("m.fits.gz") or filename.endswith("pol.fits.gz") \
               or filename.endswith("m.fits") or filename.endswith("pol.fits"):
                header = getheader(filepath.rstrip('\n'),0)
                objname = header[objectkey]
                if objname == objectname :
                    file_paths.append(filepath)

    file_paths.sort()
    return file_paths
######################

######################
def generateListWithRVs(directory, objectname, objectkey="OBJECT"):
    
    """
        Generates a list of m.fits.gz file names in a given directory
        and also gives a second list of rv2.gz files, containing
        source radial velocities which are associated to the
        spectra.
        
        Parameters
        ----------
        directory : directory path
        objectname : object name
        objectkey: header keyword for object name
        
        Returns
        -------
        mfiles, rvfiles: two lists of file paths
        """
    
    
    pattern = directory + '*.m.fits.gz'

    files = glob.glob(pattern)

    mfiles = []
    rvfiles = []
    
    for file in files :
        header = getheader(file,0)
        objname = header[objectkey]
        
        if objname.replace(' ','') == objectname.replace(' ','') :
            rvfilepath = file.strip('m.fits.gz') + '.rv2.gz'
            if os.path.exists(rvfilepath) :
                rvfiles.append(rvfilepath)
                mfiles.append(file)

    return mfiles, rvfiles
######################

######################
def wlrange(wlstr, spc):
    
    """
        Get wavelength range from input command line string string
        
        Parameters
        ----------
        wlstr : input range in string format: "wl0 wlf". E.g. "500 550"
        spc   : Spectrum() class
        
        Returns
        -------
        wavelength floats: wl0, wlf
        """

    wl0=0.0
    wlf=0.0
    
    if wlstr :
        wl0 = float(wlstr.split()[0])
        wlf = float(wlstr.split()[1])
    else :
        wl0 = spc.wl[0]
        wlf = spc.wl[-1]
    
    return wl0, wlf
######################

######################
def continuum(wl, flux, binsize = 30, overlap = 20, sigmaclip = 2.5, window = 3, linesmask=None) :

    lines = []
    if linesmask :
        lines = loadwlranges(linesmask)
 
    nbins = len(wl)/binsize
        
    wlbin = []
    fluxbin = []
    
    for i in range(nbins):
            
        idx0 = i*binsize - overlap
        idxf = (i+1)*binsize + overlap
            
        if idx0 < 0 : idx0=0
        if idxf > len(wl) : idxf=len(wl)-1
            
        wltmp = wl[idx0:idxf]
        fluxtmp = flux[idx0:idxf]
            
        wlmean = np.mean(wltmp)
        
        medflux = np.median(fluxtmp)
        medfluxdev = np.median(np.abs(fluxtmp - medflux))
        filtermask = np.where(np.logical_and(fluxtmp > medflux, fluxtmp < medflux + sigmaclip*medfluxdev))
        
        withinLine = False
        for line in lines :
            if wlmean > line[0] and wlmean < line[1] :
                withinLine = True
                break
        if not withinLine :
            wlbin.append(wlmean)
            #fluxbin.append(np.max(fluxtmp[filtermask]))
            fluxbin.append(np.median(fluxtmp[filtermask]))

    newwlbin = []
    newfluxbin = []
        
    newwlbin.append(wl[0])
    newfluxbin.append(fluxbin[0])
        
    for i in range(len(wlbin)):
        idx0 = i - window
        idxf = i + 1 + window
            
        if idx0 < 0 : idx0=0
        if idxf > nbins : idxf=nbins-1
            
        slope, intercept, r_value, p_value, std_err = stats.linregress(wlbin[idx0:idxf], fluxbin[idx0:idxf])
            
        newwlbin.append(wlbin[i])
        newfluxbin.append(intercept + slope*wlbin[i])

    #plt.plot(newwlbin,newfluxbin,marker="o")

    s = UnivariateSpline(newwlbin, newfluxbin, s=0)
    continuum = s(wl)
    return continuum
######################

######################
def loadSourceRVs(rvfilelist):
    """
        Load a list of RV files, read RV values and return list of values
        
        Parameters
        ----------
        rvfilelist : input list of RV files
        
        Returns
        -------
        rvs: list of RV values
        """
    rvs = []
    
    for file in rvfilelist :
        with gzip.open(file, 'rb') as f:
            file_content = f.read()
        rv = float(file_content.split('\n')[6].split(' ')[7])
        rvs.append(rv)
 
    return rvs
######################

######################
def loadwlranges(wlrangesfile):
    """
        Load a list of wavelength ranges from file
        
        Parameters
        ----------
        wlrangesfile : input file with list of wavelength ranges
        
        Returns
        -------
        wlranges: list of wavelength ranges values
        """
    wlranges = []

    file = open(wlrangesfile, "r")
    for line in file:
        wl0wlf = [float(line.split(' ')[0]),float(line.split(' ')[1])]
        wlranges.append(wl0wlf)
    
    return wlranges
######################

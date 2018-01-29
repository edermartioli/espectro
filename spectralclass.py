# -*- coding: utf-8 -*-

"""
Spectral Classes
---------------------------
Created on Mar 10 2017
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import os
import numpy as np
from astropy.io import fits
import espectrolib
from scipy import constants
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy import stats
from scipy.signal import medfilt
import copy

########## SPECTRUM CLASS ############
class Spectrum :
    'Common base class for a spectrum'

    def __init__(self, Filename, Spectype="raw", polar=False, Telluric=False, Helio=False, SourceRV=0.0, Sort=True):

        """
        Create a Spectrum object.
        
        Parameters
        ----------
        filename : string
            File to read the spectrum from.

        Examples
        --------
        >>> spc = Spectrum("spectrum.m.fits.gz")
        """
        self.sourceRV = SourceRV
        
        self.telluric = Telluric
        self.helio = Helio
        self.spectype = Spectype
        self.filepath = Filename
        basename = espectrolib.getbasename(self.filepath)
        
        self.id = basename[0]
        self.filename = os.path.basename(self.filepath)
        
        self.hascontinuum = False
        
        try :
            if self.filepath.endswith(".m.fits.gz") or self.filepath.endswith(".m.fits"):
                self.wl,self.flux,self.fluxerr=self.loadSpectrumFromMFITS(self.filepath, self.spectype, self.telluric, self.helio, sort=Sort)
            elif self.filepath.endswith(".pol.fits.gz") or self.filepath.endswith(".pol.fits") :
                self.wl,self.flux,self.fluxerr=self.loadSpectrumFromPOLFITS(self.filepath, self.spectype, polar, self.telluric, self.helio, sort=Sort)
            elif self.filepath.endswith(".s.fits.gz") or self.filepath.endswith(".s.fits") :
                self.wl,self.flux,self.fluxerr=self.loadSpectrumFromSFITS(self.filepath, sort=Sort)
            else :
                print "Error: file type not supported for input spectrum: ",self.filepath
                exit()
        except :
            print "Error: could not open file: ",self.filepath
            exit()

    #--- Function to load spectrum from .s.fits file
    def loadSpectrumFromSFITS(self, fitsfilename, sort=False):
        
        hdu = fits.open(fitsfilename)
        
        self.header = hdu[0].header
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
        self.HJDTT = hdu[0].header['HJDTT']
        self.exptime = hdu[0].header['EXPTIME']

        wl = hdu[0].data[0]
        flux = hdu[0].data[1]
        fluxerr = hdu[0].data[2]
        
        if sort :
            indices = wl.argsort()
            return wl[indices], flux[indices], fluxerr[indices]
        else :
            return wl, flux, fluxerr

    #--- Function to load spectrum from .m.fits.gz file
    def loadSpectrumFromMFITS(self,fitsfilename, spectype="raw", telluric=False, helio=False, sort=True):

        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        
        self.header = hdu[0].header
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
        self.HJDTT = hdu[0].header['HJDTT']
        self.exptime = hdu[0].header['EXPTIME']

        self.orders = hdu[0].data[0]

        if spectype == "raw" :
            influx = hdu[0].data[8]
            influxvar = hdu[0].data[9]
        elif spectype == "norm" :
            influx = hdu[0].data[10]
            influxvar = hdu[0].data[11]
        elif spectype == "fcal" :
            influx = hdu[0].data[12]
            influxvar = hdu[0].data[13]
        else :
            print "Error: spectrum type not recognized: ",spectype
            exit()

        influxerr = np.sqrt(influxvar)

        if telluric :
            wltmp = hdu[0].data[5]
        else :
            wltmp = hdu[0].data[4]

        if helio :
            wltmp += hdu[0].data[6]

        wltmp *= (1.0 - self.sourceRV*1000.0/constants.c)

        if sort :
            indices = wltmp.argsort()
            wl = wltmp[indices]
            flux = influx[indices]
            fluxerr = influxerr[indices]
            self.orders = self.orders[indices]
        else :
            wl = wltmp
            flux = influx
            fluxerr = influxerr

        return wl,flux,fluxerr
    #------------

    #--- Function to load spectrum from .pol.fits.gz file
    def loadSpectrumFromPOLFITS(self,fitsfilename, spectype="raw", polar=False, telluric=False, helio=False, sort=True):
    
        wl,flux = [],[]
        hdu = fits.open(fitsfilename)
        self.header = hdu[0].header
        
        self.instrument = hdu[0].header['INSTRUME']
        self.object = hdu[0].header['OBJECT']
        self.HJDTT = hdu[0].header['HJDTT']
        self.exptime = hdu[0].header['EXPTIME']

        self.orders = hdu[0].data[0]

        if polar :
            if spectype == "raw" or spectype == "fcal":
                influx = hdu[0].data[13]
                influxvar = hdu[0].data[14]
            elif spectype == "norm":
                influx = hdu[0].data[15]
                influxvar = hdu[0].data[16]
            else :
                print "Error: spectrum type not recognized: ",spectype
                exit()
        else :
            if spectype == "raw" :
                influx = hdu[0].data[7]
                influxvar = hdu[0].data[8]
            elif spectype == "norm" :
                influx = hdu[0].data[9]
                influxvar = hdu[0].data[10]
            elif spectype == "fcal" :
                influx = hdu[0].data[11]
                influxvar = hdu[0].data[12]
            else :
                print "Error: spectrum type not recognized: ",spectype
                exit()

        influxerr = np.sqrt(influxvar)

        if telluric :
            wltmp = hdu[0].data[4]
        else :
            wltmp = hdu[0].data[3]

        if helio :
            wltmp += hdu[0].data[5]

        wltmp *= (1.0 - self.sourceRV*1000.0/constants.c)
        
        if sort :
            indices = wltmp.argsort()
            wl = wltmp[indices]
            flux = influx[indices]
            fluxerr = influxerr[indices]
            self.orders = self.orders[indices]
        else :
            wl = wltmp
            flux = influx
            fluxerr = influxerr

        return wl, flux, fluxerr
    #------------
    
    #--- resampling spectrum
    def resampling(self, wlsampling, wl0, wlf) :
        npoints = int((wlf-wl0)/wlsampling)
        wl_new = np.linspace(wl0, wlf, npoints)
        flux_new = np.interp(wl_new, self.wl, self.flux)
        self.wl = wl_new
        self.flux = flux_new
    #------------

    #--- resample spectrum
    def resample(self, wl_new) :
        
        flux_new = copy.copy(self.flux)
        fluxvar_new = copy.copy(self.fluxerr)
        
        minorder,maxorder = self.getMinMaxOders()
        for order in range(minorder, maxorder+1) :
            ordermask = np.where(self.orders == order)
            flux_new[ordermask] = np.interp(wl_new[ordermask], self.wl[ordermask], self.flux[ordermask])
            fluxvar_new[ordermask] = np.interp(wl_new[ordermask], self.wl[ordermask], self.fluxerr[ordermask]*self.fluxerr[ordermask])

        self.wl = wl_new
        self.flux = flux_new
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------

    #--- apply median filter
    def applyMedianFilter(self, window) :
        self.flux = medfilt(self.flux, kernel_size=window)
    #--------------------------

    #--- Print spectrum information
    def info(self) :
        print "#**************************"
        print "#Info for spectrum: ",self.filename
        print "#Object:", self.object
        print "#Instrument:",self.instrument
        print "#wl0,wlf =",self.wl[0],",",self.wl[-1],"nm"
        sampling = (self.wl[-1] - self.wl[0])/float(len(self.wl))
        print "#sampling =",sampling," nm/pixel"
        print "#Wave telluric correction =",self.telluric
        print "#Wave heliocentric correction =",self.helio
        print "#Flux type =",self.spectype
        print "#HJD (TT) =",self.HJDTT
        print "#EXPTIME = ",self.exptime, "s"
        print "#<F> =",self.flux.mean(),"+-",self.flux.std()
        print "#**************************"
    #------------

    #--- Print spectrum data
    def printdata(self, printerrors=True) :
        
        if printerrors :
            for i in range(len(self.wl)) :
                print self.wl[i],self.flux[i],self.fluxerr[i]
        else :
            for i in range(len(self.wl)) :
                print self.wl[i],self.flux[i]
    #------------

    #--- Extract spectral range
    def extractChunk(self, wl0, wlf, printdata=False) :
        mask = np.where(np.logical_and(self.wl > wl0, self.wl < wlf))

        if printdata :
            wlout = self.wl[mask]
            fluxout = self.flux[mask]
            fluxerrout = self.fluxerr[mask]
            for i in range(len(wlout)) :
                print wlout[i], fluxout[i], fluxerrout[i]

        return self.wl[mask],self.flux[mask],self.fluxerr[mask]
    #------------

    #--- Get time from img header
    def getTimeHJDTT(self) :
        header = fits.getheader(self.filepath,0)
        return header['HJDTT']
    #------------

    #--- Mask data
    def maskdata(self, lines, width) :
        for line in lines :
            wl0 = line - width
            wlf = line + width
            mask = np.where(np.logical_or(self.wl < wl0, self.wl > wlf))
            self.wl = self.wl[mask]
            self.flux = self.flux[mask]
            self.fluxerr = self.fluxerr[mask]
    #------------

    #--- bin spectrum
    def binning(self, rvsampling_kps, wl0=0.0, wlf=0.0, median=False) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        fluxvar = self.fluxerr*self.fluxerr
        
        bins = []
        wl = wl0
        while wl <= wlf :
            bins.append(wl)
            wl *= (1.0 + rvsampling_kps*1000/constants.c)
        bins = np.array(bins)
        
        digitized = np.digitize(self.wl, bins)

        wl_new = []
        flux_new = []
        fluxvar_new = []

        for i in range(1, len(bins)):
            if len(self.wl[digitized == i]) :
                try:
                    wl_new.append(self.wl[digitized == i].mean())
                    if median :
                        flux_new.append(np.median(self.flux[digitized == i]))
                    else :
                        flux_new.append(self.flux[digitized == i].mean())
                    fluxvar_new.append(self.flux[digitized == i].std())
                except :
                    continue

        self.wl = np.array(wl_new)
        self.flux = np.array(flux_new)
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------
    
    #--- bin spectrum for a single order
    def binOrder(self, order, rvsampling_kps, wl0=0.0, wlf=0.0, median=False) :

        ordermask = np.where(self.orders == order)

        if wl0 == 0.0:
            wl0 = self.wl[ordermask][0]
        if wlf == 0.0:
            wlf = self.wl[ordermask][-1]
        
        fluxvar = self.fluxerr[ordermask]*self.fluxerr[ordermask]

        bins = []
        wl = wl0
        while wl <= wlf :
            bins.append(wl)
            wl *= (1.0 + rvsampling_kps*1000/constants.c)
        
        npbins = np.array(bins)
        
        digitized = np.digitize(self.wl[ordermask], npbins)
        
        wl_new = []
        flux_new = []
        fluxvar_new = []
        orders_new = []
        
        for i in range(1, len(npbins)):
            if len(self.wl[ordermask][digitized == i]) :
                try:
                    wl_new.append(self.wl[ordermask][digitized == i].mean())
                    if median :
                        flux_new.append(np.median(self.flux[ordermask][digitized == i]))
                    else :
                        flux_new.append(self.flux[ordermask][digitized == i].mean())
                    fluxvar_new.append(self.flux[ordermask][digitized == i].std())
                    orders_new.append(float(order))
                except :
                    continue
        
        self.wl[ordermask] = np.array(wl_new)
        self.flux[ordermask] = np.array(flux_new)
        self.fluxerr[ordermask] = np.sqrt(fluxvar_new)
        self.orders[ordermask] = np.array(orders_new)
    #--------------------------

    #--- bin spectrum
    def binningOrderByOrder(self, rvsampling_kps, median=False) :
        minorder,maxorder = self.getMinMaxOders()
        for order in range(minorder, maxorder+1) :
            self.binOrder(order, rvsampling_kps=rvsampling_kps, median=median)
    #--------------------------

    #--- Extract order
    def extractOrder(self, order) :
        ordermask = np.where(self.orders == order)
        wlout = copy.copy(self.wl[ordermask])
        fluxout = copy.copy(self.flux[ordermask])
        fluxerrout = copy.copy(self.fluxerr[ordermask])
        return wlout,fluxout,fluxerrout
    #--------------------------

    #--- Extract order continuum
    def extractOrderContinuum(self, order) :
        ordermask = np.where(self.orders == order)
        if self.hascontinuum :
            return self.continuum[ordermask]
        else :
            return None
    #--------------------------

    #--- Calculate continuum for a single order
    def ordercontinuum(self, order, binsize=30, overlap=10, linesmask=None) :
        #--- Extract order
        orderwl, orderflux, orderfluxerr = self.extractOrder(order)
        ordercontinuum = espectrolib.continuum(orderwl, orderflux, binsize = binsize, overlap = overlap, linesmask=linesmask)
        return ordercontinuum
    #--------------------------

    #--- Calculate continuum order by order
    def continuumOrderByOrder(self, binsize = 30, overlap = 10, linesmask=None) :
        
        self.continuum = copy.copy(self.flux)
        
        minorder,maxorder = self.getMinMaxOders()
        for order in range(minorder, maxorder+1) :
            continuumflux = self.ordercontinuum(order, binsize, overlap, linesmask=linesmask)
            ordermask = np.where(self.orders == order)
            self.continuum[ordermask] = continuumflux

        self.hascontinuum = True
    #--------------------------

    #--- filter by orders
    def getMinMaxOders(self) :
        
        minorder = np.min(self.orders)
        maxorder = np.max(self.orders)

        return int(minorder), int(maxorder)
    #--------------------------

    #--- Set all fluxes to 0.0
    def setAllFluxesToZero(self) :
        self.flux = np.zeros_like(self.flux)
        self.fluxerr = np.zeros_like(self.flux)
    #--------------------------
    
    #--- Sort spectrum by wavelength
    def sortByWavelength(self) :
        indices = (self.wl).argsort()
        self.wl = self.wl[indices]
        self.flux = self.flux[indices]
        self.fluxerr = self.fluxerr [indices]
        self.orders = self.orders[indices]
    #--------------------------
    
    #--- apply mask (keep only data within the mask)
    def applyMask(self, wavemaskfile) :
        
        file = open(wavemaskfile, "r")

        for line in file:
            o = float(line.split(' ')[0])
            wl0 = float(line.split(' ')[1])
            wlf = float(line.split(' ')[2])
            
            logi_wl = np.logical_or(self.wl < wl0, self.wl > wlf)
            logi_order =  np.logical_and(self.orders>o-1,self.orders<o+1)
            cutmask = np.where(np.logical_not(np.logical_and(logi_wl,logi_order)))
            self.wl = self.wl[cutmask]
            self.flux = self.flux[cutmask]
            self.fluxerr = self.fluxerr[cutmask]
            self.orders = self.orders[cutmask]
            if self.hascontinuum :
                self.continuum = self.continuum[cutmask]
    #--------------------------

    #--- save spectrum to output file
    def saveToFile(self, output, format='fits') :

        if format == 'fits' :
            self.header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            self.header['ORIGDATA'] = self.id
            self.header['TELLCORR'] = self.telluric
            self.header['HERVCORR'] = self.helio
            self.header['SPECTYPE'] = self.spectype
            self.header['SOURCERV'] = self.sourceRV

            self.header.set('COL1','Wavelength', 'wavelength (nm)')
            self.header.set('COL2','Flux', 'flux (electron)')
            self.header.set('COL3','FluxError', 'flux error (electron)')

            self.header.remove('COL4')
            self.header.remove('COL5')
            self.header.remove('COL6')
            self.header.remove('COL7')
            self.header.remove('COL8')
            self.header.remove('COL9')
            self.header.remove('COL10')
            self.header.remove('COL11')
            self.header.remove('COL12')
            self.header.remove('COL13')
            self.header.remove('COL14')
            
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.flux)
            scidata.append(self.fluxerr)
            
            if os.path.exists(output) :
                os.remove(output)
            fits.writeto(output, scidata, self.header)
    
        elif format == 'ascii' :
            
            '''
            from astropy.io import ascii
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.flux)
            ascii.write(scidata, output)
            '''
            
            f = open(output, 'w')
            strdata = ''
            for i in range(len(self.wl)) :
                wlstr = "%.5f" % self.wl[i]
                fluxstr = "%.5f" % self.flux[i]
                strdata += wlstr + ' ' + fluxstr + '\n'
            f.write(strdata)
            f.close()
            
        else :
            print 'file format ',format, ' not supported'
            exit()
    #--------------------------

    #--- calculate continuum
    def continuum(self, binsize = 200, overlap = 100, sigmaclip = 3.0, window = 3) :

        nbins = len(self.wl)/binsize

        wlbin = []
        fluxbin = []

        for i in range(nbins):
    
            idx0 = i*binsize - overlap
            idxf = (i+1)*binsize + overlap
    
            if idx0 < 0 : idx0=0
            if idxf > len(self.wl) : idxf=len(self.wl)-1
    
            wltmp = self.wl[idx0:idxf]
            fluxtmp = self.flux[idx0:idxf]
    
            wlmean = np.mean(wltmp)
    
            wlbin.append(wlmean)
            medflux = np.median(fluxtmp)
            medfluxdev = np.median(np.abs(fluxtmp - medflux))
            filtermask = np.where(np.logical_and(fluxtmp > medflux, fluxtmp < medflux + sigmaclip*medfluxdev))
    
            fluxbin.append(np.max(fluxtmp[filtermask]))

        newwlbin = []
        newfluxbin = []

        newwlbin.append(self.wl[0])
        newfluxbin.append(fluxbin[0])

        for i in range(nbins):
            idx0 = i - window
            idxf = i + 1 + window
    
            if idx0 < 0 : idx0=0
            if idxf > nbins : idxf=nbins-1
    
            slope, intercept, r_value, p_value, std_err = stats.linregress(wlbin[idx0:idxf], fluxbin[idx0:idxf])
    
            newwlbin.append(wlbin[i])
            newfluxbin.append(intercept + slope*wlbin[i])


        s = UnivariateSpline(newwlbin, newfluxbin, s=1)
        self.continuum = s(self.wl)
        self.hascontinuum = True
    #--------------------------

    #--- normalize by continuum
    def normalizeByContinuum(self, constant=0.0) :
        try :
            self.flux = (self.flux + constant) / (self.continuum + constant)
        except :
            print "Could not perfom continuum normalization."
    #--------------------------

    #--- normalize order by order
    def normalizeOrderByOrder(self, bkgbin = 20) :
        minorder, maxorder = self.getMinMaxOders()
        
        for order in range(minorder, maxorder+1) :
        
            ordermask = np.where(self.orders == order)
            
            orderflux = copy.copy(self.flux[ordermask])
            sortedorderflux = sorted(orderflux)
            bkg = np.mean(sortedorderflux[bkgbin:2*bkgbin])
        
            self.flux[ordermask] = copy.copy((orderflux - bkg)/np.median(orderflux - bkg))
    #--------------------------

    #--- save continuum spectrum to output file
    def saveContinuumToFile(self, output, format='fits') :
    
        if format == 'fits' and self.hascontinuum :
            self.header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            self.header['ORIGDATA'] = self.id
            
            self.header.set('COL1','Wavelength', 'wavelength (nm)')
            self.header.set('COL2','Flux', 'flux (electron)')
            self.header.set('COL3','FluxError', 'flux error (electron)')
            self.header.set('COL4','Order', 'order number')

            self.header.remove('COL5')
            self.header.remove('COL6')
            self.header.remove('COL7')
            self.header.remove('COL8')
            self.header.remove('COL9')
            self.header.remove('COL10')
            self.header.remove('COL11')
            self.header.remove('COL12')
            self.header.remove('COL13')
            self.header.remove('COL14')
            
            scidata = []
            scidata.append(self.wl)
            
            scidata.append(self.continuum)
            scidata.append(self.fluxerr)
            scidata.append(self.orders)
            
            if os.path.exists(output) :
                os.remove(output)
            fits.writeto(output, scidata, self.header)
        else :
            if self.hascontinuum :
                print 'file format ',format, ' not supported'
            else:
                print 'No continuum has been calculated'
            exit()
    #--------------------------

########## SPECTRUM CHUNK CLASS ############
class SpectrumChunk :
    'Common base class for a spectrum chunk'

    def __init__(self, Wl, Flux, Fluxerr):
        self.wl,self.flux,self.fluxerr = Wl, Flux, Fluxerr

    def removeBackground(self, backgroundsize, wl0=0.0, wlf=0.0) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]
        wlranges = [(float(wl0),float(wl0)+backgroundsize),(float(wlf)-backgroundsize,float(wlf))]
        self.wl,self.flux = espectrolib.removeBackground(self.wl,self.flux,wlranges)


    def snrfilter(self,snrcut) :
        snr = self.flux / self.fluxerr 
        mask = np.where(snr > snrcut)
        self.wl,self.flux,self.fluxerr = self.wl[mask],self.flux[mask],self.fluxerr[mask]

    #--- resampling spectrum
    def resampling(self, wlsampling, wl0=0.0, wlf=0.0) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        npoints = int((wlf-wl0)/wlsampling)
        wl_new = np.linspace(wl0, wlf, npoints)
        flux_new = np.interp(wl_new, self.wl, self.flux)
        fluxvar = self.fluxerr*self.fluxerr
        fluxvar_new = np.interp(wl_new, self.wl, fluxvar)

        self.wl = wl_new
        self.flux = flux_new
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------

    #--- bin spectrum
    def binning(self, wlsampling, wl0=0.0, wlf=0.0, median=False) :
        if wl0 == 0.0:
            wl0 = self.wl[0]
        if wlf == 0.0:
            wlf = self.wl[-1]

        fluxvar = self.fluxerr*self.fluxerr

        npoints = int((wlf-wl0)/wlsampling)
        bins = np.linspace(wl0, wlf, npoints)
        digitized = np.digitize(self.wl, bins)

        wl_new = [self.wl[digitized == i].mean() for i in range(1, len(bins))]

        if median :
            flux_new = [np.median(self.flux[digitized == i]) for i in range(1, len(bins))]
        else :
            flux_new = [self.flux[digitized == i].mean() for i in range(1, len(bins))]

        fluxvar_new = [self.flux[digitized == i].std() for i in range(1, len(bins))]

        self.wl = wl_new
        self.flux = flux_new
        self.fluxerr = np.sqrt(fluxvar_new)
    #--------------------------

    #----------- FFT filtering -----------
    def fft_filter(self, step=0.01, samp=5):
        """
            Perform Fast Fourier Transform filtering of flux data
        
            Parameters
            ----------
            step : step to define frequency sampling
            samp : number of sample points 
        """
        fourier = np.fft.fft(self.flux)
        n = len(fourier)
        freq = np.fft.fftfreq(n, d=step)
        iy = np.zeros(n, dtype=complex)
    
        for j in xrange(n):
            if -samp < freq[j] < samp:
                iy[j] = fourier[j]
            else:
                iy[j] = 0
        self.flux = np.real(np.fft.ifft(iy, n))
    #-----------

    #----------- Fit Gaussian -----------
    def fitgaussian(self, guess) :

        popt, pcov = curve_fit(espectrolib.gaussfunc, self.wl, self.flux, p0=guess)

        self.model = espectrolib.gaussfunc(self.wl, *popt)

        wlcen = []
        amp = []
        sigma = []

        for i in range(0, len(popt), 3) :
            wlcen.append(popt[i])
            amp.append(popt[i+1])
            sigma.append(popt[i+2])
 
        return wlcen, amp, sigma
    #-----------

    #----------- get chunk spectrum -----------
    def getSpectrum(self) :
        return self.wl, self.flux, self.fluxerr
    #-----------

    #----------- print data and model -----------
    def printdataWithModel(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i], self.model[i]
    #-----------

    #----------- print data -----------
    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i], self.flux[i], self.fluxerr[i]
    #-----------

    #--- save spectrum to output file
    def saveToFile(self, output, format='fits') :
        if format == 'fits' :
            scidata = []
            scidata.append(self.wl)
            scidata.append(self.flux)
            scidata.append(self.fluxerr)

            hdu = fits.PrimaryHDU(scidata)
            hdulist = fits.HDUList([hdu])
            
            hdulist[0].header.set('COL1','Wavelength', 'wavelength (nm)')
            hdulist[0].header.set('COL2','Flux', 'flux (electron)')
            hdulist[0].header.set('COL3','FluxError', 'flux error (electron)')
            hdulist[0].header['COMMENT'] = "Data processed by ESPECTRO 1.0"
            if os.path.exists(output) :
                os.remove(output)
            hdulist.writeto(output)
        else :
            print 'file format ',format, ' not supported'
            exit()
    #--------------------------

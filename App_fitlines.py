# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Oct 31 2017
    
    Description: Fit gaussian to a set of lines
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    python $PATH/App_fitlines.py --input=spectrum.s.fits
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

from scipy import constants
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from spectralclass import Spectrum
from spectralclass import SpectrumChunk

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum",type='string', default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_fitlines.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input

lines = [402.619, 412.082, 447.148, 471.315, 492.193, 501.568, 587.562, 667.815, 706.519]

#lines = [447.148, 471.315, 492.193, 501.568, 587.562, 667.815, 706.519]

delta_rv_kps = 400

gaussian = lambda x,a,c,w: a * np.exp(-((x - c)/w)**2)
gaussfunc = np.vectorize(gaussian)

spc = Spectrum(options.input, Sort=False)

shift = 0.0
dshift = 0.015

for line in lines :
    wl_line = line

    delta_wl = 1000.0 * delta_rv_kps * wl_line / constants.c

    wl0, wlf = wl_line-delta_wl, wl_line+delta_wl
    
    wl,flux,fluxerr = spc.extractChunk(float(wl0), float(wlf))
    chunk = SpectrumChunk(wl,flux,fluxerr)
    chunk.removeBackground(0.1)
    init_guess = [wl_line,-0.1, 0.3]
    wlcen, amp, sigma = chunk.fitgaussian(init_guess)
    delta = wlcen[0] - wl_line
    
    print "line=",wl_line, wlcen[0], delta

    linemodel = gaussfunc(wl,amp[0],wlcen[0],sigma[0])
    
    plt.plot(wl-wlcen[0],flux+shift, marker=".")
    plt.plot(wl-wlcen[0],linemodel+shift)

    shift += dshift

plt.show()





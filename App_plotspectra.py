# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
    Created on Mar 29 2017
    
    Description: Plot product spectrum 
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_plotspectra.py --input='spectrum1.s.fits spectrum1.s.fits'
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """
from optparse import OptionParser
import os,sys

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectrum file(s)",type='string', default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_plotspectra.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum file(s): ', options.input

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

inputimages = options.input.split()

for img in inputimages:
    hdu = fits.open(img)
    wl = hdu[0].data[0]
    flux = hdu[0].data[1]
    plt.plot(wl,flux)

plt.show()

# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
    #!/usr/bin/python
    #!/opt/anaconda/bin/python
    Created on Jan 24 2018
    
    Description: Convert a 2-column (wl flux) .txt spectrum to IRAF FITS file
    
    @author: Eder Martioli
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python $PATH/App_exportToIRAFFITS.py --input=spectrum.txt --output=spectrum.fits 
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
from pyraf import iraf

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input 2-column .txt spectrum file",type='string', default="")
parser.add_option("-o", "--output", dest="output", help="Output IRAF FITS spectrum file",type='string', default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with App_extract.py -h "; sys.exit(1);

if options.verbose:
    print 'Input spectrum: ', options.input
    print 'Output spectrum: ', options.output

iraf.rspectext(options.input, options.output, dtype='interp')
#iraf.splot(options.output)



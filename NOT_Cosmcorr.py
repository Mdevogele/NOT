#!/usr/bin/env python

""" RP_CosmCorr - Correct for Cosmic rays
    v1.0: 2018-10-02, mdevogele@lowell.edu
"""

import argparse, shlex

#import RP_Toolbox as RP
import numpy as np

#from SP_CheckInstrument import CheckInstrument

from astropy.io import fits

import cosmic

def Cosmic(filename):
    
#    telescope, obsparam = CheckInstrument([filename[0]])
    for idx, elem in enumerate(filename): 
        hdulist = fits.open(elem)
        data=hdulist[1].data
        c = cosmic.cosmicsimage(data,satlevel=100000)
        c.run(maxiter = 1)
        hdulist[1].data = c.cleanarray
        hdulist.writeto(elem.replace('.fits','').replace('_Procc','') + '_' + '_CosmCorr' + '.fits')
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Create a log file')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    filenames = args.images    
    
    # call run_the_pipeline only on filenames
    Cosmic(filenames)
    pass

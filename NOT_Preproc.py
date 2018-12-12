#!/usr/bin/env python

""" NOT_Preproc - Apply the bias to science data
    v1.0: 2018-10-02, mdevogele@lowell.edu        
"""
import argparse

from astropy.io import fits
import numpy as np

def Preproc(filenames,MasterBias,Verbose,Suffix,MasterFlat):

    
    if MasterBias:
        if Verbose: 
            print('Opening master bias file: ' + str(MasterBias))            
        hdulist = fits.open(MasterBias)
        Bias_data = hdulist[1].data  

    if MasterFlat:
        if Verbose: 
            print('Opening master flat file: ' + str(MasterFlat))            
        hdulist = fits.open(MasterBias)
        Flat_data = hdulist[1].data
        
        
    for idx, elem in enumerate(filenames):
        hdulist = fits.open(elem)
        data = hdulist[1].data[100:300,100:300]

        if MasterBias and not MasterFlat:
            data = data - Bias_data
        if not MasterBias and MasterFlat:
            data = data/Flat_data 
        if MasterBias and MasterFlat:
            data = (data-Bias_data)/Flat_data 

        hdulist[1].data = data
        hdulist.writeto(elem.split('.')[0] + '_' +  Suffix + '.fits',overwrite = True)
        if Verbose: 
            print('{} \t {} processed'.format(idx+1,elem))


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Preprocessing of the NOT polarimetric files')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
    parser.add_argument('-s',
                        help='Suffix to add to processed files',
                        default='Procc')
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Default to None if no bias',
                        default=None)
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    
    parser.add_argument('-f',
                        help='Name of the master flat to use \n || || Default to None if no flat',
                        default=None)    

    args = parser.parse_args()
    Suffix = args.s
    MasterBias = args.b
    Verbose = args.v
    filenames = args.images  
    MasterFlat = args.f

    
    print(filenames)
    
    Preproc(filenames,MasterBias,Verbose,Suffix,MasterFlat)
    pass
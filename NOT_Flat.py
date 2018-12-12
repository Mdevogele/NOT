#!/usr/bin/env python

""" NOT_Flat - Function to create master flat for NOT polarimetric data
    v1.0: 2018-10-02, mdevogele@lowell.edu
"""

import argparse

import numpy as np

from astropy.io import fits


def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


class Result(object):
    def Indiv_to_Master(self, method, Flat):
        Master = getattr(self, method, lambda x: getattr(self, 'Default')(Flat))(Flat)
        return Master
 
    def Median(self,Flat):
        Master = np.median(Flat,axis=0)/np.median(Flat)
        return Master
 
    def Mean(self,Flat):
        Master = np.mean(Flat,axis=0)/np.mean(Flat)
        return Master

    def Default(self, Flat):
        Master = np.median(Flat,axis=0)/np.median(Flat)
        print("Invalid method, use of the median as default")
        return Master 


def Create_Flat(filenames,MasterName,Verbose,Method,Bias,Bin):
    
    
    if Bias:
        if Verbose: 
            print('Opening master bias file')            
        hdulist = fits.open(Bias)
        Bias_data = hdulist[1].data    


    Flat = []
    if Verbose:
        print('Beginning flats processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))
            hdulist = fits.open(elem)
            data = hdulist[1].data[750:1500, 640:1390] # Polarimetric images are cropped
            data = rebin(data,(int(np.shape(data)[0]/Bin),int(np.shape(data)[1]/Bin)))
            data = data[100:300,100:300]     # Cropping the data to avoid the vignetting in science frames  
            Flat.append(data - Bias_data)
        print('Creating the master flat')
    else:
        for idx,elem in enumerate(filenames):
            hdulist = fits.open(elem)
            data = hdulist[1].data[750:1500, 640:1390] # Polarimetric images are cropped
            data = rebin(data,(int(np.shape(data)[0]/Bin),int(np.shape(data)[1]/Bin)))
            data = data[100:300,100:300]     # Cropping the data to avoid the vignetting in science frames  
            Flat.append(data - Bias_data)
    
    
    Res = Result()
    MasterFlat = Res.Indiv_to_Master(Method, Flat)   
    
    print(len(MasterFlat))
    
    hdulist[1].data = MasterFlat

    hdulist.writeto(MasterName, overwrite = True)
    hdulist.close()     
        

    if Verbose:
        print('Master flat save to {}'.format(MasterName))
        hdulist = fits.open(MasterName)
        data = hdulist[1].data
        print('Statistics of the Master flat')
        print('Mean: {} \t Median: {} \t std: {}'.format(np.mean(data), np.median(data), np.std(data)))
        print('End of flat processing')
        hdulist.close()
        
    




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master bias')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    
    parser.add_argument('-o',
                        help='Prefix of the master flats files',
                        default='MasterFlat.fits')  
    
    parser.add_argument('-m',
                        help='Method to use to compute the master bias: Mean, Median',
                        default='Median') 
    
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Default to None (no bias)',
                        default=None)
    
    parser.add_argument('-bin',
                        help='Binning factor',
                        default=2) 
    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()

    Verbose = args.v
    MasterName = args.o
    filenames = args.images    
    Method = args.m 
    Bias = args.b
    Bin = int(args.bin) 

    
    Create_Flat(filenames,MasterName,Verbose,Method,Bias,Bin)
    pass
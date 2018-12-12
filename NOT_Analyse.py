#!/usr/bin/env python

""" CAPS_Anal - Extract best q and u from curve of growth
    v1.0: 2018-05-05, mdevogele@lowell.edu
"""

import matplotlib.pyplot as plt
import argparse
import numpy as np
import operator
from scipy.optimize import curve_fit
from astroquery.jplhorizons import Horizons

def Anal(filenames,Plot,Target,Aperture):

    Res = []
    JD = []
    Alpha = []
    PlAng = []
    All = []
    retarder =[]
    Stoke_Final = []
    for elem in filenames:
        stokes = []
        with open(elem,'r') as f:
            for elem in f.readlines():
                print(elem)
                stokes.append(float(elem.replace('\n',' ').replace('\t',' ').split()[2]))
            
        Stoke_Final.append(stokes)
        JD.append(float(elem.replace('\n',' ').replace('\t',' ').split()[0]))
        retarder.append(float(elem.replace('\n',' ').replace('\t',' ').split()[1]))

        obj = Horizons(id=Target, location='950', epochs=float(elem.replace('\n',' ').replace('\t',' ').split()[0]))
        eph = obj.ephemerides()
        Alpha.append(eph['alpha'][0])
        PlAng.append(eph['sunTargetPA'][0])  
        All.append(stokes)
        if Plot:    
            plt.plot(stokes)
            
            
#        Res.append(np.median(stokes[Aperture]))  


    All = np.array(All)
    print(np.sqrt((1./2*(All[0,:]-All[2,:]))**2+(1./2*(All[1,:]+All[3,:]))**2))
    plt.plot(np.sqrt((1./2*(All[0,:]-All[2,:]))**2+(1./2*(All[1,:]-All[3,:]))**2))
    
    plt.show() 

    if not Aperture:
        Aperture = int(raw_input("What aperture to you want to use? "))

    print()
  
    Res = All[:,Aperture]
    
    with open('result.txt', 'w') as f:
        for idx,elem in enumerate(Res):
            f.write(str(retarder[idx]) + '\t' + str(JD[idx]) + '\t' + str(Alpha[idx]) + '\t' + str(PlAng[idx]) + '\t' + str(elem) + '\n')
        
    

if __name__ == '__main__':
    
    
    # define command line arguments
    parser = argparse.ArgumentParser(description='manual target identification')
    parser.add_argument('-aper', default = None, help='Aperture to use to compute the polarization')
    parser.add_argument('-plot', action="store_true")
    parser.add_argument('-object', help='Name of the target for retrieving alpha and scaterring plane angle values',
                        default=False) 
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    
    filenames = args.images
    aperture = args.aper
    Plot = args.plot
    Target = args.object.replace('_',' ')


    Anal(filenames,Plot,Target,aperture)


    pass
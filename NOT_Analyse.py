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


def EXP(x, *p):
    S0,S1,S2 = p
    bb = S0+S1*np.exp(-S2*x)
    
    return bb

def EXP2(x, *p):
    S0,S1,S2,S3,S4 = p
    bb = S0+S1*np.exp(-S2*x)+S3*np.exp(-S4*x)
    
    return bb

def Anal(filenames,Auto,Plot,Target):

    Res1 = []
    Res2 = []
    Res = []
    JD = []
    Alpha = []
    PlAng = []
    All = []
    retarder =[]
    for elem in filenames:
        stokes = []
        with open(elem,'r') as f:
            for elem in f.readlines():
                print(elem)
                stokes.append(float(elem.replace('\n',' ').replace('\t',' ').split()[2]))
            
        stokes = np.array(stokes)
        JD.append(float(elem.replace('\n',' ').replace('\t',' ').split()[0]))
        retarder.append(float(elem.replace('\n',' ').replace('\t',' ').split()[1]))

        Res.append(np.median(stokes[0:3]))
        obj = Horizons(id=Target, location='950', epochs=float(elem.replace('\n',' ').replace('\t',' ').split()[0]))
        eph = obj.ephemerides()
        Alpha.append(eph['alpha'][0])
        PlAng.append(eph['sunTargetPA'][0])  
        All.append(stokes)
        if Plot:    
#            plt.figure()
            plt.plot(stokes)
    All = np.array(All)
    print(np.sqrt((1./2*(All[0,:]-All[2,:]))**2+(1./2*(All[1,:]+All[3,:]))**2))
    plt.plot(np.sqrt((1./2*(All[0,:]-All[2,:]))**2+(1./2*(All[1,:]-All[3,:]))**2))

        
#        xs = range(1,18)
#        p0 = [np.median(stokes),0,1]
#        print(p0)
#        coeff, var_matrix = curve_fit(EXP,xs,stokes,p0,maxfev=100000)
#        print(coeff)
#        fit = EXP(np.array(xs), *coeff)
#        if Plot:
#            plt.figure()
#            plt.plot(stokes)
#            plt.plot(fit)
#            plt.show()  
#        Res1.append(coeff[0])
#        
#        p0 = [np.median(stokes),0,1,0,1]
#        coeff, var_matrix = curve_fit(EXP2,xs,stokes,p0,maxfev=100000)
#        print(coeff)
#        fit = EXP2(np.array(xs), *coeff)
#        Res2.append(coeff[0])
        
#    Res = (np.array(Res1)+np.array(Res2))/2
        
    with open('result.txt', 'w') as f:
        for idx,elem in enumerate(Res):
            f.write(str(retarder[idx]) + '\t' + str(JD[idx]) + '\t' + str(Alpha[idx]) + '\t' + str(PlAng[idx]) + '\t' + str(elem) + '\n')
        
    plt.show()        
    
#    if Auto:
#        min_index, min_value = min(enumerate(np.std(SS,axis=0)), key=operator.itemgetter(1))
#        print('%.5f +- %.5f' % (np.median(SS[:,min_index]), np.std(SS[:,min_index])/np.sqrt(len(SS[:,min_index]))))
#    
#    with open('Best_' + filenames[0],'w') as f:
#        f.write('%.5f +- %.5f' % (np.median(SS[:,min_index]), np.std(SS[:,min_index])/np.sqrt(len(SS[:,min_index]))))
#        
#    print(SS[:,Rad])    


if __name__ == '__main__':
    
    
    # define command line arguments
    parser = argparse.ArgumentParser(description='manual target identification')
    parser.add_argument('-auto', action="store_true")
    parser.add_argument('-plot', action="store_true")
    parser.add_argument('-object', help='Name of the target for retrieving alpha and scaterring plane angle values',
                        default=False) 
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    
    filenames = args.images
    Auto = args.auto
    Plot = args.plot
    Target = args.object.replace('_',' ')


    Anal(filenames,Auto,Plot,Target)


    pass
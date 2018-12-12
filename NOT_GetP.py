#!/usr/bin/env python

""" CAPS_GetPola - Convert individual observation to q, u, P 
    v1.0: 2018-10-05, mdevogele@lowell.edu
"""

import matplotlib.pyplot as plt
import argparse
import numpy as np
import operator
from scipy.optimize import curve_fit
from astroquery.jplhorizons import Horizons
import NOT_Toolbox as tb


def Analyse(filenames):

    
    Retarder = []
    JD = []
    Alpha = []
    PlAng = []
    P = []
    
    with open(filenames[0],'r') as f:
        for elem in f.readlines():
            Retarder.append(float(elem.replace('\n',' ').replace('\t',' ').split()[0]))
            JD.append(float(elem.replace('\n',' ').replace('\t',' ').split()[1]))
            Alpha.append(float(elem.replace('\n',' ').replace('\t',' ').split()[2]))
            PlAng.append(float(elem.replace('\n',' ').replace('\t',' ').split()[3]))
            P.append(float(elem.replace('\n',' ').replace('\t',' ').split()[4]))            

    print(PlAng)
    Old_Ret = 999
    Same = False
    Retarder = np.array(Retarder)
    Q1=[]
    Q1_Ind=[]
    Q1_Ind.append(list(np.where(Retarder==0)[0]))
    Q1_Ind.append(list(np.where(Retarder==180)[0]))
    Q1_Ind = sum(Q1_Ind,[])
    
    Q2=[]
    Q2_Ind=[]
    Q2_Ind.append(list(np.where(Retarder==45)[0]))
    Q2_Ind.append(list(np.where(Retarder==225)[0]))
    Q2_Ind = sum(Q2_Ind,[])

    Q3=[]
    Q3_Ind=[]
    Q3_Ind.append(list(np.where(Retarder==90)[0]))
    Q3_Ind.append(list(np.where(Retarder==270)[0]))
    Q3_Ind = sum(Q3_Ind,[])


    Q4=[]
    Q4_Ind=[]
    Q4_Ind.append(list(np.where(Retarder==135)[0]))
    Q4_Ind.append(list(np.where(Retarder==315)[0]))
    Q4_Ind = sum(Q4_Ind,[])


    U1=[]
    U1_Ind=[]
    U1_Ind.append(list(np.where(Retarder==22.5)[0]))
    U1_Ind.append(list(np.where(Retarder==202)[0]))
    U1_Ind = sum(U1_Ind,[])
    
    U2=[]
    U2_Ind=[]
    U2_Ind.append(list(np.where(Retarder==67.5)[0]))
    U2_Ind.append(list(np.where(Retarder==247)[0]))
    U2_Ind = sum(U2_Ind,[])

    U3=[]
    U3_Ind=[]
    U3_Ind.append(list(np.where(Retarder==112.5)[0]))
    U3_Ind.append(list(np.where(Retarder==292)[0]))
    U3_Ind = sum(U3_Ind,[])


    U4=[]
    U4_Ind=[]
    U4_Ind.append(list(np.where(Retarder==157.5)[0]))
    U4_Ind.append(list(np.where(Retarder==337)[0]))
    U4_Ind = sum(U4_Ind,[])


    JD_F =[]
    A_F = []
    PS_F =[]
    P_F = [] 

    Series_Q1 = tb.Get_Consecutive(Q1_Ind)
    Q1 = []
    JDD_Q1 = []
    alpha_Q1 = [] 
    PSA_Q1 = [] 
    for elem in Series_Q1:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        Q1.append(np.mean(PP))
        JDD_Q1.append(np.mean(JJD))
        alpha_Q1.append(np.mean(AA))
        PSA_Q1.append(np.mean(PsAn))
 

    Series_U1 = tb.Get_Consecutive(U1_Ind)
    U1 = []
    JDD_U1 = []
    alpha_U1 = [] 
    PSA_U1 = [] 
    for elem in Series_U1:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        U1.append(np.mean(PP))
        JDD_U1.append(np.mean(JJD))
        alpha_U1.append(np.mean(AA))
        PSA_U1.append(np.mean(PsAn))
        
    P_F.append(np.sqrt(np.array(Q1)**2+np.array(U1)**2))
    JD_F.append(np.mean([JDD_Q1,JDD_U1],axis=0))
    A_F.append(np.mean([alpha_Q1,alpha_U1],axis=0))    
    PS_F.append(np.mean([PSA_Q1,PSA_U1],axis=0))    
           
    Series_Q2 = tb.Get_Consecutive(Q2_Ind)
    Q2 = []
    JDD_Q2 = []
    alpha_Q2 = [] 
    PSA_Q2 = [] 
    for elem in Series_Q2:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        Q2.append(np.mean(PP))
        JDD_Q2.append(np.mean(JJD))
        alpha_Q2.append(np.mean(AA))
        PSA_Q2.append(np.mean(PsAn))
 

    Series_U2 = tb.Get_Consecutive(U2_Ind)
    U2 = []
    JDD_U2 = []
    alpha_U2 = [] 
    PSA_U2 = [] 
    for elem in Series_U2:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        U2.append(np.mean(PP))
        JDD_U2.append(np.mean(JJD))
        alpha_U2.append(np.mean(AA))
        PSA_U2.append(np.mean(PsAn))
        
    P_F.append(np.sqrt(np.array(Q2)**2+np.array(U2)**2))
    JD_F.append(np.mean([JDD_Q2,JDD_U2],axis=0))
    A_F.append(np.mean([alpha_Q2,alpha_U2],axis=0))    
    PS_F.append(np.mean([PSA_Q2,PSA_U2],axis=0))  
    
    
    QQ1 = np.array(Q1)
    UU1 = np.array(U1)
    QQ2 = np.array(Q2)
    UU2 = np.array(U2)
    
    QQ = 1./2*(QQ1-QQ2)
    UU = 1./2*(UU1-UU2)
   
    P_Fin = []
    JD_Fin = []
    A_Fin = []
    PS_Fin = []
    Q_Fin = []
    U_Fin = []
    P_Fin.append(np.sqrt(QQ**2+UU**2))
    JD_Fin.append(np.mean(np.array(JD_F),axis=0))
    A_Fin.append(np.mean(np.array(A_F),axis=0))
    PS_Fin.append((np.mean(np.array(PS_F),axis=0)))
    Q_Fin.append(QQ)
    U_Fin.append(UU)    





    JD_F =[]
    A_F = []
    PS_F =[]
    P_F = [] 
    
    Series_Q3 = tb.Get_Consecutive(Q3_Ind)
    Q1 = []
    JDD_Q1 = []
    alpha_Q1 = [] 
    PSA_Q1 = [] 
    for elem in Series_Q3:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        Q1.append(np.mean(PP))
        JDD_Q1.append(np.mean(JJD))
        alpha_Q1.append(np.mean(AA))
        PSA_Q1.append(np.mean(PsAn))
 

    Series_U3 = tb.Get_Consecutive(U3_Ind)
    U1 = []
    JDD_U1 = []
    alpha_U1 = [] 
    PSA_U1 = [] 
    for elem in Series_U3:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        U1.append(np.mean(PP))
        JDD_U1.append(np.mean(JJD))
        alpha_U1.append(np.mean(AA))
        PSA_U1.append(np.mean(PsAn))
        
    P_F.append(np.sqrt(np.array(Q1)**2+np.array(U1)**2))
    JD_F.append(np.mean([JDD_Q1,JDD_U1],axis=0))
    A_F.append(np.mean([alpha_Q1,alpha_U1],axis=0))    
    PS_F.append(np.mean([PSA_Q1,PSA_U1],axis=0))  


    Series_Q4 = tb.Get_Consecutive(Q4_Ind)
    Q2 = []
    JDD_Q2 = []
    alpha_Q2 = [] 
    PSA_Q2 = [] 
    for elem in Series_Q4:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        Q2.append(np.mean(PP))
        JDD_Q2.append(np.mean(JJD))
        alpha_Q2.append(np.mean(AA))
        PSA_Q2.append(np.mean(PsAn))
 

    Series_U4 = tb.Get_Consecutive(U4_Ind)
    U2 = []
    JDD_U2 = []
    alpha_U2 = [] 
    PSA_U2 = [] 
    for elem in Series_U2:
        PP = []
        JJD = []
        AA = []
        PsAn = []
        for elem2 in elem:
            PP.append(P[elem2[1]])
            JJD.append(JD[elem2[1]])
            AA.append(Alpha[elem2[1]])
            PsAn.append(PlAng[elem2[1]])
        U2.append(np.mean(PP))
        JDD_U2.append(np.mean(JJD))
        alpha_U2.append(np.mean(AA))
        PSA_U2.append(np.mean(PsAn))
        
    P_F.append(np.sqrt(np.array(Q1)**2+np.array(U1)**2))
    JD_F.append(np.mean([JDD_Q1,JDD_U1],axis=0))
    A_F.append(np.mean([alpha_Q1,alpha_U1],axis=0))    
    PS_F.append(np.mean([PSA_Q1,PSA_U1],axis=0))  


    QQ1 = np.array(Q1)
    UU1 = np.array(U1)
    QQ2 = np.array(Q2)
    UU2 = np.array(U2)
    
    QQ = 1./2*(QQ1-QQ2)
    UU = 1./2*(UU1-UU2)
    
    P_Fin.append(np.sqrt(QQ**2+UU**2))
    JD_Fin.append(np.mean(np.array(JD_F),axis=0))
    A_Fin.append(np.mean(np.array(A_F),axis=0))
    PS_Fin.append((np.mean(np.array(PS_F),axis=0)))
    Q_Fin.append(QQ)
    U_Fin.append(UU)  


    P_Final = np.array(P_Fin).flatten()         
    JD_Final = np.array(JD_Fin).flatten()             
    A_Final = np.array(A_Fin).flatten()         
    PS_Final = np.array(PS_Fin).flatten() 
    Q_Final = np.array(Q_Fin).flatten() 
    U_Final = np.array(U_Fin).flatten()
    
    Angle = 90
    Pr =  (np.sin((2*(PS_Final+Angle))*np.pi/180)*(U_Final)+np.cos((2*(PS_Final+Angle))*np.pi/180)*(Q_Final));
    Pr2 = (-np.sin((2*(PS_Final+Angle))*np.pi/180)*(-Q_Final)+np.cos((2*(PS_Final+Angle))*np.pi/180)*(-U_Final))
    
    print(P_Final)
    f = open('Final_Result.txt','w')
    f.write('JD \t Alpha \t Q \t U \t P \t Pr \n')
    for idx,elem in enumerate(P_Final):        
        f.write(str(JD_Final[idx]) + '\t' + str(A_Final[idx]) + '\t' + str(Q_Final[idx]) + '\t' + str(U_Final[idx]) + '\t' + str(P_Final[idx]) + '\t' + str(Pr[idx]) + '\n')
    
    
    
    

if __name__ == '__main__':
    
    
    # define command line arguments
    parser = argparse.ArgumentParser(description='manual target identification')
    parser.add_argument('-auto', action="store_true")
    parser.add_argument('-plot', action="store_true")
    parser.add_argument('-object', help='Name of the target for retrieving alpha and scaterring plane angle values',
                        default=False)
    parser.add_argument('-retarder', help='Possition of the retarder information on the filename',
                        default=3)    
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    
    filenames = args.images
    Auto = args.auto
    Plot = args.plot
    Retarder = args.retarder


    Analyse(filenames)


    pass
#!/usr/bin/env python

""" CAPS_Phot - Perform the wavelength calibration
    v1.0: 2018-04-19, mdevogele@lowell.edu
"""

from __future__ import print_function
from __future__ import division

from past.utils import old_div
import os, sys
import numpy
from tkinter import *
from PIL import Image
from PIL import ImageTk
import argparse
from astropy.io import fits
from scipy.ndimage import interpolation as interp
from astropy.time import Time

from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from photutils import make_source_mask

from photutils import centroid_com, centroid_1dg, centroid_2dg
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

import time



class Clicker(object):

    def __init__(self, master, zoom, filelist):
        self.top = master
        self.files = filelist
        self.zoom = zoom
        self.target_index = [None for i in range(len(self.files))]
        self.interp_index = [None for i in range(len(self.files))]
        self.index = 0
        self.images = []
        self.ldac   = []
        self.mjd    = []
        self.retarder = []

        
        self.redcircle = []
        
        self.JD = []

        ### load image data
        print('please wait, loading images...', end=' ')
        sys.stdout.flush()

        self.read_all_fits(self.files)

        print('done!')

        # create title bar
        self.title = Label(text='%s (%d/%d)' %
                           (self.images[0],
                            self.index+1, len(self.files)))
        self.title.pack()


        # select first image
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(master, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)    

        # create position indicators:
        # green: sources, yellow: inter/extrapolated, red: manually
        # selected
        self.green_circles = []

        self.redcircle_id = self.canvas.create_oval(-100, -100, -100,
                                                    -100, outline='red',
                                                    width=2)
        self.yellowcircle_id = self.canvas.create_oval(-100, -100, -100,
                                                    -100, outline='orange',
                                                       width=1)

        # frame counter variable
        self.evar = IntVar()
        self.evar.set(1)

        self.canvas.pack(side='top', expand=1)

        # display image
        self.nextframe()

        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 3>", self.right_click)


    def left_click(self, event):
        """ select source """
        x, y = old_div(event.x,self.zoom), old_div(event.y,self.zoom)
        
        
        
        self.Select_Target(self.data[self.index],x,y,5)


    def key(self, event):
        """ keyboard events """
        if event.char == 'a':
            # previous frame
            self.nextframe(-1)
        elif event.char == 'd':
            # next frame
            self.nextframe(1)
        elif event.char == 'q':
            # quit
            self.top.quit()
        elif event.char == 'z':
            self.Auto_Detect()
        
        elif event.char == 'x':
            self.Show_Targets()
        elif event.char == 'p':
            self.Photometrie()
            print('Processing DONE')

    def right_click(self, event):
        """ next frame """
        self.nextframe(1)

    def read_all_fits(self, filenames, zoom=0.5):
        """ read in all image data, scale images """
        self.data=[]
        self.Tx = []
        self.Ty = []
        for idx, filename in enumerate(filenames):
            if idx > 0:
                print('\b\b\b\b%3d' % (idx+1), end=' ')
            else:
                print('%3d' % (idx+1), end=' ')
            sys.stdout.flush()

            self.Tx.append([1385.0146174630074, 1389.1632183587542, 1408.5932392652608, 1416.3528818011403])
            self.Ty.append([1385.0146174630074, 1389.1632183587542, 1408.5932392652608, 1416.3528818011403])
            ## read image data
            hdulist = fits.open(filename, ignore_missing_end=True)
            imgdat = hdulist[1].data
            Date = hdulist[0].header['DATE-AVG']
            t = Time(Date)
            self.JD.append(t.jd)
            self.data.append(imgdat)
            self.retarder.append(hdulist[0].header['FARETANG'])

            median = numpy.nanmedian(imgdat[50:150,50:150])
            std    = numpy.nanstd(imgdat[50:150,50:150])

            imgdat = old_div(numpy.clip(imgdat, median-0.5*std,
                                median+2*std),(old_div(2.5*std,256)))
            imgdat = imgdat - numpy.min(imgdat)
            imgdat = interp.zoom(imgdat, self.zoom)

            self.images.append(Image.fromarray(imgdat))
            
    def nextframe(self,i=1, imgnum=-1):
        """ display frame using iterator i"""

        if imgnum == -1:
            self.index += i
        else:
            self.index = imgnum - 1
        if self.index >= len(self.files):
            self.index = 0
        elif self.index < 0:
            self.index = len(self.files) - 1
        filename = self.files[self.index]
        if not os.path.exists(filename):
            print("Unable to find %s" % filename)
            self.top.quit()
        self.evar.set(self.index+1)

        self.title.configure(text='%s (%d/%d)' %
                           (os.path.basename(filename),
                            self.index+1, len(self.files)))

        im = self.images[self.index]

        self.tkimage.paste(im)
        
        #delete previous circles
        for elem in self.redcircle:
            self.canvas.delete(elem)
            
    def Select_Target(self,image,x,y,box):
        
        X = []
        Y = []
        
        #delete previous circles
        for elem in self.redcircle:
            self.canvas.delete(elem)
    
        D = image[int(y)-box:int(y)+box,int(x)-box:int(x)+box]
                
        xn,yn = centroid_com(D-numpy.median(D))
        xn = xn+int(x)-box
        yn = yn+int(y)-box

        print('Center of the target: x = {} y = {}'.format(xn,yn))

        
        X.append(xn)
        Y.append(yn)
        
        NX = X[0]-35
        NY = Y[0]
        D = image[int(NY)-box:int(NY)+box,int(NX)-box:int(NX)+box]
        xn,yn = centroid_com(D-numpy.median(D))
        xn = xn+int(NX)-box
        yn = yn+int(NY)-box        
        
        X.append(xn)
        Y.append(yn)
 
    
        
        for xe,ye in zip(X,Y):
            self.redcircle.append(self.canvas.create_oval(xe*self.zoom-10, ye*self.zoom-10, xe*self.zoom+10,
                                                ye*self.zoom+10, outline='blue',
                                                    width=2))        
            self.redcircle.append(self.canvas.create_oval(xe*self.zoom-1, ye*self.zoom-1, xe*self.zoom+1,
                                                ye*self.zoom+1, outline='red',
                                                    width=2))    
        self.Tx[self.index] = X 
        self.Ty[self.index] = Y

    def Auto_Detect(self):
        
        print('Dectection of the target')
        old_index = self.index
        for idx, elem in enumerate(self.data):
            self.nextframe(1)
            self.Select_Target(self.data[self.index],self.Tx[old_index][0],self.Ty[old_index][0],10)
            old_index = self.index
        print('Dectection done')
        
    def Show_Targets(self):
        self.nextframe(1)


        for xe,ye in zip(self.Tx[self.index],self.Ty[self.index]):
            self.redcircle.append(self.canvas.create_oval(xe*self.zoom-10, ye*self.zoom-10, xe*self.zoom+10,
                                                ye*self.zoom+10, outline='blue',
                                                    width=2))        
            self.redcircle.append(self.canvas.create_oval(xe*self.zoom-1, ye*self.zoom-1, xe*self.zoom+1,
                                                ye*self.zoom+1, outline='red',
                                                    width=2)) 
        

    def Photometrie(self):  
        q = []
        for elem in self.data:
            phot = []
            for x,y in zip(self.Tx[self.index],self.Ty[self.index]):
                radii = range(3,25)
                positions = [(x, y)]
                apertures = [CircularAperture(positions, r=r) for r in radii]
                annulus_apertures = CircularAnnulus(positions, r_in=12., r_out=20.)
                phot_table = aperture_photometry(self.data[self.index], apertures)
                phot_background = aperture_photometry(self.data[self.index], annulus_apertures)
                bkg_mean = (phot_background['aperture_sum'] / annulus_apertures.area())
                Int= []
                Area = []
                Int_bkg1 = []

                               
                for num in range(3,25):
                    Int.append(phot_table[0][num])
                    Area.append(apertures[num-3].area())
                    Int_bkg1.append(phot_table[0][num] - apertures[num-3].area()*bkg_mean)
    
            
                phot.append(Int_bkg1)

            filename = self.files[self.index]
            print('Processing file: ' + str(filename))
            phot = numpy.array(phot)
            qq = (phot[0]-phot[1])/(phot[0]+phot[1])
            with open(filename +'.txt', 'w') as f:
                for elem3 in qq:
                    for elem2 in elem3:
                        f.write(str(self.JD[self.index]) + '\t' + str(self.retarder[self.index]) + '\t' + str(elem2) + '\n')
            q.append(qq.flatten())
            self.nextframe(1)

        q = numpy.array(q)
        with open('q', 'w') as f:
            for elem in q:
                for elem2 in elem:
                    f.write(str(elem2) + '\n')

if __name__ == '__main__':
    
    
    # define command line arguments
    parser = argparse.ArgumentParser(description='manual target identification')
    parser.add_argument('-zoom', help='image zoom factor', default=3)
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    zoom = float(args.zoom)
    filenames = args.images

    root = Tk()
    app = Clicker(root, zoom, filenames)
    root.mainloop()

    pass



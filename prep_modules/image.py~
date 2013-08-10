#!/usr/bin/python
# -*- coding:  cp1251 -*-
# Module to convert nMgy to ADU and some other things

#http://www.sdss.org/dr5/algorithms/fluxcal.html
#http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html

import sys
import shutil
import pyfits
import numpy as np
import math
import os

# http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
def prep_ima(survey,gal,new_gal, gain):
	if survey=="SDSS_dr8":
		# This function re-writes pixel values given in NMgy to ADU
		hdulist = pyfits.open(gal, do_not_scale_image_data=True, mode='update')
		img = hdulist[0].data
		(dimy,dimx) = img.shape
		cimg = np.zeros(shape=(dimy,dimx))
		nrowc = dimy
		calib = hdulist[1].data
		for i in range(dimy):
			cimg[i] = calib 

		#img_new = img / cimg
		#print scidata2 # image in DN
		#simg = sdss_prep_sky(gal)
		#dn = img/cimg + simg # now we have a fits-image with ADU with unsubtracted SKY
		#dn_err = sqrt(dn/gain + darkVariance) # Error of the image in dn	

		dn = img/cimg
		shutil.copy(gal,new_gal) 
		hdulist1 = pyfits.open(new_gal, do_not_scale_image_data=True, mode='update')
		img_new = hdulist1[0].data
		for i in range(dimy):
			for k in range(dimx):
				img_new[i,k] = dn[i,k]  #new fits-file img_new with ADU	
		hdulist1.flush()
		os.remove(gal)


def add_to_header(gal,EXPTIME,m0,GAIN,NCOMBINE,pix2sec,ron):
        """ The function to update the fits header adding some keywords to it. Theses are required for Galfitif you do not have a weight image!"""
	hdulist = pyfits.open(gal, do_not_scale_image_data=True, mode='update')
	prihdr = hdulist[0].header

	prihdr.update('EXPTIME',EXPTIME)#,before='DATE')
	if 'GAIN' in prihdr: z=1
	else:
  		prihdr.update('GAIN',GAIN)
	if 'NCOMBINE' in prihdr: z=1
	else:
  		prihdr.update('NCOMBINE',NCOMBINE)
	if 'RDNOISE' in prihdr: z=1
	else:
  		prihdr.update('RDNOISE',ron)
	prihdr.update('m0',m0)
	hdulist.flush()

def pa(x):
        """ The function which will bring position angle 
         measured by sextrator in the range -90 and 90"""		
        if(float(x)>=0 and float(x)<=180.0): 
            pos_ang = float(x) - 90.0 #position angle
        if(float(x)<0 and float(x)>=-180.0):
            pos_ang = 90.0 - abs(float(x))  #position angle
        if(float(x)>180 and float(x)<=360.0):
            pos_ang = float(x) - 360.0 + 90.0 #position angle
        if(float(x)>=-360 and float(x)<-180.0):
            pos_ang = float(x) + 360.0 - 90.0 #position angle	
        return pos_ang

def pa_arr(x):
        """ The function which will bring position angle 
         measured by sextrator in the range -90 and 90"""
	pos_ang = []
	for k in range(len(x)):		
		if float(x[k])>=0 and float(x[k])<=180.0: 
			pos_ang.append(float(x[k]) - 90.0)
		if float(x[k])<0 and float(x[k])>=-180.0:
			pos_ang.append(90.0 - abs(float(x[k])))
		if float(x[k])>180 and float(x[k])<=360.0:
			pos_ang.append(float(x[k]) - 360.0 + 90.0)
		if float(x[k])>=-360 and float(x[k])<-180.0:
			pos_ang.append(float(x[k]) + 360.0 - 90.0)	
        return pos_ang	 

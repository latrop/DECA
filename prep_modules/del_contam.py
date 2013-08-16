#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import pyfits

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')


tmp_out = sys.stdout

def fill_cont():
	# Function to fill masked areas with model pixel values
	try:
		x,y = loadtxt('badpix.txt', usecols=[0,1], unpack=True)
		x = x - 1
		y = y - 1

		hdulist = pyfits.open('out.fits', do_not_scale_image_data=True)
		model = hdulist[2].data

		(ny,nx) = model.shape
		shutil.copy('galaxy_clean.fits','galaxy_wt_cont.fits') 
		hdulist1 = pyfits.open('galaxy_wt_cont.fits', do_not_scale_image_data=True, mode='update')
		gal = hdulist1[0].data

		import sky_est
		sky_est.imcopy_func('out.fits[2]','model.fits',1,1,nx,ny)
		sky_est.imcopy_func('out.fits[3]','resid_wt_cont.fits',1,1,nx,ny)
		hdulist2 = pyfits.open('resid_wt_cont.fits', do_not_scale_image_data=True, mode='update')
		resid = hdulist2[0].data

		'''
		for i in y:
			for k in x:
				gal[i,k] = model[i,k]
		'''
		for k in range(len(x)):
			gal[y[k],x[k]] = model[y[k],x[k]]
			resid[y[k],x[k]] = 0.
		
		
		hdulist1.flush()
		hdulist2.flush()

	except:
		print 'No badpix.txt file is found!'
		shutil.copy('galaxy_clean.fits','galaxy_wt_cont.fits') 
		import sky_est
		#sky_est.imcopy_func('out.fits[2]','model.fits',1,1,nx,ny)
		hdulist = pyfits.open('out.fits', do_not_scale_image_data=True)
		model = hdulist[2].data

		(ny,nx) = model.shape
		sky_est.imcopy_func('out.fits[3]','resid_wt_cont.fits',1,1,nx,ny)

	#sky_est.imarith_func('galaxy_wt_cont.fits','-','model.fits','resid_wt_cont.fits')

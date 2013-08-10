#!/usr/bin/python
from pyraf import iraf
from scipy import interpolate
import random
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
import fileinput
import pyfits
import re
#from kapteyn import wcs
from scipy import special
from scipy.odr.odrpack import *

tmp_out = sys.stdout


file_gal_clean = 'galaxy_clean.fits'
file_gal_clean1 = 'galaxy_clean1.fits'
file_badpix = 'badpix.txt'


def mask(galfit_out,xc,yc,a,b):
	shutil.copy(file_gal_clean,file_gal_clean1) 
	hdulist3 = pyfits.open(file_gal_clean1, do_not_scale_image_data=True, mode='update')
	img3 = hdulist3[0].data


	hdulist = pyfits.open(galfit_out)
	prihdr = hdulist[3].header
	scidata = hdulist[3].data
	dimx = prihdr['NAXIS1']
	dimy = prihdr['NAXIS2']
	
	'''
	I = []
	x = []
	y = []

	for i in range(dimy):
		for k in range(dimx):
			if (k-int(xc+0.5))**2/(a*a) + (i-int(yc+0.5))**2/(b*b)<=1.:
				I.append(scidata[i,k])
				x.append(k)
				y.append(i)
	'''
	Imed = np.median(img3)
	#Istd = np.std(img3)
	#print Imed,Istd

	for i in range(dimy):
		for k in range(dimx):
			if fabs(scidata[i,k])>3.*Imed:
				img3[i,k] = 0.

	hdulist3.flush()

	f = open(file_badpix, "w") 
	sys.stdout = f

	for i in range(dimy):
		for k in range(dimx):
			if img3[i,k]==0.:
				print "%i %i" % (k+1,i+1)
	#exit()

	sys.stdout = tmp_out
	f.close()	
		



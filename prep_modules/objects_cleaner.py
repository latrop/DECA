#!/usr/bin/python
# The script to clean the image from contaminents like stars
#
#from pyraf import iraf
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

def r_cut(scidata,xc,y,pix2sec,m0,rmin,rmax,plot=0):
	#*** The cut along the major axis (x-axis) ***
	# INPUT:
	# 	xc - pixel coordinate of the center [pix]
	# 	y - the cut at this point (y-axis) [pix]
	# 	rmin, rmax - borders where to investigate [pix]
	# 	plot = 1,2 if you want to plot this cut (DN- and mag- /arcsec^2)
	# OUTPUT:
	#	x - shifted scaled coordinates to the center [pix]
	# 	x1 - unshifted scaled coordinates [pix]
	#	Imajor - Surface Brightness [DN]

	ny,nx = scidata.shape
	Imajor = []
	x = []
	x1 = []
	for k in range(0, nx):
		if fabs(k-xc+0.5)<rmax and fabs(k-xc+0.5)>rmin:
			Imajor.append(scidata[y, k])
			x.append(k- xc+0.5)
			x1.append(k+0.5)

	if plot==1:
		# Plotting
		plt.plot(x, Imajor,'r-',color='black', lw=1)
		xlim(-nx/2,nx/2)
		#ylim(floor(min(Iminor)),ceil(max(Iminor)+np.mean(Iminor)/4.))
		plt.xlabel(r'r (pix)', fontsize=30)
		plt.ylabel(r'Intensity (DN)', fontsize=30)
	if plot==2:
		# Plotting
		plt.plot(x, -2.5*log10(Imajor)+m0 + 5.*log10(pix2sec),'r-',color='black', lw=1)
		xlim(-nx/2.,nx/2.)
		#ylim(floor(min(Iminor)),ceil(max(Iminor)+np.mean(Iminor)/4.))
		plt.xlabel(r'r (arcsec)', fontsize=30)
		plt.ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=30)
	return np.array(x),np.array(x1),np.array(Imajor)

def cleaner_stars(gal,xc,yc,pix2sec,m0,rmin,rmax):
	hdulist = pyfits.open(gal)
	scidata = hdulist[0].data
	ny,nx = scidata.shape

	y = yc
	x,x1,I = r_cut(scidata,xc,y,pix2sec,m0,rmin,rmax,plot=2)

cleaner_stars('3.fits',1024,744,0.396,26.,0.,800.)

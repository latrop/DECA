#!/usr/bin/python
#*** THE SCRIPT TO CREATE A PICTURE OF THE SHAPE OF THE BULGE (AS X-STRUCTURES ETC.) ***
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import pyfits
import re
tmp_out = sys.stdout

def median_filter(gal,gal_new):
	shutil.copy(gal,gal_new) 
	hdulist = pyfits.open(gal_new, do_not_scale_image_data=True, mode='update')
	scidata = hdulist[0].data
	ny,nx = scidata.shape
	R = 5

	f = open(r"image.txt", "w")
	sys.stdout = f
	print '%i\t%i\t%i' % (nx,ny,R)
	for k in range(ny):
		for i in range(nx):
			print '%.3f' % (scidata[k,i])
	sys.stdout = tmp_out
	f.close()

	subprocess.call("./median", shell=True)
	I_new = loadtxt('image_new.txt', usecols=[0], unpack=True)

	for k in range(ny):
		for i in range(nx):
			scidata[k,i] = I_new[nx*k+i]
	hdulist.flush()
		

#!/usr/bin/python
# -*- coding:  cp1251 -*-

#*** Common modules ***
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
import time
import pyfits
import re
from scipy import special
#from pyraf import iraf

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')


import setup

N_min = setup.N_min
coeff_disk_backgr = setup.coeff_disk_backgr
coeff_backgr = setup.coeff_backgr
# We are analyzing SB of the azimuthal profile or the major cut profile.

def Idisk(r,m0d,h,m0,pix2sec):
	return 10**(0.4*( m0 + 5.*log10(pix2sec) - m0d - 1.0857*r/h ))

def Idisk_edge(r,m0d,h,m0,pix2sec):
	I0d = 10**(0.4*(m0 + 5.*log10(pix2sec) - m0d)) 
	return I0d * fabs(r)/h * sp.special.kn(1,fabs(r)/h)

def rmax_find(r,I,h,m0d,m0,pix2sec,NOISE):
	n = 0
	Rmax = []


	'''
	plt.plot(r,m0 + 5.*log10(pix2sec) - 2.5*log10(Idisk(r,m0d,h,m0,pix2sec)))
	plt.plot(r,m0 + 5.*log10(pix2sec) - 2.5*log10(I))
	plt.show()

	exit()
	'''
	#print Idisk(r,m0d,h,m0,pix2sec)
	#print coeff_backgr*NOISE
	#exit()
	for k in range(len(r)):
		if r[k]>2.*h and (fabs(I[k]-Idisk(r[k],m0d,h,m0,pix2sec))>coeff_disk_backgr*NOISE or I[k]<coeff_backgr*NOISE and r[k]<6.*h):
			#print 'YES'
			n = n + 1
			if n == 1:	rmax = r[k]
			if n == N_min:	Rmax.append(rmax)
		else:
			n = 0
	#print Rmax
	RMAX = max(Rmax)
	#try:	RMAX = max(Rmax)
	#except:
	#print 'Could not find the true cutout of the disk. We will use RMAX=5h!' 
	#	RMAX = 5.*h
	print 'RMAX=%.3f' % (RMAX)	
	return RMAX

def rmax_find_edge(r,I,h,m0d,m0,pix2sec,NOISE):
	n = 0
	Rmax = []

	for k in range(len(r)):
		if r[k]>2.*h and (fabs(I[k]-Idisk_edge(r[k],m0d,h,m0,pix2sec))>coeff_disk_backgr*NOISE or I[k]<coeff_backgr*NOISE):
			#print 'YES'
			n = n + 1
			if n == 1:	rmax = r[k]
			if n == N_min:	Rmax.append(rmax)
		else:
			n = 0
	try:	RMAX = max(Rmax)
	except:
		print 'Could not find the true cutout of the disk. We will use RMAX=5h!' 
		RMAX = 4.*h	
	return RMAX



def rmax_d25(r,m0d,h,m0,pis2sec):	# only for g-band
	RMAX = 0.
	d25 = 25.0 - 0.4 
	for k in range(len(r)):
		if m0 + 5.*log10(pix2sec) - 2.5*log10(Idisk(r[k],m0d,h,m0,pix2sec))>=d25:
			RMAX = r[k]

	if RMAX == 0.:
		print 'Could not find the true cutout of the disk. We will use RMAX=5h!' 
		RMAX = 5.*h

	return RMAX	
		
	
	

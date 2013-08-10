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

#******** Functions for calculating m0,mags and errors ********#
# _dr7 - for dr<=7
# _dr8 - for dr>=8


def m0_dr7(exptime,aa,kk,airmass,pix2sec):
	# Calculated for arcsec^2!
	return 2.5*log10(exptime) - (aa+kk*airmass) #5.*log10(pix2sec)

def m0_dr8(pix2sec):
	# Calculated for arcsec^2!
	return 22.5 #+ 5.*log10(pix2sec)

def m0_dr8_dn(k):
	# Calculated for arcsec^2!
	return 22.5 - 2.5*log10(k)		

def mag_dr7(DN,exptime,aa,kk,airmass,pix2sec):
	# Calculated for arcsec^2!
	return -2.5*log10(DN) + 5.*log10(pix2sec) + 2.5*log10(exptime) - (aa+kk*airmass)

def DN_err_dr7(DN,skyErr,GAIN,darkVariance,pix2sec):
	# Calculated for arcsec^2!
	# DN = DN(sky) + DN(object) 
	return sqrt( DN/((pix2sec**2.)*GAIN) + 1./(pix2sec**2.) *(darkVariance + skyErr) )

def mag_err_dr7(DN,skyErr,GAIN,darkVariance):
	# Calculated for arcsec^2!
	return 2.5/ln(10.) * DN_err_dr7(DN,skyErr,GAIN,darkVariance,pix2sec)/DN

def mag_dr8(DN,pix2sec):
	return 22.5 - 2.5*log10(DN) + 5.*log10(pix2sec)

def DN_err_dr8(DN, GAIN,darkVariance,pix2sec):
	# Calculated for arcsec^2!
	# DN = DN(sky) + DN(object) 
	return sqrt( DN/((pix2sec**2.)*GAIN) + 1./(pix2sec**2.)*darkVariance )	

def mag_err_dr8(DN,GAIN,darkVariance,pix2sec):
	# Calculated for arcsec^2!
	return 2.5/ln(10.) * DN_err_dr8(DN, GAIN,darkVariance,pix2sec)/DN

# For UKIDSS
def header_extr(gal_image):
	hdulist = pyfits.open(gal_image)#, do_not_scale_image_data=True, mode='update')
	prihdr = hdulist[1].header
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']
	pix2sec = prihdr['PIXLSIZE']
	GAIN = prihdr['GAIN']
	read_out_noise = prihdr['READNOIS']
	sky_level = prihdr['SKYLEVEL']
	sky_noise = prihdr['SKYNOISE']
	fwhm = prihdr['SEEING'] * pix2sec # now in arcsec
	m0 = prihdr['MAGZPT']
	A = prihdr['EXTINCT']
	prihdr0 = hdulist[0].header
	exptime = prihdr0['EXP_TIME']
	NCOMBINE = prihdr0['NEXP']
	#del prihdr['CTYPE1']
	#del prihdr['CTYPE2']
	#hdulist.flush()
	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE

#!/usr/bin/python
# -*- coding:  cp1251 -*-
# SCRIPT TO ROTATE CUT OUT THE IMAGE AND CHANGE UNCORRECT COORDINATES

import random as random_number
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

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

import sky_est
import setup

extr_coeff = setup.extr_coeff

a_coeff = 15. # coefficient for extraction along the x-axis
b_coeff = 15. # coefficient for extraction along the y-axis

def run_imstat(input):
 iraf.images(_doprint=0)
 imtype = iraf.images.imutil
 all_exist = 1 # boolean flag
 for image in input:
  if not iraf.imaccess(image):
   all_exist = 0
   print "Error: can’t open", image
 if not all_exist:
  return
 iraf.imstat(",".join(input))


#run_imstat(["field.fits[0]"])

#iraf.images.imgeom.rotate("field.fits[0]","galaxy.fits",rotation=-32.4,xin=1029.427, yin=929.542)
#iraf.imcopy("galaxy.fits[690:1260,630:880]","galaxy1.fits")

def extr_box(xc,yc,a,b,theta,rmin,nx,ny):
	theta=radians(-theta)
	t=arange(-math.pi,math.pi,math.pi/100.)
	x = xc + a*cos(theta)*cos(t) + b*sin(theta)*sin(t)
	y = yc - a*sin(theta)*cos(t) + b*cos(theta)*sin(t)
	xmin = min(x)
	ymin = min(y)
	xmax = max(x)
	ymax = max(y)
	if xmin>xmax:
		xmin1 = xmin
		xmin = xmax
		xmax = xmin1
	if ymin>ymax:
		ymin1 = ymin
		ymin = ymax
		ymax = ymin1
	'''
	if rmin<sqrt((xmax-xmin)**2+(ymax-ymin)**2)/2.:
		xmin = xc - rmin/sqrt(2.)*a/b
		xmax = xc + rmin/sqrt(2.)*a/b
		ymin = yc - rmin/sqrt(2.)
		ymax = yc + rmin/sqrt(2.)
	'''


	if xc>=nx/2:
		if xc+a*fabs(cos(theta))>=nx:
			print 'here 1'
			xmax = nx-2
			xmin = xc-a*fabs(cos(theta))
			if xmin<=0:
				xmin = 2
	if xc<nx/2:
		if xc<=a*fabs(cos(theta)):
			print 'here 2'
			xmin = 2
			xmax = xc+a*fabs(cos(theta))
			if xmax>=nx:
				xmax = nx-2
	if yc>=ny/2:
		if yc+a*fabs(sin(theta))>=ny:
			print 'here 3'
			ymax = ny-2
			ymin = yc-a*fabs(sin(theta))
			if ymin<=0:
				ymin = 2
	if yc<ny/2:
		if yc<=a*fabs(sin(theta)):
			print 'here 4'
			ymin = 2
			ymax = yc + a*fabs(sin(theta)) 
			if ymax>=ny:
				ymax = ny - 2
	#exit()
	return int(xmin),int(ymin),int(xmax),int(ymax)

def rotate(survey,file_in,file_out,angle,xc,yc,nx,ny):
	if survey == "SDSS_dr8":
		file_in = file_in+'[0]' # for cube data file (SDSS frames)
	if survey == "UKIDSS":
		file_in = file_in+'[1]' # for cube data file (UKIDSS frames)
	if os.path.exists(file_out)==True:
		os.remove(file_out)

	#iraf.images.imgeom.rotate(file_in,file_out,rotation=-angle,xin=xc, yin=yc,xout=nx/2,yout=ny/2)
	sky_est.rotate_func(file_in,file_out,-angle,xc,yc,nx/2.,ny/2.)
'''
def extract(survey,file_in,file_out,a_image,b_image):
	hdulist = pyfits.open(file_in)
	prihdr = hdulist[0].header
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']
	xc = int(nx/2)
	yc = int(ny/2)
	dimx2 = a_image*a_coeff 
	dimy2 = b_image*b_coeff 
	x1 = int(xc - dimx2)
	y1 = int(yc - dimy2)
	x2 = int(xc + dimx2)
	y2 = int(yc + dimy2)
	if survey == "SDSS":
		file_in = file_in #+'[0]' # for cube data file (SDSS frames)
	if os.path.exists(file_out)==True:
		os.remove(file_out)
	#file_in = file_in + '['+str(x1)+':'+str(x2)+','+str(y1)+':'+str(y2)+']'
	#iraf.imcopy(file_in,file_out)
	sky_est.imcopy_func(file_in,file_out,x1,y1,x2,y2)
'''
def extract1(survey,file_in,file_out,xc,yc,a,b,theta):

	hdulist = pyfits.open(file_in)
	prihdr = hdulist[0].header
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']
	rmin = min([xc,nx-xc,ny,ny-yc])
	xmin,ymin,xmax,ymax = extr_box(xc,yc,extr_coeff*a,extr_coeff*b,theta,rmin,nx,ny)
	if xmin<0:	xmin=1
	if ymin<0:	ymin=1
	if xmax>nx:	xmax=nx
	if ymax>ny:	ymax=ny
	#if survey == "SDSS_dr8":
	#	file_in = file_in +'[0]' # for cube data file (SDSS frames)
	if os.path.exists(file_out)==True:
		os.remove(file_out)
	'''
	from pyraf import iraf
	file_in = file_in + '['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+']'
	print file_in
	iraf.imcopy(file_in,file_out)
	'''
	sky_est.imcopy_func(file_in,file_out,xmin,ymin,xmax,ymax)

def coord_convert(image,RA,DEC):
	#http://www.astro.rug.nl/software/kapteyn/wcs.html
	hdulist = pyfits.open(image)
	proj3 = wcs.Projection(hdulist[0].header) # create Projection object
	world = (RA, DEC)
	pixel = proj3.topixel(world)
	return pixel[0],pixel[1]





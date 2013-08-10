#!/usr/bin/python
# Module to find sky background
# Also there are some useful functions to operate with images: rotation, extraction, arithmetical expressions.
#!/usr/bin/python
# -*- coding:  cp1251 -*-
# RMS(signal) = Stdev(signal) if the mean signal is 0


import random as random_number
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
import fileinput
import pyfits
import re
tmp_out = sys.stdout

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

import sextr

import setup
star_detect = setup.star_detect
#way = setup.sky_find 


def imarith_func(file1,oper,file2,file_out):
	# Function iraf.imarith
	f_res = open("imari.cl", "w") 
	sys.stdout = f_res
	print "# Script: imarith"
	print "images"
	print "imutil"
	if type(file2)==float:
		print "imarith %s %s %f %s" % (file1,oper,file2,file_out)
	else:
		print "imarith %s %s %s %s" % (file1,oper,file2,file_out)
	print "logout"	
	sys.stdout = tmp_out
	f_res.close()
	os.chmod(r"imari.cl",0777)
	subprocess.call("cl < imari.cl -o", shell=True)
	os.remove('imari.cl')

def imsurfit_func(file1,file2):
	# Function iraf.imsurfit
	f_res = open("imsur.cl", "w") 
	sys.stdout = f_res
	print "# Script: imsurfit"
	print "images"
	print "imfit"
	print "imsurfit %s %s 3 3 function=leg" % (file1,file2)
	print "logout"	
	sys.stdout = tmp_out
	f_res.close()
	os.chmod(r"imsur.cl",0777)
	subprocess.call("cl < imsur.cl -o", shell=True)
	os.remove('imsur.cl')



def rotate_func(file_in,file_out,angle,xc,yc,x_out,y_out):
	# Function imgeom.rotate
	f_res = open("rot.cl", "w") 
	sys.stdout = f_res
	print "# Script: rotate"
	print "images"
	print "imgeom"
	print "rotate %s %s %.2f xin=%.2f yin=%.2f xout=%.2f yout=%.2f interp=poly3 " % (file_in,file_out,angle,xc,yc,x_out,y_out)
	print "logout"	
	sys.stdout = tmp_out
	f_res.close()
	os.chmod(r"rot.cl",0777)
	subprocess.call("cl < rot.cl -o", shell=True)
	os.remove('rot.cl')


def imcopy_func(file_in,file_out,xmin,ymin,xmax,ymax):
	# Function iraf.imcopy
	f_res = open("imcop.cl", "w") 
	sys.stdout = f_res
	print "# Script: imcopy"
	print "images"
	print "imutil"
	print "imcopy %s[%i:%i,%i:%i] %s" % (file_in,xmin,xmax,ymin,ymax,file_out)
	print "logout"	
	sys.stdout = tmp_out
	f_res.close()
	os.chmod(r"imcop.cl",0777)
	subprocess.call("cl < imcop.cl -o", shell=True)
	os.remove('imcop.cl')

def sky1():
	# From the SExtractor output file: CLASSTAR > 0.
	backgr,class_star= loadtxt('field.cat', usecols=[16,17], unpack=True, skiprows = 17)
	sky = []
	for k in range(len(class_star)):
		if class_star[k]>star_detect:
			sky.append(backgr[k])
	sky_level = np.median(sky)
	return sky_level

def Variance():
	hdulist = pyfits.open('fon.fits') # open a FITS file
	scidata = hdulist[0].data	
	return (np.std(scidata))**2

def sky2():
	#iraf.imarith('field.fits','-','objects.fits','backgr1.fits')
	imarith_func('field.fits','-','objects.fits','backgr1.fits')
	hdulist = pyfits.open('backgr1.fits') # open a FITS file
	scidata = hdulist[0].data
	sky_level = np.median(scidata)
	#iraf.imarith('backgr1.fits','-',sky_level,'backgr1.fits')
	imarith_func('backgr1.fits','-',sky_level,'backgr1.fits')
	hdulist1 = pyfits.open('backgr1.fits') # open a FITS file
	scidata1 = hdulist1[0].data
	variance = (np.std(scidata1))**2
	os.remove(r"backgr1.fits") 
	return sky_level,variance

def sky3():
	hdulist = pyfits.open('field.fits') # open a FITS file
	scidata = hdulist[0].data
	ny,nx = scidata.shape	
	hdulist1 = pyfits.open('segm.fits') # open a FITS file
	scidata1 = hdulist1[0].data
	
	I = []
	for i in range(0,ny):
		for k in range(0,nx):
			if scidata1[i,k]==0:
				I.append(scidata[i,k])
	sky_level = np.median(I)
	return sky_level,(np.std(I-sky_level))**2

def sky4(survey):
	#http://www.ast.cam.ac.uk/~rgm/scratch/e.eglez/iraf/imsurfit.txt
	#iraf.imarith('field.fits','-','objects.fits','backgr1.fits')
	if survey == 'UKIDSS':
		imarith_func('field.fits[1]','-','objects.fits','backgr1.fits')	
	elif survey == 'SDSS_dr8':
		imarith_func('field.fits[0]','-','objects.fits','backgr1.fits')	
	else:
		imarith_func('field.fits','-','objects.fits','backgr1.fits')


	#iraf.imsurfit('backgr1.fits','backgr2.fits',xorder="3",yorder="3",function="leg")
	imsurfit_func('backgr1.fits','backgr2.fits')
	hdulist = pyfits.open('backgr2.fits') # open a FITS file
	I = hdulist[0].data	
	sky_level = np.median(I)
	#iraf.imarith('backgr1.fits','-',sky_level,'backgr1.fits')
	imarith_func('backgr1.fits','-',sky_level,'backgr1.fits')
	hdulist1 = pyfits.open('backgr1.fits') # open a FITS file
	scidata1 = hdulist1[0].data
	variance = (np.std(scidata1))**2
	os.remove(r"backgr1.fits") 
	os.remove(r"backgr2.fits") 
	return sky_level,variance


def sky_aper(survey,xc,yc,kron_r,ellip):
	#iraf.imarith('field.fits','-','objects.fits','backgr1.fits')
	sextr.gal_clean_ima1(survey,'segm.fits','field.fits','backgr1.fits',xc,yc)


	step = 1.0
	#print xc,yc
	#exit()

	f = open("ell_sky.cl", "w") 
	sys.stdout = f
	print "stsdas"
	print "stsdas"
	print "analysis"
	print "isophote"
	print "geompar.linear=yes"
	print "geompar.step=%.1f" % (step)
	print "geompar.minsma=%.1f" % (0.)
	print "geompar.maxsma=%.1f" % (2.5*kron_r)
	print "geompar.x0=%.3f" % (xc)
	print "geompar.y0=%.3f" % (yc)
	print "geompar.pa=%.3f" % (90.)
	#print "geompar.ellip=%.3f" % (ellip)
	print "ellipse backgr1.fits sky.tab"
	print "tprint sky.tab pwidth=600 plength=1000 > ell_sky.txt"
	print "!rm -f sky.tab"
	print "logout"




	sys.stdout = tmp_out
	f.close()

	try:
		os.chmod(r"ell_sky.cl",0777)
		subprocess.call("cl < ell_sky.cl -o", shell=True)
	

		sma,inten,var = loadtxt('ell_sky.txt', usecols=[1,2,4], unpack=True, skiprows = 6, dtype='str')
		#exit()
		for k in range(len(sma)):
			if sma[k]=='INDEF': sma[k]=0.
			if inten[k]=='INDEF': inten[k]=0.
			if var[k]=='INDEF': var[k]=0.

		sma = np.array(sma,dtype='float')
		inten = np.array(inten,dtype='float')
		var = np.array(var,dtype='float')

		I_med = []
		variance = []
		for k in range(len(sma)):
			if sma[k]>2.*kron_r:
				I_med.append(inten[k])
				variance.append(var[k])
	
		os.remove("ell_sky.cl")
		os.remove("ell_sky.txt")
		os.remove("backgr1.fits")
		print I_med
		return mean(I_med),mean((np.array(variance))**2)
	except:
		return 99999.,99999.


def sky_random(survey,kron_r,xc_gal,yc_gal):
	# THE MAIN FUNCTION TO FIND RELIABLE SKY BACKGROUND
	from random import randint

	# See Martin_Navarro et al. (2012), http://arxiv.org/pdf/1208.2893v2.pdf
	if survey=='UKIDSS':
		imarith_func('field.fits[1]','-','objects.fits','backgr1.fits')
	elif survey=='SDSS_dr8':
		imarith_func('field.fits[1]','-','objects.fits','backgr1.fits')
	else:
		imarith_func('field.fits','-','objects.fits','backgr1.fits')

	hdulist = pyfits.open('backgr1.fits') # open a FITS file
	scidata = hdulist[0].data	
	ny,nx = np.shape(scidata)
	I_med = []
	I_std = []
	ii = 0
	while ii<2:
		for k in range(10000):
			I = []
			repeat = 0
			while repeat == 0:
				xc = randint(10,nx-11)
				yc = randint(10,ny-11)
				if sqrt((xc-xc_gal)**2 + (yc-yc_gal)**2)>2.*kron_r:
					repeat = 1
				else:	repeat = 0
			x = arange(xc-10,xc+10,1)
			y = arange(yc-10,yc+10,1)
			for m in y:
				for n in x:
					I.append(scidata[m,n])

			
			I_med.append(np.mean(I))	#I_med.append(np.median(I))
			I_std.append(np.std(I))
			del I

		ii = ii + 1
	os.remove("backgr1.fits")
	I_med1 = []
	I_std1 = []
	#print median(I_med),(median(I_std))**2

	#print len(I_med)
	med = mean(I_med)	#med = median(I_med)
	stdev = std(I_med)
	med_std = mean(I_std)	#med_std = median(I_std)
	stdev_std = std(I_std)
	#print stdev
	#exit()

	for k in range(len(I_med)):
		if fabs(I_med[k]-med)<stdev and fabs(I_std[k]-med_std)<stdev_std:
			I_med1.append(I_med[k])
			I_std1.append(I_std[k])

	return mean(I_med1),(mean(I_std1))**2
	#return median(I_med1),(median(I_std1))**2


def sky_box(survey,sma,smb,xc_gal,yc_gal,theta):
	from random import randint
	import changer

	# As in Martin_Navarro et al. (2012)
	if survey=='UKIDSS':
		#iraf.imarith('field[1].fits','-','objects.fits','backgr1.fits')
		imarith_func('field.fits[1]','-','objects.fits','backgr1.fits')

	elif survey=='SDSS_dr8':
		#iraf.imarith('field[0].fits','-','objects.fits','backgr1.fits')
		imarith_func('field.fits[0]','-','objects.fits','backgr1.fits')
	else:
		#iraf.imarith('field.fits','-','objects.fits','backgr1.fits')
		imarith_func('field.fits','-','objects.fits','backgr1.fits')


	hdulist = pyfits.open('backgr1.fits') # open a FITS file
	scidata = hdulist[0].data	
	ny,nx = np.shape(scidata)
	I_med = []
	I_std = []


	prihdr = hdulist[0].header
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']
	rmin = min([xc_gal,nx-xc_gal,ny,ny-yc_gal])
	xmin,ymin,xmax,ymax = changer.extr_box(xc_gal,yc_gal,2.5*sma,2.5*smb,theta,rmin,nx,ny)
	print xmin,ymin,xmax,ymax
	#print theta
	#exit()
	

	step = 1
	N = 50
	

	for k in range(N):
		x1 = []
		I = []
		y1 = arange(ymin,ymax,1)
		for k in range(len(y1)):
			if y1[k]<ny:
				x1.append(xmin)
				I.append(scidata[y1[k],xmin])

		x2 = arange(xmin,xmax,1)
		y2 = []
		for k in range(len(x2)):
			if x2[k]<nx:
				y2.append(ymax)	
				I.append(scidata[ymax,x2[k]])

		x3 = []
		y3 = arange(ymax,ymin,-1)
		for k in range(len(y3)):
			if y3[k]<ny:
				x3.append(xmax)
				I.append(scidata[y3[k],xmax])

		x4 = arange(xmax,xmin,-1)
		y4 = []
		for k in range(len(x4)):
			if x4[k]<nx:
				y4.append(ymin)
				I.append(scidata[ymin,x4[k]])
		I_med.append(median(I))
		I_std.append(std(I))
		'''
		plt.plot(x1,y1,'o',markersize=1.)
		plt.plot(x2,y2,'o',markersize=1.)
		plt.plot(x3,y3,'o',markersize=1.)
		plt.plot(x4,y4,'o',markersize=1.)
		'''

		xmin = xmin - step; ymin = ymin - step; xmax = xmax + step; ymax = ymax + step
	#plt.show()
	os.remove("backgr1.fits")
	return median(I_med),(median(I_std))**2

		
	



def mu_crit(r,I,sigma,m0):
	I1 = I - sigma
	I2 = I + sigma
	for k in range(len(r)):
		if fabs(2.5*log10(I2)-2.5*log10(I1))<0.2:
			mu = -2.5*log10(I[k]) + m0
			rmax = r[k]
	return mu,rmax



def sky_estim(survey,xc,yc,kron_r,ellip,theta,way):
	#print way
	#exit()
	#global way
	hdulist = pyfits.open('field.fits') # open a FITS file
	if survey=='UKIDSS':
		scidata = hdulist[1].data
		ny,nx = scidata.shape
		cube=0
	else:
		scidata = hdulist[0].data
		ny,nx = scidata.shape	
		cube=1		
	if nx<2.*kron_r or ny<2.*kron_r:
		way=4
	if way==1:
		sky_level = sky1()
		variance = Variance()
	if way==2:
		sky_level,variance = sky2()
	if way==3:
		sky_level,variance = sky3()
	if way==4:
		sky_level,variance = sky4(survey)
	if way==5:
		sky_level,variance = sky_random(survey,kron_r,xc,yc)
	if way==6:
		sky_level,variance = sky_aper(survey,xc,yc,kron_r,ellip)
		#if sky_level==99999. or variance==99999.:
		#	way=7
	if way==7:
		sky_level,variance = sky_box(survey,kron_r,kron_r*(1.-ellip),xc,yc,theta)

	#print 'sky_level= %.3f' % (sky_level)
	#print 'variance= %.3f' % (variance)
	return sky_level,variance




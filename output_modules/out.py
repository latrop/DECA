#!/usr/bin/python
# Script to plot profiles of edge-on galaxy
import random as random_number
import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
from os.path import exists
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

import ell

# SETUP FOR OUTPUT
fsize=27	# font size
fig_size = (10,10)
m_crit_add = 0.0
levels_step = 0.5
sample_pts = 500




def isophotes(number,file_input,xc,yc,m0,pix2sec,m_crit,rmax,zmax,h=0.,z0=0.,re=0.):
	m_crit = m_crit + m_crit_add
	m0 = m0 + 5.*log10(pix2sec)
	rmax = rmax*pix2sec
	zmax = zmax*pix2sec
	z0 = z0*pix2sec
	h = h*pix2sec
	re = re*pix2sec

	# ISOPHOTES map
	hdulist = pyfits.open(file_input)
	scidata = hdulist[0].data
	ny,nx = np.shape(scidata)

	f = plt.figure(0,figsize=(fig_size[0],fig_size[0]*zmax/rmax))
	ax1 = subplot(111)
	x = []
	y = []
	I = []
	I_crit = 10**(0.4*(m0-m_crit))

	for k in range(ny):
		for i in range(nx):
			if scidata[k,i]>I_crit:
				x.append((i + 0.5-xc)*pix2sec)
				y.append((k + 0.5-yc)*pix2sec)
				I.append(scidata[k,i])
	x = np.array(x)
	y = np.array(y)
	I = np.array(I)

	m = -2.5*log10(I) + m0 
	m = np.around(m,decimals=1)

	con_levels = np.arange(min(m),m_crit,levels_step)

	xi = np.linspace(-rmax,rmax,sample_pts)
	yi = np.linspace(-zmax,zmax,sample_pts)

	Ii = griddata(x,y,m,xi,yi)

	ax1.contour(xi,yi,Ii,con_levels,linewidths=1,colors='black')
	ax1.set_xlim(rmax,-rmax,-1)
	ax1.set_ylim(-zmax,zmax)
	if z0==0.:
		ax1.set_ylabel(r'y (arcsec)', fontsize=fsize)
	else:
		ax1.set_ylabel(r'z (arcsec)', fontsize=fsize)		
	ax1.set_xlabel(r'r (arcsec)', fontsize=fsize)

	if z0>0.:
		ax12 = ax1.twinx()
		ax12.set_ylabel(r'z (z$_{0}$)', fontsize=fsize)
		ax12.set_ylim(-rmax/z0,rmax/z0)
		ax12.set_yticks([0,2,4,6])
		ax12 = ax1.twiny()
		ax12.set_xlabel(r'r (h)', fontsize=fsize)
		ax12.set_xlim(-rmax/h,rmax/h)
	elif h>0. and z0==0.:
		ax12 = ax1.twinx()
		ax12.set_ylabel(r'y (h)', fontsize=fsize)
		ax12.set_ylim(-rmax/h,rmax/h)
		#ax12.set_yticks([0,2,4,6])
		ax12 = ax1.twiny()
		ax12.set_xlabel(r'r (h)', fontsize=fsize)
		ax12.set_xlim(-rmax/h,rmax/h)

	elif re>0. and z0==0.:
		ax12 = ax1.twinx()
		ax12.set_ylabel(r'y (r$_{e}$)', fontsize=fsize)
		ax12.set_ylim(-rmax/re,rmax/re)
		#ax12.set_yticks([0,2,4,6])
		ax12 = ax1.twiny()
		ax12.set_xlabel(r'r (r$_{e}$)', fontsize=fsize)
		ax12.set_xlim(-rmax/re,rmax/re)
	plt.savefig("isophotes.png", transparent = False)
	shutil.move('isophotes.png','./pics/%i/isophotes.png' % (number))


def crea_ell(file_in,xc,yc,step,minsma,maxsma):
	f = open("ell.cl", "w") 
	sys.stdout = f

	print "# Script: ellipse"
	print "!rm -f ellipse.txt"
	print "stsdas"
	#print "fitsio"
	#print "set imtype=hhh"
	#print "strfits %s \"\" model.hhh xdim=yes oldiraf=yes force=yes" % (file_in)
	#print "bye"
	print "analysis"
	print "isophote"
	print "geompar.step=%.1f" % (step) 
	print "geompar.minsma=%.1f" % (minsma)
	print "geompar.maxsma=%.1f" % (maxsma)
	print "geompar.x0=%.3f" % (xc)
	print "geompar.y0=%.3f" % (yc)
	print "ellipse %s model.tab" % (file_in)
	print "tprint model.tab pwidth=600 plength=1000 > ellipse.txt"
	print "!rm -f model.tab"
	print "logout"

	sys.stdout = tmp_out
	f.close()


def profiles_galfit_extract(xc_ima,yc_ima,theta,model,z0=0.,h=0.):
	# Along the major axis:
	# FITS Cube reading:

	if os.path.exists('profileima.fits')==True:	os.remove(r"profileima.fits")
	if os.path.exists('zima.fits')==True:	os.remove(r"zima.fits")
	if os.path.exists('profilemodel.fits')==True:	os.remove(r"profilemodel.fits")
	if os.path.exists('profiledisk.fits')==True:	os.remove(r"profiledisk.fits")
	if os.path.exists('profilebulge.fits')==True:	os.remove(r"profilebulge.fits")
	if os.path.exists('profileresid.fits')==True:	os.remove(r"profileresid.fits")
	if os.path.exists('zmodel.fits')==True:	os.remove(r"zmodel.fits")

	if os.path.exists('disk.fits')==True:	os.remove(r"disk.fits")
	if os.path.exists('bulge.fits')==True:	os.remove(r"bulge.fits")
	if os.path.exists('model.fits')==True:	os.remove(r"model.fits")
	if os.path.exists('resid.fits')==True:	os.remove(r"resid.fits")

	if os.path.exists('disk.imh')==True:	os.remove(r"disk.imh")
	if os.path.exists('bulge.imh')==True:	os.remove(r"bulge.imh")
	if os.path.exists('model.imh')==True:	os.remove(r"model.imh")
	if os.path.exists('resid.imh')==True:	os.remove(r"resid.imh")

	if theta<0.:	theta = 360+theta
	print 'THETA=%.3f' % (theta)

	f = open("profiles.cl", "w") 
	sys.stdout = f
	print "# Script: profiles"
	print "set imtype = \"imh\""
	#print "rfits galfit_out.fits 1 galaxy.imh"
	#print "wfits galaxy.imh galaxy_galfit.fits"
	if model=='edge+ser':
		print "rfits subcomps.fits 2 disk.imh"
		print "wfits disk.imh disk.fits"
		print "rfits subcomps.fits 3 bulge.imh"
		print "wfits bulge.imh bulge.fits"
		print "rfits galfit_out.fits 2 model.imh"
		print "wfits model.imh model.fits"
		print "rfits galfit_out.fits 3 resid.imh"
		print "wfits resid.imh resid.fits"
		print "set imtype = \"fits\""
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"zima.fits\" out_type=\"image\"" % (xc_ima+1.5*h,yc_ima,theta+90.)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"profilemodel.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"zmodel.fits\" out_type=\"image\"" % (xc_ima+1.5*h,yc_ima,theta+90.)
		print "pvector disk.fits xc=%i yc=%i theta=%.1f vec_output= \"profiledisk.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta+90.)
		print "pvector bulge.fits xc=%i yc=%i theta=%.1f vec_output= \"profilebulge.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector resid.fits xc=%i yc=%i theta=%.1f vec_output= \"profileresid.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "logout"

	if model=='edge':
		print "rfits subcomps.fits 2 disk.imh"
		print "wfits disk.imh disk.fits"
		print "rfits galfit_out.fits 2 model.imh"
		print "wfits model.imh model.fits"
		print "rfits galfit_out.fits 3 resid.imh"
		print "wfits resid.imh resid.fits"
		print "set imtype = \"fits\""
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"zima.fits\" out_type=\"image\"" % (xc_ima+1.5*h,yc_ima,theta+90.)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"profilemodel.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"zmodel.fits\" out_type=\"image\"" % (xc_ima+1.5*h,yc_ima,theta+90.)
		print "pvector disk.fits xc=%i yc=%i theta=%.1f vec_output= \"profiledisk.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta+90.)
		print "pvector resid.fits xc=%i yc=%i theta=%.1f vec_output= \"profileresid.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "logout"

	if model=='exp+ser':
		print "rfits subcomps.fits 2 disk.imh"
		print "wfits disk.imh disk.fits"
		print "rfits subcomps.fits 3 bulge.imh"
		print "wfits bulge.imh bulge.fits"
		print "rfits galfit_out.fits 2 model.imh"
		print "wfits model.imh model.fits"
		print "rfits galfit_out.fits 3 resid.imh"
		print "wfits resid.imh resid.fits"
		print "set imtype = \"fits\""
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"profilemodel.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector disk.fits xc=%i yc=%i theta=%.1f vec_output= \"profiledisk.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector bulge.fits xc=%i yc=%i theta=%.1f vec_output= \"profilebulge.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector resid.fits xc=%i yc=%i theta=%.1f vec_output= \"profileresid.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "logout"

	if model=='ser':
		print "rfits subcomps.fits 2 bulge.imh"
		print "wfits bulge.imh bulge.fits"
		print "rfits galfit_out.fits 2 model.imh"
		print "wfits model.imh model.fits"
		print "rfits galfit_out.fits 3 resid.imh"
		print "wfits resid.imh resid.fits"
		print "set imtype = \"fits\""
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"profilemodel.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector bulge.fits xc=%i yc=%i theta=%.1f vec_output= \"profilebulge.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector resid.fits xc=%i yc=%i theta=%.1f vec_output= \"profileresid.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "logout"

	if model=='exp':
		print "rfits subcomps.fits 2 disk.imh"
		print "wfits disk.imh disk.fits"
		print "rfits galfit_out.fits 2 model.imh"
		print "wfits model.imh model.fits"
		print "rfits galfit_out.fits 3 resid.imh"
		print "wfits resid.imh resid.fits"
		print "set imtype = \"fits\""
		print "pvector galaxy_clean.fits xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector model.fits xc=%i yc=%i theta=%.1f vec_output= \"profilemodel.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector disk.fits xc=%i yc=%i theta=%.1f vec_output= \"profiledisk.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "pvector resid.fits xc=%i yc=%i theta=%.1f vec_output= \"profileresid.fits\" out_type=\"image\"" % (xc_ima,yc_ima+z0/2.,theta)
		print "logout"

	sys.stdout = tmp_out
	f.close()
	subprocess.call("cl < profiles.cl -o", shell=True)


def appear(symb):
	size = 0
	line = 0
	if symb==0:
		symb = 'o'
		size = 3
		col = 'sienna'
		name = 'residual'

	if symb==1:
		symb = 'o'
		size = 7
		col = 'black'
		name = 'galaxy'

	if symb==2:
		symb = 'o'
		size = 3
		col = 'magenta'
		name = '2nd disk'

	if symb==3:
		symb = 'r-'
		line = 3
		col = 'red'
		name = 'model'

	if symb==4:
		symb = '-.'
		line = 4
		col = 'lime'
		name = 'disk'

	if symb==5:
		symb = 'r--'
		line = 4
		col = 'blue'
		name = 'bulge'

	if symb==6:
		symb = ':'
		line = 3
		col = 'turquoise'
		name = 'bar'

	if symb==7:
		symb = ':'
		line = 3
		col = 'lime'
		name = 'AGN'


	return symb,size,line,col,name


def appear1(symb):
	if symb==1:
		symb = 'o'
		size = 7
		col = 'black'
		name = 'galaxy'

	if symb==2:
		symb = 'o'
		size = 2
		col = 'magenta'
		name = '2nd disk'

	if symb==3:
		symb = 's'
		size = 5
		col = 'red'
		name = 'model'

	if symb==4:
		symb = '^'
		size = 3
		col = 'lime'
		name = 'disk'

	if symb==5:
		symb = 'D'
		size = 3
		col = 'blue'
		name = 'bulge'

	if symb==6:
		symb = 'o'
		size = 2
		col = 'turquoise'
		name = 'bar'

	if symb==7:
		symb = 'o'
		size = 2
		col = 'lime'
		name = 'AGN'


	return symb,size,col,name





def major_prof(m0,pix2sec,rmax,m_crit,model):
	m0 = m0 + 5.*log10(pix2sec)

	f = plt.figure(1,figsize=(fig_size[0],fig_size[1]))
	gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1],sharex=ax1)

	def read_1dprofile(ima_fits,m0,pix2sec,symb):
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		xc = list(scidata).index(max(scidata))

		Imajor = []  # the cut along the major axis
		r = []
		
		for k in range(0, nx):
			#if scidata[k]>0.:
			Imajor.append(scidata[k]) # In case if there is sky subtraction in these pixels
			r.append(pix2sec*(k - xc + 0.5))		# add 0.5pix shift (BUDDA)
		majorMag = m0 - 2.5*log10(Imajor)  # Poghson formula

		symb,size,line,col,name = appear(symb)
		if line==0:	ax1.plot(r,majorMag,symb,color=col,markersize=size, label = name)
		if size==0:	ax1.plot(r,majorMag,symb,color=col,lw=line, label = name)
		return r,majorMag
	r,majorMag1 = read_1dprofile("profileima.fits",m0,pix2sec,1)
	ax1.set_xlim(-rmax*pix2sec,rmax*pix2sec)
	r,majorMag2 = read_1dprofile("profilemodel.fits",m0,pix2sec,3)
	delta = majorMag1-majorMag2
	delta1 = []
	majorMag11 = []
	r1 = [] 
	for k in range(len(r)):
		if delta[k]>0.:
			majorMag11.append(majorMag1[k])
			delta1.append(delta[k])
			r1.append(r[k])
	ax1.set_ylim(m_crit,min(majorMag11)-1.0)

	ax2.plot(r1,delta1,'o',color='black',markersize=7)

	if model=='exp' or model=='edge' or model=='exp+ser' or model=='edge+ser':	r,majorMag = read_1dprofile("profiledisk.fits",m0,pix2sec,4)
	if model=='ser' or model=='exp+ser' or model=='edge+ser':	r,majorMag = read_1dprofile("profilebulge.fits",m0,pix2sec,5)
	ax1.legend(loc=1, numpoints=1)

	ax1.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax2.set_xlabel(r'r (arcsec)', fontsize=fsize)
	ax2.set_ylabel(r'$\Delta \mu$', fontsize=fsize)
	#ax2.set_ylim(-max(fabs(delta1)),max(fabs(delta1)))
	ax2.set_ylim(median(delta1)-std(delta1),median(delta1)+std(delta1))
	ax2.axhline(median(delta1),linestyle='--',color='black',lw=2)
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

	plt.savefig("major_prof.png", transparent = False)

def minor_prof(m0,pix2sec,rmax,m_crit,model):
	m0 = m0 + 5.*log10(pix2sec)

	f = plt.figure(2,figsize=(fig_size[0],fig_size[1]))
	gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1],sharex=ax1)

	def read_1dprofile(ima_fits,m0,pix2sec,symb):
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		xc = list(scidata).index(max(scidata))

		Imajor = []  # the cut along the major axis
		r = []
		
		for k in range(0, nx):
			#if scidata[k]>0.:
			Imajor.append(scidata[k]) # In case if there is sky subtraction in these pixels
			r.append(pix2sec*(k - xc + 0.5))		# add 0.5pix shift (BUDDA)
		majorMag = m0 - 2.5*log10(Imajor)  # Poghson formula

		symb,size,line,col,name = appear(symb)
		if line==0:	ax1.plot(r,majorMag,symb,color=col,markersize=size, label = name)
		if size==0:	ax1.plot(r,majorMag,symb,color=col,lw=line, label = name)
		return r,majorMag
	r,majorMag1 = read_1dprofile("profileima.fits",m0,pix2sec,1)
	ax1.set_xlim(-rmax*pix2sec,rmax*pix2sec)
	r,majorMag2 = read_1dprofile("profilemodel.fits",m0,pix2sec,3)
	delta = majorMag1-majorMag2
	delta1 = []
	majorMag11 = []
	r1 = [] 
	for k in range(len(r)):
		if delta[k]>0.:
			majorMag11.append(majorMag1[k])
			delta1.append(delta[k])
			r1.append(r[k])
	ax1.set_ylim(m_crit,min(majorMag11)-1.0)

	ax2.plot(r1,delta1,'o',color='black',markersize=7)

	if model=='exp' or model=='edge' or model=='exp+ser' or model=='edge+ser':	r,majorMag = read_1dprofile("profiledisk.fits",m0,pix2sec,4)
	if model=='ser' or model=='exp+ser' or model=='edge+ser':	r,majorMag = read_1dprofile("profilebulge.fits",m0,pix2sec,5)
	ax1.legend(loc=1, numpoints=1)

	ax1.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax2.set_xlabel(r'r (arcsec)', fontsize=fsize)
	ax2.set_ylabel(r'$\Delta \mu$', fontsize=fsize)
	ax2.set_ylim(median(delta1)-std(delta1),median(delta1)+std(delta1))
	ax2.axhline(median(delta1),linestyle='--',color='black',lw=2)
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

	plt.savefig("minor_prof.png", transparent = False)


def azim_prof(xc,yc,step,minsma,maxsma,m0,pix2sec):
	m0 = m0 + 5.*log10(pix2sec)
	sma,inten,inten_err = loadtxt('ellipse.txt', usecols=[1,2,3], unpack=True, skiprows = 6)
	r = sma*pix2sec
	mag1 = m0 - 2.5*log10(inten)
	mag_err = fabs((2.5/log(10.0)) * inten_err/inten)


	f = plt.figure(3,figsize=(fig_size[0],fig_size[1]))
	gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1],sharex=ax1)

	ax1.errorbar(r,mag1,yerr=mag_err,fmt='o',color='black', markersize=7, label='galaxy')

	#ax1.errorbar(r,mag1,yerr=mag_err,fmt='o',color='red', markersize=7, label='galaxy')
	ax1.set_ylim(max(mag1)+max(mag_err)+0.3,min(mag1)-max(mag_err)-1.)

	if os.path.exists('ellipse.txt')==True:
		os.remove(r"ellipse.txt")
	if os.path.exists('ell.cl')==True: 
		os.remove(r"ell.cl") 
	crea_ell("model.fits",xc,yc,step,minsma,maxsma)
	os.chmod(r"ell.cl",0777)
	subprocess.call("cl < ell.cl -o", shell=True)

	sma,inten,inten_err = loadtxt('ellipse.txt', usecols=[1,2,3], unpack=True, skiprows = 6)
	r = sma*pix2sec
	mag2 = m0 - 2.5*log10(inten)
	#mag_err = fabs((2.5/log(10.0)) * inten_err/inten)
	ax1.plot(r,mag2,'r-',color='red',lw=4, label = 'model')

	#print mag_err
	
	ax1.legend(loc=1, numpoints=1)
	ax2.set_xlabel(r'r (arcsec)', fontsize=fsize)
	ax1.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax2.set_ylabel(r'$\Delta \mu$', fontsize=fsize)
	
	try:
		delta = mag1-mag2

		r1 = []
		delta1 = []

		for k in range(len(delta)):
			if math.isnan(float(delta[k]))==False:
				r1.append(r[k])
				delta1.append(delta[k])


		ax2.plot(r1,delta1,'o',color='black',markersize=7)
		#ax2.set_ylim(median(delta1)-std(delta1),median(delta1)+std(delta1))
		ax2.set_ylim(-0.5,0.5)
		ax2.axhline(median(delta1),linestyle='--',color='black',lw=2)
	except:
		print 'There is a problem with azimuthal profile!'
	ax1.set_ylim(max(mag1),min(mag1)-0.5)
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

	plt.savefig("azim_prof.png", transparent = False)


def radius_prof(xc,yc,m0,pix2sec,rmax):
	def ima_plotter(ima_fits,xc,yc,m0,pix2sec,symb):
		f = plt.figure(4,figsize=(3,3))
		msize = 5.
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		ny,nx = scidata.shape
		I = []
		r = []
		for k in range (0,ny):
			for i in range (0,nx):
				I.append(scidata[k, i])	
				r.append( sqrt( (i-xc+0.5)**2 + (k-yc+0.5)**2 ) )
		r = np.array(r)*pix2sec
		mag = m0 - 2.5*log10(I)

		if symb==1:
			mag1 = []
			r1 = []
			for k in range(len(mag)):
				if mag[k]>0.:
					mag1.append(mag[k])
					r1.append(r[k])
			del mag
			del r
			mag = mag1
			r = r1
		if symb==1:		
			symb,size,col,name = appear1(symb)
			size = size/3.
			plt.plot(r,mag,'ro',markeredgecolor=col,markerfacecolor = 'white',markersize=0.4, label = name, alpha=1)
		else:
			symb,size,col,name = appear1(symb)
			plt.plot(r,mag,symb,markeredgecolor=col,markerfacecolor = col,markersize=0.4, label = name, alpha=1)
		return nx,ny,mag
	m0 = m0 + 5.*log10(pix2sec)
	#f = plt.figure(4,figsize=(3,3))

	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	minorLocator   = MultipleLocator(5)

	xsize = 9
	ysize = 9
	fsize = 9

	xticks_length = 3.
	xticks_width = 0.5

	#nx,ny,majorMag = ima_plotter("galaxy_clean.fits",xc,yc,m0,pix2sec,1)

	xlim(0,rmax*pix2sec)

	#ylim(np.median(majorMag),min(majorMag)-1.)

	nx,ny,majorMag = ima_plotter("galaxy_clean.fits",xc,yc,m0,pix2sec,1)
	plt.yticks(size=ysize)
	plt.xticks(size=xsize)
	plt.tick_params(length=xticks_length, width=xticks_width)
	plt.minorticks_on()
	plt.tick_params(length=xticks_length-5, width=xticks_width,which='minor',direction='out')
	plt.tight_layout(pad=1.0, h_pad=None, w_pad=None, rect=[0.03,0.03,1.01, 1.04])


	ylim(np.median(majorMag),min(majorMag)-1.)
	nx,ny,majorMag = ima_plotter("model.fits",xc,yc,m0,pix2sec,3)
	try:	
		nx,ny,majorMag = ima_plotter("disk.fits",xc,yc,m0,pix2sec,4)
	except IOError:
		print 'There is no disk!'
	try:
		nx,ny,majorMag = ima_plotter("bulge.fits",xc,yc,m0,pix2sec,5)
	except IOError:
		print 'There is no bulge!'


	legend(loc=1, numpoints=1,prop={'size':fsize})
	plt.xlabel(r'r (arcsec)', fontsize=fsize)
	plt.ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	plt.savefig("radius_prof.png",dpi=600, transparent = False)
	plt.clf()


def vert_prof(m0,pix2sec,rmax,m_crit):
	m0 = m0 + 5.*log10(pix2sec)

	f = plt.figure(5,figsize=(fig_size[0],fig_size[1]))
	gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1],sharex=ax1)

	def read_1dprofile(ima_fits,m0,pix2sec,symb):
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		yc = list(scidata).index(max(scidata))

		Imajor = []  # the cut along the major axis
		r = []
		
		for k in range(0, nx):
			#if scidata[k]>0.:
			Imajor.append(scidata[k]) # In case if there is sky subtraction in these pixels
			r.append(pix2sec*(k - yc + 0.5))		# add 0.5pix shift (BUDDA)
		majorMag = m0 - 2.5*log10(Imajor)  # Poghson formula

		symb,size,line,col,name = appear(symb)
		if line==0:	ax1.plot(r,majorMag,symb,color=col,markersize=size, label = name)
		if size==0:	ax1.plot(r,majorMag,symb,color=col,lw=line, label = name)
		return r,majorMag
	r,majorMag1 = read_1dprofile("zima.fits",m0,pix2sec,1)
	ax1.set_xlim(-rmax*pix2sec,rmax*pix2sec)
	r,majorMag2 = read_1dprofile("zmodel.fits",m0,pix2sec,3)
	delta = majorMag1-majorMag2
	delta1 = []
	majorMag11 = []
	r1 = [] 
	for k in range(len(r)):
		if delta[k]>0.:
			majorMag11.append(majorMag1[k])
			delta1.append(delta[k])
			r1.append(r[k])
	ax1.set_ylim(m_crit,min(majorMag11)-1.0)

	ax2.plot(r1,delta1,'o',color='black',markersize=7)
	ax1.legend(loc=1, numpoints=1)

	ax1.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax2.set_xlabel(r'r (arcsec)', fontsize=fsize)
	ax2.set_ylabel(r'$\Delta \mu$', fontsize=fsize)
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

	plt.savefig("vert_prof.png", transparent = False)	




#*******************************************************
#************ ALL PROFILES TOGETHER ********************
def edge_prof(file_input,xc,yc,m0,pix2sec,rmax,m_crit,zmax,z0,h):
	rmax = rmax*2.5
	zmax = zmax*2.5

	#m_crit = m_crit + 0.6
	m0 = m0 + 5.*log10(pix2sec)
	rmax = rmax*pix2sec
	zmax = zmax*pix2sec
	z0 = z0*pix2sec
	h = h*pix2sec
	xsize = 10
	ysize = int(xsize*z0/h)

	f = plt.figure(6,figsize=(fig_size[0],fig_size[1]))
	#gs = gridspec.GridSpec(3, 2,height_ratios=[1,3,1],width_ratios=[2,1])
	gs = gridspec.GridSpec(3, 2,height_ratios=[1,rmax/zmax,1],width_ratios=[2,1])


	ax1 = plt.subplot(gs[0,0])	# Isophote map
	ax2 = plt.subplot(gs[1,0],sharex=ax1)	# Major Axis Profile
	ax3 = plt.subplot(gs[1,1])	# Minor Axis Profile
	ax4 = plt.subplot(gs[2,0],sharex=ax1)	# Delta Major
	ax5 = plt.subplot(gs[2,1],sharex=ax3)	# Delta Minor


	# ISOPHOTE map
	hdulist = pyfits.open(file_input)
	scidata = hdulist[0].data
	ny,nx = np.shape(scidata)

	x = []
	y = []
	I = []
	I_crit = 10**(0.4*(m0-m_crit))

	for k in range(ny):
		for i in range(nx):
			if scidata[k,i]>I_crit:
				x.append((i + 0.5-xc)*pix2sec)
				y.append((k + 0.5-yc)*pix2sec)
				I.append(scidata[k,i])
	x = -np.array(x)
	y = np.array(y)
	I = np.array(I)

	m = -2.5*log10(I) + m0
	m = np.around(m,decimals=1)


	sample_pts = 500
	con_levels = np.arange(min(m),m_crit,0.5)
	#levels = np.arange(16,20,0.1)

	xi = np.linspace(-rmax,rmax,sample_pts)
	yi = np.linspace(-zmax,zmax,sample_pts)

	Ii = griddata(x,y,m,xi,yi)
	#print rmax

	#f = plt.figure()
	#ax1 = f.add_subplot(111)
	ax1.contour(xi,yi,Ii,con_levels,linewidths=1,colors='black')
	ax1.set_xlim(-rmax,rmax)
	ax1.set_ylim(-zmax,zmax)
	ax1.set_ylabel(r'z (arcsec)', fontsize=fsize)
	ax12 = ax1.twinx()
	ax12.set_ylabel(r'z (z$_0$)', fontsize=fsize)
	ax12.set_ylim(-zmax*pix2sec/z0,zmax*pix2sec/z0)
	ax12.set_yticks([0,2,4,6])
	ax12 = ax1.twiny()
	ax12.set_xlabel(r'r (h)', fontsize=fsize)
	ax12.set_xlim(-rmax*pix2sec/h,rmax*pix2sec/h)
	ax12.set_xticks([-6,-4,-2,0,2,4,6])

	# MINOR AXIS PROFILE
	def read_1dprofile(ima_fits,m0,pix2sec,symb):
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		xc = list(scidata).index(max(scidata))

		Imajor = []  # the cut along the major axis
		r = []
		
		for k in range(0, nx):
			#if scidata[k]>0.:
			Imajor.append(scidata[k]) # In case if there is sky subtraction in these pixels
			r.append(pix2sec*(k - xc + 0.5))		# add 0.5pix shift (BUDDA)
		majorMag = m0 - 2.5*log10(Imajor)  # Poghson formula

		symb,size,line,col,name = appear(symb)
		if line==0:	ax3.plot(r,majorMag,symb,color=col,markersize=size, label = name)
		if size==0:	ax3.plot(r,majorMag,symb,color=col,lw=line, label = name)
		return r,majorMag
	r,majorMag1 = read_1dprofile("profileima.fits",m0,pix2sec,1)
	ax2.set_xlim(-rmax,rmax)
	ax5.set_xlim(-rmax,rmax)
	r,majorMag2 = read_1dprofile("profilemodel.fits",m0,pix2sec,3)
	delta = majorMag1-majorMag2
	delta1 = []
	majorMag11 = []
	r1 = [] 
	for k in range(len(r)):
		if delta[k]>0.:
			majorMag11.append(majorMag1[k])
			delta1.append(delta[k])
			r1.append(r[k])
	#ax2.set_ylim(m_crit+1.,min(majorMag11)-1.0)

	ax5.plot(r1,delta1,'o',color='black',markersize=4)
	#try:	
	#	r,majorMag = read_1dprofile("profiledisk.fits",m0,pix2sec,4)
	#except:
	#	print 'There is no disk!'
	#try:
	#	r,majorMag = read_1dprofile("profilebulge.fits",m0,pix2sec,5)
	#except:
	#	print 'There is no bulge!'
	ax2.legend(loc=1, numpoints=1)

	ax2.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax5.set_xlabel(r'z (arcsec)', fontsize=fsize)
	ax4.set_ylabel(r'$\Delta \mu$', fontsize=fsize)



	#MAJOR AXIS PROFILE
	def read_1dprofile(ima_fits,m0,pix2sec,symb):
		hdulist = pyfits.open(ima_fits)
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		yc = list(scidata).index(max(scidata))

		Imajor = []  # the cut along the major axis
		r = []
		
		for k in range(0, nx):
			#if scidata[k]>0.:
			Imajor.append(scidata[k]) # In case if there is sky subtraction in these pixels
			r.append(pix2sec*(k - yc + 0.5))		# add 0.5pix shift (BUDDA)
		majorMag = m0 - 2.5*log10(Imajor)  # Poghson formula

		symb,size,line,col,name = appear(symb)
		if line==0:	ax2.plot(r,majorMag,symb,color=col,markersize=size, label = name)
		if size==0:	ax2.plot(r,majorMag,symb,color=col,lw=line, label = name)
		return r,majorMag
	r,majorMag1 = read_1dprofile("zima.fits",m0,pix2sec,1)
	ax3.set_xlim(-zmax,zmax)
	r,majorMag2 = read_1dprofile("zmodel.fits",m0,pix2sec,3)

	delta = majorMag1-majorMag2
	delta2 = []
	majorMag12 = []
	r1 = [] 
	for k in range(len(r)):
		if delta[k]>0.:
			majorMag12.append(majorMag1[k])
			delta2.append(delta[k])
			r1.append(r[k])
	try:
		ax3.set_ylim(m_crit+1.,min([min(majorMag11),min(majorMag12)])-1.0)
		ax2.set_ylim(m_crit+1.,min([min(majorMag11),min(majorMag12)])-1.0)
		#ax2.set_ylim(m_crit+1.,min(majorMag11)-1.0)
	except:
		ax3.set_ylim(-0.5,0.5)
		ax2.set_ylim(-0.5,0.5)		

	ax4.plot(r1,delta2,'o',color='black',markersize=4)

	try:	
		r,majorMag = read_1dprofile("profiledisk.fits",m0,pix2sec,4)
	except:
		print 'There is no disk!'
	try:
		r,majorMag = read_1dprofile("profilebulge.fits",m0,pix2sec,5)
	except:
		print 'There is no bulge!'
	#ax4.set_ylim(min([min(delta1),min(delta2)]),max([max(delta1),max(delta2)]))
	#ax5.set_ylim(min([min(delta1),min(delta2)]),max([max(delta1),max(delta2)]))

	#ax4.set_ylim(-max([max(fabs(delta1)),max(fabs(delta2))]),max([max(fabs(delta1)),max(fabs(delta2))]))
	#ax5.set_ylim(-max([max(fabs(delta1)),max(fabs(delta2))]),max([max(fabs(delta1)),max(fabs(delta2))]))

	ax4.set_ylim(-0.5,0.5)
	ax5.set_ylim(-0.5,0.5)

	#ax3.legend(loc=1, numpoints=1)
	ax32 = ax3.twiny()
	ax32.set_xlabel(r'z (z$_0$)', fontsize=fsize)
	ax32.set_xlim(-zmax*pix2sec/z0,zmax*pix2sec/z0)
	ax32.set_xticks([0,2,4,6])

	#ax3.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	ax4.set_xlabel(r'r (arcsec)', fontsize=fsize)
	#ax4.set_ylabel(r'$\Delta \mu$', fontsize=fsize)
	plt.setp(ax5.get_yticklabels(), visible=False)
	plt.setp(ax3.get_yticklabels(), visible=False)
	ax5.axhline(0.,linestyle='--',color='black',lw=2)
	ax4.axhline(0.,linestyle='--',color='black',lw=2)

	plt.setp([a.get_xticklabels() for a in f.axes[:0]], visible=False)
	plt.setp([a.get_yticklabels() for a in f.axes[:0]], visible=False)
	f.subplots_adjust(hspace=0)
	f.subplots_adjust(wspace=0)
	plt.savefig("edge_prof.png", transparent = False)	


def pictures(model):
	try:
		import aplpy
	except:
		print "AplPy module is not found!"
	def read_fits(gal): 
		# The function to read the fits file.
		hdulist = pyfits.open(gal)
		prihdr = hdulist[2].header
		scidata = hdulist[2].data
		return scidata
	def image(gal,name,number,ratio):
		dx = 5.
		dy = 5./ratio
		'''
		if dy>5.:
			dy = 5.
			print 'I am here!!!!!!^&^&E'
		'''
		try:
			fig = aplpy.FITSFigure(gal, hdu=number, figsize=(dx, dy))
			fig.show_grayscale(stretch='log',vmid=-10,vmin=-5, invert=True)
			fig.show_contour(gal, colors='red')
			#fig.tick_labels.set_font(size='small')
			fig.axis_labels.hide() # to hide axis labels
			fig.tick_labels.hide() # to hide tick labels
			#fig.add_label(0.2, 0.9, name, relative=True,size=28)
			fig.add_scalebar(0.002777)
			fig.scalebar.set_label(' 10" ')
			pic = name + '.png'
			savefig(pic,bbox_inches=0)
		except:
			print "AplPy module is not found!"

	galfit_out = 'galfit_out.fits'
	scidata=read_fits(galfit_out)
	[ny,nx] = scidata.shape

	l=200.
	ratio=float(nx)/float(ny)
	print ratio

	#ratio = 3
	dimy = l/(3.*ratio)
	dimx=dimy*ratio


	try:
		import ds9
		ds9.ds9_image(galfit_out)
	except:
		image(galfit_out,'galaxy',1,ratio)
		image(galfit_out,'model',2,ratio)
		image(galfit_out,'residual',3,ratio)

	galfit_out = 'subcomps.fits'
	scidata=read_fits(galfit_out)
	[ny,nx] = scidata.shape
	l=200.
	ratio=float(nx)/float(ny)
	#ratio = 3
	dimy = l/(3.*ratio)
	dimx=dimy*ratio

	if model=='exp' or model=='edge':
		image(galfit_out,'disk',2,ratio)

	if model=='ser':
		image(galfit_out,'bulge',2,ratio)

	if model=='exp+ser' or model=='edge+ser':
		image(galfit_out,'disk',2,ratio)
		image(galfit_out,'bulge',3,ratio)




#radius_prof(695,718,27.92,0.396,265.)
#azim_prof(695,718,0.2,1.,265.,27.92,0.396)
##pictures('exp+ser')
#isophotes('galaxy_clean.fits',695,718,27.92,0.396,21.89,265.,265.*1387/1432,h=77.)

'''
profiles_galfit_extract(695,718,0.,'exp+ser',z0=0.,h=0.)
major_prof(27.92,0.396,265.,21.89,'exp+ser')

profiles_galfit_extract(695,718,90.,'exp+ser',z0=0.,h=0.)
minor_prof(27.92,0.396,265.,21.89,'exp+ser')
'''

#gadotti_style_profile(498,350,26.0,0.396,150)
##profiles_galfit_extract(100,50,0,4.17,18.)
#profiles_plot(26.,0.396,70,23.5)
#z_profile(26.,0.396,70,23.5)

#prof_tog('galaxy_clean.fits',100,50,26.,0.396,80,22.0,30,4.17,18.)

#step = 0.1
#minsma = 1.0
#azim_aver_plot(100,50,step,minsma,26.,0.396)

#gadotti_style_profile(100,50,26.,0.396,70.)








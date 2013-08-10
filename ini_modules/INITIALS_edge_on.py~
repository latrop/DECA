#!/usr/bin/python
# Module to make vertical cuts for rotated edge-on galaxy
# You should first find rmin and rmax to do this
# -*- coding:  cp1251 -*-
# Program to rotate and extract fits-files


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

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')


import setup

#import superellipse

rmin_bulge_coeff = setup.rmin_bulge
Ccr = setup.Ccr
El_cr = setup.El_cr
#*************** CUTS ALONG MAJOR AND MINOR AXES ********************
def z_cut(scidata,yc,x,m0,zmin1,zmin2,zmax,backgr,plot=0):
	#*** The cut along the minor axis (y-axis) ***
	# INPUT:
	# 	yc - pixel coordinate of the center [pix]
	# 	x - the cut at this point (x-axis) [pix]
	# 	zmin1>=<0.,zmin2>=<0., zmax>0. - borders where to investigate [pix], zmin1<=zmin2<zmax
	# 	plot = 1 if you want to plot this cut
	# OUTPUT:
	#	z - shifted scaled coordinates to the center [pix]
	# 	z1 - unshifted scaled coordinates [pix]
	#	Iminor - Surface brightness [DN]
 
	ny,nx = scidata.shape
	Iminor = []
	z = []
	z1 = []
	
	for k in range(0, ny):		
		if x<nx:
			if fabs(k- yc+0.5)<zmax and (k- yc+0.5<zmin1 or k- yc+0.5>zmin2) and scidata[k, x]>1*backgr:
				Iminor.append(scidata[k, x])
				z.append(k- yc+1.)
				z1.append(k+1.)
	# Earlier it was:
	'''
	for k in range(0, ny):		
		if fabs(k- yc+0.5)<zmax and (k- yc+0.5<zmin1 or k- yc+0.5>zmin2) and scidata[k, x]>1*backgr:
			Iminor.append(scidata[k, x])
			z.append(k- yc+1.)
			z1.append(k+1.)
	'''


	if plot==1:
		# Plotting
		plt.plot(z, Iminor,'r-',color='black', lw=3)
		xlim(-ny/2,ny/2)
		ylim(floor(min(Iminor)),ceil(max(Iminor)+np.mean(Iminor)/4.))
		plt.xlabel(r'r (pix)', fontsize=30)
		plt.ylabel(r'Intensity (DN)', fontsize=30)
	return np.array(z),np.array(z1),np.array(Iminor)

def r_cut(scidata,xc,y,pix2sec,m0,rmin,rmax,backgr,plot=0):
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
		if fabs(k-xc+0.5)<rmax and fabs(k-xc+0.5)>rmin and scidata[y, k]>1.*backgr:
			Imajor.append(scidata[y, k])
			x.append(k- xc+1.)
			x1.append(k+1.)

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



#*************** FUNCTIONS TO FIT ********************
def line(B, x):
	#*** For non edge-on disk SB in mag/arcsec^2 (radial) ***
	# B[0]=h
	# B[1]=m0d
	x = np.array(x)
	return (1.0857/B[0])*fabs(x) + B[1]

def disk_edge_soph(B, z):
	#*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
	# B[0] = I0d
	# B[1] = n
	# B[2] = z0 
	# B[3] = z_c
	z = np.array(z)	
	return B[0] * (1.0 / np.cosh(B[1]*fabs(z-B[3])/(2.*B[2])) ) **(2./B[1])

def disk_edge(B,z):
	#*** For edge-on disk SB in mag/arcsec^2 (along z-axis).  Sech^2-law. ***
	# B[0] = I0d
	# B[1] = z0 
	# B[2] = z_c
	z = np.array(z)	
	return B[0] * ( 1. / (np.cosh(fabs(z-B[2])/B[1]))**2   )

def disk_edge_r(B,x):
	#*** For edge-on disk SB in mag/arcsec^2 (along r-axis). ***
	# B[0] = I0d
	# B[1] = h
	# B[2] = x_c
	x = np.array(x)
	return B[0] * fabs(x-B[2])/B[1] * sp.special.kn(1,fabs(x-B[2])/B[1])

def bulge_func(B, x):
	return B[0]*np.exp( -(1.9992*B[2] - 0.3271)*((fabs(x-B[3])/B[1])**(1./B[2]) - 1.) )

def bulge_f(B, x):
	return B[0]*np.exp( -(1.9992*B[2] - 0.3271)*((fabs(x)/B[1])**(1./B[2]) - 1.) )

def bul_disk_f(B, z):
	z = np.array(z)	
	return B[0] * ( 1. / (np.cosh(fabs(z)/B[1]))**2   ) + B[2]*np.exp( -(1.9992*B[4] - 0.3271)*((fabs(z)/B[3])**(1./B[4]) - 1.) )

def bul_disk_e(B, r):
	r = np.array(r)	
	return B[0] * fabs(r)/B[1] * sp.special.kn(1,fabs(r)/B[1]) + B[2]*np.exp( -(1.9992*B[4] - 0.3271)*((fabs(r)/B[3])**(1./B[4]) - 1.) )

#*************** FITTING FUNCTIONS  ********************
'''
def fitting_extinction(z,I,eta0,k0,zd,zs,h,r):
	def extinction(B,z):
		# B[0] = eta0
		# B[1] = k0
		# B[2] = zd 
		# B[3] = zs
		return ( (B[0]/B[1]) * ( np.cosh(fabs(z)/B[2]) / np.cosh(fabs(z)/B[3])) ** 2) * (1. - np.exp(-2.*B[1]*r*sp.special.kn(1,r/h)/((np.cosh(fabs(z)/B[2])) ** 2) ) )
	dataToFit = RealData(z, I)
	mo = Model(extinction)
	fitting = ODR(dataToFit, mo, [eta0,k0,zd,zs])
	fitting.set_job()
	fit_result = fitting.run()
	eta0 = fit_result.beta[0]
	eta0_e = fit_result.sd_beta[0]
	k0 = fit_result.beta[1]
	k0_e = fit_result.sd_beta[1]
	zd = fit_result.beta[2]
	zd_e = fit_result.sd_beta[2]
	zs = fit_result.beta[3]
	zs_e = fit_result.sd_beta[3]
	return eta0,eta0_e,k0,k0_e,zd,zd_e,zs,zs_e
def edge_ext(z,r,h,eta0,k0,zd,zs):
	# B[0] = eta0
	# B[1] = k0
	# B[2] = zd 
	# B[3] = zs
	return ( (eta0/k0) * (np.cosh(fabs(z)/zd ) / np.cosh(fabs(z)/zs)) ** 2) * (1. - np.exp(-2.*k0*r*sp.special.kn(1,r/h)/((np.cosh(fabs(z)/zd )) ** 2) ))

def eta(ro,z,eta0,hs,zs):
	return eta0*np.exp(-ro/hs)/((np.cosh(fabs(z)/zs)) ** 2)

def k(ro,z,k0,hd,zd):
	return k0*np.exp(-ro/hd)/((np.cosh(fabs(z)/zd)) ** 2)

def full_fitting(galaxy_image,xc,yc,rtr,rmin,rmax,zmin,zmax,pix2sec,m0):
	# scales are in pixels!
	xc = np.ceil(xc-0.5)
	yc = np.ceil(yc-0.5)
	xmin = xc + int(rmin)
	xmax = xc + int(rmax)

	Z0 = []
	xx = []
	I0 = []
	# STEP: in the region outside the bulge and dust
	for x in arange(xmax-2*rmax,xmax,5):
		if (x<xmax and x > xmin):
		#if (x<xmax and x > xmin) or (x<xmin-2*rmin and x>xmax-2*rmax):
		#if (x<xmin-2*rmin and x>xmax-2*rmax):
			z,I,ny = z_cut(galaxy_image,yc,x,pix2sec,m0,zmin,zmax)
			I0d,I0d_e,z0,z0_e,z_c = fitting(z,I,max(I),(rtr/20.)*pix2sec,yc*pix2sec)
			I0.append(I0d)
			Z0.append(z0)
			xx.append(x*pix2sec)
			print x
	#plt.show()
	plt.plot(z-z_c, Iz(z-z_c,z0,I0d),'r-',color='blue', lw=3)
	hdulist = pyfits.open(galaxy_image)
	scidata = hdulist[0].data
	I0d = scidata[yc,xc]
	z0 = np.median(Z0)
	z0_err = np.std(Z0)
	print 'z0 = %.2f +/- %.2f [arcsec]' % (z0,z0_err) 

	I0new = []
	xxnew = []
	Ilim = np.median(I0) + 1.*np.std(I0)
	for k in range(len(I0)):
		if I0[k]<Ilim:
			I0new.append(I0[k])
			xxnew.append(xx[k])

	I0d,I0d_e,h,h_err,x_c = fitting_radius(xxnew,I0new,2500.,108.,xc*pix2sec)
	#plt.figure(1)
	#plt.plot(xxnew,I0new,'o',color='blue',markersize=5)
	#plt.show()
	m0d = -2.5*log10(I0d) + m0 
	m0d_err = 2.5*I0d_e/(log(10.)*I0d)
	print 'm0d = %.2f +/- %.2f [mag/arcsec^-2],\t h = %.2f +/- %.2f [arcsec]' % (m0d,m0d_err,h,h_err)

	#plt.figure(1)
	#plt.plot(xxnew-x_c,-2.5*log10(I0new) + m0 + 5*log10(pix2sec),'o',color='blue',markersize=5)
	#plt.plot(xxnew-x_c,-2.5*log10(Iedge(fabs(xxnew-x_c),h,I0d)) + m0 + 5*log10(pix2sec),'r-',color='red',lw=4)
	#plt.show()
	#print I0d
	# STEP: Extinction inclusion
	tau = []
	zdd = []
	zss = []
	for k in range(len(xx)):
		#k0 = 0.16
		k0 = 0.01
		r = fabs(xx[k]-x_c)
		eta0,eta0_e,k0,k0_e,zd,zd_e,zs,zs_e = fitting_extinction(z,I,82600.*pix2sec**2,k0,8.,11.,150,r)
		#tau.append(2.*k0*r*sp.special.kn(1,r/h)/((np.cosh(fabs(z)/zd)) ** 2))
		tau.append(2.*k0*z)
		zdd.append(zd)
		zss.append(zs)
		#print xx[k]
		plt.plot(z-z_c, edge_ext(fabs(z-z_c),r,150,eta0/(2*h),k0,zd,zs),'r-',color='red', lw=4)
	#print zdd,zss
	#print eta0
	#print median(zdd),median(zss),median(tau)
	#plt.show()
'''
  

def fitting_soph(z,I,I0d,n,z0,z_c,plot=0):
	# *** To fit z-cuts with sophisticated law ***
	# All is given in arcsec
	z_c1 = z_c
	dataToFit = RealData(z, I)
	mo = Model(disk_edge_soph)
	fitting = ODR(dataToFit, mo, [I0d,n,z0,z_c])
	fitting.set_job()
	fit_result = fitting.run()
	I0d = fit_result.beta[0]
	I0d_e = fit_result.sd_beta[0]
	n = fit_result.beta[1]
	n_e = fit_result.sd_beta[1]
	z0 = fit_result.beta[2]
	z0_e = fit_result.sd_beta[2]
	z_c = fit_result.beta[3]

	if plot==1:
		plt.plot(z-z_c, I,'r-',color='black', lw=3)
		plt.plot(z-z_c, disk_edge_soph([I0d,n,z0,z_c], z),'r--',color='blue', lw=1)
		plt.plot(z-z_c, fabs(I - disk_edge_soph([I0d,n,z0,z_c], z)),'r--',color='red', lw=1)
	return I0d,I0d_e,n,n_e,z0,z0_e,z_c1

def fitting(z,I,I0d,z0,z_c,plot=0):
	# *** To fit z-cuts with sech^2 law ***
	# All is given in arcsec
	z_c1 = z_c
	z = np.array(z)
	dataToFit = RealData(z, I)
	mo = Model(disk_edge)
	fitting = ODR(dataToFit, mo, [I0d,z0,z_c])
	fitting.set_job()
	fit_result = fitting.run()
	I0d = fit_result.beta[0]
	I0d_e = fit_result.sd_beta[0]
	z0 = fit_result.beta[1]
	z0_e = fit_result.sd_beta[1]
	z_c = fit_result.beta[2]

	if plot==1:
		plt.plot(z-z_c1, I,'r-',color='black', lw=3)
		plt.plot(z-z_c1, disk_edge([I0d,z0,z_c1],z),'r--',color='blue', lw=1)
		plt.plot(z-z_c1, fabs(I - disk_edge([I0d,z0,z_c1],z)),'r--',color='red', lw=1)
	return I0d,I0d_e,z0,z0_e,z_c

def fitting_radius(x,I,I0d,h,x_c):
	# *** To fit r-cuts of edge-on galaxies ***
	# All is given in arcsec
	dataToFit = RealData(x, I)
	mo = Model(disk_edge_r)
	fitting = ODR(dataToFit, mo, [I0d,h,x_c])
	fitting.set_job()
	fit_result = fitting.run()
	I0d = fit_result.beta[0]
	I0d_e = fit_result.sd_beta[0]
	h = fit_result.beta[1]
	h_e = fit_result.sd_beta[1]
	x_c = fit_result.beta[2]
	return I0d,I0d_e,h,h_e,x_c


#************* OTHER USEFUL FUNCTIONS **************
def z0_est(z,I):
	# Function to estimate z0 for sech^2 law
	I_z0 = max(I)/2.3811 
	tck = interpolate.splrep(z,I,s=0)
	z_interp = np.arange(min(z),max(z),fabs(max(z)-min(z))/1000.)
	I_interp = interpolate.splev(z_interp,tck,der=0)

	'''		
	plt.plot(z_interp,I_interp,'*')
	plt.plot(z,I,'o')
	plt.show()
	exit()
	'''
	for k in range(len(z_interp)):
		if I_interp[k]>I_z0:
			z0 = fabs(z_interp[k])
	return z0



def bulge_profile_minor(scidata,xc,z_c,zmin1,zmin2,zmax,m0d,z0,m0,pix2sec,backgr_dev,rmin_bulge,dust):
	ny,nx = scidata.shape
	dust = -1


	z,z1,I = z_cut(scidata,z_c,xc,m0,ny/2.,-ny/2.,zmax,backgr_dev,plot=0)


	z_new = []
	I_new = []
	for k in range(len(z)):
		if dust==-1 and z[k]>=0.:
			z_new.append(z[k])
			I_new.append(I[k])
		if dust==0 and z[k]>=0.:
			z_new.append(z[k])
			I_new.append(I[k])
		if dust==1 and z[k]<=0.:
			z_new.append(-z[k])
			I_new.append(I[k])
		
	tck = interpolate.splrep(z_new,I_new,s=0)
	z_interp = np.arange(min(z_new),max(z_new),fabs(max(z_new)-min(z_new))/1000.)
	I_interp = interpolate.splev(z_interp,tck,der=0)
	z_new = z_interp
	I_new = I_interp


	I_res = I_new	# - disk_edge([10**(0.4*(m0-m0d)),z0,z_c],np.array(z_new)+z_c)

	z_new1 = []
	I_new1 = []



	for k in range(len(I_res)):
		if I_res[k]>0. and z_new[k]>rmin_bulge and I_res[k]>1.*backgr_dev: # and  z_new[k]<zmax:
			z_new1.append(z_new[k])
			I_new1.append(I_res[k])

	z_new2 = []
	I_new2 = []

	I_res = I_new - disk_edge([10**(0.4*(m0-m0d+5.*log10(pix2sec))),z0,z_c],np.array(z_new)+z_c)
	for k in range(len(I_res)):
		if I_res[k]>0. and I_res[k]>1.*backgr_dev: # and  z_new[k]<zmax:
			z_new2.append(z_new[k])
			I_new2.append(I_res[k])


	'''
	plt.plot(z_new,-2.5*log10(disk_edge([10**(0.4*(m0-m0d)),z0,z_c],np.array(z_new)+z_c))+m0+5.*log10(pix2sec))
	plt.plot(z_new1,-2.5*log10(I_new1)+5.*log10(pix2sec)+m0,'*',markersize=15)
	plt.plot(z_new,-2.5*log10(I_new)+m0+5.*log10(pix2sec),'s')
	plt.plot(z,-2.5*log10(I)+m0+5.*log10(pix2sec),'o')

	plt.show()
	exit()
	'''

	
	return z_new1,I_new1,z_new2,I_new2



def bulge_profile_major(scidata,xc,z_c,zmin1,zmin2,zmax,m0d,h,m0,pix2sec,backgr_dev,rmin_bulge,dust):
	ny,nx = scidata.shape
	dust = -1

	x,x1,I = r_cut(scidata,xc-1.,z_c,pix2sec,m0,0.,nx/2.,backgr_dev,plot=0)



	# Averaging:
	rr = []
	II = []
	for k in range(len(x)):
		if x[k]<0.:
			rr.append(int(floor(x[k])))
			II.append(I[k])
		if x[k]>0.:
			rr.append(int(ceil(x[k])))
			II.append(I[k])
	rrr = []
	III = []
	k_left = rr.index(-1)
	k_right = rr.index(1)
	for xx in arange(k_right,len(rr),1):		
		for yy in arange(0,k_left,1):
			if fabs(rr[xx])==fabs(rr[yy]):
				rrr.append(rr[xx])
				III.append((II[xx]+II[yy])/2.)

	tck = interpolate.splrep(rrr,III,s=0)
	r_interp = np.arange(min(rrr),max(rrr),fabs(max(rrr)-min(rrr))/1000.)
	I_interp = interpolate.splev(r_interp,tck,der=0)
	rrr = r_interp
	III = I_interp


	I_res = III - disk_edge_r([10**(0.4*(m0-m0d+5.*log10(pix2sec))),h,0.],rrr)

	r_new2=[]
	I_new2=[]
	for k in range(len(I_res)):
		if I_res[k]>0. and I_res[k]>1*backgr_dev and rrr[k]>rmin_bulge and rrr[k]<5.*h:
			r_new2.append(rrr[k])
			I_new2.append(I_res[k])

	I_res = III	# - disk_edge_r([10**(0.4*(m0-m0d)),h,0.],rrr)

	r_new1=[]
	I_new1=[]
	for k in range(len(III)):
		if III[k]>0. and rrr[k]>rmin_bulge and III[k]>1*backgr_dev and rrr[k]<5.*h:
			r_new1.append(rrr[k])
			I_new1.append(III[k])

	'''
	plt.plot(x,-2.5*log10(I)+m0+5.*log10(pix2sec),'o')
	plt.plot(rrr,-2.5*log10(III)+m0+5.*log10(pix2sec),'o')
	plt.plot(rrr,-2.5*log10(I_res)+m0+5.*log10(pix2sec),'*')
	plt.plot(r_new1,-2.5*log10(I_new1)+m0+5.*log10(pix2sec),'s',markersize=10)
	plt.show()
	'''
	return r_new1,I_new1,r_new2,I_new2

def bulge_fitting1(r,I,r1,I1,m0,pix2sec,ellip,C31,m0d,z0,h):
	L_bulge = trapz(I1, r1)
	S = 0.

	for k in range(len(r1)):
		if S<=L_bulge/2.:
			S = (r1[k+1]-r1[k])*(I1[k]+I1[k+1])/2.+S
			reb = r1[k]
			Ieb = I1[k]


	meb = -2.5*log10(Ieb)+m0+5.*log10(pix2sec)# + 0.5	# !!! ADD 1 mag !!!
	if ellip>El_cr and C31<Ccr:
		n = 2.0
	else:
		n = 3.4

	meb = meb + 1.5
	reb = reb*2.

	'''
	dataToFit = RealData(r, I)
	mo = Model(bulge_f)
	fitting = ODR(dataToFit, mo, [10**(0.4*(m0-meb+5.*log10(pix2sec))),reb,n])
	fitting.set_job()
	fit_result = fitting.run()
	a = fit_result.beta[0]
	b = fit_result.beta[1]
	c = fit_result.beta[2]
	eps = fabs(median(fit_result.eps))
	#print -2.5*log10(a)+m0+5.*log10(pix2sec),b,c

	print '\n'
	print 'm0d = %.2f [mag arcsec^-2]' % (m0d) 
	print 'z0 = %.2f [arcsec]' % (z0*0.396) 
	print 'meb = %.2f [mag arcsec^-2]' % (-2.5*log10(a)+m0+5.*log10(pix2sec)) 
	print 'reb = %.2f [arcsec]' % (b*pix2sec) 
	print 'n = %.2f ' % (c) 
	return m0d,z0,-2.5*log10(a)+m0+5.*log10(pix2sec),b,c
	'''
	I0d = 10**(0.4*(m0-m0d+5.*log10(pix2sec)))

	dataToFit = RealData(r, I)
	mo = Model(bul_disk_f)
	fitting = ODR(dataToFit, mo, [I0d,z0,10**(0.4*(m0-meb+5.*log10(pix2sec))),reb,n])
	fitting.set_job()
	fit_result = fitting.run()
	a = fit_result.beta[0]
	b = fit_result.beta[1]
	c = fit_result.beta[2]
	d = fit_result.beta[3]
	e = fit_result.beta[4]

	if (e<=0.5 or e>=5.) and ellip>El_cr:
		e = 2.0
	elif (e<=0.5 or e>=5.) and ellip<El_cr:
		e = 3.4

	print '\n'
	print -2.5*log10(a)+m0+5.*log10(pix2sec)-2.5*log10(b/h)
	print 'm0d = %.2f [mag arcsec^-2]' % (-2.5*log10(a)+m0+5.*log10(pix2sec)-2.5*log10(b/h)) 
	print 'z0 = %.2f [arcsec]' % (b*0.396) 
	print 'meb = %.2f [mag arcsec^-2]' % (-2.5*log10(c)+m0+5.*log10(pix2sec)) 
	print 'reb = %.2f [arcsec]' % (d*pix2sec) 
	print 'n = %.2f' % (e) 
	

	#print -2.5*log10(a)+m0+5.*log10(pix2sec),b,-2.5*log10(c)+m0+5.*log10(pix2sec),d,e
	return m0d,z0,-2.5*log10(c)+m0+5.*log10(pix2sec),d,e




def bulge_fitting2(r,I,r1,I1,m0,pix2sec,ellip,C31,m0d,h,z0):
	'''
	plt.figure(99)
	plt.plot(np.array(r1)*pix2sec,I1)
	plt.show()
	'''
	L_bulge = trapz(I1, r1)
	S = 0.

	for k in range(len(r1)):
		if S<=L_bulge/2.:
			S = (r1[k+1]-r1[k])*(I1[k]+I1[k+1])/2.+S
			reb = r1[k]
			Ieb = I1[k]


	meb = -2.5*log10(Ieb)+m0+5.*log10(pix2sec) # + 1.	# !!! ADD 1 mag !!!
	if ellip>El_cr and C31<Ccr:
		n = 2.0
	else:
		n = 3.4

	#print m0d,h*pix2sec,meb,reb*pix2sec,'\n'

	meb = meb + 1.
	reb = reb*2.
	dataToFit = RealData(r, I)
	mo = Model(bul_disk_e)
	fitting = ODR(dataToFit, mo, [10**(0.4*(m0-m0d+5.*log10(pix2sec))),h,10**(0.4*(m0-meb+5.*log10(pix2sec))),reb,n])
	fitting.set_job()
	fit_result = fitting.run()
	a = fit_result.beta[0]
	b = fit_result.beta[1]
	c = fit_result.beta[2]
	d = fit_result.beta[3]
	e = fit_result.beta[4]
	'''

	dataToFit = RealData(r1, I1)
	mo = Model(bulge_f)
	fitting = ODR(dataToFit, mo, [10**(0.4*(m0-meb+5.*log10(pix2sec))),reb,n])
	fitting.set_job()
	fit_result = fitting.run()

	c = fit_result.beta[0]
	d = fit_result.beta[1]
	e = fit_result.beta[2]


	print '\n'
	print 'm0d = %.2f [mag arcsec^-2]' % (m0d) 
	print 'h = %.2f [arcsec]' % (h*pix2sec) 
	print 'meb = %.2f [mag arcsec^-2]' % (-2.5*log10(c)+m0+5.*log10(pix2sec)) 
	print 'reb = %.2f [arcsec]' % (d*pix2sec) 
	print 'n = %.2f' % (e) 
	'''

	if (e<=0.5 or e>=5.) and ellip>El_cr:
		e = 2.0
	elif (e<=0.5 or e>=5.) and ellip<El_cr:
		e = 3.4


	print '\n'
	print 'm0d = %.2f [mag arcsec^-2]' % (-2.5*log10(a)+m0+5.*log10(pix2sec)-2.5*log10(z0/b)) 
	print 'h = %.2f [arcsec]' % (b*pix2sec) 
	print 'meb = %.2f [mag arcsec^-2]' % (-2.5*log10(c)+m0+5.*log10(pix2sec)) 
	print 'reb = %.2f [arcsec]' % (d*pix2sec) 
	print 'n = %.2f' % (e) 
	
	#return m0d,h,-2.5*log10(c)+m0+5.*log10(pix2sec),d,e
	return -2.5*log10(a)+m0+5.*log10(pix2sec),b,-2.5*log10(c)+m0+5.*log10(pix2sec),d,e

def bulge(scidata,xc,z_c,zmin1,zmin2,zmax,m0d,z0,m0,pix2sec,backgr_dev):
	z,z1,I = z_cut(scidata,z_c,xc,m0,zmin1,zmin2,zmax,backgr_dev,plot=0)
	I_res = I - disk_edge([10**(0.4*(m0-m0d)),z0,z_c],z1)


	I_bulge1 = []
	I_bulge2 = []

	#plt.plot(z,-2.5*log10(I)+m0+5.*log10(pix2sec),'o')
	#plt.plot(z,-2.5*log10(disk_edge([10**(0.4*(m0-m0d)),z0,z_c],z1))+m0+5.*log10(pix2sec),'o')
	#plt.plot(z,-2.5*log10(I_res)+5.*log10(pix2sec)+m0,'o')
	#ylim(24,16)
	#plt.show()

	z2 = []
	del z
	z = []

	for k in range(len(I_res)):
		if I_res[k]>0.:
			I_bulge1.append(I_res[k])
			z.append(z1[k]-z_c)

	L_bulge = fabs(trapz(I_bulge1, fabs(z)))

	S = 0.
	#print L_bulge

	for k in range(len(z)):
		if fabs(z[k])==min(fabs(z)):
			zmin0 = k

	for k in range(int((len(z)-1)/2)):
		if S<=L_bulge/2.:
			S = fabs(z[zmin0+k]-z[zmin0-k])*(I_bulge1[zmin0+k]+I_bulge1[zmin0-k])/2.+S
			reb = ( fabs(z[zmin0+k]) + fabs(z[zmin0-k]) ) / 2.
			Ieb = I_bulge1[zmin0+k]
	
	for k in range(len(I_bulge1)):
		if (z[k]<zmin1 or z[k]>zmin2)  and fabs(z[k])<zmax:
			I_bulge2.append(I_bulge1[k])
			z2.append(z[k])
	#print reb, -2.5*log10(Ieb)+m0+5.*log10(pix2sec)

	#exit() 
	dataToFit = RealData(z2, I_bulge2)
	mo = Model(bulge_func)
	fitting = ODR(dataToFit, mo, [Ieb,reb,3.,0.])
	fitting.set_job()
	fit_result = fitting.run()
	a = fit_result.beta[0]
	b = fit_result.beta[1]
	c = fit_result.beta[2]
	d = fit_result.beta[3]
	eps = fabs(median(fit_result.eps))
	if c>5.:	
		a = Ieb
		b = reb
		c = 3.
	#plt.plot(z2,I_bulge2,'o',markersize=5,color='red')
	#plt.plot(z2,bulge_func([a,b,c,d], z2), 'r-', lw=2, color='red')
	#print -2.5*log10(a)+m0,b/pix2sec,c,d
	#plt.show()
	return -2.5*log10(a)+m0+5.*log10(pix2sec),b,c	


def bulge_wod(scidata,nx,xc,z_c,m0d,h,m0,pix2sec,backgr_dev):
	rmin = 1.			#!!!
	x,x1,I = r_cut(scidata,xc,z_c,pix2sec,m0,rmin,nx/2.,backgr_dev,plot=0)
	I_res = I - disk_edge_r([10**(0.4*(m0-m0d)),h,xc],x1)
	I_bulge1 = []
	I_bulge2 = []

	#plt.plot(z,-2.5*log10(I)+m0+5.*log10(pix2sec))
	#plt.plot(z,-2.5*log10(disk_edge([10**(0.4*(m0-m0d)),z0,z_c],z1))+m0+5.*log10(pix2sec))
	#plt.plot(z,-2.5*log10(I_res)+5.*log10(pix2sec))
	#plt.show()

	x2 = []
	del x
	x = []

	for k in range(len(I_res)):
		if I_res[k]>0.:
			I_bulge1.append(I_res[k])
			x.append(x1[k]-xc)

	L_bulge = fabs(trapz(I_bulge1, fabs(x)))

	S = 0.
	#print L_bulge

	for k in range(len(x)):
		if fabs(x[k])==min(fabs(x)):
			xmin0 = k

	for k in range(int((len(x)-1)/2)):
		if S<=L_bulge/2.:
			S = fabs(x[xmin0+k]-x[xmin0-k])*(I_bulge1[xmin0+k]+I_bulge1[xmin0-k])/2.+S
			reb = ( fabs(x[xmin0+k]) + fabs(x[xmin0-k]) ) / 2.
			Ieb = I_bulge1[xmin0+k]
	
	#print reb, -2.5*log10(Ieb)+m0+5.*log10(pix2sec)

	#exit() 
	dataToFit = RealData(x, I_bulge1)
	mo = Model(bulge_func)
	fitting = ODR(dataToFit, mo, [Ieb,reb,3.,0.])
	fitting.set_job()
	fit_result = fitting.run()
	a = fit_result.beta[0]
	b = fit_result.beta[1]
	c = fit_result.beta[2]
	d = fit_result.beta[3]
	eps = fabs(median(fit_result.eps))
	if c>5.:	
		a = Ieb
		b = reb
		c = 3.
	#plt.plot(z2,I_bulge2,'o',markersize=5,color='red')
	#plt.plot(z2,bulge_func([a,b,c,d], z2), 'r-', lw=2, color='red')
	#print -2.5*log10(a)+m0+5.*log10(pix2sec),b,c,d
	#plt.show()
	return -2.5*log10(a)+m0+5.*log10(pix2sec),b,c



def main_edge(file_gal_clean,xc,yc,R_kron,a_image,b_image,backgr_dev,pix2sec,m0,fwhm,model,ellip,C31,NIR):
	#print file_gal_clean,xc,yc,R_kron,a_image,b_image,backgr_dev,pix2sec,m0,fwhm,model,ellip
	# *** FUNCTION TO FIND INITIAL PARAMETERS OF EDGE-ON DISK GALAXIES ***
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	# INPUT:
	#	 xc,yc - pixel coordinates of the center from SExtractor catalogue
	#	 R_kron - KRON_RADIUS from SExtractor catalogue
	#	 a_image, b_image - A_IMAGE, B_IMAGE from SExtractor catalogue
	#	 pix2sec - pixel scale [arcsec/pix]
	#	 m0 - magnitude Zeropoint for the Poghson formula m = m0 + 2.5*log10(I) + 5.*log10(pix2sec), where I - surface brightness
	# OUTPUT:
	#	 (x_c,y_c,m0d,h,z0,rtr,meb,reb,n,incl)
	#	 Where x_c,y_c - refined pixel coordinates of the galaxy center 
	#	 m0d - central SB of the disk [mag arcsec^2]; h,z0,rtr - in pixels
	# 	 meb - effective SB of the bulge [mag arcsec^2]; reb - in pixels
	#	 incl - estimated inclination angle if found [degrees]

	sma = R_kron*a_image	#semimajor radius 
	smi = R_kron*b_image	#semiminor radius

	rmin_bulge = rmin_bulge_coeff*fwhm/pix2sec
	#rmin_bulge = rmin_bulge_coeff/pix2sec

	hdulist = pyfits.open(file_gal_clean)
	scidata = hdulist[0].data
	ny,nx = scidata.shape


	# STEP 1: Estimation of the center y-coordinate to compare with the SExtractor's one
	zc = []
	z00 = []
	s = [-4,-3,3,4]
	for k in s:			# !!!
		x = int(xc + sma*k/6)
		z,z1,I = z_cut(scidata,yc,x,m0,0.,0.,ny/2,backgr_dev,plot=0)
		I0d,I0d_e,z01,z0_e1,z_c = fitting(z1,I,max(I),z0_est(z,I),yc,plot=0)
		if z01>0.:
			zc.append(z_c)
			z00.append(z01)

	z_c = np.median(np.array(zc))
	z0 = np.min(np.array(z00))
	print 'z0 = %.2f +/- %.2f [arcsec]' % (z0*0.396,0.5) 

	


	#print z0*0.396
	#plt.show()
	#exit()

	del zc

	if fabs(yc-z_c)>2.:	yc = z_c
	z_c = yc
	#plt.show()
	#print yc
	# Now we have yc=[pix] and z_c=[arcsec]


	# STEP 2: Find the location of the dust lane (up the found center or down it) - zmin1, zmin2 and zmax
	z,z1,I = z_cut(scidata,yc,np.ceil(xc-0.5),m0,0.,0.,ny/2,backgr_dev,plot=0)
	I0d,I0d_e,z01,z0_e1,z_c1 = fitting(z1,I,max(I),z0,z_c,plot=0)

	I_res = fabs(I - disk_edge([I0d,z0,z_c],z1))

	z_Imax = z1[list(I).index(max(I))]-z_c
	z_Ires = z1[list(I_res).index(max(I_res))]-z_c

	if z_Imax<z_Ires and z_Ires>0:
		print 'Dust lane is upper the center'
		delta_dust = fabs(z_Ires - z_Imax)
		dust = 1
		zmin1 = z_Imax
		zmin2 = zmin1 + delta_dust 
	elif z_Imax>z_Ires and z_Ires<0:
		print 'Dust lane is lower the center'
		delta_dust = fabs(z_Imax - z_Ires)
		dust = -1
		zmin2 = z_Imax
		zmin1 = zmin2 - delta_dust 
	elif fabs(z_Imax)<2. and fabs(z_Ires)<2.:
		print 'Dust lane is at the center or is not detected'
		dust = 0
		delta_dust = fabs(2.*(z0/4.))		# !!!
		zmin1 = - delta_dust/2.
		zmin2 = + delta_dust/2.
	else:
		dust = 0
		delta_dust = 0.		# !!!
		zmin1 = 0.
		zmin2 = 0.

	#print 'zmin',zmin1,zmin2
	

	'''
	def func(x):
		return disk_edge([I0d,z0,0.],x) - 3.*backgr_dev		# !!! 3sigma
	zmax = sp.optimize.bisect(func,0,ny/2)
	print zmax
	'''
	zmax = smi+smi/3.		#!!! Other way
	# Now we have zmin1,zmin2 and zmax [arcsec]!
	#print zmin1,zmin2,zmax


	# STEP 3: Averaging the radial cuts
	i = 0
	#print np.ceil(z_c - 0.5)
	s = [-z0,-z0/2,z0/2,z0]

	for m in s:
			y = int(np.ceil(z_c - 0.5 + m))
			if i==0:
				x,x1,I = r_cut(scidata,xc,y,pix2sec,m0,0,a_image*R_kron,-1000.,plot=0)
				xx = x1
				II = I
				i = i + 1		
			else:
				x,x1,I = r_cut(scidata,xc,y,pix2sec,m0,0,a_image*R_kron,-1000.,plot=0)
				II = np.vstack([II,I])
	#plt.show()

	size=II.size
	ny,nx=II.shape
	#print nx
	III = np.reshape(II,size,order='F').reshape((nx,ny))	
	#print III[0]
	Intensity = []
	for k in range(nx):
		Intensity.append(np.median(III[k]))

	xc = xx[list(Intensity).index(max(Intensity))]	
	#plt.plot(xx-xc, m0-2.5*log10(Intensity),'r-',color='red', lw=3)
	#plt.plot(xx-xc, Intensity,'r-',color='red', lw=3)
	#plt.show()
	#print xc

	rmin = smi

	xxx= []
	III = []
	for k in range(len(xx)):
		if xx[k]>xc-smi and xx[k]<xc+smi:
			xxx.append(xx[k])
			III.append(Intensity[k])


	#plt.show()

	def Igau(B,x):
		#B[0]=I0
		#B[1]=fwhm
		#B[2]=xc
		return B[0]*np.exp(-fabs(x-B[2])**2/(2.*(B[1]/2.35482)**2))
	dataToFit = RealData(xxx, III)
	mo = Model(disk_edge)

	fitt = ODR(dataToFit, mo, [max(III),z0,xc])
	fitt.set_job()
	fit_result = fitt.run()
	A = fit_result.beta[0]
	B = fit_result.beta[1]
	xc_new = fit_result.beta[2]
	xc = xc_new
	#print xc_new
	#plt.plot(xxx-xc_new, III,'r-',color='red', lw=3)
	#plt.plot(xxx-xc_new, disk_edge([A,B,xc_new],xxx),'r-',color='blue', lw=3)
	#plt.show()
	# Now we have xc



	# STEP 4: Finding h
	del II
	r = []
	Int = []
	#print rmin*pix2sec
	#exit()
	for k in range(len(xx)):
		if Intensity[k]>3.*backgr_dev and fabs(xx[k]-xc)>rmin:
			r.append(fabs(xx[k]-xc))
			Int.append(Intensity[k])

	dataToFit = RealData(np.array(r)*pix2sec, m0-2.5*log10(Int) + 5.*log10(pix2sec))
	mo = Model(line)
	fitt = ODR(dataToFit, mo, [1.,1.])
	fitt.set_job()
	fit_result = fitt.run()
	h = fit_result.beta[0]
	h_e = fit_result.sd_beta[0]
	m0d = fit_result.beta[1]
	print 'h = %.2f +/- %.2f [arcsec]' % (h,h_e) 
	#plt.plot(np.array(xx-xc)*pix2sec, m0-2.5*log10(Intensity)+ 5.*log10(pix2sec))
	#plt.plot(r,line([h,m0d],r))
	#plt.show()
	h = h/pix2sec 	# now h is in pixels

	h = h/1.27	# !!! Trend!!!

	#exit()
	# Now we have h

	try:
		for k in range(len(r)):
			if r[k]>sma/2. and fabs( 10**(0.4*(m0-line([h,m0d],r[k]))) - Int[k] )>=1.*backgr_dev:	#!!!
				rmax = r[k]
	except:
		rmax = R_kron*a_image

	rmax = R_kron*a_image		#!!!!

	#print 'rtr = %.2f [arcsec]' % (rmax) 
	#plt.show()

	npix = int(np.floor(rmax-rmin))
	rmax_new = (int(np.floor(rmax)) - int(npix/2))

	'''
	# STEP 5: Finding z0 in the region outside the bulge and dust
	step = 15
	I0 = []
	Z0 = []
	xx = []


	npix = int(np.floor(rmax-rmin))
	rmax_new = (int(np.floor(rmax)) - int(npix/2))
	x_good = []
	s = [int(xc-3.*h),int(xc-2.*h),int(xc-h),int(xc+h),int(xc+2.*h),int(xc+3.*h)]

	for x in s:
		#if (x<(xc-rmin) and x > (xc-rmax_new)) or (x>(xc+rmin) and x < (xc+rmax)):
		x_good.append(x)

	i=0
	Z0 = []
	for k in range(len(x_good)):
		x = x_good[k]
		z,z1,I = z_cut(scidata,yc,x,m0,zmin1,zmin2,zmax,backgr_dev,plot=0)
		good_pix = 0
		for m in range(len(I)):
			if I[m]!=0:	good_pix = good_pix + 1
		if good_pix==len(I):
			I0d,I0d_e,z0,z0_e,z_c = fitting(z1,I,max(I),h/5.,z_c)
			Z0.append(z0)
			i = i + 1
	z0 = np.median(Z0)
	z0_err = np.std(Z0)
	print 'z0 = %.2f +/- %.2f [arcsec]' % (z0,z0_err) 
	del x_good
	del x
	del y
	exit()
	'''
	# STEP 6: Finding m0d of the disk in the region outside the bulge and dust
	m0d = []
	x_good = []
	y_good = []
	step1 = int(ceil(2*zmax/(100)))
	step2 = int(ceil(2*rmax/(100)))
	ymax_l = int(yc-zmax)
	ymax_r = int(yc+zmax)
	xmax_l = int(xc-rmax_new)
	xmax_r = int(xc+rmax_new)
	ymin1 = int(yc+zmin1)
	ymin2 = int(yc+zmin2)
	xmin1 = int(xc-rmin)
	xmin2 = int(xc+rmin)
	#print ymax_l,ymin1,ymin2,ymax_r
	#print xmax_l,xmin1,xmin2,xmax_r

	#exit()

	for y in [int(yc-z0),int(yc-z0/2),int(yc+z0/2),int(yc+z0)]:
		for x in [int(xc-3.*h),int(xc-2.*h),int(xc-h),int(xc+h),int(xc+2.*h),int(xc+3.*h)]:
				#if (y<yc+zmin1 or y>yc+zmin2) and (x<xc-rmin or x>xc+rmin) and scidata[y,x]!=0:
				x_good.append(x)
				y_good.append(y)
				I_obs = scidata[y,x] 
				r = fabs(xc - x)
				z = fabs(yc - y)
				L_s = I_obs/(2.*r*sp.special.kn(1,r/h)/((np.cosh(z/z0))**2))
				I0 = 2.*L_s*h
				m0d.append(-2.5*log10(I0) + m0 )

	m0_d = np.median(m0d)+5.*log10(pix2sec)
	m0_d_e = np.std(m0d)

	print 'S0d = %.2f +/- %.2f [mag arcsec^-2]' % (m0_d-2.5*log10(z0/h),m0_d_e) 
	#exit()

	if model=='edge+ser':
		# STEP 7: Bulge parameters
		#if delta_dust <= fabs(2.*(z0/4.)):
		#meb,reb,n = bulge(scidata,xc,yc,zmin1,zmin2,zmax,m0_d-5.*log10(pix2sec),z0,m0,pix2sec,backgr_dev)
		#NIR = 0
		#print 'zmax=%.3f' % (zmax)

		if NIR==0:
			try:
				r,I,r1,I1, = bulge_profile_minor(scidata,xc,z_c,zmin1,zmin2,zmax,m0_d-5.*log10(pix2sec),z0,m0,pix2sec,backgr_dev,rmin_bulge,dust)
				m0d,z0,meb,reb,n = bulge_fitting1(r,I,r1,I1,m0,pix2sec,ellip,C31,m0_d,z0,h)
			except:
				r,I,r1,I1 = bulge_profile_major(scidata,xc,z_c,zmin1,zmin2,zmax,m0_d-5.*log10(pix2sec),h,m0,pix2sec,backgr_dev,rmin_bulge,dust)
				m0d,h,meb,reb,n = bulge_fitting2(r,I,r1,I1,m0,pix2sec,ellip,C31,m0_d,h,z0)


		else:
			try:
				print 'yes!'
				r,I,r1,I1 = bulge_profile_major(scidata,xc,z_c,zmin1,zmin2,zmax,m0_d-5.*log10(pix2sec),h,m0,pix2sec,backgr_dev,rmin_bulge,dust)
				m0d,h,meb,reb,n = bulge_fitting2(r,I,r1,I1,m0,pix2sec,ellip,C31,m0_d,h,z0)
			except:
				print 'YES!'
				r,I,r1,I1, = bulge_profile_minor(scidata,xc,z_c,zmin1,zmin2,zmax,m0_d-5.*log10(pix2sec),z0,m0,pix2sec,backgr_dev,rmin_bulge,dust)
				m0d,z0,meb,reb,n = bulge_fitting1(r,I,r1,I1,m0,pix2sec,ellip,C31,m0_d,z0,h)

	else:
		meb = 99999.
		reb = 99999.
		n = 99999. 
		m0d = m0_d
	#print 'meb = %.2f [mag arcsec^-2], reb = %.2f [arcsec], n = %.2f' % (meb,reb,n) 
	if dust==-1 or dust==1:	incl = 90.-degrees(np.arctan(delta_dust/(rmax*2.)))
	else:	incl = 90.

	cd=-1.
	import ell
	try:
		ellb,bb4 = ell.read_ell('ellipse.txt',reb,1.)
		if bb4>0.01 and bb4<1.:
			cb = 0.5
		else:
			cb = 0.

	except:
		ellb = 0.
		cb = 0.


	if model=='edge':
		ellb=99999.
		cb=99999.

	elld=99999.

	import rmax_finder
	try:
		rmax = rmax_find_edge(r,I,h,m0d,m0,pix2sec,NOISE)
	except:
		rmax = 4.*h
	
	#zmin1 = int(round(yc - fabs(zmin1)))
	#zmin2 = int(round(yc + fabs(zmin2))) 

	zmin1 = int(math.floor(yc - fabs(zmin1)))-1
	zmin2 = int(math.ceil(yc + fabs(zmin2)))+1 

	#return m0d,99999.,h,z0,elld,cd,rmax,meb,99999.,reb,n,ellb,cb,xc,yc,incl
	return m0d,99999.,h,z0,elld,cd,rmax,meb,99999.,reb,n,ellb,cb,xc,yc,incl,zmin1,zmin2


#xc,yc,m0d,h,z0,rtr,meb,reb,n,incl = main_edge('galaxy_clean.fits',2190,430,3.5,517,38.,10.,0.4,24.04)
#xc,yc,m0d,h,z0,rtr,meb,reb,n,incl = main_edge('galaxy_clean.fits',99,49,3.5,28.,8.,0.,0.396,26.0)
#print xc,yc,m0d,h,z0,rtr,meb,reb,n,incl
#print xc,yc,m0d,h,z0,rtr,meb,reb,n,incl
#main_edge('galaxy_clean.fits',300.,300.,4.,12.,4.95,7.78,0.396,27.9,1.244,'edge+ser',4.06)




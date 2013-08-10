#!/usr/bin/python
# Script to find warp parameters
# -*- coding:  cp1251 -*-

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
from kapteyn import wcs
from scipy import special
from scipy.odr.odrpack import *

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')


fsize = 25
fig_size = (10,10)

def left_warp(file_image,xmin,ymin,xmax,ymax):
	# (x1,y1) - bottom left corner
	# (x2,y2) - top right corner
	import sky_est

	sky_est.imcopy_func(file_image,'left.fits',xmin,ymin,xmax,ymax)

	hdulist = pyfits.open('left.fits')
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']

	x = []; y = []

	backgr = median(scidata)

	y1 = []; Imax = []
	yyy = []
	for k in range(nx):
		I1 = []
		for i in range(ny):
			I1.append(scidata[i,k])
			yyy.append(i)
		#y1.append(I1.index(max(I1)))
		#Imax.append(max(I1))
		yy,II = y_int(yyy,I1)
		y1.append(yy)
		Imax.append(II)

	for k in range(nx):
		if Imax[k]>backgr:
			#I.append(Imax[k])
			x.append(k)
			y.append(y1[k])

	#print np.array(x),np.array(y)
	#exit()

	return xmin+np.array(x),ymin+np.array(y)



def right_warp(x,y,xc,yc,rmax,side,number):
	X_l1,Y_l1,X_r1,Y_r1,X_l2,Y_l2,X_r2,Y_r2 = loadtxt('warps_coords.txt', usecols=[0,1,2,3,4,5,6,7], unpack=True, skiprows = 1)
	try:
		x_l1=X_l1[number-1]; y_l1=Y_l1[number-1]; x_r1=X_r1[number-1]; y_r1=Y_r1[number-1]
		x_l2=X_l2[number-1]; y_l2=Y_l2[number-1]; x_r2=X_r2[number-1]; y_r2=Y_r2[number-1]
	except:
		x_l1=X_l1; y_l1=Y_l1; x_r1=X_r1; y_r1=Y_r1
		x_l2=X_l2; y_l2=Y_l2; x_r2=X_r2; y_r2=Y_r2

	if side=='left':
		#bb = fabs((xc-x_l2)/(x_l2-x_l1))
		bb = fabs((x_l2-min(x))/(x_l1-x_l2))
		print 'bb_l=%.1f' % (bb)


		#x_l=110;y_l=39
		#x_r=165;y_r=63

		x1,y1 = left_warp('galaxy_clean.fits',x_l1,y_l1,x_r1,y_r1)


		x1 = list(fabs(np.array(x1)-xc))
		y1 = list(np.array(y1)-yc)

		x_left2 = x1
		y_left2 = []
		x_left3 = sorted(x_left2, key=float)
		x_left2 = x_left3
		#print x_left2

		for k in range(len(x_left2)):
			y_left2.append(y1[x1.index(x_left2[k])])

		x1 = x_left2
		y1 = y_left2
		
	else:
		bb = fabs((x_r1-min(x))/(x_r2-x_r1))
		print 'bb_r=%.1f' % (bb)
		#x_l=462;y_l=29
		#x_r=505;y_r=50

		x1,y1 = left_warp('galaxy_clean.fits',x_l2,y_l2,x_r2,y_r2)

		x1 = list(fabs(np.array(x1)-xc))
		y1 = list(np.array(y1)-yc)	



	tck = interpolate.splrep(np.array(x1),np.array(y1),s=0)
	x_interp = np.arange(min(x1),max(x1),fabs(x1[0]-x1[1])/bb)
	y_interp = interpolate.splev(x_interp,tck,der=0)

	'''
	med = mean(y_interp)
	stand = std(y_interp)

	x2 = []
	y2 = []
	for k in range(len(y_interp)):
		if fabs(y_interp[k]-med)<1.*stand:
			x2.append(x_interp[k])
			y2.append(y_interp[k])			
	'''
	x2 = x_interp
	y2 = y_interp

	x_2 = []
	y_2 = []
	for k in range(len(x)):
		if x[k]>0. and x[k]<min(x1):
			x_2.append(x[k])
			y_2.append(y[k])

	xx =  x_2 + list(x2)
	yy =  y_2 + list(y2)
	return xx,yy

def fitting(x,y,slope):
		#plt.plot(x,y)
		#plt.show()
		dataToFit = RealData(np.array(x), np.array(y))
		#print x,y
		mo = Model(warp_fit)
		fitt = ODR(dataToFit, mo, [30.,slope])
		fitt.set_job()
		fit_result = fitt.run()
		A = fit_result.beta[0]
		C = fit_result.beta[1]
		return A,C

def Heavyside(x,A):
	y = []
	for k in range(len(x)):
		if fabs(x[k])<fabs(A):
			y.append(0.)
		elif fabs(x[k])>=fabs(A):
			y.append(1.)
	return np.array(y)

def warp_fit(B, x):
	#B[0]=A
	#B[1]=C
	return  Heavyside(x,B[0])*B[1]*(fabs(x)-B[0])

def y_int(y,I):
	Isum = sum(fabs(I))
	I1 = 0.
	k = 0
	#print Isum
	while I1<=Isum/2.:
		I1 = I1 + fabs(I[k])
		y1 = y[k]
		k = k + 1
		#print y1,I1
	return y1,I1

def warp_pars(scidata,ff,x,y,xc,yc,pix2sec,m0,backgr,h,z0,rmax,number):
	#plt.figure(27)
	#plt.plot(x,y,'o')
	#plt.show()
	#exit()


	omega = 1./(max(x)**3.) * sum(x*y)
	omega_r = 0.
	omega_l = 0.	
	for k in range(len(x)):
		if x[k]>=0.:
			omega_r = omega_r + x[k]*y[k]
		if x[k]<0.:
			omega_l = omega_l + x[k]*y[k]
	omega_l = omega_l/(4.*(max(-x)**3.))
	omega_r = omega_r/(4.*(max(x)**3.))
	alfa_s = fabs(omega_r-omega_l)/(omega_r+omega_l)

	r_last_right = max(x)
	r_last_left = max(-x)
	y_last_right = y[list(x).index(max(list(x)))]
	y_last_left = y[list(-x).index(max(list(-x)))]

	if y_last_right>=0.:	slope_right = 0.5
 	if y_last_right<0.:	slope_right = -0.5
	if y_last_left>=0.:	slope_left = 0.5
 	if y_last_left<0.:	slope_left = -0.5


	psi = max( degrees(arctan(y_last_right/r_last_right)),degrees(arctan(y_last_left/r_last_left))) 
	print psi




	x_left1 = []; y_left1 = []; x_right1 = []; y_right1 = []
	for k in range(len(x)):
		if x[k]<0.:
			x_left1.append(-x[k])
			y_left1.append(y[k])
		else:
			x_right1.append(x[k])
			y_right1.append(y[k])
	#plt.figure(25)
	#plt.plot(x_left1,y_left1,'o')
	#plt.show()
	#exit()

	#for k in range(len(x_right1)):
	#	print x_right1[k],y_right1[k]

	#exit()

	x_left2 = x_left1
	y_left2 = []
	x_left3 = sorted(x_left2, key=float)
	x_left2 = x_left3
	#print x_left2



	
	for k in range(len(x_left2)):
		y_left2.append(y_left1[x_left1.index(x_left2[k])])
	'''
	for k in range(len(x_right1)):
		print x_left2[k],y_left2[k]
	'''

	x_left,y_left = right_warp(x_left2,y_left2,xc,yc,rmax,'left',number)
	x_right,y_right = right_warp(x_right1,y_right1,xc,yc,rmax,'right',number)
	#print x_left,y_left
	#exit()


	x_left = np.array(x_left)
	y_left = np.array(y_left)
	x_right = np.array(x_right)
	y_right = np.array(y_right)

	plt.figure(25)
	plt.plot(x_right,y_right,'o')
	plt.plot(x_left,y_left,'o')
	plt.plot(x_left2,y_left2,'x')
	plt.plot(x_right1,y_right1,'x')
	plt.savefig("warps.png", transparent = False)

	#exit()
	#yy = warp_fit(x_left,[70,0.7])
	#exit()	



	A_left,C_left = fitting(fabs(x_left),y_left,slope_left)
	print A_left,-C_left

	A_right,C_right = fitting(x_right,y_right,slope_right)
	print A_right,C_right



	'''
	plt.figure(5)
	plt.plot(x_right,warp_fit([A_right,C_right], x_right))
	plt.plot(x_right,y_right)
	plt.plot(-x_left,warp_fit([A_left,C_left], -x_left))
	plt.plot(-x_left,y_left)
	plt.show()
	'''


	#plt.figure(12)
	ny,nx = scidata.shape
	zmin = -4.*z0*pix2sec #min(y*pix2sec)-1. 
	zmax = 4.*z0*pix2sec #max(y*pix2sec)+1.
	rmin = -4.*h*pix2sec #min(x*pix2sec)-1.
	rmax = 4.*h*pix2sec #max(x*pix2sec)+1.


	f = plt.figure(50,figsize=(fig_size[0],fig_size[0]*0.1+4.5))
	# Center yc(r)
	#plt.axvline(rmin,linestyle='--',color='blue')
	#plt.axvline(-rmin,linestyle='--',color='blue')
	plt.xlabel(r'r (arcsec)', fontsize=fsize)
	plt.ylabel(r'z (arcsec)', fontsize=fsize)
	plt.plot(x*pix2sec,y*pix2sec,'o',color='black')
	plt.plot(x_right*pix2sec,warp_fit([A_right,C_right], x_right)*pix2sec,lw=2,color='red')
	plt.plot(-x_left*pix2sec,warp_fit([A_left,C_left], -x_left)*pix2sec,lw=2,color='red')



	ylim(zmin,zmax)
	xlim(rmin,rmax)


	xx = []
	yy = []
	II = []


	for k in range(ny):
		for i in range(nx):
			if scidata[k,i]>backgr:
				xx.append((i + 0.5-xc)*pix2sec)
				yy.append((k + 0.5-yc)*pix2sec)
				II.append(scidata[k,i])
	xx = np.array(xx)
	yy = np.array(yy)
	II = np.array(II)

	m = -2.5*log10(II) + m0 + 5.*log10(pix2sec)
	m_crit = -2.5*log10(backgr) + m0 + 5.*log10(pix2sec)
	m = np.around(m,decimals=1)


	sample_pts = 500
	con_levels = np.arange(min(m),m_crit,0.5)
	#levels = np.arange(16,20,0.1)

	xi = np.linspace(rmin,rmax,sample_pts)
	yi = np.linspace(zmin,zmax,sample_pts)

	Ii = griddata(xx,yy,m,xi,yi)
	#print rmax

	#f = plt.figure()
	#ax1 = f.add_subplot(111)
	plt.contour(xi,yi,Ii,con_levels,linewidths=1,colors='black')
	#ax1.set_xlim(-rmax,rmax)
	#ax1.set_ylim(-zmax,zmax)
	#ax1.set_ylabel(r'z (arcsec)', fontsize=fsize)

	plt.savefig("yc_r.png", transparent = False)



	print 'Warp(s) parameters:'
	print '\tomega=%.2f, omega_r=%2.f, omega_l=%.2f, alfa_s=%.2f,r_last_right=%.2f,r_last_left=%.2f,y_last_right=%.2f,y_last_left=%.2f,psi=%.2f,A_right=%.2f,A_left=%.2f,C_right=%.2f,C_left=%.2f' % (omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,A_left,C_right,C_left)

	print >> ff, 'Warp(s) parameters:'
	print >> ff, '\tomega=%.2f, omega_r=%2.f, omega_l=%.2f, alfa_s=%.2f,r_last_right=%.2f,r_last_left=%.2f,y_last_right=%.2f,y_last_left=%.2f,psi=%.2f,A_right=%.2f,A_left=%.2f,C_right=%.2f,C_left=%.2f' % (omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,A_left,C_right,C_left)

	return omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left
'''
x,y = loadtxt('warp.txt', usecols=[0,1], unpack=True)
omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left = warp_pars(x,y)
print omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left
'''


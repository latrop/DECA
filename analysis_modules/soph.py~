#!/usr/bin/python
# Module to make vertical cuts for rotated edge-on galaxy
# You should first find rmin and rmax to do this
#!/usr/bin/python
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
import INITIALS_edge_on
import median_fil
import plots
import warp
import ell

file_out_edge = setup.file_out_edge

tmp_out = sys.stdout

#rmin_bulge_coeff = setup.rmin_bulge_coeff

fsize = 25

def line(B, x):
	#*** For non edge-on disk SB in mag/arcsec^2 (radial) ***
	# B[0]=h
	# B[1]=m0d
	x = np.array(x)
	return (1.0857/B[0])*fabs(x) + B[1]

def line1(B, x):
	#*** For non edge-on disk SB in mag/arcsec^2 (radial) ***
	# B[0]=h
	# B[1]=m0d
	# B[2]=x_c
	x = np.array(x)
	return (1.0857/B[0])*fabs(x-B[2]) + B[1]

def disk_edge_r(B,x):
	#*** For edge-on disk SB in mag/arcsec^2 (along r-axis). ***
	# B[0] = I0d
	# B[1] = h
	# B[2] = x_c
	x = np.array(x)
	return B[0] * fabs(x-B[2])/B[1] * sp.special.kn(1,fabs(x-B[2])/B[1])

def Iedge_i(r,h,I0d):
	# Radial distribution for the edge-on disc
	r = fabs(r)
	I = I0d * r/h * sp.special.kn(1,r/h)
	return I

def disk_edge_soph(B, z):
	#*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
	# B[0] = I0d
	# B[1] = n
	# B[2] = z0 
	# B[3] = z_c
	z = np.array(z)	
	return B[0] * (1.0 / np.cosh(B[1]*fabs(z-B[3])/(2.*B[2])) ) **(2./B[1])

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

def fitting(r,mag,plot=0):
	# All is given in arcsec
	r = np.array(r)
	dataToFit = RealData(r, mag)
	mo = Model(line)
	fitting = ODR(dataToFit, mo, [1,1])
	fitting.set_job()
	fit_result = fitting.run()
	m0d = fit_result.beta[0]
	m0d_e = fit_result.sd_beta[0]
	h = fit_result.beta[1]
	h_e = fit_result.sd_beta[1]

	if plot==1:
		plt.plot(r, line([m0d,h],r),'r--',color='blue', lw=2)
	return m0d,h



def breaks(r,I,rmin,rmax,m0,pix2sec):
	plt.plot(r,-2.5*log10(I)+m0+5.*log10(pix2sec),'o',color='black',markersize=5)
	nx = rmax - rmin
	rmin = rmin + nx/6.
	rmax = rmax - nx/6.
	nx = rmax - rmin

	r1 = []
	mag1 = []
	r2 = []
	mag2 = []
	r_tr = []
	mag_tr = []
	for k in range(len(r)):
		if r[k]>rmin and r[k]<rmin+nx/4.:
			r1.append(r[k]*pix2sec)
			mag1.append(-2.5*log10(I[k])+m0+5.*log10(pix2sec))
	
	m0d1,h1 = fitting(r1,mag1,plot=1)

	for k in range(len(r)):
		if r[k]>rmax-nx/4. and r[k]<rmax:
			r2.append(r[k]*pix2sec)
			mag2.append(-2.5*log10(I[k])+m0+5.*log10(pix2sec))

	m0d2,h2 = fitting(r2,mag2,plot=1)

	r_B1 = (m0d2-m0d1)/(1.0857*(1./h1-1./h2))
	plt.axvline(r_B1,linestyle='--',color='black')

	if r_B1>rmax:
		r3 = []
		mag3 = []
		# 2 breaks!
		for k in range(len(r)):
			if r[k]>rmin + 3.*nx/4. and r[k]<rmax-3.*nx/4.:
				r3.append(r[k]*pix2sec)
				mag3.append(-2.5*log10(I[k])+m0+5.*log10(pix2sec))
			
		m0d3,h3 = fitting(r3,mag3,plot=1)
		r_B1 = (m0d3-m0d1)/(1.0857*(1./h1-1./h3))
		r_B2 = (m0d2-m0d3)/(1.0857*(1./h3-1./h2))
		plt.axvline(r_B1,linestyle='--',color='red')
		plt.axvline(r_B2,linestyle='--',color='blue')
		
	for k in range(len(r)):
		if r[k]>rmax-5. and r[k]<rmax:			#!!!
			r_tr.append(r[k]*pix2sec)
			mag_tr.append(-2.5*log10(I[k])+m0+5.*log10(pix2sec))

	m0d_tr,h_tr = fitting(r_tr,mag_tr,plot=1)
	r_Btr = (m0d_tr-m0d2)/(1.0857*(1./h2-1./h_tr))
	plt.axvline(r_Btr,linestyle='--',color='black')
	plt.show()

def slit(scidata,xc,yc,z0,rmax,m0,pix2sec):
	ny,nx = np.shape(scidata)
	#rmax = nx/2
	I_eo = []
	r = []
	for x in range(int(xc-rmax),int(xc+rmax)):
		z,z1,I = INITIALS_edge_on.z_cut(scidata,yc,x,m0,z0/2,-z0/2,ny/2,backgr,plot=0)
		I_eo.append(mean(I))
		r.append(x)
	xc_new = r[I_eo.index(max(I_eo))]
	r = np.array(r)-xc_new
	#plt.plot(r,-2.5*log10(I_eo)+m0+5.*log10(pix2sec))
	#plt.show()
	#exit()



	I_right = []
	r_right = []
	for k in range(len(r)):
		if r[k]>=0:
			I_right.append(I_eo[k])
			r_right.append(r[k])

	#plt.plot(r_right,-2.5*log10(I_right)+m0+5.*log10(pix2sec))
	#plt.show()
	#exit()

	b = 1.03
	x = 0
	n = 0
	x_uzl = []
	while x<rmax:
		x = x + b**n
		x_uzl.append(x)
		n = n + 1


	r_bin_right=[]
	I_bin_right=[]
	for i in range(len(x_uzl)-1):
		I = 0
		n = 0
		for k in range(len(r_right)):
			if r_right[k]<=x_uzl[i+1] and r_right[k]>=x_uzl[i]:
				 I = I + I_right[k]
				 n = n + 1
		if n==0:
			r_bin_right.append(1)
			I_bin_right.append(I_right[k])
		else:
			r_bin_right.append(x_uzl[i]+(x_uzl[i+1]-x_uzl[i])/2.)
			I_bin_right.append(I/n)
	#plt.plot(r_bin_right,-2.5*log10(I_bin_right)+m0+5.*log10(pix2sec), 'o', color='red',markersize=5)


	I_left = []
	r_left = []
	for k in range(len(r)):
		if r[k]<0:
			I_left.append(I_eo[k])
			r_left.append(fabs(r[k]))

	#plt.plot(r_left,-2.5*log10(I_left)+m0+5.*log10(pix2sec))
	#plt.show()
	#exit()

	b = 1.03
	x = 0
	n = 0
	x_uzl = []
	while x<rmax:
		x = x + b**n
		x_uzl.append(x)
		n = n + 1


	r_bin_left=[]
	I_bin_left=[]
	for i in range(len(x_uzl)-1):
		I = 0
		n = 0
		for k in range(len(r_left)):
			if r_left[k]<=x_uzl[i+1] and r_left[k]>=x_uzl[i]:
				 I = I + I_left[k]
				 n = n + 1
		if n==0:
			r_bin_left.append(1)
			I_bin_left.append(I_left[k])
		else:
			r_bin_left.append(x_uzl[i]+(x_uzl[i+1]-x_uzl[i])/2.)
			I_bin_left.append(I/n)
	#plt.plot(r_bin_left,-2.5*log10(I_bin_left)+m0+5.*log10(pix2sec), 'o', color='green',markersize=3)

	rr = []
	Int = []
	for k in range(len(r_bin_left)):
		for m in range(len(r_bin_right)):
			if r_bin_left[k]==r_bin_right[m]:
				Int.append((I_bin_left[k]+I_bin_right[m])/2.)
				rr.append(r_bin_right[m])

	plt.plot(rr,-2.5*log10(Int)+m0+5.*log10(pix2sec), 'v', color='black',markersize=6)

	#plt.show()
	return rr,Int

def z0_est(z,I):
	# Function to estimate z0 for sech^2 law
	I_z0 = max(I)/2.3811 
	for k in range(len(z)):
		if I[k]>I_z0:
			z0 = fabs(z[k])
	return z0




def cuts_soph(scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z01,h1,m0d1,m0,pix2sec,backgr):
	ny,nx = np.shape(scidata)
	I_eo = []
	r = []
	#z0 = []
	#z_c = []
	I0d = []
	x_new = []
	n = []
	rmax = 3.*h1

	zc = []
	z00 = []

	s = range(int(xc-rmax),int(xc-rmin)) + range(int(xc+rmin),int(xc+rmax))

	for x in s:
		z,z1,I = INITIALS_edge_on.z_cut(scidata,yc,x,m0,0.,0.,ny/2,backgr,plot=0)
		try:
			I0d1,I0d_e1,n1,n_e1,z0,z0_e1,z_c = INITIALS_edge_on.fitting_soph(z1,I,max(I),1.,z01,z1[list(I).index(max(I))],plot=0)
			if I0d1>0. and z0<100 and z0>0.:
				z00.append(z0)
				zc.append(z_c)
				I0d.append(I0d1)
				x_new.append(x)
				n.append(n1)
		except:
			no = 1
	z0 = np.median(z00)
	yc = np.median(zc)
	N = np.median(n)


	plt.figure(15)
	# Height z0(r)
	plt.plot((np.array(x_new)-xc)*pix2sec,pix2sec*np.array(z00),'o',color='black')
	#plt.plot((np.array(x_new_bul)-xc)*pix2sec,pix2sec*np.array(z0_bul),'o',color='black')
	plt.axvline(rmin*pix2sec,linestyle='--',color='blue')
	plt.axvline(-rmin*pix2sec,linestyle='--',color='blue')
	plt.xlabel(r'r (arcsec)', fontsize=fsize)
	plt.ylabel(r'z$_0$ (arcsec)', fontsize=fsize)
	plt.savefig("z0_r_sph.png", transparent = False)

	plt.figure(16)
	# Center yc(r)
	plt.plot((np.array(x_new)-xc)*pix2sec,zc,'o',color='black')
	#plt.plot(np.array(x_new_bul)-xc,z_c_bul,'o',color='black')
	plt.axvline(rmin,linestyle='--',color='blue')
	plt.axvline(-rmin,linestyle='--',color='blue')
	plt.xlabel(r'r (pix)', fontsize=fsize)
	plt.ylabel(r'y$_c$ (pix)', fontsize=fsize)
	plt.savefig("yc_r_sph.png", transparent = False)

	plt.figure(17)
	# Center n(r)
	plt.plot((np.array(x_new)-xc)*pix2sec,n,'o',color='black')
	#plt.plot(np.array(x_new_bul)-xc,n_bul,'o',color='black')
	plt.axvline(rmin,linestyle='--',color='blue')
	plt.axvline(-rmin,linestyle='--',color='blue')
	plt.xlabel(r'r (pix)', fontsize=fsize)
	plt.ylabel(r'n$_{disk}$', fontsize=fsize)
	plt.ylim(0,5)
	plt.savefig("n_r_sph.png", transparent = False)

	return yc,z0,N







def cuts(mod,scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z01,h1,m0d1,m0,pix2sec,backgr):
	'''	
	if m0d1==nan:
		m0d1 = -2.5*log10(max(scidata)) + m0 + 5.*log10(pix2sec)
		print 'ok'
	else:
		print 'not ok'
	'''
	ny,nx = np.shape(scidata)
	I_eo = []
	r = []
	#z0 = []
	#z_c = []
	I0d = []
	x_new = []
	rmax = 4.*h1

	zc = []
	z00 = []

	s = range(int(xc-rmax),int(xc-rmin)) + range(int(xc+rmin),int(xc+rmax))
	#s = [int(xc-3.*h1)]
	for x in s:			# !!!
		#x = int(xc + rmax*k/6)
		z,z1,I = INITIALS_edge_on.z_cut(scidata,yc,x,m0,0.,0.,ny/2,backgr,plot=0)
		try:
			I0d1,I0d_e,z0,z0_e1,z_c = INITIALS_edge_on.fitting(z1,I,max(I),z01,z1[list(I).index(max(I))],plot=0)
			if I0d1>0. and z0<100 and z0>0.:
				zc.append(z_c)
				z00.append(z0)
				I0d.append(I0d1)
				x_new.append(x)
		except:
			no = 0

	z0 = np.median(np.array(z00))
	yc = np.median(np.array(zc))

	#print 'z0 = %.2f +/- %.2f [arcsec]' % (z0*0.396,0.5)
	#print median(zc)





	#m0d1 = m0d1 + 2.5*log10(z01/h1)
	#print m0d1
	#exit()

	r = np.array(x_new) - xc

	rr = []
	I0dr = []


	'''
	plt.figure(43)
	plt.plot(r,-2.5*log10(I0d)+m0,'o',color='black')
	I0 = Iedge_i(r,h1,10**(0.4*(m0-m0d1+5.*log10(pix2sec))))

	plt.plot(r,-2.5*log10(I0)+m0,'s',color='red',markersize=3)
	plt.show()
	exit()
	'''

	if mod==0:
		for k in range(len(r)):
			if fabs(I0d[k]-Iedge_i(r[k],h1,10**(0.4*(m0-m0d1+5.*log10(pix2sec)))))<3.*backgr:
				rr.append(r[k])
				I0dr.append(I0d[k])
	else:
		for k in range(len(r)):
			if fabs(z00[k]-z0)<std(z00) and fabs(zc[k]-yc)<std(zc) and fabs(I0d[k]-Iedge_i(r[k],h1,10**(0.4*(m0-m0d1+5.*log10(pix2sec)))))<10.*backgr:
				rr.append(r[k])
				I0dr.append(I0d[k])		


	#plt.figure(30)
	#plt.plot(np.array(rr)*pix2sec,-2.5*log10(I0dr)+m0+5.*log10(pix2sec),'o',color='black')
	#plt.show()
	#exit()
	dataToFit = RealData(rr, I0dr)
	mo = Model(disk_edge_r)
	fitting = ODR(dataToFit, mo, [10**(0.4*(m0-m0d1+5.*log10(pix2sec))),h1,0.])
	fitting.set_job()
	fit_result = fitting.run()
	m0d = -2.5*log10(fit_result.beta[0])+m0+5.*log10(pix2sec)
	h = fit_result.beta[1]
	xc = fit_result.beta[2] + xc



	z0_bul = []
	z_c_bul = []
	I0d_bul = []
	x_new_bul = []

	s1 = range(int(xc-rmin),int(xc+rmin))
	for x in s1:
		z,z1,I = INITIALS_edge_on.z_cut(scidata,yc,x,m0,0.,0.,ny/2,backgr,plot=0)
		I0d1,I0d_e,z01,z0_e1,z_c = INITIALS_edge_on.fitting(z1,I,max(I),z01,z1[list(I).index(max(I))],plot=0)
		z_c_bul.append(z_c)
		z0_bul.append(z01)
		I0d_bul.append(I0d1)
		x_new_bul.append(x)

	Z0_bul = np.median(z0_bul)
	yc_bul = np.median(z_c_bul)
	r_bul = np.array(x_new_bul) - xc


	plt.figure(33)
	# Profile SB(r)
	plt.plot(r*pix2sec,-2.5*log10(I0d)+m0+5.*log10(pix2sec),'o',color='black')
	#plt.plot(r_bul*pix2sec,-2.5*log10(I0d_bul)+m0+5.*log10(pix2sec),'o',color='black')
	plt.plot(r*pix2sec,-2.5*log10(Iedge_i(r*pix2sec,h*pix2sec,fit_result.beta[0]))+m0+5.*log10(pix2sec),'s',color='red',markersize=3)
	plt.axvline(rmin*pix2sec,linestyle='--',color='blue')
	plt.axvline(-rmin*pix2sec,linestyle='--',color='blue')
	plt.xlabel(r'r (arcsec)', fontsize=fsize)
	plt.ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	#xlim(-rmax*pix2sec,rmax*pix2sec)
	#ylim(min(majorMag)+11,min(majorMag)-1.)
	ymin, ymax = ylim()
	ylim( ymax, ymin )
	plt.savefig("edge_sph_prof.png", transparent = False)

	print h*0.396,m0d,xc
	#exit()

	plt.figure(34)
	# Height z0(r)
	plt.plot((np.array(x_new)-xc)*pix2sec,pix2sec*np.array(z00),'o',color='black')
	plt.plot((np.array(x_new_bul)-xc)*pix2sec,pix2sec*np.array(z0_bul),'o',color='black')
	plt.axvline(rmin*pix2sec,linestyle='--',color='blue')
	plt.axvline(-rmin*pix2sec,linestyle='--',color='blue')
	plt.xlabel(r'r (arcsec)', fontsize=fsize)
	plt.ylabel(r'z$_0$ (arcsec)', fontsize=fsize)
	plt.savefig("z0_r.png", transparent = False)
	'''
	plt.figure(12)
	# Center yc(r)
	#plt.plot(np.array(x_new)-xc,zc,'o',color='black')
	#plt.plot(np.array(x_new_bul)-xc,z_c_bul,'o',color='black')
	plt.axvline(rmin,linestyle='--',color='blue')
	plt.axvline(-rmin,linestyle='--',color='blue')
	plt.xlabel(r'r (pix)', fontsize=fsize)
	plt.ylabel(r'y$_c$ (pix)', fontsize=fsize)
	'''
	

	zz = []
	xx = []

	#print I0d


	for k in range(len(zc)):
		if I0d[k]>3.*backgr and fabs(z00[k]-z0)<std(z00) and fabs(zc[k]-yc)<3.*z0: #and fabs(zc[k]-yc)<1.*z0:
			xx.append(x_new[k]-xc)
			zz.append(zc[k]-yc)
	#plt.plot(xx,zz,'o',color='red')
	
	#plt.savefig("yc_r.png", transparent = False)
	#exit()

	#return xc,yc,m0d,h,z0,np.array(x_new)-xc,np.array(zc)-yc
	#print m0d,h,z0
	#exit()

	'''
	plt.figure(100)
	plt.plot(xx,zz,'o')
	plt.show()
	exit()
	'''
	return xc,yc,m0d,h,z0,np.array(xx),np.array(zz)

def Iser(r,meb,reb,n,m0,pix2sec):
	# Sersic function: for bulge and bar
	Ieb = 10**(0.4*(m0+5.*log10(pix2sec)-meb))
	nu = 1.9992*n - 0.3271
	Iser = Ieb * np.exp( -nu * ( (fabs(r)/reb)**(1.0/n) - 1.0 ) )
	return Iser

def Iedge(r,m0d,h,m0,pix2sec):
	# Radial distribution for the edge-on disc
	I0d = 10**(0.4*(m0+5.*log10(pix2sec)-m0d))
	r = fabs(r)
	I = I0d * r/h * sp.special.kn(1,r/h)

	return I

def rmax_bulge(meb,reb,n,backgr,m0,pix2sec):
	def func1(r):
		return Iser(r,meb,reb,n,m0,pix2sec) - 2.*backgr
	rmax = sp.optimize.bisect(func1,0,reb*100.)
	return rmax

def rmax_disk(m0d,h,backgr,m0,pix2sec):
	def func1(r):
		return Iedge(r,m0d,h,m0,pix2sec) - 2.*backgr
	rmax = sp.optimize.bisect(func1,1.*h,h*100.)
	return rmax


def main(ff,gal,model,number,numb,meb,reb,n,backgr,m0,pix2sec,FWHM,xc,yc,m0d,h,z0,smi,rmax,status): #,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,m0,pix2sec):
 if model=='edge':
 	ell_bul=99999.;b4_bul=99999.

 if model=='ser':
 	ell_disk=99999.;b4_disk=99999.

 if status!=1:
	print 'I am here 1'
	#sub_out = './pics/'+str(number)+'/sub_'+str(numb)+'.fits'
	#shutil.copy(sub_out,'subcomps.fits')

	if model=='exp+ser' or model=='ser' or model=='edge+ser':
		# BULGE
		#os.remove('disk.fits')
		#os.remove('bulge.fits')
		#os.chmod(r"disk.fits",0777)
		#os.chmod(r"galaxy_clean.fits",0777)
		#iraf.imcopy('subcomps.fits'+'[2]','disk.fits')
		#exit()

		#iraf.imarith(operand1="galaxy_clean.fits",op="-",operand2="disk.fits",result="bulgeRES.fits")
		f = open("resid.cl", "w") 
		sys.stdout = f
		print "#Soph"
		print "images"
		print "imutil"
		print "imarith galaxy_clean.fits - bulge.fits diskRES.fits"
		print "imarith galaxy_clean.fits - disk.fits bulgeRES.fits"
		print "logout"	
		sys.stdout = tmp_out
		f.close()
		os.chmod(r"resid.cl",0777)
		subprocess.call("cl < resid.cl -o", shell=True)




		#iraf.imarith('galaxy_clean.fits','-','disk.fits','bulgeRES.fits')
		#os.remove('disk.fits')

		#1. Median filtered image of the bulge:
		hdulist = pyfits.open('bulgeRES.fits')
		scidata = hdulist[0].data
		ny,nx = scidata.shape
		median_fil.median_filter('bulgeRES.fits','gal_med.fits')
		plots.image('gal_med.fits','gal_med',0,float(nx)/float(ny))

		#2. Ellipse task for it
		rmax = rmax_bulge(meb,reb,n,backgr,m0,pix2sec)
		step = 1.0
		minsma = 1.0
		maxsma = rmax
		#print bcolors.OKBLUE+ 'Ellipse fitting ...' + bcolors.ENDC
		#print >> ff, 'Ellipse fitting ...'
		try:	
			ell.main_ell('bulgeRES.fits',xc,yc,step,minsma,maxsma,m0,pix2sec,FWHM,rmax,'bulge_ell.txt')
			ell_bul,b4_bul = ell.read_ell('bulge_ell.txt',reb,step)
		except:
			step = 1.
			ell_bul,b4_bul = ell.read_ell('ellipse.txt',reb,step)
			
		print ell_bul,b4_bul


	if model=='edge+ser' or model=='edge':
		# DISK
		# A. Subtract bulge from galaxy image => disk	
		#iraf.imcopy('subcomps.fits'+'[3]','bulge.fits')
		#iraf.imarith(operand1='galaxy_clean.fits',op="-",operand2='bulge.fits',result='diskRES.fits')
		#iraf.imarith('galaxy_clean.fits','-','bulge.fits','diskRES.fits')
		#os.remove('bulge.fits')

		try:
			#1. Add the high peak at the disk center to run Ellipse task
		 	hdulist1 = pyfits.open('diskRES.fits', do_not_scale_image_data=True, mode='update')
			img_new = hdulist1[0].data
			img_new[int(yc),int(xc)] = 10**(0.4*(m0 - m0d + 2. + 5.*log10(pix2sec)))
			hdulist1.flush()

			#2. Ellipse task for it
			rmax = rmax_disk(m0d,h,backgr,m0,pix2sec)
			step = 1.0
			minsma = 1.0
			maxsma = rmax
			#print bcolors.OKBLUE+ 'Ellipse fitting ...' + bcolors.ENDC
			#print >> ff, 'Ellipse fitting ...'
			ell.main_ell('diskRES.fits',xc,yc,step,minsma,maxsma,m0,pix2sec,FWHM,rmax,'disk_ell.txt')
			ell_disk,b4_disk = read_ell('disk_ell.txt',rmax,step)
			print ell_disk,b4_disk
		except:
			ell_disk=99999.;b4_disk=99999.


		rmin = smi
		zmin1 = 0.
		zmin2 = 0.
		zmax = smi+smi/3.


		#3. Photometric cuts
		hdulist = pyfits.open('diskRES.fits')
		scidata = hdulist[0].data

		xc,yc,m0d,h,z0,xx,yy = cuts(0,scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,h,m0d,m0,pix2sec,backgr)

		yc_s,z0_s,n_s = cuts_soph(scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,h,m0d,m0,pix2sec,backgr)
		#r,I = slit(scidata,xc,yc,z0,rmax,m0,pix2sec)
		#breaks(r,I,rmin,rmax,m0,pix2sec)
		#4. Warps of disks
		omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left = warp.warp_pars(scidata,ff,xx,yy,xc,yc,pix2sec,m0,backgr,h,z0,rmax,number)



 else:
	print 'I am here 2'
	#1. Median filtered image of the galaxy:
	hdulist = pyfits.open('galaxy_clean.fits')
	scidata = hdulist[0].data
	ny,nx = scidata.shape
	#median_fil.median_filter('galaxy_clean.fits','gal_med.fits')
	#plots.image('gal_med.fits','gal_med',0,float(nx)/float(ny))


	if model=='exp+ser' or model=='ser' or model=='edge+ser':
		#2. Bulge isophotes
		step = 1.
		ell_bul,b4_bul = ell.read_ell('ellipse.txt',reb,step)
		print ell_bul,b4_bul


	if model=='edge+ser' or model=='edge':
		ell_disk,b4_disk = 99999.,99999.
		# DISK
		step = 1.
		#ell_disk,b4_disk = ell.read_ell('ellipse.txt',rmax,step)
		#if ell_disk==nan or b4_disk==nan:
			
		print ell_disk,b4_disk


		rmin = smi
		zmin1 = 0.
		zmin2 = 0.
		zmax = smi+smi/3.


		#3. Photometric cuts
		hdulist = pyfits.open('galaxy_clean.fits')
		scidata = hdulist[0].data
		#print xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,h,m0d,m0,pix2sec,backgr
		#exit()
		xc,yc,m0d,h,z0,xx,yy = cuts(1,scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,h,m0d,m0,pix2sec,backgr)
		yc_s,z0_s,n_s = cuts_soph(scidata,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,h,m0d,m0,pix2sec,backgr)
		#r,I = slit(scidata,xc,yc,z0,rmax,m0,pix2sec)
		#breaks(r,I,rmin,rmax,m0,pix2sec)
		#4. Warps of disks

		omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left = warp.warp_pars(scidata,ff,xx,yy,xc,yc,pix2sec,m0,backgr,h,z0,rmax,number)


 if not os.path.exists("%s" % (file_out_edge)):
	ffff = open(file_out_edge,'w')
	print >>ffff, '#\tm0d\th\tz0\telld\tb4d\tellb\tb4bul\tn_disk\tomega\tomega_l\tomega_r\talfa_s\tr_last_r\tr_last_l\ty_last_r\ty_last_l\tpsi\tA_r\t,C_r\tA_l\tC_l'
	print >>ffff, '\tmas\tarcsec\tarcsec\t\t\t\t\t\t\t\t\t\t\t\t\t\tarcsec\t\tarcsec'
	ffff.close()


 ffff = open(file_out_edge,'a')
 print >>ffff, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (number,m0d-2.5*log10(z0/h),h*pix2sec,z0*pix2sec,ell_disk,b4_disk,ell_bul,b4_bul,n_s,omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right*pix2sec,C_right,A_left*pix2sec,C_left)
 ffff.close()


 print 'Sech2 distribution:'
 print '\tm0d=%.2f (mag arcsec-2), h=%2.f (arcsec), z0=%.2f (arcsec)' % (m0d-2.5*log10(z0/h),h*pix2sec,z0*pix2sec)
 print 'N-law distribution:'
 print ' n=%.1f' % (n_s)

 print '\n Disk and Bulge Isophotes:'
 print '\te_disk=%.2f, b4_disk=%.2f' % (ell_disk,b4_disk)
 print '\te_bulge=%.2f, b4_bulge=%.2f' % (ell_bul,b4_bul)
 return 0


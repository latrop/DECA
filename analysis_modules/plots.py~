#!/usr/bin/python
# -*- coding:  cp1251 -*-
import fpdf
#import aplpy
from fpdf import FPDF
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
#http://code.google.com/p/pyfpdf/
#http://aplpy.github.com/

pix2sec = 1.
#model = 'exp+ser'
model = 'edge+ser'
m0_d = 18.23
m0_d_e = 0.12
h_d = 1.23
h_d_e = 0.12
#z0 = 0.23
z0=0.54
z0_e = 0.02
ell_d = 0.74
ell_d_e = 0.02
PA_d = 43.2
PA_d_e = 0.5
c_d = 0.
me_bul = 19.32
me_bul_e = 0.32
re_bul = 1.23
re_bul_e = 0.21
n_bul = 1.23
n_bul_e = 0.13
ell_bul = 0.38
ell_bul_e = 0.02
PA_bul = 23.2
PA_bul_e = 0.2
c_bul = 0. 
m_gal = 11.42
M_gal = -20.32
BT = 0.35
chi2 = 1.261
scale = 'arcsec'
if scale=='kpc':
	if (z0<=0. or z0==0) and h_d>0 and re_bul>0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f kpc, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,ell_d,ell_d_e,PA_d,PA_d_e,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,PA_bul,PA_bul_e,m_gal,M_gal,BT,chi2)


	if (z0<=0. or z0==0) and h_d>0 and re_bul<=0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f kpc, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n  Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,ell_d,ell_d_e,PA_d,PA_d_e,m_gal,M_gal,BT,chi2)

	if z0>0. and h_d>0 and re_bul>0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f kpc, z0=%.2f +/- %.2f kpc\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,z0,z0_e,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,m_gal,M_gal,BT,chi2)

	if z0>0. and h_d>0 and re_bul<=0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f kpc, z0=%.2f +/- %.2f kpc\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,z0,z0_e,m_gal,M_gal,BT,chi2)

	if h_d<=0 and re_bul>0:
		text = 'Model: %s\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,PA_bul,PA_bul_e,m_gal,M_gal,BT,chi2)



if scale=='arcsec':
	if (z0<=0. or z0==0) and h_d>0 and re_bul>0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,ell_d,ell_d_e,PA_d,PA_d_e,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,PA_bul,PA_bul_e,m_gal,M_gal,BT,chi2)


	if (z0<=0. or z0==0) and h_d>0 and re_bul<=0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n  Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,ell_d,ell_d_e,PA_d,PA_d_e,m_gal,M_gal,BT,chi2)

	if z0>0. and h_d>0 and re_bul>0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, z0=%.2f +/- %.2f arcsec\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, n=%.2f +/- %.2f, e=%.2f +/- %.2f\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,z0,z0_e,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,m_gal,M_gal,BT,chi2)

	if z0>0. and h_d>0 and re_bul<=0:
		text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, z0=%.2f +/- %.2f arcsec\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d,m0_d_e,h_d,h_d_e,z0,z0_e,m_gal,M_gal,BT,chi2)

	if h_d<=0 and re_bul>0:
		text = 'Model: %s\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,me_bul,me_bul_e,re_bul,re_bul_e,n_bul,n_bul_e,ell_bul,ell_bul_e,PA_bul,PA_bul_e,m_gal,M_gal,BT,chi2)









def read_fits(gal): 
# The function to read the fits file.
	hdulist = pyfits.open(gal)
	prihdr = hdulist[2].header
	scidata = hdulist[2].data
	return scidata

def image(gal,name,number,ratio):
	dx = 5.
	dy = 5./ratio
	try:
		fig = aplpy.FITSFigure(gal, hdu=number, figsize=(dx, dy))
		fig.show_grayscale(stretch='log',vmid=-10,vmin=-9, invert=True)
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
'''
gal = 'galfit_out.fits'
name = 'NGC4143'

scidata=read_fits(gal)
#print scidata.size # number of pixels at all
#print scidata.shape # (dimx,dimy)
[ny,nx] = scidata.shape
#print nx,ny
l=200.
ratio=float(nx)/float(ny)
#ratio = 3
dimy = l/(3.*ratio)
dimx=dimy*ratio
'''
'''
image(gal,'galaxy',1)
image(gal,'model',2)
image(gal,'residual',3)
'''
'''
pdf=FPDF()
pdf.add_page()
pdf.set_font('Times','B',18)
pdf.set_text_color(220,50,50)
pdf.cell(170,10,name,0,1,'C')
pdf.set_text_color(0,0,0)
pdf.set_font('Times','B',14)
pdf.write(5,text)


pdf.image('galaxy.png',x=5,y=50.,h=dimy,w=dimx)
pdf.image('model.png',x=75,y=50.,h=dimy,w=dimx)
pdf.image('residual.png',x=145,y=50.,h=dimy,w=dimx)
'''
'''
pdf.image('ell.png',x=5,y=112.,w=80.,h=190.)
#pdf.image('major_prof.png',x=70,y=110.,w=70.)
pdf.image('azim_prof.png',x=100,y=115.,w=90.)
pdf.image('radius_prof.png',x=100,y=205,w=90.)
#pdf.cell(100,100,text,0,1,'B')
'''
'''
pdf.image('ell.png',x=5,y=70.,w=65.,h=160.)
pdf.image('edge_prof.png',x=75,y=85,w=130.)
pdf.output('plots.pdf','F')
'''

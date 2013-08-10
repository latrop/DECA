#!/usr/bin/python
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
from matplotlib.mlab import griddata
import csv


fsize=25

hdulist = pyfits.open('galaxy_clean.fits')
scidata = hdulist[0].data
#prihdr = hdulist[0].header
ny,nx = np.shape(scidata)

x = []
y = []
I = []
I_crit = 10.
xc = 100.
yc = 50.
for k in range(ny):
	for i in range(nx):
		if scidata[k,i]>I_crit:
			x.append(i + 0.5-xc)
			y.append(k + 0.5-yc)
			I.append(scidata[k,i])
x = np.array(x)
y = np.array(y)
I = np.array(I)

m0 = 26.
pix2sec = 0.396
m = -2.5*log10(I) + m0 + 5.*log10(pix2sec)
m = np.around(m,decimals=1)
#print m[1]
#exit()

sample_pts = 500
con_levels = np.arange(16,20,0.5)
levels = np.arange(16,20,0.1)
'''
x = data['x']
xmin = x.min()
xmax = x.max()

y = data['y']
ymin = y.min()
ymax = y.max()

z = data['z']
'''
m1 = []
x1 = []
y1 = []

for k in range(len(m)):
	if m[k]<20 and m[k]>19.5:
		m1.append(m[k])
		x1.append(x[k])
		y1.append(y[k])

xi = np.linspace(min(x),max(x),sample_pts)
yi = np.linspace(min(y),max(y),sample_pts)

Ii = griddata(x,y,m,xi,yi)


f = plt.figure()
ax1 = f.add_subplot(111)
ax1.contour(xi,yi,Ii,con_levels,linewidths=1,colors='black')
#plt.contour(Ii,levels,linewidths=1)
plt.scatter(x1,y1,c=m1,s=20)
ax1.set_xlim(min(x),max(x))
ax1.set_ylim(min(y),max(y))
ax1.set_ylabel(r'z (arcsec)', fontsize=fsize)
ax2 = ax1.twinx()
ax2.set_ylabel(r'z (z$_0$)', fontsize=fsize)
ax2.set_ylim(-3,3)
ax2 = ax1.twiny()
ax2.set_xlabel(r'r (h)', fontsize=fsize)
ax2.set_xlim(-4,4)
plt.show()

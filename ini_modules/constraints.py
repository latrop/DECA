#!/usr/bin/python
# Script to create constraints file for GALFIT
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
#from pyraf import iraf

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

tmp_out = sys.stdout


def constr_spirals(model,trunc,fwhm,pix2sec,h0,r_in,r_out,ang_rot,incl,R10,n_iter):
	if (model=='exp+ser') and trunc==0:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    1	  	   9	    0.5 to 0.75	# Soft constraint: Constrains' 
		print '    1	  	   4	  %.3f to %.3f  # Soft constraint: Constrains' % (0.5*h0, 1.5*h0) 
		print '    1	  	   10	     0 to  180	# Soft constraint: Constrains'
		print '    1	  	   R1	  %.3f to  %.3f	# Soft constraint: Constrains' % (0.25*h0, r_in) 
		print '    1	  	   R2    %.3f to  %.3f	# Soft constraint: Constrains' % (r_in+(r_out-r_in)/3, 3.5*h0) 
		print '    1	  	   R3	    %.3f to %.3f	# Soft constraint: Constrains' % (ang_rot-50.,ang_rot+50.) 
		print '    1	  	   R9	    %.3f to %.3f	# Soft constraint: Constrains' % (incl-10.,incl+10.) 
		if n_iter==3:
			print '    1	  	   R10	    %.3f to %.3f	# Soft constraint: Constrains' % (R10-30.,R10+30.) 


		print '    2              n        0.5 to 8    # Sersic n'
		print '    2              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    2              c0       0 to 1.5    # c0'
		print '    2	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    2	  	   y	     -2  2	# Soft constraint: Constrains'

	else:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    1	  	   9	    0.5 to 0.75	# Soft constraint: Constrains' 
		print '    1	  	   4	  %.3f to %.3f  # Soft constraint: Constrains' % (0.5*h0, 1.5*h0) 
		print '    1	  	   10	     0 to  180	# Soft constraint: Constrains'
		print '    1	  	   R1	  %.3f to  %.3f	# Soft constraint: Constrains' % (0.25*h0, r_in) 
		print '    1	  	   R2    %.3f to  %.3f	# Soft constraint: Constrains' % (r_in+(r_out-r_in)/3, 3.5*h0) 
		print '    1	  	   R3	    %.3f to %.3f	# Soft constraint: Constrains' % (ang_rot-50.,ang_rot+50.) 
		print '    1	  	   R9	    %.3f to %.3f	# Soft constraint: Constrains' % (incl-10.,incl+10.) 
		if n_iter==3:
			print '    1	  	   R10	    %.3f to %.3f	# Soft constraint: Constrains' % (R10-30.,R10+30.) 


def constr_spirals_fix(model,trunc,fwhm,pix2sec,h0,r_in,r_out,ang_rot,incl,R10,n_iter):
	if (model=='exp+ser') and trunc==0:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    1	  	   9	    0.5 to 0.75	# Soft constraint: Constrains' 
		print '    1	  	   4	  %.3f to %.3f  # Soft constraint: Constrains' % (0.5*h0, 1.5*h0) 
		print '    1	  	   10	     0 to  180	# Soft constraint: Constrains'
		print '    1	  	   R1	  %.3f to  %.3f	# Soft constraint: Constrains' % (0.25*h0, r_in) 
		print '    1	  	   R2    %.3f to  %.3f	# Soft constraint: Constrains' % (r_in+(r_out-r_in)/3, 3.5*h0) 
		print '    1	  	   R3	    %.3f to %.3f	# Soft constraint: Constrains' % (ang_rot-50.,ang_rot+50.) 
		print '    1	  	   R9	    %.3f to %.3f	# Soft constraint: Constrains' % (incl-10.,incl+10.) 
		if n_iter==3:
			print '    1	  	   R10	    %.3f to %.3f	# Soft constraint: Constrains' % (R10-30.,R10+30.) 

		print '    2              n        0.5 to 8    # Sersic n'
		print '    2              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    2              c0       0 to 1.5    # c0'
		print '    2	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    2	  	   y	     -2  2	# Soft constraint: Constrains'

		print '    1/2	  	   re	     1.3   4	# Soft constraint: Constrains'
		print '    1-2	  	   mag	     -1.5   1.5	# Soft constraint: Constrains'
		print '    1	  	   mag	     -0.5   0.5	# Soft constraint: Constrains'
		print '    2	  	   mag	     -0.5   0.5	# Soft constraint: Constrains'

	else:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    1	  	   9	    0.5 to 0.75	# Soft constraint: Constrains' 
		print '    1	  	   4	  %.3f to %.3f  # Soft constraint: Constrains' % (0.5*h0, 1.5*h0) 
		print '    1	  	   10	     0 to  180	# Soft constraint: Constrains'
		print '    1	  	   R1	  %.3f to  %.3f	# Soft constraint: Constrains' % (0.25*h0, r_in) 
		print '    1	  	   R2    %.3f to  %.3f	# Soft constraint: Constrains' % (r_in+(r_out-r_in)/3, 3.5*h0) 
		print '    1	  	   R3	    %.3f to %.3f	# Soft constraint: Constrains' % (ang_rot-50.,ang_rot+50.) 
		print '    1	  	   R9	    %.3f to %.3f	# Soft constraint: Constrains' % (incl-10.,incl+10.) 
		if n_iter==3:
			print '    1	  	   R10	    %.3f to %.3f	# Soft constraint: Constrains' % (R10-30.,R10+30.) 



def constr_bul(model,trunc,fwhm,pix2sec,h0):
	if (model=='edge+ser' or model=='exp+ser') and trunc==0:
		print '    1	  	   x	     -3  3	# Soft constraint: Constrains' 
		print '    1	  	   y	     -3  3	# Soft constraint: Constrains'  
		print '    2              n        0.5 to 6    # Sersic n'
		print '    2              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    2              c0       0 to 2    # c0'
		print '    2	  	   x	     -3  3	# Soft constraint: Constrains' 
		print '    2	  	   y	     -3  3	# Soft constraint: Constrains'
		print '    2	  	   q	     0.4 to  1	# Soft constraint: Constrains'
		#print '    1/2             re         1  100    # Soft constraint: Constrains'
		if model=='exp+ser':
			print '    1	  	   q	     0.4 to  1	# Soft constraint: Constrains'
			print '    1/2	  	   re	     1.3   4	# Soft constraint: Constrains'
			print '    1-2	  	   mag	     -1.5   1.5	# Soft constraint: Constrains'

	if (model=='edge+ser' or model=='exp+ser') and trunc==1:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'   
		print '    3              n        0.5 to 8    # Sersic n'
		print '    3              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    3              c0       0 to 2    # c0'
		print '    3	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    3	  	   y	     -2  2	# Soft constraint: Constrains'


	if (model=='edge+ser' or model=='exp+ser') and trunc==2:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    4              n        0.5 to 8    # Sersic n'
		print '    4              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    4	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    4	  	   y	     -2  2	# Soft constraint: Constrains' 


	if (model=='edge+ser' or model=='exp+ser') and trunc==3:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    5              n        0.5 to 8    # Sersic n'
		print '    5              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    5	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    5	  	   y	     -2  2	# Soft constraint: Constrains' 


	if model=='ser':
		print '    1              n        0.5 to 8    # Sersic n'
		print '    1              re       %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec)
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'


def constr_bul_fix(model,trunc,fwhm,pix2sec,h0,rmax):
	if (model=='edge+ser' or model=='exp+ser') and trunc==0:
		print '    1	  	   x	     -3  3	# Soft constraint: Constrains' 
		print '    1	  	   y	     -3  3	# Soft constraint: Constrains'  
		print '    2              n        0.5 to 6    # Sersic n'
		print '    2              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    2              c0       0 to 2    # c0'
		print '    2	  	   x	     -3  3	# Soft constraint: Constrains' 
		print '    2	  	   y	     -3  3	# Soft constraint: Constrains'
		print '    2	  	   q	     0.4 to  1	# Soft constraint: Constrains'
		#print '    1/2             re         1  100    # Soft constraint: Constrains'
		if model=='exp+ser':
			print '    1	  	   q	     0.4 to  1	# Soft constraint: Constrains'
			print '    2	  	   re	     %.3f to %.3f	# Soft constraint: Constrains' % (h0/6.,h0)
			print '    1	  	   re	     0 to %.3f	# Soft constraint: Constrains' % (rmax/2.5)
			#print '    1-2	  	   mag	     -1.5   1.5	# Soft constraint: Constrains'
			print '    1	  	   mag	     -0.5   0.5	# Soft constraint: Constrains'
			print '    2	  	   mag	     -0.5   0.5	# Soft constraint: Constrains'

	if (model=='edge+ser' or model=='exp+ser') and trunc==1:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'   
		print '    3              n        0.5 to 8    # Sersic n'
		print '    3              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    3              c0       0 to 2    # c0'
		print '    3	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    3	  	   y	     -2  2	# Soft constraint: Constrains'


	if (model=='edge+ser' or model=='exp+ser') and trunc==2:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    4              n        0.5 to 8    # Sersic n'
		print '    4              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    4	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    4	  	   y	     -2  2	# Soft constraint: Constrains' 


	if (model=='edge+ser' or model=='exp+ser') and trunc==3:
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  
		print '    5              n        0.5 to 8    # Sersic n'
		print '    5              re       %.3f to %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec,5.*h0)
		print '    5	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    5	  	   y	     -2  2	# Soft constraint: Constrains' 


	if model=='ser':
		print '    1              n        0.5 to 8    # Sersic n'
		print '    1              re       %.3f    # reb>0.8*fwhm(PSF)' % (0.8*fwhm/pix2sec)
		print '    1	  	   x	     -2  2	# Soft constraint: Constrains' 
		print '    1	  	   y	     -2  2	# Soft constraint: Constrains'  


def constr1(file_constraints,model,trunc,fwhm,pix2sec,h0,r_in,r_out,ang_rot,incl,R10=0.,n_iter=0,fix=0):
	f = open(file_constraints, 'w') 
	sys.stdout = f

	trunc = 0
	if fix==0:
		constr_spirals(model,trunc,fwhm,pix2sec,h0,r_in,r_out,ang_rot,incl,R10,n_iter)
	else:
		constr_spirals_fix(model,trunc,fwhm,pix2sec,h0,r_in,r_out,ang_rot,incl,R10,n_iter)

	sys.stdout = tmp_out
	f.close()


def constr(file_constraints,model,trunc,fwhm,pix2sec,h0,rmax, fix=0):
	f = open(file_constraints, 'w') 
	sys.stdout = f

	if fix==0:
		constr_bul(model,trunc,fwhm,pix2sec,h0)
	else:
		constr_bul_fix(model,trunc,fwhm,pix2sec,h0,rmax)		

	sys.stdout = tmp_out
	f.close()

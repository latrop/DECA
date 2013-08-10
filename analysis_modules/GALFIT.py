#!/usr/bin/python
# -*- coding:  cp1251 -*-
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
#import subprocess
import signal

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')
import setup



tmp_out = sys.stdout

file_weight = 'weight.fits'
file_constraints = 'constr.txt'
file_galfit_input = 'input.txt'
file_field = 'field.fits'
file_psf_field = 'psField.fit'	# for SDSS
file_field1 = 'field1.fits'
file_field2 = 'field2.fits'
file_gal = 'galaxy.fits'
file_sex_field = 'field.cat'
file_sex_galaxy = 'galaxy.cat'
file_gal_clean = 'galaxy_clean.fits'
file_gal_clean_new = 'galaxy_clean1.fits'
file_segm = 'segm.fits'
file_badpix = 'badpix.txt'
file_inter = 'interact.cat'
file_psf = 'psf.fits'
file_psf_model = 'psf_model.fits'
file_nonpar = 'nonpar.txt'

file_galfit_outimage = 'galfit_out.fits'

import subprocess as sub
import threading
'''
class RunCmd(threading.Thread):
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
        self.p = sub.Popen(self.cmd)
        self.p.wait()

    def Run(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            self.p.terminate()
            self.join()
'''

class RunCmd(threading.Thread):
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
        self.p = sub.Popen(self.cmd)
        self.p.wait()
	return  

    def Run(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            self.p.terminate()
            self.join()
	return self.p.wait()


# ********************************* Simulation with GALFIT *********************************
def sym_galfit(nx,ny,xc,yc,I0d=0.,h=0.,z0=0.,rmind=0.,PAd=0.,elld=0.,cd=-1.,Mb=0.,reb=0.,n=0.,ellb=0.,PAb=0.,cb=0.,rmaxb=0.,sky=0.0,gain=9.9,ncombine=1):
	rtr=nx
	if reb>0.: rmaxb=nx

	m0 = 20.4885

	f = open(r"modelIN.txt", "w") 
	sys.stdout = f

	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) none                # Input data image (FITS file)"
	print "B) galaxy.fits         # Output data image block"
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) none                # Bad pixel mask (FITS image or ASCII coord list)"
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx, ny)
	print "I) 0    0              # Size of the convolution box (x y)"
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) 1.000  1.000        # Plate scale (dx dy)    [arcsec per pixel]"
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"


	if I0d>0.: 
	 m0d = -2.5*log10(I0d)
	 print "# Edge-on disk function\n"
 	 print "0) edgedisk            # Object type"
 	 print "1) %.3f  %.3f  0 0     # position x, y        [pixel]" % (xc,yc)
 	 print "3) %.3f       0       # central surface brightness  [mag/arcsec^2]" % (m0d)
 	 print "4) %.3f       0       # disk scale-height    [Pixels]" % (z0)
	 print "5) %.3f       0       # disk scale-length    [Pixels]" % (h)
	 print "10) 90.0         0       # position angle (PA)  [Degrees: Up=0, Left=90]"
 	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n"
	 print "C0) %.3f       0       # disk scale-height    [Pixels]" % (cd)


	if reb>0.:
	 print "# Sersic function\n"
	 print " 0) sersic             # Object type"
	 print " 1) %.3f  %.3f  0 0    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       0       # total magnitude"  % (Mb)   
	 print " 4) %.3f       0       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       0       # Sersic exponent (deVauc=4, expdisk=1)" % (n) 
	 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) 90.0       0       # position angle (PA)  [Degrees: Up=0, Left=90]"
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

 	print "\n================================================================================"                   

	sys.stdout = tmp_out
	f.close()
	
	path = str(i)
	os.mkdir(path)
	os.chmod(r"modelIN.txt",0777)
	subprocess.call("galfit modelIN.txt", shell=True)
	os.chmod(path,0777)




# ********************************* Decomposition with GALFIT *********************************
def dec_galfit(find_psf,m0,pix2sec,nx,ny,xc,yc,m0d,Md,h,z0,elld,cd,rtr,meb,Mb,reb,n,ellb,cb,sky_level=0.,PAd=90.,PAb=90.):
	if os.path.exists('warps_coords.txt')==True and z0!=99999. and setup.warps_cut==0:
		X_l1,Y_l1,X_r1,Y_r1,X_l2,Y_l2,X_r2,Y_r2 = loadtxt('warps_coords.txt', usecols=[0,1,2,3,4,5,6,7], unpack=True, skiprows = 1)
		try:
			x_l1=X_l1[number-1]; y_l1=Y_l1[number-1]; x_r1=X_r1[number-1]; y_r1=Y_r1[number-1]
			x_l2=X_l2[number-1]; y_l2=Y_l2[number-1]; x_r2=X_r2[number-1]; y_r2=Y_r2[number-1]
		except:
			x_l1=X_l1; y_l1=Y_l1; x_r1=X_r1; y_r1=Y_r1
			x_l2=X_l2; y_l2=Y_l2; x_r2=X_r2; y_r2=Y_r2

		xmin = int(x_r1); xmax = int(x_l2)
	else:
		xmin = 1
		xmax = nx

	print 'mod!!!',m0d
	global file_psf
	if os.path.exists(file_weight)==False:
		file_sigma = 'none'
	else:
		file_sigma = file_weight
	if find_psf==10:
		file_psf = 'none'


	if reb>0. and reb!=99999.: rmaxb=nx
	#if rtr==99999.:	rtr=nx
	f = open(file_galfit_input, "w") 
	sys.stdout = f


	
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                  # Input data image (FITS file)" % (file_gal_clean)
	print "B) %s                  # Output data image block" % (file_galfit_outimage)
	print "C) %s                # Sigma image name (made from data if blank or none)" % (file_sigma)
	print "D) %s                # Input PSF image and (optional) diffusion kernel" % (file_psf)
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (file_badpix)
	print "G) %s                # File with parameter constraints (ASCII file)" % (file_constraints)
	print "H) %i    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (xmin, xmax, ny)
	print "I) %i    %i              # Size of the convolution box (x y)" % (setup.box_psf,setup.box_psf)
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2sec,pix2sec)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"


	if h!=99999. and z0!=99999. and m0d!=99999. and Md==99999.: 
	 #m0d = -2.5*log10(I0d)
	 print "# Edge-on disk function\n"
 	 print "0) edgedisk            # Object type"
 	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
 	 print "3) %.3f       1       # central surface brightness  [mag/arcsec^2]" % (m0d)
 	 print "4) %.3f       1       # disk scale-height    [Pixels]" % (z0)
	 print "5) %.3f       1       # disk scale-length    [Pixels]" % (h)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
 	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n"
	 print "C0) %.3f       0       # diskiness/boxiness" % (-1.0)
	 if rtr!=99999.:
		 print "To) 2                  # Outer truncation by component 2"
		 print "T0) radial	 #  truncation" 
		 #print "T3) %.3f       0      #  Trancation radius-break beginning" % (rmax_d)
		 print "T4) %.3f       1      #  Trancation radius-break end" % (rtr)
		 print "T5) %.3f      1          #  Softening length (1 normal flux) [pixels]" % (5.)


	if h!=99999. and z0!=99999. and m0d==99999. and Md!=99999.: 
	 #m0d = -2.5*log10(I0d)
	 print "# Edge-on disk function\n"
 	 print "0) edgedisk1            # Object type"
 	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
 	 print "3) %.3f       1       # central surface brightness  [mag/arcsec^2]" % (Md)
 	 print "4) %.3f       1       # disk scale-height    [Pixels]" % (z0)
	 print "5) %.3f       1       # disk scale-length    [Pixels]" % (h)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
 	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n"
	 print "C0) %.3f       1       # diskiness/boxiness" % (cd)
	 '''
	 print "To) 2                  # Outer truncation by component 2"
	 print "T0) radial	 #  truncation" 
	 #print "T3) %.3f       0      #  Trancation radius-break beginning" % (rmax_d)
	 print "T4) %.3f       1      #  Trancation radius-break end" % (rtr)
	 print "T5) %.3f      1          #  Softening length (1 normal flux) [pixels]" % (5.)
	 '''

	if h!=99999. and z0==99999. and m0d!=99999. and Md==99999.:
	 print "# Exponential disk function\n"
	 print "0) expdisk1            # Object type"
	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (m0d)   
	 print " 4) %.3f       1       #     R_s              [Pixels]" % (h)
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-elld)  
	 print " 10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

	if h!=99999. and z0==99999. and m0d==99999. and Md!=99999.:
	 print "# Exponential disk function\n"
	 print "0) expdisk            # Object type"
	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (Md)   
	 print " 4) %.3f       1       #     R_s              [Pixels]" % (h)
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-elld)  
	 print " 10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"


	if reb!=99999. and meb!=99999. and Mb==99999.:
	 print "# Sersic function\n"
	 print " 0) sersic2             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (meb)   
	 print " 4) %.3f       1       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       1       # Sersic exponent (deVauc=4, expdisk=1)" % (n)
	 if setup.fix_ell==1: 
		 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-ellb)
	 else:
		 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 if cb!=99999. and setup.fix_c==1:
	 	print "C0) %.3f       1       # diskiness/boxiness" % (cb)

	if reb!=99999. and meb==99999. and Mb!=99999.:
	 print "# Sersic function\n"
	 print " 0) sersic             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (Mb)   
	 print " 4) %.3f       1       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       1       # Sersic exponent (deVauc=4, expdisk=1)" % (n) 
	 if setup.fix_ell==1: 
		 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-ellb)
	 else:
		 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 #if cb!=99999.:
	 #	print "C0) %.3f       1       # diskiness/boxiness" % (cb)

	if sky_level!=99999.:
	 print "# sky\n"
	 print " 0) sky                    #  object type"
	 print " 1) %.3f      1          #  sky background at center of fitting region [ADUs]" % (sky_level)
	 print " 2) 0.0000      0          #  dsky/dx (sky gradient in x)"
	 print " 3) 0.0000      0          #  dsky/dy (sky gradient in y)"
	 print " Z) 0                 #  Skip this model in output image?  (yes=1, no=0)"

 	print "\n================================================================================"                   

	sys.stdout = tmp_out
	f.close()
	

	#p = subprocess.Popen('galfit input.txt',  shell=True, close_fds=False)
	#p.kill
	#status = p.wait() 

	status = RunCmd(["galfit", "input.txt"], setup.TIME).Run()
	'''
	if os.path.exists('galfit.01')==False:
		status = 1
	else:
		status = 0
	'''
	#print status
	#exit()

	if int(status)==-6:	status=0

	if status!=0:
		return status
	else:
		sub.call("galfit -o3 galfit.01", shell=True, close_fds=False)
		return status





################################## ==========SPIRALS==========#################################
# ********************************* Decomposition with GALFIT *********************************
def dec_galfit_spirals(find_psf,m0,pix2sec,nx,ny,xc,yc,m0d0,h0,PAd,optimalInnerRadius,diskR25,angularRotation,galaxyIncl,pa,rtr,meb,reb,n,ellb,cb,sky_level=0.,fix=1,PAb=90., q_d=0.6): # q_d was added as last parameter (Sergey)
	#print 'mod!!!',m0d
	global file_psf
	if os.path.exists(file_weight)==False:
		file_sigma = 'none'
	else:
		file_sigma = file_weight
	if find_psf==10:
		file_psf = 'none'


	if reb>0. and reb!=99999.: rmaxb=nx
	#if rtr==99999.:	rtr=nx
	f = open(file_galfit_input, "w") 
	sys.stdout = f


	
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                  # Input data image (FITS file)" % (file_gal_clean)
	print "B) %s                  # Output data image block" % (file_galfit_outimage)
	print "C) %s                # Sigma image name (made from data if blank or none)" % (file_sigma)
	print "D) %s                # Input PSF image and (optional) diffusion kernel" % (file_psf)
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (file_badpix)
	print "G) %s                # File with parameter constraints (ASCII file)" % (file_constraints)
	print "H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx, ny)
	print "I) %i    %i              # Size of the convolution box (x y)" % (50,50)
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2sec,pix2sec)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"


	if h0!=99999.  and m0d0!=99999.:
	 print "# Exponential disk function\n"
	 print "0) expdisk1                     # Object type"
	 print "1)   %7s  %7s   1 1    # position x, y        [pixel]" % (xc, yc)
	 print "3)   %7s             1      # total magnitude" % (m0d0)
	 print "4)   %7s             1      # disk scale-length    [Pixels]" % (h0)
	 print "9)     %5.3f             1      # Axis ratio (b/a) -- width of spirals" % (q_d)
	 print "10)  %7.2f             1      # Spirals position angle" % (pa)
	 print "R0)      log                    # PA rotation func. (tanh, sqrt, log, linear, none)"
	 print "R1)  %7i             1      # inner (bar radius)" % optimalInnerRadius
	 print "R2)  %7i             1      # outer radius" % diskR25
	 print "R3)  %7.1f             1      # Cumul. coord.rotation [degrees]" % (angularRotation)
	 print "R4)        1             0      # winding scale"
	 print "R9)  %7.1f             1      # inclination of entire galaxy" % (galaxyIncl)
	 print "R10) %7.1f             1      # position angle (of galaxy)" % (PAd-90)
	 #    print "C0)        0             1      # Fourier mode."
	 print "Z)         0                    #  Skip this model in output image?  (yes=1, no=0)"




	if reb!=99999. and meb!=99999. and fix==1:
	 print "# Sersic function\n"
	 print " 0) sersic2             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (meb)   
	 print " 4) %.3f       1       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       1       # Sersic exponent (deVauc=4, expdisk=1)" % (n)
	 if setup.fix_ell==1: 
		 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-ellb)
	 else:
		 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 if cb!=99999.:
	 	print "C0) %.3f       1       # diskiness/boxiness" % (cb)

	if reb!=99999. and meb!=99999. and fix==0:
	 print "# Sersic function\n"
	 print " 0) sersic2             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       0       # total magnitude"  % (meb)   
	 print " 4) %.3f       0       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       0       # Sersic exponent (deVauc=4, expdisk=1)" % (n)
	 if setup.fix_ell==1: 
		 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)
	 else:
		 print " 9) %.3f       0       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 if cb!=99999. and cb!=0.:
	 	print "C0) %.3f       0       # diskiness/boxiness" % (cb)



	if sky_level!=99999.:
	 print "# sky\n"
	 print " 0) sky                    #  object type"
	 print " 1) %.3f      1          #  sky background at center of fitting region [ADUs]" % (sky_level)
	 print " 2) 0.0000      0          #  dsky/dx (sky gradient in x)"
	 print " 3) 0.0000      0          #  dsky/dy (sky gradient in y)"
	 print " Z) 0                 #  Skip this model in output image?  (yes=1, no=0)"

 	print "\n================================================================================"                   

	sys.stdout = tmp_out
	f.close()
	

	#p = subprocess.Popen('galfit input.txt',  shell=True, close_fds=False)
	#p.kill
	#status = p.wait() 

	status = RunCmd(["galfit", "input.txt"], setup.TIME).Run()


	#print status
	#exit()

	if int(status)==-6:	status=0

	if status!=0:
		return status
	else:
		sub.call("galfit -o3 galfit.01", shell=True, close_fds=False)
		return status



#########################################################################################################
#################################################### END ################################################










# ********************************* Decomposition with GALFIT *********************************
def dec_galfit_object(find_psf,file_in,file_out,m0,pix2sec,nx,ny,xc,yc,m0d,Md,h,z0,elld,cd,rtr,meb,Mb,reb,n,ellb,cb,sky_level=0.,PAd=90.,PAb=90.):
	global file_psf
	if os.path.exists(file_weight)==False:
		file_sigma = 'none'
	else:
		file_sigma = file_weight
	if find_psf==10:
		file_psf = 'none'
		

	if reb>0. and reb!=99999.: rmaxb=nx
	if rtr==99999.:	rtr=nx
	f = open(file_galfit_input, "w") 
	sys.stdout = f

	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                  # Input data image (FITS file)" % (file_in)
	print "B) %s                  # Output data image block" % (file_out)
	print "C) %s                # Sigma image name (made from data if blank or none)" % (file_sigma)
	print "D) %s                # Input PSF image and (optional) diffusion kernel" % (file_psf)
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (file_badpix)
	print "G) %s                # File with parameter constraints (ASCII file)" % (file_constraints)
	print "H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx, ny)
	print "I) %i    %i              # Size of the convolution box (x y)" % (50,50)
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2sec,pix2sec)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"


	if h!=99999. and z0!=99999. and m0d!=99999. and Md==99999.: 
	 #m0d = -2.5*log10(I0d)
	 print "# Edge-on disk function\n"
 	 print "0) edgedisk            # Object type"
 	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
 	 print "3) %.3f       1       # central surface brightness  [mag/arcsec^2]" % (m0d)
 	 print "4) %.3f       1       # disk scale-height    [Pixels]" % (z0)
	 print "5) %.3f       1       # disk scale-length    [Pixels]" % (h)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
 	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n"
	 print "C0) %.3f       1       # diskiness/boxiness" % (cd)
	 print "To) 2                  # Outer truncation by component 2"
	 #print "T0) radial	 #  truncation" 
	 #print "T3) %.3f       0      #  Trancation radius-break beginning" % (rmax_d)
	 print "T4) %.3f       1      #  Trancation radius-break end" % (rtr)
	 print "T5) %.3f      1          #  Softening length (1 normal flux) [pixels]" % (5.)

	if h!=99999. and z0!=99999. and m0d==99999. and Md!=99999.: 
	 #m0d = -2.5*log10(I0d)
	 print "# Edge-on disk function\n"
 	 print "0) edgedisk1            # Object type"
 	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
 	 print "3) %.3f       1       # central surface brightness  [mag/arcsec^2]" % (Md)
 	 print "4) %.3f       1       # disk scale-height    [Pixels]" % (z0)
	 print "5) %.3f       1       # disk scale-length    [Pixels]" % (h)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
 	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n"
	 print "C0) %.3f       1       # diskiness/boxiness" % (cd)
	 print "To) 2                  # Outer truncation by component 2"
	 print "T0) radial	 #  truncation" 
	 #print "T3) %.3f       0      #  Trancation radius-break beginning" % (rmax_d)
	 print "T4) %.3f       1      #  Trancation radius-break end" % (rtr)
	 print "T5) %.3f      1          #  Softening length (1 normal flux) [pixels]" % (5.)


	if h!=99999. and z0==99999. and m0d!=99999. and Md==99999.:
	 print "# Exponential disk function\n"
	 print "0) expdisk1            # Object type"
	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (m0d)   
	 print " 4) %.3f       1       #     R_s              [Pixels]" % (h)
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-elld)  
	 print " 10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

	if h!=99999. and z0==99999. and m0d==99999. and Md!=99999.:
	 print "# Exponential disk function\n"
	 print "0) expdisk            # Object type"
	 print "1) %.3f  %.3f  1 1     # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (Md)   
	 print " 4) %.3f       1       #     R_s              [Pixels]" % (h)
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-elld)  
	 print " 10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAd)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"


	if reb!=99999. and meb!=99999. and Mb==99999.:
	 print "# Sersic function\n"
	 print " 0) sersic2             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (meb)   
	 print " 4) %.3f       1       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       1       # Sersic exponent (deVauc=4, expdisk=1)" % (n) 
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 #if cb!=99999.:
	 #	print "C0) %.3f       1       # diskiness/boxiness" % (cb)

	if reb!=99999. and meb==99999. and Mb!=99999.:
	 print "# Sersic function\n"
	 print " 0) sersic             # Object type"
	 print " 1) %.3f  %.3f  1 1    # position x, y        [pixel]" % (xc,yc)
	 print " 3) %.3f       1       # total magnitude"  % (Mb)   
	 print " 4) %.3f       1       #     R_e              [Pixels]" % (reb)
	 print " 5) %.3f       1       # Sersic exponent (deVauc=4, expdisk=1)" % (n) 
	 print " 9) %.3f       1       # axis ratio (b/a)" % (1.-ellb)  
	 print "10) %.3f        1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PAb)
	 print " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
	 #if cb!=99999.:
	 #	print "C0) %.3f       1       # diskiness/boxiness" % (cb)

	if sky_level!=0.:
	 print "# sky\n"
	 print " 0) sky                    #  object type"
	 print " 1) %.3f      1          #  sky background at center of fitting region [ADUs]" % (sky_level)
	 print " 2) 0.0000      0          #  dsky/dx (sky gradient in x)"
	 print " 3) 0.0000      0          #  dsky/dy (sky gradient in y)"
	 print " Z) 0                 #  Skip this model in output image?  (yes=1, no=0)"

 	print "\n================================================================================"                   

	sys.stdout = tmp_out
	f.close()

	#1
	#p = subprocess.Popen('galfit input.txt',  shell=True, close_fds=False)
	#p.kill
	#status = p.wait()

	#2
	#RunCmd(["galfit", "input.txt"], 60).Run()

	#3



	if status!=0:
		return status
	else:
		subprocess.call("galfit -o3 galfit.01", shell=True, close_fds=False)
		return status

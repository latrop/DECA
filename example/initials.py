#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Exponential disk
m0_d = 19.38		# Observed central SB of the disk: [mag arcsec^-2]	(should be >0 if magDisk=0)
magDisk = 99999.	# Apparent magnitude of the disk:  [mag]		(should be >0 if m0_d=0)
h_d = 16.		# Exponential scale of the disk: [arcsec]
z0 = 99999.		# Vertical scale of the disk: [arcsec]			(=0 for not edge-on disks)
ell_d = 0.		# Ellipticity of the disk				(=0 for edge-on disks)
c_d = 0.		# Ellipse index of the disk				(=0 for face-ons and =-1 for edge-ons) 
r_tr = 99999.		# Minimum radius: [arcsec]


#Sersic Bulge
me_bul = 22.		# Effective SB of the Sersic model: [mag arcsec^-2]	(should be >0 if magBulge=0) 
magBul = 99999.		# Effective radius of the Sersic model: [mag]		(should be >0 if me_bul=0)
re_bul = 15.		# Effective radius of the Sersic model: [arcsec] 
n_bul = 1.2		# Sersic index
ell_bul = 0.3		# Ellipticity
c_bul = 0.		# Ellipse index						(=0 for classical bulges and >0 for pseudo bulges)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

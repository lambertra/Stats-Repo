# Script that bins up object positions using elliptical bins
# Plots bin radius vs surface density as well as the expected surface density based on
# given params: ellipticity, position angle, half-light radius, and 
# background surface density. Area of CCD and del_RA and del_Dec are also needed
# for the binning.

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import gammainc
from scipy.special import gamma
from matplotlib.patches import Ellipse

def Main(argv):
	if len(sys.argv) != 2:
		print '''Can't read input.
		python e_r_bins.py input_file'''
		exit(1)
	else: file = sys.argv[1]
	
	print "Are your ellipticity and position angle correct?"
	
	# Leave sersic index free or not. 0 = Free and 1 = Fixed
	nswitch = 1
	
	#Get params
	##Position angle starting on positive y axis and moving counter-clockwise
	if nswitch == 0:
		ecc,pa,r_h,n,g_x,g_y = np.loadtxt('params.txt',unpack=True)
		E_bg = 6.818
	elif nswitch == 1:
		ecc,pa,r_h,g_x,g_y = np.loadtxt('params_nfix.txt',unpack=True)
		n = 1 
		E_bg = 6.818
	
	area = (11001*.11/60)**2 #Area of the CCD in arcminutes

	
	r_e = r_h
	#r_e = r_h/1.68  #Effective radius. Only need to change 1.68 to something else if 
					# n!=1. The term you're looking for can be approximated as 
					# b_n = 1.9992*n - .3271 as long as n is between 0.5 and 10
	print "Ell = {0:0.2f} PA = {1:1.2f}".format(ecc,pa)
	
	x,y = np.loadtxt(file,unpack=True) #Coordinates of objects with respect to 
									   # reported galaxy position. Coordinates should be
									   # in arcminutes.
# 	x=(x-8.997500)*60 * np.cos(51.559722/180*np.pi)
# 	y=(y-51.559722)*60
#  	x,y = np.loadtxt('LacI_rgb.dat',unpack=True)
#  	x = (x - 344.55)*np.cos((41.29)/180*np.pi)*60
#  	y = (y - 41.29)*60

	
	#####Matching Dave's Plot####
		#Dave's MLE Galaxy coord's
		#((344.55525-344.55 )*np.cos((41.2981944-41.29)/180*np.pi),41.2981944-41.29)#
		

	#g_x,g_y = ((344.55525-344.55 )*np.cos((41.29)/180*np.pi)*60,(41.2981944-41.29)*60)

# 	draw_x,draw_y = np.loadtxt('/Users/ryalambe/Downloads/good_draw_laci.txt',unpack=True)
# 	draw_x = 60*(draw_x - 344.55)*np.cos(41.29/180*np.pi)
# 	draw_y = 60*(draw_y - 41.29)
# 
# 	d_ellip_r = np.sqrt(
#                       (  (1 - ecc) ** -1 
#                       * ((draw_x - g_x) * np.cos(pa) 
#                        - (draw_y - g_y) * np.sin(pa)))**2
#                       + ((draw_x - g_x) * np.sin(pa)
#                         +(draw_y - g_y) * np.cos(pa))**2
#                     )
# 	avgfunc = np.sum(np.exp(-d_ellip_r/r_e))/len(d_ellip_r)
# 	surfdensint = avgfunc*area
# 	S0 = (len(x) - E_bg*area)/surfdensint
	
	#############################

	#g_x,g_y = (0.2500768701790329,0.45622963522871207) #Del_RA and Del_Dec respectively.
										#Assumes x and y coordinates are with respect to
										#the galaxy center as the origin
	N = len(x)
	ellip_r = np.sqrt(
                      (  (1 - ecc) ** -1 
                      * ((x - g_x) * np.cos(pa) 
                       - (y - g_y) * np.sin(pa)))**2
                      + ((x - g_x) * np.sin(pa)
                        +(y - g_y) * np.cos(pa))**2
                    )

#Create bins based on elliptical radii of points. You can uncomment the bins line
#to base the histogram on your most distant point, but the script does not account for
#elliptical bins so large that areas of them fall outside the CCD space.
	bins = np.linspace(0.,10,21)
	bins_r = np.zeros(len(bins)-1)  #Array for bin radius
	density = np.zeros(len(bins)-1) #Array for bin surface density
	err = np.zeros(len(bins)-1)		#Array for uncertainties.


	for i in xrange(1,len(bins)):
		area_e = (1-ecc)*(bins[i]**2 - bins[i-1]**2) * np.pi
		density[i-1] = (len(ellip_r[(ellip_r <= bins[i]) & (ellip_r>bins[i-1])]) 
						/ area_e ) 
		err[i-1] = (np.sqrt(len(ellip_r[(ellip_r <= bins[i]) & (ellip_r>bins[i-1])]))
					/area_e)
		
		bins_r[i-1] = (bins[i] + bins[i-1])/2

	#bins_r = np.loadtxt('/Users/ryalambe/Downloads/surf1d_area.txt',usecols=(0))
	fig = plt.figure()
	ax = plt.subplot()
	plt.errorbar(bins_r,density,yerr = err,fmt='kD', capsize=0.1,markeredgecolor='k',
				 marker='D',fillstyle= 'none',ls = 'none',fs = None)

# Uses an exponential sersic profile to estimate the surface density at each radius
# 	N_star = N - E_bg*area
# 	S0 = N_star/2/np.pi/r_e**2/(1-ecc) * 1.68**2
	b_n = 1.9992*n - .3271
	max_rad = np.max(ellip_r)

# 	area = np.pi * (1-ecc) * max_rad**2
	S0 = ((N - E_bg * area) * b_n**(2*n) 
                 / (r_e**2 * 2 * np.pi * n * (1 - ecc))
                 / (
                      gammainc(2*n,b_n*(max_rad / r_e)**(1./n)) 
                    * gamma(2*n)))
	

####Was using to estimate expected number of points in 1 - 1.75 r_e range
####Not relevant to code goals
# 	rad = np.linspace(1*r_e,1.75*r_e,30)
# 	surf_ = S0 * np.exp(-b_n * (rad / r_e)**(1./n)) + E_bg
# 	num = np.zeros(surf_.shape)
# 	for i in xrange(len(surf_)):
# 		if i!=0:
# 			num[i] = surf_[i]*(1-ecc)*np.pi*(rad[i]**2 - rad[i-1]**2)
# 	print np.sum(num)
			
	if nswitch==0: 
		label = 'Sersic profile'
	else: label = 'Exponential profile'
	plt.plot(np.linspace(0,10,1000),
				S0*np.exp(-b_n*(np.linspace(0,10,1000)/r_e)**(1./n)) + E_bg,'k-',label=label)
# 	plt.vlines([r_e,1.75*r_e],0,plt.gca().get_ylim()[1])
	plt.hlines(E_bg,0.1,10,linestyles='dashed',label='Background Surface Density')
	plt.xscale('Log')
	plt.yscale('Log')
	plt.xlabel('Elliptical Radius [arcmin]',size = 10)
	plt.ylabel('Surface Density [arcmin$^{-2}$]', size=10)
	plt.legend(loc='center left',frameon=False)
	min,max = plt.gca().get_ylim()
# 	plt.vlines(3,0,max,color='red',linestyle='dashed')
# 	ax.text(.1666666,12.2957,'''Per I fit with n = 1''')
	plt.gca().set_xlim([0.1,10])
	plt.gca().set_ylim([min,max])
	plt.show()

	
	plt.plot(x,y,'k.')
	plt.gca().add_artist(Ellipse((g_x,g_y),2*r_e,2*r_e*(1-ecc),90-pa/np.pi*180, fill=False,
						color = 'red',ls='-'))
# 	for i in xrange(1,len(bins)):
# 		e_ = Ellipse((g_x,g_y),2*i,2*i*(1-ecc),90-pa/np.pi*180, fill=False,
# 						color = 'red',ls='-')
# 		plt.gca().add_artist(e_)
	plt.show()
Main(sys.argv)
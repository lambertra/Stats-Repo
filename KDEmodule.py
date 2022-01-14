####A script that uses Kernel Density Estimation to locate over-densities within a 
####set of data. Requires a comparison data set; usually a uniformly generated
####random set of data.
import numpy as np
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import sys
from matplotlib import colors
from matplotlib import rc
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
import matplotlib.patheffects as path_effects
import sig_finder
import astropy.wcs as wcs
from astropy.io import fits



###Ugly way of priming the variables from the command line
if len(sys.argv) == 3:
	data = sys.argv[1]
	coor = sys.argv[2]
	bd = 'cv_ml'
elif len(sys.argv) == 4:
	data = sys.argv[1]
	coor = sys.argv[2]
	bd = np.array([np.float(sys.argv[3]),np.float(sys.argv[3])])
elif len(sys.argv) == 5:
	data = sys.argv[1]
	coor = sys.argv[2]
	bd = np.array([np.float(sys.argv[3]),np.float(sys.argv[4])])

else: print 'python KDEmodule.py Data_file coordinates(RADEC or Pix) \
bandwidth_x bandwidth_y (optional)'; exit(1)


#Set coordinates 
while coor !='RADEC' and coor != 'Pix':
	radec = ['RADEC','radec','ra','dec','r']
	pix = ['PIX','pix','Pix','pixel','p']
	if coor in radec: coor = 'RADEC'
	elif coor in pix: coor = 'Pix'
	else:
		coor = raw_input('Position coordinates in radec or pixel: ')


data = np.loadtxt(data,usecols=(0,1))	

if len(data[0,:]) > 2:
	data = data[:,:2]

###########Changeable variables########
Nx = 150
Ny = 150
trials = 500
p_s = 0.258 #pixel scale
r_num = 413 ##Expected number of background points
#######################


#######Create grid and background simulations##########
xg,yg = np.loadtxt('gal_center.txt',unpack=True)
xaxis,yaxis=np.loadtxt('CCD_axes.txt',unpack=True)
xmin,xmax=(0,xaxis)
ymin,ymax=(0,yaxis)
rand = np.random.uniform(0,1,(trials*r_num,2))


##Convert CCD coordinates to RADEC
##Have to do this before making grid
if coor == 'RADEC':

	hdulist=fits.open('n4406_v_gcs_2010_masked.fits',ignore_missing_end=True)
	w = wcs.WCS(hdulist[0].header,hdulist)
	xg,yg = w.all_pix2world(xg,yg,1,ra_dec_order=True)
	ra1,dec1 = w.all_pix2world(1,1,1,ra_dec_order=True)
	ra2,dec2 = w.all_pix2world(xaxis,yaxis,1,ra_dec_order=True)

	if ra1 > ra2: xmin = ra2; xmax = ra1
	else: xmin = ra1; xmax = ra2
	if dec1 > dec2: ymax = dec1; ymin = dec2
	else: ymax = dec2; ymin = dec1
	
# 	xmin = 344.376		#LacI coordinates
# 	xmax = 344.823
# 	ymax = 41.461
# 	ymin = 41.127
	
	rand[:,0] = (rand[:,0]*(xmax-xmin)) + xmin 
	rand[:,1] = (rand[:,1]*(ymax-ymin)) + ymin
	bd = np.array([200. * .26 / 3600,200*.26/3600])
    


Xgrid=np.vstack(map(np.ravel,
					np.meshgrid(np.linspace(float(xmin+(xmax-xmin)/Nx/2.),
						     				float(xmax-abs(xmax-xmin)/Nx/2.),Nx),
								np.linspace(float(ymin+(ymax-ymin)/Ny/2.),
											float(ymax-abs(ymax-ymin)/Ny/2.),Ny)))).T

##Convert pixel coordinates to arcminutes and center on galaxy	
if coor == 'Pix':
	data[:,0] = (data[:,0] - xg) * (p_s / 60)
	data[:,1] = (data[:,1] - yg) * (p_s / 60)
	Xgrid[:,0] = (Xgrid[:,0] - xg) * (p_s / 60)
	Xgrid[:,1] = (Xgrid[:,1] - yg) * (p_s / 60)
	rand[:,0] = (rand[:,0]*xmax - xg) * (p_s / 60)
	rand[:,1] = (rand[:,1]*ymax - yg) * (p_s / 60)
###

print 'Running KDE on galaxy data'
kde = KDEMultivariate(data=data,var_type='cc',bw=bd)

if coor == 'RADEC': 
	print "Bandwidth size used is: {0:f}\' x {1:f}\'".format(kde.bw[0]*60,kde.bw[1]*60)
else: print "Bandwidth size used is: {0:f}\' x {1:f}\'".format(kde.bw[0],kde.bw[1])

pdf = kde.pdf(Xgrid)


##Section to determine significance
##One day will add a switch on command line for it.
# bg_pdf = []
# print 'Running KDE on simulated background data'
# print 'This may take some time...'
# for i in xrange(trials):
# 
# 	kde = KDEMultivariate(data=rand[i*r_num:(i+1)*r_num,:],var_type='cc',bw = bd)
# 	# kde.fit(rand[i*r_num:(i+1)*r_num,:])
# 	bg_pdf.append(kde.pdf(Xgrid))
# 	if ((i*10) % trials) == 0:
# 		print '{0:d}% finished'.format(100*i/trials)
# 
# bg_pdf = np.array(bg_pdf).reshape((Ny,Nx,trials)) * float(r_num) / len(data[:,0])
# 
# cdf = sig_finder.Sig_Finder(bg_pdf.flatten(),pdf)
# cdf = cdf.reshape((Ny,Nx))

####Taking the average background subtracted KDE
bg_ = []
tots = 10
for i in xrange(tots):
	mask = np.random.choice(xrange(len(data[:,0])),r_num,replace=False)
	mask_array = np.zeros(len(data[:,0]))
	mask_array[mask] = 1
	mask_array = np.vstack((mask_array,mask_array)).T
	new_data = np.ma.masked_array(data,mask=mask_array)
	kde = KDEMultivariate(data=new_data,var_type='cc',bw=bd)
	bg_.append(kde.pdf(Xgrid))
bg_ = np.array(bg_).reshape((tots,Ny,Nx)) 

bg_pdf = np.zeros((Ny,Nx))
for i in xrange(Ny):
	for j in xrange(Nx):
		bg_pdf[i,j] = np.mean(bg_[:,i,j,])
		
# pdf = bg_pdf.reshape((Ny,Nx))

pdf = bg_pdf

dens = pdf * (len(data[:,0]) - 413)
dens = dens.reshape((Ny,Nx))


if coor == 'RADEC':
	dens = dens / 3600 
   ###Density was GCCs/deg^2. Now GCCs/arcminute^2


	

###It's possible to undo the stacking and raveling, but it's easier on readability 
###to just remake the meshgrid
Xgrid = np.meshgrid(np.linspace(float(xmin+(xmax-xmin)/Nx/2.),
								float(xmax-abs(xmax-xmin)/Nx/2.),Nx),
					np.linspace(float(ymin+(ymax-ymin)/Ny/2.),
								float(ymax-abs(ymax-ymin)/Ny/2.),Ny))


if coor == 'Pix':	
	Xgrid[0] = (Xgrid[0] - xg) * (p_s / 60)
	Xgrid[1] = (Xgrid[1] - yg) * (p_s / 60)
	xmin = (xmin - xg) * (p_s / 60)
	ymin = (ymin - yg) * (p_s / 60)
	xmax = (xmax - xg) * (p_s / 60)
	ymax = (ymax - yg) * (p_s / 60)

#####Figure Section#########

##set to latex font
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex = True)

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)
levels = np.linspace(np.min(dens.flatten()),np.max(dens.flatten()),12)
# levels = np.linspace(0,np.max(dens.flatten()),6)


img = plt.imshow(dens,extent=(xmin,xmax,ymin,ymax),origin = 'lower', 
					interpolation='nearest')

cs = plt.contour(Xgrid[0],Xgrid[1],dens,levels[2:],colors = 'black',zorder=3)

##Contour of the significance line. Isn't used unless the bg_pdf section is uncommented.
# sig_level = [.99]
# cs2 = plt.contour(Xgrid[0],Xgrid[1],cdf,sig_level,colors = '#d55e00',linewidths = 2.0)

if coor == 'RADEC':
	plt.xlabel('RA [Degrees]',size = 15)
	plt.ylabel('Dec [Degrees]',size = 15)
	ax.invert_xaxis()

if coor == 'Pix':
	plt.xlabel(r'X [arcmin]',size = 15)
	plt.ylabel(r'Y [arcmin]',size = 15)


cbar = plt.colorbar(img,format = '%1.2f',fraction=0.046, pad=0.04)
cbar.set_label(r'GCC Surface Density [GCC / arcmin$^{2}]$',size=12)




# plt.plot(data[:,0],data[:,1],'r.')
plt.show()
plt.clf()
plt.close()
# 	plt.savefig('n7331_{0:1.2f}_{1:1.2f}.png'.format(bd[0],bd[1]))






########M86 and M84 Model Stuff

################
#Subtracting from "models" like the ref wants
# data1 = np.loadtxt('N4406_gen.txt')
# pdf1 = []
# pdf2 = []
# m_pdf1 = np.zeros((Ny,Nx))
# m_pdf2 = np.zeros((Ny,Nx))
# for i in xrange(100):
# 	kde = KDEMultivariate(data=data1[i*1250:(i+1)*1250,:],var_type='cc',bw=bd)
# 	pdf1.append(kde.pdf(Xgrid).reshape((Ny,Nx)))
# 
# # pdf -= m_pdf
# data2 = np.loadtxt('N4374_gen.txt')
# for i in xrange(100):
# 	kde = KDEMultivariate(data=data2[i*1250:(i+1)*1250,:],var_type='cc',bw=bd)
# 	pdf2.append(kde.pdf(Xgrid).reshape((Ny,Nx)))
# 	
# for i in xrange(Ny):
# 	for j in xrange(Nx):
# 		m_pdf1[i,j] = np.mean(np.array(pdf1)[:,i,j])
# 		m_pdf2[i,j] = np.mean(np.array(pdf2)[:,i,j])

# m_dens1 = m_pdf1 * (len(data1[:,0])/100 ) * .42
# m_dens2 = m_pdf2 * (len(data2[:,0])/100 ) * .3
# 
# dens -= m_dens1
# dens -= m_dens2
###Adding contours for both M86 and M84 KDE
# data = np.loadtxt('N4406_gen.txt')
# data[:,0] = data[:,0] - xg * (p_s / 60)
# data[:,1] = data[:,1] - yg * (p_s / 60)
# kde = KDEMultivariate(data=data,var_type='cc',bw=bd)
# m_pdf = kde.pdf(Xgrid)
# m_pdf = m_pdf.reshape((Ny,Nx))
# m_pdf[m_pdf > np.max(m_pdf)/6.] = 0
# levels1 = np.linspace(np.min(m_pdf.flatten()),np.max(m_pdf.flatten()),18)
# 
# data = np.loadtxt('N4374_gen.txt')
# data[:,0] = data[:,0] - xg * (p_s / 60)
# data[:,1] = data[:,1] - yg * (p_s / 60)
# kde = KDEMultivariate(data=data,var_type='cc',bw=bd)
# m_pdf1 = kde.pdf(Xgrid)
# m_pdf1 = m_pdf1.reshape((Ny,Nx))
# m_pdf1[m_pdf1 > np.max(m_pdf1)/6.] = 0
# levels2 = np.linspace(np.min(m_pdf.flatten()),np.max(m_pdf.flatten()),18)
#################
###Piece of script for marking galaxies in N4406 Field
# m86 = [186.54946504566234,12.94611906928599, -224]
# m84 = [186.265597,12.886983, 1017]
# n4387 = [186.423668,12.810515, 565]
# n4388 = [186.444780,12.662086, 2524]
# n4402 = [186.531890,13.113320, 232]
# n4425 = [186.805551,12.734723, 1908]
# ic3303 = [186.313342,12.714618,-188]
# n4406b = [186.562774,12.964110,1101]	
# gal_data = np.vstack((m86,m84,n4387,n4388,n4402,n4425,n4406b,ic3303))
# plt.plot(gal_data[:,0],gal_data[:,1],marker = 'P',
# 			 mec='black',mfc='white',ms = 8,ls='none',zorder=10)
# coords = [[186.589,12.9731,'A'],[186.421,12.9024,'B'],[186.667,12.828,'C']]
# fs = 18
# for i in xrange(len(coords[0])):
#  ###White on black text effect for easy reading
# 	txt = plt.text(coords[i][0],coords[i][1],coords[i][2],family = 'serif',fontsize=fs,
# 				weight = 'medium',color = 'w')
# 	txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
###Piece of script used to mask portions of the imshow plot
# dens[dens < np.min(dens)/4.] = 0 #This line has to go before the first imshow.
#.
#.
#.
# mask = np.zeros((150,150,4))
# mask[dens==0] = [0,0,0,1]
# for i in xrange(71,76):
# 	for j in xrange(72,76):
# 		mask[i,j] = [0,0,0,1]
# plt.imshow(mask,extent=(xmin,xmax,ymin,ymax),origin = 'lower', 
# 					interpolation='nearest')

###Making the velocity scatter plot with galaxy positions
# class MidpointNormalize(colors.Normalize):
#     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#         self.midpoint = midpoint
#         colors.Normalize.__init__(self, vmin, vmax, clip)
# 
#     def __call__(self, value, clip=None):
#         # I'm ignoring masked values and all kinds of edge cases to make a
#         # simple example...
#         x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#         return np.ma.masked_array(np.interp(value, x, y))
# 
# vel_data = np.loadtxt('all_gcc_vel_with_err.dat')
# ra = vel_data[:,0]
# dec = vel_data[:,1]
# vel = vel_data[:,2]
# # # ra = (ra - xg) 
# # # dec = (dec - yg) 
# 
# norm = MidpointNormalize(midpoint = 0)
# 
# sc = plt.scatter(ra,dec,c=vel,cmap='RdBu_r',norm=norm,
# 				# vmin = -1* np.max(vel),vmax=np.max(vel),
# 				edgecolors='k',zorder=8)
# g_sc = plt.scatter(gal_data[:,0],gal_data[:,1],c = gal_data[:,2], s = 169,
# 					cmap = 'RdBu_r',norm=norm,marker= '*',edgecolors = 'k',
# 					linewidths = 1.2, zorder=7)
# cbar = plt.colorbar(sc,fraction=0.045, pad=0.04)
# print cbar.get_ticks()
# cbar.ax.set_yticklabels(['<-500','0','500','1000','1500','2000','>2500'])
# # cbar.ax.set_yticklabels(['<-800.', '-600.', '-400.', '-200.', 
# # '0.', '200.', '400.', '600.', '>800.'])
# # cbar.ax.set_yticklabels(['<0', '500', '1000', '1500', '2000', '>2500'])
# cbar.set_label('GC Radial Velocity [km/s]')
# 
# 
# fig.subplots_adjust(right = 0.78)
# fig.subplots_adjust(top = 0.8)
# x_axis = ax.get_xlim()
# # ax.set_xlim(x_axis[1],x_axis[0])
# plt.show()
# plt.clf()
# plt.close()


#############
###Regions plotter for Stephen's work
# def Regions():
#    regions = np.genfromtxt('Mr38sP1regionsfile_edit.txt',delimiter=',',filling_values=0)
# 
# 
#    #Regions we care about: 1-9
#    R_Dict = {
#    1:regions[0,:][regions[0,:]!=0].reshape(len(regions[0,:][regions[0,:]!=0])/2,2),
#    2:regions[1,:][regions[1,:]!=0].reshape(len(regions[1,:][regions[1,:]!=0])/2,2),
#    3:regions[2,:][regions[2,:]!=0].reshape(len(regions[2,:][regions[2,:]!=0])/2,2),
#    4:regions[3,:][regions[3,:]!=0].reshape(len(regions[3,:][regions[3,:]!=0])/2,2),
#    5:regions[4,:][regions[4,:]!=0].reshape(len(regions[4,:][regions[4,:]!=0])/2,2),
#    6:regions[5,:][regions[5,:]!=0].reshape(len(regions[5,:][regions[5,:]!=0])/2,2),
#    7:regions[6,:][regions[6,:]!=0].reshape(len(regions[6,:][regions[6,:]!=0])/2,2),
#    8:regions[7,:][regions[7,:]!=0].reshape(len(regions[7,:][regions[7,:]!=0])/2,2),
#    9:regions[8,:][regions[8,:]!=0].reshape(len(regions[8,:][regions[8,:]!=0])/2,2)}
# 
#    hdulist=fits.open('Mr38sP1un21212ms.fits',ignore_missing_end=True)
#    w = wcs.WCS(hdulist[0].header,hdulist)
#   
#    polygons = []
#    for i in xrange(9):
#        wcs_coord = w.all_pix2world(R_Dict[i+1],1,ra_dec_order=True)
#        dummy = np.array([wcs_coord[-1,:],wcs_coord[0,:]])
#        plt.plot(wcs_coord[:,0],wcs_coord[:,1],color='#f0e442',ls='-',zorder=4)
#        plt.plot(dummy[:,0],dummy[:,1],color='#f0e442',ls='-',zorder=4)
#            
# 
# Regions()
# plt.xlim((xmax,xmin))
# plt.ylim((ymin,ymax))
# plt.show()



#################
#Uses the biweight location and scale estimators                 
# def Biweight_Loc_est(vel):
#     #vel = np.loadtxt(input_file,usecols=(2))
#     ERR = 10000
# #     vel = np.loadtxt('normal_dis_velocities.txt')
# #     vel = np.random.normal(-200,300,10000)
# 
#     median = np.median(vel)
#     loop_count = 0
#     while ERR >= 0.0001:
#        
#         m_diff = abs(median - vel)
#         mad = np.median(m_diff)
#     
#         mu = (vel - median) / 6. / mad
#       
#         bw_loc= (  median 
#                  + np.sum((vel-median)*(1-mu**2)**2)
#                  / np.sum((1-mu**2)**2))
#        
#     
#         mad = mad/0.6745
#         bw_sca = ( np.sqrt(len(vel)
#                            *  np.sum((vel - median)**2 * (1-((vel-median)/9./mad)**2)**4))
#                   /abs(np.sum((1-((vel-median)/9./mad)**2)*(1 - 5*((vel-median)/9./mad)**2))))
#     
#         
#         ERR = abs(median - bw_loc)
#         median = bw_loc
#         loop_count+=1
#         if loop_count == 1:
#             print bw_loc,bw_sca
#         if loop_count > 100:
#             print 'Loop exceeded maximum number of iterations. Convergence too slow \
#             or not possible.'
#             break
#         
#    	return(bw_loc,bw_sca)
# 
# ######################
# def Vel_Dis_Calc(input_file):
#     m86 = [186.54946504566234,12.94611906928599]
#     m84 = [186.265597,12.886983]
#     emp_Re = [6.4, 6.6]
#     ra,dec,vel,err = np.loadtxt(input_file,unpack=True)
#     
#     rad_wrt_m86 = np.sqrt((ra - m86[0])**2 + (dec - m86[1])**2) * 60
#     rad_wrt_m84 = np.sqrt((ra - m84[0])**2 + (dec - m84[1])**2) * 60
#     
#     m86_gcs = (rad_wrt_m86 < emp_Re[0])
#     m84_gcs = (rad_wrt_m84 < emp_Re[1])
#     
#     m86_w_avg,m86_w_std = Biweight_Loc_est(vel[m86_gcs])
#     m84_w_avg,m84_w_std = Biweight_Loc_est(vel[m84_gcs])
# # #     m84_weights = 1./err[m84_gcs]**2
# # #     m86_weights = 1./err[m86_gcs]**2
# # #      
# # #     m84_w_avg = np.average(vel[m84_gcs],weights = m84_weights)
# # #     m86_w_avg = np.average(vel[m86_gcs],weights = m86_weights)
# # #     
# # #     m84_w_var = (   np.sum(m84_weights * (vel[m84_gcs] - m84_w_avg)**2) 
# # #                  / (np.sum(m84_weights)) ) 
# # #     m86_w_var = (   np.sum(m86_weights * (vel[m86_gcs] - m86_w_avg)**2)
# # #                  / (np.sum(m86_weights)) )
# # #     
# # #     m84_w_std = np.sqrt(m84_w_var)
# # #     m86_w_std = np.sqrt(m86_w_var)
# # 
#     m86gc_std = np.zeros(len(ra))
#     m84gc_std = np.zeros(len(ra))
#     bothgc_std = np.zeros(len(ra))
#     for i in xrange(len(ra)):
#        m86gc_std[i] = abs(vel[i] - m86_w_avg) / m86_w_std
#        m84gc_std[i] = abs(vel[i] - m84_w_avg) / m84_w_std
#        if rad_wrt_m86[i] <= rad_wrt_m84[i]:
#            bothgc_std[i] = abs(vel[i] - m86_w_avg) / m86_w_std
#        else: 
#            bothgc_std[i] = abs(vel[i] - m84_w_avg) / m84_w_std
# # 
#     print 'M86 Stats:'
#     print 'weighted average: {0:f}  weighted stand. dev.: {1:f}'.format(m86_w_avg,m86_w_std)
#     print 'Total GCs: {0:d}'.format(len(vel[m86_gcs]))
#     
#     print 'M84 Stats:'
#     print 'weighted average: {0:f}  weighted stand. dev.: {1:f}'.format(m84_w_avg,m84_w_std)
#     print 'Total GCs: {0:d}'.format(len(vel[m84_gcs]))
#     print np.max(bothgc_std)
#   	
#     for i in xrange(len(ra)):
#   	    print '{0:3.3f} {1:2.3f} {2:1.3f}'.format(ra[i],dec[i],bothgc_std[i])
#   	
#     ax = plt.subplot(111)
# #     ax.invert_xaxis()
#     sc = plt.scatter(ra,dec,c=bothgc_std,cmap = 'Reds',
#                 edgecolors='k',zorder=8)
#     cbar = plt.colorbar(sc,fraction=0.045, pad=0.04)
#     cbar.set_label('$\sigma_{V}$ from Average Velocity of Closest Massive Elliptical')
#     plt.show()
# # 
# Vel_Dis_Calc('all_gcc_vel_with_err.dat')










# if coor == 'RADEC':
# 	
# 	m86 = [186.54946504566234,12.94611906928599]
# 	m84 = [186.265597,12.886983]
# 	n4387 = [186.423668,12.810515]
# 	n4388 = [186.444780,12.662086]
# 	n4402 = [186.531890,13.113320]
# 	n4425 = [186.805551,12.734723]
# 	
# 	m84_c = plt.Circle((m84[0],m84[1]),6.6/60.,color = 'y',ls='-',fill=False)
# 	m86_c = plt.Circle((xg,yg),6.4/60.,color = 'y',ls='-',fill=False)
# 	n4387_c = plt.Circle((n4387[0],n4387[1]),0.5/60.,color = 'k',ls='-',fill=False)
# 	n4388_c = plt.Circle((n4388[0],n4388[1]),1.2/60.,color = 'k',ls='-',fill=False)
# 	n4402_c = plt.Circle((n4402[0],n4402[1]),4.4/60.,color = 'k',ls='-',fill=False)
# 	n4425_c = plt.Circle((n4425[0],n4425[1]),2.0/60.,color = 'k',ls='-',fill=False)
# 	plt.plot(data[:,0],data[:,1],'r.')
# 	ax.add_artist(m84_c)
# 	ax.add_artist(m86_c)
# 	ax.add_artist(n4387_c)
# 	ax.add_artist(n4388_c)
# 	ax.add_artist(n4402_c)
# 	ax.add_artist(n4425_c)
# 	m86 = [186.54946504566234,12.94611906928599, -224]
# 	m84 = [186.265597,12.886983, 1017]
# 	n4387 = [186.423668,12.810515, 565]
# 	n4388 = [186.444780,12.662086, 2524]
# 	n4402 = [186.531890,13.113320, 232]
# 	n4425 = [186.805551,12.734723, 1908]
# 	
# 	ra,dec,vel,err = np.loadtxt('all_gcc_vel_with_err.dat',unpack=True)
# 	gal_data = np.vstack((m86,m84,n4387,n4388,n4402,n4425))
# 	radii = np.zeros((len(ra),6))
# 	for i in xrange(len(ra)):
# 		radii[i] = np.sqrt((ra[i] - gal_data[:,0])**2 + (dec[i] - gal_data[:,1])**2)
# 	argsort = np.argsort(radii,axis = -1)
# 	sig = np.zeros(len(ra))
# 	for i in xrange(len(ra)):
# 		gal_vel = gal_data[argsort[i,0], 2]
# 		if i <= 92:
# 			if vel[i] < 400. : 
# 				sig[i] = np.abs(gal_vel - vel[i]) / 172.
# 			else: sig[i] = np.abs(gal_vel - vel[i]) / 314.
#         
# 		else: sig[i] = np.abs(gal_vel - vel[i]) / 292.
# 
# 	sc = plt.scatter(ra,dec,c=sig,vmax = 5,cmap='Reds', edgecolors='k')
# 	plt.colorbar(sc,fraction=0.045, pad=0.04)
# # 	ax.invert_xaxis()
# 	plt.show()
# plt.plot(.6,.3,marker='*',markeredgecolor='k',markerfacecolor='w') #M86
# txt = plt.text(.6,1.1,'M86',family = 'serif', weight = 'medium',color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
#                        
# plt.plot(-15,-12.9,marker='+',markeredgecolor='w',markerfacecolor='w') #N4425
# txt = plt.text(-16,-14.7,'NGC 4425',family = 'serif', weight = 'medium',color='w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# 
# plt.plot(6.6,-17.1,marker='+',markeredgecolor='w',markerfacecolor='w') #N4388
# txt = plt.text(6,-16.2,'NGC 4388',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
#                       
# plt.plot(7.2,-8.1,marker='+',markeredgecolor='w',markerfacecolor='w') #N4387
# txt = plt.text(6.6,-9.9,'NGC 4387',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# 
# plt.plot(1.8,9.5,marker='+',markeredgecolor='w',markerfacecolor='w') #N4402
# txt = plt.text(1.2,10.3,'NGC 4402',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# 
# plt.plot(15.8,-3.8,marker='+',markeredgecolor='k',markerfacecolor='w') #M84
# txt = plt.text(15.2,-5.6,'M84',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# 
# 
# txt = plt.arrow(-8,6,0,5,edgecolor = 'k',facecolor = 'w')
# # txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
# #                       path_effects.Normal()])
# txt = plt.text(-8.5,12,'N',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# 
# txt = plt.arrow(-8,6,-5,0,edgecolor = 'k',facecolor = 'w')
# # txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
# #                       path_effects.Normal()])
# txt = plt.text(-15,5.5,'E',family = 'serif', weight = 'medium', color = 'w')
# txt.set_path_effects([path_effects.Stroke(linewidth = 3, foreground = 'black'),
#                       path_effects.Normal()])
# plt.show()
#########
# vmin = np.min(dens)
# vmax = np.max(dens)
# data = np.loadtxt('red_n4406.dat')
# xg,yg = np.loadtxt('gal_center.txt',unpack=True)
# xaxis,yaxis=np.loadtxt('CCD_axes.txt',unpack=True)
# xmin,xmax=(0,xaxis)
# ymin,ymax=(0,yaxis)
# if len(data[0,:]) > 2:
# 	data = data[:,:2]
# Xgrid=np.vstack(map(np.ravel,
# 					np.meshgrid(np.linspace(float(xmin+(xmax-xmin)/Nx/2.),
# 						     				float(xmax-abs(xmax-xmin)/Nx/2.),Nx),
# 								np.linspace(float(ymin+(ymax-ymin)/Ny/2.),
# 											float(ymax-abs(ymax-ymin)/Ny/2.),Ny)))).T
# 
# ##Convert pixel coordinates to arcminutes and center on galaxy	
# if coor == 'Pix':
# 	data[:,0] = (data[:,0] - xg) * (p_s / 60)
# 	data[:,1] = (data[:,1] - yg) * (p_s / 60)
# 	Xgrid[:,0] = (Xgrid[:,0] - xg) * (p_s / 60)
# 	Xgrid[:,1] = (Xgrid[:,1] - yg) * (p_s / 60)
# 	rand[:,0] = (rand[:,0]*xmax - xg) * (p_s / 60)
# 	rand[:,1] = (rand[:,1]*ymax - yg) * (p_s / 60)
# 	bd = 200. * p_s / 60
# 
# print 'Running KDE on galaxy data'
# kde = KDE(bandwidth = bd)#grid.best_estimator_
# kde.fit(data) 
# pdf = np.exp(kde.score_samples(Xgrid))
# 
# pdf = pdf.reshape((Ny,Nx))
# dens = pdf * len(data[:,0])
# 
# ###It's possible to undo the stacking and raveling, but it's easier on readability 
# ###to just remake the meshgrid
# Xgrid = np.meshgrid(np.linspace(float(xmin+(xmax-xmin)/Nx/2.),
# 								float(xmax-abs(xmax-xmin)/Nx/2.),Nx),
# 					np.linspace(float(ymin+(ymax-ymin)/Ny/2.),
# 								float(ymax-abs(ymax-ymin)/Ny/2.),Ny))
# 
# ####Conversion from Pixel to Arcmin
# if coor == 'Pix':	
# 	Xgrid[0] = (Xgrid[0] - xg) * (p_s / 60)
# 	Xgrid[1] = (Xgrid[1] - yg) * (p_s / 60)
# 	xmin = (xmin - xg) * (p_s / 60)
# 	ymin = (ymin - yg) * (p_s / 60)
# 	xmax = (xmax - xg) * (p_s / 60)
# 	ymax = (ymax - yg) * (p_s / 60)
# ax = plt.subplot(212)	
# img = plt.imshow(dens,extent=(xmin,xmax,ymin,ymax),origin = 'lower', 
# 					cmap = 'Reds',interpolation='nearest',vmin = vmin, vmax = vmax)
# levels = np.linspace(np.min(dens.flatten()),np.max(dens.flatten()),12)
# cs = plt.contour(Xgrid[0],Xgrid[1],dens,levels,cmap = 'gray',zorder=3)
# 
# cbar = plt.colorbar(img,format = '%1.2f',fraction=0.046, pad=0.04)
# cbar.set_label(r'Red GCC Surface Density',size=12)
# if coor == 'Pix':
# 	plt.xlabel(r'X [arcmin]',size = 15)
# 	plt.ylabel(r'Y [arcmin]',size = 15)
# plt.show()
# plt.clf()
# plt.close()

# fig.subplots_adjust(right = 0.78)
# fig.subplots_adjust(top = 0.8)
# 
# plt.show()
# plt.clf()
# plt.close()


# ####Making the colormap	
# cdict = {'red':   [(0.0,   1.0, 1.0),
#                    (0.25, 0.75, 0.75),
#                    (0.75, 0.25, 0.25),
#                    (1.0,  0.0, 0.0)],
# 
#          'green': [(0.0,  1.0, 1.0),
#                    (0.25, 1.0, 1.0),
#                    (0.75, 0.75, 0.75),
#                    (0.85, 0.55, 0.55),
#                    (0.95, 0.25, 0.25),
#                    (1.0,  0.0, 0.0)],
# 
#          'blue':  [(0.0,   1.0, 1.0),
#                    (0.25, 0.25, 0.25),
#                    (0.50, 0.15,0.15), 
#                    (1.0,  0.0, 0.0)]}
# 
# cmap=colors.LinearSegmentedColormap('Upgrade',cdict,N=256,gamma=1.0)
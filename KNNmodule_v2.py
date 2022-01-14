import numpy as np
from astroML.density_estimation import KNeighborsDensity as knn
import sys
import matplotlib.pyplot as plt

if len(sys.argv)==4:
	data = sys.argv[1]
	rand = sys.argv[2]
	method = sys.argv[3]

else:
	print "python KNNmodule.py data_file random_data_file grid/point"
	exit(2)
if method != 'grid' and method != 'point':
	print 'Density estimation method not understood. Choices are \'grid\' and \'point\''
	exit(1)


#####Bookkeeping variables that may need to be changed####
Nx=50
Ny=50
trials= 5000
############

##Xgrid reads across the x axis first before moving up one row
xaxis,yaxis=np.loadtxt('CCD_axes.txt',unpack=True)
xmin,xmax=(0,xaxis)
ymin,ymax=(0,yaxis)
Xgrid=np.vstack(map(np.ravel,
					np.meshgrid(np.linspace(float(abs(xmin-xmax)/Nx/2.),
						     				float(xmax-abs(xmax-xmin)/Nx/2.),Nx),
								np.linspace(float(abs(ymax-ymin)/Ny/2.),
											float(ymax-abs(ymax-ymin)/Ny/2.),Ny)))).T
#################


real_data=np.loadtxt(data)
rand_data=np.loadtxt(rand)
prob = real_data[:,2]
real_data = np.delete(real_data,2,1)


####Determine the radius of each pixel so that we can compare pixels with similar radii
# try:
# 	gal_x,gal_y = np.loadtxt('gal_center.txt',unpack = True)
# 	
# except IOError: 
# 	print '	!!!!WARNING!!!!: Galaxy center text file not available'
# 	print '			 Using the middle of the CCD'
# 	gal_x = xaxis/2.
# 	gal_y = yaxis/2.
# 
# pix_rad = np.sqrt( (float((xmax-xmin)/Nx))**2 + (float((ymax-ymin)/Ny))**2 ) / 2.
# 
# radial_grid = np.zeros((Ny*Nx))
# for i in xrange(Ny*Nx):
# 	radial_grid[i] = np.sqrt((Xgrid[i,0] - gal_x)**2 + (Xgrid[i,1] - gal_y)**2)
#############



#####Manually reshaping the array since nothing else works. Gross.
counter=0
reshaped_data=np.zeros((trials,len(rand_data[:,0])/trials,len(rand_data[0,:])))

for i in xrange(1,trials+1):
	for j in xrange(1,len(rand_data[:,0])/trials+1):
		for k in xrange(len(rand_data[0,:]/trials)):
			reshaped_data[i-1,j-1,k]=rand_data[counter,k]
		counter+=1
###########		

###KNN Loop####
for kth in range(5,11):
	poisson=[]
	for i in xrange(trials+1): ##Trials+1 to get all of the random data and the real data
		if i == 0:
			if method == 'grid':
				density=knn('bayesian',kth).fit(real_data).eval(Xgrid)
			else: density = knn('bayesian',kth).fit(real_data).eval(real_data)
		else: 
			if method == 'grid':
				rand_density=knn('bayesian',kth).fit(reshaped_data[i-1]).eval(Xgrid)
														###i-1 so that it doesn't go
														###out of bounds
			else: 
				rand_density = knn('bayesian',kth).fit(reshaped_data[i-1]).eval(
														reshaped_data[i-1])
			
			poisson.append(rand_density)
	poisson = np.array(poisson)
	flat_poisson=poisson.flatten()

	cdf = []
	for i in xrange(len(density)):
		cdf.append(len(flat_poisson[density[i] > flat_poisson]) / float(len(flat_poisson)))
	cdf = np.array(cdf)


###Significance for contour plots
	for l in xrange(10):
		cut_den = density[cdf >= 0.95]
		cut_cdf = cdf[cdf >= 0.95]
		argsort = np.argsort(cut_cdf)
		expected_density = cut_den[argsort][0]
		num_of_points = int(round(xaxis*yaxis*expected_density))
		rand_sim = np.random.uniform(0,1,(trials*num_of_points,2))
		rand_sim[:,0] *= xaxis
		rand_sim[:,1] *= yaxis
		poisson = []

	
		counter=0
		reshaped_data=np.zeros((trials,len(rand_sim[:,0])/trials,len(rand_sim[0,:])))

		for i in xrange(1,trials+1):
			for j in xrange(1,len(rand_sim[:,0])/trials+1):
				for k in xrange(len(rand_sim[0,:]/trials)):
					reshaped_data[i-1,j-1,k]=rand_sim[counter,k]
				counter+=1
	
		for i in xrange(trials):
			if method == 'grid':
				rand_density=knn('bayesian',kth).fit(reshaped_data[i]).eval(Xgrid)
														###i-1 so that it doesn't go
														###out of bounds
			else: 
				rand_density = knn('bayesian',kth).fit(reshaped_data[i]).eval(
													reshaped_data[i])
	
			poisson.append(rand_density)
	
		poisson = np.array(poisson)
		flat_poisson=poisson.flatten()
	
		for i in xrange(len(cdf[cdf >= 0.95])):
			new_percent = len(flat_poisson[density[i] > flat_poisson]) / float(len(flat_poisson))
			cdf[cdf >= 0.95][i] = 1 - (1-cdf[cdf>=0.95][i])*(1-new_percent)
		
	
####Print to output files
	if kth>=10:
		name='knn_result_{0:d}.txt'.format(kth)
	else: name = 'knn_result_0{0:d}.txt'.format(kth)
	print 'Printing into {0:20s}'.format(name)
	f=open(name,'w+')
	print >> f,'#Nx  Ny'
	print >> f, "#%d %d" %(Nx,Ny)
	print >> f, "#Result   CDF	Significance"
	for i in xrange(len(density)):
		print >> f, '{0:16.14f} {1:16.14f}'.format(density[i],cdf[i])
		
	f.close()

###Kuzma significance criterion###	
# 	significance = []	
# 	for i in xrange(len(cdf)):
# 		significance.append(0)
# 		if cdf[i] >= 0.95:
# 			point_rad = np.sqrt( ( real_data[:,0] - Xgrid[i,0] ) ** 2 
# 								+( real_data[:,1] - Xgrid[i,1] ) ** 2 )
# 					
# 			sort_index = np.argsort(point_rad,axis=0)
# 			###Sorts the indices by radius so that I can get the sorted matrix by 
# 			###Plugging sort_index into the array. ex: point_rad[sort_index]
# 			
# 			od_points = real_data[sort_index[:kth]]
# 			mask = np.ones(point_rad.shape,dtype=bool)
# 			mask[sort_index[:kth]] = 0
# 			prob_pool = prob[mask]
# 
# 			# avg_od_prob = np.mean(prob[sort_index[:kth]])
# # 			avg_prob_pool = np.mean(prob_pool)
# # 			std_prob_pool = np.std(prob_pool)
# # 			significance.append((avg_od_prob - avg_prob_pool) / std_prob_pool)
# 			number_of_high_weight = []
# 			for j in xrange(5000):
# 				rand = np.random.randint(0,len(prob_pool),kth)
# 				number_of_high_weight.append(len(rand[prob_pool[rand] >= 0.2]))
# 			
# 			significance.append(( len(od_points[prob[sort_index[:kth]] >= 0.2]) -
# 							 	  np.mean(number_of_high_weight) ) 
# 							 	/ np.std(number_of_high_weight))
# 			
# 		else: significance.append(0)
# 	significance = np.array(significance)
# 	significance[significance < 0] = 0
	#significance[significance!=0] += abs(np.min(significance))			
#	significance = np.zeros(Nx*Ny)


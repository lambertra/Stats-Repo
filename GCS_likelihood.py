# A python script that takes the spatial coordinates of a data set and 
# uses Maximum-likelihood Estimation to determine the best-fitting parameters
# for a Sersic model (or exponential) and write them to an output file

import numpy as np
import scipy.optimize as op
import emcee
import matplotlib.pyplot as plt
import corner
from scipy.special import gammainc
from scipy.special import gamma
import sys
from statsmodels.nonparametric.kernel_density import KDEMultivariate


####switch determines whether n is a free param or not. 0=yes, 1=no
nswitch=1

###Priors
ecc = 0.45
pa = np.pi / 2.
r_e = 3.
n = .6
E_bg = 6.818
#6.818 = LacI, 2.64 = CasIII, 2.88/3.2 = PerI


pix_scale = 0.11


#Sersic function
def likelihood(p,f,nswitch):
    if nswitch == 1:
        ecc,pa,r_e,g_x,g_y= p	#Galaxy position, eccentricity, position
        x,y,N,area,max_rad,E_bg,n = f	#angle, half-light radius, background density
    else:
        ecc,pa,r_e,n,g_x,g_y= p
        x,y,N,area,max_rad,E_bg = f       

    
#     rand = np.random.choice(len(x),size=(int(round(area*E_bg))),replace=False)
# 
#     mask = np.zeros(len(x))
#     mask[rand] = 1
#     x = np.ma.masked_array(x,mask)
#     y = np.ma.masked_array(y,mask)
#     N = len(x) - int(round(area*E_bg))
# Eqn from Martin, De Jong, Rix 2008
    
    ellip_r = np.sqrt(
                      (  (1 - ecc) ** -1 
                       *((x - g_x) * np.cos(pa) 
                       - (y - g_y) * np.sin(pa)))**2
                      + ((x - g_x) * np.sin(pa)
                        +(y - g_y) * np.cos(pa))**2
                     )
 	
 	
    b_n = 1.9992*n - .3271
    ##Normalization term comes from Graham and Driver, 2005
    norm_term = ((N - E_bg * area) * b_n**(2*n) 
                 / (r_e**2 * 2 * np.pi * n * (1 - ecc))
                 / (gammainc(2*n,b_n*(max_rad / r_e)**(1./n))
                 *  gamma(2*n))
                )


    fcn = (
     		np.sum(np.log(norm_term 
                        * np.exp(-b_n * ((ellip_r/r_e)**(1./n)))
                        + E_bg
                         )
                 )
           )


########David Sand's LL fcn
#     N_ = N - area * E_bg
#     S0 = N_ / (2 * np.pi * (r_e / b_n)**2 * (1 - ecc))
#     fcn = np.sum(np.log(S0 * np.exp(-b_n* (ellip_r / r_e)) + E_bg))
    

#     b_n = 1.9992*n - .3271
# #     N = N - E_bg*area
#     fcn = (-2 * N * np.log(r_e) 
#           + N * np.log(b_n**(2*n) / n 
#           			  / (gammainc(2*n,b_n * (max_rad/r_e)**(1./n)) * gamma(2*n))
#           			  )
#           - N * np.log(2 * np.pi)
#           - N * np.log(1-ecc)
#           - b_n * np.sum(ellip_r ** (1/n)) / r_e ** (1/n)
# #           + np.log(E_bg)
#           )

    if ( 
        np.min(x) < g_x < np.max(x) and np.min(y) < g_y < np.max(y) and 
        0 <= ecc < 1   and 0 <= pa < np.pi and# 0 < r_e < 10
       .3 < n < 10    and 0 < r_e < 10 
        #and 0 <= E_bg < N/area
        ):
        return(fcn)
    else: return(-np.inf)
    


#Priors - Uninformative      
def prior(p,f,nswitch):
    if nswitch == 1:
        ecc,pa,r_e,g_x,g_y= p
        x,y,N,area,max_rad,E_bg,n = f
    else:
        ecc,pa,r_e,n,g_x,g_y= p
        x,y,N,area,max_rad,E_bg = f
    if ( 
        np.min(x) < g_x < np.max(x) and np.min(y) < g_y < np.max(y) and 
        0 <= ecc < 1   and 0 <= pa < np.pi and
       .3 < n < 10    and 0 < r_e < 10 
       #and 0 <= E_bg < N/area
        ):
        return(0.0)
    else: return(-np.inf)


def prob(p,f,nswitch):
    probability = prior(p,f,nswitch)
    if not np.isfinite(probability):
        return(-np.inf)
    return(probability + likelihood(p,f,nswitch))



if len(sys.argv)==2:
    file = sys.argv[1]
else:
    print 'python GCS_likelihood.py data_file'
    exit(1)	
  
x,y = np.loadtxt(file,usecols=(0,1),unpack=True)

N = len(x) 


g_x,g_y = np.loadtxt('gal_center.txt',unpack=True)



######## E_bg variables #####################
max_rad = np.max(np.sqrt((x-g_x)**2 + (y-g_y)**2))

xax,yax = np.loadtxt('CCD_axes.txt',unpack=True)
area = (yax * xax) * (pix_scale / 60) ** 2 


if nswitch == 1:
    p = np.array([ecc,pa,r_e,g_x,g_y])
    f = np.array([x,y,N,area,max_rad,E_bg,n])
else: 
    p = np.array([ecc,pa,r_e,n,g_x,g_y])
    f = np.array([x,y,N,area,max_rad,E_bg])


nll = lambda *args: -likelihood(*args)

if nswitch == 1:
    result = op.minimize(nll,p,args = (f,nswitch), bounds = [(0,0.99),(0,np.pi),
    											(0.1,9.9),(np.min(x),np.max(x)),
    											(np.min(y),np.max(y))])
else: result = op.minimize(nll,p,args = (f,nswitch), bounds = [(0,0.99),(0,np.pi),
									(0.1,9.9),(0.5,10),(np.min(x),np.max(x)),
									(np.min(y),np.max(y))])


if nswitch == 1:
    ecc,pa,r_e,g_x,g_y = result["x"]
    print ecc,pa/np.pi*180,r_e,E_bg		##PA is measured clockwise if the x axis isn't 
else:									##flipped to account for RA decreasing
    ecc,pa,r_e,n,g_x,g_y = result["x"]  ##as the plot moves to the right.
    print ecc,pa/np.pi*180,r_e,n,E_bg,g_x,g_y


ndim,nwalkers = len(p),100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers,ndim,prob,args=(f,nswitch))
sampler.run_mcmc(pos,10000)

samples = sampler.chain[:,3000:,:].reshape((-1,ndim))

for i in xrange(len(p)):
	print (np.percentile(samples[:,i],50),
		   np.percentile(samples[:,i],50)-np.percentile(samples[:,i],16),
		   np.percentile(samples[:,i],84) - np.percentile(samples[:,i],50))








if nswitch == 1:
    fig = corner.corner(samples,labels = ["ecc","pa","r_e","del_ra","del_dec"])
else: fig = corner.corner(samples,labels = ["ecc","pa","r_e","n","del_ra","del_dec"])

plt.show()
if nswitch == 1:
    f = open('params_nfix.txt','w+')
    print >> f, '#ell, pa, r_h, del_ra,del_dec'
else: f = open('params.txt','w+'); print >>f, '#ell, pa, r_h, n, del_ra, del_dec, E_bg'
for i in xrange(len(p)):
    print >> f, np.median(samples[:,i])
f.close()


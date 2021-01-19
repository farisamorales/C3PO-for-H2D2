


#----- System Information for testing this program.
if __name__ == '__main__':
	import sys
	print "\n\nName of Script: ", sys.argv[0]	
	print "Number of arguments: ", len(sys.argv)
	print "The arguments are: " , str(sys.argv), "\r\n"


#----- Declare primary function 'flxin'
def flxin(Tstar, Rstar, Gsize, Rdist):
	
	#print "Gsize is:", Gsize, "\n\n"
	
	from scipy import integrate, pi, e, log10
	import numpy as np	
	import pickle
	from astropy import units as u, constants as const,analytic_functions as func
	
        Tstar_K = Tstar*u.K
        Rstar_m = Rstar* 1.4959787066e11*u.AU
	Gsize_m = Gsize* 1.e-6 *u.m
	Rdist_m = Rdist* 1.4959787066e11*u.AU

#----- Test to be sure that the proper function parameters have been passed from main code.
	#if not (type(Gsize) == int or type(Gsize) == float) or Gsize == 0:
	#	print "Invalid parameter entry, inflx.flxin, Gsize, at index", index
 	#if not (type(Tstar) == int or type(Tstar) == float) or Tstar == 0:
	#	print "Invalid parameter entry, Tstar, at index", index 
 	#if not (type(Rstar) == int or type(Rstar) == float) or Rstar == 0:
	#	print "Invalid parameter entry, Rstar, at index", index
 	#if not (type(Rdist) == int or type(Rdist) == float) or Rdist == 0:
	#	print "Invalid parameter entry, Tstar, at index", index	

#----- Read in wavelengths and emissivities from data pickle.
	#waves = np.array([],dtype=np.float128)
	waves = pickle.load(open('suwav.txt','rb'))
	eps= 0.9 #Can start practicing reading in from astrosil tables? 

#----- Make astropy constants easier to write.
	c = const.c		# c: Speed of light in vacuum	Units: m/s
	k = const.k_B		# k: Boltzman Constant 		Units: J/K
	h = const.h		# h: Planck Constant		Units: m2 kg / s
	S = const.sigma_sb	# sigma: Stefan-Boltzman const	Units: W/(m2*K4)

	c1 = (2.*h/c**2).decompose()
	c2 = (h/k).decompose()

#----- Begin to Calculate BB curve for Stellar Radiation
	
	freq = np.array([u.micron.to(u.Hz, waves, equivalencies=u.spectral())*u.Hz],dtype=np.float128)
	val = np.array([((freq*c2/Tstar).decompose())/u.K], dtype=np.float128)	#Divide by K to strip Unit

	flux = ((c1 * freq**3 / (np.exp(val)-1.)).decompose()/(u.kg/u.s**2))/u.s**3

#----- Integrating the BB Curve over Freq 
	newflux = flux * freq * eps
	intflux = sorted(newflux[0,0,:])
	intfreq = sorted(log10(freq[0,:]))

	Lin = integrate.simps(intflux,intfreq)

#----- Calculate Energy In!
	
	Ein = ((Gsize_m/(2.*Rdist_m))**2)*(4*pi*Rstar_m**2) * Lin


	lam = np.array([u.micron.to(u.m,waves)])

#------ Plotting Flux over Wavelength for Analyzing program.
        if __name__== "__main__":
		L = np.squeeze(lam[0,:])
		I = np.squeeze(flux[0,0,:])
                Lp = L.astype(np.float32)*10**6 
		Ip = I.astype(np.float32) *10**26
		catch = np.where(Ip)
		import matplotlib.pyplot as plt
		plt.xlim(10**-1,10**3)
		plt.ylim(10**15,10**19)
                plt.loglog(Lp,Ip)
		plt.title('Flx/Lam SED IN')

                plt.show()


#----- Additional Testing Zone
	tempEin = Ein/(u.m)**2
	return tempEin

if __name__ == '__main__':

	print "Value returned by flxin = ", flxin(6e3,1,1,2)

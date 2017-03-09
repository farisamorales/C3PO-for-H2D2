#--------------------MAJOR NOTES REGARDING THIS CODE CURRENTLY
# I've been noodling with the astropy units, at some points in this code you'll find that I randomly divide the units out in order to have a dimensionless value. Moreover, the final value returned in the flxout function is in m2, where it should be in flx units. 

# If you run this program it will try to bring up a log plotted SED for the grain. It won't return the value for Energy in (Ein) until you close the graph window.

# eps (emissivity) is currently set to .9. The code will eventually have to deal with an array of Cody's epsilon values for different materials... for now .9




#----- System Information for testing this program.
if __name__ == '_main__':
	import sys
	print "\n\nName of Script: ", sys.argv[0]
	print "Number of arguments: ", len(sys.argv)
	print "The arguments are: " , str(sys.argv), "\r\n"


#----- Declare primary function 'flxout'
def flxout(Gsize, Tout):

#----- Test to be proper parameters have been passed.
#	if not (type(Gsize) == int or type(Gsize) == float) or Gsize == 0:
#		print ("Enter valid value for Gsize in microns (int, float, non-zero)")
#	if not (type(Tout) == int or type(Tout) == float) or Tout == 0:
 #               print ("Enter valid value for Tout in Kelvin (int, float, non-zero)")


#----- Import additional modules and functions
        from scipy import integrate, pi, e, log10
        import numpy as np
        import pickle
        from astropy import units as u, constants as const,analytic_functions as func

        Tout_K = Tout*u.K
        Gsize_m = Gsize* 1.e-6 *u.m

#----- Read in wavelengths and emissivities from data pickle.
        waves = pickle.load(open('suwav.txt','rb'))

	eps = .9

#----- Make astropy constants less cumbersome.
        c = const.c             # c: Speed of light in vacuum   Units: m/s
        k = const.k_B           # k: Boltzman Constant          Units: J/K
        h = const.h             # h: Planck Constant            Units: m2 kg / s

        c1 = (2.*h/c**2).decompose()
        c2 = (h/k).decompose()
    
#----- Begin to Calculate BB curve for Stellar Radiation

        freq = np.array([u.micron.to(u.Hz, waves, equivalencies=u.spectral())*u.Hz],dtype=np.float128)
        val = np.array([((freq*c2/Tout).decompose())/u.K], dtype=np.float128)  #Divide by K to strip Unit

        flux = ((c1 * freq**3 / (np.exp(val)-1.)).decompose()/(u.kg/u.s**2))/u.s**3

#----- Integrating the BB Curve over Freq 
        newflux = flux * freq * eps
        intflux = sorted(newflux[0,0,:])
        intfreq = sorted(log10(freq[0,:]))

        Lout = integrate.simps(intflux,intfreq)

#----- Calculate Energy Out!
	Eout = (4.*pi*Gsize_m**2) * Lout


#----- TEST ZONE

	lam = np.array([u.micron.to(u.m,waves)])

        if __name__== "__main__":
                L = np.squeeze(lam[0,:])
                I = np.squeeze(flux[0,0,:])
                Lp = L.astype(np.float32)
                Ip = I.astype(np.float32)
                catch = np.where(Ip)
                
		import matplotlib.pyplot as plt
                plt.loglog(Lp,Ip)
		plt.xlim(10**-7,1)
                plt.ylim(10**-18)
		plt.title('Flx/lam SED OUT')

                plt.show()
	tempEout = Eout/(u.m)**2
	return tempEout

if __name__ == "__main__":
	print "Value returned by flxout =", flxout(1,200)

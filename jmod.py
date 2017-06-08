


#******************************************************************************
#---- jmod is currently in essense a small library of functions specifically designed to be importable by other modules in the overall C3PO code. To use jmod as intended one need simply import it into their code, then call jmod.tempout. tempout being the main function of the module will automatically call the other subsequent functions. tempout takes two parameters, Rstar(radius of the star in AU.), and Tstar(temperature of the star in K)
#******************************************************************************
from scipy import integrate, pi, e, log10, interpolate, constants as const
import numpy as np
from kmod import listifier, whiteSpaceParser
import pickle
import collections
#----- Listifier and whitespaceparser are from Cody's Kmod. 



#----- declare global scipy constants. I assume many of the constants will have been declared before this module is actually called, but because I can't know how they'll be globally named I simply restate them here specifically for this module.

c = const.c             # c: Speed of light in vacuum   Units: m/s
c_um = c * 1.e6         # c: Speed of light in vacuum   Units: micron/s
k = const.k             # k: Boltzman Constant          Units: J/K
h = const.h             # h: Planck Constant            Units: m2 kg / s
S = const.sigma         # sigma: Stefan-Boltzman const  Units: W/(m2*K4)

c1 = (2.*h/c**2)
c2 = (h/k)

#here I define tempout as the first step function in the module, all other functions are defined within tempout to maintain global naming conventions and ease o data passage between the seperate functions.

def tempout(Rstar,Tstar):

#readEmiss should return a tuple consisting of Lambda (List of wavelengths) and Emiss (List of lists) first level of Emiss returns grain sizes, second level is a list of emissivities corresponding with a certain wavelength in Lambda ( Emiss[*][i] <=> Lambda[i-1] )


    def readEmiss(FILE):
        FILE.seek(0)
        Lamda = listifier(whiteSpaceParser(FILE.readline()))
        Grainsizes = []
        Emiss = []
        for ele in FILE:
            temp = listifier(whiteSpaceParser(ele))

# the first element of each line is the grain size
	    Grainsizes.append(temp[0])
# drop the first element; the rest of the line is emissivity values
	    Emiss.append(temp[1:])

        return Lamda, Grainsizes, Emiss

#flxin takes four parameters Tstar(K), Rstar(AU), Gsize(mic),Rdist(AU). Then, using the list of wavelengths and grainlists returned from reademiss, it calculates the flux headed into each individual grain.

    def flxin(Tstar, Rstar, Gsize, Rdist):
            Rstar_m = Rstar* 6.957e8          #R_Sun to m
            Gsize_m = Gsize* 1.e6             #um to m
            Rdist_m = Rdist* 1.4959787066e11  #au to m


#----- Find corresponding emiss and wave lists from this particular csv file.
            waves = emisTable['wave']       #microns
	    epslist=[1]
            for i, ele in enumerate (grainsizes):
		    if ele == Gsize:
                            epslist = emisTable['emis'][i]
            eps=epslist

#----- Calculate Flux in
            freq = np.array(c_um/waves,dtype=np.float64)	    
            val = np.array(freq*c2/Tstar, dtype=np.float64)
            flux = (c1 * freq**3 / (np.exp(val)-1.))

#----- Integrating the BB Curve over Freq 
            epsflux = flux * eps
	    #BB: for blackbody temperatures, uncomment the next line
            #epsflux = flux

# the frequency array goes the wrong direction, so add a minus sign here
            Lin = -integrate.simps(epsflux,freq)

	    if 0:
	      LinBB = -integrate.simps(flux,freq)
	      print 'LinBB',LinBB

	      #  THIS IS A GOOD THING TO CHECK:
	      print 'grain size, optical albedo',Gsize,(LinBB - Lin)/LinBB

	      #  THIS IS A GOOD TEST OF THE INTEGRATOR ACCURACY
              print 'T star real, T calculated',Tstar,(np.pi*LinBB/S)**.25

#----- Calculate Energy In!
            Ein = ((Gsize_m/(2.*Rdist_m))**2)*(4*pi*Rstar_m**2) * Lin
            return Ein

#flxout uses individual grain sizes from the list of grainsizes populated in reademiss, and an arbitrary Tout(K) to calculate the energy leaving a grain.
    def flxout(Gsize, Tout):
            Gsize_m = Gsize* 1.e6       #mic to m

#----- Find Corresponding eps and lam lists for this loop.
            waves = emisTable['wave']       #microns
            #print 'wavelength range',waves[0][0],waves[0][-1]
            
	    epslist=[1]
            for i, ele in enumerate (grainsizes):
                    if ele == Gsize:
                            epslist = emisTable['emis'][i]
            eps = epslist

            freq = np.array(c_um/waves, dtype=np.float64)
            val = np.array(freq*c2/Tout, dtype=np.float64)
            flux = (c1 * freq**3 / (np.exp(val)-1.))    
            
#----- Integrating the BB Curve over Freq 
            newflux = flux * eps
	    #BB: for blackbody temperatures, uncomment the next line
            #newflux = flux
	    
            intflux = newflux
            intfreq = freq
# the frequency array goes the wrong direction, so add a minus sign here
            Lout = -integrate.simps(intflux,intfreq)
	    #BB: if emis=1, this is a good check for integrator
            #print 'T dust real, T calculated',Tout,(np.pi*Lout/S)**.25

#----- Calculate Energy Out!
            Eout = (4.*pi*Gsize_m**2) * Lout

            return Eout

    
#----- Read in the emissivity list we'd like to use for this run.

    (wave,size,emis) = readEmiss(open('Emissivities[DI;AS;50-50].csv'))
    #(wave,size,emis) = readEmiss(open('Emissivities[0.1;1.0;10.0].csv'))

    emisTable={}
    emisTable['wave']=np.array(wave)
    emisTable['grainsize']=np.array(size)
    emisTable['emis']=np.array(emis)
        
    outtable = {}

#----- Declare lists and test values for interpolation loops.
    grainsizes=emisTable['grainsize']

# Consider a range of orbital locations
    radii = np.logspace(0,2,5)

# Set a range of temperatures for flux_in/flux_out calculation
#  (if temperature goes outside this range, grain temp becomes NaN)
    Tlist = np.logspace(0,3,20)

    print 'Rstar,Tstar',Rstar,Tstar

    templist = []
    Tdusts= []

#Begin out put loops
    for r in radii:
            print "\nrdist:", r

	    Tdusts.append([]) 
	    
            for g in grainsizes: 
                    valin = flxin(Tstar, Rstar, g, r)
                    valsout = []
	            
                    for T in Tlist:
                            tempvalout = flxout(g,T)
                            valsout.append(tempvalout)

		    #Interpolating in logspace for more accuracy
                    logTg = np.interp(np.log10(valin), np.log10(valsout), np.log10(Tlist))
		    Tg = 10.**logTg
                    Tdusts[-1].append(Tg)
		    templist.append(Tg)

#################################TESTING TESTER#########################################
		    try:
			    Tdusts.index(Tg) #then what?
		    except ValueError: #flip except else? Check online to format.
			    outtable[Tg] = {'grainsize':g,'radius':r}
	    	    else:
			    exit( 'DUPLICATE TEMPERATURE FOUND: THIS WILL NOT BE ADDED TO THE DICTIONARY BUT WILL BE ADDED TO THE LIST CANNOT DO, BRUH!')

##################################TESTIN OVER###########################################

		   outtable[Tg] = {'grainsize':g,'radius':r} 
    return templist,grainsizes,radii, outtable

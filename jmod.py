

from scipy import integrate, pi, e, log10, interpolate, constants as const
import numpy as np
from kmod import listifier, whiteSpaceParser

#----- declare global scipy constants.
c = const.c             # c: Speed of light in vacuum   Units: m/s
k = const.k           # k: Boltzman Constant          Units: J/K
h = const.h             # h: Planck Constant            Units: m2 kg / s
S = const.sigma         # sigma: Stefan-Boltzman const  Units: W/(m2*K4)

c1 = (2.*h/c**2)
c2 = (h/k)

def readEmiss(FILE):
	FILE.seek(0)
	Lamda = listifier(whiteSpaceParser(FILE.readline()))
	Emiss = []
        
	for ele in FILE:
		temp = listifier(whiteSpaceParser(ele))
		Emiss.append(temp)
	return Lamda, Emiss

def flxin(Tstar, Rstar, Gsize, Rdist):

        Tstar_K = Tstar
        Rstar_m = Rstar* 1.4959787066e11
        Gsize_m = Gsize* 1.e6 
        Rdist_m = Rdist* 1.4959787066e11

#----- Find corresponding emiss and wave lists from this particular csv file.
        tag = []
        eps = []
        waves = np.array([codyin[0]])
        count = 1
	Gsize = 10 #THIS IS A TEST VALUE- WILL NEED TO BE SURE IF STATEMENT BELOW FINDS ACTUAL VALUES
        while (count < len(codyin[1])):
                tag.append(count-1)
                count += 1
        for i, ele in enumerate(tag):
                if codyin[1][i-1][0] == Gsize:
                        epslist = codyin[1][i-1]
        for i, ele in enumerate(waves):
                eps.append(epslist[i+1])

#----- Calculate Flux in
        freq = np.array([c/waves],dtype=np.float64)
        val = np.array([freq*c2/Tstar], dtype=np.float64)
	flux = (c1 * freq**3 / (np.exp(val)-1.))

#----- Integrating the BB Curve over Freq 
        newflux = flux * freq * eps
        intflux = sorted(newflux[0,0,:])
        intfreq = sorted(log10(freq[0,0,:]))

        Lin = integrate.simps(intflux,intfreq)

#----- Calculate Energy In!
        Ein = ((Gsize_m/(2.*Rdist_m))**2)*(4*pi*Rstar_m**2) * Lin[0]
	return Ein

def flxout(Gsize, Tout):
        Tout_K = Tout
        Gsize_m = Gsize* 1.e6

#----- Find Corresponding eps and lam lists for this loop.
        tag = []
        eps = []
        waves = np.array([codyin[0]])
        count = 1
        Gsize = 10 #THIS IS A TEST VALUE- WILL NEED TO BE SURE IF STATEMENT BELOW FINDS ACTUAL VALUES


        while (count < len(codyin[1])):
                tag.append(count-1)
                count += 1
        for i, ele in enumerate(tag):
                if codyin[1][i-1][0] == Gsize:
                        epslist = codyin[1][i-1]
        for i, ele in enumerate(waves):
                eps.append(epslist[i+1])

        freq = np.array([c/waves],dtype=np.float64)
        val = np.array([freq*c2/Tstar], dtype=np.float64)
        flux = (c1 * freq**3 / (np.exp(val)-1.))

#----- Integrating the BB Curve over Freq 
        newflux = flux * freq * eps
        intflux = sorted(newflux[0,0,:])
        intfreq = sorted(log10(freq[0,0,:]))

        Lout = integrate.simps(intflux,intfreq)

#----- Calculate Energy Out!
        Eout = (4.*pi*Gsize_m**2) * Lout[0]
	return Eout

#------------------------------------------------------------------------------------------

#----- Read in the emissivity list we'd like to use for this run.
codyin = readEmiss(open('Emissivities[0.1;1.0;10.0].csv'))

#----- Declare lists and test values for interpolation loops.
valsout = []
graintemp = []

if __name__ == "__main__":
        radii = np.linspace(.1,40.,10)
        Tout = np.linspace(100,300,10)
        Rdist =1.
        Rstar = 1.
        Tstar = 6. *10**6

#This is a hard coded test value
grainlist = np.linspace(.001,1000,10)

for i, ele in enumerate(radii):
        rdist = ele
        print "rdist:", rdist, "\n\n"

        for j, eleg in enumerate(grainlist):
                valin = flxin(Tstar, Rstar, eleg, rdist)

                for k, elek in enumerate(Tout):
                        tempout = flxout(eleg,Tout[k])
                        valsout.append(tempout)
                Tg = (eleg, valin, np.interp(valin, Tout, valsout))
                graintemp.append(Tg)
                del valsout[:]
print Tg

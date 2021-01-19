
from scipy import interpolate
import matplotlib.pyplot as plt
import pickle
import inflx
import outflx
import numpy as np

valsout = []
graintemp = []
#radii =  np.linspace(1.,7.5,5)
radii = [2.3]
 
#Tout = radii *0.333333333333333*100
#Tout = np.linspace(100,300,10)
Tout = [200]

Rdist =1.
Rstar = 1.
Tstar = 6. *10**6
index = 0
count =0

#print np.logspace(100,300,50)
#print "Tout :", "\n", Tout, "\n\n"
#print "radii:", radii, "\n\n"
#print "Tout:", Tout, "\n\n"
#print "Temps:", temps, "\n\n"
grainlist = pickle.load(open('su_grain.txt','rb'))
#grains = np.array(grainlist,)
#print "Type - grains is",type(grains)


#grains = np.array([grainlist],)
#print"grains", grains

for i, ele in enumerate(radii):
	rdist = ele
	print "rdist:", rdist, "\n\n"

	for j, eleg in enumerate(grainlist):	
		valin = inflx.flxin(Tstar, Rstar, eleg, rdist)
		print valin
	
		for k, elek in enumerate(Tout):
			count += 1
			#print "\n\n", count
			tempout = outflx.flxout(eleg,elek)
			#print tempout
			valsout.append(outflx.flxout(eleg,elek))
			#print elek

		#print valsout
		#print "\n\n", len(valsout), len(Tout)
		
		#interpfunc = interpolate.interp1d(valsout, Tout, kind='linear')
		#Tg = interpfunc(valin)
		Tg = np.interp(valin, Tout, valsout)
		graintemp.append(Tg)
		del valsout[:]

print graintemp

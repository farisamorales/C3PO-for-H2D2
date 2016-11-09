# NAME:
#	FLXIN()
#
# DATE	December 2009 - Rewritten in Python November 2016
#
# PURPOSE:
#	To calculate the Energy or flux going INTO the grain
#	Use conservation of energy to determine the final grain Temperature
#	since there is no analytical way to solve for T in Planck function.
#
# CALLING SEQUENCE: 
#       Ein = FLXIN(Tstar,Rstar,Gsize,rdist)
#
# INPUT PARAMETERS:
#	TSTAR	Temperature of star in degree K
#	RSTAR	Radius of star in AU
#	GSIZE	Size of grain in um (single value)
#	RDIST	Distance from grain to star in AU
#
# OUTPUT PARAMETERS:
#       EIN - Scalar or vector giving the flux (i.e. Area*Fnu)
#               in Watts for a specified grain size, location, and
#		stellar properties.



#**********************************************************************************************************************
# NOTES:  There are a few points of the code that still need to be translated and a couple of snags that I've run into and have yet to unsnag.

# 1.) The section that reads in values from tables - We'll likely use the pickle approach to either load in CSV tables, or have Star Devourer arrange these CSVs into list arrays before my module interacts with them. The code currently doesn't read in values from tables. The variable 'waves' has been hard coded to hold a short list of values that have been specifically chosen to avoid floating point overflow, which is still an issue. (See issue #2) and eps (epsilon) has been set temporarily to 1. 
# 2.) The second bit of code that hasn't been translated is the machine precision section. Cody and I have been having a bit of trouble with that section. It's proven difficult to find a python analogue, moreover the IDL code there is a tad confusing. I still don't have a clear means of overcomming this obstacle but will continue to research and see if I can't find a way over that doesn't change the overall code, or math too dramatically.
# 3.) The third issue is in regards to ordering arrays. IDL has a nifty keyword that one can simply switch on when usng the int_tabulated function. The closest python analogue I could find was scypi.integrate.simps(), which does not have an integrated keyword for ordering arrays... However as I write now I'm remembering that it might have a keyword that deals with randomly ordered tables. Since the return value is just a floating point value it might not matter whether the parameter arrays for scypi.integrate.simps() are ordered, so long as their elements at each index actually correspond. I'll look into this.
#4.) Finally, I simply haven't coded the simple check at the beginning of the function that should test to make sure that appropriate parameters have been passed into the function. Comming soon.
# -JB
#**********************************************************************************************************************



#check that parameters have been passed for Gsize, Tstar, Rstar, Rdist					****Unfinished***** See Issue 4.)
def flxin(Gsize, Tstar, Rstar, Rdist):

	from scipy import integrate, pi, e, log10
	from numpy import sort
	
	#Need to test to see if parameters have been passed.
	 
#*********************************************************************************************************************
# RESTORE, FILENAME='graindata_SuEmiss.idl'	# for grain emiss data 
 # ----- Going to have to import from this file, or equivelent file for python code.			****Unfinished***** See Issue 1.)

 #waves = Gdata.SuWave				# um

 #eps = Gdata.SuEmiss[where(Gdata.SuGsize EQ Gsize),*] #emiss for specific Gsize
#*********************************************************************************************************************

	#########TEMPORARY DO NOT FORGET.... DUDE SRSLY#################	Added temporary values for Waves and EPS for testing. 1.) -JB
	waves = list(sort([10e10,30e10,20e10,50e10,40e10,70e10,60e10,]))
	eps = 1 

#Define Constants.
	c = 2.99792458e8				# m/s 
	k = 1.3806503e-23				# m2 kg s-2 K-1 
	h = 6.626068e-34				# m2 kg / s 
	sigma 	= 5.670400e-8				# W m^-2 K^-4 


#constants appropriate to MKS units. 
	c1 = 2.*pi*h/(c*c)				# 4.6322752e-50 
	c2 = 4.7992370e-11				# h/k (-F)


# Calculate BB curve for Stellar Radiation, given Tstar

	for i, ele in enumerate(waves):
		waves[i] = ele*1e-6                     # um to m 

  	freq = []

	for i, ele in enumerate(waves):
		freq.append(c/waves[i])			# Hz 
    
  	val =  []

	for i, ele in enumerate(freq):
		val.append (freq[i]*c2/Tstar)  		# Val = (freq*h)/(k*Tstar) 
		

	

# mstr = machar(double = (size(val,/type) EQ 5) ) 	# See Issue 2.)      
#good = where( val LT alog(mstr.xmax), Ngood )   	
#if ( Ngood GT 0 ) then  $				 
#flux[ good ] =  C1 *freq[good]^3 / ( exp( val[good])-1. )		


	flux = []
	for i, ele in enumerate (freq): 
		flux.append(c1*freq[i]**3 / (e**val[i]-1))	#J/s/m2/Hz [J/m2]





# Multiply BB (Bnu) times Emissivity 
#_______________________________________________


	for i, ele in enumerate(flux):		# Might be able to "push up" into preceding for-loop. -JB/CK
		flux[i] = ele * eps * freq[i]	# for better integration -FM	 [gives energy -JB]


# integrate BB*e curve over freq (from w)	  
#_______________________________________________

  			
  						 							
	for i, ele in enumerate(freq):		# This was logged in the idl code to avoid accuracy loss due to large intervals between descrete
		freq[i] = log10(ele)		# values. However, python seems to give exact smooth integration via integrate.simps, do we still 							need to use the Log scale? -JB/CK

	Lin = integrate.simps(flux,freq)	# integrate flux over all waves=w, [W/m^2/s] -FM	 MAY need to sort later on. See Issue 3.) -JB


# Calculate the Energy Intersected by grain	
#_______________________________________________

	AU2m = 1.4959787066e11			# meters
  #Rsun = 6.95508d8				# meters

	rdist_m = Rdist * AU2m			# m, radial location
	a_m = Gsize * 1.e-6			# m, grainsize
	Rstar_m = Rstar * AU2m			# m, Rstar
	Ein = ((a_m/(2.*rdist_m))**2) * (4*pi*Rstar_m**2) * Lin


	print('Gsize, rdist: ',Gsize, Rdist)

	print('Ein: ',Ein)			#E absorbed

	Fin = (a_m**2*pi)*(4*pi*Rstar_m**2*sigma*Tstar**4)/(4*pi*rdist_m**2)  #Incident E

	print('Compare w/: ', Fin)

	print( 'Ratio: ', Ein/Fin)


 	return Ein					# Watts
						# solar constant 1400 W/m2
						# sun emitts 3.839d26 Watts

if __name__== "__main__":			# For testing set all parameter values to 1 and print the outcome.
	x = flxin(1,1,1,1)
	print (x)
	
	
 

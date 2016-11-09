function flxout, Gsize, Tout


#
# NAME:
#	FLXOUT()
#
# DATE	December 2009
#
# PURPOSE:
#	To calculate the Energy or flux coming OUT of a grain of size Gsize
#	Use conservation of energy to determine the final grain Temperature
#	since there is no analytical way to solve for T in Planck function.
#
# CALLING SEQUENCE: 
#       Eout = FLXOUT(Gsize,Tout)
#
# INPUT PARAMETERS:
#	GSIZE	Size of grain in um (single value)
#	TOUT	Temperature of grain, K (single value)
#
# OUTPUT PARAMETERS:
#       Eout - Scalar or vector giving the Energy flux (i.e. Area*Fnu)
#               in Watts for a specified grainsize, and temperature.
#
#
#******************************************************************************************************************************************
# THIS CODE IS STILL IN A STRANGE IDL/PYTHON FRANKENSTEIN STATE: DON'T TRY TO MAKE TOO MUCH SENSE OF IT.
#******************************************************************************************************************************************

 On_error,2

  if ( N_elements( Gsize ) NE 1 ) then $
      read,'Enter grain size (um) ', Gsize

  if ( N_elements( Tout ) NE 1 ) then $
      read,'Enter distance of grain to the star (AU) ', rdist



 RESTORE, FILENAME='graindata_SuEmiss.idl'	# for grain emiss data 

 waves = Gdata.SuWave				# um

 eps = Gdata.SuEmiss[where(Gdata.SuGsize EQ Gsize),*] #emissivities for specific Gsize

 	c = 2.99792458d8				# m/s
 	k = 1.3806503d-23				# m2 kg s-2 K-1
 	h = 6.626068d-34				# m2 kg / s


#constants appropriate to MKS units.
 c1 = 2.*!DPI*h/(c*c)
 c2 = h/k


# Calculate BB curve for grain radiation, given grain size and Temp
#_______________________________________________

  flux = waves*0.				# array same size as waves

  w = waves * 1.E-6                             # um to m

  freq = c/w					# Hz
    
  val =  freq*c2/Tout  

  mstr = machar(double = (size(val,/type) EQ 5) )  #Get machine precision      
  good = where( val LT alog(mstr.xmax), Ngood )    #Avoid floating underflow

  if ( Ngood GT 0 ) then  $
      flux[ good ] =  C1 *freq[good]^3 / ( exp( val[good])-1. )		#J/s/m2/Hz




# Multiply BB (Bnu) times Emissivity
#_______________________________________________

  flux = eps * flux				# emiss array and flux array must be the
						# same size.  This is for 1 grain size
						# fed into routine


# integrate BB*e curve over freq (from w)	;
#_______________________________________________

  newflux= flux * freq				# for better integration
  newfreq = alog(freq)				# since int_tabulated integrates linearly

  Lout = int_tabulated(newfreq,newflux,/sort)	#integrate flux over all waves=w, [W/m^2]

  #oplot,freq,flux				; to see area under curve




# Energy radiated by grain over surface area
#_______________________________________________

  a_m = Gsize * 1.d-6				# m, grainsize

  Eout = (4.*!dpi*a_m^2) * Lout

  #print, 't,Lout,Eout',tout,Lout,Eout

  return, Eout					# Watts


  end 

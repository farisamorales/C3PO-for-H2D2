function flxin, Tstar, Rstar, Gsize, rdist


;
; NAME:
;	FLXIN()
;
; DATE	December 2009
;
; PURPOSE:
;	To calculate the Energy or flux going INTO the grain
;	Use conservation of energy to determine the final grain Temperature
;	since there is no analytical way to solve for T in Planck function.
;
; CALLING SEQUENCE: 
;       Ein = FLXIN(Tstar,Rstar,Gsize,rdist)
;
; INPUT PARAMETERS:
;	TSTAR	Temperature of star in degree K
;	RSTAR	Radius of star in AU
;	GSIZE	Size of grain in um (single value)
;	RDIST	Distance from grain to star in AU
;
; OUTPUT PARAMETERS:
;       EIN - Scalar or vector giving the flux (i.e. Area*Fnu)
;               in Watts for a specified grain size, location, and
;		stellar properties.
;
;
;

 On_error,2

  if ( N_elements( Gsize ) NE 1 ) then $
      read,'Enter Temp of the star ', Gsize

  if ( N_elements( Tstar ) NE 1 ) then $
      read,'Enter Temp of the star ', Tstar

  if ( N_elements( Rstar ) NE 1 ) then $
      read,'Enter Radius of the star ', Rstar

  if ( N_elements( rdist ) NE 1 ) then $
      read,'Enter distance of grain to the star ', rdist



 RESTORE, FILENAME='graindata_SuEmiss.idl'	; for grain emiss data 


 waves = Gdata.SuWave				; um


 eps = Gdata.SuEmiss[where(Gdata.SuGsize EQ Gsize),*] ;emiss for specific Gsize

 c = 2.99792458d8				; m/s
 k = 1.3806503d-23				; m2 kg s-2 K-1
 h = 6.626068d-34				; m2 kg / s
sigma 	= 5.670400d-8				; W m^-2 K^-4


;constants appropriate to MKS units.
; c1 = 4.6322752e-50				; 2.*!DPI*h/(c*c)
 c1 = 2.*!DPI*h/(c*c)				; 4.6322752e-50
 c2 = 4.7992370d-11				; h/k


; Calculate BB curve for Stellar Radiation, given Tstar
;_______________________________________________

 flux = waves*0.				; array same size as waves

  w = waves * 1.E-6                             ; um to m

  freq = c/w					; Hz
    
  val =  freq*c2/Tstar  			; print, val

  mstr = machar(double = (size(val,/type) EQ 5) )  ;Get machine precision     ((Check for Underflow)) 
  good = where( val LT alog(mstr.xmax), Ngood )    ;Avoid floating underflow

  if ( Ngood GT 0 ) then  $
      flux[ good ] =  C1 *freq[good]^3 / ( exp( val[good])-1. )		;J/s/m2/Hz

;      flux =  C1 *freq^3 / ( exp(val)-1. )		;J/s/m2/Hz



; Multiply BB (Bnu) times Emissivity
;_______________________________________________

  flux = eps * flux				; emiss array and flux array must be the
						; same size.  This is for 1 grain size
						; fed into routine


; integrate BB*e curve over freq (from w)	;
;_______________________________________________

  newflux= flux * freq				; for better integration
  newfreq = alog(freq)				; since int_tabulated integrates linearly

  Lin = int_tabulated(newfreq,newflux,/sort)	;integrate flux over all waves=w, [W/m^2]


; Calculate the Energy Intersected by grain	;
;_______________________________________________

  AU2m = 1.4959787066d11			; meters
  ;Rsun = 6.95508d8				; meters

  rdist_m = rdist * AU2m			; m, radial location
  a_m = Gsize * 1.d-6				; m, grainsize
  Rstar_m = Rstar * AU2m			; m, Rstar

  Ein = ((a_m/(2.*rdist_m))^2) * (4*!dpi*Rstar_m^2) * Lin

;  E_SolarConst = ((1/(2.*rdist_m))^2)/!dpi * (4*!dpi*Rstar_m^2) * Lin
;  print,'Solar Constant (W/m2): ',E_SolarConst

;  E_Lsun = (4*!dpi*Rstar_m^2) *  Lin
;  print,'L-Sun (W): ',E_Lsun


						;Check: compare with Lstar/4piR^2=4

  print,''

  print, 'Gsize, rdist: ',Gsize, rdist

  print, 'Ein: ',Ein						;E absorbed

  Fin = (a_m^2*!dpi)*(4*!dpi*Rstar_m^2*sigma*Tstar^4)/(4*!dpi*rdist_m^2)  ;Incident E

  print, 'Compare w/: ', Fin

  print, 'Ratio: ', Ein/Fin


  return, Ein					; Watts
						; solar constant 1400 W/m2
						; sun emitts 3.839d26 Watts
  end 

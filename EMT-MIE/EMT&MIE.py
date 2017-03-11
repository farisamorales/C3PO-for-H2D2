'''
	Combines Real and Imaginary optical constants using Effective Medium Theory for SCIENCE
	Then calculates emissivites based on the optical constants derived- also for SCIENCE
	
	Author: Dr. Farisa Morales (original IDL code) transcribed and modified by Cody King
	
	TODO:
		-Revert to declaring appropriate values imaginary in the EMT function to degenerate imaginator function
		-Change to usage of numpy readers(?)
		-Generalize variable names to accomodate other substances?
		-Make stylistic changes as per suggestion
'''

from numpy import interp, real, imag, pi
#interp(x-coord of desired interpolates, x-coord of input data, y-coord of input data)
#Returns a numpy array object; this can be quickly converted to a standard python list
#I don't know if it's necessary (probably not) but I do it anyway just in case
from kmod import whiteSpaceParser, getStart, getCol, listifier, columnizer
from cmath import sqrt
from math import log10
from os import times

TIME1 = times()[0]

print('Initializing...')

#--------------------
#- Global Variables -
#--------------------
#For whether or not we will be making our own polluted matrix prior to re-inclusion; 
#if true we do not make our own
preDMat = True

#For whether or not we will be running MIE theory on the optical constants calculated
#here. The value of this constant should be self-explanitory to what it means.
RUNMIE = True

#For whether or not we will be writing new lists to memory. If true, checks if we want to
#write the results of interpolations, EMT, and if running MIE, the results of MIE theory
WRITE = True
if WRITE:
	WRITE_INT = False
	WRITE_EMT = False
	if RUNMIE:
		WRITE_MIE = True

#For whether or not we will be generating graphs of the results of the intpolations,
#EMT, and if running MIE, the results of MIE theory.
GRAPH = True
if GRAPH:
	GRAPH_EMT = False
	if RUNMIE:
		GRAPH_MIE = True
	import matplotlib.pyplot as plt

print('Orbital Variables:\n\tpreDMat-' + str(preDMat) + '\n\tRUNMIE-' + str(RUNMIE) + '\n\tWRITE-' + str(WRITE) + '\n\tGRAPH-' + str(GRAPH) + '\nSuborbital Variables:')
if WRITE:
	print('\tWRITE_INT- ' + str(WRITE_INT) + '\n\tWRITE_EMT- ' + str(WRITE_EMT))
	if RUNMIE:
		print('\tWRITE_MIE- ' + str(WRITE_MIE))
if GRAPH:
	print('\tGRAPH_EMT- ' + str(GRAPH_EMT))
	if RUNMIE:
		print('\tGRAPH_MIE- ' + str(GRAPH_MIE))

#---------------------------------------------
#- File directories for the necessary tables -
#---------------------------------------------

#File directory for written tables
directory = r'C:\Users\Cody\Desktop\School Documents\Physics\Independent Study\DielectricEffect'

print('File Directory:\n' + directory)

if not preDMat:
	#Amorphous Carbon Optical Constants
	dir_AC = directory + r'\AmorphCarbOptConst.dat'
	
	#Water Optical Constants
	dir_WA = directory + r'\waterOpticalConstants.dat'

else:
   #Dirty Ice Optical Constants
	dir_DI = directory + r'\DirtyIceOptConst.dat'
	
#AstroSil Optical Constants
dir_AS = directory+ r'\AstroSilOptConst.dat'

#Desired Wavelengths
dir_W = directory + r'\wav_set.dat'

#Desired Grain Sizes; current table is currently unreadable by my modules so it needs to manually
#be modified; sorry!
dir_G = directory + r'\grain_size_grid_suvsil.list'

#------------------------------------------------------
#- Short functions which simplify (maybe) some things -
#------------------------------------------------------
#----------------------------------------------------------------
#Short function to trim lists based on designated mins and maxes
def trimmer(LIST, mini, maxi):
	tempList = []
	for ele in LIST:
		if ele >= mini and ele <= maxi:
			tempList.append(ele)
	return tempList
#End trimming function
#-------------------------------------------------------------
#Short function to write to a file data from a dictionary.
#Writes to a csv file the keys as a header followed by the values
#associated with each.
#Python doesn't care if we want the string read-in as a raw string
#if it sees a \' at the end of a string it complains about end-of-line errors!
def writerD(DICT, FILENAME):
	tempFile = open(directory + '\\' + FILENAME, 'w')
	for ele in DICT.keys():
		if ele is not 'NAME':
			tempFile.write(str(ele) + ',')
	tempFile.write('\n')
	row = 0
	while row < len(DICT[DICT.keys()[0]]):
		for key in DICT.keys():
			if key is not 'NAME':
				tempFile.write(str((DICT[key][row])) + ',')
		tempFile.write('\n')
		row += 1
	tempFile.close()
#End short writing function
#--------------------------------------------------------
#Short writing function specifically for the emissivities.
#Writes to a csv file a header of wavelengths followed by
#rows of emissivities with the first value in each row being
#the associated grain radius.
#Done as dictionaries are (seemingly) randomly ordered when it
#comes to numerical keys
def writerE(HEADER, ARRAYARRAY, FILENAME):
	tempFile = open(directory + '\\' + FILENAME, 'w')
	tempFile.write(str(None) + ',')	#Include header explaining table and content and generating parameters?
	for ele in HEADER:
		tempFile.write(str(ele) + ',')
	tempFile.write('\n')
	for ARRAY in ARRAYARRAY:
		for ele in ARRAY:
			tempFile.write(str(ele) + ',')
		tempFile.write('\n')
	tempFile.close()
#End short writing function
#------------------------------------------------------------
#Short function to change values in a list to complex numbers
def imaginator(LIST, demaginate=False):
	if demaginate:
		for i, ele in enumerate(LIST):
			LIST[i] = float(imag(ele))
	else:
		for i, ele in enumerate(LIST):
			LIST[i] = 1j * ele
	return LIST
#No, it wasn't necessary to write a function for this too but I'm doing it anyway
#------------------------------------------------------------------------------------------
#Short function to calculate the average optical constants based on Effective Medium Theory
#NM: Refractive Index of Matrix, KM: Extinction Coefficient of Matrix
#NI: Refractive Index of Inclusion, KI: Extinction Coefficient of Inclusion
#V: Ratio of pollutant volume to total volume
def EMT(NM, KM, NI, KI, V=.5):
	M, I = (NM + KM) ** 2, (NI + KI) ** 2
	F = (I - M) / (I + 2.0 * M)
	AMOC = sqrt(M * (1 + ((3.0 * V * F) / (1 - V * F))))
	return float(real(AMOC)), float(imag(AMOC))
#End REALLY short EMT function

#-----------------------------------------------
#- Writing data from tables in files to memory -
#-----------------------------------------------

print('Writing required files to memory...')

if not preDMat:
	#--------------------
	#- Amorphous Carbon -
	#--------------------
	#Open the Amorphous Carbon table in binary mode
	AC_Tab = open(dir_AC, 'rb')
	
	#Start the reader at the appropriate postion
	AC_Tab.seek(getStart(AC_Tab))

	#Do the thing! Seperate the values!
	#Assign to mList the list of lists returned from Columnizer which took a single
	#list returned by listifier (as well as the number of columns in the table
	#found using the getCol module) which itself took a tuple returned by
	#whiteSpaceParser which consisted of the modified contents of the read file
	#in addition to the index location of whitespaces in that file.
	#For more information, check kmod.py
	mList = columnizer(listifier(whiteSpaceParser(AC_Tab.read())), getCol(AC_Tab))
	AmorphCarb = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[1], 'Extinction Coefficient':mList[2], 'NAME':'Amorphous Carbon'}

	#Close the file and delete variables no longer needed
	AC_Tab.close()
	del(AC_Tab, mList, dir_AC)
	
	#-------------
	#- Water Ice -
	#-------------
	#Similar proceedure for Water
	WA_Tab = open(dir_WA, 'rb')
	
	WA_Tab.seek(getStart(WA_Tab))
	
	mList = columnizer(listifier(whiteSpaceParser(WA_Tab.read())), getCol(WA_Tab))
	
	Water = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[1], 'Extinction Coefficient':mList[2], 'NAME':'Water Ice'}

	WA_Tab.close()

	del(WA_Tab, mList, dir_WA)

	#----------------------
	#- Itermediary Matrix -
	#----------------------
	#For graphing purposes
	IMM = {'Wavelength(Microns)':[], 'Refractive Index':[], 'Extinction Coefficient':[], 'NAME':'IMM'}
	
	
else:
	#-------------
	#- Dirty Ice -
	#-------------
	#Similar proceedure for Dirty Ice
	DI_Tab = open(dir_DI, 'rb')
	
	DI_Tab.seek(getStart(DI_Tab))

	mList = columnizer(listifier(whiteSpaceParser(DI_Tab.read())), getCol(DI_Tab))

	DirtyIce = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[1], 'Extinction Coefficient':mList[2], 'NAME':'Dirty Ice'}

	DI_Tab.close()

	del(DI_Tab, mList, dir_DI)
	
#------------
#- AstroSil -
#------------
#Similar proceedure for AstroSil
AS_Tab = open(dir_AS, 'rb')

AS_Tab.seek(getStart(AS_Tab))

mList = columnizer(listifier(whiteSpaceParser(AS_Tab.read())), getCol(AS_Tab))

#Values must be put in reverse order for ease of interpolation later
#Reversing can't be done in a dictionary assignment because it returns
#a None type for the value for some reason
mList[0].reverse()
mList[3].reverse()
mList[4].reverse()
AstroSil = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[3], 'Extinction Coefficient':mList[4], 'NAME':'Astro Silicates'}

#Apparently the author of the AstroSil table created the table with the real parts decremented by 1?
#Some of the values in the list now incremented by 1 show floating point arithmetic relics
#I'm pretty confident that's just a visual problem
for i, ele in enumerate(AstroSil['Refractive Index']):
	AstroSil['Refractive Index'][i] = ele + 1.0

AS_Tab.close()

del(AS_Tab, mList, dir_AS, i, ele)

#-----------------
#- Desired Waves -
#-----------------
#Similar proceedure for the desired Wavelengths
W_Tab = open(dir_W, 'rb')

W_Tab.seek(getStart(W_Tab))

mList = columnizer(listifier(whiteSpaceParser(W_Tab.read())), getCol(W_Tab))

Waves = {'Wavelength(Microns)':mList[0]}

W_Tab.close()

del(W_Tab, mList, dir_W)

#-----------------------
#- Desired Grain Sizes -
#-----------------------
#A litte different proceedure for the desired Grain Sizes
G_Tab = open(dir_G, 'rb')

G_Tab.seek(getStart(G_Tab))

mList = columnizer(listifier(whiteSpaceParser(G_Tab.read())), getCol(G_Tab))

Waves.update({'SIZE':mList[0]})

G_Tab.close()

del(G_Tab, mList, dir_G)

#-----------------------------------------------------------
#- Trimming Wave lists for useful Wavelengths; the prequel - 
#-----------------------------------------------------------
#New Wavelengths based on Primary Matrix Material; determines Wavelength for all others.
#As this is needed for interpolation it must be done before the interpolation
#but assignment to other dictionaries must be done after.

print('Trimming Wave list...')

if not preDMat:
	Waves.update({'Wavelength(Microns)':trimmer(Waves['Wavelength(Microns)'], min(Water['Wavelength(Microns)']), max(Water['Wavelength(Microns)']))})

else:
	Waves.update({'Wavelength(Microns)':trimmer(Waves['Wavelength(Microns)'], min(DirtyIce['Wavelength(Microns)']), max(DirtyIce['Wavelength(Microns)']))})

#---------------------------------------------------------------------------
#- Interpolate using SciPy's (NumPy's?) interp function for missing values -
#---------------------------------------------------------------------------
#Update the dictionaries with these interpolations, changing complex coefficients
#into complex type. For some reason I get erroneous optical constants from EMT if
#I change these into complex numbers piecewise during EMT, I still can't figure
#out why

print('Interpolating for missing values and updating the appropriate dictionaries...')

if not preDMat:
	#Update Amorphous Carbon Dict with interpolates for reals
	AmorphCarb.update({'Refractive Index':list(interp(Waves['Wavelength(Microns)'], AmorphCarb['Wavelength(Microns)'], AmorphCarb['Refractive Index']))})

	#Update Amorphous Carbon Dict with interpolates for complex
	AmorphCarb.update({'Extinction Coefficient':imaginator(list(interp(Waves['Wavelength(Microns)'], AmorphCarb['Wavelength(Microns)'], AmorphCarb['Extinction Coefficient'])))})

	#Update Water Dict with interpolates for reals
	Water.update({'Refractive Index':list(interp(Waves['Wavelength(Microns)'], Water['Wavelength(Microns)'], Water['Refractive Index']))})

	#Update Water Dict with interpolates for complex
	Water.update({'Extinction Coefficient':imaginator(list(interp(Waves['Wavelength(Microns)'], Water['Wavelength(Microns)'], Water['Extinction Coefficient'])))})

else:
	#Update DirtyIce Dict with interpolates for reals
	DirtyIce.update({'Refractive Index':list(interp(Waves['Wavelength(Microns)'], DirtyIce['Wavelength(Microns)'], DirtyIce['Refractive Index']))})

	#Update DirtyIce Dict with interpolates for complex
	DirtyIce.update({'Extinction Coefficient':imaginator(list(interp(Waves['Wavelength(Microns)'], DirtyIce['Wavelength(Microns)'], DirtyIce['Extinction Coefficient'])))})


#Update AstroSil Dict with interpolates for reals
AstroSil.update({'Refractive Index':list(interp(Waves['Wavelength(Microns)'], AstroSil['Wavelength(Microns)'], AstroSil['Refractive Index']))})

#Update AstroSil Dict with interpolates for complex
AstroSil.update({'Extinction Coefficient':imaginator(list(interp(Waves['Wavelength(Microns)'], AstroSil['Wavelength(Microns)'], AstroSil['Extinction Coefficient'])))})

#----------------------------------------------------------
#- Trimming Wave lists for useful Wavelengths; the sequel - 
#----------------------------------------------------------

print('Doin\' some more trimming...')

if not preDMat:
	#New Wavelenghts for Amorphous Carbon
	AmorphCarb.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})

	#New Wavelenghts for Water
	Water.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})
	
	#New Wavelenghts for the intermediary matrix
	IMM.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})

else:
	#New Wavelenghts for DirtyIce
	DirtyIce.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})
	
#New Wavelenghts for AstroSil
AstroSil.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})

#--------------------------------------------------------------
#- Computing Optical Constants for Inclusion-Matrix Particles -
#--------------------------------------------------------------
#Finally doing part of what this program is named after! Yay!
#Initialize the Optical Constants dictionary

print('Implimenting EMT to compute Optical Constants...')

IMPOptConst = {'Wavelength(Microns)':Waves['Wavelength(Microns)'], 'Refractive Index':[], 'Extinction Coefficient':[], 'NAME':'Inclusion Matrix Particles'}

#Iterate through designated reals and utilize EMT to generate new
#optical constants

if not preDMat:
	for i, ele in enumerate(Water['Refractive Index']):
	
		#Utilize EMT for Amorphous Carbon and Water
		RI, EC = EMT(ele, Water['Extinction Coefficient'][i], AmorphCarb['Refractive Index'][i], AmorphCarb['Extinction Coefficient'][i])

		#Append the appropriate values to the IMM dictionary
		IMM['Refractive Index'].append(RI)
		if EC < 0:
			EC = 1.0e-11
		IMM['Extinction Coefficient'].append(EC)

		#Utilize the results from the last utilization of EMT for this run of EMT
		#now with AstroSil
		RI, EC = EMT(IMM['Refractive Index'][i], IMM['Extinction Coefficient'][i], AstroSil['Refractive Index'][i], AstroSil['Extinction Coefficient'][i])
	
		#Append the appropriate values to the IMP dictionary
		IMPOptConst['Refractive Index'].append(RI)
		if EC < 0:
			EC = 1.0e-11
		IMPOptConst['Extinction Coefficient'].append(EC)
			
else:
	for i, ele in enumerate(DirtyIce['Refractive Index']):
		
		#Utilize EMT with the Predetermined Matrix and Inclusion
		RI, EC = EMT(ele, DirtyIce['Extinction Coefficient'][i], AstroSil['Refractive Index'][i], AstroSil['Extinction Coefficient'][i])
		
		#Append the appropriate values to the IMP dictionary
		IMPOptConst['Refractive Index'].append(RI)
		if EC < 0:
			EC = 1.0e-11
		IMPOptConst['Extinction Coefficient'].append(EC) 

del(i, ele, RI, EC)

if WRITE:

	if WRITE_INT:
		print('Writing Interpolate Data to disk...')
		if not preDMat:
			writerD(AmorphCarb, 'AmorphCarbInterp2.csv')
			writerD(Water, 'WaterInterp2.csv')	
		else:
			writerD(DirtyIce, 'DirtyIceInterp2.csv')
			writerD(AstroSil, 'AstroSilInterp2.csv')
	if WRITE_EMT:
		print('Writing EMT Data to disk...')
		writerD(IMPOptConst, 'IMPOptConst2.csv')
	
if GRAPH and GRAPH_EMT:
#------------------------------------------------
#- Modify complex coefficients to be real again -
#------------------------------------------------
#Must be done else MatPlotLib complains when attempting to plot the values

	print('De-imaginating some values for graphing purposes...')

	if not preDMat:
		imaginator(AmorphCarb['Extinction Coefficient'], True)
		imaginator(Water['Extinction Coefficient'], True)

	else:
		imaginator(DirtyIce['Extinction Coefficient'], True)
	
	imaginator(AstroSil['Extinction Coefficient'], True)
	
#--------------------------------------------------
#- Plot appropriate values if designated to do so -
#--------------------------------------------------
	def flotEMT(*flargs):
		plt.xlabel('Wavelength(Microns)')
		plt.ylabel('Exitinction Coefficient')
		plt.title('Extinction Coefficient Vs. Wavelength')
		plt.xscale('log')
		plt.yscale('log')
		plt.grid(True)
		for Exti in flargs:
			plt.plot(Exti['Wavelength(Microns)'], Exti['Extinction Coefficient'], label = Exti['NAME'])
		plt.legend()
		plt.show()
		
		plt.xlabel('Wavelength(Microns)')
		plt.ylabel('Refractive Index')
		plt.title('Refractive Index Vs. Wavelength')
		plt.xscale('log')
		plt.grid(True)
		for Refracti in flargs:
			plt.plot(Refracti['Wavelength(Microns)'], Refracti['Refractive Index'], label = Refracti['NAME'])
		plt.legend(loc = 2)
		plt.show()

	if not preDMat:		  
		flotEMT(Water, AmorphCarb, IMM, AstroSil, IMPOptConst)
	else:
		flotEMT(DirtyIce, AstroSil, IMPOptConst)

if RUNMIE:
#--------------
#- MIE Theory -
#--------------
#Each refractive index value is associated with a particular Wavelength.
#Inside an iteration of the values of the refractive index, we iterate
#through the grain sizes. We call bhmie with the size parameter associated
#with that particular grain size with that particular Wavelength, along
#with the sum of the refractive index and extinction coefficient associated
#with that Wavelength. The resulting Emiss list is a list of lists; each
#sublist is the emissivity of the grains of the various sizes associated
#with the Wavelength of incident light.

#Seems to give me something new everytime its run without changing anything!
#Delete old dictionaries to free up memory for the Emiss array

	print('Utilizing MIE Theory to calculate Emissivity of the IMPs.\nDepending on the number of Wavelengths and Grain Sizes,\nthis can take qite a bit of time... [>1hr]')

	if not preDMat:
		del(AmorphCarb, Water, IMM)
	else:
		del(DirtyIce)
	del(AstroSil)
	
	from bhmie_herbert_kaiser_july2012_GEOFF_EDITION import bhmie
	Emiss = []
	Time0 = times()[0]
	sizey = len(Waves['SIZE'])
	for iter, radius in enumerate(Waves['SIZE']):
		print('Emissivites for grains of radius: ' + str(radius) + ' microns. ' + str(iter + 1) + ' of ' + str(sizey) + '.\nTime for last radius: ' + str(times()[0] - Time0) + ' seconds')
		Time0 = times()[0]
		Emiss.append([radius])
		for i, lamda in enumerate(Waves['Wavelength(Microns)']):
			RI = IMPOptConst['Refractive Index'][i]
			if type(IMPOptConst['Extinction Coefficient'][i]) is not type(1j):
				EC = IMPOptConst['Extinction Coefficient'][i]*1j
			else:
				EC = IMPOptConst['Extinction Coefficient'][i]
			qabs = bhmie((2*pi*radius)/lamda, RI + EC, 2)
			if (qabs <= 0):
				qabs = 1.0e-11
			Emiss[iter].append(qabs)

if WRITE and WRITE_MIE:
	print('Writing MIE Data to disk...')
	writerE(Waves['Wavelength(Microns)'], Emiss, 'Emissivities.csv')

if GRAPH and GRAPH_MIE:
	def flotMIE(EMISS):
		plt.xlabel('Wavelength(Microns)')
		plt.ylabel('Emissivity')
		plt.title('Emissivity Vs. Wavelength')
		plt.xscale('log')
		plt.yscale('log')
		plt.grid(True)
		for ARRAY in EMISS:
			if (log10(ARRAY[0]) == int(log10(ARRAY[0]))):
				plt.plot(Waves['Wavelength(Microns)'], ARRAY[1:], label = ARRAY[0])
		plt.legend(loc = 0)
		plt.show()

	flotMIE(Emiss)

print('Run Time: ' + str(times()[0] - TIME1) + ' seconds')
raw_input('Press ENTER to exit')

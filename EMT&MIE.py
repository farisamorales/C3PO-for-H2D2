'''
    Combines Real and Imaginary optical constants using Effective Medium Theory for SCIENCE
    Then calculates emissivites based on the optical constants derived- also for SCIENCE
    
    Author: Dr. Farisa Morales (original IDL code) transcribed and modified by Cody King
    
    TODO:
    	-Deimpliment the 'as's from the importations
	-Remove superflous & misleading interper function
	-Remove the unecessary dictUpdate function
	-Rename Kate's waves to something more general
    	-Impliment MIE theory
	-
'''

#import matplotlib.pyplot as plt
from numpy import interp, real as rl, imag as im
from kmod import whiteSpaceParser, getStart as GS, getCol as GC, listifier as li, columnizer as clm
from cmath import sqrt as sqrt

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
#---------------------------------------------------
#Short function to update dictionaries with new data because I'm /this/ lazy
def dictUpdate(dictName, LIST, KEYNAME):
    dictName.update({KEYNAME:LIST})
#End dictionary updater
#-------------------------------------------------------------------------------
#Short function to ease doing /this/ again and again and again and again and...
#interp(x-coord of desired interpolates, x-coord of input data, y-coord of input data)
#returns a numpy array object; should be converted to standard python list obj
def interper(x, y):
    tempList = list(interp(Kates['WAV'], x, y))
    return tempList
#End short function for ease of not doing whatever /this/ was
#-------------------------------------------------------------
#Short function to write to a file data from a dicitonary
#Python doesn't care if we want the string read-in as a raw string
#if it sees a \' at the end of a string it complains about end-of-line errors!
def writer(DICT, FILENAME):
    tempFile = open(directory + '\\' + FILENAME, 'w')
    for ele in DICT.keys():
        tempFile.write(str(ele) + ' ')
    tempFile.write('\n')
    row = 0
    while row < len(DICT[DICT.keys()[0]]):
        for key in DICT.keys():
            tempFile.write(str((DICT[key][row])) + ' ')
        tempFile.write('\n')
        row += 1
    tempFile.close()
#End short writing function
#------------------------------------------------------------
#Short function to change values in a list to complex numbers
def imaginator(LIST, demaginate=False):
    if demaginate:
        for i, ele in enumerate(LIST):
            LIST[i] = float(im(ele))
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
    M, I = NM + KM, NI + KI
    F = (I ** 2 - M ** 2) / (I + 2.0 * M)
    AMOC = sqrt(M * (1 + ((3.0 * V * F) / (1 - V * F))))
    return float(rl(AMOC)), float(im(AMOC))
#End short EMT function
#---------------------------------------------
#- File directories for the necessary tables -
#---------------------------------------------

#File directory for wirtten tables
directory = r'C:\Users\Cody\Desktop\School Documents\Physics\Independent Study\DielectricEffect'

#Amorphous Carbon Optical Constants
dir_AC = directory + r'\AmorphCarbOptConst.dat'

#AstroSil Optical Constants
dir_AS = directory+ r'\AstroSilOptConst.dat'

#Dirty Ice Optical Constants
dir_DI = directory + r'\DirtyIceOptConst.dat'

#Water Optical Constants
dir_WA = directory + r'\waterOpticalConstants.dat'

#Kate's Waves
dir_KW = directory + r'\wav_set.dat'

#Kate's Grain Sizes; table is currently unreadable by my modules so it needs to manually
#be modified; sorry!
dir_KG = directory + r'\grain_size_grid_suvsil.list'


#-----------------------------------------------
#- Writing data from tables in files to memory -
#-----------------------------------------------

#--------------------
#- Amorphous Carbon -
#--------------------

#Open the Amorphous Carbon table in binary mode
AC_Tab = open(dir_AC, 'rb')

#Start the reader at the appropriate postion
AC_Tab.seek(GS(AC_Tab))

#Do the thing! Seperate the values!
mList = clm(li(wsp(AC_Tab.read())), GC(AC_Tab))
AmorphCarb = {'WAV':mList[0], 'N':mList[1], 'K':mList[2]}

#Close the file and delete variables no longer needed
AC_Tab.close()
del(AC_Tab, mList, dir_AC)

#------------
#- AstroSil -
#------------
#Similar proceedure for AstroSil
AS_Tab = open(dir_AS, 'rb')

AS_Tab.seek(GS(AS_Tab))

mList = clm(li(wsp(AS_Tab.read())), GC(AS_Tab))

#Values must be put in reverse order for ease of interpolation later
#Reversing can't be done in a dictionary assignment because it returns
#a None type for the value for some reason
mList[0].reverse()
mList[3].reverse()
mList[4].reverse()
AstroSil = {'WAV':mList[0], 'N':mList[3], 'K':mList[4]}

#Apparently the author of the AstroSil table created the table with the real parts decremented by 1?
#Some of the values in the list now incremented by 1 show floating point arithmetic relics
#I'm pretty confident that's just a visual problem
for i, ele in enumerate(AstroSil['N']):
		AstroSil['N'][i] = ele + 1.0

AS_Tab.close()

del(AS_Tab, mList, dir_AS, i, ele)

#-------------
#- Dirty Ice -
#-------------
#Similar proceedure for Dirty Ice
DI_Tab = open(dir_DI, 'rb')

DI_Tab.seek(GS(DI_Tab))

mList = clm(li(wsp(DI_Tab.read())), GC(DI_Tab))

DirtyIce = {'WAV':mList[0], 'N':mList[1], 'K':mList[2]}

DI_Tab.close()

del(DI_Tab, mList, dir_DI)

#---------
#- Water -
#---------
#Similar proceedure for Water
WA_Tab = open(dir_WA, 'rb')

WA_Tab.seek(GS(WA_Tab))

mList = clm(li(wsp(WA_Tab.read())), GC(WA_Tab))

Water = {'WAV':mList[0], 'N':mList[1], 'K':mList[2]}

WA_Tab.close()

del(WA_Tab, mList, dir_WA)

#----------------
#- Kate's Waves -
#----------------
#Similar proceedure for Kate's Waves
KW_Tab = open(dir_KW, 'rb')

KW_Tab.seek(GS(KW_Tab))

mList = clm(li(wsp(KW_Tab.read())), GC(KW_Tab))

Kates = {'WAV':mList[0]}

KW_Tab.close()

del(KW_Tab, mList, dir_KW)

#----------------------
#- Kate's Grain Sizes -
#----------------------
#A litte different proceedure for Kate's Grain Sizes
KG_Tab = open(dir_KG, 'rb')

KG_Tab.seek(GS(KG_Tab))

mList = clm(li(wsp(KG_Tab.read())), GC(KG_Tab))

dictUpdate(Kates, mList[0], 'SIZE')

KG_Tab.close()

del(KG_Tab, mList, dir_KG)

#-----------------------------------------------------------
#- Trimming Wave lists for useful wavelengths; the prequel - 
#-----------------------------------------------------------
#New wavelengths for Kates based on DirtyIce; determines wavelength for all others.
#As this is needed for interpolation it must be done before the interpolation
#but assignment to other dictionaries must be done after
dictUpdate(Kates, trimmer(Kates['WAV'], min(DirtyIce['WAV']), max(DirtyIce['WAV'])), 'WAV')

#---------------------------------------------------------------------------
#- Interpolate using SciPy's (NumPy's?) interp function for missing values -
#---------------------------------------------------------------------------
#Update the dictionaries with these interpolations, changing complex coefficients
#into complex type. For some reason I get erroneous optical constants from EMT if
#I change these into complex numbers piecewise during EMT, I still can't figure
#out why

#Update Amorphous Carbon Dict with interpolates for reals
dictUpdate(AmorphCarb, interper(AmorphCarb['WAV'], AmorphCarb['N']), 'N')

#Update Amorphous Carbon Dict with interpolates for complex
dictUpdate(AmorphCarb, imaginator(interper(AmorphCarb['WAV'], AmorphCarb['K'])), 'K')

#Update AstroSil Dict with interpolates for reals
dictUpdate(AstroSil, interper(AstroSil['WAV'], AstroSil['N']), 'N')

#Update AstroSil Dict with interpolates for complex
dictUpdate(AstroSil, imaginator(interper(AstroSil['WAV'], AstroSil['K'])), 'K')

#Update DirtyIce Dict with interpolates for reals
dictUpdate(DirtyIce, interper(DirtyIce['WAV'], DirtyIce['N']), 'N')

#Update DirtyIce Dict with interpolates for complex
#dictUpdate(DirtyIce, imaginator(interper(DirtyIce['WAV'], DirtyIce['K'])), 'K')
dictUpdate(DirtyIce, interper(DirtyIce['WAV'], DirtyIce['K']), 'K')

#Update Water Dict with interpolates for reals
dictUpdate(Water, interper(Water['WAV'], Water['N']), 'N')

#Update Water Dict with interpolates for complex
dictUpdate(Water, imaginator(interper(Water['WAV'], Water['K'])), 'K')

#----------------------------------------------------------
#- Trimming Wave lists for useful wavelengths; the sequel - 
#----------------------------------------------------------
#New wavelengths for Amorphous Carbon
dictUpdate(AmorphCarb, Kates['WAV'], 'WAV')

#New wavelengths for AstroSil
dictUpdate(AstroSil, Kates['WAV'], 'WAV')

#New wavelengths for Water
dictUpdate(Water, Kates['WAV'], 'WAV')

#New wavelengths for DirtyIce
dictUpdate(DirtyIce, Kates['WAV'], 'WAV')

#--------------------------------------------------------------
#- Computing Optical Constants for Inclusion-Matrix Particles -
#--------------------------------------------------------------
#Finally doing part of what this program is named after! Yay!
#Initialize the Optical Constants dictionary
IMPOptConst = {'WAV':Kates['WAV'], 'N':[], 'K':[]}

#Iterate through the Water reals and utilize EMT to generate new
#optical constants
for i, ele in enumerate(Water['N']):
    
    #Utilize EMT for Amorphous Carbon and Water
    N, K = EMT(ele, Water['K'][i], AmorphCarb['N'][i], AmorphCarb['K'][i])
    
    #Utilize the results from the last utilization of EMT for this run of EMT
    #now with AstroSil
    N, K = EMT(N, K, AstroSil['N'][i], AstroSil['K'][i])
    
    #Append the appropriate values to the IMP dictionary
    IMPOptConst['N'].append(N)
    if K < 0:
        K = 1e-11
    IMPOptConst['K'].append(K)

del(i, ele, N, K)

#------------------------------------------------
#- Modify complex coefficients to be real again -
#------------------------------------------------
#Must be done else MatPlotLib complains when attempting to plot the values
imaginator(AmorphCarb['K'], True)

imaginator(AstroSil['K'], True)

#imaginator(DirtyIce['K'], True)

imaginator(Water['K'], True)

'''
def flot(DICT, y):
    if y == 'K':
        plt.ylabel('Exitinction Coefficient')
        plt.title('Extinction Coefficient Vs. Wavelength')
        plt.yscale('log')
    elif y == 'N':
        plt.ylabel('Refractive Index')
        plt.title('Refractive Index Vs. Wavelength')
    plt.plot(DICT['WAV'], DICT[y])
    plt.xscale('log')
    plt.xlabel('Wavelength (microns)')    
    plt.grid(True)
    plt.show()

flot(Water, 'N')
flot(IMPOptConst, 'N')
flot(AstroSil, 'N')
flot(AmorphCarb, 'N')
flot(DirtyIce, 'N')
'''

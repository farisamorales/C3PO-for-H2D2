from numpy import interp, pi
from kmod import whiteSpaceParser, getStart, getCol, listifier, columnizer
from os import times

directory = r'C:\Users\Cody\Desktop\School Documents\Physics\Independent Study\DielectricEffect'
dir_DI = directory + r'\DirtyIceOptConst.dat'
dir_AS = directory+ r'\AstroSilOptConst.dat'
dir_W = directory + r'\wav_set.dat'
dir_G = directory + r'\grain_size_grid_suvsil.list'

def trimmer(LIST, mini, maxi):
	tempList = []
	for ele in LIST:
		if ele >= mini and ele <= maxi:
			tempList.append(ele)
	return tempList

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

DI_Tab = open(dir_DI, 'rb')
DI_Tab.seek(getStart(DI_Tab))
mList = columnizer(listifier(whiteSpaceParser(DI_Tab.read())), getCol(DI_Tab))
DirtyIce = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[1], 'Extinction Coefficient':mList[2], 'NAME':'Dirty Ice'}
DI_Tab.close()
del(DI_Tab, mList, dir_DI)

AS_Tab = open(dir_AS, 'rb')
AS_Tab.seek(getStart(AS_Tab))
mList = columnizer(listifier(whiteSpaceParser(AS_Tab.read())), getCol(AS_Tab))
mList[0].reverse()
mList[3].reverse()
mList[4].reverse()
AstroSil = {'Wavelength(Microns)':mList[0], 'Refractive Index':mList[3], 'Extinction Coefficient':mList[4], 'NAME':'Astro Silicates'}
for i, ele in enumerate(AstroSil['Refractive Index']):
	AstroSil['Refractive Index'][i] = ele + 1.0
AS_Tab.close()
del(AS_Tab, mList, dir_AS, i, ele)
W_Tab = open(dir_W, 'rb')
W_Tab.seek(getStart(W_Tab))
mList = columnizer(listifier(whiteSpaceParser(W_Tab.read())), getCol(W_Tab))
Waves = {'Wavelength(Microns)':mList[0]}
W_Tab.close()
del(W_Tab, mList, dir_W)

G_Tab = open(dir_G, 'rb')
G_Tab.seek(getStart(G_Tab))
mList = columnizer(listifier(whiteSpaceParser(G_Tab.read())), getCol(G_Tab))
Waves.update({'SIZE':mList[0]})
G_Tab.close()
del(G_Tab, mList, dir_G)

Waves.update({'Wavelength(Microns)':trimmer(Waves['Wavelength(Microns)'], min(DirtyIce['Wavelength(Microns)']), max(DirtyIce['Wavelength(Microns)']))})

AstroSil.update({'Refractive Index':list(interp(Waves['Wavelength(Microns)'], AstroSil['Wavelength(Microns)'], AstroSil['Refractive Index']))})
AstroSil.update({'Extinction Coefficient':list(interp(Waves['Wavelength(Microns)'], AstroSil['Wavelength(Microns)'], AstroSil['Extinction Coefficient']))})
AstroSil.update({'Wavelength(Microns)':Waves['Wavelength(Microns)']})

from bhmie_herbert_kaiser_july2012_GEOFF_EDITION import bhmie
Emiss = []
Time0 = times()[0]
sizey = len(Waves['SIZE'])
for iter, radius in enumerate(Waves['SIZE']):
	print('Emissivites for grains of radius: ' + str(radius) + ' microns. ' + str(iter + 1) + ' of ' + str(sizey) + '.\nTime for last radius: ' + str(times()[0] - Time0) + ' seconds')
	Time0 = times()[0]
	Emiss.append([radius])
	for i, lamda in enumerate(Waves['Wavelength(Microns)']):
		RI = AstroSil['Refractive Index'][i]
		if type(AstroSil['Extinction Coefficient'][i]) is not type(1j):
			EC = AstroSil['Extinction Coefficient'][i]*1j
		else:
			EC = AstroSil['Extinction Coefficient'][i]
		qabs = bhmie((2*pi*radius)/lamda, RI + EC, 2)
		if (qabs <= 0):
			qabs = 1.0e-11
		Emiss[iter].append(qabs)

writerE(Waves['Wavelength(Microns)'], Emiss, 'Emissivities[AstroSil].csv')
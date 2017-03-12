from kmod import listifier, whiteSpaceParser
import matplotlib.pyplot as plt
from math import log10

DEFAULTDIR = r'C:\Users\Cody\Desktop\School Documents\Physics\Independent Study\DielectricEffect\testfiles\\'

directory = raw_input('Please enter file directory: ')

fileName = raw_input('Please enter FULL file name: ')

if directory == '':
	directory = DEFAULTDIR

if fileName == '':
	exit()

ffile = open(directory + fileName, 'rb')

def readEmiss(FILE):
	FILE.seek(0)
	Lamda = listifier(whiteSpaceParser(FILE.readline()))
	Emiss = []
	for ele in FILE:
		Emiss.append(listifier(whiteSpaceParser(ele)))
	return Lamda, Emiss

Lamda, Emiss = readEmiss(ffile)

def flotMIE(EMISS):
		plt.xlabel('Wavelength(Microns)')
		plt.ylabel('Emissivity')
		plt.title('Emissivity Vs. Wavelength')
		plt.xscale('log')
		plt.yscale('log')
		plt.grid(True)
		for ARRAY in EMISS:
			if (log10(ARRAY[0]) == int(log10(ARRAY[0]))):
				plt.plot(Lamda, ARRAY[1:], label = ARRAY[0])
		plt.legend(loc = 0)
		plt.show()

flotMIE(Emiss)

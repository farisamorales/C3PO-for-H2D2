try:
	from kmod import listifier, whiteSpaceParser
	import matplotlib.pyplot as plt
	from math import log10
except ImportError:
	print('One or more required modules could not be found:\n\tkmod.listifier, kmod.whiteSpaceParser\n\tmatplotlib.pyplot\n\tmath.log10')
	raw_input('Press any key to exit')
	exit()

DEFAULTDIR = r'C:\Users\Cody\Desktop\School Documents\Physics\Independent Study\DielectricEffect\testfiles'
print('Current default directory: ' + DEFAULTDIR)
directory = raw_input('Enter the file directory (press Enter for default): ')

fileName = raw_input('Enter the designator (e.x. DI;AS;50-50); type HELP for help: ')

if directory == '':
	directory = DEFAULTDIR

if fileName.lower() == 'help':
	print('Designators take the form MATRIX[;INCLUSION;VOL_MATRIX-VOL_INCLUSION]\nWhere MATRIX or INCLUSION is the two letter symbol for\nthe matrix/inclusion substance (e.x. DI for DirtyIce)\nVOL_MATRIX is a number between 0 and 100\nAnd VOL_INCLUSION is 100 - VOL_MATRIX')
	fileName = raw_input('Enter the designator; type HELP for more help: ')
	if fileName.lower() == 'help':
		print('Current substances and there associated symbols:\n\tDirtyIce - DI\n\tWaterIce - WA\n\tAmorphousCarbon - AC\n\tAstroSil - AS')
		fileName = raw_input('Enter the designator; don\'t type HELP for even more help: ')
		if fileName.lower() == 'help':
			print('You are beyond help')
			raw_input('Press any key to exit')
			exit()

try:
	ffile = open(directory + '\Emissivities[' + fileName + '].csv', 'rb')
except IOError:
	print(directory + '\Emissivities[' + fileName + '].csv' + ' does not exist')
	raw_input('Press any key to exit')
	exit()
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

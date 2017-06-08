import pickle

infile = open('temptable.50-50','rb') 	#open the file temptable for reading('rb')


(temperatures, grainsizes, radii, Tdict) = pickle.load(infile)	#Load in the 'pickled' file. This file returns a tuple (single item consisting of multiple values of potentially varying data types. In this case the first item returned is a list (here assigned to jlist). The second is a dictionary (jtable).

#example 1: Retreive grainsize and radius for every temperature
#for T in temperatures:
#	print Tdict[T]

#temperatures is a list of temperatures, grainsizes naturally is a list of grainsizes. Tdict is a dictionary where the key is the temperature and the values are the radius and grainsize used to calculate that temperature. Temperature was calculated by interpolating between flux in and flux out of a particular grain. Emissivity is accounted for in the flux calculations as a function of grain size and Cody's grain properties.

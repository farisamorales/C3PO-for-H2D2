import pickle

infile = open('temptable.50-50','rb')

				# + cats 'filename:' with infile
(jlist, jtable) = pickle.load(infile)

if __name__ == '__main__':
	for T in jtable:
		print jtable[T]

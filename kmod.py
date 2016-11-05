'''
	Various modules with specific but limited usability.
'''

'''
	Takes a loosly formated string and reformats it to have only spaces as its white-
	space characters. Returns said string as well as a list of index positions in that
	string where a whitespace had been separating elements in that string.
	
	Author: Cody "This is the third time I've changed what it does" King
	LastUpdated: October 31, 2016
	Python Version: 2.7.12
'''
#--------------------------------------------------------------------------------------------------	
def whiteSpaceParser(master):

	#To hold the index position of whitespaces
	whiteSpaceLocal = []

	#Appends a space at the end of master to prevent future index problems
	#in whiteSpaceLocal provided one isn't already there
	if not master.endswith(' '):
		master = master + ' '
	if not master.startswith(' '):
		master = ' ' + master

	#Deletes all tabs, newlines, and carriage returns (\r) replacing them with whitespaces
	master=master.replace('\t',' ')
	master=master.replace('\n',' ')
	master=master.replace('\r',' ')

	#Deletes any redundant whitespaces by replacing all double whitespace with a single one
	#loops through master until all of them are gone
	while master.find('  ') != -1:
		master=master.replace('  ',' ')
	
	#Iterates through the master string and appends to whiteSpaceLocal the index
	#position of a found whitespace
	for i, ele in enumerate(master):
		if ele == ' ':
			whiteSpaceLocal.append(i)

	#Deletes all the whitespaces
	master=master.translate(None,' ')

	#Deleting all the whitespaces causes master to shrink and reposition where
	#the whitespaces where by an amount equal to the index of whiteSpaceLocal
	#We compensate with the following
	for i, val in enumerate(whiteSpaceLocal):
		whiteSpaceLocal[i]=val-i
	
	#Returns a tuple of the reformated string and a list of the index positions of the whitespaces
	#that where in that string
	return master, whiteSpaceLocal
#--------------------------------------------------------------------------------------------------	
'''
	Takes in a tuple formatted by whiteSpaceParser and finds all the numerical values
	in that string; returning a list made of those values

	Author: Cody "I made this one too" King
	LastUpdated: October 14, 2016
	Python Version: 2.7.12
'''
#--------------------------------------------------------------------------------------------------	
def listifier(tuply):
	i, masterList = 0, []
	while i < len(tuply[1]) - 1:
		#print(i)   #Testing stuff to see if the module was getting hung-up
		try:
			#Tries to append to masterList a float made from the substring in tuply[0] (a string)
			#designated by the boundaries tuply[1][i] and tuply[1][i+1] (both ints in the list tuply[1])
			masterList.append(float(tuply[0][tuply[1][i]:tuply[1][i+1]]))
		except ValueError:
			#If the substring isn't floatable we ignore it
			pass
		finally:
			i += 1
	del(tuply, i)
	return masterList
#--------------------------------------------------------------------------------------------------	
'''
	Takes in a master list and returns a list of lists from
	the elements of the input list based on the number of columns
	specified.
	
	Author: Cody "WE NEED MORE MODULES" King
	LastUpdated: October 14, 2016
	Python Version: 2.7.12
'''
#--------------------------------------------------------------------------------------------------	
def columnizer(mList, col):

	#i serves as the iterator through the columns
	i, theList = 0, []
	while i < col:
		theList.append([])
		
		#j serves as the iterator through the indices; must be reset to zero after every run
		j = 0
		while j < len(mList):
		
			#If the difference between the index (j) and the current column (i) is evenly
			#divisible by the number of columns (col) then the value of that index belongs
			#in that column
			if (j-i)%col == 0:
				theList[i].append(mList[j])
			j += 1
		i += 1
	del(i, j, mList, col)
	return theList	
#--------------------------------------------------------------------------------------------------	
'''
	For getStart and getCol
	Functions which get info about numerical table files
	
	Author: Cody "Why did I make this so complicated" King
	LastUpdated: October 31, 2016
	Python Version: 2.7.12
'''

'''
	Takes in a file object linked to a loosly formatted numerical table and
	returns a long of the index of where the table starts in that file.
	This file object should be opened in binary mode ('rb') to avoid reading
	issues which occur on the windows OS.
'''
#--------------------------------------------------------------------------------------------------	
def getStart(tFile):

	#Defining the boolean variables and a list to contain the start of each line
	first, second, lineIndex=False, False, [0L]
	
	#Reset the position of the reader incase the file has been read previously
	tFile.seek(0L)
	
	#The idea is if two subsequent lines of the file only contain elements that can be
	#converted to floats, then its likely that those lines are the start of the table.
	#This while loop will find it
	while (not first) and (not second):
		'''
		Reads in one line from the tFile object and calls wsp to generate a tuple of
		the parsed line (a string) and a list of where whitespaces had been slicing up
		that string. Then appends to lineIndex where the reader is (which is the beginning
		of the next line).
		'''
		toTest=tFile.readline()
		if toTest == '':
			print('Start of table not found')
			return None
		toTest = whiteSpaceParser(toTest)
		lineIndex.append(tFile.tell())

		#Debugging code
		#print(toTest)
		#print(lineIndex)
		#print(first, second)

		'''
		This while loop takes the toTest tuple and tests the elements in the string portion
		of said tuple using the list of index positions in the list portion. Starting from both
		first and second being false, it finds the first line in which all elements in the string
		portion are floatable and sets first to true before trying the next line. If the next line
		contains some non-floatable elements it resets first to false and starts the next line. Once
		both first and second are true, it's likely we've found the start of the table so we take
		lineIndex and return the second to last value in that list; this is the line position in
		tFile in which the table starts.
		'''
		i = 0
		while i < len(toTest[1]) - 1:
		
			#Checking to see if the elements are floatable
			try:
				float(toTest[0][toTest[1][i]:toTest[1][i+1]])
				i += 1
			
			#If not we start over from the next line
			except ValueError:
				if first:
					first, second=False, False
				break
			if first:
				second=True
			else:
				first=True
		continue
	#If we've made it out of the loop then we've found the starting position (probably)
	del(tFile, first, second, i)
	return lineIndex[-2]
#--------------------------------------------------------------------------------------------------		
'''
	Takes in a file object linked to a numerical table and uses getStart to find the start of the table
	then finds the number of columns designated in this table.
'''
#--------------------------------------------------------------------------------------------------	
def getCol(tFile):
	tFile.seek(getStart(tFile))
	toTest = whiteSpaceParser(tFile.readline())
	del(tFile)
	return len(toTest[1]) - 1
#--------------------------------------------------------------------------------------------------	
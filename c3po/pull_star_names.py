import os


fdir = "C:\\Users\\justi\\C3PO-for-H2D2\\c3po\\Data\\StarFiles\\"

files = os.listdir(fdir)

starNames = []
for i in range(len(files)):
    starNames.append( files[i].rstrip("_stitched.txt") )

print( fdir.rstrip("StarFiles\\") )

with open("Data\\totalStarNames.txt", "w+") as f:
    for name in starNames:
        f.write(name + '\n')

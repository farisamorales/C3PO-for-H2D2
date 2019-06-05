#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 00:12:15 2017

@author: marlydmejia
"""

import csv
import os.path

userhome = os.path.expanduser('~')
csvfile= os.path.join(userhome,"Desktop/C3PO_for_H2D2",'complete_simbad_photometric_.csv')
database1_file= open(csvfile, "rU")
a = csv.DictReader(database1_file)
#database1_file.close

csvfile= os.path.join(userhome,"Desktop/C3PO_for_H2D2", 'Herschel_diskList_Detections.csv')
database2_file= open(csvfile, "rU")
b = csv.DictReader(database2_file)
#database2_file.close

csvfile= os.path.join(userhome,"Desktop/C3PO_for_H2D2", 'GAIA_Official_precise.csv')
database3_file= open(csvfile, "rU")
c = csv.DictReader(database3_file)
#database3_file.close

csvfile= os.path.join(userhome,"Desktop/C3PO_for_H2D2", 'gator_allwise.tbl')
database4_file= open(csvfile, "rU")
d = csv.DictReader(database4_file)
#database4_file.close

csvfile= os.path.join(userhome,"Desktop/C3PO_for_H2D2", 'Identifiers_Lister.csv')
database5_file= open(csvfile, "rU")
e = csv.DictReader(database5_file)
#database5_file.close

Tmass= [row for row in a]     #these
Herschel= [row for row in b]  #those
Gaia= [row for row in c]      #more
Allwise= [row for row in d]   #cont
Namers= [row for row in e]    #component


#updating everything to the master Namers dictionary.
#Starts by adding 2MASS data
for component in Namers:
    for these in Tmass:
        if component['HDname'].strip()== these['typed ident'].strip():
                component.update(these)
                print component
       # elif component['HIPname'].strip()== these['typed ident'].strip():
        #        component.update(these)

         #       print component

"""    
#to that, we add Herschel                
for component in Namers:
    for those in Herschel:
        if component['HDname'].strip()== those['name'].strip():
                component.update(those)
        elif component['HIPname'].strip()== those['name'].strip():
                component.update(those)   
        elif component['HDname'].strip()== those['HIP name'].strip():
                component.update(those)
        elif component['HIPname'].strip()== those['HIP name'].strip():
                component.update(those)

                
                
#next, we add Allwise by matching the RA and DEC with the least uncertainties                
for component in Namers:
    for cont in Allwise:
        if abs(component['RAval'].strip()- cont['ra'].strip())<0.05:
                component.update(cont)
        elif abs(component['DECval'].strip()- cont['dec'].strip())<0.05:
                component.update(cont)          
        
"""        
        


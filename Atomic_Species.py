#!/usr/bin/env python2
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
# AUTHOR:
#            BEKKALI HAMZA 
# EXPLANATION:
'''
-----------------------------
>$  Atomic_Species.py  Mo Se Te bi
...
pseudo_dir = '/home/hp/QE/SSSP_v1.1/SSSP_efficiency/'
ATOMIC_SPECIES
  Mo     95.95   Mo_ONCV_PBE-1.0.oncvpsp.upf
  Se     78.971  Se_pbe_v1.uspp.F.UPF
  Te     126.7   Te_pbe_v1.uspp.F.UPF
  bi     208.98  Bi_pbe_v1.uspp.F.UPF
-----------------------------
'''
#____________________________________________________________________________________

from sys   import stdin
from numpy import *
import numpy as np
import subprocess
import os.path
import shutil
import glob
import math
import sys
import os
import sys
#
from  itertools import islice
folder_path = '/home/hp/QE/SSSP_v1.1/SSSP_efficiency/'  # path where the pesudo potentials located 
# Dictionary of all elements matched with their atomic masses.
elements_dict = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}
# List of all elements to allow for easy inclusion testing.
elements_list2 = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA',\
                 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI',\
                 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE',\
                 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO',\
                 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE',\
                 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM',\
                 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF',\
                 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB',\
                 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U',\
                 'NP', 'PU', 'AM', 'CM', 'BK', 'CT', 'ES', 'FM', 'MD', 'NO',\
                 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN',\
                 'NH', 'FL', 'MC', 'LV', 'TS', 'OG']

elements = sys.argv[1:]
extrain=[]
el=[]  
def extn(extrain):
    global elements 
    print(" you are here ")
    if isinstance(extrain,list) and len(extrain)>0:
        elements=extrain
        return elements
    else: 
        print(";:")
        elements = sys.argv[1:]
        return elements
 
def extn_out(list_out,psudodir):
    return list_out,psudodir
 
if len(elements) == 0:
    # Print an error message and usage instructions
    print("Error: Usage: python Masse_atomic.py [elemnt1] [element2] [element3]")
    print("       Usage:        Masse_atomic.py  Mo S Se ... almost one element")
    sys.exit(1)
 
# Check if there are three arguments
# Split the arguments into separate variables
elements_tags=elements
# For each tag name, search for files that start with the tag name
data_list = [] 
for element in elements_tags:
    element_file = None  # initialize variable to hold file name for current element
    if element.upper() in elements_list2:
          print("Mass Atomic of  {:>5} = {:<13}".format(element.capitalize(),elements_dict[element.upper()]))
          mass_atomic=elements_dict[element.upper()]
    else:
          print(element.capitalize(),"Element not found !!")
#######
    for file_name in os.listdir(folder_path):
        if (file_name.startswith(element.capitalize()) or file_name.startswith(element.lower())) and (file_name.find(".") == len(element) or file_name.find("_") == len(element)): # condition to check if file_name starts with the element name and has either a period or underscore after
            element_file = file_name
            break
    data_list.append((element, mass_atomic, element_file))    
    if element_file is not None:
        print("\t\t   {} -> {}".format(element.capitalize(),element_file))
    else:
        print("No pseudo found here for element ", element)
        print(" here some suggestions : ") 
        if file_name.startswith(str(element.capitalize())):
            print(str(file))
        if file_name.startswith(str(element.lower())):
            print(str(file))

#--print card "ATOMIC_SPECIES" / QE input    
 
print ("\n-----------------------------\n")
print ("pseudo_dir = '{}'".format(folder_path))
print("ATOMIC_SPECIES")
for element  in data_list:
     print("  {:<5}{:^10}{:>13}".format(element[0],element[1],element[2]))
print ("-----------------------------")
print("\n\n !!>: Be aware and ensure that the pseudonum of the elements is what you intended.")

print ("I hope that helps!")
	
extn_out(data_list,folder_path)



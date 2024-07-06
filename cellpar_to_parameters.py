#!/usr/bin/env python3
#####################################################
# @@ Bekk@14 
#          THIS PROGRAM HELP ME TO EXTRACT the Crystalo parametres a b c cosbc cosac cosab           
#                  FROM CELL_PARAMETERS Matrix 
######################################################
import subprocess
import os.path
import shutil
import glob
import math
import sys
import os
from  numpy import *
import numpy as np

# Create a 3x3 matrix filled with zeros
CP = np.zeros((3, 3))

# Prompt the user to input each line for the matrix
for i in range(3):
    line = input(f"Enter the values for row {i+1} (e.g., 0.000 0.000 0.000): ")
    values = line.split()
    for j in range(3):
        CP[i, j] = float(values[j])


print('_='*26)
print()
print("CELL_PARAMETERS")
for row in CP:
    print(f"{row[0]:16.10f} {row[1]:16.10f} {row[2]:16.10f}")

print('_='*26)
celldm1=1.0
a1 = celldm1 * sqrt(CP[0,0]**2 + CP[0,1]**2 + CP[0,2]**2)
a2 = celldm1 * sqrt(CP[1,0]**2 + CP[1,1]**2 + CP[1,2]**2)
a3 = celldm1 * sqrt(CP[2,0]**2 + CP[2,1]**2 + CP[2,2]**2)
alpha = degrees(math.acos((CP[1,0]*CP[2,0]+CP[1,1]*CP[2,1]+CP[1,2]*CP[2,2])*celldm1**2/(a2*a3)))
beta  = degrees(math.acos((CP[0,0]*CP[2,0]+CP[0,1]*CP[2,1]+CP[0,2]*CP[2,2])*celldm1**2/(a1*a3)))
gamma = degrees(math.acos((CP[0,0]*CP[1,0]+CP[0,1]*CP[1,1]+CP[0,2]*CP[1,2])*celldm1**2/(a1*a2)))


print('_='*26)
cell_parameters = np.array([a1,a2,a3,alpha,beta,gamma])
# Assuming the values of bohrToA and cell_parameters are defined elsewhere
bohrToA =1.0  #0.52917720859

print(f'a = {bohrToA*cell_parameters[0]:>16.11f}         alpha   = {cell_parameters[3]:>5.2f}  (cosbc)')
print(f'b = {bohrToA*cell_parameters[1]:>16.11f}         beta    = {cell_parameters[4]:>5.2f}  (cosac) ')
print(f'c = {bohrToA*cell_parameters[2]:>16.11f}         gamma   = {cell_parameters[5]:>5.5f} (cosab)')

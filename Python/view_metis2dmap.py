# ********************************************************************
# Copyright (c) 2021   Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the PSUADE team.
# All rights reserved.
#
# Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
# disclaimer, contact information and the GNU Lesser General Public 
# License.
#
# This is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License (as published by 
# the Free Software Foundation) version 2.1 dated February 1999.
#
# This software is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF 
# MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 
# USA
# *******************************************************************
# This file contains Python codes that read metisMeshPlot.m (with a
# few lines deleted) created by PSUADE and display the 2D METIS
# partitioning plot using matplotlib
# *******************************************************************
import os
import sys
import numpy as np 
try:
   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import axes3d, Axes3D
   MATPLOTLIB = True
except ImportError:
   print('INFO: matplotlib unavailable')
   MATPLOTLIB = False
   sys.exit(0)

print('************************************************************')
print('A file created by metis_2dmap (metisMeshPlot.m) is needed.')
print('NOTE: you must delete the first line and the last few lines')
print('      for this code to work.')
print('------------------------------------------------------------')
if sys.version_info[0] < 3:
   fname = raw_input('File created by metis_2dmap : ')
else:
   fname = input('File created by metis_2dmap : ')

fdata = np.genfromtxt(fname)
fig, ax = plt.subplots()
plt.imshow(fdata)
plt.xlabel('X1',fontsize=14,fontweight='bold')
plt.ylabel('X2',fontsize=14,fontweight='bold')
plt.title('METIS Partitioning',fontsize=14,fontweight='bold')
for axis in ['top', 'bottom', 'left', 'right']:
   ax.spines[axis].set_linewidth(4)

plt.savefig('psu_metisPartition2d.png')
print('Metis partitioning plot is now in psu_metisPartition2d.png')
plt.show()


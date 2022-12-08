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
# This file contains Python codes that read MCMCPostSample created by
# PSUADE and display the posterior histogram and heat plots
# *******************************************************************
import os
import sys
try:
   import numpy as np 
except:
   print('ERROR: numpy not available.')
   sys.exit(1)
try:
   import matplotlib.pyplot as plt
   import matplotlib.mlab as mlab
except:
   print('ERROR: matplotlib not available.')
   sys.exit(0)

# *******************************************************************
# Plot parameters that users may change
# -------------------------------------------------------------------
nhist = 20
numXTicks = 3
numYTicks = 2
# *******************************************************************

print('************************************************************')
print('This script takes a MCMC posterior (MCMCPostSample) file and')
print('plots the parameter posteriors.')
print('------------------------------------------------------------')

if sys.version_info[0] < 3:
   mcmcfile = raw_input('Enter name of MCMC posterior sample file : ')
else:
   mcmcfile = input('Enter name of MCMC posterior sample file : ')
if not os.path.isfile(mcmcfile):
   print('ERROR: mcmc posterior file not found')
   sys.exit(0)

# read the posterior sample file
try:
   with open(mcmcfile, 'r') as infile:
      mcmcdata = infile.readlines();
except:
   print('ERROR: mcmc posterior file not found')
   sys.exit(0)

if len(mcmcdata) < 2:
   print('ERROR: Invalid posterior sample file format')
   print('       Number of lines in file = ' + str(len(mcmcdata)))
   sys.exit(0)

# find read in the posterior sample
# line 1 should have 'PSUADE_BEGIN'
lineIn = mcmcdata[0]
cols = lineIn.split()
if len(cols) != 1 and cols[0] != 'PSUADE_BEGIN':
   print('ERROR: Invalid posterior sample file format')
   print('       The first line must be PSUADE_BEGIN') 
   print('Your line 1 = ' + lineIn)
   sys.exit(0)

# line 2 should have <nSamples nInps>
lineIn = mcmcdata[1]
cols = lineIn.split()
if len(cols) != 2:
   print('ERROR: Invalid posterior sample file format')
   print('       The second line must be <nSamples nInputs>')
   print('Your line 2 = ' + lineIn)
   sys.exit(0)
try:
   nSamp = int(cols[0])
   nInps = int(cols[1])
except:
   print('ERROR: Invalid posterior sample file format')
   print('       The second line must be <nSamples nInputs>')
   print('Your line 2 = ' + lineIn)
   sys.exit(0)
if nSamp <= 0 or nInps <= 0:
   print('ERROR: nSamples or nInputs <= 0')
   print('       nSamples = ' + str(nSamp))
   print('       nInputs  = ' + str(nInps))
   sys.exit(0)

# line 3 should have variable names
lineIn = mcmcdata[2]
vnames = lineIn.split()
if len(vnames) != nInps+1:
   print('ERROR: Invalid posterior sample file format')
   print('       The third line must be <# variable names>')
   print('Your line 3 = ' + lineIn)
   sys.exit(0)

# line 4 on should have the posterior sample
postSample = np.zeros([nInps, nSamp])
for ii in range(nSamp):
   lineIn = mcmcdata[ii+3]
   cols = lineIn.split()
   if len(cols) < nInps+1:
      print('ERROR: Invalid posterior sample file format')
      print('       Offending line ' + str(ii+4) + ' = ' + lineIn) 
      sys.exit(0)
   try:
      ind = int(cols[0])
      if ind != ii+1:
         print('ERROR: Invalid posterior sample file format')
         print('       Offending line ' + str(ii+4) + ' = ' + lineIn) 
         sys.exit(0)
   except:
      print('ERROR: Invalid posterior sample file format')
      print('       Offending line ' + str(ii+4) + ' = ' + lineIn) 
      sys.exit(0)

   for jj in range(nInps):
      try:
         postSample[jj][ii] = float(cols[jj+1])
      except:
         print('ERROR: Invalid posterior sample file format')
         print('       Offending line ' + str(ii+4) + ' = ' + lineIn) 
         sys.exit(0)

# get upper and lower bounds
lbs = []
ubs = []
for ii in range(nInps):
   print('Posterior sample input ' + str(ii) + ' has a range of ')
   print('          min = ' + str(np.min(postSample[ii])))
   print('          max = ' + str(np.max(postSample[ii])))
   if sys.version_info[0] < 3:
      bndStr = raw_input('Enter lower bound for input ' + str(ii+1) + ' : ')
   else:
      bndStr = input('Enter lower bound for input ' + str(ii+1) + ' : ')
   try:
      lbs.append(float(bndStr))
   except:
      print('ERROR: input bound entered incorrectly.')
      sys.exit(1)
   if sys.version_info[0] < 3:
      bndStr = raw_input('Enter upper bound for input ' + str(ii+1) + ' : ')
   else:
      bndStr = input('Enter upper bound for input ' + str(ii+1) + ' : ')
   try:
      ubs.append(float(bndStr))
   except:
      print('ERROR: input bound entered incorrectly.')
      sys.exit(1)
   if lbs[ii] >= ubs[ii]:
      print('ERROR: input lower bound >= upper bound.')
      sys.exit(1)

# now plot the histogram 
fig, ax = plt.subplots()
for ii in range(nInps):
   plt.subplot(nInps,nInps,ii*nInps+ii+1)
   counts, bins = np.histogram(postSample[ii],nhist)
   probs = 1.0 * counts / len(postSample[ii])
   plt.hist(bins[:-1], bins, weights=probs)
   mu  = np.mean(postSample[ii])
   sig = np.std(postSample[ii])
   print('Input = ' + vnames[ii+1] + ' : mean, stdev = ' + str(mu) + \
         ', ' + str(sig))

   # to be replaced by scipy's norm.pdf function
   #normalFit = mlab.normpdf(bins, mu, sig)
   #plt.plot(bins, normalFit, 'r--')
   plt.xlabel(vnames[ii+1], fontsize=8,fontweight='bold')
   swidth = '%7.2e' % (bins[1]-bins[0])
   #pstr = 'Probability/(width=' + swidth + ')'
   pstr = 'Probability'
   plt.ylabel(pstr, fontsize=8,fontweight='bold')
   # should use true prior ranges, not these
   #xmax = np.max(postSample[ii])
   #xmin = np.min(postSample[ii])
   xmax = ubs[ii]
   xmin = lbs[ii]
   xstep = (xmax-xmin)/numXTicks
   plt.xticks(np.arange(xmin,xmax+0.01*xstep,xstep))
   #plt.xticks([],[])
   ymax = np.max(probs)
   ymin = 0
   ystep = (ymax-ymin)/numYTicks
   plt.yticks(np.arange(ymin,ymax+0.01*ystep,ystep))
   #plt.yticks([],[])
   plt.xlim([xmin,xmax])
   plt.ylim([0,ymax])
   plt.grid(b=True)
   for axis in ['top', 'bottom', 'left', 'right']:
      ax.spines[axis].set_linewidth(4)

# now plot 2D histogram 
for ii in range(nInps):
   for jj in range(nInps):
      if jj > ii:
         X = postSample[ii]
         Y = postSample[jj]
         plt.subplot(nInps,nInps,ii*nInps+jj+1)
         plt.hist2d(Y,X,bins=50)
         #xmax = np.max(postSample[jj])
         #xmin = np.min(postSample[jj])
         xmax = ubs[jj]
         xmin = lbs[jj]
         xstep = (xmax-xmin)/numXTicks
         plt.xticks(np.arange(xmin,xmax+0.01*xstep,xstep))
         #plt.xticks([],[])
         #ymax = np.max(postSample[ii])
         #ymin = np.min(postSample[ii])
         ymax = ubs[ii]
         ymin = lbs[ii]
         ystep = (ymax-ymin)/numXTicks
         plt.yticks(np.arange(ymin,ymax+0.01*ystep,ystep))
         #plt.yticks([],[])
         plt.xlim([xmin, xmax])
         plt.ylim([ymin, ymax])
         plt.grid(b=True)
         for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(4)
         #plt.colorbar()
         #plt.xlabel(vnames[jj+1])
         #plt.ylabel(vnames[ii+1])
tStr = 'Individual and Pairwise Posteriors'
plt.suptitle(tStr, fontsize=15,fontweight='bold')

savefile = 'psu_mcmcpost.png'
plt.savefig(savefile, bbox_inches='tight')
print('Posterior distribution plot is in psu_mcmcpost.png')
print('Dashed red line in histograms: normal distribution approximation')
#plt.show()
#plt.axis('off')

# look to see if there are more information
# 1. search for marker ('MLE')
ii = nSamp + 3
mleMarker = -1
while ii < len(mcmcdata):
   lineIn = mcmcdata[ii]
   cols = lineIn.split()
   if len(cols) >= 1 and cols[0] == 'MLE':
      mleMarker = ii
      break
   else:
      ii = ii + 1

# Ask for more plots
askStr = 'n'
if mleMarker != -1:
   if sys.version_info[0] < 3:
      askStr = raw_input('Generate parity plot for experiments ? (y or n) ')
   else:
      askStr = input('Generate parity plot for experiments ? (y or n) ')

# Search for individual terms in the posteriors
if askStr == 'y':

# 2. find the number of outputs
   nOuts = 0
   ii = mleMarker + 1
   while ii < len(mcmcdata): 
      lineIn = mcmcdata[ii]
      cols = lineIn.split()
      col1 = cols[0]
      ind = int(cols[0])
      if ind == 1:
         nOuts = nOuts + 1
         ii = ii + 1
      else:
         break
# 3. find the number of experiments
   ii = mleMarker + 1
   nExp = 0
   while ii < len(mcmcdata): 
      lineIn = mcmcdata[ii]
      cols = lineIn.split()
      if len(cols) < 5:
         break
      else:
         nExp = nExp + 1
         ii = ii + 1
   nExp = int(nExp / nOuts)

# 4. now read individual likelihood
   likelihoods = np.zeros([nOuts, nExp])
   expData = np.zeros([nOuts, nExp])
   predData = np.zeros([nOuts, nExp])
   offset = mleMarker + 1
   for ii in range(nExp):
      for jj in range(nOuts):
         lineIn = mcmcdata[offset+ii*nOuts+jj]
         cols = lineIn.split()
         likelihoods[jj][ii] = float(cols[4])
         predData[jj][ii] = float(cols[1])
         expData[jj][ii] = float(cols[2])

# plot
   nplt1D = np.sqrt(nOuts)
   pltNR = int(nplt1D) 
   pltNC = int(nplt1D)
   if pltNR*pltNC < nOuts:
      pltNC = pltNC + 1
   if pltNR*pltNC < nOuts:
      pltNR = pltNR + 1
   fig, ax = plt.subplots()
   for ii in range(pltNR):
      for jj in range(pltNC):
         ind = ii * pltNC + jj
         if ind < nOuts:
            plt.subplot(pltNR,pltNC,ind+1)
            plt.scatter(expData[ind],predData[ind],s=20,marker='*')
            xmax = np.max(expData[ind])
            if np.max(predData[ind]) > xmax:
               xmax = np.max(predData[ind])
            xmin = np.min(expData[ind])
            if np.min(predData[ind]) < xmin:
               xmin = np.min(predData[ind])
            plt.xlim([xmin, xmax])
            plt.ylim([xmin, xmax])
            xstep = (xmax-xmin)/3
            plt.xticks(np.arange(xmin,xmax+0.01*xstep,xstep),fontsize=13,fontweight='bold')
            plt.yticks(np.arange(xmin,xmax+0.01*xstep,xstep),fontsize=13,fontweight='bold')
            plt.plot([xmin,xmax],[xmin,xmax],'r-',linewidth=2)
            plt.xlabel('Experiment',fontsize=13,fontweight='bold')
            plt.ylabel('Prediction',fontsize=13,fontweight='bold')
            plt.grid(b=True)
            # thicken the border lines
            for axis in ['top', 'bottom', 'left', 'right']:
               ax.spines[axis].set_linewidth(6)
            tStr = 'Output ' + str(ind+1)
            plt.title(tStr,fontsize=14,fontweight='bold')
   savefile = 'psu_mcmcparity.png'
   plt.savefig(savefile, bbox_inches='tight')
plt.show()
sys.exit(0)

# finally plot the data 
#plt.figure(2)
fig, ax = plt.subplots()
nhist = 20
for ii in range(nOuts):
   plt.subplot(1,nOuts,ii+1)
   counts, bins = np.histogram(likelihoods[ii],nhist)
   probs = 1.0 * counts / len(likelihoods[ii])
   plt.hist(bins[:-1], bins, weights=probs)
   xStr = 'Output ' + str(ii+1)
   plt.xlabel(xStr, fontsize=12,fontweight='bold')
   yStr = 'Probability'
   plt.ylabel(yStr,fontsize=12,fontweight='bold')
   xmax = np.max(likelihoods[ii])
   xmin = np.min(likelihoods[ii])
   xstep = (xmax-xmin)/2
   plt.xticks(np.arange(xmin,xmax+0.01*xstep,xstep))
   ymax = np.max(probs)
   ymin = 0
   ystep = (ymax-ymin)/2
   plt.yticks(np.arange(ymin,ymax+0.01*ystep,ystep))
   plt.xlim([xmin,xmax])
   plt.ylim([0,ymax])
   plt.grid(b=True)
   for axis in ['top', 'bottom', 'left', 'right']:
      ax.spines[axis].set_linewidth(4)
tStr = 'Distribution of Likelihood for Each Output'
plt.suptitle(tStr, fontsize=15,fontweight='bold')



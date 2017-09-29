// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC. 
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the multivariate normal distribution
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Globals.h"
#include "PsuadeUtil.h"
#include "PDFMVNormal.h"
#include "PDFNormal.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFMVNormal::PDFMVNormal(Vector &means, Matrix &corMat)
{
   int    ii, jj, length, status;
   double sigmaI, sigmaJ, ddata, *localSds;;

   if (means.length() <= 0)
   {
      printf("PDFMVNormal ERROR: number of inputs <= 0.\n");
      exit(1);
   }
   if (means.length() != corMat.ncols())
   {
      printf("PDFMVNormal ERROR: mismatch mean and cov dimensions.\n");
      printf("                   means has dimension %d.\n",means.length());
      printf("                   cov   has dimension %d.\n",corMat.ncols());
      exit(1);
   }
   length = means.length();
   for (ii = 0; ii < length; ii++)
   {
      ddata = corMat.getEntry(ii,ii);
      if (ddata <= 0)
      {
         printf("PDFMVNormal ERROR: stdev should be > 0 (%d,%e).\n",
                ii+1,ddata);
         exit(1);
      }
   }

   means_.load(means);
   localSds = new double[length];
   for (ii = 0; ii < length; ii++)
      localSds[ii] = corMat.getEntry(ii,ii);
   covMat_.setDim(length,length);
   for (ii = 0; ii < length; ii++)
   {
      sigmaI = localSds[ii];
      for (jj = 0; jj < length; jj++)
      {
         sigmaJ = localSds[jj];
         if (ii == jj) ddata = 1.0;
         else          ddata = corMat.getEntry(ii,jj);
         ddata *= sigmaI * sigmaJ;
         covMat_.setEntry(ii,jj,ddata);
      }
   }
   if (psAnaExpertMode_ == 1)
   {
      printf("Covariance Matrix:\n");
      for (ii = 0; ii < length; ii++)
   {
         for (jj = 0; jj < length; jj++)
            printf("%12.4e ", covMat_.getEntry(ii,jj));
         printf("\n");
      }
   }

   status = covMat_.CholDecompose();
   if (status != 0)
   {
      printf("PDFMVNormal ERROR: covariance matrix not positive definite.\n");
      exit(1);
   }
   delete [] localSds;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFMVNormal::~PDFMVNormal()
{
}

// ************************************************************************
// draw random samples 
// ------------------------------------------------------------------------
int PDFMVNormal::genSample(int length, Vector &vecOut, Vector &lower,
                           Vector &upper)
{
   int    ii, jj, nInputs, count, flag, total;
   double iZero=0.0, One=1.0, *localData1, *localData2;
   double low=-4, high=4;
   Vector localVec;
   PDFNormal *normalPtr;

   if (length <= 0)
   {
      printf("PDFMVNormal genSample ERROR - length <= 0.\n");
      exit(1);
   }
   nInputs = means_.length();
   if (vecOut.length()/length != nInputs)
   {
      printf("PDFMVNormal genSample ERROR - system error.\n");
      exit(1);
   }

   normalPtr  = new PDFNormal(iZero, One);
   localData1 = new double[length];
   localData2 = new double[length*nInputs];
   count = total = 0;
   while (count < length)
   {
      for (jj = 0; jj < nInputs; jj++)
      {
         normalPtr->genSample(length,localData1,low,high);
         for (ii = 0; ii < length; ii++)
            localData2[ii*nInputs+jj] = localData1[ii];
      }
      for (ii = 0; ii < length; ii++)
      {
         localVec.load(nInputs, &(localData2[ii*nInputs]));
         covMat_.CholMatvec(localVec);
         for (jj = 0; jj < nInputs; jj++)
            vecOut[count*nInputs+jj] = means_[jj] + localVec[jj];
         flag = 0;
         for (jj = 0; jj < nInputs; jj++)
            if (vecOut[count*nInputs+jj] < lower[jj] ||
                vecOut[count*nInputs+jj] > upper[jj]) {flag = 1; break;}
         if (flag == 0) count++;
         total++;
         if (count >= length) break;
      }
      if (total > length*1000)
      {
         printf("PDFMVNormal ERROR : cannot generate enough sample points\n");
         printf("                    within the given ranges. Check your \n");
         printf("                    parameter lower and upper bounds.\n");
         for (jj = 0; jj < nInputs; jj++)
            printf("Input bounds = %e %e \n", lower[jj], upper[jj]);
         exit(1);
      }
   }
   delete [] localData1; 
   delete [] localData2; 
   delete normalPtr;
   return 0;
}


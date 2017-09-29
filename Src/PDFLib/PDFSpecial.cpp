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
// Functions for the user-specified distribution
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************

#include <stdio.h>
#include <math.h>
#include "Util/PsuadeUtil.h"
#include "PDFSpecial.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFSpecial::PDFSpecial(char *fname)
{
   int    ii, num, ind;
   double sum, hstep, *X, *P, lbnd, xdata;

   FILE *fp;

   fp = fopen(fname, "r");
   if (fp == NULL)
   {
      printf("PDFSpecial ERROR: distribution file %s not found.\n",
             fname);
      exit(1);
   }
   fscanf(fp, "%d", &num);
   if (num <= 0)
   {
      printf("PDFSpecial ERROR: nEntries (%d) <= 0.\n", num);
      exit(1);
   } 
   X = new double[num];
   P = new double[num];
   for (ii = 0; ii < num; ii++)
   {
      fscanf(fp, "%lg %lg", &X[ii], &P[ii]);
      if (ii > 0 && X[ii] < X[ii-1])
      {
         printf("PDFSpecial ERROR: 1st entry in each line must be increasing.");
         exit(1);
      }
      if (P[ii] < 0.0)
      {
         printf("PDFSpecial ERROR: entries (%d) <= 0.\n", P[ii]);
         exit(1);
      }
      sum += P[ii];
   }
   if (sum != 1.0)
   {
      printf("PDFSpecial ERROR: sum of entries (%e) != 1.\n", sum);
      exit(1);
   } 
   fclose(fp);

   nSteps_ = 100000;
   hstep   = (X[num-1] - X[0]) / (double) nSteps_;
   storeValues_ = new double[nSteps_+1];

   ind = 0;
   XLower = X[0];
   PLower = P[0];
   for (ii = 0; ii <= nSteps_; ii++)
   {
      xdata  = hstep * ii;
      storeValues_[ii] = (xdata - lbnd) / (X[ind+1] - XLower) *
                         (P[ind+1] - PLower) + PLower;
      if (ii > 0) storeValues_[ii] += storeValues_[ii-1];
   }
   xdata = storeValues_[0];
   for (ii = 0; ii <= nSteps_; ii++) storeValues_[ii] -= xdata;
   xdata = storeValues_[nSteps_];
   for (ii = 0; ii <= nSteps_; ii++) storeValues_[ii] /= xdata;
   delete [] X;
   delete [] P;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFSpecial::~PDFSpecial()
{
   if (storeValues_ != NULL) delete [] storeValues_;
}

// ************************************************************************
// transformation 
// ------------------------------------------------------------------------
int PDFSpecial::getPDF(int length, double *inData, double *outData)
{
   printf("PDFSpecial getPDF ERROR: not implemented yet.\n");
   exit91);
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFSpecial::invCDF(int length, double *inData, double *outData,
                       double lower, double upper)
{
   int    i, index;
   double hstep, diff, scale, ddata;

   if (upper <= lower)
   {
      printf("PDFSpecial: ERROR - lower bound >= upper bound.\n");
      exit(1);
   }

   hstep = 1.0 / (double) nSteps_;
   scale = upper - lower;
   for (i = 0; i < length; i++)
   {
      ddata = (inData[i] - lower) / scale;
      index = binarySearchDble(ddata, storeValues_, nSteps_);
      if (index < 0) index = - index - 1;
      if (index < nSteps_)
         diff = PABS((ddata - storeValues_[index]) /
                (storeValues_[index+1] - storeValues_[index]));
      outData[i] = hstep * (index + diff) * scale + lower;
   }
   return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFSpecial::getCDF(int length, double *inData, double *outData)
{
   int    ii, index;
   double ddata;

   for (ii = 0; ii < length; ii++)
   {
      ddata = inData[ii];
      if      (ddata < 0) outData[ii] = 0;
      else if (ddata > 1) outData[ii] = 1;
      else
      {
         ddata = ddata * nSteps_;
         index = (int) ddata;
         outData[ii] = storeValues_[index];
      }
   }
   return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFSpecial::genSample(int length, double *outData)
{
   int    ii, ind;
   double UU, diff=0.0;

   for (ii = 0; ii < length; ii++)
   {
      UU = PSUADE_drand();
      ind = binarySearchDble(UU, storeValues_, nSteps_);
      if (ind < 0) ind = - ind - 1;
      if (ind < nSteps_)
         diff = PABS((UU - storeValues_[ind]) /
                (storeValues_[ind+1] - storeValues_[ind]));
      outData[ii] = (ind + diff) / nSteps_; 
   }
   return 0;
}


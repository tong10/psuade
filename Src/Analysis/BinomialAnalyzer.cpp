// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class BinomialAnalyzer (binomial test)
// Reference: "Hypothesis Testing Procedures for Eliminating Variables of 
//            Negligible Impact in Large Scale Computer Simulations"
//            by Jason Lenderman
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "BinomialAnalyzer.h"
#include "Util/sysdef.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
BinomialAnalyzer::BinomialAnalyzer() : Analyzer()
{
   setName("BINOMIAL");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BinomialAnalyzer::~BinomialAnalyzer()
{
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double BinomialAnalyzer::analyze(aData &adata)
{
   int     nSamples, ss, nBelow, printLevel, info;
   double  *Y, p0=0.5, thresh, *Ybin, *BinomialCDF, typeI;

   printLevel = adata.printLevel_;
   nSamples   = adata.nSamples_;
   Y          = adata.sampleOutputs_;
   thresh     = adata.analysisThreshold_;

   if (nSamples <= 1)
   {
      printf("BinomialAnalyzer INFO: not meaningful to ");
      printf("do this when nSamples <= 1.\n");
      return PSUADE_UNDEFINED;
   } 
   if (Y == NULL)
   {
      printf("BinomialAnalyzer ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   } 
   info = 0;
   for (ss = 0; ss < nSamples; ss++)
      if (Y[ss] == PSUADE_UNDEFINED) info++;
   if (info > 0)
   {
      printf("BinomialAnalyzer ERROR: Some outputs are undefined.\n");
      printf("                        Prune them first before analyze.\n");
      return PSUADE_UNDEFINED;
   }
   
   Ybin = new double[nSamples];
   for (ss = 0; ss < nSamples; ss++)
   {
      if (PABS(Y[ss]) <= thresh) Ybin[ss] = 1.0;
      else                       Ybin[ss] = 0.0;
   }

   nBelow = 0;
   for (ss = 0; ss < nSamples; ss++) if (Ybin[ss] == 1.0) nBelow++;
   if (printLevel > 2)
      printf("   ====> BinomialAnalyzer: number below threshold = %d (%d)\n",
             nBelow, nSamples);

   BinomialCDF = setupBinomialCDF(nSamples, p0);
   for (ss = 0; ss <= nSamples; ss++)
      BinomialCDF[ss] = 1.0 - BinomialCDF[ss];

   typeI = BinomialCDF[nBelow];

   delete [] Ybin;
   delete [] BinomialCDF;
   return typeI;
}

// *************************************************************************
// create a binomial cumulative density function
// -------------------------------------------------------------------------
double *BinomialAnalyzer::setupBinomialCDF(int n, double p0)
{
   int    k;
   double *retValues, nFact, kFact, nkFact;

   nFact = factorial(n);
   retValues = new double[n+1];
   for (k = 0; k <= n; k++)
   {
      kFact  = factorial(k);
      nkFact = factorial(n-k);
      retValues[k] = pow(p0,1.0*k) * pow(1-p0,1.0*(n-k));
      retValues[k] *= (nFact / (kFact * nkFact));
      if (k > 0) retValues[k] += retValues[k-1];
   }
   return retValues;
}

// *************************************************************************
// Compute factorial
// -------------------------------------------------------------------------
double BinomialAnalyzer::factorial(int n)
{
   int    ii;
   double fact;

   fact = 1.0;
   if (n == 0) return fact;
   for (ii = 2; ii <= n; ii++) fact *= (double) ii;
   return fact;
}


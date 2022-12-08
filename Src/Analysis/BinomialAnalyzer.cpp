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
#include <algorithm>

#include "BinomialAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
BinomialAnalyzer::BinomialAnalyzer() : Analyzer(), nSamples_(0), 
                                       nBelow_(0), typeI_(0)
{
  setName("BINOMIAL");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BinomialAnalyzer::~BinomialAnalyzer()
{
  VecYBins_.clean();
  VecBinomialCDF_.clean();
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double BinomialAnalyzer::analyze(aData &adata)
{
  int    ss, info;
  double p0=0.5;

  //**/ ---------------------------------------------------------------
  //**/ extract data
  //**/ ---------------------------------------------------------------
  int printLevel = adata.printLevel_;
  nSamples_ = adata.nSamples_;
  double *Y = adata.sampleOutputs_;
  double thresh = adata.analysisThreshold_;

  //**/ ---------------------------------------------------------------
  //**/ error checking and clean up
  //**/ ---------------------------------------------------------------
  if (nSamples_ <= 1)
  {
    printOutTS(PL_ERROR, "BinomialAnalyzer INFO: not meaningful to ");
    printOutTS(PL_ERROR, "do this when nSamples_ <= 1.\n");
    return PSUADE_UNDEFINED;
  } 
  if (Y == NULL)
  {
    printOutTS(PL_ERROR, "BinomialAnalyzer ERROR: no data.\n");
    return PSUADE_UNDEFINED;
  } 
  info = 0;
  for (ss = 0; ss < nSamples_; ss++)
    if (Y[ss] == PSUADE_UNDEFINED) info++;
  if (info > 0)
  {
    printOutTS(PL_ERROR,
         "BinomialAnalyzer ERROR: Some outputs are undefined.\n");
    printOutTS(PL_ERROR,
         "                        Prune them first before analyze.\n");
    return PSUADE_UNDEFINED;
  }

  //**/ ---------------------------------------------------------------
  //**/ clean up
  //**/ ---------------------------------------------------------------
  VecYBins_.clean();
  VecBinomialCDF_.clean();
   
  //**/ ---------------------------------------------------------------
  //**/ fill VecYBins array (1 - not important, 0 - significant)
  //**/ ---------------------------------------------------------------
  VecYBins_.setLength(nSamples_);
  for (ss = 0; ss < nSamples_; ss++)
  {
    if (PABS(Y[ss]) <= thresh) VecYBins_[ss] = 1.0;
    else                       VecYBins_[ss] = 0.0;
  }

  //**/ ---------------------------------------------------------------
  //**/ count the number of values that are below the threshold
  //**/ ---------------------------------------------------------------
  nBelow_ = 0;
  for (ss = 0; ss < nSamples_; ss++) 
    if (VecYBins_[ss] == 1.0) nBelow_++;
  if (printLevel > 2)
    printOutTS(PL_INFO,
         "   ====> BinomialAnalyzer: number below threshold = %d (%d)\n",
         nBelow_, nSamples_);

  //**/ ---------------------------------------------------------------
  //**/ create binomial cdf (VecBinomialCDF_[i] - probabilities that 
  //**/     there are more than i values below the threshold given p0)
  //**/ If p0 is smaller, VecBinomialCDF_[i] will be smaller (smaller
  //**/ probabilities that there are more than i values below threshold)
  //**/ ---------------------------------------------------------------
  setupBinomialCDF(nSamples_, p0);
  for (ss = 0; ss <= nSamples_; ss++)
    VecBinomialCDF_[ss] = 1.0 - VecBinomialCDF_[ss];

  //**/ ---------------------------------------------------------------
  //**/ Given the observation that there are nBelow_ values below the
  //**/ threshold, typeI_ is the error of making the assumption that
  //**/ p0 is the probability that an observation is below threshold.
  //**/ (if p0 is smaller, then this error is smaller)
  //**/ ---------------------------------------------------------------
  typeI_ = VecBinomialCDF_[nBelow_];

  return typeI_;
}

// *************************************************************************
// create a binomial cumulative density function
// -------------------------------------------------------------------------
double *BinomialAnalyzer::setupBinomialCDF(int n, double p0)
{
  int    k;
  double nFact, kFact, nkFact;

  nFact = factorial(n);
  VecBinomialCDF_.setLength(n+1);
  for (k = 0; k <= n; k++)
  {
    kFact  = factorial(k);
    nkFact = factorial(n-k);
    VecBinomialCDF_[k] = pow(p0,1.0*k) * pow(1-p0,1.0*(n-k));
    VecBinomialCDF_[k] *= (nFact / (kFact * nkFact));
    if (k > 0) VecBinomialCDF_[k] += VecBinomialCDF_[k-1];
  }
  return 0;
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

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
BinomialAnalyzer& BinomialAnalyzer::operator=(const BinomialAnalyzer &)
{
   printOutTS(PL_ERROR, 
        "BinomialAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// Getters for analysis results
// ------------------------------------------------------------------------
int BinomialAnalyzer::get_nSamples()
{
  return nSamples_;
}
double *BinomialAnalyzer::get_Ybin()
{
  psVector vecT;
  vecT = VecYBins_;
  return vecT.takeDVector();
}
int BinomialAnalyzer::get_nBelow()
{
  return nBelow_;
}
double *BinomialAnalyzer::get_BinomialCDF()
{
  psVector vecT;
  vecT = VecBinomialCDF_;
  return vecT.takeDVector();
}
double BinomialAnalyzer::get_typeI()
{
   return typeI_;
}


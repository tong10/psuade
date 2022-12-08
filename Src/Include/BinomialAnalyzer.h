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
// Definition for the class BinomialAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __BINOMIALANALYZERH__
#define __BINOMIALANALYZERH__

#include "psVector.h"
#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
class BinomialAnalyzer : public Analyzer
{
   int    nSamples_;
   int    nBelow_;
   double typeI_;
   psVector VecYBins_;
   psVector VecBinomialCDF_;

public:

   //**/ Constructor
   BinomialAnalyzer();

   //**/ Destructor
   ~BinomialAnalyzer();

   //**/ Perform analysis
   //**/ @param adata - all data needed for analysis
   double analyze(aData &adata);

   //**/ assign operator
   //**/ @param analyzer
   BinomialAnalyzer& operator=(const BinomialAnalyzer &analyzer);

   //**/ Perform set up a binomial cumulative density function
   //**/ @param n  - number of outcomes
   //**/ @param p0 - probability of outcomes 
   double *setupBinomialCDF(int n, double p0);

   //**/ compute factorial
   //**/ @param n  - input integer
   double factorial(int n);

   /** Getters for analysis results */
   int    get_nSamples();
   double *get_Ybin();
   int    get_nBelow();
   double *get_BinomialCDF();
   double get_typeI();

};

#endif // __BINOMIALANALYZERH__


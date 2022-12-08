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
// Definitions for the class FASTAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#ifndef __FASTANALYZERH__
#define __FASTANALYZERH__

#include "psVector.h"
#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class FASTAnalyzer : public Analyzer
{
private:

   int    nInputs_;
   int    M_;
   double FASTvariance_;
   psVector VecFourierCoefs_;

public:

   //**/ Constructor
   FASTAnalyzer();

   //**/ Destructor
   ~FASTAnalyzer();

   //**/ Perform analysis
   //**/ @param adata - all data needed for analysis
   double analyze(aData &adata);

   //**/ assign operator
   //**/ @param analyzer
   FASTAnalyzer& operator=(const FASTAnalyzer &analyzer);

   //**/ calculate frequencies
   //**/ @param nInputs  - number of inputs
   //**/ @param nSamples - number of sample points
   //**/ @param omegas   - frequencies
   int calculateOmegas(int nInputs, int nSamples, int *omegas);

   //**/ calculate Fourier coefficients
   //**/ @param nSamples    - number of sample points
   //**/ @param nInputs     - number of inputs
   //**/ @param Y           - function values
   //**/ @param coefs       - Fourier coefficients
   //**/ @param outputLevel - print level
   int computeCoefficents(int nSamples, int nInputs, double *Y, 
                          psVector &, int outputLevel);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_M();
   double *get_fourierCoefs();
   double get_FASTvariance();
};

#endif // __FASTANALYZERH__


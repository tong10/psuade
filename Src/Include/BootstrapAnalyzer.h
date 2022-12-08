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
// Definition for the class BootstrapAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#ifndef __BOOTSTRAPANALYZERH__
#define __BOOTSTRAPANALYZERH__

#include "psVector.h"
#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class BootstrapAnalyzer : public Analyzer
{
   int    nSteps_;
   double mean_;
   double stdev_;
   psVector VecStoredVals_;

public:

   //**/ Constructor
   BootstrapAnalyzer();

   //**/ Copy Constructor created by Bill Oliver
   BootstrapAnalyzer(const BootstrapAnalyzer &);

   //**/ Destructor
   ~BootstrapAnalyzer();

   //**/ Perform analysis
   //**/ @param adata - all data needed for analysis
   double analyze(aData &adata);

   //**/ assign operator
   //**/ @param analyzer 
   BootstrapAnalyzer& operator=(const BootstrapAnalyzer &analyzer);

   //**/ Compute sample mean and variance
   //**/ @param nInputs       - number of input variables
   //**/ @param nSamples      - number of sample points
   //**/ @param nOutputs      - number of output variables
   //**/ @param sampleOutputs - sample outputs 
   //**/ @param mean          - pointer to sample mean
   //**/ @param variance      - pointer to sample variance
   //**/ @param outputID      - which output to analyze
   //**/ @param flag          - print flag
   int computeMeanVariance(int nSamples, int nOutputs, 
                           double *sampleOutputs, double *mean,
                           double *variance,int outputID, int flag);

   //**/ Compute bootstrap sample means
   //**/ @param nSamples - number of sample points
   //**/ @param Y        - sample data 
   //**/ @param ntimes   - number of means to compute
   //**/ @param bmeans   - bootstrap means
   int genBootstrap(int nSamples, double *Y, int ntimes, double *bmeans);

   //**/ Compute Jackknife statistics
   //**/ @param nSamples - number of sample points
   //**/ @param Y        - sample data 
   double genJackknife(int nSamples, double *Y);

   //**/ set up normal PDF 
   //**/ @param mean  - mean of the normal distribution
   //**/ @param stdev - standard deviation of the normal distribution
   int setupNormalCDF(double mean, double stdev);

   //**/ perform inverse CDF lookup
   //**/ @param value - probability to be looked up for the variable data
   double normalCDFInv(double &value);

   /** Getters for analysis results */
   int get_nSteps();
   double get_mean();
   double get_stdev();
   double *get_storedValues();

};

#endif // __BOOTSTRAPANALYZERH__


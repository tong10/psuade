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
// Definition for the class MomentAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************

#ifndef __MOMENTANALYZERH__
#define __MOMENTANALYZERH__

#include "psVector.h"
#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class MomentAnalyzer : public Analyzer
{
   int nSamples_;
   int nGroups_;
   int nInputs_;
   int nOutputs_;
   psVector VecStats_;

public:

   //**/ Constructor
   MomentAnalyzer();

   //**/ Destructor
   ~MomentAnalyzer();

   //**/ Perform analysis
   //**/ @param nInps - number of inputs
   //**/ @param nSamp - number of samples
   //**/ @param ilbs  - input lower bounds
   //**/ @param iubs  - input upper bounds
   //**/ @param sInps - sample inputs
   //**/ @param sOuts - sample outputs
   void analyze(int nInps, int nSamp, double *ilbs,
                double *iubs, double *sInps, double *sOuts);

   //**/ Perform analysis
   //**/ @param adata - all data needed for analysis
   double analyze(aData &adata);

   //**/ assign operator
   //**/ @param analyzer
   MomentAnalyzer& operator=(const MomentAnalyzer &analyzer);

   //**/ Perform analysis on merging 2 groups into 1
   //**/ @param nInputs - number of input variables
   //**/ @param nOutputs - number of output variables
   //**/ @param nSamples - number of sample points
   //**/ @param nSubSamples - number of sample points in a group
   //**/ @param sampleOutputs - sample outputs 
   //**/ @param outputID - which output to analyze
   int analyzeMore(int nInputs, int nOutputs, int nSamples, 
                   int nSubSamples, double *sampleOutputs, int outputID);

   //**/ Compute sample mean and variance
   //**/ @param nInputs - number of input variables
   //**/ @param nOutputs - number of output variables
   //**/ @param nSamples - number of sample points
   //**/ @param sampleOutputs - sample outputs 
   //**/ @param mean - pointer to sample mean
   //**/ @param variance - pointer to sample variance
   //**/ @param outputID - which output to analyze
   int computeMeanVariance(int nInputs, int nOutputs, int nSamples,
                           double * sampleOutputs, double *mean,
                           double *variance,int outputID);

   //**/ Compute sample skewness and kurtosis
   //**/ @param nInputs - number of input variables
   //**/ @param nOutputs - number of output variables
   //**/ @param nSamples - number of sample points
   //**/ @param sampleOutputs - sample outputs 
   //**/ @param stdev - pointer to sample standard deviation
   //**/ @param skewness - pointer to sample skewness
   //**/ @param kurtosis - pointer to sample kurtosis
   //**/ @param outputID - which output to analyze
   int computeSkewnessKurtosis(int nInputs, int nOutputs, int nSamples,
                      double * sampleOutputs, double stdev,
                      double *skewness, double *kurtosis,int outputID);

   //**/ create matlab/scilab distribution plot
   //**/ @param fname    - output file name
   //**/ @param nInputs  - number of input variables
   //**/ @param nOutputs - number of output variables
   //**/ @param nSamples - number of sample points
   //**/ @param outID    - which output to analyze
   //**/ @param outVals  - output values
   //**/ @param outNames - output names
   int genPlot(char *fname, int nInputs, int nOutputs, int nSamples, 
               int outID, double *Y, char **outNames);

   /** Getters for analysis results */
   int get_nSamples();
   int get_nGroups();
   int get_nInputs();
   int get_nOutputs();
   double get_mean();
   double get_stdev();
   double get_skewness();
   double get_kurtosis();
};

#endif // __MOMENTANALYZERH__


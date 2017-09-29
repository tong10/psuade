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
// Definition for the class MCMCAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#ifndef __MCMCANALYZERH__
#define __MCMCANALYZERH__

#include "Analyzer.h"
#include "pData.h"

// ************************************************************************
// class definition
// ************************************************************************
class MCMCAnalyzer : public Analyzer
{
private:
   int    mode_;
   int    nInputs_;
   int    nOutputs_;
   double *means_;           //length is nInputs_
   double *sigmas_;          //length is nInputs_
   double *mostLikelyInput_; //length is nInputs_
   double *mostLikelyOutput_;

public:

   MCMCAnalyzer();

   ~MCMCAnalyzer();

   double analyze(aData &adata);

   MCMCAnalyzer& operator=(const MCMCAnalyzer &analyzer);

   double genMatlabFile(int nInputs, double *lower, double *upper,
                        double *ranges, int nPlots, int *plotIndices,
                        int nbins, int **bins, int ****bins2, 
                        pData &pData, int nChains, int chainCnt,
                        double ***XChains, int *chainStatus);

   int setParams(int nParams, char **params);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_nOutputs();
   double *get_means();
   double *get_sigmas();
   double *get_mostLikelyInput();
   double *get_mostLikelyOutput();
};

#endif // __MCMCANALYZERH__


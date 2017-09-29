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
#include "FuncApprox.h"

// ************************************************************************
// class definition
// ************************************************************************
class MCMCAnalyzer : public Analyzer
{
private:
   int    mode_;
   int    bfmode_;
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

   double analyze_bf(aData &adata);

   MCMCAnalyzer& operator=(const MCMCAnalyzer &analyzer);

   double genMatlabFile(int nInputs, double *lower, double *upper,
                        double *ranges, int nPlots, int *plotIndices,
                        int nbins, int **pbins, int ****pbins2, 
                        int **bins, int ****bins2, 
                        pData &pData, int nChains, int chainCnt,
                        double ***XChains, int *chainStatus, 
                        double *Xmax, double Ymin);

   int    genPostLikelihood(int nInputs, double *lower,
                 double *upper, double *XRange, int numChains,
                 int chainCnt, double ***XChains, int *chainStatus,
                 int chainLimit, int *rsIndices, double *rsValues,
                 int *designParams, int dnInputs, int dnSample, 
                 double *dSamInputs, FuncApprox **faPtrs, 
                 FuncApprox **faPtrs1, int nOutputs, double *discOutputs,
                 double *discFuncConstantMeans, double *dSamMeans,
                 double *dSamStdevs);

   double readSpecFile(int nInputs, int nOutputs, int *dnSamp,
                       int *dnInps, int **dParams, double **dSamIns,
                       double **dMeans, double **dStds, int printLevel);

   int setParams(int nParams, char **params);

   int checkConvergence(int num, double *means, double *stds, int leng);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_nOutputs();
   double *get_means();
   double *get_sigmas();
   double *get_mostLikelyInput();
   double *get_mostLikelyOutput();
};

#endif // __MCMCANALYZERH__


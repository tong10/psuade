// ************************************************************************
// Copyright (c) 2015   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// DASSI is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// DASSI is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// KPCA functions 
// DATE   : 2015
// ************************************************************************
#ifndef __MetropolisHastingsMCMC__
#define __MetropolisHastingsMCMC__

#include "psMatrix.h"
#include "psVector.h"
#include "KPCA.h"

class MetropolisHastingsMCMC : public KPCA 
{
protected : 
  int  burnInSize_, MCMCMaxSamples_, MCMCNDim_, MCMCRandSeed_;
  int  breakPtFlag_, likelihoodEvaluationCount_;
  char PosteriorFile_[10000];
  char likelihoodFunction_[10000];
  psMatrix MatMcmcChain_;
  psVector VecIG_;

public :

  MetropolisHastingsMCMC();
  void initialize(int nInps, int rand_seed,int MCMC_nSamples);
  psVector genSampleFromPrior(psVector &x_current,double prop_scale);
  int  setInitialGuess(psVector vecIG);
  int  setPosteriorFile(char *);
  int  setLikelihoodFunction(char *);
  int  setBreakPoint();
  void runMCMC();

private:
  double computeLikelihood(psVector &vecXi,double scale);
  double computePosterior(psVector &eta, double scale);
};

#endif 


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
// Definitions for the class DNN
// AUTHOR : Charles Tong
// DATE   : 2019
// ************************************************************************

#ifndef __DNNH__
#define __DNNH__

#include <stdio.h>
#include "psVector.h"
#include "FuncApprox.h"

// ************************************************************************
// class definition
// ************************************************************************
class DNN : public FuncApprox
{
  int nLevels_; /* level 0: input layer; level nLevels_: output layer */
  psIVector VecActivationFcns_;
  psIVector VecNumNodes_;
  int batchSize_;
  psMatrix **MatWs_;
  psVector **VecBs_;
  psMatrix **MatZs_;
  psMatrix **MatAs_;
  psMatrix **MatdWs_;
  psMatrix **MatdZs_;
  psVector **VecdBs_;
  int      optScheme_;/* optimization method */
  double   alpha0_;   /* alpha0 in adaptive learning rate */
  double   alpha_;
  double   beta_;     /* beta in momentum equation */
  double   beta1_;    /* beta1 in momentum (Adam optimization */
  double   beta2_;    /* beta2 in RMSprop (Adam optimization) */
  double   epsilon_;  /* for Adam */
  int      regularizationOn_;
  int      maxIter_;
  psVector VecDropOutRatios_;
  
public:

   /** constructor */
   DNN(int nInputs, int nSamples);

   /** destructor */
   ~DNN();

   int initialize(double *, double *);
   int genNDGridData(double *, double *, int *, double **, double **);
   int gen1DGridData(double *, double *, int, double *, int *, 
                     double **, double **);
   int gen2DGridData(double *, double *, int, int, double *, int *, 
                     double **, double **);
   int gen3DGridData(double *, double *, int, int, int, double *, int *, 
                     double **, double **);
   int gen4DGridData(double *, double *, int, int, int, int, double *, 
                     int *, double **, double **);
   double evaluatePoint(double *);
   double evaluatePoint(int, double *, double *);
   double evaluatePointFuzzy(double *, double &);
   double evaluatePointFuzzy(int, double *, double *, double *);

   int  train(psVector, psVector);
   int  optimizeLBFGS(psVector, psVector);
   int  optimizeGradDescent(psVector, psVector);
   int  propagateForward(psVector);
   int  propagateBackward(psVector);
   void showDetails();
};

#endif

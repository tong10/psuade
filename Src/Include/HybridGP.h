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
// Definition for the class Hybrid homogeneous GP
// AUTHOR : CHARLES TONG
// DATE   : 2019
// ************************************************************************
#ifndef __HYBRIDGPH__
#define __HYBRIDGPH__

#include <stdio.h>
#include "FuncApprox.h"
#include "psMatrix.h"
#include "psVector.h"
#include "HomLegendreRegression.h"

// ************************************************************************
// class definition
// ************************************************************************
class HybridGP : public FuncApprox 
{
public:
   int      polyOrder_;
   int      sampleIncrement_;
   int      hasQuantiles_;
   double   expPower_;
   psVector VecHyperPs_;
   psVector VecCInvY_;
   psVector VecLegendreCoefs_;
   psVector VecLegendreStds_;
   psVector VecXDataN_;
   psVector VecYDataN_;
   psVector VecXDistances_;
   psVector VecXLows_;
   psVector VecXHighs_;
   psMatrix CMatrix_;
   HomLegendreRegression *homRegressor_;

public:

   HybridGP(int, int);
   ~HybridGP();

   int    initialize(double*,double*);
   int    genNDGridData(double*,double*,int*,double**,double**);
   int    gen1DGridData(double*,double *,int,double*, 
                        int *, double **, double **);
   int    gen2DGridData(double*,double *,int,int,double*, 
                        int *, double **, double **);
   int    gen3DGridData(double*,double *,int,int,int,double*, 
                        int *, double **, double **);
   int    gen4DGridData(double*,double *,int,int,int,int,double*, 
                        int *, double **, double **);

   double evaluatePoint(double *);
   double evaluatePoint(int, double *, double *);
   double evaluatePointFuzzy(double *, double &);
   double evaluatePointFuzzy(int, double *, double *, double *);
   double computeGradients(psVector, psVector &, int &);

private:
   int    train();
   int    computeDistances();
   int    interpolate(int, double *, double *, double *);
   void   constructCMatrix(psMatrix &, psVector);
};

#endif // __HYBRIDGPH__


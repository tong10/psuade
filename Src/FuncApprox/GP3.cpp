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
// Functions for the class GP3
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#ifdef WINDOWS
#include <windows.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "GP3.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "Sampling.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// ************************************************************************
// external functions
// ------------------------------------------------------------------------
extern "C" {
  void dpotrf_(char *, int *, double *, int *, int *);
  void dpotrs_(char *, int *, int *, double *, int *, double *,int *,int *);
#ifdef HAVE_LBFGS
#include "../../External/L-BFGS-B-C/src/lbfgsb.h"
#endif
}

// ************************************************************************
// Constructor for object class GP3
// ------------------------------------------------------------------------
GP3::GP3(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  optLinTerm_ = 0;
  faID_ = PSUADE_RS_GP2;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GP3::~GP3()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int GP3::initialize(double *XIn, double *YIn)
{
  double *dVec;
  XDataN_.setLength(nSamples_*nInputs_);
  dVec = XDataN_.getDVector();
  initInputScaling(XIn, dVec, 1);
  YData_.setLength(nSamples_);
  dVec = YData_.getDVector();
  initOutputScaling(YIn, dVec);
   
  if (outputLevel_ > 1) printf("GP3 training begins....\n");
  train();
  if (outputLevel_ > 1) printf("GP3 training completed.\n");
  if (psRSCodeGen_ == 1) genCode();
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int GP3::genNDGridData(double *XIn, double *YIn, int *N, double **X2, 
                      double **Y2)
{
  int    totPts;
  double *XX, *YY;

  initialize(XIn, YIn);
  if ((*N) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
  genNDGrid(N, &XX);
  if ((*N) == 0) return 0;
  totPts = (*N);

  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  YY = new double[totPts];
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*X2) = XX;
  (*Y2) = YY;
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int GP3::gen1DGridData(double *XIn, double *YIn, int ind1,
                       double *settings, int *n, double **X2, double **Y2)
{
  int    totPts, ii, kk;
  double *XX, *YY, HX;

  initialize(XIn, YIn);
  if ((*n) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*X2) = new double[totPts];
  XX = new double[totPts*nInputs_];
  checkAllocate(XX, "XX in GP3::gen1DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
     for (kk = 0; kk < nInputs_; kk++) 
        XX[ii*nInputs_+kk] = settings[kk]; 
     XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
     (*X2)[ii] = HX * ii + lowerBounds_[ind1];
  }
   
  YY = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int GP3::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                       double *settings, int *n, double **X2, double **Y2)
{
  int    totPts, ii, jj, kk, index;
  double *XX, *YY, *HX;

  initialize(XIn, YIn);
  if ((*n) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[2];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  XX = new double[totPts*nInputs_];
  (*X2) = new double[2*totPts];
  checkAllocate(*X2, "X2 in GP3::gen2DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      index = ii * nPtsPerDim_ + jj;
      for (kk = 0; kk < nInputs_; kk++) 
        XX[index*nInputs_+kk] = settings[kk]; 
      XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
      XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
      (*X2)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
      (*X2)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
    }
  }
    
  YY = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int GP3::gen3DGridData(double *XIn, double *YIn,int ind1,int ind2,int ind3,
                       double *settings, int *n, double **X2, double **Y2)
{
  int    totPts, ii, jj, kk, ll, index;
  double *XX, *YY, *HX;

  initialize(XIn, YIn);
  if ((*n) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[3];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  XX = new double[totPts*nInputs_];
  (*X2) = new double[3*totPts];
  checkAllocate(*X2, "X2 in GP3::gen3DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        for (kk = 0; kk < nInputs_; kk++) 
          XX[index*nInputs_+kk] = settings[kk]; 
        XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
        XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
        XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
        (*X2)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
        (*X2)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
        (*X2)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
      }
    }
  }
    
  YY = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int GP3::gen4DGridData(double *XIn,double *YIn,int ind1,int ind2, int ind3,
                       int ind4, double *settings, int *n, double **X2, 
                       double **Y2)
{
  int    totPts, ii, jj, kk, ll, mm, index;
  double *XX, *YY, *HX;

  initialize(XIn, YIn);
  if ((*n) == -999 || X2 == NULL || Y2 == NULL) return 0;
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[4];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  XX = new double[totPts*nInputs_];
  (*X2) = new double[4*totPts];
  checkAllocate(*X2, "X2 in GP3::gen4DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        for (mm = 0; mm < nPtsPerDim_; mm++) 
        {
          index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                  jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
          for (kk = 0; kk < nInputs_; kk++) 
            XX[index*nInputs_+kk] = settings[kk]; 
          XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
          XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
          XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
          XX[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
          (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
          (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
          (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
          (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
        }
      }
    }
  }
    
  YY = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double GP3::evaluatePoint(double *X)
{
  double Y=0.0;
  int    ii, iOne=1;
  double *XX;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  XX = new double[nInputs_];
  checkAllocate(XX, "XX in GP3::evaluatePoint");
  for (ii = 0; ii < nInputs_; ii++)
    XX[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(iOne, XX, &Y, NULL);
  Y = Y * YStd_ + YMean_;
  delete [] XX;
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double GP3::evaluatePoint(int npts, double *X, double *Y)
{
  int    ii, jj;
  double *XX;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  XX = new double[npts*nInputs_];
  checkAllocate(XX, "XX in GP3::evaluatePoint");
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      XX[jj*nInputs_+ii] = (X[jj*nInputs_+ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(npts, XX, Y, NULL);
  for (jj = 0; jj < npts; jj++)
     Y[jj] = Y[jj] * YStd_ + YMean_;
  delete [] XX;
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double GP3::evaluatePointFuzzy(double *X, double &std)
{
  int    ii, iOne=1;
  double *XX, Y=0.0;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  XX = new double[nInputs_];
  checkAllocate(XX, "XX in GP3::evaluatePointFuzzy");
  for (ii = 0; ii < nInputs_; ii++)
    XX[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(iOne, XX, &Y, &std);
  Y = Y * YStd_ + YMean_;
  if (std < 0) printf("GP3 ERROR: variance (%e) < 0\n",std);
  else         std = sqrt(std) * YStd_;
  delete [] XX;
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double GP3::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystds)
{
  int    ii, jj;
  double *XX;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  XX = new double[npts*nInputs_];
  checkAllocate(XX, "XX in GP3::evaluatePointFuzzy");
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      XX[jj*nInputs_+ii] = (X[jj*nInputs_+ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(npts, XX, Y, Ystds);
  for (int ii = 0; ii < npts; ii++)
  {
    Y[ii] = Y[ii] * YStd_ + YMean_;
    if (Ystds[ii] < 0) printf("GP3 ERROR: variance (%e) < 0\n", Ystds[ii]);
    else               Ystds[ii] = sqrt(Ystds[ii]) * YStd_;
  }
  delete [] XX;
  return 0.0;
}

// ************************************************************************
// train to get lengthScales and sigF2 and sigN2
// ------------------------------------------------------------------------
int GP3::train()
{
  int     ii, jj, kk, its, count, status, newnInps, nhypers;
  int     nSamp=1000, *SS, fail, nHist, maxHist=10;
  double  minVal=1e35, dmax, dmin, dtemp, *work;
  double  *PVals, *lBounds, *uBounds, *Grads, FValue, *XS, *YS, *history;
  double  sqrt2pi=sqrt(2*3.1415928);
  Sampling *sampler;
  psVector YT;

  nhypers = nInputs_ + 4;
  if (optLinTerm_) nhypers += nInputs_;

  hyperparameters_.setLength(nhypers);
  for (ii = 0; ii < nInputs_; ii++)
  {
    dmax = -PSUADE_UNDEFINED;
    dmin =  PSUADE_UNDEFINED;
    for (kk = 0; kk < nSamples_; kk++)
    {
      if (XDataN_[kk*nInputs_+ii] > dmax) dmax = XDataN_[kk*nInputs_+ii];
      if (XDataN_[kk*nInputs_+ii] < dmin) 
        dmin = XDataN_[kk*nInputs_+ii];
    }
    hyperparameters_[ii] = 0.25 * (dmax - dmin);
    if (hyperparameters_[ii] == 0)
    {
      printf("GP3 ERROR: Input %d is a constant.\n", ii+1);
      printf("           Prune this input first (use idelete).\n");
      exit(1);
    }
    hyperparameters_[ii] = log(hyperparameters_[ii]);
  }
  dtemp = 0.0;
  for (kk = 0; kk < nSamples_; kk++) 
    if (PABS(YData_[kk]) > dtemp) dtemp = PABS(YData_[kk]);
  hyperparameters_[nInputs_] = dtemp;
    
  hyperparameters_[nInputs_+1] = log(1e-12);
  hyperparameters_[nInputs_+2] = 0;
  hyperparameters_[nInputs_+3] = log(1e-10);
  count = nInputs_ + 4;
  if (optLinTerm_) 
  {
    for (ii = 0; ii < nInputs_; ii++) 
      hyperparameters_[count++] = 0.0;
  }
  if (outputLevel_ > 0) 
  {
    printf("GP3 initial hyperparameter values:\n");
    for (ii = 0; ii < hyperparameters_.length(); ii++) 
      printf("   Initial hyperparameters %2d  = %e\n",ii+1,
             hyperparameters_[ii]);
  }

  int  optimizeFlag=1;
  char pString[1000], winput[1000];
  if (psRSExpertMode_ == 1 && isScreenDumpModeOn() == 1)
  {
    sprintf(pString, "Use user-provided hyperparameters? (y or n) ");
    getString(pString, winput);
    if (winput[0] == 'y')
    {
      for (ii = 0; ii < hyperparameters_.length(); ii++)
      {
        sprintf(pString,"Enter hyperparameter %d : ",ii+1);
        hyperparameters_[ii] = getDouble(pString);
        if (hyperparameters_[ii] <= 0.0)
          printf("WARNING: hyperparameter <= 0 may not be valid.\n");
      }
      optimizeFlag = 0;
    }
  }

  status = computeDistances();
  if (status != 0)
  {
    printf("GP3 INFO: since there are repeated sample points.\n");
    printf("    GP will continue without optimizing the length\n");
    printf("    scales but instead set them to 1.\n");
    printf("    This may not give a good quality response surface.\n");
    printf("    So if you want to prune the sample and do it again,\n");
    printf("    terminate now. I am giving you 30 seconds to decide.\n");
    optimizeFlag = 0;
    for (ii = 0; ii < nInputs_; ii++) hyperparameters_[ii] = 1.0;
    for (ii = nInputs_; ii < hyperparameters_.length(); ii++)
      hyperparameters_[ii] = 0.0;
#ifdef WINDOWS
    Sleep(20000);
#else
    sleep(20);
#endif
    printf("    10 more seconds.\n");
#ifdef WINDOWS
    Sleep(10000);
#else
    sleep(10);
#endif
  }

  XLows_.setLength(hyperparameters_.length());
  XHighs_.setLength(hyperparameters_.length());
  for (ii = 0; ii < hyperparameters_.length(); ii++) 
  {
    XLows_[ii]  = +PSUADE_UNDEFINED;
    XHighs_[ii] = -PSUADE_UNDEFINED;
  }

  if (optimizeFlag == 1)
  {
    nhypers = hyperparameters_.length();
    PVals = new double[nhypers];
    lBounds = new double[nhypers];
    uBounds = new double[nhypers];
    for (ii = 0; ii < nInputs_; ii++) 
    {
      lBounds[ii] = -4.0;
      uBounds[ii] =  4.0;
    }
    lBounds[nInputs_] = -3.0;
    uBounds[nInputs_] = 10.0;
    lBounds[nInputs_+1] = -28;
    uBounds[nInputs_+1] = - 3;
    lBounds[nInputs_+2] = -0.5;
    uBounds[nInputs_+2] =  0.5;
    lBounds[nInputs_+3] = -24;
    uBounds[nInputs_+3] = - 2;
    count = nInputs_ + 4;
    if (optLinTerm_) 
    {
      for (ii = 0; ii < nInputs_; ii++) 
      {
        lBounds[count] = - 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
        uBounds[count++] = + 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
      }
    }
    if (outputLevel_ > 0) 
    {
      for (ii = 0; ii < nhypers; ii++)
        printf("Hyperparameter bounds %3d = %12.4e %12.4e\n",ii+1,
               lBounds[ii],uBounds[ii]);
    }
    nSamp = 100;
    if (nhypers > 51)
      sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else
      sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nhypers, lBounds, uBounds);
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 0);
    sampler->initialize(0);
    SS = new int[nSamp];
    XS = new double[nSamp*nhypers];
    YS = new double[nSamp];
    sampler->getSamples(nSamp, nhypers, 1, XS, YS, SS);
    delete [] SS;
    delete sampler;
    for (ii = 0; ii < nhypers; ii++) XS[ii] = hyperparameters_[ii];
 
    integer nInps, iprint=0, itask, *task=&itask, lsave[4], isave[44];
    integer *iwork, nCorr=5, *nbds, csave[60], nLBFGS=1;
    double  factr, pgtol, dsave[29];

    nInps = nhypers;
    Grads = new double[nInps];
    for (ii = 0; ii < nInps; ii++) Grads[ii] = 0.0;
    nbds = new integer[nInps];
    for (ii = 0; ii < nInps; ii++) nbds[ii] = 2;
    factr = 1e7;
    pgtol = 1e-8;
    kk = (2 * nCorr + 5) * nInps + 11 * nCorr * nCorr + 8 * nCorr;
    work  = new double[kk];
    iwork = new integer[3*nInps];
    its   = 0;
    history = new double[maxHist];

    for (ii = 0; ii < nLBFGS; ii++)
    {
      if (outputLevel_ > 0 && nLBFGS > 1 && ii > 0) 
        printf("GP3 LBFGS sample %d (%ld)\n",ii+1, nLBFGS);
      for (jj = 0; jj < nInps; jj++) PVals[jj] = XS[ii*nInps+jj];
      if (outputLevel_ > 0 && ii > 0) 
        for (jj = 0; jj < nInps; jj++) 
          printf("   Hyperparameter (init) %2d = %e\n",jj+1,PVals[jj]);
      nHist = 0;
      fail = 0;
      its = 0;
      *task = (integer) START;
      while (1 && its < 1000)
      {
        its++;
        setulb(&nInps, &nCorr, PVals, lBounds, uBounds, nbds, &FValue,
              Grads, &factr, &pgtol, work, iwork, task, &iprint, csave,
              lsave, isave, dsave);
        if (IS_FG(*task))
        {
          if (outputLevel_ > 2) 
            for (kk = 0; kk < nInps; kk++)
              printf("     Hyperparameter %5d = %e\n",kk+1,PVals[kk]);

#if 1
          FValue = computeGradients(PVals, Grads, status);
          if (FValue == 1e35) 
          {
            fail = 1;
            break;
          }
          dtemp = 0.0;
          for (kk = 0; kk < nInps; kk++)
            dtemp += pow(Grads[kk], 2.0);
          dtemp = sqrt(dtemp);
#else
          status = 0;
          FValue = computeLikelihood(PVals);
          if (FValue == 1e35) 
          {
            fail = 1;
            break;
          }
          dtemp = 0.0;
          for (kk = 0; kk < nInps; kk++)
          {
            PVals[kk] += 1.0e-11;
            Grads[kk] = computeLikelihood(PVals);
            Grads[kk] = 1e11 * (Grads[kk] - FValue);
            PVals[kk] -= 1.0e-11;
            dtemp += pow(Grads[kk], 2.0);
          }
          dtemp = sqrt(dtemp);
#endif
          if (outputLevel_ > 1) 
            printf("   LBFGS Current FValue = %e (its=%d, |grad|=%e)\n",
                   FValue,its,dtemp);
          if (outputLevel_ > 2) 
            for (kk = 0; kk < nInps; kk++)
              printf("         Gradient %5d = %e\n",kk+1,Grads[kk]);
          if (FValue < minVal && status == 0)
          {
            for (jj = 0; jj < nhypers; jj++) hyperparameters_[jj] = PVals[jj];
            minVal = FValue;
          }
          if (nHist < maxHist && FValue != 1e35) history[nHist++] = FValue;
          else if (FValue != 1e35)
          {
            for (kk = 0; kk < 4; kk++) history[kk] = history[kk+1];
            history[4] = FValue;
          }
          if (nHist >= maxHist && FValue != 1e35)
          {
            dmin = 0.0;
            for (kk = 0; kk < maxHist; kk++) dmin += history[kk];
            dmin *= 0.2;
            dmax = 0.0;
            for (kk = 0; kk < maxHist; kk++) dmax += pow(history[kk]-dmin, 2.0);
            dmax = sqrt(dmax/9.0);
            if (dmin == 0) dmin = 1;
            dmax = dmin = 1;
            if (PABS(dmax/dmin) < 1e-6)
            {
              *task = STOP_ITER;
              if (outputLevel_ > 1) 
              {
                printf("INFO: PSUADE issues a stop (converged?))\n");
                printf("INFO: current iteration   = %d\n", its);
                printf("INFO: current best FValue = %e (%e)\n",minVal,
                       PABS(dmax/dmin));
                printf("INFO: current grad norm   = %e\n",dtemp);
              }
            }
          }
          if (isave[33] >= 1000) 
          {
            *task = STOP_ITER;
            if (outputLevel_ > 1) 
            {
              printf("INFO: PSUADE issues a stop (> 1000 iterations)\n");
              printf("INFO: current best FValue = %e\n", minVal);
            }
          }
        }
        else if (*task == NEW_X)
        {
          if (isave[33] >= 1000) 
          {
            *task = STOP_ITER;
            if (outputLevel_ > 1) 
              printf("INFO: PSUADE issues a stop (> 1000 iterations)\n");
          }
        }
        else
        {
          if (outputLevel_ > 1) 
          {
            printf("INFO: LBFGS issues a stop (its = %d)\n",its);
            printf("INFO: current best FValue = %e\n", minVal);
          }
          break;
        }
      }
      if (outputLevel_ > 1) 
      {
        for (jj = 0; jj < hyperparameters_.length(); jj++) 
          printf("   Hyperparameter (final) %2d = %e\n",jj+1,
                 hyperparameters_[jj]);
        printf("   Sample %5d: Fvalue (final) = %e\n", ii+1, minVal);
      }
      if (ii == (nLBFGS-1) && (fail == 1)) nLBFGS++;
    }
    delete [] work;
    delete [] iwork;
    delete [] nbds;
    delete [] Grads;
    delete [] history;
    delete [] PVals;
    delete [] lBounds;
    delete [] uBounds;
    delete [] XS;
    delete [] YS;

    if (outputLevel_ > 0) 
    {
      printf("GP3 final hyperparameter values:\n");
      for (ii = 0; ii < hyperparameters_.length(); ii++) 
        printf("   Hyperparameter %2d = %e\n",ii+1,
               hyperparameters_[ii]);
    }
  }

  constructCMatrix(CMatrix_, hyperparameters_.getDVector());
  CMatrix_.LUDecompose();
  CInvY_.setLength(nSamples_);
  YT.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    YT[jj] = YData_[jj] - hyperparameters_[nInputs_+2];
  CMatrix_.LUSolve(YT, CInvY_);

  //Cannot do because inverse seems to have problems 3/2017
  //double likelihood = 0.0;
  //for (jj = 0; jj < nSamples_; jj++) likelihood += CInvY_[jj] * YT[jj];
  //likelihood *= 0.5;
  //for (jj = 0; jj < nSamples_; jj++) 
  //{
  //  dtemp = CMatrix.getEntry(jj, jj);
  //  likelihood += log(dtemp*sqrt2pi+1e-20);
  //}
  //if (outputLevel_ > 0) 
  //  printf("GP3 Best likelihood = %e\n", likelihood);
  if (psGMMode_ == 1)
  { 
    psMatrix CMatrix;
    constructCMatrix(CMatrix, hyperparameters_.getDVector());
    FILE *fp = fopen("Cmat.m", "w");
    fprintf(fp, "A = [\n");
    for (ii = 0; ii < nSamples_; ii++) 
    {
      for (jj = 0; jj < nSamples_; jj++) 
        fprintf(fp, "%16.8e ", CMatrix.getEntry(ii,jj));
      fprintf(fp, "\n");
    }
    fprintf(fp,"];\n");
    fclose(fp);
    printf("*** Final Cmatrix given in Cmat.m\n");
  }
  if (outputLevel_ > 1 && optimizeFlag == 1)
  {
    printf("Optimization summary: \n");
    for (ii = 0; ii < nInputs_+4; ii++) 
      printf("Input %2d search bounds = %12.4e %12.4e\n",ii+1,
             XLows_[ii], XHighs_[ii]);
  }
  return 0;
}

// ************************************************************************
// compute log likelihood value
// ------------------------------------------------------------------------
double GP3::computeLikelihood(double *params)
{
  int    jj, status;
  double likelihood, sqrt2pi=sqrt(2*3.1415928), ddata;
  psMatrix CMatrix;
  psVector YVec, YT;

  if (outputLevel_ > 3) 
  {
    printf("computeLikelihood:\n");
    for (jj = 0; jj < hyperparameters_.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, params[jj]);
  }
  constructCMatrix(CMatrix, params);

  status = CMatrix.CholDecompose();
  if (status != 0)
  {
    printf("GP3 ERROR (1): failed in linear system factorization.\n");
    for (jj = 0; jj < hyperparameters_.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, params[jj]);
    return 1e35;
  }
  YVec.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    YVec[jj] = YData_[jj] - params[nInputs_+2];
  YT.setLength(nSamples_);
  status = CMatrix.CholSolve(YVec, YT);
  if (status != 0)
  {
    printf("GP3 ERROR (1): failed in linear system solve.\n");
    return 1e35;
  }

  likelihood = 0.0;
  for (jj = 0; jj < nSamples_; jj++) 
    likelihood += YT[jj] * YVec[jj];
  likelihood *= 0.5;
  for (jj = 0; jj < nSamples_; jj++) 
  {
    ddata = CMatrix.getEntry(jj, jj);
    likelihood += log(ddata*sqrt2pi+1e-50);
  }
  if (outputLevel_ > 3) 
    printf("   computeLikelihood: Likelihood = %e\n", likelihood);
  return likelihood;
}

// ************************************************************************
// compute gradients of log likelihood 
// ------------------------------------------------------------------------
double GP3::computeGradients(double *params, double *grads, int &retstat)
{
  int    jj, kk, ii, ll, status, count;
  double ddata, *TMat,likelihood,sqrt2pi=sqrt(2*3.1415928),dist;
  double fudge;
  psMatrix CMatrix, CInverse, CPrime, CT;
  psVector YVec, CInvY, YT;

  if (outputLevel_ > 2) 
  {
    printf("computeGradients:\n");
    for (jj = 0; jj < hyperparameters_.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, params[jj]);
  }
  constructCMatrix(CMatrix, params);

  for (ii = 0; ii < hyperparameters_.length(); ii++)
  {
    if (params[ii] > XHighs_[ii]) XHighs_[ii] = params[ii]; 
    if (params[ii] < XLows_[ii])  XLows_[ii]  = params[ii]; 
  }

  status = CMatrix.computeInverse(CInverse);
  if (status != 0)
  {
    printf("GP3 ERROR (1): failed in matrix factorization (%d).\n",
           status);
    for (jj = 0; jj < hyperparameters_.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, params[jj]);
    for (jj = 0; jj < hyperparameters_.length(); jj++) grads[jj] = 100;
    return 1e35;
  }
  YVec.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    YVec[jj] = YData_[jj] - params[nInputs_+2];
  CInvY.setLength(nSamples_);
  CInverse.matvec(YVec, CInvY, 0);

  YT.setLength(nSamples_);
  TMat = new double[nSamples_*nSamples_];
  for (ll = 0; ll < nInputs_; ll++)
  {
    count = 0;
    for (jj = 0; jj < nSamples_; jj++)
    {
      TMat[jj*nSamples_+jj] = 0;
      for (kk = jj+1; kk < nSamples_; kk++)
      {
        dist = 0.0;
        for (ii = 0; ii < nInputs_; ii++)
          dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                  exp(2.0*params[ii]);
        dist *= 0.5;
        dist = exp(params[nInputs_]) * exp(-dist);
        ddata = pow(XDistances_[count*nInputs_+ll]/
                exp(2.0*params[ll]),2.0);
        ddata = ddata * exp(2.0*params[ll]);
        dist *= ddata;
        TMat[jj*nSamples_+kk] = dist;
        TMat[kk*nSamples_+jj] = dist;
        count++;
      }
    }
    CPrime.load(nSamples_, nSamples_, TMat);
    CPrime.matvec(CInvY, YT, 0);
    ddata = 0.0;
    for (jj = 0; jj < nSamples_; jj++) ddata += YT[jj] * CInvY[jj];
    grads[ll] = - ddata;
    ddata = 0.0;
    for (ii = 0; ii < nSamples_; ii++)
      for (jj = 0; jj < nSamples_; jj++) 
        ddata += CPrime.getEntry(ii, jj) * CInverse.getEntry(ii, jj);
    grads[ll] += ddata;
    grads[ll] *= 0.5;
  }

  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    TMat[jj*nSamples_+jj] = 1;
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                exp(2.0*params[ii]);
      dist = exp(-0.5 * dist);
      if (dist < 1.0e-50) dist = 0;
      TMat[jj*nSamples_+kk] = dist;
      TMat[kk*nSamples_+jj] = dist;
      count++;
    }
  }
  CPrime.load(nSamples_, nSamples_, TMat);
  CPrime.matvec(CInvY, YT, 0);
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += YT[jj] * CInvY[jj];
  grads[nInputs_] = - ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) 
    for (jj = 0; jj < nSamples_; jj++) 
      ddata += CInverse.getEntry(ii, jj) * CPrime.getEntry(ii, jj);
  grads[nInputs_] += ddata;
  grads[nInputs_] *= 0.5 * exp(params[nInputs_]);

  fudge = sqrt(1.0 * nSamples_);
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj];
  grads[nInputs_+1] = - ddata * ddata - (fudge - 1.0) * ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++)
    for (jj = 0; jj < nSamples_; jj++) 
      ddata += CInverse.getEntry(ii, jj);
  grads[nInputs_+1] += ddata;
  ddata = exp(params[nInputs_+1]);
  grads[nInputs_+1] *= 0.5 * ddata;

  //*/ compute: - 1' * C' * C^{-1} (Y - mu)
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj];
  grads[nInputs_+2] = - ddata;

  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj] * CInvY[jj];
  grads[nInputs_+3] = - ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) 
    ddata += CInverse.getEntry(ii, ii);
  grads[nInputs_+3] += ddata;
  grads[nInputs_+3] *= 0.5 * exp(params[nInputs_+3]);

  int errCount=0;
  CMatrix.eigenSolve(CT, YT, 1);
  likelihood = 0.0;
  for (jj = 0; jj < nSamples_; jj++) 
  {
    ddata = YT[jj];
    if (ddata <= 0) errCount++;
    if (ddata != 0)
      likelihood += log(sqrt(PABS(ddata))*sqrt2pi+1e-50);
  }
  retstat = 0;
  if (errCount > 0)
  {
    retstat = 1;
    if (outputLevel_ > 0)
    {
      printf("GP3 WARNING: CMatrix non-positive definite (%d)\n",
             errCount);
      if (outputLevel_ > 2)
      {
        for (jj = 0; jj < nInputs_; jj++)
          printf("   hyperparameter %2d = %e\n",jj+1,params[jj]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+1,params[nInputs_]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+2,params[nInputs_+1]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+3,params[nInputs_+2]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+4,params[nInputs_+3]);
      }
    }
    if (psGMMode_ == 1)
    { 
      printf("INFO: non-PD covariance matrix, terminate in GM mode.\n");
      printf("      diagnostics information is stored in errmat.m\n");
      FILE *fp = fopen("errmat.m", "w");
      fprintf(fp, "n = %d;\n", nSamples_);
      fprintf(fp, "f = %e;\n", sqrt(1.0*nSamples_));
      fprintf(fp, "A = [\n");
      for (ii = 0; ii < nSamples_; ii++) 
      {
        for (jj = 0; jj < nSamples_; jj++) 
          fprintf(fp, "%16.8e ", CMatrix.getEntry(ii,jj));
        fprintf(fp, "\n");
      }
      fprintf(fp, "];\n");
      fprintf(fp, "L = [\n");
      for (ii = 0; ii < nSamples_*(nSamples_-1)/2*nInputs_; ii++) 
        fprintf(fp, "%24.16e\n", XDistances_[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "H = [\n");
      for (jj = 0; jj < nInputs_; jj++)
        fprintf(fp, "%e\n", exp(2.0*params[jj]));
      fprintf(fp, "];\n");
      fprintf(fp, "%e\n", exp(params[nInputs_]));
      fprintf(fp, "%e\n", exp(params[nInputs_+1]));
      fprintf(fp, "%e\n", params[nInputs_+2]);
      fprintf(fp, "%e\n", exp(params[nInputs_+3]));
      fprintf(fp, "count = 0;\n");
      fprintf(fp, "for ii = 1 : n\n");
      fprintf(fp, "  B(ii,ii) = H(4) + ll * H(5) + H(7);\n");
      fprintf(fp, "  for jj = ii+1 : n \n");
      fprintf(fp, "    dist = 0;\n");
      fprintf(fp, "    for kk = 1 : 3\n");
      fprintf(fp, "      dist = dist + L(count*3+kk)^2 / H(kk);\n");
      fprintf(fp, "    end;\n");
      fprintf(fp, "    ddata = H(4) * exp(-0.5*dist);\n");
      fprintf(fp, "    B(ii,jj) = ddata + H(5);\n");
      fprintf(fp, "    B(jj,ii) = ddata + H(5);\n");
      fprintf(fp, "    count = count + 1;\n");
      fprintf(fp, "  end;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "Bmin = min(eig(B))\n");
      fclose(fp);
      exit(1);
    }
  }
  for (jj = 0; jj < nSamples_; jj++) 
    likelihood += 0.5 * (CInvY[jj] * YVec[jj]);
  if (outputLevel_ > 2) 
    printf("   computeGradients: Likelihood = %e\n", likelihood);
  delete [] TMat;
  return likelihood;
}

// ************************************************************************
// compute pairwise distances 
// ------------------------------------------------------------------------
int GP3::computeDistances()
{
  int    ii, jj, kk, count, error=0;
  double *LDists, dist;

  LDists = new double[(nSamples_*(nSamples_-1)/2)*nInputs_];
  checkAllocate(LDists, "LDists in GP3::computeDistances");
  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        LDists[count*nInputs_+ii] = XDataN_[jj*nInputs_+ii] -
                                    XDataN_[kk*nInputs_+ii];
        if (LDists[count*nInputs_+ii] < 0)
           LDists[count*nInputs_+ii] = - LDists[count*nInputs_+ii];
        dist += pow(LDists[count*nInputs_+ii], 2.0);
      }
      if (dist == 0.0)
      {
        printf("GP3 ERROR: repeated sample points.\n");
        printf("           Prune repeated points and re-run.\n");
        printf("Sample %d : (with sample point %d)\n", kk+1, jj+1);
        for (ii = 0; ii < nInputs_; ii++)
          printf("   Input %d : %e\n",ii+1,
                 XDataN_[kk*nInputs_+ii]*XStds_[ii] +XMeans_[ii]);
        error = 1;
      }
      count++;
    }
  }
  XDistances_.load(count*nInputs_, LDists);
  return error;
}

// ************************************************************************
// interpolation 
// ------------------------------------------------------------------------
int GP3::interpolate(int npts, double *XX, double *Y, double *Ystds)
{
  int    ii, kk, nn, iOne=1, status, offset;
  double dist, expn, ddata;
  psVector YT, YT2;

  YT.setLength(nSamples_);
  YT2.setLength(nSamples_);
  Y[0] = 0;
  offset = nInputs_ + 4;
  if (optLinTerm_) offset += nInputs_;
  for (nn = 0; nn < npts; nn++)
  {
    for (kk = 0; kk < nSamples_; kk++)
    {
      expn = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        dist = XDataN_[kk*nInputs_+ii] - XX[nn*nInputs_+ii];
        expn += dist * dist / exp(2.0*hyperparameters_[ii]);
      }
      expn *= 0.5;
      YT[kk] = exp(hyperparameters_[nInputs_]) * exp(-expn);
    }
    Y[nn] = hyperparameters_[nInputs_+2];
    for (kk = 0; kk < nSamples_; kk++)
      Y[nn] += YT[kk] * CInvY_[kk];
    if (Ystds != NULL)
    {
      Ystds[nn] = 0.0;
#if 1
      ddata = exp(hyperparameters_[nInputs_]) + 
              exp(hyperparameters_[nInputs_+1]) +
              exp(hyperparameters_[nInputs_+3]);
      CMatrix_.LUSolve(YT, YT2);
      for (kk = 0; kk < nSamples_; kk++)
        ddata -= YT[kk] * YT2[kk];
      Ystds[nn] = ddata;
#endif
    }
  }
  return 0;
}

// ************************************************************************
// construct C matrix
// ------------------------------------------------------------------------
void GP3::constructCMatrix(psMatrix &CMatrix, double *params)
{
  int    jj, kk, ii, count, status, offset;
  double dist, *TMat, likelihood, sqrt2pi=sqrt(2*3.1415928);
  double ddata, fudge;
  psMatrix CPrime, CT;
  psVector YVec, CInvY, YT;

  if (outputLevel_ > 3) 
  {
    printf("Construct CMatrix:\n");
    for (jj = 0; jj < hyperparameters_.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, params[jj]);
  }
  TMat = new double[nSamples_*nSamples_];
  fudge = sqrt(nSamples_);
  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    TMat[jj*nSamples_+jj] = exp(params[nInputs_]) + 
                fudge*exp(params[nInputs_+1]) + exp(params[nInputs_+3]);
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                exp(2.0*params[ii]);
      dist *= 0.5;
      dist = exp(params[nInputs_]) * exp(-dist);
      if (dist < 1.0e-50) dist = 0;
      TMat[jj*nSamples_+kk] = dist + exp(params[nInputs_+1]);
      TMat[kk*nSamples_+jj] = dist + exp(params[nInputs_+1]);
      count++;
    }
  }
  CMatrix.load(nSamples_, nSamples_, TMat);
  delete [] TMat;
  return;
}

// ************************************************************************
// generate C code
// ------------------------------------------------------------------------
void GP3::genCode()
{
  int ii, jj;
  FILE *fp = fopen("psuade_rs.info", "w");
  if (psRSCodeGen_ == 0) return;
  if (fp == NULL) return;
  fprintf(fp,"This file contains information to re-construct GP\n");
  fprintf(fp,"response surface offline. Follow the steps below:\n");
  fprintf(fp,"1. Search for the keywords 'SPLIT HERE' in this file.\n");
  fprintf(fp,"2. Store the lines below keywords into main.c\n");
  fprintf(fp,"3. Add the to-be-evaluated points inside main() in main.c\n");
  fprintf(fp,"   (or, embed all except main() in your own program).\n");
  fprintf(fp,"4. Compile main.c (cc -o main main.c -lm) and run\n");
  fprintf(fp,"   (Do not change the psuade_rs.info file.\n");
  fprintf(fp,"Note: if std dev is desired, uncomment the corresponding\n");
  fprintf(fp,"      section and compile with: cc main.c -llapack -lm\n");
  fprintf(fp,"\n");
  fprintf(fp,"PSUADE_BEGIN\n");
  fprintf(fp, "%d %d\n", nSamples_, nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    fprintf(fp, "%24.16e %24.16e\n", XMeans_[ii], XStds_[ii]);
  fprintf(fp, "%24.16e %24.16e\n", YMean_, YStd_);
  fprintf(fp, "%d\n", hyperparameters_.length());
  for (jj = 0; jj < hyperparameters_.length(); jj++)
    fprintf(fp, "%24.16e \n", hyperparameters_[jj]);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, "%24.16e ", XDataN_[jj*nInputs_+ii]);
    fprintf(fp, "\n");
  }
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, "%24.16e\n", CInvY_[jj]);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nSamples_; ii++)
      fprintf(fp, "%24.16e ", CMatrix_.getEntry(jj,ii));
    fprintf(fp, "\n");
  }
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, "%d\n", CMatrix_.pivots_[jj]);
  fprintf(fp,"==================== SPLIT HERE =====================\n");
  fprintf(fp,"/* ***********************************************/ \n");
  fprintf(fp,"/* GP interpolator from PSUADE.                  */ \n");
  fprintf(fp,"/* To estimate prediction uncertainty uncomment, */ \n");
  fprintf(fp,"/* dgetrs and the corresponding code segment.    */ \n");
  fprintf(fp,"/* ==============================================*/ \n");
  fprintf(fp,"#include <math.h>\n");
  fprintf(fp,"#include <stdlib.h>\n");
  fprintf(fp,"#include <stdio.h>\n");
  fprintf(fp,"/*void dgetrs_(char*,int*,int*,double*,int*,int*,\n");
  fprintf(fp,"               double*,int*,int*);*/\n");
  fprintf(fp,"int initialize();\n");
  fprintf(fp,"int finalize();\n");
  fprintf(fp,"int interpolate(int,double*,double*,double*);\n");
  fprintf(fp,"main(int argc, char **argv) {\n");
  fprintf(fp,"  int    i, iOne=1, nInps;\n");
  fprintf(fp,"  double X[%d], Y, Std;\n",nInputs_);
  fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
  fprintf(fp,"  if (argc < 2) {\n");
  fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
  fprintf(fp,"  if (fIn == NULL) {\n");
  fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
  fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
  fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
  fprintf(fp,"    exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",
          nInputs_);
  fprintf(fp,"  fclose(fIn);\n");
  fprintf(fp,"  initialize();\n");
  fprintf(fp,"  interpolate(iOne, X, &Y, &Std);\n");
  fprintf(fp,"  printf(\"Y = %%e (stdev = %%e)\\n\", Y, Std);\n");
  fprintf(fp,"  finalize();\n");
  fprintf(fp,"  if (argc == 3) {\n");
  fprintf(fp,"    fOut = fopen(argv[2], \"w\");\n");
  fprintf(fp,"    if (fOut == NULL) {\n");
  fprintf(fp,"      printf(\"ERROR: cannot open output file.\\n\");\n");
  fprintf(fp,"      exit(1);\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"    fprintf(fOut,\" %%e\\n\", Y);\n");
  fprintf(fp,"  fclose(fOut);}\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"/* Regression interpolation function         */\n");
  fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
  fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
  fprintf(fp," * ... */\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int    nInps, nSamples, nParams, *pivs=NULL;\n");
  fprintf(fp,"double *XMeans=NULL,*XStds=NULL,YMean,YStd,*Thetas=NULL;\n");
  fprintf(fp,"double *CMat=NULL,*CInvY=NULL, *XN=NULL; \n");
  fprintf(fp,"int interpolate(int npts,double *X,double *Y,\n");
  fprintf(fp,"                double *YStds) {\n");
  fprintf(fp,"  int    ss, ii, jj, kk;\n");
  fprintf(fp,"  double expn, *xt, *yt, *zt, dist, ddata;\n");
  fprintf(fp,"  xt = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  yt = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  zt = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  for (ss = 0; ss < npts; ss++) {\n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++)\n");
  fprintf(fp,"      xt[ii] = (X[ss*nInps+ii]-XMeans[ii])/XStds[ii];\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++) {\n");
  fprintf(fp,"      expn = 0.0;\n");
  fprintf(fp,"      for (ii = 0; ii < nInps; ii++) {\n");
  fprintf(fp,"        dist = XN[kk*nInps+ii] - xt[ii];\n");
  fprintf(fp,"        expn += dist*dist/exp(2.0*Thetas[ii]);}\n");
  fprintf(fp,"      expn *= 0.5;\n");
  fprintf(fp,"      yt[kk] = exp(Thetas[nInps]) * exp(-expn);}\n");
  fprintf(fp,"    Y[ss] = Thetas[nInps+2];\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++)\n");
  fprintf(fp,"      Y[ss] += yt[kk] * CInvY[kk];\n");
  fprintf(fp,"    Y[ss] = Y[ss] * YStd + YMean;\n");
  fprintf(fp,"    YStds[ss] = 0.0;\n");
  fprintf(fp,"    /* ==== if need to compute std dev. =====\n");
  fprintf(fp,"    ddata = exp(Thetas[nInps])+exp(Thetas[nInps+1])+");
  fprintf(fp,"exp(Thetas[nInps+3]);\n");
  fprintf(fp,"    LUSolve(yt, zt);\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++)\n");
  fprintf(fp,"      ddata -= yt[kk] * zt[kk];\n");
  fprintf(fp,"    YStds[ss] = sqrt(ddata); */\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  free(xt);\n");
  fprintf(fp,"  free(yt);\n");
  fprintf(fp,"  free(zt);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int LUSolve(double *X,double *Y){\n");
  fprintf(fp,"  int ii, jj, iOne=1, status, M;\n");
  fprintf(fp,"  char   trans='N';\n");
  fprintf(fp,"  M = nSamples;\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++) Y[ii] = X[ii];\n");
  fprintf(fp,"  /*dgetrs_(&trans,&M,&iOne,CMat,&M,pivs,Y,&M,&status);*/\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int initialize() {\n");
  fprintf(fp,"  int    ii, jj;\n");
  fprintf(fp,"  double ddata;\n");
  fprintf(fp,"  char   line[1001], word[1001];\n");
  fprintf(fp,"  FILE *fp = fopen(\"psuade_rs.info\", \"r\");\n");
  fprintf(fp,"  if (fp == NULL){\n");
  fprintf(fp,"     printf(\"Data file (psuade_rs.info) not found.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  while (1) {\n");
  fprintf(fp,"    fgets(line, 500, fp);\n");
  fprintf(fp,"    sscanf(line, \"%%s\",word);\n");
  fprintf(fp,"    if (!strcmp(word, \"PSUADE_BEGIN\")) break;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fp, \"%%d %%d\", &nSamples, &nInps);\n");
  fprintf(fp,"  XMeans = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  XStds  = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  Thetas = (double *) malloc((nInps+4)*sizeof(double));\n");
  fprintf(fp,"  for (ii = 0; ii < nInps; ii++) {\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"    XMeans[ii] = ddata;\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"    XStds[ii] = ddata;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fp, \"%%lg %%lg\", &YMean, &YStd);\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nParams);\n");
  fprintf(fp,"  for (ii = 0; ii < nParams; ii++) {\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"    Thetas[ii] = ddata;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  XN = (double *) malloc(nSamples*nInps*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++) \n");
  fprintf(fp,"      fscanf(fp, \"%%lg\", &XN[jj*nInps+ii]);\n");
  fprintf(fp,"  CInvY = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &CInvY[ii]);\n");
  fprintf(fp,"  CMat=(double*) malloc(nSamples*nSamples*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"      fscanf(fp, \"%%lg \", &CMat[jj+ii*nSamples]);\n");
  fprintf(fp,"  pivs = (int *) malloc(nSamples*sizeof(int));\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"    fscanf(fp, \"%%d\", &pivs[ii]);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int finalize() {\n");
  fprintf(fp,"  if (XMeans != NULL) free(XMeans);\n");
  fprintf(fp,"  if (XStds  != NULL) free(XStds);\n");
  fprintf(fp,"  if (Thetas != NULL) free(Thetas);\n");
  fprintf(fp,"  if (CInvY  != NULL) free(CInvY);\n");
  fprintf(fp,"  if (CMat   != NULL) free(CMat);\n");
  fprintf(fp,"  if (XN     != NULL) free(XN);\n");
  fprintf(fp,"  if (pivs   != NULL) free(pivs);\n");
  fprintf(fp,"}\n");
  fclose(fp);
  printf("FILE psuade_rs.info contains the GP interpolator.\n");
}


















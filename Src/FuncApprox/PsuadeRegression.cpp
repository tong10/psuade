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
// Functions for the class PsuadeRegression
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PsuadeRegression.h"
#include "sysdef.h"
#include "PrintingTS.h"
#include "Sampling.h"

// ************************************************************************
// External function 
// ------------------------------------------------------------------------
extern "C" {
  void newuoa_(int *,int *,double *,double *,double *,int *,int *,double*);
}

// ************************************************************************
// Resident function to be called by newuoa 
// ------------------------------------------------------------------------
PsuadeRegression *PsRegr_=NULL;
psVector PsRegr_OptX;
double   PsRegr_OptY;
psVector PsRegr_SampleX;
psVector PsRegr_SampleY;

#ifdef __cplusplus
extern "C"
{
#endif
  void *psrnewuoaevalfunc_(int *nInps, double *XValues, double *YValue)
  {
    if (PsRegr_ == NULL)
    {
      printf("PsuadeRegression ERROR: no regression object in ");
      printf("function evaluation\n");
      exit(1);
    }
    psVector VecP;
    VecP.load(*nInps, XValues);
    (*YValue) = PsRegr_->optFunction(VecP,PsRegr_SampleX,PsRegr_SampleY);
    if ((*YValue) < PsRegr_OptY)
    {
      PsRegr_OptY = *YValue;
      PsRegr_OptX = VecP;
    }
    return NULL;
  }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
PsuadeRegression::PsuadeRegression(int nInputs,int nSamples):
                                   FuncApprox(nInputs,nSamples)
{
  faID_ = PSUADE_RS_LOCAL;

  //**/ ===============================================================
  // display banner and additonal information
  //**/ ===============================================================
  if (psConfig_.InteractiveIsOn() && outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*                Psuade Regression Analysis\n");
    printf("* To use this regression, a user needs to go into the\n");
    printf("* PsuadeRegression.cpp file and insert their function\n");
    printf("* (to replace userFunction and optFunction).\n");
    printEquals(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PsuadeRegression::~PsuadeRegression()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PsuadeRegression::initialize(double *XIn, double *YIn)
{
  //**/ These values are for reference only
  int nParams=3;
  VecCoeffs_.setLength(nParams);
  //VecCoeffs_[0] = -8.8817 * 1.807;
  //VecCoeffs_[1] = -8.8817;
  //VecCoeffs_[2] = -8.8817 * 1.735;
  PsRegr_SampleX.load(nInputs_*nSamples_,XIn);
  PsRegr_SampleY.load(nSamples_,YIn);

  int doOptimization = 1;

  if (doOptimization > 0)
  {
    //**/ set up bounds for generating multi-start points)
    psVector vecLB, vecUB;
    vecLB.setLength(nParams);
    vecUB.setLength(nParams);
    vecLB[0] = -30; vecUB[0] = 0;
    vecLB[1] = -16; vecUB[1] = 0;
    vecLB[2] =  10; vecUB[2] = 20;
    if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
    {
      for (int ii = 0; ii < nParams; ii++)
        printf("   Parameter bounds %3d = %12.4e %12.4e\n",ii+1,
               vecLB[ii],vecUB[ii]);
    }
    //**/ create multi-starting points
    int nSamp = 5;
    Sampling *sampler = NULL;
    if (nParams > 51)
         sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nParams,vecLB.getDVector(),
                            vecUB.getDVector());
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 0);
    sampler->initialize(0);
    psVector  vecXS, vecYS;
    psIVector vecIS;
    vecIS.setLength(nSamp);
    vecXS.setLength(nSamp*nParams);
    vecYS.setLength(nSamp);
    sampler->getSamples(nSamp,nParams,1,vecXS.getDVector(),
                        vecYS.getDVector(), vecIS.getIVector());
    delete sampler;
    for (int ii = 0; ii < nParams; ii++) 
      vecXS[ii] = 0.5 * (vecLB[ii] + vecUB[ii]);

    //**/ set up for newuoa
    int    nPts, pLevel=9999, maxfun=1000, kk1, kk2;
    double tol=1e-5, rhobeg, rhoend, dtemp;
    psVector vecW;
    nPts = (nParams + 1) * (nParams + 2) / 2;
    kk1 = (nPts+13) * (nPts+nParams) + 3*nParams*(nParams+3)/2;
    kk2 = (nPts+5)*(nPts+nParams)+3*nInputs_*(nParams+5)/2+1;
    if (kk1 > kk2) vecW.setLength(kk1);
    else           vecW.setLength(kk2);
    rhobeg = vecUB[0] - vecLB[0];
    for (int ii = 1; ii < nParams; ii++)
    {
      dtemp = vecUB[ii] - vecLB[ii];
      if (dtemp < rhobeg) rhobeg = dtemp;
    }
    rhobeg *= 0.5;
    rhoend = rhobeg * tol;
 
    //**/ run newuoa
    int its = 0, jj;

    PsRegr_OptY = 1e35;
    if (PsRegr_ != NULL) PsRegr_ = NULL;
    PsRegr_ = new PsuadeRegression(nInputs_, nSamples_);
    if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
      printf("PsuadeRegression NEWUOA optimization (%d):\n",nSamp);

    psVector vecPVals;
    vecPVals.setLength(nParams);
    for (int ss = 0; ss < nSamp; ss++)
    {
      if (psConfig_.InteractiveIsOn() && outputLevel_ > 3 &&
          nSamp > 1)
        printf("PsuadeRegression multi-start %d (%d)\n",ss+1, nSamp);
      for (int ii = 0; ii < nParams; ii++) 
        vecPVals[ii] = vecXS[ss*nParams+ii];
      if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
      {
        printf("   PsuadeRegression initial parameter values:\n");
        for (int ii = 0; ii < vecPVals.length(); ii++)
          printf("   Parameter %2d = %e\n",ii+1,vecPVals[ii]);
      }
#ifdef HAVE_NEWUOA
      newuoa_(&nParams,&nPts,vecPVals.getDVector(),&rhobeg,&rhoend,
              &pLevel, &maxfun, vecW.getDVector());
#else
      printf("ERROR: NEWUOA optimizer unavailable.\n");
      exit(1);
#endif
      if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
      {
        printf("   PsuadeRegression best parameter values so far:\n");
        for (int ii = 0; ii < PsRegr_OptX.length(); ii++)
          printf("   Parameter %2d = %e\n",ii+1,
                 PsRegr_OptX[ii]);
        printf("   Best objective function so far = %e\n",PsRegr_OptY);
      }
    }
    VecCoeffs_ = PsRegr_OptX;
    delete PsRegr_;
    PsRegr_ = NULL;
    if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
    {
      printf("PsuadeRegression final best parameter values:\n");
      for (int ii = 0; ii < VecCoeffs_.length(); ii++)
        printf("   Parameter %2d = %24.16e\n",ii+1,
               VecCoeffs_[ii]);
      printf("   Final best objective function = %e\n",PsRegr_OptY);
    }
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int PsuadeRegression::genNDGridData(double *X, double *Y, int *N2,
                                    double **X2, double **Y2)
{
  int totPts, ss;
  psVector vecYOut;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if ((*N2) <= 0)
  {
    printf("PsuadeRegression::genNDGridData - ERROR detected (N <= 0).\n");
    (*N2) = 0;
    return -1;
  }

  //**/ ===============================================================
  //**/ return if there is no request to create lattice points
  //**/ ===============================================================
  if ((*N2) == -999) return 0;

  //**/ ===============================================================
  //**/ generating regular grid data
  //**/ ===============================================================
  genNDGrid(N2, X2);
  if ((*N2) <= 0) return 0;
  totPts = (*N2);

  //**/ ===============================================================
  //**/ allocate storage for the data points and generate them
  //**/ ===============================================================
  vecYOut.setLength(totPts);
  (*Y2) = vecYOut.takeDVector();
  for (ss = 0; ss < totPts; ss++)
    (*Y2)[ss] = evaluatePoint(&((*X2)[ss*nInputs_]));

  return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int PsuadeRegression::gen1DGridData(double *X, double *Y, int ind1,
                                    double *settings, int *NN, 
                                    double **XX, double **YY)
{
  int    totPts, mm, nn;
  double HX;
  psVector vecXT, vecXOut, vecYOut;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if ((*NN) <= 0)
  {
    printf("PsuadeRegression::gen1DGridData - ERROR detected (N <= 0).\n");
    (*NN) = 0;
    return -1;
  }

  //**/ ===============================================================
  //**/ set up for generating regular grid data
  //**/ ===============================================================
  totPts = nPtsPerDim_;
  HX = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1); 

  //**/ ===============================================================
  //**/ allocate storage for and then generate the data points
  //**/ ===============================================================
  (*NN) = totPts;
  vecXOut.setLength(totPts);
  vecYOut.setLength(totPts);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    vecXT[ind1] = HX * mm + VecLBs_[ind1];
    (*XX)[mm] = vecXT[ind1];
    (*YY)[mm] = evaluatePoint(vecXT.getDVector());
  }
  return 0;
}

// ************************************************************************
// Generate 2D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int PsuadeRegression::gen2DGridData(double *X, double *Y, int ind1,
                                    int ind2, double *settings, int *NN, 
                                    double **XX, double **YY)
{
  int totPts, mm, nn, index;
  psVector vecXT, vecXOut, vecYOut, vecHX;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if ((*NN) <= 0)
  {
    printf("PsuadeRegression::gen2DGridData - ERROR detected (N <= 0).\n");
    (*NN) = 0;
    return -1;
  }

  //**/ ===============================================================
  //**/ set up for generating regular grid data
  //**/ ===============================================================
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 

  //**/ ===============================================================
  //**/ allocate storage for and then generate the data points
  //**/ ===============================================================
  (*NN) = totPts;
  vecXOut.setLength(totPts*2);
  vecYOut.setLength(totPts);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      index = mm * nPtsPerDim_ + nn;
      vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
      vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
      (*XX)[index*2]   = vecXT[ind1];
      (*XX)[index*2+1] = vecXT[ind2];
      (*YY)[index] = evaluatePoint(vecXT.getDVector());
    }
  }
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PsuadeRegression::gen3DGridData(double *X, double *Y, int ind1,
                                    int ind2, int ind3, double *settings, 
                                    int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, index;
  psVector vecXT, vecXOut, vecYOut, vecHX;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if ((*NN) <= 0)
  {
    printf("PsuadeRegression::gen3DGridData - ERROR detected (N <= 0).\n");
    (*NN) = 0;
    return -1;
  }

  //**/ ===============================================================
  //**/ set up for generating regular grid data
  //**/ ===============================================================
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3])/(nPtsPerDim_ - 1); 

  //**/ ===============================================================
  //**/ allocate storage for and then generate the data points
  //**/ ===============================================================
  (*NN) = totPts;
  vecXOut.setLength(totPts*3);
  vecYOut.setLength(totPts);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
        vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
        vecXT[ind3] = vecHX[2] * pp + VecLBs_[ind3];
        (*XX)[index*3]   = vecXT[ind1];
        (*XX)[index*3+1] = vecXT[ind2];
        (*XX)[index*3+2] = vecXT[ind3];
        (*YY)[index] = evaluatePoint(vecXT.getDVector());
      }
    }
  }
  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PsuadeRegression::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                                  int ind3, int ind4, double *settings, 
                                  int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, qq, index;
  psVector vecXT, vecXOut, vecYOut, vecHX;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if ((*NN) <= 0)
  {
    printf("PsuadeRegression::gen4DGridData - ERROR detected (N <= 0).\n");
    (*NN) = 0;
    return -1;
  }

  //**/ ===============================================================
  //**/ set up for generating regular grid data
  //**/ ===============================================================
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3])/(nPtsPerDim_ - 1); 
  vecHX[3] = (VecUBs_[ind4] - VecLBs_[ind4])/(nPtsPerDim_ - 1); 

  //**/ ===============================================================
  //**/ allocate storage for and then generate the data points
  //**/ ===============================================================
  (*NN) = totPts;
  vecXOut.setLength(totPts*4);
  vecYOut.setLength(totPts);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        for (qq = 0; qq < nPtsPerDim_; qq++)
        {
          index = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                  nn*nPtsPerDim_*nPtsPerDim_ + pp * nPtsPerDim_ + qq;
          vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
          vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
          vecXT[ind3] = vecHX[2] * pp + VecLBs_[ind3];
          vecXT[ind4] = vecHX[3] * qq + VecLBs_[ind4];
          (*XX)[index*4]   = vecXT[ind1];
          (*XX)[index*4+1] = vecXT[ind2];
          (*XX)[index*4+2] = vecXT[ind3];
          (*XX)[index*4+3] = vecXT[ind4];
          (*YY)[index] = evaluatePoint(vecXT.getDVector());
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double PsuadeRegression::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y;
  userFunction(nInputs_, X, iOne, &Y);
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double PsuadeRegression::evaluatePoint(int npts, double *X, double *Y)
{
  int ii, iOne=1;
  for (ii = 0 ; ii < npts; ii++)
     userFunction(nInputs_, &X[nInputs_*ii], iOne, &Y[ii]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double PsuadeRegression::evaluatePointFuzzy(double *X, double &stdev)
{
  int    iOne=1;
  double Y;
  evaluatePoint(iOne, X, &Y);
  stdev = 0.0;
  return Y;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double PsuadeRegression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                            double *Ystd)
{
  evaluatePoint(npts, X, Y);
  for (int ii = 0 ; ii < npts; ii++) Ystd[ii] = 0;
  return 0.0;
}

// ************************************************************************
// User function 
// ------------------------------------------------------------------------
int PsuadeRegression::userFunction(int nInputs, double *inputs, 
                                   int nOutputs, double *outputs)
{
  if (nOutputs != 1)
  {
    printf("PsuadeRegression user function ERROR: nOutputs != 1\n");
    exit(1);
  }
  outputs[0] = 2.0*VecCoeffs_[0]/(inputs[1] * inputs[2]) + 
               VecCoeffs_[1] * pow(inputs[0],2.0) + VecCoeffs_[2];
  if (outputs[0] < 0) outputs[0] = 0;
  return 0;
}

// ************************************************************************
// function for parameter optimization
// ------------------------------------------------------------------------
double PsuadeRegression::optFunction(psVector vecP, psVector sampleX,
                                     psVector sampleY)
{
  int nSamp = sampleY.length();
  int nInps = sampleX.length() / nSamp;
  double errNorm = 0, ddata;
  for (int ss = 0; ss < nSamp; ss++)
  {
    ddata = 2.0 * vecP[0]/(sampleX[ss*nInps+1] * sampleX[ss*nInps+2]) + 
            vecP[1] * pow(sampleX[ss*nInps],2.0) + vecP[2];
    if (ddata < 0) ddata = 0;
    errNorm += pow(ddata - sampleY[ss], 2.0);
  }
  errNorm = sqrt(errNorm);
  return errNorm;
}


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
// Functions for the class PWLinear
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "sysdef.h"
#include "pData.h"
#include "PsuadeUtil.h"
#include "PWLinear.h"

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
PWLinear::PWLinear(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
  faID_ = PSUADE_RS_PWL;
  nROMs_ = nSamples;
  threshold_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PWLinear::~PWLinear()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PWLinear::initialize(double *X, double *Y)
{
  int     totPts, mm, iOne=1;
  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("PWLinear::genNDGridData ERROR - invalid argument.\n");
    exit(1);
  }
  if (nSamples_ <= nInputs_)
  {
    printf("PWLinear::genNDGridData INFO - not enough points.\n");
    return 0;
  }

  //**/ ---------------------------------------------------------------
  //**/ return if there is no request to create lattice points
  //**/ ---------------------------------------------------------------
  if (nInputs_ > 12)
  {
    printf("PWLinear::genNDGridData INFO - nInputs > 12.\n");
    printf("          No lattice points generated.\n");
    return 0;
  }

  //**/ ---------------------------------------------------------------
  //**/ create PWL
  //**/ ---------------------------------------------------------------
  if (VecROMStorePts_.length() == 0) setParams(iOne, NULL);
  if (psConfig_.RSCodeGenIsOn())
    printf("PWLinear INFO: response surface stand-alone code not available.\n");
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int PWLinear::genNDGridData(double *X, double *Y, int *NN, double **XX,
                            double **YY)
{
  int totPts, mm, iOne=1;

  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("PWLinear::genNDGridData ERROR - invalid argument.\n");
    exit(1);
  }
  if (nSamples_ <= nInputs_)
  {
    printf("PWLinear::genNDGridData INFO - not enough points.\n");
    return 0;
  }

  //**/ ---------------------------------------------------------------
  //**/ return if there is no request to create lattice points
  //**/ ---------------------------------------------------------------
  if (nInputs_ > 12)
  {
    printf("PWLinear::genNDGridData INFO - nInputs > 12.\n");
    printf("          No lattice points generated.\n");
    (*NN) = 0;
    (*XX) = 0;
    (*YY) = 0;
    return 0;
  }

  //**/ ---------------------------------------------------------------
  //**/ create PWL
  //**/ ---------------------------------------------------------------
  if (VecROMStorePts_.length() == 0) setParams(iOne, NULL);
  if ((*NN) == -999) return 0;

  //**/ ---------------------------------------------------------------
  //**/ generating regular grid data
  //**/ ---------------------------------------------------------------
  genNDGrid(NN, XX);
  if ((*NN) == 0) return 0;
  totPts = (*NN);

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for the data points and generate them
  //**/ ---------------------------------------------------------------
  psVector vecYOut;
  vecYOut.setLength(totPts);
  (*YY) = vecYOut.takeDVector();
  for (mm = 0; mm < totPts; mm++)
    (*YY)[mm] = evaluatePoint(&((*XX)[mm*nInputs_]));
  return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen1DGridData(double *X, double *Y, int ind1,
                            double *settings, int *NN,
                            double **XX, double **YY)
{
  int    totPts, mm, nn;
  double HX;
  psVector vecXT;

  //**/ ---------------------------------------------------------------
  //**/ create response surface
  //**/ ---------------------------------------------------------------
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_;
  HX = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1);

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  psVector vecXOut, vecYOut;
  vecXOut.setLength(totPts*2);
  (*XX) = vecXOut.takeDVector();
  vecYOut.setLength(totPts);
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
// Generate 2D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                            double *settings, int *NN, double **XX, 
                            double **YY)
{
  int totPts, mm, nn, index;
  psVector vecHX, vecXT;

  //**/ ---------------------------------------------------------------
  //**/ create response surface
  //**/ ---------------------------------------------------------------
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1);
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2]) / (nPtsPerDim_ - 1);

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  psVector vecXOut, vecYOut;
  vecXOut.setLength(totPts*2);
  (*XX) = vecXOut.takeDVector();
  vecYOut.setLength(totPts);
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
int PWLinear::gen3DGridData(double *X, double *Y, int ind1, int ind2,
                            int ind3, double *settings, int *NN,
                            double **XX, double **YY)
{
  int totPts, mm, nn, pp, index;
  psVector vecHX, vecXT;

  //**/ ---------------------------------------------------------------
  //**/ create response surface
  //**/ ---------------------------------------------------------------
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1);
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2]) / (nPtsPerDim_ - 1);
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3]) / (nPtsPerDim_ - 1);

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  psVector vecXOut, vecYOut;
  vecXOut.setLength(totPts*3);
  (*XX) = vecXOut.takeDVector();
  vecYOut.setLength(totPts);
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
int PWLinear::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                            int ind3, int ind4, double *settings, int *NN,
                            double **XX, double **YY)
{
  int totPts, mm, nn, pp, qq, index;
  psVector vecHX, vecXT;

  //**/ ---------------------------------------------------------------
  //**/ create response surface
  //**/ ---------------------------------------------------------------
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1);
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2]) / (nPtsPerDim_ - 1);
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3]) / (nPtsPerDim_ - 1);
  vecHX[3] = (VecUBs_[ind4] - VecLBs_[ind4]) / (nPtsPerDim_ - 1);

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  psVector vecXOut, vecYOut;
  vecXOut.setLength(totPts*4);
  (*XX) = vecXOut.takeDVector();
  vecYOut.setLength(totPts);
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
double PWLinear::evaluatePoint(double *X)
{
  int    iR, ii, count;
  double Y, Ytemp, dist2, diff, accum, dtemp, thresh, minDist2;

  thresh = threshold_ * 1.05;
  while (thresh < 2.0 * threshold_)
  {
    count = 0;
    accum = 0.0;
    Y = 0.0;
    minDist2 = PSUADE_UNDEFINED;
    for (iR = 0; iR < nROMs_; iR++)
    {
      Ytemp = VecROMStoreVals_[iR];
      dist2 = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        diff = X[ii] - VecROMStorePts_[iR*nInputs_+ii]; 
        Ytemp += VecROMStoreGrads_[iR*nInputs_+ii] * diff;
        dist2 += diff * diff;
      }
      dist2 = 0.5 * dist2 * VecROMStoreEigens_[iR];
      if (dist2 <= 1.0e-12)
      { 
        Y = Ytemp;
        return Y;
      }
      if (dist2 < thresh)
      {
        dtemp = (thresh - dist2) / dist2;
        if (dtemp < 0.0) dtemp = 0.0;
        if (dtemp > 0.0)
        {
          Y += Ytemp * dtemp;
          accum += dtemp;
          count++;
        }
      }
      if (dist2 < minDist2) minDist2 = dist2;
    }
    if (count == 0)
    {
      printf("PWLinear ERROR: no ROV found, min dist2 = %e > %e.\n",
             minDist2, thresh);
      printf("PWLinear ERROR: no ROV found, reduce threshold to %e.\n",
             thresh);
      thresh = 1.1 * thresh;
    }
    else 
    {
      Y /= accum;
      return Y;
    }
  }
  printf("PWLinear ERROR: no ROV found for point, return 0.\n");
  printf("Data point is: \n");
  for (ii = 0; ii < nInputs_; ii++)
    printf("X %3d = %e\n", ii+1, X[ii]);
  return 0.0;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double PWLinear::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviations
// ------------------------------------------------------------------------
double PWLinear::evaluatePointFuzzy(double *X, double &std)

{
  double Y;
  Y = evaluatePoint(X);
  std = 0.0;
  return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double PWLinear::evaluatePointFuzzy(int npts, double *X, double *Y,
                                    double *Ystd)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePointFuzzy(&X[kk*nInputs_], Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double PWLinear::setParams(int targc, char **targv)
{
  int  status, iR, ii, nOutputs, kk;
  char filename[500], lineIn[500];
  PsuadeData *psIO = NULL;
  pData      pPtr, pInputs, pOutputs; 

  if (targc == 1)
  {
    printf("PWLinear:: enter data file name : ");
    scanf("%s", filename);
    fgets(lineIn, 500, stdin);
    psIO = new PsuadeData();
    status = psIO->readPsuadeFile(filename);
    if (status != 0)
    {
      printf("ERROR: Problem reading file %s.\n", filename);
      exit(1);
    }
    psIO->getParameter("input_ninputs", pPtr);
    if (pPtr.intData_ != nInputs_)
    {
      printf("nInputs in file %s does not match with current file.\n",
             filename);
      exit(1);
    }
    psIO->getParameter("output_noutputs", pPtr);
    nOutputs = nInputs_ + nInputs_ * (nInputs_ + 1) / 2 + 2;
    if (pPtr.intData_ != nOutputs)
    {
      printf("nOutputs in file %s does not match with current file.\n",
             filename);
      exit(1);
    }
    psIO->getParameter("method_nsamples", pPtr);
    nROMs_ = pPtr.intData_;
    if (nROMs_ <= 0)
    {
      printf("ERROR: nROMs in file %s <= 0\n", filename);
      exit(1);
    }
    psIO->getParameter("input_sample", pInputs);
    VecROMStorePts_.setLength(nROMs_*nInputs_);
    for (iR = 0; iR < nROMs_*nInputs_; iR++)
      VecROMStorePts_[iR] = pInputs.dbleArray_[iR];
    psIO->getParameter("output_sample", pOutputs);
    VecROMStoreVals_.setLength(nROMs_);
    VecROMStoreGrads_.setLength(nROMs_*nInputs_);
    for (iR = 0; iR < nROMs_; iR++)
    {
      VecROMStoreVals_[iR] = pOutputs.dbleArray_[iR*nOutputs];
      for (ii = 0; ii < nInputs_; ii++)
        VecROMStoreGrads_[iR*nInputs_+ii] =
                              pOutputs.dbleArray_[iR*nOutputs+ii+1];
    }
    VecROMStoreEigens_.setLength(nROMs_);
    for (iR = 0; iR < nROMs_; iR++)
      VecROMStoreEigens_[iR] = pOutputs.dbleArray_[(iR+1)*nOutputs-1];
    psIO->getParameter("ana_threshold", pPtr);
    threshold_ = pPtr.dbleData_;
    if (threshold_ <= 0.0 || threshold_ >= 1.0)
    {
      printf("PWLinear ERROR: invalid analysis threshold %e.\n",threshold_);
    }
  }
  if (targc == 6)
  {
    printf("PWLinear: getting ROV information.\n");
    nROMs_ = *(int *) targv[0];
    nInputs_ = *(int *) targv[1];
    kk = *(int *) targv[2];
    double *sampleInputs = (double *) targv[3];
    double *sampleOutputs = (double *) targv[4];
    threshold_ = * (double *) targv[5];
    nOutputs = nInputs_ + nInputs_ * (nInputs_ + 1) / 2 + 2;
    if (kk != nOutputs)
    {
      printf("PWLinear ERROR: nOutputs mismatch (%d != %d).\n",
             kk, nOutputs);
      exit(1);
    }
    VecROMStorePts_.setLength(nROMs_*nInputs_);
    for (iR = 0; iR < nROMs_*nInputs_; iR++)
      VecROMStorePts_[iR] = sampleInputs[iR];
    VecROMStoreVals_.setLength(nROMs_);
    VecROMStoreGrads_.setLength(nROMs_*nInputs_);
    for (iR = 0; iR < nROMs_; iR++)
    {
      VecROMStoreVals_[iR] = sampleOutputs[iR*nOutputs];
      for (ii = 0; ii < nInputs_; ii++)
        VecROMStoreGrads_[iR*nInputs_+ii] =
                                 sampleOutputs[iR*nOutputs+ii+1];
    }
    VecROMStoreEigens_.setLength(nROMs_);
    for (iR = 0; iR < nROMs_; iR++)
      VecROMStoreEigens_[iR] = sampleOutputs[(iR+1)*nOutputs-1];
    if (threshold_ <= 0.0 || threshold_ >= 1.0)
    {
      printf("PWLinear ERROR: invalid threshold %e.\n",threshold_);
      exit(1);
    }
  }
  if (psIO != NULL) delete psIO;
  return 0.0;
}


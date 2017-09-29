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
// Functions for the class MGP3
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "GP3.h"
#include "MGP3.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include "Psuade.h"
#include "Sampling.h"
#include "MainEffectAnalyzer.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class MGP3
// ------------------------------------------------------------------------
MGP3::MGP3(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  char *strPtr, winput[1000], equal[100];
  faID_ = PSUADE_RS_MGP2;
  nPartitions_ = 0;
  rsPtrs_ = NULL;
  boxes_ = NULL;
  partSize_ = 500;
  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MGP_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &partSize_);
      if (partSize_ < 100) partSize_ = 500;
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  double ddata = log(1.0*nPartitions_) / log(2.0);
  int   idata = (int) ddata;
  if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MGP3::~MGP3()
{
  int ii;
  if (nPartitions_ > 0 && rsPtrs_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++)
      if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
    delete rsPtrs_;
    rsPtrs_ = NULL;
  }
  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++)
    {
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
    }
    delete [] boxes_;
  }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int MGP3::initialize(double *XIn, double *YIn)
{
  int    ii, jj, ss, incr, *indices, nSubs, index, samCnt;
  double *dVec, *vces, *XX, *YY, diff, var, ddata, *ubs, *lbs;
  MainEffectAnalyzer *me;
  if (nPartitions_ > 0 && rsPtrs_ != NULL)
  {
    for (int ii = 0; ii < nPartitions_; ii++)
      if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
    delete rsPtrs_;
    rsPtrs_ = NULL;
  }
  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++)
    {
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
    }
    delete [] boxes_;
  }

  XData_.setLength(nSamples_*nInputs_);
  for (ii = 0; ii < nSamples_*nInputs_; ii++) XData_[ii] = XIn[ii];
  YData_.setLength(nSamples_);
  for (ii = 0; ii < nSamples_; ii++) YData_[ii] = YIn[ii];
   

  me = new MainEffectAnalyzer();
  turnPrintTSOff();
  vces = new double[nInputs_];
  XX = XData_.getDVector();
  YY = YData_.getDVector();
  me->computeVCECrude(nInputs_,nSamples_,XX,YY,lowerBounds_,
                      upperBounds_, var, vces);
  delete me;
  indices = new int[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) indices[ii] = ii;
  sortDbleList2a(nInputs_, vces, indices);
  if (outputLevel_ > 1)
  {
    for (ii = 0; ii < nInputs_; ii++)
      printf("VCE %d = %e\n", indices[ii], vces[ii]);
  }

  nSubs = int(log(1.0*nPartitions_ + 1.0e-9) / log(2.0));
  boxes_ = new MGP3_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MGP3_Box();
    boxes_[ii]->lBounds_ = new double[nInputs_];
    boxes_[ii]->uBounds_ = new double[nInputs_];
    for (jj = 0; jj < nInputs_; jj++)
    {
      boxes_[ii]->lBounds_[jj] = lowerBounds_[jj];
      boxes_[ii]->uBounds_[jj] = upperBounds_[jj];
    }
  }
  for (ii = 0; ii < nSubs; ii++)
  {
    index = indices[nInputs_-1];
    if (outputLevel_ > 0)
    {
      for (jj = 0; jj < nInputs_; jj++)
        printf("vce %3d = %e\n", indices[jj]+1, vces[jj]);
      printf("Selected input for partition = %d\n", index+1);
    }
    incr = 1 << (nSubs - ii - 1);
    for (jj = 0; jj < nPartitions_; jj++)
    {
      if (((jj / incr) % 2) == 0)
      {
        boxes_[jj]->uBounds_[index] = 0.5 *
          (boxes_[jj]->uBounds_[index] + boxes_[jj]->lBounds_[index]);
      }
      else
      {
        boxes_[jj]->lBounds_[index] = 0.5 *
          (boxes_[jj]->uBounds_[index] + boxes_[jj]->lBounds_[index]);
      }
    }
    vces[nInputs_-1] *= 0.25;
    sortDbleList2a(nInputs_, vces, indices);
  }
  if (outputLevel_ > 3)
  {
    for (jj = 0; jj < nPartitions_; jj++)
    {
      printf("Partition %d:\n", jj);
        for (ii = 0; ii < nInputs_; ii++)
          printf("Input %2d = %12.4e %12.4e\n",ii+1,
            boxes_[jj]->lBounds_[ii], boxes_[jj]->uBounds_[ii]);
    }
  }
  double dcheck1=0, dcheck2=0;
  for (jj = 0; jj < nPartitions_; jj++)
  {
    ddata = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      ddata *= (boxes_[jj]->uBounds_[ii]-boxes_[jj]->lBounds_[ii]);
    dcheck2 += ddata;
  }
  dcheck1 = 1;
  for (ii = 0; ii < nInputs_; ii++)
    dcheck1 *= (upperBounds_[ii] - lowerBounds_[ii]);
  if (outputLevel_ > 2)
    printf("Coverage check: %e (sum) ?= %e (orig)\n",dcheck1,dcheck2);

  if (outputLevel_ > 1) printf("MGP3 training begins....\n");
  int total=0;
  XX = new double[nSamples_*nInputs_];
  YY = new double[nSamples_];
  rsPtrs_ = new GP3*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    lbs = boxes_[ii]->lBounds_;
    ubs = boxes_[ii]->uBounds_;
    samCnt = 0;
    for (ss = 0; ss < nSamples_; ss++)
    {
      for (jj = 0; jj < nInputs_; jj++)
      {
        diff = 0.05 * (ubs[jj] - lbs[jj]);
        ddata = XData_[ss*nInputs_+jj];
        if (ddata < lbs[jj]-diff || ddata > ubs[jj]+diff) break;
      }
      if (jj == nInputs_)
      {
        for (jj = 0; jj < nInputs_; jj++)
           XX[samCnt*nInputs_+jj] = XData_[ss*nInputs_+jj];
        YY[samCnt] = YData_[ss];
        samCnt++;
      }
    }
    if (outputLevel_ > 0)
    {
      printf("MGP3: number of partitions = %d\n", nPartitions_);
      printf("Partition %d has %d sample points.\n",ii+1,samCnt);
    }
    rsPtrs_[ii] = new GP3(nInputs_, samCnt);
    rsPtrs_[ii]->setOutputLevel(0);
    rsPtrs_[ii]->initialize(XX, YY);
    total += samCnt;
  }
  if (outputLevel_ >= 0)
    printf("Total number of sample points in all partitions = %d\n",total);
  delete [] vces;
  delete [] indices;
  delete [] XX;
  delete [] YY;
  turnPrintTSOn();

  if (outputLevel_ > 1) printf("MGP3 training completed.\n");
  if (psRSCodeGen_ == 1) 
    printf("MGP3 INFO: response surface stand-alone code not available.\n");
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MGP3::genNDGridData(double *XIn, double *YIn, int *N, double **X2, 
                        double **Y2)
{
  int    totPts;
  double *XX, *YY;

  initialize(XIn, YIn);
  if ((*N) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
  genNDGrid(N, &XX);
  if ((*N) == 0) return 0;
  totPts = (*N);

  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  YY = new double[totPts];
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*X2) = XX;
  (*Y2) = YY;
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int MGP3::gen1DGridData(double *XIn, double *YIn, int ind1,
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
  checkAllocate(XX, "XX in MGP3::gen1DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
     for (kk = 0; kk < nInputs_; kk++) 
        XX[ii*nInputs_+kk] = settings[kk]; 
     XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
     (*X2)[ii] = HX * ii + lowerBounds_[ind1];
  }
   
  YY = new double[totPts];
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int MGP3::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
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
  checkAllocate(*X2, "X2 in MGP3::gen2DGrid");
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
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int MGP3::gen3DGridData(double *XIn,double *YIn,int ind1,int ind2,int ind3,
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
  checkAllocate(*X2, "X2 in MGP3::gen3DGrid");
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
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int MGP3::gen4DGridData(double *XIn,double *YIn,int ind1,int ind2,int ind3,
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
  checkAllocate(*X2, "X2 in MGP3::gen4DGrid");
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
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, XX, YY);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*n) = totPts;
  (*Y2) = YY;
  delete [] XX;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MGP3::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y=0.0;
  interpolate(iOne, X, &Y, NULL);
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double MGP3::evaluatePoint(int npts, double *X, double *Y)
{
  interpolate(npts, X, Y, NULL);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double MGP3::evaluatePointFuzzy(double *X, double &std)
{
  int    iOne=1;
  double Y=0.0;
  interpolate(iOne, X, &Y, &std);
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double MGP3::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystds)
{
  interpolate(npts, X, Y, Ystds);
  return 0.0;
}

// ************************************************************************
// interpolation 
// ------------------------------------------------------------------------
int MGP3::interpolate(int npts, double *X, double *Y, double *Ystds)
{
  int    ii, pp, ss, count, highFlag, hasStds=0;
  double Yt, ddata, diff, *iRanges, *lbs, *ubs;

  if (Ystds != NULL) hasStds = 1;
  iRanges = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++)
     iRanges[ii] = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
  for (ss = 0; ss < npts; ss++)
  {
    count = 0;
    Yt = 0.0;
    for (pp = 0; pp < nPartitions_; pp++)
    {
      lbs = boxes_[pp]->lBounds_;
      ubs = boxes_[pp]->uBounds_;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (ubs[ii] == upperBounds_[ii]) highFlag = 1;
        else                             highFlag = 0;
        ddata = X[ss*nInputs_+ii];
        diff = 0.05 * (ubs[ii] - lbs[ii]);
        if (highFlag == 0)
        {
          if (ddata < (lbs[ii]-diff) || ddata >= (ubs[ii]+diff)) break;
        }
        else
        {
          if (ddata < (lbs[ii]-diff) || ddata > (ubs[ii]+diff)) break;
        }
      }
      if (ii == nInputs_)
      {
        if (hasStds) 
          Yt += rsPtrs_[pp]->evaluatePointFuzzy(&X[ss*nInputs_],Ystds[ss]);
        else
          Yt += rsPtrs_[pp]->evaluatePoint(&X[ss*nInputs_]);
        count++;
      }
    }
    if (count == 0)
    {
      printf("MGP3 evaluate ERROR: sample point outside range.\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %d = %e (in [%e, %e]?)\n",ii+1,
           X[ss*nInputs_+ii],lowerBounds_[ii],upperBounds_[ii]);
      return 0.0;
    }
    Y[ss] = Yt / (double) count;
  }
  delete [] iRanges;
  return 0;
}


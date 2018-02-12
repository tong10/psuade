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
// Functions for the class MMars (for large data set)
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MMars.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include "Psuade.h"
#include "MainEffectAnalyzer.h"

// ************************************************************************
// Constructor for object class Multi-domain Mars
// ------------------------------------------------------------------------
MMars::MMars(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    idata, ii;
  double ddata;
  char   pString[501], *strPtr, equal[100], winput[5000];

  boxes_ = NULL;
  partSize_ = 6000;

  // display banner and additonal information
  if (isScreenDumpModeOn() == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*   Multi-Multivariate Regression Function (MMARS) Analysis\n");
    printOutTS(PL_INFO,"* Set printlevel to 1-4 to see details.\n");
    printOutTS(PL_INFO,"* Turn on rs_expert mode to make changes.\n");
    printEquals(PL_INFO, 0);
  }
  
  // Xd contains the amount of overlaps between partitions
  Xd_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) Xd_[ii] = 0.05;
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("You can improve smoothness across partitions by allowing\n");
    printf("overlaps. The recommended overlap is 0.1 (or 10%%).\n");
    sprintf(pString, "Enter the degree of overlap (0 - 0.4) : ");
    ddata = getDouble(pString);
    if (ddata < 0 || ddata > 0.4)
    {
      ddata = 0.1;
      printf("ERROR: Degree of overlap should be > 0 and <= 0.4.\n");
      printf("INFO:  Degree of overlap set to default = 0.1\n");
    }
    for (ii = 0; ii < nInputs_; ii++) Xd_[ii] = ddata;
    printf("You can decide the sample size of each partition.\n");
    printf("Larger sample size per partition will take more setup time.\n");
    printf("The default is 1000 (will have more if there is overlap).\n");
    sprintf(pString, "Enter the partition sample size (1000 - 10000) : ");
    partSize_ = getInt(500, 20000, pString);
  }

  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MMARS_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &idata);
      if (idata >= 1000) partSize_ = idata;
      else 
      {
        printf("MMars INFO: config parameter setting not done.\n");
        printf("            max_samples_per_group %d too small.\n",idata);
        printf("            max_samples_per_group should be >= 1000.\n");
      }
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  ddata = log(1.0*nPartitions_) / log(2.0);
  idata = (int) ddata;
  if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
  printf("MMars: number of partitions = %d\n", nPartitions_);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MMars::~MMars()
{
  if (Xd_ != NULL) delete [] Xd_;
  if (boxes_ != NULL)
  {
    for (int ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii] != NULL)
      {
        if (boxes_[ii]->marsPtr_ != NULL) delete boxes_[ii]->marsPtr_;
        if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
        if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      }
    }
    delete [] boxes_;
  }
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MMars::initialize(double *X, double *Y)
{
  int    ii, jj, ss, nSubs, index, *indices, incr, samCnt;
  double range, *XX, *YY, *lbs, *ubs, *vces, *coefs, var=0,ddata, diff;

  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii] != NULL)
      {
        if (boxes_[ii]->marsPtr_ != NULL) delete boxes_[ii]->marsPtr_;
        if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
        if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      }
    }
    delete [] boxes_;
  }
  boxes_ = NULL;

  if (lowerBounds_ == NULL)
  {
    printOutTS(PL_ERROR,
         "MMars initialize ERROR: sample bounds not set yet.\n");
    return -1;
  }

  MainEffectAnalyzer *me = new MainEffectAnalyzer();
  turnPrintTSOff();
  vces = new double[nInputs_];
  me->computeVCECrude(nInputs_,nSamples_,X,Y,lowerBounds_,upperBounds_, 
                      var, vces);
  delete me;
  indices = new int[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) indices[ii] = ii;
  ddata = 0;
  for (ii = 0; ii < nInputs_; ii++) ddata += vces[ii];
  for (ii = 0; ii < nInputs_; ii++) vces[ii] /= ddata;
  sortDbleList2a(nInputs_, vces, indices);
  if (outputLevel_ > 1) 
  {
    for (ii = 0; ii < nInputs_; ii++) 
      printf("VCE %4d = %12.4e\n", indices[ii], vces[ii]);
  }

  nSubs = int(log(1.0*nPartitions_ + 1.0e-9) / log(2.0));
  boxes_ = new MMars_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MMars_Box();
    boxes_[ii]->lBounds_ = new double[nInputs_];
    boxes_[ii]->uBounds_ = new double[nInputs_];
    boxes_[ii]->marsPtr_ = NULL;
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
      printf("Selected input for partition %4dd = %4dd, VCE = %12.4e\n",
             ii+1, index+1, vces[nInputs_-1]);
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
  printf("MMars: Partition coverage check: %e (sum) ?= %e (orig)\n", 
         dcheck1, dcheck2);

  int total=0;
  YY = new double[nSamples_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    lbs = boxes_[ii]->lBounds_;
    ubs = boxes_[ii]->uBounds_;
    XX = new double[nSamples_*nInputs_];
    samCnt = 0;
    for (ss = 0; ss < nSamples_; ss++)
    {
      for (jj = 0; jj < nInputs_; jj++)
      {
        diff = Xd_[jj] * (ubs[jj] - lbs[jj]);
        ddata = X[ss*nInputs_+jj];
        if (ddata < lbs[jj]-diff || ddata > ubs[jj]+diff) break;
      } 
      if (jj == nInputs_)
      {
        for (jj = 0; jj < nInputs_; jj++)
           XX[samCnt*nInputs_+jj] = X[ss*nInputs_+jj];
        YY[samCnt] = Y[ss];
        samCnt++;
      }
    }
    if (outputLevel_ >= 0) 
      printf("Partition %d has %d sample points.\n",ii+1,samCnt);
    if (samCnt == 0)
    {
      printf("MMars INFO: some partition has no sample points.\n");
      boxes_[ii]->marsPtr_ = NULL;
    }
    else
    {
      boxes_[ii]->marsPtr_ = new Mars(nInputs_, samCnt);
      boxes_[ii]->marsPtr_->initialize(XX,YY);
    }
    boxes_[ii]->nSamples_ = samCnt;
    total += samCnt;
  }
  if (outputLevel_ > 0) 
  {
    printf("Sample size = %d\n",nSamples_);
    printf("Total sample sizes from all partitions = %d\n",total);
    printf("INFO: Total from all partitions may be larger than original\n");
    printf("      sample due to overlap. If the total is too large so\n");
    printf("      that it is close to the original size, partitioning\n");
    printf("      is not worthwhile -> you may want to reduce overlap.\n");
  }
  delete [] vces;
  delete [] indices;
  delete [] YY;
  turnPrintTSOn();
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::genNDGridData(double *X,double *Y,int *N2,double **X2,double **Y2)
{
  int totPts;

  initialize(X,Y);

  if ((*N2) == -999) return 0;
 
  genNDGrid(N2, X2);
  if ((*N2) == 0) return 0;
  totPts = (*N2);

  (*Y2) = new double[totPts];
  checkAllocate(*Y2, "Y2 in MMars::genNDGridData");
  evaluatePoint(totPts, *X2, *Y2);

  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                       int *N2, double **X2, double **Y2)
{
  int    ii, ss, totPts;
  double *XT, *XX, *YY, HX;

  initialize(X,Y);

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*X2) = new double[totPts];
  (*Y2) = new double[totPts];
  (*N2) = totPts;
  XX = (*X2);
  YY = (*Y2);

  XT = new double[totPts*nInputs_];
  checkAllocate(XT, "XT in MMars::gen1DGridData");
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
   
  for (ss = 0; ss < totPts; ss++) 
  {
    XT[ss*nInputs_+ind1]  = HX * ss + lowerBounds_[ind1];
    XX[ss] = HX * ss + lowerBounds_[ind1];
    YY[ss] = 0.0;
  }

  evaluatePoint(totPts, XT, YY);

  delete [] XT;
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *N2, double **X2, double **Y2)
{
  int    ii, ss, jj, index, totPts;
  double *XT, *XX, *YY, *HX;
 
  initialize(X,Y);

  totPts = nPtsPerDim_ * nPtsPerDim_;
  HX = new double[2];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*X2) = new double[2*totPts];
  (*Y2) = new double[totPts];
  (*N2) = totPts;
  XX = (*X2);
  YY = (*Y2);

  XT = new double[totPts*nInputs_];
  checkAllocate(XT, "XT in MMars::gen2DGridData");
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      index = ii * nPtsPerDim_ + jj;
      XT[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
      XT[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
      XX[index*2]   = HX[0] * ii + lowerBounds_[ind1];
      XX[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
    }
  }

  evaluatePoint(totPts, XT, YY);

  delete [] XT;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int MMars::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
                       double *settings, int *N2, double **X2, double **Y2)
{
  int    ii, ss, jj, ll, index, totPts;
  double *XT, *XX, *YY, *HX;

  initialize(X,Y);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX = new double[3];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*X2) = new double[3*totPts];
  (*Y2) = new double[totPts];
  (*N2) = totPts;
  XX = (*X2);
  YY = (*Y2);

  XT = new double[totPts*nInputs_];
  checkAllocate(XT, "XT in MMars::gen3DGridData");
  for (ss = 0; ss < totPts; ss++)
    for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii];

  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      for (ll = 0; ll < nPtsPerDim_; ll++)
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
        XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
        XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
        XX[index*3]   = HX[0] * ii + lowerBounds_[ind1];
        XX[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
        XX[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
      }
    }
  }

  evaluatePoint(totPts, XT, YY);

  delete [] XT;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int MMars::gen4DGridData(double *X,double *Y, int ind1, int ind2, int ind3, 
                       int ind4, double *settings, int *N2, double **X2, 
                       double **Y2)
{
  int    ii, ss, jj, ll, mm, index, totPts;
  double *XT, *XX, *YY, *HX;

  initialize(X,Y);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX = new double[4];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*X2) = new double[4*totPts];
  (*Y2) = new double[totPts];
  (*N2) = totPts;
  XX = (*X2);
  YY = (*Y2);

  XT = new double[totPts*nInputs_];
  checkAllocate(XT, "XT in MMars::gen4DGridData");
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      for (ll = 0; ll < nPtsPerDim_; ll++)
      {
        for (mm = 0; mm < nPtsPerDim_; mm++)
        {
          index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                  jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
          XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
          XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
          XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
          XT[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
          XX[index*4]   = HX[0] * ii + lowerBounds_[ind1];
          XX[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
          XX[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
          XX[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
        }
      }
    }
  }

  evaluatePoint(totPts, XT, YY);

  delete [] XT;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Evaluate a point
// ------------------------------------------------------------------------
double MMars::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y;
  evaluatePoint(1, X, &Y);
  return Y;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MMars::evaluatePoint(int nPts, double *X, double *Y)
{
  int    ss, pp, ii, nSamp, count, highFlag;
  double diff, Yt, ddata, *lbs, *ubs;

  for (ss = 0; ss < nPts; ss++) 
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
        diff = Xd_[ii] * (ubs[ii] - lbs[ii]);
        if (highFlag == 0)
        {
          if (ddata < (lbs[ii]-diff) || ddata >= (ubs[ii]+diff)) break;
        }
        else 
        {
          if (ddata < (lbs[ii]-diff) || ddata > (ubs[ii]+diff)) break;
        }
      } 
      if (ii == nInputs_ && boxes_[pp]->nSamples_ > 0)
      {
        if (boxes_[pp]->marsPtr_ != NULL)
        {
          Yt += boxes_[pp]->marsPtr_->evaluatePoint(&X[ss*nInputs_]);
          count++;
        }
      } 
    }
    if (count == 0)
    {
      printf("MMars evaluate WARNING: sample point not in any partition.\n");
      printf("INFO: this may happen during cross validation.\n");
      printf("INFO: prediction - average of all partitions (maybe wrong).\n");
      Yt = 0;
      for (pp = 0; pp < nPartitions_; pp++)
      {
        if (boxes_[pp]->marsPtr_ != NULL)
        {
          Yt += boxes_[pp]->marsPtr_->evaluatePoint(&X[ss*nInputs_]);
          count++;
        }
      }
    }
    Y[ss] = Yt / (double) count;
  }
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double MMars::evaluatePointFuzzy(double *X, double &std)
{
  int    iOne=1;
  double Y=0.0;
  evaluatePoint(iOne, X, &Y);
  std = 0.0;
  return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double MMars::evaluatePointFuzzy(int npts, double *X, double *Y, double *Ystd)
{
  evaluatePoint(npts, X, Y);
  for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
  return 0.0;
}


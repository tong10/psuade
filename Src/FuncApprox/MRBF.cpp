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
// Functions for the class RBF
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MRBF.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "MainEffectAnalyzer.h"
#include "PrintingTS.h"

extern "C" {
#if 0
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetrs_(char *,int *,int*,double*,int*,int*,double*,int*,int*);
#endif
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#define PS_RBF1

// ************************************************************************
// Constructor for object class RBF
// ------------------------------------------------------------------------
MRBF::MRBF(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    idata, ii;
  double ddata;
  char   pString[501], *strPtr, equal[100], winput[5000];

  type_ = 0;

  XNormalized_ = NULL;
  YNormalized_ = NULL;
  svdThresh_ = 1e-15;
  gaussScale_ = 1;
  boxes_ = NULL;
  partSize_ = 500;

  // display banner and additonal information
  if (isScreenDumpModeOn() == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
         "*   Multiple Radial Basis Function (RBF) Analysis\n");
    printOutTS(PL_INFO,"* Set printlevel to 1-4 to see RBF details.\n");
    printOutTS(PL_INFO,"* Default kernel    = multi-quadratic \n");
    printOutTS(PL_INFO,
         "* Default threshold = 1.0e-15 (for SVD truncation)\n");
    printOutTS(PL_INFO,"* Turn on rs_expert mode to make changes.\n");
    printEquals(PL_INFO, 0);
  }
  
  Xd_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) Xd_[ii] = 0.05;
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("In the following you have the option to select the kernel. \n");
    printf("0. multi-quadratic\n");
    printf("1. inverse multi-quadratic\n");
    printf("2. Gaussian\n");
    printf("3. thin plate spline\n");
    sprintf(pString,"Enter your choice (0 - 3) : ");
    type_ = getInt(0, 3, pString);
    if (type_ == 2)
    {
      sprintf(pString,
         "Enter scaling factor for Gaussian kernel (default=1) : ");
      gaussScale_ = getDouble(pString);
    }
    printOutTS(PL_INFO,
         "The RBF matrix to be constructed may be near-singular.\n");
    printOutTS(PL_INFO,
         "Currently, singular values < max(svd)*1e-15 are truncated.\n");
    printOutTS(PL_INFO,
         "You have the option to change this threshold (1e-15).\n");
    printOutTS(PL_INFO,
         "NOTE: truncating singular values may lead to erroneous results.\n");
    sprintf(pString, "Enter new threshold for SVD (> 0 but << 1) : ");
    svdThresh_ = getDouble(pString);
    printf("You can improve smoothness across partitions by allowing\n");
    printf("overlaps. The recommended overlap is 0.1 (or 10%%).\n");
    sprintf(pString, "Enter the degree of overlap (0 - 0.4) : ");
    ddata = getDouble(pString);
    if (ddata < 0 || ddata > 0.4)
    {
      ddata = 0.1;
      printf("Degree of overlap set to default = 0.1\n");
    }
    for (ii = 0; ii < nInputs_; ii++) Xd_[ii] = ddata;
    printf("You can decide the sample size of each partition.\n");
    printf("Larger sample size per partition will take more setup time.\n");
    printf("The default is 500 (will have more if there is overlap).\n");
    sprintf(pString, "Enter the partition sample size (500 - 2000) : ");
    partSize_ = getInt(500, 20000, pString);
  }

  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MRBF_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &idata);
      if (idata >= 100) partSize_ = idata;
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  ddata = log(1.0*nPartitions_) / log(2.0);
  idata = (int) ddata;
  if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
  printf("MRBF: number of partitions = %d\n", nPartitions_);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MRBF::~MRBF()
{
  if (XNormalized_ != NULL) delete [] XNormalized_;
  if (YNormalized_ != NULL) delete [] YNormalized_;
  if (Xd_          != NULL) delete [] Xd_;
  if (boxes_ != NULL)
  {
    for (int ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii]->regCoeffs_ != NULL) delete [] boxes_[ii]->regCoeffs_;
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      if (boxes_[ii]->X_ != NULL) delete [] boxes_[ii]->X_;
    }
    delete [] boxes_;
  }
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MRBF::initialize(double *X, double *Y)
{
  int    ii, jj, ss, nSubs, index, *indices, incr, samCnt;
  double range, *XX, *YY, *lbs, *ubs, *vces, *coefs, var=0,ddata, diff;

  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii]->regCoeffs_ != NULL) delete [] boxes_[ii]->regCoeffs_;
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      if (boxes_[ii]->X_ != NULL) delete [] boxes_[ii]->X_;
    }
    delete [] boxes_;
  }
  boxes_ = NULL;
  if (XNormalized_ != NULL) delete [] XNormalized_;
  if (YNormalized_ != NULL) delete [] YNormalized_;
  XNormalized_ = NULL;
  YNormalized_ = NULL;

  if (lowerBounds_ == NULL)
  {
    printOutTS(PL_ERROR,
         "RBF initialize ERROR: sample bounds not set yet.\n");
    return -1;
  }
  XNormalized_ = new double[nSamples_*nInputs_];
  checkAllocate(XNormalized_, "XNormalized_ in RBF::initialize");
  for (ii = 0; ii < nInputs_; ii++)
  {
    range = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
    for (ss = 0; ss < nSamples_; ss++)
      XNormalized_[ss*nInputs_+ii] = 
         (X[ss*nInputs_+ii] - lowerBounds_[ii]) * range;
  }
  YNormalized_ = new double[nSamples_];
  checkAllocate(YNormalized_, "YNormalized_ in RBF::initialize");
  initOutputScaling(Y, YNormalized_);
  for (ii = 0; ii < nSamples_; ii++) YNormalized_[ii] = Y[ii] - YMean_;

  MainEffectAnalyzer *me = new MainEffectAnalyzer();
  turnPrintTSOff();
  vces = new double[nInputs_];
  me->computeVCECrude(nInputs_,nSamples_,X,Y,lowerBounds_,upperBounds_, 
                      var, vces);
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
  boxes_ = new MRBF_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MRBF_Box();
    boxes_[ii]->lBounds_ = new double[nInputs_];
    boxes_[ii]->uBounds_ = new double[nInputs_];
    boxes_[ii]->regCoeffs_ = NULL;
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
  printf("Coverage check: %e (sum) ?= %e (orig)\n", dcheck1, dcheck2);

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
           XX[samCnt*nInputs_+jj] = XNormalized_[ss*nInputs_+jj];
        YY[samCnt] = YNormalized_[ss];
        samCnt++;
      }
    }
    if (outputLevel_ >= 0) 
      printf("Partition %d has %d sample points.\n",ii+1,samCnt);
    if (samCnt == 0)
    {
      printf("MRBF INFO: some partition has no sample points.\n");
      boxes_[ii]->regCoeffs_ = NULL;
    }
    else
    {
      initialize(samCnt, XX, YY, &coefs);
      boxes_[ii]->regCoeffs_ = coefs;
    }
    boxes_[ii]->X_ = XX;
    boxes_[ii]->nSamples_ = samCnt;
    total += samCnt;
  }
  if (outputLevel_ >= 0) 
    printf("Total number of sample points in all partitions = %d\n",total);
  delete [] vces;
  delete [] indices;
  delete [] YY;
  turnPrintTSOn();
  return 0;
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MRBF::initialize(int nSamples, double *X, double *Y, double **coefs)
{
  int    ii, kk, ss, ss2, nSamp1;
  double *Dmat, ddata, *regCoeffs;

#ifdef PS_RBF1
  nSamp1 = nSamples + 1;
#else
  nSamp1 = nSamples;
#endif
  Dmat = new double[nSamp1*nSamp1];
  checkAllocate(Dmat, "Dmat in RBF::initialize");
  switch(type_) 
  {
    case 0: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = multi-quadratic\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        Dmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          Dmat[ss*nSamp1+ss2] = 
                   Dmat[ss2*nSamp1+ss] = sqrt(ddata+1.0);
        }
      }
      break;

    case 1: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = inverse multi-quadratic\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        Dmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          Dmat[ss*nSamp1+ss2] = 
                Dmat[ss2*nSamp1+ss] = 1.0/sqrt(ddata+1.0);
        }
      }
      break;

    case 2: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = Gaussian\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        Dmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          Dmat[ss*nSamp1+ss2] = 
            Dmat[ss2*nSamp1+ss] = exp(-gaussScale_*ddata/2.0);
        }
      }
      break;

    case 3: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = thin plate spline\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        Dmat[ss*nSamp1+ss] = 0.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          Dmat[ss*nSamp1+ss2] = 
             Dmat[ss2*nSamp1+ss] = (ddata+1.0)*log(sqrt(ddata+1.0));
        }
      }
      break;
  }
#ifdef PS_RBF1
  for (ss = 0; ss < nSamples; ss++)
    Dmat[ss*nSamp1+nSamples] = Dmat[nSamples*nSamp1+ss] = 1.0; 
  Dmat[nSamples*nSamp1+nSamples] = 0.0;
#endif

  int    info, cnt=0;
  int    wlen = 5 * nSamp1;
  char   jobu = 'A', jobvt = 'A';
  double *SS = new double[nSamp1];
  double *UU = new double[nSamp1*nSamp1];
  double *VV = new double[nSamp1*nSamp1];
  double *WW = new double[wlen];
  checkAllocate(WW, "WW in RBF::initialize");
  dgesvd_(&jobu,&jobvt,&nSamp1,&nSamp1,Dmat,&nSamp1,SS,UU,&nSamp1,VV,
          &nSamp1,WW, &wlen,&info);
  if (info != 0) 
  {
    printOutTS(PL_WARN,"RBF ERROR: dgesvd returns error %d.\n",info);
    delete [] SS;
    delete [] UU;
    delete [] VV;
    delete [] WW;
    delete [] Dmat;
    return -1;
  }
  regCoeffs = new double[nSamp1];
  checkAllocate(regCoeffs, "regCoeffs in MRBF::initialize");
  for (ii = 0; ii < nSamples; ii++) regCoeffs[ii] = Y[ii];
#ifdef PS_RBF1
  regCoeffs[nSamples] = 0.0;
#endif
  for (ss = 1; ss < nSamp1; ss++)
  {
    if (SS[ss]/SS[0] < svdThresh_)
    {
      SS[ss] = 0;
      cnt++;
    }
  }
  if (cnt > 0 && psInteractive_ == 1 && outputLevel_ > 0) 
  {
    printOutTS(PL_WARN,
         "WARNING: RBF matrix is near-singular. Small singular values\n");
    printOutTS(PL_WARN,
         "         (%d out of %d) are truncated.\n",cnt,nSamp1);
    printOutTS(PL_WARN,"         Approximation may be inaccurate.\n");
  }
  for (ss = 0; ss < nSamp1; ss++)
  {
    WW[ss] = 0.0;
    for (ss2 = 0; ss2 < nSamp1; ss2++)
      WW[ss] += UU[ss*nSamp1+ss2] * regCoeffs[ss2];
  }
  for (ss = 0; ss < nSamp1; ss++) 
  {
    if (SS[ss] != 0) WW[ss] /= SS[ss];
    else             WW[ss] = 0;
  }
  for (ss = 0; ss < nSamp1; ss++)
  {
    regCoeffs[ss] = 0.0;
    for (ss2 = 0; ss2 < nSamp1; ss2++) 
      regCoeffs[ss] += VV[ss*nSamp1+ss2] * WW[ss2];
  }
  (*coefs) = regCoeffs;
  delete [] SS;
  delete [] UU;
  delete [] VV;
  delete [] WW;
  delete [] Dmat;
  return 0;
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MRBF::genNDGridData(double *X,double *Y,int *N2,double **X2,double **Y2)
{
   int totPts;

   initialize(X,Y);

   if ((*N2) == -999) return 0;
  
   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   (*Y2) = new double[totPts];
   checkAllocate(*Y2, "Y2 in MRBF::genNDGridData");
   evaluatePoint(totPts, *X2, *Y2);

   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MRBF::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
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
   checkAllocate(XT, "XT in MRBF::gen1DGridData");
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
int MRBF::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
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
   checkAllocate(XT, "XT in MRBF::gen2DGridData");
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
int MRBF::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
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
   checkAllocate(XT, "XT in MRBF::gen3DGridData");
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
int MRBF::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
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
   checkAllocate(XT, "XT in MRBF::gen4DGridData");
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
double MRBF::evaluatePoint(double *X)
{
   int    iOne=1;
   double Y;
   evaluatePoint(1, X, &Y);
   return Y;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MRBF::evaluatePoint(int nPts, double *X, double *Y)
{
  int    ss, ss2, ii, pp, nSamp, count, highFlag;
  double dist, diff, Yt, ddata, *iRanges, *lbs, *ubs, *XP;

  iRanges = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++)
     iRanges[ii] = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
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
        nSamp = boxes_[pp]->nSamples_;
        XP = boxes_[pp]->X_;
        for (ss2 = 0; ss2 < nSamp; ss2++) 
        {
          dist = 0.0;
          for (ii = 0; ii < nInputs_; ii++) 
          {
            ddata = X[ss*nInputs_+ii];
            ddata = (ddata - lowerBounds_[ii]) * iRanges[ii];
            ddata -= XP[ss2*nInputs_+ii];
            dist += ddata * ddata;
          }
          switch (type_)
          {
            case 0: dist = sqrt(dist + 1.0); break;
            case 1: dist = 1.0/sqrt(dist + 1.0); break;
            case 2: dist = exp(-0.5*dist*gaussScale_); break;
            case 3: dist = (dist+1)*log(sqrt(dist+1)); break;
          }
          Yt += dist * boxes_[pp]->regCoeffs_[ss2];
        }
        Yt += boxes_[pp]->regCoeffs_[nSamp];
        count++;
      } 
    }
    if (count == 0)
    {
      printf("MRBF evaluate ERROR: sample point outside range.\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %d = %e (in [%e, %e]?)\n",ii+1,
           X[ss*nInputs_+ii],lowerBounds_[ii],upperBounds_[ii]);
      printf("Sample prediction will be set to Ymean: may be wrong.\n");
      Yt = 0;
      count = 1;
    }
    Yt /= (double) count;
    Y[ss] = Yt + YMean_;
  }
  delete [] iRanges;
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double MRBF::evaluatePointFuzzy(double *X, double &std)
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
double MRBF::evaluatePointFuzzy(int npts, double *X, double *Y, double *Ystd)
{
   evaluatePoint(npts, X, Y);
   for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
   return 0.0;
}


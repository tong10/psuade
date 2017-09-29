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
// Functions for the class Regression
// (Up to fourth order only. Maybe not be stable for higher order)
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Psuade.h"
#include "Regression.h"
#include "PsuadeUtil.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C"
{
  void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
              double *, int *, double *, int *, double *, int *, int *);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
Regression::Regression(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
  faID_  = PSUADE_RS_REGR4;
  order_ = 4;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Regression::~Regression()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Regression::initialize(double *X, double *Y)
{
  int status;
 
  status = analyze(X, Y);
  if (status != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    return -1;
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int Regression::genNDGridData(double *X, double *Y, int *NN, double **XX, 
                              double **YY)
{
  int mm, totPts;

  if (initialize(X,Y) != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    (*NN) = 0;
    return -1;
  }

  if ((*NN) == -999) return 0;

  genNDGrid(NN, XX);
  if ((*NN) == 0) return 0;
  totPts = (*NN);

  (*YY) = new double[totPts];
  checkAllocate(*YY, "YY in Regression::genNDGrid");
  (*NN) = totPts;
  for (mm = 0; mm < totPts; mm++)
    (*YY)[mm] = evaluatePoint(&((*XX)[mm*nInputs_]));
  return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen1DGridData(double *X, double *Y, int ind1,
                              double *settings, int *NN, 
                              double **XX, double **YY)
{
  int    totPts, mm, nn;
  double HX, *Xloc;

  if (initialize(X,Y) != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    (*NN) = 0;
    return -1;
  }

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in Regression::gen1DGrid");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
   
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    Xloc[ind1] = HX * mm + lowerBounds_[ind1];
    (*XX)[mm] = Xloc[ind1];
    (*YY)[mm] = evaluatePoint(Xloc);
  }

  delete [] Xloc;
  return 0;
}

// ************************************************************************
// Generate 2D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen2DGridData(double *X, double *Y, int ind1,
                              int ind2, double *settings, int *NN, 
                              double **XX, double **YY)
{
  int    totPts, mm, nn, ind;
  double *HX, *Xloc;

  if (initialize(X,Y) != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    (*NN) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[2];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts * 2];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in Regression::gen2DGrid");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
   
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      ind = mm * nPtsPerDim_ + nn;
      Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
      Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
      (*XX)[ind*2]   = Xloc[ind1];
      (*XX)[ind*2+1] = Xloc[ind2];
      (*YY)[ind] = evaluatePoint(Xloc);
    }
  }

  delete [] Xloc;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen3DGridData(double *X, double *Y, int ind1,
                              int ind2, int ind3, double *settings, 
                              int *NN, double **XX, double **YY)
{
  int    totPts, mm, nn, pp, ind;
  double *HX, *Xloc;

  if (initialize(X,Y) != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    (*NN) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[3];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts * 3];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in Regression::gen3DGrid");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        ind = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
        Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
        Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
        (*XX)[ind*3]   = Xloc[ind1];
        (*XX)[ind*3+1] = Xloc[ind2];
        (*XX)[ind*3+2] = Xloc[ind3];
        (*YY)[ind] = evaluatePoint(Xloc);
      }
    }
  }

  delete [] Xloc;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int Regression::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                              int ind3, int ind4, double *settings, 
                              int *NN, double **XX, double **YY)
{
  int    totPts, mm, nn, pp, qq, ind;
  double *HX, *Xloc;

  if (initialize(X,Y) != 0)
  {
    printf("Regression: ERROR detected in regression analysis.\n");
    (*NN) = 0;
    return -1;
  }
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[4];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts * 4];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in Regression::gen4DGrid");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        for (qq = 0; qq < nPtsPerDim_; qq++)
        {
          ind = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                nn*nPtsPerDim_*nPtsPerDim_ + pp*nPtsPerDim_ + qq;
          Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
          Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
          Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
          Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
          (*XX)[ind*4]   = Xloc[ind1];
          (*XX)[ind*4+1] = Xloc[ind2];
          (*XX)[ind*4+2] = Xloc[ind3];
          (*XX)[ind*4+3] = Xloc[ind4];
          (*YY)[ind] = evaluatePoint(Xloc);
        }
      }
    }
  }

  delete [] Xloc;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double Regression::evaluatePoint(double *X)
{
  int    mm, nn, pp, qq, offset;
  double Xdata, Xdata2, Xdata3, Xdata4, Y;

  if (invCovMat_.nrows() <= 0)
  {
    printf("Regression ERROR: initialization has not been done.\n");
    exit(1);
  }

  Y = regCoeffs_[0];
  offset = 1;
  if (order_ >= 1)
  {
    for (mm = 0; mm < nInputs_; mm++)
    {
       Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
       Y += regCoeffs_[mm+1] * Xdata;
       offset++;
    }
  }
  if (order_ >= 2)
  {
    for (mm = 0; mm < nInputs_; mm++)
    {
      Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
      for (nn = mm; nn < nInputs_; nn++)
      {
        Xdata2 = (X[nn] - XMeans_[nn]) / XStds_[nn]; 
        Y += (regCoeffs_[offset++] * Xdata * Xdata2);
      }
    }
  }
  if (order_ >= 3)
  {
    for (mm = 0; mm < nInputs_; mm++)
    {
      Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
      for (nn = mm; nn < nInputs_; nn++)
      {
        Xdata2 = (X[nn] - XMeans_[nn]) / XStds_[nn]; 
        for (pp = nn; pp < nInputs_; pp++)
        {
           Xdata3 = (X[pp] - XMeans_[pp]) / XStds_[pp]; 
           Y += regCoeffs_[offset++] * Xdata * Xdata2 * Xdata3;
        }
      }
    }
  }
  if (order_ >= 4)
  {
    for (mm = 0; mm < nInputs_; mm++)
    {
      Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
      for (nn = mm; nn < nInputs_; nn++)
      {
        Xdata2 = (X[nn] - XMeans_[nn]) / XStds_[nn]; 
        for (pp = nn; pp < nInputs_; pp++)
        {
          Xdata3 = (X[pp] - XMeans_[pp]) / XStds_[pp]; 
          for (qq = pp; qq < nInputs_; qq++)
          {
            Xdata4 = (X[qq] - XMeans_[qq]) / XStds_[qq]; 
            Y += regCoeffs_[offset++] * Xdata * Xdata2 * Xdata3 * Xdata4;
          }
        }
      }
    }
  }
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Regression::evaluatePoint(int npts, double *X, double *Y)
{
  if (invCovMat_.nrows() <= 0)
  {
    printf("Regression ERROR: initialization has not been done.\n");
    exit(1);
  }
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double Regression::evaluatePointFuzzy(double *X, double &std)
{
  int    mm, nn, pp, qq, offset;
  double accum, *Xs, dtmp;
  double Xdata, Xdata2, Xdata3;

  if (invCovMat_.nrows() <= 0)
  {
    printf("Regression ERROR: initialization has not been done.\n");
    exit(1);
  }

  mm = invCovMat_.nrows();
  Xs = new double[mm];
  checkAllocate(Xs, "Xs in Regression::evaluatePointFuzzy");
  Xs[0] = 1.0;
  for (mm = 0; mm < nInputs_; mm++)
    Xs[mm+1] = (X[mm] - XMeans_[mm]) / XStds_[mm]; 

  accum = regCoeffs_[0];
  offset = 1;
  if (order_ >= 1)
  {
    for (mm = 1; mm <= nInputs_; mm++) 
    {
      accum += regCoeffs_[mm] * Xs[mm];
      offset++;
    }
    if (order_ >= 2)
    {
      for (mm = 1; mm <= nInputs_; mm++)
      {
        Xdata = Xs[mm];
        for (nn = mm; nn <= nInputs_; nn++)
        {
          accum += (regCoeffs_[offset] * Xdata * Xs[nn]);
          Xs[offset] = Xdata * Xs[nn];
          offset++;
        }
      }
    }
    if (order_ >= 3)
    {
      for (mm = 1; mm <= nInputs_; mm++)
      {
        Xdata = Xs[mm];
        for (nn = mm; nn <= nInputs_; nn++)
        {
          Xdata2 = Xs[nn];
          for (pp = nn; pp <= nInputs_; pp++)
          {
            accum += (regCoeffs_[offset] * Xdata * Xdata2 * Xs[pp]);
            Xs[offset] = Xdata * Xdata2 * Xs[pp];
            offset++;
          }
        }
      }
    }
    if (order_ >= 4)
    {
      for (mm = 1; mm <= nInputs_; mm++)
      {
        Xdata = Xs[mm];
        for (nn = mm; nn <= nInputs_; nn++)
        {
          Xdata2 = Xs[nn];
          for (pp = nn; pp <= nInputs_; pp++)
          {
            Xdata3 = Xs[pp];
            for (qq = pp; qq <= nInputs_; qq++)
            {
              accum += (regCoeffs_[offset]*Xdata*Xdata2*Xdata3*Xs[qq]);
              Xs[offset] = Xdata * Xdata2 * Xdata3 * Xs[qq];
              offset++;
            }
          }
        }
      }
    }
  }
  accum = accum * YStd_ + YMean_;

  std = 0.0;
  for (mm = 0; mm < offset; mm++)
  {
    dtmp = 0.0;
    for (nn = 0; nn < offset; nn++)
      dtmp += invCovMat_.getEntry(mm,nn) * Xs[nn];
    std += dtmp * Xs[mm];
  }
  std = sqrt(std) * YStd_;
  delete [] Xs;
  return accum;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double Regression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                      double *Ystd)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePointFuzzy(&(X[kk*nInputs_]), Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Regression::setParams(int targc, char **targv)
{
  order_ = *(int *) targv[0];
  if (order_ <= 0 || order_ > 4) order_ = 0;
  switch (order_)
  {
    case 0: faID_ = 0;               break;
    case 1: faID_ = PSUADE_RS_REGR1; break;
    case 2: faID_ = PSUADE_RS_REGR2; break;
    case 3: faID_ = PSUADE_RS_REGR3; break;
    case 4: faID_ = PSUADE_RS_REGR4; break;
  }
  return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int Regression::analyze(double *Xin, double *Y)
{
  int    M, N, ii, mm, nn, nn2, nn3, nn4, info, wlen, ind, last, NRevised;
  double *B, *XX, SSresid, SStotal, R2, var, *Bstd, *X;
  double esum, ymax, *txArray, *AA, *SS, *UU, *VV, *WW;
  char   jobu  = 'S', jobvt = 'A';
  char   pString[1000], response[1000];
  FILE   *fp;
  psMatrix eigMatT;
  psVector eigVals;

  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("Regression ERROR: nInputs or nSamples <= 0.\n");
    exit( 1 );
  } 
   
  txArray = new double[nSamples_];
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (mm = 0; mm < nSamples_; mm++)
      txArray[mm] = Xin[mm*nInputs_+ii];
    sortDbleList(nSamples_, txArray);
    last = 1;
    for (mm = 1; mm < nSamples_; mm++)
    {
      if (txArray[mm] != txArray[last-1])
      {
        txArray[last] = txArray[mm];
        last++;
      }
    }
    if (order_ >= last)
    {
      order_ = last - 1;
      printf("* Regression: order reduced to %d - not enough levels.\n",
             order_);
    }
  }
  delete [] txArray;

  if (outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    if (order_ == 0)
      printf("*               Constant Regression Analysis\n");
    else if (order_ == 1)
      printf("*               Linear Regression Analysis\n");
    else if (order_ == 2)
      printf("*             Quadratic Regression Analysis\n");
    else if (order_ == 3)
      printf("*               cubic Regression Analysis\n");
    else if (order_ == 4)
      printf("*               quartic Regression Analysis\n");
    printf("* R-squared gives a measure of the goodness of the model.\n");
    printf("* R-squared should be close to 1 if it is a good model.\n");
    printf("* TURN ON rs_expert mode to output regression matrix.\n");
    printf("* TURN ON rs_expert mode to output regression function.\n");
    printf("* SET print level to 4 to output data standard deviations.\n");
    printDashes(PL_INFO, 0);
    printf("* Suggestion: if your parameter ranges are too high, SCALE\n");
    printf("*             them first using 'irerange' command in PSUADE\n");
    printf("*             command line mode.\n");
    printEquals(PL_INFO, 0);
  }

  X = new double[nSamples_*nInputs_];
  checkAllocate(X, "X in Regression::analyze");
  if (psMasterMode_ == 1) 
  {
    printf("* Regression INFO: scaling turned off.\n");
    printf("*            To turn on scaling, use rs_expert mode.\n");
    initInputScaling(Xin, X, 0);
  }
  else initInputScaling(Xin, X, 1);

  N = loadXMatrix(X, &XX); 
  if (N == 0) return -1;
  if (N > nSamples_ && order_ >= 0)
  {
    if (XX != NULL) delete [] XX;
    printf("* Regression ERROR: sample too small for order %d.\n",order_);
    printf("                    Try lower order or larger sample.\n");
    return -1;
  }
  if (N > nSamples_)
  {
    if (XX != NULL) delete [] XX;
    printf("* Regression: sample size too small (%d > %d).\n",N,nSamples_);
    return -1;
  }
  M = nSamples_;

  wlen = 5 * M;
  AA = new double[M*N];
  UU = new double[M*N];
  SS = new double[N];
  VV = new double[N*N];
  WW = new double[wlen];
  checkAllocate(WW, "WW in Regression::analyze");

  for (mm = 0; mm < M; mm++) 
    for (nn = 0; nn < N; nn++) 
      AA[mm+nn*M] = sqrt(weights_[mm]) * XX[mm+nn*M];

  if (psMasterMode_ == 1)
  {
    sprintf(pString, "Store regression matrix? (y or n) ");
    getString(pString, response);
    if (response[0] == 'y')
    {
      fp = fopen("regression_matrix.m", "w");
      if(fp == NULL)
      {
         printf("fopen returned NULL in file %s line %d, exiting\n",
                __FILE__, __LINE__);
         exit(1);
      }
      fprintf(fp, "%% the sample matrix where svd is computed\n");
      fprintf(fp, "%% the last column is the right hand side\n");
      fprintf(fp, "%% B is the vector of coefficients\n");
      fprintf(fp, "AA = [\n");
      for (mm = 0; mm < M; mm++) 
      {
        for (nn = 0; nn < N; nn++) 
          fprintf(fp, "%16.6e ", AA[mm+nn*M]);
        fprintf(fp, "%16.6e \n",Y[mm]);
      }
      fprintf(fp, "];\n");
      fprintf(fp, "A = AA(:,1:%d);\n", N);
      fprintf(fp, "Y = AA(:,%d);\n", N+1);
      fprintf(fp, "B = A \\ Y;\n");
      fclose(fp);
      printf("Regression matrix now in regression_matrix.m\n");
    }
  }

  if (outputLevel_ > 3) printf("Running SVD ...\n"); 
  dgesvd_(&jobu, &jobvt, &M, &N, AA, &M, SS, UU, &M, VV, &N, WW,
          &wlen, &info);
  if (outputLevel_ > 3) 
    printf("SVD completed: status = %d (should be 0).\n",info); 

  if (info != 0)
  {
    printf("* Regression ERROR: dgesvd returns a nonzero (%d).\n",info);
    printf("* Regression terminates further processing.\n");
    printf("* To diagnose problem, re-run with rs_expert on.\n");
    delete [] XX;
    delete [] AA;
    delete [] UU;
    delete [] SS;
    delete [] VV;
    delete [] WW;
    return -1;
  }

  mm = 0;
  for (nn = 0; nn < N; nn++) if (SS[nn] < 0) mm++;
  if (mm > 0)
  {
    printf("* Regression WARNING: some of the singular values\n");
    printf("*            are < 0. May spell trouble but will.\n");
    printf("*            proceed anyway (%d).\n",mm);
  }
  if (SS[0] == 0.0) NRevised = 0;
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++) 
      if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8) NRevised--;
  }
  if (NRevised < N)
  {
    printf("* Regression WARNING: true rank of sample = %d (need %d)\n",
           NRevised, N);
    printf("* This can be due to the quality of the sample.\n");
  }
  if (psMasterMode_ == 1)
  {
    printf("Regression: matrix singular values\n");
    printf("The VERY small ones may cause poor numerical accuracy,\n");
    printf("but not keeping them may ruin the approximation power.\n");
    printf("So, select them judiciously.\n");
    for (nn = 0; nn < N; nn++) 
      printf("Singular value %5d = %e\n", nn+1, SS[nn]);
    sprintf(pString, "How many to keep (1 - %d, 0 - all) ? ", N); 
    NRevised = getInt(0,N,pString);
    if (NRevised == 0) NRevised = N;
    for (nn = NRevised; nn < N; nn++) SS[nn] = 0.0;
  }
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++) 
    {
      if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8)
      {
        SS[nn] = 0.0;
        NRevised--;
      }
    }
    if (NRevised != N) 
      printf("Regression INFO: %d singular values have been truncated.\n",
             N-NRevised);
  }

  for (mm = 0; mm < N; mm++) 
  {
    WW[mm] = 0.0;
    for (nn = 0; nn < M; nn++) 
      WW[mm] += UU[mm*M+nn] * sqrt(weights_[nn]) * Y[nn]; 
  }
  for (nn = 0; nn < NRevised; nn++) WW[nn] /= SS[nn];
  for (nn = NRevised; nn < N; nn++) WW[nn] = 0.0;
  B = new double[N];
  checkAllocate(WW, "B in Regression::analyze");
  for (mm = 0; mm < N; mm++) 
  {
    B[mm] = 0.0;
    for (nn = 0; nn < N; nn++) B[mm] += VV[mm*N+nn] * WW[nn]; 
  }

  eigMatT.load(N, N, VV);
  eigVals.load(N, SS);
  for (nn = 0; nn < N; nn++) eigVals[nn] = pow(eigVals[nn], 2.0);
  delete [] AA;
  delete [] SS;
  delete [] UU;
  delete [] VV;

  fp = NULL;
  if (psMasterMode_ == 1)
  {
    fp = fopen("regression_err_splots.m", "w");
    if(fp == NULL)
    {
      printf("fopen returned NULL in file %s line %d, exiting\n",
             __FILE__, __LINE__);
      exit(1);
    }
    fprintf(fp, "%% This file contains scatter plots of\n");
    fprintf(fp, "%% errors with respect to each input.\n");
    fprintf(fp, "%% This is useful to see which input fits the worst.\n");
    fprintf(fp,"X = [\n");
    for (mm = 0; mm < nSamples_; mm++)
    {
      for (nn = 0; nn < nInputs_; nn++) 
        fprintf(fp,"%16.8e ", X[mm*nInputs_+nn]);
      fprintf(fp,"\n");
    }
    fprintf(fp,"];\n");
    fprintf(fp,"R = [\n");
  }

  esum = ymax = 0.0;
  for (mm = 0; mm < M; mm++)
  {
    WW[mm] = 0.0;
    for (nn = 0; nn < N; nn++) WW[mm] += XX[mm+nn*nSamples_] * B[nn];
    WW[mm] -= Y[mm];
    esum = esum + WW[mm] * WW[mm] * weights_[mm];
    if (fp != NULL) 
      fprintf(fp, "%6d %24.16e\n", mm+1, WW[mm]*sqrt(weights_[mm]));
    if (PABS(Y[mm]) > ymax) ymax = PABS(Y[mm]);
  }
  esum /= (double) nSamples_;
  esum = sqrt(esum);
  if (outputLevel_ > 1)
    printf("* Regression: interpolation rms error = %11.4e (Ymax=%9.2e)\n",
           esum, ymax); 
  if (fp != NULL) 
  {
    fprintf(fp,"];\n");
    fprintf(fp,"for ii = 1 : size(X,2)\n");
    fprintf(fp,"   plot(X(:,ii), R(:,2), 'x')\n");
    fprintf(fp,"   title('Residual Plot')\n");
    fprintf(fp,"   disp(['Residual versus input ' int2str(ii)])\n");
    fprintf(fp,"   disp('Press Enter to continue')\n");
    fprintf(fp,"   pause\n");
    fprintf(fp,"end\n");
    fclose(fp);
    printf("FILE regression_err_splots.m contains fitting errors.\n");
    printf("     (resubstitution test) with respect to each input.\n");
  }

  computeSS(N, XX, Y, B, SSresid, SStotal);
  R2 = 1.0;
  if (SStotal != 0.0) R2  = 1.0 - SSresid / SStotal;
  if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
  else               var = 0.0;
  if (var < 0)
  { 
    if (PABS(var) > 1.0e-12)
    {
      printf("Regression WARNING: variance < 0.\n");
      printf("           Temporarily absolutize var (may have problems).\n");
      var = PABS(var);
    }
    else var = 0;
  }
  regCoeffs_.load(N, B);

  computeCoeffVariance(eigMatT, eigVals, var);
  Bstd = new double[N];
  for (ii = 0; ii < N; ii++)
    Bstd[ii] = sqrt(invCovMat_.getEntry(ii,ii));

  if (outputLevel_ >= 0)
  {
    printRC(N, B, Bstd, XX, Y);
    printf("* Regression R-square = %12.4e (SSresid,SStotal=%10.2e,%10.2e)\n",
           R2, SSresid, SStotal);
    if ((M - N - 1) > 0)
      printf("* adjusted   R-square = %12.4e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
    if (outputLevel_ > 1) printSRC(X, B, SStotal);
  }
  fp = NULL;
  if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.info", "w");
  if (fp != NULL)
  {
    fprintf(fp,"/* *************************************/\n");
    fprintf(fp,"/* Regression interpolator from PSUADE.*/\n");
    fprintf(fp,"/* ====================================*/\n");
    fprintf(fp,"/* This file contains information for interpolation\n");
    fprintf(fp,"   using response surface. Follow the steps below:\n");
    fprintf(fp,"   1. move this file to *.c file (e.g. main.c)\n");
    fprintf(fp,"   2. Compile main.c (cc -o main main.c -lm) \n");
    fprintf(fp,"   3. run: main input output\n");
    fprintf(fp,"          where input has the number of inputs and\n");
    fprintf(fp,"          the input values\n");
    fprintf(fp,"*/\n");
    fprintf(fp,"/* ====================================*/\n");
    fprintf(fp,"#include <math.h>\n");
    fprintf(fp,"#include <stdlib.h>\n");
    fprintf(fp,"#include <stdio.h>\n");
    fprintf(fp,"int interpolate(int,double*,double*,double*);\n");
    fprintf(fp,"main(int argc, char **argv) {\n");
    fprintf(fp,"  int    i, iOne=1, nInps;\n");
    fprintf(fp,"  double X[%d], Y, S;\n",nInputs_);
    fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
    fprintf(fp,"  if (argc != 3) {\n");
    fprintf(fp,"    printf(\"ERROR: not enough argument.\\n\");\n");
    fprintf(fp,"    exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
    fprintf(fp,"  if (fIn == NULL) {\n");
    fprintf(fp,"    printf(\"ERROR: cannot open input file.\\n\");\n");
    fprintf(fp,"    exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
    fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
    fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
    fprintf(fp,"    exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",
            nInputs_);
    fprintf(fp,"  fclose(fIn);\n");
    fprintf(fp,"  interpolate(iOne, X, &Y, &S);\n");
    fprintf(fp,"  printf(\"Y = %%e\\n\", Y);\n");
    fprintf(fp,"  printf(\"S = %%e\\n\", S);\n");
    fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
    fprintf(fp,"  if (fOut == NULL) {\n");
    fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
    fprintf(fp,"     exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
    fprintf(fp,"  fclose(fOut);\n");
    fprintf(fp,"}\n\n");
    fprintf(fp,"/* *************************************/\n");
    fprintf(fp,"/*  Regression interpolation function  */\n");
    fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
    fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
    fprintf(fp," * ... */\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"static double\n");
    fprintf(fp,"regCoefs[%d] = \n", N);
    fprintf(fp,"{\n");
    for (mm = 0; mm < N; mm++)
      fprintf(fp," %24.16e,\n", regCoeffs_[mm]);
    fprintf(fp,"};\n");
    fprintf(fp,"static double invCovMat[%d][%d] = \n", N, N);
    fprintf(fp,"{\n");
    for (mm = 0; mm < N; mm++)
    {
      fprintf(fp," { %24.16e", invCovMat_.getEntry(mm,0));
      for (nn = 1; nn < N; nn++)
        fprintf(fp,", %24.16e", invCovMat_.getEntry(mm,nn));
      fprintf(fp," },\n");
    }
    fprintf(fp,"};\n");
    fprintf(fp,"static int N=%d;\n",N);
    fprintf(fp,"/* ====================================*/\n");
    fprintf(fp,"int interpolate(int npts,double *X,double *Y,double *S){\n");
    fprintf(fp,"  int    ii, jj, kk, nInps=%d;\n",nInputs_);
    fprintf(fp,"  double y, *x, *x2, std, dtmp;\n");
    fprintf(fp,"  x2 = (double *) malloc(%d * sizeof(double));\n",
            N);
    fprintf(fp,"  for (ii = 0; ii < npts; ii++) {\n");
    fprintf(fp,"    x = &X[ii * %d];\n", nInputs_);
    fprintf(fp,"    y = regCoefs[0];\n");
    fprintf(fp,"    x2[0] = 1.0;\n");
    for (nn = 0; nn < nInputs_; nn++)
    {
      if (XMeans_[nn] != 0.0 || XStds_[nn] != 1.0)
        fprintf(fp,"    x2[%d] = (x[%d] - %e) / %e;\n", nn+1, nn,
                XMeans_[nn], XStds_[nn]);
      else
        fprintf(fp,"    x2[%d] = x[%d];\n", nn+1, nn);
    }
    for (nn = 1; nn <= nInputs_; nn++)
    {
      fprintf(fp,"    y += regCoefs[%d] * x2[%d];\n", nn, nn);
    }
    if (order_ >= 2)
    {
      ind = nInputs_ + 1;
      for (nn = 1; nn <= nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <= nInputs_; nn2++)
        {
          fprintf(fp, "    y += regCoefs[%d] * x2[%d] * x2[%d];\n",ind,
                  nn,nn2);
          fprintf(fp, "    x2[%d] = x2[%d] * x2[%d];\n",ind, nn,nn2);
          ind++;
        }
      }
    }
    if (order_ >= 3)
    {
      for (nn = 1; nn <= nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <= nInputs_; nn2++)
        {
          for (nn3 = nn2; nn3 <= nInputs_; nn3++)
          {
            fprintf(fp, "    y += regCoefs[%d] * x2[%d] * ",ind,nn);
            fprintf(fp, "x2[%d] * x2[%d];\n", nn2, nn3);
            fprintf(fp, "    x2[%d]=x2[%d]*x2[%d]*x2[%d];\n",ind,nn,nn2,nn3);
            ind++;
          }
        }
      }
    }
    if (order_ >= 4)
    {
      for (nn = 1; nn <= nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <= nInputs_; nn2++)
        {
          for (nn3 = nn2; nn3 <= nInputs_; nn3++)
          {
            for (nn4 = nn3; nn4 <= nInputs_; nn4++)
            {
              fprintf(fp, "    y += regCoefs[%d] * x2[%d] * ",ind,nn);
              fprintf(fp, "x2[%d] * x2[%d] * x2[%d];\n",nn2,nn3,nn4);
              fprintf(fp, "    x2[%d]=x2[%d]*x2[%d]*x2[%d]*x2[%d];\n",
                      ind,nn,nn2,nn3,nn4);
              ind++;
            }
          }
        }
      }
    }
    fprintf(fp,"    Y[ii] = y * %e + %e;\n",YStd_, YMean_);
    fprintf(fp,"    std = 0.0;\n");
    fprintf(fp,"    for (jj = 0; jj < N; jj++) {\n");
    fprintf(fp,"      dtmp = 0.0;\n");
    fprintf(fp,"      for (kk = 0; kk < N; kk++)\n");
    fprintf(fp,"        dtmp += invCovMat[jj][kk] * x2[kk];\n");
    fprintf(fp,"      std += dtmp * x2[jj];\n");
    fprintf(fp,"    }\n");
    fprintf(fp,"    std = sqrt(std);\n");
    fprintf(fp,"    S[ii] = std;\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  free(x2);\n");
    fprintf(fp,"  return 0;\n");
    fprintf(fp,"}\n\n");
    fprintf(fp,"/* ==============================================*/\n");
    printf("FILE psuade_rs.info contains the final polynomial\n");
    printf("     functional form.\n");
    fclose(fp);
  }

  fp = NULL;
  if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.py", "w");
  if (fp != NULL)
  {
    fwriteRSPythonHeader(fp);
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"# Regression interpolation\n");
    fprintf(fp,"#==================================================\n");
    fwriteRSPythonCommon(fp);
    fprintf(fp,"regCoefs = [\n");
    for (mm = 0; mm < N-1; mm++)
      fprintf(fp," %24.16e,\n", regCoeffs_[mm]);
    fprintf(fp," %24.16e ]\n", regCoeffs_[N-1]);
    fprintf(fp,"invCovMat = [\n");
    for (mm = 0; mm < N; mm++)
    {
      fprintf(fp," [ %24.16e", invCovMat_.getEntry(mm,0));
      for (nn = 1; nn < N; nn++)
        fprintf(fp,", %24.16e", invCovMat_.getEntry(mm,nn));
      fprintf(fp," ],\n");
    }
    fprintf(fp,"]\n");
    fprintf(fp,"###################################################\n");
    fprintf(fp,"# Regression interpolation function  \n");
    fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
    fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
    fprintf(fp,"# ... \n");
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"def interpolate(X): \n");
    fprintf(fp,"  nSamp = int(len(X) / %d)\n",nInputs_);
    fprintf(fp,"  Xt = %d * [0.0]\n", nInputs_);
    fprintf(fp,"  X2 = %d * [0.0]\n", N);
    fprintf(fp,"  Ys = 2 * nSamp * [0.0]\n");
    fprintf(fp,"  for ss in range(nSamp) : \n");
    fprintf(fp,"    for ii in range(%d) : \n", nInputs_);
    fprintf(fp,"      Xt[ii] = X[ss*%d+ii]\n",nInputs_);
    fprintf(fp,"    X2[0] = 1.0;\n");
    for (nn = 0; nn < nInputs_; nn++)
    {
      if (XMeans_[nn] != 0.0 || XStds_[nn] != 1.0)
           fprintf(fp,"    X2[%d] = (Xt[%d] - %e) / %e;\n", 
                   nn+1, nn, XMeans_[nn], XStds_[nn]);
      else fprintf(fp,"    X2[%d] = Xt[%d];\n", nn+1, nn);
    }
    fprintf(fp,"    Y = regCoefs[0]\n");
    if (order_ >= 1)
    {
      for (nn = 1; nn <= nInputs_; nn++)
        fprintf(fp,"    Y = Y + regCoefs[%d] * X2[%d];\n", nn, nn);
    }
    if (order_ >= 2)
    {
      ind = nInputs_ + 1;
      for (nn = 1; nn <= nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <=nInputs_; nn2++)
        {
          fprintf(fp, "    Y=Y+regCoefs[%d]*X2[%d]*X2[%d]\n",ind,nn,nn2);
          fprintf(fp, "    X2[%d]=X2[%d]*X2[%d]\n",ind,nn,nn2);
          ind++;
        }
      }
    }
    if (order_ >= 3)
    {
      for (nn = 1; nn <=nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <=nInputs_; nn2++)
        {
          for (nn3 = nn2; nn3 <= nInputs_; nn3++)
          {
            fprintf(fp, "    Y=Y+regCoefs[%d]*X2[%d]*X2[%d]*X2[%d]\n",
                    ind,nn,nn2,nn3);
            fprintf(fp, "    X2[%d]=X2[%d]*X2[%d]*X2[%d]\n",ind,nn,nn2,nn3);
            ind++;
          }
        }
      }
    }
    if (order_ >= 4)
    {
      for (nn = 1; nn <= nInputs_; nn++)
      {
        for (nn2 = nn; nn2 <= nInputs_; nn2++)
        {
          for (nn3 = nn2; nn3 <= nInputs_; nn3++)
          {
            for (nn4 = nn3; nn4 <= nInputs_; nn4++)
            {
              fprintf(fp, "    Y=Y+regCoefs[%d]*X2[%d]*X2[%d]",ind,nn,nn2);
              fprintf(fp, "*X2[%d]*X2[%d]\n", nn3, nn4);
              fprintf(fp, "    X2[%d]=X2[%d]*X2[%d]*X2[%d]*X2[nn4]\n",
                      ind,nn,nn2,nn3,nn4);
              ind++;
            }
          }
        }
      }
    }
    fprintf(fp,"    Ys[ss*2] = Y * %e + %e\n",YStd_, YMean_);
    fprintf(fp,"    std = 0.0\n");
    fprintf(fp,"    for jj in range(%d): \n", N);
    fprintf(fp,"      dtmp = 0.0\n");
    fprintf(fp,"      for kk in range(%d): \n", N);
    fprintf(fp,"        dtmp = dtmp + invCovMat[jj][kk] * X2[kk]\n");
    fprintf(fp,"      std = std + dtmp * X2[jj]\n");
    fprintf(fp,"    Ys[ss*2+1] = math.sqrt(std)\n");
    fprintf(fp,"  return Ys\n");
    fprintf(fp,"###################################################\n");
    fprintf(fp,"# main program\n");
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"infileName  = sys.argv[1]\n");
    fprintf(fp,"outfileName = sys.argv[2]\n");
    fprintf(fp,"inputs = getInputData(infileName)\n");
    fprintf(fp,"outputs = interpolate(inputs)\n");
    fprintf(fp,"genOutputFile(outfileName, outputs)\n");
    fprintf(fp,"###################################################\n");
    printf("FILE psuade_rs.py contains the final polynomial\n");
    printf("     functional form.\n");
    fclose(fp);
  }

  delete [] XX;
  delete [] WW;
  delete [] B;
  delete [] Bstd;
  return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int Regression::loadXMatrix(double *X, double **XXOut)
{
  int    M, N=0, mm, nn, nn2, nn3, nn4, ind, N2;
  double *XX=NULL;

  (*XXOut) = NULL;
  if (order_ > 4) return 0;

  M = nSamples_;
  if (order_ >= 0) N = 1;
  if (order_ >= 1)
  {
    N2 = N;
    N += nInputs_;
    if (nSamples_ < N)
    {
      N = N2;
      order_ = 0;
      printf("Regression INFO: order reduced to 0.\n");
    }
  }
  if (order_ >= 2)
  {
    N2 = N;
    N += nInputs_ * (nInputs_ + 1) / 2;
    if (nSamples_ < N)
    {
       N = N2;
       order_ = 1;
       printf("Regression INFO: order reduced to 1.\n");
    }
  }
  if (order_ >= 3)
  {
    N2 = N;
    for (nn = 0; nn < nInputs_; nn++)
      for (nn2 = nn; nn2 < nInputs_; nn2++)
        for (nn3 = nn2; nn3 < nInputs_; nn3++) N++;
    if (nSamples_ < N)
    {
      N = N2;
      order_ = 2;
      printf("Regression INFO: order reduced to 2.\n");
    }
  }
  if (order_ >= 4)
  {
    N2 = N;
    for (nn = 0; nn < nInputs_; nn++)
      for (nn2 = nn; nn2 < nInputs_; nn2++)
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
          for (nn4 = nn3; nn4 < nInputs_; nn4++) N++;
    if (nSamples_ < N)
    {
      N = N2;
      order_ = 3;
      printf("Regression INFO: order reduced to 3.\n");
    }
  }
  if (N > M) return N;
  XX = new double[M*N];
  checkAllocate(XX, " Regression XX");
  if (order_ >= 0)
  {
    for (mm = 0; mm < M; mm++) XX[mm] = 1.0;
  }
  if (order_ >= 1)
  {
    for (mm = 0; mm < M; mm++)
    {
      XX[mm] = 1.0;
      for (nn = 0; nn < nInputs_; nn++)
        XX[M*(nn+1)+mm] = X[mm*nInputs_+nn];
    }
  }
  if (order_ >= 2)
  {
    ind = nInputs_ + 1;
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        for (mm = 0; mm < M; mm++)
          XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2];
        ind++;
      }
    }
  }
  if (order_ >= 3)
  {
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (mm = 0; mm < M; mm++)
            XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2] * 
                           X[mm*nInputs_+nn3];
          ind++;
        }
      }
    }
  }
  if (order_ >= 4)
  {
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
             for (mm = 0; mm < M; mm++)
                XX[M*ind+mm] = X[mm*nInputs_+nn] * X[mm*nInputs_+nn2] * 
                               X[mm*nInputs_+nn3] * X[mm*nInputs_+nn4];
             ind++;
          }
        }
      }
    }
  }
  (*XXOut) = XX;
  return N;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int Regression::computeSS(int N, double *XX, double *Y,
                          double *B, double &SSresid, double &SStotal)
{
  int    nn, mm;
  double rdata, ymean, SSreg, ddata, SSresidCheck;

  SSresid = SSresidCheck = SStotal = SSreg = ymean = 0.0;
  for (mm = 0; mm < nSamples_; mm++) ymean += sqrt(weights_[mm]) * Y[mm];
  ymean /= (double) nSamples_;
  for (mm = 0; mm < nSamples_; mm++)
  {
    ddata = 0.0;
    for (nn = 0; nn < N; nn++) ddata += (XX[mm+nn*nSamples_] * B[nn]);
    rdata = Y[mm] - ddata;
    SSresidCheck += rdata * rdata * weights_[mm];
    SSresid += rdata * Y[mm] * weights_[mm];
    SSreg += (ddata - ymean) * (ddata - ymean);
  }
  for (mm = 0; mm < nSamples_; mm++)
    SStotal += weights_[mm] * (Y[mm] - ymean) * (Y[mm] - ymean);
  if (outputLevel_ > 0)
  {
    printf("* Regression: SStot  = %24.16e\n", SStotal);
    printf("* Regression: SSreg  = %24.16e\n", SSreg);
    printf("* Regression: SSres  = %24.16e\n", SSresid);
    printf("* Regression: SSres  = %24.16e (true)\n", SSresidCheck);
  }
  SSresid = SSresidCheck;
  if (outputLevel_ > 0 && nSamples_ != N)
  {
    printf("* Regression: eps(Y) = %24.16e\n",
           SSresidCheck/(nSamples_-N));
  }

  return 0;
}

// *************************************************************************
// compute coefficient variances (diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int Regression::computeCoeffVariance(psMatrix &eigMatT, psVector &eigVals,
                                     double var)
{
  int      ii, jj, nRows;
  double   invEig, dtmp;
  psMatrix tMat;

  nRows = eigMatT.nrows();
  tMat.setDim(nRows, nRows);

  for (ii = 0; ii < nRows; ii++)
  {
    invEig = eigVals[ii];
    if (invEig != 0.0) invEig = 1.0 / invEig;
    for (jj = 0; jj < nRows; jj++)
    {
      dtmp = invEig * eigMatT.getEntry(ii,jj) * var;
      tMat.setEntry(jj, ii, dtmp);
    }
  }
  eigMatT.matmult(tMat, invCovMat_);
  return 0;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int Regression::printRC(int N,double *B,double *Bvar,double *XX,double *Y)
{
  int    nn1, ii, ind, nn2, nn3, nn4;
  double coef, Bmax;
  char   fname[200];
  FILE   *fp;

  if (order_ < 0 || order_ > 4) return 0;
  printEquals(PL_INFO, 0);
  printf("*** Note: these coefficients may not be true coefficients due\n");
  printf("***       to sample matrix scaling (i.e., they may be scaled).\n");
  printDashes(PL_INFO, 0);
  printf("*            ");
  for (ii = 1; ii < order_; ii++) printf("    ");
  printf("  coefficient   std. error   t-value\n");
  printDashes(PL_INFO, 0);
  if (PABS(Bvar[0]) < 1.0e-15) coef = 0.0;
  else                         coef = B[0] / Bvar[0]; 
  printf("* Constant  ");
  for (ii = 1; ii < order_; ii++) printf("    ");
  printf("= %12.4e %12.4e %12.4e\n", B[0], Bvar[0], coef);
  if (order_ >= 1)
  {
    for (nn1 = 1; nn1 <= nInputs_; nn1++)
    {
      if (PABS(Bvar[nn1]) < 1.0e-15) coef = 0.0;
      else                           coef = B[nn1] / Bvar[nn1]; 
      {
        printf("* Input %3d ", nn1);
        for (ii = 1; ii < order_; ii++) printf("    ");
        printf("= %12.4e %12.4e %12.4e\n", B[nn1], Bvar[nn1], coef);
      }
      strcpy(fname, "dataVariance1");
    }
  }
  if (order_ >= 2)
  {
    Bmax = 0.0;
    for (nn1 = nInputs_+1; nn1 < N; nn1++)
      if (PABS(B[nn1]) > Bmax) Bmax = PABS(B[nn1]); 
    ind = nInputs_ + 1;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        if (PABS(B[ind]) > 1.0e-6 * Bmax) 
        {
          if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
          else coef = B[ind] / Bvar[ind]; 
          {
            printf("* Input %3d %3d ", nn1+1, nn2+1);
            for (ii = 2; ii < order_; ii++) printf("    ");
            printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],coef); 
          }
        }
        ind++;
      }
    }
    strcpy(fname, "dataVariance2");
  }
  if (order_ >= 3)
  {
    Bmax = 0.0;
    for (nn1 = ind; nn1 < N; nn1++)
      if (PABS(B[nn1]) > Bmax) Bmax = PABS(B[nn1]); 
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          if (PABS(B[ind]) > 1.0e-6 * Bmax) 
          {
            if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
            else coef = B[ind] / Bvar[ind]; 
            {
              printf("* Input %3d %3d %3d ", nn1+1, nn2+1, nn3+1);
              for (ii = 3; ii < order_; ii++) printf("    ");
              printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],
                     coef);
            }
          }
          ind++;
        }
      }
    }
    strcpy(fname, "dataVariance3");
  }
  if (order_ >= 4)
  {
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            if (PABS(B[ind]) > 1.0e-6 * Bmax) 
            {
              if (PABS(Bvar[ind]) < 1.0e-15) coef = 0.0;
              else coef = B[ind] / Bvar[ind]; 
              {
                printf("* Input %3d %3d %3d %3d ",nn1+1,nn2+1,nn3+1,
                       nn4+1);
                for (ii = 4; ii < order_; ii++) printf("    ");
                printf("= %12.4e %12.4e %12.4e\n",B[ind],Bvar[ind],
                       coef); 
              }
            }
            ind++;
          }
        }
      }
    }
    strcpy(fname, "dataVariance4");
  }
  printDashes(PL_INFO, 0);

  if (order_ >= 0 && order_ <= 4 && psMasterMode_ == 1)
  {
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      printf("fopen returned NULL in file %s line %d, exiting\n",
              __FILE__, __LINE__);
      exit(1);
    }
    fprintf(fp, "data number     data         standard error\n");
    for (ii = 0; ii < nSamples_; ii++)
    {
      coef = 0.0;
      for (nn1 = 0; nn1 < N; nn1++) 
        coef += PABS(XX[nn1*nSamples_+ii]*Bvar[nn1]);
      fprintf(fp,"%7d         %12.4e %12.4e\n",ii+1,Y[ii],sqrt(coef));
    }
    fclose(fp);
    printf("FILE %s contains output data standard errors.\n",fname);
  }
  return 0;
}

// *************************************************************************
// print standardized regression coefficients
// -------------------------------------------------------------------------
int Regression::printSRC(double *X, double *B, double SStotal)
{
  int    nn, mm, ind, ii, nn2, itemp, nn3, nn4, *iArray;
  double denom, xmean, coef, Bmax, coef1, coef2, xmean1, xmean2;
  double xmean3, coef3, xmean4, coef4, *B2;

  if (order_ < 0 || order_ > 4) return 0;

  printAsterisks(PL_INFO, 0);
  printf("*       Standardized Regression Coefficients (SRC)\n");
  printf("* When R-square is acceptable (order assumption holds), the\n");
  printf("* absolute values of SRCs provide variable importance.\n"); 
  printEquals(PL_INFO, 0);
  printf("* based on nSamples = %d\n", nSamples_);

  B2 = new double[nSamples_];
  checkAllocate(B2, " Regression B2");
  if (order_ >= 1)
  {
    denom = sqrt(SStotal / (double) (nSamples_ - 1));
    Bmax = 0.0;
    for (nn = 0; nn < nInputs_; nn++)
    {
      xmean = 0.0;
      for (mm = 0; mm < nSamples_; mm++) 
        xmean += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
      xmean /= (double) nSamples_;
      coef = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
      {
        coef += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean) * 
                (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean);
      }
      coef = sqrt(coef / (double) (nSamples_ - 1)) / denom;
      printf("* Input %3d ", nn+1);
      for (ii = 1; ii < order_; ii++) printf("    ");
      if ((B[nn+1]*coef) > Bmax) Bmax = B[nn+1]*coef;
      printf("= %12.4e\n", B[nn+1]*coef);
      B2[nn] = PABS(B[nn+1] * coef);
    }
    printDashes(PL_INFO, 0);
    printf("*    ordered list of SRCs\n");
    printDashes(PL_INFO, 0);
    iArray = new int[nInputs_];
    for (nn = 0; nn < nInputs_; nn++) iArray[nn] = nn;
    sortDbleList2a(nInputs_, B2, iArray); 
    for (nn = nInputs_-1; nn >= 0; nn--)
    {
      printf("* Input %3d ", iArray[nn]+1);
      for (ii = 1; ii < order_; ii++) printf("    ");
      printf("= %12.4e\n", B2[nn]);
    }
    delete [] iArray; 
  }
  if (order_ >= 2)
  {
    ind = nInputs_ + 1;
    for (nn = 0; nn < nInputs_; nn++)
    {
      xmean1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++) 
        xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
      xmean1 /= (double) nSamples_;
      coef1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                 (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
      coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        xmean2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
        xmean2 /= (double) nSamples_;
        coef2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                   (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
        coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
        B2[ind] = B[ind] * coef1 * coef2 / denom;
        if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
        ind++;
      }
    }
    ind = nInputs_ + 1;
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        if (PABS(B2[ind]) > 1.0e-3 * Bmax)
        {
          printf("* Input %3d,%3d ", nn+1, nn2+1);
          for (ii = 2; ii < order_; ii++) printf("    ");
          printf("= %12.4e\n", B2[ind]);
        }
        ind++;
      }
    }
  }
  if (order_ >= 3)
  {
    itemp = ind;
    for (nn = 0; nn < nInputs_; nn++)
    {
      xmean1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
      xmean1 /= (double) nSamples_;
      coef1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                 (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
      coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        xmean2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
        xmean2 /= (double) nSamples_;
        coef2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                   (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
        coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          xmean3 = 0.0;
          for (mm = 0; mm < nSamples_; mm++)
            xmean3 += X[mm*nInputs_+nn3] * sqrt(weights_[mm]);
          xmean3 /= (double) nSamples_;
          coef3 = 0.0;
          for (mm = 0; mm < nSamples_; mm++)
            coef3 += (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3) * 
                     (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3);
          coef3 = sqrt(coef3 / (double) (nSamples_ - 1));
          B2[ind] = B[ind] * coef1 * coef2 *coef3 / denom;
          if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
          ind++;
        }
      }
    }
    ind = itemp;
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          if (PABS(B2[ind]) > 1.0e-3 * Bmax)
          {
            printf("* Input %3d,%3d,%3d ", nn+1, nn2+1, nn3+1);
            for (ii = 3; ii < order_; ii++) printf("    ");
            printf("= %12.4e\n", B2[ind]);
          }
          ind++;
        }
      }
    }
  }
  if (order_ >= 4)
  {
    itemp = ind;
    for (nn = 0; nn < nInputs_; nn++)
    {
      xmean1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        xmean1 += X[mm*nInputs_+nn] * sqrt(weights_[mm]);
      xmean1 /= (double) nSamples_;
      coef1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        coef1 += (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1) * 
                 (sqrt(weights_[mm])*X[mm*nInputs_+nn] - xmean1);
      coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        xmean2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          xmean2 += X[mm*nInputs_+nn2] * sqrt(weights_[mm]);
        xmean2 /= (double) nSamples_;
        coef2 = 0.0;
        for (mm = 0; mm < nSamples_; mm++)
          coef2 += (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2) * 
                   (sqrt(weights_[mm])*X[mm*nInputs_+nn2] - xmean2);
        coef2 = sqrt(coef2 / (double) (nSamples_ - 1));
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          xmean3 = 0.0;
          for (mm = 0; mm < nSamples_; mm++)
            xmean3 += X[mm*nInputs_+nn3] * sqrt(weights_[mm]);
          xmean3 /= (double) nSamples_;
          coef3 = 0.0;
          for (mm = 0; mm < nSamples_; mm++)
            coef3 += (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3) * 
                     (sqrt(weights_[mm])*X[mm*nInputs_+nn3] - xmean3);
          coef3 = sqrt(coef3 / (double) (nSamples_ - 1));
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            xmean4 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
              xmean3 += X[mm*nInputs_+nn4] * sqrt(weights_[mm]);
            xmean4 /= (double) nSamples_;
            coef4 = 0.0;
            for (mm = 0; mm < nSamples_; mm++)
              coef4 += (sqrt(weights_[mm])*X[mm*nInputs_+nn4] - xmean4) * 
                       (sqrt(weights_[mm])*X[mm*nInputs_+nn4] - xmean4);
            coef4 = sqrt(coef4 / (double) (nSamples_ - 1));
            B2[ind] = B[ind] * coef1 * coef2 *coef3 *coef4 / denom;
            if (PABS(B2[ind]) > Bmax) Bmax = PABS(B2[ind]);
            ind++;
          }
        }
      }
    }
    ind = itemp;
    for (nn = 0; nn < nInputs_; nn++)
    {
      for (nn2 = nn; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            if (PABS(B2[ind]) > 1.0e-3 * Bmax)
            {
              printf("* Input %3d,%3d,%3d,%3d ",nn+1,nn2+1,nn3+1,nn4+1);
              for (ii = 4; ii < order_; ii++) printf("    ");
              printf("= %12.4e\n",B2[ind]);
            }
            ind++;
          }
        }
      }
    }
  }
  delete [] B2;
  printAsterisks(PL_INFO, 0);
  return 0;
}

// *************************************************************************
// print coefficients 
// -------------------------------------------------------------------------
int Regression::printCoefs(int N, double *B)
{
  int    ii, jj, nn1, nn2, nn3, nn4, **indexTable, ptr2, ptr3, ptr4, cnt;
  int    mm1, mm2, kk, *indices;
  double *trueCoefs, Bmax, ddata;

  indexTable = new int*[N];
  for (ii = 0; ii < N; ii++) indexTable[ii] = new int[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) indexTable[0][ii] = 0; 
  for (ii = 0; ii < nInputs_; ii++) 
  {
    for (jj = 0; jj < nInputs_; jj++) 
      indexTable[ii+1][jj] = 0; 
    indexTable[ii+1][ii] = 1; 
  }
  ptr2 = ptr3 = nInputs_ + 1;
  if (order_ >= 2)
  {
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (jj = 0; jj < nInputs_; jj++) 
          indexTable[ptr3][jj] = 0; 
        indexTable[ptr3][nn1]++; 
        indexTable[ptr3][nn2]++; 
        ptr3++;
      }
    }
  }
  ptr4 = ptr3;
  if (order_ >= 3)
  {
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (jj = 0; jj < nInputs_; jj++) 
            indexTable[ptr4][jj] = 0; 
          indexTable[ptr4][nn1]++; 
          indexTable[ptr4][nn2]++; 
          indexTable[ptr4][nn3]++; 
          ptr4++;
        }
      }
    }
  }
  cnt = ptr4;
  if (order_ >= 4)
  {
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            for (jj = 0; jj < nInputs_; jj++) 
              indexTable[cnt][jj] = 0; 
            indexTable[cnt][nn1]++; 
            indexTable[cnt][nn2]++; 
            indexTable[cnt][nn3]++; 
            indexTable[cnt][nn4]++; 
            cnt++;
          }
        }
      }
    }
  }
  for (ii = 0; ii < N; ii++)
  {
    for (jj = 0; jj < nInputs_; jj++)
      printf("%d ", indexTable[ii][jj]);
    printf("\n");
  }
    
  printEquals(PL_INFO, 0);
  printf("*** Note: these coefficients are true coefficients.\n");
  printDashes(PL_INFO, 0);
  printf("*            ");
  for (ii = 1; ii < order_; ii++) printf("    ");
  printf("  coefficient\n");
  printDashes(PL_INFO, 0);

  trueCoefs = new double[N];
  for (ii = 0; ii < N; ii++) trueCoefs[ii] = 0.0;
  trueCoefs[0] = B[0];
  indices = new int[nInputs_];
  if (order_ >= 1)
  {
    for (nn1 = 1; nn1 <= nInputs_; nn1++)
    {
      trueCoefs[nn1] = B[nn1] / XStds_[nn1-1];
      trueCoefs[0]  -= B[nn1] * XMeans_[nn1-1] / XStds_[nn1-1];
    }
  }
  if (order_ >= 2)
  {
    cnt = ptr2;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        trueCoefs[cnt] = B[cnt] / (XStds_[nn1] * XStds_[nn2]);
        trueCoefs[nn1+1] -= B[cnt]*XMeans_[nn2]/
                            (XStds_[nn1]*XStds_[nn2]);
        trueCoefs[nn2+1] -= B[cnt]*XMeans_[nn1]/
                            (XStds_[nn1]*XStds_[nn2]);
        trueCoefs[0] += B[cnt]*(XMeans_[nn1]*XMeans_[nn2])/
                               (XStds_[nn1]*XStds_[nn2]);
        cnt++;
      }
    }
  }
  if (order_ >= 3)
  {
    cnt = ptr3;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          ddata = B[cnt];
          ddata /= XStds_[nn1];
          ddata /= XStds_[nn2];
          ddata /= XStds_[nn3];
          trueCoefs[cnt] = ddata;

          ddata = XMeans_[nn1] / XStds_[nn1];
          ddata *= XMeans_[nn2] / XStds_[nn2];
          ddata *= XMeans_[nn3] / XStds_[nn3];
          trueCoefs[0] -= B[cnt] * ddata;

          ddata = 1.0 / XStds_[nn1];
          ddata *= XMeans_[nn2] / XStds_[nn2];
          ddata *= XMeans_[nn3] / XStds_[nn3];
          trueCoefs[nn1+1] += B[cnt] * ddata;

          ddata = 1.0 / XStds_[nn2];
          ddata *= XMeans_[nn1] / XStds_[nn1];
          ddata *= XMeans_[nn3] / XStds_[nn3];
          trueCoefs[nn2+1] += B[cnt] * ddata;

          ddata = 1.0 / XStds_[nn3];
          ddata *= XMeans_[nn1] / XStds_[nn1];
          ddata *= XMeans_[nn2] / XStds_[nn2];
          trueCoefs[nn3+1] += B[cnt] * ddata;

          ii = ptr2;
          for (mm1 = 0; mm1 < nInputs_; mm1++) indices[mm1] = 0;
          indices[nn1]++;
          indices[nn2]++;
          indices[nn3]++;
          
          for (mm1 = 0; mm1 < nInputs_; mm1++)
          {
             for (mm2 = mm1; mm2 < nInputs_; mm2++)
             {
                ddata = 1.0;
                if (indices[nn1] >= indexTable[ii][nn1] &&
                    indices[nn2] >= indexTable[ii][nn2])
                {
                   ddata = 1.0 / XStds_[nn1];
                   ddata /= XStds_[nn2];
                   ddata  *= XMeans_[nn3] / XStds_[nn3];
                   trueCoefs[ii] -= B[cnt] * ddata;
                   printf("(a) ii = %d\n",ii+1);
                }
                else if (indices[nn1] >= indexTable[ii][nn1] &&
                         indices[nn3] >= indexTable[ii][nn3])
                {
                   ddata = 1.0 / XStds_[nn1];
                   ddata /= XStds_[nn3];
                   ddata  *= XMeans_[nn2] / XStds_[nn2];
                   trueCoefs[ii] -= B[cnt] * ddata;
                   printf("(b) ii = %d\n",ii+1);
                }
                else if (indices[nn2] >= indexTable[ii][nn2] &&
                         indices[nn3] >= indexTable[ii][nn3])
                {
                   ddata = 1.0 / XStds_[nn2];
                   ddata /= XStds_[nn3];
                   ddata  *= XMeans_[nn1] / XStds_[nn1];
                   trueCoefs[ii] -= B[cnt] * ddata;
                   printf("(c) ii = %d\n",ii+1);
                }
                ii++;
             }
          }
          cnt++;
        }
      }
    }
  }
#if 0
  if (order_ >= 4)
  {
    cnt = ptr4;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            ddata = B[cnt];
            ddata /= XStds_[nn1];
            ddata /= XStds_[nn2];
            ddata /= XStds_[nn3];
            ddata /= XStds_[nn4];
            trueCoefs[cnt] = ddata;

            ddata = XMeans_[nn1] / XStds_[nn1];
            ddata *= XMeans_[nn2] / XStds_[nn2];
            ddata *= XMeans_[nn3] / XStds_[nn3];
            ddata *= XMeans_[nn4] / XStds_[nn4];
            trueCoefs[0] += B[cnt] * ddata;

            if (indexTable[cnt][0] = 2)
            {
               ddata = 1.0 / XStds_[nn1];
               ddata *= XMeans_[nn2] / XStds_[nn2];
               ddata *= XMeans_[nn3] / XStds_[nn3];
               ddata *= XMeans_[nn4] / XStds_[nn4];
               trueCoefs[nn1+1] -= B[cnt] * ddata;
            }

            if (indexTable[cnt][1] = 1)
            {
               ddata = 1.0 / XStds_[nn2];
               ddata *= XMeans_[nn1] / XStds_[nn1];
               ddata *= XMeans_[nn3] / XStds_[nn3];
               ddata *= XMeans_[nn4] / XStds_[nn4];
               trueCoefs[nn2+1] -= B[cnt] * ddata;
            }

            if (indexTable[cnt][2] == 1)
            {
               ddata = 1.0 / XStds_[nn3];
               ddata *= XMeans_[nn1] / XStds_[nn1];
               ddata *= XMeans_[nn2] / XStds_[nn2];
               ddata *= XMeans_[nn4] / XStds_[nn4];
               trueCoefs[nn3+1] -= B[cnt] * ddata;
            }

            if (indexTable[cnt][3] == 1)
            {
               ddata = 1.0 / XStds_[nn4];
               ddata *= XMeans_[nn1] / XStds_[nn1];
               ddata *= XMeans_[nn2] / XStds_[nn2];
               ddata *= XMeans_[nn3] / XStds_[nn3];
               trueCoefs[nn4+1] -= B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn1 != nn2) && (indexTable[ii][nn1] >= 1 && 
                   indexTable[ii][nn2] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn3] / XStds_[nn3];
               ddata  *= XMeans_[nn4] / XStds_[nn4];
               ddata  /= XStds_[nn1];
               ddata  /= XStds_[nn2];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn1 != nn3) && (indexTable[ii][nn1] == 1 && 
                   indexTable[ii][nn3] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn2] / XStds_[nn2];
               ddata  *= XMeans_[nn4] / XStds_[nn4];
               ddata  /= XStds_[nn1];
               ddata  /= XStds_[nn3];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn1 != nn4) && (indexTable[ii][nn1] >= 1 && 
                   indexTable[ii][nn4] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn2] / XStds_[nn2];
               ddata  *= XMeans_[nn3] / XStds_[nn3];
               ddata  /= XStds_[nn1];
               ddata  /= XStds_[nn4];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn2 != nn3) && (indexTable[ii][nn2] == 1 && 
                   indexTable[ii][nn3] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn1] / XStds_[nn1];
               ddata  *= XMeans_[nn4] / XStds_[nn4];
               ddata  /= XStds_[nn2];
               ddata  /= XStds_[nn3];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn2 != nn4) && (indexTable[ii][nn2] == 1 && 
                   indexTable[ii][nn4] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn1] / XStds_[nn1];
               ddata  *= XMeans_[nn3] / XStds_[nn3];
               ddata  /= XStds_[nn2];
               ddata  /= XStds_[nn4];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr2;
            while (ii < ptr3)
            {
               if ((nn3 != nn4) && (indexTable[ii][nn3] == 1 && 
                   indexTable[ii][nn4] == 1))
                  break;
               ii++;
            }
            if (ii < ptr3)
            {
               ddata  = XMeans_[nn1] / XStds_[nn1];
               ddata  *= XMeans_[nn2] / XStds_[nn2];
               ddata  /= XStds_[nn3];
               ddata  /= XStds_[nn4];
               trueCoefs[ii] += B[cnt] * ddata;
            }

            ii = ptr3;
            while (ii < N)
            {
               if (indexTable[ii][nn1] >= 1 && 
                   indexTable[ii][nn2] >= 1 && 
                   indexTable[ii][nn3] >= 1) 
                  break;
               ii++;
            }
            ddata  = XMeans_[nn4] / XStds_[nn4];
            ddata  /= XStds_[nn1];
            ddata  /= XStds_[nn2];
            ddata  /= XStds_[nn3];
            trueCoefs[ii] -= B[cnt] * ddata;

            ii = ptr3;
            while (ii < N)
            {
               if (indexTable[ii][nn1] >= 1 && 
                   indexTable[ii][nn2] >= 1 && 
                   indexTable[ii][nn4] >= 1) 
                  break;
               ii++;
            }
            ddata  = XMeans_[nn3] / XStds_[nn3];
            ddata  /= XStds_[nn1];
            ddata  /= XStds_[nn2];
            ddata  /= XStds_[nn4];
            trueCoefs[ii] -= B[cnt] * ddata;

            ii = ptr3;
            while (ii < N)
            {
               if (indexTable[ii][nn1] >= 1 && 
                   indexTable[ii][nn3] >= 1 && 
                   indexTable[ii][nn4] >= 1) 
                  break;
               ii++;
            }
            ddata  = XMeans_[nn2] / XStds_[nn2];
            ddata  /= XStds_[nn1];
            ddata  /= XStds_[nn3];
            ddata  /= XStds_[nn4];
            trueCoefs[ii] -= B[cnt] * ddata;

            ii = ptr3;
            while (ii < N)
            {
               if (indexTable[ii][nn2] >= 1 && 
                   indexTable[ii][nn3] >= 1 && 
                   indexTable[ii][nn4] >= 1) 
                  break;
               ii++;
            }
            ddata  = XMeans_[nn1] / XStds_[nn1];
            ddata  /= XStds_[nn2];
            ddata  /= XStds_[nn3];
            ddata  /= XStds_[nn4];
            trueCoefs[ii] -= B[cnt] * ddata;
            cnt++;
          }
        }
      }
    }
  }
#endif
  Bmax = trueCoefs[0];
  for (nn1 = 1; nn1 < N; nn1++)
    if (PABS(trueCoefs[nn1]) > Bmax) 
      Bmax = PABS(trueCoefs[nn1]);
  if (Bmax == 0) Bmax = 1;
  for (nn1 = 0; nn1 < N; nn1++)
    if (PABS(trueCoefs[nn1]/Bmax) < 1.0e-8)
      trueCoefs[nn1] = 0;

  printf("* Constant  ");
  for (ii = 1; ii < order_; ii++) printf("    ");
  printf("= %16.8e \n", trueCoefs[0]);
  if (order_ >= 1)
  {
    for (nn1 = 1; nn1 <= nInputs_; nn1++)
    {
      if (trueCoefs[nn1] != 0)
      {
        printf("* Input %3d ", nn1);
        for (ii = 1; ii < order_; ii++) printf("    ");
        printf("= %16.8e \n", trueCoefs[nn1]);
      }
    }
  }
  if (order_ >= 2)
  {
    cnt = nInputs_ + 1;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        if (trueCoefs[cnt] != 0)
        {
          printf("* Input %3d %3d ", nn1+1, nn2+1);
          for (ii = 2; ii < order_; ii++) printf("    ");
            printf("= %16.8e \n",trueCoefs[cnt]);
        }
        cnt++;
      }
    }
  }
  if (order_ >= 3)
  {
    cnt = ptr3;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          if (trueCoefs[cnt] != 0)
          {
            printf("* Input %3d %3d %3d ", nn1+1, nn2+1, nn3+1);
            for (ii = 3; ii < order_; ii++) printf("    ");
              printf("= %16.8e \n",trueCoefs[cnt]);
          }
          cnt++;
        }
      }
    }
  }
  if (order_ >= 4)
  {
    cnt = ptr4;
    for (nn1 = 0; nn1 < nInputs_; nn1++)
    {
      for (nn2 = nn1; nn2 < nInputs_; nn2++)
      {
        for (nn3 = nn2; nn3 < nInputs_; nn3++)
        {
          for (nn4 = nn3; nn4 < nInputs_; nn4++)
          {
            if (trueCoefs[cnt] != 0)
            {
              printf("* Input %3d %3d %3d %3d ",nn1+1,nn2+1,
                     nn3+1, nn4+1);
              for (ii = 4; ii < order_; ii++) printf("    ");
                printf("= %16.8e \n",trueCoefs[cnt]);
            }
            cnt++;
          }
        }
      }
    }
  }
  printDashes(PL_INFO, 0);
  delete [] trueCoefs;
  for (ii = 0; ii < N; ii++) delete [] indexTable[ii];
  delete [] indexTable;
  return 0;
}


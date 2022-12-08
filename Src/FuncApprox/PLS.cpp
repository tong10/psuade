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
// Functions for the class PLS
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Psuade.h"
#include "PLS.h"
#include "PDFManager.h"
#include "PsuadeUtil.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C"
{
   void dgels_(char *, int *, int *, int *, double *, int *,
               double *, int *, double *, int *, int *);
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetri_(int *, double *, int *, int *, double*, int *, int *);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
PLS::PLS(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
  faID_ = PSUADE_RS_PLS;
  numTerms_  = 0;
  initialized_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PLS::~PLS()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PLS::initialize(double *X, double *Y)
{
  //**/ ---------------------------------------------------------------
  //**/ print header
  //**/ ---------------------------------------------------------------
  if (psConfig_.InteractiveIsOn() && outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*     Partial Least Squares Linear Regression Analysis\n");
    printf("* R-squared gives a measure of the goodness of the model.\n");
    printf("* R-squared should be close to 1 if it is a good model.\n");
    printf("* TURN ON rs_expert mode to output orthonormal matrix.\n");
    printDashes(PL_INFO, 0);
  }

  //**/ ---------------------------------------------------------------
  //**/ preliminary error checking
  //**/ ---------------------------------------------------------------
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("PLS initialize ERROR: consult PSUADE developers.\n");
    exit( 1 );
  } 
  return analyze(X,Y);
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int PLS::genNDGridData(double *X, double *Y, int *NN, double **XX, 
                       double **YY)
{
  int mm, totPts, status=0;
  psVector vecYOut;

  //**/ ---------------------------------------------------------------
  //**/ initialization
  //**/ ---------------------------------------------------------------
  if (initialized_ == 0)
  {
    status = initialize(X,Y);
    if (status != 0)
    {
      printf("PLS: ERROR detected in regression analysis.\n");
      (*NN) = 0;
      return -1;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ return if there is no request to create lattice points
  //**/ ---------------------------------------------------------------
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
  vecYOut.setLength(totPts);
  (*YY) = vecYOut.takeDVector();
  (*NN) = totPts;
  for (mm = 0; mm < totPts; mm++)
    (*YY)[mm] = evaluatePoint(&((*XX)[mm*nInputs_]));
  return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                       int *NN, double **XX, double **YY)
{
  int mm, nn, status=0;

  //**/ ---------------------------------------------------------------
  //**/ initialization
  //**/ ---------------------------------------------------------------
  if (initialized_ == 0)
  {
    status = initialize(X,Y);
    if (status != 0)
    {
      printf("PLS: ERROR detected in regression analysis.\n");
      (*NN) = 0;
      return -1;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  int totPts = nPtsPerDim_;
  double HX = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1); 

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  psVector vecXOut, vecYOut, vecXT;
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
// Generate 2D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *NN, double **XX, double **YY)
{
  int mm, nn, ind, status=0;
  psVector vecHX, vecXT, vecXOut, vecYOut;

  //**/ ---------------------------------------------------------------
  //**/ initialization
  //**/ ---------------------------------------------------------------
  if (initialized_ == 0)
  {
    status = initialize(X,Y);
    if (status != 0)
    {
      printf("PLS: ERROR detected in regression analysis.\n");
      (*NN) = 0;
      return -1;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  int totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  vecXOut.setLength(totPts * 2);
  vecYOut.setLength(totPts);
  vecXT.setLength(nInputs_);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      ind = mm * nPtsPerDim_ + nn;
      vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
      vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
      (*XX)[ind*2]   = vecXT[ind1];
      (*XX)[ind*2+1] = vecXT[ind2];
      (*YY)[ind] = evaluatePoint(vecXT.getDVector());
    }
  }
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
                       double *settings, int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, ind, status=0;
  psVector vecHX, vecXT, vecXOut, vecYOut;

  //**/ ---------------------------------------------------------------
  //**/ initialization
  //**/ ---------------------------------------------------------------
  if (initialized_ == 0)
  {
    status = initialize(X,Y);
    if (status != 0)
    {
      printf("PLS: ERROR detected in regression analysis.\n");
      (*NN) = 0;
      return -1;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3])/(nPtsPerDim_ - 1); 

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  vecXOut.setLength(totPts * 3);
  vecYOut.setLength(totPts);
  vecXT.setLength(nInputs_);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        ind = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
        vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
        vecXT[ind3] = vecHX[2] * pp + VecLBs_[ind3];
        (*XX)[ind*3]   = vecXT[ind1];
        (*XX)[ind*3+1] = vecXT[ind2];
        (*XX)[ind*3+2] = vecXT[ind3];
        (*YY)[ind] = evaluatePoint(vecXT.getDVector());
      }
    }
  }
  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                              int ind3, int ind4, double *settings, 
                              int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, qq, ind, status=0;
  psVector vecHX, vecXT, vecXOut, vecYOut;

  //**/ ---------------------------------------------------------------
  //**/ initialization
  //**/ ---------------------------------------------------------------
  if (initialized_ == 0)
  {
    status = initialize(X,Y);
    if (status != 0)
    {
      printf("PLS: ERROR detected in regression analysis.\n");
      (*NN) = 0;
      return -1;
    }
  }
 
  //**/ ---------------------------------------------------------------
  //**/ set up for generating regular grid data
  //**/ ---------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3])/(nPtsPerDim_ - 1); 
  vecHX[3] = (VecUBs_[ind4] - VecLBs_[ind4])/(nPtsPerDim_ - 1); 

  //**/ ---------------------------------------------------------------
  //**/ allocate storage for and then generate the data points
  //**/ ---------------------------------------------------------------
  (*NN) = totPts;
  vecXOut.setLength(totPts * 4);
  vecYOut.setLength(totPts);
  vecXT.setLength(nInputs_);
  (*XX) = vecXOut.takeDVector();
  (*YY) = vecYOut.takeDVector();
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
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
          vecXT[ind1] = vecHX[0] * mm + VecLBs_[ind1];
          vecXT[ind2] = vecHX[1] * nn + VecLBs_[ind2];
          vecXT[ind3] = vecHX[2] * pp + VecLBs_[ind3];
          vecXT[ind4] = vecHX[3] * qq + VecLBs_[ind4];
          (*XX)[ind*4]   = vecXT[ind1];
          (*XX)[ind*4+1] = vecXT[ind2];
          (*XX)[ind*4+2] = vecXT[ind3];
          (*XX)[ind*4+3] = vecXT[ind4];
          (*YY)[ind] = evaluatePoint(vecXT.getDVector());
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double PLS::evaluatePoint(double *X)
{
  int    mm;
  double Xdata, Y=0.0;

  if (VecRegCoeffs_.length() == 0)
  {
    printf("PLS ERROR: need to call initialize first.\n");
    return 0.0;
  }
  for (mm = 0; mm < nInputs_; mm++)
  {
    Xdata = (X[mm] - VecXMeans_[mm]) / VecXStds_[mm]; 
    Y += VecRegCoeffs_[mm] * Xdata;
  }
  Y = Y * YStd_ + YMean_;
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double PLS::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double PLS::evaluatePointFuzzy(double *X, double &std)
{
  int    mm, nn;
  double YOut, coef, Xdata, dtmp;

  //**/ if no coefficient, flag error
  if (VecRegCoeffs_.length() == 0)
  {
    printf("PLS ERROR: initialize has not been called.\n");
    exit(1);
  }
 
  //**/ compute the mean
  YOut = 0.0;
  for (mm = 0; mm < nInputs_; mm++)
  {
    coef = VecRegCoeffs_[mm];
    Xdata = (X[mm] - VecXMeans_[mm]) / VecXStds_[mm]; 
    YOut += coef * Xdata;
  }
  YOut = YOut * YStd_ + YMean_;

  //**/ compute the standard deviation
  std = 0.0;
  for (mm = 0; mm < nInputs_; mm++)
  {
    dtmp = 0.0;
    for (nn = 0; nn < nInputs_; nn++)
    {
      Xdata = (X[nn] - VecXMeans_[nn]) / VecXStds_[nn]; 
      dtmp += invCovMat_.getEntry(mm,nn) * Xdata;
    }
    Xdata = (X[mm] - VecXMeans_[mm]) / VecXStds_[mm]; 
    std += dtmp * Xdata;
  }
  std = sqrt(std) * YStd_;
  return YOut;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double PLS::evaluatePointFuzzy(int npts, double *X, double *Y,double *Ystd)
{
  //**/ if no coefficient, flag error
  if (VecRegCoeffs_.length() == 0)
  {
     printf("PLS ERROR: initialize has not been called.\n");
     exit(1);
  }
  //**/ evaluate one point at a time
  for (int kk = 0; kk < npts; kk++)
     Y[kk] = evaluatePointFuzzy(&(X[kk*nInputs_]), Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int PLS::analyze(double *Xin, double *Y)
{
  psVector VecX, VecY;
  VecX.load(nSamples_*nInputs_, Xin);
  VecY.load(nSamples_, Y);
  return analyze(VecX, VecY);
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int PLS::analyze(psVector &vecXIn, psVector &vecYIn)
{
  int    ii, jj, kk, mm, nn, info, status=0;
  double ddata, snorm, esum, ymax, SSresid, SStotal, R2, var;
   
  //**/ =================================================================
  //**/ initial error checking
  //**/ =================================================================
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("UserRegression ERROR: nInputs or nSamples <= 0.\n");
    exit( 1 );
  }

  //**/ =================================================================
  //**/ display banner 
  //**/ =================================================================
  if (psConfig_.InteractiveIsOn() && outputLevel_ >= 1)
  {
    printAsterisks(PL_INFO, 0);
    printf("*          User Regression Analysis\n");
    printf("* R-square gives a measure of the goodness of the model.\n");
    printf("* R-square should be close to 1 if it is a good model.\n");
    printEquals(PL_INFO, 0);
  }

  //**/ ---------------------------------------------------------------
  //**/ optional scaling of the sample matrix
  //**/ ---------------------------------------------------------------
  psVector vecXT, vecYT;
  vecXT.setLength(nSamples_ * nInputs_);
  initInputScaling(vecXIn.getDVector(), vecXT.getDVector(), 1);
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) ddata += vecYIn[ii];
  YMean_ = ddata / (double) nSamples_;
  for (ii = 0; ii < nSamples_; ii++) vecYIn[ii] = vecYIn[ii] - YMean_;
  YStd_ = 1.0;

  //**/ ---------------------------------------------------------------
  //**/ compute S = X' * Y
  //**/ ---------------------------------------------------------------
  psVector vecS;
  vecS.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
  {
    ddata = 0.0;
    for (jj = 0; jj < nSamples_; jj++)
      ddata += vecXT[jj*nInputs_+ii] * vecYIn[jj];
    vecS[ii] = ddata;
  }
  for (ii = 0; ii < nSamples_; ii++) vecYIn[ii] += YMean_;

  //**/ ---------------------------------------------------------------
  //**/ loop 
  //**/ ---------------------------------------------------------------
  snorm = 0.0;
  for (jj = 0; jj < nInputs_; jj++) snorm += vecS[jj] * vecS[jj]; 
  snorm = sqrt(snorm);
  if (snorm == 0.0) 
  {
    printf("PLS ERROR: null correlation matrix.\n");
    status = 1;
    return status;
  }
  psVector vecR, vecT, vecP, vecD;
  vecR.setLength(nInputs_*nInputs_);
  vecT.setLength(nSamples_*nInputs_);
  vecP.setLength(nInputs_*nInputs_);
  vecD.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
  {
    //**/ r = S / norm(S)
    for (jj = 0; jj < nInputs_; jj++)
      vecR[ii*nInputs_+jj] = vecS[jj] / snorm; 
    //**/ t = X * r
    for (kk = 0; kk < nSamples_; kk++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nInputs_; jj++)
        ddata += vecXT[jj+kk*nInputs_] * vecR[ii*nInputs_+jj];
      vecT[ii*nSamples_+kk] = ddata;
    }
    //**/ p = X' * t
    for (kk = 0; kk < nInputs_; kk++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
        ddata += vecXT[kk+jj*nInputs_] * vecT[ii*nSamples_+jj];
        vecP[ii*nInputs_+kk] = ddata;
    }
    //**/ orthonormalize p
    for (jj = 0; jj < ii; jj++)
    {
      ddata = 0.0;
      for (kk = 0; kk < nInputs_; kk++)
        ddata = ddata + vecP[ii*nInputs_+kk] * vecP[jj*nInputs_+kk];
      for (kk = 0; kk < nInputs_; kk++)
        vecP[ii*nInputs_+kk] -= ddata * vecP[jj*nInputs_+kk];
    }
    ddata = 0.0;
    for (jj = 0; jj < nInputs_; jj++)
      ddata += vecP[ii*nInputs_+jj] * vecP[ii*nInputs_+jj]; 
    ddata = sqrt(ddata);
    if (ddata == 0.0) 
    {
      printf("PLS ERROR: rank deficient P matrix.\n");
      status = 1;
      break;
    }
    for (jj = 0; jj < nInputs_; jj++) vecP[ii*nInputs_+jj] /= ddata; 
    //**/ normalize T
    ddata = 0.0;
    for (jj = 0; jj < nSamples_; jj++)
      ddata += vecT[ii*nSamples_+jj] * vecT[ii*nSamples_+jj]; 
    ddata = sqrt(ddata);
    if (ddata == 0.0) 
    {
      printf("PLS ERROR: null T vector.\n");
      status = 1;
      break;
    }
    for (jj = 0; jj < nSamples_; jj++) vecT[ii*nSamples_+jj] /= ddata; 
    //**/ normalize R wrt T
    for (jj = 0; jj < nInputs_; jj++) vecR[ii*nInputs_+jj] /= ddata; 
    vecD[ii] = ddata;
    //**/ deflate S
    ddata = 0.0; 
    for (jj = 0; jj < nInputs_; jj++) 
      ddata += vecP[ii*nInputs_+jj] * vecS[jj];
    for (jj = 0; jj < nInputs_; jj++) 
      vecS[jj] = vecS[jj] - vecP[ii*nInputs_+jj] * ddata;

    //**/ compute norm of S
    snorm = 0.0;
    for (jj = 0; jj < nInputs_; jj++) snorm += vecS[jj] * vecS[jj]; 
    snorm = sqrt(snorm);
    if (snorm < 1.0e-8) break;
  }
  if (ii < nInputs_) ii++;
  numTerms_ = ii;
  //**/ At the end, X = T * D * P' and B = R T' Y
  //**/ or T is the basis
  snorm = 0.0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      for (kk = 0; kk < numTerms_; kk++)
        ddata = vecT[kk*nSamples_+jj] * vecP[kk*nInputs_+ii] * vecD[kk];
      ddata = ddata - vecXT[jj*nInputs_+ii];
      snorm += ddata * ddata;
    }
  }
  snorm = sqrt(snorm);
  if (psConfig_.InteractiveIsOn() && outputLevel_ >= 1)
    printf("PLS residual norm (||X-TDP^t||) = %e\n", snorm);
  if (status == 1) 
  {
    printf("PLS ERROR: failed.\n");
    status = 1;
    return status;
  }

  //**/ ---------------------------------------------------------------
  //**/ compute B
  //**/ ---------------------------------------------------------------
  if (psConfig_.InteractiveIsOn() && outputLevel_ > 1)
  {
    printf("Basis set T (so that X = T D P') is:\n");
    for (jj = 0; jj < nSamples_; jj++)
    {
      for (ii = 0; ii < numTerms_; ii++)
        printf("%12.4e ", vecT[ii*nSamples_+jj]);
      printf("\n");
    }
    printf("Given new x, transform to reduced space via x*P*D\n");
    printf("where x is 1xn\n");
    printf("      P is nxp  (p << n is the reduced input size)\n");
    printf("      D is pxp\n");
  }
  psVector vecB;
  vecB.setLength(nInputs_);
  //**/ P = T' * Y
  for (ii = 0; ii < numTerms_; ii++)
  {
    ddata = 0.0;
    for (jj = 0; jj < nSamples_; jj++)
      ddata += vecT[jj+ii*nSamples_] * vecYIn[jj];
    vecP[ii] = ddata;
  }
  //**/ B = R * P
  for (ii = 0; ii < nInputs_; ii++)
  {
    ddata = 0.0;
    for (jj = 0; jj < nInputs_; jj++)
      ddata += vecR[jj*nInputs_+ii] * vecP[jj];
    vecB[ii] = ddata;
  }
  VecRegCoeffs_ = vecB;
  initialized_ = 1;
   
  //**/ ---------------------------------------------------------------
  //**/ compute residual
  //**/ ---------------------------------------------------------------
  psVector vecW;
  vecW.setLength(nSamples_);
  esum = ymax = 0.0;
  for (ii = 0; ii < nSamples_; ii++)
  {
    vecW[ii] = 0.0;
    for (jj = 0; jj < nInputs_; jj++) 
      vecW[ii] += vecXT[ii*nInputs_+jj] * vecB[jj];
    vecW[ii] -= vecYIn[ii];
    esum = esum + vecW[ii] * vecW[ii];
    if (PABS(vecYIn[ii]) > ymax) ymax = PABS(vecYIn[ii]);
  }
  esum /= (double) nSamples_;
  esum = sqrt(esum);
  if (psConfig_.InteractiveIsOn() && outputLevel_ > 1)
    printf("* PLS: rms of interpolation error = %11.4e (Ymax=%9.2e)\n",
           esum, ymax); 

  //**/ ------------- form compute SS statistics ------------------
  computeSS(vecXT, vecYIn, vecB, SSresid, SStotal);
  R2 = 1.0;
  //**/ 1 - SSE / SST = (SST - SSE) / SST 
  //**/ MSE = SSE / (nSamples - N)
  //**/ MSR = SSR = SST - SSE
  //**/ F   = MSR / MSE = (SST - SSE) / SSE * (nSamples - N)
  if (SStotal != 0.0) R2  = 1.0 - SSresid / SStotal;
  if (nSamples_ > numTerms_)
       var = SSresid / (double) (nSamples_ - numTerms_);
  else var = 0.0;
  if (var < 0)
  { 
    if (PABS(var) > 1.0e-12) printf("PLS WARNING: var < 0.\n");
    else var = 0;
  }

  //**/ ------------- find variance of each coefficient ------------------
  psVector vecA, eigVals;
  psMatrix MatA, MatU, MatV, eigMatT;
  vecA.setLength(nSamples_*nInputs_);
  for (mm = 0; mm < nSamples_; mm++)
    for (nn = 0; nn < nInputs_; nn++)
      vecA[mm+nn*nSamples_] = vecXT[mm*nInputs_+nn];
  MatA.load(nSamples_, nInputs_, vecA.getDVector());
  if (psConfig_.InteractiveIsOn() && outputLevel_ > 3)
    printf("Running SVD ...\n");
  status = MatA.computeSVD(MatU, vecS, MatV);
  if (status != 0)
  {
    printf("* PLS Info: dgesvd returns a nonzero (%d).\n",status);
    printf("* PLS terminates further processing.\n");
    return -1;
  }
  eigMatT.load(nInputs_, nInputs_, MatV.getMatrix1D());
  eigVals.load(vecS);
  for (nn = 0; nn < nInputs_; nn++) eigVals[nn] = pow(eigVals[nn], 2.0);
  computeCoeffVariance(eigMatT, eigVals, var);
  VecRegStdevs_.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++)
    VecRegStdevs_[nn] = sqrt(invCovMat_.getEntry(nn,nn));

  //**/ --------- print out regression coefficients --------------
  if (psConfig_.InteractiveIsOn() && outputLevel_ > 0)
  {
    if (outputLevel_ > 0) printRC(vecB, VecRegStdevs_);
    printf("* PLS R-square = %12.4e (SSresid, SStotal =%10.2e,%10.2e)\n",
           R2, SSresid, SStotal);
    if ((nSamples_ - numTerms_ - 1) > 0)
      printf("* Adjusted R2  = %12.4e\n",
             1.0-(1.0-R2)*((nSamples_-1)/(nSamples_-numTerms_-1)));
    if (SSresid > 0)
      printf("* F-statistics = %12.4e (>6 for fit to be significant)\n",
             (SStotal-SSresid)/SSresid*(nSamples_-numTerms_));
    printEquals(PL_INFO, 0);
  }
  return 0;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int PLS::computeSS(psVector vecX, psVector vecY, psVector vecB, 
                   double &SSresid, double &SStotal)
{
  int    N, nn, mm;
  double rdata, ymean, SSreg, ddata;

  SSresid = SStotal = SSreg = ymean = 0.0;
  for (mm = 0; mm < nSamples_; mm++) ymean += vecY[mm];
  ymean /= (double) nSamples_;
  N = vecB.length();
  for (mm = 0; mm < nSamples_; mm++)
  {
    ddata = 0.0;
    for (nn = 0; nn < N; nn++) ddata += (vecX[mm*N+nn] * vecB[nn]);
    ddata += YMean_;
    rdata = vecY[mm] - ddata;
    SSresid += rdata * rdata;
    SSreg += (ddata - ymean) * (ddata - ymean);
  }
  for (mm = 0; mm < nSamples_; mm++)
    SStotal += (vecY[mm] - ymean) * (vecY[mm] - ymean);

  if (psConfig_.InteractiveIsOn() && outputLevel_ > 0)
  {
    //**/ SST
    printf("* PLS: SStot  = %16.8e\n", SStotal);
    //**/ SSE
    printf("* PLS: SSres  = %16.8e\n", SSresid);
    printf("* PLS: SSreg  = %16.8e\n", SSreg);
  }
  if (psConfig_.InteractiveIsOn() && outputLevel_ > 0 &&
      nSamples_ != N)
    printf("* PLS: MSE of residual = %16.8e\n",SSresid/(nSamples_-N));
  return 0;
}

// *************************************************************************
// compute variance matrix
// -------------------------------------------------------------------------
int PLS::computeCoeffVariance(psMatrix eigMat, psVector eigVals, double var)
{
  int      ii, jj, nRows;
  double   invEig, dtmp;
  psMatrix tMat;

  nRows = eigMat.nrows();
  tMat.setDim(nRows, nRows);

  //**/ compute Sigma^{-2} * V^T 
  for (ii = 0; ii < nRows; ii++)
  {
    invEig = eigVals[ii];
    if (invEig != 0.0) invEig = 1.0 / invEig;
    for (jj = 0; jj < nRows; jj++)
    {
      dtmp = invEig * eigMat.getEntry(ii,jj) * var;
      tMat.setEntry(jj, ii, dtmp);
    }
  }

  //**/ compute V * Sigma^{-2} * V^T 
  eigMat.matmult(tMat, invCovMat_);
  return 0;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int PLS::printRC(psVector vecB, psVector vecBStd)
{
  int    nn1;
  double coef, scaled;

  printEquals(PL_INFO, 0);
  scaled = 0.0;
  for (nn1 = 0; nn1 < nInputs_; nn1++) scaled += PABS(VecXMeans_[nn1]);
  if (scaled == 0)
  {
    for (nn1 = 0; nn1 < nInputs_; nn1++) if (VecXStds_[nn1] != 1) break;
    if (nn1 != nInputs_) scaled = 1;
  }
  if (scaled != 0)
  {
    printf("* NOTE: The coefficients below have been scaled.\n");
    printf("*       See above for input scaling information.\n");
  }
  printDashes(PL_INFO, 0);
  printf("*            ");
  printf("  coefficient   std. error   t-value\n");
  printDashes(PL_INFO, 0);
  //**/ ------ print out the scaled linear coefficients -----------
  printf("* Constant   = %12.4e\n", YMean_);
  for (nn1 = 0; nn1 < nInputs_; nn1++)
  {
    if (PABS(vecBStd[nn1]) < 1.0e-15) coef = 0.0;
    else                              coef = vecB[nn1] / vecBStd[nn1]; 
    //**/ if (PABS(coef) > 1.0)
    {
      printf("* Input  %3d ", nn1+1);
      printf("= %12.4e %12.4e %12.4e\n",vecB[nn1],vecBStd[nn1],coef);
    }
  }
  printEquals(PL_INFO, 0);
  return 0;
}


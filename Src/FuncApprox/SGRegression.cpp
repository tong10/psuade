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
// Functions for the class SparseGridRegression
// AUTHOR : CHARLES TONG
// DATE   : 2011
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SGRegression.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Globals.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
SparseGridRegression::SparseGridRegression(int nInputs,int nSamples):
                                           FuncApprox(nInputs,nSamples)
{
  int    ii, kk;
  double ddata, ddata2;
  FILE   *fp;

  faID_     = PSUADE_RS_REGSG;
  numPerms_ = 0;

  if (outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*       Sparse Grid Regression Analysis\n");
    printDashes(PL_INFO, 0);
    printf("* Note: This function looks for a ps_sparse_grid_info file\n");
    printf("*       for regression information.\n");
    printf("The file should be in the format:\n");
    printf(" line 1 : <nSamples> <nInputs> <pOrder>\n");
    printf(" line 2 : 1 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf(" line 3 : 2 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf("          ...\n");
    printEquals(PL_INFO, 0);
  }
  fp = fopen("ps_sparse_grid_info", "r");
  if (fp == NULL)
  {
    printf("SparseGridRegression ERROR: ps_sparse_grid_info file not found.\n");
    printf("This file is used to specify information needed for \n");
    printf("sparse grid regression. The file has the format:\n");
    printf(" line 1 : <nSamples> <nInputs> <pOrder>\n");
    printf(" line 2 : 1 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf(" line 3 : 2 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf("          ...\n");
    exit(1);
  }
  fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &pOrder_);
  printf("SparseGridRegression INFO: polynomial order = %d.\n", pOrder_);
  if (nSamples != nSamples_ || nInputs != nInputs_)
  {
    printf("SparseGridRegression ERROR: nSamples or nInputs does not match.\n");
    printf("                            SparseGrid is rigid in its sample.\n");
    fclose(fp);
    return;
  }
  sampleInputs_.setLength(nSamples_*nInputs_);
  sampleWeights_.setLength(nSamples_);
  for (ii = 0; ii < nSamples_; ii++)
  {
    for (kk = 0; kk < nInputs_; kk++)
    {
      fscanf(fp, "%lg", &ddata);
      sampleInputs_[ii*nInputs_+kk] = ddata;
    }
    fscanf(fp, "%lg", &ddata);
    sampleWeights_[ii] = ddata;
  }
  VecLBs_.setLength(nInputs_);
  VecUBs_.setLength(nInputs_);
  for (kk = 0; kk < nInputs_; kk++)
  {
    fscanf(fp, "%lg %lg", &ddata, &ddata2);
    VecLBs_[kk] = ddata;
    VecUBs_[kk] = ddata2;
    if (VecLBs_[kk] >= VecUBs_[kk])
    {
      printf("SparseGridRegression ERROR: invalid input bounds.\n");
      printf("       lbound, ubound (input %d) = %e %e\n", kk+1, 
             VecLBs_[kk], VecUBs_[kk]);
      exit(1);
    }
  }
  fclose(fp); 
  GenPermutations();
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SparseGridRegression::~SparseGridRegression()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int SparseGridRegression::initialize(double *X, double *Y)
{
  int totPts, ss;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if (sampleInputs_.length() == 0)
  {
    printf("SparseGridRegression::initialize ERROR - invalid sample.\n");
    return -1;
  }
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::initialize ERROR - invalid argument.\n");
    exit(1);
  } 
   
  //**/ ===============================================================
  //**/ analyze the data
  //**/ ===============================================================
  analyze(X, Y);
  if (psConfig_.RSCodeGenIsOn())
  {
    printf("SparseGridRegression INFO: response surface stand-alone ");
    printf("code not available.\n");
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int SparseGridRegression::genNDGridData(double *X, double *Y, int *N2,
                                        double **XOut, double **YOut)
{
  int totPts, ss;
  psVector vecYOut;

  //**/ ===============================================================
  //**/ error checking
  //**/ ===============================================================
  if (sampleInputs_.length() == 0)
  {
    printf("SparseGridRegression::genNDGridData ERROR - invalid sample.\n");
    return -1;
  }
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::genNDGridData ERROR - invalid argument.\n");
    exit(1);
  } 
   
  //**/ ===============================================================
  //**/ analyze the data
  //**/ ===============================================================
  analyze(X, Y);

  //**/ ===============================================================
  //**/ return if there is no request to create lattice points
  //**/ ===============================================================
  if ((*N2) == -999) return 0;

  //**/ ===============================================================
  //**/ generating regular grid data
  //**/ ===============================================================
  genNDGrid(N2, XOut);
  if ((*N2) == 0) return 0;
  totPts = (*N2);

  //**/ ===============================================================
  //**/ allocate storage for the data points and generate them
  //**/ ===============================================================
  vecYOut.setLength(totPts);
  (*YOut) = vecYOut.takeDVector();
  for (ss = 0; ss < totPts; ss++)
    (*YOut)[ss] = evaluatePoint(&((*XOut)[ss*nInputs_]));
  return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int SparseGridRegression::gen1DGridData(double *X, double *Y, int ind1,
                                        double *settings, int *NN, 
                                        double **XX, double **YY)
{
  int    totPts, mm, nn;
  double HX;
  psVector vecYOut, vecXT, vecXOut;
  
  //**/ ===============================================================
  //**/ create response surface
  //**/ ===============================================================
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

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
int SparseGridRegression::gen2DGridData(double *X, double *Y, int ind1,
                                      int ind2, double *settings, int *NN, 
                                      double **XX, double **YY)
{
  int totPts, mm, nn, index;
  psVector vecYOut, vecXT, vecXOut, vecHX;

  //**/ ===============================================================
  //**/ create response surface
  //**/ ===============================================================
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

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
int SparseGridRegression::gen3DGridData(double *X, double *Y, int ind1,
                                  int ind2, int ind3, double *settings, 
                                  int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, index;
  psVector vecYOut, vecXT, vecXOut, vecHX;

  //**/ ===============================================================
  //**/ create response surface
  //**/ ===============================================================
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

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
int SparseGridRegression::gen4DGridData(double *X, double *Y, int ind1, 
                             int ind2, int ind3, int ind4, double *settings, 
                             int *NN, double **XX, double **YY)
{
  int totPts, mm, nn, pp, qq, index;
  psVector vecYOut, vecXT, vecXOut, vecHX;

  //**/ ===============================================================
  //**/ create response surface
  //**/ ===============================================================
  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

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
double SparseGridRegression::evaluatePoint(double *X)
{
  int    ii, kk;
  double Y, multiplier, ddata;
  psMatrix MatLTable;

  if (VecRegCoefs_.length() == 0) return 0.0;
  Y = 0.0;
  MatLTable.setFormat(PS_MAT2D);
  MatLTable.setDim(nInputs_, pOrder_+1);
  double **LTable = MatLTable.getMatrix2D();
  Y = 0.0;
  for (kk = 0; kk < numPerms_; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = (X[ii] - VecLBs_[ii])/(VecUBs_[ii] - VecLBs_[ii])*2-1;
      EvalLegendrePolynomials(ddata, LTable[ii]);
    }
    multiplier = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      multiplier *= LTable[ii][MatPcePerms_.getEntry(kk,ii)];
    Y += VecRegCoefs_[kk] * multiplier;
  }
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePointFuzzy(double *X, double &std)
{
  printf("SparseGridRegression INFO: not implemented yet.\n");
  std = 0.0;
  return 0.0;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePointFuzzy(int npts, double *X, 
                                                double *Y, double *Ystd)
{
   printf("SparseGridRegression INFO: not implemented yet.\n");
   for (int kk = 0; kk < npts; kk++) Y[kk] = 0.0;
   return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int SparseGridRegression::analyze(double *X, double *Y)
{
  int    ii, jj, kk, ll;
  double ddata, ddata2, wt;
  psVector vecCoefs;

  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::analyze ERROR - invalid arguments.\n");
    exit(1);
  } 
  if (VecLBs_.length() == 0)
  {
    printf("SparseGridRegression::analyze ERROR - input bounds not set.\n");
    exit(1);
  } 
   
  //**/ ---------------------------------------------------------------
  //**/ check sample
  //**/ ---------------------------------------------------------------
  for (ii = 0; ii < nSamples_; ii++)
  {
    for (kk = 0; kk < nInputs_; kk++)
    {
      ddata = sampleInputs_[ii*nInputs_+kk];
      ddata = ddata * (VecUBs_[kk] - VecLBs_[kk]) + VecLBs_[kk];
      if (PABS(ddata - X[ii*nInputs_+kk]) > 1.0e-13)
      {
        printf("SparseGridRegression::analyze ERROR - sample mismatch\n");
        exit(1);
      }
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ computing the coefficients 
  //**/ ---------------------------------------------------------------
  VecRegCoefs_.setLength(numPerms_);
  psMatrix MatLTables;
  MatLTables.setFormat(PS_MAT2D);
  MatLTables.setDim(nInputs_, pOrder_+1);
  double **LTables = MatLTables.getMatrix2D();
  vecCoefs.setLength(numPerms_);
  for (ii = 0; ii < numPerms_; ii++) vecCoefs[ii] = 0.0;

  for (ii = 0; ii < nSamples_; ii++)
  {
    for (jj = 0; jj < nInputs_; jj++)
    {
      ddata = 2.0 * sampleInputs_[ii*nInputs_+jj] - 1.0;
      EvalLegendrePolynomials(ddata, LTables[jj]);
    }
    for (kk = 0; kk < numPerms_; kk++)
    {
      ddata = 1.0;
      for (ll = 0; ll < nInputs_; ll++)
        ddata *= LTables[ll][MatPcePerms_.getEntry(kk,ll)];
      ddata2 = ddata * ddata;
      wt = sampleWeights_[ii];
      VecRegCoefs_[kk] += (ddata * Y[ii] * wt);
      vecCoefs[kk] += (ddata2 * wt);
    }
  } 
  //**/ finally divide numerators by denominators
  if (outputLevel_ > 0)
  {
    printf("Legendre polynomial functional forms: \n");
    printf("X normalized to Z in [-1 1] (e.g. X in [0,1] -> Z=2X-1)\n");
    printf("P_0(Z) = 1\n");
    printf("P_1(Z) = Z\n");
    printf("P_{n+1} = 1/(n+1) {(2n + 1) Z P_n(Z) - n P_{n-1}(Z)}\n");
    printEquals(PL_INFO, 0);
  }
  for (kk = 0; kk < numPerms_; kk++)
  {
    if (vecCoefs[kk] == 0.0) 
         printf("ERROR in SparseGridRegression: divide by 0.\n");
    else VecRegCoefs_[kk] /= vecCoefs[kk];
    if (outputLevel_ > 0)
    {
      printf("Legendre polynomial (");
      for (jj = 0; jj < nInputs_; jj++)
        printf("%d ", MatPcePerms_.getEntry(kk,jj));
      ddata = VecRegCoefs_[kk];
      if (PABS(ddata) < 1.0e-13) ddata = 0.0;
      printf(") coefficient = %e\n", ddata);
    }
  }
  if (outputLevel_ > 0) printAsterisks(PL_INFO, 0);
  return 0;
}

// *************************************************************************
// generate all combinations of a multivariate Legendre expansion
// This code is a direct translation from Burkardt's matlab code)
// -------------------------------------------------------------------------
int SparseGridRegression::GenPermutations()
{
  int  ii, kk, orderTmp, rvTmp;

  //**/ compute number of permutations
  numPerms_ = computeNumPCEPermutations(nInputs_, pOrder_);

  //**/ construct the permutations
  MatPcePerms_.setFormat(PS_MAT2D);
  MatPcePerms_.setDim(numPerms_, nInputs_);
  int **matPerms = MatPcePerms_.getIMatrix2D();

  numPerms_ = 0;
  for (kk = 0; kk <= pOrder_; kk++)
  {
    orderTmp = kk;
    rvTmp = 0;
    matPerms[numPerms_][0] = orderTmp;
    for (ii = 1; ii < nInputs_; ii++) matPerms[numPerms_][ii] = 0;
    while (matPerms[numPerms_][nInputs_-1] != kk)
    {
      numPerms_++;
      for (ii = 0; ii < nInputs_; ii++)
        matPerms[numPerms_][ii] = matPerms[numPerms_-1][ii];
      if (orderTmp > 1) rvTmp = 1;
      else              rvTmp++;
      matPerms[numPerms_][rvTmp-1] = 0;
      orderTmp = matPerms[numPerms_-1][rvTmp-1];
      matPerms[numPerms_][0] = orderTmp - 1;
      matPerms[numPerms_][rvTmp] = matPerms[numPerms_-1][rvTmp] + 1;
    }
    numPerms_++;
  }
  return 0;
}

// *************************************************************************
// Purpose: evaluate 1D Legendre polynomials (normalized)
// -------------------------------------------------------------------------
int SparseGridRegression::EvalLegendrePolynomials(double X, double *LTable)
{
  int    ii;
  LTable[0] = 1.0;
  if (pOrder_ >= 1)
  {
    LTable[1] = X;
    for (ii = 2; ii <= pOrder_; ii++)
      LTable[ii] = ((2 * ii - 1) * X * LTable[ii-1] -
                    (ii - 1) * LTable[ii-2]) / ii;
  }
  //**/ normalize
  //**/for (ii = 0; ii <= pOrder_; ii++) LTable[ii] *= sqrt(0.5+ii);
  return 0;
}


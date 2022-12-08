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
// Functions for the class Splines
// AUTHOR : CHARLES TONG
// DATE   : 2013
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Splines.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "PsuadeConfig.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

extern "C" 
{
  void spline1d(int,double*,double*,int,double*,double*);
  void spline2d(int,int,double*,double*,double*,int,
                double*,double*,double*);
  void spline3d(int,int,int,double*,double*,double*,double*,int,
                double*,double*,double*,double*);
}

// ************************************************************************
// Constructor for object class Mars
// ------------------------------------------------------------------------
Splines::Splines(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  double ddata;

  //**/ ----------------------------------------------------------------
  //**/ set identifier
  //**/ ----------------------------------------------------------------
  faID_ = PSUADE_RS_SPLINES;

  //**/ ----------------------------------------------------------------
  //**/ only valid for 1, 2, or 3 inputs
  //**/ ----------------------------------------------------------------
  if (nSamples_ <= 0)
  {
    printf("Splines ERROR: nSamples <= 0.\n");
    exit(1);
  }
  if (nInputs <= 0 || nInputs > 3)
  {
    printf("Splines ERROR: this method currently only supports 1, 2,\n");
    printf("               or 3 inputs on a regular grid.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Splines::~Splines()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Splines::initialize(double *XIn, double *YIn)
{
  genNDGridData(XIn, YIn, NULL, NULL, NULL);
  if (psConfig_.RSCodeGenIsOn())
    printf("Splines INFO: RS stand-alone code not available.\n");
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::genNDGridData(double *XIn, double *YIn, int *, double **, 
                           double **)
{
  //**/ ----------------------------------------------------------------
  //**/ checking sample
  //**/ ----------------------------------------------------------------
  if (nInputs_ == 1)
  {
    int ii;
    for (ii = 1; ii < nSamples_; ii++)
    {
      if (XIn[ii] <= XIn[ii-1])
      {
        printf("Splines initialize ERROR: sample not factorial (1).\n");
        printf("        Suggestion: re-order the inputs first.\n");
        exit(1);
      }
    }
    n1_ = nSamples_;
    n2_ = 1;
    n3_ = 1;
  }
  else if (nInputs_ == 2)
  {
    int ii, jj;
    //**/ find n1_
    for (ii = 1; ii < nSamples_; ii++)
      if (XIn[2*ii+1] != XIn[2*(ii-1)+1]) break;
    n1_ = ii;
    //**/ check that nSamples is divisible by n1_
    n2_ = nSamples_ / n1_;
    if (n1_ * n2_ != nSamples_)
    {
      printf("Splines initialize ERROR: sample not factorial (2a).\n");
      printf("   Note: factorial should be laid out in order ");
      printf("beginning with input 1.\n");
      exit(1);
    }
    n3_ = 1;
    //**/ check that the first n1_ elements of input 1 are increasing
    for (ii = 1; ii < n1_; ii++)
    {
      if (XIn[2*ii] <= XIn[2*(ii-1)])
      {
        printf("Splines initialize ERROR: sample not factorial (2b).\n");
        exit(1);
      }
    }
    //**/ check that the first n1_ element of input 1 are repeating
    for (ii = 1; ii < n2_; ii++)
    {
      for (jj = 0; jj < n1_; jj++)
      {
        if (XIn[2*(n1_*ii+jj)] != XIn[2*jj])
        {
          printf("Splines initialize ERROR: sample not factorial (2c).\n");
          exit(1);
        }
      }
    }
    //**/ check that every n1_ element of input 2 are increasing
    for (ii = 1; ii < n2_; ii++)
    {
      if (XIn[2*ii*n1_+1] <= XIn[2*(ii-1)*n1_+1])
      {
        printf("Splines initialize ERROR: sample not factorial (2d).\n");
        exit(1);
      }
    }
    //**/ check that every n2_ group of input 2 are the same
    for (ii = 0; ii < n2_; ii++)
    {
      for (jj = 1; jj < n1_; jj++)
      {
        if (XIn[2*(ii*n1_+jj)+1] != XIn[2*(ii*n1_+jj-1)+1])
        {
          printf("Splines initialize ERROR: sample not factorial (2e).\n");
          exit(1);
        }
      }
    }
    printf("Splines: number of distinct values for input 1 = %d\n", n1_);
    printf("Splines: number of distinct values for input 2 = %d\n", n2_);
  } 
  else if (nInputs_ == 3)
  {
    int ii, jj, kk;
    psVector vecInp1, vecInp2, vecInp3;
    vecInp1.setLength(nSamples_);
    vecInp2.setLength(nSamples_);
    vecInp3.setLength(nSamples_);
    for (ii = 0; ii < nSamples_; ii++) vecInp1[ii] = XIn[3*ii];
    for (ii = 0; ii < nSamples_; ii++) vecInp2[ii] = XIn[3*ii+1];
    for (ii = 0; ii < nSamples_; ii++) vecInp3[ii] = XIn[3*ii+2];
    //**/ find n1_
    for (ii = 1; ii < nSamples_; ii++)
      if (vecInp2[ii] != vecInp2[ii-1]) break;
    n1_ = ii;
    printf("Splines: number of distinct values for input 1 = %d\n", n1_);
    //**/ find n2_
    for (ii = 1; ii < nSamples_; ii++)
      if (vecInp3[ii] != vecInp3[ii-1]) break;
    n2_ = ii / n1_;
    n3_ = nSamples_ / n1_ / n2_;
    printf("Splines: number of distinct values for input 2 = %d\n", n2_);
    printf("Splines: number of distinct values for input 3 = %d\n", n3_);
    if (n1_ * n2_ * n3_ != nSamples_)
    {
      printf("Splines initialize ERROR: sample not factorial (3a).\n");
      printf("   Note: factorial should be laid out in order ");
      printf("beginning with input 1.\n");
      exit(1);
    }
    //**/ check that the first n1_ elements of input 1 are increasing
    for (ii = 1; ii < n1_; ii++)
    {
      if (vecInp1[ii] <= vecInp1[ii-1])
      {
        printf("Splines initialize ERROR: sample not factorial (3b).\n");
        exit(1);
      }
    }
    //**/ check that the first n1_ elements of input 1 are repeating
    for (kk = 1; kk < n2_*n3_; kk++)
    {
      for (ii = 0; ii < n1_; ii++)
      {
        if (vecInp1[kk*n1_+ii] != vecInp1[ii])
        {
          printf("Splines initialize ERROR: sample not factorial (3c).\n");
          exit(1);
        }
      }
    }
    //**/ check that the first n2_ elements of input 2 are increasing
    for (ii = 1; ii < n2_; ii++)
    {
      if (vecInp2[ii*n1_] <= vecInp2[(ii-1)*n1_])
      {
        printf("Splines initialize ERROR: sample not factorial (3d).\n");
        exit(1);
      }
    }
    //**/ check that each group of n1_ elements of input 2 are the same
    for (kk = 0; kk < n2_*n3_; kk++)
    {
      for (ii = 1; ii < n1_; ii++)
      {
        if (vecInp2[kk*n1_+ii] != vecInp2[kk*n1_])
        {
          printf("Splines initialize ERROR: sample not factorial (3e).\n");
          exit(1);
        }
      }
    }
    //**/ check that the first n2_ elements of input 2 are repeating
    for (kk = 1; kk < n3_; kk++)
    {
      for (ii = 0; ii < n1_*n2_; ii++)
      {
        if (vecInp2[kk*n1_*n2_+ii] != vecInp2[ii])
        {
          printf("Splines initialize ERROR: sample not factorial (3f).\n");
          exit(1);
        }
      }
    }
    //**/ check that the first n3_ elements of input 3 are increasing
    for (ii = 1; ii < n3_; ii++)
    {
      if (vecInp3[ii*n1_*n2_] <= vecInp3[(ii-1)*n1_*n2_])
      {
        printf("Splines initialize ERROR: sample not factorial (3g).\n");
        exit(1);
      }
    }
    //**/ check that the n3_ elements of input 3 are repeating
    for (kk = 0; kk < n3_; kk++)
    {
      for (ii = 1; ii < n1_*n2_; ii++)
      {
        if (vecInp3[kk*n1_*n2_+ii] != vecInp3[kk*n1_*n2_])
        {
          printf("Splines initialize ERROR: sample not factorial (3h).\n");
          exit(1);
        }
      }
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ no preprocessing needed
  //**/ ----------------------------------------------------------------
  vecSamIns_.setLength(n1_+n2_+n3_);
  int ii;
  for (ii = 0; ii < n1_; ii++) vecSamIns_[ii] = XIn[ii*nInputs_];
  if (nInputs_ >= 2)
  {
    for (ii = 0; ii < n2_; ii++)
      vecSamIns_[n1_+ii] = XIn[nInputs_*ii*n1_+1];
  }
  if (nInputs_ == 3)
  {
    for (ii = 0; ii < n3_; ii++)
      vecSamIns_[n1_+n2_+ii] = XIn[nInputs_*ii*n1_*n2_+2];
  }
  vecSamOut_.setLength(nSamples_);
  for (ii = 0; ii < nSamples_; ii++) vecSamOut_[ii] = YIn[ii];
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::gen1DGridData(double *XIn,double *YIn,int ind1,double *setting,
                           int *NOut,double **XOut,double **YOut)
{
  int    ii;
  double HX;
  psVector vecX1, vecX2, vecX3, vecYOut;
  //**/ ----------------------------------------------------------------
  //**/ error checking
  //**/ ----------------------------------------------------------------
  if (ind1 < 0 || ind1 > 2)
  {
    printf("Splines gen1DGrid ERROR: invalid input.\n");
    printf("                         Input 1 = %d\n", ind1+1);
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  genNDGridData(XIn, YIn, NOut, XOut, YOut);
  if ((*NOut) == -999) return 0;

  //**/ ----------------------------------------------------------------
  //**/ generate regular grid data in X1 and X2
  //**/ ----------------------------------------------------------------
  HX = (VecUBs_[ind1] - VecLBs_[ind1]) / (nPtsPerDim_ - 1); 
  vecX1.setLength(nPtsPerDim_);
  vecX2.setLength(nPtsPerDim_);
  vecX3.setLength(nPtsPerDim_);
  vecYOut.setLength(nPtsPerDim_);
  if (ind1 == 0)
  {
    for (ii = 0; ii < nPtsPerDim_; ii++)
    {
      vecX1[ii] = HX * ii + VecLBs_[ind1];
      if (nInputs_ > 1) vecX2[ii] = setting[1];
      if (nInputs_ > 2) vecX3[ii] = setting[2];
    }
  }
  else if (ind1 == 1)
  {
    for (ii = 0; ii < nPtsPerDim_; ii++)
    {
      vecX2[ii] = HX * ii + VecLBs_[ind1];
      vecX1[ii] = setting[0];
      if (nInputs_ > 2) vecX3[ii] = setting[2];
    }
  }
  else
  {
    for (ii = 0; ii < nPtsPerDim_; ii++)
    {
      vecX3[ii] = HX * ii + VecLBs_[ind1];
      vecX1[ii] = setting[0];
      vecX2[ii] = setting[1];
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ interpolate
  //**/ ----------------------------------------------------------------
  double *samIns = vecSamIns_.getDVector();
  if (nInputs_ == 1) 
    spline1d(nSamples_, XIn, YIn, nPtsPerDim_, vecX1.getDVector(), 
             vecYOut.getDVector()); 
  else if (nInputs_ == 2)
  {
    spline2d(n1_, n2_, samIns, &samIns[n1_],YIn,nPtsPerDim_,
             vecX1.getDVector(),vecX2.getDVector(),vecYOut.getDVector()); 
  }
  else
  {
    spline3d(n1_,n2_,n3_,samIns,&samIns[n1_],&samIns[n1_+n2_],YIn,
             nPtsPerDim_,vecX1.getDVector(),vecX2.getDVector(),
             vecX3.getDVector(),vecYOut.getDVector()); 
  }

  //**/ ----------------------------------------------------------------
  //**/ clean up
  //**/ ----------------------------------------------------------------
  if (ind1 == 0)
  {
    (*XOut) = vecX1.takeDVector();
  }
  else if (ind1 == 1)
  {
    (*XOut) = vecX2.takeDVector();
  }
  else if (ind1 == 2)
  {
    (*XOut) = vecX3.takeDVector();
  }
  (*YOut) = vecYOut.takeDVector();
  (*NOut) = nPtsPerDim_;
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           double *settings,int *NOut,double **XOut, 
                           double **YOut)
{
  int totPts, ii, jj, index;
  psVector vecHX, vecX1, vecX2, vecX3, vecXOut, vecYOut;

  //**/ ----------------------------------------------------------------
  //**/ error checking
  //**/ ----------------------------------------------------------------
  if (nInputs_ < 2)
  {
    printf("Splines gen2DGrid ERROR: nInputs < 2.\n");
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  if (ind1 < 0 || ind1 > 2 || ind2 < 0 || ind2 > 2 || ind1 >= nInputs_ ||
      ind2 >= nInputs_)
  {
    printf("Splines gen2DGrid ERROR: invalid input.\n");
    printf("                         Input 1 = %d\n", ind1+1);
    printf("                         Input 2 = %d\n", ind2+1);
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  if (ind1 == ind2)
  {
    printf("Splines gen2DGrid ERROR: the two inputs should be distinct.\n");
    printf("                         Input 1 = %d\n", ind1+1);
    printf("                         Input 2 = %d\n", ind2+1);
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  genNDGridData(XIn, YIn, NOut, XOut, YOut);
  if ((*NOut) == -999) return 0;

  //**/ ----------------------------------------------------------------
  //**/ set up regular grid
  //**/ ----------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1);
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1);

  vecX1.setLength(totPts);
  vecX2.setLength(totPts);
  vecX3.setLength(totPts);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      index = ii * nPtsPerDim_ + jj;
      vecX1[index] = settings[0];
      vecX2[index] = settings[1];
      vecX3[index] = settings[2];
      if (ind1 == 0)
        vecX1[index] = vecHX[0] * ii + VecLBs_[ind1];
      else if (ind1 == 1)
        vecX2[index] = vecHX[0] * ii + VecLBs_[ind1];
      else if (ind1 == 2)
        vecX3[index] = vecHX[0] * ii + VecLBs_[ind1];

      if (ind2 == 0)
        vecX1[index] = vecHX[1] * jj + VecLBs_[ind2];
      else if (ind2 == 1)
        vecX2[index] = vecHX[1] * jj + VecLBs_[ind2];
      else if (ind2 == 2)
        vecX3[index] = vecHX[1] * jj + VecLBs_[ind2];
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ generate spline interpolations
  //**/ ----------------------------------------------------------------
  double *samIns = vecSamIns_.getDVector();
  vecYOut.setLength(totPts);
  if (nInputs_ == 2)
  {
    spline2d(n1_,n2_,samIns, &samIns[n1_],YIn,totPts,vecX1.getDVector(),
             vecX2.getDVector(), vecYOut.getDVector()); 
  }
  else
  {
    spline3d(n1_,n2_,n3_,samIns,&samIns[n1_],&samIns[n1_+n2_],YIn,
             totPts,vecX1.getDVector(),vecX2.getDVector(),
             vecX3.getDVector(),vecYOut.getDVector()); 
  }

  //**/ ----------------------------------------------------------------
  //**/ return values
  //**/ ----------------------------------------------------------------
  vecXOut.setLength(2*totPts);
  (*XOut) = vecXOut.takeDVector();
  (*YOut) = vecYOut.takeDVector();
  (*NOut) = totPts;
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      index = ii * nPtsPerDim_ + jj;
      (*XOut)[2*index]   = vecHX[0] * ii + VecLBs_[ind1];
      (*XOut)[2*index+1] = vecHX[1] * jj + VecLBs_[ind2];
    }
  }
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Splines::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           int ind3, double *settings, int *NOut, 
                           double **XOut, double **YOut)
{
  int totPts, ii, jj, ll, index;
  psVector vecHX, vecX1, vecX2, vecX3, vecXOut, vecYOut;

  //**/ ----------------------------------------------------------------
  //**/ error checking
  //**/ ----------------------------------------------------------------
  if (nInputs_ < 3)
  {
    printf("Splines gen3DGrid ERROR: nInputs < 3.\n");
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  if (ind1 < 0 || ind1 > 2 || ind2 < 0 || ind2 > 2 || ind3 < 0 || ind3 > 2)
  {
    printf("Splines gen3DGrid ERROR: invalid input.\n");
    printf("                         Input 1 = %d\n", ind1+1);
    printf("                         Input 2 = %d\n", ind2+1);
    printf("                         Input 3 = %d\n", ind3+1);
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  if (ind1 == ind2 || ind2 == ind3 || ind1 == ind3)
  {
    printf("Splines gen3DGrid ERROR: the three inputs should be distinct.\n");
    printf("                         Input 1 = %d\n", ind1+1);
    printf("                         Input 2 = %d\n", ind2+1);
    printf("                         Input 3 = %d\n", ind3+1);
    (*NOut) = 0;
    (*XOut) = NULL;
    (*YOut) = NULL;
    return -1;
  }
  genNDGridData(XIn, YIn, NOut, XOut, YOut);
  if ((*NOut) == -999) return 0;

  //**/ ----------------------------------------------------------------
  //**/ set up regular grid
  //**/ ----------------------------------------------------------------
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (VecUBs_[ind1] - VecLBs_[ind1])/(nPtsPerDim_ - 1);
  vecHX[1] = (VecUBs_[ind2] - VecLBs_[ind2])/(nPtsPerDim_ - 1);
  vecHX[2] = (VecUBs_[ind3] - VecLBs_[ind3])/(nPtsPerDim_ - 1);

  vecX1.setLength(totPts);
  vecX2.setLength(totPts);
  vecX3.setLength(totPts);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      for (ll = 0; ll < nPtsPerDim_; ll++)
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        vecX1[index] = settings[0];
        vecX2[index] = settings[1];
        vecX3[index] = settings[2];
        if (ind1 == 0)
          vecX1[index] = vecHX[0] * ii + VecLBs_[ind1];
        else if (ind1 == 1)
          vecX2[index] = vecHX[0] * ii + VecLBs_[ind1];
        else if (ind1 == 2)
          vecX3[index] = vecHX[0] * ii + VecLBs_[ind1];

        if (ind2 == 0)
          vecX1[index] = vecHX[1] * jj + VecLBs_[ind2];
        else if (ind2 == 1)
          vecX2[index] = vecHX[1] * jj + VecLBs_[ind2];
        else if (ind2 == 2)
          vecX3[index] = vecHX[1] * jj + VecLBs_[ind2];

        if (ind3 == 0)
          vecX1[index] = vecHX[2] * ll + VecLBs_[ind3];
        else if (ind3 == 1)
          vecX2[index] = vecHX[2] * ll + VecLBs_[ind3];
        else if (ind3 == 2)
          vecX3[index] = vecHX[2] * ll + VecLBs_[ind3];
      }
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ generate spline interpolations
  //**/ ----------------------------------------------------------------
  vecYOut.setLength(totPts);
  double *samIns = vecSamIns_.getDVector();
  spline3d(n1_,n2_,n3_,samIns,&samIns[n1_],&samIns[n1_+n2_],YIn,totPts,
           vecX1.getDVector(),vecX2.getDVector(),vecX3.getDVector(),
           vecYOut.getDVector()); 

  //**/ ----------------------------------------------------------------
  //**/ return the outputs
  //**/ ----------------------------------------------------------------
  vecXOut.setLength(3*totPts);
  for (ii = 0; ii < nPtsPerDim_; ii++)
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      for (ll = 0; ll < nPtsPerDim_; ll++)
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        vecXOut[index*3]   = vecHX[0] * ii + VecLBs_[ind1];
        vecXOut[index*3+1] = vecHX[1] * jj + VecLBs_[ind2];
        vecXOut[index*3+2] = vecHX[2] * ll + VecLBs_[ind3];
      }
    }
  }
  (*XOut) = vecXOut.takeDVector();
  (*YOut) = vecYOut.takeDVector();
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Splines::gen4DGridData(double *, double *, int, int, int, int, double *, 
                           int *NOut, double **XOut, double **YOut)
{
  printf("Splines ERROR: not relevant since nInputs <= 2.\n");
  (*NOut) = 0;
  (*XOut) = NULL;
  (*YOut) = NULL;
  return -1;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Splines::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y=0.0, *samIns;
 
  samIns = vecSamIns_.getDVector();
  if (nInputs_ == 1)
    spline1d(nSamples_,samIns,vecSamOut_.getDVector(),iOne, X, &Y); 
  else if (nInputs_ == 2)
  {
    spline2d(n1_,n2_,samIns,&samIns[n1_],vecSamOut_.getDVector(),
             iOne,X,&X[1],&Y);
  }
  else 
  {
    spline3d(n1_,n2_,n3_,samIns,&samIns[n1_],&samIns[n1_+n2_],
             vecSamOut_.getDVector(),iOne,X,&X[1],&X[2],&Y);
  }
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Splines::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[nInputs_*kk]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double Splines::evaluatePointFuzzy(double *X, double &std)
{
  double Y = evaluatePoint(X);
  std = 0.0;
  return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double Splines::evaluatePointFuzzy(int npts,double *X,double *Y,double *Ysd)
{
  evaluatePoint(npts, X, Y);
  for (int kk = 0; kk < npts; kk++) Ysd[kk] = 0.0;
  return 0.0;
}


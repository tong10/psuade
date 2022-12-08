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
// Functions for the class PDFBase
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "DataIO/PsuadeData.h"
#include "PDFBase.h"
#include "PDFNormal.h"
#include "PDFLogNormal.h"
#include "PDFTriangle.h"
#include "PDFBeta.h"
#include "psVector.h"

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
PDFBase::PDFBase()
{
}
   
// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFBase::~PDFBase()
{
}

// ************************************************************************
// forward function 
// ------------------------------------------------------------------------
int PDFBase::getPDF(int, double *, double *)
{
  return 0;
}

// ************************************************************************
// look up cumulative distribution function 
// ------------------------------------------------------------------------
int PDFBase::getCDF(int, double *, double *)
{
  printf("genCDF not available\n");
  return 0;
}

// ************************************************************************
// look up cumulative distribution function 
// ------------------------------------------------------------------------
int PDFBase::invCDF(int, double *, double *)
{
  printf("invCDF not available\n");
  return 0;
}

// ************************************************************************
// generate sample
// ------------------------------------------------------------------------
int PDFBase::genSample(int, double *, double *, double *) 
{
  return 0;
}

// ************************************************************************
// get the mean 
// ------------------------------------------------------------------------
double PDFBase::getMean()
{
  return 0;
}

// ************************************************************************
// get means (for group) 
// ------------------------------------------------------------------------
int PDFBase::getMeans(double *)
{
  return 0;
}

// ************************************************************************
// set internal parameter
// ------------------------------------------------------------------------
int PDFBase::setParam(char *)
{
   return 0;
}

// ************************************************************************
// transform from uniform to normal distribution
// ------------------------------------------------------------------------
                                                                                
extern "C" 
{
int PDFTransform(PsuadeData *psuadeIO, int nSamples, int nInputs, 
                 double *sampleData)
{
  int      ii, ss;
  pData    pPtr, pLowerB, pUpperB, pPDFs, pMeans, pStdevs;
  psVector vecData1, vecData2;
                                                                               
  //**/ ----------------------------------------------------------------
  //**/ fetch input parameters
  //**/ ----------------------------------------------------------------
  psuadeIO->getParameter("input_lbounds", pLowerB);
  double *LBounds = pLowerB.dbleArray_;
  psuadeIO->getParameter("input_ubounds", pUpperB);
  double *UBounds = pUpperB.dbleArray_;
  psuadeIO->getParameter("input_pdfs", pPDFs);
  int *inputPDFs = pPDFs.intArray_;
  psuadeIO->getParameter("input_means", pMeans);
  double *inputMeans = pMeans.dbleArray_;
  psuadeIO->getParameter("input_stdevs", pStdevs);
  double *inputStdevs = pStdevs.dbleArray_;
                                                                               
  vecData1.setLength(nSamples);
  vecData2.setLength(nSamples);
  PDFBase **PDFPtrs = new PDFBase*[nInputs];
                                                                               
  //**/ ----------------------------------------------------------------
  //**/ construct the corresponding PDF objects
  //**/ ----------------------------------------------------------------
  for (ii = 0; ii < nInputs; ii++)
  {
    if (inputPDFs[ii] == 1)
         PDFPtrs[ii] = (PDFBase *) new PDFNormal(inputMeans[ii],
                                            inputStdevs[ii]);
    else if (inputPDFs[ii] == 2)
         PDFPtrs[ii] = (PDFBase *) new PDFLogNormal(inputMeans[ii],
                                            inputStdevs[ii]);
    else if (inputPDFs[ii] == 3)
         PDFPtrs[ii] = (PDFBase *) new PDFTriangle(inputMeans[ii],
                                            inputStdevs[ii]);
    else if (inputPDFs[ii] == 4)
         PDFPtrs[ii] = (PDFBase *) new PDFBeta(inputMeans[ii],
                                               inputStdevs[ii]);
    else PDFPtrs[ii] = NULL;
  }
                                                                               
  //**/ ----------------------------------------------------------------
  //**/ perform transformation
  //**/ ----------------------------------------------------------------
  for (ii = 0; ii < nInputs; ii++)
  {
    if (PDFPtrs[ii] != NULL)
    {
      for (ss = 0; ss < nSamples; ss++)
        vecData1[ss] = sampleData[ss*nInputs+ii];
      PDFPtrs[ii]->invCDF(nSamples,vecData1.getDVector(),
                          vecData1.getDVector());
      for (ss = 0; ss < nSamples; ss++)
      {
        sampleData[ss*nInputs+ii] = vecData2[ss];
        if (sampleData[ss*nInputs+ii] < LBounds[ii])
          sampleData[ss*nInputs+ii] = LBounds[ii];
        if (sampleData[ss*nInputs+ii] > UBounds[ii])
          sampleData[ss*nInputs+ii] = UBounds[ii];
      }
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ clean up
  //**/ ----------------------------------------------------------------
  for (ii = 0; ii < nInputs; ii++)
     if (PDFPtrs[ii] != NULL) delete PDFPtrs[ii];
  delete [] PDFPtrs;
  return 0;
}
}

// ************************************************************************
// transform from uniform to normal distribution for one input
// ------------------------------------------------------------------------
                                                                                
extern "C" 
{
int PDFTransformSingle(int ptype, int nSamples, double *sampleData, 
                       double lbound, double ubound, double mean, 
                       double stdev)
{
  int ss;
  PDFBase *PDFPtrs=NULL;
  psVector vecData;

  //**/ ----------------------------------------------------------------
  //**/ memory allocation
  //**/ ----------------------------------------------------------------
  vecData.setLength(nSamples);
                                                                               
  //**/ ----------------------------------------------------------------
  //**/ construct the corresponding PDF objects
  //**/ ----------------------------------------------------------------
  if (ptype == 1)
       PDFPtrs = (PDFBase *) new PDFNormal(mean, stdev);
  else if (ptype == 2)
       PDFPtrs = (PDFBase *) new PDFLogNormal(mean, stdev);
  else if (ptype == 3)
       PDFPtrs = (PDFBase *) new PDFTriangle(mean, stdev);
  else if (ptype == 4)
       PDFPtrs = (PDFBase *) new PDFBeta(mean, stdev);
  else PDFPtrs = NULL;
                                                                               
  //**/ ----------------------------------------------------------------
  //**/ perform transformation
  //**/ ----------------------------------------------------------------
  if (PDFPtrs != NULL)
  {
    PDFPtrs->invCDF(nSamples,sampleData,vecData.getDVector());
    for (ss = 0; ss < nSamples; ss++) 
    {
      sampleData[ss] = vecData[ss];
      if (sampleData[ss] < lbound) sampleData[ss] = lbound;
      if (sampleData[ss] > ubound) sampleData[ss] = ubound;
    }
  }

  //**/ ----------------------------------------------------------------
  //**/ clean up
  //**/ ----------------------------------------------------------------
  if (PDFPtrs != NULL) delete PDFPtrs;
  return 0;
}
}


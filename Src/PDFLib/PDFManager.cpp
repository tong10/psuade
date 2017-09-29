// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC. 
// Produced at the Lawrence Livermore National Laboratory.
// Written by the Charles Tong <tong10@llnl.gov>.
//
// CODE UCRL-CODE-235523 All rights reserved.
//
// This file is part of PSUADE.
//
// Please see the GPL.pdf file for the copyright notice,
//
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License (as published by the 
// Free Software Foundation) version 2.1 dated February 1991.
//
// This program is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
// or FITNESS FOR A PARTICULAR PURPOSE.  See the terms and conditions of the 
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class PDFManager
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "PDFManager.h"
#include "DataIO/PsuadeData.h"
#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Util/Matrix.h"
#include "Util/Vector.h"
#include "PDFNormal.h"
#include "PDFLogNormal.h"
#include "PDFMVNormal.h"
#include "PDFMVLogNormal.h"
#include "PDFTriangle.h"
#include "PDFBeta.h"
#include "PDFWeibull.h"
#include "PDFGamma.h"
#include "PDFExponential.h"

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
PDFManager::PDFManager()
{
   nInputs_ = 0;
   pdfMap_ = NULL;
   PDFptrs_ = NULL;
   PDFMVNormalPtr_ = NULL;
   PDFMVLogNormalPtr_ = NULL;
   PSIOptr_ = NULL;
   nGNormal_ = 0;
   gnormalInputs_ = NULL;
   nGLognormal_ = 0;
   glognormalInputs_ = NULL;
   printLevel_ = 0;
}
   
// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFManager::~PDFManager()
{
   if (pdfMap_ != NULL) delete [] pdfMap_;
   if (PDFptrs_ != NULL)
   {
      for (int ii = 0; ii < nInputs_; ii++)
         if (PDFptrs_[ii] != NULL) delete PDFptrs_[ii];
      delete [] PDFptrs_;
   }
   if (PDFMVNormalPtr_ != NULL) delete PDFMVNormalPtr_;
   if (PDFMVLogNormalPtr_ != NULL) delete PDFMVLogNormalPtr_;
   if (gnormalInputs_ != NULL) delete [] gnormalInputs_;
   if (glognormalInputs_ != NULL) delete [] glognormalInputs_;
}

// ************************************************************************
// initialize function 
// ------------------------------------------------------------------------
int PDFManager::initialize(PsuadeData *psData)
{
   int     ii, nInputs, *inputPDFs, count;
   double  *inputMeans, *inputStdevs;
   pData   pPtr, pPDFs, pMeans, pStdevs, pCorMat;
   Matrix  *corMatrix;

   if (pdfMap_ != NULL) delete [] pdfMap_;
   if (PDFptrs_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (PDFptrs_[ii] != NULL) delete PDFptrs_[ii];
      delete [] PDFptrs_;
   }
   if (PDFMVNormalPtr_ != NULL) delete PDFMVNormalPtr_;
   if (PDFMVLogNormalPtr_ != NULL) delete PDFMVLogNormalPtr_;
   if (gnormalInputs_ != NULL) delete [] gnormalInputs_;
   if (glognormalInputs_ != NULL) delete [] glognormalInputs_;
   nInputs_ = nGNormal_ = nGLognormal_ = 0;
   pdfMap_ = NULL;
   PDFptrs_ = NULL;
   PDFMVNormalPtr_ = NULL;
   PDFMVLogNormalPtr_ = NULL;
   PSIOptr_ = NULL;
   gnormalInputs_ = glognormalInputs_ = NULL;

   if (psData == NULL)
   {
      printf("PDFManager ERROR: no PsuadeData object given.\n");
      exit(1);
   } 

   PSIOptr_ = psData;
   psData->getParameter("ana_diagnostics", pPtr);
   printLevel_ = pPtr.intData_;
   psData->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psData->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psData->getParameter("input_means", pMeans);
   inputMeans = pMeans.dbleArray_;
   psData->getParameter("input_stdevs", pStdevs);
   inputStdevs = pStdevs.dbleArray_; 
   psData->getParameter("input_cor_matrix", pCorMat);
   corMatrix = (Matrix *) pCorMat.psObject_; 
   pCorMat.psObject_ = NULL; 
   count = 0;
   for (ii = 0; ii < nInputs; ii++) count += inputPDFs[ii];
   if (count > 0) 
      printf("PDFManager INFO: some distributions are not UNIFORM.\n");
   initialize(nInputs, inputPDFs, inputMeans, inputStdevs, *corMatrix);
   return 0;
}

// ************************************************************************
// another initialize function 
// ------------------------------------------------------------------------
int PDFManager::initialize(int nInputs, int *inputPDFs, 
                           double *inputMeans, double *inputStdevs,
                           Matrix &corMatrix)
{
   int     ii, jj, errFlag, normalFlag, lognormalFlag;
   double  *means, mean, stdev;
   Matrix  corMat;
   Vector  meanVec;

   nInputs_ = nInputs;
   corMat_.load(corMatrix);
   if (corMatrix.nrows() != corMatrix.ncols())
   {
      printf("PDFManager ERROR: correlation matrix is not square.\n");
      exit(1);
   }
   if (corMatrix.nrows() != nInputs)
   {
      printf("PDFManager ERROR: correlation matrix dimension mismatch.\n"); 
      exit(1);
   }

   errFlag = 0;
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (jj = 0; jj < nInputs_; jj++)
      {
         //*/ first check that only same distributions are allowed to 
         //*/ have correlation != 0
         if (inputPDFs[ii] != inputPDFs[jj] && 
             corMat_.getEntry(ii,jj) != 0.0)
         {
            printf("PDFManager ERROR: correlation != 0 for different PDFs.\n");
            errFlag = 1;
         }
         //*/ check that only normal and lognormal distributions are 
         //*/ allowed to have correlation != 0
         if (inputPDFs[ii] != PSUADE_PDF_LOGNORMAL && 
             inputPDFs[ii] != PSUADE_PDF_NORMAL && ii != jj && 
             corMat_.getEntry(ii,jj) != 0.0)
         {
            printf("PDFManager ERROR: correlation != 0 for non-normal PDFs.\n");
            errFlag = 1;
         }
         //*/ check that correlation cannot be less than -1 or > 1
         if (corMat_.getEntry(ii,jj) < -1 || corMat_.getEntry(ii,jj) > 1)
         {
            printf("PDFManager ERROR: correlation should be in [-1,1].\n");
            errFlag = 1;
         }
         //*/ check that covariance matrix is symmetric
         if (corMat_.getEntry(ii,jj) != corMat_.getEntry(jj,ii))
         {
            printf("PDFManager ERROR : correlation matrix not symmetric.\n");
            errFlag = 1;
         }
      }
   }
   assert(errFlag == 0);         

   nGNormal_ = nGLognormal_ = 0;
   gnormalInputs_ = new int[nInputs_];
   glognormalInputs_ = new int[nInputs_];
   pdfMap_ = new int[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      pdfMap_[ii] = inputPDFs[ii];
      if (inputPDFs[ii] == PSUADE_PDF_NORMAL)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            if (ii != jj && corMat_.getEntry(ii,jj) != 0.0)
            {
               pdfMap_[ii] = 1000 + PSUADE_PDF_NORMAL;
               gnormalInputs_[nGNormal_++] = ii;
            }
         }
      }
      if (inputPDFs[ii] == PSUADE_PDF_LOGNORMAL)
      {
         for (jj = 0; jj < nInputs_; jj++)
         {
            if (ii != jj && corMat_.getEntry(ii,jj) != 0.0)
            {
               pdfMap_[ii] = 1000 + PSUADE_PDF_LOGNORMAL;
               glognormalInputs_[nGLognormal_++] = ii;
            }
         }
      }
   } 

   if (corMat_.nrows() != nInputs) corMat_.setDim(nInputs, nInputs);
   means_.setLength(nInputs_);
   for (ii = 0; ii < nInputs_; ii++)
   {
      stdev = inputStdevs[ii]; 
      corMat_.setEntry(ii, ii, stdev);
      means_[ii] = inputMeans[ii];
   }

   PDFptrs_ = new PDFBase*[nInputs_];
   normalFlag = lognormalFlag = 0;
   means = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      mean = inputMeans[ii]; 
      stdev = inputStdevs[ii]; 
      PDFptrs_[ii] = NULL;
      if (pdfMap_[ii] == PSUADE_PDF_NORMAL)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFNormal(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_LOGNORMAL)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFLogNormal(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_TRIANGLE)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFTriangle(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_BETA)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFBeta(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_WEIBULL)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFWeibull(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_GAMMA)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFGamma(mean, stdev);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_EXPONENTIAL)
      {
         PDFptrs_[ii] = (PDFBase *) new PDFExponential(mean, stdev);
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         if (normalFlag == 0)
         {
            corMat.submatrix(corMat_, nGNormal_, gnormalInputs_);
            for (jj = 0; jj < nGNormal_; jj++)
               means[jj] = inputMeans[gnormalInputs_[jj]];
            meanVec.load(nGNormal_, means);
            PDFMVNormalPtr_ = new PDFMVNormal(meanVec, corMat);
            normalFlag = 1;
         }
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         if (lognormalFlag == 0)
         {
            corMat.submatrix(corMat_, nGLognormal_, glognormalInputs_);
            for (jj = 0; jj < nGLognormal_; jj++)
               means[jj] = inputMeans[glognormalInputs_[jj]];
            meanVec.load(nGLognormal_, means);
            PDFMVLogNormalPtr_ = new PDFMVLogNormal(meanVec, corMat);
            lognormalFlag = 1;
         }
      }
   }
   delete [] means;
   return 0;
}

// ************************************************************************
// forward transformation 
// ------------------------------------------------------------------------
int PDFManager::getPDF(int nSamples, Vector &vecIn, Vector &vecOut,
                       Vector &vecLower, Vector &vecUpper)
{
   int    ss, ii, jj, nInputs;
   double *localData1, *localData2, *inData, *outData;
   Vector vecLocalIn, vecLocalOut;

   if (nInputs_ <= 0)
   {
      printf("PDFManager evaluate ERROR: nInputs <= 0.\n");
      exit(1);
   }
   nInputs = vecLower.length();
   if (nInputs_ != nInputs)
   {
      printf("PDFManager evaluate ERROR: nInputs mismatch.\n");
      exit(1);
   }
   if (nInputs_ != vecUpper.length() || nInputs_ != vecLower.length())
   {
      printf("PDFManager evaluate ERROR: vecBounds length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecIn.length()/nInputs)
   {
      printf("PDFManager evaluate ERROR: vecIn length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecOut.length()/nInputs)
   {
      printf("PDFManager evaluate ERROR: vecOut length mismatch.\n");
      exit(1);
   }

   inData = vecIn.getDVector();
   outData = vecOut.getDVector();
   if (PDFptrs_ == NULL)
   {
      for (ii = 0; ii < nSamples; ii++)
         for (jj = 0; jj < nInputs_; jj++)
            outData[ii*nInputs_+jj] = 1.0 / (vecUpper[jj] - vecLower[jj]);
      return 0;
   }

   localData1 = new double[nSamples*nInputs_];
   localData2 = new double[nSamples*nInputs_];
   for (ii = 0; ii < nSamples*nInputs_; ii++) localData2[ii] = 0.0;
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (pdfMap_[ii] == PSUADE_PDF_UNIFORM)
      {
         if (printLevel_ > 2) printf("Input %d: uniform distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = 1.0 / (vecUpper[ii] - vecLower[ii]);
      }
      else if (pdfMap_[ii] == PSUADE_PDF_NORMAL)
      {
         if (printLevel_ > 2) printf("Input %d: normal distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_LOGNORMAL)
      {
         if (printLevel_ > 2) printf("Input %d: lognormal distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_TRIANGLE)
      {
         if (printLevel_ > 2) printf("Input %d: triangle distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_BETA)
      {
         if (printLevel_ > 2) printf("Input %d: beta distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_WEIBULL)
      {
         if (printLevel_ > 2) printf("Input %d: Weibull distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_GAMMA)
      {
         if (printLevel_ > 2) printf("Input %d: gamma distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == PSUADE_PDF_EXPONENTIAL)
      {
         if (printLevel_ > 2) printf("Input %d: exponential distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = inData[ss*nInputs+ii];
         PDFptrs_[ii]->getPDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         printf("ERROR: getPDF for multivariate normal not available.\n");
         exit(1);
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         printf("ERROR: getPDF for multivariate lognormal not available.\n");
         exit(1);
      }
   }
   delete [] localData1;
   delete [] localData2;
   return 0;
}

// ************************************************************************
// generate sample
// ------------------------------------------------------------------------
int PDFManager::genSample(int nSamples, Vector &vecOut, 
                          Vector &vecLower, Vector &vecUpper)
{
   int    ss, ii, jj, normalFlag, lognormalFlag, ind, nInputs, *iVec;
   double *localData, *outData, range, low, high;
   Vector vecLocalOut, vLower, vUpper;

   if (nInputs_ <= 0)
   {
      printf("PDFManager genSample ERROR: nInputs <= 0.\n");
      exit(1);
   }
   nInputs = vecLower.length();
   if (nInputs_ != nInputs)
   {
      printf("PDFManager genSample ERROR: nInputs mismatch.\n");
      exit(1);
   }
   if (nInputs_ != vecUpper.length() || nInputs_ != vecLower.length())
   {
      printf("PDFManager genSample ERROR: vecBounds length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecOut.length()/nInputs_)
   {
      printf("PDFManager genSample ERROR: vecOut length mismatch.\n");
      exit(1);
   }

   outData    = vecOut.getDVector();
   normalFlag = lognormalFlag = 0;
   localData  = new double[nSamples*nInputs];
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (pdfMap_[ii] == PSUADE_PDF_UNIFORM)
      {
         range = vecUpper[ii] - vecLower[ii];
         iVec = new int[nSamples]; 
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = PSUADE_drand() * range + vecLower[ii];
         delete [] iVec;
      }
      if (pdfMap_[ii] == PSUADE_PDF_TRIANGLE ||
          pdfMap_[ii] == PSUADE_PDF_BETA ||
          pdfMap_[ii] == PSUADE_PDF_WEIBULL || 
          pdfMap_[ii] == PSUADE_PDF_NORMAL ||
          pdfMap_[ii] == PSUADE_PDF_LOGNORMAL ||
          pdfMap_[ii] == PSUADE_PDF_GAMMA ||
          pdfMap_[ii] == PSUADE_PDF_EXPONENTIAL)
      {
         low = vecLower[ii];
         high = vecUpper[ii];
         PDFptrs_[ii]->genSample(nSamples,localData,low,high);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData[ss];
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         if (normalFlag == 0)
         {
            vLower.setLength(nGNormal_);
            vUpper.setLength(nGNormal_);
            vecLocalOut.setLength(nSamples*nGNormal_);
            for (jj = 0; jj < nGNormal_; jj++)
            {
               ind = gnormalInputs_[jj];
               vLower[jj] = vecLower[ind];
               vUpper[jj] = vecUpper[ind];
            }
            PDFMVNormalPtr_->genSample(nSamples,vecLocalOut,vLower,vUpper);
            for (jj = 0; jj < nGNormal_; jj++)
            {
               ind = gnormalInputs_[jj];
               for (ss = 0; ss < nSamples; ss++)
               {
                  outData[ss*nInputs+ind] = vecLocalOut[ss*nGNormal_+jj];
                  if (outData[jj*nInputs_+ind] > vecUpper[ind])
                     outData[jj*nInputs_+ind] = vecUpper[ind];
                  if (outData[jj*nInputs_+ind] < vecLower[ind])
                     outData[jj*nInputs_+ind] = vecLower[ind];
               }
            }
            normalFlag = 1;
         }
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         if (lognormalFlag == 0)
         {
            vLower.setLength(nGLognormal_);
            vUpper.setLength(nGLognormal_);
            vecLocalOut.setLength(nSamples*nGLognormal_);
            for (jj = 0; jj < nGLognormal_; jj++)
            {
               ind = glognormalInputs_[jj];
               vLower[jj] = vecLower[ind];
               vUpper[jj] = vecUpper[ind];
            }
            PDFMVLogNormalPtr_->genSample(nSamples,vecLocalOut,vLower,vUpper);
            for (jj = 0; jj < nGLognormal_; jj++)
            {
               ind = glognormalInputs_[jj];
               for (ss = 0; ss < nSamples; ss++)
               {
                  outData[ss*nInputs+ind] = vecLocalOut[ss*nGLognormal_+jj];
                  if (outData[jj*nInputs_+ind] > vecUpper[ind])
                     outData[jj*nInputs_+ind] = vecUpper[ind];
                  if (outData[jj*nInputs_+ind] < vecLower[ind])
                     outData[jj*nInputs_+ind] = vecLower[ind];
               }
            }
            lognormalFlag = 1;
         }
      }
   }
   delete [] localData;
   return 0;
}

// ************************************************************************
// generate sample (and put it back to the PsuadeData object)
// ------------------------------------------------------------------------
int PDFManager::genSample()
{
   int    ss, ii, nSamples, nInputs, nOutputs, *sampleStates;
   double *sampleInputs, *sampleOutputs, *iLowerB, *iUpperB;
   pData  pPtr, pLowerB, pUpperB;
   Vector vecLocalIn, vecLocalOut, vecUpper, vecLower;

   if (printLevel_ > 2) printf("PDFManager genSample\n");
   if (PSIOptr_ == NULL)
   {
      printf("PDFManager ERROR: no PsuadeData object given.\n");
      exit(1);
   } 

   if (PDFptrs_ == NULL) return 0;

   PSIOptr_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   PSIOptr_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   PSIOptr_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   PSIOptr_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   PSIOptr_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   sampleInputs = new double[nSamples*nInputs];
   sampleOutputs = new double[nSamples*nOutputs];
   sampleStates = new int[nSamples];

   vecUpper.load(nInputs_, iUpperB);
   vecLower.load(nInputs_, iLowerB);
   vecLocalOut.setLength(nSamples*nInputs_);
   genSample(nSamples, vecLocalOut, vecLower, vecUpper);
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ii = 0; ii < nInputs; ii++)
         sampleInputs[ss*nInputs+ii] = vecLocalOut[ss*nInputs+ii];
      for (ii = 0; ii < nOutputs; ii++)
         sampleOutputs[ss*nOutputs+ii] = PSUADE_UNDEFINED;
      sampleStates[ss] = 0;
   }
   PSIOptr_->updateInputSection(nSamples,nInputs_,NULL,NULL,NULL,
                                sampleInputs,NULL);
   PSIOptr_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                 sampleStates,NULL);
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   return 0;
}

// ************************************************************************
// look up the cumulative density function
// ------------------------------------------------------------------------
int PDFManager::invCDF(int nSamples, Vector &vecIn, Vector &vecOut, 
                       Vector &vecLower, Vector &vecUpper)
{
   int    ss, ii, nInputs;
   double *localData1, *localData2, *inData, *outData;
   Vector vecLocalIn, vecLocalOut;

   if (nInputs_ <= 0)
   {
      printf("PDFManager invCDF ERROR: nInputs <= 0.\n");
      exit(1);
   }
   nInputs = vecLower.length();
   if (nInputs_ != nInputs)
   {
      printf("PDFManager invCDF ERROR: nInputs mismatch.\n");
      exit(1);
   }
   if (nInputs_ != vecUpper.length() || nInputs_ != vecLower.length())
   {
      printf("PDFManager invCDF ERROR: vecBounds length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecIn.length()/nInputs_)
   {
      printf("PDFManager invCDF ERROR: vecIn length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecOut.length()/nInputs_)
   {
      printf("PDFManager invCDF ERROR: vecOut length mismatch.\n");
      exit(1);
   }

   inData = vecIn.getDVector();
   outData = vecOut.getDVector();
   if (PDFptrs_ == NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = inData[ss*nInputs_+ii];
      }
      return 0;
   }

   localData1 = new double[nSamples*nInputs];
   localData2 = new double[nSamples*nInputs];
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (pdfMap_[ii] == PSUADE_PDF_UNIFORM)
      {
         if (printLevel_ > 2) printf("Input %d: uniform distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = vecIn[ss*nInputs+ii];
      }
      if (pdfMap_[ii] == PSUADE_PDF_NORMAL)
      {
         if (printLevel_ > 2) printf("Input %d: normal distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_LOGNORMAL)
      {
         if (printLevel_ > 2) printf("Input %d: lognormal distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_TRIANGLE)
      {
         if (printLevel_ > 2) printf("Input %d: triangle distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_BETA)
      {
         if (printLevel_ > 2) printf("Input %d: beta distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_WEIBULL)
      {
         if (printLevel_ > 2) printf("Input %d: Weibull distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_GAMMA)
      {
         if (printLevel_ > 2) printf("Input %d: gamma distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      if (pdfMap_[ii] == PSUADE_PDF_EXPONENTIAL)
      {
         if (printLevel_ > 2) printf("Input %d: exponential distribution.\n",ii+1);
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->invCDF(nSamples,localData1,localData2,
                              vecLower[ii], vecUpper[ii]);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         printf("invCDF for multivariate normal is not available.\n");
         exit(1);
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         printf("invCDF for multivariate lognormal is not available.\n");
         exit(1);
      }
   }
   delete [] localData1;
   delete [] localData2;
   return 0;
}

// ************************************************************************
// look up the cumulative density function
// ------------------------------------------------------------------------
int PDFManager::getCDF(int nSamples, Vector &vecIn, Vector &vecOut, 
                       Vector &vecLower, Vector &vecUpper)
{
   int    ss, ii, nInputs;
   double *localData1, *localData2, *inData, *outData;
   double range;
   Vector vecLocalIn, vecLocalOut;

   if (nInputs_ <= 0)
   {
      printf("PDFManager getCDF ERROR: nInputs <= 0.\n");
      exit(1);
   }
   nInputs = vecLower.length();
   if (nInputs_ != nInputs)
   {
      printf("PDFManager getCDF ERROR: nInputs mismatch.\n");
      exit(1);
   }
   if (nSamples != vecIn.length()/nInputs_)
   {
      printf("PDFManager getCDF ERROR: vecIn length mismatch.\n");
      exit(1);
   }
   if (nSamples != vecOut.length()/nInputs_)
   {
      printf("PDFManager getCDF ERROR: vecOut length mismatch.\n");
      exit(1);
   }

   inData = vecIn.getDVector();
   outData = vecOut.getDVector();
   if (PDFptrs_ == NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = inData[ss*nInputs_+ii];
      }
      return 0;
   }

   localData1 = new double[nSamples*nInputs];
   localData2 = new double[nSamples*nInputs];
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (pdfMap_[ii] == PSUADE_PDF_UNIFORM)
      {
         range = vecUpper[ii] - vecLower[ii];
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs+ii] = (vecIn[ss*nInputs+ii]-vecLower[ii])/range;
      }
      if (pdfMap_[ii] == PSUADE_PDF_TRIANGLE ||
          pdfMap_[ii] == PSUADE_PDF_BETA ||
          pdfMap_[ii] == PSUADE_PDF_WEIBULL ||
          pdfMap_[ii] == PSUADE_PDF_NORMAL ||
          pdfMap_[ii] == PSUADE_PDF_LOGNORMAL ||
          pdfMap_[ii] == PSUADE_PDF_GAMMA ||
          pdfMap_[ii] == PSUADE_PDF_EXPONENTIAL)
      {
         for (ss = 0; ss < nSamples; ss++)
            localData1[ss] = vecIn[ss*nInputs_+ii];
         PDFptrs_[ii]->getCDF(nSamples,localData1,localData2);
         for (ss = 0; ss < nSamples; ss++)
            outData[ss*nInputs_+ii] = localData2[ss];
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_NORMAL)
      {
         printf("getCDF for multivariate normal is not available.\n");
         exit(1);
      }
      else if (pdfMap_[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
         printf("getCDF for multivariate lognormal is not available.\n");
         exit(1);
      }
   }
   delete [] localData1;
   delete [] localData2;
   return 0;
}


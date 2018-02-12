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
// Functions for the class Sampling
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sysdef.h"
#include "sysdef.h"
#include "Sampling.h"
#include "MCSampling.h"
#include "FactorialSampling.h"
#include "LHSampling.h"
#include "OASampling.h"
#include "OALHSampling.h"
#include "MOATSampling.h"
#include "SobolSampling.h"
#include "LPtauSampling.h"
#include "MetisSampling.h"
#include "FASTSampling.h"
#include "SFASTSampling.h"
#include "BoxBehnkenSampling.h"
#include "PlackettBurmanSampling.h"
#include "FractFactSampling.h"
#include "CentralCompositeSampling.h"
#include "LSASampling.h"
#include "UserMetisSampling.h"
#include "GMOATSampling.h"
#include "GMetisSampling.h"
#include "SparseGridSampling.h"
#include "DiscreteSampling.h"
#include "RFractFactSampling.h"
#include "PrintingTS.h"
#include "PsuadeUtil.h"

// ------------------------------------------------------------------------
// local includes : The OTHER way to do sampling...
// ------------------------------------------------------------------------
#include "PDFManager.h"
#include "PDFBase.h"

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
Sampling::Sampling()
{
   printLevel_    = 0;
   samplingID_    = -1;
   nSamples_      = 0;
   nInputs_       = 0;
   nOutputs_      = 0;
   randomize_     = 0;
   nReplications_ = 1;
   lowerBounds_   = NULL;
   upperBounds_   = NULL;
   sampleMatrix_  = NULL;
   sampleOutput_  = NULL;
   sampleStates_  = NULL;
   symTable_      = NULL;
   inputSettings_ = NULL;
   inputNumSettings_ = NULL;
}

// ************************************************************************
// Copy Constructor created by Bill Oliver 
// ------------------------------------------------------------------------
Sampling::Sampling(const Sampling & smp)
{
   printLevel_ = smp.printLevel_;
   samplingID_ = smp.samplingID_;
   nSamples_ = smp.nSamples_;
   nInputs_ = smp.nInputs_;
   nOutputs_ = smp.nOutputs_;
   randomize_ = smp.randomize_;
   nReplications_ = smp.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = smp.lowerBounds_[i];
      upperBounds_[i] = smp.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = smp.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = smp.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = smp.sampleStates_[i];
   // symTable_,inputSettings_,and inputNumSettings_ are currently not 
   // being allocated
} 

// ************************************************************************
// Overload operator=  by Bill Oliver 
// ------------------------------------------------------------------------
Sampling & Sampling::operator=(const Sampling & smp)
{
   if(this == & smp) return *this;

   printLevel_ = smp.printLevel_;
   samplingID_ = smp.samplingID_;
   nSamples_ = smp.nSamples_;
   nInputs_ = smp.nInputs_;
   nOutputs_ = smp.nOutputs_;
   randomize_ = smp.randomize_;
   nReplications_ = smp.nReplications_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0; i < nInputs_; i++)
   {
      lowerBounds_[i] = smp.lowerBounds_[i];
      upperBounds_[i] = smp.upperBounds_[i];
   }
   sampleMatrix_ = new double*[nSamples_];
   for (int i = 0; i < nSamples_; i++)
   {
      sampleMatrix_[i] = new double[nInputs_];
      for(int j = 0; j < nInputs_; j++)
         sampleMatrix_[i][j] = smp.sampleMatrix_[i][j];
   }
   sampleOutput_ = new double[nSamples_*nOutputs_];
   for (int i = 0; i < nSamples_*nOutputs_; i++)
      sampleOutput_[i] = smp.sampleOutput_[i];
   sampleStates_ = new int[nSamples_];
   for (int i = 0; i < nSamples_; i++)
      sampleStates_[i] = smp.sampleStates_[i];

   return *this;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
Sampling::~Sampling()
{
   int inputID;

   if (lowerBounds_ != NULL) delete [] lowerBounds_;
   if (upperBounds_ != NULL) delete [] upperBounds_;
   lowerBounds_ = upperBounds_ = NULL;
   deleteSampleData();
   if (symTable_ != NULL) delete [] symTable_;
   if (inputSettings_ != NULL)
   {
      for (inputID = 0; inputID < nInputs_; inputID++)
         if (inputSettings_[inputID] != NULL)
            delete [] inputSettings_[inputID];
      delete [] inputSettings_;
   }
   if (inputNumSettings_ != NULL) delete [] inputNumSettings_;
}

// ************************************************************************
// This function fully constructs a sampling from the psuadeIO_
// input.  All the common functionality it done in the Sampleing
// class, subclasses should overwrite this function.
// Currently this function does not support adaptive sampling or refinement, 
// although it probably could.  (The psuadeIO object has all the info for it)
// @param psuadeIO_ : A Psuade Data file, should have Input/Output variables 
// and Method defined.
// ------------------------------------------------------------------------
int Sampling::doSampling(PsuadeData* psuadeIO_) 
{
   pData  pPtr, pSymTable, pLowerB, pUpperB, pOutputs, pInputs;
   pData  pStates, pPDFs, pSettings;

   int outputLevel_ = 0;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;

   int nInputs = 0;
   int nOutputs = 0;
   double *iLowerB, *iUpperB = NULL;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   int *symTable = NULL;
   psuadeIO_->getParameter("input_symtable", pSymTable);
   symTable = pSymTable.intArray_;

   int    randomize, nReps, nSamples = 0;
   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("method_randomflag", pPtr);
   randomize = pPtr.intData_;
   psuadeIO_->getParameter("method_nreplications", pPtr);
   nReps = pPtr.intData_;  

   int    *inputPDFs, iSum, usePDFs, nRefines, samplingMethod = 0;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samplingMethod = pPtr.intData_;
   psuadeIO_->getParameter("input_use_input_pdfs", pPtr);
   usePDFs = pPtr.intData_;
   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefines = pPtr.intData_;

   int haveSettings = 0;
   psuadeIO_->getParameter("input_settings", pSettings);
   haveSettings = pSettings.nInts_;

   psuadeIO_->resetSamples();

   setPrintLevel(outputLevel_);
   setInputBounds(nInputs, iLowerB, iUpperB);
   setOutputParams(nOutputs);
   setInputParams(nInputs, NULL, NULL, symTable);
   setSamplingParams(nSamples, nReps, randomize);
   if (haveSettings) 
   {
      setInputParams(nInputs, pSettings.intArray_,
                    pSettings.dbleArray2D_, NULL);
   }

   for (int ii = 0; ii < nInputs; ii++) iSum += inputPDFs[ii];
   bool doingPDFs = false;
   if (iSum > 0 && usePDFs) doingPDFs = true;
  
   if (doingPDFs) 
   {
      if (nRefines > 0)
      {
         printf("PSUADE ERROR: both sample refinement and PDFs requested.\n");
         printf("       Sample refinement requires that all the input\n");
         printf("       probability distribution be uniform. You need to\n");
         printf("       either reset the number of refinements to be zero,\n");
         printf("       or not take out the probability distributions.\n");
         exit(1);
       }

       printf("INFO: Creating a sample where some uncertain parameters\n");
       printf("      are not uniformly distributed.\n");
       PDFManager *pdfman;
       pdfman = new PDFManager();
       pdfman->initialize(psuadeIO_);
       pdfman->genSample();
       delete pdfman;
    
       pData pPtr;
       double *sampleInputs, *sampleOutputs = NULL;
       int *sampleStates = NULL;

       psuadeIO_->getParameter("input_sample", pPtr);
       sampleInputs = pPtr.dbleArray_;
       psuadeIO_->getParameter("output_states", pPtr);
       sampleStates = pPtr.intArray_;
       psuadeIO_->getParameter("output_sample", pPtr);
       sampleOutputs = pPtr.dbleArray_;

       initialize(1);
       loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                   sampleOutputs, sampleStates);
   }
   else
   {
      printf("INFO: Creating a sample assuming all uncertain parameters\n");
      printf("      are uniformly distributed.\n");

      double *sampleInputs, *sampleOutputs = NULL;
      int *sampleStates = NULL;

      initialize(0);
      nSamples = getNumSamples(); /* may have been modified */
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                         sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                  sampleInputs,NULL,NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                   sampleStates,NULL); 
      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
   }
   return 0;
}

// ************************************************************************
// get print level
// ------------------------------------------------------------------------
int Sampling::setPrintLevel(int level)
{
   printLevel_ = level;
   return 0;
}

// ************************************************************************
// setInputBounds : set up the input lower and upper bounds
// ------------------------------------------------------------------------
int Sampling::setInputBounds(int nInputs, double *lower, double *upper)
{
   if (nInputs_ != 0 && nInputs != nInputs_)
   {
      printf("Sampling::setInputBounds ERROR : nInputs mismatch.\n");
      exit(1);
   }
   nInputs_ = nInputs;
   if (lowerBounds_  == NULL) delete [] lowerBounds_;
   if (upperBounds_  == NULL) delete [] upperBounds_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int jj = 0; jj < nInputs_; jj++) 
   {
      lowerBounds_[jj] = lower[jj];
      upperBounds_[jj] = upper[jj];
      if (lower[jj] >= upper[jj]) 
      {
         printf("Sampling::setInputBounds ERROR : lbound >= ubound.\n");
         printf("          Input %d : %e %e\n", jj+1, lower[jj], upper[jj]);
         exit(1);
      }    
   }
   return 0;
}

// ************************************************************************
// setSamplingParams : set up other sampling parameters
// ------------------------------------------------------------------------
int Sampling::setSamplingParams(int nSamples, int nReps, int randomize)
{
   if (nSamples_ != 0 && nSamples_ != nSamples)
   {
      printf("Sampling::setSamplingParams ERROR : nSamples mismatch.\n");
      printf("          nSamples = %d versus %d\n", nSamples, nSamples_);
      exit(1);
   }
   if (sampleMatrix_ != NULL)
   {
      for (int ii = 0; ii < nSamples_; ii++)
         if (sampleMatrix_[ii] != NULL)
            delete [] sampleMatrix_[ii];
      delete [] sampleMatrix_;
      sampleMatrix_ = NULL;
   }
   nSamples_ = nSamples;
   if (nReps > 0) nReplications_ = nReps;
   if (randomize > 0) randomize_ = randomize;
   return 0;
}

// ************************************************************************
// setInputParams : set up the input settings
// ------------------------------------------------------------------------
int Sampling::setInputParams(int, int *, double **, int *)
{
   return -1;
}

// ************************************************************************
// setOutputParams : set up the output settings
// ------------------------------------------------------------------------
int Sampling::setOutputParams(int nOutputs)
{
   if (nOutputs_ != 0 && nOutputs_ != nOutputs)
   {
      printf("Sampling::setOutputParams ERROR : nOutputs mismatch.\n");
      exit(1);
   }
   nOutputs_ = nOutputs;
   return 0;
}

// ************************************************************************
// load samples
// ------------------------------------------------------------------------
int Sampling::loadSamples(int nSamples, int nInputs, int nOutputs,
                       double *samples, double *outputs, int *states)
{
   int ii, jj, offset;

   deleteSampleData();
   if (nSamples != nSamples_) 
   {
      printf("Sampling::loadSamples WARNING: mismatched nSamples. \n");
      printf("          nSamples %d reset to incoming nSamples %d. \n",
             nSamples_, nSamples);
      nSamples_ = nSamples;
   }
   if (nInputs != nInputs_) 
   {
      printf("Sampling::loadSamples WARNING: invalid nInputs (%d,%d).\n",
             nInputs, nInputs_);
      exit(1);
   }
   if (nOutputs != nOutputs_) 
   {
      printf("Sampling::loadSamples WARNING: invalid nOutputs. \n");
      exit(1);
   }

   allocSampleData();
   for (ii = 0; ii < nSamples_; ii++) 
   {
      offset = ii * nInputs_;
      for (jj = 0; jj < nInputs_; jj++) 
         sampleMatrix_[ii][jj] = samples[offset+jj];
      sampleStates_[ii] = states[ii];
   }
   for (ii = 0; ii < nSamples_*nOutputs_; ii++)
      sampleOutput_[ii] = outputs[ii];
   return 0;
}

// ************************************************************************
// get number of samples
// ------------------------------------------------------------------------
int Sampling::getNumSamples()
{
   return nSamples_;
}

// ************************************************************************
// get number of inputs
// ------------------------------------------------------------------------
int Sampling::getNumInputs()
{
   return nInputs_;
}

// ************************************************************************
// get number of outputs
// ------------------------------------------------------------------------
int Sampling::getNumOutputs()
{
   return nOutputs_;
}

// ************************************************************************
// get all sample points 
// ------------------------------------------------------------------------
int Sampling::getSamples(int nSamples, int nInputs, int nOutputs,
                         double *sampData, double *outputs, int *states)
{
   int ii, jj, offset, offset2;

   if (nSamples != nSamples_) 
   {
      printf("Sampling::getSamples ERROR : invalid nSamples. \n");
      exit(1);
   }
   if (nInputs != nInputs_) 
   {
      printf("Sampling::getSamples ERROR : invalid nInputs. \n");
      exit(1);
   }
   if (nOutputs_ != nOutputs) 
   {
      printf("Sampling::getSamples ERROR : invalid nOutputs. \n");
      exit(1);
   }
   if (sampleMatrix_ == NULL)
   {
      printf("Sampling::getSamples ERROR : samples not ready. \n");
      exit(1);
   }
   if (sampleOutput_ == NULL && outputs != NULL)
   {
      printf("Sampling::getSamples ERROR : sample output not ready. \n");
      exit(1);
   }
   if (sampleStates_ == NULL && states != NULL)
   {
      printf("Sampling::getSamples ERROR : sample states not ready. \n");
      exit(1);
   }

   for (ii = 0; ii < nSamples; ii++) 
   {
      offset = ii * nInputs;
      offset2 = ii * nOutputs;
      for (jj = 0; jj < nInputs; jj++) 
         sampData[offset+jj] = sampleMatrix_[ii][jj];
      if (outputs != NULL)
      {
         for (jj = 0; jj < nOutputs; jj++) 
            outputs[offset2+jj] = sampleOutput_[offset2+jj];
      }
      if (states != NULL) states[ii] = sampleStates_[ii];
   }
   return nSamples;
}

// ************************************************************************
// get all sample points (this is for library call - no state information)
// ------------------------------------------------------------------------
int Sampling::getSamples(int nSamples, int nInputs, int nOutputs,
                         double *sampData, double *outputs)
{
   int ii, jj, offset, offset2;

   if (nSamples != nSamples_) 
   {
      printf("Sampling::getSamples ERROR : invalid nSamples. \n");
      exit(1);
   }
   if (nInputs != nInputs_) 
   {
      printf("Sampling::getSamples ERROR : invalid nInputs. \n");
      exit(1);
   }
   if (nOutputs_ != nOutputs) 
   {
      printf("Sampling::getSamples ERROR : invalid nOutputs. \n");
      exit(1);
   }
   if (sampleMatrix_ == NULL)
   {
      printf("Sampling::getSamples ERROR : samples not ready. \n");
      exit(1);
   }
   if (sampleOutput_ == NULL && outputs != NULL)
   {
      printf("Sampling::getSamples ERROR : sample output not ready. \n");
      exit(1);
   }

   for (ii = 0; ii < nSamples; ii++) 
   {
      offset = ii * nInputs;
      offset2 = ii * nOutputs;
      for (jj = 0; jj < nInputs; jj++) 
         sampData[offset+jj] = sampleMatrix_[ii][jj];
      if (outputs != NULL)
      {
         for (jj = 0; jj < nOutputs; jj++) 
            outputs[offset2+jj] = sampleOutput_[offset2+jj];
      }
   }
   return nSamples;
}

// ************************************************************************
// store sample outputs
// ------------------------------------------------------------------------
int Sampling::storeResult(int sampleID, int nOutputs, double *outputs, 
                          int *state)
{
   int jj;

   if (sampleID < 0 || sampleID >= nSamples_) 
   {
      printf("Sampling::storeResult ERROR : invalid sample ID. \n");
      exit(1);
   }
   if (nOutputs != nOutputs_) 
   {
      printf("Sampling::storeResult ERROR : invalid nOutputs. \n");
      exit(1);
   }
   if (sampleOutput_ == NULL || sampleStates_ == NULL)
   {
      printf("Sampling::storeResult ERROR : no storage. \n");
      exit(1);
   }

   for (jj = 0; jj < nOutputs_; jj++)
      sampleOutput_[sampleID*nOutputs_+jj] = outputs[jj];
   sampleStates_[sampleID] = 1;
   (*state) = 1;
   return 0;
}

// ************************************************************************
// perform refinement 
// ------------------------------------------------------------------------
int Sampling::refine(int, int, double, int, double *)
{
   return -1;
}

// ************************************************************************
// perform repair 
// ------------------------------------------------------------------------
int Sampling::repair(char *, int)
{
   return -1;
}

// ************************************************************************
// allocate sample matrix, output, and states
// ------------------------------------------------------------------------
int Sampling::allocSampleData()
{
   int ii;
   sampleMatrix_ = new double*[nSamples_];
   for (ii = 0; ii < nSamples_; ii++)
      sampleMatrix_[ii] = new double[nInputs_];
   sampleOutput_ = new double[nSamples_*nOutputs_];
   sampleStates_ = new int[nSamples_];
   for (ii = 0; ii < nSamples_*nOutputs_; ii++)
      sampleOutput_[ii] = PSUADE_UNDEFINED;
   for (ii = 0; ii < nSamples_; ii++)
      sampleStates_[ii] = 0;
   return 0;
}

// ************************************************************************
// clean the sample matrix, output, and status
// ------------------------------------------------------------------------
int Sampling::deleteSampleData()
{
   if (sampleMatrix_ != NULL)
   {
      for (int ii = 0; ii < nSamples_; ii++)
         if (sampleMatrix_[ii] != NULL)
            delete [] sampleMatrix_[ii];
      delete [] sampleMatrix_;
      sampleMatrix_ = NULL;
   }
   if (sampleOutput_ != NULL) delete [] sampleOutput_;
   if (sampleStates_ != NULL) delete [] sampleStates_;
   sampleOutput_ = NULL;
   sampleStates_ = NULL;
   return 0;
}

// ************************************************************************
// set parameter
// ------------------------------------------------------------------------
int Sampling::setParam(char *)
{
   return -1;
}

// ************************************************************************
// create a sampling scheme from identifier
// ------------------------------------------------------------------------
Sampling *SamplingCreateFromID(int samplingMethod)
{
   Sampling *sampler;
   char     sparam[101];

   switch (samplingMethod)
   {
      case PSUADE_SAMP_MC: 
           sampler = (Sampling *) new MCSampling();
           break;
      case PSUADE_SAMP_FACT:
           sampler = (Sampling *) new FactorialSampling();
           break;
      case PSUADE_SAMP_LHS:
           sampler = (Sampling *) new LHSampling();
           break;
      case PSUADE_SAMP_OA:
           sampler = (Sampling *) new OASampling();
           break;
      case PSUADE_SAMP_OALH:
           sampler = (Sampling *) new OALHSampling();
           break;
      case PSUADE_SAMP_MOAT:
           sampler = (Sampling *) new MOATSampling();
           break;
      case PSUADE_SAMP_SOBOL:
           sampler = (Sampling *) new SobolSampling();
           break;
      case PSUADE_SAMP_LPTAU:
           sampler = (Sampling *) new LPtauSampling();
           break;
      case PSUADE_SAMP_METIS:
           sampler = (Sampling *) new MetisSampling();
           break;
      case PSUADE_SAMP_FAST:
           sampler = (Sampling *) new FASTSampling();
           break;
      case PSUADE_SAMP_BBD:
           sampler = (Sampling *) new BoxBehnkenSampling();
           break;
      case PSUADE_SAMP_PBD:
           sampler = (Sampling *) new PlackettBurmanSampling();
           break;
      case PSUADE_SAMP_FF4:
           sampler = (Sampling *) new FractFactSampling();
           strcpy(sparam, "setResolution 4");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_FF5:
           sampler = (Sampling *) new FractFactSampling();
           strcpy(sparam, "setResolution 5");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCI4:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 4");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCI5:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 5");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCIF:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 0");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCF4:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 4");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCF5:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 5");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCFF:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 0");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCC4:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 4");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 2");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCC5:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 5");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 2");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCCF:
           sampler = (Sampling *) new CentralCompositeSampling();
           strcpy(sparam, "setResolution 0");
           sampler->setParam(sparam);
           strcpy(sparam, "setScheme 2");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_SFAST:
           sampler = (Sampling *) new SFASTSampling();
           break;
      case PSUADE_SAMP_UMETIS:
           sampler = (Sampling *) new UserMetisSampling();
           break;
      case PSUADE_SAMP_GMOAT:
           sampler = (Sampling *) new GMOATSampling();
           break;
      case PSUADE_SAMP_GMETIS:
           sampler = (Sampling *) new GMetisSampling();
           break;
      case PSUADE_SAMP_SG:
           sampler = (Sampling *) new SparseGridSampling();
           break;
      case PSUADE_SAMP_DISCRETE:
           sampler = (Sampling *) new DiscreteSampling();
           break;
      case PSUADE_SAMP_LSA:
           sampler = (Sampling *) new LSASampling();
           break;
      case PSUADE_SAMP_RFF4:
           sampler = (Sampling *) new RFractFactSampling();
           strcpy(sparam, "setResolution 4");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_RFF5:
           sampler = (Sampling *) new RFractFactSampling();
           strcpy(sparam, "setResolution 5");
           sampler->setParam(sparam);
           break;
      default :
           printf("PSUADE::SamplingCreateFromID ERROR - ");
           printf("invalid sampling method (%d).\n",samplingMethod);
           exit(1);
           break;
   }
   return sampler;
}

// ************************************************************************
// destroy a sampling object 
// ------------------------------------------------------------------------
int SamplingDestroy(Sampling *sampler)
{
  MCSampling               *MCSampler;
  FactorialSampling        *FactSampler;
  LHSampling               *LHSampler;
  OASampling               *OASampler;
  OALHSampling             *OALHSampler;
  MOATSampling             *MOATSampler;
  SobolSampling            *SobolSampler;
  LPtauSampling            *LPtauSampler;
  MetisSampling            *MetisSampler;
  FASTSampling             *FastSampler;
  SFASTSampling            *SFastSampler;
  BoxBehnkenSampling       *BBDSampler;
  PlackettBurmanSampling   *PBDSampler;
  FractFactSampling        *FFSampler;
  CentralCompositeSampling *CCDSampler;
  UserMetisSampling        *UMetisSampler;
  GMOATSampling            *GMOATSampler;
  GMetisSampling           *GMetisSampler;
  SparseGridSampling       *SparseGridSampler;
  DiscreteSampling         *DiscreteSampler;
  LSASampling              *LSASampler;
  RFractFactSampling       *RFFSampler;

  switch (sampler->samplingID_)
  {
    case PSUADE_SAMP_MC: 
         MCSampler = (MCSampling *) sampler;
         delete MCSampler;
         break;
    case PSUADE_SAMP_FACT:
         FactSampler = (FactorialSampling *) sampler;
         delete FactSampler;
         break;
    case PSUADE_SAMP_LHS:
         LHSampler = (LHSampling *) sampler;
         delete LHSampler;
         break;
    case PSUADE_SAMP_OA:
         OASampler = (OASampling *) sampler;
         delete OASampler;
         break;
    case PSUADE_SAMP_OALH:
         OALHSampler = (OALHSampling *) sampler;
         delete OALHSampler;
         break;
    case PSUADE_SAMP_MOAT:
          MOATSampler = (MOATSampling *) sampler;
          delete MOATSampler;
          break;
    case PSUADE_SAMP_SOBOL:
         SobolSampler = (SobolSampling *) sampler;
         delete SobolSampler;
         break;
    case PSUADE_SAMP_LPTAU:
         LPtauSampler = (LPtauSampling *) sampler;
         delete LPtauSampler;
         break;
    case PSUADE_SAMP_METIS:
         MetisSampler = (MetisSampling *) sampler;
         delete MetisSampler;
         break;
    case PSUADE_SAMP_FAST:
         FastSampler = (FASTSampling *) sampler;
         delete FastSampler;
         break;
    case PSUADE_SAMP_BBD:
         BBDSampler = (BoxBehnkenSampling *) sampler;
         delete BBDSampler;
         break;
    case PSUADE_SAMP_PBD:
         PBDSampler = (PlackettBurmanSampling *) sampler;
         delete PBDSampler;
         break;
    case PSUADE_SAMP_FF4:
         FFSampler = (FractFactSampling *) sampler;
         delete FFSampler;
         break;
    case PSUADE_SAMP_FF5:
         FFSampler = (FractFactSampling *) sampler;
         delete FFSampler;
         break;
    case PSUADE_SAMP_CCI4:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCI5:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCIF:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCF4:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCF5:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCFF:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCC4:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCC5:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_CCCF:
         CCDSampler = (CentralCompositeSampling *) sampler;
         delete CCDSampler;
         break;
    case PSUADE_SAMP_SFAST:
         SFastSampler = (SFASTSampling *) sampler;
         delete SFastSampler;
         break;
    case PSUADE_SAMP_UMETIS:
         UMetisSampler = (UserMetisSampling *) sampler;
         delete UMetisSampler;
         break;
    case PSUADE_SAMP_GMOAT:
         GMOATSampler = (GMOATSampling *) sampler;
         delete GMOATSampler;
         break;
    case PSUADE_SAMP_GMETIS:
         GMetisSampler = (GMetisSampling *) sampler;
         delete GMetisSampler;
         break;
    case PSUADE_SAMP_SG:
         SparseGridSampler = (SparseGridSampling *) sampler;
         delete SparseGridSampler;
         break;
    case PSUADE_SAMP_DISCRETE:
         DiscreteSampler = (DiscreteSampling *) sampler;
         delete DiscreteSampler;
         break;
    case PSUADE_SAMP_LSA:
         LSASampler = (LSASampling *) sampler;
         delete LSASampler;
         break;
    case PSUADE_SAMP_RFF4:
         RFFSampler = (RFractFactSampling *) sampler;
         delete RFFSampler;  
         break;
    case PSUADE_SAMP_RFF5:
         RFFSampler = (RFractFactSampling *) sampler;
         delete RFFSampler;
         break;
  }
  return 0;
}

// ************************************************************************
// compute sample quality
// ------------------------------------------------------------------------
int SamplingQuality(int nSamples, int nInputs, double *sampleInputs,
                    double *lbnds, double *ubnds, double *quality)
{
  int    ss1, ss2, ii, numPerms;
  double distance, minDistance, dtemp, minMinDist, avgDist, *darray;

  minMinDist =  1.0e35;
  for (ss1 = 0; ss1 < nSamples; ss1++)
  {
    minDistance = 1.0e35;
    for (ss2 = 0; ss2 < nSamples; ss2++)
    {
      if (ss1 != ss2)
      {
        distance = 0.0;
        for (ii = 0; ii < nInputs; ii++)
        {
          dtemp = sampleInputs[ss1*nInputs+ii] - 
                  sampleInputs[ss2*nInputs+ii];
          distance += pow(dtemp/(ubnds[ii]-lbnds[ii]), 2.0);
        }
        distance = sqrt(distance);
        if (distance < minDistance) minDistance = distance;
      }
    }
    if (minDistance < minMinDist) minMinDist = minDistance;
  }
  quality[0] = minMinDist;

  numPerms = 1;
  for (ii = 0; ii < nInputs; ii++) numPerms *= 2;
  darray = new double[nInputs];
  for (ss1 = 0; ss1 < numPerms; ss1++)
  {
    for (ii = 0; ii < nInputs; ii++) 
    {
      if (ss1 & (1 << ii)) darray[ii] = 1.0;
      else                 darray[ii] = 0.0;
    }
    minDistance = 1.0e35;
    for (ss2 = 0; ss2 < nSamples; ss2++)
    {
      distance = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
        dtemp = (sampleInputs[ss2*nInputs+ii] - lbnds[ii]) /
                (ubnds[ii] - lbnds[ii]);
        dtemp -= darray[ii];
        distance += dtemp * dtemp;
      }
      distance = sqrt(distance);
      if (distance < minDistance) minDistance = distance;
    } 
    avgDist += minDistance;
  }
  quality[1] = minDistance / (double) numPerms;
 
  avgDist =  0.0;
  for (ss1 = 0; ss1 < nSamples; ss1++)
  {
    for (ss2 = 0; ss2 < nSamples; ss2++)
    {
      if (ss1 != ss2)
      {
        distance = 0.0;
        for (ii = 0; ii < nInputs; ii++)
        {
          dtemp = sampleInputs[ss1*nInputs+ii] - 
                  sampleInputs[ss2*nInputs+ii];
          distance += pow(dtemp/(ubnds[ii]-lbnds[ii]), 2.0);
        }
        distance = sqrt(distance);
        avgDist += distance;
      }
    }
  }
  avgDist /= (1.0 * (nSamples) * (nSamples - 1.0));
  quality[2] = avgDist;
  delete [] darray;
  return 0;
}

// ************************************************************************
// ask user for a sampling method 
// ------------------------------------------------------------------------
int getSamplingMethod(char *pString)
{
  printDashes(PL_INFO, 0);
  printf("Available sampling methods: \n");
  printDashes(PL_INFO, 0);
  printf("1. Monte Carlo sampling\n");
  printf("2. Full Factorial sampling\n", PSUADE_SAMP_FACT);
  printf("3. Latin hypercube sampling\n", PSUADE_SAMP_LHS);
  printf("4. Orthogonal array sampling\n", PSUADE_SAMP_OA);
  printf("5. Morris sampling\n", PSUADE_SAMP_MOAT);
  printf("6. Quasi Monte Carlo (LPTAU) sampling\n",PSUADE_SAMP_LPTAU);
  printf("7. Full space-filling (METIS) sampling\n",PSUADE_SAMP_METIS);
  printf("8. Sparse grid sampling\n",PSUADE_SAMP_SG);
  int samType = getInt(1, 8, pString);
  if      (samType == 1) samType = PSUADE_SAMP_MC;
  else if (samType == 2) samType = PSUADE_SAMP_FACT;
  else if (samType == 3) samType = PSUADE_SAMP_LHS;
  else if (samType == 4) samType = PSUADE_SAMP_OA;
  else if (samType == 5) samType = PSUADE_SAMP_MOAT;
  else if (samType == 6) samType = PSUADE_SAMP_LPTAU;
  else if (samType == 7) samType = PSUADE_SAMP_METIS;
  else if (samType == 8) samType = PSUADE_SAMP_SG;
  return samType;
}



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
#include "Util/sysdef.h"
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

   if (nSamples != nSamples_) 
   {
      printf("Sampling::loadSamples WARNING: mismatched nSamples. \n");
      printf("          nSamples reset to incoming nSamples. \n");
      nSamples_ = nSamples;
   }
   if (nInputs != nInputs_) 
   {
      printf("Sampling::loadSamples WARNING: invalid nInputs. \n");
      exit(1);
   }
   if (nOutputs != nOutputs_) 
   {
      printf("Sampling::loadSamples WARNING: invalid nOutputs. \n");
      exit(1);
   }
   deleteSampleData();

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
int Sampling::setParam(string)
{
   return -1;
}

// ************************************************************************
// create a sampling scheme from identifier
// ------------------------------------------------------------------------
Sampling *SamplingCreateFromID(int samplingMethod)
{
   Sampling *sampler;
   string   sparam;

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
           sparam.clear();
           sparam.append("setResolution 4");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_FF5:
           sampler = (Sampling *) new FractFactSampling();
           sparam.clear();
           sparam.append("setResolution 5");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCI4:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 4");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCI5:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 5");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCIF:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 0");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 0");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCF4:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 4");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCF5:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 5");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCFF:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 0");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 1");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCC4:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 4");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 2");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCC5:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 5");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 2");
           sampler->setParam(sparam);
           break;
      case PSUADE_SAMP_CCCF:
           sampler = (Sampling *) new CentralCompositeSampling();
           sparam.clear();
           sparam.append("setResolution 0");
           sampler->setParam(sparam);
           sparam.clear();
           sparam.append("setScheme 2");
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
   }
   return 0;
}

// ************************************************************************
// compute sample quality
// ------------------------------------------------------------------------
double SamplingQuality(int nSamples, int nInputs, double *sampleInputs)
{
   int    ss1, ss2, ii;
   double distance, minDistance, minMaxDist, dtemp;

   minMaxDist = -1.0e35;
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
               distance += dtemp * dtemp;
            }
            distance = sqrt(distance);
            if (distance < minDistance) minDistance = distance;
         }
      }
      if (minDistance > minMaxDist) minMaxDist = minDistance;
   }
   return minMaxDist;
}


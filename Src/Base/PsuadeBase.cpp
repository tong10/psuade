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
// Functions for the class PsuadeBase
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "PsuadeBase.h"
#include "dtype.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"

// ------------------------------------------------------------------------
// local includes : sampling methods and input distributions
// ------------------------------------------------------------------------
#include "PDFManager.h"
#include "PDFBase.h"

// ------------------------------------------------------------------------
// local includes : function approximator
// ------------------------------------------------------------------------
#include "FuncApprox.h"

// ------------------------------------------------------------------------
// local includes : optimizers
// ------------------------------------------------------------------------
#include "Optimizers/Optimizer.h"

// ------------------------------------------------------------------------
// local includes : python
// ------------------------------------------------------------------------
#ifdef HAVE_PYTHON
#include "Python.h"
#endif

// ------------------------------------------------------------------------
// local defines 
// ------------------------------------------------------------------------
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PsuadeBase::PsuadeBase()
{
   outputLevel_  = 0;
   sampler_      = NULL;
   psuadeIO_     = NULL;
   useRSModel_   = 0;
   optimalCount_ = 0;
   optimalXData_ = NULL;
   optimalYData_ = NULL;
   optInitXData_ = NULL;
   optInitYData_ = NULL;
   jobsCompleted = -1;
   psuade_stop   = 0;

#ifdef HAVE_PYTHON
   update_gui   = NULL;
   yesno_dialog = NULL;
#endif
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PsuadeBase::~PsuadeBase()
{
   cleanUp();

#ifdef HAVE_PYTHON
   if (update_gui != NULL) Py_DECREF(update_gui);
   update_gui = NULL;
   if (yesno_dialog != NULL) Py_DECREF(yesno_dialog);
   yesno_dialog = NULL;
#endif
}

// ************************************************************************
// get data and parameters from PsuadeData object
// ------------------------------------------------------------------------
int PsuadeBase::getInputFromFile(const char *fname, const char *psuadeio_filename)
{
   int    status, nInputs, nOutputs, *sampleStates, samplingMethod;
   int    randomize, nReps, *symTable, nSamples, nRefines, flag, ii;
   int    *inputPDFs, iSum, usePDFs;
   double *sampleInputs, *iLowerB, *iUpperB, *sampleOutputs, haveSettings;
   char   cString[500];
   string sparam;
   pData  pPtr, pSymTable, pLowerB, pUpperB, pOutputs, pInputs, pStates;
   pData  pSettings, pPDFs, pOutputs2, pInputs2, pStates2;
   PDFManager *pdfman;

   cleanUp();
   if (psuadeIO_ != NULL) delete psuadeIO_;
   psuadeIO_ = new PsuadeData();
   status = psuadeIO_->readPsuadeFile(fname);
   if (status != 0) return 1;
   if (psuadeio_filename != NULL)
      psuadeIO_->readPsuadeIO(psuadeio_filename);

   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;
   if (outputLevel_ >= 4) printf("PSUADE::getInputFromFile begins.\n");
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("input_symtable", pSymTable);
   symTable = pSymTable.intArray_;
   psuadeIO_->getParameter("input_settings", pSettings);
   haveSettings = pSettings.nInts_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samplingMethod = pPtr.intData_;
   psuadeIO_->getParameter("method_nreplications", pPtr);
   nReps = pPtr.intData_;
   psuadeIO_->getParameter("method_randomflag", pPtr);
   randomize = pPtr.intData_;
   psuadeIO_->getParameter("input_use_input_pdfs", pPtr);
   usePDFs = pPtr.intData_;
   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;

   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefines = pPtr.intData_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;
   psuadeIO_->getParameter("output_sample", pOutputs);
   sampleOutputs = pOutputs.dbleArray_;

   if (sampleOutputs != NULL)
   {
      iSum = 0;
      for (ii = 0; ii < nInputs; ii++) iSum += inputPDFs[ii];
      if (iSum > 0 && nRefines > 0 && usePDFs == 1)
      {
         printf("PSUADE ERROR: you have requested sample refinement.\n");
         printf("       Sample refinement requires that all the input\n");
         printf("       probability distributions be uniform.\n");
         exit(1);
      }
      psuadeIO_->getParameter("input_sample", pInputs);
      sampleInputs = pInputs.dbleArray_;
      psuadeIO_->getParameter("output_states", pStates);
      sampleStates = pStates.intArray_;
      if (nRefines == 0 && samplingMethod == -1) 
         samplingMethod = PSUADE_SAMP_MC;
      sampler_ = (Sampling *) SamplingCreateFromID(samplingMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setInputParams(nInputs, NULL, NULL, symTable);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, nReps, randomize);
      sampler_->setParam(sparam);
      sampler_->initialize(1);
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                            sampleOutputs, sampleStates);
   }
   else
   {
     sampler_ = (Sampling *) SamplingCreateFromID(samplingMethod);
     sampler_->doSampling(psuadeIO_);
   }
   psuadeIO_->writePsuadeFile(NULL,1);

   if (outputLevel_ >= 4) printf("PSUADE::getInputFromFile ends.\n");
   return 0;
}

// ************************************************************************
// run integrated design and analysis
// ------------------------------------------------------------------------
int PsuadeBase::run() throw(Psuade_Stop_Exception)
{
   int   refineType, anaMethod;
   char  *appName, inString[200];
   FILE  *fp;
   pData pAppFiles, pPtr;

   if (psuadeIO_ == NULL)
   {
      printf("PSUADE::run ERROR - no PsuadeData object.\n");
      exit(1);
   }

   psuadeIO_->getParameter("app_files", pAppFiles);
   appName = pAppFiles.strArray_[0];
   fp = fopen(appName, "r");
   if (fp != NULL)
   {
      fscanf(fp, "%10c", inString);
      if (!strncmp(inString, "PSUADE_IO",9)) useRSModel_ = 1;
      fclose(fp);
   }

   psuade_stop = 0;

   psuadeIO_->getParameter("method_refine_type", pPtr);
   refineType = pPtr.intData_;
   psuadeIO_->getParameter("ana_method", pPtr);
   anaMethod = pPtr.intData_;

   if (refineType == 1 && anaMethod == PSUADE_ANA_ONESIGMA) 
        runAdaptiveErrBased0();
   else if (anaMethod == PSUADE_ANA_ARSM)
        runAdaptiveNN();
   else if (anaMethod == PSUADE_ANA_ARSMMB)
        runAdaptiveErrBased1();
   else if (anaMethod == PSUADE_ANA_ARSMMBBS)
        runAdaptiveErrBasedG();
   else if (anaMethod == PSUADE_ANA_REL)
        runAdaptivePRA();
   else if (anaMethod == PSUADE_ANA_AOPT)
        runAdaptiveOpt();
   else if (anaMethod == PSUADE_ANA_GLSA)
        runAdaptiveGradBased();
   else runUniform();

   if (psuadeIO_ != NULL) psuadeIO_->processOutputData();
   return 0;
}

// ************************************************************************
// interactive session
// ------------------------------------------------------------------------
int PsuadeBase::analyzeInteractive()
{
  setPrintLevelTS(3);
  interpretInteractive();
  return 0;
}

// ************************************************************************
// ************************************************************************

// ************************************************************************
// run the sample points on a single processor
// ------------------------------------------------------------------------
int PsuadeBase::runUniform()
{
   int     nInputs, nSamples, nOutputs, oldNSamples, sampleID, refineFlag=1;
   int     *sampleStates, status, nRefinements;
   int     iR, sampleGraphics, maxJobWaitTime, maxParallelJobs;
   int     parallelJobCount, count, runType, iteration;
   int     *refineNSamples, nReUsed, ii, askFlag=1;
   int     minJobWaitTime, maxState, mm, analysisMethod, jobsCompletedLast;
   int     launchInterval, saveFrequency, launchOnly, noAnalysis=0;
   int     limitedJobCount, optSwitch=0, randomize, refineRatio;
   long    nlong=0, *sampleStatesLong=NULL;
   double  *sampleInputs, *sampleOutputs, *tempOutputs, *tempInputs;
   double  *iLowerB, *iUpperB, analysisThreshold;
   double  *optData[4], refineThreshold, *sampleErrors;
   char    systemCommand[200], *appDriver, winput[500];
   FILE    *fp;
   pData   pPtr, pLowerB, pUpperB, pAppFiles, pStates;
   pData   pInpData, pOutData;
   FunctionInterface *funcIO=NULL;

   if (outputLevel_ >= 4) printf("PSUADE::run begins.\n");
   if (psuadeIO_ == NULL)
   {
      printf("PSUADE run: ERROR - no PsuadeData object.\n");
      exit(1);
   }

   psuadeIO_->getParameter("app_runtype", pPtr);
   runType = pPtr.intData_;
   if (runType & 8)
   {
      printf("PSUADE run: INFO - GENERATE INPUT FILE ONLY MODE.\n");
      funcIO = createFunctionInterfaceSimplified(psuadeIO_);
      psuadeIO_->getParameter("input_sample", pInpData);
      sampleInputs = pInpData.dbleArray_;
      psuadeIO_->getParameter("output_sample", pOutData);
      sampleOutputs = pOutData.dbleArray_;
      psuadeIO_->getParameter("output_states", pStates);
      sampleStates = pStates.intArray_;
      psuadeIO_->getParameter("method_nsamples", pPtr);
      nSamples = pPtr.intData_;
      psuadeIO_->getParameter("input_ninputs", pPtr);
      nInputs = pPtr.intData_;
      psuadeIO_->getParameter("output_noutputs", pPtr);
      nOutputs = pPtr.intData_;
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      for (sampleID = 0; sampleID < nSamples; sampleID++) 
         status = funcIO->evaluate(sampleID, nInputs,
                          &sampleInputs[sampleID*nInputs], nOutputs, 
                          &sampleOutputs[sampleID*nOutputs],1);   
      delete funcIO;
      return 0;
   }

   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_randomflag", pPtr);
   randomize = pPtr.intData_;
   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;
   psuadeIO_->getParameter("ana_graphicsflag", pPtr);
   sampleGraphics = pPtr.intData_;
   psuadeIO_->getParameter("ana_method", pPtr);
   analysisMethod = pPtr.intData_;
   psuadeIO_->getParameter("app_files", pAppFiles);
   appDriver = pAppFiles.strArray_[0];
   psuadeIO_->getParameter("ana_threshold", pPtr);
   analysisThreshold = pPtr.dbleData_;
   psuadeIO_->getParameter("ana_opt_switch", pPtr);
   optSwitch = pPtr.intData_;

   if ((runType & 2) && ((runType & 4) == 0))
   {
      launchOnly = 1;
      maxParallelJobs = 100000;
      if (outputLevel_ > 0) 
      {
         printf("PSUADE run: launch_only, max parallel jobs set to 100000.\n");
         printf("            launch_interval has been set to %d seconds\n",
                launchInterval);
      }
      outputLevel_ = 3;
   }
   else if (runType & 4)
   {
      launchOnly = 1;
      if (outputLevel_ > 0) 
      {
         printf("PSUADE run: limited_launch_only mode, max parallel jobs = %d.\n",
                maxParallelJobs);
         printf("            launch_interval has been set to %d seconds\n",
                launchInterval);
      }
   }
   else launchOnly = 0;

   if (!strcmp(appDriver, "NONE"))
   {
      psuadeIO_->getParameter("input_sample", pInpData);
      sampleInputs = pInpData.dbleArray_;
      psuadeIO_->getParameter("output_sample", pOutData);
      sampleOutputs = pOutData.dbleArray_;
      psuadeIO_->getParameter("output_states", pStates);
      sampleStates = pStates.intArray_;
      psuadeIO_->getParameter("method_nsamples", pPtr);
      nSamples = pPtr.intData_;
      psuadeIO_->getParameter("input_ninputs", pPtr);
      nInputs = pPtr.intData_;
      psuadeIO_->getParameter("output_noutputs", pPtr);
      nOutputs = pPtr.intData_;
      count = 0;
      for (sampleID = 0; sampleID < nSamples; sampleID++)
         if (sampleStates[sampleID] == 1) count++;
      if (count < nSamples)
      {
         if (optSwitch == 1)
         {
            if (outputLevel_ > 1) 
            {
               printf("PSUADE INFO: no driver given. The sample points\n");  
               printf("             will be used for optimization only.\n");  
            }
            for (sampleID = 0; sampleID < nSamples; sampleID++)
            {
               if (sampleStates[sampleID] == 0)
               {
                  for (ii = 0; ii < nOutputs; ii++)
                     sampleOutputs[sampleID*nOutputs+ii] = PSUADE_UNDEFINED;
                  sampleStates[sampleID] = 1;
               }
            }
         }
         else printf("PSUADE WARNING: no driver given.\n");  
         nRefinements = 0;
         noAnalysis = 1;
      }
      funcIO = createFunctionInterfaceSimplified(psuadeIO_);
   }
   else
   {
      psuadeIO_->getParameter("input_sample", pInpData);
      sampleInputs = pInpData.dbleArray_;
      psuadeIO_->getParameter("output_sample", pOutData);
      sampleOutputs = pOutData.dbleArray_;
      psuadeIO_->getParameter("output_states", pStates);
      sampleStates = pStates.intArray_;
      psuadeIO_->getParameter("method_nsamples", pPtr);
      nSamples = pPtr.intData_;
      psuadeIO_->getParameter("input_ninputs", pPtr);
      nInputs = pPtr.intData_;
      psuadeIO_->getParameter("output_noutputs", pPtr);
      nOutputs = pPtr.intData_;
      count = 0;
      for (sampleID = 0; sampleID < nSamples; sampleID++)
         if (sampleStates[sampleID] == 1) count++;
      if (count < nSamples || nRefinements > 0)
      {
         if (outputLevel_ > 0) 
            printf("PSUADE run: creating interface to user driver.\n");
         funcIO = createFunctionInterface(psuadeIO_);
      }
      else funcIO = createFunctionInterfaceSimplified(psuadeIO_);

      psuadeIO_->getParameter("ana_opt_switch", pPtr);
   }

   if (useRSModel_)
   {
      maxParallelJobs = 1;
      minJobWaitTime  = 0;
      maxJobWaitTime  = 0;
   }

   analysisManager_.clearAnalyzers();
   analysisManager_.setup(psuadeIO_);

   if ((sampleGraphics & 2) != 0) sampleGraphics = 1;
   else                           sampleGraphics = 0;
   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);
   funcIO->setOutputLevel(outputLevel_);

   if (outputLevel_ > 0) 
   {
      printf("PSUADE run: output level = %d\n", outputLevel_);
      printf("PSUADE run: max parallel jobs = %d\n", maxParallelJobs);
      printf("PSUADE run: max job wait time = %d seconds\n", maxJobWaitTime);
      printf("PSUADE run: min job wait time = %d seconds\n", minJobWaitTime);
      printf("PSUADE run: launch interval   = %d seconds\n", launchInterval);
      printf("PSUADE run: save frequency    = every %d runs\n", saveFrequency);
      printf("NOTE: if evaluation should be fast but is slow, check save frequency.\n");
      printf("Note: use psuade_pmachine to dynamically change max jobs.\n");
      printf("Note: use psuade_stop to terminate gracefully.\n");
   }
   refineNSamples = new int[nRefinements+1];
   for (iR = 0; iR < nRefinements+1; iR++)
   {
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 0) 
      {
         printEquals(0);
         if (nRefinements > 0)
            printf("PSUADE run: refinement %d(out of %d), nSamples = %d\n",
                   iR, nRefinements, nSamples);
         else
            printf("PSUADE run: running sample, nSamples = %d \n", nSamples);
      }

      refineNSamples[iR] = nSamples;
      psuadeIO_->getParameter("input_sample", pInpData);
      sampleInputs = pInpData.dbleArray_;
      psuadeIO_->getParameter("output_sample", pOutData);
      sampleOutputs = pOutData.dbleArray_;
      psuadeIO_->getParameter("output_states", pStates);
      sampleStates = pStates.intArray_;
      psuadeIO_->getParameter("method_nsamples", pPtr);
      nSamples = pPtr.intData_;
      psuadeIO_->getParameter("input_ninputs", pPtr);
      nInputs = pPtr.intData_;
      psuadeIO_->getParameter("output_noutputs", pPtr);
      nOutputs = pPtr.intData_;

      jobsCompleted = 0;
      for (sampleID = 0; sampleID < nSamples; sampleID++) 
      {
         if (sampleStates[sampleID] == 1) jobsCompleted++;
         else                             sampleStates[sampleID] = 0;
      }

      if ((noAnalysis == 0) && (sampleGraphics > 0) && (nInputs >= 2))
      {
         nlong = nSamples;
         sampleStatesLong = new long[nSamples];
         for (sampleID = 0; sampleID < nSamples; sampleID++) 
            sampleStatesLong[sampleID] = (long) sampleStates[sampleID]; 
         tempInputs = new double[nSamples];
         tempOutputs = new double[nSamples];
         for (sampleID = 0; sampleID < nSamples; sampleID++) 
         {
            tempInputs[sampleID] = sampleInputs[sampleID*nInputs];
            tempOutputs[sampleID] = sampleInputs[sampleID*nInputs+1];
         }
         Plotbegin(iLowerB[0],iUpperB[0],iLowerB[1], iUpperB[1]);
      }

      parallelJobCount = limitedJobCount = iteration = 0;
      while (noAnalysis == 0 && jobsCompleted < nSamples)
      {
         jobsCompletedLast = jobsCompleted;
         iteration++;
         for (sampleID = 0; sampleID < nSamples; sampleID++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
               psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                              sampleStates,NULL); 
	       psuadeIO_->writePsuadeFile(NULL,0);
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[sampleID] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               if ((minJobWaitTime > 0) && (runType & 1) == 0)
                    status = compareSamples(sampleID,nSamples,nInputs,
                                            sampleInputs, sampleStates);
               else status = -1;

               if (status < 0) /* repeated sample not found */
               {
                  if ((runType & 4))
                  {
                     if (limitedJobCount < maxParallelJobs)
                     {
                        if (outputLevel_ > 2)
                           printf("Limited launch: job = %6d\n",sampleID+1);
                        status = funcIO->evaluate(sampleID,nInputs,
                                 &sampleInputs[sampleID*nInputs], nOutputs, 
                                 &sampleOutputs[sampleID*nOutputs],0);   
                        limitedJobCount++;
                     }
                     else
                     {
                        if (outputLevel_ > 2)
                           printf("Limited launch: creating jobfile %6d\n",
                                  sampleID+1);
                        status = funcIO->evaluate(sampleID,nInputs,
                                 &sampleInputs[sampleID*nInputs], nOutputs, 
                                 &sampleOutputs[sampleID*nOutputs],1);   
                     }
                  }
                  else
                  {
                     if (outputLevel_ > 2)
                        printf("Launch: job = %6d\n",sampleID+1);
                     status = funcIO->evaluate(sampleID,nInputs,
                                &sampleInputs[sampleID*nInputs], nOutputs, 
                                &sampleOutputs[sampleID*nOutputs],0);   
                     if ((runType & 2) == 0) parallelJobCount++;
                  }
               }
               else /* repeat found, copy the results or wait */
               {
                  if (outputLevel_ > 5)
                  {
                     for (mm = 0; mm < nInputs; mm++)
                        printf("(%4d,%4d) : %12.4e %12.4e\n",status+1,
                               sampleID+1,sampleInputs[status*nInputs+mm],
                               sampleInputs[sampleID*nInputs+mm]);
                     printf("PSUADE run: sample %d repeated with %d.\n",
                            sampleID+1,status+1);
                  }
                  if (sampleStates[status] == 1)
                  {
                     for (ii = 0; ii < nOutputs; ii++)
                        sampleOutputs[sampleID*nOutputs+ii] =
                            sampleOutputs[status*nOutputs+ii];
                     status = 0;
                  }
                  else status = -(status + 1);
               }

               if (status == 0) 
               {
                  sampler_->storeResult(sampleID, nOutputs,
                                   &(sampleOutputs[sampleID*nOutputs]),
                                   &(sampleStates[sampleID]));
                  if (parallelJobCount > 0) parallelJobCount--;
                  if (outputLevel_ > 3)
                  {
                     printf("Completed: job = %6d\n", sampleID+1);
                     for (mm = 0; mm < nInputs; mm++)
                        printf("\t\tinput data %3d = %e\n", mm+1,
                               sampleInputs[sampleID*nInputs+mm]);
                     for (mm = 0; mm < nOutputs; mm++)
                        printf("\t\toutput data %2d = %e\n",
                               mm+1, sampleOutputs[sampleID*nOutputs+mm]);
                  }
                  jobsCompleted++;
                  if (((jobsCompleted-jobsCompletedLast) % saveFrequency) == 
                      (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                                    sampleStates,NULL); 
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }

                  if (outputLevel_ > 0)
                  {
                     if (((sampleID+1) % 100) == 0) 
                       printf("\nSample point %6d completed (out of %d).\n",
                               sampleID+1, nSamples);
                     else if ((sampleID % 11) == 1) printf(".");
                     fflush(stdout);
                  }
                  if ((sampleID % 1) == 0) 
                  {
                     fp = fopen("psuade_stop","r");
                     if (fp != NULL)
                     {
                        fclose(fp);
                        printf("\nTermination signal detected (1).\n");
                        psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                                       sampleStates,NULL); 
                        psuadeIO_->writePsuadeFile(NULL,0);
			throw Psuade_Stop_Exception();
                     }
                  }
               }
               else // status > 0 or launchOnly = 1
               {
                  sampleStates[sampleID] = status; /* running or waiting */
                  if ((status > 0) && (launchOnly == 1)) 
                  {
                     sprintf(systemCommand, "sleep %d", launchInterval);
                     system(systemCommand);
                  }
               }
            }
            else if (sampleStates[sampleID] >= 2)
            {
               status = funcIO->evaluate(sampleID,nInputs,
                              &sampleInputs[sampleID*nInputs], nOutputs, 
                              &sampleOutputs[sampleID*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(sampleID, nOutputs,
                                   &(sampleOutputs[sampleID*nOutputs]),
                                   &(sampleStates[sampleID]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else 
               {
                  sampleStates[sampleID]++;
                  if (outputLevel_ > 0) 
                     printf("Waiting for Job %d to complete (status = %d)\n",
                            sampleID+1, sampleStates[sampleID]);
               } 
               if (minJobWaitTime > 0 &&
                   (sampleStates[sampleID]-2)*minJobWaitTime > maxJobWaitTime)
               {
                  sampleStates[sampleID] = 0;
                  parallelJobCount--;
                  sampleID--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE run: sample %6d to be restarted.\n",
                            sampleID+1);
               }
            }

            if ((sampleGraphics > 0) && (nInputs >= 2) && (nlong != 0))
            {
               sampleStatesLong[sampleID] = (long) sampleStates[sampleID];
               PlotSamples2D(nlong,tempInputs,tempOutputs,sampleStatesLong);
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         if (jobsCompletedLast < jobsCompleted)
         {
            fillInSamples(nSamples,nOutputs,sampleOutputs,sampleStates);
            psuadeIO_->writePsuadeFile(NULL,0);
         }

         maxState = 0;
         for (sampleID = 0; sampleID < nSamples; sampleID++)
            if (sampleStates[sampleID] > maxState)
               maxState = sampleStates[sampleID];
         if ((maxState > 100) && (minJobWaitTime <= 0))
         {
            minJobWaitTime  = 1;
            if (outputLevel_ > 0) 
               printf("PSUADE run: min job wait time re-set to %d.\n",
                      minJobWaitTime);
         }
         if ((maxState > 200) && (minJobWaitTime >  0))
         {
            minJobWaitTime *= 2;
            if (outputLevel_ > 0) 
               printf("PSUADE run: min job wait time re-set to %d.\n",
                      minJobWaitTime);
         }

         if (launchOnly) break;

         if ((jobsCompleted < nSamples) && (minJobWaitTime > 0))
         {
            if (outputLevel_ > 0) 
               printf("PSUADE run: sleep for %d seconds.\n",minJobWaitTime); 
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
            if (minJobWaitTime > 100)
            {
               fp = fopen("psuadeStatus", "a");
               if (fp != NULL)
               {
                  fprintf(fp,"\nRefinement %d, iteration %d: \n",iR,iteration);
                  for (ii = 0; ii < nSamples; ii++)
                     fprintf(fp, "Sample %7d = %d\n",ii+1,sampleStates[ii]);
                  fclose(fp);
               }
               printf("PSUADE sample state information in psuadeStatus.\n");
            }
         }

         fp = fopen("psuade_pmachine","r");
         if (fp != NULL)
         {
            fscanf(fp, "%d", &status);
            if ((status > 0) && (status < 200))
            {
               maxParallelJobs = status;
               if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
            }
            printf("PSUADE run: psuade_pmachine found, max jobs reset to %d.\n",
                   maxParallelJobs);
            fclose(fp);
         }

         fp = fopen("psuade_restart","r");
         if (fp != NULL)
         {
            fscanf(fp, "%d", &count);
            for (ii = 0; ii < count; ii++)
            {
               fscanf(fp, "%d", &sampleID);
               if (sampleID >= 1 && sampleID <= nSamples) 
               {
                  sampleID--;
                  if (sampleStates[sampleID] == 1) 
                  {
                     printf("Sample %d completed, no need to restart.\n",
                            sampleID+1);
                  }
                  else
                  {
                     if (sampleStates[sampleID] != 0 && parallelJobCount > 0)
                        parallelJobCount--;
                     sampleStates[sampleID] = 0; 
                     printf("Sample %d to be restarted.\n", sampleID+1);
                  }
               }
            }
            fclose(fp);
         }

         fp = fopen("psuade_stop","r");
         if (fp != NULL)
         {
            fclose(fp);
            printf("\nTermination signal detected (2).\n");
	    psuadeIO_->writePsuadeFile(NULL,0);
	    throw Psuade_Stop_Exception();
         }
         if (outputLevel_ > 0) 
         {
            printf("\nPSUADE run: jobs completed = %d(out of %d)\n",
                   jobsCompleted, nSamples);
            if (jobsCompleted == 0 && askFlag == 1 && iteration > 10 &&
                outputLevel_ >= 4)
            {
#ifdef HAVE_PYTHON
               if (yesno_dialog != NULL) {
                  PyObject* temp;
	          temp = PyObject_CallFunction( yesno_dialog, "ss",
		     "No Progress", "No progress has been made.\n"
		     "Something may be wrong.\n\n"
                     "Do you want to proceed anyway?" );
                  if (temp != NULL && PyObject_Not(temp) == 1) {
		     Py_DECREF(temp);
		     throw Psuade_Stop_Exception();
		  }
	          else if (temp != NULL) Py_DECREF(temp);
               }
#else
	       printf("The following message is for precautionary purpose:\n");
	       printf("No progress has been made. Something may be wrong if\n");
	       printf("your job should take less than a second to run, and\n");
	       printf("in that case you should answer no to the following\n");
               printf("question. Otherwise, answer yes.\n");
               printf("You can turn off this checking by setting diagnostics\n");
               printf("level to less than 4. You can avoid this complaint by\n");
               printf("setting the minimum wait time to be your estimated\n");
               printf("time per run.\n");
               printf("Do you want to proceed? (y or n) ");
               scanf("%s", winput);
               if (winput[0] == 'n') throw Psuade_Stop_Exception();
#endif
               askFlag = 0;
            }
         }
      }

      if ((noAnalysis == 0) && (sampleGraphics > 0) && (nInputs >= 2))
      {
         delete [] sampleStatesLong;
         delete [] tempInputs;
         delete [] tempOutputs;
         Plotend();
      }

      if (launchOnly && (jobsCompleted < nSamples))
      {
         if (outputLevel_ > 0) 
            printf("PSUADE run: terminate - launch only mode.\n");
         break;
      }

      if (noAnalysis == 0)
      {
         refineFlag = analysisManager_.analyze(psuadeIO_,iR+1,refineNSamples,-1);
      }
      else refineFlag = 1;

      if (noAnalysis == 0) pgPlotResponseSurface();

      status = 0;
      if (optSwitch == 1)
      {
         if (optInitXData_ != NULL) delete [] optInitXData_;
         if (optInitYData_ != NULL) delete [] optInitYData_;
         if (optimalYData_ != NULL) delete [] optimalYData_;
         if (optimalXData_ != NULL) delete [] optimalXData_;
         if (optimalYData_ != NULL) delete [] optimalYData_;
         optimalCount_ = 0;
         optData[0] = NULL;;
         optData[1] = NULL;;
         optData[2] = NULL;
         optData[3] = NULL;
         status = OptimizerSearch(psuadeIO_, funcIO, optData, &optimalCount_);
         optInitXData_ = optData[0];
         optInitYData_ = optData[1];
         optimalXData_ = optData[2];
         optimalYData_ = optData[3];
      }
      if (status > 0) refineFlag = 0;

      if ((refineFlag == 1) && (nRefinements > 0) && (iR < nRefinements))
      {
         sampleErrors = analysisManager_.getSampleErrors();
         refineRatio = 2;
         refineThreshold = analysisThreshold;
         sampler_->refine(refineRatio,randomize,refineThreshold,nSamples,
                          sampleErrors);


         oldNSamples = nSamples;
         nSamples = sampler_->getNumSamples();
         if (nSamples == oldNSamples)
         {
             refineFlag = 0;
             break;
         }
         sampleInputs  = new double[nSamples*nInputs];
         sampleOutputs = new double[nSamples*nOutputs];
         sampleStates  = new int[nSamples];
            
         sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                              sampleOutputs, sampleStates);
         nReUsed = 0;
         for (sampleID = 0; sampleID < nSamples; sampleID++)
            if (sampleStates[sampleID] == 1) nReUsed++;
         printf("PSUADE run: refinement - number reused = %d(out of %d)\n", 
                nReUsed, oldNSamples);

         psuadeIO_->updateMethodSection(-1,nSamples,-1,nRefinements-iR-1,-1);
         psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                       sampleInputs,NULL);
         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
      }
      if (outputLevel_ > -1) 
      {
          if (nRefinements > 0) 
             printf("PSUADE run: refinements completed = %d (out of %d)\n",
                    iR, nRefinements);
          printEquals(0);
      }
      if (refineFlag == 0) break;
   }

   delete funcIO;
   delete [] refineNSamples;
   if (outputLevel_ >= 4) printf("PSUADE run: exiting...\n");
   return 0;
}

// ************************************************************************
// run the special adaptive mode (GP or ANN with METIS)
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveErrBased0()
{
   int    faType, samMethod, nInputs, nOutputs, nRefinements, ss, ii, iSum;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, iR = 0, nSamples, *refineNSamples;
   int    loopFlag, currNJobs, *sampleStates, nJobsDiff;
   int    parallelJobCount, status, maxState, jobsCompletedLast, maxInd;
   int    refineFlag, initFlag, count, usedFlag=0, refineRatio, randomize;
   int    *inputPDFs;
   double *iLowerB, *iUpperB, *sampleErrors;
   double *sampleInputs, *sampleOutputs, maxError, refineThreshold;
   char   systemCommand[100];
   pData  pPtr, pLowerB, pUpperB, pPDFs;
   FunctionInterface *funcIO;

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ii = 0; ii < nInputs; ii++) iSum += inputPDFs[ii];
   if (iSum)
   {
      printf("PSUADE ERROR: adaptive(0) does not currenlty allow \n");
      printf("       non-uniform probability distribution functions\n");
      printf("       to be used in INPUT section.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   printf("PSUADE adaptive(0): GP, Kriging, or ANN with METIS.\n");
   psuadeIO_->getParameter("ana_rstype", pPtr);
   faType = pPtr.intData_;
   if (faType != PSUADE_RS_GP1   && faType != PSUADE_RS_ANN &&
       faType != PSUADE_RS_MARSB && faType != PSUADE_RS_KR)
   {
      printf("PSUADE adaptive(0) INFO: response surface should be one of\n");
      printf("                         GP, ANN, Kriging, or MarsBag.\n");
      printf("                         Default to Kriging.\n");
      faType = PSUADE_RS_KR;
      psuadeIO_->updateAnalysisSection(-1, -1, faType, -1, -1, -1);
   }
   psuadeIO_->getParameter("method_sampling", pPtr);
   samMethod = pPtr.intData_;
   if (samMethod != PSUADE_SAMP_METIS && samMethod != PSUADE_SAMP_FACT)
   {
      printf("PSUADE adaptive(0): sampling defaulted to METIS.\n");
      samMethod = PSUADE_SAMP_METIS; 
      psuadeIO_->updateMethodSection(samMethod,-1,-1,-1,-1);
      initFlag = 1;
      if (sampler_ != NULL) SamplingDestroy(sampler_);
      sampler_ = NULL;
   }
   else initFlag = 0;

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;

   if (initFlag == 1)
   {
      sampler_ = (Sampling *) SamplingCreateFromID(samMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, 1, 1);
      sampler_->initialize(0);
      nSamples = sampler_->getNumSamples(); /* may have been modified */
      psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
   }
   else
   {
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      count = 0;
      for (ss = 0; ss < nSamples; ss++)
         if (sampleStates[ss] == -1) count++;
      if (count > 0) usedFlag = 1;
      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
   }

   funcIO = createFunctionInterface(psuadeIO_);

   analysisManager_.clearAnalyzers();
   analysisManager_.setup(psuadeIO_);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   refineRatio = 2;
   randomize = 1;
   refineThreshold = 1.0;
   refineNSamples = new int[nRefinements+1];
   refineNSamples[0] = nSamples;
   if (usedFlag == 0)
   {
      for (iR = 0; iR < nRefinements; iR++)
      {
         sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
         nSamples = sampler_->getNumSamples();
         refineNSamples[iR+1] = nSamples;
         if (nSamples > 2000) 
         {
            printf("PSUADE adaptive(0): max nSamples clamped to %d\n",
                   nSamples);
            break;
         }
      }
      nRefinements = iR;
      psuadeIO_->updateMethodSection(-1,nSamples,-1,nRefinements,-1);
   }

   sampleInputs  = new double[nSamples * nInputs];
   sampleOutputs = new double[nSamples * nOutputs];
   sampleStates  = new int[nSamples];
   sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                        sampleOutputs, sampleStates);
   if (usedFlag == 0)
   {
      for (ss = 0; ss < refineNSamples[0]; ss++) sampleStates[ss] = 0; 
      for (ss = refineNSamples[0]; ss < nSamples; ss++) sampleStates[ss]=-1; 
   }
   psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                 sampleInputs,NULL); 
   psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                  sampleStates,NULL);
   psuadeIO_->writePsuadeFile(NULL,0);

   loopFlag = 1;
   jobsCompleted = 0;
   for (ss = 0; ss < nSamples; ss++) 
      if (sampleStates[ss] == 1) jobsCompleted++; 
   while (loopFlag)
   {
      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) 
         if (sampleStates[ss] == 0) currNJobs++;

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
               delete [] refineNSamples;
               delete [] sampleInputs;
               delete [] sampleOutputs;
               delete [] sampleStates;
               delete funcIO;
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptive(0): sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("\nPSUADE adaptive(0): jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      refineFlag = analysisManager_.analyze(psuadeIO_,iR+1,refineNSamples,-1);
      if (refineFlag == 0) break;

      sampleErrors = analysisManager_.getSampleErrors();
      maxInd = -1;
      maxError = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         if (sampleStates[ss] < 0 && sampleErrors[ss] > maxError)
         {
            maxError = sampleErrors[ss];
            maxInd   = ss;
         }
      }
      if (maxInd >= 0)
      {
         sampleStates[maxInd] = 0;
         if (outputLevel_ > 1)
         {
            printf("PSUADE adaptive(0): new data point (error = %12.4e) = \n",
                   sampleErrors[maxInd]);
            for (ii = 0; ii < nInputs; ii++)
               printf("\tX%3d = %12.4e\n",ii+1,sampleInputs[maxInd*nInputs+ii]);
         }
      }
   }

   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   delete funcIO;
   delete [] refineNSamples;
   return 0;
}

// ************************************************************************
// run the special adaptive mode for response surface analysis (METIS)
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveNN()
{
   int    samMethod, initFlag, refineLevel, nInputs, nOutputs, nRefinements;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, ss, nSamples, refineRatio, randomize, iSum, *inputPDFs;
   int    loopFlag, currNJobs, *sampleStates, nJobsDiff, ii, jj, marsMode=0;
   int    parallelJobCount, status, maxState, jobsCompletedLast, refineType;
   int    lastNSamples, rstype, nPtsPerDim=64, length, refineSize;
   int    tstNSamples=0, tstNInputs, tstNOutputs, ivar1, ivar2, numMars=100;
   int    marsNSamples;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs, outData, errAvg;
   double refineThreshold=1.0, errMax, dtemp;
   double anaThreshold, errRMS, *tstSamInputs, *tstSamOutputs;
   double totalSum, *tstOutputs, **marsDataX, **marsDataY;
   char   systemCommand[100], cString[100], winput[500];
   char   *targv[6];
   pData  pPtr, pLowerB, pUpperB, pPDFs, pTstInputs, pTstOutputs;
   string sparam;
   FILE   *fp;
   FuncApprox        *faPtr=NULL;
   FunctionInterface *funcIO=NULL;
   PsuadeData        *tstIO=NULL;

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ss = 0; ss < nInputs; ss++) iSum += inputPDFs[ss];
   if (iSum)
   {
      printf("PSUADE ERROR: adaptiveRSM does not currently allow \n");
      printf("       non-uniform probability distributions to be\n");
      printf("       defined in the INPUT SECTION.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   printAsterisks(0);
   printf("PSUADE adaptiveNN: Use Metis sampling.\n");
   printf("Adaptive sampling based on discrepany of the output\n");
   printf("with its neighbor's output (select maximum).\n");
   printDashes(0);
   printf("Note: Turn on rs_expert mode to set RS parameters.\n");
   printf("      Turn on outputLevel (>0) to get error histogram.\n");
   printEquals(0);
   psuadeIO_->getParameter("method_refine_type", pPtr);
   refineType = pPtr.intData_;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samMethod = pPtr.intData_;
   if (samMethod != PSUADE_SAMP_METIS)
   {
      printf("PSUADE adaptiveNN: sampling defaulted to METIS.\n");
      samMethod = PSUADE_SAMP_METIS; 
      psuadeIO_->updateMethodSection(samMethod,-1,-1,-1,-1);
      initFlag = 1;
      if (sampler_ != NULL) SamplingDestroy(sampler_);
      sampler_ = NULL;
   }
   else initFlag = 0;

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   if (nOutputs > 1)
   {
      printf("PSUADE adaptiveNN: nOutputs should be 1.\n");
      return 0;
   }
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;
   psuadeIO_->getParameter("ana_rstype", pPtr);
   rstype = pPtr.intData_;
   psuadeIO_->getParameter("ana_threshold", pPtr);
   anaThreshold = pPtr.dbleData_;

   printf("You may test the quality of the response surface using a test\n");
   printf("sample (in psuadeData format).\n");
   sprintf(cString, "Use a test sample ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y')
   {
      sprintf(cString, "Enter the test sample file name : ");
      getString(cString, winput);
      ss = strlen(winput);
      winput[ss-1] = '\0';
      tstIO = new PsuadeData();
      status = tstIO->readPsuadeFile(winput);
      if (status == -1) exit(1);
      tstIO->getParameter("method_nsamples", pPtr);
      tstNSamples = pPtr.intData_;
      tstIO->getParameter("input_ninputs", pPtr);
      tstNInputs = pPtr.intData_;
      tstIO->getParameter("output_noutputs", pPtr);
      tstNOutputs = pPtr.intData_;
      if (tstNInputs != nInputs)
      {
         printf("PSUADE adaptiveNN ERROR : test sample nInputs != %d\n",
                nInputs);
         delete tstIO;
         return 0;
      }
      if (tstNOutputs > 1)
      {
         printf("PSUADE adaptiveNN: test sample nOutputs != 1.\n");
         printf("       INFO: use 'write' in interactive model to select\n");
         printf("             1 output only and re-run this.\n");
         delete tstIO;
         return 0;
      }
      tstIO->getParameter("input_sample", pTstInputs);
      tstSamInputs = pTstInputs.dbleArray_;
      tstIO->getParameter("output_sample", pTstOutputs);
      tstSamOutputs = pTstOutputs.dbleArray_;
   }
   if (psRSExpertMode_ == 1)
   {
      if (rstype == PSUADE_RS_MARSB)
      {
         sprintf(cString, "Number of MARS (default = 100, >2, <502) = ");
         numMars = getInt(3, 501, cString);
         sprintf(cString, "Use mean (0) or median (1) of MarsBag : ");
         marsMode = getInt(0, 1, cString);
      }
   }

   if (initFlag == 1)
   {
      sampler_ = (Sampling *) SamplingCreateFromID(samMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, 1, 1);
      sampler_->initialize(0);
      nSamples = sampler_->getNumSamples();
      psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
   }
   sparam.clear();
   if (refineType == 0)
   {
      sparam.append("setUniformRefinement");
      refineSize = 100000;
   }
   else sparam.append("setAdaptiveRefinementBasedOnOutputs");
   sampler_->setParam(sparam);
   if (refineSize > 0)
   {
      sparam.clear();
      sprintf(cString, "setRefineSize %d", refineSize);
      sparam.append(cString);
      sampler_->setParam(sparam);
      printf("PSUADE adaptiveNN: refineSize = %d\n",refineSize);
   }
   refineRatio = 2;
   randomize = 1;

   funcIO = createFunctionInterface(psuadeIO_);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   nSamples = sampler_->getNumSamples();
   marsDataX = new double*[numMars];
   marsDataY = new double*[numMars];
   length = (nSamples+nRefinements*refineSize);
   for (ii = 0; ii < numMars; ii++)
   {
      marsDataX[ii] = new double[length*nInputs];
      marsDataY[ii] = new double[length];
   }

   loopFlag = 1;
   jobsCompleted = 0;
   refineLevel = 0;
   nSamples = -1;
   marsNSamples = 0;

   while (loopFlag)
   {
      lastNSamples = nSamples;
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 1)
      {
         printAsterisks(0);
         printf("PSUADE adaptiveNN: current level    = %d (of %d)\n",
                refineLevel, nRefinements);
         printf("PSUADE adaptiveNN: current nSamples = %d\n",nSamples);
      }
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                    sampleInputs,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                     sampleStates,NULL);
      psuadeIO_->writePsuadeFile(NULL,0);

      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss]==0) currNJobs++;
      if (psAnalysisInteractive_ == 1 && jobsCompleted < currNJobs)
      {
         printf("A sample has been created in the psuadeData file. \n");
         printf("At this point you can choose to run your model with\n");
         printf("this sample via psuade or by yourself (if the model\n");
         printf("is expensive to run, you want to choose the latter).\n");
         sprintf(systemCommand, "Run the model yourself (y or n) ? ");
         getString(systemCommand, cString);
         if (cString[0] == 'y')
         {
            printf("You have chosen to run the sample yourself.\n");
            printf("The following steps are for ARSM refinements:\n");
            printf("(1) Rename psuadeData to something else. \n");
            printf("(2) Comment out num_refinements in this file.\n");
            printf("(3) Change sampling method to be MC in this file.\n");
            printf("(4) Comment out analysis method = ARSM in this file.\n");
            printf("(5) Run the sample in this file and collect outputs.\n");
            printf("(6) Restore num_refinements in this file.\n");
            printf("(7) Restore sampling method = METIS in this file.\n");
            printf("(8) Restore analysis method = ARSM in this file.\n");
            printf("(9) Finally, restart psuade with this file.\n");
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            delete funcIO;
            for (ii = 0; ii < numMars; ii++)
            {
               delete [] marsDataX[ii];
               delete [] marsDataY[ii];
            }
            delete [] marsDataX;
            delete [] marsDataY;
            return 0;
         }
      }

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptiveNN: sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("\nPSUADE adaptiveNN: jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      errRMS = anaThreshold;
      if (lastNSamples > 0)
      {
         printf("PSUADE adaptiveNN: response surface analysis.\n");
         printf("   construct response surface with %d sample points.\n",
                lastNSamples);
         printf("   test response surface with previous %d sample points.\n",
                nSamples-lastNSamples);
         faPtr = genFA(rstype, nInputs, nOutputs, lastNSamples);
         if (faPtr == NULL)
         {
            printf("ERROR: cannot create function approximator.\n");
            exit(1);
         }
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);
         if (psRSExpertMode_ == 0 && 
             (rstype == PSUADE_RS_MARS || rstype == PSUADE_RS_MARSB))
         {
            strcpy(cString, "mars_params");
            targv[0] = (char *) cString;
            ivar1 = lastNSamples;
            targv[1] = (char *) &ivar1;
            ivar2 = 2 * nInputs / 3 + 1;
            targv[2] = (char *) &ivar2;
            faPtr->setParams(3, targv);
            if (rstype == PSUADE_RS_MARSB)
            {
               strcpy(cString, "num_mars");
               targv[0] = (char *) cString;
               targv[1] = (char *) &numMars;
               faPtr->setParams(2, targv);
               if (marsMode == 1)
               {
                  strcpy(cString, "median");
                  targv[0] = (char *) cString;
                  faPtr->setParams(1, targv);
               }
            }
         }
         length = -999;
         faPtr->genNDGridData(sampleInputs,sampleOutputs,&length,NULL,NULL);
         errMax = errAvg = errRMS = 0.0;
         totalSum = 0.0;
         for (ss = lastNSamples; ss < nSamples; ss++)
         {
            outData = faPtr->evaluatePoint(&sampleInputs[ss*nInputs]);
            if (outputLevel_ > 3)
            {
               printf("Data %5d : predicted = %12.4e, (actual) = %12.4e,",
                      ss, outData, sampleOutputs[ss]);
               printf(" diff = %12.4e\n", outData - sampleOutputs[ss]);
            }
            totalSum += PABS(sampleOutputs[ss]);
            dtemp   = outData - sampleOutputs[ss];
            errAvg += dtemp;
            errRMS += (dtemp * dtemp);
            dtemp   = PABS(dtemp);
            errMax = (dtemp > errMax) ? dtemp : errMax;
         }
         totalSum /= (double) (nSamples - lastNSamples);;
         errRMS    = sqrt(errRMS/(nSamples-lastNSamples));
         errAvg    = errAvg / (nSamples-lastNSamples);
         printf("     response surface unscaled max error = %e\n",errMax);
         printf("     response surface   scaled max error = %e\n",
                errMax/totalSum);
         printf("     response surface unscaled rms error = %e\n",errRMS);
         printf("     response surface   scaled rms error = %e\n",
                errRMS/totalSum);
         printf("     response surface unscaled avg error = %e\n",errAvg);
         printf("     response surface   scaled avg error = %e\n",
                errAvg/totalSum);
         delete faPtr;
      }
      if (tstNSamples > 0)
      {
         printEquals(0);
         faPtr = genFA(rstype, nInputs, nOutputs, nSamples);
         if(faPtr == NULL){
	   printf("function genFA returned NULL in file %s line %d, exiting\n", __FILE__, __LINE__);
           exit(1);
         }
         faPtr->setNPtsPerDim(nPtsPerDim);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(outputLevel_);

         if (psRSExpertMode_ == 0 && 
             (rstype == PSUADE_RS_MARS || rstype == PSUADE_RS_MARSB))
         {
            strcpy(cString, "mars_params");
            targv[0] = (char *) cString;
            targv[1] = (char *) &nSamples;
            ivar2 = 2 * nInputs / 3 + 1;
            targv[2] = (char *) &ivar2;
            faPtr->setParams(3, targv);
            if (rstype == PSUADE_RS_MARSB)
            {
               strcpy(cString, "num_mars");
               targv[0] = (char *) cString;
               targv[1] = (char *) &numMars;
               faPtr->setParams(2, targv);
               if (marsMode == 1)
               {
                  strcpy(cString, "median");
                  targv[0] = (char *) cString;
                  faPtr->setParams(1, targv);
               }
               for (ii = 0; ii < numMars; ii++)
               {
                  for (ss = 0; ss < nSamples-marsNSamples; ss++)
                  {
                     if (marsNSamples == 0)
                          ivar1 = PSUADE_rand() % 
                                  (nSamples-marsNSamples);
                     else ivar1 = ss;
                     ivar2 = ivar1 + marsNSamples;
                     for (jj = 0; jj < nInputs; jj++)
                     {
                        marsDataX[ii][(marsNSamples+ss)*nInputs+jj] =
                            sampleInputs[ivar2*nInputs+jj];
                     }
                     marsDataY[ii][marsNSamples+ss] = sampleOutputs[ivar2];
                  }
               }
               marsNSamples = nSamples;
               strcpy(cString, "mars_sample");
               targv[0] = (char *) cString;
               targv[2] = (char *) &marsNSamples;
               for (ii = 0; ii < numMars; ii++)
               {
                  targv[1] = (char *) &ii;
                  targv[3] = (char *) marsDataX[ii];
                  targv[4] = (char *) marsDataY[ii];
                  faPtr->setParams(5, targv);
               }
            }
         }
         length = -999;
         faPtr->genNDGridData(sampleInputs,sampleOutputs,&length,NULL,NULL);
         tstOutputs = new double[tstNSamples];
         faPtr->evaluatePoint(tstNSamples, tstSamInputs, tstOutputs);
         totalSum = errMax = errAvg = errRMS =0.0;
         for (ss = 0; ss < tstNSamples; ss++)
         {
            totalSum += PABS(tstOutputs[ss]);
            dtemp   = tstOutputs[ss] - tstSamOutputs[ss];
            errAvg += dtemp;
            errRMS  += dtemp * dtemp;
            dtemp   = PABS(dtemp);
            errMax  = (dtemp > errMax) ? dtemp : errMax;
         }
         errRMS = sqrt(errRMS / tstNSamples);
         totalSum /= (double) tstNSamples;
         errAvg  = errAvg / tstNSamples;
         printf("     test sample RS unscaled max error = %e\n",errMax);
         printf("     test sample RS   scaled max error = %e\n",
                errMax/totalSum);
         printf("     test sample RS unscaled rms error = %e\n",errRMS);
         printf("     test sample RS   scaled rms error = %e\n",
                errRMS/totalSum);
         printf("     test sample RS unscaled avg error = %e\n",errAvg);
         printf("     test sample RS   scaled avg error = %e\n",
                errAvg / totalSum);
         if (errRMS < anaThreshold || refineLevel >= nRefinements)
         {
            sprintf(cString, "arsm_nn_err.m");
            fp = fopen(cString, "w");
            if (fp != NULL)
            {
               fprintf(fp, "%% inputs, true outputs, predicted outputs\n");
               fprintf(fp, "clf\n");
               fprintf(fp, "A = [\n");
               for (ss = 0; ss < tstNSamples; ss++)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     fprintf(fp, "%e ", tstSamInputs[ss*nInputs+ii]);
                  fprintf(fp, "%e %e\n", tstSamOutputs[ss], tstOutputs[ss]);
               }
               fprintf(fp, "];\n");
               fwritePlotCLF(fp);
               fprintf(fp, "m = %d;\n", nInputs);
               fprintf(fp, "X1 = A(:,1);\n");
               fprintf(fp, "X2 = A(:,2);\n");
               fprintf(fp, "Y  = A(:,m+1) - A(:,m+2);\n");
               fprintf(fp, "subplot(2,2,1)\n");
               fprintf(fp, "plot3(X1,X2,Y,'*','markerSize',13')\n");
               fwritePlotXLabel(fp, "Input 1");
               fwritePlotYLabel(fp, "Input 2");
               fwritePlotZLabel(fp, "Output");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Errors w.r.t. Input 1 and 2");
               fprintf(fp, "subplot(2,2,2)\n");
               fprintf(fp, "plot(X1,Y,'*','markerSize',13')\n");
               fwritePlotXLabel(fp, "Input 1");
               fwritePlotYLabel(fp, "Output");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Output vs Input 1");
               fprintf(fp, "subplot(2,2,4)\n");
               fprintf(fp, "plot(X2,Y,'*','markerSize',13')\n");
               fwritePlotXLabel(fp, "Input 2");
               fwritePlotYLabel(fp, "Output");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Output vs Input 2");
               fprintf(fp, "AA = [\n");
               for (ss = 0; ss < nSamples; ss++)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     fprintf(fp, "%e ", sampleInputs[ss*nInputs+ii]);
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
               fprintf(fp, "XX1 = AA(:,1);\n");
               fprintf(fp, "XX2 = AA(:,2);\n");
               fprintf(fp, "subplot(2,2,3)\n");
               fprintf(fp, "plot(XX1,XX2,'*','markerSize',13')\n");
               fwritePlotXLabel(fp, "Input 1");
               fwritePlotYLabel(fp, "Input 2");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Input 1 vs 2 (original sample)");
               printf("PSUADE adaptiveNN: error plots are in %s.\n",
                      cString);
               fclose(fp);
            }
         }
         delete [] tstOutputs;
         delete faPtr;
      }

      refineLevel++;
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);
      if (errRMS < anaThreshold)
      {
         printf("PSUADE adaptiveNN: threshold reached (using unscaled rms).\n");
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         break;
      }
      if (psAnalysisInteractive_ == 1)
      {
         sprintf(systemCommand, "Do you want to quit now (y or n) ? ");
         getString(systemCommand, cString);
         if (cString[0] == 'y')
         {
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            break;
         }
      }
      if (refineLevel > nRefinements)
      {
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         break;
      }

      sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
      psuadeIO_->updateMethodSection(-1,-1,-1,nRefinements-1,-1);
      psuadeIO_->writePsuadeFile(NULL,0);

      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
   }

   delete funcIO;
   if (tstIO != NULL) delete tstIO;
   for (ii = 0; ii < numMars; ii++)
   {
      delete [] marsDataX[ii];
      delete [] marsDataY[ii];
   }
   delete [] marsDataX;
   delete [] marsDataY;
   return 0;
}

// ************************************************************************
// run the special adaptive mode for response surface analysis (METIS)
// using MARS with bagging as response surface
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveErrBased1()
{
   int    samMethod, initFlag, refineLevel, nInputs, nOutputs, nRefinements;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, nSamples, refineRatio, randomize, iSum, *inputPDFs;
   int    loopFlag, currNJobs, *sampleStates, nJobsDiff, *samStates2=NULL;
   int    parallelJobCount, status, maxState, jobsCompletedLast, refineType;
   int    rstype, nPtsPerDim=64, length, refineSize;
   int    nSamples2, *iArray, useRandomPts=0, auxNSamples=0, ivar1, ivar2;
   int    auxNInputs, auxNOutputs, tstNSamples=0, tstNInputs, tstNOutputs;
   int    ii, jj, ss, numMars=100, marsMode=0, marsNSamples;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs;
   double refineThreshold=1.0, dtemp, **marsDataX, **marsDataY;
   double anaThreshold, *samInputs2=NULL, *samOutputs2=NULL;
   double *samStds2=NULL, totalSum=0.0, errMax=0.0, errAvg=0.0, errL2=0.0;
   double *auxSamInputs=NULL, *auxSamOutputs=NULL;
   double *tstSamInputs=NULL, *tstSamOutputs=NULL, *tstOutputs;
   char   systemCommand[100], cString[100], winput[100], *targv[6];
   FILE   *fp;
   pData  pPtr, pLowerB, pUpperB, pPDFs, pAuxInputs, pAuxOutputs;
   pData  pTstInputs, pTstOutputs;
   string sparam;
   FuncApprox *faPtr=NULL;
   Sampling   *samplerAux=NULL;
   PsuadeData *auxIO=NULL, *tstIO=NULL;
   FunctionInterface *funcIO=NULL;

   printAsterisks(0);
   printf("PSUADE adaptive(1): Use Metis with MarsBagging/GP/Kriging.\n");
   printf("Adaptive sampling based on predicted errors from\n");
   printf("the MARS-with-bagging/GP/Kriging response surface.\n");
   printEquals(0);
   printf("To be able to exercise more control of this method,\n");
   printf("you can turn on interactive and/or expert modes in\n");
   printf("the input file.\n");
   printf("Turn on outputLevel (>0) to get error histogram.\n");
   printf("Turn on interactive to have stepwise control.\n");
   printEquals(0);

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ss = 0; ss < nInputs; ss++) iSum += inputPDFs[ss];
   if (iSum)
   {
      printf("PSUADE ERROR: adaptive(1) does not currently allow\n");
      printf("       non-uniform probability distribution to be\n");
      printf("       defined in the INPUT SECTION.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   printf("This random strategy allows comparing this adaptivity\n");
   printf("strategy with just using random sample points.\n");
   printf("To run this adaptivity scheme, answer 'n' below.\n");
   sprintf(cString, "Use random points for refinement ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y') useRandomPts = 1;

   psuadeIO_->getParameter("method_refine_type", pPtr);
   refineType = pPtr.intData_;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samMethod = pPtr.intData_;
   if (samMethod != PSUADE_SAMP_METIS)
   {
      printf("PSUADE adaptive(1): sampling defaulted to METIS.\n");
      samMethod = PSUADE_SAMP_METIS; 
      psuadeIO_->updateMethodSection(samMethod,-1,-1,-1,-1);
      initFlag = 1;
      if (sampler_ != NULL) SamplingDestroy(sampler_);
      sampler_ = NULL;
   }
   else initFlag = 0;

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   if (nOutputs > 1)
   {
      printf("PSUADE adaptive(1): nOutputs should be 1.\n");
      printf("       INFO: use 'write' in interactive model to select\n");
      printf("             1 output only and re-run this.\n");
      return 0;
   }
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;

   printf("You may test the quality of the response surface using a test\n");
   printf("sample (in psuadeData format).\n");
   sprintf(cString, "Use a test sample ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y')
   {
      sprintf(cString, "Enter the test sample file name : ");
      getString(cString, winput);
      ss = strlen(winput);
      winput[ss-1] = '\0';
      tstIO = new PsuadeData();
      status = tstIO->readPsuadeFile(winput);
      if (status == -1) exit(1);
      tstIO->getParameter("method_nsamples", pPtr);
      tstNSamples = pPtr.intData_;
      tstIO->getParameter("input_ninputs", pPtr);
      tstNInputs = pPtr.intData_;
      tstIO->getParameter("output_noutputs", pPtr);
      tstNOutputs = pPtr.intData_;
      if (tstNInputs != nInputs)
      {
         printf("PSUADE adaptive(1) ERROR : test sample nInputs != %d\n",
                nInputs);
         return 0;
      }
      if (tstNOutputs > 1)
      {
         printf("PSUADE adaptive(1): test sample nOutputs != 1.\n");
         printf("       INFO: use 'write' in interactive model to select\n");
         printf("             1 output only and re-run this.\n");
         return 0;
      }
      tstIO->getParameter("input_sample", pTstInputs);
      tstSamInputs = pTstInputs.dbleArray_;
      tstIO->getParameter("output_sample", pTstOutputs);
      tstSamOutputs = pTstOutputs.dbleArray_;
   }

   printf("You may add to the base sample an auxiliary sample which\n");
   printf("covers the corners (for example, factorial or fractional\n");
   printf("factorial).\n");
   sprintf(cString, "Add an auxiliary sample ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y')
   {
      sprintf(cString, "Enter auxiliary sample file name : ");
      getString(cString, winput);
      ss = strlen(winput);
      winput[ss-1] = '\0';
      auxIO = new PsuadeData();
      status = auxIO->readPsuadeFile(winput);
      if (status == -1) exit(1);
      auxIO->getParameter("method_nsamples", pPtr);
      auxNSamples = pPtr.intData_;
      auxIO->getParameter("input_ninputs", pPtr);
      auxNInputs = pPtr.intData_;
      auxIO->getParameter("output_noutputs", pPtr);
      auxNOutputs = pPtr.intData_;
      if (auxNInputs != nInputs)
      {
         printf("PSUADE adaptive(1) ERROR : auxiliary nInputs != %d\n",
                nInputs);
         return 0;
      }
      if (auxNOutputs > 1)
      {
         printf("PSUADE adaptive(1): auxiliary sample nOutputs != 1.\n");
         printf("       INFO: use 'write' in interactive model to select\n");
         printf("             1 output only and re-run this.\n");
         return 0;
      }
      auxIO->getParameter("input_sample", pAuxInputs);
      auxSamInputs = pAuxInputs.dbleArray_;
      samInputs2 = auxSamInputs;
      length = (nSamples+auxNSamples+nRefinements*refineSize)*nInputs;
      auxSamInputs = new double[length];
      for (ss = 0; ss < auxNSamples*nInputs; ss++)
         auxSamInputs[ss] = samInputs2[ss];
      auxIO->getParameter("output_sample", pAuxOutputs);
      samOutputs2 = pAuxOutputs.dbleArray_;
      length = (nSamples+auxNSamples+nRefinements*refineSize)*nOutputs;
      auxSamOutputs = new double[length];
      for (ss = 0; ss < auxNSamples*nInputs; ss++)
         auxSamOutputs[ss] = samOutputs2[ss];
   }

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;
   psuadeIO_->getParameter("ana_rstype", pPtr);
   rstype = pPtr.intData_;
   if (rstype != PSUADE_RS_MARSB && rstype != PSUADE_RS_GP1 &&
       rstype != PSUADE_RS_TGP   && rstype != PSUADE_RS_KR)
   {
      printf("PSUADE adaptive(1) INFO: RS type defaulted to MarsBag.\n");
      rstype = PSUADE_RS_MARSB; 
   }
   if (psRSExpertMode_ == 1)
   {
      if (rstype == PSUADE_RS_MARSB)
      {
         sprintf(cString, "Number of MARS (default = 100, >2, <502) = ");
         numMars = getInt(3, 501, cString);
         sprintf(cString, "Use mean (0) or median (1) of MarsBag : ");
         marsMode = getInt(0, 1, cString);
      }
   }
   psuadeIO_->getParameter("ana_threshold", pPtr);
   anaThreshold = pPtr.dbleData_;

   if (initFlag == 1)
   {
      sampler_ = (Sampling *) SamplingCreateFromID(samMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, 1, 1);
      sampler_->initialize(0);
      nSamples = sampler_->getNumSamples();
      psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
   }
   sparam.clear();
   if (refineType == 0)
   {
      sparam.append("setUniformRefinement");
      refineSize = 100000;
   }
   else sparam.append("setAdaptiveRefinementBasedOnErrors");
   sampler_->setParam(sparam);
   if (refineSize > 0)
   {
      sparam.clear();
      sprintf(cString, "setRefineSize %d", refineSize);
      sparam.append(cString);
      sampler_->setParam(sparam);
      printf("PSUADE adaptive(1): refineSize = %d\n",refineSize);
   }
   refineRatio = 2;
   randomize = 1;

   funcIO = createFunctionInterface(psuadeIO_);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   nSamples = sampler_->getNumSamples();
   marsDataX = new double*[numMars];
   marsDataY = new double*[numMars];
   length = (nSamples+auxNSamples+nRefinements*refineSize);
   for (ii = 0; ii < numMars; ii++)
   {
      marsDataX[ii] = new double[length*nInputs];
      marsDataY[ii] = new double[length];
   }
   for (ss = 0; ss < auxNSamples*nInputs; ss++)
      marsDataX[ii][ss] = auxSamInputs[ss];
   for (ss = 0; ss < auxNSamples; ss++)
      marsDataY[ii][ss] = auxSamOutputs[ss];

   loopFlag = 1;
   jobsCompleted = 0;
   refineLevel = 0;
   marsNSamples = auxNSamples;

   while (loopFlag)
   {
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 1)
      {
         printAsterisks(0);
         printf("PSUADE adaptive(1): current level    = %d (of %d)\n",
                refineLevel, nRefinements);
         printf("PSUADE adaptive(1): current nSamples = %d\n",
                nSamples);
      }
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                    sampleInputs,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                     sampleStates,NULL);
      psuadeIO_->writePsuadeFile(NULL,0);

      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss]==0) currNJobs++;
      if (psAnalysisInteractive_ == 1 && jobsCompleted < currNJobs)
      {
         printf("A sample has been created in the psuadeData file. \n");
         printf("At this point you can choose to run your model with\n");
         printf("this sample via psuade or by yourself (if the model\n");
         printf("is expensive to run, you want to choose the latter).\n");
         sprintf(systemCommand, "Run the model yourself (y or n) ? ");
         getString(systemCommand, cString);
         if (cString[0] == 'y')
         {
            printf("You have chosen to run the sample yourself.\n");
            printf("The following are the steps to continue ARSM refinements:\n");
            printf("(1) Rename psuadeData to something else (e.g. psData). \n");
            printf("(2) Run the sample in this file and collect outputs.\n");
            printf("(3) Replace the outputs in psData (do not change others).\n");
            printf("(4) Finally, restart psuade with this file (psData).\n");
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            delete funcIO;
            return 0;
         }
      }

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) 
            {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptive(1): sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("PSUADE adaptive(1): jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      samplerAux = (Sampling *) SamplingCreateFromID(samMethod);
      samplerAux->setInputBounds(nInputs, iLowerB, iUpperB);
      samplerAux->setOutputParams(nOutputs);
      samplerAux->setSamplingParams(nSamples, -1, -1);
      samplerAux->initialize(1);
      samplerAux->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                              sampleOutputs, sampleStates);
      sparam.clear();
      sparam.append("setUniformRefinement");
      samplerAux->setParam(sparam);
      samplerAux->refine(refineRatio,randomize,0,0,NULL);
      nSamples2 = samplerAux->getNumSamples();
      if (nSamples2 != 2 * nSamples)
      {
         printf("PSUADE adaptive(1): something is wrong. Consult developers.\n");
         printf("                    refined sample size != 2 * original size\n");
         delete samplerAux;
         delete sampler_;
         sampler_ = NULL;
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         delete funcIO;
         return 0;
      }
      samInputs2  = new double[nSamples2 * nInputs];
      samOutputs2 = new double[nSamples2 * nOutputs];
      samStates2  = new int[nSamples2];
      samStds2    = new double[nSamples2];
      samplerAux->getSamples(nSamples2, nInputs, nOutputs, samInputs2, 
                             samOutputs2, samStates2);

      printAsterisks(0);
      printf("PSUADE adaptive(1): response surface analysis.\n");
      printf("PSUADE adaptive(1): current level     = %d (of %d)\n",
             refineLevel, nRefinements);
      printf("                    training nSamples = %d\n",
             nSamples + auxNSamples);
      printf("                    test set nSamples = %d\n",
             nSamples2-nSamples);

      faPtr = genFA(rstype, nInputs, nOutputs, nSamples+auxNSamples);
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel_);
      length = -999;

      for (ii = 0; ii < numMars; ii++)
      {
         for (ss = 0; ss < nSamples-marsNSamples+auxNSamples; ss++)
         {
            if (marsNSamples == auxNSamples)
                 ivar1 = PSUADE_rand() % (nSamples-marsNSamples+auxNSamples);
            else ivar1 = ss;
            for (jj = 0; jj < nInputs; jj++)
               marsDataX[ii][(marsNSamples+ss)*nInputs+jj] =
                   sampleInputs[(ivar1+marsNSamples-auxNSamples)*nInputs+jj];
            marsDataY[ii][marsNSamples+ss] = 
                   sampleOutputs[ivar1+marsNSamples-auxNSamples];
         }
      }
      marsNSamples = nSamples + auxNSamples;

      if (psRSExpertMode_ == 0 && rstype == PSUADE_RS_MARSB)
      {
         strcpy(cString, "mars_params");
         targv[0] = (char *) cString;
         ivar1 = (nSamples + auxNSamples);
         targv[1] = (char *) &ivar1;
         ivar2 = 2 * nInputs / 3 + 1;
         targv[2] = (char *) &ivar2;
         faPtr->setParams(3, targv);
         strcpy(cString, "num_mars");
         targv[0] = (char *) cString;
         targv[1] = (char *) &numMars;
         faPtr->setParams(2, targv);
         if (marsMode == 1)
         {
            strcpy(cString, "median");
            targv[0] = (char *) cString;
            faPtr->setParams(1, targv);
         }
         strcpy(cString, "mars_sample");
         targv[0] = (char *) cString;
         targv[2] = (char *) &marsNSamples;
         for (ii = 0; ii < numMars; ii++)
         {
            targv[1] = (char *) &ii;
            targv[3] = (char *) marsDataX[ii];
            targv[4] = (char *) marsDataY[ii];
            faPtr->setParams(5, targv);
         }
      }
      faPtr->genNDGridData(sampleInputs,sampleOutputs,&length,NULL,NULL);

      faPtr->evaluatePointFuzzy(nSamples2-nSamples, 
                         &samInputs2[nInputs*nSamples], 
                         &samOutputs2[nSamples], &samStds2[nSamples]);
      if (tstNSamples == 0)
      {
         // extract the magnitude of each sample points
         totalSum = errMax = errAvg = errL2 =0.0;
         for (ss = nSamples; ss < nSamples2; ss++)
         {
            totalSum += PABS(samOutputs2[ss]);
            errMax = (PABS(samStds2[ss]) > errMax) ? PABS(samStds2[ss]) : errMax;
            errAvg += samStds2[ss];
            errL2  += pow(samStds2[ss], 2.0e0);
         }
         errL2 = sqrt(errL2 / (nSamples2 - nSamples));
         totalSum /= (nSamples2 - nSamples);
         errAvg  = sqrt(errAvg/(nSamples2-nSamples));
         printf("     response surface unscaled max error = %e\n",errMax);
         printf("     response surface   scaled max error = %e\n",errMax/totalSum);
         printf("     response surface unscaled rms error = %e\n",errL2);
         printf("     response surface   scaled rms error = %e\n",errL2/totalSum);
         printf("     response surface unscaled avg error = %e\n",errAvg);
         printf("     response surface   scaled avg error = %e\n",
                errAvg/totalSum);
      }
      if (tstNSamples > 0)
      {
         printEquals(0);
         tstOutputs = new double[tstNSamples];
         faPtr->evaluatePoint(tstNSamples, tstSamInputs, tstOutputs);
         totalSum = errMax = errAvg = errL2 =0.0;
         for (ss = 0; ss < tstNSamples; ss++)
         {
            totalSum += PABS(tstOutputs[ss]);
            dtemp = tstOutputs[ss] - tstSamOutputs[ss];
            errAvg += dtemp;
            dtemp = PABS(dtemp);
            errL2  += dtemp * dtemp;
            errMax = (dtemp > errMax) ? dtemp : errMax;
         }
         errL2 = sqrt(errL2 / tstNSamples);
         totalSum /= (double) tstNSamples;
         errAvg  = errAvg / tstNSamples;
         printf("     test sample RS unscaled max error = %e\n",errMax);
         printf("     test sample RS   scaled max error = %e\n",errMax/totalSum);
         printf("     test sample RS unscaled rms error = %e\n",errL2);
         printf("     test sample RS   scaled rms error = %e\n",errL2/totalSum);
         printf("     test sample RS unscaled avg error = %e\n",errAvg);
         printf("     test sample RS   scaled avg error = %e\n",
                errAvg/totalSum);
         if (outputLevel_ > 0)
         {
            sprintf(cString, "arsm_marsb_err.m");
            fp = fopen(cString, "w");
            fprintf(fp, "%% inputs, true outputs, predicted outputs\n");
            fprintf(fp, "clf\n");
            fprintf(fp, "A = [\n");
            for (ss = 0; ss < tstNSamples; ss++)
            {
               for (ii = 0; ii < nInputs; ii++)
                  fprintf(fp, "%e ", tstSamInputs[ss*nInputs+ii]);
               fprintf(fp, "%e %e\n", tstSamOutputs[ss], tstOutputs[ss]);
            }
            fprintf(fp, "];\n");
            fwritePlotCLF(fp);
            fprintf(fp, "m = %d;\n", nInputs);
            fprintf(fp, "X1 = A(:,1);\n");
            fprintf(fp, "X2 = A(:,2);\n");
            fprintf(fp, "Y  = A(:,m+1) - A(:,m+2);\n");
            fprintf(fp, "subplot(2,2,1)\n");
            fprintf(fp, "plot3(X1,X2,Y,'*','markerSize',13')\n");
            fwritePlotXLabel(fp, "Input 1");
            fwritePlotYLabel(fp, "Input 2");
            fwritePlotZLabel(fp, "Output");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Errors w.r.t. Input 1 and 2");
            fprintf(fp, "subplot(2,2,2)\n");
            fprintf(fp, "plot(X1,Y,'*','markerSize',13')\n");
            fwritePlotXLabel(fp, "Input 1");
            fwritePlotYLabel(fp, "Output");
            fwritePlotAxes(fp); 
            fwritePlotTitle(fp, "Output vs Input 1");
            fprintf(fp, "subplot(2,2,4)\n");
            fprintf(fp, "plot(X2,Y,'*','markerSize',13')\n");
            fwritePlotXLabel(fp, "Input 2");
            fwritePlotYLabel(fp, "Output");
            fwritePlotAxes(fp); 
            fwritePlotTitle(fp, "Output vs Input 2");
            fprintf(fp, "AA = [\n");
            for (ss = 0; ss < nSamples; ss++)
            {
               for (ii = 0; ii < nInputs; ii++)
                  fprintf(fp, "%e ", sampleInputs[ss*nInputs+ii]);
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "XX1 = AA(:,1);\n");
            fprintf(fp, "XX2 = AA(:,2);\n");
            fprintf(fp, "subplot(2,2,3)\n");
            fprintf(fp, "plot(XX1,XX2,'*','markerSize',13')\n");
            fwritePlotXLabel(fp, "Input 1");
            fwritePlotYLabel(fp, "Input 2");
            fwritePlotAxes(fp); 
            fwritePlotTitle(fp, "Input 1 vs 2 (original sample)");
            printf("PSUADE adaptive(1): error plots are in %s.\n",
                   cString);
            fclose(fp);
         }
         delete [] tstOutputs;
      }
      printEquals(0);

      refineLevel++;
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);
      if (errL2 < anaThreshold)
      {
         printf("PSUADE adaptive(1): threshold reached.\n");
         printf("                    unscaled rms = %en", errL2);
         printf("                    threshold    = %e\n", anaThreshold);
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         delete [] samInputs2;
         delete [] samOutputs2;
         delete [] samStates2;
         delete [] samStds2;
         delete faPtr;
         break;
      }
      if (psAnalysisInteractive_ == 1)
      {
         cString[0] = 'n';
         if (cString[0] == 'y')
         {
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            delete [] samInputs2;
            delete [] samOutputs2;
            delete [] samStates2;
            delete [] samStds2;
            delete faPtr;
            break;
         }
      }
      if (refineLevel > nRefinements)
      {
         printf("PSUADE adaptive(1): number of refinements %d reached.\n", 
                nRefinements);
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         delete [] samInputs2;
         delete [] samOutputs2;
         delete [] samStates2;
         delete [] samStds2;
         delete faPtr;
         break;
      }

      if (useRandomPts == 1)
      {
         iArray = new int[nSamples];
         generateRandomIvector(nSamples, iArray);
         for (ss = nSamples-1; ss >= nSamples-refineSize; ss--) 
            samStds2[iArray[ss]] = 1.0;
         for (ss = 0;  ss < nSamples-refineSize; ss++) 
            samStds2[iArray[ss]] = 0.0;
         delete [] iArray;
      }
      else
      {
         for (ss = 0; ss < nSamples; ss++) 
            samStds2[ss] = PABS(samStds2[ss+nSamples]);
      }
      sampler_->refine(refineRatio,randomize,refineThreshold,nSamples,samStds2);
      psuadeIO_->updateMethodSection(-1,-1,-1,nRefinements-1,-1);
      psuadeIO_->writePsuadeFile(NULL,0);

      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
      delete [] samInputs2;
      delete [] samOutputs2;
      delete [] samStates2;
      delete [] samStds2;
      delete samplerAux;
      delete faPtr;
   }

   delete funcIO;
   if (auxIO != NULL)
   {
      if (auxSamInputs  != NULL) delete [] auxSamInputs;
      if (auxSamOutputs != NULL) delete [] auxSamOutputs;
      delete auxIO;
   }
   for (ii = 0; ii < numMars; ii++)
   {
      delete [] marsDataX[ii];
      delete [] marsDataY[ii];
   }
   delete [] marsDataX;
   delete [] marsDataY;
   if (tstIO != NULL) delete tstIO;
   return 0;
}

// ************************************************************************
// run the special adaptive mode for response surface analysis (GMETIS)
// using MARS with bagging as response surface (bootstrap)
// This method uses arbitrary initial set of sample points.
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveErrBasedG()
{
   int    samMethod, refineLevel, nInputs, nOutputs, nRefinements;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, ss, nSamples, refineRatio, randomize, iSum, *inputPDFs;
   int    loopFlag, currNJobs, *sampleStates, nJobsDiff, *samStates2=NULL;
   int    parallelJobCount, status, maxState, jobsCompletedLast;
   int    rstype, nPtsPerDim=64, length, refineSize, marsMode=0;
   int    nSamples2, *iArray, auxNSamples = 0, auxNInputs, auxNOutputs;
   int    tstNSamples = 0, tstNInputs, tstNOutputs, ivar1, ivar2, numMars=100;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs;
   double refineThreshold=1.0, dtemp;
   double anaThreshold, *samInputs2=NULL, *samOutputs2=NULL;
   double *samStds2=NULL, totalSum=0.0, errMax=0.0, errAvg=0.0, errL2=0.0;
   double *auxSamInputs=NULL, *auxSamOutputs=NULL;
   double *tstSamInputs=NULL, *tstSamOutputs=NULL, *tstOutputs;
   char   systemCommand[100], cString[100], winput[500];
   char   *targv[3];
   FILE   *fp;
   pData  pPtr, pLowerB, pUpperB, pPDFs, pInpData, pOutData, pStates;
   pData  pAuxInputs, pAuxOutputs, pTstInputs, pTstOutputs;
   string sparam;
   FuncApprox *faPtr=NULL;
   Sampling   *samplerAux=NULL;
   PsuadeData *auxIO=NULL, *tstIO=NULL;
   FunctionInterface *funcIO=NULL;

   printAsterisks(0);
   printf("PSUADE adaptive(G): Use GMetis with MarsBagging.\n");
   printf("Adaptive sampling based on predicted errors from the\n");
   printf("MARS-with-bagging response surfaces (differs from RSMMB\n");
   printf("by allowing initial sample to be arbitrary instead of Metis).\n");
   printDashes(0);
   printf("To run this method, an initial sample can be provided.\n");
   printf("The sample should be somewhat space-filling.\n");
   printf("To use an initial sample, insert it in PSUADE_IO section,\n");
   printf("use ARSMG with desired number of refinements, and set\n");
   printf("refinement size.\n");
   printEquals(0);
   printf("To be able to exercise more control of this method,\n");
   printf("you can turn on interactive and/or expert modes in\n");
   printf("the input file.\n");
   printf("Turn on outputLevel (>0) to get error histogram.\n");
   printEquals(0);

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ss = 0; ss < nInputs; ss++) iSum += inputPDFs[ss];
   if (iSum)
   {
      printf("PSUADE adaptive(G) ERROR: does not currently allow \n");
      printf("       non-uniform probability distribution to be\n");
      printf("       defined in the INPUT SECTION.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   if (nOutputs > 1)
   {
      printf("PSUADE adaptive(G): nOutputs should be 1.\n");
      printf("       INFO: use 'write' in interactive model to select\n");
      printf("             1 output only and re-run this.\n");
      return 0;
   }
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;

   printf("You may test the quality of the response surface using a test\n");
   printf("sample (in psuadeData format).\n");
   sprintf(cString, "Use a test sample ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y')
   {
      sprintf(cString, "Enter the test sample file name : ");
      getString(cString, winput);
      tstIO = new PsuadeData();
      status = tstIO->readPsuadeFile(winput);
      tstIO->getParameter("method_nsamples", pPtr);
      tstNSamples = pPtr.intData_;
      tstIO->getParameter("input_ninputs", pPtr);
      tstNInputs = pPtr.intData_;
      tstIO->getParameter("output_noutputs", pPtr);
      tstNOutputs = pPtr.intData_;
      if (tstNInputs != nInputs)
      {
         printf("PSUADE adaptive(G) ERROR : test sample nInputs != %d\n",
                nInputs);
         return 0;
      }
      if (tstNOutputs > 1)
      {
         printf("PSUADE adaptive(G): test sample nOutputs != 1.\n");
         printf("       INFO: use 'write' in interactive model to select\n");
         printf("             1 output only and re-run this.\n");
         return 0;
      }
      tstIO->getParameter("input_sample", pTstInputs);
      tstSamInputs = pTstInputs.dbleArray_;
      tstIO->getParameter("output_sample", pTstOutputs);
      tstSamOutputs = pTstOutputs.dbleArray_;
   }

   printf("You may add to the base sample an auxiliary sample which\n");
   printf("covers the corners (for example, factorial or fractional\n");
   printf("factorial).\n");
   sprintf(cString, "Add an auxiliary sample ? (y or n) ");
   getString(cString, winput);
   if (winput[0] == 'y')
   {
      sprintf(cString, "Enter auxiliary sample file name : ");
      getString(cString, winput);
      ss = strlen(winput);
      winput[ss-1] = '\0';
      auxIO = new PsuadeData();
      status = auxIO->readPsuadeFile(winput);
      if (status == -1) exit(1);
      auxIO->getParameter("method_nsamples", pPtr);
      auxNSamples = pPtr.intData_;
      auxIO->getParameter("input_ninputs", pPtr);
      auxNInputs = pPtr.intData_;
      auxIO->getParameter("output_noutputs", pPtr);
      auxNOutputs = pPtr.intData_;
      if (auxNInputs != nInputs)
      {
         printf("PSUADE adaptive(G) ERROR : auxiliary nInputs != %d\n",
                nInputs);
         return 0;
      }
      if (auxNOutputs > 1)
      {
         printf("PSUADE adaptive(G): auxiliary sample nOutputs != 1.\n");
         printf("       INFO: use 'write' in interactive model to select\n");
         printf("             1 output only and re-run this.\n");
         return 0;
      }
      auxIO->getParameter("input_sample", pAuxInputs);
      samInputs2 = pAuxInputs.dbleArray_;
      length = (nSamples+auxNSamples+nRefinements*refineSize)*nInputs;
      auxSamInputs = new double[length];
      for (ss = 0; ss < auxNSamples*nInputs; ss++)
         auxSamInputs[ss] = samInputs2[ss];
      auxIO->getParameter("output_sample", pAuxOutputs);
      samOutputs2 = pAuxOutputs.dbleArray_;
      length = (nSamples+auxNSamples+nRefinements*refineSize)*nOutputs;
      auxSamOutputs = new double[length];
      for (ss = 0; ss < auxNSamples*nInputs; ss++)
         auxSamOutputs[ss] = samOutputs2[ss];
   }
   if (psRSExpertMode_ == 1)
   {
      sprintf(cString, "Number of MARS (default = 100, > 2, < 502) = ");
      numMars = getInt(3, 501, cString);
      sprintf(cString, "Use mean (0) or median (1) of MarSBag : ");
      marsMode = getInt(0, 1, cString);
   }

   psuadeIO_->getParameter("input_sample", pInpData);
   sampleInputs = pInpData.dbleArray_;
   psuadeIO_->getParameter("output_sample", pOutData);
   sampleOutputs = pOutData.dbleArray_;
   psuadeIO_->getParameter("output_states", pStates);
   sampleStates = pStates.intArray_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;
   psuadeIO_->getParameter("ana_rstype", pPtr);
   rstype = pPtr.intData_;
   if (rstype != PSUADE_RS_MARSB)
   {
      printf("PSUADE adaptive(G): RS type defaulted to MarsBag.\n");
      rstype = PSUADE_RS_MARSB; 
   }
   psuadeIO_->getParameter("ana_threshold", pPtr);
   anaThreshold = pPtr.dbleData_;
   refineRatio = 2;
   randomize = 1;

   if (sampler_ != NULL) SamplingDestroy(sampler_);
   samMethod = PSUADE_SAMP_GMETIS;
   sampler_ = (Sampling *) SamplingCreateFromID(samMethod);
   sampler_->setPrintLevel(outputLevel_);
   sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
   sampler_->setOutputParams(nOutputs);
   //sampler_->setSamplingParams(nSamples, 1, 1);
   sparam.clear();
   sprintf(cString, "reset");
   sparam.append(cString);
   sampler_->setParam(sparam);
   sparam.clear();
   sparam.append("setAdaptiveRefinementBasedOnErrors");
   sampler_->setParam(sparam);
   if (refineSize > 0)
   {
      sparam.clear();
      sprintf(cString, "setRefineSize %d", refineSize);
      sparam.append(cString);
      sampler_->setParam(sparam);
      printf("PSUADE adaptive(G): refineSize = %d\n",refineSize);
   }
   sampler_->initialize(1);
   sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                         sampleOutputs, sampleStates);
 
   funcIO = createFunctionInterface(psuadeIO_);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   loopFlag = 1;
   jobsCompleted = 0;
   refineLevel = 0;
   nSamples = -1;

   while (loopFlag)
   {
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 1)
      {
         printEquals(0);
         printf("PSUADE adaptive(G): current level    = %d (of %d)\n",
                refineLevel, nRefinements);
         printf("PSUADE adaptive(G): current nSamples = %d\n",nSamples);
      }
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                    sampleInputs,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                     sampleStates,NULL);
      psuadeIO_->writePsuadeFile(NULL,0);

      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss]==0) currNJobs++;
      if (psAnalysisInteractive_ == 1 && jobsCompleted < currNJobs)
      {
         printf("A sample has been created in the psuadeData file. \n");
         printf("At this point you can choose to run your model with\n");
         printf("this sample via psuade or by yourself (if the model\n");
         printf("is expensive to run, you want to choose the latter).\n");
         sprintf(systemCommand, "Run the model yourself (y or n) ? ");
         getString(systemCommand, cString);
         if (cString[0] == 'y')
         {
            printf("You have chosen to run the sample yourself.\n");
            printf("The following are steps to continue ARSM refinement:\n");
            printf("(1) Rename psuadeData to something else (e.g. psData). \n");
            printf("(2) Run the sample in this file and collect outputs.\n");
            printf("(3) Replace the outputs in psData (outputs only).\n");
            printf("(4) Finally, restart psuade with this file (psData).\n");
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            delete funcIO;
            return 0;
         }
      }

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptive(G): sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("PSUADE adaptive(G): jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      printf("Perform uniform refinement of the current sample.\n");
      samplerAux = (Sampling *) SamplingCreateFromID(samMethod);
      samplerAux->setInputBounds(nInputs, iLowerB, iUpperB);
      samplerAux->setOutputParams(nOutputs);
      samplerAux->setSamplingParams(nSamples, -1, -1);
      samplerAux->initialize(1);
      samplerAux->loadSamples(nSamples, nInputs, nOutputs, sampleInputs,
                              sampleOutputs, sampleStates);
      sparam.clear();
      sparam.append("setUniformRefinement");
      samplerAux->setParam(sparam);
      samplerAux->refine(refineRatio,randomize,0,0,NULL);
      nSamples2 = samplerAux->getNumSamples();
      samInputs2  = new double[nSamples2 * nInputs];
      samOutputs2 = new double[nSamples2 * nOutputs];
      samStates2  = new int[nSamples2];
      samStds2    = new double[nSamples2];
      samplerAux->getSamples(nSamples2, nInputs, nOutputs, samInputs2, 
                             samOutputs2, samStates2);

      printAsterisks(0);
      printf("PSUADE adaptive(G): response surface analysis.\n");
      printf("PSUADE adaptive(G): current level     = %d (of %d)\n",
             refineLevel, nRefinements);
      printf("                    training nSamples = %d\n",
             nSamples+auxNSamples);
      printf("                    test set nSamples = %d\n",
             nSamples2-nSamples);

      faPtr = genFA(rstype, nInputs, nOutputs, nSamples+auxNSamples);
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel_);
      if (psRSExpertMode_ == 0)
      {
         strcpy(cString, "mars_params");
         targv[0] = (char *) cString;
         ivar1 = (nSamples + auxNSamples);
         targv[1] = (char *) &ivar1;
         ivar2 = 2 * nInputs / 3 + 1;
         printf("Set degree of interaction = 3\n");
         ivar2 = 3;
         targv[2] = (char *) &ivar2;
         faPtr->setParams(3, targv);
         strcpy(cString, "num_mars");
         targv[0] = (char *) cString;
         targv[1] = (char *) &numMars;
         faPtr->setParams(2, targv);
         if (marsMode == 1)
         {
            strcpy(cString, "median");
            targv[0] = (char *) cString;
            faPtr->setParams(1, targv);
         }
      }
      length = -999;
      if (auxIO != NULL)
      {
         for (ss = 0; ss < nSamples*nInputs; ss++)
            auxSamInputs[auxNSamples*nInputs+ss] = sampleInputs[ss];
         for (ss = 0; ss < nSamples; ss++)
            auxSamOutputs[auxNSamples+ss] = sampleOutputs[ss];
         faPtr->genNDGridData(auxSamInputs,auxSamOutputs,&length,NULL,NULL);
      }
      else faPtr->genNDGridData(sampleInputs,sampleOutputs,&length,NULL,NULL);

      faPtr->evaluatePointFuzzy(nSamples2-nSamples, 
                                &samInputs2[nInputs*nSamples], 
                                &samOutputs2[nSamples], &samStds2[nSamples]);

      // extract the magnitude of each sample points
      totalSum = errMax = errAvg = errL2 =0.0;
      for (ss = nSamples; ss < nSamples2; ss++)
      {
         totalSum += PABS(samOutputs2[ss]);
         errMax = (samStds2[ss] > errMax) ? samStds2[ss] : errMax;
         errAvg += samStds2[ss];
         errL2  += pow(samStds2[ss], 2.0e0);
      }
      errL2 = sqrt(errL2 / (nSamples2 - nSamples));
      totalSum /= (nSamples2 - nSamples);
      errAvg  = errAvg / (nSamples2 - nSamples);
      printf("     response surface unscaled max error = %e\n",errMax);
      printf("     response surface   scaled max error = %e\n",errMax/totalSum);
      printf("     response surface unscaled rms error = %e\n",errL2);
      printf("     response surface   scaled rms error = %e\n",errL2/totalSum);
      printf("     response surface unscaled avg error = %e\n",errAvg);
      printf("     response surface   scaled avg error = %e\n",
             errAvg/totalSum);
      if (tstNSamples > 0)
      {
         printEquals(0);
         tstOutputs = new double[tstNSamples];
         faPtr->evaluatePoint(tstNSamples, tstSamInputs, tstOutputs);
         totalSum = errMax = errAvg = errL2 =0.0;
         for (ss = 0; ss < tstNSamples; ss++)
         {
            totalSum += PABS(tstOutputs[ss]);
            dtemp = tstOutputs[ss] - tstSamOutputs[ss];
            errAvg += dtemp;
            dtemp = PABS(dtemp);
            errL2  += dtemp * dtemp;
            errMax = (dtemp > errMax) ? dtemp : errMax;
         }
         errL2 = sqrt(errL2 / tstNSamples);
         totalSum /= (double) tstNSamples;
         errAvg  = errAvg / tstNSamples;
         printf("     test sample RS unscaled max error = %e\n",errMax);
         printf("     test sample RS   scaled max error = %e\n",errMax/totalSum);
         printf("     test sample RS unscaled rms error = %e\n",errL2);
         printf("     test sample RS   scaled rms error = %e\n",errL2/totalSum);
         printf("     test sample RS unscaled avg error = %e\n",errAvg);
         printf("     test sample RS   scaled avg error = %e\n",
                errAvg/totalSum);
         if (outputLevel_ > 0)
         {
            sprintf(cString, "rsag_err_hist_%d.m", refineLevel);
            fp = fopen(cString, "w");
            if (fp != NULL)
            {
               fprintf(fp, "%% first column: true outputs\n");
               fprintf(fp, "A = [\n");
               for (ss = 0; ss < tstNSamples; ss++)
                  fprintf(fp, "%e %e\n", tstSamOutputs[ss], tstOutputs[ss]);
               fprintf(fp, "];\n");
               fwritePlotCLF(fp);
               fprintf(fp, "hist(A(:,1)-A(:,2), 10)\n");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp, "Distribution of Prediction Errors");
               fwritePlotXLabel(fp, "Output Value");
               sprintf(winput, "Count (total = %d)", tstNSamples);
               fwritePlotYLabel(fp, winput);
               fprintf(fp, "E = hist(A(:,1)-A(:,2);\n");
               fprintf(fp, "mean = sum(E)/length(E)\n");
               fprintf(fp, "std  = sqrt(sum((E - m) .* (E - m)) / length(E))\n");
               fclose(fp);
               printf("PSUADE adaptive(G): error distribution plot is in %s.\n", 
                      cString);
            }
         }
         delete [] tstOutputs;
      }

      refineLevel++;
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);
      if (errL2 < anaThreshold)
      {
         printf("PSUADE adapive(G): threshold reached (based on unscaled L2).\n");
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         delete [] samInputs2;
         delete [] samOutputs2;
         delete [] samStates2;
         delete [] samStds2;
         break;
      }
      if (psAnalysisInteractive_ == 1)
      {
         sprintf(systemCommand, "Do you want to quit now (y or n) ? ");
         getString(systemCommand, cString);
         if (cString[0] == 'y')
         {
            delete [] sampleInputs;
            delete [] sampleOutputs;
            delete [] sampleStates;
            delete [] samInputs2;
            delete [] samOutputs2;
            delete [] samStates2;
            delete [] samStds2;
            break;
         }
      }
      if (refineLevel > nRefinements)
      {
         printf("PSUADE adapive(G): number of refinements %d reached.\n", 
                nRefinements);
         delete [] sampleInputs;
         delete [] sampleOutputs;
         delete [] sampleStates;
         delete [] samInputs2;
         delete [] samOutputs2;
         delete [] samStates2;
         delete [] samStds2;
         break;
      }

      iArray = new int[nSamples2-nSamples];
      for (ss = 0; ss < nSamples2-nSamples; ss++) iArray[ss] = ss;
      for (ss = 0; ss < nSamples2-nSamples; ss++) 
         samStds2[ss] = PABS(samStds2[ss+nSamples]);
      sortDbleList2a(nSamples2-nSamples, samStds2, iArray);
      for (ss=nSamples2-nSamples-1;ss >= nSamples2-nSamples-refineSize; ss--) 
         samStds2[iArray[ss]] = 1.0;
      for (ss = 0;  ss < nSamples2-nSamples-refineSize; ss++) 
         samStds2[iArray[ss]] = 0.0;
      sampler_->refine(refineRatio,randomize,refineThreshold,
                       nSamples2-nSamples,samStds2);
      psuadeIO_->updateMethodSection(-1,-1,-1,nRefinements-1,-1);
      psuadeIO_->writePsuadeFile(NULL,0);
      delete [] iArray;

      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
      delete [] samInputs2;
      delete [] samOutputs2;
      delete [] samStates2;
      delete [] samStds2;
      delete samplerAux;
      delete faPtr;
   }

   delete funcIO;
   if (auxIO != NULL)
   {
      if (auxSamInputs  != NULL) delete [] auxSamInputs;
      if (auxSamOutputs != NULL) delete [] auxSamOutputs;
      delete auxIO;
   }
   if (tstIO != NULL) delete tstIO;
   return 0;
}

// ************************************************************************
// run the special adaptive mode for risk analysis (with METIS)
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptivePRA()
{
   int    samMethod, initFlag, refineLevel, nInputs, nOutputs, nRefinements;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, nSamples, refineRatio, randomize, count, askFlag=0;
   int    ss, loopFlag, currNJobs, *sampleStates, nJobsDiff;
   int    parallelJobCount, status, maxState, jobsCompletedLast, curVol;
   int    totVol, refineType, numSuccess=0, refineSize;
   int    nUniform, iSum, *inputPDFs;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs, *tempY;
   double refineThreshold=1.0, anaThreshold, curVal=1.0, lastVal, relMax;
   double relMin, *sampErrors, ddata;
   char   systemCommand[100], cString[100];
   pData  pPtr, pLowerB, pUpperB, pPDFs;
   string sparam;
   FunctionInterface *funcIO=NULL;

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ss = 0; ss < nInputs; ss++) iSum += inputPDFs[ss];
   if (iSum)
   {
      printf("PSUADE adaptivePRA ERROR: does not currently allow \n");
      printf("       non-uniform probability distribution to be\n");
      printf("       defined in the INPUT SECTION.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   psuadeIO_->getParameter("method_refine_type", pPtr);
   refineType = pPtr.intData_;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samMethod = pPtr.intData_;
   if (samMethod != PSUADE_SAMP_METIS)
   {
      printf("PSUADE adaptivePRA: sampling defaulted to METIS.\n");
      samMethod = PSUADE_SAMP_METIS; 
      psuadeIO_->updateMethodSection(samMethod,-1,-1,-1,-1);
      initFlag = 1;
      if (sampler_ != NULL) SamplingDestroy(sampler_);
      sampler_ = NULL;
   }
   else initFlag = 0;

   psuadeIO_->getParameter("ana_threshold", pPtr);
   anaThreshold = pPtr.dbleData_;
   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   if (nOutputs > 1)
   {
      printf("PSUADE adaptivePRA: nOutputs should be 1.\n");
      return 0;
   }
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;

   if (initFlag == 1)
   {
      sampler_ = (Sampling *) SamplingCreateFromID(samMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, 1, 1);
      sampler_->initialize(0);
      nSamples = sampler_->getNumSamples();
      psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
   }
   sparam.clear();
   if (refineType == 0)
   {
      sparam.append("setUniformRefinement");
      refineSize = 100000;
   }
   else sparam.append("setAdaptiveRefinementBasedOnOutputs");
   sampler_->setParam(sparam);
   if (refineSize > 0)
   {
      sparam.clear();
      sprintf(cString, "setRefineSize %d", refineSize);
      sparam.append(cString);
      sampler_->setParam(sparam);
      printf("PSUADE adaptivePRA: refineSize = %d\n",refineSize);
   }
   refineRatio = 2;
   randomize = 1;

   funcIO = createFunctionInterface(psuadeIO_);

   printf("PSUADE adaptivePRA: initial nSamples = %d\n", nSamples);
   printf("PSUADE adaptivePRA: number of refinements = %d\n", nRefinements);
   printf("This function finds the first division between fail/no fail\n");
   printf("and focuses the adaptive sampling on the interface. Hence, it\n");
   printf("is important that the initial nSamples is sufficient large to\n");
   printf("uncover all fail/no fail interfaces. Alternatively, you can\n");
   printf("enforce the number of uniform refinements before adaptive\n");
   printf("sampling is applied. You have the opportunity to do it now.\n");
   sprintf(cString,"Number of initial uniform refinements : (0 - %d): ",
           nRefinements);
   nUniform = getInt(0, nRefinements, cString);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   relMax = relMin = 0.0;
   loopFlag = 1;
   jobsCompleted = 0;
   refineLevel = 0;
   nSamples = -1;

   while (loopFlag)
   {
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 1)
      {
         printAsterisks(0);
         printf("PSUADE adaptivePRA: current level    = %d (of %d)\n",
                refineLevel, nRefinements);
         printf("PSUADE adaptivePRA: current nSamples = %d\n",nSamples);
      }
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                    sampleInputs,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                     sampleStates,NULL);
      psuadeIO_->writePsuadeFile(NULL,0);

      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss]==0) currNJobs++;

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptivePRA: sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("\nPSUADE adaptivePRA: jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      refineLevel++;
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);

      tempY = NULL;
      for (ss = 0; ss < nSamples; ss++)
      {
         if (sampleOutputs[ss] != 0 && sampleOutputs[ss] != 1)
            break;
      }

      if (ss != nSamples)
      {
         printf("PSUADE adaptivePRA: non 0/1 outputs (transform to 0/1).\n");
         if (askFlag == 0)
         {
            while (relMin >= relMax)
            {
               sprintf(cString,"Upper bound for safe region (-999 if none): ");
               relMax = getDouble(cString);
               if (relMax == -999) relMax = 1.0e10;
               sprintf(cString,"Lower bound for safe region (-999 if none): ");
               relMin = getDouble(cString);
               if (relMin == -999) relMin = -1.0e10;
               if (relMin >= relMax) printf("INVALID bounds.\n");
            } 
            askFlag = 1;
         }
         tempY = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++)
         {
            if (sampleOutputs[ss] < relMin || sampleOutputs[ss] > relMax)
                 tempY[ss] = 0;
            else tempY[ss] = 1;
         }
         sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                               tempY, sampleStates);
      }

      sparam.clear();
      sparam.append("calVolumes");
      curVol = sampler_->setParam(sparam);
      sparam.clear();
      sparam.append("totalVolumes");
      totVol = sampler_->setParam(sparam);
      lastVal = curVal;
      curVal = 1.0 * curVol / totVol;
      if (refineLevel == 1) lastVal = curVal;
      printf("Current reliability = %8.5f%% (%d/%d*100)\n",
                100*curVal,curVol,totVol);
      if (refineLevel > 1 && curVal != 1.0 && curVal != 0.0)
         printf("Convergence check = %e (%e), trial %d (of 2)\n",
                PABS(2.0*(curVal-lastVal)/(curVal+lastVal)),anaThreshold,
                numSuccess+1);
      if (refineLevel > 1 && curVal != 1.0 && curVal != 0.0 &&
          PABS(2.0*(lastVal-curVal)/(lastVal+curVal)) < anaThreshold)
      {
         if (curVal != 1.0)
         {
            if (numSuccess <= 0) numSuccess++;
            else
            {
               printf("Convergence check: %e <? %e\n", 
                      PABS(2.0*(lastVal-curVal)/(lastVal+curVal)), anaThreshold);
               printf("Convergence achieved: reliability = %8.5f%%\n",
                      100*curVal);
               break;
            }
         }
      }
      else numSuccess = 0;

      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);
      delete [] tempY;

      if (refineLevel > nRefinements) break;

      if (refineLevel <= nUniform)
      {
         printf("PSUADE adaptivePRA INFO: user requested uniform refinement.\n");
         sparam.clear();
         sparam.append("setUniformRefinement");
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
      }
      else if (relMax != relMin && (curVal == 0.0 || curVal == 1.0))
      {
         printf("INFO: outputs != 0/1, reliability = 0/1 => error adaptive.\n");
         sampErrors = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++) 
         {
            ddata = sampleOutputs[ss];
            if      (ddata == relMin) ddata = PSUADE_UNDEFINED;
            else if (ddata == relMax) ddata = PSUADE_UNDEFINED;
            else if (PABS(ddata-relMax) < PABS(ddata-relMin))
                                      ddata = 1.0 / PABS(ddata - relMax);
            else                      ddata = 1.0 / PABS(ddata - relMin);
            sampErrors[ss] = ddata;
         }
         sparam.clear();
         sparam.append("setAdaptiveRefinementBasedOnErrors");
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,
                          nSamples,sampErrors);
         delete [] sampErrors;
      }
      else if (relMax != relMin && curVal != 0.0 && curVal != 1.0)
      {
         printf("INFO: outputs != 0/1, reliability != 0/1 => adaptive.\n");
         tempY = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++)
         {
            if (sampleOutputs[ss] < relMin || sampleOutputs[ss] > relMax)
                 tempY[ss] = 0;
            else tempY[ss] = 1;
         }
         sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                               tempY, sampleStates);
         sparam.clear();
         sparam.append("setAdaptiveRefinementBasedOnOutputs");
         sampler_->setParam(sparam);
         sparam.clear();
         sprintf(cString, "setRefineSize 100000");
         sparam.append(cString);
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
         sparam.clear();
         sprintf(cString, "setRefineSize %d", refineSize);
         sparam.append(cString);
         sampler_->setParam(sparam);
         delete [] tempY;
         tempY = sampleOutputs;
         delete [] sampleInputs;
         delete [] sampleStates;
         count = nSamples;
         nSamples = sampler_->getNumSamples();
         sampleInputs  = new double[nSamples * nInputs];
         sampleOutputs = new double[nSamples * nOutputs];
         sampleStates  = new int[nSamples];
         sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                              sampleOutputs, sampleStates);
         for (ss = 0; ss < count; ss++) sampleOutputs[ss] = tempY[ss];
         sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                               sampleOutputs, sampleStates);
         delete [] tempY;
      }
      else if (relMax == relMin && (curVal == 0.0 || curVal == 1.0))
      {
         printf("INFO: outputs = 0/1, reliability = 0 or 1 => uniform\n");
         printf("      refinement (refinement size not limited)\n");
         sparam.clear();
         sparam.append("setUniformRefinement");
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
      }
      else
      {
         printf("INFO: outputs = 0/1, reliability != 0 or 1 => adaptive\n");
         printf("      refinement (refinement size not limited)\n");
         sparam.clear();
         sparam.append("setAdaptiveRefinementBasedOnOutputs");
         sampler_->setParam(sparam);
         sparam.clear();
         sprintf(cString, "setRefineSize 100000");
         sparam.append(cString);
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,
                          nSamples,sampErrors);
      }

      psuadeIO_->updateMethodSection(-1,-1,-1,nRefinements-1,-1);
      psuadeIO_->writePsuadeFile(NULL,0);

      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
   }

   delete funcIO;
   return 0;
}

// ************************************************************************
// run the special adaptive mode for optimization (with METIS)
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveOpt()
{
   int    samMethod, initFlag, refineLevel, nInputs, nOutputs, nRefinements;
   int    maxParallelJobs, maxJobWaitTime, minJobWaitTime, launchInterval;
   int    saveFrequency, nSamples, refineRatio, randomize;
   int    ss, loopFlag, currNJobs, *sampleStates, nJobsDiff;
   int    parallelJobCount, status, maxState, jobsCompletedLast;
   int    ii, refineType, refineSize;
   int    nUniform, iSum, *inputPDFs, nMins=1, *sortList;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs, *tempY;
   double refineThreshold=1.0, *sampErrors, dmin;
   char   systemCommand[100], cString[100];
   pData  pPtr, pLowerB, pUpperB, pPDFs;
   string sparam;
   FunctionInterface *funcIO=NULL;

   psuadeIO_->getParameter("input_pdfs", pPDFs);
   inputPDFs = pPDFs.intArray_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   iSum = 0;
   for (ss = 0; ss < nInputs; ss++) iSum += inputPDFs[ss];
   if (iSum)
   {
      printf("PSUADE ERROR: adaptiveOpt does not currently allow \n");
      printf("       non-uniform probability distribution to be \n");
      printf("       defined in the INPUT SECTION.\n");
      printf("       Please fix it and then run again.\n");
      exit(1);
   }

   psuadeIO_->getParameter("method_refine_type", pPtr);
   refineType = pPtr.intData_;
   psuadeIO_->getParameter("method_sampling", pPtr);
   samMethod = pPtr.intData_;
   if (samMethod != PSUADE_SAMP_METIS)
   {
      printf("PSUADE adaptiveOpt: sampling defaulted to METIS.\n");
      samMethod = PSUADE_SAMP_METIS; 
      psuadeIO_->updateMethodSection(samMethod,-1,-1,-1,-1);
      initFlag = 1;
      if (sampler_ != NULL) SamplingDestroy(sampler_);
      sampler_ = NULL;
   }
   else initFlag = 0;

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   if (nOutputs > 1)
   {
      printf("PSUADE adaptiveOpt: nOutputs should be 1.\n");
      return 0;
   }
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefinements = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;

   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("app_minjobwaittime", pPtr);
   minJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_maxjobwaittime", pPtr);
   maxJobWaitTime = pPtr.intData_;
   psuadeIO_->getParameter("app_launchinterval", pPtr);
   launchInterval = pPtr.intData_;
   psuadeIO_->getParameter("app_savefrequency", pPtr);
   saveFrequency = pPtr.intData_;
   psuadeIO_->getParameter("ana_diagnostics", pPtr);
   outputLevel_ = pPtr.intData_;

   if (initFlag == 1)
   {
      sampler_ = (Sampling *) SamplingCreateFromID(samMethod); 
      sampler_->setPrintLevel(outputLevel_);
      sampler_->setInputBounds(nInputs, iLowerB, iUpperB);
      sampler_->setOutputParams(nOutputs);
      sampler_->setSamplingParams(nSamples, 1, 1);
      sampler_->initialize(0);
      nSamples = sampler_->getNumSamples();
      psuadeIO_->updateMethodSection(-1,nSamples,-1,-1,-1);
   }
   sparam.clear();
   if (refineType == 0)
   {
      sparam.append("setUniformRefinement");
      refineSize = 100000;
   }
   else sparam.append("setAdaptiveRefinementBasedOnOutputs");
   sampler_->setParam(sparam);
   if (refineSize > 0)
   {
      sparam.clear();
      sprintf(cString, "setRefineSize %d", refineSize);
      sparam.append(cString);
      sampler_->setParam(sparam);
      printf("PSUADE adaptiveOpt: refineSize = %d\n",refineSize);
   }
   refineRatio = 2;
   randomize = 1;

   funcIO = createFunctionInterface(psuadeIO_);

   printf("PSUADE adaptiveOpt: initial nSamples = %d\n", nSamples);
   printf("PSUADE adaptiveOpt: number of refinements = %d\n", nRefinements);
   printf("This function focuses the adaptive sampling on the regions\n");
   printf("where most of the small objective function values are located.\n");
   printf("It is important that the initial nSamples is sufficient large\n");
   printf("to uncover all potential troughs. Alternatively, you can\n");
   printf("enforce the number of uniform refinements before adaptive\n");
   printf("sampling is applied. You have the opportunity to do it now.\n");
   sprintf(cString,"Number of initial uniform refinements : (0 - %d): ",
           nRefinements);
   nUniform = getInt(0, nRefinements, cString);
   sprintf(cString,"Number of minimum to track : (1 - 10): ");
   nMins = getInt(1, 10, cString);

   if (maxParallelJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   loopFlag = 1;
   jobsCompleted = 0;
   refineLevel = 0;
   nSamples = -1;

   while (loopFlag)
   {
      nSamples = sampler_->getNumSamples();
      if (outputLevel_ > 1)
      {
         printAsterisks(0);
         printf("PSUADE adaptiveOpt: current level    = %d (of %d)\n",
                refineLevel, nRefinements);
         printf("PSUADE adaptiveOpt: current nSamples = %d\n",nSamples);
      }
      sampleInputs  = new double[nSamples * nInputs];
      sampleOutputs = new double[nSamples * nOutputs];
      sampleStates  = new int[nSamples];
      sampler_->getSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                           sampleOutputs, sampleStates);
      psuadeIO_->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                    sampleInputs,NULL); 
      psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                     sampleStates,NULL);
      psuadeIO_->writePsuadeFile(NULL,0);

      currNJobs = jobsCompleted;
      for (ss = 0; ss < nSamples; ss++) if (sampleStates[ss]==0) currNJobs++;

      parallelJobCount = 0;
      while (jobsCompleted < currNJobs)
      {
         jobsCompletedLast = jobsCompleted;
         for (ss = 0; ss < nSamples; ss++)
         {
#ifdef HAVE_PYTHON
	    PyObject* temp;
	    if (update_gui != NULL) {
	       temp = PyObject_CallObject(update_gui, NULL);
	       if (temp != NULL) Py_DECREF(temp);
	    }
#endif
	    if (psuade_stop == 1)
	    {
	       psuadeIO_->writePsuadeFile(NULL,0);
               delete [] sampleStates;
               delete [] sampleOutputs;
               delete [] sampleInputs;
               delete funcIO;
	       throw Psuade_Stop_Exception();
	    }

            if ((sampleStates[ss] == 0) && 
                (parallelJobCount < maxParallelJobs))
            {
               status = funcIO->evaluate(ss,nInputs,
                            &sampleInputs[ss*nInputs], nOutputs, 
                            &sampleOutputs[ss*nOutputs],0);   

               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                             &(sampleOutputs[ss*nOutputs]),
                             &(sampleStates[ss]));
                  jobsCompleted++;
                  nJobsDiff = jobsCompleted - jobsCompletedLast;
                  if ((nJobsDiff % saveFrequency) == (saveFrequency-1))
                  {
                     psuadeIO_->updateOutputSection(nSamples,nOutputs,
                                      sampleOutputs,sampleStates,NULL);
                     psuadeIO_->writePsuadeFile(NULL,0);
                  }
                  if (outputLevel_ > 0)
                  {
                     if ((jobsCompleted % 11) == 1) printf(".");
                     fflush(stdout);
                  }
               }
               else
               {
                  sampleStates[ss] = status;
                  parallelJobCount++;
               }
            }
            else if (sampleStates[ss] >= 2)
            {
               status = funcIO->evaluate(ss,nInputs,
                              &sampleInputs[ss*nInputs], nOutputs, 
                              &sampleOutputs[ss*nOutputs],2);   
               if (status == 0) 
               {
                  sampler_->storeResult(ss, nOutputs,
                                   &(sampleOutputs[ss*nOutputs]),
                                   &(sampleStates[ss]));
                  jobsCompleted++;
                  parallelJobCount--;
               }
               else sampleStates[ss]++;

               if (minJobWaitTime > 0 &&
                   (sampleStates[ss]-2)*minJobWaitTime>maxJobWaitTime)
               {
                  sampleStates[ss] = 0;
                  parallelJobCount--;
                  ss--; /* roll back to the sample to be restarted */
                  if (outputLevel_ > 0) 
                     printf("PSUADE adaptiveOpt: sample %6d to be restarted.\n",
                            ss+1);
               }
            }
         }

         psuadeIO_->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                        sampleStates,NULL);
         psuadeIO_->writePsuadeFile(NULL,0);

         maxState = 0;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > maxState) maxState = sampleStates[ss];
         if ((maxState > 100) && (minJobWaitTime <= 0)) minJobWaitTime *= 2;
         for (ss = 0; ss < nSamples; ss++)
            if (sampleStates[ss] > 10) sampleStates[ss] /= 2;

         if ((jobsCompleted < currNJobs) && (minJobWaitTime > 0))
         {
            sprintf(systemCommand, "sleep %d", minJobWaitTime);
            system(systemCommand);
         }
         if (outputLevel_ > 0) 
            printf("\nPSUADE adaptiveOpt: jobs completed = %d(of %d)\n",
                   jobsCompleted, currNJobs);
      }

      refineLevel++;
      sampler_->loadSamples(nSamples, nInputs, nOutputs, sampleInputs, 
                            sampleOutputs, sampleStates);

      tempY = new double[nSamples];
      sortList = new int[nSamples];
      for (ss = 0; ss < nSamples; ss++) tempY[ss] = sampleOutputs[ss];
      for (ss = 0; ss < nSamples; ss++) sortList[ss] = ss;
      sortDbleList2a(nSamples, tempY, sortList);
      for (ss = 0; ss < nMins; ss++)
         printf("PSUADE adaptiveOpt: min %2d = %e\n", ss+1, tempY[ss]);
      printf("PSUADE adaptiveOpt: min at \n");
      for (ii = 0; ii < nInputs; ii++)
         printf("   Input %2d = %e\n", ii+1, 
                sampleInputs[sortList[0]*nInputs+ii]);

      delete [] tempY;
      delete [] sortList;

      if (refineLevel > nRefinements) break;

      if (refineLevel <= nUniform)
      {
         printf("INFO: user requested uniform refinement.\n");
         sparam.clear();
         sparam.append("setUniformRefinement");
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,0,NULL);
      }
      else
      {
         sampErrors = new double[nSamples];
         dmin = PSUADE_UNDEFINED;
         for (ss = 0; ss < nSamples; ss++) 
            if (sampleOutputs[ss] < dmin) dmin = sampleOutputs[ss];
         for (ss = 0; ss < nSamples; ss++)
         {
            sampErrors[ss] = sampleOutputs[ss] - dmin;
            if (sampleOutputs[ss] == 0.0) sampErrors[ss] = PSUADE_UNDEFINED;
            else                          sampErrors[ss] = 1.0 / sampErrors[ss];
         }
         sparam.clear();
         sparam.append("setAdaptiveRefinementBasedOnErrors");
         sampler_->setParam(sparam);
         sampler_->refine(refineRatio,randomize,refineThreshold,
                          nSamples,sampErrors);
         delete [] sampErrors;
      }

      psuadeIO_->updateMethodSection(-1,-1,-1,nRefinements-1,-1);
      psuadeIO_->writePsuadeFile(NULL,0);

      delete [] sampleInputs;
      delete [] sampleOutputs;
      delete [] sampleStates;
   }

   delete funcIO;
   delete [] sampleStates;
   delete [] sampleInputs;
   delete [] sampleOutputs;
   return 0;
}

// ************************************************************************
// display response surface using Pgplot
// ------------------------------------------------------------------------
void PsuadeBase::pgPlotResponseSurface()
{
   int        loopFlag, iplot1, iplot2, iplot3, numIncrements, jplot;
   int        ind1, ind2, ind3, faLeng, minIndex, ii, ss, graphics;
   int        nPtsPerDim=64, threshFlag, nInputs, nOutputs, nSamples;
   long       nlong;
   double     *inputSettings, *faYIn, spacing, *faXOut, *faYOut;
   double     Ymin, Ymax, **faXOut2, threshL, threshU;
   double     *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs;
   char       pString[100], sInput[500];
   FuncApprox *faPtr;
   pData      pPtr, pLowerB, pUpperB, pInpData, pOutData; 

   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;
   psuadeIO_->getParameter("input_sample", pInpData);
   sampleInputs = pInpData.dbleArray_;
   psuadeIO_->getParameter("output_sample", pOutData);
   sampleOutputs = pOutData.dbleArray_;
   psuadeIO_->getParameter("ana_graphicsflag", pPtr);
   graphics = pPtr.intData_;

   if ((graphics & 1) != 0 && (nInputs > 1))
   {
      faPtr = genFAInteractive(psuadeIO_, 0);
      if (faPtr == NULL)
      {
         printf("INFO: cannot create function approximator.\n");
         return;
      }
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel_);
      loopFlag = 1;
      inputSettings = new double[nInputs];
      while (loopFlag == 1)
      {
         if (nInputs == 2)
         {
            iplot1 = 0;
            iplot2 = 1;
            numIncrements = 1;
         }
         else
         {
            sprintf(pString, "Enter the first axis number (1 - %d) ",nInputs);
            iplot1 = getInt(1, nInputs, pString);
            iplot1--;
            iplot2 = iplot1;
            sprintf(pString, "Enter the second axis number (1 - %d), not %d : ",
                    nInputs, iplot1+1);
            while (iplot2 < 0 || iplot2 >= nInputs || iplot1 == iplot2)
            {
               iplot2 = getInt(1, nInputs, pString);
               iplot2--;
               if ((iplot2 < 0) || (iplot2 >= nInputs) || (iplot1 == iplot2))
                  printf("ERROR : wrong input number %d.\n", iplot2+1);
            }
            if (nInputs == 3)
            {
               iplot3 = 0;
               if ((iplot3 == iplot1) || (iplot3 == iplot2)) iplot3 = 1;
               if ((iplot3 == iplot1) || (iplot3 == iplot2)) iplot3 = 2;
            }
            else
            {
               iplot3 = iplot2;
               sprintf(pString, "Enter the time axis number (1 - %d), not %d %d : ",
                       nInputs, iplot1+1, iplot2+1);
               while (iplot3 < 0 || iplot3 >= nInputs || 
                      iplot3 == iplot1 || iplot3 == iplot2)
               {
                  iplot3 = getInt(1, nInputs, pString);
                  iplot3--;
                  if ((iplot3 < 0) || (iplot3 >= nInputs) || 
                      (iplot3 == iplot1) || (iplot3 == iplot2))
                     printf("ERROR : wrong input number %d.\n", iplot3+1);
               }
            }
            for (ii = 0; ii < nInputs; ii++)
            {
               if ((ii != iplot1) && (ii != iplot2) && (ii != iplot3))
               {
                  inputSettings[ii] = iLowerB[ii] - 1.0;
                  sprintf(pString,"Enter nominal value of input %d : ",ii+1);
                  while (inputSettings[ii] < iLowerB[ii] ||
                         inputSettings[ii] > iUpperB[ii])
                     inputSettings[ii] = getDouble(pString);
               }
            }
            numIncrements = 10;
         }
         jplot = 0;
         if (nOutputs > 1)
         {
            sprintf(pString, "Enter the output number for z axis (1 - %d) ",
                    nOutputs);
            jplot = getInt(1, nOutputs, pString);
            jplot--;
         }
         if (nInputs > 2) 
            printf("\nPSUADE: plot input %d in 10 increment in time.\n", 
                   iplot3); 
         faYIn = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++) 
            faYIn[ss] = sampleOutputs[ss*nOutputs+jplot];

         Plotbegin(iLowerB[iplot1],iUpperB[iplot1],iLowerB[iplot2],
                   iUpperB[iplot2]);
         if (nInputs > 2) spacing = (iUpperB[iplot3] - iLowerB[iplot3]) * 0.1;
         for (ind1 = 0; ind1 < numIncrements; ind1++)
         {
            if (nInputs > 2)
            {
               inputSettings[iplot3] = iLowerB[iplot3] + ind1 * spacing;
               printf("Displaying graph %d out of 10.\n", ind1+1);
            }
            faPtr->gen2DGridData(sampleInputs,faYIn, iplot1, iplot2, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
            Ymin = Ymax = faYOut[0];
            minIndex = 0;
            for (ind2 = 1; ind2 < faLeng; ind2++)
            {
               if (faYOut[ind2] > Ymax) Ymax = faYOut[ind2];
               if (faYOut[ind2] < Ymin) 
               {
                  Ymin = faYOut[ind2];
                  minIndex = ind2;
               }
            }
            printf("PSUADE: X at YMIN = %16.8e %16.8e\n", 
                   faXOut[2*minIndex], faXOut[2*minIndex+1]);
            printf("PSUADE: YMIN, YMAX = %16.8e %16.8e\n", Ymin, Ymax);
            if ((Ymax - Ymin) < 1.0E-6)
            {
               printf("The surface is almost flat ===> \n");
               printf("   lift one corner to enable proper plotting.\n");
               faYOut[faLeng-1] = 1.001 * faYOut[faLeng-1];
            }
            threshFlag = 0;
            sprintf(pString,"Set lower threshold for plotting ? (y or n) ");
            getString(pString, sInput);
            if (sInput[0] == 'y')
            {
               threshFlag |= 1;
               sprintf(pString, "Enter lower threshold : ");
               threshL = getDouble(pString);
            }
            sprintf(pString,"Set upper threshold for plotting ? (y or n) ");
            getString(pString, sInput);
            if (sInput[0] == 'y')
            {
               threshFlag |= 2;
               sprintf(pString,"Enter upper threshold : ");
               threshU = getDouble(pString);
            }
            if ((threshFlag & 1) != 0)
            {
               for (ind2 = 0; ind2 < faLeng; ind2++)
                  if (faYOut[ind2] < threshL) 
                     faYOut[ind2] = Ymin - 0.2 * PABS(Ymin);
            }
            if ((threshFlag & 2) != 0)
            {
               for (ind2 = 0; ind2 < faLeng; ind2++)
                  if (faYOut[ind2] > threshU) 
                        faYOut[ind2] = Ymin - 0.2 * PABS(Ymin);
            }
            faXOut2 = new double*[2];
            for (ind2 = 0; ind2 < 2; ind2++)
               faXOut2[ind2] = new double[faLeng];
            ind3 = 0;
            for (ind2 = 0; ind2 < faLeng; ind2++)
            {
               faXOut2[0][ind2] = faXOut[ind3++];
               faXOut2[1][ind2] = faXOut[ind3++];
            }
            nlong = faLeng;
            Plot3d(nlong, faXOut2[0], faXOut2[1], faYOut);
            delete [] faXOut;
            delete [] faYOut;
            delete [] faXOut2[0];
            delete [] faXOut2[1];
            delete [] faXOut2;
            sprintf(pString,"Enter 1 to continue, otherwise to terminate.\n");
            if (ind1 < (numIncrements-1)) 
            {
               getString(pString, sInput);
               if (sInput[0] != '1') break;
            }
         }
         delete [] faYIn;
         Plotend();
         sprintf(pString,"\nEnter 1 to continue, otherwise to terminate.\n");
         getString(pString, sInput);
         if (strcmp(sInput, "1")) loopFlag = 1;
         else                     loopFlag = 0;
      }
      delete [] inputSettings;
      delete faPtr;
   }
}

// ************************************************************************
// deallocate variables 
// ------------------------------------------------------------------------
void PsuadeBase::cleanUp()
{
   if (sampler_  != NULL) SamplingDestroy(sampler_); 
   if (psuadeIO_ != NULL) delete psuadeIO_;
   if (optimalXData_  != NULL) delete optimalXData_; 
   if (optimalYData_  != NULL) delete optimalYData_; 
   sampler_ = NULL;
   psuadeIO_ = NULL;
   optimalXData_ = NULL;
   optimalYData_ = NULL;
   optimalCount_ = 0;
}

// ************************************************************************
// fill in samples (for duplicates)
// ------------------------------------------------------------------------
int PsuadeBase::fillInSamples(int nSamples, int nOutputs, 
                              double *sampleOutputs, int *sampleStates)
{
   int sampleID, index, ii;

   for (sampleID = 0; sampleID < nSamples; sampleID++)
   {
      if (sampleStates[sampleID] < 0)
      {
         index = -(sampleStates[sampleID] + 1);
         if (sampleStates[index] == 1)
         {
            for (ii = 0; ii < nOutputs; ii++)
               sampleOutputs[sampleID*nOutputs+ii] =  
                  sampleOutputs[index*nOutputs+ii];  
            sampleStates[sampleID] = 1;  
         }
      }
   } 
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
PsuadeBase& PsuadeBase::operator=(const PsuadeBase &)
{
   printf("PsuadeBase operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


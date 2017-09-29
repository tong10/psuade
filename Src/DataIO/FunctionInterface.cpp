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
// Functions for FunctionInterface 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
//#include <stream.h>
#include "dtype.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "FuncApprox.h"
#include "pData.h"
#include "PDFBase.h"
#include "PsuadeData.h"


// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
FunctionInterface::FunctionInterface()
{ 
   nInputs_ = 0;
   nOutputs_ = 0;
   inputNames_ = NULL;
   outputNames_ = NULL;
   strcpy(appDriver_, "true");
   strcpy(optDriver_, "true");
   strcpy(auxOptDriver_, "true");
   strcpy(appInputTemplate_, "psuadeApps_ct.in");
   strcpy(appOutputTemplate_, "psuadeApps_ct.out");
   executionMode_ = 0;
   launchInterval_ = 0;
   appOptFlag_ = 0;
   useRSModel_ = 0;
   rsPtrs_ = NULL;
   rsIndices_ = NULL;
   rsValues_ = NULL;
   rsRanFlag_ = 0;
   printLevel_ = 0;
   nInps_ = 0;
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
FunctionInterface::FunctionInterface(const FunctionInterface & fi)
{
   nInputs_ = fi.nInputs_;
   nOutputs_ = fi.nOutputs_;
   strcpy(appDriver_, fi.appDriver_);
   strcpy(optDriver_, fi.optDriver_);
   strcpy(auxOptDriver_, fi.auxOptDriver_);
   strcpy(appInputTemplate_, fi.appInputTemplate_);
   strcpy(appOutputTemplate_, fi.appOutputTemplate_);
   rsRanFlag_ = fi.rsRanFlag_;
   printLevel_ = fi.printLevel_;
   nInps_ = fi.nInps_;

   //Now for the dynamic memory

   // inputNames_
   inputNames_ = new char*[nInputs_];
   for(int i = 0; i < nInputs_; i++)
   {
      inputNames_[i] = new char[80];
      strcpy(inputNames_[i], fi.inputNames_[i]);
   }
  
   //  outputNames_
   outputNames_ = new char *[nOutputs_];
   for(int i = 0; i < nOutputs_; i++)
   {
      outputNames_[i] = new char[80];
      strcpy(outputNames_[i], fi.outputNames_[i]);
   }

   // rsPtrs_
   rsPtrs_ = new FuncApprox*[nOutputs_];
   // FuncApprox is a class whose copy constructor has been created and the
   // operator= has been defined so we can just do a copy
   for(int i = 0; i < nOutputs_; i++) rsPtrs_[i] = fi.rsPtrs_[i];

   //rsIndices_ and rsValues_
   rsIndices_ = new int[nInps_];
   rsValues_ = new double[nInps_];
   for(int i = 0; i < nInps_; i++)
   {
      rsIndices_[i] = fi.rsIndices_[i];
      rsValues_[i] = fi.rsValues_[i];
   }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FunctionInterface::~FunctionInterface()
{ 
   if (inputNames_ != NULL)
   {
      for (int i = 0; i < nInputs_; i++)
         if (inputNames_[i] != NULL) delete [] inputNames_[i];
      delete [] inputNames_;
   }
   if (outputNames_ != NULL)
   {
      for (int j = 0; j < nOutputs_; j++)
         if (outputNames_[j] != NULL) delete [] outputNames_[j];
      delete [] outputNames_;
   }
   if (rsPtrs_ != NULL)
   {
      for (int k = 0; k < nOutputs_; k++)
         if (rsPtrs_[k] != NULL) delete rsPtrs_[k];
      delete [] rsPtrs_;
   }
   if (rsIndices_ != NULL)
   {
      delete [] rsIndices_;
      rsIndices_ = NULL;
   }
   if (rsValues_ != NULL)
   {
      delete [] rsValues_;
      rsValues_ = NULL;
   }
}

// ************************************************************************
// load input data 
// ------------------------------------------------------------------------
int FunctionInterface::loadInputData(int nInputs, char **names)
{ 
   if (nInputs <= 0)
   {
      printf("FunctionInterface LoadInput ERROR: nInputs <= 0\n");
      return 1; 
   }
   nInputs_ = nInputs;
   inputNames_ = new char*[nInputs_];
   for (int i = 0; i < nInputs; i++)
   {
      inputNames_[i] = new char[80];
      if (names != NULL && names[i] != NULL)
         strcpy(inputNames_[i], names[i]);
      else sprintf(inputNames_[i], "X%d", i);
   }
   return 0;
}

// ************************************************************************
// load output data 
// ------------------------------------------------------------------------
int FunctionInterface::loadOutputData(int nOutputs, char **names)
{ 
   if (nOutputs <= 0)
   {
      printf("FunctionInterface loadOutput ERROR: nOutputs <= 0\n");
      return 1; 
   }
   nOutputs_ = nOutputs;
   outputNames_ = new char*[nOutputs_];
   for (int i = 0; i < nOutputs_; i++)
   {
      outputNames_[i] = new char[80];
      if (names != NULL && names[i] != NULL)
         strcpy(outputNames_[i], names[i]);
      else sprintf(outputNames_[i], "Y%d", i);
   }
   return 0;
}

// ************************************************************************
// load the Function information 
// ------------------------------------------------------------------------
int FunctionInterface::loadFunctionData(int length, char **names)
{ 
   int  ii, kk, nInps, cnt, status;
   char inString[200], fname[200];
   FILE *fp;
   PsuadeData *psIO;
   pData pPtr;

   if (length != 5)
   {
      printf("FunctionInterface loadFunction ERROR: length != 5\n");
      return 1; 
   }

   if (rsPtrs_ != NULL)
   {
      for (ii = 0; ii < nOutputs_; ii++)
         if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
      delete [] rsPtrs_;
      rsPtrs_ = NULL;
   }
   if (rsIndices_ != NULL)
   {
      delete [] rsIndices_;
      rsIndices_ = NULL;
   }
   if (rsValues_ != NULL)
   {
      delete [] rsValues_;
      rsValues_ = NULL;
   }

   strcpy(appDriver_, names[0]);
   if (strcmp(names[1], "NONE")) strcpy(appInputTemplate_, names[1]);
   if (strcmp(names[2], "NONE")) strcpy(appOutputTemplate_, names[2]);
   strcpy(optDriver_, names[3]);
   strcpy(auxOptDriver_, names[4]);
   if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
        strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 2;
      return 0; 
   }
   if ( strcmp(appDriver_,"PSUADE_LOCAL") &&
       !strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 3;
      return 0; 
   }
   if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
       !strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 4;
      return 0; 
   }

   useRSModel_ = 0;
   if      (appOptFlag_ == 0) strcpy(fname, appDriver_);
   else if (appOptFlag_ == 1) strcpy(fname, optDriver_);
   if      (appOptFlag_ == 2) strcpy(fname, auxOptDriver_);
   fp = fopen(fname, "r");
   if (fp != NULL)
   {
      fscanf(fp, "%10c", inString);
      if (!strncmp(inString, "PSUADE_IO",9)) useRSModel_ = 1;
      fclose(fp);
      if (useRSModel_ == 1)
      {
         psIO = new PsuadeData();
         psIO->setOutputLevel(0);
         status = psIO->readPsuadeFile(fname);
         if (status < 0) exit(1);
         if (status > 0)
         {
            printf("ERROR : cannot read file %d in PSUADE format\n",fname);
            exit(1);
         }
         psIO->getParameter("input_ninputs", pPtr);
         kk = pPtr.intData_;
         if (kk < nInputs_ || kk > nInputs_)
         {
            printf("WARNING: nInputs in driver not the same as nInputs in input file.\n");
            printf("** In order for this to work, you need to specify\n");
            printf("** a RS index file inside the driver (data) file\n");
            printf("(in ANALYSIS section, analyzer rs_index_file = <file>).\n");
            printf("Now, Checking for rs index file (1).\n");
         }
         psIO->getParameter("ana_rsindexfile", pPtr);
         if (kk != nInputs_ && !strcmp(pPtr.strArray_[0], "NONE"))
         {
            printf("ERROR: RS index file not specified.\n");
            printf("       nInputs in original .in file = %d\n", nInputs_);
            printf("       nInputs in driver data  file = %d\n", kk);
            exit(1);
         }
         if (strcmp(pPtr.strArray_[0], "NONE"))
         {
            fp = fopen(pPtr.strArray_[0], "r");
            if (fp == NULL)
            {
               printf("ERROR: rs index file %s not found.\n",pPtr.strArray_[0]);
               exit(1);
            }
            else
            {
               fscanf(fp,"%d", &nInps);
               if (kk != nInps)
               {
                  printf("ERROR: invalid nInputs in rs index file.\n");
                  printf("  data format should be (num: input number from calling file): \n");
                  printf("  line 1: nInputs in RS driver data file\n");
                  printf("  line 2: 1 <num> <default value if num == 0>\n");
                  printf("  line 3: 2 <num> <0 if num != 0>\n");
                  printf("  line 4: 3 <num> <default value if num == 0>\n");
                  printf("  line 5: 4 <num> <0 if num != 0>\n");
                  printf("  ...\n");
                  printf("  line nInputs+1: nInputs <num> <0 if num != 0>\n");
                  exit(1);
               }
               rsIndices_ = new int[nInps];
               rsValues_ = new double[nInps];
               cnt = 0;
               for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
               for (ii = 0; ii < nInps; ii++)
               {
                  fscanf(fp, "%d", &kk);
                  if (kk != ii+1)
                  {
                     printf("ERROR: first index in rsIndexFile = %d (must be %d]).\n",
                            rsIndices_[ii], ii+1);
                     printf("  data format should be (num: from calling file): \n");
                     printf("  line 1: nInputs in RS driver data file\n");
                     printf("  line 2: 1 <num> <default value if num == 0>\n");
                     printf("  line 3: 2 <num> <0 if num != 0>\n");
                     printf("  line 4: 3 <num> <default value if num == 0>\n");
                     printf("  line 5: 4 <num> <0 if num != 0>\n");
                     printf("  ...\n");
                     printf("  line nInputs+1: nInputs <num> <0 if num != 0>\n");
                     exit(1);
                  } 
                  fscanf(fp, "%d", &rsIndices_[ii]);
                  if (rsIndices_[ii] == 0)
                     printf("INFO: input %3d inactive\n",ii+1);
                  if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInps)
                  {
                     printf("INFO: input %3d = %d not valid\n",ii+1,rsIndices_[ii]);
                     exit(1);
                  }
                  rsIndices_[ii]--;
                  fscanf(fp, "%lg", &rsValues_[ii]);
                  if (rsIndices_[ii] >= -1) cnt++;
                  printf("RS_index %5d = %5d %16.8e\n",ii+1,rsIndices_[ii]+1,
                         rsValues_[ii]);
               }
               if (cnt != nInps)
               {
                  printf("ERROR: invalid number of active inputs.\n");
                  printf("  data format should be (num: from calling file): \n");
                  printf("  line 1: nInputs in RS driver data file\n");
                  printf("  line 2: 1 <num> <default value if num == 0>\n");
                  printf("  line 3: 2 <num> <0 if num != 0>\n");
                  printf("  line 4: 3 <num> <default value if num == 0>\n");
                  printf("  line 5: 4 <num> <0 if num != 0>\n");
                  printf("  ...\n");
                  printf("  line nInputs+1: nInputs <num> <0 if num != 0>\n");
                  exit(1);
               }
               fclose(fp);
            }
         }
         delete psIO;
         rsPtrs_ = new FuncApprox*[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            printf("Creating response surface for output %d\n", ii+1);
            rsPtrs_[ii] = (FuncApprox *) genFAFromFile(fname,ii);
            if (rsPtrs_[ii] == NULL)
            {
               printf("FunctionInterface ERROR: no RS model given.\n");
               useRSModel_ = 0;
               exit(1);
            }
         }
      }
   }
   return 0;
}

// ************************************************************************
// set execution mode to be asynchronous
// ------------------------------------------------------------------------
int FunctionInterface::setAsynchronousMode()
{
   executionMode_ = 1;
   return 0;
}

// ************************************************************************
// set execution mode to be synchronous
// ------------------------------------------------------------------------
int FunctionInterface::setSynchronousMode()
{
   executionMode_ = 0;
   return 0;
}

// ************************************************************************
// set stochastic response surface mode (1:minus 3 sigma, 2: plus 3 sigma)
// ------------------------------------------------------------------------
int FunctionInterface::setStochasticRSMode(int indata)
{
   rsRanFlag_ = indata;
   return 0;
}

// ************************************************************************
// set time interval between asynchronous jobs
// ------------------------------------------------------------------------
int FunctionInterface::setLaunchInterval(int interval)
{
   launchInterval_ = interval;
   if (launchInterval_ < 0) launchInterval_ = 0;
   if (launchInterval_ > 100000) launchInterval_ = 0;
   return 0;
}

// ************************************************************************
// set print level
// ------------------------------------------------------------------------
int FunctionInterface::setOutputLevel(int level)
{
   printLevel_ = level;
   if (printLevel_ < 0)  printLevel_ = 0;
   if (printLevel_ > 10) printLevel_ = 0;
   return 0;
}

// ************************************************************************
// set which driver to use
// ------------------------------------------------------------------------
int FunctionInterface::setDriver(int which)
{
   int  ii, kk, nInps, cnt, status;
   char inString[200], fname[200];
   FILE *fp;
   PsuadeData *psIO;
   pData pPtr;

   if (rsPtrs_ != NULL)
   {
      for (ii = 0; ii < nOutputs_; ii++)
         if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
      delete [] rsPtrs_;
      rsPtrs_ = NULL;
   }
   if (rsIndices_ != NULL)
   {
      delete [] rsIndices_;
      rsIndices_ = NULL;
   }
   if (rsValues_ != NULL)
   {
      delete [] rsValues_;
      rsValues_ = NULL;
   }

   if      (which >= 2) appOptFlag_ = 2;
   else if (which == 1) appOptFlag_ = 1;
   else                 appOptFlag_ = 0;
   useRSModel_ = 0;
   if      (appOptFlag_ == 0) strcpy(fname, appDriver_);
   else if (appOptFlag_ == 1) strcpy(fname, optDriver_);
   if      (appOptFlag_ == 2) strcpy(fname, auxOptDriver_);
   if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
        strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 2;
      return 0; 
   }
   if ( strcmp(appDriver_,"PSUADE_LOCAL") &&
       !strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 3;
      return 0; 
   }
   if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
       !strcmp(optDriver_,"PSUADE_LOCAL"))
   {
      useRSModel_ = 4;
      return 0; 
   }

   fp = fopen(fname, "r");
   if (fp != NULL)
   {
      fscanf(fp, "%10c", inString);
      if (!strncmp(inString, "PSUADE_IO",9)) useRSModel_ = 1;
      if (useRSModel_ == 1)
      {
         psIO = new PsuadeData();
         psIO->setOutputLevel(0);
         status = psIO->readPsuadeFile(fname);
         if (status < 0) exit(1);
         if (status > 0)
         {
            printf("ERROR : cannot read file %d in PSUADE format\n",fname);
            exit(1);
         }
         psIO->getParameter("input_ninputs", pPtr);
         kk = pPtr.intData_;
         if (kk < nInputs_)
         {
            printf("ERROR: nInputs in driver < nInputs in input file.\n");
            printf("** nInputs in driver should be >= nInputs in input file\n");
            exit(1);
         }
         else if (kk >= nInputs_)
         {
            if (kk > nInputs_)
            {
               printf("WARNING: nInputs in driver > nInputs in input file.\n");
               printf("** In order for this to work, you need to specify\n");
               printf("** a rs index file in the driver data file.\n");
            }
            psIO->getParameter("ana_rsindexfile", pPtr);
            if (kk != nInputs_ && !strcmp(pPtr.strArray_[0], "NONE"))
            {
               printf("ERROR: rs index file not specified.\n");
               printf("       nInputs in original .in file = %d\n", nInputs_);
               printf("       nInputs in driver data  file = %d\n", kk);
               exit(1);
            }
            if (strcmp(pPtr.strArray_[0], "NONE"))
            {
               fp = fopen(pPtr.strArray_[0], "r");
               if (fp == NULL)
               {
                  printf("ERROR: rs index file %s not found.\n",
                         pPtr.strArray_[0]);
                  exit(1);
               }
               else
               {
                  printf("INFO: rs index file %s found.\n", pPtr.strArray_[0]);
                  fscanf(fp,"%d", &nInps);
                  nInps_ = nInps;
                  if (kk != nInps)
                  {
                     printf("ERROR: invalid nInputs in rs index file.\n");
                     printf("  data format should be: \n");
                     printf("  line 1: nInputs in rs data (driver) file\n");
                     printf("  line 2: 1 <if selected>\n");
                     printf("  line 3: 2 <if selected>\n");
                     printf("  line 4: 0 <if not selected>\n");
                     printf("  line 5: 4 <if selected>\n");
                     printf("  ...\n");
                     exit(1);
                  }
                  rsIndices_ = new int[nInps];
                  rsValues_ = new double[nInps];
                  cnt = 0;
                  for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
                  for (ii = 0; ii < nInps; ii++)
                  {
                     fscanf(fp, "%d", &kk);
                     if (kk != ii+1)
                     {
                        printf("ERROR: first index in rsIndexFile = %d (must be %d]).\n",
                               rsIndices_[ii], ii+1);
                        printf("  data format should be: \n");
                        printf("  line 1: nInputs in rs data (driver) file\n");
                        printf("  line 2: 1 <num> <default if num == 0>\n");
                        printf("  line 3: 2 <num> <0 if num != 0>\n");
                        printf("  line 4: 3 <num> <default if num == 0>\n");
                        printf("  line 5: 4 <num> <0 if num != 0>\n");
                        printf("  ...\n");
                        exit(1);
                     } 
                     fscanf(fp, "%d", &rsIndices_[ii]);
                     if (rsIndices_[ii] == 0)
                          printf("INFO: input %3d inactive\n",ii+1);
                     if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInps)
                     {
                        printf("INFO: input %3d = %d not valid\n",ii+1,
                               rsIndices_[ii]);
                        exit(1);
                     }
                     rsIndices_[ii]--;
                     fscanf(fp,"%lg", &rsValues_[ii]);
                     if (rsIndices_[ii] >= 0) cnt++;
                     printf("RS_index %5d = %5d %16.8e\n",ii+1,rsIndices_[ii]+1,
                            rsValues_[ii]);
                  }
                  if (cnt != nInputs_)
                  {
                     printf("ERROR: invalid number of active inputs.\n");
                     printf("  data format should be: \n");
                     printf("  line 1: nInputs in rs data (driver) file\n");
                     printf("  line 2: 1 <num> <default if num == 0>\n");
                     printf("  line 3: 2 <num> <0 if num != 0>\n");
                     printf("  line 4: 3 <num> <default if num == 0>\n");
                     printf("  line 5: 4 <num> <0 if num != 0>\n");
                     printf("  ...\n");
                     exit(1);
                  }
               }
            }
         }
         delete psIO;
         rsPtrs_ = new FuncApprox*[nOutputs_];
         for (ii = 0; ii < nOutputs_; ii++)
         {
            printf("Creating response surface for output %d.\n", ii+1);
            rsPtrs_[ii] = (FuncApprox *) genFAFromFile(fname,ii);
            if (rsPtrs_[ii] == NULL)
            {
               printf("FunctionInterface setDriver ERROR: no RS model.\n");
               useRSModel_ = 0;
               exit(1);
            }
         }
      }
   }
   if(fp != NULL) fclose(fp);
   return 0;
}

// ************************************************************************
// run application
// ------------------------------------------------------------------------
// flag = 0, this job has never been launched before 
// flag = 1, just create input file and return (for limited launch)
// flag = 2, check the status of this job
// ------------------------------------------------------------------------
int FunctionInterface::evaluate(int sampleID,int nInputs,double *inputs, 
                                int nOutputs, double *outputs, int flag)
{
   int    ii, outputCount, status, length, outputType;
   double value, *myInputs, stdev;
   char   lineIn[500], command[500];
   char   outfile[500], infile[500];
   FILE   *fp, *fIn, *fOut;

   if (nInputs_ != nInputs || nOutputs_ != nOutputs)
   {
      printf("FunctionInterface ERROR: nInputs/nOutputs mismatch.\n");
      exit(1);
   }

   if ((useRSModel_ == 0 || (useRSModel_ >= 2 && useRSModel_ <= 4)) && 
       rsPtrs_ == NULL)
   {
      if (appOptFlag_ == 0)
      {
	  status = createInputFile(sampleID+1,inputs);
          if (flag == 1 && (status == 0)) return 2;
      }
      else status = 1;

      if (appOptFlag_ == 0)
      {
         length = strlen(appInputTemplate_);
         for (ii = length-1; ii >= 0; ii--)
            if (appInputTemplate_[ii] == '/') break; 
         sprintf(infile, "%s.%d", &(appInputTemplate_[ii+1]), sampleID+1);
      }
      else sprintf(infile, "psuadeOpt.in.%d", sampleID+1);

      if (appOptFlag_ == 0)
      {
	 length = strlen(appOutputTemplate_);
	 for (ii = length-1; ii >= 0; ii--)
	    if (appOutputTemplate_[ii] == '/') break; 
	 sprintf(outfile, "%s.%d", &(appOutputTemplate_[ii+1]), sampleID+1);
      }
      else sprintf(outfile, "psuadeOpt.out.%d", sampleID+1);

      if ((fIn=fopen(outfile, "r")) == NULL) 
      {
         if (flag == 2) return 2;

         if (appOptFlag_ == 0 && (useRSModel_ == 2 || useRSModel_ == 4))
         {
            psLocalFunction(nInputs, inputs, nOutputs, outputs);
            return 0;
         }
         if (appOptFlag_ == 1 && (useRSModel_ == 3 || useRSModel_ == 4))
         {
            psLocalFunction(nInputs, inputs, nOutputs, outputs);
            return 0;
         }

         if (status == 1)
         {
            fOut = fopen(infile, "w");
            if (fOut == NULL)
            {
               printf("FunctionInterface ERROR: cannot open %s file\n", 
                      infile);
               exit(1);
            }
            fprintf(fOut, "%d\n", nInputs);
            for (ii = 0; ii < nInputs; ii++)
               fprintf(fOut, "%20.12e\n", inputs[ii]);
            fclose(fOut);
            if (flag == 1) return 2;
         }

         if (appOptFlag_ == 0)
         {
            if ((!strcmp(appDriver_, "NONE") || !strcmp(appDriver_, "true"))) 
            {
               printf("FunctionInterface ERROR: app driver not found.\n");
               exit(1);
            }
            fp = fopen(appDriver_, "r");
            if (fp == NULL)
            {
               printf("FunctionInterface ERROR: app driver %s not found.\n",
                      appDriver_);
               exit(1);
            }
            else fclose(fp);
         }
         if (appOptFlag_ == 1)
         {
            if ((!strcmp(optDriver_, "NONE") || !strcmp(optDriver_, "true")))
            {
               printf("FunctionInterface ERROR: opt driver not found.\n");
               exit(1);
            }
            fp = fopen(optDriver_, "r");
            if (fp == NULL)
            {
               printf("FunctionInterface ERROR: opt driver %s not found.\n",
                      optDriver_);
               exit(1);
            }
            else fclose(fp);
         }
         if (appOptFlag_ == 2)
         {
            if ((!strcmp(auxOptDriver_, "NONE") || 
                 !strcmp(auxOptDriver_, "true")))
            {
               printf("FunctionInterface ERROR: aux opt driver not found.\n");
               exit(1);
            }
            fp = fopen(auxOptDriver_, "r");
            if (fp == NULL)
            {
               printf("FunctionInterface ERROR: aux opt driver %s not found.\n",
                      auxOptDriver_);
               exit(1);
            }
            else fclose(fp);
         }

         if (appOptFlag_ == 0)
         {
            if (status == 0 || status == 1)
            {
               if (executionMode_ == 1)
                  sprintf(command, "\"%s\" %s %s %d %d %d&",appDriver_, infile,
                          outfile, sampleID, flag, printLevel_);
               else
                  sprintf(command, "\"%s\" %s %s %d %d %d", appDriver_, infile,
                          outfile, sampleID, flag, printLevel_);
            }
         }
         else if (appOptFlag_ == 1)
         {
            if (executionMode_ == 1)
               sprintf(command, "%s %s %s %d %d %d&", optDriver_, infile, outfile,
                       sampleID+1, flag, printLevel_);
            else
               sprintf(command, "%s %s %s %d %d %d", optDriver_, infile, outfile,
                       sampleID+1, flag, printLevel_);
         }
         else if (appOptFlag_ == 2)
         {
            sprintf(command, "%s %s %s %d %d %d", auxOptDriver_, infile, 
                    outfile, sampleID+1, flag, printLevel_);
         }
         if ((appOptFlag_ == 0 && (status == 0 || status == 1)) ||
             (appOptFlag_ > 0))
         {
            if (strstr((const char*) appDriver_, "rm ") != NULL ||
                strstr((const char*) appDriver_, "mv ") != NULL ||
                strstr((const char*) appDriver_, " -f ") != NULL ||
                strstr((const char*) appDriver_, "/bin/") != NULL) 
            {
               printf("FunctionInterface::evaluate ERROR: \n");
               printf("\t\t for security reason do not use rm in executable.\n");
               exit(1);
            }
            system(command);   
            if (executionMode_ == 1 && launchInterval_ > 0)
            {
#ifdef WINDOWS
                  Sleep(1000 * launchInterval_);
#else
		  sleep(launchInterval_);
#endif
            }
         }
      }
      else if (flag == 0)
      {
         printf("WARNING: Output file %s exists before it is run.\n",outfile);
         printf("WARNING: PSUADE will use this output file.\n");
         printf("WARNING: If this is a mistake, stop PSUADE and clean up.\n");
         fclose(fIn);
      }
      else fclose(fIn);

      if (executionMode_ == 0)
      {
	
         if (launchInterval_ == 0)
	   launchInterval_ = 1;

         while ((fIn=fopen(outfile, "r")) == NULL)
         {
            if (printLevel_ > 2)
            {
               printf("Waiting for Job %d to complete.\n", sampleID+1);
               printf("If you run the simulator yourself, use the inputs\n");
               if (appOptFlag_ == 0)
               {
                  printf("from psuadeApps_ct.in.%d for your simulation and\n", 
                         sampleID+1);
                  printf("write the outputs to psuadeApps_ct.out.%d\n",
                         sampleID+1);
               }
               else if (appOptFlag_ == 1)
               {
                  printf("from psuadeOpt.in.%d for your simulation and\n", 
                         sampleID+1);
                  printf("write the outputs to psuadeOpt.out.%d\n",
                         sampleID+1);
               }
            }
#ifdef WINDOWS
	    Sleep(1000 * launchInterval_);
#else
	    sleep(launchInterval_);
#endif

         }
         fclose(fIn);
      }
      if (executionMode_ == 1)
      {
         if ((fIn=fopen(outfile, "r")) == NULL) 
         {
            return 2;
         }
         fclose(fIn);
      }
      length = 0;
      for (ii = 0; ii < 100000; ii++) length = (length + ii) %  32768;

      outputType = 0;

      fIn = fopen(outfile, "r");
      if(fIn == NULL) return 2;
      outputCount = 0;
      if (outputType == 0)
      {
         outputCount = 0;
         for (ii = 0; ii < nOutputs_; ii++)
         {
            fgets(lineIn, 100, fIn);
            sscanf(lineIn, "%lg", &(outputs[outputCount]));
            outputCount++;
            if (feof(fIn) == 1) break; 
         }
      }
      else
      {
#if 0
         char space=' ';
         while (fgets(lineIn, 100, fIn) != NULL)
         {
            for (ii = 0; ii < nOutputs_; ii++)
            {
               const char * strPtr=strstr((const char *) lineIn,
                                          (const char*)outputNames_[ii]);
               if (strPtr != NULL)
               {
                  if (strPtr != lineIn && strPtr[-1] != space) strPtr--;
                  sscanf(strPtr, "%s", lineIn);
                  if (!strcmp(lineIn, (const char*)outputNames_[ii]))
                  {
                     getPatternData(lineIn, outputNames_[ii], &value);
                     outputs[ii] = value;
                     outputCount++;
                     if (feof(fIn) == 1) break; 
                  }
               }
            }
            if (feof(fIn) == 1) break; 
         }
#else
         printf("FunctionInterface: output_template option not currently supported.\n");
         exit(1);
#endif
      }
      fclose(fIn);

      if (outputCount != nOutputs) 
      {
         printf("\t\t output file %s found but nOutputs mismatch (%d,%d).\n", 
                outfile, outputCount, nOutputs);
         printf("\t\t (check your output format).\n");
         return 2;
      }

      if (strcmp(infile, "*"))
      {
         unlink(infile);
      }
      if (strcmp(outfile, "*"))
      {
         unlink(outfile);
      }
   }
   else if (useRSModel_ == 1 && rsPtrs_ != NULL)
   {
      if (rsIndices_ != NULL)
      {
         myInputs = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            if (rsIndices_[ii] < 0) myInputs[ii] = rsValues_[ii];
            else                    myInputs[ii] = inputs[rsIndices_[ii]];
         }
         if (rsRanFlag_ == 0)
         {
            for (ii = 0; ii < nOutputs; ii++)
               outputs[ii] = rsPtrs_[ii]->evaluatePoint(myInputs);
         }
         else
         {
            for (ii = 0; ii < nOutputs; ii++)
            {
               outputs[ii] = rsPtrs_[ii]->evaluatePointFuzzy(myInputs,stdev);
               value = 3.0 * stdev;
               if (rsRanFlag_ == 1) outputs[ii] -= value;
               else                 outputs[ii] += value;
            } 
         }
         delete [] myInputs;
      }
      else
      {
         if (rsRanFlag_ == 0)
         {
            for (ii = 0; ii < nOutputs; ii++)
               outputs[ii] = rsPtrs_[ii]->evaluatePoint(inputs);
         }
         else
         {
            for (ii = 0; ii < nOutputs; ii++)
            {
               outputs[ii] = rsPtrs_[ii]->evaluatePointFuzzy(inputs,stdev);
               value = 3.0 * stdev;
               if (rsRanFlag_ == 1) outputs[ii] -= value;
               else                 outputs[ii] += value;
            }
         }
      }
   }
   else
   {
      printf("FunctionInterface ERROR: evaluate error.\n");
      if (useRSModel_ == 1)
         printf("       Did you forget to declare opt_driver?\n");
      exit(1);
   } 
   return 0;
}

// ************************************************************************
// get number of input variables
// ------------------------------------------------------------------------
int FunctionInterface::getNumInputs()
{
   return nInputs_;
}

// ************************************************************************
// get number of output variables
// ------------------------------------------------------------------------
int FunctionInterface::getNumOutputs()
{
   return nOutputs_;
}

// ************************************************************************
// get current driver code
// ------------------------------------------------------------------------
int FunctionInterface::getDriver()
{
   return appOptFlag_;
}

// ************************************************************************
// get names of input variables
// ------------------------------------------------------------------------
char **FunctionInterface::getInputNames()
{
   return inputNames_;
}

// ************************************************************************
// get names of output variables
// ------------------------------------------------------------------------
char **FunctionInterface::getOutputNames()
{
   return outputNames_;
}

// ************************************************************************
// get application driver
// ------------------------------------------------------------------------
char *FunctionInterface::getApplicationDriver()
{
   return appDriver_;
}

// ************************************************************************
// get optimization driver
// ------------------------------------------------------------------------
char *FunctionInterface::getOptimizationDriver()
{
   return optDriver_;
}

// ************************************************************************
// get auxiliary optimization driver
// ------------------------------------------------------------------------
char *FunctionInterface::getAuxOptimizationDriver()
{
   return auxOptDriver_;
}

// ************************************************************************
// ************************************************************************
// create application input file from input template
// ------------------------------------------------------------------------
int FunctionInterface::createInputFile(int sampleID, double *inputs)
{
   if (strcmp(appInputTemplate_, "psuadeApps_ct.in")) 
   {
     printf("FunctionInterface: user provided input template feature is\n");
     printf("                   no longer supported.\n");
     exit(1);
   }
   return 1;
#if 0
   int  i, stat, length, *trackArray;
   char infile[500], lineIn[500], lineOut[500], lineTemp[500];
   FILE *fIn, *fOut;

   if ((fIn=fopen(appInputTemplate_, "r")) == NULL) 
   {
      printf("FunctionInterface createInputFile ERROR: \n");
      printf("\t\t input template %s does not exist.\n",
             appInputTemplate_);
      exit(1);
   }
   length = strlen(appInputTemplate_);
   for (i = length-1; i >= 0; i--)
      if (appInputTemplate_[i] == '/') break; 


   sprintf(infile, "%s.%d", &(appInputTemplate_[i+1]), sampleID);
   if ((fOut=fopen(infile, "w")) == NULL) 
   {
      printf("FunctionInterface createInputFile ERROR: \n");
      printf("\t\t cannot open input file %s.\n", infile);
      //fclose(fOut);
      exit(1);
   }
   trackArray = new int[nInputs_];
   for (i = 0; i < nInputs_; i++) trackArray[i] = 0;
   while ((fgets(lineIn, 100, fIn) != NULL) && (feof(fIn) == 0))
   {
      strcpy(lineTemp, lineIn);
      for (i = 0; i < nInputs_; i++)
      {
         if (strstr((const char *) lineTemp, (const char*) inputNames_[i]) != NULL) 
         {
            trackArray[i] = 1;
            stat = substitutePattern(lineTemp,lineOut,inputNames_[i],
                                     inputs[i]);
            if ( stat != 0 ) strcpy(lineOut, lineTemp);
            else             strcpy(lineTemp, lineOut);
         } 
      } 
      fprintf(fOut, "%s", lineTemp);
   }
   stat = 0;
   for (i = 0; i < nInputs_; i++)
   {
      if (trackArray[i] == 0)
      {
         printf("FunctionInterface createInputFile ERROR: input \n");
         printf("                   name %s not found in template file.\n",
                inputNames_[i]);
         stat++;
      }
   }
   delete [] trackArray;
   fclose(fIn);
   fclose(fOut);
   if (stat > 0) exit(1);
   return 0;
#endif
}

// ************************************************************************
// Given an input string inputString and a pattern, take out the 
// pattern from the string
// ------------------------------------------------------------------------
int FunctionInterface::removePattern(char *, char *, char *)
{
   return 0;
}

// ************************************************************************
// Given an input string inputString and a pattern, replace the 
// pattern with the given value.
// ------------------------------------------------------------------------
int FunctionInterface::substitutePattern(char *, char *, char *, double)
{
   // Bill Oliver added code to check string lengths
   return 0;
}

// ************************************************************************
// Given an input string inputString and a pattern, find the pattern
// and get the value
// ------------------------------------------------------------------------
int FunctionInterface::getPatternData(char *, char *, double *)
{
   return 0;
}

// ************************************************************************
// Local function
// ------------------------------------------------------------------------
int FunctionInterface::psLocalFunction(int nInputs, double *inputs,
                                       int nOutputs, double *outputs)
{
#if 0
   double W[209];
   double T[209];
   double C[209];
   double O[209];
   W[ 0] =   3.0000e-01; T[ 0] =   2.9815e+02; C[0] =   1.0000e-01; O[0] =   1.0333e+00;
   W[ 1] =   3.0000e-01; T[ 1] =   2.9815e+02; C[1] =   2.1000e-01; O[1] =   1.0534e+00;
   W[ 2] =   3.0000e-01; T[ 2] =   2.9815e+02; C[2] =   3.2000e-01; O[2] =   1.0756e+00;
   W[ 3] =   3.0000e-01; T[ 3] =   2.9815e+02; C[3] =   4.4000e-01; O[3] =   1.0964e+00;
   W[ 4] =   3.0000e-01; T[ 4] =   2.9815e+02; C[4] =   5.6000e-01; O[4] =   1.1142e+00;
   W[ 5] =   3.0000e-01; T[ 5] =   3.1315e+02; C[5] =   1.0000e-01; O[5] =   1.0253e+00;
   W[ 6] =   3.0000e-01; T[ 6] =   3.1315e+02; C[6] =   2.1000e-01; O[6] =   1.0464e+00;
   W[ 7] =   3.0000e-01; T[ 7] =   3.1315e+02; C[7] =   3.2000e-01; O[7] =   1.0669e+00;
   W[ 8] =   3.0000e-01; T[ 8] =   3.1315e+02; C[8] =   4.4000e-01; O[8] =   1.0891e+00;
   W[ 9] =   3.0000e-01; T[ 9] =   3.1315e+02; C[9] =   5.6000e-01; O[9] =   1.1068e+00;
   W[10] =   3.0000e-01; T[10] =   3.2315e+02; C[10] =  1.0000e-01; O[10] =  1.0196e+00;
   W[11] =   3.0000e-01; T[11] =   3.2315e+02; C[11] =  2.1000e-01; O[11] =  1.0412e+00;
   W[12] =   3.0000e-01; T[12] =   3.2315e+02; C[12] =  3.2000e-01; O[12] =  1.0613e+00;
   W[13] =   3.0000e-01; T[13] =   3.2315e+02; C[13] =  4.4000e-01; O[13] =  1.0838e+00;
   W[14] =   3.0000e-01; T[14] =   3.2315e+02; C[14] =  5.6000e-01; O[14] =  1.1014e+00;
   W[15] =   3.0000e-01; T[15] =   3.3315e+02; C[15] =  1.0000e-01; O[15] =  1.0138e+00;
   W[16] =   3.0000e-01; T[16] =   3.3315e+02; C[16] =  2.1000e-01; O[16] =  1.0356e+00;
   W[17] =   3.0000e-01; T[17] =   3.3315e+02; C[17] =  3.2000e-01; O[17] =  1.0556e+00;
   W[18] =   3.0000e-01; T[18] =   3.3315e+02; C[18] =  4.4000e-01; O[18] =  1.0782e+00;
   W[19] =   3.0000e-01; T[19] =   3.3315e+02; C[19] =  5.6000e-01; O[19] =  1.0957e+00;
   W[20] =   3.0000e-01; T[20] =   3.4315e+02; C[20] =  1.0000e-01; O[20] =  1.0076e+00;
   W[21] =   3.0000e-01; T[21] =   3.4315e+02; C[21] =  2.1000e-01; O[21] =  1.0297e+00;
   W[22] =   3.0000e-01; T[22] =   3.4315e+02; C[22] =  3.2000e-01; O[22] =  1.0496e+00;
   W[23] =   3.0000e-01; T[23] =   3.4315e+02; C[23] =  4.4000e-01; O[23] =  1.0723e+00;
   W[24] =   3.0000e-01; T[24] =   3.4315e+02; C[24] =  5.6000e-01; O[24] =  1.0887e+00;
   W[25] =   3.0000e-01; T[25] =   3.5315e+02; C[25] =  1.0000e-01; O[25] =  1.0002e+00;
   W[26] =   3.0000e-01; T[26] =   3.5315e+02; C[26] =  2.1000e-01; O[26] =  1.0234e+00;
   W[27] =   3.0000e-01; T[27] =   3.5315e+02; C[27] =  3.2000e-01; O[27] =  1.0434e+00;
   W[28] =   3.0000e-01; T[28] =   3.5315e+02; C[28] =  4.4000e-01; O[28] =  1.0660e+00;
   W[29] =   3.0000e-01; T[29] =   3.5315e+02; C[29] =  5.6000e-01; O[29] =  1.0812e+00;
   W[30] =   4.0000e-01; T[30] =   2.9815e+02; C[30] =  1.0000e-01; O[30] =  1.0376e+00;
   W[31] =   4.0000e-01; T[31] =   2.9815e+02; C[31] =  2.1000e-01; O[31] =  1.0627e+00;
   W[32] =   4.0000e-01; T[32] =   2.9815e+02; C[32] =  3.3000e-01; O[32] =  1.0945e+00;
   W[33] =   4.0000e-01; T[33] =   2.9815e+02; C[33] =  4.5000e-01; O[33] =  1.1296e+00;
   W[34] =   4.0000e-01; T[34] =   3.1315e+02; C[34] =  1.0000e-01; O[34] =  1.0295e+00;
   W[35] =   4.0000e-01; T[35] =   3.1315e+02; C[35] =  2.1000e-01; O[35] =  1.0547e+00;
   W[36] =   4.0000e-01; T[36] =   3.1315e+02; C[36] =  3.3000e-01; O[36] =  1.0867e+00;
   W[37] =   4.0000e-01; T[37] =   3.1315e+02; C[37] =  4.5000e-01; O[37] =  1.1199e+00;
   W[38] =   4.0000e-01; T[38] =   3.2315e+02; C[38] =  1.0000e-01; O[38] =  1.0237e+00;
   W[39] =   4.0000e-01; T[39] =   3.2315e+02; C[39] =  2.1000e-01; O[39] =  1.0490e+00;
   W[40] =   4.0000e-01; T[40] =   3.2315e+02; C[40] =  3.3000e-01; O[40] =  1.0811e+00;
   W[41] =   4.0000e-01; T[41] =   3.2315e+02; C[41] =  4.5000e-01; O[41] =  1.1138e+00;
   W[42] =   4.0000e-01; T[42] =   3.3315e+02; C[42] =  1.0000e-01; O[42] =  1.0178e+00;
   W[43] =   4.0000e-01; T[43] =   3.3315e+02; C[43] =  2.1000e-01; O[43] =  1.0430e+00;
   W[44] =   4.0000e-01; T[44] =   3.3315e+02; C[44] =  3.3000e-01; O[44] =  1.0752e+00;
   W[45] =   4.0000e-01; T[45] =   3.3315e+02; C[45] =  4.5000e-01; O[45] =  1.1087e+00;
   W[46] =   4.0000e-01; T[46] =   3.4315e+02; C[46] =  1.0000e-01; O[46] =  1.0110e+00;
   W[47] =   4.0000e-01; T[47] =   3.4315e+02; C[47] =  2.1000e-01; O[47] =  1.0367e+00;
   W[48] =   4.0000e-01; T[48] =   3.4315e+02; C[48] =  3.3000e-01; O[48] =  1.0686e+00;
   W[49] =   4.0000e-01; T[49] =   3.4315e+02; C[49] =  4.5000e-01; O[49] =  1.1032e+00;
   W[50] =   4.0000e-01; T[50] =   3.5315e+02; C[50] =  1.0000e-01; O[50] =  1.0048e+00;
   W[51] =   4.0000e-01; T[51] =   3.5315e+02; C[51] =  2.1000e-01; O[51] =  1.0292e+00;
   W[52] =   4.0000e-01; T[52] =   3.5315e+02; C[52] =  3.3000e-01; O[52] =  1.0626e+00;
   W[53] =   4.0000e-01; T[53] =   3.5315e+02; C[53] =  4.5000e-01; O[53] =  1.0963e+00;
   W[54] =   2.0000e-01; T[54] =   3.0315e+02; C[54] =  0.0000e+00; O[54] =  1.0036e+00;
   W[55] =   2.0000e-01; T[55] =   3.0315e+02; C[55] =  1.0000e-01; O[55] =  1.0191e+00;
   W[56] =   2.0000e-01; T[56] =   3.0315e+02; C[56] =  2.0000e-01; O[56] =  1.0324e+00;
   W[57] =   2.0000e-01; T[57] =   3.0315e+02; C[57] =  3.0000e-01; O[57] =  1.0467e+00;
   W[58] =   2.0000e-01; T[58] =   3.0315e+02; C[58] =  4.0000e-01; O[58] =  1.0611e+00;
   W[59] =   2.0000e-01; T[59] =   3.0315e+02; C[59] =  5.0000e-01; O[59] =  1.0744e+00;
   W[60] =   2.0000e-01; T[60] =   3.1315e+02; C[60] =  0.0000e+00; O[60] =  9.9940e-01;
   W[61] =   2.0000e-01; T[61] =   3.1315e+02; C[61] =  1.0000e-01; O[61] =  1.0148e+00;
   W[62] =   2.0000e-01; T[62] =   3.1315e+02; C[62] =  2.0000e-01; O[62] =  1.0281e+00;
   W[63] =   2.0000e-01; T[63] =   3.1315e+02; C[63] =  3.0000e-01; O[63] =  1.0424e+00;
   W[64] =   2.0000e-01; T[64] =   3.1315e+02; C[64] =  4.0000e-01; O[64] =  1.0567e+00;
   W[65] =   2.0000e-01; T[65] =   3.1315e+02; C[65] =  5.0000e-01; O[65] =  1.0698e+00;
   W[66] =   2.0000e-01; T[66] =   3.2315e+02; C[66] =   0.0000e+00; O[66] =   9.9460e-01;
   W[67] =   2.0000e-01; T[67] =   3.2315e+02; C[67] =   1.0000e-01; O[67] =   1.0102e+00;
   W[68] =   2.0000e-01; T[68] =   3.2315e+02; C[68] =   2.0000e-01; O[68] =   1.0233e+00;
   W[69] =   2.0000e-01; T[69] =   3.2315e+02; C[69] =   3.0000e-01; O[69] =   1.0376e+00;
   W[70] =   2.0000e-01; T[70] =   3.2315e+02; C[70] =   4.0000e-01; O[70] =   1.0518e+00;
   W[71] =   2.0000e-01; T[71] =   3.2315e+02; C[71] =   5.0000e-01; O[71] =   1.0649e+00;
   W[72] =   2.0000e-01; T[72] =   3.3315e+02; C[72] =   0.0000e+00; O[72] =   9.8910e-01;
   W[73] =   2.0000e-01; T[73] =   3.3315e+02; C[73] =   1.0000e-01; O[73] =   1.0048e+00;
   W[74] =   2.0000e-01; T[74] =   3.3315e+02; C[74] =   2.0000e-01; O[74] =   1.0180e+00;
   W[75] =   2.0000e-01; T[75] =   3.3315e+02; C[75] =   3.0000e-01; O[75] =   1.0323e+00;
   W[76] =   2.0000e-01; T[76] =   3.3315e+02; C[76] =   4.0000e-01; O[76] =   1.0464e+00;
   W[77] =   2.0000e-01; T[77] =   3.3315e+02; C[77] =   5.0000e-01; O[77] =   1.0594e+00;
   W[78] =   3.0000e-01; T[78] =   3.0315e+02; C[78] =   0.0000e+00; O[78] =   1.0082e+00;
   W[79] =   3.0000e-01; T[79] =   3.0315e+02; C[79] =   1.0000e-01; O[79] =   1.0278e+00;
   W[80] =   3.0000e-01; T[80] =   3.0315e+02; C[80] =   2.0000e-01; O[80] =   1.0475e+00;
   W[81] =   3.0000e-01; T[81] =   3.0315e+02; C[81] =   3.0000e-01; O[81] =   1.0661e+00;
   W[82] =   3.0000e-01; T[82] =   3.0315e+02; C[82] =   4.0000e-01; O[82] =   1.0863e+00;
   W[83] =   3.0000e-01; T[83] =   3.0315e+02; C[83] =   5.0000e-01; O[83] =   1.1060e+00;
   W[84] =   3.0000e-01; T[84] =   3.1315e+02; C[84] =   0.0000e+00; O[84] =   1.0033e+00;
   W[85] =   3.0000e-01; T[85] =   3.1315e+02; C[85] =   1.0000e-01; O[85] =   1.0229e+00;
   W[86] =   3.0000e-01; T[86] =   3.1315e+02; C[86] =   2.0000e-01; O[86] =   1.0427e+00;
   W[87] =   3.0000e-01; T[87] =   3.1315e+02; C[87] =   3.0000e-01; O[87] =   1.0614e+00;
   W[88] =   3.0000e-01; T[88] =   3.1315e+02; C[88] =   4.0000e-01; O[88] =   1.0815e+00;
   W[89] =   3.0000e-01; T[89] =   3.1315e+02; C[89] =   5.0000e-01; O[89] =   1.1010e+00;
   W[90] =   3.0000e-01; T[90] =   3.2315e+02; C[90] =   0.0000e+00; O[90] =   9.9810e-01;
   W[91] =   3.0000e-01; T[91] =   3.2315e+02; C[91] =   1.0000e-01; O[91] =   1.0178e+00;
   W[92] =   3.0000e-01; T[92] =   3.2315e+02; C[92] =   2.0000e-01; O[92] =   1.0376e+00;
   W[93] =   3.0000e-01; T[93] =   3.2315e+02; C[93] =   3.0000e-01; O[93] =   1.0563e+00;
   W[94] =   3.0000e-01; T[94] =   3.2315e+02; C[94] =   4.0000e-01; O[94] =   1.0763e+00;
   W[95] =   3.0000e-01; T[95] =   3.2315e+02; C[95] =   5.0000e-01; O[95] =   1.0958e+00;
   W[96] =   3.0000e-01; T[96] =   3.3315e+02; C[96] =   0.0000e+00; O[96] =   9.9230e-01;
   W[97] =   3.0000e-01; T[97] =   3.3315e+02; C[97] =   1.0000e-01; O[97] =   1.0120e+00;
   W[98] =   3.0000e-01; T[98] =   3.3315e+02; C[98] =   2.0000e-01; O[98] =   1.0320e+00;
   W[99] =   3.0000e-01; T[99] =   3.3315e+02; C[99] =   3.0000e-01; O[99] =   1.0507e+00;
   W[100] =  3.0000e-01; T[100] =  3.3315e+02; C[100] =  4.0000e-01; O[100] =  1.0707e+00;
   W[101] =  3.0000e-01; T[101] =  3.3315e+02; C[101] =  5.0000e-01; O[101] =  1.0901e+00;
   W[102] =  4.0000e-01; T[102] =  3.0315e+02; C[102] =  0.0000e+00; O[102] =  1.0133e+00;
   W[103] =  4.0000e-01; T[103] =  3.0315e+02; C[103] =  1.0000e-01; O[103] =  1.0378e+00;
   W[104] =  4.0000e-01; T[104] =  3.0315e+02; C[104] =  2.0000e-01; O[104] =  1.0628e+00;
   W[105] =  4.0000e-01; T[105] =  3.0315e+02; C[105] =  3.0000e-01; O[105] =  1.0876e+00;
   W[106] =  4.0000e-01; T[106] =  3.0315e+02; C[106] =  4.0000e-01; O[106] =  1.1139e+00;
   W[107] =  4.0000e-01; T[107] =  3.0315e+02; C[107] =  5.0000e-01; O[107] =  1.1396e+00;
   W[108] =  4.0000e-01; T[108] =  3.1315e+02; C[108] =  0.0000e+00; O[108] =  1.0078e+00;
   W[109] =  4.0000e-01; T[109] =  3.1315e+02; C[109] =  1.0000e-01; O[109] =  1.0325e+00;
   W[110] =  4.0000e-01; T[110] =  3.1315e+02; C[110] =  2.0000e-01; O[110] =  1.0576e+00;
   W[111] =  4.0000e-01; T[111] =  3.1315e+02; C[111] =  3.0000e-01; O[111] =  1.0825e+00;
   W[112] =  4.0000e-01; T[112] =  3.1315e+02; C[112] =  4.0000e-01; O[112] =  1.1086e+00;
   W[113] =  4.0000e-01; T[113] =  3.1315e+02; C[113] =  5.0000e-01; O[113] =  1.1344e+00;
   W[114] =  4.0000e-01; T[114] =  3.2315e+02; C[114] =  0.0000e+00; O[114] =  1.0021e+00;
   W[115] =  4.0000e-01; T[115] =  3.2315e+02; C[115] =  1.0000e-01; O[115] =  1.0269e+00;
   W[116] =  4.0000e-01; T[116] =  3.2315e+02; C[116] =  2.0000e-01; O[116] =  1.0521e+00;
   W[117] =  4.0000e-01; T[117] =  3.2315e+02; C[117] =  3.0000e-01; O[117] =  1.0771e+00;
   W[118] =  4.0000e-01; T[118] =  3.2315e+02; C[118] =  4.0000e-01; O[118] =  1.1032e+00;
   W[119] =  4.0000e-01; T[119] =  3.2315e+02; C[119] =  5.0000e-01; O[119] =  1.1290e+00;
   W[120] =  4.0000e-01; T[120] =  3.3315e+02; C[120] =  0.0000e+00; O[120] =  9.9570e-01;
   W[121] =  4.0000e-01; T[121] =  3.3315e+02; C[121] =  1.0000e-01; O[121] =  1.0208e+00;
   W[122] =  4.0000e-01; T[122] =  3.3315e+02; C[122] =  2.0000e-01; O[122] =  1.0462e+00;
   W[123] =  4.0000e-01; T[123] =  3.3315e+02; C[123] =  3.0000e-01; O[123] =  1.0713e+00;
   W[124] =  4.0000e-01; T[124] =  3.3315e+02; C[124] =  4.0000e-01; O[124] =  1.0975e+00;
   W[125] =  4.0000e-01; T[125] =  3.3315e+02; C[125] =  5.0000e-01; O[125] =  1.1232e+00;
   W[126] =  2.0000e-01; T[126] =  2.9815e+02; C[126] =  0.0000e+00; O[126] =  1.0053e+00;
   W[127] =  2.0000e-01; T[127] =  2.9815e+02; C[127] =  1.0000e-01; O[127] =  1.0188e+00;
   W[128] =  2.0000e-01; T[128] =  2.9815e+02; C[128] =  2.0000e-01; O[128] =  1.0327e+00;
   W[129] =  2.0000e-01; T[129] =  2.9815e+02; C[129] =  3.0000e-01; O[129] =  1.0476e+00;
   W[130] =  2.0000e-01; T[130] =  2.9815e+02; C[130] =  4.0000e-01; O[130] =  1.0640e+00;
   W[131] =  2.0000e-01; T[131] =  2.9815e+02; C[131] =  5.0000e-01; O[131] =  1.0800e+00;
   W[132] =  2.0000e-01; T[132] =  3.1315e+02; C[132] =  0.0000e+00; O[132] =  9.9910e-01;
   W[133] =  2.0000e-01; T[133] =  3.1315e+02; C[133] =  1.0000e-01; O[133] =  1.0125e+00;
   W[134] =  2.0000e-01; T[134] =  3.1315e+02; C[134] =  2.0000e-01; O[134] =  1.0264e+00;
   W[135] =  2.0000e-01; T[135] =  3.1315e+02; C[135] =  3.0000e-01; O[135] =  1.0413e+00;
   W[136] =  2.0000e-01; T[136] =  3.1315e+02; C[136] =  4.0000e-01; O[136] =  1.0579e+00;
   W[137] =  2.0000e-01; T[137] =  3.1315e+02; C[137] =  5.0000e-01; O[137] =  1.0735e+00;
   W[138] =  2.0000e-01; T[138] =  3.2315e+02; C[138] =  0.0000e+00; O[138] =  9.9430e-01;
   W[139] =  2.0000e-01; T[139] =  3.2315e+02; C[139] =  1.0000e-01; O[139] =  1.0076e+00;
   W[140] =  2.0000e-01; T[140] =  3.2315e+02; C[140] =  2.0000e-01; O[140] =  1.0215e+00;
   W[141] =  2.0000e-01; T[141] =  3.2315e+02; C[141] =  3.0000e-01; O[141] =  1.0364e+00;
   W[142] =  2.0000e-01; T[142] =  3.2315e+02; C[142] =  4.0000e-01; O[142] =  1.0530e+00;
   W[143] =  2.0000e-01; T[143] =  3.2315e+02; C[143] =  5.0000e-01; O[143] =  1.0680e+00;
   W[144] =  2.0000e-01; T[144] =  3.4315e+02; C[144] =  0.0000e+00; O[144] =  9.8300e-01;
   W[145] =  2.0000e-01; T[145] =  3.4315e+02; C[145] =  1.0000e-01; O[145] =  9.9650e-01;
   W[146] =  2.0000e-01; T[146] =  3.4315e+02; C[146] =  2.0000e-01; O[146] =  1.0105e+00;
   W[147] =  2.0000e-01; T[147] =  3.4315e+02; C[147] =  3.0000e-01; O[147] =  1.0254e+00;
   W[148] =  2.0000e-01; T[148] =  3.4315e+02; C[148] =  4.0000e-01; O[148] =  1.0419e+00;
   W[149] =  2.0000e-01; T[149] =  3.4315e+02; C[149] =  5.0000e-01; O[149] =  1.0570e+00;
   W[150] =  2.0000e-01; T[150] =  3.5315e+02; C[150] =  0.0000e+00; O[150] =  9.7660e-01;
   W[151] =  2.0000e-01; T[151] =  3.5315e+02; C[151] =  1.0000e-01; O[151] =  9.9020e-01;
   W[152] =  2.0000e-01; T[152] =  3.5315e+02; C[152] =  2.0000e-01; O[152] =  1.0043e+00;
   W[153] =  2.0000e-01; T[153] =  3.5315e+02; C[153] =  3.0000e-01; O[153] =  1.0192e+00;
   W[154] =  2.0000e-01; T[154] =  3.5315e+02; C[154] =  4.0000e-01; O[154] =  1.0360e+00;
   W[155] =  3.0000e-01; T[155] =  2.9815e+02; C[155] =  0.0000e+00; O[155] =  1.0106e+00;
   W[156] =  3.0000e-01; T[156] =  2.9815e+02; C[156] =  1.0000e-01; O[156] =  1.0280e+00;
   W[157] =  3.0000e-01; T[157] =  2.9815e+02; C[157] =  2.0000e-01; O[157] =  1.0480e+00;
   W[158] =  3.0000e-01; T[158] =  2.9815e+02; C[158] =  3.0000e-01; O[158] =  1.0700e+00;
   W[159] =  3.0000e-01; T[159] =  2.9815e+02; C[159] =  4.0000e-01; O[159] =  1.0957e+00;
   W[160] =  3.0000e-01; T[160] =  2.9815e+02; C[160] =  5.0000e-01; O[160] =  1.1211e+00;
   W[161] =  3.0000e-01; T[161] =  3.1315e+02; C[161] =  0.0000e+00; O[161] =  1.0034e+00;
   W[162] =  3.0000e-01; T[162] =  3.1315e+02; C[162] =  1.0000e-01; O[162] =  1.0210e+00;
   W[163] =  3.0000e-01; T[163] =  3.1315e+02; C[163] =  2.0000e-01; O[163] =  1.0410e+00;
   W[164] =  3.0000e-01; T[164] =  3.1315e+02; C[164] =  3.0000e-01; O[164] =  1.0629e+00;
   W[165] =  3.0000e-01; T[165] =  3.1315e+02; C[165] =  4.0000e-01; O[165] =  1.0885e+00;
   W[166] =  3.0000e-01; T[166] =  3.1315e+02; C[166] =  5.0000e-01; O[166] =  1.1140e+00;
   W[167] =  3.0000e-01; T[167] =  3.2315e+02; C[167] =  0.0000e+00; O[167] =  9.9810e-01;
   W[168] =  3.0000e-01; T[168] =  3.2315e+02; C[168] =  1.0000e-01; O[168] =  1.0160e+00;
   W[169] =  3.0000e-01; T[169] =  3.2315e+02; C[169] =  2.0000e-01; O[169] =  1.0355e+00;
   W[170] =  3.0000e-01; T[170] =  3.2315e+02; C[170] =  3.0000e-01; O[170] =  1.0580e+00;
   W[171] =  3.0000e-01; T[171] =  3.2315e+02; C[171] =  4.0000e-01; O[171] =  1.0830e+00;
   W[172] =  3.0000e-01; T[172] =  3.2315e+02; C[172] =  5.0000e-01; O[172] =  1.1080e+00;
   W[173] =  3.0000e-01; T[173] =  3.4315e+02; C[173] =  0.0000e+00; O[173] =  9.8580e-01;
   W[174] =  3.0000e-01; T[174] =  3.4315e+02; C[174] =  1.0000e-01; O[174] =  1.0040e+00;
   W[175] =  3.0000e-01; T[175] =  3.4315e+02; C[175] =  2.0000e-01; O[175] =  1.0240e+00;
   W[176] =  3.0000e-01; T[176] =  3.4315e+02; C[176] =  3.0000e-01; O[176] =  1.0464e+00;
   W[177] =  3.0000e-01; T[177] =  3.4315e+02; C[177] =  4.0000e-01; O[177] =  1.0719e+00;
   W[178] =  3.0000e-01; T[178] =  3.5315e+02; C[178] =  0.0000e+00; O[178] =  9.7940e-01;
   W[179] =  3.0000e-01; T[179] =  3.5315e+02; C[179] =  1.0000e-01; O[179] =  9.9700e-01;
   W[180] =  3.0000e-01; T[180] =  3.5315e+02; C[180] =  2.0000e-01; O[180] =  1.0176e+00;
   W[181] =  3.0000e-01; T[181] =  3.5315e+02; C[181] =  3.0000e-01; O[181] =  1.0402e+00;
   W[182] =  3.0000e-01; T[182] =  3.5315e+02; C[182] =  4.0000e-01; O[182] =  1.0660e+00;
   W[183] =  4.0000e-01; T[183] =  2.9815e+02; C[183] =  0.0000e+00; O[183] =  1.0158e+00;
   W[184] =  4.0000e-01; T[184] =  2.9815e+02; C[184] =  1.0000e-01; O[184] =  1.0380e+00;
   W[185] =  4.0000e-01; T[185] =  2.9815e+02; C[185] =  2.0000e-01; O[185] =  1.0630e+00;
   W[186] =  4.0000e-01; T[186] =  2.9815e+02; C[186] =  3.0000e-01; O[186] =  1.0930e+00;
   W[187] =  4.0000e-01; T[187] =  2.9815e+02; C[187] =  4.0000e-01; O[187] =  1.1285e+00;
   W[188] =  4.0000e-01; T[188] =  2.9815e+02; C[188] =  5.0000e-01; O[188] =  1.1597e+00;
   W[189] =  4.0000e-01; T[189] =  3.1315e+02; C[189] =  0.0000e+00; O[189] =  1.0077e+00;
   W[190] =  4.0000e-01; T[190] =  3.1315e+02; C[190] =  1.0000e-01; O[190] =  1.0300e+00;
   W[191] =  4.0000e-01; T[191] =  3.1315e+02; C[191] =  2.0000e-01; O[191] =  1.0550e+00;
   W[192] =  4.0000e-01; T[192] =  3.1315e+02; C[192] =  3.0000e-01; O[192] =  1.0850e+00;
   W[193] =  4.0000e-01; T[193] =  3.1315e+02; C[193] =  4.0000e-01; O[193] =  1.1210e+00;
   W[194] =  4.0000e-01; T[194] =  3.2315e+02; C[194] =  0.0000e+00; O[194] =  1.0018e+00;
   W[195] =  4.0000e-01; T[195] =  3.2315e+02; C[195] =  1.0000e-01; O[195] =  1.0240e+00;
   W[196] =  4.0000e-01; T[196] =  3.2315e+02; C[196] =  2.0000e-01; O[196] =  1.0490e+00;
   W[197] =  4.0000e-01; T[197] =  3.2315e+02; C[197] =  3.0000e-01; O[197] =  1.0797e+00;
   W[198] =  4.0000e-01; T[198] =  3.2315e+02; C[198] =  4.0000e-01; O[198] =  1.1150e+00;
   W[199] =  4.0000e-01; T[199] =  3.4315e+02; C[199] =  0.0000e+00; O[199] =  9.8890e-01;
   W[200] =  4.0000e-01; T[200] =  3.4315e+02; C[200] =  1.0000e-01; O[200] =  1.0120e+00;
   W[201] =  4.0000e-01; T[201] =  3.4315e+02; C[201] =  2.0000e-01; O[201] =  1.0370e+00;
   W[202] =  4.0000e-01; T[202] =  3.4315e+02; C[202] =  3.0000e-01; O[202] =  1.0680e+00;
   W[203] =  4.0000e-01; T[203] =  3.4315e+02; C[203] =  4.0000e-01; O[203] =  1.1040e+00;
   W[204] =  4.0000e-01; T[204] =  3.5315e+02; C[204] =  0.0000e+00; O[204] =  9.8190e-01;
   W[205] =  4.0000e-01; T[205] =  3.5315e+02; C[205] =  1.0000e-01; O[205] =  1.0050e+00;
   W[206] =  4.0000e-01; T[206] =  3.5315e+02; C[206] =  2.0000e-01; O[206] =  1.0310e+00;
   W[207] =  4.0000e-01; T[207] =  3.5315e+02; C[207] =  3.0000e-01; O[207] =  1.0620e+00;
   W[208] =  4.0000e-01; T[208] =  3.5315e+02; C[208] =  4.0000e-01; O[208] =  1.0977e+00;
   if (nInputs != 6)
   {
      printf("psLocalFunction ERROR: nInputs != 5\n");
      exit(1);
   }
   if (nOutputs != 1)
   {
      printf("psLocalFunction ERROR: nOutputs != 1\n"); 
      exit(1);
   }
   double MWco2 = 44.01;
   double MWh2o = 18.01528;
   double MWmea = 61.08;
   double err = 0.0;
   double CO2L, XMEA, TEMP, R;
   double XCO2, XH2O, MWsin, VMEA, VH2O, RHO;
   double a1 = -3.2484e-6;
   double a2 = 0.00165;
   double a3 = 0.793;
   double A = inputs[0];
   double B = inputs[1];
   double CC = inputs[2];
   double D = inputs[3];
   double E = inputs[4];
   for (int ii = 0; ii < 209; ii++)
   {
      CO2L = C[ii];
      R    = W[ii];
      TEMP = T[ii];
      XMEA = 1.0 / (1 + CO2L + (MWmea/MWh2o)*(1-R)/R);
      XCO2 = CO2L * XMEA;
      XH2O = 1 - XMEA - XCO2;
      MWsin = MWmea * XMEA + MWco2 * XCO2 + MWh2o * XH2O;
      VMEA = MWmea / (-5.35162e-7 * TEMP * TEMP - 4.51417e-4 * TEMP + 1.19451);
      VH2O = MWh2o / (a1*TEMP*TEMP+a2*TEMP+a3);
      RHO = MWsin/(XMEA*VMEA+XH2O*VH2O+A*XCO2+(B+CC*XMEA)*XMEA*XH2O+(D+E*XMEA)*XMEA*XCO2);
      err += pow(RHO-O[ii],2.0);
   }
   outputs[0] = err;
#endif
#if 1
   if (nInputs != 7)
   {
      printf("psLocalFunction ERROR: nInputs != 7\n");
      exit(1);
   }
   if (nOutputs != 1)
   {
      printf("psLocalFunction ERROR: nOutputs != 1\n");
      exit(1);
   }
   int nW = 3;
   int nA = 6;
   int nT = 5;
   double *Table = new double[nW * nA * nT];
   Table[0]=1.7; Table[1] =1.8;Table[2] =1.9;Table[3] =1.9;Table[4] =2.1;Table[5] =2.2;
   Table[6]=1.18; Table[7] =1.3;Table[8] =1.3;Table[9] =1.3;Table[10]=1.4;Table[11]=1.6;
   Table[12]=0.95;Table[13]=1.0;Table[14]=1.0;Table[15]=1.1;Table[16]=1.2;Table[17]=1.3;
   Table[18]=0.67;Table[19]=0.7;Table[20]=0.7;Table[21]=0.8;Table[22]=0.8;Table[23]=0.9;
   Table[24]=0.58;Table[25]=0.6;Table[26]=0.6;Table[27]=0.7;Table[28]=0.7;Table[29]=0.8;

   Table[30]=2.48;Table[31]=2.6;Table[32]=2.9;Table[33]=3.1;Table[34]=3.5;Table[35]=3.9;
   Table[36]=1.67;Table[37]=1.7;Table[38]=2.0;Table[39]=2.0;Table[40]=2.4;Table[41]=2.7;
   Table[42]=1.33;Table[43]=1.4;Table[44]=1.6;Table[45]=1.6;Table[46]=1.9;Table[47]=2.1;
   Table[48]=0.92;Table[49]=0.9;Table[50]=1.1;Table[51]=1.1;Table[52]=1.3;Table[53]=1.5;
   Table[54]=0.77;Table[55]=0.8;Table[56]=0.9;Table[57]=0.9;Table[58]=1.1;Table[59]=1.3;

   Table[60]=3.58;Table[61]=4.0;Table[62]=4.6;Table[63]=5.1;Table[64]=6.0;Table[65]=7.0;
   Table[66]=2.28;Table[67]=2.5;Table[68]=3.0;Table[69]=3.3;Table[70]=4.0;Table[71]=4.6;
   Table[72]=1.75;Table[73]=2.0;Table[74]=2.3;Table[75]=2.6;Table[76]=3.1;Table[77]=3.8;
   Table[78]=1.14;Table[79]=1.3;Table[80]=1.5;Table[81]=1.7;Table[82]=2.0;Table[83]=2.3;
   Table[84]=0.95;Table[85]=1.1;Table[86]=1.3;Table[87]=1.4;Table[88]=1.7;Table[89]=1.9;

   /* start processing */
   double *Ts = new double[nT];
   Ts = (double *) malloc(nT * sizeof(double));
   Ts[0] = 25 + 273.15;
   Ts[1] = 40 + 273.15;
   Ts[2] = 50 + 273.15;
   Ts[3] = 70 + 273.15;
   Ts[4] = 80 + 273.15;
   double *As = new double[nA];
   As[0] = 0.0;
   As[1] = 0.1;
   As[2] = 0.2;
   As[3] = 0.3;
   As[4] = 0.4;
   As[5] = 0.5;
   double *Ws = new double[nW];
   Ws[0] = 20;
   Ws[1] = 30;
   Ws[2] = 40;

   double MEA, TEM, CO2;
   double A   = inputs[0];
   double B   = inputs[1];
   double C   = inputs[2];
   double D   = inputs[3];
   double E   = inputs[4];
   double F   = inputs[5];
   double G   = inputs[6];
   double muH2O, dtemp, err=0;
   int    index;
   
   for (int kk = 0; kk < nW; kk++)
   {
      MEA = Ws[kk];
      for (int ii = 0; ii < nT; ii++)
      {
         TEM = Ts[ii];
         muH2O = 1.002 * pow(10.0, 1.3272*(293.15-TEM-0.001053*(TEM-293.15)*(TEM-293.15))/(TEM-168.15));
         for (int jj = 0; jj < nA; jj++)
         {
            CO2 = As[jj];
            dtemp = (A * MEA + B) * TEM + C * MEA + D;
            dtemp = muH2O * exp(dtemp * (CO2 * (E * MEA + F * TEM + G) + 1) * MEA / (TEM * TEM));
            index = kk * nT * nA + ii * nA + jj;
            err = err + pow(dtemp-Table[index],2.0);
         }
      }
   }
   outputs[0] = err;
#endif
#if 0
   int    ii;
   double A, B, V, Y1, Y2, Y3, Y;
   double R0, US, R1, R2, W, E0, Pa, Pb, EdaErr;
   double Yt, Yb, PJ, C, VS[3], Eds[3], vs, Ed0;
   double EdErr, Edd, sigma;

   if (nInputs != 8)
   {
      printf(":psLocalFunction ERROR: nInputs != 8\n");
      exit(1);
   }
   V  = inputs[0];
   A  = inputs[1];
   B  = inputs[2];
   E0 = inputs[3];
   US = inputs[4];
   R1 = inputs[5];
   R2 = inputs[6];
   W  = inputs[7];
   R0 = 1.905;
                                                                                
   VS[0] = 2.20000000e+00; Eds[0] = 4.75e-02;
   VS[1] = 4.40000000e+00; Eds[1] = 5.36e-02;
   VS[2] = 7.20000000e+00; Eds[2] = 5.61e-02;

   Yt = A*(1-W/(R1*V))*exp(-R1*V)+B*(1-W/(R2*V))*exp(-R2*V)+W*E0/V;
   Yb = 1 - W*(1-V)/(2*V);
   PJ = Yt/Yb;
   C = pow(V,1+W)*(PJ - A*exp(-R1*V) - B*exp(-R2*V));
                                                                                
   Pa = R0 * US * US * (1.0 - V);
   Pb = A * exp(-R1*V) + B * exp(-R2*V) + C / pow(V, 1+W);
                                                                                
   Y1 = (Pa - Pb) * 100.0 / Pb;
   Y2 = 0.0;
   for (ii = 0; ii < 3; ii++)
   {
      Ed0 = Eds[ii];
      vs  = VS[ii];
      Edd = E0-A/R1*exp(-vs*R1)-B/R2*exp(-vs*R2)-C/(W*pow(vs,W));
      EdErr = (Edd-Ed0)* 100/Ed0;
      Y2 += pow(EdErr,2.0);
   }
   Y2 /= 3.0;
   Y3 = (sqrt(1.0/R0*(A*R1*exp(-R1*V)+B*R2*exp(-R2*V)+(1+W)*C/pow(V,2+W)))-US)*
        100/US;
   sigma = 5.375e8;
   outputs[0] = 10000.0 * Y1 * Y1 + 0.0001 * Y2;
#endif
   return 0;
}

// ************************************************************************
// create an FunctionInterface instantiation
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterface(PsuadeData *psuadeIO)
{
   int   nInputs, nOutputs, nAppFiles=5;
   char  **inputNames, **outputNames;
   pData pPtr, pINames, pONames, pAppFiles;
   FunctionInterface *funcIO;

   funcIO = new FunctionInterface();
   psuadeIO->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO->getParameter("input_names", pINames);
   inputNames = pINames.strArray_;
   psuadeIO->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO->getParameter("output_names", pONames);
   outputNames = pONames.strArray_;
   psuadeIO->getParameter("app_files", pAppFiles);
   funcIO->loadInputData(nInputs, inputNames);
   funcIO->loadOutputData(nOutputs,outputNames);
   funcIO->loadFunctionData(nAppFiles, pAppFiles.strArray_);
   return funcIO;
}

// ************************************************************************
// create an FunctionInterface instantiation without setting up the
// driver (in case it is a response surface, this will save the time to
// set it up)
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterfaceSimplified(PsuadeData *psuadeIO)
{
   int   nInputs, nOutputs, nAppFiles=5;
   char  **inputNames, **outputNames;
   pData pPtr, pINames, pONames, pAppFiles;
   FunctionInterface *funcIO;

   // ----------------------------------------------------------------
   // set up the FunctionInterface (for sample runs)
   // ----------------------------------------------------------------
   funcIO = new FunctionInterface();
   assert(psuadeIO->getParameter("input_ninputs", pPtr) == 0);
   nInputs = pPtr.intData_;
   assert(psuadeIO->getParameter("input_names", pINames) == 0);
   inputNames = pINames.strArray_;
   assert(psuadeIO->getParameter("output_noutputs", pPtr) == 0);
   nOutputs = pPtr.intData_;
   assert(psuadeIO->getParameter("output_names", pONames) == 0);
   outputNames = pONames.strArray_;
   assert(psuadeIO->getParameter("app_files", pAppFiles) == 0);
   funcIO->loadInputData(nInputs, inputNames);
   funcIO->loadOutputData(nOutputs,outputNames);
   strcpy(pAppFiles.strArray_[0], "NULL");
   funcIO->loadFunctionData(nAppFiles, pAppFiles.strArray_);
   return funcIO;
}

// ************************************************************************
// create an FunctionInterface instantiation using data file for
// response surface
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterfaceGivenAppDriver(int nInputs,
                         int nOutputs, char *fname)
{
   int   ii, nAppFiles=5;
   char  **inputNames, **outputNames, **appFiles;
   FunctionInterface *funcIO;

   appFiles = new char*[nAppFiles];
   for (ii = 0; ii < nAppFiles; ii++)
   {
      appFiles[ii] = new char[500];
      strcpy(appFiles[ii], "NULL");
   }
   strcpy(appFiles[0], fname);

   inputNames = new char*[nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      inputNames[ii] = new char[500];
      strcpy(inputNames[ii], "XX");
   }

   outputNames = new char*[nOutputs];
   for (ii = 0; ii < nOutputs; ii++)
   {
      outputNames[ii] = new char[500];
      strcpy(outputNames[ii], "YY");
   }

   funcIO = new FunctionInterface();
   funcIO->loadInputData(nInputs, inputNames);
   funcIO->loadOutputData(nOutputs,outputNames);
   funcIO->loadFunctionData(nAppFiles, appFiles);

   for (ii = 0; ii < nAppFiles; ii++) delete [] appFiles[ii];
   delete [] appFiles;
   for (ii = 0; ii < nInputs; ii++) delete [] inputNames[ii];
   delete [] inputNames;
   for (ii = 0; ii < nOutputs; ii++) delete [] outputNames[ii];
   delete [] outputNames;

   return funcIO;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FunctionInterface& FunctionInterface::operator=(const FunctionInterface &)
{
   printf("FunctionInterface operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


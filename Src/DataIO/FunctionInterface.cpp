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
   int  ii, kk, nInps, cnt;
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
         psIO->readPsuadeFile(fname);
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
   int  ii, kk, nInps, cnt;
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
         psIO->readPsuadeFile(fname);
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
                  sprintf(command, "%s %s %s %d %d %d&",appDriver_, infile,
                          outfile, sampleID, flag, printLevel_);
               else
                  sprintf(command, "%s %s %s %d %d %d", appDriver_, infile,
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
               sprintf(command, "sleep %d", launchInterval_);
               system(command);   
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
            sprintf(command, "sleep 1");
         else
            sprintf(command, "sleep %d", launchInterval_);
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
            system(command);
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
   int    ii;
   double A, B, V, Y1, Y2, Y3, Y;
   double R0, US, R1, R2, W, E0, Pa, Pb, EdaErr;
   double Yt, Yb, PJ, C, VS[3], Eds[3], vs, Ed0;
   double EdErr, Edd, sigma;

   if (nInputs != 8)
   {
      printf(":psLocalFunction ERROR: nInputs != 4\n");
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

#if 0
   int    i, j;
   double X[25], Y[25], Z[25], A = 0.25, B = 0.75, error;

   if (nInputs != 2)
   {
      printf(":psLocalFunction ERROR: nInputs != 4\n");
      exit(1);
   }
   for (i = 0; i < 5; i++)
   {
      for (j = 0; j < 5; j++)
      {
         X[i*5+j] = i * 0.25;
         Y[i*5+j] = j * 0.25;
      }
   }
   for (i = 0; i < 25; i++)
      Z[i] = A * X[i] * X[i] + B * Y[i] * Y[i];

   outputs[0] = 0.0;
   for (i = 0; i < 25; i++)
   {
      error = Z[i] - inputs[0] * X[i] * X[i] - inputs[1] * Y[i] * Y[i];
      outputs[0] += error * error;
   }
   outputs[0] /= 25.0;
   outputs[0] = exp(-0.5*outputs[0]/pow(0.104,2.0));
#endif

#if 0
   int    i, j;
   double X[25], Y[25], Z[25], A = 0.25, B = 0.75, error;
   double A2 = 0.75, B2 = 0.25;

   if (nInputs != 4)
   {
      printf(":psLocalFunction ERROR: nInputs != 4\n");
      exit(1);
   }
   for (i = 0; i < 5; i++)
   {
      for (j = 0; j < 5; j++)
      {
         X[i*5+j] = i * 0.25;
         Y[i*5+j] = j * 0.25;
      }
   }
   for (i = 0; i < 25; i++)
      Z[i] = A  * X[i] * X[i] + B  * Y[i] * Y[i] +
             A2 * X[i] * X[i] + B2 * Y[i] * Y[i];

   outputs[0] = 0.0;
   for (i = 0; i < 25; i++)
   {
      error = Z[i] - inputs[0] * X[i] * X[i] - inputs[1] * Y[i] * Y[i] -
              inputs[2] * X[i] * X[i] - inputs[3] * Y[i] * Y[i];
      outputs[0] += error * error;
   }
   outputs[0] /= 25.0;
   outputs[0] = exp(-0.5*outputs[0]/pow(0.076,2.0));
#endif

#if 0
   double pi=3.14159;
   outputs[0] = sin(2.0*pi*inputs[0]) * sin(2.0*pi*inputs[1]);
#endif
#if 1
   /* Table 5 : t=25-80, w=20%-100%, load=0 */

   int nW = 3, ii, jj;
   int nA = 6;
   int nT = 5;
   int    index;
   double xwMEA, T, T2, muH2O, muH2O_20, CO2Load, muX, muBlend, ddata;
   muH2O_20 = 1.002; /* mPa-s at 20C */
   double *Table = (double *) malloc(nW * nA * nT*sizeof(double));
   Table[0] =1.7; Table[1] =1.8;Table[2] =1.9;Table[3] =1.9;Table[4] =2.1;Table[5] =2.2;
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
   double *Ts = (double *) malloc(nT * sizeof(double));
   Ts[0] = 25;
   Ts[1] = 40;
   Ts[2] = 50;
   Ts[3] = 70;
   Ts[4] = 80;
   double *As = (double *) malloc(nA * sizeof(double));
   As[0] = 0.0;
   As[1] = 0.1;
   As[2] = 0.2;
   As[3] = 0.3;
   As[4] = 0.4;
   As[5] = 0.5;
   double *Ws = (double *) malloc(nW * sizeof(double));
   Ws[0] = 20;
   Ws[1] = 30;
   Ws[2] = 40;
   double A = 0;
   double B = 0;
   double C = inputs[0];
   double D = inputs[1];
   double E = inputs[2];
   double F = inputs[3];
   double G = inputs[4];
   double rms = 0.0;
   double nrms = 0.0;
   double avg = 0.0;
   double navg = 0.0;
   for (int kk = 0; kk < nW; kk++)
   {
      xwMEA = Ws[kk];
      for (ii = 0; ii < nT; ii++)
      {
         T = Ts[ii];
         T2 = 1.0 / ((T+273.15) * (T+273.15));
         muH2O = muH2O_20 *
                 pow(10.0, (1.3272*(20-T)-0.001053*(T-20)*(T-20))/(T+105));
         for (jj = 0; jj < nA; jj++)
         {
            CO2Load = As[jj];
            muBlend = (A*xwMEA + B) * (T+273.15) + (C*xwMEA + D);
            muBlend = muBlend *(CO2Load*(E*xwMEA + F*(T+273.15) + G)+1)*xwMEA;
            muBlend = muBlend * T2;
            muBlend = exp(muBlend);
            muX = muBlend * muH2O;
            index = kk * nT * nA + ii * nA + jj;
            ddata = muX - Table[index];
            rms += ddata * ddata;
            ddata = ddata / Table[index];
            nrms += ddata * ddata;
            if (muX >= Table[index]) avg += muX - Table[index];
            else                     avg += Table[index] - muX;
            if (ddata >= 0) navg += ddata;
            else            navg -= ddata;
         }
      }
   }
   rms /= (nW * nA * nT);
   rms = sqrt(rms);
   nrms /= (nW * nA * nT);
   rms = sqrt(nrms);
   avg /= (nW * nA * nT);
   navg /= (nW * nA * nT);
   free(Ts);
   free(As);
   free(Ws);
   free(Table);
   outputs[0] = avg;
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


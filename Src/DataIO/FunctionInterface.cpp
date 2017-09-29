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
#include "dtype.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "FuncApprox.h"
#include "pData.h"
#include "PDFBase.h"
#include "PsuadeData.h"
#include "Psuade.h"

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
  strcpy(ensembleDriver_, "true");
  strcpy(ensembleOptDriver_, "true");
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
  rsnInps_ = 0;
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
FunctionInterface::FunctionInterface(const FunctionInterface & fi)
{
  int ii;

  nInputs_ = fi.nInputs_;
  nOutputs_ = fi.nOutputs_;
  strcpy(appDriver_, fi.appDriver_);
  strcpy(optDriver_, fi.optDriver_);
  strcpy(auxOptDriver_, fi.auxOptDriver_);
  strcpy(ensembleDriver_, fi.ensembleDriver_);
  strcpy(ensembleOptDriver_, fi.ensembleOptDriver_);
  strcpy(appInputTemplate_, fi.appInputTemplate_);
  strcpy(appOutputTemplate_, fi.appOutputTemplate_);
  rsRanFlag_ = fi.rsRanFlag_;
  printLevel_ = fi.printLevel_;
  rsnInps_ = fi.rsnInps_;

  inputNames_ = new char*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++)
  {
    inputNames_[ii] = new char[80];
    strcpy(inputNames_[ii], fi.inputNames_[ii]);
  }
 
  outputNames_ = new char *[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++)
  {
    outputNames_[ii] = new char[80];
    strcpy(outputNames_[ii], fi.outputNames_[ii]);
  }

  rsPtrs_ = new FuncApprox*[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) rsPtrs_[ii] = fi.rsPtrs_[ii];

  rsIndices_ = new int[rsnInps_];
  rsValues_ = new double[rsnInps_];
  for (ii = 0; ii < rsnInps_; ii++)
  {
    rsIndices_[ii] = fi.rsIndices_[ii];
    rsValues_[ii] = fi.rsValues_[ii];
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
  int    ii, kk, nInps, status;
  char   inString[200], fname[200], *cString, pString[1001];
  char   winput1[1001], winput2[1001], winput3[1001];
  double ddata;
  FILE   *fp;
  PsuadeData *psIO;
  pData pPtr;

  if (length < 5)
  {
    printf("FunctionInterface loadFunction ERROR: length < 5\n");
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
  if (length >= 6) strcpy(ensembleDriver_, names[5]);
  if (length >= 7) strcpy(ensembleOptDriver_, names[6]);

  if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
       strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 2;
    return 0; 
  }
  if (strcmp(appDriver_,"PSUADE_LOCAL") &&
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
      if (status != 0)
      {
        printf("ERROR : cannot read file %s in PSUADE format\n",fname);
        exit(1);
      }
      psIO->getParameter("input_ninputs", pPtr);
      int rsModelnInps = pPtr.intData_;
      if (rsModelnInps < nInputs_ || rsModelnInps > nInputs_)
      {
        printAsterisks(PL_INFO,0);
        printf("FunctionInterface: setting up RS driver from %s.\n",
               fname);
        printEquals(PL_INFO,0);
        printf("WARNING: nInputs in RS driver does not match nInputs");
        printf(" in original psuade file.\n");
        printf("   nInputs in original psuade file = %d\n",nInputs_);
        printf("   nInputs in RS driver data file  = %d\n",rsModelnInps);
        psIO->getParameter("ana_rsindexfile", pPtr);
        if (!strcmp(pPtr.strArray_[0], "NONE"))
        {
          psIO->getParameter("input_names", pPtr);
          if (psConfig_ != NULL && pPtr.strArray_ != NULL)
          {
            cString = psConfig_->getParameter("num_fixed");
            if (cString != NULL)
            {
              sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
              if (nInps > 0)
              {
                rsIndices_ = new int[rsModelnInps];
                rsValues_ = new double[rsModelnInps];
                for (ii = 0; ii < rsModelnInps; ii++) 
                  rsIndices_[ii] = ii;
                for (ii = 0; ii < nInps; ii++)
                {
                  sprintf(pString,"fixed-%d", ii+1);
                  cString = psConfig_->getParameter(pString);
                  if (cString != NULL)
                  {
                    sscanf(cString,"%s %s %s %lg",winput1,winput2,
                           winput3, &ddata);
                    for (kk = 0; kk < rsModelnInps; kk++)
                    {
                      if (!strcmp(winput2,pPtr.strArray_[kk]))
                      {
                        rsIndices_[kk] = -1;
                        rsValues_[kk] = ddata;
                        printf("Input %4d fixed at %e\n",kk+1,ddata);
                      }
                    }
                    if (kk == rsModelnInps) break;
                  }
                }
                if (ii == nInps)
                {
                  delete [] rsIndices_;
                  delete [] rsValues_;
                  rsIndices_ = NULL;
                  rsValues_ = NULL;
                }
                else rsnInps_ = rsModelnInps;
              }
            }
          }
          if (rsIndices_ == NULL)
          { 
            printf("ERROR: nInputs mismatch and missing rs_index_file.\n");
            printf("   nInputs in original psuade file = %d\n",nInputs_);
            printf("   nInputs in RS driver data file  = %d\n",
                   rsModelnInps);
            printf("ADVICE: Put rs_index_file in %s or in the",fname);
            printf(" original psuade file.\n");
            exit(1);
          }
        }
        else
        {
          printf("WARNING: rs_index_file found in the RS driver file.\n");
          fp = fopen(pPtr.strArray_[0], "r");
          if (fp == NULL)
          {
            printf("ERROR: missing rs_index_file %s in current folder.\n",
                   pPtr.strArray_[0]);
            exit(1);
          }
          else
          {
            fscanf(fp,"%d", &nInps);
            rsnInps_ = nInps;
            if (rsModelnInps != nInps)
            {
              printf("ERROR: invalid nInputs in rs index file.\n");
              printf("  It has to match nInputs in the RS data file.\n");
              printf("  nInputs read     = %d\n", nInps);
              printf("  nInputs expected = %d\n", rsModelnInps);
              printf("  Data format in rs index file should be: \n");
              printf("  line 1: nInputs in RS driver data file\n");
              printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
              printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
              printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
              printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
              printf("  ...\n");
              exit(1);
            }
            rsIndices_ = new int[nInps];
            rsValues_ = new double[nInps];
            for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
            for (ii = 0; ii < nInps; ii++)
            {
              fscanf(fp, "%d", &kk);
              if (kk != ii+1)
              {
                printf("ERROR: first index in rs index file = %d.\n",
                       rsIndices_[ii]);
                printf("       Must be equal to %d.\n",ii+1);
                printf("  Data format in rs index file should be: \n");
                printf("  line 1: nInputs in RS driver data file\n");
                printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
                printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
                printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
                printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
                printf("  ...\n");
                exit(1);
              } 
              fscanf(fp, "%d", &rsIndices_[ii]);
              if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInputs_)
              {
                printf("ERROR: input %3d = %d not valid\n",ii+1,
                       rsIndices_[ii]);
                printf("       Need to be between 1 and %d\n",nInputs_);
                exit(1);
              }
              rsIndices_[ii]--;
              fscanf(fp, "%lg", &rsValues_[ii]);
              if (rsIndices_[ii] == -1)
                printf("   RS Input %d inactive, fixed at %16.8e\n",
                       ii+1, rsValues_[ii]);
              else
              {
                printf("   RS Input %d   active, mapped to",ii+1);
                printf(" PSUADE input %d\n", rsIndices_[ii]+1);
              }
            }
            fclose(fp);
          }
        }
      }
      psIO->getParameter("input_names", pPtr);
      if (psConfig_ != NULL && rsIndices_ == NULL && pPtr.strArray_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
          if (nInps > 0)
          {
            rsIndices_ = new int[rsModelnInps];
            rsValues_ = new double[rsModelnInps];
            for (ii = 0; ii < rsModelnInps; ii++) rsIndices_[ii] = ii;
            for (ii = 0; ii < nInps; ii++)
            {
              sprintf(pString,"fixed-%d", ii+1);
              cString = psConfig_->getParameter(pString);
              if (cString != NULL)
              {
                sscanf(cString,"%s %s %s %lg",winput1,winput2,
                       winput3, &ddata);
                for (kk = 0; kk < rsModelnInps; kk++)
                {
                  if (!strcmp(winput2,pPtr.strArray_[kk]))
                  {
                    rsIndices_[kk] = -1;
                    rsValues_[kk] = ddata;
                    break;
                  }
                }
                if (kk == rsModelnInps) break;
              }
            }
            if (ii != nInps)
            {
              printf("WARNING: Config info on fixed variables not used.\n");
              delete [] rsIndices_;
              delete [] rsValues_;
              rsIndices_ = NULL;
              rsValues_ = NULL;
            }
            else rsnInps_ = rsModelnInps;
          }
          for (kk = 0; kk < rsModelnInps; kk++)
          {
            if (rsIndices_[kk] == -1) 
              printf("Input %4d fixed at %e\n",kk+1,rsValues_[kk]);
          }
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
  int    ii, kk, nInps, status;
  char   inString[200], fname[200], winput1[1001], winput2[1001];
  char   winput3[1001], *cString, pString[1001];
  double ddata;
  FILE   *fp;
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
    fclose(fp);
    if (useRSModel_ == 1)
    {
      psIO = new PsuadeData();
      psIO->setOutputLevel(0);
      status = psIO->readPsuadeFile(fname);
      if (status != 0)
      {
        printf("ERROR : cannot read file %s in PSUADE format\n",fname);
        exit(1);
      }
      psIO->getParameter("input_ninputs", pPtr);
      int rsModelnInps = pPtr.intData_;
      if (rsModelnInps < nInputs_ || rsModelnInps > nInputs_)
      {
        printAsterisks(PL_INFO,0);
        printf("FunctionInterface: setting up RS driver from %s.\n",
               fname);
        printEquals(PL_INFO,0);
        printf("WARNING: nInputs in RS driver does not match nInputs");
        printf(" in original psuade file.\n");
        printf("   nInputs in original psuade file = %d\n",nInputs_);
        printf("   nInputs in RS driver data file  = %d\n",
               rsModelnInps);
        psIO->getParameter("ana_rsindexfile", pPtr);
        if (!strcmp(pPtr.strArray_[0], "NONE"))
        {
          printf("ERROR: nInputs mismatch and missing rs_index_file.\n");
          printf("   nInputs in original psuade file = %d\n",nInputs_);
          printf("   nInputs in RS driver data file  = %d\n",
                 rsModelnInps);
          printf("ADVICE: Put rs_index_file in %s or in the",fname);
          printf(" original psuade file.\n");
          exit(1);
        }
        printf("WARNING: rs_index_file found in the RS driver file.\n");
        fp = fopen(pPtr.strArray_[0], "r");
        if (fp == NULL)
        {
          printf("ERROR: missing rs_index_file %s in current folder.\n",
                 pPtr.strArray_[0]);
          exit(1);
        }
        else
        {
          fscanf(fp,"%d", &nInps);
          rsnInps_ = nInps;
          if (rsModelnInps != nInps)
          {
            printf("ERROR: invalid nInputs in rs index file.\n");
            printf("  nInputs read     = %d\n", nInps);
            printf("  nInputs expected = %d\n", rsModelnInps);
            printf("  Data format in rs index file should be: \n");
            printf("  line 1: nInputs in RS driver data file\n");
            printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
            printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
            printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
            printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
            printf("  ...\n");
            exit(1);
          }
          rsIndices_ = new int[nInps];
          rsValues_ = new double[nInps];
          for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
          for (ii = 0; ii < nInps; ii++)
          {
            fscanf(fp, "%d", &kk);
            if (kk != ii+1)
            {
              printf("ERROR: first index in rs index file = %d.\n",
                     rsIndices_[ii]);
              printf("       Must be equal to %d.\n",ii+1);
              printf("  Data format in rs index file should be: \n");
              printf("  line 1: nInputs in RS driver data file\n");
              printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
              printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
              printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
              printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
              printf("  ...\n");
              exit(1);
            }
            fscanf(fp, "%d", &rsIndices_[ii]);
            if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInputs_)
            {
              printf("INFO: input %3d = %d not valid\n",ii+1,
                     rsIndices_[ii]);
              printf("       Need to be between 1 and %d\n",nInputs_);
              exit(1);
            }
            rsIndices_[ii]--;
            fscanf(fp, "%lg", &rsValues_[ii]);
            if (rsIndices_[ii] == -1)
              printf("   RS Input %d inactive, fixed at %16.8e\n",
                     ii+1, rsValues_[ii]);
            else
            {
              printf("   RS Input %d   active, mapped to",ii+1);
              printf(" PSUADE input %d\n", rsIndices_[ii]+1);
            }
          }
          fclose(fp);
        }
      }
      psIO->getParameter("input_names", pPtr);
      if (psConfig_ != NULL && rsIndices_ == NULL && pPtr.strArray_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
          if (nInps > 0)
          {
            rsIndices_ = new int[rsModelnInps];
            rsValues_ = new double[rsModelnInps];
            for (ii = 0; ii < rsModelnInps; ii++) rsIndices_[ii] = ii;
            for (ii = 0; ii < nInps; ii++)
            {
              sprintf(pString,"fixed-%d", ii+1);
              cString = psConfig_->getParameter(pString);
              if (cString != NULL)
              {
                sscanf(cString,"%s %s %s %lg",winput1,winput2,
                       winput3, &ddata);
                for (kk = 0; kk < rsModelnInps; kk++)
                {
                  if (!strcmp(winput2,pPtr.strArray_[kk]))
                  {
                    rsIndices_[kk] = -1;
                    rsValues_[kk] = ddata;
                    break;
                  }
                }
                if (kk == rsModelnInps) break;
              }
            }
            if (ii != nInps)
            {
              printf("WARNING: Config info on fixed variables not used.\n");
              delete [] rsIndices_;
              delete [] rsValues_;
              rsIndices_ = NULL;
              rsValues_ = NULL;
            }
            else rsnInps_ = nInps;
          }
          for (kk = 0; kk < rsModelnInps; kk++)
          {
            if (rsIndices_[kk] == -1) 
              printf("Input %4d fixed at %e\n",kk+1,rsValues_[kk]);
          }
        }
      }
      delete psIO;
      rsPtrs_ = new FuncApprox*[nOutputs_];
      for (ii = 0; ii < nOutputs_; ii++)
      {
        if (printLevel_ > 3)
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
  int    ii, outputCount, length, nfixed;
  double value, *myInputs, stdev;
  char   lineIn[500], command[500], winput1[500], winput2[500];
  char   outfile[500], infile[500], *cString;
  FILE   *fp, *fIn, *fOut;

  if (nInputs_ != nInputs || nOutputs_ != nOutputs)
  {
    printf("FunctionInterface evaluate ERROR: nInputs/nOutputs mismatch.\n");
    printf("   nInputs  = %d versus %d (local)\n", nInputs, nInputs_);
    printf("   nOutputs = %d versus %d (local)\n", nOutputs, nOutputs_);
    exit(1);
  }

  if ((useRSModel_ == 0 || (useRSModel_ >= 2 && useRSModel_ <= 4)) && 
      rsPtrs_ == NULL)
  {
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

      fOut = fopen(infile, "w");
      if (fOut == NULL)
      {
         printf("FunctionInterface ERROR: cannot open %s file\n",infile);
         exit(1);
      }
      fprintf(fOut, "%d\n", nInputs);
      for (ii = 0; ii < nInputs; ii++)
        fprintf(fOut, "%20.12e\n", inputs[ii]);
      if (psConfig_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nfixed);
          fprintf(fOut,"num_fixed = %d\n", nfixed);
          ii = 0;
          while (ii < nfixed)
          {
            ii++;
            sprintf(winput1, "fixed-%d=",ii);
            cString = psConfig_->getParameter(winput1);
            if (cString != NULL)
            {
              sscanf(cString, "%s %s %lg", winput1, winput2, &value);
              fprintf(fOut,"fixed %d %s = %24.16e\n",ii,winput2,value);
            }
          }
        }
      }
      fclose(fOut);
      if (flag == 1) return 2;

      if (appOptFlag_ == 0)
      {
        if ((!strcmp(appDriver_, "NONE") || !strcmp(appDriver_, "true"))) 
        {
          printf("FunctionInterface ERROR: app driver not found.\n");
          exit(1);
        }
      }
      if (appOptFlag_ == 1)
      {
        if ((!strcmp(optDriver_, "NONE") || !strcmp(optDriver_, "true")))
        {
          printf("FunctionInterface ERROR: opt driver not found.\n");
          exit(1);
        }
      }
      if (appOptFlag_ == 2)
      {
        if ((!strcmp(auxOptDriver_, "NONE") || 
             !strcmp(auxOptDriver_, "true")))
        {
          printf("FunctionInterface ERROR: aux opt driver not found.\n");
          exit(1);
        }
      }

      if (appOptFlag_ == 0)
      {
        if (executionMode_ == 1)
          sprintf(command, "\"%s\" %s %s %d %d %d&",appDriver_,infile,
                  outfile, sampleID, flag, printLevel_);
        else
          sprintf(command, "\"%s\" %s %s %d %d %d", appDriver_,infile,
                   outfile, sampleID, flag, printLevel_);
        if (strstr((const char*) appDriver_, "rm ") != NULL ||
            strstr((const char*) appDriver_, "mv ") != NULL ||
            strstr((const char*) appDriver_, " -f ") != NULL ||
            strstr((const char*) appDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
      }
      else if (appOptFlag_ == 1)
      {
        if (executionMode_ == 1)
          sprintf(command, "%s %s %s %d %d %d&",optDriver_,infile,outfile,
                  sampleID+1, flag, printLevel_);
        else
          sprintf(command, "%s %s %s %d %d %d",optDriver_,infile,outfile,
                  sampleID+1, flag, printLevel_);
        if (strstr((const char*) optDriver_, "rm ") != NULL ||
            strstr((const char*) optDriver_, "mv ") != NULL ||
            strstr((const char*) optDriver_, " -f ") != NULL ||
            strstr((const char*) optDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
      }
      else if (appOptFlag_ == 2)
      {
        sprintf(command, "%s %s %s %d %d %d", auxOptDriver_, infile, 
                outfile, sampleID+1, flag, printLevel_);
        if (strstr((const char*) auxOptDriver_, "rm ") != NULL ||
            strstr((const char*) auxOptDriver_, "mv ") != NULL ||
            strstr((const char*) auxOptDriver_, " -f ") != NULL ||
            strstr((const char*) auxOptDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
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
      //if (launchInterval_ == 0) launchInterval_ = 1;

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
            printf("write the outputs to psuadeOpt.out.%d\n",sampleID+1);
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

    fIn = fopen(outfile, "r");
    if (fIn == NULL) return 2;
    outputCount = 0;
    for (ii = 0; ii < nOutputs_; ii++)
    {
      fgets(lineIn, 100, fIn);
      sscanf(lineIn, "%lg", &(outputs[outputCount]));
      outputCount++;
      if (feof(fIn) == 1) break; 
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
      myInputs = new double[rsnInps_];
      for (ii = 0; ii < rsnInps_; ii++)
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
      printf("       Did you forget to declare driver?\n");
    exit(1);
  } 
  return 0;
}

// ************************************************************************
// run ensemble simulation
// ------------------------------------------------------------------------
int FunctionInterface::ensembleEvaluate(int nSamp,int nInputs,double *inputs, 
                               int nOutputs, double *outputs, int ID)
{
  int    ii, ss, outputCount, nfixed;
  double value;
  char   outfile[500], infile[500], lineIn[5000], command[500];
  char   winput1[500], winput2[500], *cString;
  FILE   *fp, *fIn, *fOut;

  if (nInputs_ != nInputs || nOutputs_ != nOutputs)
  {
    printf("FunctionInterface ERROR: nInputs/nOutputs mismatch.\n");
    printf("   nInputs  = %d versus %d (local)\n", nInputs, nInputs_);
    printf("   nOutputs = %d versus %d (local)\n", nOutputs, nOutputs_);
    printf("NOTE: ensembleEvaluate expects driver is an actual simulator");
    printf(" and not RS data file.\n");
    printf("      Therefore, it does not take a rs_index_file.\n");
    exit(1);
  }

  sprintf(infile, "psuadeEval.in.%d", ID);
  sprintf(outfile, "psuadeEval.out.%d", ID);

  if ((fIn=fopen(outfile, "r")) == NULL) 
  {
    if (appOptFlag_ == 0)
    {
      if (!strcmp(ensembleDriver_,"PSUADE_LOCAL"))
      {
        psEnsembleLocalFunction(nSamp,nInputs,inputs,nOutputs,outputs);
        return 0;
      } 
    }
    else
    {
      if (!strcmp(ensembleOptDriver_,"PSUADE_LOCAL"))
      {
        psEnsembleLocalFunction(nSamp,nInputs,inputs,nOutputs,outputs);
        return 0;
      } 
    }

    fOut = fopen(infile, "w");
    if (fOut == NULL)
    {
      printf("FunctionInterface ERROR: cannot open %s file\n",infile);
      exit(1);
    }
    fprintf(fOut, "%d %d\n", nSamp, nInputs);
    for (ss = 0; ss < nSamp; ss++)
    {
      for (ii = 0; ii < nInputs; ii++)
        fprintf(fOut, "%20.12e ", inputs[ss*nInputs+ii]);
      fprintf(fOut, "\n");
    }
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("num_fixed");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nfixed);
        fprintf(fOut,"num_fixed = %d\n", nfixed);
        ii = 0;
        while (ii < nfixed)
        {
          ii++;
          sprintf(winput1, "fixed-%d=",ii);
          cString = psConfig_->getParameter(winput1);
          if (cString != NULL)
          {
            sscanf(cString, "%s %s %lg", winput1, winput2, &value);
            fprintf(fOut,"fixed %d %s = %24.16e\n", ii, winput2, value);
          }
        }
      }
    }
    fclose(fOut);

    if (appOptFlag_ == 0)
    {
      if ((!strcmp(ensembleDriver_, "NONE") || 
           !strcmp(ensembleDriver_, "true")))
      {
        printf("FunctionInterface ERROR: ensemble driver not set.\n");
        exit(1);
      }
      fp = fopen(ensembleDriver_, "r");
      if (fp == NULL)
      {
        printf("FunctionInterface ERROR: ensemble driver %s not found.\n",
               ensembleDriver_);
        exit(1);
      }
      else fclose(fp);
      sprintf(command, "%s %s %s %d %d", ensembleDriver_,infile,outfile,ID, 
              printLevel_);
    }
    else
    {
      if ((!strcmp(ensembleOptDriver_, "NONE") || 
           !strcmp(ensembleOptDriver_, "true")))
      {
        printf("FunctionInterface ERROR: ensemble opt driver not set.\n");
        exit(1);
      }
      fp = fopen(ensembleOptDriver_, "r");
      if (fp == NULL)
      {
        printf("FunctionInterface ERROR: ensemble opt driver %s not found\n",
               ensembleOptDriver_);
        exit(1);
      }
      else fclose(fp);
      sprintf(command, "%s %s %s %d %d", ensembleOptDriver_,infile,outfile,
              ID, printLevel_);
    }

    system(command);   
  }
  else 
  {
    printf("WARNING: Output file %s exists before it is run.\n",outfile);
    printf("WARNING: PSUADE will use this output file.\n");
    printf("WARNING: If this is a mistake, stop PSUADE and clean up.\n");
    fclose(fIn);
  }

  while ((fIn=fopen(outfile, "r")) == NULL)
  {
    if (printLevel_ > 2)
    {
      printf("FunctionInterface: waiting for Job %d to complete.\n",ID);
#ifdef WINDOWS
      Sleep(1000 * launchInterval_);
#else
      sleep(launchInterval_);
#endif
    }
  }

  fIn = fopen(outfile, "r");
  if (fIn == NULL)
  {
    printf("FunctionInterface ERROR: output file %s not found.\n",outfile);
    exit(1);
  }
  outputCount = 0;
  for (ss = 0; ss < nSamp; ss++)
  {
    fgets(lineIn, 4000, fIn);
    for (ii = 0; ii < nOutputs_; ii++)
      sscanf(lineIn, "%lg", &(outputs[ss*nOutputs_+ii]));
    outputCount++;
    if (feof(fIn) == 1) break; 
  }
  fclose(fIn);

  if (outputCount != nSamp) 
  {
    printf("FunctionInterface ERROR: output file %s found but with\n",
           outfile);
    printf("                         insufficient data.\n");
    printf("Advice: Check the output format of your aux opt driver.\n");
    printf("        It should have %d lines each with %d output data.\n",
           nSamp, nOutputs_);
    exit(1);
  }

  if (strcmp(infile, "*")) unlink(infile);
  if (strcmp(outfile, "*")) unlink(outfile);
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
// Local function
// ------------------------------------------------------------------------
int FunctionInterface::psLocalFunction(int nInputs, double *inputs,
                                       int nOutputs, double *outputs)
{
  int    ss, ii;
  double X1, X2, X3, X4, D, W1, W2, W3, T1, T2, Y, D1, D2, D3, D4, W4;
  double alpha=5, beta=1, delta=-10, gamma3, gamma4;
  static int initializeProblem=2;

  //if (printLevel_ > 0)
  //   printf("FunctionInterface local function: Dowling's toy problem.\n");
  if (initializeProblem == -1)
  {
    printf("Available test problem :\n");
    printf(" 1. Dowling's toy problem - function in single-stage OUU.\n");
    printf(" 2. Dowling's toy problem - optimizer in two-stage OUU.\n");
    printf("Which problem to use ? (1 or 2) ");
    scanf("%d", &initializeProblem);
    if (initializeProblem != 2) initializeProblem = 1;
    if (initializeProblem == 1)
         printf("Dowling's toy function for single-stage OUU selected.\n");
    else printf("Dowling's toy optimizer for two-stage OUU selected.\n");
  }
  if (initializeProblem == 1)
  {
    D  = inputs[0];
    X1 = inputs[1];
    X2 = inputs[2];
    W1 = inputs[3];
    W2 = inputs[4];
    W3 = inputs[5];
    Y  = pow(X1 - D + W1, 2.0) + (1 + W1 * W1) * pow(X2 - D + W2, 2.0) +
         (W1 + W3) * X1 + (W2 + W3) * X2 + pow(W2*W2+(2+W3*W3)*D, 2.0);
    outputs[0] = Y;
  }
  else
  {
    if (nInputs != 12)
    {
      printf("psEnsembleLocalFunction ERROR: nInputs do not match.\n");
      exit(1);
    }
    D1 = inputs[0];
    D2 = inputs[1];
    D3 = inputs[2];
    D4 = inputs[3];
    X1 = inputs[4];
    X2 = inputs[5];
    X3 = inputs[6];
    X4 = inputs[7];
    W1 = inputs[8];
    W2 = inputs[9];
    W3 = inputs[10];
    W4 = inputs[11];
    X1 = - W1;
    gamma3 = 1 + W3 * W3;
    gamma4 = 1 + W4 * W4;
    X2 = - (delta + beta * D2 * gamma4 + beta * gamma4 * W2) / (beta*gamma4);
    X3 = - (delta * W3 + D3 * gamma3 * gamma4 + gamma3 * gamma4 * W3) / 
            (gamma3 * gamma4);
    T1 = delta * gamma3 + beta * delta *W3*W3 + beta * gamma3 * gamma4 * W2;
    T2 = beta * gamma3 * gamma4 * W3 * W3 + beta * D2 * gamma3 * gamma4;
    T1 = T1 + T2;
    T2 = beta * D4 * gamma3 * gamma4;
    T1 = T1 - T2;
    T2 = beta * delta * gamma3 * gamma4 + beta * D3 * gamma3 * gamma4 * W3;
    X4 = (T1 + T2) / (beta * gamma3 * gamma4 * gamma4);

    Y = (X1 + W1) * (X1 + W1);
    Y = Y + beta * pow(X2 + W2 + D2, 2.0);
    Y = Y + (1 + W3 * W3) * pow(X3 + D3 + W3, 2.0);
    T1 = pow(D4 + X2 + W3 * X3 + X4 * (1 + W4 * W4), 2.0);
    Y = Y + 1.0/(1+W4*W4) * T1;
    Y = Y - 2.0 * delta * X4;
    Y = Y + alpha * pow(D1+W1, 2.0);
    Y = Y + (10 - alpha) * D1 * D1;
    Y = Y + (10 - beta) * D2 * D2;
    Y = Y + 3 * D3 * D3;
    Y = Y + D4 * D4 * sqrt(1+W3*W3+W4*W4);
    outputs[0] = Y;
  }
  return 0;
#if 0
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
  Table[0]=1.7; Table[1] =1.8;Table[2] =1.9;Table[3] =1.9;Table[4] =2.1;
  Table[5] =2.2;
  Table[6]=1.18; Table[7] =1.3;Table[8] =1.3;Table[9] =1.3;Table[10]=1.4;
  Table[11]=1.6;
  Table[12]=0.95;Table[13]=1.0;Table[14]=1.0;Table[15]=1.1;Table[16]=1.2;
  Table[17]=1.3;
  Table[18]=0.67;Table[19]=0.7;Table[20]=0.7;Table[21]=0.8;Table[22]=0.8;
  Table[23]=0.9;
  Table[24]=0.58;Table[25]=0.6;Table[26]=0.6;Table[27]=0.7;Table[28]=0.7;
  Table[29]=0.8;

  Table[30]=2.48;Table[31]=2.6;Table[32]=2.9;Table[33]=3.1;Table[34]=3.5;
  Table[35]=3.9;
  Table[36]=1.67;Table[37]=1.7;Table[38]=2.0;Table[39]=2.0;Table[40]=2.4;
  Table[41]=2.7;
  Table[42]=1.33;Table[43]=1.4;Table[44]=1.6;Table[45]=1.6;Table[46]=1.9;
  Table[47]=2.1;
  Table[48]=0.92;Table[49]=0.9;Table[50]=1.1;Table[51]=1.1;Table[52]=1.3;
  Table[53]=1.5;
  Table[54]=0.77;Table[55]=0.8;Table[56]=0.9;Table[57]=0.9;Table[58]=1.1;
  Table[59]=1.3;

  Table[60]=3.58;Table[61]=4.0;Table[62]=4.6;Table[63]=5.1;Table[64]=6.0;
  Table[65]=7.0;
  Table[66]=2.28;Table[67]=2.5;Table[68]=3.0;Table[69]=3.3;Table[70]=4.0;
  Table[71]=4.6;
  Table[72]=1.75;Table[73]=2.0;Table[74]=2.3;Table[75]=2.6;Table[76]=3.1;
  Table[77]=3.8;
  Table[78]=1.14;Table[79]=1.3;Table[80]=1.5;Table[81]=1.7;Table[82]=2.0;
  Table[83]=2.3;
  Table[84]=0.95;Table[85]=1.1;Table[86]=1.3;Table[87]=1.4;Table[88]=1.7;
  Table[89]=1.9;

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
      muH2O = 1.002 * pow(10.0, 
       1.3272*(293.15-TEM-0.001053*(TEM-293.15)*(TEM-293.15))/(TEM-168.15));
      for (int jj = 0; jj < nA; jj++)
      {
        CO2 = As[jj];
        dtemp = (A * MEA + B) * TEM + C * MEA + D;
        dtemp = muH2O*exp(dtemp*(CO2*(E*MEA+F*TEM+G)+1)*MEA/(TEM*TEM));
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
// Local function
// ------------------------------------------------------------------------
int FunctionInterface::psEnsembleLocalFunction(int nSamples, int nInputs, 
                              double *inputs, int nOutputs, double *outputs)
{
  int    ss, ii;
  double X1, X2, X3, X4, D, W1, W2, W3, T1, T2, Y, D1, D2, D3, D4, W4;
  double alpha=5, beta=1, delta=-10, gamma3, gamma4;
  static int initializeProblem=2;

  if (printLevel_ > 0)
     printf("FunctionInterface local function: Dowling's toy problem.\n");
  if (initializeProblem == -1)
  {
    printf("Available test problem :\n");
    printf(" 1. Dowling's toy problem - function in single-stage OUU.\n");
    printf(" 2. Dowling's toy problem - optimizer in two-stage OUU.\n");
    printf("Which problem to use ? (1 or 2) ");
    scanf("%d", &initializeProblem);
    if (initializeProblem != 2) initializeProblem = 1;
    if (initializeProblem == 1)
         printf("Dowling's toy function for single-stage OUU selected.\n");
    else printf("Dowling's toy optimizer for two-stage OUU selected.\n");
  }
  if (initializeProblem == 1)
  {
    for (ss = 0; ss < nSamples; ss++)
    {
      D  = inputs[ss*nInputs];
      X1 = inputs[ss*nInputs+1];
      X2 = inputs[ss*nInputs+2];
      W1 = inputs[ss*nInputs+3];
      W2 = inputs[ss*nInputs+4];
      W3 = inputs[ss*nInputs+5];
      Y = pow(X1 - D + W1, 2.0) + (1 + W1 * W1) * pow(X2 - D + W2, 2.0) +
          (W1 + W3) * X1 + (W2 + W3) * X2 + pow(W2*W2+(2+W3*W3)*D, 2.0);
      outputs[ss] = Y;
    }
  }
  else
  {
    if (nInputs != 12)
    {
      printf("psEnsembleLocalFunction ERROR: nInputs do not match.\n");
      exit(1);
    }
    for (ss = 0; ss < nSamples; ss++)
    {
      D1  = inputs[ss*nInputs];
      D2 = inputs[ss*nInputs+1];
      D3 = inputs[ss*nInputs+2];
      D4 = inputs[ss*nInputs+3];
      X1 = inputs[ss*nInputs+4];
      X2 = inputs[ss*nInputs+5];
      X3 = inputs[ss*nInputs+6];
      X4 = inputs[ss*nInputs+7];
      W1 = inputs[ss*nInputs+8];
      W2 = inputs[ss*nInputs+9];
      W3 = inputs[ss*nInputs+10];
      W4 = inputs[ss*nInputs+11];
      X1 = - W1;
      gamma3 = 1 + W3 * W3;
      gamma4 = 1 + W4 * W4;
      X2 = - (delta + beta * D2 * gamma4 + beta * gamma4 * W2)/(beta*gamma4);
      X3 = - (delta * W3 + D3 * gamma3 * gamma4 + gamma3 * gamma4 * W3) / 
           (gamma3 * gamma4);
      T1 = delta * gamma3 + beta*delta*W3*W3 + beta * gamma3 * gamma4 * W2;
      T2 = beta * gamma3 * gamma4 * W3 * W3 + beta * D2 * gamma3 * gamma4;
      T1 = T1 + T2;
      T2 = beta * D4 * gamma3 * gamma4;
      T1 = T1 - T2;
      T2 = beta * delta * gamma3 * gamma4 + beta*D3 * gamma3 * gamma4 * W3;
      X4 = (T1 + T2) / (beta * gamma3 * gamma4 * gamma4);

      Y = (X1 + W1) * (X1 + W1);
      Y = Y + beta * pow(X2 + W2 + D2, 2.0);
      Y = Y + (1 + W3 * W3) * pow(X3 + D3 + W3, 2.0);
      T1 = pow(D4 + X2 + W3 * X3 + X4 * (1 + W4 * W4), 2.0);
      Y = Y + 1.0/(1+W4*W4) * T1;
      Y = Y - 2.0 * delta * X4;
      Y = Y + alpha * pow(D1+W1, 2.0);
      Y = Y + (10 - alpha) * D1 * D1;
      Y = Y + (10 - beta) * D2 * D2;
      Y = Y + 3 * D3 * D3;
      Y = Y + D4 * D4 * sqrt(1+W3*W3+W4*W4);
      outputs[ss] = Y;
    }
  }
  return 0;
}

// ************************************************************************
// create an FunctionInterface instantiation
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterface(PsuadeData *psuadeIO)
{
  int   nInputs, nOutputs;
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
  funcIO->loadFunctionData(pAppFiles.nStrings_, pAppFiles.strArray_);
  return funcIO;
}

// ************************************************************************
// create an FunctionInterface instantiation without setting up the
// driver (in case it is a response surface, this will save the time to
// set it up)
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterfaceSimplified(PsuadeData *psuadeIO)
{
  int   nInputs, nOutputs;
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
  funcIO->loadFunctionData(pAppFiles.nStrings_, pAppFiles.strArray_);
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


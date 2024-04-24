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
// Functions for the class PsuadeBase specifically for interactive mode
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "PsuadeBase.h"
#include "dtype.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"
#include "psMatrix.h"
#include "psVector.h"
#include "ProbMatrix.h"
#include "mData.h"

// ------------------------------------------------------------------------
// local includes : function approximator and others
// ------------------------------------------------------------------------
#include "FuncApprox.h"
#include "AnalysisManager.h"
#include "TSIAnalyzer.h"
#include "SobolAnalyzer.h"
#include "MCMCAnalyzer.h"
#include "PDFManager.h"
#include "Sampling.h"
#include "FunctionInterface.h"
#include "PsuadeData.h"
#include "Optimizer.h"
#include "PsuadeSession.h"
#include "PrintingTS.h"
#include "PDFHistogram.h"
#include "HMCMCAnalyzer.h"
#include "KSDensity.h"
#include "GP3.h"
#include "FactorialSampling.h"
#include "KPCA.h"
#include "KPCAInterface.h"

// ------------------------------------------------------------------------
// local defines 
// ------------------------------------------------------------------------
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::interpretInteractive()
{
  int    ss, ii, jj, kk, ll, it, ind, ind2, iInd, iInd1, iInd2, sInd, oInd;
  int    nPlots, iplot1, iplot2, iplot3, jplot=-1, oplot1, oplot2, iOne=1;
  int    status, inputID, outputID, count, scriptMode=0;
  int    nPtsPerDim=64, faType, faLeng, faFlag;
  double Ymax, Ymin, Xmin, Xmax, ddata, dtemp, thresh, threshL, threshU;
  char   command[5001], dataFile[5001], winput[50001], lineIn[50001];
  char   cString[10001],pString[50001],lineIn2[10001], *targv[8], **names;
  FILE   *fp, *scriptFp=NULL;
  FuncApprox **faPtrsRsEval=NULL;
  pData  pINames, pLower, pUpper, pONames, pInputs, pOutputs, pStates;
  pData  pPtr;
  psVector  vecInps,vecOuts,vecUBs,vecLBs,vecXT,vecYT,vecWT,vecVT;
  psIVector vecST, vecIT, vecTags;

  //**/ ----------------------------------------------------------------
  // Launch the command interpreter
  //**/ ----------------------------------------------------------------
  printf("PSUADE - A Problem Solving environment for \n");
  printf("         Uncertainty Analysis and Design Exploration (%d.%d.%d)\n",
         psuade_VERSION_MAJOR, psuade_VERSION_MINOR, psuade_VERSION_PATCH);
  printf("(for help, enter <help>)\n");
  printEquals(PL_INFO, 0);
  //**/ allow only a limited number of contiguous invalid commands
  int commandCnt = 0;
  //**/ command status that can be obtained by using getstatus
  int cmdStatus = 0;
  //**/ capture information of the current session
  PsuadeSession *currSession = new PsuadeSession();

  //**/ continue to receive and interpret command until 'quit'
  while (1)
  {
    //**/ =============================================================
    //**/ permit users to run a script in command line mode
    //**/ Users enter the name of a script, and psuade
    //**/ searches to see if the 'script' exists and
    //**/ if so, it will execute the commands in the script
    //**/ =============================================================
    if (scriptMode == 1 && scriptFp != NULL)
    {
      fgets(lineIn, 5000, scriptFp);
      if (feof(scriptFp) != 0)
      {
        fclose(scriptFp);
        scriptMode = 0;
      }
      for (ii = 0; ii < 100; ii++) command[ii] = '\0';
      sscanf(lineIn, "%s", command);
      printf("script> %s\n", command);
    }
    else
    {
      printf("psuade> ");
      for (ii = 0; ii < 500; ii++) lineIn[ii] = '\0';
      fgets(lineIn,5000,stdin); 
      winput[0] = '\0';
      pString[0] = '\0';
      command[0] = '\0';
      sscanf(lineIn, "%s", command);
    }
    //**/ =============================================================
    //**/ this prevents an infinite loop (which can happen when psuade
    //**/ is executing a script with `quit' as the last line) that may 
    //**/ create a very large file of garbage
    //**/ =============================================================
    if (!strcmp(command, "\0"))
    {
      commandCnt++;
      if (commandCnt >= 10)
      {
        printf("Enter carriage return > 10 times ==> terminate.\n");
        return 0;
      }
    }
    else commandCnt = 0;

    //**/ =============================================================
    // +++ help menus
    //**// ============================================================
    if (!strcmp(command, "help") || !strcmp(command, "h"))
    {
      strcpy(winput, "\0");
      sscanf(lineIn, "%s %s %s", command, winput, pString);
      displayHelpMenu(winput, pString);
      cmdStatus = 0;
    }

    //**/ =============================================================
    // +++ rename a file (especially rename psuadeData)
    //**// ============================================================
    else if (!strcmp(command, "rename"))
    {
      sscanf(lineIn,"%s %s %s",command,winput,cString);
      if (!strcmp(winput, "-h"))
      {
        printf("rename: renames a file in the current directory.\n");
        printf("Syntax: rename <file1> <file2>\n");
        printf("where <file1> and <file2> are source/destination files.\n");
        printf("\nThis command enables file name change without ");
        printf("leaving the PSUADE\n");
        printf("command line mode.\n");
        continue;
      }
      fp = fopen(winput,"r");
      if (fp == NULL) 
      {
        printf("ERROR: file %s not found.\n", winput);
        cmdStatus = 1;
        continue;
      }
      fclose(fp);
      if (rename(winput, cString) == 0) 
      {
        printf("%s has been renamed to %s\n", winput, cString);
        cmdStatus = 0;
      }
      else 
      {
        printf("ERROR: in renaming %s to %s\n", winput, cString);
        cmdStatus = 1;
      }
    }

    //**/ =============================================================
    // +++ run a script
    //**/ =============================================================
    else if (!strcmp(command, "run"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("run: run a PSUADE input file in batch mode.\n");
        printf("Syntax: run <file>\n");
        printf("where <file> is a PSUADE input file.\n");
        printf("This command allows you to run PSUADE batch files ");
        printf("in command line mode.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn,"%s %s",command,dataFile);
      kk = strlen(dataFile);
      //**/ if there is no argument, run from resident setup
      if (kk <= 0)
      {
        run();
        cmdStatus = 0;
        continue;
      }
      //**/ if there is argument, check whether the file exists
      if ((fp=fopen(dataFile,"r")) == NULL)
      {
        printf("file %s not found.\n", dataFile);
        printf("Syntax: run <file>.\n");
        printf("where <file> is a PSUADE input file.\n");
        cmdStatus = 1;
        continue;
      }
      else fclose(fp);
      //**/ if file exists, get sample information from file
      status = getInputFromFile(dataFile);

      //**/ if everything is fine, run the sample
      if (status != 0) 
      { 
        printf("ERROR: run file not valid : check file format.\n");
        cmdStatus = 1;
        continue;
      }
      run();
      cmdStatus = 0;
    }

    //**/ =============================================================
    // RS commands that are in another file (PsuadeRSInterpreter.cpp)
    //**/ =============================================================
    else if (!strcmp(command,"rsua")  || !strcmp(command,"rsua2") ||
         !strcmp(command,"rsuab2") || !strcmp(command,"rs_ua2") ||
         !strcmp(command,"rsuab2") || !strcmp(command,"rsuab") ||
         !strcmp(command,"rsuap")  || !strcmp(command,"rssobol1") || 
         !strcmp(command,"rssobol2") || !strcmp(command,"rssoboltsi") || 
         !strcmp(command,"rssobolg") || !strcmp(command,"rs_qsa") ||
         !strcmp(command,"rs1") || !strcmp(command,"rs1s") ||
         !strcmp(command,"rs2") || !strcmp(command,"rs3") ||
         !strcmp(command,"rs3m") || !strcmp(command,"rs4") ||
         !strcmp(command,"rssd") || !strcmp(command,"rssd_ua") ||
         !strcmp(command,"rsi2") || !strcmp(command,"rsi3") ||
         !strcmp(command,"rsi3m") || !strcmp(command,"rsmeb") || 
         !strcmp(command,"rsieb") || !strcmp(command,"rssobol1b") || 
         !strcmp(command,"rssobol2b") || !strcmp(command,"rssoboltsib"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        cmdStatus = RSBasedAnalysis(lineIn, NULL);
        continue;
      }

      //**/ create a Psuade session and call RS functions
      PsuadeSession *newSession = new PsuadeSession();
      newSession->outputLevel_ = outputLevel_;
      newSession->nInputs_ = nInputs_;
      newSession->nOutputs_ = nOutputs_;
      newSession->nSamples_ = nSamples_;
      if (nInputs_ > 0 && nSamples_ > 0)
      {
        newSession->vecSamInputs_.load(nInputs_*nSamples_, 
                                       VecSamInputs_.getDVector());
      }
      if (nOutputs_ > 0 && nSamples_ > 0)
      {
        newSession->vecSamOutputs_.load(nOutputs_*nSamples_, 
                                        VecSamOutputs_.getDVector());
        newSession->vecSamStates_.load(nSamples_,
                                       VecSamStates_.getIVector());
      }
      if (nInputs_ > 0)
      {
        newSession->vecInpLBounds_.load(nInputs_,
                                        VecILowerBs_.getDVector());
        newSession->vecInpUBounds_.load(nInputs_,
                                        VecIUpperBs_.getDVector());
        newSession->inputNames_  = StrInpNames_;
        newSession->outputNames_ = StrOutNames_;
      }
      newSession->psuadeIO_ = psuadeIO_;
      cmdStatus = RSBasedAnalysis(lineIn, newSession);
      delete newSession;
      newSession = NULL;
    }

    //**/ =============================================================
    // Input/output commands
    //**/ =============================================================
    //**/ =============================================================
    //**/ -------------------------------------------------------------
    // +++ clear workspace +++
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "clear"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      { 
        printf("clear: clear PSUADE's workspace.\n");
        printf("Syntax: clear <filename>.\n");
        continue;
      }
      cleanUp();
      if (currSession != NULL) delete currSession;
      currSession = NULL;
      cmdStatus = 0;
      printf("INFO: PSUADE's workspace is now empty.\n");
      fflush(stdout);
    }

    //**/ -------------------------------------------------------------
    // +++ load/read and loadp +++
    //**/ load sample from a psuade file. 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "load") || !strcmp(command, "loadp") ||
             !strcmp(command, "read"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      { 
        printf("load/read: load PSUADE data from a file to workspace\n");
        printf("Syntax: load (or read) <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      if (!strcmp(command,"loadp")) 
           strcpy(dataFile, "psuadeData");
      else sscanf(lineIn,"%s %s",command,dataFile);
      if (!strcmp(dataFile,psOutputFilename_))
      {
        printAsterisks(PL_INFO, 0);
        printf("INFO: You are loading a file with the same name ");
        printf("as the default output\n");
        printf("      file %s, which is periodically overwritten ",
               psOutputFilename_);
        printf("during some\n");
        printf("      of PSUADE's internal operations. As such, ");
        printf("you are advised to\n");
        printf("      rename your 'psuadeData' file to something else (");
        printf("You can use\n");
        printf("      the output_file command to change the default ");
        printf("output filename.)\n");
        printAsterisks(PL_INFO, 0);
      }
      if (!checkFileExists(dataFile))
      {
        printf("Syntax: load (or read) <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        cmdStatus = 1;
        continue;
      }
      cleanUp();
      //**/ probably should not be reset, since users may want to
      //**/ reload but keep the settings
      //psConfig_.reset();
      psuadeIO_ = new PsuadeData();
      psuadeIO_->setOutputLevel(0);
      status = psuadeIO_->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
               dataFile);
        cleanUp();
        cmdStatus = 1;
        continue;
      }
      getSampleFromPsuadeIO();
      if (VecSamInputs_.length() <= 0 || VecSamOutputs_.length() <= 0)
      {
        printf("WARNING: no sample matrix or output found.\n");
        nSamples_ = 0;
      }
      else
      {
        printf("load complete : nSamples = %d\n", nSamples_);
        printf("                nInputs  = %d\n", nInputs_);
        printf("                nOutputs = %d\n", nOutputs_);
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
      }
      currSession->psuadeIO_ = psuadeIO_;
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      outputLevel_ = pPtr.intData_;
      cmdStatus = 0;

      fflush(stdout);
      checkResidentSampleOutputs(); 
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ loadmore  +++
    //**/ load sample from a file and add to the existing data set
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "loadmore") || 
             !strcmp(command, "readmore")) 
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      { 
        printf("loadmore: ADD another PSUADE sample to workspace\n");
        printf("          (it is equivalent to sample concatenation)\n");
        printf("Syntax: loadmore <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psOutputFilename_))
      {
        printf("WARNING: you are loading a file with the same name ");
        printf("as the default\n");
        printf("         output file %s\n",psOutputFilename_);
        printf("This default file may be overwritten during a ");
        printf("PSUADE session so you\n");
        printf("may want to change the default output filename ");
        printf("with the 'output_file'\n");
        printf("command.\n");
      }
      if (!checkFileExists(dataFile))
      {
        printf("ERROR: file %s not found.\n", dataFile);
        printf("Syntax: loadmore <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        cmdStatus = 1;
        continue;
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->setOutputLevel(0);
      status = ioPtr->readPsuadeFile(dataFile);
      int flag = 0;
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
               dataFile);
        cmdStatus = 1;
        continue;
      }
      //**/ get sample input information
      //**/ if current sample is not empty, append new sample
      if (VecSamInputs_.length() == 0)
      {
        printf("ERROR: Nothing has been loaded. Use 'load' instead\n");
        cmdStatus = 1;
        continue;
      }
      //**/ check nInputs consistency
      ioPtr->getParameter("input_ninputs", pPtr);
      int nInps = pPtr.intData_;
      if (nInps != nInputs_)
      {
        printf("ERROR: nInputs are different.\n");
        printf("       incoming nInputs = %d\n", ind);
        printf("       expected nInputs = %d\n", nInputs_);
        cmdStatus = 1;
        continue;
      }
      //**/ check nOutputs consistency
      ioPtr->getParameter("output_noutputs",pPtr);
      int nOuts = pPtr.intData_;
      if (nOuts != nOutputs_)
      {
        printf("ERROR: nOutputs are different.\n");
        printf("       incoming nOutputs = %d\n", ind);
        printf("       expected nOutputs = %d\n", nOutputs_);
        printf("INFO: local data set not changed.\n");
        cmdStatus = 1;
        continue;
      }
      //**/ get input bounds
      ioPtr->getParameter("input_lbounds", pLower);
      ioPtr->getParameter("input_ubounds", pUpper);
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (pLower.dbleArray_[ii] < VecILowerBs_[ii]) 
          VecILowerBs_[ii] = pLower.dbleArray_[ii];
        if (pUpper.dbleArray_[ii] > VecIUpperBs_[ii]) 
          VecIUpperBs_[ii] = pUpper.dbleArray_[ii];
      }
      names = NULL;
      //**/ get sample 
      ioPtr->getParameter("method_nsamples", pPtr);
      int nSamp = pPtr.intData_;
      vecXT = VecSamInputs_;
      vecYT = VecSamOutputs_;
      vecST = VecSamStates_;
      VecSamInputs_.setLength((nSamples_+nSamp)*nInputs_);
      VecSamOutputs_.setLength((nSamples_+nSamp)*nOutputs_);
      VecSamStates_.setLength(nSamples_+nSamp);
      for (ii = 0; ii < nSamples_*nInputs_; ii++)
        VecSamInputs_[ii] = vecXT[ii];
      for (ii = 0; ii < nSamples_*nOutputs_; ii++)
        VecSamOutputs_[ii] = vecYT[ii];
      for (ii = 0; ii < nSamples_; ii++)
        VecSamStates_[ii] = vecST[ii];
      ioPtr->getParameter("input_sample", pPtr);
      for (ii = 0; ii < nSamp*nInputs_; ii++)
        VecSamInputs_[nSamples_*nInputs_+ii] = pPtr.dbleArray_[ii];
      ioPtr->getParameter("output_sample", pPtr);
      for (ii = 0; ii < nSamp*nOutputs_; ii++)
        VecSamOutputs_[nSamples_*nOutputs_+ii] = pPtr.dbleArray_[ii];
      ioPtr->getParameter("output_states", pPtr);
      for (ii = 0; ii < nSamp; ii++)
        VecSamStates_[nSamples_+ii] = pPtr.intArray_[ii];

      //**/ not needed any more (5/19/08)
      //**/ PDFTransform(psuadeIO_,nSamp,nInputs_,VecSamInputs_);
      nSamples_ += nSamp;
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                    NULL,NULL,VecSamInputs_.getDVector(),NULL,NULL,
                    NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                         VecSamOutputs_.getDVector(), 
                         VecSamStates_.getIVector(), NULL); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("loadmore complete : nSamples = %d\n", nSamples_);
      printf("                    nInputs  = %d\n", nInputs_);
      printf("                    nOutputs = %d\n", nOutputs_);
      cmdStatus = 0;
      fflush(stdout);
      VecDataReg_.clean();
      delete ioPtr;
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ write +++
    //**/ write the data to an output file - this is needed if PDF ~=U
    //**/ and users want to have the actual input data
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      { 
        printf("write: save the resident sample to ");
        printf("another PSUADE data file.\n");
        printf("Syntax: write <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, psOutputFilename_);
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n");
        printf("          name %s.\n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: no file to write to.\n");
        printf("Syntax: write <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        cmdStatus = 1;
      }
      else
      {
        if (nOutputs_ > 1)
        {
          sprintf(pString,
            "Save only one output (n - save all outputs)? (y or n) "); 
          getString(pString, winput);
        }
        else
        {
          winput[0] = 'n';
          outputID = 0;
        }
        if (winput[0] == 'y')
        {
          sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
          outputID = getInt(1, nOutputs_, pString);
          outputID--;
          for (sInd = 0; sInd < nSamples_; sInd++)
            VecSamOutputs_[sInd] = 
                   VecSamOutputs_[sInd*nOutputs_+outputID];
          pONames.clean();
          psuadeIO_->getParameter("output_names", pONames);
          names = pONames.strArray_;
          strcpy(names[0], names[outputID]);
          nOutputs_ = 1;
          psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
        }
        else names = NULL;
        updateSampleInputBounds();
        psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                     VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                     VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                     VecInpPDFs_.getIVector(),VecInpMeans_.getDVector(),
                     VecInpStds_.getDVector(),inputCMat_);
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                         VecSamOutputs_.getDVector(), 
                         VecSamStates_.getIVector(), names);
        if (winput[0] == 'y') pONames.clean();
        psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,
                                       nReplications_, -1,-1);
        psuadeIO_->writePsuadeFile(dataFile,0);
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ nwrite +++ 
    //**/ write the unevaluated sample points data to an output file 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "nwrite"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("nwrite: save the unevaluated points to ");
        printf("a PSUADE data file\n");
        printf("NOTE: sample points with undefined outputs or ");
        printf("status=0 are considered\n");
        printf("      unevaluated points.\n");
        printf("Syntax: nwrite <filename>.\n");
        printf("where <filename> is a PSUADE data file.\n");
        continue;
      }

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, psOutputFilename_);
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n");
        printf("          name %s.\n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        printf("Syntax: nwrite <filename>.\n");
        cmdStatus = 1;
      }
      else
      {
        vecXT.setLength(nInputs_*nSamples_);
        vecYT.setLength(nOutputs_*nSamples_);
        vecST.setLength(nSamples_);
        count = 0;
        for (ii = 0; ii < nSamples_; ii++)
        {
          for (jj = 0; jj < nOutputs_; jj++)
            if (VecSamOutputs_[ii*nOutputs_+jj] == PSUADE_UNDEFINED)
              VecSamStates_[ii] = 0;
          if (VecSamStates_[ii] == 0)
          {
            for (jj = 0; jj < nInputs_; jj++)
              vecXT[count*nInputs_+jj] = VecSamInputs_[ii*nInputs_+jj];
            for (jj = 0; jj < nOutputs_; jj++)
              vecYT[count*nOutputs_+jj] = 
                      VecSamOutputs_[ii*nOutputs_+jj]; 
            vecST[count++] = VecSamStates_[ii]; 
          }
        }
        if (count == 0)
        {
          printf("INFO: no unevaluated sample points\n");
          printf("      ==> no file generated.\n");
          cmdStatus = 1;
        }
        else
        {
          PsuadeData *ioPtr = new PsuadeData();
          ioPtr->updateInputSection(count, nInputs_, NULL, 
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   vecXT.getDVector(), StrInpNames_.getStrings(), 
                   NULL,NULL,NULL,NULL); 
          ioPtr->updateOutputSection(count,nOutputs_,vecYT.getDVector(),
                       vecST.getIVector(), StrOutNames_.getStrings());
          ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
          ioPtr->writePsuadeFile(dataFile, 0);
          delete ioPtr;
          printf("INFO: unevaluated sample has been written to file %s.\n", 
                 dataFile);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ iread +++
    //**/ read the data from an input file in input only format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iread"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iread: read sample inputs from a file written ");
        printf("using iwrite or a file\n");
        printf("       with the following format:\n");
        printf("line 1: PSUADE_BEGIN (optional)\n");
        printf("line 2: <nSamples> <nInputs>\n");
        printf("line 3: (optional) '#' in column 1, then input names\n");
        printf("line 4: 1 <sample point 1 inputs> \n");
        printf("line 5: 2 <sample point 2 inputs> \n");
        printf(".......\n");
        printf("line n: PSUADE_END (optional)\n\n");
        printf("Syntax: iread <filename>.\n");
        printf("where <filename> is the name of the data file:\n\n");
        printf("NOTE: iread can be used to read a posterior sample ");
        printf("(MCMCPostSample)\n");
        printf("created by rsmcmc.\n");
        continue;
      }
      //**/ check file
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!checkFileExists(dataFile))
      {
        printf("Syntax: iread <filename>.\n");
        printf("where <filename> is the name of the data file:\n");
        cmdStatus = 1;
        continue;
      }
      cleanUp();
      //**/ open file and check key words
      psuadeIO_ = new PsuadeData();
      psuadeIO_->setOutputLevel(0);
      fp = fopen(dataFile,"r");
      fscanf(fp, "%s", winput);
      if (strcmp(winput, "PSUADE_BEGIN"))
      {
        printf("INFO: Optional PSUADE_BEGIN keyword is not found.\n");
        printf("      Make sure the first line has nSamples ");
        printf("and nInputs.\n");
        fclose(fp);
        fp = fopen(dataFile,"r");
      }
      else fgets(lineIn, 5000, fp);
      fscanf(fp, "%d %d", &nSamples_, &nInputs_);
      nOutputs_ = 1;
      if (nSamples_ <= 0 || nInputs_ <= 0)
      {
        printf("iread ERROR: some sample parameters <= 0.\n");
        printf("      nSamples read = %d\n", nSamples_);
        printf("      nInputs  read = %d\n", nInputs_);
        printf("NOTE: the format for the second line should be: \n");
        printf("  <nSamples> <nInputs> \n");
        delete psuadeIO_;
        psuadeIO_ = NULL;
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      //**/ detect comment line
      int fileHasInpNames = 0;
      fgets(pString, 5000, fp);
      while (1)
      {
        kk = getc(fp);
        if (kk == '#')
        {
          fileHasInpNames = 1;
          StrInpNames_.setNumStrings(nInputs_);
          for (ii = 0; ii < nInputs_; ii++)
          {
            fscanf(fp,"%s", pString);
            StrInpNames_.loadOneString(ii, pString);
          }
          fgets(pString, 5000, fp);
        }
        else
        {
          ungetc(kk, fp);
          break;
        }
      }
      //**/ read in sample data
      VecSamInputs_.setLength(nSamples_*nInputs_);
      VecSamOutputs_.setLength(nSamples_*nOutputs_);
      VecSamStates_.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++) VecSamStates_[ii] = 1;
      for (ii = 0; ii < nSamples_; ii++)
      {
        fscanf(fp, "%d", &jj);
        if ((ii+1) != jj)
        {
          printf("iread ERROR: sample index mismatch.\n");
          printf("         sample number read     = %d\n", jj);
          printf("         sample number expected = %d\n", ii+1);
          printf("INFO: the first field must be the sample number.\n");
          cmdStatus = 1;
          break;
        }
        for (jj = 0; jj < nInputs_; jj++)
        {
          fscanf(fp,"%lg", &ddata);
          VecSamInputs_[ii*nInputs_+jj] = ddata;
        }
        VecSamOutputs_[ii] = PSUADE_UNDEFINED;
        VecSamStates_[ii] = 0;
      }
      if (ii != nSamples_) 
      {
        VecSamInputs_.clean();
        VecSamOutputs_.clean();
        VecSamStates_.clean();
        delete psuadeIO_;
        psuadeIO_ = NULL;
        StrInpNames_.clean();
        fclose(fp);
        continue;
      }
      fclose(fp);

      //**/ update INPUT section 
      StrInpNames_.setNumStrings(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString, "X%d", ii+1);
        StrInpNames_.loadOneString(ii, pString);
      }
      VecILowerBs_.setLength(nInputs_);
      VecIUpperBs_.setLength(nInputs_);
      psVector VecTW;
      VecTW.setLength(nSamples_);
      double ireadMax, ireadMin;
      for (ii = 0; ii < nInputs_; ii++)
      {
        for (jj = 0; jj < nSamples_; jj++)
          VecTW[jj] = VecSamInputs_[jj*nInputs_+ii];
        ireadMax = ireadMin = VecTW[0];
        for (jj = 1; jj < nSamples_; jj++)
        {
          if (VecTW[jj] > ireadMax) ireadMax = VecTW[jj];
          if (VecTW[jj] < ireadMin) ireadMin = VecTW[jj];
        } 
        VecILowerBs_[ii] = ireadMin;
        VecIUpperBs_[ii] = ireadMax;
        if (VecILowerBs_[ii] == VecIUpperBs_[ii]) 
          VecIUpperBs_[ii] += 1.0e-12;
      }
      VecInpPDFs_.setLength(nInputs_);
      VecInpMeans_.setLength(nInputs_);
      VecInpStds_.setLength(nInputs_);
      if (inputCMat_ != NULL) delete inputCMat_;
      inputCMat_ = new psMatrix();
      inputCMat_->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        VecInpPDFs_[ii] = 0;
        VecInpMeans_[ii] = 0;
        VecInpStds_[ii] = 0;
        inputCMat_->setEntry(ii,ii,1.0);
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                  VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                  VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                  VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                  VecInpStds_.getDVector(), inputCMat_); 

      //**/ update OUTPUT section 
      StrOutNames_.setNumStrings(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        sprintf(pString, "Y%d", ii+1);
        StrOutNames_.loadOneString(ii, pString);
      }
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(),
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, nSamples_, 1, 0, 0);
      psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
      nReplications_ = 1;
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("iread: data have been read.\n");
      printf("nSamples = %d\n", nSamples_);
      printf("nInputs  = %d\n", nInputs_);
      printf("nOutputs = %d (set to 1)\n", nOutputs_);
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iwrite +++ 
    //**/ write the data to an output file in input only format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iwrite") || !strcmp(command, "iwrite2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iwrite: save only the sample inputs to an ASCII file\n");
        printf("Syntax: iwrite2 <filename>  OR\n");
        printf("        iwrite  <filename>  (PSUADE_BEGIN not written)\n");
        printf("where <filename> is the name of the target data file.\n");
        printf("\nThe target file will have the following format: \n");
        printf("line 1: PSUADE_BEGIN\n");
        printf("line 2: <nSamples> <nInputs>\n");
        printf("line 3: # input parameter names\n");
        printf("line 4: 1 <sample point 1 inputs> \n");
        printf("line 5: 2 <sample point 2 inputs> \n");
        printf(".......\n");
        printf("line n: PSUADE_END\n\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, psOutputFilename_);
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n"); 
        printf("          name %s. \n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
        cmdStatus = 1;
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        cmdStatus = 1;
      }
      else
      {
        if (VecSamInputs_.length() == 0) 
        {
          printf("ERROR: no input to write.\n");
          cmdStatus = 1;
        }
        else
        {
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open file %s.\n", dataFile);
            cmdStatus = 1;
            continue;
          }
          if (!strcmp(command, "iwrite2")) fprintf(fp, "PSUADE_BEGIN\n");
          fprintf(fp, "%d %d\n", nSamples_, nInputs_);
          for (ii = 0; ii < nSamples_; ii++)
          {
            fprintf(fp, "%d ", ii+1);
            for (jj = 0; jj < nInputs_; jj++)
              fprintf(fp,"%24.16e ",VecSamInputs_[ii*nInputs_+jj]); 
            fprintf(fp, "\n");
          }
          if (!strcmp(command, "iwrite2"))
          {
            fprintf(fp, "PSUADE_END\n");
            fprintf(fp, "# ");
            for (ii = 0; ii < nInputs_; ii++)
              fprintf(fp, "%s ", StrInpNames_[ii]);
          }
          fclose(fp);
          printf("iwrite: sample inputs has been written to %s.\n",
                 dataFile);
          cmdStatus = 0;
        }
      }
    }
    
    //**/ -------------------------------------------------------------
    // +++ owrite +++
    //**/ write the data to an output file in output only format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "owrite") || !strcmp(command, "owrite2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("owrite: save only the sample outputs to an ASCII file\n");
        printf("Syntax: owrite2 <filename>  OR\n");
        printf("        owrite  <filename>  (PSUADE_BEGIN not written)\n");
        printf("where <filename> is the name of the target data file.\n");
        printf("\nThe target file will have the following format: \n");
        printf("line 1: PSUADE_BEGIN\n");
        printf("line 2: <nSamples> <nOutputs>\n");
        printf("line 3: 1 <sample point 1 outputs> \n");
        printf("line 4: 2 <sample point 2 outputs> \n");
        printf(".......\n");
        printf("line n: PSUADE_END\n\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        cmdStatus = 1;
      }
      else
      {
        if (VecSamOutputs_.length() == 0)
        {
          printf("ERROR: No output to write.\n");
          cmdStatus = 1;
        }
        else
        {
          if (nOutputs_ == 1) outputID = 0;
          else
          {
            sprintf(pString, "Output number (1 - %d, 0 for all): ",
                    nOutputs_);
            outputID = getInt(0, nOutputs_, pString);
            outputID--;
          }
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open file %s.\n", dataFile);
            cmdStatus = 1;
            continue;
          }
          if (!strcmp(command, "owrite2")) fprintf(fp,"PSUADE_BEGIN\n");
          if (outputID == -1)
               fprintf(fp, "%d %d\n", nSamples_, nOutputs_);
          else fprintf(fp, "%d 1\n", nSamples_);
          for (ii = 0; ii < nSamples_; ii++)
          {
            fprintf(fp, "%d ", ii+1);
            if (outputID == -1)
            {
              for (jj = 0; jj < nOutputs_; jj++)
                fprintf(fp,"%24.16e ",VecSamOutputs_[ii*nOutputs_+jj]); 
              fprintf(fp, "\n");
            }
            else
              fprintf(fp,"%24.16e\n",
                      VecSamOutputs_[ii*nOutputs_+outputID]);
          }
          if (!strcmp(command, "owrite2")) fprintf(fp,"PSUADE_END\n");
          fclose(fp);
          printf("owrite: sample outputs has been written to %s.\n",
                 dataFile);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ read_std +++ 
    //**/ read data in standard format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "read_std"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("read_std: read a data file in the following format:\n");
        printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
        printf("line 2: <sample 1 inputs> < sample 1 outputs>\n");
        printf("line 3: <sample 2 inputs> < sample 2 outputs>\n");
        printf(".......\n\n");
        printf("Syntax: read_std <filename>.\n");
        printf("where <filename> is a data file in standard format.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: please specify a file to read from.\n");
        printf("Syntax: read_std <filename>.\n");
        printf("where <filename> is a data file in standard format.\n");
        cmdStatus = 1;
        continue;
      }
      if (!checkFileExists(dataFile))
      {
        printf("Syntax: read_std <filename>.\n");
        printf("where <filename> is a data file in standard format.\n");
        cmdStatus = 1;
        continue;
      }
      cleanUp();
      fp = fopen(dataFile, "r");
      kk = '#';
      while (kk == '#')
      {
        kk = getc(fp);
        if (kk != '#') 
        {
          ungetc (kk, fp);
          break;
        }
        else fgets(pString, 50000, fp);
      }
      fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &nOutputs_);
      if (nSamples_ <= 0) 
      {
        printf("ERROR: nSamples <= 0 (%d)\n", nSamples_);
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ <= 0) 
      {
        printf("ERROR: nInputs <= 0 (%d)\n", nInputs_);
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      int flag = 0;
      if (nOutputs_ <= 0) 
      {
        printf("INFO: nOutputs = 0 ==> create a dummy output\n");
        nOutputs_ = 1;
        flag = 1;
      }
      printf("nSamples = %d\n", nSamples_);
      printf("nInputs  = %d\n", nInputs_);
      printf("nOutputs = %d\n", nOutputs_);
      VecSamInputs_.setLength(nSamples_*nInputs_);
      VecSamOutputs_.setLength(nSamples_*nOutputs_);
      VecSamStates_.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
        {
          fscanf(fp, "%lg", &ddata);
          VecSamInputs_[ii*nInputs_+jj] = ddata; 
        }
        if (flag == 0)
        {
          for (jj = 0; jj < nOutputs_; jj++)
          {
            fscanf(fp,"%lg", &ddata);
            VecSamOutputs_[ii*nOutputs_+jj] = ddata;
          }
        }
        else VecSamOutputs_[ii] = PSUADE_UNDEFINED;
        VecSamStates_[ii] = 1;
        for (jj = 0; jj < nOutputs_; jj++)
          if (VecSamOutputs_[ii*nOutputs_+jj] > 0.5*PSUADE_UNDEFINED)
            VecSamStates_[ii] = 0;
      }
      fclose(fp);
      StrInpNames_.setNumStrings(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString, "X%d", ii+1);
        StrInpNames_.loadOneString(ii, pString);
      }
      StrOutNames_.setNumStrings(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        sprintf(pString, "Y%d", ii+1);
        StrOutNames_.loadOneString(ii, pString);
      }
      VecILowerBs_.setLength(nInputs_);
      VecIUpperBs_.setLength(nInputs_);
      for (jj = 0; jj < nInputs_; jj++)
      {
        VecILowerBs_[jj] = VecSamInputs_[jj];
        VecIUpperBs_[jj] = VecSamInputs_[jj];
        for (ii = 1; ii < nSamples_; ii++)
        {
          if (VecSamInputs_[ii*nInputs_+jj] < VecILowerBs_[jj])
            VecILowerBs_[jj] = VecSamInputs_[ii*nInputs_+jj];
          if (VecSamInputs_[ii*nInputs_+jj] > VecIUpperBs_[jj])
            VecIUpperBs_[jj] = VecSamInputs_[ii*nInputs_+jj];
        }
        if (VecILowerBs_[jj] == VecIUpperBs_[jj])
          VecIUpperBs_[jj] = VecILowerBs_[jj] + 1e-12; 
      }
      SamplingMethod_ = PSUADE_SAMP_MC;
      VecInpPDFs_.setLength(nInputs_);
      VecInpMeans_.setLength(nInputs_);
      VecInpStds_.setLength(nInputs_);
      if (inputCMat_ != NULL) delete inputCMat_;
      inputCMat_  = new psMatrix();
      inputCMat_->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        VecInpPDFs_[ii] = 0;
        VecInpMeans_[ii] = 0;
        VecInpStds_[ii] = 0;
        inputCMat_->setEntry(ii,ii,1.0);
      }
      psuadeIO_ = new PsuadeData();
      psuadeIO_->setOutputLevel(0);
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                  VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                  VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                  VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                  VecInpStds_.getDVector(), inputCMat_); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                  VecSamOutputs_.getDVector(), 
                  VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      nReplications_ = 1;
      psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,
                                     nReplications_,0,0);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
      checkResidentSampleOutputs(); 
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ write_std +++ 
    //**/ write the data to an output file in standard format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write_std"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("write_std: write a data file in the standard format.\n");
        printf("           (run <read_std -h> to see the standard ");
        printf("format.)\n");
        printf("Syntax: write_std <filename>.\n");
        printf("where <filename> is the name of the target data file.\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n"); 
        printf("          name %s. \n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
        cmdStatus = 1;
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        printf("Syntax: write_std <filename>.\n");
        cmdStatus = 1;
      }
      else
      {
        if (VecSamInputs_.length() == 0) 
        {
          printf("ERROR: no inputs to write.\n");
          cmdStatus = 1;
        }
        else
        {
          if (VecSamOutputs_.length() == 0) outputID = -9;
          else
          {
            if (nOutputs_ == 1) outputID = 0;
             else
             {
               sprintf(pString,"Enter output ID (1 - %d, 0 for all) : ",
                       nOutputs_);
               outputID = getInt(0, nOutputs_, pString);
               outputID--;
             }
          }
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open file %s.\n", dataFile);
            continue;
          }
          if (outputID == -9) 
            fprintf(fp, "%d %d\n", nSamples_, nInputs_);
          else if (outputID == -1) 
            fprintf(fp, "%d %d %d\n", nSamples_, nInputs_, nOutputs_);
          else fprintf(fp, "%d %d 1\n", nSamples_, nInputs_);
          for (ii = 0; ii < nSamples_; ii++)
          {
            for (jj = 0; jj < nInputs_; jj++)
              fprintf(fp, "%24.16e ", VecSamInputs_[ii*nInputs_+jj]); 
            if (outputID == -9) fprintf(fp, "\n");
            else if (outputID == -1)
            {
              for (jj = 0; jj < nOutputs_; jj++)
                fprintf(fp,"%24.16e ",VecSamOutputs_[ii*nOutputs_+jj]);
              fprintf(fp,"\n");
            }
            else 
            {
              fprintf(fp,"%24.16e\n",
                      VecSamOutputs_[ii*nOutputs_+outputID]);
            }
          }
          fclose(fp);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ read_csv +++ 
    //**/ read data in CSV format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "read_csv"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("read_csv: read a data file in CSV format.\n");
        printf("Syntax: read_csv <filename>.\n");
        printf("where <filename> is in the following CSV format.\n");
        printf("line 1 : (optional) variable name list\n");
        printf("line 2 : val1, val2, ... valn\n");
        printf(".....\n");
        printf("line m : val1, val2, ... valn\n\n");
        printf("That is, Line 2 specifies input/output variable names.\n");
        printf("  Line 3 has the input/output values of sample point 1\n");
        printf("  Line 4 has the input/output values of sample point 2\n");
        printf("  And so on.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: please specify a file to read from.\n");
        printf("Syntax: read_csv <filename>.\n");
        printf("where <filename> is in the following CSV format.\n");
        printf("line 1 : (optional) variable name list\n");
        printf("line 2 : val1, val2, ... valn\n");
        printf(".....\n");
        printf("line m : val1, val2, ... valn\n\n");
        printf("That is, Line 2 specifies input/output variable names.\n");
        printf("  Line 3 has the input/output values of sample point 1\n");
        printf("  Line 4 has the input/output values of sample point 2\n");
        printf("  And so on.\n");
        continue;
      }
      if ((fp=fopen(dataFile,"r")) == NULL)
      {
        printf("ERROR: file %s not found.\n", dataFile);
        printf("Syntax: read_csv <filename>.\n");
        printf("where <filename> is to be in the following CSV format.\n");
        printf("line 1 : (optional) variable list separated by commas)\n");
        printf("line 2 : val1, val2, ... valn\n");
        printf(".....\n");
        printf("line m : val1, val2, ... valn\n\n");
        printf("That is, Line 2 specifies input/output variable names.\n");
        printf("  Line 3 has the input/output values of sample point 1\n");
        printf("  Line 4 has the input/output values of sample point 2\n");
        printf("  And so on.\n");
        cmdStatus = 1;
        continue;
      }
      cleanUp();

      //**/ read the CSV file
      int    nInOuts;
      double *tempX;
      char   **inames, **onames;
      status = read_csv(dataFile,&nSamples_,&nInputs_,&tempX,&nOutputs_, 
                        &inames,&onames);
      if (status < 0 || nSamples_ <= 0 || nInputs_ <= 0) 
      {
        printf("ERROR: not a valid csv file (%d,%d,%d).\n",status,
               nSamples_,nInputs_);
        cmdStatus = 1;
        continue;
      }
      nInOuts = nInputs_ + nOutputs_;
      int flag = 0;
      if (nOutputs_ == 0)
      {
        sprintf(pString,"Enter number of input variables (1 - %d) : ",
                nInputs_);
        kk = getInt(1, nInputs_, pString);
        if (nInputs_ == kk)
        {
          printf("INFO: nOutputs = 0 ==> create a dummy output\n");
          nOutputs_ = 1;
          flag = 1;
          nInOuts = nInputs_;
        } 
        else
        {
          nInOuts = nInputs_;
          nOutputs_ = nInputs_ - kk;
          nInputs_  = kk;
          flag = 2;
        }
      } 
      printf("nSamples = %d\n", nSamples_);
      printf("nInputs  = %d\n", nInputs_);
      printf("nOutputs = %d\n", nOutputs_);
      VecSamInputs_.setLength(nSamples_*nInputs_);
      VecSamOutputs_.setLength(nSamples_*nOutputs_);
      VecSamStates_.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
          VecSamInputs_[ii*nInputs_+jj] = tempX[ii*nInOuts+jj]; 
        if (flag == 0 || flag == 2)
        {
          for (jj = 0; jj < nOutputs_; jj++)
            VecSamOutputs_[ii*nOutputs_+jj] = 
                      tempX[ii*nInOuts+nInputs_+jj];
        }
        else VecSamOutputs_[ii] = PSUADE_UNDEFINED;
        VecSamStates_[ii] = 1;
        for (jj = 0; jj < nOutputs_; jj++)
          if (VecSamOutputs_[ii*nOutputs_+jj] > 0.5*PSUADE_UNDEFINED)
            VecSamStates_[ii] = 0;
      }
      free(tempX);
      StrInpNames_.setNumStrings(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (inames != NULL && inames[ii] != NULL)
        {
          if (inames[ii][0] != '\0')
            StrInpNames_.loadOneString(ii, inames[ii]);
          else
          {
            sprintf(pString, "X%d", ii+1);
            StrInpNames_.loadOneString(ii, pString);
          }
          free(inames[ii]);
        }
        else 
        {
          sprintf(pString, "X%d", ii+1);
          StrInpNames_.loadOneString(ii, pString);
        }
      }
      StrOutNames_.setNumStrings(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        if (flag == 0 && onames != NULL && onames[ii] != NULL)
        {
          if (onames[ii][0] != '\0')
            StrOutNames_.loadOneString(ii, onames[ii]);
          else
          {
            sprintf(pString,"Y%d",ii+1);
            StrOutNames_.loadOneString(ii, pString);
          }
          free(onames[ii]);
        }
        else if (flag == 2 && inames != NULL && 
                 inames[nInputs_+ii] != NULL)
        {
          if (inames[nInputs_+ii][0] != '\0')
            StrOutNames_.loadOneString(ii, inames[nInputs_+ii]);
          else
          {
            sprintf(pString,"Y%d",ii+1);
            StrOutNames_.loadOneString(ii, pString);
          }
          free(inames[nInputs_+ii]);
        }
        else 
        {
          sprintf(pString,"Y%d",ii+1);
          StrOutNames_.loadOneString(ii, pString);
        }
      }
      if (inames != NULL) free(inames);
      if (onames != NULL) free(onames);
      VecILowerBs_.setLength(nInputs_);
      VecIUpperBs_.setLength(nInputs_);
      for (jj = 0; jj < nInputs_; jj++)
      {
        VecILowerBs_[jj] = VecSamInputs_[jj];
        VecIUpperBs_[jj] = VecSamInputs_[jj];
        for (ii = 1; ii < nSamples_; ii++)
        {
          if (VecSamInputs_[ii*nInputs_+jj] < VecILowerBs_[jj])
            VecILowerBs_[jj] = VecSamInputs_[ii*nInputs_+jj];
          if (VecSamInputs_[ii*nInputs_+jj] > VecIUpperBs_[jj])
            VecIUpperBs_[jj] = VecSamInputs_[ii*nInputs_+jj];
        }
        if (VecILowerBs_[jj] == VecIUpperBs_[jj]) 
          VecIUpperBs_[jj] += 1.0e-12;
      }
      SamplingMethod_ = PSUADE_SAMP_MC;
      VecInpPDFs_.setLength(nInputs_);
      VecInpMeans_.setLength(nInputs_);
      VecInpStds_.setLength(nInputs_);
      if (inputCMat_ != NULL) delete inputCMat_;
      inputCMat_  = new psMatrix();
      inputCMat_->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        VecInpPDFs_[ii] = 0;
        VecInpMeans_[ii] = 0;
        VecInpStds_[ii] = 0;
        inputCMat_->setEntry(ii,ii,1.0);
      }
      psuadeIO_ = new PsuadeData();
      psuadeIO_->setOutputLevel(0);
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                   VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                   VecInpStds_.getDVector(), inputCMat_); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      nReplications_ = 1;
      psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,
                                     nReplications_,0,0);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
      checkResidentSampleOutputs(); 
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ write_csv +++ 
    //**/ write data in CSV format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write_csv"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("write_csv: write a data file in CSV format.\n");
        printf("Syntax: write_csv <filename>.\n");
        printf("where <filename> will be in the following CSV format.\n");
        printf("line 1 : variable list separated by commas\n");
        printf("line 2 : val1, val2, ... valn\n");
        printf(".....\n");
        printf("line m : val1, val2, ... valn\n\n");
        printf("That is, Line 2 specifies input/output variable names.\n");
        printf("  Line 3 has the input/output values of sample point 1\n");
        printf("  Line 4 has the input/output values of sample point 2\n");
        printf("  And so on.\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: please specify a file to write to.\n");
        printf("Syntax: write_csv <filename>.\n");
        printf("where <filename> will be in the following CSV format.\n");
        printf("line 1 : variable list separated by commas\n");
        printf("line 2 : val1, val2, ... valn\n");
        printf(".....\n");
        printf("line m : val1, val2, ... valn\n\n");
        printf("That is, Line 2 specifies input/output variable names.\n");
        printf("  Line 3 has the input/output values of sample point 1\n");
        printf("  Line 4 has the input/output values of sample point 2\n");
        printf("  And so on.\n");
        cmdStatus = 1;
        continue;
      }
      if ((fp=fopen(dataFile,"w")) == NULL)
      {
        printf("write_csv ERROR: file %s cannot be opened.\n", dataFile);
        cmdStatus = 1;
        continue;
      }
      status = write_csv(dataFile,nSamples_,nInputs_,
                    VecSamInputs_.getDVector(),nOutputs_,
                    VecSamOutputs_.getDVector(),
                    StrInpNames_.getStrings(),StrOutNames_.getStrings()); 
      if (status != 0)
      {
        printf("write_csv ERROR: write not successful.\n");
        cmdStatus = 1;
      }
      else cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ write_matlab +++ 
    //**/ write the data to an output file in matlab format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write_matlab"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("write_matlab: write a data file in the matlab format.\n");
        printf("Syntax: write_matlab <filename>.\n");
        printf("where <filename> is the name of the target data file.\n\n");
        printf("INFO: the target file will have 2 matrices: X ");
        printf("has the sample inputs\n");
        printf("      and Y has the sample outputs.\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n"); 
        printf("          name %s. \n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
        cmdStatus = 1;
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        printf("Syntax: write_matlab <filename>.\n");
        cmdStatus = 1;
      }
      else
      {
        if (VecSamInputs_.length() == 0) 
        {
          printf("ERROR: no inputs to write.\n");
          cmdStatus = 1;
        }
        else
        {
          if (VecSamOutputs_.length() == 0)
            printf("INFO: no outputs to write (only inputs).\n");
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open file %s.\n", dataFile);
            cmdStatus = 1;
            continue;
          }
          fprintf(fp, "%% X - input matrix\n");
          fprintf(fp, "%% Y - output matrix\n");
          fprintf(fp, "X = [\n");
          for (ii = 0; ii < nSamples_; ii++)
          {
            for (jj = 0; jj < nInputs_; jj++)
              fprintf(fp, "%24.16e ", VecSamInputs_[ii*nInputs_+jj]); 
            fprintf(fp,"\n");
          }
          fprintf(fp,"];\n");
          if (VecSamOutputs_.length() > 0)
          {
            fprintf(fp, "Y = [\n");
            for (ii = 0; ii < nSamples_; ii++)
            {
              for (jj = 0; jj < nOutputs_; jj++)
                fprintf(fp,"%24.16e ",VecSamOutputs_[ii*nOutputs_+jj]);
              fprintf(fp,"\n");
            }
            fprintf(fp,"];\n");
          }
          fclose(fp);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ read_xls +++  
    //**/ read the data in Excel format and convert to psuade form
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "read_xls"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("read_xls: read a data file in Excel format.\n");
        printf("Syntax: read_xls <filename>.\n");
        printf("where <filename> is a data file in xls format:\n");
        printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
        printf("line 2: (optional) input and output names\n");
        printf("line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
        printf("line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
        printf(".......\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn,"%s %s",command,dataFile);
      if (!checkFileExists(dataFile))
      {
        printf("Syntax: read_xls <filename>.\n");
        printf("where <filename> is a data file in xls format:\n");
        cmdStatus = 1;
        continue;
      }
      cleanUp();

      //**/ read sample data 
      psuadeIO_ = new PsuadeData();
      psuadeIO_->setOutputLevel(0);
      fp = fopen(dataFile,"r");
      fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &nOutputs_);
      printf("nSamples = %d\n", nSamples_);
      printf("nInputs  = %d\n", nInputs_);
      printf("nOutputs = %d\n", nOutputs_);
      if (nSamples_ <= 0 || nInputs_ <= 0 || nOutputs_ < 0)
      {
        printf("read_xls ERROR: some sample parameters <= 0.\n");
        printf("NOTE: the first line should be: \n");
        printf("  <nSamples> <nInputs> <nOutputs>\n");
        fclose(fp);
        nSamples_ = nInputs_ = nOutputs_ = 0;
        cmdStatus = 1;
        continue;
      }
      int flag = 0;
      if (nOutputs_ == 0)
      {
        flag = 1;
        printf("INFO: No output given. Create a new output.\n");
        nOutputs_ = 1;
      }
      VecSamInputs_.setLength(nSamples_*nInputs_);
      VecSamOutputs_.setLength(nSamples_*nOutputs_);
      VecSamStates_.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++) VecSamStates_[ii] = 1;
      VecILowerBs_.setLength(nInputs_);
      VecIUpperBs_.setLength(nInputs_);
      //**/ detect comment line
      fgets(pString, 50000, fp);
      fgets(pString, 50000, fp);
      int hasComment = 0;
      if (pString[0] == '#') hasComment = 1;
      else
      {
        ii = 0;
        while ((pString[ii] == ' '  || pString[ii] == '\t') && 
                pString[ii] != '\r' && pString[ii] != '\n' && 
          ii < 5000) ii++;
        if (pString[ii] == '-') ii++;
        kk = (int) pString[ii] - '0';
        if (kk != 1) hasComment = 1;
        if (kk == 1 && (pString[ii+1] == ' ' || pString[ii+1] == '\t'))
             hasComment = 0;
        else hasComment = 1;
      }
      fclose(fp);
      fp = fopen(dataFile,"r");
      fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &kk);
      if (hasComment == 1)
      {
        StrInpNames_.setNumStrings(nInputs_);
        for (ii = 0; ii < nInputs_; ii++)
        {
          fscanf(fp,"%s", pString);
          StrInpNames_.loadOneString(ii, pString);
        }
        StrOutNames_.setNumStrings(nOutputs_);
        if (flag == 1) 
        {
          strcpy(pString, "Y");
          StrOutNames_.loadOneString(0, pString);
        }
        else
        {
          for (ii = 0; ii < nOutputs_; ii++)
          {
            fscanf(fp,"%s", pString);
            StrOutNames_.loadOneString(ii, pString);
          }
        }
        fgets(pString, 5000, fp);
      }
      for (ii = 0; ii < nSamples_; ii++)
      {
        fscanf(fp, "%d", &jj);
        if ((ii+1) != jj)
        {
          printf("read_xls ERROR: sample index mismatch.\n");
          printf("         sample number read     = %d\n", jj);
          printf("         sample number expected = %d\n", ii+1);
          printf("INFO: the first field must be the sample number.\n");
          cmdStatus = 1;
          break;
        }
        for (jj = 0; jj < nInputs_; jj++)
        {
          fscanf(fp,"%lg", &ddata);
          VecSamInputs_[ii*nInputs_+jj] = ddata; 
        }
        if (flag == 0)
        {
          for (jj = 0; jj < nOutputs_; jj++)
          {
            fscanf(fp,"%lg", &ddata);
            VecSamOutputs_[ii*nOutputs_+jj] = ddata; 
          }
        }
        else
        {
          VecSamOutputs_[ii] = PSUADE_UNDEFINED;
          VecSamStates_[ii] = 0;
        }
      }
      fclose(fp);
      if (ii != nSamples_)
      {
        VecSamInputs_.clean();
        VecSamOutputs_.clean();
        VecSamStates_.clean();
        VecILowerBs_.clean();
        VecIUpperBs_.clean();
        StrInpNames_.clean();
        StrOutNames_.clean();
        nInputs_ = nSamples_ = nOutputs_ = 0;
        continue;
      }
      psVector VecTW;
      VecTW.setLength(nSamples_);
      double xlsMin, xlsMax;
      for (ii = 0; ii < nInputs_; ii++)
      {
        for (jj = 0; jj < nSamples_; jj++) 
          VecTW[jj] = VecSamInputs_[jj*nInputs_+ii];
        xlsMin = xlsMax = VecTW[0];
        for (jj = 1; jj < nSamples_; jj++) 
        {
          if (VecTW[jj] > xlsMax) xlsMax = VecTW[jj];
          if (VecTW[jj] < xlsMin) xlsMin = VecTW[jj];
        }
        VecILowerBs_[ii] = xlsMin;
        VecIUpperBs_[ii] = xlsMax;
        if (VecILowerBs_[ii] == VecIUpperBs_[ii]) 
          VecIUpperBs_[ii] += 1.0e-12;
      }
      VecInpPDFs_.setLength(nInputs_);
      VecInpMeans_.setLength(nInputs_);
      VecInpStds_.setLength(nInputs_);
      if (inputCMat_ != NULL) delete inputCMat_;
      inputCMat_  = new psMatrix();
      inputCMat_->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        VecInpPDFs_[ii] = 0;
        VecInpMeans_[ii] = 0;
        VecInpStds_[ii] = 0;
        inputCMat_->setEntry(ii,ii,1.0);
      }
      StrInpNames_.setNumStrings(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString, "X%d", ii+1);
        StrInpNames_.loadOneString(ii, pString);
      }
      StrOutNames_.setNumStrings(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        sprintf(pString, "Y%d", ii+1);
        StrOutNames_.loadOneString(ii, pString);
      } 
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                   VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                   VecInpStds_.getDVector(),inputCMat_); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),StrOutNames_.getStrings());
      psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, nSamples_,1,0,0);
      psuadeIO_->updateAnalysisSection(0, 0, 0, 0, 0, 0);
      nReplications_ = 1;
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("Excel data have been read.\n");
      printf("Use 'write' to store the data in PSUADE format.\n");
      cmdStatus = 0;
      checkResidentSampleOutputs(); 
      VecTagArray_.clean();
      //**/ check for large sample only, since repeated points is a
      //**/ problem mainly for response surfaces
      if (nSamples_ < 10000)
      {
        status = checkRepeatedSamplePts(nSamples_,nInputs_,
                                      VecSamInputs_.getDVector());
        if (status != 0)
        {
          printf("WARNING: Repeated sample points have been ");
          printf("detected. This may cause\n");
          printf("         problems for some commands. So beware ");
          printf("(see if spurge or\n");
          printf("         rm_dup is needed).\n");
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ write_xls +++ 
    //**/ write the data to an output file in Excel format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write_xls"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("write_xls: write to a data file in Excel format\n");
        printf("           Use <read_xls -h> to see the Excel\n");
        printf("           format.\n");
        printf("Syntax: write_xls <filename>.\n");
        printf("where <filename> is the name of the target data file.\n");
        printf("\nThe target file will have the following format:\n");
        printf("line 1: <nSamples> <nInputs> <nOutputs>\n");
        printf("line 2: # input and output names\n");
        printf("line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
        printf("line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
        printf(".......\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n"); 
        printf("          name %s. \n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
        cmdStatus = 1;
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        printf("Syntax: write_xls <filename>.\n");
        cmdStatus = 1;
      }
      else
      {
        if (VecSamInputs_.length() == 0 || VecSamOutputs_.length() == 0)
        {
          printf("ERROR: no inputs or outputs to write.\n");
          cmdStatus = 1;
        }
        else
        {
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open file %s.\n", dataFile);
            cmdStatus = 1;
            continue;
          }
          fprintf(fp, "%d \t %d \t %d\n",nSamples_,nInputs_,nOutputs_);
          fprintf(fp, "#Sample \t ");
          for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp, "%s \t ", StrInpNames_[ii]);
          for (ii = 0; ii < nOutputs_-1; ii++)
            fprintf(fp, "%s \t ", StrOutNames_[ii]);
          fprintf(fp, "%s \n", StrOutNames_[nOutputs_-1]);
          for (ii = 0; ii < nSamples_; ii++)
          {
            fprintf(fp, "%d", ii+1);
            for (jj = 0; jj < nInputs_; jj++)
              fprintf(fp,"\t%24.16e", VecSamInputs_[ii*nInputs_+jj]); 
            for (jj = 0; jj < nOutputs_; jj++)
              fprintf(fp,"\t%24.16e", VecSamOutputs_[ii*nOutputs_+jj]);
            fprintf(fp, "\n");
          }
          fclose(fp);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ write_ultra 
    //**/ write the data to an output file in ULTRA format 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "write_ultra"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("write_ultra: write to a data file in the ULTRA data\n");
        printf("             format (suitable for 2 input only).\n");
        printf("Syntax: write_ultra <filename>.\n");
        printf("where <filename> is the name of the target data file.\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      if (!strcmp(dataFile,psInputFilename_))
      {
        printf("WARNING : output file name should not be the ");
        printf("same as the input file\n"); 
        printf("          name %s. \n", psInputFilename_);
        printf("NOTE: Command not executed.\n");
        cmdStatus = 1;
      }
      else if (!strcmp(dataFile, "\0"))
      {
        printf("ERROR: need to specify a file to write to.\n");
        printf("Syntax: write_ultra <filename>.\n");
        printf("where <filename> will be a data file in ultra format.\n");
        cmdStatus = 1;
        continue;
      }
      else
      {
        if (VecSamInputs_.length() == 0 || VecSamOutputs_.length() == 0)
        {
          printf("ERROR: no inputs or outputs to save.\n");
          cmdStatus = 1;
        }
        else if (nInputs_ != 2)
        {
          printf("ERROR: command available only for nInputs=2.\n");
          cmdStatus = 1;
        }
        else
        {
          if (nOutputs_ == 1) outputID = 0;
          else
          {
            sprintf(pString, "Enter output number (1 - %d) : ",nOutputs_);
            outputID = getInt(1, nOutputs_, pString);
            outputID--;
          }
          psVector  VecTX, VecTY, VecTW, VecTT, VecTI;
          VecTX.setLength(nSamples_);
          VecTY.setLength(nSamples_);
          VecTW.setLength(nSamples_);
          VecTT.setLength(nSamples_);
          VecTI.setLength(nSamples_);
          for (ii = 0; ii < nSamples_; ii++)
          {
            VecTX[ii] = VecSamInputs_[ii*nInputs_]; 
            VecTW[ii] = VecSamInputs_[ii*nInputs_+1]; 
            VecTY[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
            VecTI[ii] = (double) ii;
          }
          sortDbleList2(nSamples_, VecTW.getDVector(), 
                        VecTT.getDVector());
          for (ii = 0; ii < nSamples_; ii++) VecTT[ii] = VecTX[ii];
          for (ii = 0; ii < nSamples_; ii++)
          {
            ind2 = (int) VecTI[ii];
            VecTX[ii] = VecTT[ind2];
          }
          for (ii = 0; ii < nSamples_; ii++) VecTT[ii] = VecTY[ii];
          for (ii = 0; ii < nSamples_; ii++)
          {
            ind2 = (int) VecTI[ii];
            VecTY[ii] = VecTT[ind2];
          }
          ind = 0;
          double *TX = VecTX.getDVector();
          double *TY = VecTY.getDVector();
          for (ii = 1; ii < nSamples_; ii++)
          {
            if (VecTW[ii] != VecTW[ii-1])
            {
              sortDbleList2(ii-ind, &TX[ind], &TY[ind]);
              ind = ii;
            }
          } 
          sortDbleList2(nSamples_-ind, &TX[ind], &TY[ind]);
          fp = fopen(dataFile, "w");
          if (fp == NULL)
          {
            printf("ERROR: cannot open output file %s\n", dataFile);
            cmdStatus = 1;
            continue;
          }
          ind = 0;
          ind2 = 1;
          for (ii = 1; ii < nSamples_; ii++)
          {
            if (VecTW[ii] != VecTW[ii-1])
            {
              fprintf(fp, "#var2_%d\n", ind2);
              for (jj = ind; jj < ii; jj++)
                fprintf(fp, "%24.16e %24.16e\n",VecTX[jj],VecTY[jj]);
              ind = ii;
              ind2++;
            }
          } 
          fprintf(fp, "#var2_%d\n", ind2);
          for (jj = ind; jj < nSamples_; jj++)
            fprintf(fp, "%14.6e %24.16e\n",VecTX[jj],VecTY[jj]);
          fclose(fp);
          printf("ultra file created in %s\n", dataFile);
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ oupdate +++ 
    //**/ merge in the output of identical sample point from another file
    //**/ ----------------------------------------------------------------
    else if (!strcmp(command, "oupdate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oupdate: update the loaded sample outputs from ");
        printf("another PSUADE data\n");
        printf("         file (i.e. the current sample outputs ");
        printf("will be replaced with\n");
        printf("         the evaluated outputs from the provided ");
        printf("sample file). New\n");
        printf("         points in the second sample will be ignored.\n");
        printf("Syntax: oupdate <filename> (merge from <filename>).\n");
        printf("where <filename> is a data file in PSUADE data format.\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("Syntax: oupdate <filename>).\n");
        printf("where <filename> is a data file in PSUADE data format.\n");
        cmdStatus = 1;
        continue;
      }
      fclose(fp);
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
               dataFile);
        cmdStatus = 1;
        continue;
      }
      ioPtr->getParameter("input_ninputs", pPtr);
      ind = pPtr.intData_;
      int flag = 1;
      if (ind != nInputs_) flag = 0;
      if (flag == 1)
      {
        pINames.clean();
        pLower.clean();
        pUpper.clean();
      }
      if (flag == 1)
      {
        ioPtr->getParameter("output_noutputs", pPtr);
        ind = pPtr.intData_;
        if (ind != nOutputs_) flag = 0;
      }
      if (flag == 1)
      {
        pONames.clean();
        ioPtr->getParameter("output_names", pONames);
        names = pONames.strArray_;
        for (ii = 0; ii < ind; ii++)
          if (strcmp(StrOutNames_[ii],names[ii])) flag = 0;
        pONames.clean();
      }
      if (flag == 0) 
      {
        printf("ERROR: invalid data file.\n");
        cmdStatus = 1;
      }
      else
      {
        pInputs.clean();
        pOutputs.clean();
        pStates.clean();
        ioPtr->getParameter("input_sample",pInputs);
        ioPtr->getParameter("output_sample",pOutputs);
        ioPtr->getParameter("output_states",pStates);
        ioPtr->getParameter("method_nsamples", pPtr);
        int nSamp2 = pPtr.intData_;
        count = 0;
        vecTags.setLength(nSamp2);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          //**/ update regardless of sample state
          //if (VecSamStates_[sInd] == 1) continue;
          kk = sInd * nInputs_;
          for (ii = 0; ii < nSamp2; ii++)
          {
            if (pStates.intArray_[ii] != 1) 
            {
              vecTags[ii] = -1;
              continue;
            }
            jj = ii * nInputs_;
            for (iInd = 0; iInd < nInputs_; iInd++)
              if (VecSamInputs_[kk+iInd] != pInputs.dbleArray_[jj+iInd])
                break;
            if (iInd == nInputs_)
            {
              kk = sInd * nOutputs_;
              jj = ii * nOutputs_;
              for (oInd = 0; oInd < nOutputs_; oInd++)
                VecSamOutputs_[kk+oInd] = pOutputs.dbleArray_[jj+oInd];
              VecSamStates_[sInd] = 1;
              printf("   Matched: sample %d <-- sample %d (%d)\n",
                     sInd+1,ii+1,nSamp2);
              count++;
              vecTags[ii]++;
              break;
            }
          }
        }
        printf("Number of sample points updated = %d\n", count);
        for (ii = 0; ii < nSamp2; ii++)
        {
          if (vecTags[ii] > 0) 
            printf("Sample %d has been used %d times\n",ii+1,vecTags[ii]);
          else if (vecTags[ii] < 0) 
            printf("Sample %d is invalid (has undefined values)\n",ii+1);
          else 
            printf("Sample %d has not been used\n",ii+1);
        }
        pInputs.clean();
        pOutputs.clean();
        pStates.clean();
      }
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(),
                   VecSamStates_.getIVector(), NULL); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      fflush(stdout);
      if (ioPtr != NULL) delete ioPtr;
      cmdStatus = 0;
      checkResidentSampleOutputs(); 
    }

    //**/ -------------------------------------------------------------
    // +++ moat_adjust +++ 
    //**/ adjust an MOAT sample file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moat_adjust"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moat_adjust: adjust a MOAT sample with new input values\n");
        printf("Syntax: moat_adjust <fname>.\n");
        printf("where <fname> is a MOAT file (created by moatgen)\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      status = 0;
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("Syntax: moat_adjust <fname>.\n");
        printf("where <fname> is a MOAT file (created by moatgen)\n");
        status = 1;
        cmdStatus = 1;
      }
      else fclose(fp);

      if (SamplingMethod_ != PSUADE_SAMP_MOAT)
      {
        printf("Cannot adjust: the sampling method is not MOAT.\n");
        status = 1;
        cmdStatus = 1;
        continue;
      }
      Sampling *sampPtr =
           (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MOAT);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, 
                 VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(nSamples_,nSamples_/(nInputs_+1),0);
      sampPtr->initialize(1);
      sampPtr->loadSamples(nSamples_,nInputs_,nOutputs_,
                   VecSamInputs_.getDVector(), 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector());
      sampPtr->repair(dataFile, 0);
      sampPtr->getSamples(nSamples_,nInputs_,nOutputs_,
                VecSamInputs_.getDVector(), 
                VecSamOutputs_.getDVector(),
                VecSamStates_.getIVector());
      delete sampPtr;
      sampPtr = NULL;
      fflush(stdout);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gmoat_adjust +++ 
    //**/ adjust an GMOAT sample file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gmoat_adjust"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gmoat_adjust: adjust a MOAT sample with new values\n");
        printf("NOTE: gmoat_adjust differs from moat_adjust ");
        printf("in that it uses GMOAT\n");
        printf("      instead of MOAT to create the initial sample.\n");
        printf("Syntax: gmoat_adjust <filename>.\n");
        printf("where <filename> is a MOAT file (created by moatgen)\n");
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("When the input space is not a hyper-rectangle, ");
      printf("the standard MOAT\n");
      printf("method will have difficulties (sample points ");
      printf("may not be inside the\n");
      printf("feasible region). This command circumvents this ");
      printf("difficulties by\n");
      printf("adjusting the loaded GMOAT sample with valid ");
      printf("GMOAT sample points\n");
      printf("created by gmoatgen.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      status = 0;
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("Syntax: gmoat_adjust <file>.\n");
        printf("where <file> is a MOAT adjust file (created by moatgen)\n");
        status = 1;
        cmdStatus = 1;
      }
      else fclose(fp);

      if (SamplingMethod_ != PSUADE_SAMP_GMOAT)
      {
        printf("Cannot adjust: the sampling method is not GMOAT.\n");
        status = 1;
        cmdStatus = 1;
        continue;
      }
      Sampling *sampPtr =
              (Sampling *) SamplingCreateFromID(PSUADE_SAMP_GMOAT);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(nSamples_,nSamples_/(nInputs_+1),0);
      sampPtr->initialize(1);
      sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, 
               VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
               VecSamStates_.getIVector());
      sampPtr->repair(dataFile, 0);
      sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
               VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
               VecSamStates_.getIVector());
      delete sampPtr;
      sampPtr = NULL;
      cmdStatus = 0;
      fflush(stdout);
    }

    //**/ -------------------------------------------------------------
    // +++ moatgen +++ 
    //**/ generate MOAT adjust file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moatgen"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moatgen: create a Morris adjust file with 1 constraint\n");
        printf("Syntax: moatgen\n");
        printf("NOTE: A PSUADE datafile should have been ");
        printf("loaded before using this\n");
        printf("      command. The data will be used as ");
        printf("a response surface to find\n");
        printf("      feasible points in creating a MOAT ");
        printf("adjust sample.\n");
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("When the input space is not a hyper-rectangle, ");
      printf("the standard MOAT\n");
      printf("method will have difficulties (sample points ");
      printf("may not be inside the\n");
      printf("feasible region). This command circumvents this ");
      printf("difficulties by\n");
      printf("creating MOAT sample points that live inside ");
      printf("the constrained input\n");
      printf("space. To do so, it needs to know how the ");
      printf("input space is constrained\n");
      printf("(with inequality constraints) via evaluating ");
      printf("and bounding the outputs\n");
      printf("at new MOAT sample points using the response ");
      printf("surface constructed from\n");
      printf("the loaded sample.");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }

      printf("Constructing a RS from the loaded sample ");
      printf("for setting up constraints...\n");
      faFlag = 3;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) {printf("ERROR detected in RS.\n"); continue;}
      faPtr->setOutputLevel(outputLevel_);
      psuadeIO_->getParameter("ana_outputid", pPtr);
      outputID = pPtr.intData_;
      Ymax = - 1.0e35;
      Ymin =   1.0e35;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
          Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
          Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
      }
      printf("Obtaining constraint bounds from user (a MOAT ");
      printf("candidate sample point X\n");
      printf("is feasible if Ymin <= F(X) <= Ymax) where F(X) ");
      printf("is from response\n");
      printf("surface evaluation.\n");
      sprintf(pString,
              "Enter the lower bound constraint (Ymin=%e) : ",Ymin);
      threshL = getDouble(pString);
      sprintf(pString,
              "Enter the upper bound constraint (Ymax=%e) : ",Ymax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        cmdStatus = 1;
        continue;
      }

      printf("P is the MOAT resolution (normally equal to 4) ");
      printf("but can be otherwise.\n");
      int nPaths = 5000;
      sprintf(pString, "Enter P (resolution: try 4-10) : ");
      int currP = getInt(4, 10, pString);
      sprintf(pString, "Enter the number of trials (> 100) : ");
      int nTrials = getInt(101, 10000000, pString);
      psMatrix matMOATSample;
      matMOATSample.setFormat(PS_MAT2D);
      matMOATSample.setDim(nPaths*(nInputs_+1), nInputs_);
      double **moatSample = matMOATSample.getMatrix2D();
      psVector  VecTW;
      psIVector VecTI;
      VecTW.setLength(nInputs_);
      VecTI.setLength(nInputs_);
      count = 0;
      double sLo, sHi, currX1, currX2, currY1, currY2;
      double filterRange = threshU - threshL;
      int    trial=0;
      for (ii = 0; ii < nPaths; ii++)
      {
        trial = 0; 
        while (trial < nTrials)
        {
          iInd = count;
          trial++;
          for (jj = 0; jj < nInputs_; jj++)
          {
            ind = PSUADE_rand() % currP;
            dtemp = ind * (VecIUpperBs_[jj] - VecILowerBs_[jj]) / 
                   (currP - 1.0);
            moatSample[iInd][jj] = VecILowerBs_[jj] + dtemp;
          }
          dtemp = faPtr->evaluatePoint(moatSample[iInd]);
          if (dtemp < threshL || dtemp > threshU) continue;

          //**/ randomize the selection of next input
          generateRandomIvector(nInputs_, VecTI.getIVector());

          //**/ construct the path
          for (jj = 0; jj < nInputs_; jj++)
          {
            iInd++;
            //**/ always start from last point
            for (kk = 0; kk < nInputs_; kk++)
              moatSample[iInd][kk] = moatSample[iInd-1][kk];

            //**/ pick the next input
            ind2 = VecTI[jj];

            //**/ see if both end within feasible region
            ddata = moatSample[iInd][ind2]; 
            moatSample[iInd][ind2] = VecILowerBs_[ind2];
            currY1 = faPtr->evaluatePoint(moatSample[iInd]);
            moatSample[iInd][ind2] = VecIUpperBs_[ind2];
            currY2 = faPtr->evaluatePoint(moatSample[iInd]);
            currX1 = VecILowerBs_[ind2];
            currX2 = VecIUpperBs_[ind2];
            moatSample[iInd][ind2] = ddata;

            //**/ if not, the search range is empty
            if (currY2 >= threshU && currY1 >= threshU)
              currX1 = currX2 = 0.0;
            else if (currY2 <= threshL && currY1 <= threshL)
              currX1 = currX2 = 0.0;
            else if (currY2 > currY1)
            {
              //**/ binary search
              if (currY2 <= threshU) currX2 = VecIUpperBs_[ind2];
              else
              {
                sLo = VecILowerBs_[ind2];
                sHi = VecIUpperBs_[ind2];
                while (PABS((currY2-threshU)/filterRange)>0.0001)
                {
                  moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                  currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                  if (currY2 > threshU) sHi = 0.5 * (sLo + sHi);
                  else                  sLo = 0.5 * (sLo + sHi);
                }
                currX2 = moatSample[iInd][ind2];
              }
              if (currY1 >= threshL) currX1 = VecILowerBs_[ind2];
              else
              {
                sLo = VecILowerBs_[ind2];
                sHi = VecIUpperBs_[ind2];
                while (PABS((currY1-threshL)/filterRange)>0.0001)
                {
                  moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                  currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                  if (currY1 < threshL) sLo = 0.5 * (sLo + sHi);
                  else                  sHi = 0.5 * (sLo + sHi);
                }
                currX1 = moatSample[iInd][ind2];
              }
            }
            else
            {
              if (currY1 <= threshU) currX1 = VecILowerBs_[ind2];
              else
              {
                sLo = VecILowerBs_[ind2];
                sHi = VecIUpperBs_[ind2];
                while (PABS((currY1-threshU)/filterRange)>0.0001)
                {
                  moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                  currY1 = faPtr->evaluatePoint(moatSample[iInd]);
                  if (currY1 > threshU) sLo = 0.5 * (sLo + sHi);
                  else                  sHi = 0.5 * (sLo + sHi);
                }
                currX1 = moatSample[iInd][ind2];
              }
              if (currY2 >= threshL) currX2 = VecIUpperBs_[ind2];
              else
              {
                sLo = VecILowerBs_[ind2];
                sHi = VecIUpperBs_[ind2];
                while (PABS((currY2-threshL)/filterRange)>0.0001)
                {
                  moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                  currY2 = faPtr->evaluatePoint(moatSample[iInd]);
                  if (currY2 < threshL) sHi = 0.5 * (sLo + sHi);
                  else                  sLo = 0.5 * (sLo + sHi);
                }
                currX2 = moatSample[iInd][ind2];
              }
            }
            if (PABS(currX2-currX1)<
                0.1*(VecIUpperBs_[ind2]-VecILowerBs_[ind2])) 
              break;
            moatSample[iInd][ind2] = ddata;
            VecTW[ind2] = PABS(currX2 - currX1) / (currP-1.0);
            moatSample[iInd][ind2] += VecTW[ind2];
            if (moatSample[iInd][ind2] > VecIUpperBs_[ind2])
                 dtemp = threshL - 1.0;
            else dtemp = faPtr->evaluatePoint(moatSample[iInd]);
            if (dtemp < threshL || dtemp > threshU) 
            {
              moatSample[iInd][ind2] -= 2.0 * VecTW[ind2];
              if (moatSample[iInd][ind2] < VecILowerBs_[ind2])
                break;
              dtemp = faPtr->evaluatePoint(moatSample[iInd]);
              if (dtemp < threshL || dtemp > threshU) 
                break;
            }
          }
          if (jj == nInputs_) 
          {
            count += (nInputs_ + 1);
            if (outputLevel_ > 2)
              printf("moatgen: path %d (out of %d) found.\n", ii+1,
                     nPaths);
            break; 
          }
          else
          {
            if (outputLevel_ > 2)
              printf("Current path fails (%d out of max %d).\n",
                     trial, nTrials); 
          }
        }
        if (trial >= nTrials)
        {
          printf("ERROR: moatgen fails to find all possible paths.\n");
          printf("Suggestion: try again with a larger P than %d.\n",
                 currP);
          cmdStatus = 1;
          break;
        }
      }
      for (ii = 0; ii < nPaths; ii++)
      {
        dtemp = faPtr->evaluatePoint(moatSample[ii]);
        if (dtemp < threshL || dtemp > threshU)
          printf("moatgen: sample %d fails final test (%e <? %e <? %e)\n",
                 ii, threshL, dtemp, threshU);
      }
      if (trial >= nTrials) continue;

      fp = fopen("MOAT_adjust_file", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file MOAT_adjust_file.\n");
        cmdStatus = 1;
        continue;
      }
      fprintf(fp, "BEGIN\n");
      fprintf(fp, "%d %d\n", nPaths*(nInputs_+1), nInputs_);
      for (ii = 0; ii < nInputs_; ii++) fprintf(fp, "%d ", ii+1);
      fprintf(fp, "\n");
      for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
          fprintf(fp, "%e ", moatSample[ii][jj]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "END\n");
      fclose(fp);
      count = nPaths * (nInputs_ + 1); 
      vecWT.setLength(count*nInputs_);
      vecST.setLength(count);
      for (ii = 0; ii < count; ii++) vecST[ii] = 1;
      for (ii = 0; ii < count; ii++)
        for (jj = 0; jj < nInputs_; jj++)
          vecWT[ii*nInputs_+jj] = moatSample[ii][jj];
      printf("moatgen: check for repeated sample points.\n");
      for (ii = 0; ii < count; ii++)
      {
        status = compareSamples(ii,count,nInputs_,vecWT.getDVector(), 
                                vecST.getIVector());
        if (status >= 0)
          printf("moatgen check: sample %d and %d are identical.\n",
                 ii+1,status+1);
      }
      printf("moatgen: adjust file created in MOAT_adjust_file. ");
      printf("Make sure to change\n");
      printf("         the input indices to match the indices ");
      printf("in the original MOAT\n");
      printf("         file when used with gmoat_adjust.\n");
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ moatgen2 +++ 
    //**/ generate MOAT adjust file with multiple constraints
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moatgen2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moatgen2: create a MOAT adjust file with ");
        printf("multiple constraints\n");
        printf("          (moatgen: for single constraint)\n");
        printf("Syntax: moatgen (no argument needed)\n");
        printf("NOTE: A PSUADE datafile should have been loaded ");
        printf("prior to using this\n");
        printf("      command to define input bounds for MOAT ");
        printf("search. Users will be\n");
        printf("      prompted for multiple data file to build ");
        printf("response surfaces for\n");
        printf("      carving feasible regions in guiding MOAT ");
        printf("path search.\n");
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command is similar to moatgen but is ");
      printf("generalized for multiple\n");
      printf("constraints.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      status = 0;
      sprintf(pString,
         "How many constraint data files are there (1-10)? ");
      int nFiles = getInt(1, 10, pString);
      faFlag = 2;
      FuncApprox **faPtrs = new FuncApprox*[nFiles];
      psVector vecThreshLs, vecThreshUs;
      vecThreshLs.setLength(nFiles);
      vecThreshUs.setLength(nFiles);
      PsuadeData *ioPtr = NULL;
      for (kk = 0; kk < nFiles; kk++)
      {
        sprintf(pString,
                "Enter name of constraint file #%d : ",kk+1);
        getString(pString, winput);
        kk = strlen(winput);
        winput[kk-1] = '\0';
        ioPtr = new PsuadeData;
        status = ioPtr->readPsuadeFile(winput);
        if (status != 0)
        {
          printf("moatgen2 FILE READ ERROR: file = %s\n", winput);
          cmdStatus = 1;
          break;
        }
        ioPtr->getParameter("input_ninputs", pPtr);
        jj = pPtr.intData_;
        if (jj != nInputs_)
        {
          printf("moatgen2 ERROR: nInputs mismatch.\n");
          printf("         local nInputs = %d.\n", nInputs_);
          printf("         file  nInputs = %d.\n", jj);
          status = cmdStatus = 1;
          break;
        }
        pLower.clean();
        ioPtr->getParameter("input_lbounds", pLower);
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (VecILowerBs_[ii] != pLower.dbleArray_[ii])
          {
            printf("moatgen2 ERROR: lower bound mismatch (input %d)\n", 
                   ii+1);
            status = cmdStatus = 1;
            break;
          }
        }
        pUpper.clean();
        ioPtr->getParameter("input_ubounds", pUpper);
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (VecIUpperBs_[ii] != pUpper.dbleArray_[ii])
          {
            printf("moatgen2 ERROR: upper bound mismatch (input %d)\n",
                   ii+1);
            status = cmdStatus = 1;
            break;
          }
        }
        faPtrs[kk] = genFAInteractive(ioPtr, faFlag);
        if (faPtrs[kk] == NULL) 
        {
          printf("ERROR detected in RS.\n");
          status = cmdStatus = 1;
          break;
        }
        faPtrs[kk]->setOutputLevel(outputLevel_);
        sprintf(pString,"Constraint %d lower bound : ",kk+1);
        vecThreshLs[kk] = getDouble(pString);
        sprintf(pString,"Constraint %d upper bound : ",kk+1);
        vecThreshUs[kk] = getDouble(pString);
        if (vecThreshLs[kk] >= vecThreshUs[kk])
        {
          printf("ERROR: lower bound >= upper bound.\n");
          status = 1;
          cmdStatus = 1;
          delete ioPtr;
          break;
        }
        delete ioPtr;
      }
      if (status != 0) continue;
      int nPaths = 5000;
      printf("P is the MOAT resolution (normally = 4).\n");
      sprintf(pString, "Enter P (resolution: try 4-10) : ");
      int currP = getInt(4, 10, pString);
      printf("INFO: Searching for a MOAT path may come up ");
      printf("empty due to complex\n");
      printf("      constraints. So we allow a number of ");
      printf("trials before giving up.\n"); 
      sprintf(pString, 
              "Enter the maximum number of trials (> 100) : ");
      int nTrials = getInt(101, 10000000, pString);
      psMatrix matMOATSample;
      matMOATSample.setFormat(PS_MAT2D);
      matMOATSample.setDim(nPaths*(nInputs_+1), nInputs_);
      double **moatSample = matMOATSample.getMatrix2D();
      vecWT.setLength(nInputs_);
      vecIT.setLength(nInputs_);
      count = 0;
      double sLo,sHi,filterRange;
      double currX1, currX2, currY1, currY2, currX1F, currX2F;
      int    trial=0;
      for (ii = 0; ii < nPaths; ii++)
      {
        trial = 0; 
        while (trial < nTrials)
        {
          iInd = count;
          trial++;
          //**/ find a random initial point for all inputs 
          for (jj = 0; jj < nInputs_; jj++)
          {
            ind = PSUADE_rand() % currP;
            dtemp = ind*(VecIUpperBs_[jj]-VecILowerBs_[jj])/(currP-1);
            moatSample[iInd][jj] = VecILowerBs_[jj] + dtemp;
          }
          //**/ see if the point satisfies the multiple constraints
          //**/ (note: constrain the outputs so RS is needed
          for (kk = 0; kk < nFiles; kk++)
          {
            dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
            if (dtemp < vecThreshLs[ii] || dtemp > vecThreshUs[ii]) 
              break;
          }
          if (kk != nFiles) continue;

          //**/ randomize the selection of next input
          generateRandomIvector(nInputs_, vecIT.getIVector());

          //**/ construct the path
          for (jj = 0; jj < nInputs_; jj++)
          {
            iInd++;
            //**/ always start from last point
            for (kk = 0; kk < nInputs_; kk++)
              moatSample[iInd][kk] = moatSample[iInd-1][kk];

            //**/ pick the next input
            ind2 = vecIT[jj];

            //**/ see if both end within feasible region
            ddata = moatSample[iInd][ind2]; 
            currX1F = - PSUADE_UNDEFINED;
            currX2F =   PSUADE_UNDEFINED;
            for (kk = 0; kk < nFiles; kk++)
            {
              filterRange = vecThreshUs[kk] - vecThreshLs[kk];
              moatSample[iInd][ind2] = VecILowerBs_[ind2];
              currY1 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
              moatSample[iInd][ind2] = VecIUpperBs_[ind2];
              currY2 = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
              currX1 = VecILowerBs_[ind2];
              currX2 = VecIUpperBs_[ind2];

              //**/ if not, the search range is empty
              if (currY2 >= vecThreshUs[kk] && currY1 >= vecThreshUs[kk])
                 currX1 = currX2 = 0.0;
              else if (currY2 <= vecThreshLs[kk] && 
                       currY1 <= vecThreshLs[kk])
                 currX1 = currX2 = 0.0;
              else if (currY2 > currY1)
              {
                //**/ binary search
                if (currY2 <= vecThreshUs[kk]) 
                  currX2 = VecIUpperBs_[ind2];
                else
                {
                  sLo = VecILowerBs_[ind2];
                  sHi = VecIUpperBs_[ind2];
                  while (PABS((currY2-vecThreshUs[kk])/filterRange)>1e-4)
                  {
                    moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                    currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                    if (currY2 > vecThreshUs[kk]) sHi = 0.5 * (sLo+sHi);
                    else                          sLo = 0.5 * (sLo+sHi);
                  }
                  currX2 = moatSample[iInd][ind2];
                }
                if (currY1 >= vecThreshLs[kk]) currX1 = VecILowerBs_[ind2];
                else
                {
                  sLo = VecILowerBs_[ind2];
                  sHi = VecIUpperBs_[ind2];
                  while (PABS((currY1-vecThreshLs[kk])/filterRange)>1e-4)
                  {
                    moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                    currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                    if (currY1 < vecThreshLs[kk]) sLo = 0.5 * (sLo+sHi);
                    else                          sHi = 0.5 * (sLo+sHi);
                  }
                  currX1 = moatSample[iInd][ind2];
                }
              }
              else
              {
                if (currY1 <= vecThreshUs[kk]) currX1 = VecILowerBs_[ind2];
                else
                {
                  sLo = VecILowerBs_[ind2];
                  sHi = VecIUpperBs_[ind2];
                  while (PABS((currY1-vecThreshUs[kk])/filterRange)>1e-4)
                  {
                    moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                    currY1=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                    if (currY1 > vecThreshUs[kk]) sLo = 0.5 * (sLo+sHi);
                    else                          sHi = 0.5 * (sLo+sHi);
                  }
                  currX1 = moatSample[iInd][ind2];
                }
                if (currY2 >= vecThreshLs[kk]) currX2 = VecIUpperBs_[ind2];
                else
                {
                  sLo = VecILowerBs_[ind2];
                  sHi = VecIUpperBs_[ind2];
                  while (PABS((currY2-vecThreshLs[kk])/filterRange)>1e-4)
                  {
                    moatSample[iInd][ind2] = 0.5 * (sLo + sHi);
                    currY2=faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                    if (currY2 < vecThreshLs[kk]) sHi = 0.5 * (sLo+sHi);
                    else                          sLo = 0.5 * (sLo+sHi);
                  }
                  currX2 = moatSample[iInd][ind2];
                }
              }
              if (PABS(currX2-currX1) < 
                  0.1*(VecIUpperBs_[ind2]-VecILowerBs_[ind2])) 
                break;
              if (currX1 > currX1F) currX1F = currX1;
              if (currX2 < currX2F) currX2F = currX2;
            }
            moatSample[iInd][ind2] = ddata;
            if (kk != nFiles) break;
            vecWT[ind2] = PABS(currX2F - currX1F) / (currP-1.0);
            moatSample[iInd][ind2] += vecWT[ind2];
            if (moatSample[iInd][ind2] > VecIUpperBs_[ind2])
              dtemp = vecThreshLs[kk] - 1.0;
            else
            {
              for (kk = 0; kk < nFiles; kk++)
              {
                moatSample[iInd][ind2] = ddata;
                moatSample[iInd][ind2] += vecWT[ind2];
                dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                if (dtemp < vecThreshLs[kk] || dtemp > vecThreshUs[kk]) 
                {
                  moatSample[iInd][ind2] -= 2.0 * vecWT[ind2];
                  if (moatSample[iInd][ind2] < VecILowerBs_[ind2])
                    break;
                  dtemp = faPtrs[kk]->evaluatePoint(moatSample[iInd]);
                  if (dtemp < vecThreshLs[kk] || dtemp > vecThreshUs[kk]) 
                    break;
                }
              }
              moatSample[iInd][ind2] = ddata;
              if (kk != nFiles) break;
            }
          }
          if (jj == nInputs_) 
          {
            count += (nInputs_ + 1);
            if (outputLevel_ > 2)
              printf("moatgen2: path %d (out of %d) found.\n", ii+1,
                     nPaths);
            break; 
          }
          else
          {
            if (outputLevel_ > 2)
              printf("Current path fails (%d out of max %d).\n",
                     trial, nTrials); 
          }
        }
        if (trial >= nTrials)
        {
          printf("moatgen2 fails to find all possible paths.\n");
          printf("Suggestion: try again with a larger P than %d.\n", 
                 currP);
          break;
        }
      }
      for (ii = 0; ii < nPaths; ii++)
      {
        for (kk = 0; kk < nFiles; kk++)
        {
          dtemp = faPtrs[kk]->evaluatePoint(moatSample[ii]);
          if (dtemp < vecThreshLs[kk] || dtemp > vecThreshUs[kk])
            printf("moatgen2: sample %d fails final test (%e <? %e <? %e)\n",
                   ii, vecThreshLs[kk], dtemp, vecThreshUs[kk]);
        }
      }
      for (kk = 0; kk < nFiles; kk++) delete faPtrs[kk];
      delete [] faPtrs;
      faPtrs = NULL;
      if (trial >= nTrials) continue;
      fp = fopen("MOAT_adjust_file", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file MOAT_adjust_file.\n");
        continue;
      }
      fprintf(fp, "BEGIN\n");
      fprintf(fp, "%d %d\n", nPaths*(nInputs_+1), nInputs_);
      for (ii = 0; ii < nInputs_; ii++) fprintf(fp, "%d ", ii+1);
      fprintf(fp, "\n");
      for (ii = 0; ii < nPaths*(nInputs_+1); ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
          fprintf(fp, "%e ", moatSample[ii][jj]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "END\n");
      fclose(fp);
      count = nPaths * (nInputs_ + 1); 
      vecWT.setLength(count*nInputs_);
      vecST.setLength(count);
      for (ii = 0; ii < count; ii++) vecST[ii] = 1;
      for (ii = 0; ii < count; ii++)
        for (jj = 0; jj < nInputs_; jj++)
          vecWT[ii*nInputs_+jj] = moatSample[ii][jj];
      printf("moatgen2: check for repeated sample points.\n");
      for (ii = 0; ii < count; ii++)
      {
        status = compareSamples(ii,count,nInputs_, vecWT.getDVector(), 
                                vecST.getIVector());
        if (status >= 0)
          printf("moatgen2 check: sample %d and %d are identical.\n",
                 ii+1,status+1);
      }
      printf("moatgen2: adjust file created in MOAT_adjust_file.\n");
      printf("          Make sure to change the input indices ");
      printf("to match the indices\n");
      printf("          in the original MOAT file when used with ");
      printf("gmoat_adjust.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ moat_concatenate +++
    //**/ concatenate 2 MOAT samples with 2 different input sets
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moat_concat"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moat_concat: concatenate 2 MOAT samples (with ");
        printf("different input sets)\n");
        printf("Syntax: moat_concat <file>\n");
        printf("NOTE: A PSUADE MOAT datafile must have been ");
        printf("loaded before using this\n");
        printf("      command. \n");
      }

      printAsterisks(PL_INFO, 0);
      printf("This command is useful if you have created 2 ");
      printf("MOAT samples using 2\n");
      printf("distinct sets of inputs, and you would like ");
      printf("to concatenate them to\n");
      printf("form a single MOAT sample for the 2 input sets.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        continue;
      }

      //**/ check the second MOAT file
      strcpy(dataFile, "\0");
      sscanf(lineIn,"%s %s",command,dataFile);
      if ((fp=fopen(dataFile,"r")) == NULL)
      {
        printf("file %s not found.\n", dataFile);
        printf("Syntax: moat_concat <file>.\n");
        printf("where <file> is a PSUADE data file.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ check and read the second MOAT file
      psuadeIO_->getParameter("method_sampling", pPtr);
      if (pPtr.intData_ != PSUADE_SAMP_MOAT)
      {
        printf("ERROR: loaded sample is not MOAT. \n");
        cmdStatus = 1;
        continue;
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->setOutputLevel(0);
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: when reading second MOAT sample.\n");
        cmdStatus = 1;
        continue;
      }
      ioPtr->getParameter("method_sampling", pPtr);
      if (pPtr.intData_ != PSUADE_SAMP_MOAT)
      {
        printf("ERROR: second sample is not MOAT. \n");
        cmdStatus = 1;
        continue;
      }

      //**/ check number of replications for both
      nReplications_ = nSamples_ / (nInputs_ + 1);
      ioPtr->getParameter("input_ninputs", pPtr);
      ind = pPtr.intData_;
      ioPtr->getParameter("method_nsamples", pPtr);
      count = pPtr.intData_;
      ll = count / (ind + 1);
      if (nReplications_ != ll)
      {
        printf("ERROR: different number of replications.\n");
        printf("       num_replcations for sample 1 = %d\n",
               nReplications_);
        printf("       num_replcations for sample 2 = %d\n",ll);
        cmdStatus = 1;
        continue;
      }

      //**/ add new input names
      StrInpNames_.addMoreStrings(ind);
      pINames.clean();
      psuadeIO_->getParameter("input_names", pINames);
      names = pINames.strArray_;
      for (ii = nInputs_; ii < nInputs_+ind; ii++)
        StrInpNames_.loadOneString(ii, names[ii-nInputs_]);

      //**/ add to inputs, outputs, and status
      vecXT = VecSamInputs_;
      VecSamOutputs_.clean();
      VecSamStates_.clean();
      kk = nReplications_ * (nInputs_ + ind + 1) * (nInputs_ + ind);
      VecSamInputs_.setLength(kk);
      kk = nReplications_ * (nInputs_ + ind + 1) * nOutputs_;
      VecSamOutputs_.setLength(kk);
      for (ii = 0; ii < kk; ii++) 
        VecSamOutputs_[ii] = PSUADE_UNDEFINED;
      kk = nReplications_ * (nInputs_ + ind + 1);
      VecSamStates_.setLength(kk);
      for (ii = 0; ii < kk; ii++) VecSamStates_[ii] = 0;
      for (ii = 0; ii < nReplications_; ii++)
      {
        for (jj = 0; jj <= nInputs_; jj++)
        {
          for (kk = 0; kk < nInputs_; kk++)
          {
            ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
            VecSamInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                vecXT[ii*(nInputs_+1)*nInputs_+jj*nInputs_+kk];
          }
        }
        for (jj = nInputs_+1; jj < nInputs_+ind+1; jj++)
        {
          for (kk = 0; kk < nInputs_; kk++)
          {
            ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
            VecSamInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                     VecSamInputs_[ind2+nInputs_*(nInputs_+ind)+kk]; 
          }
        }
      }
      ioPtr->getParameter("input_sample", pPtr);
      for (ii = 0; ii < nReplications_; ii++)
      {
        for (jj = 0; jj <= nInputs_; jj++)
        {
          for (kk = nInputs_; kk < nInputs_+ind; kk++)
          {
            ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
            VecSamInputs_[ind2+jj*(nInputs_+ind)+kk] = 
                        pPtr.dbleArray_[ii*(ind+1)*ind+kk-nInputs_];
          }
        }
        for (jj = nInputs_+1; jj < nInputs_+ind+1; jj++)
        {
          for (kk = nInputs_; kk < nInputs_+ind; kk++)
          {
            ind2 = ii * (nInputs_ + ind + 1) * (nInputs_ + ind);
            VecSamInputs_[ind2+jj*(nInputs_+ind)+kk] = 
              pPtr.dbleArray_[ii*(ind+1)*ind+(jj-nInputs_)*ind+kk-nInputs_];
          }
        }
      }
      pPtr.clean();

      //**/ update input lower/upper bounds
      pLower.clean();
      ioPtr->getParameter("input_lbounds", pLower);
      psVector VecTV = VecILowerBs_;
      VecILowerBs_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecILowerBs_[ii] = VecTV[ii];
      for (ii = nInputs_; ii < nInputs_+ind; ii++)
        VecILowerBs_[ii] = pLower.dbleArray_[ii-nInputs_];
      pLower.clean();

      pUpper.clean();
      ioPtr->getParameter("input_ubounds", pUpper);
      VecTV = VecIUpperBs_;
      VecIUpperBs_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecIUpperBs_[ii] = VecTV[ii];
      for (ii = nInputs_; ii < nInputs_+ind; ii++)
        VecIUpperBs_[ii] = pUpper.dbleArray_[ii-nInputs_];
      pUpper.clean();

      delete ioPtr;
      nSamples_ = (nInputs_ + ind + 1) * nReplications_;

      //**/ update input distributions
      psIVector VecTI = VecInpPDFs_;
      VecInpPDFs_.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) VecInpPDFs_[ii] = VecTI[ii];
      for (ii = nInputs_; ii < nInputs_+ind; ii++) VecInpPDFs_[ii] = 0;

      psVector VecTW = VecInpMeans_;
      VecInpMeans_.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) VecInpMeans_[ii] = VecTW[ii];
      for (ii = nInputs_; ii < nInputs_+ind; ii++) VecInpMeans_[ii] = 0;

      VecTW = VecInpStds_;
      VecInpStds_.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) VecInpStds_[ii] = VecTW[ii];
      for (ii = nInputs_; ii < nInputs_+ind; ii++) VecInpStds_[ii] = 0;

      if (inputCMat_ != NULL)
      {
        psMatrix *tmpMat = new psMatrix();
        tmpMat->setDim(nInputs_+ind,nInputs_+ind);
        for (ii = 0; ii < nInputs_; ii++)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            ddata = inputCMat_->getEntry(ii,jj);
            tmpMat->setEntry(ii,jj,ddata);
          }
        }
        for (ii = nInputs_; ii < nInputs_+ind; ii++)
          tmpMat->setEntry(ii,ii,1.0);
        delete inputCMat_;
        inputCMat_ = tmpMat;
      }

      //**/ update IO
      nInputs_ += ind;
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                VecSamInputs_.getDVector(),StrInpNames_.getStrings(), 
                VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                VecInpStds_.getDVector(), inputCMat_); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                VecSamOutputs_.getDVector(), 
                VecSamStates_.getIVector(), NULL); 
      psuadeIO_->updateMethodSection(-1, nSamples_, 1, -1, -1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("The two samples have been concatenated.\n");
      printf("The new sample has nInputs = %d\n", nInputs_);
      printf("                  nSamples = %d\n", nSamples_);
      printf("                  nOutputs = %d\n", nOutputs_);
      printf("Use 'write' to store the expanded sample to a file.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iadd +++ 
    //**/ put additional inputs from another file which has the same
    //**/ number of sample points
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iadd"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iadd: add more inputs to the loaded ");
        printf("sample from another PSUADE file\n");
        printf("Syntax: iadd <filename>\n");
        printf("where <filename> should be a PSUADE ");
        printf("data file containing additional\n");
        printf("inputs.\n");
        printf("NOTE: The sample size in the data file should ");
        printf("be the same as that in\n");
        printf("      the loaded sample.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command appends one input to the loaded sample. ");
      printf("The values of\n");
      printf("this input is to be obtained from a file in ");
      printf("PSUADE format.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }

      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("       Syntax: iadd <file>\n");
        printf("       where <file> should be a PSUADE data file.\n");
        cmdStatus = 1;
        continue;
      }
      fclose(fp);

      //**/ read the sample file
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
               dataFile);
        cmdStatus = 1;
        continue;
      }

      //**/ check that sample sizes are the same as that of 
      //**/ the loaded sample
      ioPtr->getParameter("method_nsamples", pPtr);
      ind = pPtr.intData_;
      if (ind != nSamples_)
      {
        printf("ERROR: nSamples not the same in both sets of data.\n");
        printf("       nSamples in local memory = %d\n", nSamples_);
        printf("       nSamples from file       = %d\n", ind);
        cmdStatus = 1;
        continue;
      }
      ioPtr->getParameter("input_ninputs", pPtr);
      ind = pPtr.intData_;

      //**/ add new input lower and upper bounds
      psVector VecTW = VecILowerBs_;
      ioPtr->getParameter("input_lbounds", pLower);
      VecILowerBs_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecILowerBs_[ii] = VecTW[ii];
      for (ii = 0; ii < ind; ii++) 
        VecILowerBs_[nInputs_+ii] = pLower.dbleArray_[ii];
      pLower.clean();
             
      VecTW = VecIUpperBs_;
      ioPtr->getParameter("input_ubounds", pUpper);
      VecIUpperBs_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecIUpperBs_[ii] = VecTW[ii];
      for (ii = 0; ii < ind; ii++) 
        VecIUpperBs_[nInputs_+ii] = pUpper.dbleArray_[ii];
      pUpper.clean();

      //**/ add new input pdf type
      psIVector VecTI = VecInpPDFs_;
      pData pPDFs;
      ioPtr->getParameter("input_pdfs", pPDFs);
      VecInpPDFs_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecInpPDFs_[ii] = VecTI[ii];
      for (ii = 0; ii < ind; ii++)
        VecInpPDFs_[nInputs_+ii] = pPDFs.intArray_[ii];
      pPDFs.clean();

      //**/ add new input distribution means
      pData pMeans;
      ioPtr->getParameter("input_means", pMeans);
      VecTW = VecInpMeans_;
      VecInpMeans_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecInpMeans_[ii] = VecTW[ii];
      for (ii = 0; ii < ind; ii++)
        VecInpMeans_[nInputs_+ii] = pMeans.dbleArray_[ii];
      pMeans.clean();

      //**/ add new input distribution std dev
      pData pStds;
      ioPtr->getParameter("input_stdevs", pStds);
      VecTW = VecInpStds_;
      VecInpStds_.setLength(nInputs_+ind);
      for (ii = 0; ii < nInputs_; ii++) VecInpStds_[ii] = VecTW[ii];
      for (ii = 0; ii < ind; ii++)
        VecInpStds_[nInputs_+ii] = pStds.dbleArray_[ii];
      pStds.clean();

      //**/ add new input distribution sample files
      if (SamPDFFiles_ != NULL)
      {
        names = SamPDFFiles_;
        SamPDFFiles_ = new char*[nInputs_+ind];
        for (ii = 0; ii < nInputs_; ii++)
          SamPDFFiles_[ii] = names[ii];
        for (ii = 0; ii < ind; ii++)
        {
          SamPDFFiles_[nInputs_+ii] = new char[1000];
          strcpy(SamPDFFiles_[nInputs_+ii], "NONE");
        }
      }

      //**/ update correlation matrix from new inputs
      if (inputCMat_ != NULL)
      {
        ioPtr->getParameter("input_cor_matrix", pPtr);
        psMatrix *tmpMat1 = (psMatrix *) pPtr.psObject_;
        psMatrix *tmpMat2 = new psMatrix();
        tmpMat2->setDim(nInputs_+ind,nInputs_+ind);
        for (ii = 0; ii < nInputs_; ii++)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            ddata = inputCMat_->getEntry(ii,jj);
            tmpMat2->setEntry(ii,jj,ddata);
          }
        }
        for (ii = nInputs_; ii < nInputs_+ind; ii++)
        {
          for (jj = nInputs_; jj < nInputs_+ind; jj++)
          {
            ddata = tmpMat1->getEntry(ii-nInputs_,jj-nInputs_);
            tmpMat2->setEntry(ii,jj,ddata);
          }
        }
        delete inputCMat_;
        inputCMat_ = tmpMat2;
      }

      //**/ update input names
      StrInpNames_.addMoreStrings(ind);
      pINames.clean();
      ioPtr->getParameter("input_names", pINames);
      for (ii = 0; ii < ind; ii++)
        StrInpNames_.loadOneString(nInputs_+ii,pINames.strArray_[ii]);
      pINames.clean();

      //**/ update sample inputs
      ioPtr->getParameter("input_sample", pInputs);
      double *dPtrW = pInputs.dbleArray_;
      psVector vecXT = VecSamInputs_;
      VecSamInputs_.setLength(nSamples_*(nInputs_+ind));
      kk = 0;
      for (ii = 0; ii < nSamples_; ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
          VecSamInputs_[kk++] = vecXT[ii*nInputs_+jj];
        for (jj = 0; jj < ind; jj++)
          VecSamInputs_[kk++] = dPtrW[ii*ind+jj];
      }
      nInputs_ += ind;
      pInputs.clean();

      //**/ update input section 
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
            VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
            VecSamInputs_.getDVector(),StrInpNames_.getStrings(), 
            VecInpPDFs_.getIVector(),VecInpMeans_.getDVector(),
            VecInpStds_.getDVector(), inputCMat_); 
      cmdStatus = 0;
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      fflush(stdout);
      if (ioPtr != NULL) delete ioPtr;
      printf("iadd completed. Use 'write' to store the ");
      printf("modified sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ oadd +++ 
    //**/ put additional outputs from another file which has the same
    //**/ sample points
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oadd"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oadd: add more outputs to the loaded sample ");
        printf("from another PSUADE file\n");
        printf("Syntax: oadd <filename>\n");
        printf("where <filename> should be a PSUADE data file ");
        printf("containing additional\n");
        printf("outputs.\n\n");
        printf("NOTE: The sample size in the data file can be ");
        printf("different from that in\n");
        printf("      the loaded sample, in which case ");
        printf("the outputs in the sample file\n");
        printf("      are interpolated on the loaded sample ");
        printf("using response surfaces.");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command appends one output to the loaded sample. ");
      printf("The values\n");
      printf("of this output are obtained from another file in ");
      printf("PSUADE format.\n");
      printf("In case the sample in the second file is not the same ");
      printf("as the loaded\n");
      printf("sample, a response surface will be built using the ");
      printf("second sample and\n");
      printf("interpolated onto the loaded sample.\n");
      printf("Also, nInputs in the second file can be different ");
      printf("from (but has to be\n");
      printf("less than) nInputs in the loaded sample, in which ");
      printf("case input matching\n");
      printf("information will be needed.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      //**/ check whether a sample has been loaded
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ check that the data file exists
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("       Syntax: oadd <file>\n");
        printf("       where <file> should be a PSUADE data file.\n");
        cmdStatus = 1;
        continue;
      }
      else fclose(fp);

      //**/ read the data file
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
                 dataFile);
        cmdStatus = 1;
        continue;
      }

      //**/ get data file information
      ioPtr->getParameter("method_nsamples", pPtr);
      int nSamps2 = pPtr.intData_;
      ioPtr->getParameter("input_ninputs", pPtr);
      int nInps2 = pPtr.intData_;
      ioPtr->getParameter("output_noutputs", pPtr);
      int nOuts2 = pPtr.intData_;
      pData pInpPtr, pOutPtr;
      ioPtr->getParameter("input_sample", pInpPtr);
      double *samInps2 = pInpPtr.dbleArray_;
      ioPtr->getParameter("output_sample", pOutPtr);
      double *samOuts2 = pOutPtr.dbleArray_;
   
      //**/ see if the two sample inputs are the same
      int needRS=0;
      if (nInputs_ != nInps2)
      {
        printf("ERROR: nInputs in the data file != %d\n",nInputs_);
        printf("If you just want to modify the outputs and ");
        printf("not inputs, use oreplace.\n");
        delete ioPtr;
        cmdStatus = 1;
        continue;
      }
      if (nSamps2 != nSamples_) needRS = 1;
      else
      {
        for (ss = 0; ss < nSamples_; ss++)
          for (ii = 0; ii < nInputs_; ii++)
            if (VecSamInputs_[ss*nInputs_+ii] != 
                samInps2[ss*nInputs_+ii])
              break;
        if (ss != nSamples_ && ii != nInputs_) needRS = 1;
      }

      //**/ if the samples are the same, just add more outputs
      if (needRS == 0)
      {
        psVector vecYT;
        vecYT = VecSamOutputs_;
        VecSamOutputs_.setLength(nSamples_ * (nOutputs_ + nOuts2));
        for (ss = 0; ss < nSamples_; ss++)
        {
          for (ii = 0; ii < nOutputs_; ii++)
            VecSamOutputs_[ss*(nOutputs_+nOuts2)+ii] = 
              vecYT[ss*nOutputs_+ii]; 
          for (ii = 0; ii < nOuts2; ii++)
            VecSamOutputs_[ss*(nOutputs_+nOuts2)+nOutputs_+ii] = 
              samOuts2[ss*nOuts2+ii]; 
        }
        StrOutNames_.addMoreStrings(nOuts2);
        pData pONames;
        ioPtr->getParameter("output_names", pONames);
        for (ii = 0; ii < nOuts2; ii++)
          StrOutNames_.loadOneString(nOutputs_+ii,
                                     pONames.strArray_[ii]);
        pONames.clean();
        nOutputs_ += nOuts2; 
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
               VecSamOutputs_.getDVector(),VecSamStates_.getIVector(),
               StrOutNames_.getStrings());
        cmdStatus = 0;
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        delete ioPtr;
        printf("oadd successful: use write to store the ");
        printf("modified sample.\n");
        continue;
      }
      else if (nOuts2 != 1)
      {
        printf("ERROR: Since the loaded and incoming sample ");
        printf("inputs are different,\n");
        printf("       interpolation is needed. In this case, ");
        printf("the data sample must\n");
        printf("       have nOutputs=1.\n");
        printf("       But, the incoming sample has nOutputs = %d\n", 
               nOuts2);
        delete ioPtr;
        cmdStatus = 1;
        continue;
      }

      //**/ otherwise, interpolate
      psIVector vecInds;
      if (nInputs_ > nInps2)
      {
        printf("The incoming sample has less inputs than the ");
        printf("loaded sample. As such,\n");
        printf("you need to specify which inputs in the local ");
        printf("sample correspond to\n");
        printf("which inputs in the data file so response surface ");
        printf("interpolation makes\n");
        printf("sense.\n");
        vecInds.setLength(nInps2);
        for (ii = 0; ii < nInps2; ii++)
        {
          sprintf(pString,
             "Which local input corresponds to data file sample input %d ? ",
             ii+1);
          kk = getInt(1, nInputs_, pString);
          vecInds[ii] = kk - 1;
        }
      }
 
      //**/ create response surface for data file sample and evaluate
      //**/ the sample in local memory 
      printf("** Creating response surface for the second sample.\n");
      faType = -1;
      faFlag = 3;
      FuncApprox *faPtr = genFAInteractive(ioPtr, faFlag);
      if (faPtr == NULL)
      {
        printf("ERROR: cannot create response surface.\n");
        cmdStatus = 1;
        delete ioPtr;
        continue;
      }
      faPtr->setOutputLevel(outputLevel_);
      psVector vecXT, vecYT;
      if (nInputs_ == nInps2) vecXT = VecSamInputs_;
      else
      {
        vecXT.setLength(nInps2*nSamples_);
        for (ss = 0; ss < nSamples_; ss++)
          for (ii = 0; ii < nInps2; ii++)
            vecXT[ss*nInps2+ii] = 
                 VecSamInputs_[ss*nInputs_+vecInds[ii]];
      }
      vecYT.setLength(nSamples_);
      faPtr->evaluatePoint(nSamples_, vecXT.getDVector(),
                           vecYT.getDVector());

      //**/ add a new output to VecSamOutputs_
      vecXT = VecSamOutputs_;
      VecSamOutputs_.setLength(nSamples_ * (nOutputs_ + 1));
      for (ss = 0; ss < nSamples_; ss++)
      {
        for (ii = 0; ii < nOutputs_; ii++)
          VecSamOutputs_[ss*(nOutputs_+1)+ii] = vecXT[ss*nOutputs_+ii]; 
        VecSamOutputs_[ss*(nOutputs_+1)+nOutputs_] = vecYT[ss]; 
      }
      StrOutNames_.addMoreStrings(nOuts2);
      pData pONames;
      ioPtr->getParameter("output_names", pONames);
      StrOutNames_.loadOneString(nOutputs_,pONames.strArray_[0]);
      pONames.clean();
      nOutputs_ += nOuts2; 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                    VecSamOutputs_.getDVector(),VecSamStates_.getIVector(),
                    StrOutNames_.getStrings());
      cmdStatus = 0;
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
      delete ioPtr;
      delete faPtr;
      printf("oadd completed. Use 'write' to store the ");
      printf("modified sample.\n");
      continue;
    }

    //**/ =============================================================
    // commands that manipulate local memory
    //**/ =============================================================
    //**/ -------------------------------------------------------------
    // +++ iadd1 +++ 
    //**/ add an input 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iadd1"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iadd: add one more input to the existing sample\n");
        printf("      The values of the new input may be set ");
        printf("to a constant, or they\n");
        printf("      may be drawn randomly between 0 and 1.\n");
        printf("Syntax: iadd1 \n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command appends one input to the loaded sample. ");
      printf("The values of\n");
      printf("this input are generated uniformly in [0,1], or ");
      printf("set to a constant.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }

      printf("How to initialize this additional input ? \n");
      printf("1. Set it to a constant, or\n");
      printf("2. Set it to be uniformly distributed in [0,1].\n");
      sprintf(pString,"Select an option (1 or 2) : ");
      int setOption = getInt(1, 2, pString);
      double setVal;
      if (setOption == 1)
      {
        sprintf(pString,"What is the constant ? ");
        setVal = getDouble(pString);
      }

      //**/ adjust the lower/upper bounds for the additional input
      psVector VecTW = VecILowerBs_;
      VecILowerBs_.setLength(nInputs_+1);
      for (ii = 0; ii < nInputs_; ii++) VecILowerBs_[ii] = VecTW[ii];
      if (setOption == 1) VecILowerBs_[nInputs_] = setVal - 1e-6;
      else                VecILowerBs_[nInputs_] = 0;

      VecTW = VecIUpperBs_;
      VecIUpperBs_.setLength(nInputs_+1);
      for (ii = 0; ii < nInputs_; ii++) VecIUpperBs_[ii] = VecTW[ii];
      if (setOption == 1) VecIUpperBs_[nInputs_] = setVal + 1e-6;
      else                VecIUpperBs_[nInputs_] = 1;

      //**/ adjust the PDF information for the additional input
      if (VecInpPDFs_.length() > 0)
      {
        psIVector VecTI = VecInpPDFs_;
        VecInpPDFs_.setLength(nInputs_+1);
        for (jj = 0; jj < nInputs_; jj++) VecInpPDFs_[jj] = VecTI[jj];
        VecInpPDFs_[nInputs_] = 0;
      }
      if (VecInpMeans_.length() > 0)
      {
        psVector VecTW = VecInpMeans_;
        VecInpMeans_.setLength(nInputs_+1);
        for (jj = 0; jj < nInputs_; jj++) VecInpMeans_[jj] = VecTW[jj];
        VecInpMeans_[nInputs_] = 0.0;
      }
      if (VecInpStds_.length() > 0)
      {
        psVector VecTW = VecInpStds_;
        VecInpStds_.setLength(nInputs_+1);
        for (jj = 0; jj < nInputs_; jj++) VecInpStds_[jj] = VecTW[jj];
        VecInpStds_[nInputs_] = 0.0;
      }
      if (inputCMat_ != NULL)
      {
        psMatrix *tmpMat2 = new psMatrix();
        tmpMat2->setDim(nInputs_+1,nInputs_+1);
        for (ii = 0; ii < nInputs_; ii++)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            ddata = inputCMat_->getEntry(ii,jj);
            tmpMat2->setEntry(ii,jj,ddata);
          }
        }
        ddata = 1.0;
        tmpMat2->setEntry(nInputs_,nInputs_,ddata);
        delete inputCMat_;
        inputCMat_ = tmpMat2;
      }
      if (SamPDFFiles_ != NULL)
      {
        names = SamPDFFiles_;
        SamPDFFiles_ = new char*[nInputs_+1];
        for (ii = 0; ii < nInputs_; ii++)
          SamPDFFiles_[ii] = names[ii];
        SamPDFFiles_[nInputs_] = new char[1000];
        strcpy(SamPDFFiles_[nInputs_], "NONE");
      }

      //**/ add one more input names
      StrInpNames_.addMoreStrings(iOne);
      sprintf(pString, "X%d", nInputs_+1);
      StrInpNames_.loadOneString(nInputs_, pString);

      //**/ add more sample inputs 
      vecXT = VecSamInputs_;
      VecSamInputs_.setLength(nSamples_*(nInputs_+1));
      for (ii = 0; ii < nSamples_; ii++)
      {
        for (jj = 0; jj < nInputs_; jj++)
          VecSamInputs_[ii*(nInputs_+1)+jj] = vecXT[ii*nInputs_+jj];
        if (setOption == 1)
          VecSamInputs_[ii*(nInputs_+1)+nInputs_] = setVal;
        else
          VecSamInputs_[ii*(nInputs_+1)+nInputs_] = PSUADE_drand();
      }
      nInputs_++;

      //**/ update the input section
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                   VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                   VecInpStds_.getDVector(), inputCMat_); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      fflush(stdout);
      printf("iadd1 : one input added.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oadd1 +++
    //**/ add an output 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oadd1"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oadd1: add one output to the existing sample.\n");
        printf("       The values of the new output are set ");
        printf("to be UNDEFINED.\n");
        printf("Syntax: oadd1 \n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command appends one output to the loaded sample. ");
      printf("The values of\n");
      printf("this output will be set to UNDEFINED (1e35)\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }

      if (VecSamOutputs_.length() == 0)
      {
        printf("ERROR: No existing sample outputs in the ");
        printf("loaded sample.\n");
        printf("       This command is not to be executed.\n");
        cmdStatus = 1;
        continue;
      }
      else
      {
        StrOutNames_.addMoreStrings(iOne);
        sprintf(pString, "Y%d", nOutputs_+1);
        StrOutNames_.loadOneString(nOutputs_, pString);
        vecYT = VecSamOutputs_;
        VecSamOutputs_.setLength(nSamples_*(nOutputs_+1));
        for (ii = 0; ii < nSamples_; ii++)
        {
          for (jj = 0; jj < nOutputs_; jj++)
            VecSamOutputs_[ii*(nOutputs_+1)+jj] = vecYT[ii*nOutputs_+jj];
          VecSamOutputs_[ii*(nOutputs_+1)+nOutputs_] = PSUADE_UNDEFINED;
        }
        nOutputs_++;
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
             VecSamOutputs_.getDVector(),VecSamStates_.getIVector(), 
             StrOutNames_.getStrings()); 
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        printf("oadd1 : 1 output added.\n");
        fflush(stdout);
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ ireplace +++ 
    //**/ replace one existing input with an input from a given file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ireplace"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ireplace: replace an input in the existing sample ");
        printf("from another sample\n");
        printf("          file in ASCII format.\n");
        printf("Syntax: ireplace <filename>\n");
        printf("where <filename> should be in the format given below: \n");
        printf("line 1: PSUADE_BEGIN (optional)\n");
        printf("line 2: <nSamples>\n");
        printf("line 3: input value for sample point 1\n");
        printf("line 4: input value for sample point 2\n");
        printf(".......\n");
        printf("last line: PSUADE_END\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      status = 0;
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("Use <ireplace -h> to see command syntax\n");
        status = 1;
        cmdStatus = 1;
        continue;
      }
      else fclose(fp);
      if (status == 0)
      {
        fp = fopen(dataFile, "r");
        fscanf(fp, "%s", winput);
        if (strcmp(winput, "PSUADE_BEGIN"))
        {
          printf("INFO: optional PSUADE_BEGIN keyword not found.\n");
          printf("      First line should have nSamples\n");
          fclose(fp);
          fp = fopen(dataFile, "r");
        }
        fscanf(fp, "%d", &kk);
        if (kk != nSamples_)
        {
          fclose(fp);
          printf("ERROR: File and local data do not match.\n");
          printf("     size of sample in local memory = %d\n",nSamples_);
          printf("     nSamples from external file    = %d\n",kk);
          status = 1;
          cmdStatus = 1;
          continue;
        }
      }
      if (status == 0)
      {
        inputID = 0;
        if (nInputs_ > 1)
        {
          sprintf(pString,"Enter input number to replace (1 - %d) : ",
                  nInputs_);
          inputID = getInt(1, nInputs_, pString);
          inputID--;
        }
        for (ii = 0; ii < nSamples_; ii++)
        {
          fscanf(fp, "%lg", &ddata);
          VecSamInputs_[ii*nInputs_+inputID] = ddata;
        }
        Xmin = PSUADE_UNDEFINED;
        Xmax = - Xmin;
        for (ii = 0; ii < nSamples_; ii++)
        {
          if (VecSamInputs_[ii*nInputs_+inputID] < Xmin)
            Xmin = VecSamInputs_[ii*nInputs_+inputID];
          else if (VecSamInputs_[ii*nInputs_+inputID] > Xmax)
            Xmax = VecSamInputs_[ii*nInputs_+inputID];
        }
        VecILowerBs_[inputID] = Xmin; 
        VecIUpperBs_[inputID] = Xmax; 
        VecInpPDFs_[inputID] = 0;
        VecInpMeans_[inputID] = 0;
        VecInpStds_[inputID] = 0;
        if (inputCMat_ != NULL)
        {
          for (ii = 0; ii < nInputs_; ii++)
            inputCMat_->setEntry(inputID,ii,0);
          for (ii = 0; ii < nInputs_; ii++)
            inputCMat_->setEntry(ii,inputID,0);
          inputCMat_->setEntry(inputID,inputID,1);
        }

        //**/ add one more input name and then update input section
        sprintf(pString, "X%d", inputID+1);
        StrInpNames_.loadOneString(inputID, pString);
        psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                     VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                     VecSamInputs_.getDVector(),StrInpNames_.getStrings(), 
                     VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                     VecInpStds_.getDVector(), inputCMat_); 
        fscanf(fp, "%s", winput);
        if (strcmp(winput, "PSUADE_END"))
          printf("INFO: the optional PSUADE_END keyword not found.\n");

        fclose(fp);
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        cmdStatus = 0;
        VecTagArray_.clean();
      }
      fflush(stdout);
    }

    //**/ -------------------------------------------------------------
    // +++ oreplace +++ 
    //**/ replace all existing output with outputs from a given file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oreplace"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oreplace: replace all outputs in the existing ");
        printf("sample from another\n");
        printf("          another sample file in ASCII format.\n");
        printf("Syntax: oreplace <file>\n");
        printf("where <file> should be in the format given below: \n");
        printf("line 1: PSUADE_BEGIN (optional)\n");
        printf("line 2: <nSamples> <nOutputs>\n");
        printf("line 3: output value for sample point 1\n");
        printf("line 4: output value for sample point 2\n");
        printf(".......\n");
        printf("last line: PSUADE_END\n");
        continue;
      }
      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: Need to load a sample first.\n");
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sscanf(lineIn, "%s %s", command, dataFile);
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: File %s not found.\n", dataFile);
        printf("Syntax: oreplace <sample file>\n");
        cmdStatus = 1;
        continue;
      }
      fclose(fp);
      int nOut2;
      fp = fopen(dataFile, "r");
      fscanf(fp, "%s", winput);
      if (strcmp(winput, "PSUADE_BEGIN"))
      {
        printf("INFO: The optional PSUADE_BEGIN keyword is not found\n");
        printf("      Make sure that the first line has ");
        printf("nSamples and nOutputs.\n");
        fclose(fp);
        fp = fopen(dataFile, "r");
      }
      fscanf(fp, "%d %d", &kk, &nOut2);
      if (kk != nSamples_ || nOut2 <= 0)
      {
        printf("ERROR: File and local parameters do not match.\n");
        printf("       nSamples (local) = %d\n", nSamples_);
        printf("       nSamples (file)  = %d\n", kk);
        printf("==> nSamples (local) should be equal ");
        printf("to nSamples (file).\n");
        printf("Advice: Check your file format. Make sure ");
        printf("that the first line has\n");
        printf("        nSamples and nOutputs.\n");
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      printf("This command only checks that nSamples of the ");
      printf("output data file matches\n");
      printf("that of the loaded sample, so make sure the ");
      printf("sample data are correctly\n");
      printf("ordered.\n");
      printf("       nSamples (local) = %d\n", nSamples_);
      printf("       nOutputs (file)  = %d\n", nOut2);
      if (nOutputs_ != nOut2)
      {
        nOutputs_ = nOut2;
        StrOutNames_.setNumStrings(nOutputs_);
        for (ii = 0; ii < nOutputs_; ii++)
        {
          sprintf(pString, "Y%d", ii+1);
          StrOutNames_.loadOneString(ii, pString);
        }
        VecSamOutputs_.clean();
      }
      if (VecSamOutputs_.length() == 0)
        VecSamOutputs_.setLength(nSamples_*nOutputs_); 
      for (ii = 0; ii < nSamples_*nOutputs_; ii++)
      {
        fscanf(fp, "%lg", &ddata);
        VecSamOutputs_[ii] = ddata;
      }
      for (ii = 0; ii < nSamples_; ii++) VecSamStates_[ii] = 1;
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
            VecSamOutputs_.getDVector(),
            VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      fscanf(fp, "%s", winput);
      if (strcmp(winput, "PSUADE_END"))
        printf("INFO: The optional PSUADE_END keyword is not found.\n");

      fclose(fp);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
      VecTagArray_.clean();
      checkResidentSampleOutputs(); 
      fflush(stdout);
    }

    //**/ -------------------------------------------------------------
    // +++ ssplit (splitsample)
    //**/ split the file into 2 separate data file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command,"splitsample") || !strcmp(command,"ssplit"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ssplit: split the data set into 2 subsets.\n");
        printf("Syntax: ssplit (no argument needed)\n");
        printf("NOTE: the splitted samples will be stored in ");
        printf("psuadeSample1 and\n");
        printf("      psuadeSample2.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command splits the loaded sample inputs ");
      printf("into 2 sub-samples.\n");
      printf("A few splitting options are available.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      printf("The current sample size is %d.\n", nSamples_);
      sprintf(pString,"Sample size of the first set? (1 - %d) ",
              nSamples_-1);
      int splitCnt = getInt(1, nSamples_, pString);
      printf("Select from the following options: \n");
      printf("1. random draw from the sample\n");
      printf("2. draw the first %d sample points\n",splitCnt);
      printf("3. draw every other %d sample points\n",nSamples_/splitCnt);
      printf("4. draw the first m points from every chunk of M points\n");
      sprintf(pString,"your choice: (1 - 4) ");
      int splitChoice = getInt(1, nSamples_, pString);
      if (splitChoice == 3 && (nSamples_/splitCnt*splitCnt != nSamples_))
      {
        printf("ERROR: nSamples not a multiple of draw size.\n");
        cmdStatus = 1;
        continue;
      }
      if (splitChoice == 3) 
      {
        sprintf(pString,"Starting sample number: (1 - %d) ",
                nSamples_/splitCnt);
        ll = getInt(1, nSamples_/splitCnt, pString);
        ll--;
      }
      int chunkSize=1, mStep;
      if (splitChoice == 4) 
      {
        sprintf(pString,"Enter m = ");
        chunkSize = getInt(1, nSamples_, pString);
        mStep = nSamples_ * chunkSize / splitCnt;
        printf("m = %d\n", chunkSize);
        printf("M = %d\n", mStep);
      }
      vecXT.setLength(splitCnt*nInputs_);
      vecYT.setLength(splitCnt*nOutputs_);
      vecST.setLength(splitCnt);
      vecTags.setLength(nSamples_);
      for (ii = 0; ii < splitCnt; ii++)
      {
        if (splitChoice == 1)
        {
          ind = PSUADE_rand() % nSamples_;
          ind2 = 0;
          while (vecTags[ind] == 1 && ind2 < 1000)
          {
            ind = PSUADE_rand() % nSamples_;
            ind2++;
          }
          if (vecTags[ind] == 1)
            for (ind = 0; ind < nSamples_; ind++)
              if (vecTags[ind] == 0) break;
          if (vecTags[ind] == 1)
          {
            printf("ERROR: cannot split data. \n");
            cmdStatus = 1;
            break;
          }
        } 
        else if (splitChoice == 2) ind = ii;
        else if (splitChoice == 3) ind = ii*nSamples_/splitCnt+ll;
        else if (splitChoice == 4)
          ind = ii / chunkSize * mStep + ii % chunkSize;
        for (jj = 0; jj < nInputs_; jj++)
          vecXT[ii*nInputs_+jj] = VecSamInputs_[ind*nInputs_+jj];
        for (jj = 0; jj < nOutputs_; jj++)
          vecYT[ii*nOutputs_+jj] = VecSamOutputs_[ind*nOutputs_+jj];
        vecTags[ind] = 1;
        vecST[ii] = VecSamStates_[ind];
      }
      ind = 0;
      for (ii = 0; ii < nSamples_; ii++)
      {
        if (vecTags[ii] == 0)
        {
          for (jj = 0; jj < nInputs_; jj++)
            VecSamInputs_[ind*nInputs_+jj] =
                          VecSamInputs_[ii*nInputs_+jj];
          for (jj = 0; jj < nOutputs_; jj++)
            VecSamOutputs_[ind*nOutputs_+jj] = 
                          VecSamOutputs_[ii*nOutputs_+jj];
          VecSamStates_[ind] = VecSamStates_[ii]; 
          ind++;
        }
      }
      psuadeIO_->updateInputSection(splitCnt,nInputs_,NULL,NULL,NULL,
                      vecXT.getDVector(),NULL,NULL,NULL,NULL,NULL);
      psuadeIO_->updateOutputSection(splitCnt,nOutputs_,vecYT.getDVector(),
                        vecST.getIVector(),StrOutNames_.getStrings()); 
      psuadeIO_->updateMethodSection(-1,splitCnt,-1,-1,-1);
      psuadeIO_->writePsuadeFile("psuadeSample1",0);

      psuadeIO_->updateInputSection(nSamples_-splitCnt,nInputs_,NULL,
                       NULL, NULL,VecSamInputs_.getDVector(),NULL, 
                       NULL, NULL, NULL, NULL);
      psuadeIO_->updateOutputSection(nSamples_-splitCnt,nOutputs_,
                   VecSamOutputs_.getDVector(),
                   VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      psuadeIO_->updateMethodSection(-1,nSamples_-splitCnt,-1,-1,-1);
      psuadeIO_->writePsuadeFile("psuadeSample2",0);
      if (currSession != NULL) delete currSession;
      currSession = NULL;

      printf("The 2 data files are in psuadeSample1 and psuadeSample2.\n");
      printf("The loaded sample in the workspace have been erased ");
      printf("so you need to\n");
      printf("re-load to restore the original sample.\n");
      fflush(stdout);
      VecSamInputs_.clean();
      VecSamOutputs_.clean();
      VecSamStates_.clean();
      VecILowerBs_.clean();
      VecIUpperBs_.clean();
      VecInpPDFs_.clean();
      VecInpMeans_.clean();
      VecInpStds_.clean();
      StrInpNames_.clean();
      StrOutNames_.clean();
      nSamples_ = 0;
      nInputs_ = 0;
      nOutputs_ = 0;
      cmdStatus = 0;
      VecTagArray_.clean();
    }

    //**/ =============================================================
    // UQ/SA commands
    //**/ =============================================================
    //**/ -------------------------------------------------------------
    // +++ uq 
    //**/ uncertainty analysis + output to a matlab file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ua"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ua: uncertainty analysis (compute moments).\n");
        printf("Syntax: ua (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command computes the basic statistics (mean, ");
      printf("standard deviation,\n");
      printf("skewness, and kurtosis) of a selected output ");
      printf("in the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(0, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_MOMENT;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sem 
    //**/ compute standard of mean
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sem"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sem: standard error of mean analysis.\n");
        printf("Syntax: sem (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command computes directly from the ");
      printf("loaded sample the standard\n");
      printf("error of mean using bootstrapping (Note: ");
      printf("the standard error of mean\n");
      printf("from 'ua' uses a formula and not bootstrapping).\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(0, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_BSTRAP;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ca 
    //**/ correlation analysis 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ca"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ca: correlation analysis\n");
        printf("Syntax: ca (after data have been loaded)\n");
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command computes correlation information ");
      printf("between inputs and a\n");
      printf("selected output of the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ > 1)
	 printf("INFO: This analysis currently is for nInputs = 1\n");

      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_CORRELATION;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ anova 
    //**/ analysis of variation
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "anova"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("anova: analysis of variation\n");
        printf("Syntax: anova (after data have been loaded)\n");
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command performs analysis of variation (ANOVA) ");
      printf("on a selected\n");
      printf("output of the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_ANOVA;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ moat 
    //**/ Morris' std deviation versus revised mean plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moat"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moat: Morris screening analysis\n");
        printf("Syntax: moat (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs a Morris-one-at-a-time ");
      printf("(MOAT) analysis on a\n");
      printf("selected output of the loaded sample. The loaded ");
      printf("sample needs to be\n");
      printf("a Morris sample. The MOAT method is for coarse ");
      printf("parameter ranking.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_MOAT;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ moatmo 
    //**/ Morris' modified means for multiple outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "moatmo"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("moatmo: Morris screening analysis for multiple\n");
        printf("        outputs simultaneously\n");
        printf("Syntax: moatmo (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs a Morris-one-at-a-time ");
      printf("(MOAT) analysis on ALL\n");
      printf("outputs of the loaded sample. The multi-output ");
      printf("MOAT heat map gives a\n");
      printf("global view of parameter importance for all ");
      printf("outputs. The loaded sample\n");
      printf("needs to be a Morris sample. The MOAT method ");
      printf("is for coarse parameter\n");
      printf("ranking.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || VecSamOutputs_.length() == 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ <= 1)
      {
        printf("INFO: only one output -> use moat instead.\n");
        cmdStatus = 1;
        continue;
      } 
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      int diagSave = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      int analysisMethod = PSUADE_ANA_MOAT;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      vecWT.setLength(nInputs_*nOutputs_);
      if (psConfig_.AnaExpertModeIsOn())
      {
        printf("INFO: Analysis expert mode will be turned off.\n");
        printf("      This command will only create modified means ");
        printf("plots. If screen or\n");
        printf("      scatter plots are also needed, use 'moat' ");
        printf("for each output.\n");
        psConfig_.AnaExpertModeSaveAndReset();
        psConfig_.AnaExpertModeOff();
      }
      for (ii = 0; ii < nOutputs_; ii++)
      {
        anaManager->analyze(psuadeIO_, 0, NULL, ii);
        pData *pdata = psuadeIO_->getAuxData();
        if (pdata->nDbles_ == nInputs_)
        {
          for (jj = 0; jj < nInputs_; jj++) 
            vecWT[ii*nInputs_+jj] = pdata->dbleArray_[jj];
          pdata->clean();
        }
      }
      psuadeIO_->updateAnalysisSection(-1,-1,-1,diagSave,-1, -1);
      psConfig_.AnaExpertModeRestore();
      delete anaManager;
      fp = NULL;
      if (plotMatlab()) fp = fopen("matlabmoatmo.m", "w");
      if (fp != NULL)
      {
        sprintf(pString,"This file contains Morris modified Means");
        fwriteComment(fp, pString);
        fprintf(fp, "A = [\n");
        for (ii = 0; ii < nInputs_; ii++)
        {
          for (jj = 0; jj < nOutputs_; jj++)
            fprintf(fp,"%16.8e ", vecWT[jj*nInputs_+ii]);
          fprintf(fp, "\n");
        }
        fprintf(fp, "];\n");
        fprintf(fp, "A2 = A * inv(diag(max(A)));\n");

        fprintf(fp, "  XStr = {");
        for (ii = 0; ii < nInputs_-1; ii++)
        {
          if (StrInpNames_[ii] != NULL)
               fprintf(fp,"'%s',",StrInpNames_[ii]);
          else fprintf(fp,"'X%d',",ii+1);
        }
        if (StrInpNames_[nInputs_-1] != NULL)
             fprintf(fp,"'%s'};\n",StrInpNames_[nInputs_-1]);
        else fprintf(fp,"'X%d'};\n",nInputs_);

        fprintf(fp, "  YStr = {");
        for (ii = 0; ii < nOutputs_-1; ii++) fprintf(fp,"'%d',",ii+1);
        fprintf(fp,"'%d'};\n",nOutputs_);
        fwriteHold(fp, 0);
        fprintf(fp,"nn = %d;\n", nInputs_);
        fprintf(fp,"mm = %d;\n", nOutputs_);
        fprintf(fp,"X = 0.5 : 1 : nn-0.5;\n");
        fprintf(fp,"Y = 0.5 : 1 : mm-0.5;\n");
        fprintf(fp, "imagesc(X,Y,A2');\n");
        fprintf(fp, "axis([0 nn 0 mm])\n");
        fwriteHold(fp, 1);
        fprintf(fp,"XX = 0 : nn;\n");
        fprintf(fp,"for ii = 1 : mm-1\n");
        fprintf(fp,"  plot(XX,ii*ones(nn+1,1),'linewidth',2)\n");
        fprintf(fp,"end;\n");
        fprintf(fp,"YY = 0 : mm;\n");
        fprintf(fp,"for ii = 1 : nn-1\n");
        fprintf(fp,"  plot(ii*ones(mm+1,1),YY,'linewidth',2)\n");
        fprintf(fp,"end;\n");
        fprintf(fp,"set(gca,'XTickLabel',[]);\n");
        fprintf(fp,"set(gca,'YTickLabel',[]);\n");
        fprintf(fp,"th=text((1:nn)-0.5,repmat(mm+0.05*mm,nn,1),XStr,");
        fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
        fprintf(fp,"set(th, 'fontsize', 12)\n");
        fprintf(fp,"set(th, 'fontweight', 'bold')\n");
        fprintf(fp,"th=text(repmat(-0.05,mm,1),(1:mm)-0.5,YStr,");
        fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
        fprintf(fp,"set(th, 'fontsize', 12)\n");
        fprintf(fp,"set(th, 'fontweight', 'bold')\n");
        fwritePlotAxes(fp);
        fprintf(fp,"colorbar\n");
        fwritePlotTitle(fp,"Morris Relative Importance Measure");
        fwritePlotXLabel(fp, "Inputs");
        fwritePlotYLabel(fp, "Outputs");
        fprintf(fp,"disp('The colors denote relative magnitudes')\n");
        fwriteHold(fp, 0);
        fclose(fp);
        printf("Morris plot file = matlabmoatmo.m\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ ff 
    //**/ fractional factorial first and second order effects
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ff"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ff: Fractional factorial screening analysis\n");
        printf("Syntax: ff (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs parameter screening using ");
      printf("factorial/fractional\n");
      printf("factorial sampling designs.\n");
      printf("Thus, the loaded sample has to be a factorial ");
      printf("or fractional factorial\n");
      printf("sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_FF;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ lsa 
    //**/ local sensitivity analysis fof first order effects
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "lsa"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("lsa: local sensitivity analysis\n");
        printf("Syntax: lsa (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs local gradient-based linear ");
      printf("sensitivity analysis\n");
      printf("on the selected output. Thus, the loaded sample ");
      printf("has to be a Plackett-\n");
      printf("Burman, Box-Behnken, or LSA design.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_LSA;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ mars_sa 
    //**/ MARS screening of parameters
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mars_sa"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mars_sa: MARS-based sensitivity analysis\n");
        printf("Syntax: mars_sa (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs COARSE global sensitivity ");
      printf("analysis based on the\n");
      printf("result of training the loaded sample with the ");
      printf("MARS response surface.\n");
      printf("The loaded sample may be any space-filling ");
      printf("design.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString,"MARS (0) or MARS with bagging (1) ? ");
      kk = getInt(0, 1, pString);
      if (kk == 0) faType = PSUADE_RS_MARS;
      else         faType = PSUADE_RS_MARSB;
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      if (faPtr != NULL)
      {
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        if (faType == PSUADE_RS_MARSB) faPtr->setOutputLevel(4);
        vecYT.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++)
          vecYT[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
        status = faPtr->initialize(VecSamInputs_.getDVector(),
                                   vecYT.getDVector());
        if (faType == PSUADE_RS_MARS)
        {
          strcpy(pString, "rank");
          targv[0] = (char *) pString;
          faPtr->setParams(1, targv);
        }
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ gp_sa 
    //**/ GP screening of parameters
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gp_sa"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gp_sa: Gaussian Process-based sensitivity analysis\n");
        printf("Syntax: gp_sa (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs COARSE global sensitivity ");
      printf("analysis based on the\n");
      printf("result of training the loaded sample with ");
      printf("the Gaussian process/Kriging\n");
      printf("model (that is, magnitudes of the hyperparameters). The ");
      printf("loaded sample\n");
      printf("may be any space-filling design.\n");
      printf("NOTE: this command does not work on samples with ");
      printf("replicated points.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      faType = -1;
      printf("Which Gaussian process ? \n");
#ifdef HAVE_TPROS
      printf("1. MacKay's Tpros\n");
#endif
      printf("2. Tong's GP\n");
      printf("3. Kriging\n");
      sprintf(pString, "Enter number (1, 2, or 3) = ");
      faType = getInt(1, 3, pString);
#ifdef HAVE_TPROS
      if (faType == 1) faType = PSUADE_RS_GP1;
#endif
      if (faType == 2) faType = PSUADE_RS_GP3;
      if (faType == 3) faType = PSUADE_RS_KR;
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      if (faPtr != NULL)
      {
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        //**/ use the slowest and more accurate Kriging
        //**/if (faType == PSUADE_RS_KR)
        //**/{
        //**/   strcpy(pString, "setMode2");
        //**/   targv[0] = (char *) pString;
        //**/   faPtr->setParams(1, targv);
        //**/}
        //**/ temporary turn off rs expert mode
        psConfig_.RSExpertModeSaveAndReset();
        psVector VecTY;
        VecTY.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++)
          VecTY[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
        status = faPtr->initialize(VecSamInputs_.getDVector(),
                                   VecTY.getDVector());
        strcpy(pString, "rank");
        targv[0] = (char *) pString;
        faPtr->setParams(1, targv);
        psConfig_.RSExpertModeRestore();
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ sot_sa 
    //**/ sum-of-trees screening of parameters
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sot_sa"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sot_sa: Sum-of-trees-based sensitivity analysis\n");
        printf("Syntax: sot_sa (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs COARSE global sensitivity ");
      printf("analysis based on the\n");
      printf("result of training the loaded sample with ");
      printf("the sum-of-trees model (the\n");
      printf("number of splittings). The sample may be any ");
      printf("space-filling design.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      faType = PSUADE_RS_SOTS;
      FuncApprox *faPtr  = genFA(faType, nInputs_, iOne, nSamples_);
      if (faPtr == NULL)
      {
        printf("ERROR: cannot create response surface.\n");
        cmdStatus = 1;
        continue;
      }
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      vecYT.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
        vecYT[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
      status = faPtr->initialize(VecSamInputs_.getDVector(),
                                 vecYT.getDVector());
      strcpy(pString, "mode0");
      targv[0] = (char *) pString;
      faPtr->setParams(1, targv);
      strcpy(pString, "rank");
      targv[0] = (char *) pString;
      ddata = faPtr->setParams(1, targv);
      delete faPtr;
      faPtr = NULL;
      printf("sot_sa score (sum of all std dev) = %e\n", ddata);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ me 
    //**/ quantitative main effect analysis + output to a matlab file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "me"))
    {
      printAsterisks(PL_INFO, 0);
      printf("This command takes the 'large' (tens of thousands or more) ");
      printf("sample that\n");
      printf("has been loaded and then computes the approximate main ");
      printf("effects (Sobol'\n");
      printf("first-order indices for all inputs).  It can handle ");
      printf("uncorrelated and\n");
      printf("correlated inputs (in the form of joint multivariate ");
      printf("distributions or\n");
      printf("inequality constraints, which are expected to have been ");
      printf("embedded in\n");
      printf("loaded sample.)\n");
      printf("This command operates directly on the sample, meaning ");
      printf("that no response\n");
      printf("surface is used in the process (as opposed to 'rssobol1'), ");
      printf("although\n");
      printf("the sample may have been evaluated using a response surface ");
      printf("before\n");
      printf("it is loaded.\n");
      printf("This command performs best with replicated Latin ");
      printf("hypercube samples,\n");
      printf("although it works for any random or quasi-random samples, ");
      printf("albeit a\n");
      printf("a little less accurate. Larger samples (tens to ");
      printf("hundred of thousands)\n");
      printf("should give more accurate results. For small samples ");
      printf("hundreds), the\n");
      printf("alternative is 'rssobol1' or 'rssobol1b', which ");
      printf("internally generates\n");
      printf("large samples and computes main effects via ");
      printf("response surfaces.\n");
      printf("NOTE: internal parameters for this command can be ");
      printf("changed by first\n");
      printf("      turning on ana_expert mode before calling ");
      printf("this command.\n");
      printDashes(PL_INFO, 0);
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("me: main effect analysis (variance-based)\n");
        printf("Syntax: me (after data have been loaded)\n");
        continue;
      }
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_ME;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      psuadeIO_->getAuxData()->clean(); 
      delete anaManager;
      //**/ auxData has the results too, but currently it is not used,
      //**/ so just clean it
      psuadeIO_->getAuxData()->clean();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ie 
    //**/ quantitative interaction effect analysis 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ie"))
    {
      printAsterisks(PL_INFO, 0);
      printf("This command takes the large (tens of thousands or more) ");
      printf("sample that\n");
      printf("has been loaded and computes the approximate 2-parameter ");
      printf("effects for\n");
      printf("all pairs of inputs.  It can handle uncorrelated and ");
      printf("correlated inputs\n");
      printf("alike (in the form of joint multivariate ");
      printf("distributions or in the\n");
      printf("form of inequality constraints that have been embedded ");
      printf("in the sample). \n");
      printf("NOTE: additional inequality constraints may be imposed ");
      printf("through the\n");
      printf("      'rs_consraint' feature defined in the ANALYSIS ");
      printf("section of your\n");
      printf("      sample file.\n");
      printf("This command operates directly on the sample, meaning ");
      printf("that no response\n"); 
      printf("surface will be used in the process, although the sample ");
      printf("may have been\n");
      printf("evaluated using a response surface before it is loaded.\n");
      printf("This command performs best with replicated orthogonal ");
      printf("array samples,\n");
      printf("although it works for any random or quasi-random samples, ");
      printf("albert a\n");
      printf("little less accurate. Larger samples (tens to hundreds of");
      printf(" thousands)\n");
      printf("should give more accurate results.\n");
      printf("For small samples (hundreds to thousands), the alternative ");
      printf("is to use\n");
      printf("'rssobol2' or 'rssobol2b', which internally creates ");
      printf("large samples and\n");
      printf("then compute the interaction effects using some ");
      printf("user-selected response\n");
      printf("surfaces.\n");
      printf("NOTE: internal parameters for this command can be ");
      printf("changed by first\n");
      printf("      turning on ana_expert mode before calling this command.\n");
      printDashes(PL_INFO, 0);
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ie: 2-way interaction effect analysis (variance-based)\n");
        printf("Syntax: ie (after data have been loaded)\n");
        continue;
      }
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ <= 2)
      {
        printf("INFO: There is no point doing this for nInputs <= 2\n");
        printf("      since two-parameter effect = total variance\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_IE;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      //**/ auxData has the results too, but currently it is not used,
      //**/ so just clean it
      psuadeIO_->getAuxData()->clean();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ tsi 
    //**/ quantitative total sensitivity analysis 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "tsi"))
    {
      printAsterisks(PL_INFO, 0);
      printf("This command takes the 'large' (tens of thousands or more) ");
      printf("sample that\n");
      printf("has been loaded and then computes the approximate total ");
      printf("effects (Sobol'\n");
      printf("total-order indices for all inputs).  It can handle ");
      printf("uncorrelated and\n");
      printf("correlated inputs (in the form of joint multivariate ");
      printf("distributions or\n");
      printf("inequality constraints, which are expected to have been ");
      printf("embedded in\n");
      printf("loaded sample.)\n");
      printf("This command operates directly on the sample, meaning ");
      printf("that no response\n");
      printf("surface is used in the process (as opposed to 'rssoboltsi'), ");
      printf("although\n");
      printf("the sample may have been evaluated using a response surface ");
      printf("before\n");
      printf("it is loaded. This command is intended for samples with ");
      printf("a relatively\n");
      printf("small number of inputs (<= 10).\n");
      printf("This command works for any space-filling samples. ");
      printf("Larger samples (tens\n");
      printf("to hundreds of thousands) should give more accurate ");
      printf("results. For small\n");
      printf("samples (hundreds), the alternative is 'rssoboltsi1' ");
      printf("or 'rssobol1b',\n");
      printf("which internally creates large samples and computes total ");
      printf("effects via\n");
      printf("response surfaces.\n");
      printf("NOTE: internal parameters for this command can be ");
      printf("changed by first\n");
      printf("      turning on ana_expert mode before calling this command.\n");
      printDashes(PL_INFO, 0);
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("tsi: total sensitivity analysis (variance-based)\n");
        printf("     (suitable for raw data)\n");
        printf("Syntax: tsi (after data have been loaded)\n");
        continue;
      }
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ >= 20 || nSamples_ < 50*nInputs_) 
      {
        printf("This command is not recommended for small sample or\n");
        printf("large number of inputs.\n");
        printf("Use at most 20 inputs.\n");
        printf("Need at least %d sample points.\n",50*nInputs_);
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      TSIAnalyzer *tsiAnalyzer = new TSIAnalyzer();
      aData aPtr;
      aPtr.printLevel_ = outputLevel_;
      aPtr.nSamples_ = nSamples_;
      aPtr.nInputs_ = nInputs_;
      aPtr.nOutputs_ = nOutputs_;
      aPtr.sampleInputs_ = VecSamInputs_.getDVector();
      aPtr.sampleOutputs_ = VecSamOutputs_.getDVector();
      aPtr.sampleStates_ = VecSamStates_.getIVector();
      aPtr.iLowerB_ = VecILowerBs_.getDVector();
      aPtr.iUpperB_ = VecIUpperBs_.getDVector();
      aPtr.outputID_ = outputID;
      aPtr.ioPtr_ = psuadeIO_;
      tsiAnalyzer->analyze(aPtr);
      delete tsiAnalyzer;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sobol 
    //**/ sensitivity analysis using Sobol sampling and method
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sobol"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("tsi: Sobol sensitivity analysis (variance-based)\n");
        printf("     (suitable for raw data)\n");
        printf("Syntax: sobol (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs variance-based global ");
      printf("sensitivity analysis using\n");
      printf("the Sobol' method (first- and total-Sobol' indices). ");
      printf("The loaded sample\n");
      printf("must be a Sobol' sample.");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nSamples_ < 50*nInputs_) 
      {
        printf("This command is not recommended for small sample\n");
        printf("Need at least %d sample points.\n",50*nInputs_);
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      SobolAnalyzer *sobolAnalyzer = new SobolAnalyzer();
      aData aPtr;
      aPtr.printLevel_ = outputLevel_;
      aPtr.nSamples_ = nSamples_;
      aPtr.nInputs_ = nInputs_;
      aPtr.nOutputs_ = nOutputs_;
      aPtr.sampleInputs_ = VecSamInputs_.getDVector();
      aPtr.sampleOutputs_ = VecSamOutputs_.getDVector();
      aPtr.sampleStates_ = VecSamStates_.getIVector();
      aPtr.iLowerB_ = VecILowerBs_.getDVector();
      aPtr.iUpperB_ = VecIUpperBs_.getDVector();
      aPtr.outputID_ = outputID;
      aPtr.ioPtr_ = psuadeIO_;
      sobolAnalyzer->analyze(aPtr);
      delete sobolAnalyzer;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ fast 
    //**/ Fourier Amplitude Sensitivity Test
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "fast"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("fast: Fourier Amplitude Sensitivity Test to\n");
        printf("      compute first order sensitivity indices.\n");
        printf("Syntax: fast (after data have been loaded)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs variance-based global ");
      printf("sensitivity analysis using\n");
      printf("the Fourier Amplitude Sensitivity Test (or FAST) ");
      printf("method. The loaded\n");
      printf("sample must be a FAST sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_FAST;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ splot 
    //**/ generate matlab scatter plot file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "splot"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("splot: create scatter plots (against each parameter).\n");
        printf("Syntax: splot (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      printf("This command creates %d scatter plots (selected output ",
             nInputs_);
      printf("vs each input.)\n"); 
      sprintf(pString, "Select which output (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      //**/ write to scilab file
      if (plotScilab())
      {
        fp = fopen("scilabsp.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabsp.sci.\n");
          cmdStatus = 1;
          continue;
        }
      }
      else
      {
        //**/ write to matlab file
        fp = fopen("matlabsp.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabsp.m.\n");
          cmdStatus = 1;
          continue;
        }
      }
      sprintf(pString," plotMode=0  : plot all in a single plot");
      fwriteComment(fp, pString);
      sprintf(pString," plotMode!=0 : plot one at a time");
      fwriteComment(fp, pString);
      fprintf(fp, "plotMode=0;\n");
      fprintf(fp, "Y = [\n");
      for (sInd = 0; sInd < nSamples_; sInd++)
        fprintf(fp, "%24.16e\n",VecSamOutputs_[sInd*nOutputs_+outputID]);
      fprintf(fp, "];\n");
      for (iInd = 0; iInd < nInputs_; iInd++)
      {
        fprintf(fp, "X%d = [\n", iInd+1);
        for (sInd = 0; sInd < nSamples_; sInd++)
          fprintf(fp, "%24.16e\n",VecSamInputs_[sInd*nInputs_+iInd]);
        fprintf(fp, "];\n");
      }
      fprintf(fp, "S = [\n");
      for (sInd = 0; sInd < nSamples_; sInd++)
        fprintf(fp, "%d\n",VecSamStates_[sInd]);
      fprintf(fp, "];\n");
      if (plotScilab())
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
        {
          fwritePlotCLF(fp);
          fprintf(fp, "drawlater\n");
          fprintf(fp, "plot(X%d,Y,'.','markersize',10)\n",iInd+1);
          fprintf(fp, "a = gca();\n");
          fprintf(fp, "a.children.children.mark_foreground = 5;\n");
          sprintf(winput, "%s vs %s", StrOutNames_[outputID],
                  StrInpNames_[iInd]);
          fwritePlotTitle(fp, winput);
          fwritePlotAxes(fp);
          fwritePlotXLabel(fp, StrInpNames_[iInd]);
          fwritePlotYLabel(fp, StrOutNames_[outputID]);
          fprintf(fp, "drawnow\n");
          if (iInd < nInputs_-1) 
          {
            fprintf(fp, "disp(\'Press enter to advance,\')\n");
            fprintf(fp, "halt;\n");
          }
        }
        printf("scilabsp.sci is now available for scatter plots.\n");
      }
      else
      {
        fprintf(fp, "fs=6;\n");
        fwritePlotCLF(fp);
        kk = ((int) pow(1.0*nInputs_-0.1, 0.5)) + 1;
        ll = kk;
        while ((ll - 1) * kk >= nInputs_) ll--; 
        for (iInd = 0; iInd < nInputs_; iInd++)
        {
          fprintf(fp,"if plotMode == 0\n");
          fprintf(fp,"subplot(%d,%d,%d)\n",kk,ll,iInd+1);
          fprintf(fp,"else\n");
          if (iInd > 0)
          {
            fprintf(fp,"pause\n");
            fprintf(fp,"disp('Press enter to continue')\n");
          }
          fwritePlotCLF(fp);
          fprintf(fp,"end;\n");
          fprintf(fp,"iset = find(S == 0);\n");
          fprintf(fp,"plot(X%d(iset),Y(iset),'rX','markersize',fs)\n",
                  iInd+1);
          fprintf(fp,"hold on\n");
          fprintf(fp,"iset = find(S == 1);\n");
          fprintf(fp,"plot(X%d(iset),Y(iset),'b*','markersize',fs)\n",
                  iInd+1);
          fprintf(fp,"hold off\n");
          fprintf(fp,"axis([%24.16e %24.16e min(Y) max(Y)])\n",
                  VecILowerBs_[iInd], VecIUpperBs_[iInd]);
          sprintf(winput, "%s vs %s", StrOutNames_[outputID],
                  StrInpNames_[iInd]);
          fwritePlotTitle(fp, winput);
          fwritePlotAxes(fp);
          fwritePlotXLabel(fp, StrInpNames_[iInd]);
          fwritePlotYLabel(fp, StrOutNames_[outputID]);
        }
        printf("matlabsp.m is now available for scatter plots.\n");
      }
      fclose(fp);    
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rawi2
    //**/ generate intersection surfaces for multiple outputs for 
    //**/ display with matlab (using raw instead of RS data)
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rawi2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rawi2: similar to rsi2 except with sample (not RS) data\n");
        printf("Syntax: rawi2 (no argument needed).\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 2)
      {
        printf("ERROR: rawi2 requires 2 inputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ < 2)
      {
        printf("ERROR: rawi2 requires 2 or more outputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: rawi2 is currently not available for scilab.\n");
        cmdStatus = 1;
        continue;
      }
      dtemp = pow(1.0*nSamples_, 0.5) + 1.0e-12;
      nPtsPerDim = (int) dtemp;
      if (nPtsPerDim * nPtsPerDim != nSamples_)
      {
        printf("rawi2 error: nSamples must be a square.\n");
        cmdStatus = 1;
      }
      else
      { 
        //**/ ask users to specify the output set
        int rsiNOutputs = 2;
        sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs_);
        rsiNOutputs = getInt(2, nOutputs_, pString);

        //**/ get the collection of output set
        psIVector vecRsiSet;
        vecRsiSet.setLength(rsiNOutputs);
        if (rsiNOutputs == nOutputs_)
        {
          for (ii = 0; ii < rsiNOutputs; ii++) vecRsiSet[ii] = ii;
        }
        else
        {
          for (ii = 0; ii < rsiNOutputs; ii++)
          {
            sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                    ii+1, nOutputs_);
            vecRsiSet[ii] = getInt(1, nOutputs_, pString);
            vecRsiSet[ii]--;
          }
        }
        psIMatrix matRSI;
        matRSI.setFormat(PS_MAT2D); 
        matRSI.setDim(nPtsPerDim, nPtsPerDim);
        int **rsiMatrix = matRSI.getIMatrix2D();
        for (ii = 0; ii < nPtsPerDim; ii++)
          for (jj = 0; jj < nPtsPerDim; jj++)
            rsiMatrix[ii][jj] = rsiNOutputs;

        fp = fopen("matlabrawi2.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabrawi2.m.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        psVector vecFaYOut;
        vecFaYOut.setLength(nSamples_);

        //**/ filter
        for (ii = 0; ii < rsiNOutputs; ii++)
        {
          jplot = vecRsiSet[ii];
          for (sInd = 0; sInd < nSamples_; sInd++)
            vecFaYOut[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];

          Ymin = vecFaYOut[0];
          for (sInd = 1; sInd < nSamples_; sInd++)
            if (vecFaYOut[sInd] < Ymin) Ymin = vecFaYOut[sInd];
          Ymax = vecFaYOut[0];
          for (sInd = 1; sInd < nSamples_; sInd++)
            if (vecFaYOut[sInd] > Ymax) Ymax = vecFaYOut[sInd];

          printf("Ymin and Ymax = %e %e\n", Ymin, Ymax);
          sprintf(pString,
             "Enter the lower threshold for output %d (min = %16.8e) : ",
             jplot, Ymin);
          threshL = getDouble(pString);
          sprintf(pString,
             "Enter the upper threshold for output %d (max = %16.8e) : ",
             jplot, Ymax);
          threshU = getDouble(pString);

          for (sInd = 0; sInd < nSamples_; sInd++)
          {
            ind  = sInd % nPtsPerDim;
            ind2 = sInd / nPtsPerDim;
            if (vecFaYOut[sInd] < threshL) rsiMatrix[ind][ind2]--;
            if (vecFaYOut[sInd] > threshU) rsiMatrix[ind][ind2]--;
          }
        }

        //**/ write data to a matlab file
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd+=nPtsPerDim)
          fprintf(fp, "%e\n", VecSamInputs_[sInd*2]);
        fprintf(fp, "];\n");
        fprintf(fp, "Y = [\n");
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
          fprintf(fp, "%e\n", VecSamInputs_[sInd*2+1]);
        fprintf(fp, "];\n");
        fprintf(fp, "A = [\n");
        count = 0;
        for (ii = 0;  ii < nPtsPerDim; ii++)
          for (jj = 0;  jj < nPtsPerDim; jj++)
            if (rsiMatrix[jj][ii] == 0) count++;
        if (count == nPtsPerDim*nPtsPerDim)
        {
          for (ii = 0;  ii < nPtsPerDim; ii++)
            for (jj = 0;  jj < nPtsPerDim; jj++) fprintf(fp, "0\n");
        }
        else
        {
          for (ii = 0;  ii < nPtsPerDim; ii++)
          {
            for (jj = 0;  jj < nPtsPerDim; jj++)
              if (rsiMatrix[jj][ii] == 0) fprintf(fp, "NaN\n");
                else fprintf(fp, "%d\n", rsiMatrix[jj][ii]);
          }
        }
        fprintf(fp, "];\n");
        fprintf(fp, "A = reshape(A,%d,%d);\n",nPtsPerDim, nPtsPerDim);
        fprintf(fp, "A(%d,%d) = %e;\n", nPtsPerDim, nPtsPerDim, 
                (double) rsiNOutputs);
        fprintf(fp, "contourf(X,Y,A)\n");
        fprintf(fp, "axis([%e %e %e %e])\n",VecILowerBs_[0],
                VecIUpperBs_[0],VecILowerBs_[1],VecIUpperBs_[1]);
        fwritePlotAxes(fp);
        fwritePlotXLabel(fp, StrInpNames_[0]);
        fwritePlotYLabel(fp, StrInpNames_[1]);
        fwritePlotTitle(fp, "Intersection Contour");
        fprintf(fp, "colorbar\n");
        fprintf(fp, "colormap(cool)\n");
        fclose(fp);
        printf("matlabrawi2.m is now available for plotting.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ rawi3 
    //**/ generate 3D response surface and write the grid data to file
    //**/ for display with matlab using raw data
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rawi3"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rawi3: similar to rsi3 except with sample (not RS) data\n");
        printf("Syntax: rawi3 (no argument needed).\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 3)
      {
        printf("ERROR: rawi3 requires 3 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ < 2)
      {
        printf("ERROR: rawi3 requires 2 or more outputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: rawi3 is currently not available for scilab.\n");
        continue;
      }
      //**/ set up the function approximator
      dtemp = pow(1.0*nSamples_, 0.333333) + 0.1;
      if (nPtsPerDim*nPtsPerDim*nPtsPerDim != nSamples_)
      {
        printf("rawi3 error: nSamples must be an integer 3rd power.\n");
      }
      iplot1 = 0; iplot2 = 1; iplot3 = 2;

      //**/ ask users to specify the output set
      int rsiNOutputs = 2;
      sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs_);
      rsiNOutputs = getInt(2, nOutputs_, pString);

      //**/ get the collection of output set
      psIVector vecRsiSet;
      vecRsiSet.setLength(rsiNOutputs);
      if (rsiNOutputs == nOutputs_)
      {
        for (ii = 0; ii < rsiNOutputs; ii++) vecRsiSet[ii] = ii;
      }
      else
      {
        for (ii = 0; ii < rsiNOutputs; ii++)
        {
          sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                  ii+1, nOutputs_);
          vecRsiSet[ii] = getInt(1, nOutputs_, pString);
          vecRsiSet[ii]--;
        }
      }

      //**/ allocate space for the intersection set
      double ***rsi3Matrix = new double**[nPtsPerDim];
      for (ii = 0; ii < nPtsPerDim; ii++)
      {
        rsi3Matrix[ii] = new double*[nPtsPerDim];
        for (jj = 0; jj < nPtsPerDim; jj++)
        {
          rsi3Matrix[ii][jj] = new double[nPtsPerDim];
          for (kk = 0; kk < nPtsPerDim; kk++)
            rsi3Matrix[ii][jj][kk] = rsiNOutputs;
        }
      }

      //**/ generate the intersection set
      psVector vecFaYOut;
      vecFaYOut.setLength(nSamples_);
      for (ii = 0; ii < rsiNOutputs; ii++)
      {
        jplot = vecRsiSet[ii];
        for (sInd = 0; sInd < nSamples_; sInd++)
          vecFaYOut[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];

        Ymax = vecFaYOut[0];
        Ymin = vecFaYOut[0];
        for (sInd = 1; sInd < nSamples_; sInd++)
          if (vecFaYOut[sInd] < Ymin) Ymin = vecFaYOut[sInd];
        for (sInd = 1; sInd < nSamples_; sInd++)
          if (vecFaYOut[sInd] > Ymax) Ymax = vecFaYOut[sInd];

        printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot,
               Ymin, Ymax);
        sprintf(pString,"Enter the lower threshold (min = %e) : ", Ymin);
        threshL = getDouble(pString);
        sprintf(pString,"Enter the upper threshold (max = %e) : ", Ymax);
        threshU = getDouble(pString);

        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          ind  = (sInd % (nPtsPerDim * nPtsPerDim)) % nPtsPerDim;
          ind2 = (sInd % (nPtsPerDim * nPtsPerDim)) / nPtsPerDim;
          kk   = sInd / (nPtsPerDim * nPtsPerDim);
          if (vecFaYOut[sInd] < threshL) rsi3Matrix[ind][ind2][kk]--;
          if (vecFaYOut[sInd] > threshU) rsi3Matrix[ind][ind2][kk]--;
        }
      }
      for (ii = 0; ii < nPtsPerDim; ii++)
      {
        for (jj = 0; jj < nPtsPerDim; jj++)
        {
          if (rsi3Matrix[ii][jj][kk] == 0.0)
            rsi3Matrix[ii][jj][kk] = 0.5;
        }
      }

      //**/ generate and write response surface data
      fp = fopen("matlabrawi3.m", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file matlabrawi3.m.\n");
        continue;
      }
      fwritePlotCLF(fp);
      fprintf(fp, "xlo = %e; \n", VecILowerBs_[1]);
      fprintf(fp, "xhi = %e; \n", VecIUpperBs_[1]);
      fprintf(fp, "ylo = %e; \n", VecILowerBs_[0]);
      fprintf(fp, "yhi = %e; \n", VecIUpperBs_[0]);
      fprintf(fp, "zlo = %e; \n", VecILowerBs_[2]);
      fprintf(fp, "zhi = %e; \n", VecIUpperBs_[2]);
      fprintf(fp, "X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
      fprintf(fp, "Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
      fprintf(fp, "Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
      fprintf(fp, "V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
      for (jj = 0; jj < nPtsPerDim; jj++)
      {
        fprintf(fp, "Y(:,:,%d) = [\n", jj + 1);
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
        {
          for (ii = 0; ii < nPtsPerDim; ii++)
          {
            ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
            fprintf(fp, "%e ", VecSamInputs_[ind*3]);
          }
          fprintf(fp, "\n");
        }
        fprintf(fp, "];\n");
        fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
        {
          for (ii = 0; ii < nPtsPerDim; ii++)
          {
            ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
            fprintf(fp, "%e ", VecSamInputs_[ind*3+1]);
          }
          fprintf(fp, "\n");
        }
        fprintf(fp, "];\n");
        fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
        {
          for (ii = 0; ii < nPtsPerDim; ii++)
          {
            ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
            fprintf(fp, "%e ", VecSamInputs_[ind*3+2]);
          }
          fprintf(fp, "\n");
        }
        fprintf(fp, "];\n");
      }
      count = 0;
      for (jj = 0; jj < nPtsPerDim; jj++)
      {
        fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
        {
          for (ii = 0; ii < nPtsPerDim; ii++)
          {
            ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
            fprintf(fp, "%e ", rsi3Matrix[jj][ii][sInd]);
            if (rsi3Matrix[jj][ii][sInd] == 0.5) count++;
          }
          fprintf(fp, "\n");
        }
        fprintf(fp, "];\n");
      }
      if (count == nPtsPerDim*nPtsPerDim*nPtsPerDim)
      {
        fprintf(fp, "V(1,1,1)=0;\n");
        fprintf(fp, "V(%d,%d,%d)=1;\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
      }

      fprintf(fp, "xt = [%e:%e:%e];\n", VecILowerBs_[1],
              (VecIUpperBs_[1]-VecILowerBs_[1])*0.01, VecIUpperBs_[1]);
      fprintf(fp, "yt = [%e:%e:%e];\n", VecILowerBs_[iplot1],
              (VecIUpperBs_[0]-VecILowerBs_[0])*0.01, VecIUpperBs_[0]);
      fprintf(fp, "zt = [%e:%e:%e];\n", VecILowerBs_[2],
              (VecIUpperBs_[2]-VecILowerBs_[2])*0.01, VecIUpperBs_[2]);
      fwritePlotCLF(fp);
      fprintf(fp, "isoval = 0.5;\n");
      fprintf(fp, "h = patch(isosurface(X,Y,Z,V,isoval),... \n");
      fprintf(fp, "          'FaceColor', 'blue', ... \n");
      fprintf(fp, "          'EdgeColor', 'none', ... \n");
      fprintf(fp, "          'AmbientStrength', 0.2, ... \n");
      fprintf(fp, "          'SpecularStrength', 0.7, ... \n");
      fprintf(fp, "          'DiffuseStrength', 0.4);\n");
      fprintf(fp, "isonormals(X,Y,Z,V,h);\n");
      fprintf(fp, "patch(isocaps(X,Y,Z,V,isoval), ...\n");
      fprintf(fp, "      'FaceColor', 'interp', ... \n");
      fprintf(fp, "      'EdgeColor', 'none'); \n");
      fprintf(fp, "axis([xlo xhi ylo yhi zlo zhi])\n");
      fprintf(fp, "daspect([%e,%e,%e])\n",VecIUpperBs_[1]-VecILowerBs_[1],
              VecIUpperBs_[0]-VecILowerBs_[0], 
              VecIUpperBs_[2]-VecILowerBs_[2]);
      fprintf(fp, "colormap('default'); colorbar\n");
      fprintf(fp, "%%axis tight\n");
      fprintf(fp, "view(3) \n");
      fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
      fprintf(fp, "box on\n");
      fprintf(fp, "grid on\n");
      fprintf(fp, "lighting phong\n");
      fwritePlotAxes(fp);
      fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
              StrInpNames_[1]);
      fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
              StrInpNames_[0]);
      fprintf(fp, "zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
              StrInpNames_[2]);
      fprintf(fp, "title('Intersection Contour','FontWeight',");
      fprintf(fp, "'bold','FontSize',12)\n");
      fprintf(fp, "colorbar\n");
      fprintf(fp, "colormap(cool)\n");
      fclose(fp);
      printf("matlabrawi3.m is now available for plotting.\n");
      for (ii = 0; ii < nPtsPerDim; ii++) 
      {
        for (jj = 0; jj < nPtsPerDim; jj++) 
          delete [] rsi3Matrix[ii][jj];
        delete [] rsi3Matrix[ii];
      }
      delete [] rsi3Matrix;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rspairs
    //**/ generate response surface of all 2-input pairs and write the
    //**/ grid data to file for display with matlab
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rspairs"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rspairs: generate RS of all 2-input pairs\n");
        printf("Syntax: rspairs (no argument needed).\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 2)
      {
        printf("ERROR: rspairs requires 2 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: rspairs is currently not available for scilab.\n");
        continue;
      }
      //**/ set up the function approximator
      nPtsPerDim = 64;
      sprintf(pString, "Grid resolution ? (32 - 128) ");
      nPtsPerDim = getInt(32, 128, pString);
      faFlag = 1;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);

      //**/ get which inputs to plot
      nPlots = 0;
      sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
              nInputs_);
      nPlots = getInt(1, nInputs_, pString);
      psIVector vecPlotInds;
      vecPlotInds.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) vecPlotInds[ii] = ii;
      if (nPlots < nInputs_)
      {
        for (ii = 0; ii < nPlots; ii++)
        {
          sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                  nInputs_);
          vecPlotInds[ii] = getInt(1, nInputs_, pString);
          vecPlotInds[ii]--;
        }
      }
      psVector vecInpVals;
      vecInpVals.setLength(nInputs_);
      if (nInputs_ > 2)
      {
        sprintf(pString,
                "Set other nominal values at mid point ? (y or n) ");
        getString(pString, command);
      }
      if (command[0] == 'y')
      {
        for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
          vecInpVals[iInd1] = 
              0.5*(VecILowerBs_[iInd1]+VecIUpperBs_[iInd1]);
      }
      else
      {
        printf("Enter data file for nominal values. Format: \n");
        printf("PSUADE_BEGIN (optional)\n");
        printf("<numInputs>\n");
        printf("1   <data>\n");
        printf("2   <data>\n");
        printf("..  <data>\n");
        printf("PSUADE_END (optional)\n");
        printf("Data file name : ");
        scanf("%s", dataFile);
        fgets(lineIn,5000,stdin); 
        fp = fopen(dataFile, "r");
        if (fp == NULL)
        {
          printf("ERROR: data file %s not found.\n", dataFile);
          continue;
        }
        else
        {
          fscanf(fp, "%s", winput);
          if (strcmp(winput, "PSUADE_BEGIN"))
          {
            printf("INFO: the optional PSUADE_BEGIN keyword is not found.\n");
            printf("      First line should have nInputs\n");
            fclose(fp);
            fp = fopen(dataFile, "r");
          }
          else
          {
            fscanf(fp, "%d", &kk);
            if (kk != nInputs_)
            {
              printf("ERROR: input size does not match nInputs.\n");
              fclose(fp);
              cmdStatus = 1;
              continue;
            }
            for (ii = 0; ii < nInputs_; ii++)
            {
              fscanf(fp, "%d", &ind);
              if (ind != (ii+1))
              {
                printf("ERROR: input index mismatch (%d,%d)\n",
                       ii+1,ind);
                break;
              }
              fscanf(fp, "%lg", &ddata);
              vecInpVals[ii] = ddata;
            }
            if (ii != nInputs_)
            {
              fclose(fp);
              continue;
            }
            fscanf(fp, "%s", winput);
            fscanf(fp, "%s", winput);
            fclose(fp);
            if (strcmp(winput, "PSUADE_END"))
              printf("INFO: the optional PSUADE_END keyword is not found.\n");
          }
        }
      }

      //**/ open matlab file
      fp = fopen("matlabrspairs.m", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file matlabrspairs.m.\n");
        cmdStatus = 1;
        continue;
      }
      fwritePlotCLF(fp);

      jplot = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      jplot = getInt(1, nOutputs_, pString);
      jplot--;
      Ymax = Ymin = VecSamOutputs_[0*nOutputs_+jplot];
      for (sInd = 1; sInd < nSamples_; sInd++)
      {
       if (VecSamOutputs_[sInd*nOutputs_+jplot] > Ymax) 
         Ymax = VecSamOutputs_[sInd*nOutputs_+jplot];
       if (VecSamOutputs_[sInd*nOutputs_+jplot] < Ymin) 
         Ymin = VecSamOutputs_[sInd*nOutputs_+jplot];
      }
      printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
      sprintf(pString,"Set lower threshold ? (y or n) ");
      getString(pString, command);
      if (command[0] == 'y')
      {
        sprintf(pString,"Enter the lower threshold (min = %e) : ", 
                Ymin);
        thresh = getDouble(pString);
        fprintf(fp, "Ymin = %e;\n", thresh);
      }
      else fprintf(fp, "Ymin = -1.0e35;\n");
      sprintf(pString,"Set upper threshold ? (y or n) : ");
      getString(pString, command);
      if (command[0] == 'y')
      {
        sprintf(pString,"Enter the upper threshold (max = %e) : ", 
                Ymax);
        thresh = getDouble(pString);
        fprintf(fp, "Ymax = %e;\n", thresh);
      }
      else fprintf(fp, "Ymax = 1.0e35;\n");

      psVector vecFaYIn;
      vecFaYIn.setLength(nSamples_);
      double diagMax = - PSUADE_UNDEFINED;
      double diagMin =   PSUADE_UNDEFINED;
      double *faXOut = NULL;
      double *faYOut = NULL;
      for (ii = 0; ii < nPlots; ii++)
      {
        iplot2 = vecPlotInds[ii];
        for (jj = 0; jj <= ii; jj++)
        {
          iplot1 = vecPlotInds[jj];
          if (iplot1 == iplot2)
          {
            for (sInd = 0; sInd < nSamples_; sInd++)
              vecFaYIn[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];
            faPtr->gen1DGridData(VecSamInputs_.getDVector(),
                      vecFaYIn.getDVector(),iplot2,
                      vecInpVals.getDVector(), &faLeng,&faXOut,&faYOut);
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
            {
              fprintf(fp, "%e\n", faYOut[sInd]);
              if (faYOut[sInd] > diagMax) diagMax = faYOut[sInd];
              if (faYOut[sInd] < diagMin) diagMin = faYOut[sInd];
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd++)
              fprintf(fp, "%e\n", faXOut[sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "B = A;\n");
            fprintf(fp, "[ia,aa] = find(A<Ymin);\n");
            fprintf(fp, "for ii = 1 : length(ia)\n");
            fprintf(fp, "   B(ia(ii)) = NaN;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "n1 = length(ia);\n");
            fprintf(fp, "[ia,aa] = find(A>Ymax);\n");
            fprintf(fp, "for ii = 1 : length(ia)\n");
            fprintf(fp, "   B(ia(ii)) = NaN;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "n2 = length(ia);\n");
            fprintf(fp, "if (n1 + n2 == %d)\n",faLeng);
            fprintf(fp, "   B(1) = 0;\n");
            fprintf(fp, "   B(%d) = 1;\n",faLeng);
            fprintf(fp, "end;\n");
            fprintf(fp, "subplot(%d,%d,%d), ",
                    nPlots,nPlots,ii*nPlots+jj+1);
            fprintf(fp, "plot(X,B,'LineWidth',2.0)\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, StrInpNames_[iplot1]);
            fwritePlotYLabel(fp, StrOutNames_[jplot]);
            delete [] faXOut;
            delete [] faYOut;
            continue;
          }

          //**/ generate 2D data
          for (sInd = 0; sInd < nSamples_; sInd++)
            vecFaYIn[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];
          faPtr->gen2DGridData(VecSamInputs_.getDVector(),
                     vecFaYIn.getDVector(),iplot1,iplot2,
                     vecInpVals.getDVector(),&faLeng,&faXOut,&faYOut);

          //**/ write to matlab file
          fprintf(fp, "A = [\n");
          for (sInd = 0; sInd < faLeng; sInd++)
            fprintf(fp, "%e\n", faYOut[sInd]);
          fprintf(fp, "];\n");
          fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
          fprintf(fp, "X = [\n");
          for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
            fprintf(fp, "%e\n", faXOut[sInd*2]);
          fprintf(fp, "];\n");
          fprintf(fp, "Y = [\n");
          for (sInd = 0; sInd < nPtsPerDim; sInd++)
            fprintf(fp, "%e\n", faXOut[sInd*2+1]);
          fprintf(fp, "];\n");
          fprintf(fp, "B = A;\n");
          fprintf(fp, "[ia,ja,aa] = find(A<Ymin);\n");
          fprintf(fp, "for ii = 1 : length(ia)\n");
          fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
          fprintf(fp, "end;\n");
          fprintf(fp, "n1 = length(ia);\n");
          fprintf(fp, "[ia,ja,aa] = find(A>Ymax);\n");
          fprintf(fp, "for ii = 1 : length(ia)\n");
          fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
          fprintf(fp, "end;\n");
          fprintf(fp, "n2 = length(ia);\n");
          fprintf(fp, "nA = size(A,1);\n");
          fprintf(fp, "if (n1 + n2 == nA*nA)\n");
          fprintf(fp, "   B(1,1) = 0;\n");
          fprintf(fp, "   B(%d,%d) = 1;\n",nPtsPerDim,nPtsPerDim);
          fprintf(fp, "end;\n");
          fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,ii*nPlots+jj+1);
          fprintf(fp, "contourf(X,Y,B)\n");
          fwritePlotAxes(fp);
          fwritePlotXLabel(fp, StrInpNames_[iplot1]);
          fwritePlotYLabel(fp, StrInpNames_[iplot2]);
          fprintf(fp, "axis([%e %e %e %e])\n",VecILowerBs_[iplot1],
                  VecIUpperBs_[iplot1],VecILowerBs_[iplot2],
                  VecIUpperBs_[iplot2]);
          fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,jj*nPlots+ii+1);
          fprintf(fp, "surf(X,Y,B)\n");
          fwritePlotAxes(fp);
          fwritePlotXLabel(fp, StrInpNames_[iplot1]);
          fwritePlotYLabel(fp, StrInpNames_[iplot2]);
          delete [] faXOut;
          delete [] faYOut;
        }
      }
      if (diagMax - diagMin == 0) diagMax += 0.1;
      for (ii = 0; ii < nPlots; ii++)
      {
        iplot1 = vecPlotInds[ii];
        fprintf(fp, "subplot(%d,%d,%d), ",nPlots,nPlots,ii*nPlots+ii+1);
        fprintf(fp, "axis([%e %e %e %e])\n",VecILowerBs_[iplot1],
                VecIUpperBs_[iplot1],diagMin, diagMax);
      }
      fclose(fp);
      printf("matlabrspairs.m is now available for contour plots.\n");
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsipairs
    //**/ generate intersection surfaces for multiple outputs for 
    //**/ display with matlab
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsipairs"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsipairs: generate RS intersections for all input pairs\n");
        printf("Syntax: rsipairs (no argument needed).\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 2)
      {
        printf("ERROR: rsipairs requires 2 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: rsipairs is currently not available for scilab.\n");
        continue;
      }
      //**/ set up the function approximator
      nPtsPerDim = 64;
      sprintf(pString, "Grid resolution ? (32 - 128) ");
      nPtsPerDim = getInt(32, 128, pString);
      faFlag = 1;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);

      //**/ get which inputs to plot
      nPlots = 0;
      sprintf(pString, "Enter the number of inputs to plot (1 - %d) : ",
              nInputs_);
      nPlots = getInt(1, nInputs_, pString);
      psIVector vecPlotInds;
      vecPlotInds.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) vecPlotInds[ii] = ii;
      if (nPlots < nInputs_)
      {
        for (ii = 0; ii < nPlots; ii++)
        {
          sprintf(pString, "Enter the %d-th input (1 - %d) : ", ii+1,
                  nInputs_);
          vecPlotInds[ii] = getInt(1, nInputs_, pString);
          vecPlotInds[ii]--;
        }
      }
  
      //**/ ask users to specify the output set
      printf("rsipairs: will use all outputs for constraining.\n");
      int rsiNOutputs = nOutputs_;
      psIVector vecRsiSet;
      vecRsiSet.setLength(rsiNOutputs);
      for (ii = 0; ii < rsiNOutputs; ii++) vecRsiSet[ii] = ii;
      psVector vecThreshLs, vecThreshUs;
      vecThreshLs.setLength(rsiNOutputs);
      vecThreshUs.setLength(rsiNOutputs);
      for (ii = 0; ii < rsiNOutputs; ii++)
      {
        jplot = vecRsiSet[ii];
        Ymax = Ymin = VecSamOutputs_[0*nOutputs_+jplot];
        for (sInd = 1; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_+jplot] > Ymax) 
            Ymax = VecSamOutputs_[sInd*nOutputs_+jplot];
          if (VecSamOutputs_[sInd*nOutputs_+jplot] < Ymin) 
            Ymin = VecSamOutputs_[sInd*nOutputs_+jplot];
        }
        printf("Output %d: Ymin and Ymax = %e %e\n",jplot+1,Ymin,Ymax);
        sprintf(pString,
                "Enter the lower threshold (min = %16.8e) : ", Ymin);
        vecThreshLs[ii] = getDouble(pString);
        sprintf(pString,
                "Enter the upper threshold (max = %16.8e) : ", Ymax);
        vecThreshUs[ii] = getDouble(pString);
      }
      for (ii = 0; ii < rsiNOutputs; ii++)
      {
        printf("Lower and upper thresholds for output %d = %e %e\n",
               vecRsiSet[ii]+1, vecThreshLs[ii], vecThreshUs[ii]);
      }

      //**/ set up rsi data store
      psVector vecFaYIn;
      vecFaYIn.setLength(nSamples_);
      psIMatrix matRSI;
      matRSI.setFormat(PS_MAT2D);
      matRSI.setDim(nPtsPerDim, nPtsPerDim);
      int **rsiMatrix = matRSI.getIMatrix2D();
      fp = fopen("matlabrsipairs.m", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file matlabrsipairs.m.\n");
        cmdStatus = 1;
        continue;
      }
      fprintf(fp, "hold off\n");
      fwritePlotCLF(fp);

      //**/ interpolate
      double *faXOut = NULL;
      double *faYOut = NULL;
      psVector vecInpVals;
      vecInpVals.setLength(nInputs_);
      for (ii = 0; ii < nPlots; ii++)
      {
        iplot2 = vecPlotInds[ii];
        for (jj = 0; jj <= ii; jj++)
        {
          iplot1 = vecPlotInds[jj];
          for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
          {
            if (iInd1 != iplot1 && iInd1 != iplot2)
              vecInpVals[iInd1] = 0.5 *
                       (VecILowerBs_[iInd1]+VecIUpperBs_[iInd1]);
            else vecInpVals[iInd1]=1.0;
          }
          for (iInd1 = 0; iInd1 < nPtsPerDim; iInd1++)
            for (iInd2 = 0; iInd2 < nPtsPerDim; iInd2++)
              rsiMatrix[iInd1][iInd2] = rsiNOutputs;
          for (kk = 0; kk < rsiNOutputs; kk++)
          {
            jplot = vecRsiSet[kk];
            for (sInd = 0; sInd < nSamples_; sInd++)
              vecFaYIn[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];

            if (iplot1 == iplot2)
            {
              faPtr->gen1DGridData(VecSamInputs_.getDVector(),
                       vecFaYIn.getDVector(), iplot1,
                       vecInpVals.getDVector(),&faLeng,&faXOut,&faYOut);
              if (kk == 0)
              {
                fprintf(fp, "X = [\n");
                for (sInd = 0; sInd < faLeng; sInd++)
                  fprintf(fp, "%e\n", faXOut[sInd]);
                fprintf(fp, "];\n");
              }
              for (sInd = 0; sInd < faLeng; sInd++)
              {
                if (faYOut[sInd]<vecThreshLs[kk]) rsiMatrix[sInd][sInd]--;
                if (faYOut[sInd]>vecThreshUs[kk]) rsiMatrix[sInd][sInd]--;
              }
            }
            else
            {
              faPtr->gen2DGridData(VecSamInputs_.getDVector(),
                        vecFaYIn.getDVector(),iplot1,iplot2,
                        vecInpVals.getDVector(),&faLeng,&faXOut,&faYOut);
              if (kk == 0)
              {
                fprintf(fp, "X = [\n");
                for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
                  fprintf(fp, "%e\n", faXOut[sInd*2]);
                fprintf(fp, "];\n");
                fprintf(fp, "Y = [\n");
                for (sInd = 0; sInd < nPtsPerDim; sInd++)
                  fprintf(fp, "%e\n", faXOut[sInd*2+1]);
                fprintf(fp, "];\n");
              }
              for (sInd = 0; sInd < faLeng; sInd++)
              {
                ind  = sInd % nPtsPerDim;
                ind2 = sInd / nPtsPerDim;
                if (faYOut[sInd] < vecThreshLs[kk]) rsiMatrix[ind][ind2]--;
                if (faYOut[sInd] > vecThreshUs[kk]) rsiMatrix[ind][ind2]--;
              }
            }
            delete [] faXOut;
            delete [] faYOut;
          }
          if (iplot1 == iplot2)
          {
            fprintf(fp, "A = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
              fprintf(fp, "%d\n", rsiMatrix[sInd][sInd]);
            fprintf(fp, "];\n");
            fprintf(fp, "subplot(%d,%d,%d), ",
                    nPlots,nPlots,ii*nPlots+jj+1);
            fprintf(fp, "plot(X,A,'LineWidth',2.0)\n");
            fprintf(fp, "axis([%e %e 0 %d])\n",VecILowerBs_[iplot1],
                    VecIUpperBs_[iplot1],rsiNOutputs+1);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, StrInpNames_[iplot1]);
          }
          else
          {
            fprintf(fp, "A = [\n");
            for (iInd1 = 0; iInd1 < nPtsPerDim; iInd1++)
              for (iInd2 = 0; iInd2 < nPtsPerDim; iInd2++)
                fprintf(fp, "%d\n", rsiMatrix[iInd2][iInd1]);
            fprintf(fp, "];\n");
            fprintf(fp, "A = reshape(A,%d,%d);\n",nPtsPerDim,
                    nPtsPerDim);
            fprintf(fp, "A(1,1) = 0;\n");
            fprintf(fp, "A(%d,%d) = %d;\n",nPtsPerDim,nPtsPerDim,
                    rsiNOutputs); 
            fprintf(fp, "subplot(%d,%d,%d)\n",
                    nPlots,nPlots,ii*nPlots+jj+1);
            fprintf(fp, "contourf(X,Y,A)\n");
            fprintf(fp, "axis([%e %e %e %e])\n",VecILowerBs_[iplot1],
                    VecIUpperBs_[iplot1],VecILowerBs_[iplot2],
                    VecIUpperBs_[iplot2]);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, StrInpNames_[iplot1]);
            fwritePlotYLabel(fp, StrInpNames_[iplot2]);
            fprintf(fp, "subplot(%d,%d,%d)\n",
                    nPlots,nPlots,jj*nPlots+ii+1);
            fprintf(fp, "surf(X,Y,A)\n");
            fprintf(fp, "axis([%e %e %e %e 0 %d])\n",VecILowerBs_[iplot1],
                    VecIUpperBs_[iplot1],VecILowerBs_[iplot2],    
                    VecIUpperBs_[iplot2], rsiNOutputs);
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, StrInpNames_[iplot1]);
            fwritePlotYLabel(fp, StrInpNames_[iplot2]);
          }
        }
      }
      fclose(fp);
      printf("matlabrsipairs.m is now available for plotting.\n");
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rscheck, rsvalidate, rsfit 
    //**/ check quality of the response surface 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rscheck") || 
             !strcmp(command, "rstest_cv") ||
             !strcmp(command, "rsvalidate") ||
             !strcmp(command, "rsfit"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsfit: check RS quality (training errors and cross\n");
        printf("       validation errors)\n");
        printf("Syntax: rsfit or rsvalidate (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command can be used to check the quality of a ");
      printf("response surface\n");
      printf("trained on the loaded sample. The quality metrics may ");
      printf("be training\n");
      printf("errors (from resubstitution test) or cross validation ");
      printf("errors. Plots\n");
      printf("will be created for visualizing the goodness of the fit.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nSamples_ <= 1)
      {
        printf("ERROR: sample size (%d) too small.\n",nSamples_);
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      faType = -1;
      sprintf(pString, "Enter your choice ? ");
      while (faType < 0 || faType > PSUADE_NUM_RS)
      {
        writeFAInfo(outputLevel_);
        faType = getFAType(pString);
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      psIVector vecXsforms;
      vecXsforms.setLength(2);
      if (psConfig_.RSExpertModeIsOn() && psConfig_.InteractiveIsOn())
      {
        printf("Available input/output transformations :\n");
        printf("0. no transformation.\n");
        printf("1. log transformation on all the inputs.\n");
        printf("2. log transformation on all the selected output.\n");
        printf("3. log transformation on all inputs and outputs.\n");
        sprintf(pString, "Enter your choice ? ");
        int otrans = getInt(0, 3, pString);
        vecXsforms[0] = otrans & 1;
        vecXsforms[1] = otrans & 2;
      }
      else
      {
        vecXsforms[0] = 0;
        vecXsforms[1] = 0;
      }
      int analysisMethod = PSUADE_ANA_RSFA;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, faType);
      anaManager->loadLogXsformFlags(2, vecXsforms.getIVector());
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 1, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rstest_ts 
    //**/ response surface test on training set
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstest_ts"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rstest_ts: response surface test on the training set.\n");
        printf("Syntax: rstest_ts (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      faType = -1;
      sprintf(pString, "Enter your choice ? ");
      while (faType < 0 || faType > PSUADE_NUM_RS)
      {
        writeFAInfo(outputLevel_);
        faType = getFAType(pString);
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      status = faPtr->initialize(VecSamInputs_.getDVector(),
                                 vecYT.getDVector());
      vecVT.setLength(nSamples_);
      vecWT.setLength(nSamples_);
      faPtr->evaluatePointFuzzy(nSamples_,VecSamInputs_.getDVector(), 
                                vecVT.getDVector(), vecWT.getDVector());
      delete faPtr;
      faPtr = NULL;
      if (plotScilab()) fp = fopen("scilabrstestts.sci", "w");
      else              fp = fopen("matlabrstestts.m", "w");
      sprintf(pString," col 1: simulation data, col 2: rs data");
      fwriteComment(fp, pString);
      sprintf(pString," col 3: std dev");
      fwriteComment(fp, pString);
      fprintf(fp, "A = [\n");
      for (ss = 0; ss < nSamples_; ss++)
        fprintf(fp, "%e %e %e\n",vecYT[ss],vecVT[ss],vecWT[ss]);
      fprintf(fp, "];\n");
      fwritePlotCLF(fp);
      fwritePlotFigure(fp, 1);
      fprintf(fp, "subplot(1,2,1)\n");
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, A(:,1)-A(:,2), style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,2),10);\n");
        fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Error Plot (unscaled)");
      fwritePlotXLabel(fp, "Error");
      fwritePlotYLabel(fp, "Probabilities");
      fprintf(fp, "subplot(1,2,2)\n");
      fprintf(fp, "xmax = max(A(:,1));\n");
      fprintf(fp, "xmin = min(A(:,1));\n");
      fprintf(fp, "ymax = max(A(:,2));\n");
      fprintf(fp, "ymin = min(A(:,2));\n");
      fprintf(fp, "xmin = min(xmin, ymin);\n");
      fprintf(fp, "xmax = max(xmax, ymax);\n");
      fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
      fprintf(fp, "plot(A(:,1), A(:,2),'x', XX, XX)\n");
      if (plotScilab())
      {
        fprintf(fp, "a = gca();\n");
        fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
      }
      else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Interpolated versus actual data");
      fwritePlotXLabel(fp, "Actual data");
      fwritePlotYLabel(fp, "Interpolated data");
      fclose(fp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += (vecYT[ss] - vecVT[ss]);
      ddata = ddata / nSamples_;
      printf("Training set avg error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        ddata += (vecYT[ss] - vecVT[ss]);
        if (vecYT[ss] != 0.0)
        {
          ddata /= PABS(vecYT[ss]);
          dtemp = vecYT[ss];
        }
      }
      ddata = ddata / nSamples_;
      printf("Training set avg error (  scaled) = %e (base=%e)\n",
              ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += pow(vecYT[ss] - vecVT[ss], 2.0);
      ddata = sqrt(ddata / nSamples_);
      printf("Training set rms error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        ddata += pow(vecYT[ss] - vecVT[ss], 2.0);
        if (vecYT[ss] != 0.0)
        {
          ddata /= pow(vecYT[ss],2.0);
          dtemp = vecYT[ss];
        }
      }
      ddata = sqrt(ddata / nSamples_);
      printf("Training set rms error (  scaled) = %e (base=%e)\n",
             ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        if (PABS(vecVT[ss]-vecYT[ss]) > ddata)
          ddata = PABS(vecYT[ss] - vecVT[ss]);
      printf("Training set max error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (vecYT[ss] != 0 && PABS((vecVT[ss]-vecYT[ss])/vecYT[ss])>ddata)
        {
          dtemp = vecYT[ss];
          ddata = PABS((vecYT[ss]-vecVT[ss])/vecYT[ss]);
        }
      }
      printf("Training set max error (  scaled) = %e (base=%e)\n",
             ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += (vecYT[ss] - vecVT[ss]);
      ddata = ddata / nSamples_;
      printf("Training set error mean = %e\n", ddata);
      dtemp = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        dtemp += pow(vecYT[ss] - vecVT[ss] - ddata, 2.0);
      dtemp = sqrt(dtemp / (nSamples_ - 1.0));
      printf("Training set error std  = %e\n", dtemp);
      if (plotScilab())
           printf("rstest_ts: plot file in scilabrstestts.sci\n");
      else printf("rstest_ts: plot file in matlabrstestts.m\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rstest_cv 
    //**/ response surface test - random cross validation
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstest_rcv"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rstest_rcv: RS random cross validation test.\n");
        printf("Syntax: rstest_cv (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs a different cross validation ");
      printf("test from the rsfit\n");
      printf("(or rscheck) command in that, instead of partitioning ");
      printf("the sample into\n");
      printf("non-overlapping sets, it randomly selects holdout ");
      printf("sample points at\n");
      printf("each iteration and uses the remaining points to ");
      printf("construct response\n");
      printf("surfaces. As such, any point in the sample may be ");
      printf("held out multiple\n");
      printf("times, so the cross validation predictions are ");
      printf("calculated by averaging\n");
      printf("the holdout predictions and by taking the worst and ");
      printf("best cases. Also,\n");
      printf("this command accommodates an auxiliary sample (which, ");
      printf("for example,\n");
      printf("covers corners and faces) that will be used in every ");
      printf("training sample\n");
      printf("instance during cross validation for preventing ");
      printf("extrapolation. Thus,\n");  
      printf("test should be more rigorous than rsfit, but it ");
      printf("comes with additional\n");
      printf("computational cost.\n\n");
      printf("To run this command, users are to be prompted for: \n");
      printf("1. An optional auxiliary sample\n");
      printf("2. Choice of response surface\n");
      printf("3. Which sample output to perform cross validation\n");
      printf("4. Size of the holdout set M (e.g. 1/10 of sample size N)\n");
      printf("5. Number of iterations K (K * M >> N)\n");
      printEquals(PL_INFO, 0);
      printf("This analysis will use cross validation to ");
      printf("predict the sample outputs.\n");
      printf("Average, worst-case, and best-case prediction ");
      printf("errors will be analyzed\n");
      printf("and plotted.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      printf("Sometimes you may have augmented which covers ");
      printf("corners and faces to\n");
      printf("prevent extrapolation. This augmented sample ");
      printf("will be used in every\n");
      printf("training set (that is, not to be held out).\n"); 
      sprintf(pString, "Do you have an augmented sample ? (y or n) ");
      getString(pString, winput);
      PsuadeData *pIO=NULL;
      pData pX, pY;
      psVector vecXAux, vecYAux;
      int nSamAux=0;
      if (winput[0] == 'y')
      {
        sprintf(pString, "Enter name of your augment sample file : ");
        getString(pString, dataFile);
        kk = strlen(dataFile);
        dataFile[kk-1] = '\0';
        pIO = new PsuadeData();
        status = pIO->readPsuadeFile(dataFile);
        if (status != 0)
        {
          printf("ERROR: file %s either not found or in wrong format.\n",
                 dataFile);
          cmdStatus = 1;
          continue;
        }
        pIO->getParameter("input_ninputs", pPtr);
        int nInps = pPtr.intData_;
        if (nInps != nInputs_)
        {
          printf("ERROR: nInputs are different.\n");
          printf("       Augmented sample nInputs = %d\n",ind);
          printf("       Resident  sample nInputs = %d\n",nInputs_);
          cmdStatus = 1;
          delete pIO;
          continue;
        }
        pIO->getParameter("output_noutputs",pPtr);
        int nOuts = pPtr.intData_;
        if (nOuts != nOutputs_)
        {
          printf("ERROR: nOutputs are different.\n");
          printf("       Augmented sample nOutputs = %d\n",ind);
          printf("       Resident  sample nOutputs = %d\n",nOutputs_);
          cmdStatus = 1;
          delete pIO;
          continue;
        }
        pIO->getParameter("method_nsamples", pPtr);
        nSamAux = pPtr.intData_;
        pIO->getParameter("input_sample", pX);
        vecXAux.load(nInputs_ * nSamAux, pX.dbleArray_);
        pIO->getParameter("output_sample", pY);
        vecYAux.load(nOutputs_ * nSamAux, pY.dbleArray_);
        delete pIO;
      }

      faType = -1;
      sprintf(pString, "Enter your response surface choice ? ");
      while (faType < 0 || faType > PSUADE_NUM_RS)
      {
        writeFAInfo(outputLevel_);
        faType = getFAType(pString);
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString, "Size of holdout sets (2 - %d)? ", nSamples_/2);
      int hoSize = getInt(2, nSamples_/2, pString);
      kk = nSamples_ / hoSize * 2;
      sprintf(pString, "Number of iterations (>= %d)? ",kk);
      int numIter = getInt(kk, 1000000, pString);

      psVector vecXTrain, vecYTrain, vecXHO, vecYHO, vecYPred, vecYWorst;
      psVector vecYBest;
      vecXTrain.setLength((nSamAux + nSamples_) * nInputs_);
      vecYTrain.setLength(nSamAux+nSamples_);
      vecYPred.setLength(nSamples_);
      vecYWorst.setLength(nSamples_);
      vecYBest.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++) vecYBest[ss] = 1e35;
      for (ss = 0; ss < nSamples_; ss++) 
        vecYWorst[ss] = VecSamOutputs_[ss*nOutputs_+outputID];;
      vecXHO.setLength(hoSize * nInputs_);
      vecYHO.setLength(hoSize);
      psConfig_.InteractiveSaveAndReset();
      psConfig_.RSExpertModeSaveAndReset();
      psIVector vecSCnt, vecIRand;
      vecSCnt.setLength(nSamples_);
      vecIRand.setLength(nSamples_);
      int allTested = 0, iterCnt=0;
      FuncApprox *faPtr = genFA(faType,nInputs_,iOne,nSamAux+nSamples_-hoSize);
      faPtr->setBounds(VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      while (allTested == 0 || iterCnt < numIter)
      {
        printf("Processing iteration %d (out of %d)\n",iterCnt+1,numIter);
        if (iterCnt >= numIter)
        {
          printf("     Continuing: not all points have been held out.\n");
          count = 0;
          for (kk = 0; kk < nSamples_; kk++) 
            if (vecSCnt[kk] == 0) count++;
          printf("     Number of not-yet-heldout points = %d\n",count);
        }
          
        //**/ create a holdout set
        generateRandomIvector(nSamples_, vecIRand.getIVector());
        for (kk = 0; kk < hoSize; kk++) 
        {
          ind = vecIRand[kk];
          for (ii = 0; ii < nInputs_; ii++)
             vecXHO[kk*nInputs_+ii] = VecSamInputs_[ind*nInputs_+ii]; 
          vecYHO[kk] = VecSamOutputs_[ind*nOutputs_+outputID]; 
          vecSCnt[ind]++;
        }
        //**/ create a RS training set 
        for (kk = hoSize; kk < nSamples_; kk++) 
        {
          ind = vecIRand[kk];
          for (ii = 0; ii < nInputs_; ii++)
             vecXTrain[(kk-hoSize)*nInputs_+ii] = 
                     VecSamInputs_[ind*nInputs_+ii]; 
          vecYTrain[kk-hoSize] = VecSamOutputs_[ind*nOutputs_+outputID]; 
        }  
        //**/ add the auxiliary set to the RS training set 
        for (kk = 0; kk < nSamAux; kk++) 
        {
          for (ii = 0; ii < nInputs_; ii++)
             vecXTrain[(nSamples_-hoSize+kk)*nInputs_+ii] = 
                     vecXAux[kk*nInputs_+ii]; 
          vecYTrain[nSamples_-hoSize+kk] = vecYAux[kk*nOutputs_+outputID]; 
        }

        //**/ construct response surface
        status = faPtr->initialize(vecXTrain.getDVector(),
                                   vecYTrain.getDVector());

        //**/ evaluate holdout set
        faPtr->evaluatePoint(hoSize,vecXHO.getDVector(),
                             vecYHO.getDVector());

        //**/ update the predictions 
        ddata = 0;
        for (kk = 0; kk < hoSize; kk++) 
        {
          ind = vecIRand[kk];
          vecYPred[ind] += vecYHO[kk];
          dtemp = VecSamOutputs_[ind*nOutputs_+outputID];
          ddata += pow(dtemp - vecYHO[kk], 2.0);
          if (PABS(dtemp-vecYHO[kk]) > PABS(dtemp-vecYWorst[ind]))
            vecYWorst[ind] = vecYHO[kk];
          if (PABS(dtemp-vecYHO[kk]) < PABS(dtemp-vecYBest[ind]))
            vecYBest[ind] = vecYHO[kk];
        }
        printf(" ===> RMS error = %10.4e\n",sqrt(ddata/hoSize));

        //**/ check to see if all sample points have been held out
        for (kk = 0; kk < nSamples_; kk++) 
          if (vecSCnt[kk] == 0) break;
        if (kk == nSamples_) allTested = 1;

        //**/ update counter
        iterCnt++;
      }
      delete faPtr;
      //**/ take average
      for (kk = 0; kk < nSamples_; kk++) 
        vecYPred[kk] /= (double) vecSCnt[kk];

      //**/ compute errors for average case 
      printAsterisks(PL_INFO, 0);
      printf("*            Average Cross Validation Errors\n");
      printEquals(PL_INFO, 0);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += (VecSamOutputs_[ss*nOutputs_+outputID] - vecYPred[ss]);
      ddata = ddata / nSamples_;
      printf("CV avg error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = (VecSamOutputs_[ss*nOutputs_+outputID] - vecYPred[ss]);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= PABS(VecSamOutputs_[ss*nOutputs_+outputID]);
        ddata += dtemp;
      }
      ddata = ddata / nSamples_;
      printf("CV avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYPred[ss],2.0);
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYPred[ss],2.0);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= pow(VecSamOutputs_[ss*nOutputs_+outputID],2.0);
        ddata += dtemp;
      }
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        if (PABS(vecYPred[ss]-VecSamOutputs_[ss*nOutputs_+outputID]) > ddata)
          ddata = PABS(VecSamOutputs_[ss*nOutputs_+outputID] - vecYPred[ss]);
      printf("CV max error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID];
        if (dtemp != 0 && PABS((vecYPred[ss]-dtemp)/dtemp) > ddata)
          ddata = PABS((dtemp-vecYPred[ss]) / dtemp);
      }
      printf("CV max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      double dmean = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        dmean += (VecSamOutputs_[ss*nOutputs_+outputID]-vecYPred[ss]);
      dmean = dmean / nSamples_;
      printf("CV error mean = %e\n", dmean);
      double dstd = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID] - vecYPred[ss];
        dstd += pow(dtemp-dmean, 2.0);
      }
      dstd = sqrt(dstd / (nSamples_ - 1.0));
      printf("CV error std  = %e\n", dstd);

      //**/ compute errors for worst case 
      printAsterisks(PL_INFO, 0);
      printf("*            Worst Cross Validation Errors\n");
      printEquals(PL_INFO, 0);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += (VecSamOutputs_[ss*nOutputs_+outputID] - vecYWorst[ss]);
      ddata = ddata / nSamples_;
      printf("CV avg error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = (VecSamOutputs_[ss*nOutputs_+outputID] - vecYWorst[ss]);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= PABS(VecSamOutputs_[ss*nOutputs_+outputID]);
        ddata += dtemp;
      }
      ddata = ddata / nSamples_;
      printf("CV avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYWorst[ss],2.0);
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYWorst[ss],2.0);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= pow(VecSamOutputs_[ss*nOutputs_+outputID],2.0);
        ddata += dtemp;
      }
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        if (PABS(vecYWorst[ss]-VecSamOutputs_[ss*nOutputs_+outputID]) > ddata)
          ddata = PABS(VecSamOutputs_[ss*nOutputs_+outputID] - vecYWorst[ss]);
      printf("CV max error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID];
        if (dtemp != 0 && PABS((vecYWorst[ss]-dtemp)/dtemp) > ddata)
          ddata = PABS((dtemp-vecYWorst[ss]) / dtemp);
      }
      printf("CV max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      dmean = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        dmean += (VecSamOutputs_[ss*nOutputs_+outputID]-vecYWorst[ss]);
      dmean = dmean / nSamples_;
      printf("CV error mean = %e\n", dmean);
      dstd = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID] - vecYWorst[ss];
        dstd += pow(dtemp-dmean, 2.0);
      }
      dstd = sqrt(dstd / (nSamples_ - 1.0));
      printf("CV error std  = %e\n", dstd);

      //**/ compute errors for best case 
      printAsterisks(PL_INFO, 0);
      printf("*            Best Cross Validation Errors\n");
      printEquals(PL_INFO, 0);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += (VecSamOutputs_[ss*nOutputs_+outputID] - vecYBest[ss]);
      ddata = ddata / nSamples_;
      printf("CV avg error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = (VecSamOutputs_[ss*nOutputs_+outputID] - vecYBest[ss]);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= PABS(VecSamOutputs_[ss*nOutputs_+outputID]);
        ddata += dtemp;
      }
      ddata = ddata / nSamples_;
      printf("CV avg error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        ddata += pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYBest[ss],2.0);
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = pow(VecSamOutputs_[ss*nOutputs_+outputID]-vecYBest[ss],2.0);
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 0.0)
          dtemp /= pow(VecSamOutputs_[ss*nOutputs_+outputID],2.0);
        ddata += dtemp;
      }
      ddata = sqrt(ddata / nSamples_);
      printf("CV rms error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        if (PABS(vecYBest[ss]-VecSamOutputs_[ss*nOutputs_+outputID]) > ddata)
          ddata = PABS(VecSamOutputs_[ss*nOutputs_+outputID] - vecYBest[ss]);
      printf("CV max error (unscaled) = %e\n", ddata);
      ddata = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID];
        if (dtemp != 0 && PABS((vecYBest[ss]-dtemp)/dtemp) > ddata)
          ddata = PABS((dtemp-vecYBest[ss]) / dtemp);
      }
      printf("CV max error (  scaled) = %e (base=%e)\n",ddata,dtemp);
      dmean = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        dmean += (VecSamOutputs_[ss*nOutputs_+outputID]-vecYBest[ss]);
      dmean = dmean / nSamples_;
      printf("CV error mean = %e\n", dmean);
      dstd = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+outputID] - vecYBest[ss];
        dstd += pow(dtemp-dmean, 2.0);
      }
      dstd = sqrt(dstd / (nSamples_ - 1.0));
      printf("CV error std  = %e\n", dstd);
      printAsterisks(PL_INFO, 0);
      psConfig_.RSExpertModeRestore();
      psConfig_.InteractiveRestore();

      //**/ create visualization file
      if (plotScilab()) fp = fopen("scilabrstestrcv.sci", "w");
      else              fp = fopen("matlabrstestrcv.m", "w");
      sprintf(pString," col 1: simulation data, col 2: rs predicted data");
      fwriteComment(fp, pString);
      sprintf(pString," col 3: worst cases, col 4: best case");
      fwriteComment(fp, pString);
      fprintf(fp, "A = [\n");
      for (ss = 0; ss < nSamples_; ss++)
        fprintf(fp,"%e %e %e %e\n",VecSamOutputs_[ss*nOutputs_+outputID],
                vecYPred[ss], vecYWorst[ss], vecYBest[ss]);
      fprintf(fp, "];\n");
      fwritePlotCLF(fp);
      fwritePlotFigure(fp, 1);
      fprintf(fp, "subplot(1,2,1)\n");
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, A(:,1)-A(:,2), style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,2),10);\n");
        fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "CV Error Plot (unscaled)");
      fwritePlotXLabel(fp, "Error");
      fwritePlotYLabel(fp, "Probabilities");
      fprintf(fp, "subplot(1,2,2)\n");
      fprintf(fp, "xmax = max(A(:,1));\n");
      fprintf(fp, "xmin = min(A(:,1));\n");
      fprintf(fp, "ymax = max(A(:,2));\n");
      fprintf(fp, "ymin = min(A(:,2));\n");
      fprintf(fp, "xmin = min(xmin, ymin);\n");
      fprintf(fp, "xmax = max(xmax, ymax);\n");
      fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
      fprintf(fp, "plot(A(:,1), A(:,2),'x', XX, XX)\n");
      if (plotScilab())
      {
        fprintf(fp, "a = gca();\n");
        fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
      }
      else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Predicted (Means) vs actual data");
      fwritePlotXLabel(fp, "Actual data");
      fwritePlotYLabel(fp, "Predicted data (Average case)");

      fwritePlotFigure(fp, 2);
      fprintf(fp, "subplot(1,2,1)\n");
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, A(:,1)-A(:,3), style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,3),10);\n");
        fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "CV Error Plot (unscaled)");
      fwritePlotXLabel(fp, "Error");
      fwritePlotYLabel(fp, "Probabilities");
      fprintf(fp, "subplot(1,2,2)\n");
      fprintf(fp, "xmax = max(A(:,1));\n");
      fprintf(fp, "xmin = min(A(:,1));\n");
      fprintf(fp, "ymax = max(A(:,3));\n");
      fprintf(fp, "ymin = min(A(:,3));\n");
      fprintf(fp, "xmin = min(xmin, ymin);\n");
      fprintf(fp, "xmax = max(xmax, ymax);\n");
      fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
      fprintf(fp, "plot(A(:,1), A(:,3),'x', XX, XX)\n");
      if (plotScilab())
      {
        fprintf(fp, "a = gca();\n");
        fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
      }
      else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Predicted (Worst case) vs actual data");
      fwritePlotXLabel(fp, "Actual data");
      fwritePlotYLabel(fp, "Predicted data (Worst case)");

      fwritePlotFigure(fp, 3);
      fprintf(fp, "subplot(1,2,1)\n");
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, A(:,1)-A(:,4), style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,xk] = hist(A(:,1)-A(:,4),10);\n");
        fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "CV Error Plot (unscaled)");
      fwritePlotXLabel(fp, "Error");
      fwritePlotYLabel(fp, "Probabilities");
      fprintf(fp, "subplot(1,2,2)\n");
      fprintf(fp, "xmax = max(A(:,1));\n");
      fprintf(fp, "xmin = min(A(:,1));\n");
      fprintf(fp, "ymax = max(A(:,4));\n");
      fprintf(fp, "ymin = min(A(:,4));\n");
      fprintf(fp, "xmin = min(xmin, ymin);\n");
      fprintf(fp, "xmax = max(xmax, ymax);\n");
      fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
      fprintf(fp, "plot(A(:,1), A(:,4),'x', XX, XX)\n");
      if (plotScilab())
      {
        fprintf(fp, "a = gca();\n");
        fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
      }
      else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Predicted (Best case) vs actual data");
      fwritePlotXLabel(fp, "Actual data");
      fwritePlotYLabel(fp, "Predicted data (Best case)");
      fclose(fp);
      if (plotScilab())
           printf("rstest_rcv: plot file in scilabrstestrcv.sci\n");
      else printf("rstest_rcv: plot file in matlabrstestrcv.m\n");
      printf("NOTE: Average case may be better than best case due ");
      printf("to cancellation\n");
      printf("      (output can be positive or negative.)\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rstest_gt 
    //**/ response surface test - generalization test
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstest_gt"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rstest_gt: response surface generalization test.\n");
        printf("Syntax: rstest_gt (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      faType = -1;
      sprintf(pString, "Enter your choice ? ");
      while (faType < 0 || faType > PSUADE_NUM_RS)
      {
        writeFAInfo(outputLevel_);
        faType = getFAType(pString);
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      thresh = -1.0;
      while (thresh < 0.5 || thresh >= 1.0)
      {
        printf("This test partitions each input into two parts and\n");
        printf("uses one part to predict the other part. You can\n");
        printf("decide the percentage of the sample points for \n");
        printf("training (>= 50 percent).\n");
        sprintf(pString,"Enter this percentage in fraction [0.5 - 1) = ");
        thresh = getDouble(pString);
      }
      if (plotScilab()) fp = fopen("scilabrstestgt.sci", "w");
      else              fp = fopen("matlabrstestgt.m", "w");
      sprintf(pString," col 1-m: inputs, col m+1: simulation output");
      fwriteComment(fp, pString);
      sprintf(pString," col m+2: predicted output");
      fwriteComment(fp, pString);
      fwritePlotCLF(fp);

      int count2;
      psVector vecTT,vecAMeans,vecAVars;
      psIVector vecACnts;
      vecXT.setLength(nSamples_*nInputs_);
      vecVT.setLength(nSamples_*nInputs_);
      vecYT.setLength(nSamples_);
      vecWT.setLength(nSamples_);
      vecTT.setLength(nSamples_);
      vecAMeans.setLength(nInputs_*2);
      vecAVars.setLength(nInputs_*2);
      vecACnts.setLength(nInputs_*2);

      printf("Generalization test results:\n");
      for (ii = 0; ii < nInputs_; ii++)
      {
        fprintf(fp, "A%d = [\n", ii+1);
        ddata = thresh * (VecIUpperBs_[ii] - VecILowerBs_[ii]) + 
                VecILowerBs_[ii];
        count = count2 = 0;
        for (ss = 0; ss < nSamples_; ss++)
        { 
          if (VecSamInputs_[ss*nInputs_+ii] <= ddata)
          {
            for (kk = 0; kk < nInputs_; kk++)
              vecXT[count*nInputs_+kk] = VecSamInputs_[ss*nInputs_+kk];
            vecYT[count] = VecSamOutputs_[ss*nOutputs_+outputID];
            count++;
          }
          else
          {
            for (kk = 0; kk < nInputs_; kk++)
              vecVT[count2*nInputs_+kk] = VecSamInputs_[ss*nInputs_+kk];
            vecWT[count2] = VecSamOutputs_[ss*nOutputs_+outputID];
            count2++;
          }
        }
        if (count2 == 0 || count == 0)
        {
          printf("ERROR: for input %d - no training or test point.\n", 
                 ii+1);
          cmdStatus = 1;
          continue;
        }
        FuncApprox *faPtr = genFA(faType, nInputs_, iOne, count);
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        status = faPtr->initialize(vecXT.getDVector(),
                                   vecYT.getDVector());
        faPtr->evaluatePoint(count2, vecVT.getDVector(), 
                            vecTT.getDVector());
        delete faPtr;
        for (ss = 0; ss < count2; ss++)
        {
          for (kk = 0; kk < nInputs_; kk++)
            fprintf(fp, "%e ", vecVT[ss*nInputs_+kk]);
          fprintf(fp, "%e ", vecWT[ss]);
          fprintf(fp, "%e\n", vecTT[ss]);
        }
        vecAMeans[2*ii] = 0.0;
        for (ss = 0; ss < count2; ss++) 
          vecAMeans[2*ii] += (vecWT[ss] - vecTT[ss]);
        vecAMeans[2*ii] /= (double) count2;
        vecAVars[2*ii] = 0.0;
        for (ss = 0; ss < count2; ss++) 
          vecAVars[2*ii] += pow(vecWT[ss]-vecTT[ss]-vecAMeans[2*ii],2.0);
        vecAVars[2*ii] /= (double) count2;
        vecACnts[2*ii] = count2;
 
        ddata = (1.0-thresh) * (VecIUpperBs_[ii]-VecILowerBs_[ii]) + 
                VecILowerBs_[ii];
        count = count2 = 0;
        for (ss = 0; ss < nSamples_; ss++)
        { 
          if (VecSamInputs_[ss*nInputs_+ii] > ddata)
          {
            for (kk = 0; kk < nInputs_; kk++)
              vecXT[count*nInputs_+kk] = VecSamInputs_[ss*nInputs_+kk];
            vecYT[count] = VecSamOutputs_[ss*nOutputs_+outputID];
            count++;
          }
          else
          {
            for (kk = 0; kk < nInputs_; kk++)
              vecVT[count2*nInputs_+kk] = VecSamInputs_[ss*nInputs_+kk];
            vecWT[count2] = VecSamOutputs_[ss*nOutputs_+outputID];
            count2++;
          }
        }
        if (count2 == 0 || count == 0)
        {
          printf("ERROR: for input %d - no training or test point.\n",
                 ii+1);
          cmdStatus = 1;
          continue;
        }
        faPtr = genFA(faType, nInputs_, iOne, count);
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        status = faPtr->initialize(vecXT.getDVector(),
                                   vecYT.getDVector());
        faPtr->evaluatePoint(count2, vecVT.getDVector(), 
                             vecTT.getDVector());
        delete faPtr;
        faPtr = NULL;
        for (ss = 0; ss < count2; ss++)
        {
          for (kk = 0; kk < nInputs_; kk++)
            fprintf(fp, "%e ", vecVT[ss*nInputs_+kk]);
          fprintf(fp, "%e ", vecWT[ss]);
          fprintf(fp, "%e\n", vecTT[ss]);
        }
        fprintf(fp, "];\n");
        vecAMeans[2*ii+1] = 0.0;
        for (ss = 0; ss < count2; ss++) 
          vecAMeans[2*ii+1] += (vecWT[ss] - vecTT[ss]);
        vecAMeans[2*ii+1] /= (double) count2;
        vecAVars[2*ii+1] = 0.0;
        for (ss = 0; ss < count2; ss++) 
          vecAVars[2*ii+1] += 
              pow(vecWT[ss]-vecTT[ss]-vecAMeans[2*ii+1],2.0);
        vecAVars[2*ii+1] /= (double) count2;
        vecACnts[2*ii+1] = count2;
        fprintf(fp, "xmax = max(A%d(:,%d+1));\n", ii+1, nInputs_);
        fprintf(fp, "xmin = min(A%d(:,%d+1));\n", ii+1, nInputs_);
        fprintf(fp, "ymax = max(A%d(:,%d+2));\n", ii+1, nInputs_);
        fprintf(fp, "ymin = min(A%d(:,%d+2));\n", ii+1, nInputs_);
        fprintf(fp, "xmin = min(xmin, ymin);\n");
        fprintf(fp, "xmax = max(xmax, ymax);\n");
        fprintf(fp, "XX   = xmin : xmax-xmin : xmax;\n");
        fprintf(fp, "plot(A%d(:,%d+1), A%d(:,%d+2),'*', XX, XX)\n",
                ii+1, nInputs_, ii+1, nInputs_);
        if (plotScilab())
        {
          fprintf(fp, "a = gca();\n");
          fprintf(fp, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
        }
        else fprintf(fp, "axis([xmin xmax xmin xmax])\n");
        fwritePlotAxes(fp);
        sprintf(pString, "Extrapolated vs actual data (input = %d)",ii+1);
        fwritePlotTitle(fp, pString);
        fwritePlotXLabel(fp, "Actual data");
        fwritePlotYLabel(fp, "Predicted data");
        fprintf(fp, "disp('Press enter to continue to the next input')\n");
        if (ii < nInputs_-1) fprintf(fp, "pause\n");
      }
      for (ii = 0; ii < nInputs_; ii++)
      {
        printf("Input %4d: ", ii+1);
        printf("partition 1 error mean/var = %12.4e %12.4e (%d)\n",
               vecAMeans[2*ii], vecAVars[2*ii], vecACnts[2*ii]); 
        printf("Input %4d: ", ii+1);
        printf("partition 2 error mean/var = %12.4e %12.4e (%d)\n",
               vecAMeans[2*ii+1], vecAVars[2*ii+1], vecACnts[2*ii]+1); 
      }
      fclose(fp);
      if (plotScilab())
           printf("rstest_gt: plot file in scilabrstestgt.sci\n");
      else printf("rstest_gt: plot file in matlabrstestgt.m\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rstest_sd 
    //**/ response surface test - probe prediction uncertainties
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstest_sd"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rstest_sd: compute prediction errors from a response\n");
        printf("           surface (test sample outputs not used).\n");
        printf("Syntax: rstest_sd (no argument needed)\n");
        printf("This command is intended for estimating prediction\n");
        printf("uncertainties in the entire parameter space. As such\n");
        printf("the test sample should be large and space-filling.\n");
        printf("This command computes the prediction uncertainties\n");
        printf("only and no comparison is made against the truth, so\n");
        printf("the output values of the test sample is not needed.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }

      //**/ build response surface
      int faType = -1;
      printf("Select a stochastic response surface method : \n");
      printf("1. MARS with bootstrapped aggregation (MARSB)\n");
      printf("2. Gaussian process (GP)\n");
      printf("3. Kriging\n");
      printf("4. Sum of trees method\n");
      sprintf(pString,"Make your choice: ");
      faType = getInt(1,4,pString);
      switch (faType)
      {
        case 1: faType = PSUADE_RS_MARSB; break;
        case 2: faType = PSUADE_RS_GP3;   break;
        case 3: faType = PSUADE_RS_KR;    break;
        case 4: faType = PSUADE_RS_SOTS;  break;
      }
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      faPtr->setNPtsPerDim(16);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(0);

      //**/ get output information and generate sample 
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString,
         "Sample size for probing the parameter space? (10000 - 1000000) ");
      int nSamp = getInt(10000, 1000000, pString);
      psVector  vecStds;
      vecST.setLength(nSamp);
      vecInps.setLength(nSamp*nInputs_);
      vecOuts.setLength(nSamp);
      vecStds.setLength(nSamp);
      Sampling *samPtr;
      if (nInputs_ < 51)
           samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      else samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
      samPtr->setPrintLevel(0);
      samPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                             VecIUpperBs_.getDVector());
      samPtr->setOutputParams(1);
      samPtr->setSamplingParams(nSamp, -1, -1);
      samPtr->initialize(0);
      samPtr->getSamples(nSamp, nInputs_, iOne, vecInps.getDVector(),
                         vecOuts.getDVector(), vecST.getIVector());
      delete samPtr;
      samPtr = NULL;

      //**/ evaluate
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      status = faPtr->initialize(VecSamInputs_.getDVector(),
                                 vecYT.getDVector());

      faPtr->evaluatePointFuzzy(nSamp, vecInps.getDVector(),
                     vecOuts.getDVector(),vecStds.getDVector());
      delete faPtr;
      faPtr = NULL;

      //**/ generate histogram
      fp = fopen("matlabrstestsd.m", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file matlabrstestsd.m\n");
        cmdStatus = 1;
        continue;
      } 
      fprintf(fp, "Y = [ \n");
      for (jj = 0; jj < nSamp; jj++)
        fprintf(fp, "%16.8e\n", vecStds[jj]);
      fprintf(fp, "];\n");
      fprintf(fp, "[nk,xk]=hist(Y,10);\n");
      fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamp);
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Probability Distribution");
      fwritePlotXLabel(fp, "RS Std. Dev.");
      fwritePlotYLabel(fp, "Probabilities");
      fclose(fp);
      printf("rstest_sd: distribution in matlabrstestsd.m.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ svmfind 
    //**/ find best RS parameters for SVM
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "svmfind"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("svmfind: find optimal parameters for SVM surface fit\n");
        printf("Syntax: svmfind (no argument needed)\n");
        printf("NOTE: This command is for SVM kern option RBF.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
#ifndef HAVE_SVM
      printf("SVM not installed.\n");
#else
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      faType = PSUADE_RS_SVM;

      psConfig_.AnaExpertModeSaveAndReset();
      psConfig_.RSExpertModeSaveAndReset();
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      psConfig_.AnaExpertModeRestore();
      psConfig_.RSExpertModeRestore();

      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      vecYT.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
        vecYT[ii] = VecSamOutputs_[nOutputs_*ii+outputID];
      kk = 2;
      printf("Search for best gamma for RBF kernel.\n");
      printf("SVM tolerance can be adjusted for more accurate result.\n");
      printf("However, small tolerance takes more time to run.\n");
      sprintf(pString, "Enter the desired tolerance (e.g. 1e-4): ");
      double tolerance = getDouble(pString);
      double gamma1=1e-10;
      sprintf(pString,"Enter lower bound of gamma to search (1e-6 - 1e6): ");
      while (gamma1 < 1e-6 || gamma1 > 1e6) gamma1 = getDouble(pString);
      double gamma2=1e-10;
      sprintf(pString,"Enter upper bound of gamma to search (1e-6 - 1e6): ");
      while (gamma2 < 1e-6 || gamma2 > 1e6 || gamma2 < gamma1)
         gamma2 = getDouble(pString);
      Ymin  = 1.0e35;
      double gamma = gamma1;
      double sumErr, maxErr;
      while (gamma <= gamma2)
      {
        targv[0] = (char *) &gamma;
        targv[1] = (char *) &tolerance;
        targv[2] = (char *) &kk;
        faPtr->setParams(3, targv);
        faPtr->initialize(VecSamInputs_.getDVector(),vecYT.getDVector());
        sumErr = maxErr = 0.0;
        for (ii = 0; ii < nSamples_; ii++)
        {
          ddata = VecSamOutputs_[ii*nOutputs_+outputID];
          double *dPtr = VecSamInputs_.getDVector();
          dtemp = faPtr->evaluatePoint(&(dPtr[ii*nInputs_]));
          dtemp = PABS(dtemp - ddata);
          if (dtemp > maxErr) maxErr = dtemp;
          sumErr += (dtemp * dtemp);
        }
        sumErr = sqrt(sumErr);
        printf("svmfind: RBF_gamma = %e, error (L2n,max) = %e %e\n",
               gamma, sumErr, maxErr);
        if (gamma == 1.0e-3 || sumErr < Ymin)
        {
          Ymin = sumErr;
          Xmin = gamma;
        }
        gamma *= 10.0;
      }
      printf("svmfind: best RBF_gamma = %e\n", Xmin);
      targv[0] = (char *) &gamma;
      gamma = Xmin;
      faPtr->setParams(1, targv);
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
#endif
    }

    //**/ -------------------------------------------------------------
    // +++ set_rstype 
    //**/ change rstype in resident sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "set_rstype"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("set_rstype: change the response surface type in the\n");
        printf("   current session already loaded. This command is\n");
        printf("   useful for modifying rstype to do rssobol1....\n");
        printf("Syntax: set_rstype (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      faType = -1;
      while (faType < 0 || faType >= PSUADE_NUM_RS)
      {
        printf("Enter a response surface method : \n");
        int faLimit = writeFAInfo(outputLevel_);
        sprintf(pString, "Enter your choice ? ");
        faType = getInt(0, faLimit, pString);
      }
      psuadeIO_->updateAnalysisSection(-1, -1, faType, -3, -1, -1);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ set_sam_method 
    //**/ change sampling method in resident session
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "set_sam_method"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("set_sam_method: set sampling method in the current\n");
        printf("   session (after 'load'). This command is useful for\n");
        printf("   modifying psuade.in file before interactive 'run'\n");
        printf("Syntax: set_sam_method (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Select a sampling method from below: ");
      SamplingMethod_ = getSamplingMethod(pString);
      sprintf(pString,"Select a sample size: (> 10, < 1000000): ");
      count = getInt(10, 1000000, pString);
      if (sampler_ != NULL) delete sampler_;
      sampler_ = (Sampling *) SamplingCreateFromID(SamplingMethod_);
      sampler_->setPrintLevel(PL_INTERACTIVE);
      sampler_->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                               VecIUpperBs_.getDVector());
      sampler_->setOutputParams(nOutputs_);
      sampler_->setSamplingParams(count, -1, -1);
      sampler_->initialize(0);
      VecSamInputs_.clean();
      VecSamOutputs_.clean();
      VecSamStates_.clean();
      nSamples_ = sampler_->getNumSamples();
      VecSamInputs_.setLength(nInputs_*nSamples_);
      VecSamOutputs_.setLength(nOutputs_*nSamples_);
      VecSamStates_.setLength(nSamples_);
      sampler_->getSamples(nSamples_,nInputs_,nOutputs_,
                   VecSamInputs_.getDVector(),
                   VecSamOutputs_.getDVector(),
                   VecSamStates_.getIVector());
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),NULL); 
      psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,-1,-1,-1);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rstest  or rstest_hs
    //**/ check quality of the response surface using another data set 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstest") || !strcmp(command, "rstest_hs"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rstest_hs: check quality of RS with a holdout set\n");
        printf("Syntax: rstest_hs (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      strcpy(dataFile, "\0");
      sprintf(pString, "Enter test data (PSUADE format) file name : ");
      getString(pString, dataFile);
      dataFile[strlen(dataFile)-1] = '\0';
      status = 0;
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: Test data file %s not found.\n", dataFile);
        status = 1;
        cmdStatus = 1;
      }
      else fclose(fp);
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0)
      {
        printf("ERROR: file %s either not found or in wrong format.\n",
               dataFile);
        cmdStatus = 1;
        continue;
      }
      if (status == 0)
      {
        ioPtr->getParameter("output_noutputs", pPtr);
        kk = pPtr.intData_;
        if (kk > 1)
        {
          printf("ERROR: your test data file should only have 1 output\n");
          printf("       per sample point. Fix and do this again.\n");
          status = 1;
          cmdStatus = 1;
        }
      }
      delete ioPtr;
      if (status == 0)
      {
        sprintf(pString,"Enter output number (1 - %d) = ",nOutputs_);
        outputID = getInt(1, nOutputs_, pString);
        outputID--;
        faType = -1;
        int faLimit = 9;
        while (faType < 0 || faType >= PSUADE_NUM_RS)
        {
          printf("Enter a response surface method : \n");
          faLimit = writeFAInfo(outputLevel_);
          sprintf(pString, "Enter your choice ? ");
          faType = getInt(0, faLimit, pString);
        }
        psIVector vecXSForms;
        if (psConfig_.RSExpertModeIsOn())
        {
          printf("Available input/output transformations :\n");
          printf("0. no transformation.\n");
          printf("1. log transformation on all the inputs.\n");
          printf("2. log transformation on all the outputs.\n");
          printf("3. log transformation on all inputs and outputs.\n");
          sprintf(pString, "Enter your choice ? ");
          int otrans = getInt(0, 3, pString);
          vecXSForms.setLength(2);
          vecXSForms[0] = otrans & 1;
          vecXSForms[1] = otrans & 2;
        }
        else
        {
          vecXSForms.setLength(2);
          vecXSForms[0] = 0;
          vecXSForms[1] = 0;
        }
        //**/ if there is a discrepancy model, set it up
        int discFile=1, nInps, nOuts;
        FuncApprox *faPtr=NULL;
        printf("Name of the discrepancy model file? (enter n if none): ");
        scanf("%s", winput);
        fgets(lineIn,5000,stdin); 
        if (!strcmp(winput, "NONE") || winput[0] == 'n') discFile = 0;
        else
        {
          PsuadeData *ioPtr = new PsuadeData();
          status = ioPtr->readPsuadeFile(winput);
          if (status == 0)
          {
            ioPtr->getParameter("input_ninputs", pPtr);
            nInps = pPtr.intData_;
            if (nInps < nInputs_)
            {
              printf("ERROR: your sample has %d inputs but\n", nInputs_);
              printf("       your discrepancy model has %d inputs.\n",
                     nInps);
              delete ioPtr;
              ioPtr = NULL;
              continue;
            }
            ioPtr->getParameter("output_noutputs", pPtr);
            nOuts = pPtr.intData_;
            if (nOuts > 1)
            {
              printf("INFO: The discrepancy model has nOutputs > 1.\n");
              printf("      This is currently not supported.\n");
              delete ioPtr;
              ioPtr = NULL;
              continue;
            }
            faPtr = genFAInteractive(ioPtr, 3);
            delete ioPtr;
            ioPtr = NULL;
          }
          else
          {
            printf("ERROR: in reading the discrepancy model file %s.\n",
                   winput);
            discFile = 0;
            delete ioPtr;
            ioPtr = NULL;
            continue;
          }
        }
        psVector vecSO;
        vecSO.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++) 
          vecSO[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
        if (discFile == 1)
        {
          if (psConfig_.MasterModeIsOn())
          {
            fp = fopen("rstest_hs_data", "w");
            if (fp != NULL)
            {
              fprintf(fp,"%% Col 2: RS estimate, \n");
              fprintf(fp,"%% Col 3: model correction.\n");
            }
          }
          else fp = NULL; 
               
          double *dPtr = VecSamInputs_.getDVector();
          for (ii = 0; ii < nSamples_; ii++)
          {
            dtemp = 
              faPtr->evaluatePoint(&(dPtr[ii*nInputs_]));
            if (fp != NULL) 
            {
              fprintf(fp, "%6d ",ii+1);
              for (jj = 0; jj < nInputs_; jj++)
                fprintf(fp, "%12.4e ",VecSamInputs_[ii*nInputs_+jj]);
              fprintf(fp, " %12.4e %12.4e\n",vecSO[ii],dtemp);
            }
            vecSO[ii] += dtemp;
          }
          if (fp != NULL)
          {
            printf("rstest_hs_data has simulations and discrepancies.\n");
            fclose(fp);
          }
        }
        int analysisMethod = PSUADE_ANA_RSFA;
        AnalysisManager *anaManager = new AnalysisManager();
        anaManager->setup(analysisMethod, faType);
        anaManager->loadLogXsformFlags(2, vecXSForms.getIVector());
        aData aPtr;
        aPtr.printLevel_ = outputLevel_;
        aPtr.sampleInputs_ = VecSamInputs_.getDVector();
        aPtr.sampleOutputs_ = vecSO.getDVector();
        aPtr.nSamples_ = nSamples_;
        aPtr.nInputs_ = nInputs_;
        aPtr.nOutputs_ = 1;
        aPtr.iLowerB_ = VecILowerBs_.getDVector();
        aPtr.iUpperB_ = VecIUpperBs_.getDVector();
        aPtr.outputID_ = 0;
        aPtr.faType_ = faType;
        aPtr.regWgtID_ = -1;
        aPtr.ioPtr_ = NULL;
        aPtr.sampleStates_ = VecSamStates_.getIVector();
        strcpy(pString, "validate");
        targv[0] = (char *) pString;
        targv[1] = (char *) &aPtr;
        targv[2] = (char *) dataFile;
        char errFile[5001];
        if (plotScilab()) strcpy(errFile, "RSTest_hs.sci");
        else              strcpy(errFile, "RSTest_hs.m");
        targv[3] = (char *) errFile;
        kk = 4;
        anaManager->specialRequest(PSUADE_ANA_RSFA, kk, targv);
        delete anaManager;
        aPtr.sampleOutputs_ = NULL;
        if (discFile == 1)
        {
          delete faPtr;
          faPtr = NULL;
        }
        if (plotScilab())
             printf("rstest plot file RSTest_hs.sci has been created\n");
        else printf("rstest plot file RSTest_hs.m has been created\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ rstgen 
    //**/ generate sample for response surface test
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rstgen"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command creates a sample for response surface\n");
        printf("tests. Possible samples are factorial or fractional\n");
        printf("factorial to cover the corners (uniform distribution).\n");
        printf("Syntax: rstgen (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a test sample for response ");
      printf("surface analysis.\n");
      printf("Possible samples are factorial or fractional ");
      printf("factorial to cover\n");
      printf("the corners (uniform distribution).\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      printf("You have the following options:\n");
      printf("(1) If the number of inputs m is small, you can use\n");
      printf("    factorial sampling where nSamples = 2^m.\n");
      printf("(2) You can use fractional factorial sampling with\n");
      printf("    resolution V which needs nSamples=128 for m up to 11.\n");
      printf("(3) If the number of inputs m is moderate, you can use\n");
      printf("    fractional factorial sampling with resolution IV\n");
      printf("    which needs nSamples = 32 for m up to 15 and 64 for\n");
      printf("    m up to 32.\n");
      sprintf(pString, "Which option? (1 - 3) ");
      ind = getInt(1, 3, pString);
      Sampling *sampPtr = NULL;
      if (ind == 1)
      {
        sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
        count = (int) pow(2.0, 1.0*nInputs_);
        if (count > 10000)
          printf("nSamples = %d may be too large.\n", nSamples_);
      }
      else if (ind == 2)
      {
        sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF5);
        count = 128;
        if (nInputs_ > 11)
        {
          printf("nInputs = %d too large for FF5.\n", nInputs_);
          cmdStatus = 1;
          continue;
        }
      }
      else if (ind == 3)
      {
        sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
        count = 64;
        if (nInputs_ > 32)
        {
          printf("nInputs = %d too large for FF4.\n", nInputs_);
          cmdStatus = 1;
          continue;
        }
      }
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(count, -1, 0);
      sampPtr->initialize(0);
      count = sampPtr->getNumSamples();
      vecXT.setLength(count*nInputs_);
      vecYT.setLength(count*nOutputs_);
      vecST.setLength(count);
      sampPtr->getSamples(count,nInputs_,nOutputs_,vecXT.getDVector(),
                          vecYT.getDVector(),vecST.getIVector());
      delete sampPtr;
      sampPtr = NULL;

      psVector  vecIMeans, vecIStds;
      psIVector vecIPDFS;
      vecIPDFS.setLength(nInputs_);
      vecIMeans.setLength(nInputs_);
      vecIStds.setLength(nInputs_);
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        vecIPDFS[ii] = 0;
        vecIMeans[ii] = 0;
        vecIStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->updateInputSection(count, nInputs_, NULL, 
               VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
               vecXT.getDVector(), StrInpNames_.getStrings(), 
               vecIPDFS.getIVector(), vecIMeans.getDVector(), 
               vecIStds.getDVector(), iCMat);
      delete iCMat;
      ioPtr->updateOutputSection(count, nOutputs_, vecYT.getDVector(), 
                   vecST.getIVector(), StrOutNames_.getStrings()); 
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
      ioPtr->writePsuadeFile("psuade_rs_test.data",0);
      delete ioPtr;
      printf("The test sample is in file psuade_rs_test.data.\n"); 
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsint 
    //**/ numerical integration
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsint"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsint: compute volume under sample surface)\n");
        printf("Syntax: int (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      //**/ create response surface
      faFlag = 3;
      psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID,-1);
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      faPtr->setOutputLevel(outputLevel_);
      Sampling *sampPtr = NULL;
      if (nInputs_ <= 51)
           sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      else sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(1);
      kk = 1000000;
      sampPtr->setSamplingParams(kk, 1, 0);
      sampPtr->initialize(0);

      vecXT.setLength(kk*nInputs_);
      vecYT.setLength(kk);
      vecST.setLength(kk);
      sampPtr->getSamples(kk, nInputs_, 1, vecXT.getDVector(), 
                          vecYT.getDVector(), vecST.getIVector());
      faPtr->evaluatePoint(kk,vecXT.getDVector(),vecYT.getDVector());
      ddata = 0.0;
      for (ii = 0; ii < kk; ii++) ddata += vecYT[ii];
      ddata /= (double) kk;
      for (ii = 0; ii < nInputs_; ii++) 
        ddata *= (VecIUpperBs_[ii] - VecILowerBs_[ii]);
      printf("rsint: volume under the function = %16.8e.\n",ddata);
      delete faPtr;
      faPtr = NULL;
      delete sampPtr;
      sampPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsvol 
    //**/ calculate volume of constrained space
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsvol"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsvol: calculate percentage of volume in the");
        printf(" constrained space\n");
        printf("Syntax: rsvol (no argument needed)\n");
        printf("This command is useful for estimating how much of\n");
        printf("parameter space has been cut out due to constraints.\n");
        printf("It will report the percentage of the parameter space\n");
        printf("that is feasible in view of constraints.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      Ymax = - 1.0e35;
      Ymin =   1.0e35;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
          Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
          Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
      }
      sprintf(pString,"Enter the lower bound constraint (Ymin=%e) : ",
              Ymin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper bound constraint (Ymax=%e) : ",
              Ymax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        printf("       lower bound = %e\n", threshL);
        printf("       upper bound = %e\n", threshU);
        continue;
      }
      //**/ create response surface
      faFlag = 3;
      psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID,-1); 
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      faPtr->setOutputLevel(outputLevel_);
      //**/ create a large sample
      Sampling *sampPtr = NULL;
      if (nInputs_ <= 51) 
           sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      else sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(1);
      kk = 1000000;
      sampPtr->setSamplingParams(kk, 1, 0);
      sampPtr->initialize(0);
      vecXT.setLength(kk*nInputs_);
      vecYT.setLength(kk);
      vecST.setLength(kk);
      sampPtr->getSamples(kk, nInputs_, 1, vecXT.getDVector(), 
                          vecYT.getDVector(), vecST.getIVector());
      count = 0;
      double *dPtrX = vecXT.getDVector();
      for (ii = 0; ii < kk; ii++)
      {
        dtemp = faPtr->evaluatePoint(&dPtrX[ii*nInputs_]);
        if (dtemp >= threshL && dtemp <= threshU) count++; 
      } 
      printf("rsvol: percentage inside the constrained region = %5.2f%%.\n",
             100.0 * count / kk);
      delete faPtr;
      faPtr = NULL;
      delete sampPtr;
      sampPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ opca 
    //**/ output principal component analysis 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "opca"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("opca: principal component analysis on the outputs\n");
        printf("Syntax: opca (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      aData aPtr;
      aPtr.nSamples_ = nSamples_;
      aPtr.nInputs_  = nInputs_;
      aPtr.nOutputs_ = nOutputs_;
      aPtr.iLowerB_  = VecILowerBs_.getDVector();
      aPtr.iUpperB_  = VecIUpperBs_.getDVector();
      aPtr.sampleInputs_ = VecSamInputs_.getDVector();
      aPtr.sampleOutputs_ = VecSamOutputs_.getDVector();
      aPtr.sampleStates_ = VecSamStates_.getIVector();
      aPtr.printLevel_ = outputLevel_;
      AnalysisManager *anaManager = new AnalysisManager();
      int analysisMethod = PSUADE_ANA_PCA;
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      psuadeIO_->getParameter("output_noutputs", pPtr);
      kk = pPtr.intData_;
      if (kk != nOutputs_) 
      {
        nOutputs_ = kk;
        pONames.clean();
        psuadeIO_->getParameter("output_names", pONames);
        StrOutNames_.load(nOutputs_, (const char **) pONames.strArray_);
        VecSamOutputs_.clean();
        psuadeIO_->getParameter("output_sample", pPtr);
        VecSamOutputs_.load(nSamples_*nOutputs_,pPtr.dbleArray_);
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ 1stest 
    //**/ One sample test (chi-squred, distribution fitting) 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "1stest"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("1stest: perform one-sample tests (chi-squared, etc.)\n");
        printf("Syntax: lstest (no argument needed)\n");
        continue;
      }
      int analysisMethod = PSUADE_ANA_1SAMPLE;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ 2stest 
    //**/ 2 sample test (Kolmogorov Smirnov, Mann-Whitney, T-test) 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "2stest"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("2stest: two-sample tests (Kolmogorov Smirnov, etc.)\n");
        printf("Syntax: 2stest (no argument needed)\n");
        continue;
      }
      int analysisMethod = PSUADE_ANA_2SAMPLE;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod);
      delete anaManager;
      fgets(lineIn,5000,stdin); 
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iplot1
    //**/ 1D input scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot1"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot1: plot the sample points in one parameter space\n");
        printf("Syntax: iplot1 (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a scatter plot for one input. The ");
      printf("X-axis has the\n");
      printf("sample number and the Y-axis has the selected input values.\n");
      printf("NOTE: if you desire Scilab plot, call scilab before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      iplot1 = -1;
      sprintf(pString, "Enter the input number (1 - %d) : ", nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilabiplt1.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabiplt1.sci.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
              VecSamStates_[sInd] != 1)
            fprintf(fp, "%e 0\n",VecSamInputs_[sInd*nInputs_+iplot1]);
          else
            fprintf(fp, "%e 1\n",VecSamInputs_[sInd*nInputs_+iplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "n = %d;\n", nSamples_);
        fprintf(fp, "ia = find(X(:,2) == 0);\n");
        fprintf(fp, "if (length(ia) > 0)\n");
        fprintf(fp, "   plot(ia, X(ia,1),'rX','markerSize',13)\n");
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "ia = find(X(:,2) == 1);\n");
        fprintf(fp, "if (length(ia) > 0)\n");
        fprintf(fp, "   plot(ia, X(ia,1),'b*','markerSize',13)\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
        fprintf(fp, "minX = min(X(:,1)); maxX = max(X(:,1));\n");
        fprintf(fp, "a = gca();\n");
        fprintf(fp, "a.data_bounds=[0,minX;%d,maxX];\n",nSamples_);
        fwritePlotYLabel(fp, StrInpNames_[iplot1]);
        fwritePlotXLabel(fp, "Sample Number");
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "1D Input Data Plot");
        fclose(fp);    
        printf("scilabiplt1.sci is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
      else
      {
        fp = fopen("matlabiplt1.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabiplt1.m.\n");
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
              VecSamStates_[sInd] != 1)
            fprintf(fp, "%e 0\n",VecSamInputs_[sInd*nInputs_+iplot1]);
          else
            fprintf(fp, "%e 1\n",VecSamInputs_[sInd*nInputs_+iplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "iset = find(X(:,2)==0);\n");
        fprintf(fp, "plot(iset,X(iset,1),'rX','MarkerSize',13)\n");
        fprintf(fp, "hold on\n");
        fprintf(fp, "iset = find(X(:,2)==1);\n");
        fprintf(fp, "plot(iset,X(iset,1),'b*','MarkerSize',13)\n");
        fprintf(fp, "hold off\n");
        fprintf(fp, "axis([0 %d min(X(:,1)) max(X(:,1))])\n", nSamples_);
        fwritePlotYLabel(fp, StrInpNames_[iplot1]);
        fwritePlotXLabel(fp, "Sample Number");
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "1D Input Data Plot");
        fclose(fp);    
        printf("matlabiplt1.m is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ iplot2 
    //**/ 2D input scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot2: plot sample points in two-parameter space\n");
        printf("Syntax: iplot2 (no argument needed)\n");
        printf("This command is useful for examining where in the.\n");
        printf("parameter space failed runs are.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a scatter plot for 2 selected ");
      printf("inputs. The X-axis\n");
      printf("and Y-axis have the 2 selected input values.\n");
      printf("NOTE: if you desire Scilab plot, call scilab ");
      printf("before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue;

      if (nInputs_ < 2)
      {
        printf("ERROR: iplot2 requires 2 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }

      iplot1 = iplot2 = -1;
      sprintf(pString, 
              "Enter the input for x axis (1 - %d) : ",nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;
      sprintf(pString, 
              "Enter the input for y axis (1 - %d) : ",nInputs_);
      iplot2 = getInt(1, nInputs_, pString);
      iplot2--;
      if (iplot2 == iplot1)
        printf("WARNING: same input index %d for x and y axes.\n",
               iplot1+1);

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilabiplt2.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabiplt2.sci.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED)
            fprintf(fp, "%e %e 0\n",VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2]);
          else
            fprintf(fp, "%e %e 1\n",VecSamInputs_[sInd*nInputs_+iplot1],
                     VecSamInputs_[sInd*nInputs_+iplot2]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "ia = find(X(:,3) == 0);\n");
        fprintf(fp, "if (length(ia) > 0)\n");
        fprintf(fp, "   plot(X(ia,1),X(ia,2),'rX')\n");
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "ia = find(X(:,3) == 1);\n");
        fprintf(fp, "if (length(ia) > 0)\n");
        fprintf(fp, "   plot(X(ia,1),X(ia,2),'b*')\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
        fprintf(fp, "a=gca();\n");
        fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                VecILowerBs_[iplot1], VecILowerBs_[iplot2],
                VecIUpperBs_[iplot1], VecIUpperBs_[iplot2]);
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotTitle(fp, "2D Input Data Plot");
        fwritePlotAxes(fp);
        fprintf(fp,"disp('red X: failed runs.')\n");
        fclose(fp);    
        printf("scilabiplt2.sci now has input scatter plots.\n");
        cmdStatus = 0;
      }
      else
      {
        fp = fopen("matlabiplt2.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabiplt2.m.\n");
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "ranflag = 0;\n");
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_] > 0.5*PSUADE_UNDEFINED ||
              VecSamStates_[sInd] != 1)
            fprintf(fp, "%24.16e %24.16e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2]);
          else
            fprintf(fp, "%24.16e %24.16e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2]);
        }
        fprintf(fp,"];\n");
        fprintf(fp,"iset = find(X(:,3)==0);\n");
        fprintf(fp,"plot(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
        fprintf(fp,"X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
        fprintf(fp,"'markerSize',13)\n");
        fprintf(fp,"hold on\n");
        fprintf(fp,"iset = find(X(:,3)==1);\n");
        fprintf(fp,"plot(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100),");
        fprintf(fp,"X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
        fprintf(fp,"'markerSize',13)\n");
        fprintf(fp,"hold off\n");
        fprintf(fp,"axis([%e %e %e %e])\n", VecILowerBs_[iplot1], 
                VecIUpperBs_[iplot1],VecILowerBs_[iplot2], 
                VecIUpperBs_[iplot2]);
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "2D Input Data Plot");
        fprintf(fp,"disp('red X: failed runs.')\n");
        fclose(fp);    
        printf("matlabiplt2.m is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ splot2
    //**/ 3D input/output scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "splot2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("splot2: scatter plot 2-inputs/1-output in 3D\n");
        printf("Syntax: splot2 (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a scatter plot for two inputs ");
      printf("and one output.\n");
      printf("The inputs are in the X and Y axes, and the output ");
      printf("is in the Z axis.\n");
      printf("NOTE: if you desire Scilab plot, call scilab before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ < 2)
      {
        printf("ERROR: splot2 requires 2 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }

      iplot1 = iplot2 = -1;
      sprintf(pString, "X-axis input ? (1 - %d) ", nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;
      sprintf(pString, "Y-axis input ? (1 - %d) ", nInputs_);
      iplot2 = getInt(1, nInputs_, pString);
      iplot2--;
      if (iplot1 == iplot2)
        printf("WARNING: X- and Y-axes use the same input %d\n",
               iplot1+1);
      sprintf(pString, "Z-axis output ? (1 - %d) : ",nOutputs_);
      oplot1 = getInt(1, nOutputs_, pString);
      oplot1--;

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilabsp2.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabsp2.sci.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          dtemp = VecSamOutputs_[sInd*nOutputs_+oplot1];
          if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
            fprintf(fp, "%e %e %e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
          else
            fprintf(fp, "%e %e %e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "f = gcf();\n");
        fprintf(fp, "f.color_map = jetcolormap(3);\n");
        fprintf(fp, "drawlater\n");
        fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
        fprintf(fp, "if (length(ia1) > 0)\n");
        fprintf(fp, "   param3d1([X(ia1,1)' ; X(ia1,1)'],");
        fprintf(fp, "[X(ia1,2)' ; X(ia1,2)'],[X(ia1,3)' ; X(ia1,3)'])\n");
        fprintf(fp, "   e = gce();\n");
        fprintf(fp, "   e.children.mark_mode = \"on\";\n");
        fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
        fprintf(fp, "   e.children.mark_style = 10;\n");
        fprintf(fp, "   e.children.mark_size = 6;\n");
        fprintf(fp, "   for i = 1:length(e.children)\n");
        fprintf(fp, "      e.children(i).mark_foreground = 3;\n");
        fprintf(fp, "   end\n");
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
        fprintf(fp, "if (length(ia2) > 0)\n");
        fprintf(fp, "   param3d1([X(ia2,1)';X(ia2,1)'],[X(ia2,2)';");
        fprintf(fp, "X(ia2,2)'],[X(ia2,3)';X(ia2,3)'])\n");
        fprintf(fp, "   e = gce();\n");
        fprintf(fp, "   e.children.mark_mode = \"on\";\n");
        fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
        fprintf(fp, "   e.children.mark_style = 10;\n");
        fprintf(fp, "   e.children.mark_size = 6;\n");
        fprintf(fp, "   for i = 1:length(e.children)\n");
        fprintf(fp, "      e.children(i).mark_foreground = 1;\n");
        fprintf(fp, "   end\n");
        fprintf(fp, "end\n");
        fprintf(fp, "drawnow\n");
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrOutNames_[oplot1]);
        fwritePlotTitle(fp, "3D 2-Input/1-Output Data Plot");
        fwritePlotAxes(fp);
        fprintf(fp,"disp('red  *: failed runs.')\n");
        fprintf(fp,"disp('blue X: good   runs.')\n");
        fclose(fp);
        printf("scilabsp2.sci is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
      else
      {
        fp = fopen("matlabsp2.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabsp2.m.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          dtemp = VecSamOutputs_[sInd*nOutputs_+oplot1];
          if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
            fprintf(fp, "%e %e %e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
          else
            fprintf(fp, "%e %e %e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "iset = find(X(:,4)==0);\n");
        fprintf(fp, "if size(iset) > 0\n");
        fprintf(fp, "plot3(X(iset,1),X(iset,2),X(iset,3),'rX',");
        fprintf(fp, "'MarkerSize',13)\n");
        fprintf(fp, "hold on\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "iset = find(X(:,4)==1);\n");
        fprintf(fp, "if size(iset) > 0\n");
        fprintf(fp, "plot3(X(iset,1),X(iset,2),X(iset,3),'b*',");
        fprintf(fp, "'MarkerSize',13)\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "hold off\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrOutNames_[oplot1]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "3D 2-Input/1-Output Data Plot");
        fprintf(fp,"disp('red  *: failed runs.')\n");
        fprintf(fp,"disp('blue X: good   runs.')\n"); 
        fclose(fp);    
        printf("matlabsp2.m is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ splot3
    //**/ 3D 3-input/output scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "splot3"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("splot3: scatter plot 3-inputs/1-output in 3D\n");
        printf("Syntax: splot3 (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a scatter plot for 3 inputs ");
      printf("and 1 output.\n");
      printf("The inputs are in the X, Y, and Z axes, and the ");
      printf("output values are\n");
      printf("expressed by the ball sizes.\n");
      printf("NOTE: if you desire Scilab plot, call scilab ");
      printf("before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 3)
      {
        printf("ERROR: splot3 requires 3 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "X-axis input ? (1 - %d) ", nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;
      sprintf(pString, "Y-axis input ? (1 - %d) ", nInputs_);
      iplot2 = getInt(1, nInputs_, pString);
      iplot2--;
      if (iplot1 == iplot2)
        printf("WARNING: X- and Y-axes use the same input %d\n",
               iplot1+1);
      sprintf(pString, "Z-axis input ? (1 - %d) ", nInputs_);
      iplot3 = getInt(1, nInputs_, pString);
      iplot3--;
      if (iplot1 == iplot3 || iplot2 == iplot3)
        printf("WARNING: input %d is used more than once.\n",
               iplot3+1);
      sprintf(pString, "Which output ? (1 - %d) : ",nOutputs_);
      oplot1 = getInt(1, nOutputs_, pString);
      oplot1--;

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        printf("ERROR: this command is not supported on scilab.\n");
        cmdStatus = 1;
      }
      else
      {
        fp = fopen("matlabsp3.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabsp3.m.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "scale = 1;\n");
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          dtemp = VecSamOutputs_[sInd*nOutputs_+oplot1];
          if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
            fprintf(fp, "%e %e %e %e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
          else
            fprintf(fp, "%e %e %e %e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3],
                    VecSamOutputs_[sInd*nOutputs_+oplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "ymin = min(X(:,4));\n");
        fprintf(fp, "ywid = max(X(:,4)) - ymin;\n");
        fprintf(fp, "if (ywid == 0)\n");
        fprintf(fp, "  ywid = 1;\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "ywid = 1 / ywid;\n");
        fprintf(fp, "for ii = 1 : %d\n", nSamples_);
        fprintf(fp, "  if X(ii,5) == 0\n");
        fprintf(fp, "    dsize = (X(ii,4) - ymin) * ywid;\n");
        fprintf(fp, "    isize = ceil((dsize * 100 + 10)*scale);\n");
        if (plotScilab())
        {
          fprintf(fp, "    scatter3(X(ii,1),X(ii,2),X(ii,3),'r.',");
          fprintf(fp, "'MarkerSize',isize)\n");
        }
        else
        {
          fprintf(fp, "    plot3(X(ii,1),X(ii,2),X(ii,3),'r.',");
          fprintf(fp, "'MarkerSize',isize)\n");
        }
        fprintf(fp, "  else\n");
        fprintf(fp, "    dsize = (X(ii,4) - ymin) * ywid;\n");
        fprintf(fp, "    isize = ceil((dsize * 100 + 10)*scale);\n");
        if (plotScilab())
        {
          fprintf(fp, "    scatter3(X(ii,1),X(ii,2),X(ii,3),'b.',");
          fprintf(fp, "'MarkerSize',isize)\n");
        }
        else
        {
          fprintf(fp, "    plot3(X(ii,1),X(ii,2),X(ii,3),'b.',");
          fprintf(fp, "'MarkerSize',isize)\n");
        }
        fprintf(fp, "  end;\n");
        fprintf(fp, "  if ii == 1\n");
        if (plotScilab())
             fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        else fprintf(fp, "    hold on\n");
        fprintf(fp, "  end;\n");
        fprintf(fp, "end;\n");
        if (plotScilab())
             fprintf(fp, "   set(gca(),\"auto_clear\",\"on\")\n");
        else fprintf(fp, "hold off\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrInpNames_[iplot3]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "3D 3-Input/1-Output Data Plot");
        if (plotMatlab())
        {
          fprintf(fp,"disp('red  *: failed runs.')\n");
          fprintf(fp,"disp('blue X: good   runs.')\n"); 
          fprintf(fp,"disp('NOTE: dot sizes are output magnitudes.')\n");
          fprintf(fp,"disp('NOTE: use scale to adjust relative dot size')\n");
        }
        fclose(fp);    
        printf("matlabsp3.m is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ splot3m
    //**/ 3D 3-input/output scatter plot in movie mode
    //**/ (work only for factorial)
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "splot3m"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("splot3m: scatter plot 3-inputs/1-output in movie\n");
        printf("         It only works for 3-input factorial samples.\n");
        printf("Syntax: splot3m (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a movie of scatter plots for 3 inputs ");
      printf("and 1 output.\n");
      printf("Two inputs are in the X and Y axes, and the output in ");
      printf("in the Z axis.");
      printf("The third output is in the time axis.\n");
      printf("NOTE: if you desire Scilab plot, call scilab before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ != 3)
      {
        printf("ERROR: splot3m only works with 3 inputs because ");
        printf("it does not make sense\n");
        printf("       with more than 3 inputs (the other ");
        printf("inputs are hidden).\n");
        cmdStatus = 1;
        continue;
      }

      //**/ check sample is factorial
      psIVector vecDims;
      FactorialSampling *factSampler = new FactorialSampling();
      status = factSampler->checkFactorial(nSamples_,nInputs_,
                                           VecSamInputs_,vecDims);
      delete factSampler;
      if (status != 0) 
      {
        printf("splot3m ERROR: sample not factorial.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Select output ? (1 - %d) : ",nOutputs_);
      oplot1 = getInt(1, nOutputs_, pString);
      oplot1--;

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        printf("ERROR: this command is not supported on scilab.\n");
      }
      else
      {
        fp = fopen("matlabsp3m.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabsp3m.m.\n");
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          fprintf(fp, "%e %e %e %e\n", VecSamInputs_[sInd*nInputs_],
                  VecSamInputs_[sInd*nInputs_+1],
                  VecSamInputs_[sInd*nInputs_+2],
                  VecSamOutputs_[sInd*nOutputs_+oplot1]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "n1 = %d;\n", vecDims[0]);
        fprintf(fp, "n2 = %d;\n", vecDims[1]);
        fprintf(fp, "n3 = %d;\n", vecDims[2]);
        fprintf(fp, "zmin = min(X(:,4));\n");
        fprintf(fp, "zmax = max(X(:,4));\n");
        fprintf(fp, "if (zmax - zmin == 0)\n");
        fprintf(fp, "  zmax = zmax + 1;\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "x1min = min(X(:,1));\n");
        fprintf(fp, "x1max = max(X(:,1));\n");
        fprintf(fp, "x2min = min(X(:,2));\n");
        fprintf(fp, "x2max = max(X(:,2));\n");
        fprintf(fp, "for ii = 1 : n3\n");
        fprintf(fp, "  i1 = (ii - 1) * n1 * n2 + 1;\n");
        fprintf(fp, "  i2 = ii * n1 * n2;\n");
        fprintf(fp, "  plot3(X(i1:i2,1),X(i1:i2,2),X(i1:i2,4),'bx',");
        fprintf(fp, "'MarkerSize',12)\n");
        fprintf(fp, "  axis([x1min x1max x2min x2max zmin zmax]);\n");
        fwritePlotAxes(fp);
        fwritePlotXLabel(fp, StrInpNames_[0]);
        fwritePlotYLabel(fp, StrInpNames_[1]);
        fwritePlotZLabel(fp, StrOutNames_[oplot1]);
        fwritePlotTitle(fp, "4D 3-Input/1-Output Data Plot");
        fprintf(fp, "  disp(['Input3 = ', num2str(X(i1,3))])\n");
        fprintf(fp, "  pause\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "hold off\n");
        fclose(fp);    
        printf("matlabsp3m.m is now ready for input scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ iplot3
    //**/ 3D input scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot3"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot3: plot the sample points in 3-parameter space\n");
        printf("Syntax: iplot3 (no argument needed)\n");
        printf("This command is useful for examining where in the.\n");
        printf("parameter space failed runs are.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a 3D scatter plot for 3 ");
      printf("selected inputs on the\n");
      printf("X-, Y-, and Z-axes.\n");
      printf("NOTE: if you desire Scilab plot, call scilab before ");
      printf("doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 3)
      {
        printf("ERROR: iplot3 requires 3 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }

      iplot1 = iplot2 = iplot3 = -1;
      sprintf(pString,"Enter the input for x axis (1 - %d) : ",nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;
      sprintf(pString,"Enter the input for y axis (1 - %d) : ",nInputs_);
      iplot2 = getInt(1, nInputs_, pString);
      iplot2--;
      if (iplot1 == iplot2)
        printf("WARNING: input %d is used more than once.\n",
               iplot1+1);
      sprintf(pString,"Enter the input for z axis (1 - %d) : ",
              nInputs_);
      iplot3 = getInt(1, nInputs_, pString);
      iplot3--;
      if (iplot3 == iplot1 || iplot3 == iplot2)
        printf("WARNING: input %d is used more than once.\n",
               iplot3+1);

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilabiplt3.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabiplt3.sci.\n");
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          dtemp = VecSamOutputs_[sInd*nOutputs_];
          if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
            fprintf(fp, "%e %e %e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3]);
          else
            fprintf(fp, "%e %e %e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "f = gcf();\n");
        fprintf(fp, "f.color_map = jetcolormap(3);\n");
        fprintf(fp, "drawlater\n");
        fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
        fprintf(fp, "if (length(ia1) > 0)\n");
        fprintf(fp, "   param3d1([X(ia,1)' ; X(ia,1)'],");
        fprintf(fp, "[X(ia,2)' ; X(ia,2)'],[X(ia,3)' ; X(ia,3)'])\n");
        fprintf(fp, "   e = gce();\n");
        fprintf(fp, "   e.children.mark_mode = \"on\";\n");
        fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
        fprintf(fp, "   e.children.mark_style = 10;\n");
        fprintf(fp, "   e.children.mark_size = 6;\n");
        fprintf(fp, "   for i = 1:length(e.children)\n");
        fprintf(fp, "      e.children(i).mark_foreground = 3;\n");
        fprintf(fp, "   end\n");
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
        fprintf(fp, "if (length(ia2) > 0)\n");
        fprintf(fp, "   param3d1([X(ia2,1)';X(ia2,1)'],[X(ia2,2)';");
        fprintf(fp, "X(ia2,2)'],[X(ia2,3)';X(ia2,3)'])\n");
        fprintf(fp, "   e = gce();\n");
        fprintf(fp, "   e.children.mark_mode = \"on\";\n");
        fprintf(fp, "   e.children.mark_size_unit = \"point\";\n");
        fprintf(fp, "   e.children.mark_style = 10;\n");
        fprintf(fp, "   e.children.mark_size = 6;\n");
        fprintf(fp, "   for i = 1:length(e.children)\n");
        fprintf(fp, "      e.children(i).mark_foreground = 1;\n");
        fprintf(fp, "   end\n");
        fprintf(fp, "end\n");
        fprintf(fp, "drawnow\n");
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrInpNames_[iplot3]);
        fwritePlotTitle(fp, "3D Input Data Plot");
        fwritePlotAxes(fp);
        fclose(fp);    
        printf("scilabiplt3.sci now has input scatter plots.\n");
        cmdStatus = 0;
      }
      else
      {
        fp = fopen("matlabiplt3.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabiplt3.m.\n");
          continue;
        }
        sprintf(pString," ranflag: to distinguish identical points");
        fwriteComment(fp, pString);
        sprintf(pString,"     by adding a small perturbation (when on)");
        fwriteComment(fp, pString);
        fprintf(fp, "ranflag  = 0;\n");
        sprintf(pString," set cvFlag to 1 to use convex hull");
        fwriteComment(fp, pString);
        fprintf(fp, "cvFlag = 0;\n");
        fwritePlotCLF(fp);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          dtemp = VecSamOutputs_[sInd*nOutputs_];
          if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
            fprintf(fp, "%e %e %e 0\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3]);
          else
            fprintf(fp, "%e %e %e 1\n",
                    VecSamInputs_[sInd*nInputs_+iplot1],
                    VecSamInputs_[sInd*nInputs_+iplot2],
                    VecSamInputs_[sInd*nInputs_+iplot3]);
        }
        fprintf(fp,"];\n");
        fprintf(fp,"iset = find(X(:,4)==1);\n");
        fprintf(fp,"plot3(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100)");
        fprintf(fp,",X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),");
        fprintf(fp,"X(iset,3).*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
        fprintf(fp,"'MarkerSize',13)\n");
        fprintf(fp,"hold on\n");
        fprintf(fp,"iset = find(X(:,4)==0);\n");
        fprintf(fp,"plot3(X(iset,1).*(1+ranflag*rand(size(iset,1),1)/100)");
        fprintf(fp,",X(iset,2).*(1+ranflag*rand(size(iset,1),1)/100),");
        fprintf(fp,"X(iset,3).*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
        fprintf(fp,"'MarkerSize',13)\n");
        fprintf(fp,"hold off\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrInpNames_[iplot3]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "3D Input Scatter Plot");
        fprintf(fp, "if cvFlag == 1\n");
        fprintf(fp, "   figure\n");
        fprintf(fp, "   iset = find(X(:,4)==0);\n");
        fprintf(fp, "   m = size(iset,1);\n");
        fprintf(fp, "   if (m > 3)\n");
        fprintf(fp, "      XX = X(iset,1:3);\n");
        fprintf(fp, "      KK = convhulln(XX);\n");
        fprintf(fp, "      trisurf(KK,XX(:,1),XX(:,2),XX(:,3))\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrInpNames_[iplot3]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp,"3D Input Convex Hull for Good Points");
        fprintf(fp, "   else\n");
        fprintf(fp, "      disp('too few points to display')\n");
        fprintf(fp, "   end;\n");
        fprintf(fp, "end;\n");
        fclose(fp);    
        printf("matlabiplt3.m is now ready for input scatter plots.\n");
        printf("NOTE: see inside the matlab file to see ");
        printf("options to distinguish between\n");
        printf("      valid and invalid runs.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ iplot4m
    //**/ 4D input scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot4m"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot4: plot the sample points in 4-parameter space\n");
        printf("Syntax: iplot4 (no argument needed)\n");
        printf("This command is useful for examining where in the\n");
        printf("parameter space failed runs are. The 4-dimension is\n");
        printf("in time.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a movie of 3D scatter plots for ");
      printf("4 selected inputs\n");
      printf("with the fourth input in the time axis.\n");
      printf("NOTE: if you desire Scilab plot, call scilab before ");
      printf("doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 4)
      {
        printf("ERROR: iplot4m requires 4 or more inputs.\n");
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: iplot4m is currently not available for scilab.\n");
        continue;
      }

      int iplot4 = -1;
      sprintf(pString, 
              "Enter the input for x axis (1 - %d) : ",nInputs_);
      iplot1 = getInt(1, nInputs_, pString);
      iplot1--;
      sprintf(pString, 
              "Enter the input for y axis (1 - %d) : ",nInputs_);
      iplot2 = getInt(1, nInputs_, pString);
      iplot2--;
      if (iplot2 == iplot1)
        printf("WARNING: input %d is used more than once.\n",
               iplot2+1);
      sprintf(pString, 
              "Enter the input for z axis (1 - %d) : ",nInputs_);
      iplot3 = getInt(1, nInputs_, pString);
      iplot3--;
      if (iplot3 == iplot1 || iplot3 == iplot2)
        printf("WARNING: input %d is used more than once.\n",
               iplot3+1);
      sprintf(pString, 
              "Enter the input for t axis (1 - %d) : ",nInputs_);
      iplot4 = getInt(1, nInputs_, pString);
      iplot4--;
      if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
        printf("WARNING: input %d is used more than once.\n",
               iplot4+1);

      //**/ write to matlab file
      fp = fopen("matlabiplt4m.m", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file matlabiplt4m.m.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString," set cvFlag to 1 to use convex hull");
      fwriteComment(fp, pString);
      fprintf(fp, "cvFlag = 0;\n");
      sprintf(pString," change npart to change resolution of the ");
      fwriteComment(fp, pString);
      sprintf(pString," 4th dimension");
      fwriteComment(fp, pString);
      fwritePlotCLF(fp);
      fprintf(fp, "X = [\n");
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        dtemp = VecSamOutputs_[sInd*nOutputs_];
        if (dtemp > 0.5*PSUADE_UNDEFINED || VecSamStates_[sInd] != 1)
          fprintf(fp, "%e %e %e %e 0\n",
                  VecSamInputs_[sInd*nInputs_+iplot1],
                  VecSamInputs_[sInd*nInputs_+iplot2],
                  VecSamInputs_[sInd*nInputs_+iplot3],
                  VecSamInputs_[sInd*nInputs_+iplot4]);
        else
          fprintf(fp, "%e %e %e %e 1\n",
                  VecSamInputs_[sInd*nInputs_+iplot1],
                  VecSamInputs_[sInd*nInputs_+iplot2],
                  VecSamInputs_[sInd*nInputs_+iplot3],
                  VecSamInputs_[sInd*nInputs_+iplot4]);
      }
      fprintf(fp,"];\n");
      fprintf(fp,"iset = find(X(:,5)==0);\n");
      fprintf(fp,"m = size(iset,1);\n");
      fprintf(fp,"if m == 0\n");
      fprintf(fp,"   disp('No valid point.')\n");
      fprintf(fp,"else\n");
      fprintf(fp,"  if m > 100\n");
      fprintf(fp,"     npart = 10;\n");
      fprintf(fp,"  else\n");
      fprintf(fp,"     npart = m / 10;\n");
      fprintf(fp,"  end;\n");
      fprintf(fp,"  if npart == 0\n");
      fprintf(fp,"    npart = 1;\n");
      fprintf(fp,"  end;\n");
      fprintf(fp,"  xmax = max(X(:,4));\n");
      fprintf(fp,"  xmin = min(X(:,4));\n");
      fprintf(fp,"  for ii = 1 : npart\n");
      fprintf(fp,"    iset=find(X(:,4)<xmin+(xmax-xmin)/npart*ii);\n");
      fprintf(fp,"    XT = X(iset,:);\n");
      fprintf(fp,"    iset=find(XT(:,4)>=xmin+(xmax-xmin)/npart*(ii-1));\n");
      fprintf(fp,"    XV = XT(iset,:);\n");
      fprintf(fp,"    if cvFlag == 0\n");
      fprintf(fp,"      iset = find(XV(:,5)==0);\n");
      fprintf(fp,"      plot3(XV(iset,1),XV(iset,2),XV(iset,3),'r*',");
      fprintf(fp,"'MarkerSize',13)\n");
      fprintf(fp,"      hold on\n");
      fprintf(fp,"      iset = find(XV(:,5)~=0);\n");
      fprintf(fp,"      plot3(XV(iset,1),XV(iset,2),XV(iset,3),'bX',");
      fprintf(fp,"'MarkerSize',13)\n");
      fprintf(fp,"      hold off\n");
      fprintf(fp,"    else\n");
      fprintf(fp,"      iset = find(XV(:,5)==1);\n");
      fprintf(fp,"      mm = size(iset,1);\n");
      fprintf(fp,"      if (mm > 3)\n");
      fprintf(fp,"        XW = XV(iset,1:3);\n");
      fprintf(fp,"        KK = convhulln(XW);\n");
      fprintf(fp,"        trisurf(KK,XW(:,1),XW(:,2),XW(:,3))\n");
      fprintf(fp,"      else\n");
      fprintf(fp,"        disp('too few points to display')\n");
      fprintf(fp,"      end;\n");
      fprintf(fp,"    end;\n");
      fwritePlotXLabel(fp, StrInpNames_[iplot1]);
      fwritePlotYLabel(fp, StrInpNames_[iplot2]);
      fwritePlotZLabel(fp, StrInpNames_[iplot3]);
      fwritePlotAxes(fp);
      sprintf(pString,"4D Input Scatter Plot for %s)",
              StrInpNames_[iplot4]);
      fwritePlotTitle(fp, pString);
      fprintf(fp,"     disp('Press ENTER to advance')\n");
      fprintf(fp,"     pause\n");
      fprintf(fp,"   end;\n");
      fprintf(fp,"end;\n");
      fclose(fp);    
      printf("matlabiplt4m.m is now ready for input scatter plots.\n");
      printf("NOTE: see inside the matlab file to see options to \n");
      printf("      distinguish between valid and invalid runs.\n");
        cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oplot2
    //**/ 2D output scatter plot
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oplot2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oplot2: plot the sample outputs in two-parameter space.\n");
        printf("Syntax: oplot2 (no argument needed)\n");
        printf("This command is useful for examining output correlation\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a 2D scatter plot for two ");
      printf("selected outputs on the\n");
      printf("X- and Y-axes.\n");
      printf("NOTE: if you desire Scilab plot, call scilab ");
      printf("before doing this.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ < 2)
      {
        printf("ERROR: oplot2 requires 2 or more outputs.\n");
        cmdStatus = 1;
        continue;
      }

      oplot1 = oplot2 = -1;
      sprintf(pString,
              "Enter the output for x axis (1 - %d) : ",nOutputs_);
      oplot1 = getInt(1, nOutputs_, pString);
      oplot1--;
      sprintf(pString,
              "Enter the output for y axis (1 - %d) : ",nOutputs_);
      oplot2 = getInt(1, nOutputs_, pString);
      oplot2--;
      if (oplot2 == oplot1)
        printf("WARNING: output %d is used more than once.\n",
               oplot2+1);

      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilaboplt2.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilaboplt2.sci.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "Y = [\n");
        ddata = 0.5 * PSUADE_UNDEFINED;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if ((VecSamOutputs_[sInd*nOutputs_+oplot1] > ddata)||
              (VecSamOutputs_[sInd*nOutputs_+oplot2] > ddata))
            fprintf(fp, "%e %e 0\n",
                    VecSamOutputs_[sInd*nOutputs_+oplot1],
                    VecSamOutputs_[sInd*nOutputs_+oplot2]);
          else
            fprintf(fp, "%e %e 1\n",
                    VecSamOutputs_[sInd*nOutputs_+oplot1],
                    VecSamOutputs_[sInd*nOutputs_+oplot2]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "ia1 = find(X(:,4) == 0);\n");
        fprintf(fp, "ia2 = find(X(:,4) == 1);\n");
        fprintf(fp, "drawlater\n");
        fprintf(fp, "if (length(ia1) > 0)\n");
        fprintf(fp, "   plot(Y(ia1,1),Y(ia1,2),'r*')\n");
        fprintf(fp, "   set(gca(),\"auto_clear\",\"off\")\n");
        fprintf(fp, "end;\n");
        fprintf(fp, "if (length(ia2) > 0)\n");
        fprintf(fp, "   plot(Y(ia2,1),Y(ia2,2),'bX')\n");
        fprintf(fp, "end\n");
        fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
        fprintf(fp, "drawnow\n");
        fwritePlotXLabel(fp, StrOutNames_[oplot1]);
        fwritePlotYLabel(fp, StrOutNames_[oplot2]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "2D Output Scatter Plot");
        fclose(fp);    
        printf("scilaboplt2.sci is now ready for scatter plots.\n");
      }
      else
      {
        fp = fopen("matlaboplt2.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlaboplt2.m.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp, "Y = [\n");
        ddata = 0.5 * PSUADE_UNDEFINED;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if ((VecSamOutputs_[sInd*nOutputs_+oplot1] > ddata)||
              (VecSamOutputs_[sInd*nOutputs_+oplot2] > ddata))
            fprintf(fp, "%e %e 0\n",
                    VecSamOutputs_[sInd*nOutputs_+oplot1],
                    VecSamOutputs_[sInd*nOutputs_+oplot2]);
          else
            fprintf(fp, "%e %e 1\n",
                    VecSamOutputs_[sInd*nOutputs_+oplot1],
                    VecSamOutputs_[sInd*nOutputs_+oplot2]);
        }
        fprintf(fp, "];\n");
        fprintf(fp, "iset = find(Y(:,3) == 1);\n");
        fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'bX\')\n");
        fprintf(fp, "hold on\n");
        fprintf(fp, "iset = find(Y(:,3) == 0);\n");
        fprintf(fp, "plot(Y(iset,1),Y(iset,2),\'r*\')\n");
        fprintf(fp, "hold off\n");
        fwritePlotXLabel(fp, StrOutNames_[oplot1]);
        fwritePlotYLabel(fp, StrOutNames_[oplot2]);
        fwritePlotAxes(fp);
        fwritePlotTitle(fp, "2D Output Scatter Plot");
        fclose(fp);    
        printf("matlaboplt2.m is now ready for scatter plots.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ mcmc_set_option
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc_set_option"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc_set_option: set MCMC options\n");
        printf("Syntax: mcmc_set_option <option>\n");
        printf("Current options are: \n");
        printf("   MCMC_gibbs \n");
        printf("   MCMC_brute_force\n");
      }
      if ((!strcmp(winput,"MCMC_gibbs")) || 
          (!strcmp(winput,"MCMC_brute_force"))) 
        psConfig_.putParameter(winput);
      else
      {
        printf("Option not recognized: not action taken.\n");
        printf("Syntax: mcmc_set_option <option>\n");
        printf("Current options are: \n");
        printf("   MCMC_gibbs \n");
        printf("   MCMC_brute_force\n");
      }
    }

    //**/ -------------------------------------------------------------
    // +++ mcmc_dm
    //**/ discrepancy modeling for MCMC 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc_dm"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc_dm: build discrepancy model for use in MCMC\n");
        printf("Syntax: mcmc_dm (no argument needed)\n");
        printf("If you are performing statistical inference ");
        printf("(rsmcmc), there is an\n");
        printf("option in rsmcmc for including discrepancy ");
        printf("modeling. A caveat for\n");
        printf("using discrepancy modeling is that users have ");
        printf("to provide some nominal\n");
        printf("values for the calibration parameters (due to ");
        printf("identifiability). In\n");
        printf("case users do not know a good set of nominal ");
        printf("values, this command\n");
        printf("provides some suggestions. Run this command ");
        printf("without -h to see more\n");
        printf("information.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("MCMC may not work well sometimes due to systematic ");
      printf("bias induced by, \n");
      printf("for example, incomplete physics. To mitigate the ");
      printf("problem, discrepancy\n");
      printf("modeling may be used. The idea is to 'correct' the ");
      printf("simulation model by\n");
      printf("adding an 'adjustment' (or discrepancy) to it. ");
      printf("Notationally, let Y be\n");
      printf("the simulation output, F(X,theta) be the simulation ");
      printf("function, and X\n");
      printf("and theta be the design and calibration parameters, ");
      printf("respectively, then\n");
      printf("    Y = F(X,theta) + d(X)\n");
      printf("where d(X) is called the discrepancy function to ");
      printf("correct F(X,theta).\n");
      printf("So if we have experimental data Y_e(X_e) where X_e ");
      printf("are the design\n");
      printf("settings at which experiments are run, then we can ");
      printf("approximate d(X_e)\n");
      printf("by (say, X_e is an experimental sample of N_e points): \n");
      printf("    d(X_e) ~ Y_e(X_e) - F(X(X_e),theta*)\n");
      printf("So that we now have N_e points ({X_e, d(X_e)}) to ");
      printf("build a discrepancy\n");
      printf("response surface.\n");
      printf("Notice that d(X_e) depends only on the design ");
      printf("inputs X. Thus, to build\n");
      printf("a response surface d(X) from the N_e points, ");
      printf("we have to\n");
      printf("1. Fix theta at some values theta* (to ensure ");
      printf("uniqueness of d(X))\n");
      printf("2. Have sufficient experimental data (e.g. ");
      printf("size ~ 10 * size(X))\n");
      printf("3. Ensure that the experimental input covers ");
      printf("the X space well.\n");
      printf("One question is: what should theta be fixed at? ");
      printf("Ideally, theta should\n");
      printf("be set to their most likely values based on the ");
      printf("modeler's judgment.\n");
      printf("When expert judgment is not available, the ");
      printf("alternative is to search\n");
      printf("in the theta space for a good discrepancy model. ");
      printf("There may be many\n");
      printf("ways to create a good discrepancy model. Two ");
      printf("criteria are given in\n");
      printf("this command - to find calibration parameter ");
      printf("settings that minimize: \n");
      printf("1. cross validation error from the set {X_e, ");
      printf("d(X_e)} given a theta*,\n");
      printf("2. the area under d(X) within the X space.\n\n");

      printf("This command facilitates (A) the evaluation of a ");
      printf("user-specified theta*\n");
      printf("or (B) searching the theta space for a good theta* ");
      printf("model. It needs :\n");
      printf("1. an evaluated sample F(X,theta) (which should ");
      printf("have been loaded)\n");
      printf("2. a data file for {X_e,Y_e} (in rsmcmc experimental ");
      printf("file format)\n");
      printf("3. a response surface choice.\n");

      printf("This command then performs the following :\n");
      printf("a. build a response surface for the loaded ");
      printf("sample (for F(X,theta))\n");
      printf("b. draw N points from the calibration ");
      printf("parameter space (for option (B))\n");
      printf("c. for each sample  theta_k, k=1 ... N ");
      printf("(N=1 for option (B))\n");
      printf("   - create a sample input set comprising theta_k ");
      printf("and X_e ==> S_k\n");
      printf("     (X_e are design inputs of all experiments)\n");
      printf("   - evaluate S_k with RS ==> Y_k\n");
      printf("   - compute D_k = Y_e - Y_k\n");
      printf("     (D_k is the discrepancy sample for instance C_k)\n");
      printf("   - perform a left-one-out CV on D_k\n");
      printf("d. select the theta_k that gives the smallest ");
      printf("selected metric\n");
      printf("NOTE: This command is for cases when nOutputs = 1. ");
      printf("It does not make\n");
      printf("      sense for multiple outputs because each output ");
      printf("has a different\n");
      printf("      best parameter setting theta*, but in MCMC a ");
      printf("single theta* is\n");
      printf("      required for all outputs.\n"); 

      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ != 1)
      {
        printf("ERROR: nOutputs should be equal to 1.\n");
        cmdStatus = 1;
        continue;
      }
      if (checkResidentSampleOutputs()) 
      {
        cmdStatus = 1;
        continue;
      }

      //**/ ----------------------------------------------------
      //**/ first of all, get experimental data
      //**/ use MCMC analyzer's read experiment function
      //**/ check to make sure there is no repeated sample points
      //**/ ----------------------------------------------------
      MCMCAnalyzer *mcmc = new MCMCAnalyzer();
      psIVector vecDesignP;
      psMatrix  matExpInps, matExpMeans, matExpStdvs;
      int combFlag=0;
      double dstat = mcmc->readSpecFile(nInputs_, iOne, vecDesignP, 
                           matExpInps,matExpMeans,matExpStdvs,combFlag,0);
      delete mcmc;
      if (dstat != 0)
      {
        printf("ERROR: reading experimental file unsuccessful.\n");
        cmdStatus = 1;
        continue;
      }
      if (matExpMeans.ncols() != nOutputs_)
      {
        printf("ERROR: Mismatch in the number of outputs.\n");
        printf("       Experimental file has %d outputs\n",
               matExpMeans.ncols());
        printf("       The loaded sample has %d outputs\n",
               nOutputs_);
        printf("*** They have to be the same.\n");
        cmdStatus = 1;
        continue;
      }
      int dnSamples = matExpMeans.nrows();
      int dnInputs  = matExpInps.ncols();
      int nCalib    = nInputs_ - dnInputs;
      psVector vecExpInps;
      vecExpInps.setLength(dnSamples*dnInputs);
      for (ii = 0; ii < dnSamples; ii++)
        for (jj = 0; jj < dnInputs; jj++)
          vecExpInps[ii*dnInputs+jj] = matExpInps.getEntry(ii,jj); 
      if (checkRepeatedSamplePts(dnSamples,dnInputs,
                                 vecExpInps.getDVector()))
      {
        printf("ERROR: the experimental sample have ");
        printf("duplicate sample inputs. Abort.\n");
        continue;
      }

      //**/ ----------------------------------------------------
      //**/ ask for options
      //**/ ----------------------------------------------------
      printf("Do you want to : \n");
      printf("1. test your own calibration parameter values, or\n");
      printf("2. have PSUADE to search for good ");
      printf("calibration parameter values?\n");
      sprintf(pString, "Select 1 or 2 : ");
      int srchOption = getInt(1,2,pString);
 
      int dmOption = 1;
      if (srchOption == 2)
      {
        printf("There are 2 options to choose discrepancy model : \n");
        printf("1. The one that minimizes one-at-a-time CV error.\n");
        printf("2. The one that minimizes the energy of the RS.\n");
        sprintf(pString, "Select 1 or 2 : ");
        dmOption = getInt(1,2,pString);
      }

      //**/ ----------------------------------------------------
      //**/ create response surface for the original loaded sample 
      //**/ ==> faPtr
      //**/ ----------------------------------------------------
      int rstype = -1, iZero=0;
      FuncApprox *faPtr = genFA(rstype, nInputs_, iZero, nSamples_);
      faPtr->setNPtsPerDim(16);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(0);
      psVector vecSamOut1;
      vecSamOut1.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
        vecSamOut1[ii] = VecSamOutputs_[ii];
      status = faPtr->initialize(VecSamInputs_.getDVector(), 
                                 vecSamOut1.getDVector());
      if (status != 0)
      {
        printOutTS(PL_ERROR,
           "mcmc_dm ERROR: Unable to create response surface.\n");
        delete faPtr;
        continue;
      }

      //**/ ----------------------------------------------------
      //**/ create a sample for the calibration parameters
      //**/ that is, search out of multiple settings of calibration
      //**/ parameters to see which best fit discrepancy model
      //**/ (since we don't know where we should place the 
      //**/ nominal point) ==> vecTestCInps
      //**/ ----------------------------------------------------
      int nTests;
      psVector vecCLBs, vecCUBs, vecTestCInps, vecTestCOuts;
      Sampling *sampPtr = NULL;
      vecCLBs.setLength(nCalib);
      vecCUBs.setLength(nCalib);
      kk = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecDesignP[ii] == 0)
        {
          vecCLBs[kk] = VecILowerBs_[ii];
          vecCUBs[kk] = VecIUpperBs_[ii];
          kk++;
        }
      }
      if (srchOption == 1)
      {
        vecTestCInps.setLength(nCalib);
        for (ii = 0; ii < nCalib; ii++)
        {
          sprintf(pString, 
             "Enter calibration parameter %d value (in [%12.4e, %12.4e]) : ",
             ii+1,vecCLBs[ii],vecCUBs[ii]);
          vecTestCInps[ii] = getDouble(pString);
        } 
        nTests = 1;
      }
      else
      {
        printf("Choose how many samples in the uncertain ");
        printf("parameter space to explore to\n");
        printf("find the best discrepancy model. Remember: ");
        printf("the more the samples, the\n");
        printf("longer it takes for this command.\n");
        strcpy(pString, "How many samples to search ? (10-1000) ");
        nTests = getInt(10,1000,pString);
        sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
        sampPtr->setPrintLevel(0);
        sampPtr->setInputBounds(nCalib, vecCLBs.getDVector(),
                                vecCUBs.getDVector());
        sampPtr->setSamplingParams(nTests,-1,0);
        sampPtr->initialize(0);
        psVector vecTestCOuts;
        psIVector vecTestCStas;
        vecTestCInps.setLength(nTests*nCalib);
        vecTestCOuts.setLength(nTests);
        vecTestCStas.setLength(nTests);
        sampPtr->getSamples(nTests,nCalib,iZero,vecTestCInps.getDVector(),
                   vecTestCOuts.getDVector(), vecTestCStas.getIVector());
        delete sampPtr;
        sampPtr = NULL;
      }

      //**/ ----------------------------------------------------
      //**/ ask user for discrepancy response surface type 
      //**/ ==> rstype
      //**/ ----------------------------------------------------
      rstype = -1;
      FuncApprox *faPtr2 = genFA(rstype,dnInputs,iZero,dnSamples-1);
      rstype = faPtr2->getID();
      delete faPtr2;

      //**/ ----------------------------------------------------
      //**/ create multiple response surfaces (for OpenMP)
      //**/ ----------------------------------------------------
      ///**/ first allocate variables such as lower/upper bounds
      FuncApprox **faPtrs = new FuncApprox*[nTests];
      psVector vecDLBs, vecDUBs;
      vecDLBs.setLength(dnInputs);
      vecDUBs.setLength(dnInputs);

      kk = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecDesignP[ii] == 1)
        {
          vecDLBs[kk] = VecILowerBs_[ii];
          vecDUBs[kk] = VecIUpperBs_[ii];
          kk++;
        }
      }

      //**/ ----------------------------------------------------
      //**/ for option 2, generate a design for energy 
      //**/ ----------------------------------------------------
      int nSamEnergy = 2000;
      psVector vecEnergySamInps, vecEnergySamOuts;
      if (dmOption == 2)
      {
        sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
        sampPtr->setPrintLevel(0);
        sampPtr->setInputBounds(dnInputs, vecDLBs.getDVector(),
                                vecDUBs.getDVector());
        sampPtr->setSamplingParams(nSamEnergy,-1,0);
        sampPtr->initialize(0);
        vecEnergySamInps.setLength(nSamEnergy*dnInputs);
        vecEnergySamOuts.setLength(nSamEnergy);
        psIVector vecEnergySamStas;
        vecEnergySamStas.setLength(nSamEnergy);
        sampPtr->getSamples(nSamEnergy,dnInputs,iZero,
                 vecEnergySamInps.getDVector(),
                 vecEnergySamOuts.getDVector(),
                 vecEnergySamStas.getIVector());
        delete sampPtr;
        sampPtr = NULL;
      }

      //**/ ----------------------------------------------------
      //**/ set up response surfaces
      //**/ ----------------------------------------------------
      int nn, count;
      psConfig_.InteractiveSaveAndReset();
      psConfig_.RSExpertModeSaveAndReset();
      for (nn = 0; nn < nTests; nn++)
      {
        faPtrs[nn] = genFA(rstype,dnInputs,iOne,dnSamples-1);
        if (faPtrs[nn] == NULL)
        {
          printOutTS(PL_INFO,
               "mcmc_dm: genFA returned NULL in file %s line %d\n",
               __FILE__, __LINE__ );
          exit(1);
        }
        faPtrs[nn]->setNPtsPerDim(16);
        faPtrs[nn]->setBounds(vecDLBs.getDVector(),vecDUBs.getDVector());
        faPtrs[nn]->setOutputLevel(0);
      }

      //**/ set up response surface for discrepancy model for option 2
      //**/ This RS has size dnSamples instead of dnSamples for CV
      faPtr2 = NULL;
      if (dmOption == 2)
      {
        faPtr2 = genFA(rstype,dnInputs,iOne,dnSamples);
        if (faPtr2 == NULL)
        {
          printOutTS(PL_INFO,
               "mcmc_dm: genFA returned NULL in file %s line %d\n",
               __FILE__, __LINE__ );
          exit(1);
        }
        faPtr2->setNPtsPerDim(16);
        faPtr2->setBounds(vecDLBs.getDVector(),vecDUBs.getDVector());
        faPtr2->setOutputLevel(0);
      }

      //**/ allocate vectors for testing and run
      double errNorm, *expInps = matExpInps.getMatrix1D();
      psVector vecTestAllInps, vecTestAllOuts, vecDiscOuts;
      psVector vecCVInps, vecCVOuts, vecAllNorms, **vecAllDiscOuts;
      vecAllNorms.setLength(nTests);
      vecAllDiscOuts = new psVector*[nTests];
      for (nn = 0; nn < nTests; nn++)
        vecAllDiscOuts[nn] = new psVector();
#pragma omp parallel shared(vecDesignP,vecTestCInps,dnSamples,nCalib,matExpInps,faPtr,faPtr2,faPtrs,dnInputs,expInps,vecAllNorms,vecAllDiscOuts,rstype,vecEnergySamInps,vecEnergySamOuts) \
    private(nn,ii,jj,kk,vecTestAllInps,vecTestAllOuts,vecCVInps,vecCVOuts,vecDiscOuts,status,count,errNorm,ddata)
#pragma omp for
      for (nn = 0; nn < nTests; nn++)
      {
        printf("Testing uncertain input sample %d out of %d\n",
               nn+1,nTests);
        vecTestAllInps.setLength(dnSamples*nInputs_);
        vecTestAllOuts.setLength(dnSamples);

        //**/ fill in vecTestAllInps
        for (kk = 0; kk < dnSamples; kk++)
        {
          //**/ fill in the calibration test values of 
          //**/ vecTestAllInps
          jj = 0;
          for (ii = 0; ii < nInputs_; ii++)
          {
            if (vecDesignP[ii] == 0) 
            {
              vecTestAllInps[kk*nInputs_+ii] = 
                    vecTestCInps[nn*nCalib+jj]; 
              jj++;
            }
          }
          //**/ fill in the design values of vecTestAllInps 
          //**/ for all experimental samples
          jj = 0;
          for (ii = 0; ii < nInputs_; ii++)
          {
            if (vecDesignP[ii] == 1) 
            {
              vecTestAllInps[kk*nInputs_+ii] = 
                    matExpInps.getEntry(kk,jj);
              jj++;
            }
          }
        }
 
        //**/ evaluate the loaded sample response surface 
        //**/ on the experimental inputs
        status = faPtr->evaluatePoint(dnSamples,
                          vecTestAllInps.getDVector(),
                          vecTestAllOuts.getDVector());

        //**/ now compute discrepancy = expdata - simdata
        vecDiscOuts.setLength(dnSamples);
        for (kk = 0; kk < dnSamples; kk++)
           vecDiscOuts[kk] = matExpMeans.getEntry(kk,0) -
                             vecTestAllOuts[kk];

        //**/ now we have vecDiscOuts
        //**/ vecDiscInps is just matExpInps
        //**/ next we perform take-one-out CV
        vecCVInps.setLength(dnSamples*dnInputs);
        vecCVOuts.setLength(dnSamples);
        if (dmOption == 1)
        {
          errNorm = 0;  
          for (kk = 0; kk < dnSamples; kk++)
          {
            //**/ fill the vecCVInps
            count = 0;
            for (jj = 0; jj < dnSamples; jj++)
            {
              if (jj != kk)
              {
                for (ii = 0; ii < dnInputs; ii++)
                  vecCVInps[count*dnInputs+ii] =
                     matExpInps.getEntry(jj,ii);
                vecCVOuts[count] = vecDiscOuts[jj];
                count++;
              }
            }
            if (rstype != PSUADE_RS_MARS   && rstype != PSUADE_RS_MMARS &&
                rstype != PSUADE_RS_MARSB  && rstype != PSUADE_RS_GP1   && 
                rstype != PSUADE_RS_KR     && rstype != PSUADE_RS_SVM &&
                rstype != PSUADE_RS_ACOSSO && rstype != PSUADE_RS_BSSANOVA &&
                rstype != PSUADE_RS_HKR    && rstype != PSUADE_RS_GP3)
              faPtrs[nn]->initialize(vecCVInps.getDVector(), 
                                     vecCVOuts.getDVector());
            else
            {
#pragma omp critical
              faPtrs[nn]->initialize(vecCVInps.getDVector(), 
                                     vecCVOuts.getDVector());
            }
            ddata = faPtrs[nn]->evaluatePoint(&(expInps[kk*dnInputs]));
            errNorm += pow(vecDiscOuts[kk] - ddata, 2.0);
          }
          errNorm = sqrt(errNorm/dnSamples);
          printf("Current test = %4d, current error norm = %e\n",nn+1,
                 errNorm);
          vecAllNorms[nn] = errNorm;
        }
        else
        //**/ compute energy
        {
          //**/ copy design inputs 
          for (jj = 0; jj < dnSamples; jj++)
          {
            for (ii = 0; ii < dnInputs; ii++)
              vecCVInps[jj*dnInputs+ii] = matExpInps.getEntry(jj,ii);
            vecCVOuts[count] = vecDiscOuts[jj];
          }
#pragma omp critical
          faPtr2->initialize(vecCVInps.getDVector(), 
                             vecCVOuts.getDVector());
          faPtr2->evaluatePoint(nSamEnergy,vecEnergySamInps.getDVector(),
                                vecEnergySamOuts.getDVector());
          errNorm = 0;
          for (jj = 0; jj < nSamEnergy; jj++)
            errNorm += pow(vecEnergySamOuts[jj], 2.0) / nSamEnergy;
          errNorm = sqrt(errNorm);
        }
        vecAllDiscOuts[nn]->load(dnSamples,vecDiscOuts.getDVector());
      }
      psConfig_.InteractiveRestore();
      psConfig_.RSExpertModeRestore();
      for (nn = 0; nn < nTests; nn++) delete faPtrs[nn];
      if (faPtrs != NULL) delete [] faPtrs;
      if (faPtr  != NULL) delete faPtr;
      if (faPtr2 != NULL) delete faPtr2;

      //**/ search for the best one
      int    bestInd = -1;
      double bestNorm = PSUADE_UNDEFINED;
      psVector vecBestDiscOuts;
      for (nn = 0; nn < nTests; nn++)
      {
        if (vecAllNorms[nn] < bestNorm)
        {
          bestNorm = vecAllNorms[nn];
          bestInd = nn;
          vecBestDiscOuts.load(dnSamples,
                             vecAllDiscOuts[nn]->getDVector());
        }
        delete vecAllDiscOuts[nn];
        vecAllDiscOuts[nn] = NULL;
      }
      delete vecAllDiscOuts;
      if (srchOption == 2)
      {
        printf("The best calibration parameter values are ");
        printf("(test #%d, CVerr=%e):\n",bestInd+1,vecAllNorms[bestInd]);
        for (ii = 0; ii < nCalib; ii++)
          printf("Param %3d = %e\n",ii+1,vecTestCInps[bestInd*nCalib+ii]); 
        printf("You can use these settings in 'rsmcmc' as design ");
        printf("default values when\n");
        printf("discrepancy is turned on.\n");
      }
      fp = fopen("mcmc_dm_sample.std", "w");
      if (fp == NULL)
      {
        printf("mcmc_dm ERROR: unable to write to a ");
        printf("discrepancy sample file.\n");
        continue;
      }
      fprintf(fp,"%d %d 1\n",dnSamples,dnInputs);
      for (ii = 0; ii < vecBestDiscOuts.length(); ii++)
      {
        for (jj = 0; jj < dnInputs; jj++)
          fprintf(fp, "%24.16e ", matExpInps.getEntry(ii,jj));
        fprintf(fp, "%24.16e\n", vecBestDiscOuts[ii]);
      }
      fprintf(fp,"Best Calibration Parameter Values\n");
      for (ii = 0; ii < nCalib; ii++)
        fprintf(fp,"Param %3d = %e\n",ii+1,
                vecTestCInps[bestInd*nCalib+ii]); 
      fclose(fp);
      printf("INFO: The discrepancy sample based on the ");
      printf("best theta is now available\n");
      printf("      in mcmc_dm_sample.std. You will need ");
      printf("to do the following to see\n");
      printf("      if discrepancy modeling should be used ");
      printf("in your rsmcmc run: \n");
      printf("1. perform CV on this sample with your ");
      printf("selected RS method\n");
      printf("2. evaluate your loaded sample with the RS ");
      printf("created from this sample to\n");
      printf("   see if the corrections (d(X)) is acceptable.\n");
      printf("Discrepancy modeling in rsmcmc is recommended only ");
      printf("if the following\n");
      printf("conditions are met:\n");
      printf("a. the number of experimental data is at ");
      printf("least 10xn where n is the\n");
      printf("   number of design variables,\n");
      printf("b. the discrepancy RS does not look like a ");
      printf("constant function, \n");
      printf("c. CV using GP or Kriging shows good results ");
      printf("(small deviation in\n");
      printf("   parity plot), and\n");
      printf("d. the loaded sample (F(X,theta)) evaluated ");
      printf("with the discrepancy\n");
      printf("   response surface shows acceptable output values.\n"); 

      printf("Now you can choose to construct a response ");
      printf("surface from [X_e, d(X_e)]\n");
      printf("(the discrepancy sample) evaluated at the ");
      printf("best calibration parameter\n");
      printf("values given above and evaluate the ");
      printf("loaded sample with it. The result\n");
      printf("will be written a file for you to inspect.\n");
      sprintf(pString, 
              "Perform response surface evaluation ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
        //**/ build RS for discrepancy model
        rstype = -1;
        faPtr = genFA(rstype,dnInputs,iOne,dnSamples);
        if (faPtr == NULL)
        {
          printOutTS(PL_INFO,
               "mcmc_dm: genFA returned NULL in file %s line %d\n",
               __FILE__, __LINE__ );
          exit(1);
        }
        faPtr->setNPtsPerDim(16);
        faPtr->setBounds(vecDLBs.getDVector(),vecDUBs.getDVector());
        faPtr->setOutputLevel(0);
        psVector vecTInps, vecTOuts;
        vecTInps.setLength(dnSamples*dnInputs);
        vecTOuts.setLength(nSamples_);
        count = 0;
        for (kk = 0; kk < dnSamples; kk++)
        {
          for (ii = 0; ii < dnInputs; ii++)
            vecTInps[kk*dnInputs+ii] = matExpInps.getEntry(kk,ii);
        }
        faPtr->initialize(vecTInps.getDVector(), 
                          vecBestDiscOuts.getDVector());
        vecTInps.setLength(nSamples_*dnInputs);
        for (kk = 0; kk < nSamples_; kk++)
        {
          //**/ fill in the design parameter values
          jj = 0;
          for (ii = 0; ii < nInputs_; ii++)
          {
            if (vecDesignP[ii] == 1) 
            {
              vecTInps[kk*dnInputs+jj] = 
                       VecSamInputs_[kk*nInputs_+ii];
              jj++;
            }
          }
        }
        status = faPtr->evaluatePoint(nSamples_,
                  vecTInps.getDVector(),vecTOuts.getDVector());
        fp = fopen("mcmc_dm_predict.std", "w");
        if (fp == NULL)
        {
          printf("mcmc_dm ERROR: failed to crete ");
          printf("mcmc_dm_predict.std.\n");
          cmdStatus = 1;
          continue;
        }
        fprintf(fp,"%d %d 1\n", nSamples_, dnInputs);
        for (kk = 0; kk < nSamples_; kk++)
        {
          for (ii = 0; ii < dnInputs; ii++)
            fprintf(fp,"%16.8e ", vecTInps[kk*dnInputs+ii]); 
          fprintf(fp,"%16.8e\n", vecTOuts[kk]); 
        }
        fclose(fp);
        printf("INFO: mcmc_dm_predict.std now has the ");
        printf("discrepancy values for each of\n");
        printf("      the %d sample points in the loaded ",nSamples_);
        printf("sample. Now you can\n");
        printf("      examine this file to determine whether the ");
        printf("discrepancy data\n");
        printf("      (d(X_s) where X_s is the loaded sample ");
        printf("inputs) makes sense.\n");
        printf("      (Note that the d(X) function is obtained ");
        printf("from a response surface\n");
        printf("      built from Y_e(X_e)-F(X_e,theta*) where ");
        printf("theta* is shown above.)\n");
        printf("      If so, you may use the best calibration ");
        printf("parameter values above\n");
        printf("      in rsmcmc with discrepancy modeling turned on.\n");
      }
    }

    //**/ -------------------------------------------------------------
    // +++ ocov
    //**/ create a covariance matrix 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ocov"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ocov: build output covariance matrix\n");
        printf("Syntax: ocov (no argument needed)\n");
        printf("This command creates a covariance matrix ");
        printf("of all the outputs in the\n"); 
        printf("loaded sample.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command takes the outputs of the loaded ");
      printf("sample and construct a\n");
      printf("covariance matrix. There are two options: \n");
      printf("1. Use the raw sample\n");
      printf("2. Use a sample drawn from the response surfaces\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      printf("Under-development\n");
    }

    //**/ -------------------------------------------------------------
    // +++ rsmcmc
    //**/ response surface based MCMC
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsmcmc"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsmcmc: perform response surface-based inference \n");
        printf("Syntax: rsmcmc (no argument needed)\n");
        printf("This command performs statistical inference ");
        printf("using response surfaces\n");
        printf("constructed from the sample outputs.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs statistical inference ");
      printf("using response surface\n");
      printf("models from the sample outputs.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      printf("You may choose between the 3 MCMC schemes: \n");
      printf("1. Gibbs\n");
      printf("2. Metropolis-Hastings (MH)\n");
      printf("3. Brute-force (BF)\n");
      printf("A few notes about these methods: \n");
      printf("- Gibbs and MH are popular and generally ");
      printf("fast methods.\n");
      printf("- Gibbs and MH may have problems when ");
      printf("posteriors are multimodal.\n");
      printf("- BF method may be slow if large sample is used.\n");
      printf("- BF method may not be effective if the random ");
      printf("sample drawn from the\n");
      printf("  prior is not large enough to capture the shape of ");
      printf("the posterior\n");
      printf("  distribution (e.g. if the posterior distribution ");
      printf("is a narrow peak).\n");
      printf("- Therefore, BF is not recommended for ");
      printf("high-dimensional calibration.\n");
      printf("- BF may be good to capture multimodal posteriors.\n");

      //printf("2. Metropolis-Hastings\n");
      sprintf(pString, "Enter MCMC scheme (1, 2, or 3) : ");
      int mcmcScheme = getInt(1,3,pString);
      if      (mcmcScheme == 1) sprintf(pString, "MCMC_gibbs");
      else if (mcmcScheme == 2) sprintf(pString, "MCMC_mh");
      else if (mcmcScheme == 3) sprintf(pString, "MCMC_brute_force");
      psConfig_.putParameter(pString);
      int analysisMethod = PSUADE_ANA_MCMC;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 1, NULL, -1);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      psConfig_.removeParameter(pString);
      printf("INFO: Run the matlab/scilab-mcmc2 file to ");
      printf("visualize posteriors.\n");
      printf("INFO: The posterior sample can now be found in ");
      printf("'MCMCPostSample', along\n");
      printf("      with useful likelihood information:\n");
      printf("      * Information that tells you which experiment ");
      printf("may be off.\n");
      printf("      * Information that indicates whether discrepancy ");
      printf("is needed.\n");
      printf("INFO: If you have selected to use discrepancy, ");
      printf("you should see a sample\n");
      printf("      file called 'psDiscrepancyModel' (You can ");
      printf("perform a response\n");
      printf("      surface analysis on it).\n");
      printf("INFO: MCMCPostSample/psDiscrepancyModel can be used ");
      printf("for prediction\n");
      printf("      with the 'mcmc_predict' command.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ hmcmc
    //**/ hierarchical MCMC
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "hmcmc*"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("hmcmc: perform hierarchical MCMC\n");
        printf("Syntax: hmcmc (no argument needed)\n");
        printf("This command is used for single mixed effect analysis.\n");
        printf("Specifically, a family of individual members which are\n");
        printf("characterized by their means and standard deviations are\n");
        printf("divided into sub-families with each sub-family having\n");
        printf("its own mean and standard deviation.\n"); 
        continue;
      }
      HMCMCAnalyzer *hmcmc = new HMCMCAnalyzer();
      hmcmc->analyze();
      delete hmcmc;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ mcmc
    //**/ simulation based MCMC (not done yet)
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc: perform MCMC on actual simulator\n");
        printf("Syntax: mcmc (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs statistical ");
      printf("inference using actual simulator.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }
      //if (nOutputs_ != 1)
      //{
      //  printf("ERROR: mcmc does not support nOutputs > 1\n");
      //  cmdStatus = 1;
      //  continue;
      //}
      int analysisMethod = PSUADE_ANA_MCMC;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      strcpy(winput, "setsim");
      targv[0] = winput;
      anaManager->specialRequest(analysisMethod, 1, targv);
      //**/if (nOutputs_ > 1)
      //**/{
      //**/   sprintf(pString, "Which output to analyze (1 - %d) ? ",
      //**/           nOutputs_);
      //**/   outputID = getInt(1, nOutputs_, pString);
      //**/}
      //**/else outputID = 1;
      //**/outputID--;
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      int diagKeep = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      outputID = 0;
      anaManager->analyze(psuadeIO_, 1, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,diagKeep,-1,-1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ create probability plots for the posterior distributions
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc_genplot"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc_genplot: generate distribution plots for ");
        printf("posterior sample.\n");
        printf("Syntax: mcmc_genplot (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command takes the loaded sample as a MCMC ");
      printf("posterior sample and\n");
      printf("generates the 1-parameter and 2-parameter ");
      printf("probability distribution\n");
      printf("plots similar to the ones created by rsmcmc.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data has not been loaded.\n");
        cmdStatus = 1;
        continue;
      }

      McmcData mdata;
      mdata.printLevel_ = 0;
      mdata.nSamples_ = nSamples_;
      mdata.nInputs_ = nInputs_;
      mdata.nOutputs_ = nOutputs_;
      mdata.VecLowerB_ = VecILowerBs_;
      mdata.VecUpperB_ = VecIUpperBs_;
      mdata.VecSamInputs_ = VecSamInputs_;
      mdata.VecSamOutputs_ = VecSamOutputs_;
      double dOne=1.0, dZero=0;
      mdata.MatExpMeans_.setDim(1, nOutputs_);
      mdata.MatExpStds_.setDim(1, nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        mdata.MatExpMeans_.setEntry(0, ii, dZero);
        mdata.MatExpStds_.setEntry(0, ii, dOne);
      }
      mdata.StrInpNames_.load(nInputs_, 
            (const char **) StrInpNames_.getStrings());
      MCMCAnalyzer *mcmc = new MCMCAnalyzer();
      mcmc->analyzeDirect_nors(mdata);
      delete mcmc;
      cmdStatus = 0;
    }

    //**/ =============================================================
    // miscellaneous commands
    //**/ =============================================================
    //**/ -------------------------------------------------------------
    // +++ urefine
    //**/ sample refinement
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "urefine"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("urefine: refine the loaded sample uniformly\n");
        printf("Syntax: urefine (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command refines the loaded sample ");
      printf("uniformly based on the type of\n");
      printf("sampling method in the loaded sample. For ");
      printf("example, if the sampling\n");
      printf("method of the loaded sample is Latin hypercube, ");
      printf("then the refinement\n");
      printf("strategy for Latin hypercube will be applied.\n");
      printf("The refined sample should be approximately ");
      printf("twice the size of the\n");
      printf("original one (except OA).\n");
      printf("NOTE: This command operates only on a few sampling "); 
      printf("designs such as\n");
      printf("      Latin hypercube, Monte Carlo, LPTAU, ");
      printf("METIS, and orthogonal array\n");
      printf("NOTE: If the original sampling method is METIS, ");
      printf("'urefine' will require\n");
      printf("      that the psuadeMetisInfo file created ");
      printf("previously.\n"); 
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      //**/ instantiate the corresponding sampler to prepare for
      //**/ refinement
      psuadeIO_->getParameter("method_sampling", pPtr);
      SamplingMethod_ = pPtr.intData_;
      if (SamplingMethod_ != PSUADE_SAMP_MC &&
          SamplingMethod_ != PSUADE_SAMP_LPTAU &&
          SamplingMethod_ != PSUADE_SAMP_METIS &&
          SamplingMethod_ != PSUADE_SAMP_LHS &&
          SamplingMethod_ != PSUADE_SAMP_OA)
      {
        printf("ERROR: sampling method not supported.\n");
        printf("NOTE: This command operates only on a few sampling "); 
        printf("designs such as\n");
        printf("      Latin hypercube, Monte Carlo, LPTAU, ");
        printf("METIS, and orthogonal array\n");
        cmdStatus = 1;
        continue;
      } 
      psConfig_.SamExpertModeSaveAndReset();
      Sampling *sampPtr = 
         sampPtr = (Sampling *) SamplingCreateFromID(SamplingMethod_);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
      sampPtr->setOutputParams(nOutputs_);
      psuadeIO_->getParameter("method_nreplications",pPtr);
      nReplications_ = pPtr.intData_;
      sampPtr->setSamplingParams(nSamples_, nReplications_, -1);
      sampPtr->initialize(1);

      //**/ load the original sample to the sampler and refine
      sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, 
                           VecSamInputs_.getDVector(),
                           VecSamOutputs_.getDVector(), 
                           VecSamStates_.getIVector());
      status = sampPtr->refine(2, 0, 0.0, nSamples_, NULL);

      //**/ store results to psuadeIO_ 
      if (status == 0)
      {
        VecSamInputs_.clean();
        VecSamOutputs_.clean();
        VecSamStates_.clean();
        nSamples_ = sampPtr->getNumSamples();
        VecSamInputs_.setLength(nInputs_*nSamples_);
        VecSamOutputs_.setLength(nOutputs_*nSamples_);
        VecSamStates_.setLength(nSamples_);
        sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
              VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(),
              VecSamStates_.getIVector());
        psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,
              NULL,VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,
              NULL); 
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
              VecSamOutputs_.getDVector(),VecSamStates_.getIVector(),
              NULL); 
        psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
        printf("urefine successful: use write to store ");
        printf("the refined sample.\n");
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
      }
      else
      {
        printf("ERROR: urefine not successful.\n");
        cmdStatus = 1;
      }
      delete sampPtr;
      sampPtr = NULL;
      psConfig_.SamExpertModeRestore();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arefine_metis
    //**/ adaptive sample refinement
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arefine_metis"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arefine_metis: adaptively refine the loaded sample ");
        printf("based on prediction\n");
        printf("        uncertainty (std dev) at the new sample ");
        printf("points. This function\n");
        printf("        assumes the loaded sample is a METIS ");
        printf("sample (along with the\n");
        printf("        corresponding METIS information file ");
        printf("(psuadeMetisInfo).\n");
        printf("Syntax: arefine_metis (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("Given a sample which has been loaded (using ");
      printf("load or read_std), this\n");
      printf("command adaptively refine the sample near ");
      printf("either areas in the input\n");
      printf("space with the largest prediction uncertainty. ");
      printf("If the selected output\n");
      printf("is either 0 or 1, refinement will be focused ");
      printf("near the 0/1 boundary.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,
          "Select output to use for prediction estimation (1 - %d) ",
          nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      psuadeIO_->getParameter("method_sampling", pPtr);
      SamplingMethod_ = pPtr.intData_;
      if (SamplingMethod_ != PSUADE_SAMP_METIS)
      {
        printf("ERROR: adaptive refinement requires METIS sample.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ check whether sample outputs are 0/1. If so, set flag = 0
      //**/ since it uses a different method (to add points along the
      //**/ 0/1 boundary) so this one is similar to arefine
      //**/ to add new points in focused regions
      int flag = 1;
      int refineSize=0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] != 1 || 
            VecSamOutputs_[ss*nOutputs_+outputID] != 0)
        {
          flag = 0;
          break;
        }
      }
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];

      //**/ case when the loaded sample outputs are not 0/1
      if (flag == 0)
      {
        printf("This function uses the current loaded sample ");
        printf("together with previous\n");
        printf("refinement information (stored in file ");
        printf("'psuadeMetisInfo') to perform\n");
        printf("adaptive refinement. To work properly, the ");
        printf("following information are\n");
        printf("needed: \n");
        printf("(1) the original sample size (before any refinement)\n");
        printf("    NOTE: if you do not know what this mean, enter %d\n",
               nSamples_);
        printf("(2) the previous refinement file (psuadeMetisInfo)\n");
        printf("    which needs to be in the current directory.\n");
        sprintf(pString,"What is the original sample size ? ");
        int origNSamples = getInt(1, nSamples_, pString);
        sprintf(pString,
              "How many sample points to add? (1 - %d, default = %d) ",
              nSamples_, nSamples_);
        refineSize = getInt(1, nSamples_, pString);

        //**/ temporarily turn off sampling expert mode
        //**/ load sample to 2 samplers, 1 for uniform sampling (sampAux),
        //**/ and the other for adaptive sampling (sampPtr)
        //**/ The reason for the uniform sample is that it needs to give
        //**/ the prediction errors for all POSSIBLE refined points (in
        //**/ the uniformly refined mesh) for METIS refine function to
        //**/ select candidates
        psConfig_.SamExpertModeSaveAndReset();
        Sampling *sampPtr = 
            (Sampling *) SamplingCreateFromID(SamplingMethod_);
        sampPtr->setPrintLevel(outputLevel_);
        sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                                VecIUpperBs_.getDVector());
        sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
        sampPtr->setOutputParams(nOutputs_);
        sampPtr->setSamplingParams(nSamples_, -1, -1);
        sampPtr->initialize(1);
        sampPtr->loadSamples(nSamples_,nInputs_,nOutputs_,
                VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
                VecSamStates_.getIVector());
        Sampling *sampAux = SamplingCreateFromID(SamplingMethod_);
        sampAux->setPrintLevel(outputLevel_);
        sampAux->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                                VecIUpperBs_.getDVector());
        sampAux->setInputParams(nInputs_, NULL, NULL, NULL);
        sampAux->setOutputParams(iOne);
        sampAux->setSamplingParams(nSamples_, -1, -1);
        sampAux->initialize(1);
        sampAux->loadSamples(nSamples_, nInputs_, iOne, 
                  VecSamInputs_.getDVector(), vecYT.getDVector(), 
                  VecSamStates_.getIVector());

        //**/ first refine uniformly
        strcpy(pString, "setUniformRefinement");
        sampAux->setParam(pString);
        status = sampAux->refine(2, 1, 0.0, nSamples_, NULL);
        int nSamples2 = sampAux->getNumSamples();
        if (nSamples2 != 2 * nSamples_)
        {
          printf("arefine_metis ERROR: consult developers.\n");
          printf("     refined sample size != 2 * original size\n");
          printf("     May be due to too many levels of refinement.\n");
          printf("     %d versus the expected %d\n",nSamples2,2*nSamples_);
          delete sampAux;;
          sampAux = NULL;
          delete sampPtr;;
          sampPtr = NULL;
          continue;
        }

        //**/ get the uniformly refined sample
        psVector  vecSamInp2, vecSamOut2, vecSamStd2;
        psIVector vecSamSta2;
        vecSamInp2.setLength(nSamples2*nInputs_);
        vecSamOut2.setLength(nSamples2);
        vecSamStd2.setLength(nSamples2);
        vecSamSta2.setLength(nSamples2);
        double *samInp2 = vecSamInp2.getDVector();
        double *samOut2 = vecSamOut2.getDVector();
        double *samStd2 = vecSamStd2.getDVector();
        sampAux->getSamples(nSamples2,nInputs_,1,vecSamInp2.getDVector(),
                     vecSamOut2.getDVector(),vecSamSta2.getIVector());

        //**/ create a MARSB response surface
        faType = PSUADE_RS_MARSB;
        FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
        faPtr->setNPtsPerDim(32);
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(0);

        //**/ generate 100 instantiations
        int numMars = 100, ivar1;
        psMatrix matMarsX, matMarsY;
        matMarsX.setFormat(PS_MAT2D);
        matMarsY.setFormat(PS_MAT2D);
        matMarsX.setDim(numMars, nSamples_*nInputs_);
        matMarsY.setDim(numMars, nSamples_);
        double **marsX = matMarsX.getMatrix2D();
        double **marsY = matMarsY.getMatrix2D();
        for (ii = 0; ii < numMars; ii++)
        {
          for (ss = 0; ss < nSamples_; ss++)
          {
            if (ss < origNSamples)
                 ivar1 = PSUADE_rand() % origNSamples;
            else ivar1 = ss;
            for (jj = 0; jj < nInputs_; jj++)
              marsX[ii][ss*nInputs_+jj] =
                       VecSamInputs_[ivar1*nInputs_+jj];
            marsY[ii][ss] = VecSamOutputs_[ivar1*nOutputs_+outputID];
          }
        }
        strcpy(cString, "mars_params");
        int ivar2 = 2 * nInputs_ / 3 + 1;
        int ivar3 = nSamples_;
        if (ivar3 > 300) ivar3 = 300;
        targv[0] = (char *) cString;
        targv[1] = (char *) &ivar3;
        targv[2] = (char *) &ivar2;
        faPtr->setParams(3, targv);
        strcpy(cString, "num_mars");
        targv[0] = (char *) cString;
        targv[1] = (char *) &numMars;
        faPtr->setParams(2, targv);
        strcpy(cString, "mars_sample");
        targv[0] = (char *) cString;
        targv[2] = (char *) &nSamples_;
        for (ii = 0; ii < numMars; ii++)
        {
          targv[1] = (char *) &ii;
          targv[3] = (char *) marsX[ii];
          targv[4] = (char *) marsY[ii];
          faPtr->setParams(5, targv);
        }
        faPtr->initialize(VecSamInputs_.getDVector(),
                          vecYT.getDVector());

        //**/ evaluate the uniformly refined sample's errors
        //**/ (only on the new points) ==> samStds2
        //**/ on new points instead of original point because prediction
        //**/ errors on those points should be 0, but we want prediction
        //**/ errors on new points
        faPtr->evaluatePointFuzzy(nSamples2-nSamples_,
            &samInp2[nSamples_*nInputs_],&samOut2[nSamples2-nSamples_],
            &samStd2[nSamples2-nSamples_]); 
        for (ss = 0; ss < nSamples2-nSamples_; ss++) 
          vecSamStd2[ss] = PABS(vecSamStd2[ss+nSamples_]);
        if (outputLevel_ > 3)
        {
          printf("Std. dev. to be used to select refinements.\n");
          for (ss = 0; ss < nSamples_; ss++) 
            printf("Sample point %7d: stdev = %e\n",ss+1,samStd2[ss]);
        }

        //**/ now refine based on magnitudes of the errors
        strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
        sampPtr->setParam(pString);
        sprintf(pString, "setRefineSize %d", refineSize);
        sampPtr->setParam(pString);
        sampPtr->refine(2,1,1.0e-6,nSamples_,samStd2);
        VecSamOutputs_.clean();
        VecSamStates_.clean();

        //**/ fetch the new sample 
        nSamples_ = sampPtr->getNumSamples();
        VecSamInputs_.setLength(nInputs_*nSamples_);
        VecSamOutputs_.setLength(nOutputs_*nSamples_);
        VecSamStates_.setLength(nSamples_);
        sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
                  VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(),
                  VecSamStates_.getIVector());
        psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                     VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                     VecSamOutputs_.getDVector(), 
                     VecSamStates_.getIVector(),NULL); 
        psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        printf("arefine_metis successful: use write to store ");
        printf("the refined sample. Then\n");
        printf("     run the newly created sample points in ");
        printf("this refined sample (the\n");
        printf("     last %d sample points).\n",refineSize);
        delete sampPtr;
        delete sampAux;
        sampPtr = sampAux = NULL;
        delete faPtr;
        faPtr = NULL;
        psConfig_.SamExpertModeRestore();
        cmdStatus = 0;
      }
      else
      {
        printf("This function uses the current loaded sample ");
        printf("together with previous\n");
        printf("refinement information (stored in file ");
        printf("'psuadeMetisInfo') to perform\n");
        printf("adaptive refinement. The sample has been ");
        printf("detected to have outputs that\n");
        printf("are either 0 or 1, so it is assumed that ");
        printf("adaptive refinement is to be\n");
        printf("applied to the 0/1 interfaces.\n");

        //**/ temporarily turn off sampling expert mode
        psConfig_.SamExpertModeSaveAndReset();
        Sampling *sampPtr = 
              (Sampling *) SamplingCreateFromID(SamplingMethod_);
        sampPtr->setPrintLevel(PL_INTERACTIVE);
        sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                                VecIUpperBs_.getDVector());
        sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
        sampPtr->setOutputParams(nOutputs_);
        sampPtr->setSamplingParams(nSamples_, -1, -1);
        sampPtr->initialize(1);
        sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, 
             VecSamInputs_.getDVector(), VecSamOutputs_.getDVector(), 
             VecSamStates_.getIVector());

        //**/ now refine based on magnitudes of the errors
        strcpy(pString, "setAdaptiveRefinementBasedOnOutputs");
        sampPtr->setParam(pString);
        sprintf(pString, "setRefineSize %d", nSamples_);
        sampPtr->setParam(pString);
        sampPtr->refine(2,1,1.0e-6,0,NULL);

        //**/ fetch the new sample 
        kk = nSamples_;
        nSamples_ = sampPtr->getNumSamples();
        VecSamInputs_.setLength(nInputs_*nSamples_);
        VecSamOutputs_.setLength(nOutputs_*nSamples_);
        VecSamStates_.setLength(nSamples_);
        sampPtr->getSamples(nSamples_,nInputs_,nOutputs_,
                 VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
                 VecSamStates_.getIVector());
        psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                     VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                     VecSamOutputs_.getDVector(), 
                     VecSamStates_.getIVector(),NULL); 
        psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        printf("arefine_metis successful: use write to store ");
        printf("the refined sample. Then\n");
        printf("     run the newly created sample points in ");
        printf("this refined sample (the\n");
        printf("     last %d sample points).\n",refineSize);
        delete sampPtr;
        sampPtr = NULL;
        psConfig_.SamExpertModeRestore();
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ arefine
    //**/ adaptive sample refinement
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arefine") ||
             !strcmp(command, "arefine_pu"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arefine: adaptively refine the loaded sample ");
        printf("based on prediction\n");
        printf("        uncertainty (std dev) at the new sample ");
        printf("points. This function,\n");
        printf("         unlike arefine_metis, does not assume ");
        printf("the loaded sample is a\n");
        printf("         METIS sample.");
        printf("Syntax: arefine (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("Given a sample which has been loaded (using ");
      printf("load or read_std), this\n");
      printf("command adaptively refine the sample near ");
      printf("areas in the input space\n");
      printf("with the largest prediction uncertainty. This ");
      printf("command, unlike the\n");
      printf("arefine_metis command, does not assume the ");
      printf("loaded sample is a METIS\n");
      printf("sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Select output to use for arefine (1 - %d) ",
              nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] == PSUADE_UNDEFINED ||
            VecSamStates_[ss] != 1)
        {
          printf("arefine ERROR: some sample outputs are undefined.\n");
          cmdStatus = 1;
          break;
        }
      }
      if (ss != nSamples_) continue;

      printf("Enter the original sample size (before any ");
      printf("refinement) below. This\n");
      printf("information is useful to prevent clustering of ");
      printf("new sample points. If\n");
      printf("you enter 0, the current sample size will be ");
      printf("used (i.e. no favorites\n");
      printf("in selecting bootstrapped samples for MarsBagg.\n");
      sprintf(pString,"Enter the original sample size (<= %d): ",
              nSamples_);
      int origNSamples = getInt(0, nSamples_, pString);
      if (origNSamples == 0) origNSamples = nSamples_;
      SamplingMethod_ = PSUADE_SAMP_GMETIS;
      sprintf(pString,"How many sample points to add? (1 - %d) ",
              nSamples_);
      int refineSize = getInt(1, nSamples_, pString);

      //**/ temporarily turn off sampling expert mode
      //**/ load sample to 2 samplers, 1 for uniform sampling (sampAux),
      //**/ and the other for adaptive sampling (sampPtr)
      //**/ The reason for the uniform sample is that it needs to give
      //**/ the prediction errors for all POSSIBLE refined points (in
      //**/ the uniformly refined mesh) for METIS refine function to
      //**/ select candidates
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      psConfig_.SamExpertModeSaveAndReset();
      Sampling *sampPtr = 
             (Sampling *) SamplingCreateFromID(SamplingMethod_);
      sampPtr->setPrintLevel(outputLevel_);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(nSamples_, -1, -1);
      sampPtr->initialize(1);
      sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, 
                VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
                VecSamStates_.getIVector());
      Sampling *sampAux = SamplingCreateFromID(SamplingMethod_);
      sampAux->setPrintLevel(outputLevel_);
      sampAux->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampAux->setOutputParams(iOne);
      sampAux->setSamplingParams(nSamples_, -1, -1);
      sampAux->initialize(1);
      sampAux->loadSamples(nSamples_, nInputs_, iOne, 
                   VecSamInputs_.getDVector(),vecYT.getDVector(), 
                   VecSamStates_.getIVector());

      //**/ first refine uniformly
      strcpy(pString, "changeInfoName");
      sampAux->setParam(pString);
      strcpy(pString, "setUniformRefinement");
      sampAux->setParam(pString);
      status = sampAux->refine(2, 1, 0.0, nSamples_, NULL);
      int nSamples2 = sampAux->getNumSamples();
      //**/ get the uniformly refined sample
      psVector VecSamInp2, VecSamOut2, VecSamStd2;
      psIVector VecSamSta2;
      VecSamInp2.setLength(nSamples2*nInputs_);
      VecSamOut2.setLength(nSamples2);
      VecSamSta2.setLength(nSamples2);
      VecSamStd2.setLength(nSamples2);
      sampAux->getSamples(nSamples2,nInputs_,iOne,
               VecSamInp2.getDVector(), 
               VecSamOut2.getDVector(), VecSamSta2.getIVector());

      //**/ create a MARSB response surface
      faType = PSUADE_RS_MARSB;
      if (psConfig_.RSExpertModeIsOn())
      {
        faType = -1;
        printf("Select a stochastic response surface method : ");
        printf("1. MARS with bootstrapped aggregation (MARSB)\n");
        printf("2. Gaussian process (GP)\n");
        printf("3. Kriging\n");
        printf("4. Sum of trees method\n");
        sprintf(pString,"Make your choice (default: MARSB) : ");
        faType = getInt(1,4,pString);
        switch (faType)
        {
          case 1: faType = PSUADE_RS_MARSB; break;
#ifdef HAVE_TPROS
          case 2: faType = PSUADE_RS_GP1; break;
#else
          case 2: faType = PSUADE_RS_GP3; break;
#endif
          case 3: faType = PSUADE_RS_KR; break;
          case 4: faType = PSUADE_RS_SOTS; break;
        }
      }
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      faPtr->setNPtsPerDim(32);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(0);
      //**/ generate 100 instantiations
      int numMars = 100, ivar1;
      psMatrix matMarsX, matMarsY;
      matMarsX.setFormat(PS_MAT2D);
      matMarsY.setFormat(PS_MAT2D);
      matMarsX.setDim(numMars, nSamples_*nInputs_);
      matMarsY.setDim(numMars, nSamples_);
      double **marsX = matMarsX.getMatrix2D();
      double **marsY = matMarsY.getMatrix2D();
      for (ii = 0; ii < numMars; ii++)
      {
        for (ss = 0; ss < nSamples_; ss++)
        {
          if (ss < origNSamples)
               ivar1 = PSUADE_rand() % origNSamples;
          else ivar1 = ss;
          for (jj = 0; jj < nInputs_; jj++)
            marsX[ii][ss*nInputs_+jj] =
                 VecSamInputs_[ivar1*nInputs_+jj];
          marsY[ii][ss] = VecSamOutputs_[ivar1*nOutputs_+outputID];
        }
      }
      strcpy(cString, "mars_params");
      int ivar2 = 2 * nInputs_ / 3 + 1;
      int ivar3 = nSamples_;
      if (ivar3 > 300) ivar3 = 300;
      targv[0] = (char *) cString;
      targv[1] = (char *) &ivar3;
      targv[2] = (char *) &ivar2;
      faPtr->setParams(3, targv);
      strcpy(cString, "num_mars");
      targv[0] = (char *) cString;
      targv[1] = (char *) &numMars;
      faPtr->setParams(2, targv);
      strcpy(cString, "mars_sample");
      targv[0] = (char *) cString;
      targv[2] = (char *) &nSamples_;
      for (ii = 0; ii < numMars; ii++)
      {
        targv[1] = (char *) &ii;
        targv[3] = (char *) marsX[ii];
        targv[4] = (char *) marsY[ii];
        faPtr->setParams(5, targv);
      }
      faPtr->initialize(VecSamInputs_.getDVector(),
                        vecYT.getDVector());

      //**/ evaluate the uniformly refined sample's errors
      double *samInp2 = VecSamInp2.getDVector();
      double *samOut2 = VecSamOut2.getDVector();
      double *samStd2 = VecSamStd2.getDVector();
      faPtr->evaluatePointFuzzy(nSamples2-nSamples_, 
                      &samInp2[nInputs_*nSamples_], 
                      &samOut2[nSamples_], &samStd2[nSamples_]);

      for (ss = 0; ss < nSamples2-nSamples_; ss++) 
        VecSamStd2[ss] = PABS(VecSamStd2[ss+nSamples_]);
      if (outputLevel_ > 3)
      {
         printf("Std dev to be used to select refinements.\n");
         for (ss = 0; ss < nSamples2-nSamples_; ss++) 
           printf("Sample point %7d: stdev = %e\n", ss+1,VecSamStd2[ss]);
      }

      //**/ now refine based on magnitudes of the errors
      strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
      sampPtr->setParam(pString);
      sprintf(pString, "setRefineSize %d", refineSize);
      sampPtr->setParam(pString);
      if (outputLevel_ > 3)
      {
        printf("Sample data for refinement: (standard deviations)\n");
        for (ss = 0; ss < nSamples_; ss++) 
          printf("%5d: std dev = %e\n", ss+1, VecSamStd2[ss]);
      }
      sampPtr->refine(2,1,1.0e-6,nSamples_,VecSamStd2.getDVector());

      //**/ fetch the new sample 
      nSamples_ = sampPtr->getNumSamples();
      VecSamInputs_.setLength(nInputs_*nSamples_);
      VecSamOutputs_.setLength(nOutputs_*nSamples_);
      VecSamStates_.setLength(nSamples_);
      sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
                          VecSamInputs_.getDVector(),
                          VecSamOutputs_.getDVector(), 
                          VecSamStates_.getIVector());
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),NULL); 
      psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("arefine successful: use write to store ");
      printf("the refined sample. Then run\n");
      printf("     the newly created sample points in ");
      printf("this refined sample (the last\n");
      printf("     %d sample points).\n",refineSize);
      delete sampPtr;
      delete sampAux;
      sampPtr = sampAux = NULL;
      delete faPtr;
      faPtr = NULL;
      psConfig_.SamExpertModeRestore();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arefine_cv
    //**/ adaptive sample refinement based on cross validation
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arefine_cv"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arefine_cv: adaptively refine a loaded sample ");
        printf("using N-fold CV to\n");
        printf("        determine which samples (the ones with ");
        printf("largest differences\n");
        printf("        between predicted and actual values) to ");
        printf("refine.\n");
        printf("Syntax: arefine_cv (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command adaptively refines the loaded ");
      printf("sample based on results of\n"); 
      printf("N-fold (N = sample size) cross validation ");
      printf("performed on the sample.\n");
      printf("The loaded sample should have good space-filling ");
      printf("property. The steps are:\n");
      printf("1. Select which output to use for CV\n");
      printf("2. Select how many new sample points to add\n");
      printf("3. Select response surface type\n");
      printf("At completion, the new points have been appended ");
      printf("to the sample.\n");
      printf("Also, set printlevel to 5 before this command ");
      printf("to see more screendump.\n");

      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample has not been loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Select output to use for refinement (1 - %d) ",
              nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] == PSUADE_UNDEFINED ||
            VecSamStates_[ss] != 1)
        {
          printf("arefine_cv ERROR: some sample outputs are undefined.\n");
          break;
        }
      }
      if (ss != nSamples_) continue;

      psuadeIO_->getParameter("method_sampling", pPtr);
      SamplingMethod_ = PSUADE_SAMP_GMETIS;
      sprintf(pString,"How many sample points to add? (1 - %d) ",nSamples_);
      int refineSize = getInt(1, nSamples_, pString);
      //**/ temporarily turn off sampling expert mode
      //**/ load sample to 2 samplers, one for uniform sampling (sampAux), 
      //**/ and the other for adaptive sampling (sampPtr)
      psConfig_.SamExpertModeSaveAndReset();
      Sampling *sampPtr = 
           (Sampling *) SamplingCreateFromID(SamplingMethod_);
      sampPtr->setPrintLevel(outputLevel_);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(nSamples_, -1, -1);
      sampPtr->initialize(1);
      sampPtr->loadSamples(nSamples_, nInputs_, nOutputs_, 
                 VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
                 VecSamStates_.getIVector());
      Sampling *sampAux = SamplingCreateFromID(SamplingMethod_);
      sampAux->setPrintLevel(PL_INTERACTIVE);
      sampAux->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampAux->setOutputParams(iOne);
      sampAux->setSamplingParams(nSamples_, -1, -1);
      sampAux->initialize(1);
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      sampAux->loadSamples(nSamples_, nInputs_, iOne, 
                   VecSamInputs_.getDVector(), 
                   vecYT.getDVector(), VecSamStates_.getIVector());

      //**/ first refine uniformly
      strcpy(pString, "changeInfoName");
      sampAux->setParam(pString);
      strcpy(pString, "setUniformRefinement");
      sampAux->setParam(pString);
      status = sampAux->refine(2, 1, 0.0, nSamples_, NULL);
      int nSamples2 = sampAux->getNumSamples();
      //**/ get the uniformly refined sample
      psVector  vecSamInp2, vecSamOut2;
      psIVector vecSamSta2;
      vecSamInp2.setLength(nSamples2*nInputs_);
      vecSamOut2.setLength(nSamples2);
      vecSamSta2.setLength(nSamples2);
      double *samOut2  = vecSamOut2.getDVector();
      sampAux->getSamples(nSamples2,nInputs_,1,vecSamInp2.getDVector(),samOut2,
                          vecSamSta2.getIVector());
      //**/ create a response surface
      faType = PSUADE_RS_RBF;
      if (psConfig_.RSExpertModeIsOn())
      {
        faType = -1;
        sprintf(pString, "Enter your choice ? ");
        while (faType < 0 || faType > PSUADE_NUM_RS)
        {
          writeFAInfo(0);
          faType = getFAType(pString);
        }
      }
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_-1);
      faPtr->setNPtsPerDim(32);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(0);
      psVector VecTX, VecTY;
      VecTX.setLength(nSamples_*nInputs_);
      VecTY.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
      {
        count = 0;
        for (kk = 0; kk < nSamples_; kk++)
        {
          if (kk != ss)
          {
            for (jj = 0; jj < nInputs_; jj++)
              VecTX[count*nInputs_+jj] = VecSamInputs_[kk*nInputs_+jj];
            VecTY[count] = VecSamOutputs_[kk*nOutputs_+outputID];
            count++;
          }
        }
        faPtr->initialize(VecTX.getDVector(),VecTY.getDVector());

        //**/ evaluate the uniformly refined sample's errors
        double *dPtr = VecSamInputs_.getDVector();
        faPtr->evaluatePoint(iOne, &dPtr[ss*nInputs_], &samOut2[ss]);

        ddata = samOut2[ss] - VecSamOutputs_[ss*nOutputs_+outputID];
        samOut2[ss] = PABS(ddata);
      }
      if (outputLevel_ > 4)
      {
        printf("CV predicted versus actual:\n");
        for (ss = 0; ss < nSamples_; ss++) 
          printf("Sample point %7d: diff = %e\n",ss+1,samOut2[ss]);
      }

      //**/ now refine based on magnitudes of the errors
      strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
      sampPtr->setParam(pString);
      sprintf(pString, "setRefineSize %d", refineSize);
      sampPtr->setParam(pString);
      sampPtr->refine(2,1,1.0e-6,nSamples_,samOut2);

      //**/ fetch the new sample 
      nSamples_ = sampPtr->getNumSamples();
      VecSamInputs_.setLength(nInputs_*nSamples_);
      VecSamOutputs_.setLength(nOutputs_*nSamples_);
      VecSamStates_.setLength(nSamples_);
      sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
                 VecSamInputs_.getDVector(),
                 VecSamOutputs_.getDVector(),VecSamStates_.getIVector());
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),NULL); 
      psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("arefine_cv: Now use write to store the refined sample.\n");
      printf("NOTE: the new points are the last %d sample points.\n", 
             refineSize);
      delete sampPtr;
      delete sampAux;
      sampPtr = sampAux = NULL;
      delete faPtr;
      faPtr = NULL;
      psConfig_.SamExpertModeRestore();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arel_fact
    //**/ adaptive sample refinement along the threshold contour
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arel_fact"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arel_fact: adaptively refine the ");
        printf("loaded sample based on some\n"); 
        printf("threshod user-specified threshold. New ");
        printf("sample points will be created\n");
        printf("near the threshold boundaries in the ");
        printf("input space (that separate regions\n");
        printf("with output <= threshold and those with ");
        printf("output > threshold.\n");
        printf("NOTE: This command uses factorial design ");
        printf("for refinement. As such, it\n");
        printf("      is not recommended for samples ");
        printf("with nInputs large (> 3 or 4).\n");
        printf("Syntax: arel_fact (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command adaptively refines the loaded ");
      printf("sample based on some user-\n"); 
      printf("specified threshold. New sample points will be ");
      printf("created near the\n");
      printf("threshold boundaries in the input space ");
      printf("(that separate regions with\n");
      printf("output <= threshold and those with output > ");
      printf("threshold.) The steps are:\n");
      printf("1. Select which output to use for thresholding\n");
      printf("2. Specify the threshold value\n");
      printf("3. Select response surface type\n");
      printf("4. Select resolution of the m-dimensional lattice\n");
      printf("   (for use in determining where the ");
      printf("threshold boundary is).\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample has not been loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ > 3)
      {
        printf("WARNING: for nInputs > 3, too many new refined ");
        printf("points may be created\n");
        printf("         since the thresh boundary is ");
        printf("high-dimensional (> 2).\n");
        cmdStatus = 1;
        continue;
      }

      //**/ ask for which output
      sprintf(pString,"Select output to use for refinement (1 - %d) ",
              nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] == PSUADE_UNDEFINED ||
            VecSamStates_[ss] != 1)
        {
          printf("arel_fact ERROR: some sample outputs ");
          printf("are undefined or their\n");
          printf("                       states are not ready.\n");
          break;
        }
      }
      if (ss != nSamples_) continue;

      //**/ get min/max information from sample output
      double dMin=1e35, dMax=-1e-35;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] < dMin)
          dMin = VecSamOutputs_[ss*nOutputs_+outputID];
        if (VecSamOutputs_[ss*nOutputs_+outputID] > dMax)
          dMax = VecSamOutputs_[ss*nOutputs_+outputID];
      }
      printf("Sample minimum and maximum are %e and %e.\n",dMin,dMax);
      printf("The threshold should be somewhere in between ");
      printf("or not too far away.\n");
      sprintf(pString, "Please enter threshold : ");
      thresh = getDouble(pString);

      //**/ ask to see if RS is needed
      int useRS=0;
      printf("Adaptive refinement used in this command is ");
      printf("performed on a regular\n");
      printf("grid or lattice structure. If your loaded sample ");
      printf("already has a lattice\n");
      printf("structure (a factorial design), or it is an ");
      printf("expanded sample from\n");
      printf("repeatedly running arel_fact with ");
      printf("finer resolutions (factor\n");
      printf("of 2 each time), then it may be possible to ");
      printf("perform refinement on\n");
      printf("the sample directly without using response ");
      printf("surface. Otherwise, a\n");
      printf("response surface is needed for interpolation ");
      printf("at different parts of\n");
      printf("the input space. Even if the loaded sample has ");
      printf("a lattice or expanded\n");
      printf("lattice structure, you can still use a ");
      printf("response surface. So, \n");
      printf("Use response surface ? (y or n) ");
      scanf("%s", lineIn2);
      fgets(lineIn,5000,stdin); 
      if (lineIn2[0] == 'y') useRS = 1;

      //**/ see if the sample is a factorial sample
      psIVector vecResolution;
      vecResolution.setLength(nInputs_);
      printf("A lattice will be created for adaptive refinement. ");
      printf("Please enter the\n");
      printf("resolution (resolution = number of distinct ");
      printf("equally-spaced points) for\n");
      printf("each input below.\n");
      if (useRS == 0)
      {
        printf("NOTE: If your loaded sample is a factorial ");
        printf("sample of, e.g. MxN (2D),\n");
        printf("      then the resolutions you should enter ");
        printf("are M and N.\n");
      }
      printf("NOTE: If you have run arel_fact at least ");
      printf("once, the resolutions you\n");
      printf("      should enter are the ones suggested ");
      printf("upon completion of the last\n");
      printf("      arel_fact.\n");
      long lcount = 1;
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString, 
                "Enter resolution for input %d : ", ii+1);
        vecResolution[ii] = getInt(2,1000000,pString);
        lcount *= vecResolution[ii];
      }
      if (lcount > 10000000)
      {
        printf("WARNING: the resolutions may be too large ");
        printf("to be handled.\n");
        printf("         Total lattice points = %ld\n", lcount);
      }

      //**/ create a response surface
      FuncApprox *faPtr = NULL; 
      if (useRS == 1)
      {
        faType = -1;
        sprintf(pString, "Enter your choice for response surface : ");
        while (faType < 0 || faType > PSUADE_NUM_RS)
        {
          writeFAInfo(0);
          faType = getFAType(pString);
        }
        faPtr = genFA(faType, nInputs_, iOne, nSamples_);
        faPtr->setNPtsPerDim(256);
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        vecYT.setLength(nSamples_);
        for (ss = 0; ss < nSamples_; ss++)
          vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
        faPtr->initialize(VecSamInputs_.getDVector(),vecYT.getDVector());
      }

      //**/  create a nInputs-dimensional grid
      int nEvals = (int) lcount;
      FactorialSampling *factPtr = new FactorialSampling();
      factPtr->setPrintLevel(outputLevel_);
      factPtr->setInputBounds(nInputs_,VecILowerBs_.getDVector(),
                              VecIUpperBs_.getDVector());
      factPtr->setSamplingParams(nEvals, 1, 0);
      factPtr->setOutputParams(iOne);
      factPtr->setInputParams(nInputs_,vecResolution.getIVector(),
                              NULL,NULL);
      factPtr->initialize(0);
      psVector  vecEvalX, vecEvalY, vecEvalStd;
      psIVector vecEvalS;
      vecEvalX.setLength(nEvals*nInputs_);
      vecEvalY.setLength(nEvals);
      vecEvalStd.setLength(nEvals);
      vecEvalS.setLength(nEvals);
      factPtr->getSamples(nEvals,nInputs_,iOne,vecEvalX.getDVector(),
                vecEvalY.getDVector(),vecEvalS.getIVector());

      //**/ evaluate the lattice sample points 
      if (useRS == 1)
      {
        printf("Since a response surface is used in this analysis, \n");
        printf("you may want to\n");
        printf("account for response surface uncertainties (if \n");
        printf("a stochastic response\n");
        printf("surface has bee selected). One way to account \n");
        printf("for this uncertainty is\n");
        printf("to add some multiple of prediction std devs \n");
        printf("to the predictions. E.g.\n");
        printf("Prediction mean + 1 x prediction std dev ");
        printf("==> 84.0%% confidence\n");  
        printf("Prediction mean + 2 x prediction std dev ");
        printf("==> 97.5%% confidence\n");  
        printf("Prediction mean + 3 x prediction std dev ");
        printf("==> 99.5%% confidence\n");  
        sprintf(pString, 
                "Please enter this multiple factor (0, 1, 2, 3) : ");
        double multiplier = getDouble(pString);
        faPtr->evaluatePointFuzzy(nEvals, vecEvalX.getDVector(),
                             vecEvalY.getDVector(),
                             vecEvalStd.getDVector());
        //**/ subtract m*stds
        for (ss = 0; ss < nEvals; ss++)
          vecEvalY[ss] += (multiplier*vecEvalStd[ss]);
        delete faPtr;
      }
      else
      {
        //**/ load sample outputs to vecEvalX and vecEvalY
        count = 0;
        for (ss = 0; ss < nEvals; ss++) vecEvalY[ss] = PSUADE_UNDEFINED;
        psIVector vecTags;
        vecTags.setLength(nSamples_);
        for (ss = 0; ss < nEvals; ss++)
        {
          for (kk = 0; kk < nSamples_; kk++)
          {
            for (ii = 0; ii < nInputs_; ii++)
            {
              ddata = VecSamInputs_[kk*nInputs_+ii]-vecEvalX[ss*nInputs_+ii];
              if (PABS(ddata) > 1e-15) break;
            }
            if (ii == nInputs_)
            {
              //printf("Point found : %d = %d\n", ss+1, kk+1);
              //for (ii = 0; ii < nInputs_; ii++)
              //  printf("Input %d = %e\n", ii+1, vecEvalX[ss*nInputs_+ii]);
              vecEvalY[ss] = VecSamOutputs_[kk*nOutputs_+outputID];
              vecTags[kk] = 1;
              count++;
              break;
            }
          }
        }
        printf("INFO: %d of %d sample points have ", count, nSamples_);
        printf("been identified as points coinciding\n");
        printf("     with the created lattice.\n");
        if (count < nSamples_)
        {
          printf("ERROR: Not all sample points are used. ");
          printf("Possible errors are:\n"); 
          printf("  - The loaded sample is not ");
          printf("a lattice or expanded lattice.\n");
          printf("  - The resolutions you entered are not correct.\n");
          for (ii = 0; ii < nSamples_; ii++)
          {
            if (vecTags[ii] == 0) 
            {
              printf("Sample %d has not found a match on the lattice.\n",
                     ii+1);
              for (jj = 0; jj < nInputs_; jj++)
                printf("  Sample input %d = %e\n",jj+1,
                  VecSamInputs_[ii*nInputs_+jj]);
            }
          }
          exit(1);
        }
      }

      //**/ define region to search
      psIVector vecDirections;
      vecDirections.setLength(nInputs_);
      printf("To better identify the true threshold boundary, ");
      printf("please answer the\n");
      printf("following question about your knowledge ");
      printf("about the simulation model.\n");
      printf("Select an option for each input dimension ");
      printf("X_i from the following:\n");
      printf(" 1. Y(X_ij) > Y(X_ik) if X_ij > X_ik ");
      printf("(monotonically increasing)\n");
      printf(" 2. Y(X_ij) < Y(X_ik) if X_ij > X_ik ");
      printf("(monotonically decreasing)\n");
      printf(" 3. Don't know. Assume both 1 and 2 are possible.\n");
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString,"Enter option for input %d (1 - 3) : ",ii+1);
        vecDirections[ii] = getInt(1,3,pString);
      }

      //**/ search for boundaries
      psVector  vecCandX, vecCandY, vecNeighY;
      psIVector vecEvalTags, vecCnts;
      vecCandX.setLength(nInputs_*nEvals);
      vecCandY.setLength(nEvals);
      vecNeighY.setLength(nInputs_*nEvals*2);
      vecEvalTags.setLength(nEvals);

      //**/ disable undefined points
      int nCand = 0, prodM;
      for (ss = 0; ss < nEvals; ss++)
        if (vecEvalY[ss] == PSUADE_UNDEFINED) vecEvalTags[ss] = 1;
      for (ss = nEvals-1; ss >= 0; ss--)
      {
        //**/ if I am undefined, let others select me
        if (vecEvalY[ss] != PSUADE_UNDEFINED)
        {
          prodM = 1;
          for (ii = 0; ii < nInputs_; ii++)
          {
            vecCnts.setLength(nInputs_);
            if ((ss / prodM) % vecResolution[ii] > 0)
            {
              //**/ save my neighbor regardless
              vecNeighY[nCand*2*nInputs_+ii*2] = vecEvalY[ss-prodM];
              //**/ if neighbors are defined and on diferent side
              //**/ of the boundary, just turn it on
              if (vecEvalY[ss] > thresh && vecEvalY[ss-prodM] <= thresh &&
                  (vecDirections[ii] == 1 || vecDirections[ii] == 3))
                vecCnts[ii] = -1;
              //**/ if I am larger than threshold and my neighbor is
              //**/ undefined, set a tag
              if (vecEvalY[ss] > thresh && 
                  vecEvalY[ss-prodM] == PSUADE_UNDEFINED &&
                  (vecDirections[ii] == 1 || vecDirections[ii] == 3))
              {
                vecCnts[ii] = -2;
                //**/ if any point in my downstream is larger than
                //**/ threshold, then my point is not needed
                for (jj = ((ss/prodM)%vecResolution[ii])-1; jj > 1; jj--)
                  if (vecEvalY[ss-jj*prodM] != PSUADE_UNDEFINED &&
                      vecEvalY[ss-jj*prodM] > thresh) vecCnts[ii] = 0;
              }
            }
            if ((ss / prodM) % vecResolution[ii] < vecResolution[ii]-1)
            {
              //**/ save my neighbor regardless
              vecNeighY[nCand*2*nInputs_+ii*2+1] = vecEvalY[ss+prodM];
              //**/ if neighbors are defined and on diferent side
              //**/ of the boundary, just turn it on
              if (vecEvalY[ss] > thresh && vecEvalY[ss+prodM] <= thresh &&
                  (vecDirections[ii] == 2 || vecDirections[ii] == 3))
                vecCnts[ii] = 1;
              //**/ if I am larger than threshold and my neighbor is
              //**/ undefined, set a tag
              if (vecEvalY[ss] > thresh && 
                  vecEvalY[ss+prodM] == PSUADE_UNDEFINED &&
                  (vecDirections[ii] == 2 || vecDirections[ii] == 3))
              {
                vecCnts[ii] = 2;
                //**/ if any point in my upstream is larger than
                //**/ threshold, then my point is not needed
                for (jj = ((ss/prodM)%vecResolution[ii])-1; jj > 1; jj--)
                  if (vecEvalY[ss+jj*prodM] != PSUADE_UNDEFINED &&
                      vecEvalY[ss+jj*prodM] > thresh) vecCnts[ii] = 0;
              }
            }
            //**/ if some points are tagged
            if (vecCnts[ii] != 0)
            {
              vecEvalTags[ss] = 1;
              for (jj = 0; jj < nInputs_; jj++)
              {
                //**/ pick my neighbor because my neighbor is an
                //**/ undefined and I am larger than threshold
                if (ii == jj && vecCnts[ii] == -2)
                  vecCandX[nCand*nInputs_+jj] = 
                     vecEvalX[(ss-prodM)*nInputs_+jj];
                //**/ pick midpoint between me and my neighbor
                //**/ because my neighbor is <= thresh, and I am larger
                else if (ii == jj && vecCnts[ii] == -1)
                  vecCandX[nCand*nInputs_+jj] = 0.5 * 
                    (vecEvalX[ss*nInputs_+jj]+
                     vecEvalX[(ss-prodM)*nInputs_+jj]);
                //**/ pick midpoint between me and my neighbor
                //**/ because my neighbor is <= thresh, and I am larger
                else if (ii == jj && vecCnts[ii] == 1)
                  vecCandX[nCand*nInputs_+jj] = 0.5 * 
                    (vecEvalX[ss*nInputs_+jj]+
                     vecEvalX[(ss+prodM)*nInputs_+ii]);
                //**/ pick my neighbor because my neighbor is an
                //**/ undefined and I am larger than threshold
                else if (ii == jj && vecCnts[ii] == 2)
                  vecCandX[nCand*nInputs_+jj] = 
                     vecEvalX[(ss+prodM)*nInputs_+jj];
                //**/ for other inputs, keep them the sam
                else 
                  vecCandX[nCand*nInputs_+jj] = 
                     vecEvalX[ss*nInputs_+jj];
              }
              vecCandY[nCand] = vecEvalY[ss];
              //**/ search for duplicates
              for (kk = 0; kk < nCand; kk++)
              {
                for (jj = 0; jj < nInputs_; jj++)
                  if (vecCandX[kk*nInputs_+jj] != 
                      vecCandX[nCand*nInputs_+jj]) break;
                if (jj == nInputs_) break;
              }
              if (kk == nCand)
              {
                printf("Candidate %d: \n",ss+1);
                for (jj = 0; jj < nInputs_; jj++)
                  printf("  X %d = %e\n",jj+1,
                    vecCandX[nCand*nInputs_+jj]);
                printf("  Y = %e\n",vecCandY[nCand]);
                nCand++;
              }
            }
            prodM *= vecResolution[ii];
          }
        }
      }

      //**/ filled in the undefined cell for computing volume
      if (nInputs_ == 2 && vecDirections[0] == 1 && 
          vecDirections[1] == 1)
      {
        vecCnts.setLength(vecResolution[1]);
        for (ss = 0; ss < nEvals; ss+=vecResolution[0])
        {
          ind = ss / vecResolution[0];
          //**/ search from left to right to find the first
          //**/ location whereby EvalY > thresh and register
          //**/ the length of the EvalY > thresh elements
          for (kk = 0; kk < vecResolution[0]; kk++)
          {
            if (vecEvalY[ss+kk] != PSUADE_UNDEFINED &&
                vecEvalY[ss+kk] > thresh)
            {
              vecCnts[ind] = vecResolution[0] - kk;
              break;
            }
          }
        }
        //**/ search from bottom to top the first row
        //**/ whereby there is an entry > thresh
        for (ss = 0; ss < vecResolution[1]; ss++)
          if (vecCnts[ss] > 0) break;
        ss++;
        //**/ interpolate the vecCnts in between
        //**/ because there are holes
        while (ss < vecResolution[1])
        {
          //**/ search for the next vecCnts == nonzero
          for (kk = ss; kk < vecResolution[1]; kk++)
            if (vecCnts[kk] > 0) break;
          //**/ now vecCnts between (ss-1,kk) are holes
          //**/ that needs to be filled
          for (jj = ss; jj < kk; jj++)
          {
            ddata = 1.0 * (vecCnts[kk]-vecCnts[ss-1])/(kk-ss+1)*
                    (jj-ss+1) + vecCnts[ss-1];
            vecCnts[jj] = (int) ddata;
          }
          ss = kk + 1;
        }
        //**/ now fill the nonzeros and compute volume
        for (ss = 0; ss < nEvals; ss+=vecResolution[0])
        {
          ind = ss / vecResolution[0];
          for (kk = vecResolution[0]-vecCnts[ind]; 
               kk < vecResolution[0]; kk++)
            vecEvalY[ss+kk] = thresh + 1e-8;
        }
        double minVol = 0, maxVol;
        for (ss = 0; ss < vecResolution[1]; ss++)
        {
          for (kk = 0; kk < vecResolution[0]; kk++)
          {
            count = 0;
            ind = ss * vecResolution[0] + kk;
            if (vecEvalY[ind] != PSUADE_UNDEFINED &&
                vecEvalY[ind] > thresh) count++;
            ind = ss * vecResolution[0] + kk + 1;
            if (vecEvalY[ind] != PSUADE_UNDEFINED &&
                vecEvalY[ind] > thresh) count++;
            ind = (ss + 1) * vecResolution[0] + kk;
            if (vecEvalY[ind] != PSUADE_UNDEFINED &&
                vecEvalY[ind] > thresh) count++;
            ind = (ss + 1) *vecResolution[0] + kk + 1;
            if (vecEvalY[ind] != PSUADE_UNDEFINED &&
                vecEvalY[ind] > thresh) count++;
            if (count == 4) minVol += 1;
            if (count >= 1) maxVol += 1;
          }
        }
        count = (vecResolution[0] - 1) * (vecResolution[1] - 1);
        printf("Min percentage above threshold = %7.3f\n",
               100*minVol/count);
        printf("Max percentage above threshold = %7.3f\n",
               100*maxVol/count);
        printf("Based on %d cells.\n", count);
        printEquals(PL_INFO, 0);
        if (outputLevel_ > 1)
        {
          printf("2D Mesh: (# - undefined; X,O - above/below threshold)\n");
          for (jj = vecResolution[1]-1; jj >= 0; jj--)
          {
            for (ii = 0; ii < vecResolution[0]; ii++)
            {
              if (vecEvalY[ii+jj*vecResolution[0]] == PSUADE_UNDEFINED)
                printf("# ");
              else if (vecEvalY[ii+jj*vecResolution[0]] > thresh)
                printf("X ");
              else 
                printf("O ");
            }
            printf("\n");
          }
        }
      }

      //**/ now vecEvalX and vecEvalY have been filled
      if (nCand == 0) 
      {
        vecCandX.clean();
        printf("INFO: no threshold interface point has been found.\n");
        continue;
      }

      //**/ store the results
      printf("Number of threshold boundary points identified = %d\n",
             nCand);
      fp = fopen("arel_fact.std", "w");
      if (fp != NULL)
      {
        fprintf(fp,"# Each line: m inputs X, Y, Y(2 * nInputs) neighbors)\n");
        fprintf(fp,"%d %d %d\n", nCand, nInputs_, 2*nInputs_+1);
        for (ss = 0; ss < nCand; ss++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp,"%12.4e ", vecCandX[ss*nInputs_+ii]);
          fprintf(fp,"%12.4e ", vecCandY[ss]);
          for (ii = 0; ii < 2*nInputs_; ii++)
            fprintf(fp,"%12.4e ",vecNeighY[nCand*2*nInputs_+ii]);
          fprintf(fp, "\n");
        }
        fclose(fp); 
        printAsterisks(PL_INFO, 0);
        printf("INFO: threshold interface points are in ");
        printf("arel_fact.std\n");
        printf("      Number of interface points = %d\n",nCand);
        printf("NOTE: this number may be too large or too small ");
        printf("for your need.\n");
        printf("If this is too large, you have three options: \n");
        printf("1. Inspect the new points visually and make ");
        printf("your selection\n");
        printf("2. Re-run this with smaller resolution to get ");
        printf("smaller set\n");
        printf("3. Load this set into PSUADE and use odoe_mmd ");
        printf("to downselect\n");
        printf("NOTE: This last option may be computationally ");
        printf("very expensive.\n");
        if (useRS == 0)
        {
          printf("ALSO: In the next round of refinement, PSUADE ");
          printf("recommends using a finer\n");
          printf("resolution. Specifically,\n");
          for (ii = 0; ii < nInputs_; ii++)
            printf("   Input %3d resolution = %d\n",ii+1,
                   2*vecResolution[ii]-1);
        }
        else
        {
          printf("HINT: if input i is insensitive, you can keep ");
          printf("the same resolution.\n");
        }
        printAsterisks(PL_INFO, 0);
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arel_fact2
    //**/ adaptive sample refinement along the threshold contour
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arel_fact2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arel_fact2: adaptively refine the ");
        printf("loaded sample based on some user-\n");
        printf("provided threshold value. New sample ");
        printf("points will be created along the\n"); 
        printf("threshold boundaries in the input");
        printf("space (that separate regions with\n");
        printf("output <= threshold and those with ");
        printf("output > threshold.\n");
        printf("NOTE: This command only works for nInputs=2.\n");
        printf("Syntax: arel_fact2 (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command adaptively refines the loaded ");
      printf("sample based on some user-\n"); 
      printf("specified threshold. New sample points will be ");
      printf("created near the\n");
      printf("threshold boundaries in the input space ");
      printf("(that separate regions with\n");
      printf("output <= threshold and those with output > ");
      printf("threshold.) The steps are:\n");
      printf("1. Select which output to use for thresholding\n");
      printf("2. Specify the threshold value\n");
      printf("3. Select resolution of the m-dimensional lattice\n");
      printf("   (for use in determining where the ");
      printf("threshold boundary is).\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample has not been loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ !=2)
      {
        printf("ERROR: This command only works for nInputs = 2\n");
        cmdStatus = 1;
        continue;
      }

      //**/ ask for which output
      sprintf(pString,"Select output to use for refinement (1 - %d) ",
              nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] == PSUADE_UNDEFINED ||
            VecSamStates_[ss] != 1)
        {
          printf("arel2_fact ERROR: some sample outputs ");
          printf("are undefined or their\n");
          printf("                       states are not ready.\n");
          break;
        }
      }
      if (ss != nSamples_) continue;

      //**/ get min/max information from sample output
      double dMin=1e35, dMax=-1e-35;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] < dMin)
          dMin = VecSamOutputs_[ss*nOutputs_+outputID];
        if (VecSamOutputs_[ss*nOutputs_+outputID] > dMax)
          dMax = VecSamOutputs_[ss*nOutputs_+outputID];
      }
      printf("Sample minimum and maximum are %e and %e.\n",dMin,dMax);
      printf("The threshold should be somewhere in between ");
      printf("or not too far away.\n");
      sprintf(pString, "Please enter threshold : ");
      thresh = getDouble(pString);

      //**/ see if the sample is a factorial sample
      psIVector vecResolution;
      vecResolution.setLength(nInputs_);
      printf("A 2D mesh will be created for adaptive refinement. ");
      printf("Please enter the\n");
      printf("resolution (number of equally-spaced cells) ");
      printf("for each input below.\n");
      printf("NOTE: The first you run this command, the ");
      printf("sample you have loaded\n");
      printf("      should be a factorial design, e.g. MxN, ");
      printf("and the resolutions\n");
      printf("      you should enter are M and N.\n");
      printf("NOTE: Subsequent call to arel2_fact ");
      printf("should use the resolutions\n");
      printf("      suggestion upon completion of the previous ");
      printf("arel2_fact call.\n");
      long lcount = 1;
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString, 
                "Enter resolution for input %d : ", ii+1);
        vecResolution[ii] = getInt(2,1000000,pString);
        lcount *= vecResolution[ii];
      }
      if (lcount > 10000000)
      {
        printf("WARNING: the resolutions may be too large ");
        printf("to be handled.\n");
        printf("         Total lattice points = %ld\n", lcount);
      }

      //**/ create a 2D mesh
      psMatrix matMeshY;
      matMeshY.setFormat(PS_MAT2D);
      matMeshY.setDim(vecResolution[0],vecResolution[1]);
      double **meshY = matMeshY.getMatrix2D();
      psIMatrix matMeshTags;
      matMeshTags.setFormat(PS_MAT2D);
      matMeshTags.setDim(vecResolution[0],vecResolution[1]);
      int **meshTags = matMeshTags.getIMatrix2D();
      int ss1, ss2;
      for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
        for (ss2 = 0; ss2 < vecResolution[1]; ss2++)
          meshY[ss1][ss2] = PSUADE_UNDEFINED;

      //**/ load sample outputs to matMeshY
      psVector vecDX;
      vecDX.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
        vecDX[ii] = (VecIUpperBs_[ii]-VecILowerBs_[ii])/
                     vecResolution[ii];;  
      int ind1, ind2;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (outputLevel_ > 2)
        {
          printf("Sample point %d : %e %e\n",ss+1,
             VecSamInputs_[ss*nInputs_],VecSamInputs_[ss*nInputs_+1]);
        }
        ddata = VecSamInputs_[ss*nInputs_];
        ddata -= VecILowerBs_[0];
        ddata /= vecDX[0];
        ind1 = (int) (ddata+1e-15); 
        if (ind1 == vecResolution[0]) ind1--;
        ddata = VecSamInputs_[ss*nInputs_+1];
        ddata -= VecILowerBs_[1];
        ddata /= vecDX[1];
        ind2 = (int) (ddata+1e-15); 
        if (ind2 == vecResolution[1]) ind2--;
        if (outputLevel_ > 2)
          printf("   is mapped to Mesh(%d,%d)\n",ind1,ind2);
        meshY[ind1][ind2] = VecSamOutputs_[ss*nOutputs_+outputID];
        meshTags[ind1][ind2] = meshTags[ind1][ind2] + 1;
      }
      for (ss2 = vecResolution[1]-1; ss2 >= 0; ss2--)
      {
        for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
        {
          if (meshTags[ss1][ss2] > 1)
          {
            printf("ERROR: some cells have more than one points.\n");
            continue;
          }
        }
      }

      //**/ display mesh and check
      if (outputLevel_ > 3)
      {
        printf("After sample points have been mapped to ");
        printf("mesh (1 = occupied):\n");
        count = 0;
        for (ss2 = vecResolution[1]-1; ss2 >= 0; ss2--)
        {
          for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
          {
            printf("%1d ", meshTags[ss1][ss2]);
            if (meshTags[ss1][ss2] > 0) count++;
          }
          printf("\n");
        }
      }

      //**/ filled in the undefined cell for computing volume
      psIVector vecCnts;
      vecCnts.setLength(vecResolution[1]);
      //**/ scan all rows from bottom to top
      for (ss2 = 0; ss2 < vecResolution[1]; ss2++)
      {
        //**/ search from left to right to find the first
        //**/ location whereby EvalY > thresh and register
        //**/ the length of the EvalY > thresh elements
        //**/ vecCnts[ss2] holds number of > thresh elements
        //**/ in row ss2
        for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
        {
          if (meshY[ss1][ss2] != PSUADE_UNDEFINED &&
              meshY[ss1][ss2] > thresh)
          {
            vecCnts[ss2] = vecResolution[0] - ss1;
            break;
          }
        }
      }
      count = vecCnts.sum();
      if (count == 0)
      {
        printf("ERROR: threshold boundary not found.\n");
        continue;
      }

      //**/ e.g.
      //**/ 0 0 0 X X X (row 5) vecCnts[5] = 3
      //**/ 0 0 0 X X X         vecCnts[4] = 3
      //**/ 0 0 0 0 X X         vecCnts[3] = 2
      //**/ 0 0 0 0 0 0         vecCnts[2] = 0
      //**/ 0 0 0 0 0 0         vecCnts[1] = 0
      //**/ 0 0 0 0 0 0 (row 0) vecCnts[0] = 0
      //**/ search from bottom to top the first row whereby
      //**/ there is an entry > thresh (row 3 in the example)
      for (ss = 0; ss < vecResolution[1]; ss++)
        if (vecCnts[ss] > 0) break;
      ss++; /* starting from ss=4 */
      //**/ interpolate the vecCnts in between
      //**/ because there are holes
      while (ss < vecResolution[1])
      {
        //**/ search for the next vecCnts == nonzero
        //**/ vecCnts[4] = 3 so the first kk=4
        for (kk = ss; kk < vecResolution[1]; kk++)
          if (vecCnts[kk] > 0) break;
        //**/ now vecCnts between (ss-1,kk) are holes
        //**/ that needs to be filled
        //**/ for the example, there is nothing between
        //**/ (ss-1=3,4) initially, but if so, the rows
        //**/ between will be the minimum vecCnts[ss-1] 
        //**/ vecCnts[kk)
        if (vecCnts[ss-1] < vecCnts[kk])
        {
          for (jj = ss; jj < kk; jj++)
            vecCnts[jj] = vecCnts[ss-1];
        }
        else
        {
          for (jj = ss; jj < kk; jj++)
            vecCnts[jj] = vecCnts[kk];
        }
        ss = kk + 1;
      }
      //**/ now fill the nonzeros and compute volume
      double minVol = 0, maxVol=0;
      for (ss2 = 0; ss2 < vecResolution[1]; ss2++)
      {
        if (vecCnts[ss2] == 0 && ss2 > 0 && vecCnts[ss2+1] > 0)
        {
          minVol -= vecCnts[ss2+1];
          maxVol += vecCnts[ss2+1];
        }
        else
        {
          minVol += vecCnts[ss2] - 1;
          maxVol += vecCnts[ss2] + 1;
        }
        for (kk = 0; kk < vecResolution[0]-vecCnts[ss2]; kk++)
          meshY[kk][ss2] = thresh;
        for (kk = vecResolution[0]-vecCnts[ss2]; 
             kk < vecResolution[0]; kk++)
          meshY[kk][ss2] = thresh + 1e-8;
      }
      count = vecResolution[0] * vecResolution[1];
      printf("Min percentage above threshold = %7.3f\n",
             100*minVol/count);
      printf("Max percentage above threshold = %7.3f\n",
             100*maxVol/count);
      printf("Based on %d cells.\n", count);
      printEquals(PL_INFO, 0);
      if (outputLevel_ > 1)
      {
        printf("2D Mesh: (# - undefined; X,O - above/below threshold)\n");
        for (ss2 = vecResolution[1]-1; ss2 >= 0; ss2--)
        {
          for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
          {
            if (meshY[ss1][ss2] == PSUADE_UNDEFINED)
              printf("# ");
            else if (meshY[ss1][ss2] > thresh)
              printf("X ");
            else 
              printf("O ");
          }
          printf("\n");
        }
      }

      //**/ create a larger mesh and fill it
      int dim1, dim2;
      psMatrix matMeshY2;
      matMeshY2.setFormat(PS_MAT2D);
      dim1 = 2 * vecResolution[0];
      dim2 = 2 * vecResolution[1];
      matMeshY2.setDim(dim1, dim2);
      double **meshY2 = matMeshY2.getMatrix2D();
      for (ss1 = 0; ss1 < dim1; ss1++) 
        for (ss2 = 0; ss2 < dim2; ss2++) 
          meshY2[ss1][ss2] = PSUADE_UNDEFINED;

      psIMatrix matMesh2Tags;
      matMesh2Tags.setFormat(PS_MAT2D);
      matMesh2Tags.setDim(dim1,dim2);
      int **mesh2Tags = matMesh2Tags.getIMatrix2D();
      vecDX.setLength(nInputs_);
      vecDX[0] = (VecIUpperBs_[0]-VecILowerBs_[0])/dim1;
      vecDX[1] = (VecIUpperBs_[1]-VecILowerBs_[1])/dim2;
      for (ss = 0; ss < nSamples_; ss++)
      {
        ddata = VecSamInputs_[ss*nInputs_];
        ddata -= VecILowerBs_[0];
        ddata /= vecDX[0];
        ind1 = (int) (ddata+1e-15); 
        if (ind1 >= dim1) ind1 = dim1 - 1;
        ddata = VecSamInputs_[ss*nInputs_+1];
        ddata -= VecILowerBs_[1];
        ddata /= vecDX[1];
        ind2 = (int) (ddata+1e-15); 
        if (ind2 >= dim2) ind2 = dim2 - 1;
        meshY2[ind1][ind2] = VecSamOutputs_[ss*nOutputs_+outputID];
        mesh2Tags[ind1][ind2] = mesh2Tags[ind1][ind2] + 1;
      }
      if (outputLevel_ > 2)
      {
        printf("Larger mesh\n");
        for (ss2 = dim2-1; ss2 >= 0; ss2--)
        {
          for (ss1 = 0; ss1 < dim1; ss1++)
            printf("%1d ", mesh2Tags[ss1][ss2]);
          printf("\n");
        }
      }

      //**/ search for boundaries
      int nCand=0;
      double x1, x2;
      count = dim1 * dim2;
      psVector  vecCandX, vecCandY;
      vecCandX.setLength(nInputs_*count);
      vecCandY.setLength(count);
      //**/ first each row, look for cells having values
      //**/ opposite to any of its neighbors 
      for (ss2 = 0; ss2 < vecResolution[1]; ss2++)
      {
        for (ss1 = 0; ss1 < vecResolution[0]; ss1++)
        {
          count = 0;
          if (meshY[ss1][ss2] != PSUADE_UNDEFINED) 
          {
            //**/ west neighbor
            if (ss1 > 0)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1-1][ss2] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1-1][ss2] != PSUADE_UNDEFINED && 
                  meshY[ss1-1][ss2] > thresh) count++;
            }
            //**/ south neighbor
            if (ss2 > 0)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1][ss2-1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1][ss2-1] != PSUADE_UNDEFINED && 
                  meshY[ss1][ss2-1] > thresh) count++;
            }
            //**/ east
            if (ss1 < vecResolution[0]-1)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1+1][ss2] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1+1][ss2] != PSUADE_UNDEFINED && 
                  meshY[ss1+1][ss2] > thresh) count++;
            }
            //**/ north
            if (ss2 < vecResolution[1]-1)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1][ss2+1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1][ss2+1] != PSUADE_UNDEFINED && 
                  meshY[ss1][ss2+1] > thresh) count++;
            }
            //**/ southwest
            if (ss1 > 0 && ss2 > 0)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1-1][ss2-1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1-1][ss2-1] != PSUADE_UNDEFINED && 
                  meshY[ss1-1][ss2-1] > thresh) count++;
            }
            //**/ southeast
            if (ss1 < vecResolution[0]-1 && ss2 > 0)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1+1][ss2-1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1+1][ss2-1] != PSUADE_UNDEFINED && 
                  meshY[ss1+1][ss2-1] > thresh) count++;
            }
            //**/ northwest
            if (ss1 > 0 && ss2 < vecResolution[1]-1)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1-1][ss2+1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1-1][ss2+1] != PSUADE_UNDEFINED && 
                  meshY[ss1-1][ss2+1] > thresh) count++;
            }
            //**/ northeast
            if (ss1 < vecResolution[0]-1 && ss2 < vecResolution[1]-1)
            {
              if (meshY[ss1][ss2] > thresh && 
                  meshY[ss1+1][ss2+1] <= thresh) count++;
              if (meshY[ss1][ss2] <= thresh && 
                  meshY[ss1+1][ss2+1] != PSUADE_UNDEFINED && 
                  meshY[ss1+1][ss2+1] > thresh) count++;
            }
          }
          if (count > 0)
          {
            ind = ss1 * 2;
            kk  = ss2 * 2;
            for (ii = 0; ii < 2; ii++)
            {
              for (jj = 0; jj < 2; jj++)
              {
                if (mesh2Tags[ind+ii][kk+jj] == 0)
                {
                  x1 = vecDX[ii] * (ind + ii) + 0.5 * vecDX[ii];
                  x2 = vecDX[jj] * (kk + jj) + 0.5 * vecDX[jj];
                  vecCandX[nCand*nInputs_] = x1;
                  vecCandX[nCand*nInputs_+1] = x2;
                  vecCandY[nCand] = PSUADE_UNDEFINED;
                  nCand++;
                  mesh2Tags[ind+ii][kk+jj] = -999;
                }
              }
            }
          }
        }
      }

      //**/ compress the candidate set
      count = 1;      
      for (ss = 1; ss < nCand; ss++)
      {
        for (kk = 0; kk < count; kk++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            if (vecCandX[ss*nInputs_+ii] != vecCandX[kk*nInputs_+ii])
              break;
          //**/ if ii == nInputs_ ==> point ss == prior point kk
          if (ii == nInputs_) break;
        }
        //**/ if kk == ss, cannot find same point ==> unique
        if (kk == ss)
        {
          for (ii = 0; ii < nInputs_; ii++)
            vecCandX[count*nInputs_+ii] = vecCandX[ss*nInputs_+ii];
          count++;
        }
      } 
      nCand = count;
      if (outputLevel_ > 2)
      {
        for (ss = 0; ss < nCand; ss++)
        {
          printf("Candidate %d: \n",ss+1);
          for (ii = 0; ii < nInputs_; ii++)
            printf("  X %d = %e\n",ii+1,vecCandX[ss*nInputs_+ii]);
        }
      }
        
      //**/ map into next finer mesh
      if (outputLevel_ > 1)
      {
        printf("2D Large Mesh: (X - new points, 1 - already occupied)\n");
        for (ss2 = dim2-1; ss2 >= 0; ss2--)
        {
          for (ss1 = 0; ss1 < dim1; ss1++)
          {
            if      (mesh2Tags[ss1][ss2] == 0) printf("0 ");
            else if (mesh2Tags[ss1][ss2] == 1) printf("1 ");
            else    printf("X ");
          }
          printf("\n");
        }
      }

      //**/ now vecEvalX and vecEvalY have been filled
      if (nCand == 0) 
      {
        vecCandX.clean();
        printf("INFO: no threshold interface point has been found.\n");
        continue;
      }

      //**/ store the results
      printf("Number of threshold boundary points identified = %d\n",
             nCand);
      fp = fopen("arel_fact2.std", "w");
      if (fp != NULL)
      {
        fprintf(fp,"# Each line: m inputs X, Y, Y(2 * nInputs) neighbors)\n");
        fprintf(fp,"%d %d 1\n", nCand, nInputs_);
        for (ss = 0; ss < nCand; ss++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp,"%12.4e ", vecCandX[ss*nInputs_+ii]);
          fprintf(fp,"%12.4e\n", vecCandY[ss]);
        }
        fclose(fp); 
        printAsterisks(PL_INFO, 0);
        printf("INFO: threshold interface points are in ");
        printf("arel_fact2.std\n");
        printf("      Number of interface points = %d\n",nCand);
        printf("In the next refinement, you must use ");
        printf("finer resolution. Specifically,\n");
        for (ii = 0; ii < nInputs_; ii++)
          printf("   Input %3d resolution = %d\n",ii+1,
                 2*vecResolution[ii]);
        printAsterisks(PL_INFO, 0);
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arefinem_thresh
    //**/ adaptive sample refinement based on METIS
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arefinem_thresh"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arefinem_thresh: adaptively refine a sample ");
        printf("based on a user-specified\n");
        printf("threshold. This function is similar to ");
        printf("arefine1_thresh except that\n");
        printf("this command assumes the loaded sample is a ");
        printf("METIS sample (along with\n");
        printf("the corresponding METIS information file ");
        printf("(psuadeMetisInfo).\n");
        printf("Syntax: arefine_metis (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("Given a METIS sample which has been loaded ");
      printf("(using load or read_std),\n");
      printf("this command adaptively refines the sample ");
      printf("near areas where the RS\n");
      printf("prediction means are closest to (but larger than) ");
      printf("some user-given\n");
      printf("threshold.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data has not been loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,
          "Select output to use for prediction estimation (1 - %d) ",
          nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      psuadeIO_->getParameter("method_sampling", pPtr);
      SamplingMethod_ = pPtr.intData_;
      if (SamplingMethod_ != PSUADE_SAMP_METIS)
      {
        printf("ERROR: adaptive refinement requires METIS sample.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ case when the loaded sample outputs are not 0/1
      printf("This function uses the current loaded sample ");
      printf("together with previous\n");
      printf("refinement information (stored in file ");
      printf("'psuadeMetisInfo') to perform\n");
      printf("adaptive refinement. To work properly, the ");
      printf("following information are\n");
      printf("needed: \n");
      printf("(1) the original sample size (before any refinement)\n");
      printf("    NOTE: if you do not know what this mean, enter %d\n",
             nSamples_);
      printf("(2) the previous refinement file (psuadeMetisInfo)\n");
      printf("    which needs to be in the current directory.\n");
      sprintf(pString,"What is the original sample size ? ");

      //**/ get min/max information from sample output
      double dMin=1e35, dMax=-1e-35;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+outputID] < dMin)
          dMin = VecSamOutputs_[ss*nOutputs_+outputID];
        if (VecSamOutputs_[ss*nOutputs_+outputID] > dMax)
          dMax = VecSamOutputs_[ss*nOutputs_+outputID];
      }
      printf("Sample minimum and maximum are %e and %e.\n",dMin,dMax);
      printf("The threshold should be somewhere in between ");
      printf("or not too far away.\n");
      sprintf(pString, "Please enter threshold : ");
      thresh = getDouble(pString);

      int origNSamples = getInt(1, nSamples_, pString);
      sprintf(pString,
            "How many sample points to add? (1 - %d, default = %d) ",
            nSamples_, nSamples_);
      int refineSize = getInt(1, nSamples_, pString);

      //**/ temporarily turn off sampling expert mode
      //**/ load sample sampler for adaptive sampling (sampPtr)
      psConfig_.SamExpertModeSaveAndReset();
      Sampling *sampPtr = 
            (Sampling *) SamplingCreateFromID(SamplingMethod_);
      sampPtr->setPrintLevel(outputLevel_);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setInputParams(nInputs_, NULL, NULL, NULL);
      sampPtr->setOutputParams(nOutputs_);
      sampPtr->setSamplingParams(nSamples_, -1, -1);
      sampPtr->initialize(1);
      sampPtr->loadSamples(nSamples_,nInputs_,nOutputs_,
                VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
                VecSamStates_.getIVector());

      //**/ create a MARSB response surface
      faType = PSUADE_RS_MARSB;
      FuncApprox *faPtr = genFA(faType, nInputs_, iOne, nSamples_);
      faPtr->setNPtsPerDim(32);
      faPtr->setBounds(VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(0);

      //**/ generate 100 instantiations for MARSBAG
      int numMars = 100, ivar1;
      psMatrix matMarsX, matMarsY;
      matMarsX.setFormat(PS_MAT2D);
      matMarsY.setFormat(PS_MAT2D);
      matMarsX.setDim(numMars, nSamples_*nInputs_);
      matMarsY.setDim(numMars, nSamples_);
      double **marsX = matMarsX.getMatrix2D();
      double **marsY = matMarsY.getMatrix2D();
      for (ii = 0; ii < numMars; ii++)
      {
        for (ss = 0; ss < nSamples_; ss++)
        {
          if (ss < origNSamples)
               ivar1 = PSUADE_rand() % origNSamples;
          else ivar1 = ss;
          for (jj = 0; jj < nInputs_; jj++)
            marsX[ii][ss*nInputs_+jj] =
                     VecSamInputs_[ivar1*nInputs_+jj];
          marsY[ii][ss] = VecSamOutputs_[ivar1*nOutputs_+outputID];
        }
      }
      strcpy(cString, "mars_params");
      int ivar2 = 2 * nInputs_ / 3 + 1;
      int ivar3 = nSamples_;
      if (ivar3 > 300) ivar3 = 300;
      targv[0] = (char *) cString;
      targv[1] = (char *) &ivar3;
      targv[2] = (char *) &ivar2;
      faPtr->setParams(3, targv);
      strcpy(cString, "num_mars");
      targv[0] = (char *) cString;
      targv[1] = (char *) &numMars;
      faPtr->setParams(2, targv);
      strcpy(cString, "mars_sample");
      targv[0] = (char *) cString;
      targv[2] = (char *) &nSamples_;
      for (ii = 0; ii < numMars; ii++)
      {
        targv[1] = (char *) &ii;
        targv[3] = (char *) marsX[ii];
        targv[4] = (char *) marsY[ii];
        faPtr->setParams(5, targv);
      }
      vecYT.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      faPtr->initialize(VecSamInputs_.getDVector(),vecYT.getDVector());
      strcpy(cString, "SETRSPTR");
      targv[0] = (char *) cString;
      targv[1] = (char *) faPtr;
      sampPtr->setParam(2, targv);

      //**/ now refine based on magnitudes of the errors
      strcpy(pString, "setAdaptiveRefinementBasedOnErrors");
      sampPtr->setParam(pString);
      sprintf(pString, "setRefineSize %d", refineSize);
      sampPtr->setParam(pString);
      sampPtr->refine(2,0,0,nSamples_,NULL);
      VecSamOutputs_.clean();
      VecSamStates_.clean();

      //**/ fetch the new sample 
      nSamples_ = sampPtr->getNumSamples();
      VecSamInputs_.setLength(nInputs_*nSamples_);
      VecSamOutputs_.setLength(nOutputs_*nSamples_);
      VecSamStates_.setLength(nSamples_);
      sampPtr->getSamples(nSamples_, nInputs_, nOutputs_, 
                VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(),
                VecSamStates_.getIVector());
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),NULL); 
      psuadeIO_->updateMethodSection(-1, nSamples_, -1, -1, -1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("arefinem_thresh successful: use write to ");
      printf("store the refined sample.\n");
      printf("   Then run the newly created sample points ");
      printf("in this refined sample\n");
      printf("   (the last %d sample points.\n", refineSize);
      delete sampPtr;
      sampPtr = NULL;
      delete faPtr;
      faPtr = NULL;
      psConfig_.SamExpertModeRestore();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ arel_metis - adaptive reliability analysis based on METIS
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "arel_metis"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("arel_metis: adaptive reliability analysis");
        printf("Syntax: arel_metis (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs adaptive sample ");
      printf("refinement for realiability\n");
      printf("analysis. It is intended to be run repeatedly ");
      printf("to achieve a desired\n");
      printf("accuracy.\n");
      printDashes(PL_INFO, 0);
      printf("IMPORTANT: MAKE SURE\n");
      printf("- the loaded sample is a METIS sample which has\n");
      printf("- a corresponding file psuadeMetisInfo in the ");
      printf("in the working directory\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString) - 1;

      //**/ check that all outputs are 0/1 or not
      //**/ if not, transform them to 0/1 ==> vecYT
      for (ss = 0; ss < nSamples_; ss++)
      {
        if ((VecSamOutputs_[ss*nOutputs_+outputID] != 0) && 
            (VecSamOutputs_[ss*nOutputs_+outputID] != 1))
          break;
      }
      double relMin=PSUADE_UNDEFINED, relMax=-PSUADE_UNDEFINED;
      psVector vecYT;
      if (ss != nSamples_)
      {
        printf("Outputs are not 0 or 1 ==> transform to 0 or 1.\n");
        printf("To do so, need lower/upper safety margins.\n");
        while (relMin >= relMax)
        {
          printf("Safety region is defined to be inside [lower, upper].\n");
          sprintf(cString,"Upper bound for safety region (-999 if none): ");
          relMax = getDouble(cString);
          if (relMax == -999) relMax = 1.0e10;
          sprintf(cString,"Lower bound for safety region (-999 if none): ");
          relMin = getDouble(cString);
          if (relMin == -999) relMin = -1.0e10;
          if (relMin >= relMax) printf("INVALID bounds.\n");
        }
        vecYT.setLength(nSamples_);
        for (ss = 0; ss < nSamples_; ss++)
        {
          if (VecSamOutputs_[ss*nOutputs_+outputID] < relMin || 
              VecSamOutputs_[ss*nOutputs_+outputID] > relMax)
               vecYT[ss] = 0;
          else vecYT[ss] = 1;
        }
      }
      else
      {
        vecYT.setLength(nSamples_);
        for (ss = 0; ss < nSamples_; ss++)
          vecYT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      }

      //**/ create metis sampler
      int samMethod = PSUADE_SAMP_METIS;
      Sampling *sampler = SamplingCreateFromID(samMethod);
      sampler->setPrintLevel(outputLevel_);
      sampler->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampler->setOutputParams(iOne);
      sampler->setSamplingParams(nSamples_, 1, 1);
      sampler->initialize(1);
      sampler->loadSamples(nSamples_,nInputs_,iOne,
                  VecSamInputs_.getDVector(),vecYT.getDVector(),
                  VecSamStates_.getIVector());
     
      //**/ set up sampler for refinement
      char sparam[1000];
      strcpy(sparam, "calVolumes");
      int curVol = sampler->setParam(sparam);
      strcpy(sparam, "totalVolumes");
      int totVol = sampler->setParam(sparam);
      double curRel = 100.0 * curVol / totVol;
      printOutTS(PL_INFO,
       "CURRENT RELIABLITY (SUCCESS) = %8.5f%% (%d/%d*100)\n",
       curRel,curVol,totVol);
      //**/ randomize is irrelevant (it tells METIS to randomly
      //**/ select points if new points > refineSize
      //**/ threshold is irrelevant
      int refineRatio=2,randomize=0;
      double refineThresh=0;
      if (curRel == 0 || curRel == 1)
      {
        printf("Current reliability = 0 or 1 ==> uniform refinement.\n");
        strcpy(sparam, "setUniformRefinement");
        sampler->setParam(sparam);
        printf("arel_metis: refinement begins.\n");
        sampler->refine(refineRatio,randomize,refineThresh,0,NULL);
        printf("arel_metis: refinement completed.\n");
      }
      else
      {
        printf("adaptive refinement along 0/1 boundary.\n");
        strcpy(sparam, "setAdaptiveRefinementBasedOnOutputs");
        sampler->setParam(sparam);
        sprintf(sparam, "setRefineSize 100000");
        sampler->setParam(sparam);
        printf("arel_metis: refinement begins.\n");
        sampler->refine(refineRatio,randomize,refineThresh,
                         nSamples_,NULL);
        printf("arel_metis: refinement completed.\n");
      }

      //**/ get refined sample
      int nSamples = sampler->getNumSamples();
      printf("INFO: new nSamples = %d\n", nSamples);
      psVector  vecSamInps, vecSamOuts;
      psIVector vecSamStas;
      vecSamInps.setLength(nSamples*nInputs_);
      vecSamOuts.setLength(nSamples);
      vecSamStas.setLength(nSamples);
      double *sampleInputs  = vecSamInps.getDVector();
      double *sampleOutputs = vecSamOuts.getDVector();
      int    *sampleStates  = vecSamStas.getIVector();
      sampler->getSamples(nSamples,nInputs_,iOne,vecSamInps.getDVector(),
                vecSamOuts.getDVector(), vecSamStas.getIVector());

      //**/ further eliminate new points if certain criteria
      //**/ is satisfied
      if (nInputs_ == 2)
      {
        printf("If your reliability region is at upper right ");
        printf("corner of the input space\n");
        printf("(that is, when both input values are high), ");
        printf("and also the reliability\n");
        printf("region is convex upward, you can further ");
        printf("reduce the new sample points\n");
        printf("by additional analysis.\n");
        printf("Perform additional analysis ? (y or n) ");
        scanf("%s", lineIn2);
        fgets(winput,5000,stdin); 
        count = 0;
        if (lineIn2[0] == 'y')
        {
          //**/ collect all samples with output=0
          count = 0;
          for (ss = 0; ss < nSamples_; ss++)
             if (vecYT[ss] == 0) count++;
        }
        if (count > 0)
        {
          vecWT.setLength(count);
          vecIT.setLength(count);
          vecXT.setLength(count*nInputs_);
          //**/ gather samples with output = 0
          //**/ vecXT, vecIT, vecWT
          kk = 0;
          for (ss = 0; ss < nSamples_; ss++)
          {
            if (vecYT[ss] == 0)
            {  
              for (ii = 0; ii < nInputs_; ii++)
                vecXT[kk*nInputs_+ii] = VecSamInputs_[ss*nInputs_+ii];
              vecWT[kk] = VecSamInputs_[ss*nInputs_+1];
              vecIT[kk] = kk;
              kk++;
            }
          }
          //**/ sort vecWT (second input)
          sortDbleList2a(count,vecWT.getDVector(),
                         vecIT.getIVector());
          //**/ copy the re-ordered inputs back
          for (ss = 0; ss < count; ss++)
            vecXT[ss*nInputs_+1] = vecWT[ss];
          //**/ apply the reordering to the first input
          for (ss = 0; ss < count; ss++)
            vecWT[ss] = vecXT[ss*nInputs_];
          for (ss = 0; ss < count; ss++)
          {
            kk = vecIT[ss];
            vecXT[ss*nInputs_] = vecWT[kk];
          }
          //**/ message vecXT
          for (ss = count-2; ss >= 0; ss--)
          {
            if (vecXT[ss*nInputs_] < vecXT[(ss+1)*nInputs_])
              vecXT[ss*nInputs_] = vecXT[(ss+1)*nInputs_];
          }
          //**/ now vecXT has ordered 2nd inputs corresponding
          //**/ to output = 0
          //**/ now examine each of the new points
          double x1Val, x2Val;
          vecWT.setLength(nSamples-nSamples_);
          for (ss = nSamples_; ss < nSamples; ss++)
          {
            vecWT[ss-nSamples_] = PSUADE_UNDEFINED;
            x1Val = sampleInputs[ss*nInputs_];
            x2Val = sampleInputs[ss*nInputs_+1];
            //**/ if 2nd input is less than the lowest vecXT[1], 
            //**/ and if 1st input less than vecXT[0], then
            //**/ assign 0 to output
            if (x2Val < vecXT[0*nInputs_+1])
            {
              if (x1Val < vecXT[0*nInputs_])
              {
                vecWT[ss-nSamples_] = 0;
              }
            }
            //**/ if 2nd input is larger than the highest vecXT[1], 
            //**/ do nothing (will need to be evaluated) 
            else if (x2Val > vecXT[(count-1)*nInputs_+1])
            {
              vecWT[ss-nSamples_] = PSUADE_UNDEFINED;
            }
            //**/ other search in between 
            else
            {
              for (kk = 1; kk < count; kk++)
              {
                if (x2Val > vecXT[(kk-1)*nInputs_+1] &&
                    x2Val <= vecXT[kk*nInputs_+1])
                {
                  if (x1Val < vecXT[kk*nInputs_])
                  {
                    vecWT[ss-nSamples_] = 0;
                    break;
                  }
                }
              } 
            } 
          }
          for (ss = 0; ss < nSamples-nSamples_; ss++)
          {
            if (vecWT[ss] == 0)
            {
              vecSamOuts[ss+nSamples_] = vecWT[ss];
              vecSamStas[ss+nSamples_] = 1;
              printf("New point %d (1-based) assigned to 0\n",
                     ss+1);
            }
          }
        }
        //**/ continue with Y=1
        if (lineIn2[0] == 'y')
        {
          //**/ collect all samples with output=1
          count = 0;
          for (ss = 0; ss < nSamples_; ss++)
             if (vecYT[ss] == 1) count++;
        }
        if (count > 0)
        {
          vecWT.setLength(count);
          vecIT.setLength(count);
          vecXT.setLength(count*nInputs_);
          //**/ gather samples with output = 1
          //**/ vecXT, vecIT, vecWT
          kk = 0;
          for (ss = 0; ss < nSamples_; ss++)
          {
            if (vecYT[ss] == 1)
            {  
              for (ii = 0; ii < nInputs_; ii++)
                vecXT[kk*nInputs_+ii] = VecSamInputs_[ss*nInputs_+ii];
              vecWT[kk] = VecSamInputs_[ss*nInputs_+1];
              vecIT[kk] = kk;
              kk++;
            }
          }
          //**/ sort vecWT (second input)
          sortDbleList2a(count,vecWT.getDVector(),
                         vecIT.getIVector());
          //**/ copy the re-ordered inputs back
          for (ss = 0; ss < count; ss++)
            vecXT[ss*nInputs_+1] = vecWT[ss];
          //**/ apply the reordering to the first input
          for (ss = 0; ss < count; ss++)
            vecWT[ss] = vecXT[ss*nInputs_];
          for (ss = 0; ss < count; ss++)
          {
            kk = vecIT[ss];
            vecXT[ss*nInputs_] = vecWT[kk];
          }
          //**/ message vecXT
          for (ss = 1; ss < count; ss++)
          {
            if (vecXT[ss*nInputs_] > vecXT[(ss-1)*nInputs_])
              vecXT[ss*nInputs_] = vecXT[(ss-1)*nInputs_];
          }
          //**/ now vecXT has ordered 2nd inputs corresponding
          //**/ to output = 1
          //**/ now examine each of the new points
          double x1Val, x2Val;
          vecWT.setLength(nSamples-nSamples_);
          for (ss = nSamples_; ss < nSamples; ss++)
          {
            vecWT[ss-nSamples_] = PSUADE_UNDEFINED;
            x1Val = sampleInputs[ss*nInputs_];
            x2Val = sampleInputs[ss*nInputs_+1];
            //**/ if 2nd input is less than the lowest vecXT[1], 
            //**/ do nothing
            if (x2Val < vecXT[0*nInputs_+1])
            {
               vecWT[ss-nSamples_] = PSUADE_UNDEFINED;
            }
            //**/ if 2nd input is larger than the highest vecXT[1], 
            //**/ do nothing (will need to be evaluated) 
            else if (x2Val > vecXT[(count-1)*nInputs_+1])
            {
              if (x1Val > vecXT[(count-1)*nInputs_])
                vecWT[ss-nSamples_] = 1;
            }
            //**/ other search in between 
            else
            {
              for (kk = 1; kk < count; kk++)
              {
                if (x2Val > vecXT[(kk-1)*nInputs_+1] &&
                    x2Val <= vecXT[kk*nInputs_+1])
                {
                  //**/ (x2-lb) / (x2ub-x2lb) * (x1ub-x1lb)+x1lb
                  ddata = vecXT[kk*nInputs_+1] -
                          vecXT[(kk-1)*nInputs_+1];
                  dtemp = x2Val - vecXT[(kk-1)*nInputs_+1];
                  ddata = dtemp / ddata; 
                  dtemp = vecXT[kk*nInputs_]-vecXT[(kk-1)*nInputs_];
                  ddata = ddata * dtemp + vecXT[(kk-1)*nInputs_];
                  if (x1Val > ddata)
                  {
                    vecWT[ss-nSamples_] = 1;
                    break;
                  }
                }
              } 
            } 
          }
          for (ss = 0; ss < nSamples-nSamples_; ss++)
          {
            if (vecWT[ss] == 1)
            {
              vecSamOuts[ss+nSamples_] = vecWT[ss];
              vecSamStas[ss+nSamples_] = 1;
              printf("New point %d (1-based) assigned to 1\n",
                     ss+1);
            }
          }
        }
        count = 0;
        for (ss = 0; ss < nSamples; ss++)
        {
          if (vecSamOuts[ss] == PSUADE_UNDEFINED) 
          {
            count++;
            printf("True new point = %d\n", ss+1);
          }
        }
        printf("INFO: new unevaluated nSamples = %d\n", count);

        fp = fopen("arel_sample.m", "w");
        if (fp != NULL)
        {
          fprintf(fp, "A = [\n");
          for (ss = 0; ss < nSamples; ss++)
          {
            for (ii = 0; ii < nInputs_; ii++)
              fprintf(fp,"%16.8e ", vecSamInps[ss*nInputs_+ii]);
            fprintf(fp,"%16.8e\n", vecSamOuts[ss]);
          }
          fprintf(fp,"];\n");
          fwriteHold(fp, 1);
          fprintf(fp,"for ii = 1 : %d\n", nSamples_);
          fprintf(fp,"  plot(A(ii,1),A(ii,2),'b*','markersize',14)\n");
          fprintf(fp,"  if A(ii,%d) == 0\n", nInputs_+1);
          fprintf(fp,"  plot(A(ii,1),A(ii,2),'bo','markersize',14)\n");
          fprintf(fp,"  end;\n");
          fprintf(fp,"end;\n");
          fprintf(fp,"for ii = %d : %d\n", nSamples_+1,nSamples);
          fprintf(fp,"  plot(A(ii,1),A(ii,2),'r*','markersize',24)\n");
          fprintf(fp,"end;\n");
          fprintf(fp,"for ii = %d : %d\n", nSamples_+1,nSamples);
          fprintf(fp,"  if A(ii,%d) == 0\n", nInputs_+1);
          fprintf(fp,"  plot(A(ii,1),A(ii,2),'ro','markersize',24)\n");
          fprintf(fp,"  end;\n");
          fprintf(fp,"  if A(ii,%d) == 1\n", nInputs_+1);
          fprintf(fp,"  plot(X(ii,1),A(ii,2),'rs','markersize',24)\n");
          fprintf(fp,"  end;\n");
          fprintf(fp,"  if A(ii,%d) == 1e35\n", nInputs_+1);
          fprintf(fp,"  plot(A(ii,1),A(ii,2),'kd','markersize',24)\n");
          fprintf(fp,"  end;\n");
          fprintf(fp,"end;\n");
          fprintf(fp,"disp('Blue * - initial sample')\n");
          fprintf(fp,"disp('Blue */blue circle - initial sample Y = 0')\n");
          fprintf(fp,"disp('Red  * - new sample')\n");
          fprintf(fp,"disp('Red */red circle - new sample Y = 0')\n");
          fprintf(fp,"disp('Red */red square - new sample Y = 1')\n");
          fprintf(fp,"disp('Red */black diamond - new sample to be evaluated')\n");
          fwritePlotXLabel(fp, "X1");
          fwritePlotYLabel(fp, "X2");
          fwritePlotTitle(fp, "Input Scatter Plot");
          fwritePlotAxes(fp);
          fclose(fp);
          printf("arel_sample.m is available for inspection.\n");
        }
      }

      //**/ save the appended sample
      PsuadeData *psuadeIO = new PsuadeData();
      psuadeIO->updateInputSection(nSamples,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   vecSamInps.getDVector(),StrInpNames_.getStrings(),
                   NULL,NULL,NULL,NULL);
      psuadeIO->updateOutputSection(nSamples,iOne,vecSamOuts.getDVector(), 
                   vecSamStas.getIVector(),StrOutNames_.getStrings());
      psuadeIO->updateMethodSection(samMethod,nSamples,-1,-1,-1);
      strcpy(winput, "arel_sample.psu");
      psuadeIO->writePsuadeFile(winput,0);
      printf("arel_metis: expanded sample is in arel_sample.psu\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ metis_2dmap - create Matlab plot for 2D partitioning
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "metis_2dmap"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("metis_2dmap: plot 2D METIS parititoned regions\n");
        printf("Syntax: metis_2dmap (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        continue;
      }
      if (nInputs_ > 2)
      {
        printf("ERROR: this command cannot handle nInputs > 2.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a Matlab plot to show ");
      printf("the partitioned regions\n");
      printf("given in METIS-generated partitioning.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      Sampling *samPtr = SamplingCreateFromID(PSUADE_SAMP_METIS);
      samPtr->setPrintLevel(0);
      samPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                             VecIUpperBs_.getDVector());
      samPtr->setOutputParams(1);
      samPtr->setSamplingParams(nSamples_, -1, -1);
      samPtr->loadSamples(nSamples_,nInputs_,nOutputs_,
              VecSamInputs_.getDVector(),VecSamOutputs_.getDVector(), 
              VecSamStates_.getIVector());
      samPtr->initialize(1);
      strcpy(pString,"genMeshPlot");
      samPtr->setParam(pString);
    }

    //**/ -------------------------------------------------------------
    // +++ sqc
    //**/ sample quality check 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sqc"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sqc: check sample quality\n");
        printf("Syntax: sqc (no argument needed)\n");
        printf("The sample quality is measured by min-min distance\n");
        printf("between sample points and with corners.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command computes the following measures:\n");
      printf("* Min-min = minimum of minimum distances between all pairs\n");
      printf("* Avg-min = average of minimum distances between all pairs\n");
      printf("* Min-cor = avg of min dist between all points and corners.\n");
      printf("Since each input may be in different scales, if you desire ");
      printf("to give the\n");
      printf("same weights to all inputs with respect to their ranges, ");
      printf("select yes\n");
      printf("below.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      psVector VecTV;
      VecTV.setLength(3);
      winput[0] = 'a';
      while (winput[0] != 'n' && winput[0] != 'y')
      {
        printf("Scale the inputs (w.r.t. ranges)? (y or n) ");
        scanf("%s", winput);
      }
      psVector vecLs, vecHs;
      if (winput[0] == 'n') 
      {
        vecLs.setLength(nInputs_);
        vecHs.setLength(nInputs_);
        vecHs.setConstant(1.0);
      }
      else
      {
        vecLs = VecILowerBs_;
        vecHs = VecIUpperBs_;
      }
      SamplingQuality(nSamples_,nInputs_,VecSamInputs_.getDVector(),
          vecLs.getDVector(),vecHs.getDVector(),VecTV.getDVector(),0);
      printf("Min-min distance (sample-sample ) = %10.4e (large is good)\n",
             VecTV[0]);
      printf("Avg-min distance (sample pairs  ) = %10.4e (large is good)\n",
             VecTV[2]);
      if (VecTV[1] >= 0)
        printf("Min-cor distance (corner-samples) = %10.4e (small is good)\n",
               VecTV[1]);
      cmdStatus = 0;
      fgets(winput,5000,stdin); 
    }

    //**/ -------------------------------------------------------------
    // +++ ssc
    //**/ sample smoothness check 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ssc"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ssc: check sample smoothness\n");
        printf("Syntax: ssc (no argument needed)\n");
        printf("The sample smoothness is measured by taking the mean and\n");
        printf("the standard deviation of the differences (normalized by\n");
        printf("the distances) between the output of a sample point with\n");
        printf("its nearest neighbors.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
      }
      else
      {
        sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
        outputID = getInt(1, nOutputs_, pString);
        outputID--;
        sprintf(pString,
          "Check (a) against itself or (b) against another sample? (a or b) ");
        getString(pString, winput);

        if (winput[0] == 'a')  
        {      
          int nNeighbors=0;
          sprintf(pString,
             "How many nearest neighbors to use number (1 - default=100) : ");
          nNeighbors = getInt(1, 100, pString);

          psVector VecX, VecY, VecXL, VecXU, VecM, VecSD;
          VecX.load(nSamples_*nInputs_, VecSamInputs_.getDVector());
          VecY.setLength(nSamples_);
          for (ss = 0; ss < nSamples_; ss++)
            VecY[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
          VecXL.load(nInputs_, VecILowerBs_.getDVector());
          VecXU.load(nInputs_, VecIUpperBs_.getDVector());
          CheckSampleSmoothness(VecX,VecY,VecXL,VecXU,nNeighbors,VecM,
                                VecSD,outputLevel_);
          if (plotScilab()) fp = fopen("scilabssc.sci", "w");
          else              fp = fopen("matlabssc.m", "w");
          fprintf(fp,"A = [\n"); 
          for (ss = 0; ss < nSamples_; ss++) fprintf(fp,"%16.8e\n",VecM[ss]); 
          fprintf(fp,"];\n"); 
          fprintf(fp,"D = [\n"); 
          for (ss = 0; ss < nSamples_; ss++) fprintf(fp,"%16.8e\n",VecSD[ss]); 
          fprintf(fp,"];\n"); 
          fprintf(fp, "clf\n");
          fprintf(fp, "plot(A,'*','markersize',8)\n");
          fwritePlotXLabel(fp, "Sample number");
          fwritePlotYLabel(fp, "Avg Gradient Measure");
          fwritePlotTitle(fp, "Sample Smoothness Analysis");
          fwritePlotAxes(fp);
          fprintf(fp, "figure(2)\n");
          fprintf(fp, "plot(D,'*','markersize',8)\n");
          fwritePlotXLabel(fp, "Sample number");
          fwritePlotYLabel(fp, "Shortest Normalized (wrt range) Distance");
          fwritePlotTitle(fp, "Sample Smoothness Analysis 2");
          fwritePlotAxes(fp);
          fclose(fp);
          if (plotScilab()) 
            printf("INFO: sample smoothness plot is now in scilabssc.sci\n");
          else
            printf("INFO: sample smoothness plot is now in matlabssc.m\n");
        }
        else
        {
          int nNeighbors=0;
          sprintf(pString,
             "How many neighbors to use number (1 - default=100) : ");
          nNeighbors = getInt(1, 100, pString);
          sprintf(pString,"Enter name of sample file name : ");
          getString(pString, winput);
          kk = strlen(winput);
          winput[kk-1] = '\0';
          PsuadeData *ioPtr = new PsuadeData;
          status = ioPtr->readPsuadeFile(winput);
          if (status != 0)
          {
            printf("ssc FILE READ ERROR: file = %s\n", winput);
            continue;
          }
          ioPtr->getParameter("input_ninputs", pPtr);
          int ssc_ninputs = pPtr.intData_;
          if (ssc_ninputs != nInputs_)
          {
            printf("ssc ERROR: nInputs mismatch.\n");
            printf("         local nInputs = %d.\n", nInputs_);
            printf("         file  nInputs = %d.\n", jj);
            cmdStatus = 1;
            continue;
          }
          ioPtr->getParameter("method_nsamples", pPtr);
          int ssc_nsamples = pPtr.intData_;
          ioPtr->getParameter("output_noutputs", pPtr);
          int ssc_noutputs = pPtr.intData_;
          if (ssc_noutputs != 1)
          {
            printf("ssc ERROR: sample file nOutputs should be = 1\n");
            cmdStatus = 1;
            continue;
          }

          psVector VecX, VecY, VecXL, VecXU, VecM, VecSD;
          VecX.load(nSamples_*nInputs_, VecSamInputs_.getDVector());
          VecY.setLength(nSamples_);
          for (ss = 0; ss < nSamples_; ss++)
            VecY[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
          VecXL.load(nInputs_, VecILowerBs_.getDVector());
          VecXU.load(nInputs_, VecIUpperBs_.getDVector());

          psVector VecX2, VecY2;
          pData    sscPX, sscPY;
          ioPtr->getParameter("input_sample", sscPX);
          VecX2.load(ssc_ninputs * ssc_nsamples, sscPX.dbleArray_);
          ioPtr->getParameter("output_sample", sscPY);
          VecY2.load(ssc_nsamples, sscPY.dbleArray_);
          CheckSampleSmoothness2(VecX,VecY,VecXL,VecXU,VecX2,VecY2,nNeighbors,
                                 VecM,VecSD,outputLevel_);
          if (plotScilab()) fp = fopen("scilabssc.sci", "w");
          else              fp = fopen("matlabssc.m", "w");
          fprintf(fp,"A = [\n"); 
          for (ss = 0; ss < ssc_nsamples; ss++) 
            fprintf(fp,"%16.8e\n", VecM[ss]); 
          fprintf(fp,"];\n"); 
          fprintf(fp,"D = [\n"); 
          for (ss = 0; ss < ssc_nsamples; ss++) 
            fprintf(fp,"%16.8e\n", VecSD[ss]); 
          fprintf(fp,"];\n"); 
          fprintf(fp, "clf\n");
          fprintf(fp, "plot(A,'*','markersize',8)\n");
          fwritePlotXLabel(fp, "Sample number");
          fwritePlotYLabel(fp, "Avg Gradient Measure");
          fwritePlotTitle(fp, "Sample Smoothness Analysis");
          fwritePlotAxes(fp);
          fprintf(fp, "figure(2)\n");
          fprintf(fp, "plot(D,'*','markersize',8)\n");
          fwritePlotXLabel(fp, "Sample number");
          fwritePlotYLabel(fp, "Shortest Normalized Distance");
          fwritePlotTitle(fp, "Sample Smoothness Analysis 2");
          fclose(fp);
          if (plotScilab()) 
            printf("INFO: sample smoothness plot is now in scilabssc.sci\n");
          else
            printf("INFO: sample smoothness plot is now in matlabssc.m\n");
          delete ioPtr;
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ svalidate 
    //**/ validate certain sample outputs (a range)
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "svalidate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("svalidate: validate a subset of sample points (outputs)\n");
        printf("Syntax: svalidate (no argument needed)\n");
        printf("This command sets the ready flag of the selected  sample\n");
        printf("points to be 'valid' (NO OUTPUT VALUE IS MODIFIED).\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command sets the 'ready' flags of the selected ");
      printf("sample points.\n");
      printf("To do this, you need to specify the range of points ");
      printf("to be validated\n");
      printf("(that is, the beginning and ending sample points).\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", winput);
      fgets(lineIn2, 5000, stdin);
      if (winput[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString,"Enter first sample number (1 - %d) = ", nSamples_);
      ind = getInt(1, nSamples_, pString);
      ind--;
      sprintf(pString, "Enter last sample number (%d - %d) = ", 
              ind+1, nSamples_);
      ind2 = getInt(ind+1, nSamples_, pString);
      for (ii = ind; ii < ind2; ii++) VecSamStates_[ii] = 1;
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(),NULL); 
      psuadeIO_->writePsuadeFile(dataFile,0);
      cmdStatus = 0;
      printf("svalidate completed. Use 'write' to store the new sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ sinvalidate 
    //**/ invalidate all sample outputs (set to undefined)
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sinvalidate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sinvalidate: tag selected sample points to be ");
        printf("'un-ready' based on\n");
        printf("             some filtering criterion.\n");
        printf("Syntax: sinvalidate (no argument needed)\n");
        printf("This command invalidates sample points based ");
        printf("on the values of an\n");
        printf("output. If most sample points are to be ");
        printf("invalidated, we suggest\n");
        printf("you first use 'invalidate' to invalidate ");
        printf("all sample points, and\n");
        printf("then use 'validate' to restore the range of ");
        printf("desired sample points.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command invalidates sample points based on ");
      printf("sample output values.\n");
      printf("If all sample points are to be invalidated, enter 0 below.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d, 0 for all points) : ", 
              nOutputs_);
      outputID = getInt(0, nOutputs_, pString);
      outputID--;
      if (outputID == -1)
      {
        for (ii = 0; ii < nSamples_; ii++) VecSamStates_[ii] = 0;
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                     VecSamOutputs_.getDVector(), 
                     VecSamStates_.getIVector(),NULL); 
        printf("Use the write command to store the filtered sample.\n");
      }
      else
      {
        Ymax = - PSUADE_UNDEFINED;
        Ymin =   PSUADE_UNDEFINED;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
            Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
          if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
            Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
        }
        printf("INFO: Values outside the bounds are invalidated.\n");
        sprintf(pString, "Enter the filter's lower bound (Ymin=%e) : ",Ymin);
        threshL = getDouble(pString);
        sprintf(pString, "Enter the filter's upper bound (Ymax=%e) : ",Ymax);
        threshU = getDouble(pString);
        if (threshL >= threshU)
        {
          printf("ERROR: lower bound >= upper bound.\n");
          cmdStatus = 1;
          continue;
        }
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_+outputID] < threshL ||
              VecSamOutputs_[sInd*nOutputs_+outputID] > threshU)
            VecSamStates_[sInd] = 0;
        }
        psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                     VecSamOutputs_.getDVector(),
                     VecSamStates_.getIVector(), NULL);
        if (currSession != NULL) delete currSession;
        currSession = new PsuadeSession();
        psuadeIO_->getSession(currSession);
        printf("INFO: use 'spurge' to take out the invalid points.\n");
        printf("INFO: then use 'write' to store the filtered sample.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ srandomize 
    //**/ randomize the orders of the resident sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "srandomize"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("srandomize: randomize the ordering the local sample.\n");
        printf("Syntax: srandomize <n>\n");
        printf("That command can be used, for example, with rstest_cv.\n");
        printf("(that is, re-order and perform RS cross validation)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command randomly shuffles (reorders) the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      psVector  VecXT, VecYT;
      psIVector VecST, VecITmp;
      VecXT.load(nInputs_*nSamples_, VecSamInputs_.getDVector());
      VecYT.load(nOutputs_*nSamples_, VecSamOutputs_.getDVector());
      VecST.load(nSamples_, VecSamStates_.getIVector());
      VecITmp.setLength(nSamples_);
      generateRandomIvector(nSamples_, VecITmp.getIVector());
      for (ii = 0; ii < nSamples_; ii++)
      {
        ind = VecITmp[ii];
        for (jj = 0; jj < nInputs_; jj++)
          VecSamInputs_[ii*nInputs_+jj] = VecXT[ind*nInputs_+jj];
        for (jj = 0; jj < nOutputs_; jj++)
          VecSamOutputs_[ii*nOutputs_+jj] = VecYT[ind*nOutputs_+jj];
        VecSamStates_[ii] = VecST[ind];
      }
      cmdStatus = 0;
      printf("srandomize completed. Use 'write' to store the new sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ spurge 
    //**/ take out failed sample points
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "spurge"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("spurge: remove bad sample points (Output=UNDEFINED or\n");
        printf("        status != 1). There is also an option to remove\n");
        printf("        repeated sample points.)\n");
        printf("Syntax: spurge (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command removes `invalid' and  duplicate points ");
      printf("from the loaded\n");
      printf("sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,
              "Do you want duplicate sample points purged too? (y or n)");
      getString(pString, winput);
      int setCompare = 0;
      if (winput[0] == 'y') setCompare = 1;
      kk = 0;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        for (oInd = 0; oInd < nOutputs_; oInd++)
          if (VecSamOutputs_[sInd*nOutputs_+oInd] == PSUADE_UNDEFINED || 
              VecSamStates_[sInd] == 0) break; 
        if (oInd == nOutputs_)
        {
          if (setCompare == 1)
            jj = compareSamples(sInd,nSamples_,nInputs_,
                    VecSamInputs_.getDVector(),VecSamStates_.getIVector());
          else jj = -1;
          if (jj < 0 || jj > sInd)
          {
            for (iInd = 0; iInd < nInputs_; iInd++)
              VecSamInputs_[kk*nInputs_+iInd] = 
                   VecSamInputs_[sInd*nInputs_+iInd];
            for (oInd = 0; oInd < nOutputs_; oInd++)
              VecSamOutputs_[kk*nOutputs_+oInd] = 
                   VecSamOutputs_[sInd*nOutputs_+oInd];
            VecSamStates_[kk] = VecSamStates_[sInd]; 
            kk++;
          }
          else
          {
            printf("Repeated sample points: %d and %d\n",sInd+1,jj+1);
            for (iInd = 0; iInd < nInputs_; iInd++)
              printf("Inputs %4d = %12.4e %12.4e\n",iInd+1,
                     VecSamInputs_[sInd*nInputs_+iInd],
                     VecSamInputs_[jj*nInputs_+iInd]); 
            for (oInd = 0; oInd < nOutputs_; oInd++)
              printf("Outputs %3d = %12.4e %12.4e\n",oInd+1,
                     VecSamOutputs_[sInd*nOutputs_+oInd],
                     VecSamOutputs_[jj*nOutputs_+oInd]); 
          }
        }
      }
      nSamples_ = kk;
      VecSamInputs_.subvector(0, nSamples_*nInputs_-1);
      VecSamOutputs_.subvector(0, nSamples_*nOutputs_-1);
      VecSamStates_.subvector(0, nSamples_-1);
      SamplingMethod_ = PSUADE_SAMP_MC;
      nReplications_ = 1;
      psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,
                                     nReplications_,0,-1);
      psuadeIO_->updateInputSection(nSamples_, nInputs_, NULL, 
                   VecILowerBs_.getDVector(), VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL);
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("Number of sample points after spurge = %d\n", nSamples_);
      printf("spurge completed. Use write to store the reduced sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rm_dup
    //**/ remove duplicate points 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rm_dup"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rm_dup: take out or combine duplicate sample points\n");
        printf("Syntax: rm_dup (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("The loaded sample may have duplicate sample points ");
      printf("with identical\n");
      printf("sample inputs (which may be due to nondeterministic ");
      printf("simulations).\n");
      printf("Since some PSUADE methods (e.g. GP response surface) ");
      printf("do not allow\n");
      printf("duplicates, the repeated sample points may need to ");
      printf("be compressed.\n");
      printf("This command offers several schemes for sample point ");
      printf("compression.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (!checkRepeatedSamplePts(nSamples_,nInputs_,
                                  VecSamInputs_.getDVector()))
      {
        printf("INFO: the sample does not have duplicate ");
        printf("sample points. Abort.\n");
        continue;
      }

      collapseSample();
      SamplingMethod_ = PSUADE_SAMP_MC;
      nReplications_ = 1;
      psuadeIO_->updateMethodSection(SamplingMethod_,nSamples_,
                                     nReplications_,0,-1);
      psuadeIO_->updateInputSection(nSamples_, nInputs_, NULL, 
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(), NULL,NULL,NULL,NULL,NULL);
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("Number of sample points after rm_dup = %d\n", nSamples_);
      printf("rm_dup completed. Use write to store the reduced sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sopt_mmd
    //**/ optimize sample using coordinate exchange
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sopt_mmd"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sopt_mmd: permute the resident sample inputs to maximize\n");
        printf("     the inter-sample distances using the coordinate\n");
        printf("     exchange algorithm, which is especially suitable\n");
        printf("     for Latin hypercube samples.\n");
        printf("Syntax: sopt_mmd (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command uses a coordinate exchange algorithm ");
      printf("to search for a\n");
      printf("permutation (of inputs) of the resident sample that ");
      printf("maximizes the\n");
      printf("minimum distance between any 2 sample points. As such, ");
      printf("it is ideal for\n");
      printf("optimizing Latin hypercube samples without destroying ");
      printf("its property.\n");
      printDashes(PL_INFO, 0);
      if (nSamples_ > 1000)
        printf("INFO: sample size > 1000 ==> may take a very long time.\n");
      printf("Proceed? (y or n) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ put inputs and outputs away
      winput[0] = '0';
      while (winput[0] != 'n' && winput[0] != 'y')
      {
        printf("Use scaled (w.r.t. input ranges) distances? (y or n) ");
        scanf("%s", winput);
      }
      if (winput[0] == 'y')
      {
        for (ii = 0; ii < nSamples_; ii++)
        {
          for (kk = 0; kk < nInputs_; kk++)
            VecSamInputs_[ii*nInputs_+kk] /= 
                 (VecIUpperBs_[kk] - VecILowerBs_[kk]);
        }
      }
      
      //**/ perform pre-optimization quality check 
      psVector vecLs, vecHs, VecTV;
      VecTV.setLength(3);
      vecLs.setLength(nInputs_);
      vecHs.setLength(nInputs_);
      vecHs.setConstant(1.0);
      SamplingQuality(nSamples_,nInputs_,VecSamInputs_.getDVector(),
          vecLs.getDVector(),vecHs.getDVector(),VecTV.getDVector(),0);
      printEquals(PL_INFO, 0);
      printf("Sampling quality result before applying optimization: \n");
      printf("Min-min distance (sample-sample ) = %10.4e (large is good)\n",
             VecTV[0]);
      printf("Avg-min distance (sample pairs  ) = %10.4e (large is good)\n",
             VecTV[2]);
      //**/ optimize and invalidate sample outputs
      //**/ Use twice to make sure near-optimality is achieved.
      //**/ This is so because the algorithm accepts exchanges
      //**/ that gives the same minimum distance (rather than
      //**/ strictly larger) to allow searching alternate paths.
      optimizeSampleUsingCE(nSamples_,nInputs_,VecSamInputs_);
      optimizeSampleUsingCE(nSamples_,nInputs_,VecSamInputs_);
 
      //**/ post-optimization sampling quality check
      printEquals(PL_INFO, 0);
      SamplingQuality(nSamples_,nInputs_,VecSamInputs_.getDVector(),
          vecLs.getDVector(), vecHs.getDVector(), VecTV.getDVector(),0);
      printf("Sampling quality result after applying optimization: \n");
      printf("Min-min distance (sample-sample ) = %10.4e (large is good)\n",
             VecTV[0]);
      printf("Avg-min distance (sample pairs  ) = %10.4e (large is good)\n",
             VecTV[2]);
      printEquals(PL_INFO, 0);
      for (ii = 0; ii < nSamples_; ii++)
      {
        if (winput[0] == 'y')
        {
          for (kk = 0; kk < nInputs_; kk++)
            VecSamInputs_[ii*nInputs_+kk] *= 
                 (VecIUpperBs_[kk] - VecILowerBs_[kk]);
        }
        for (kk = 0; kk < nOutputs_; kk++)
          VecSamOutputs_[ii*nOutputs_+kk] = PSUADE_UNDEFINED;
        VecSamStates_[ii] = 0;
      }
      printf("sopt_mmd INFO: Optimized sample is in local memory.\n");
      printf("               Use write to store it.\n");
      printf("         NOTE: Since the order of each input has ");
      printf("been changed, the\n");
      printf("               sample outputs are no longer valid. ");
      printf("Hence, they have\n");
      printf("               been set to UNDEFINED.\n");
      fgets(winput,5000,stdin); 
      cmdStatus = 0;
      continue;
    }
        
    //**/ -------------------------------------------------------------
    // +++ convhull2d 
    //**/ find if a point is inside a convex hull 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "convhull2d"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("convhull2d: find whether a given point is inside ");
        printf("the convex hull \n");
        printf("            formed by 2 inputs in the loaded sample.\n");
        printf("Syntax: convhull2d (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command finds whether a given point is ");
      printf("inside the convex hull\n");
      printf("formed by 2 inputs in the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      double xVal, yVal;
      printf("Since this command only works for 2D, you ");
      printf("need to select 2 inputs from\n");
      printf("the %d inputs in the loaded sample.\n",nInputs_);
      sprintf(pString,"Select the first input (1 - %d) : ",nInputs_);
      ii = getInt(1,nInputs_,pString);
      jj = ii;
      sprintf(pString,"Select the second input (1 - %d, not %d) : ",
              nInputs_,ii);
      while (jj <= 0 || jj > nInputs_ || jj == ii)
        jj = getInt(1,nInputs_,pString);
      printf("Now enter the point to see if it is in the convex hull.\n");
      sprintf(pString,"Enter the value of the first input : ");
      xVal = getDouble(pString);
      sprintf(pString,"Enter the value of the second input : ");
      yVal = getDouble(pString);
      psVector vec1, vec2;
      vec1.setLength(nSamples_);
      vec2.setLength(nSamples_);
      for (kk = 0; kk < nSamples_; kk++)
      {
        vec1[kk] = VecSamInputs_[kk*nInputs_+ii-1];
        vec2[kk] = VecSamInputs_[kk*nInputs_+jj-1];
      }
      status = checkConvHull2D(xVal,yVal,vec1,vec2);
      if (status != 0)
      {
        printf("The point is either inside or on the ");
        printf("boundary of the convex hull.\n");
      }
      else
      {
        printf("The point is outside the convex hull.");
      }
    }

    //**/ -------------------------------------------------------------
    // +++ imodify 
    //**/ modify an input of a sample paoint
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "imodify"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("imodify: modify an input of a selected sample point\n");
        printf("Syntax: imodify (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command modifies one input of a selected sample point.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter sample number (1 - %d) : ", nSamples_);
      sInd = getInt(1, nSamples_, pString);
      sInd--;
      sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
      inputID = getInt(1, nInputs_, pString);
      inputID--;
      printf("Current value of input = %e\n",
             VecSamInputs_[sInd*nInputs_+inputID]);
      sprintf(pString, "Enter new value : ");
      VecSamInputs_[sInd*nInputs_+inputID] = getDouble(pString);
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                    VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ omodify 
    //**/ modify an output of a sample paoint
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "omodify"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("omodify: modify an output of a selected sample point\n");
        printf("Syntax: omodify (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command modifies one output of a selected sample point.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter sample number (1 - %d) : ", nSamples_);
      sInd = getInt(1, nSamples_, pString);
      sInd--;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      ii = getInt(1, nOutputs_, pString);
      ii--;
      printf("Current value of output = %e\n",
             VecSamOutputs_[sInd*nOutputs_+ii]);
      sprintf(pString, "Enter new value : ");
      VecSamOutputs_[sInd*nOutputs_+ii] = getDouble(pString);
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), NULL);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ifilter 
    //**/ take out infeasible sample points
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ifilter"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ifilter: take out sample points based on constraints\n");
        printf("Syntax: ifilter (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command removes some of the loaded sample ");
      printf("points based on the\n");
      printf("values of a selected input.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      inputID = 0;
      if (nInputs_ > 1)
      {
        sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
        inputID = getInt(1, nInputs_, pString);
        inputID--;
      }
      Xmin = VecSamInputs_[inputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamInputs_[sInd*nInputs_+inputID] < Xmin) 
          Xmin = VecSamInputs_[sInd*nInputs_+inputID];
      Xmax = VecSamInputs_[inputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamInputs_[sInd*nInputs_+inputID] > Xmax) 
          Xmax = VecSamInputs_[sInd*nInputs_+inputID];
      printf("Xmin and Xmax found = %e %e.\n", Xmin, Xmax);
      sprintf(pString,"Enter the lower threshold (Xmin = %e) : ",Xmin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper threshold (Xmax = %e) : ",Xmax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: Lower bound should be < upper bound.\n");
        continue;
      }
      kk = 0;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamInputs_[sInd*nInputs_+inputID] >= threshL && 
            VecSamInputs_[sInd*nInputs_+inputID] <= threshU) 
        {
          for (iInd = 0; iInd < nInputs_; iInd++)
            VecSamInputs_[kk*nInputs_+iInd] = 
                 VecSamInputs_[sInd*nInputs_+iInd];
          for (oInd = 0; oInd < nOutputs_; oInd++)
            VecSamOutputs_[kk*nOutputs_+oInd] = 
                 VecSamOutputs_[sInd*nOutputs_+oInd];
          VecSamStates_[kk] = VecSamStates_[sInd]; 
          kk++;
        }
      }
      nSamples_ = kk;
      Xmin =   PSUADE_UNDEFINED;
      Xmax = - PSUADE_UNDEFINED;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamInputs_[sInd*nInputs_+inputID] < Xmin)
          Xmin = VecSamInputs_[sInd*nInputs_+inputID];
        if (VecSamInputs_[sInd*nInputs_+inputID] > Xmax)
          Xmax = VecSamInputs_[sInd*nInputs_+inputID];
      }
      VecILowerBs_[inputID] = Xmin;
      VecIUpperBs_[inputID] = Xmax;
      nReplications_ = 1;
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                   NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      psuadeIO_->updateMethodSection(-1,nSamples_,nReplications_,0,-1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ifilter completed. Use write and load again to continue.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ofilter 
    //**/ take out sample points outside bounds
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ofilter"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ofilter: take out sample points based on constraints\n");
        printf("Syntax: ofilter (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command removes some of the loaded sample ");
      printf("points based on the\n");
      printf("values of a selected output.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      Ymax = - PSUADE_UNDEFINED;
      Ymin =   PSUADE_UNDEFINED;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
          Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
          Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
      }
      sprintf(pString,"Enter the lower threshold (Ymin=%e) : ",Ymin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper threshold (Ymax=%e) : ",Ymax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        continue;
      }
      kk = 0;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] >= threshL &&
            VecSamOutputs_[sInd*nOutputs_+outputID] <= threshU)
        {
          for (iInd = 0; iInd < nInputs_; iInd++)
            VecSamInputs_[kk*nInputs_+iInd] = 
                    VecSamInputs_[sInd*nInputs_+iInd];
          for (oInd = 0; oInd < nOutputs_; oInd++)
            VecSamOutputs_[kk*nOutputs_+oInd] = 
                    VecSamOutputs_[sInd*nOutputs_+oInd];
          VecSamStates_[kk] = VecSamStates_[sInd]; 
          kk++;
        }
      }
      nSamples_ = kk;
      nReplications_ = 1;
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,NULL,NULL,
                 VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                 NULL,NULL,NULL,NULL); 
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      psuadeIO_->updateMethodSection(-1,nSamples_,nReplications_,0,-1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ofilter completed. Use write and load again to continue.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iop 
    //**/ manipulate inputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iop"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iop: perform various operations on the inputs.\n");
        printf("Syntax: iop (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command modifies a sample input by ");
      printf("replacing it with its inverse,\n");
      printf("trigonometric transformation, logarithmic ");
      printf("transformation, etc.\n");
      printf("NOTE: to add a new input with transformation of existing ");
      printf("inputs, use\n");
      printf("      iadd2 to add an input and then run this command.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("Select from the following transformations : \n");
      printf("1. Linear transformation : \n");
      printf("   input <1> = <a> * input <2> + <b> * input <3>\n");
      printf("2. Take an input to a power: input <1> = (input <1>)^P\n");
      printf("3. Exponentiate an input : input <1> = exp(input <1>)\n");
      printf("4. Replace an input with its logarithm : ");
      printf("input <1> = log(input <1>)\n");
      printf("5. Replace an input with its sine : ");
      printf("input <1> = sin(input <1>)\n");
      printf("6. Replace an input with its cosine : ");
      printf("input <1> = cos(input <1>)\n");
      printf("7. Replace an input with its tanget : ");
      printf("input <1> = tan(input <1>)\n");
      printf("8. Replace an input with its arctanget : ");
      printf("input <1> = arctan(input <1>)\n");
      printf("9. Divide one input by another : ");
      printf("input <1> = input <2> / input <3>)\n");
      sprintf(pString, "Enter a selection : ");
      int option = getInt(1, 9, pString);

      if (option == 1)
      {
        printf("Form linear transformation of inputs : \n");
        printf("   input <1> = <a> * input <2> + <b> * input <3>\n");
        sprintf(pString, "Enter input <1> (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        sprintf(pString, "Enter the value <a>: ");
        double aVal = getDouble(pString);
        sprintf(pString, "Enter input <2> (1 - %d) : ", nInputs_);
        jj = getInt(1, nInputs_, pString);
        jj--;
        sprintf(pString, "Enter the value <b>: ");
        double bVal = getDouble(pString);
        sprintf(pString, "Enter input <3> (1 - %d) : ", nInputs_);
        kk = getInt(1, nInputs_, pString);
        kk--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                        aVal * VecSamInputs_[sInd*nInputs_+jj] + 
                        bVal * VecSamInputs_[sInd*nInputs_+kk]; 
        }
      }
      else if (option == 2)
      {
        printf("Take an input to a power: input = (input)^P\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        sprintf(pString, "Enter the power <P>: ");
        double pVal = getDouble(pString);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
            pow(VecSamInputs_[sInd*nInputs_+ii],pVal);
        }
      }
      else if (option == 3)
      {
        printf("Exponentiate an input : input = exp(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                exp(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 4)
      {
        printf("Replace an input with its logarithm : ");
        printf("input = log(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamInputs_[sInd*nInputs_+ii] <= 0)
            printf("ERROR: sample %d input %d <= 0\n",sInd+1,ii+1); 
          else
            VecSamInputs_[sInd*nInputs_+ii] = 
                log(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 5)
      {
        printf("Replace an input with its sine : ");
        printf("input = sin(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                sin(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 6)
      {
        printf("Replace an input with its cosine : ");
        printf("input  = cos(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                cos(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 7)
      {
        printf("Replace an input with its tangent : ");
        printf("input = tan(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                tan(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 8)
      {
        printf("Replace an input with its arctangent : ");
        printf("input = arctan(input)\n");
        sprintf(pString, "Enter input (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamInputs_[sInd*nInputs_+ii] = 
                atan(VecSamInputs_[sInd*nInputs_+ii]); 
        }
      }
      else if (option == 9)
      {
        printf("Division: input <1> = input <2> / input <3>\n");
        sprintf(pString, "Enter input <1> (1 - %d) : ", nInputs_);
        ii = getInt(1, nInputs_, pString);
        ii--;
        sprintf(pString, "Enter input <2> (1 - %d) : ", nInputs_);
        jj = getInt(1, nInputs_, pString);
        jj--;
        sprintf(pString, "Enter input <3> (1 - %d) : ", nInputs_);
        kk = getInt(1, nInputs_, pString);
        kk--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamInputs_[sInd*nInputs_+kk] == 0)
            printf("ERROR: sample %d input %d = 0\n",sInd+1,kk+1); 
          else
            VecSamInputs_[sInd*nInputs_+ii] = 
                        VecSamInputs_[sInd*nInputs_+jj] / 
                        VecSamInputs_[sInd*nInputs_+kk]; 
        }
      }
      if (VecILowerBs_.length() > 0)
      {
        VecILowerBs_[ii] = PSUADE_UNDEFINED;
        VecIUpperBs_[ii] = -PSUADE_UNDEFINED;
        VecInpPDFs_[ii] = 0;
        VecInpMeans_[ii] = 0;
        VecInpStds_[ii] = 0;
        for (ss = 0; ss < nSamples_; ss++)
        {
          if (VecSamInputs_[ss*nInputs_+ii] < VecILowerBs_[ii])
            VecILowerBs_[ii] = VecSamInputs_[ss*nInputs_+ii];
          if (VecSamInputs_[ss*nInputs_+ii] > VecIUpperBs_[ii])
            VecIUpperBs_[ii] = VecSamInputs_[ss*nInputs_+ii];
        }
      }
      if (inputCMat_ != NULL)
      {
        for (jj = 0; jj < nInputs_; jj++)
        {
          inputCMat_->setEntry(jj,ii,0.0);
          inputCMat_->setEntry(ii,jj,0.0);
        }
        inputCMat_->setEntry(ii,ii,1.0);
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("iop completed. Use write and load again to continue.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oop 
    //**/ manipulate outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oop"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oop: perform various operations on the outputs.\n");
        printf("Syntax: oop (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command modifies a sample output by ");
      printf("some combinations of selected\n");
      printf("outputs.\n");
      printf("NOTE: to add a new output with combination of ");
      printf("existing outputs, use\n");
      printf("      oadd2 to add an output and then run this command.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("Select from the following transformations : \n");
      printf("1. Linear transformation : \n");
      printf("   output <1> = <a> * output <2> + <b> * output <3>\n");
      printf("2. Take an output to a power: output <1> = (output <1>)^P\n");
      printf("3. Replace an input with its logarithm : ");
      printf("output = log(output)\n");
      printf("4. Division: output <1> = output <2> / output <3>\n");
      sprintf(pString, "Enter a selection : ");
      int option = getInt(1, 4, pString);

      if (option == 1)
      {
        printf("Linear transformation : \n");
        printf("   output <1> = <a> * output <2> + <b> * output <3>\n");
        sprintf(pString, "Enter output <1> (1 - %d) : ", nOutputs_);
        ii = getInt(1, nOutputs_, pString);
        ii--;
        sprintf(pString, "Enter the value <a>: ");
        double aVal = getDouble(pString);
        sprintf(pString, "Enter output <2> (1 - %d) : ", nOutputs_);
        jj = getInt(1, nOutputs_, pString);
        jj--;
        sprintf(pString, "Enter the value <b>: ");
        double bVal = getDouble(pString);
        sprintf(pString, "Enter output <3> (1 - %d) : ", nOutputs_);
        kk = getInt(1, nOutputs_, pString);
        kk--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamOutputs_[sInd*nOutputs_+ii] = 
                        aVal * VecSamOutputs_[sInd*nOutputs_+jj] + 
                        bVal * VecSamOutputs_[sInd*nOutputs_+kk]; 
        }
      }
      else if (option == 2)
      {
        printf("Take an output to a power: output = (output)^P\n");
        sprintf(pString, "Enter output (1 - %d) : ", nOutputs_);
        ii = getInt(1, nOutputs_, pString);
        ii--;
        sprintf(pString, "Enter the power <P>: ");
        double pVal = getDouble(pString);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamOutputs_[sInd*nOutputs_+ii] = 
            pow(VecSamOutputs_[sInd*nOutputs_+ii],pVal);
        }
      }
      else if (option == 3)
      {
        printf("Replace an input with its logarithm : ");
        printf("output = log(output)\n");
        sprintf(pString, "Enter output (1 - %d) : ", nOutputs_);
        ii = getInt(1, nOutputs_, pString);
        ii--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          VecSamOutputs_[sInd*nOutputs_+ii] = 
            log(VecSamOutputs_[sInd*nOutputs_+ii]);
        }
      }
      else if (option == 4)
      {
        printf("Division: output <1> = output <2> / output <3>\n");
        sprintf(pString, "Enter output <1> (1 - %d) : ", nOutputs_);
        ii = getInt(1, nOutputs_, pString);
        ii--;
        sprintf(pString, "Enter output <2> (1 - %d) : ", nOutputs_);
        jj = getInt(1, nOutputs_, pString);
        jj--;
        sprintf(pString, "Enter output <3> (1 - %d) : ", nOutputs_);
        kk = getInt(1, nOutputs_, pString);
        kk--;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (VecSamOutputs_[sInd*nOutputs_+kk] == 0)
            printf("ERROR: sample %d output %d = 0\n",sInd+1,kk+1); 
          else
            VecSamOutputs_[sInd*nOutputs_+ii] = 
                        VecSamOutputs_[sInd*nOutputs_+jj] / 
                        VecSamOutputs_[sInd*nOutputs_+kk]; 
        }
      }
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("oop completed. Use write and load again to continue.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ioop 
    //**/ manipulate inputs and outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ioop"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ioop: form linear combinations of inputs to form output\n");
        printf("Syntax: ioop (no argument needed)\n");
        printf("Sometimes one may want to combine some inputs to form\n");
        printf("a new output (e.g. create constraints and filter).\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command modifies a sample output by some ");
      printf("linear combinations of\n");
      printf("selected inputs. That is, \n");
      printf("output <1> = <a> * input <1> + <b> * input <2>\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output <1> (1 - %d) : ", nOutputs_);
      ii = getInt(1, nOutputs_, pString);
      ii--;
      sprintf(pString, "Enter the value <a>: ");
      double aVal = getDouble(pString);
      sprintf(pString, "Enter input <1> (1 - %d) : ", nInputs_);
      jj = getInt(1, nInputs_, pString);
      jj--;
      sprintf(pString, "Enter the value <b>: ");
      double bVal = getDouble(pString);
      sprintf(pString, "Enter input <2> (1 - %d) : ", nInputs_);
      kk = getInt(1, nInputs_, pString);
      kk--;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        VecSamOutputs_[sInd*nOutputs_+ii] = 
                    aVal * VecSamInputs_[sInd*nInputs_+jj] + 
                    bVal * VecSamInputs_[sInd*nInputs_+kk]; 
      }
      psuadeIO_->updateOutputSection(nSamples_, nOutputs_, 
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), StrOutNames_.getStrings());
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ioop completed. Use write and load again to continue.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ nna
    //**/ nearest neighbor analysis (for detecting outliers) 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "nna"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("nna: nearest neighbor analysis (for detecting outliers)\n");
        printf("Syntax: nna (no argument needed)\n");
        printf("For each sample point, the output will be compared to\n");
        printf("its nearest neighbor. If it is an outlier, it will be\n");
        printf("shown to have large gradient.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("nna: for each data point, find nearest neighbor and");
      printf(" plot changes in the outputs.\n");
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      //**/ write to scilab/matlab file
      if (plotScilab())
      {
        fp = fopen("scilabnna.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabnna.sci.\n");
          cmdStatus = 1;
          continue;
        }
      }
      else
      {
        fp = fopen("matlabnna.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabnna.m.\n");
          cmdStatus = 1;
          continue;
        }
      }
      sprintf(pString," nearest neighbor analysis");
      fwriteComment(fp, pString);
      sprintf(pString," for detecting outliers.");
      fwriteComment(fp, pString);
      sprintf(pString," The following plot is : ");
      fwriteComment(fp, pString);
      sprintf(pString," Y-axis: delta output / distance");
      fwriteComment(fp, pString);
      sprintf(pString," X-axis: sample number");
      fwriteComment(fp, pString);
      sprintf(pString," Column 3: nearest neighbor");
      fwriteComment(fp, pString);
      sprintf(pString," Column 4: distance with nearest neighbor");
      fwriteComment(fp, pString);
      fprintf(fp, "A = [\n");
      double minDist;
      for (ii = 0; ii < nSamples_; ii++)
      {
        minDist = PSUADE_UNDEFINED;
        ind = -1;
        for (kk = 0; kk < nSamples_; kk++)
        {
          if (ii != kk)
          {
            ddata = 0.0;
            for (jj = 0; jj < nInputs_; jj++)
            {
              dtemp = VecSamInputs_[ii*nInputs_+jj] - 
                      VecSamInputs_[kk*nInputs_+jj];
              dtemp /= (VecIUpperBs_[jj] - VecILowerBs_[jj]);
              ddata += pow(dtemp, 2.0);
            }
            if (ddata < minDist && ddata > 0.0)
            {
              minDist = ddata;
              ind = kk;
            }
          }
        }
        if (minDist > 0)
          dtemp = (VecSamOutputs_[ii*nOutputs_+outputID] -
                   VecSamOutputs_[ind*nOutputs_+outputID])/minDist;
        else dtemp = 0.0;
        fprintf(fp, "%7d %24.16e %d %e\n", ii+1, dtemp, ind+1, minDist);
      }
      fprintf(fp, "];\n");
      if (plotScilab())
           fprintf(fp, "plot(A(:,1),A(:,2),'x');\n");
      else fprintf(fp, "plot(A(:,1),A(:,2),'x','MarkerSize',12)\n");
      fwritePlotXLabel(fp, "Sample number");
      fwritePlotYLabel(fp, "Delta output / distance");
      fwritePlotTitle(fp, "Nearest Neighbor Analysis");
      fwritePlotAxes(fp);
      fclose(fp);
      if (plotScilab())
           printf("Nearest neighbor result is now in file scilabnna.sci\n");
      else printf("Nearest neighbor result is not in file matlabnna.m.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ setranseed
    //**/ set random number generator seed
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "setranseed"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("setranseed: set random number generator seed\n");
        printf("Syntax: setranseed <positive integer>\n");
        continue;
      }
      long rseed;
      sscanf(lineIn,"%s %ld",command,&rseed);
      if (rseed > 0) 
      {
        PSUADE_randInit(rseed);
        printf("New random seed = %ld\n", rseed);
        cmdStatus = 0;
      }
      else
      {
        printf("ERROR: invalid seed - no change\n");
        cmdStatus = 1;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ interface_track
    //**/ polynomial regression for interface
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "interface_track"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("interface_track: threshold interface tracking \n");
        printf("Syntax: interface_track (no argument needed)\n");
        printf("This command tracks the interface between ");
        printf("the region in the input\n");
        printf("space where the output value is equal to ");
        printf("some user-specified value.\n");
        printf("The result will be a Matlab figure and a ");
        printf("list of suggested new sample\n");
        printf("points on the interface for further exploration.\n");
        printf("This works only for up to 3 inputs.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ > 2)
      {
        printf("ERROR: currently only support nInputs <= 2.\n");
        cmdStatus = 1;
        continue;
      }
      int inp1=0, inp2=1, inp3=2;
      sprintf(pString, 
              "Select the first input (1 - %d) = ",nInputs_);
      inp1 = getInt(1, nInputs_, pString);
      inp1--;
      inp2 = -1;
      while ((inp2 < 1 || inp2 > nInputs_) && (inp2 != inp1))
      {
        sprintf(pString, 
            "Select the second input (1 - %d) = ",nInputs_);
        inp2 = getInt(1, nInputs_, pString);
        inp2--;
        if (inp2 == inp1) 
          printf("ERROR: input %d already used\n",inp2+1);
      }
      inp3 = -2;
      while ((inp3 < -1 || inp3 > nInputs_) && (inp2 != inp1) &&
             (inp3 != inp2))
      {
        sprintf(pString, 
            "Select the third input (1 - %d, or 0 if none) = ",
            nInputs_);
        inp3 = getInt(0, nInputs_, pString);
        inp3--;
        if (inp3 == inp1 || inp3 == inp2) 
          printf("ERROR: input %d already used\n",inp3+1);
      }
      sprintf(pString,"Which output to use (1 - %d) = ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int inpCnt = 0;
      if (inp1 != -1) inpCnt++;
      if (inp2 != -1) inpCnt++;
      if (inp3 != -1) inpCnt++;

      psVector vecInpSettings, vecNewLBs, vecNewUBs;
      vecInpSettings.setLength(nInputs_);
      if (nInputs_-inpCnt > 0)
      {
        sprintf(pString,
                "Set other inputs at their mid points? (y or n) ");
        getString(pString, winput);
        if (winput[0] == 'y')
        {
          for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
          {
            if (iInd1 != inp1 && iInd1 != inp2 && iInd1 != inp3)
                 vecInpSettings[iInd1] = 0.5*
                          (VecIUpperBs_[iInd1]+VecILowerBs_[iInd1]);
          }
        }
        else
        {
          for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
          {
            if (iInd1 != inp1 && iInd1 != inp2 && iInd1 != inp3)
            {
              sprintf(pString,
                    "Enter nominal value for input %d (%e - %e): ",
                    iInd1+1,VecILowerBs_[iInd1],VecIUpperBs_[iInd1]);
              vecInpSettings[iInd1] = getDouble(pString);
            }
          }
        }
      }
      vecNewLBs = VecILowerBs_;
      vecNewUBs = VecIUpperBs_;
      sprintf(pString,
              "Use a different lower/upper bounds for the inputs? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (ii == inp1 || ii == inp2)
          {
            vecNewLBs[ii] = vecNewUBs[ii];
            while (vecNewLBs[ii] >= vecNewUBs[ii])
            {
              printf("Input %d: current lower bound = %12.4e\n",ii+1,
                      VecILowerBs_[ii]);
              sprintf(pString, "New lower bound ? ");
              vecNewLBs[ii] = getDouble(pString);
              printf("Input %d: current upper bound = %12.4e\n",ii+1,
                      VecIUpperBs_[ii]);
              sprintf(pString, "New upper bound ? ");
              vecNewUBs[ii] = getDouble(pString);
              if (vecNewLBs[ii] >= vecNewUBs[ii])
                printf("ERROR: lower bound cannot be >= upper bound\\n");
            }
          }
        }
      }
           
      //**/ create response surface
      printf("Build a response surface (RS) to populate the input space.\n");
      faFlag = 1;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(vecNewLBs.getDVector(),vecNewUBs.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      psVector vecYT;
      vecYT.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++)
        vecYT[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
      faPtr->initialize(VecSamInputs_.getDVector(),vecYT.getDVector());

      //**/ create lattice 
      printf("Next, the RS is to be evaluated at selected point ");
      printf("in the input space.\n");
      int lattice1, lattice2, lattice3, totalPts;

      if (inpCnt == 2)
      {
        printf("It will be a 2D lattice of size MxN.\n");
        sprintf(pString,"What M to use ? (10 - 1000) = ");
        lattice1 = getInt(10, 1000, pString);
        sprintf(pString,"What N to use ? (10 - 1000) = ");
        lattice2 = getInt(10, 1000, pString);
        totalPts = lattice1 * lattice2;
      }
      if (inpCnt == 3)
      {
        printf("It will be a 3D lattice of size MxNxP.\n");
        sprintf(pString,"What M to use ? (10 - 100) = ");
        lattice1 = getInt(10, 100, pString);
        sprintf(pString,"What N to use ? (10 - 100) = ");
        lattice2 = getInt(10, 100, pString);
        sprintf(pString,"What P to use ? (10 - 100) = ");
        lattice3 = getInt(10, 100, pString);
        totalPts = lattice1 * lattice2 * lattice3;
      }

      //**/ response surface evaluation
      psVector vecFAX, vecFAY;
      if (inpCnt == 2)
      {
        //**/ prepare and evaluate with response surface
        vecFAX.setLength(totalPts*nInputs_);
        vecFAY.setLength(totalPts);
        int    index;
        double HX1 = (vecNewUBs[inp1]-vecNewLBs[inp1])/(lattice1 - 1);
        double HX2 = (vecNewUBs[inp2]-vecNewLBs[inp2])/(lattice2 - 1);
        for (ii = 0; ii < lattice1; ii++)
        {
          for (jj = 0; jj < lattice2; jj++)
          {
            index = ii * lattice2 + jj;
            for (kk = 0; kk < nInputs_; kk++)
            {
              if (kk == inp1)
                vecFAX[index*nInputs_+kk] = HX1*ii + vecNewLBs[inp1];
              else if (kk == inp2)
                vecFAX[index*nInputs_+kk] = HX2*jj + vecNewLBs[inp2];
              else
                vecFAX[index*nInputs_+kk] = vecInpSettings[kk];
            }
          }
        }
        printf("Please wait while generating the RS data \n");
        faPtr->evaluatePoint(totalPts,vecFAX.getDVector(),
                             vecFAY.getDVector());
        double Ymin = vecFAY[0];
        for (ii = 1; ii < totalPts; ii++)
          if (vecFAY[ii] < Ymin) Ymin = vecFAY[ii];
        double Ymax = vecFAY[0];
        for (ii = 1; ii < totalPts; ii++)
          if (vecFAY[ii] > Ymax) Ymax = vecFAY[ii];
        printf("Output Ymin and Ymax = %e %e.\n", Ymin, Ymax);
        sprintf(pString,"Enter threshold : ");
        double thresh = getDouble(pString);

        //**/ create a user-specified number of points near the 
        //**/ interface
        psVector vecCandX, vecCandY;
        vecCandX.setLength(2*totalPts);
        vecCandY.setLength(totalPts);
        int nCand = 0;
        psVector vecYY;
        vecYY.setLength(totalPts*4);
        for (ii = 0; ii < totalPts; ii++)
        {
          count = 0;
          if ((ii % lattice1) > 0)
          {
            vecYY[nCand*4+3] = vecFAY[ii-1];
            if ((vecFAY[ii-1] >= thresh && vecFAY[ii] < thresh) ||
                (vecFAY[ii-1] <= thresh && vecFAY[ii] > thresh)) 
              count++; 
          }
          if ((ii % lattice1) < lattice1-1)
          {
            vecYY[nCand*4+2] = vecFAY[ii+1];
            if ((vecFAY[ii+1] >= thresh && vecFAY[ii] < thresh) ||
                (vecFAY[ii+1] <= thresh && vecFAY[ii] > thresh))
              count++; 
          }
          if ((ii / lattice1) > 0)
          {
            vecYY[nCand*4+1] = vecFAY[ii-lattice2];
            if ((vecFAY[ii-lattice2] >= thresh && vecFAY[ii] < thresh) ||
                (vecFAY[ii-lattice2] <= thresh && vecFAY[ii] > thresh)) 
              count++; 
          }
          if ((ii / lattice1) < lattice2-1)
          {
            vecYY[nCand*4] = vecFAY[ii+lattice2];
            if ((vecFAY[ii+lattice2] >= thresh && vecFAY[ii] < thresh) ||
                (vecFAY[ii+lattice2] <= thresh && vecFAY[ii] > thresh)) 
              count++; 
          }
          if (count > 0)
          {
            printf("Candidate %d: X1 = %10.3e, X2 = %10.3e, Y = %10.3e\n",
                   ii+1,vecFAX[ii*nInputs_+inp1],vecFAX[ii*nInputs_+inp2],
                   vecFAY[ii]);
            printf("  YE = %10.3e, YW = %10.3e, YN = %10.3e, YS = %10.3e\n",
                   vecYY[nCand*4], vecYY[nCand*4+1], vecYY[nCand*4+2], 
                   vecYY[nCand*4+3]);
            vecCandX[nCand*2] = vecFAX[ii*nInputs_+inp1];
            vecCandX[nCand*2+1] = vecFAX[ii*nInputs_+inp2];
            vecCandY[nCand] = vecFAY[ii];
            nCand++;
          }
        } 
        if (nCand == 0) 
        {
          vecCandX.clean();
          printf("INFO: no threshold interface point has been found.\n");
        }
        else
        {
          fp = fopen("itrack2_data", "w");
          if (fp != NULL)
          {
            fprintf(fp,"# Each line: 2 inputs X, Y, Y(4 neighbors)\n");
            fprintf(fp,"%d 2 5\n", nCand);
            for (ii = 0; ii < nCand; ii++)
            {
              fprintf(fp,"%12.4e %12.4e %10.3e ",
                      vecCandX[ii*2],vecCandX[ii*2+1],vecCandY[ii]);
              fprintf(fp,"%10.3e %10.3e %10.3e %10.3e\n",
                      vecYY[ii*4],vecYY[ii*4+1],vecYY[ii*4+2],vecYY[ii*4+3]);
            }
            fclose(fp); 
            printAsterisks(PL_INFO, 0);
            printf("INFO: threshold interface points are in itrack2_data\n");
            printf("      Number of interface points = %d\n",nCand);
            printf("      To select a smaller subset, load this data ");
            printf("set into psuade and\n");
            printf("      use the odoe_mmd command.\n");
            printf("      Also, use rs2 to visualize response surface.\n");
            printAsterisks(PL_INFO, 0);
          }
          vecCandX.clean();
        }
      }
      if (inpCnt == 3)
      {
#if 0
        printf("Please wait while generating the RS data \n");
        faPtr->gen3DGridData(VecSamInputs_.getDVector(),
               vecYT.getDVector(),inp1,inp2,inp3,
               vecInpSettings.getDVector(),&totalPts,&faXOut,&faYOut);
        fp = fopen("matlabitrack3.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabitrack3.m.\n");
          delete [] faXOut;
          delete [] faYOut;
          delete faPtr;
          continue;
        }
        fwritePlotCLF(fp);
        fprintf(fp,"xlo = %e; \n", VecILowerBs_[inp2]);
        fprintf(fp,"xhi = %e; \n", VecIUpperBs_[inp2]);
        fprintf(fp,"ylo = %e; \n", VecILowerBs_[inp1]);
        fprintf(fp,"yhi = %e; \n", VecIUpperBs_[inp1]);
        fprintf(fp,"zlo = %e; \n", VecILowerBs_[inp3]);
        fprintf(fp,"zhi = %e; \n", VecIUpperBs_[inp3]);
        fprintf(fp,"X=zeros(%d,%d,%d);\n",lattice1,lattice2,lattice3);
        fprintf(fp,"Y=zeros(%d,%d,%d);\n",lattice1,lattice2,lattice3);
        fprintf(fp,"Z=zeros(%d,%d,%d);\n",lattice1,lattice2,lattice3);
        fprintf(fp,"V=zeros(%d,%d,%d);\n",lattice1,lattice2,lattice3);
        for (jj = 0; jj < latticeN; jj++)
        {
          fprintf(fp,"Y(:,:,%d) = [\n", jj + 1);
          for (sInd = 0; sInd < latticeN; sInd++)
          {
            for (ii = 0; ii < latticeN; ii++)
            {
              ind = sInd*latticeN*latticeN+ii*latticeN+jj;
              fprintf(fp,"%e ", faXOut[ind*3]);
            }
            fprintf(fp,"\n");
          }
          fprintf(fp, "];\n");
          fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
          for (sInd = 0; sInd < latticeN; sInd++)
          {
            for (ii = 0; ii < latticeN; ii++)
            {
              ind = sInd*latticeN*latticeN+ii*latticeN+jj;
              fprintf(fp, "%e ", faXOut[ind*3+1]);
            }
            fprintf(fp, "\n");
          }
          fprintf(fp, "];\n");
          fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
          for (sInd = 0; sInd < latticeN; sInd++)
          {
            for (ii = 0; ii < latticeN; ii++)
            {
              ind = sInd*latticeN*latticeN+ii*latticeN+jj;
              fprintf(fp, "%e ", faXOut[ind*3+2]);
            }
            fprintf(fp, "\n");
          }
          fprintf(fp, "];\n");
        }
        double Ymin = faYOut[0];
        for (sInd = 1; sInd < totalPts; sInd++)
          if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
        double Ymax = faYOut[0];
        for (sInd = 1; sInd < totalPts; sInd++)
          if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
        printf("Output Ymin and Ymax = %e %e.\n", Ymin, Ymax);
        sprintf(pString,"Enter threshold : ");
        double thresh = getDouble(pString);
        int myCnt=0;
        for (jj = 0; jj < latticeN; jj++)
        {
          fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
          for (sInd = 0; sInd < latticeN; sInd++)
          {
            for (ii = 0; ii < latticeN; ii++)
            {
              ind = sInd*latticeN*latticeN+ii*latticeN+jj;
              if (faYOut[ind] < thresh)
              {
                fprintf(fp, "%e ", Ymin);
                myCnt++;
              }
              else fprintf(fp, "%e ", faYOut[ind]);
            }
            fprintf(fp, "\n");
          }
          fprintf(fp, "];\n");
        }
        if (myCnt == latticeN*latticeN*latticeN)
        {
          fprintf(fp, "V(1,1,1)=0;\n");
          fprintf(fp, "V(%d,%d,%d)=1;\n",latticeN,latticeN,latticeN);
        }
        fprintf(fp,"xt = [%e:%e:%e];\n", VecILowerBs_[inp2],
            (VecIUpperBs_[inp2]-VecILowerBs_[inp2])*0.01, 
             VecIUpperBs_[inp2]);
        fprintf(fp,"yt = [%e:%e:%e];\n", VecILowerBs_[inp1],
            (VecIUpperBs_[inp1]-VecILowerBs_[inp1])*0.01, 
             VecIUpperBs_[inp1]);
        fprintf(fp,"zt = [%e:%e:%e];\n", VecILowerBs_[inp3],
            (VecIUpperBs_[inp3]-VecILowerBs_[inp3])*0.01, 
             VecIUpperBs_[inp3]);
        fprintf(fp,"isoval = %e;\n", thresh);
        fprintf(fp,"h = patch(isosurface(X,Y,Z,V,isoval),... \n");
        fprintf(fp,"          'FaceColor', 'blue', ... \n");
        fprintf(fp,"          'EdgeColor', 'none', ... \n");
        fprintf(fp,"          'AmbientStrength', 0.2, ... \n");
        fprintf(fp,"          'SpecularStrength', 0.7, ... \n");
        fprintf(fp,"          'DiffuseStrength', 0.4);\n");
        fprintf(fp,"isonormals(X,Y,Z,V,h);\n");
        fprintf(fp,"patch(isocaps(X,Y,Z,V,isoval), ...\n");
        fprintf(fp,"      'FaceColor', 'interp', ... \n");
        fprintf(fp,"      'EdgeColor', 'none'); \n");
        fprintf(fp,"axis([xlo xhi ylo yhi zlo zhi])\n");
        fprintf(fp,"daspect([%e,%e,%e])\n",
            VecIUpperBs_[inp2]-VecILowerBs_[inp2],
            VecIUpperBs_[inp1]-VecILowerBs_[inp1],
            VecIUpperBs_[inp3]-VecILowerBs_[inp3]);
        fprintf(fp,"   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
            StrInpNames_[inp2]);
        fprintf(fp,"   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
            StrInpNames_[inp1]);
        fprintf(fp,"   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
            StrInpNames_[inp3]);
        fprintf(fp,"   title('%s','Fontsize',12,'FontWeight','bold')\n",
            StrOutNames_[outputID]);
        fwritePlotAxes(fp);
        fprintf(fp,"colormap('default'); colorbar\n");
        fprintf(fp,"%%axis tight\n");
        fprintf(fp,"view(3) \n");
        fprintf(fp,"set(gcf,'Renderer','zbuffer')\n");
        fprintf(fp,"lighting phong\n");
        fclose(fp);
        printf("Visualization file in matlabitrack3.m\n");
#endif
      }
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ pdfconvert 
    //**/ convert data set based on distribution
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "pdfconvert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("pdfconvert: convert a sample based on its pdfs\n");
        printf("Syntax: pdfconvert (no argument needed)\n");
        printf("To use this, first create a sample using uniform\n");
        printf("distribution. Then load the sample (make sure that\n");
        printf("before you load the sample, you have modified the\n");
        printf("input section of this sample file to reflect the\n");
        printf("desired distribution. Then use this command to\n");
        printf("convert the sample to the desired distributions.\n");
        printf("Use write to store the converted sample.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command converts the loaded sample based on ");
      printf("the PDF prescription\n");
      printf("given in the loaded file.\n");
      printDashes(PL_INFO, 0);

      PDFManager *pdfman = new PDFManager();
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      kk = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1, -1);
      psuadeIO_->updateMethodSection(PSUADE_SAMP_MC, -1, -1, -1, -1);
      pdfman->initialize(psuadeIO_);
      vecInps.load(nSamples_*nInputs_, VecSamInputs_.getDVector());
      vecOuts.setLength(nSamples_*nInputs_);
      vecUBs.load(nInputs_, VecIUpperBs_.getDVector());
      vecLBs.load(nInputs_, VecILowerBs_.getDVector());
      pdfman->invCDF(nSamples_, vecInps, vecOuts, vecLBs, vecUBs);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1, -1);
      for (ii = 0; ii < nSamples_*nInputs_; ii++)
        VecSamInputs_[ii] = vecOuts[ii];

      psIVector veciPDFs;
      psVector  veciMeans, veciStds;
      veciPDFs.setLength(nInputs_);
      veciMeans.setLength(nInputs_);
      veciStds.setLength(nInputs_);
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        veciPDFs[ii] = 0;
        veciMeans[ii] = 0;
        veciStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                       NULL,NULL,VecSamInputs_.getDVector(),NULL, 
                       veciPDFs.getIVector(),veciMeans.getDVector(),
                       veciStds.getDVector(),iCMat);
      delete iCMat;
      psuadeIO_->updateAnalysisSection(-1, -1, -1, kk, -1, 0);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      delete pdfman;
      printf("INFO: The sample in local memory has been converted ");
      printf("based on the PDF\n");
      printf("      information in the INPUT section.\n");
      printf("      You can now write your converted sample to a file.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rand_draw 
    //**/ draw a random sample from the resident sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rand_draw"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rand_draw: draw a random sample from the loaded sample\n");
        printf("Syntax: rand_draw <n>\n");
        printf("That is, this command generates a bootstrapped sample\n");
        printf("from the sample which has been loaded to local memory.\n");
        printf("This command is used if your want to draw a sample\n");
        printf("from the posterior sample after Bayesian analysis.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("This command draws a sample with replacement from the ");
      printf("loaded sample.\n");
      printf("Size of the loaded sample = %d.\n",nSamples_);
      sprintf(pString,"Size of the sample to be drawn : (1-1000000) ");
      count = getInt(1, 1000000, pString);
      vecXT.setLength(nInputs_ * count);
      vecYT.setLength(nOutputs_ * count);
      vecST.setLength(count);
      for (ii = 0; ii < count; ii++)
      {
        ind = PSUADE_rand() % nSamples_;
        for (jj = 0; jj < nInputs_; jj++)
          vecXT[ii*nInputs_+jj] = VecSamInputs_[ind*nInputs_+jj];
        for (jj = 0; jj < nOutputs_; jj++)
          vecYT[ii*nOutputs_+jj] = VecSamOutputs_[ind*nOutputs_+jj];
        vecST[ii] = VecSamStates_[ind];
      }
      psIVector veciPDFs;
      psVector  veciMeans, veciStds;
      veciPDFs.setLength(nInputs_);
      veciMeans.setLength(nInputs_);
      veciStds.setLength(nInputs_);
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        veciPDFs[ii] = 0;
        veciMeans[ii] = 0;
        veciStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->updateInputSection(count, nInputs_, NULL, 
               VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
               vecXT.getDVector(), StrInpNames_.getStrings(),
               veciPDFs.getIVector(),veciMeans.getDVector(),
               veciStds.getDVector(),iCMat);
      delete iCMat;
      ioPtr->updateOutputSection(count, nOutputs_, vecYT.getDVector(), 
               vecST.getIVector(), StrOutNames_.getStrings());
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
      printf("Store random sample to : (filename) ");
      scanf("%s", dataFile);
      fgets(lineIn,5000,stdin); 
      if ((fp = fopen(dataFile, "w")) == NULL)
      {
        printf("ERROR: cannot open file %s\n", dataFile);
        cmdStatus = 1;
      }
      else
      {
        ioPtr->writePsuadeFile(dataFile,0);
        printf("The randomly drawn sample has been saved in %s.\n",dataFile);
        cmdStatus = 0;
      }    
      delete ioPtr;
    }

    //**/ -------------------------------------------------------------
    // +++ rand_drawb 
    //**/ block draw a random sample from the resident sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rand_drawb"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rand_drawb: draw a random sample from the loaded sample\n");
        printf("Syntax: rand_drawb <n>\n");
        printf("That is, this command creates a new sample ");
        printf("by randomly drawing blocks\n");
        printf("of contiguous sample points from the ");
        printf("loaded sample.  This command is\n");
        printf("useful for nondeterministic ");
        printf("simulations when the loaded sample has\n");
        printf("blocks of contiguous sample points ");
        printf("whereby each block has identical\n");
        printf("sample points but the outputs are different due ");
        printf("to stochasticity in\n"); 
        printf("simulation results.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("This command draws a sample with replacement from the ");
      printf("loaded sample.\n");
      printf("It differs from rand_draw in that it draws ");
      printf("blocks of identical and\n");
      printf("continguous sample points.\n");
      int blksize=1;
      for (ss = 1; ss < nSamples_; ss++)
      {
        for (ii = 0; ii < nInputs_; ii++)
          if (VecSamInputs_[ss*nInputs_+ii] != VecSamInputs_[ii])
            break;
        if (ii < nInputs_) break;
        blksize++;
      }
      if (blksize == 1 || (1.0*nSamples_/blksize*blksize != nSamples_))
      {
        printf("ERROR: The loaded sample does not have continguous ");
        printf("blocks of identical\n");
        printf("sample points with fixed size. So this command ");
        printf("is not applicable.\n");
        cmdStatus = 1;
        continue;
      }
      int ss2;
      for (ss = blksize; ss < nSamples_; ss+=blksize)
      {
        for (ss2 = 1; ss2 < blksize; ss2++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            if (VecSamInputs_[(ss+ss2)*nInputs_+ii] != 
                VecSamInputs_[ss*nInputs_+ii])
              break;
          if (ii < nInputs_) break;
        }
        if (ss2 != blksize) break;
      }
      if (ss != nSamples_)
      {
        printf("ERROR: The loaded sample does not have continguous ");
        printf("blocks of identical\n");
        printf("sample points with fixed size. So this command ");
        printf("is not applicable.\n");
        cmdStatus = 1;
        continue;
      } 
      printf("The loaded sample has %d blocks each with %d identical points.\n",
             nSamples_/blksize, blksize);
      printf("So the sample size to be drawn should be multiples of %d.\n",
             blksize);
      sprintf(pString,"Size of the sample to be drawn (multiples of %d) ",
              blksize);
      int samSize = getInt(blksize, 10000000, pString);
      if (samSize/blksize*blksize != samSize)
      {
        printf("Your sample size is not a multiple of %d.\n",blksize);
        samSize = (samSize / blksize + 1) * blksize;
        printf("Your sample size has been changed to %d.\n",samSize);
      }
      vecXT.setLength(nInputs_ * samSize);
      vecYT.setLength(nOutputs_ * samSize);
      vecST.setLength(samSize);
      for (ss = 0; ss < samSize/blksize; ss++)
      {
        ind = PSUADE_rand() % (nSamples_/blksize);
        for (kk = 0; kk < blksize; kk++)
        {
          for (jj = 0; jj < nInputs_; jj++)
            vecXT[(ss*blksize+kk)*nInputs_+jj] =
                  VecSamInputs_[(ind*blksize+kk)*nInputs_+jj];
          for (jj = 0; jj < nOutputs_; jj++)
            vecYT[(ss*blksize+kk)*nOutputs_+jj] =
                  VecSamOutputs_[(ind*blksize+kk)*nOutputs_+jj];
          vecST[ss*blksize+kk] = VecSamStates_[ind*blksize+kk];
        }
      }
      psIVector veciPDFs;
      psVector  veciMeans, veciStds;
      veciPDFs.setLength(nInputs_);
      veciMeans.setLength(nInputs_);
      veciStds.setLength(nInputs_);
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        veciPDFs[ii] = 0;
        veciMeans[ii] = 0;
        veciStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->updateInputSection(samSize, nInputs_, NULL, 
               VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
               vecXT.getDVector(), StrInpNames_.getStrings(),
               veciPDFs.getIVector(), veciMeans.getDVector(),
               veciStds.getDVector(),iCMat);
      delete iCMat;
      ioPtr->updateOutputSection(samSize, nOutputs_, vecYT.getDVector(),
                   vecST.getIVector(), StrOutNames_.getStrings());
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, samSize, 1, -1, -1);
      printf("Store random sample to : (filename) ");
      scanf("%s", dataFile);
      fgets(lineIn,5000,stdin);
      if ((fp = fopen(dataFile, "w")) == NULL)
      {
        printf("ERROR: cannot open file %s\n", dataFile);
        cmdStatus = 1;
      }
      else
      {
        ioPtr->writePsuadeFile(dataFile,0);
        printf("The randomly drawn sample has been saved in %s.\n",dataFile);
        cmdStatus = 0;
      }
      delete ioPtr;
    }

    //**/ -------------------------------------------------------------
    // +++ rand_draw2 
    //**/ draw a random sample from the 2 samples
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rand_draw2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rand_draw2: draw a sample randomly from 2 samples\n");
        printf("Syntax: rand_draw2 <n>\n");
        printf("That is, this command creates a new sample from 2 ");
        printf("sample files\n");
        printf("in PSUADE data format. This command is useful if, ");
        printf("e.g. you have\n");
        printf("2 MCMC posterior samples with two distinct sets of ");
        printf("inputs (say m1\n");
        printf("and m2), and you want to create a new sample with ");
        printf("m1+m2 inputs by\n");
        printf("randomly drawing and concatenating points from the 2 samples.\n");
        continue;
      }
      if (psuadeIO_ != NULL) delete psuadeIO_;

      //**/ read file names
      printf("This command randomly draws from 2 samples and ");
      printf("concatenates the inputs.\n");
      printf("Both files have to be in PSUADE data format.\n");
      sprintf(pString,"Enter name of the first sample file : ");
      getString(pString, winput);
      kk = strlen(winput);
      winput[kk-1] = '\0';
      psuadeIO_ = new PsuadeData;
      status = psuadeIO_->readPsuadeFile(winput);
      if (status != 0)
      {
        printf("rand_draw2 ERROR READING FILE %s\n", winput);
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter name of the second sample file : ");
      getString(pString, winput);
      kk = strlen(winput);
      winput[kk-1] = '\0';
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(winput);
      if (status != 0)
      {
        printf("rand_draw2 ERROR READING FILE %s\n", winput);
        cmdStatus = 1;
        continue;
      }

      //**/ read data from files ==> psuadeIO_ and ioPtr
      psuadeIO_->getParameter("input_ninputs", pPtr);
      nInputs_ = pPtr.intData_;
      psuadeIO_->getParameter("method_nsamples", pPtr);
      nSamples_ = pPtr.intData_;
      psuadeIO_->getParameter("input_sample", pPtr);
      VecSamInputs_.load(nSamples_*nInputs_, pPtr.dbleArray_);
      pPtr.clean();

      ioPtr->getParameter("input_ninputs", pPtr);
      int nInputs2 = pPtr.intData_;
      ioPtr->getParameter("method_nsamples", pPtr);
      int nSamples2 = pPtr.intData_;
      ioPtr->getParameter("input_sample", pPtr);
      double *dPtrX = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;

      //**/ get input lower and upper bounds
      VecILowerBs_.setLength(nInputs_+nInputs2);
      VecIUpperBs_.setLength(nInputs_+nInputs2);
      pLower.clean();
      psuadeIO_->getParameter("input_lbounds", pLower);
      for (ii = 0; ii < nInputs_; ii++) 
        VecILowerBs_[ii] = pLower.dbleArray_[ii];
      pLower.clean();
      ioPtr->getParameter("input_lbounds", pLower);
      for (ii = 0; ii < nInputs2; ii++) 
        VecILowerBs_[nInputs_+ii] = pLower.dbleArray_[ii];
      pLower.clean();

      pUpper.clean();
      psuadeIO_->getParameter("input_ubounds", pUpper);
      for (ii = 0; ii < nInputs_; ii++) 
        VecIUpperBs_[ii] = pUpper.dbleArray_[ii];
      pUpper.clean();
      ioPtr->getParameter("input_ubounds", pUpper);
      for (ii = 0; ii < nInputs2; ii++) 
        VecIUpperBs_[nInputs_+ii] = pUpper.dbleArray_[ii];
      pUpper.clean();

      //**/ mix and match
      sprintf(pString,"Size of the sample to be drawn : (1-2000000) ");
      count = getInt(1, 2000000, pString);
      kk = 0;
      if (count == (nSamples_ * nSamples2))
      {
        sprintf(pString,"Form sample tensor product? (y or n) ");
        getString(pString, winput);
        if (winput[0] == 'y') kk = 1;
      }
      vecWT.setLength(count*(nInputs_+nInputs2));
      vecYT.setLength(count);
      vecST.setLength(count);
      for (ii = 0; ii < count; ii++)
      {
        if (kk == 0)
        {
          ind  = PSUADE_rand() % nSamples_;
          ind2 = PSUADE_rand() % nSamples2;
        }
        else
        {
          ind  = (ii * nSamples2) % nSamples_;
          ind2 = ii % nSamples2;
        }
        for (jj = 0; jj < nInputs_; jj++)
          vecWT[ii*(nInputs_+nInputs2)+jj] = 
                VecSamInputs_[ind*nInputs_+jj];
        for (jj = 0; jj < nInputs2; jj++)
          vecWT[ii*(nInputs_+nInputs2)+nInputs_+jj] = 
                dPtrX[ind2*nInputs2+jj];
        vecYT[ii] = PSUADE_UNDEFINED;
        vecST[ii] = 0;
      } 

      //**/ compose new input name list
      StrInpNames_.setNumStrings(nInputs_+nInputs2);
      psuadeIO_->getParameter("input_names", pINames);
      for (ii = 0; ii < nInputs_; ii++)
        StrInpNames_.loadOneString(ii, pINames.strArray_[ii]);
      pINames.clean();
      ioPtr->getParameter("input_names", pINames);
      for (ii = 0; ii < nInputs2; ii++)
        StrInpNames_.loadOneString(nInputs_+ii, pINames.strArray_[ii]);
      pINames.clean();

      //**/ clean up
      delete psuadeIO_;
      delete ioPtr;
      psuadeIO_ = NULL;

      //**/ reset all input distributions to be uniform and store
      ioPtr = new PsuadeData();
      psIVector veciPDFs;
      psVector  veciMeans, veciStds;
      veciPDFs.setLength(nInputs_+nInputs2);
      veciMeans.setLength(nInputs_+nInputs2);
      veciStds.setLength(nInputs_+nInputs2);
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_+nInputs2, nInputs_+nInputs2);
      for (ii = 0; ii < nInputs_+nInputs2; ii++)
      {
        veciPDFs[ii] = 0;
        veciMeans[ii] = 0;
        veciStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      ioPtr->updateInputSection(count, nInputs_+nInputs2, NULL, 
                 VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                 vecWT.getDVector(), StrInpNames_.getStrings(),
                 veciPDFs.getIVector(), veciMeans.getDVector(), 
                 veciStds.getDVector(), iCMat);
      delete iCMat;

      //**/ set OUTPUT and METHOD fields
      nOutputs_ = 1;
      StrOutNames_.setNumStrings(iOne);
      sprintf(pString, "Y");
      StrOutNames_.loadOneString(0, pString);
      ioPtr->updateOutputSection(count,nOutputs_,vecYT.getDVector(), 
                   vecST.getIVector(),StrOutNames_.getStrings());
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);

      printf("Store random sample of size %d to : (filename) ",count);
      scanf("%s", dataFile);
      fgets(lineIn,5000,stdin); 
      if ((fp = fopen(dataFile, "w")) == NULL)
      {
        printf("ERROR: cannot open file %s\n", dataFile);
        cmdStatus = 1;
      }
      else
      {
        ioPtr->writePsuadeFile(dataFile,0);
        printf("The randomly drawn sample has been saved in %s.\n",dataFile);
        printf("This file now has %d inputs.\n",nInputs_+nInputs2);
        cmdStatus = 0;
      }    
      pPtr.dbleArray_ = dPtrX;
      pPtr.clean();
      StrInpNames_.clean();
      StrOutNames_.clean();
      VecSamInputs_.clean();
      VecSamOutputs_.clean();
      VecSamStates_.clean();
      VecILowerBs_.clean();
      VecIUpperBs_.clean();
      VecInpPDFs_.clean();
      VecInpMeans_.clean();
      VecInpStds_.clean();
      nSamples_ = 0;
      nInputs_ = 0;
      nOutputs_ = 0;
      delete ioPtr;
      ioPtr = NULL;
    }

    //**/ -------------------------------------------------------------
    // +++ gensample 
    //**/ create a sample from the PDF information in the loaded sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gensample"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gensample: create a sample from the loaded INPUT PDFs.\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a sample of a user-specified size ");
      printf("based on the\n");
      printf("PDF presciption given in the loaded file.\n");
      printDashes(PL_INFO, 0);

      sprintf(pString,"Sample size? (>= 2, <=10000000) ");
      kk = getInt(1,10000000,pString);
      PDFManager *pdfman = new PDFManager();
      pdfman->initialize(nInputs_,VecInpPDFs_.getIVector(),
                    VecInpMeans_.getDVector(),VecInpStds_.getDVector(),
                    *inputCMat_,SamPDFFiles_,
                    VecSamPDFIndices_.getIVector());
      vecLBs.load(nInputs_, VecILowerBs_.getDVector());
      vecUBs.load(nInputs_, VecIUpperBs_.getDVector());
      vecInps.setLength(kk*nInputs_);
      pdfman->genSample(kk, vecInps, vecLBs, vecUBs);
      delete pdfman;
      PsuadeData *ioPtr = new PsuadeData();
      double *dPtrX = vecInps.getDVector();
      vecYT.setLength(kk);
      for (ii = 0; ii < kk; ii++) vecYT[ii] = PSUADE_UNDEFINED; 
      vecST.setLength(kk);
      ioPtr->updateInputSection(kk,nInputs_,NULL,VecILowerBs_.getDVector(),
               VecIUpperBs_.getDVector(),dPtrX,StrInpNames_.getStrings(),
               VecInpPDFs_.getIVector(),VecInpMeans_.getDVector(),
               VecInpStds_.getDVector(), inputCMat_);

      psStrings Yname;
      if (VecSamOutputs_.length() == 0)
      {
        Yname.setNumStrings(iOne);
        strcpy(pString, "Y");
        Yname.loadOneString(0, pString);
        ioPtr->updateOutputSection(kk, iOne, vecYT.getDVector(), 
                           vecST.getIVector(), Yname.getStrings());
      }
      else 
      {
        ioPtr->updateOutputSection(kk, iOne, vecYT.getDVector(), 
                    vecST.getIVector(), StrOutNames_.getStrings());
      }
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, kk, 1, -1, -1);
      ioPtr->writePsuadeFile("psuade_sample",0);
      printf("The sample has been written into the 'psuade_sample' file.\n");
      delete ioPtr;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ cdf_lookup 
    //**/ look up the cdf of an input value
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "cdf_lookup"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("cdf_lookup: look up the cumulative probability given\n");
        printf("            a random variable value and distribution\n");
        printf("            type and parameter values.\n");
        printf("Syntax: cdf_lookup <n>\n");
        printf("Given a distribution type and its parameters, this\n");
        printf("command returns the cumulative probability when\n");
        printf("provided with a lookup random variable value.\n");
        continue;
      }
      int    gtype;
      double slbound, subound, smean, sstdev;
      psMatrix corMat;

      printAsterisks(PL_INFO, 0);
      printf("Given a PDF type, its parameter values, and the random value ");
      printf("X, this\n");
      printf("command computes the corresponding cumulative probability. ");
      printf("To do this,\n");
      printf("you need to provide the distribution parameter values ");
      printf("and the random\n");
      printf("value X.\n");
      printDashes(PL_INFO, 0);

      sprintf(pString,
        "PDF type: 1-N, 2-L, 3-T, 4-Beta, 5-Weibull, 6-Gamma, 7-Exp, 8-F ?");
      gtype = getInt(1, 8, pString);
      if (gtype == 1)
      {
        gtype = PSUADE_PDF_NORMAL;
        sprintf(pString, "PDF mean = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF std. dev. = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 2) 
      {
        gtype = PSUADE_PDF_LOGNORMAL;
        sprintf(pString, "PDF log(mean) = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF std. dev. = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 3)
      {
        gtype = PSUADE_PDF_TRIANGLE;
        sprintf(pString, "PDF center = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF half base width (assumed isosceles) = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 4)
      {
        gtype = PSUADE_PDF_BETA;
        sprintf(pString, "PDF alpha = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF beta = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 5)
      {
        gtype = PSUADE_PDF_WEIBULL;
        sprintf(pString, "PDF lambda = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF k = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 6)
      {
        gtype = PSUADE_PDF_GAMMA;
        sprintf(pString, "PDF alpha = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF beta = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 7)
      {
        gtype = PSUADE_PDF_EXPONENTIAL;
        sprintf(pString, "PDF lambda = ");
        smean = getDouble(pString);
        sstdev = 0.0;
      }
      else if (gtype == 8)
      {
        gtype = PSUADE_PDF_F;
        sprintf(pString, "PDF d1 = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF d2 = ");
        sstdev = getDouble(pString);
      }
      corMat.setDim(1,1);
      corMat.setEntry(0,0, 1.0e0);
      PDFManager *pdfman = new PDFManager();
      pdfman->initialize(1, &gtype, &smean, &sstdev, corMat, NULL, NULL);
      vecInps.setLength(iOne);
      vecOuts.setLength(iOne);
      sprintf(pString,
              "Enter parameter value to fetch cumulative probability: ");
      ddata = getDouble(pString);
      vecInps.load(1, &ddata);
      slbound = 0;
      subound = 1;
      vecUBs.load(1, &subound);
      vecLBs.load(1, &slbound);
      pdfman->getCDF(1, vecInps, vecOuts, vecLBs, vecUBs);
      printf("Cumulative probability = %e\n", vecOuts[0]);
      delete pdfman;
      pdfman = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ icdf_lookup 
    //**/ look up the value of the random variable given CDF
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "icdf_lookup"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("icdf_lookup: look up the random variable value given\n");
        printf("             a CDF between 0 and 1.\n");
        printf("Syntax: icdf_lookup <n>\n");
        printf("Given a distribution type and its parameters, this\n");
        printf("command returns the random parameter value that\n");
        printf("gives the desired CDF.\n");
        continue;
      }
      int    gtype;
      double slbound, subound, smean, sstdev;
      psMatrix corMat;

      printAsterisks(PL_INFO, 0);
      printf("Given a PDF type, its parameters, and a target CDF value, ");
      printf("this command\n");
      printf("returns the random parameter value X that gives ");
      printf("the target CDF.\n");
      printDashes(PL_INFO, 0);

      sprintf(pString,
        "PDF type: 1-N, 2-L, 3-T, 4-Beta, 5-Weibull, 6-Gamma, 7-Exp, 8-F, 9-Cauchy ? ");
      gtype = getInt(1, 9, pString);
      if (gtype == 1)
      {
        gtype = PSUADE_PDF_NORMAL;
        sprintf(pString, "PDF mean = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF std. dev. = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 2) 
      {
        gtype = PSUADE_PDF_LOGNORMAL;
        sprintf(pString, "PDF log(mean) = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF std. dev. = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 3)
      {
        gtype = PSUADE_PDF_TRIANGLE;
        sprintf(pString, "PDF center = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF half base width (assumed isosceles) = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 4)
      {
        gtype = PSUADE_PDF_BETA;
        sprintf(pString, "PDF alpha = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF beta = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 5)
      {
        gtype = PSUADE_PDF_WEIBULL;
        sprintf(pString, "PDF lambda = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF k = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 6)
      {
        gtype = PSUADE_PDF_GAMMA;
        sprintf(pString, "PDF alpha = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF beta = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 7)
      {
        gtype = PSUADE_PDF_EXPONENTIAL;
        sprintf(pString, "PDF lambda = ");
        smean = getDouble(pString);
        sstdev = 0.0;
      }
      else if (gtype == 8)
      {
        gtype = PSUADE_PDF_F;
        sprintf(pString, "PDF d1 = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF d2 = ");
        sstdev = getDouble(pString);
      }
      else if (gtype == 9)
      {
        gtype = PSUADE_PDF_CAUCHY;
        sprintf(pString, "PDF X0 = ");
        smean = getDouble(pString);
        sprintf(pString, "PDF gamma = ");
        sstdev = getDouble(pString);
      }
      corMat.setDim(1,1);
      corMat.setEntry(0,0, 1.0e0);
      PDFManager *pdfman = new PDFManager();
      pdfman->initialize(1, &gtype, &smean, &sstdev, corMat, NULL, NULL);
      vecInps.setLength(iOne);
      vecOuts.setLength(iOne);
      sprintf(pString, "Enter CDF value (0 - 1): ");
      ddata = getDouble(pString);
      vecInps.load(1, &ddata);
      slbound = -PSUADE_UNDEFINED;
      subound = +PSUADE_UNDEFINED;
      vecUBs.load(1, &subound);
      vecLBs.load(1, &slbound);
      pdfman->invCDF(1, vecInps, vecOuts, vecLBs, vecUBs);
      printf("Random variable value with CDF (%e) = %e\n",
             ddata, vecOuts[0]);
      delete pdfman;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ kde 
    //**/ Create PDF from a sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kde"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kde: create 1D PDF from from a small sample ");
        printf("using a kernel method .\n");
        printf("Syntax: kde\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command uses Use kernel density estimation to ");
      printf("create a histogram\n");
      printf("of a selected output.\n");
      printDashes(PL_INFO, 0);

      sprintf(pString, "Which output to use? (1 - %d) ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      KSDensity *ksd = new KSDensity();
      psVector vecDataSet, vecXp, vecPp;
      if (nOutputs_ == 1) 
        vecDataSet.load(nSamples_,VecSamOutputs_.getDVector());
      else
      {
        vecDataSet.setLength(nSamples_);
        for (ss = 0; ss < nSamples_; ss++)
          vecDataSet[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      }
      ksd->genDensity1D(vecDataSet, vecXp, vecPp);
      delete ksd;

      //**/ write result to file
      if (plotMatlab()) strcpy(winput, "matlabkde.m");
      else              strcpy(winput, "scilabkde.sci");              
      fp = fopen(winput,"w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file %s\n", winput);
        cmdStatus = 1;
        continue;
      }
      fprintf(fp, "X = [\n"); 
      for (ss = 0; ss < vecXp.length(); ss++)
        fprintf(fp, "%e\n", vecXp[ss]);
      fprintf(fp, "];\n"); 
      fprintf(fp, "P = [\n"); 
      for (ss = 0; ss < vecPp.length(); ss++)
        fprintf(fp, "%e\n", vecPp[ss]);
      fprintf(fp, "];\n"); 
      fprintf(fp,"plot(X,P)\n");
      sprintf(pString,"Data Values");
      fwritePlotXLabel(fp, pString);
      sprintf(pString,"Probabilities");
      fwritePlotYLabel(fp, pString);
      fwritePlotAxes(fp);
      fclose(fp);

      //**/ compute statistics
      double kdemean = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        kdemean += VecSamOutputs_[ss*nOutputs_+outputID];
      kdemean /= (double) nSamples_;
      double kdevar = 0.0;
      for (ss = 0; ss < nSamples_; ss++)
        kdevar += pow(VecSamOutputs_[ss*nOutputs_+outputID]-kdemean,2.0);
      kdevar /= (double) nSamples_;
      printDashes(PL_INFO, 0);
      printf("Statistics from original sample:\n");
      printf("Sample mean    = %12.4e\n", kdemean);
      printf("Sample std dev = %12.4e\n", sqrt(kdevar));
      kdemean = 0.0;
      for (ss = 0; ss < vecPp.length(); ss++)
        kdemean += vecXp[ss] * vecPp[ss];
      kdevar = 0.0;
      for (ss = 0; ss < vecPp.length(); ss++)
        kdevar += pow(vecXp[ss] - kdemean, 2.0) * vecPp[ss];
      printDashes(PL_INFO, 0);
      printf("Statistics after kernel density estimation:\n");
      printf("Estimated mean     = %12.4e \n", kdemean);
      printf("Estimated std dev  = %12.4e \n", sqrt(kdevar));
      printDashes(PL_INFO, 0);
      printf("Kernel density information is now in %s.\n",winput);
      printf("The file contains bin values and probabilities.\n");
      printAsterisks(PL_INFO, 0);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ kde2 
    //**/ estimate mean and sd from a sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kde2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kde2: create PDF for each unique sample point with "); 
        printf("nondeterminstic output\n");
        printf("Syntax: kde2\n");
        printf("This command first divides the loaded sample into ");
        printf("blocks of equal size\n");
        printf("and then apply kde to each block. This command ");
        printf("is intended for models\n");
        printf("with heteroskechasticity when the sample can be ");
        printf("divided into equal-size\n");
        printf("blocks whereby each block has the same input ");
        printf("values but the outputs\n");
        printf("are different due to stochastic simulations, and ");
        printf("the distribution for\n");
        printf("each block is to be characterized.\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command uses kernel density estimation to compute ");
      printf("mean/std dev\n");
      printf("of each block of sub-samples in the loaded sample ");
      printf("with identical inputs\n");
      printf("(identical due to nondeterministic simulations).\n");
      printDashes(PL_INFO, 0);

      sprintf(pString, "Which output to use? (1 - %d) ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString, "Chunk size? (1 - %d) ", nSamples_/2);
      int    chunkSize = getInt(1, nSamples_/2, pString);
      chunkSize = nSamples_ / (nSamples_ / chunkSize);
      KSDensity *ksd = new KSDensity();
      psVector vecDataSet, chunkSet, vecXp, vecPp;
      vecDataSet.setLength(nSamples_);
      for (ss = 0; ss < nSamples_; ss++)
        vecDataSet[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
      chunkSet.setLength(chunkSize);

      //**/ write results to matlab or scilab file
      fp = fopen("kde2_results","w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open file kde2_results\n");
        cmdStatus = 1;
        continue;
      }
      fprintf(fp,"# The sample has been divided into %d chunks.\n",
              nSamples_/chunkSize);
      fprintf(fp,"# Kernel density estimation is applied to each.\n");
      fprintf(fp,"# First line has number of chunks and nInputs.\n");
      fprintf(fp,"# Each subsequent line begins with the sample\n");
      fprintf(fp,"# inputs of a chunk, followed by the estimated\n");
      fprintf(fp,"# mean and standard deviation of the chunk.\n");
      fprintf(fp,"%d %d 2\n", nSamples_/chunkSize, nInputs_);
      for (ss = 0; ss < nSamples_; ss+=chunkSize)
      {
        for (kk = 0; kk < chunkSize; kk++)
          chunkSet[kk] = vecDataSet[ss+kk];
        ksd->genDensity1D(chunkSet, vecXp, vecPp);
        for (ii = 0; ii < nInputs_; ii++)
          fprintf(fp, "%e ", VecSamInputs_[ss*nInputs_+ii]);
        ddata = 0.0;
        for (kk = 0; kk < vecXp.length(); kk++)
          ddata += vecXp[kk] * vecPp[kk];
        dtemp = 0.0;
        for (kk = 0; kk < vecXp.length(); kk++)
          dtemp += pow(vecXp[kk] - ddata, 2.0) * vecPp[kk];
        printf("Estimated mean,std dev = %e %e\n", ddata, sqrt(dtemp));
        fprintf(fp, "%16.8e %16.8e\n", ddata, sqrt(dtemp));
      }
      fclose(fp);
      printDashes(PL_INFO, 0);
      printf("kde2 information is now in kde2_results\n");
      printAsterisks(PL_INFO, 0);
      delete ksd;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ output_file 
    //**/ Changes the default name of the outputFile 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "output_file"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("output_file: changes the default output filename.\n");
        printf("Syntax: output_file <filename>\n");
        continue;
      }
      psOutputFilename_ = PSUADE_strdup(winput);
      printf("Default output file has been changed to %s\n",
             psOutputFilename_);
      if (psInputFilename_ == psOutputFilename_)
      {
        printf("WARNING: The output filename is the same as the input\n");
        printf("         filename. If you save the file it will\n");
        printf("         overwrite the input file.\n");
      }
    }

    //**/ -------------------------------------------------------------
    // +++ setupguide 
    //**/ a short guide to set up
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "setupguide"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("setupGuide: show how to set up an application\n");
        continue;
      }
      setupGuide();
    }

    //**/ -------------------------------------------------------------
    // +++ genworkflow 
    //**/ generate an application workflow
    //**/ -------------------------------------------------------------
    //**/else if (!strcmp(command, "genworkflow"))
    //**/{
    //**/  sscanf(lineIn,"%s %s",command,winput);
    //**/  if (!strcmp(winput, "-h"))
    //**/  {
    //**/    printf("genworkflow: generate a UQ workfow\n");
    //**/    continue;
    //**/  }
    //**/  printf("genworkflow not implemented yet.\n");
    //**/  cmdStatus = 0;
    //**/}

    //**/ -------------------------------------------------------------
    // +++ geninputfile 
    //**/ generate an input file for psuade 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "geninputfile"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("geninputfile: create a PSUADE input file\n");
        printf("Syntax: geninputfile (no argument needed)\n");
        continue;
      }
      status = genSetup(0, dataFile);
      VecSamInputs_.clean();
      VecSamOutputs_.clean();
      VecSamStates_.clean();
      VecILowerBs_.clean();
      VecIUpperBs_.clean();
      StrInpNames_.clean();
      StrOutNames_.clean();

      if (status == 0)
      {
        printf("PSUADE can create the sample input data file for you.\n");
        sprintf(pString,"Create the file ? (y or n) ");
        getString(pString, winput);
        if (winput[0] == 'y')
        {
          getInputFromFile(dataFile);
          cleanUp();
          psuadeIO_ = new PsuadeData();
          psuadeIO_->setOutputLevel(0);
          strcpy(winput, psInputFilename_);
          status = psuadeIO_->readPsuadeFile(winput);
          if (status != 0) 
          {
            printf("ERROR: cannot read file %s\n",winput);
            cmdStatus = 1;
            continue;
          }
          psuadeIO_->getParameter("input_ninputs", pPtr);
          nInputs_ = pPtr.intData_;
          pLower.clean();
          psuadeIO_->getParameter("input_lbounds", pLower);
          VecILowerBs_.load(nInputs_, pLower.dbleArray_);
          pLower.clean();
          pUpper.clean();
          psuadeIO_->getParameter("input_ubounds", pUpper);
          VecIUpperBs_.load(nInputs_, pUpper.dbleArray_);
          pUpper.clean();

          //**/ update input names
          pINames.clean();
          psuadeIO_->getParameter("input_names", pINames);
          StrInpNames_.setNumStrings(nInputs_+1);
          for (ii = 0; ii < nInputs_; ii++)
            StrInpNames_.loadOneString(ii, pINames.strArray_[ii]);

          psuadeIO_->getParameter("output_noutputs", pPtr);
          nOutputs_ = pPtr.intData_;
          pONames.clean();
          psuadeIO_->getParameter("output_names", pONames);
          StrOutNames_.setNumStrings(nOutputs_);
          for (ii = 0; ii < nOutputs_; ii++)
            StrOutNames_.loadOneString(ii, pONames.strArray_[ii]);

          psuadeIO_->getParameter("method_sampling", pPtr);
          SamplingMethod_ = pPtr.intData_;
          psuadeIO_->getParameter("method_nsamples", pPtr);
          nSamples_ = pPtr.intData_;
          psuadeIO_->getParameter("method_nreplications",pPtr);
          nReplications_ = pPtr.intData_;
          psuadeIO_->getParameter("input_sample", pPtr);
          VecSamInputs_.load(nSamples_*nInputs_, pPtr.dbleArray_);
          psuadeIO_->getParameter("output_sample", pPtr);
          VecSamOutputs_.load(nSamples_*nOutputs_,pPtr.dbleArray_);
          psuadeIO_->getParameter("output_states", pPtr);
          VecSamStates_.load(nSamples_,pPtr.intArray_);
          pINames.clean();
          pONames.clean();
          printf("==================================================\n");
          printf("The sample matrix is now stored in %s file.\n", 
                 psOutputFilename_);
          printf("You can also use genmars command to convert the\n");
          printf("sample matrix to a row-column format.\n");
          printf("==================================================\n");
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ genbatchfile 
    //**/ generate a LLNL-specific batch file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "genbatchfile"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("genbatchfile: create a batch file for ensemble runs\n");
        printf("Syntax: genbatchfile (no argument needed)\n");
        continue;
      }
      genBatchFile(0);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gendriver 
    //**/ generate an input file for psuade 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gendriver"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gendriver: create a driver file for runing PSUADE\n");
        printf("Syntax: gendriver (no argument needed)\n");
        continue;
      }
      genDriver(0);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gendist 
    //**/ generate data based on some distribution
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gendist"))
    {
      int    ns, gtype;
      double *sData, slbound, subound, smean, sstdev;
      psMatrix corMat;

      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gendist: create a sample using selected PDFs\n");
        printf("Syntax: gendist (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a sample of a user-specified ");
      printf("size based on the\n");
      printf("PDF presciption given by user (1D sample only).\n");
      printDashes(PL_INFO, 0);

      sprintf(pString, "Enter desired sample size = ");
      ns = getInt(10, 10000000, pString);
      sprintf(pString, 
        "PDF type = (1) N (2) L (3) T (4) Beta (5) Exp (6) Weibull: ");
      gtype = getInt(1, 6, pString);
      if      (gtype == 1) gtype = PSUADE_PDF_NORMAL;
      else if (gtype == 2) gtype = PSUADE_PDF_LOGNORMAL;
      else if (gtype == 3) gtype = PSUADE_PDF_TRIANGLE;
      else if (gtype == 4) gtype = PSUADE_PDF_BETA;
      else if (gtype == 5) gtype = PSUADE_PDF_EXPONENTIAL;
      else if (gtype == 6) gtype = PSUADE_PDF_WEIBULL;
      sprintf(pString, "PDF parameter 1 (e.g. mean for N) = ");
      smean = getDouble(pString);
      if (gtype != PSUADE_PDF_EXPONENTIAL) 
      {
        sprintf(pString, "PDF parameter 2 (e.g. std dev for N) = ");
        sstdev = getDouble(pString);
      }
      else sstdev = 0.0;
      corMat.setDim(1,1);
      corMat.setEntry(0,0, 1.0e0);
      PDFManager *pdfman = new PDFManager();
      pdfman->initialize(1, &gtype, &smean, &sstdev, corMat, NULL, NULL);
      vecOuts.setLength(ns);
      subound =  PSUADE_UNDEFINED;
      slbound = -PSUADE_UNDEFINED;
      if (gtype == PSUADE_PDF_LOGNORMAL) slbound = 0;
      vecUBs.load(1, &subound);
      vecLBs.load(1, &slbound);
      pdfman->genSample(ns, vecOuts, vecLBs, vecUBs);
      delete pdfman;

      sData = vecOuts.getDVector();
      fp = fopen("sample1D", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot write to file sample1D.\n");
        continue;
      }
      fprintf(fp, "%d 1\n", ns);
      for (ii = 0; ii < ns; ii++) fprintf(fp,"%9d %e\n",ii+1,sData[ii]);
      fclose(fp);
      printf("data file created in sample1D.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ genexample 
    //**/ generate an example
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "genexample"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("genexample: create a PSUADE example (driver+input file)\n");
        printf("Syntax: genexample (no argument needed)\n");
        continue;
      }
      genDriver(1);
      genSetup(1, winput);
      printf("Now use: psuade psuade.in to run the example.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ genconfigfile 
    //**/ generate an example Psuade config file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "genconfigfile"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("genconfigfile: create a PSUADE config template file\n");
        printf("Syntax: genconfigfile (no argument needed)\n");
        printf("Config files are used to modify settings.\n");
        printf("Config files may be used to replace interactive\n");
        printf("queries from PSUADE.\n");
        continue;
      }
      sprintf(pString,"Enter the name of the configure file to write to: ");
      getString(pString, dataFile);
      dataFile[strlen(dataFile)-1] = '\0';
      if (!strcmp(dataFile, "\0")) 
      {
        printf("ERROR: invalid file name.\n");
        cmdStatus = 1;
        continue;
      }
      status = genConfigFileTemplate(dataFile);
      if (status == 0) 
      {
        printf("genconfigfile completed - file name is %s.\n", dataFile);
        cmdStatus = 0;
      }
      else 
      {
        printf("ERROR: cannot write to file %s or no filename given.\n",
               dataFile);
        cmdStatus = 1;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ setconfigoption 
    //**/ set/modify a config option in the configuration table
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "setconfigoption"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("setconfigoption: set/modify a config table option\n");
        printf("Syntax: setconfigoption\n");
        printf("The internal configuration table is used to modify\n");
        printf("settings. \n");
        continue;
      }
      printf("Enter the new config string for the table : ");
      fgets(lineIn, 50000, stdin); 
      lineIn[strlen(lineIn)-1] = '\0';
      psConfig_.putParameter(lineIn);
      printf("%s\n", lineIn);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ showconfigtable 
    //**/ set/modify a config option in the configuration table
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "showconfigtable"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("showconfigtable: show configuration table content.\n");
        printf("Syntax: showconfigtable\n");
        printf("The internal configuration table is used to modify\n");
        printf("settings. \n");
        continue;
      }
      psConfig_.print();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // resetconfig
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "resetconfig"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("resetconfig: remove all contents of PSUADE's\n");
        printf("             internal config object.\n");
        printf("Syntax: resetconfig (no argument needed)\n");
        continue;
      }
      psConfig_.reset();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ chkjobs 
    //**/ check job status and generate a report
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "chkjobs"))
    {
      int  nPatterns, choice, nJobs1, nJobs2;
      char checkFile[200], subdirName[5001], dirName[5001];
      psStrings strPatterns;

      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("chkjobs: monitor the status of ensemble runs \n");
        printf("Syntax: chkjobs (no argument needed)\n");
        printf("Use this command with '-h' to see more details.\n");
        continue;
      }
      printf("PSUADE will check for job status in the current directory.\n");
      printf("So first make sure you are in the right directory.\n");
      printf("You should have subdirectories each of which is a job.\n");
      printf("The subdirectory names should be something like: \n");
      printf("<dir_prefix>xxx where xxx is a number from 1 to njobs.\n");
      sprintf(pString,"Enter the subdirectory prefix now : ");
      getString(pString, winput);
      sscanf(winput, "%s", dirName);
      sprintf(pString,
          "Enter another level of subdirectory, if any : (or NONE) ");
      getString(pString, winput);
      sscanf(winput, "%s", subdirName);
      printf("What is considered to be failed runs?\n");
      printf("(1) that a certain file does not exist\n");
      printf("(2) both (1) and that a pattern does not exist in a file\n");
      printf("(3) that certain pattern(s)s exist in a given file\n");
      sprintf(pString,"Make your selection : ");
      choice = getInt(1, 3, pString);
      nPatterns = 0;
      if (choice == 1)
      {
        sprintf(pString,"Enter the name of file that has to exist : ");
        getString(pString, winput);
        sscanf(winput, "%s", checkFile);
      }
      else
      {
        sprintf(pString,"Enter the name of file to check patterns: ");
        getString(pString, winput);
        sscanf(winput, "%s", checkFile);
        if (choice == 2) nPatterns = 1;
        else
        {
          sprintf(pString,"How many patterns (1 - 10)? ");
          nPatterns = getInt(1, 10, pString);
        }
        strPatterns.setNumStrings(nPatterns);
        for (ii = 0; ii < nPatterns; ii++)
        {
          printf("Enter pattern %d : ",ii+1);
          fgets(pString, 500, stdin); 
          strPatterns.loadOneString(ii, pString);
        }
      }
      sprintf(pString,"Enter the first job number to be probed (1 - ?): ");
      nJobs1 = getInt(1, 1000000, pString);
      sprintf(pString,"Enter the last  job number to be probed (%d - ?): ",
              nJobs1+1);
      nJobs2 = getInt(nJobs1+1, 1000000, pString);
      FILE *fErr = fopen("relaunchJobs.py", "w");
      if (fErr == NULL)
      {
        printf("ERROR: cannot open file to store job status info.\n");
        cmdStatus = 1;
        continue;
      }
      if (fErr != NULL)
      {
        fprintf(fErr, "import os\n");
        fprintf(fErr, "import sys\n\n");
        fprintf(fErr, "jobs = [");
      }
      count = 0;
      for (ii = nJobs1; ii <= nJobs2; ii++)
      {
        if (outputLevel_ > 0) printf("Processing job %d\n", ii);
        if (choice == 1)
        {
          if (strncmp(subdirName, "NONE", 4) == 0)
             sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
          else
             sprintf(winput,"%s%d/%s/%s",dirName,ii,subdirName,checkFile);
          fp = fopen(winput, "r");
          if (fp == NULL) 
          {
            if (fErr != NULL)
            {
              if (count == 0) fprintf(fErr, "%d", ii);
              else            fprintf(fErr, ",%d", ii);
            }
            else
            {
              printf("%s%d/%s fails (1): file does not exist.\n",
                     dirName, ii, checkFile);
            } 
            count++;
          }
          else fclose(fp);
        }
        if (choice == 2)
        {
          if (strncmp(subdirName, "NONE", 4) == 0)
             sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
          else
             sprintf(winput, "%s%d/%s/%s",dirName,ii,subdirName,checkFile);
          fp = fopen(winput, "r");
          if (fp == NULL) 
          {
            if (fErr != NULL)
            {
              if (count == 0) fprintf(fErr, "%d", ii);
              else            fprintf(fErr, ",%d", ii);
            }
            else
            {
              printf("%s%d fails (2a): file does not exist.\n",
                     dirName, ii);
            }
            count++;
          }
          else 
          {
            fclose(fp);
            sprintf(command, "grep \"%s\" %s > /dev/null",strPatterns[0], 
                    winput);
            status = system(command);
            if (status != 0)
            {
              if (fErr != NULL)
              {
                if (count == 0) fprintf(fErr, "%d", ii);
                else            fprintf(fErr, ",%d", ii);
              }
              else
              {
                printf("%s%d fails (2b): file does not exist.\n",dirName, ii);
              }
              count++;
            }
          }
        }
        if (choice == 3)
        {
          if (strncmp(subdirName, "NONE", 4) == 0)
             sprintf(winput, "%s%d/%s", dirName, ii, checkFile);
          else
             sprintf(winput, "%s%d/%s/%s", dirName, ii, subdirName, 
                     checkFile);
          fp = fopen(winput, "r");
          if (fp == NULL) 
          {
            if (fErr != NULL)
            {
              if (count == 0) fprintf(fErr, "%d", ii);
              else            fprintf(fErr, ",%d", ii);
            }
            else
            {
              printf("%s%d fails (3a): file does not exist.\n",dirName,
                     ii);
            }
            count++;
          }
          else 
          {
            fclose(fp);
            for (jj = 0; jj < nPatterns; jj++)
            {
              sprintf(command, "grep \"%s\" %s > /dev/null", 
                      strPatterns[jj], winput);
              status = system(command);
              if (status == 0)
              {
                if (fErr != NULL)
                {
                  if (count == 0) fprintf(fErr, "%d", ii);
                  else            fprintf(fErr, ",%d", ii);
                }
                else
                {
                  printf("%s%d fails (3b): file does not exist.\n",
                         dirName, ii);
                }
                count++;
                break;
              }
            }
          }
        }
      }
      if (fErr != NULL)
      {
        fprintf(fErr,"]");
        fprintf(fErr,"\n");
        fprintf(fErr,"for index in jobs:\n");
        fprintf(fErr,"#  insert your clean-up procedure\n");
        fprintf(fErr,"## e.g. cmd = \"/bin/rm -r workdir.\" + str(index)\n");
        fprintf(fErr,"##      os.system(cmd)\n\n");
        fprintf(fErr,"#  insert your re-start procedure\n");
        fprintf(fErr,
             "## e.g. cmd = \"driver.py psuadeApps_ct.in.\" + str(index)\n");
        fprintf(fErr,"psuadeApps_ct.out.\" + str(index)\n");
        fprintf(fErr,"##      os.system(cmd)\n\n");
        fprintf(fErr,"#  insert your job submission procedure\n");
        fprintf(fErr,"## cmd = \"cd %s.\" + str(index) + \"; \"\n",dirName);
        fprintf(fErr,
             "## cmd = cmd + \"/usr/bin/psub batchFile.\" +str(index)\n");
        fprintf(fErr,"## os.system(cmd)\n\n");
        fclose(fErr);
        printf("The failed jobs are in the file relaunchJobs.py.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ list1 
    //**/ list one input and one output pair
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "list1"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("list1: list all points with 1 input and 1 output\n");
        printf("Syntax: list1 (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command displays all sample points with one ");
      printf("selected input and 1\n");
      printf("output, optionally sorted based on the input ");
      printf("or output values.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString,"Sort the input (y or n) ? ");
      getString(pString, winput);
      if (winput[0] == 'n') 
      {
        sprintf(pString,"Sort the output (y or n) ? ");
        getString(pString, winput);
        if (winput[0] == 'n') winput[0] = 'N'; 
        if (winput[0] == 'y') winput[0] = 'Y'; 
      }
      vecXT.setLength(nSamples_);
      vecYT.setLength(nSamples_);
      vecWT.setLength(nSamples_);
      double *dPtrX = vecXT.getDVector();
      double *dPtrY = vecYT.getDVector();
      double *dPtrW = vecWT.getDVector();
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        dPtrX[sInd] = VecSamInputs_[sInd*nInputs_+iInd];
        dPtrY[sInd] = VecSamOutputs_[sInd*nOutputs_+outputID];
        dPtrW[sInd] = (double) sInd + 1;
      }
      if (winput[0] == 'y') 
        sortDbleList3(nSamples_, dPtrX, dPtrY, dPtrW);
      else if (winput[0] == 'Y') 
        sortDbleList3(nSamples_, dPtrY, dPtrX, dPtrW);
      for (sInd = 0; sInd < nSamples_; sInd++)
        printf("%6d: Sample %7d : input = %16.8e, output = %16.8e\n",
               sInd+1, (int) dPtrW[sInd], dPtrX[sInd], dPtrY[sInd]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ list2 
    //**/ list two input and one output pair
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "list2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("list2: list all points with 2 inputs and 1 output\n");
        printf("Syntax: list2 (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command displays all sample points with two ");
      printf("selected inputs and 1\n");
      printf("output, optionally sorted based on the input ");
      printf("or output values.\n");
      printDashes(PL_INFO, 0);
 
      sprintf(pString, "Enter input number 1 (1 - %d) : ", nInputs_);
      iInd1 = getInt(1, nInputs_, pString);
      iInd1--;
      sprintf(pString, "Enter input number 2 (1 - %d) : ", nInputs_);
      iInd2 = getInt(1, nInputs_, pString);
      iInd2--;
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString,"Sort the inputs (y or n) ? ");
      getString(pString, winput);
      if (winput[0] == 'n') 
      {
        sprintf(pString,"Sort the output (y or n) ? ");
        getString(pString, winput);
        if (winput[0] == 'n') winput[0] = 'N'; 
        if (winput[0] == 'y') winput[0] = 'Y'; 
      }
      vecXT.setLength(nSamples_);
      vecYT.setLength(nSamples_);
      vecWT.setLength(nSamples_);
      vecVT.setLength(nSamples_);
      double *dPtrX = vecXT.getDVector();
      double *dPtrY = vecYT.getDVector();
      double *dPtrW = vecWT.getDVector();
      double *dPtrV = vecVT.getDVector();
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        dPtrX[sInd] = VecSamInputs_[sInd*nInputs_+iInd1];
        dPtrV[sInd] = VecSamInputs_[sInd*nInputs_+iInd2];
        dPtrY[sInd] = VecSamOutputs_[sInd*nOutputs_+outputID];
        dPtrW[sInd] = (double) sInd;
      }
      if (winput[0] == 'y') 
      {
        sortDbleList4(nSamples_, dPtrX, dPtrV, dPtrY, dPtrW);
        sInd = 0;
        count = 1; 
        while (sInd < nSamples_) 
        {
          sInd++;
          if (dPtrX[sInd] == dPtrX[sInd-1]) count++;
          else
          {
            if (count > 1)
            {
              sortDbleList4(count, &(dPtrV[sInd-count]), 
                    &(dPtrX[sInd-count]), &(dPtrW[sInd-count]),
                    &(dPtrY[sInd-count]));
            }
            count = 1;
          }
        }
        if (count > 1)
        {
          sortDbleList4(count, &(dPtrV[sInd-count]),
                 &(dPtrX[sInd-count]), &(dPtrW[sInd-count]),
                 &(dPtrY[sInd-count]));
        }
      }
      else if (winput[0] == 'Y') 
      {
        sortDbleList4(nSamples_, dPtrY, dPtrX, dPtrV, dPtrW);
      }
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        printf("%6d: Sample %7d : ",sInd+1,((int) dPtrW[sInd])+1);
        printf("inputs = (%12.4e, %12.4e), output = %12.4e\n",
               dPtrX[sInd], dPtrV[sInd], dPtrY[sInd]);
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ listall 
    //**/ list all inputs and one output 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "listall"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("listall: list all points with all inputs and 1 output\n");
        printf("Syntax: listall (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command displays all sample points with all ");
      printf("inputs and 1 selected\n");
      printf("output, optionally sorted based on the output values.\n");
      printDashes(PL_INFO, 0);

      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      sprintf(pString,"Sort the inputs (y or n) ? ");
      getString(pString, winput);
      if (winput[0] == 'n')
      {
        sprintf(pString,"Sort the output (y or n) ? ");
        getString(pString, winput);
        if (winput[0] == 'y') winput[0] = 'Y';
      }
      if (winput[0] == 'y')
      {
        psMatrix matTemp;
        matTemp.setFormat(PS_MAT2D);
        matTemp.setDim(nSamples_,nInputs_+1);
        vecIT.setLength(nSamples_);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            matTemp.setEntry(sInd, ii, VecSamInputs_[sInd*nInputs_+ii]);
          matTemp.setEntry(sInd,nInputs_,VecSamOutputs_[sInd*nOutputs_+outputID]);
        }
        sortnDbleList(nSamples_,nInputs_+1,matTemp.getMatrix2D(),
                      vecIT.getIVector());
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          printf("Sample %7d : ", sInd+1);
          for (ii = 0; ii < nInputs_; ii++)
            printf("%12.4e ", matTemp.getEntry(sInd, ii));
          printf("%12.4e\n", matTemp.getEntry(sInd, nInputs_));
        }
      }
      else if (winput[0] == 'Y')
      {
        vecIT.setLength(nSamples_);
        vecYT.setLength(nSamples_);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          vecIT[sInd] = sInd;
          vecYT[sInd] = VecSamOutputs_[sInd*nOutputs_+outputID];
        }
        sortDbleList2a(nSamples_, vecYT.getDVector(), vecIT.getIVector());
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          kk = vecIT[sInd];
          printf("%6d: Sample %7d : ", sInd+1, kk+1);
          for (ii = 0; ii < nInputs_; ii++)
            printf("%12.4e ", VecSamInputs_[kk*nInputs_+ii]);
          printf("= %12.4e\n", vecYT[sInd*nOutputs_+outputID]);
        }
      }
      else
      {
        for (ss = 0; ss < nSamples_; ss++)
        {
          printf("Sample %7d : ", ss+1);
          for (ii = 0; ii < nInputs_; ii++)
            printf("%12.4e ", VecSamInputs_[ss*nInputs_+ii]);
          printf("= %12.4e\n", VecSamOutputs_[ss*nOutputs_+outputID]);
        }
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ isort 
    //**/ sort sample by inputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "isort"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("isort: sort the sample based on a selected input\n");
        printf("Syntax: isort (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command sorts the loaded sample based on the ");
      printf("values of a selected\n");
      printf("input.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ",nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;

      vecWT.setLength(nSamples_);
      vecVT.setLength(nSamples_);
      vecIT.setLength(nSamples_);
      //**/ sort the first input 
      for (ss = 0; ss < nSamples_; ss++)
      {
        vecWT[ss] = VecSamInputs_[ss*nInputs_+iInd];
        vecIT[ss] = ss;
      }
      sortDbleList2a(nSamples_,vecWT.getDVector(),vecIT.getIVector());
      //**/ copy the re-ordered inputs back
      for (ss = 0; ss < nSamples_; ss++)
        VecSamInputs_[ss*nInputs_+iInd] = vecWT[ss];
      //**/ apply the reordering to all the other inputs
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (ii != iInd)
        {
          for (ss = 0; ss < nSamples_; ss++)
            vecWT[ss] = VecSamInputs_[ss*nInputs_+ii];
          for (ss = 0; ss < nSamples_; ss++)
          {
            kk = vecIT[ss];
            VecSamInputs_[ss*nInputs_+ii] = vecWT[kk];
          }
        }
      }
      //**/ apply the reordering to all outputs
      for (ii = 0; ii < nOutputs_; ii++)
      {
        for (ss = 0; ss < nSamples_; ss++)
          vecWT[ss] = VecSamOutputs_[ss*nOutputs_+ii];
        for (ss = 0; ss < nSamples_; ss++)
        {
          kk = vecIT[ss];
          VecSamOutputs_[ss*nOutputs_+ii] = vecWT[kk];
        }
      }
      cmdStatus = 0;
      printf("isort completed. Use 'write' to store the modified sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ osort 
    //**/ sort sample by outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "osort"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("osort: sort the sample based on outputs\n");
        printf("Syntax: osort (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command sorts the loaded sample based on ");
      printf("its output values.");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      vecWT.setLength(nSamples_);
      vecIT.setLength(nSamples_);
      //**/ sort the output 
      for (ss = 0; ss < nSamples_; ss++)
      {
        vecWT[ss] = VecSamOutputs_[ss*nOutputs_+outputID];
        vecIT[ss] = ss;
      }
      sortDbleList2a(nSamples_,vecWT.getDVector(),vecIT.getIVector());
      //**/ apply the reordering to all inputs
      for (ii = 0; ii < nInputs_; ii++)
      {
        for (ss = 0; ss < nSamples_; ss++)
          vecWT[ss] = VecSamInputs_[ss*nInputs_+ii];
        for (ss = 0; ss < nSamples_; ss++)
        {
          kk = vecIT[ss];
          VecSamInputs_[ss*nInputs_+ii] = vecWT[kk];
        }
      }
      //**/ apply the reordering to all outputs
      for (ii = 0; ii < nOutputs_; ii++)
      {
        for (ss = 0; ss < nSamples_; ss++)
          vecWT[ss] = VecSamOutputs_[ss*nOutputs_+ii];
        for (ss = 0; ss < nSamples_; ss++)
        {
          kk = vecIT[ss];
          VecSamOutputs_[ss*nOutputs_+ii] = vecWT[kk];
        }
      }
      cmdStatus = 0;
      printf("osort completed. Use 'write' to store the modified sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ disp_sample or sshow
    //**/ list one sample point 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "disp_sample") || !strcmp(command, "sshow"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("disp_sample (or sshow): display one sample point\n");
        printf("Syntax: disp_sample or sshow (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("This command shows the values of a selected sample point.\n");
      sprintf(pString, "Enter sample number (1 - %d) : ", nSamples_);
      ind = getInt(1, nSamples_, pString);
      ind--;
      printf("Sample %7d : \n", ind+1);
      for (ii = 0; ii < nInputs_; ii++)
        printf("   input  %3d = %16.8e\n", ii+1, 
               VecSamInputs_[ind*nInputs_+ii]);
      for (ii = 0; ii < nOutputs_; ii++)
        printf("   output %3d = %16.8e\n", ii+1, 
               VecSamOutputs_[ind*nOutputs_+ii]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ imax 
    //**/ find the maximum value of a selected input
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "imax"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("imax: find the sample point with maximum input value\n");
        printf("Syntax: imax (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command finds maximum value of a selected input.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      int inputID = getInt(1, nInputs_, pString);
      inputID--;
      double maxX = VecSamInputs_[inputID];
      int    maxI = 0;
      for (sInd = 1; sInd < nSamples_; sInd++)
      {
        if (VecSamInputs_[sInd*nInputs_+inputID] > maxX)
        {
          maxI = sInd;
          maxX = VecSamInputs_[sInd*nInputs_+inputID];
        }
      }
      printf("Sample %d gives maximum for input %d = %e\n",
             maxI+1,inputID+1,VecSamInputs_[maxI*nInputs_+inputID]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ omax 
    //**/ find the maximum value of a selected output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "omax"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("omax: find the sample point with maximum output value\n");
        printf("Syntax: omax (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command computes the maximum value of a selected output.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      double maxY = VecSamOutputs_[outputID];
      int    maxI = 0;
      for (sInd = 1; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > maxY)
        {
          maxI = sInd;
          maxY = VecSamOutputs_[sInd*nOutputs_+outputID];
        }
      }
      printf("Sample %d gives maximum for output %d = %e\n",
             maxI+1,outputID+1,
             VecSamOutputs_[maxI*nOutputs_+outputID]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ imin 
    //**/ find the minimum value of a selected input
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "imin"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("imax: find the sample point with minimum input value\n");
        printf("Syntax: imin (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command finds minimum value of a selected input.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      int inputID = getInt(1, nInputs_, pString);
      inputID--;
      double minX = VecSamInputs_[inputID];
      int    minI = 0;
      for (sInd = 1; sInd < nSamples_; sInd++)
      {
        if (VecSamInputs_[sInd*nInputs_+inputID] < minX)
        {
          minI = sInd;
          minX = VecSamInputs_[sInd*nInputs_+inputID];
        }
      }
      printf("Sample %d gives maximum for input %d = %e\n",
             minI+1,inputID+1,VecSamInputs_[minI*nInputs_+inputID]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ omin 
    //**/ find the minimum value of a selected output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "omin"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("omin: find the sample point with minimum output value\n");
        printf("Syntax: omin (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command finds the minimum value of a selected output.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      double minY = VecSamOutputs_[outputID];
      int    minI = 0;
      for (sInd = 1; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < minY)
        {
          minI = sInd;
          minY = VecSamOutputs_[sInd*nOutputs_+outputID];
        }
      }
      printf("Sample %d gives minimum for output %d = %e\n",
             minI+1,outputID+1,
             VecSamOutputs_[minI*nOutputs_+outputID]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ onorm 
    //**/ find 2-norm of a sample output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "onorm"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("onorm: compute the 2-norm of a sample output\n");
        printf("Syntax: onorm (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command computes the vector 2-norm of a selected output.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      double onorm = 0.0;
      for (sInd = 0; sInd < nSamples_; sInd++)
        onorm += pow(VecSamOutputs_[sInd*nOutputs_+outputID], 2.0);
      printf("2-norm of output %d = %e\n",outputID+1,sqrt(onorm));
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ osum 
    //**/ find sample sum of a sample output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "osum"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("osum: compute the sample sum of a sample output\n");
        printf("Syntax: osum (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command computes the vector sum of a selected output.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      double osum = 0.0;
      for (sInd = 0; sInd < nSamples_; sInd++)
        osum += VecSamOutputs_[sInd*nOutputs_+outputID];
      printf("Accumulative sum of output %d = %e\n", outputID+1,osum);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ inormalize 
    //**/ normalize a sample input
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "inormalize"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("inormalize: normalize a sample input\n");
        printf("Syntax: inormalize (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command re-scales a selected input ");
      printf("(restricted form of 'iscale').\n");
      printf("You are to select between the following:\n");
      printf("1. normalize by subtracting mean and divide by std dev\n");
      printf("2. normalize by scaling to [0,1]\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ",nInputs_);
      ind = getInt(1, nInputs_, pString);
      ind--;
      sprintf(pString," Which option (1 or 2) : ");
      kk = getInt(1, 2, pString);
      if (kk == 1)
      {
        double dmean = 0.0;
        for (sInd = 0; sInd < nSamples_; sInd++)
          dmean += VecSamInputs_[sInd*nInputs_+ind];
        printf("Input mean = %e\n",dmean);
        dmean = dmean / (double) nSamples_;
        double dstdv = 0.0;
        for (sInd = 0; sInd < nSamples_; sInd++)
          dstdv += pow(VecSamInputs_[sInd*nInputs_+ind]-dmean,2.0);
        dstdv = sqrt(dstdv/(nSamples_-1));
        printf("Input s.d. = %e\n",dstdv);
        if (dstdv == 0.0) 
          printf("Standard deviation = 0 ==> no normalization.\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
          VecSamInputs_[sInd*nInputs_+ind] = 
             (VecSamInputs_[sInd*nInputs_+ind] - dmean) / dstdv;
      }
      else
      {
        double ddmax=-PSUADE_UNDEFINED;
        double ddmin= PSUADE_UNDEFINED;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          ddata = VecSamInputs_[sInd*nInputs_+ind];
          if (ddata > ddmax) ddmax = ddata;
          if (ddata < ddmin) ddmin = ddata;
        }
        if (ddmax == ddmin)
          printf("INFO: all values are the same ==> no scaling.\n");
        else
        {
          for (sInd = 0; sInd < nSamples_; sInd++)
            VecSamInputs_[sInd*nInputs_+ind] =
              (VecSamInputs_[sInd*nInputs_+ind]-ddmin)/(ddmax-ddmin);
        }
      }
      if (VecILowerBs_.length() > 0)
      {
        VecILowerBs_[ind] = PSUADE_UNDEFINED;
        VecIUpperBs_[ind] = -PSUADE_UNDEFINED;
        VecInpPDFs_[ind] = 0;
        VecInpMeans_[ind] = 0;
        VecInpStds_[ind] = 0;
        for (ss = 0; ss < nSamples_; ss++)
        {
          if (VecSamInputs_[ss*nInputs_+ind] < VecILowerBs_[ind])
            VecILowerBs_[ind] = VecSamInputs_[ss*nInputs_+ind];
          if (VecSamInputs_[ss*nInputs_+ind] > VecIUpperBs_[ind])
            VecIUpperBs_[ind] = VecSamInputs_[ss*nInputs_+ind];
        }
      }
      if (inputCMat_ != NULL)
      {
        for (jj = 0; jj < nInputs_; jj++)
        {
          inputCMat_->setEntry(jj,ind,0.0);
          inputCMat_->setEntry(ind,jj,0.0);
        }
        inputCMat_->setEntry(ind,ind,1.0);
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      cmdStatus = 0;
      printf("inormalize completed. Use 'write' to store modified sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ onormalize 
    //**/ normalize a sample output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "onormalize"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("onormalize: normalize a sample output\n");
        printf("Syntax: onormalize (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command re-scales a selected output ");
      printf("(restricted form of 'oscale').\n");
      printf("You are to select between the following:\n");
      printf("1. normalize by subtracting mean and divide by std dev\n");
      printf("2. normalize by scaling to [0,1]\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      sprintf(pString," Which option (1 or 2) : ");
      kk = getInt(1, 2, pString);
      if (kk == 1)
      {
        double dmean = 0.0;
        for (sInd = 0; sInd < nSamples_; sInd++)
          dmean += VecSamOutputs_[sInd*nOutputs_+outputID];
        printf("Output mean = %e\n",dmean);
        dmean = dmean / (double) nSamples_;
        double dstdv = 0.0;
        for (sInd = 0; sInd < nSamples_; sInd++)
          dstdv += pow(VecSamOutputs_[sInd*nOutputs_+outputID]-dmean,
                       2.0);
        dstdv = sqrt(dstdv/(nSamples_-1));
        printf("Output s.d. = %e\n",dstdv);
        if (dstdv == 0.0) 
          printf("Standard deviation = 0 ==> no normalization.\n");
        for (sInd = 0; sInd < nSamples_; sInd++)
          VecSamOutputs_[sInd*nOutputs_+outputID] = 
             (VecSamOutputs_[sInd*nOutputs_+outputID]-dmean)/dstdv;
      }
      else
      {
        double ddmax=-PSUADE_UNDEFINED;
        double ddmin= PSUADE_UNDEFINED;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          ddata = VecSamOutputs_[sInd*nOutputs_+outputID];
          if (ddata > ddmax) ddmax = ddata;
          if (ddata < ddmin) ddmin = ddata;
        }
        if (ddmax == ddmin)
          printf("INFO: all values are the same ==> no scaling.\n");
        else
        {
          for (sInd = 0; sInd < nSamples_; sInd++)
            VecSamOutputs_[sInd*nOutputs_+outputID] =
              (VecSamOutputs_[sInd*nOutputs_+outputID]-ddmin)/
                (ddmax-ddmin);
        }
      }
      cmdStatus = 0;
      printf("onormalize completed. Use 'write' to store modified sample.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ idelete/ikeep
    //**/ remove/keep inputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "idelete") ||
             !strcmp(command, "ikeep"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        
        if (!strcmp(command, "idelete"))
        {
          printf("idelete: delete some inputs in the loaded sample\n");
          printf("Syntax: idelete (no argument needed)\n");
        }
        else
        {
          printf("ikeep: keep only some inputs in the loaded sample\n");
          printf("Syntax: ikeep (no argument needed)\n");
        }
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      if (!strcmp(command, "idelete"))
      {
        printf("This command deletes a subset of inputs ");
        printf("from the loaded sample.\n");
      }
      else
      {
        printf("This command keeps a subset of inputs in the ");
        printf("loaded sample and remove\n");
        printf("the rest.\n");
      }
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ == 1)
      {
        printf("You have only one input left -> no action taken.\n");
        cmdStatus = 1;
        continue;
      }
      printf("The current set of inputs are:\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %3d = %s\n", ii+1, StrInpNames_[ii]);
          
      //**/ collect information from user
      vecIT.setLength(nInputs_);
      int inpCnt = 0;
      if (!strcmp(command, "idelete"))
      {
        for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = 1;
        sprintf(pString,"How many inputs to delete? (1-%d) ",nInputs_-1);
      }
      else
      {
        for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = 0;
        sprintf(pString,"How many inputs to keep? (1-%d) ",nInputs_-1);
      }
      inpCnt = getInt(1, nInputs_-1, pString);
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      for (ii = 0; ii < inpCnt; ii++)
      {
        iInd = getInt(1, nInputs_, pString);
        if (!strcmp(command, "idelete"))
        {
          printf("You are removing input %d (%s)\n",iInd,
                 StrInpNames_[iInd-1]);
          vecIT[iInd-1] = 0;
        }
        else
        {
          printf("You are keeping input %d (%s)\n",iInd,
                 StrInpNames_[iInd-1]);
          vecIT[iInd-1] = 1;
        }
      }

      //**/ remove inputs from bounds 
      inpCnt = vecIT.sum();
      vecWT = VecILowerBs_;
      vecVT = VecIUpperBs_;
      VecILowerBs_.setLength(inpCnt);
      VecIUpperBs_.setLength(inpCnt);
      inpCnt = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecIT[ii] == 1)
        {
          VecILowerBs_[inpCnt] = vecWT[ii];
          VecIUpperBs_[inpCnt] = vecVT[ii];
          inpCnt++;
        }
      }

      //**/ remove input from sample input matrix
      vecXT = VecSamInputs_;
      VecSamInputs_.setLength(inpCnt*nSamples_);
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        kk = 0;
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (vecIT[ii] == 1)
          {
            VecSamInputs_[sInd*inpCnt+kk] = vecXT[sInd*nInputs_+ii];
            kk++;
          }
        }
      }
      vecXT.clean();

      //**/ remove input from input names
      count = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecIT[ii] == 0)
        {
          StrInpNames_.removeOneString(ii-count);
          count++;
        }
      }

      //**/ remove input from input PDF information
      vecIT = VecInpPDFs_;
      vecWT = VecInpMeans_;
      vecVT = VecInpStds_;
      VecInpPDFs_.setLength(inpCnt);
      VecInpMeans_.setLength(inpCnt);
      VecInpStds_.setLength(inpCnt);
      kk = 0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecIT[ii] == 1)
        {
          VecInpPDFs_[kk] = vecIT[ii];
          VecInpMeans_[kk] = vecWT[ii];
          VecInpStds_[kk]  = vecVT[ii];
          if (inputCMat_ != NULL)
          {
            for (jj = 0; jj < nInputs_; jj++)
            {
              ddata = inputCMat_->getEntry(ii,jj);
              inputCMat_->setEntry(kk,jj,ddata);
            }
            for (jj = 0; jj < nInputs_; jj++)
            {
              ddata = inputCMat_->getEntry(jj,ii);
              inputCMat_->setEntry(jj,kk,ddata);
            }
          }
          kk++;
        }
      }

      //**/ remove input from input correlation matrix
      if (inputCMat_ != NULL)
      {
        psMatrix *tmpMat = new psMatrix();
        tmpMat->setDim(inpCnt,inpCnt);
        for (ii = 0; ii < inpCnt; ii++)
        {
          for (jj = 0; jj < inpCnt; jj++)
          {
             ddata = inputCMat_->getEntry(ii,jj);
             tmpMat->setEntry(ii,jj,ddata);
          }
        }
        delete inputCMat_;
        inputCMat_ = tmpMat;
      }

      //**/ update data base and session
      nInputs_ = inpCnt;
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                  VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                  VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                  VecInpPDFs_.getIVector(),VecInpMeans_.getDVector(),
                  VecInpStds_.getDVector(), inputCMat_); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      if (!strcmp(command, "idelete"))
        printf("idelete completed. Use 'write' to store.\n");
      else
        printf("ikeep completed. Use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ odelete/okeep 
    //**/ remove/keep selected outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "odelete") ||
             !strcmp(command, "okeep"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        if (!strcmp(command, "odelete"))
        {
          printf("odelete: delete a subset of outputs from the ");
          printf("loaded sample\n");
          printf("Syntax: odelete (no argument needed)\n");
        }
        else
        {
          printf("okeep: keep a subset of outputs in the loaded sample\n");
          printf("Syntax: okeep (no argument needed)\n");
        }
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      if (!strcmp(command, "odelete"))
      {
        printf("This command deletes a subset of outputs ");
        printf("from the loaded sample.\n");
      }
      else
      {
        printf("This command keeps a subset of outputs in the ");
        printf("loaded sample and\n");
        printf("removes the rest.\n");
      }
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nOutputs_ <= 1)
      {
        printf("You have only one output left -> no action taken.\n");
        cmdStatus = 1;
        continue;
      }
      printf("The current set of outputs are:\n");
      for (ii = 0; ii < nOutputs_; ii++)
        if (StrOutNames_[ii] != NULL)
          printf("Output %3d = %s\n", ii+1, StrOutNames_[ii]);

      //**/ collect information
      vecIT.setLength(nOutputs_);
      int outCnt = 0;
      if (!strcmp(command, "odelete"))
      {
        for (ii = 0; ii < nOutputs_; ii++) vecIT[ii] = 1;
        sprintf(pString,"How many outputs to delete? (1-%d) ",nOutputs_-1);
      }
      else
      {
        for (ii = 0; ii < nOutputs_; ii++) vecIT[ii] = 0;
        sprintf(pString,"How many outputs to keep? (1-%d) ",nOutputs_-1);
      }
      outCnt = getInt(1, nOutputs_-1, pString);
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      for (ii = 0; ii < outCnt; ii++)
      {
        iInd = getInt(1, nOutputs_, pString);
        if (!strcmp(command, "odelete"))
        {
          printf("You are removing output %d (%s)\n",iInd,
                 StrOutNames_[iInd-1]);
          vecIT[iInd-1] = 0;
        }
        else
        {
          printf("You are keeping output %d (%s)\n",iInd,
                 StrOutNames_[iInd-1]);
          vecIT[iInd-1] = 1;
        }
      }

      //**/ remove output from input names
      for (ii = nOutputs_-1; ii >= 0; ii--)
      {
        if (vecIT[ii] == 0)
        {
          printf("Removing output %s\n",StrOutNames_[ii]);
          StrOutNames_.removeOneString(ii);
        }
      }

      //**/ remove output from sample output matrix
      if (!strcmp(command, "odelete")) 
        outCnt = nOutputs_ - outCnt;
      vecYT = VecSamOutputs_;
      VecSamOutputs_.setLength(outCnt*nSamples_);
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        kk = 0;
        for (ii = 0; ii < nOutputs_; ii++)
        {
          if (vecIT[ii] == 1)
          {
            VecSamOutputs_[sInd*outCnt+kk] = vecYT[sInd*nOutputs_+ii];
            kk++;
          }
        }
      }
      vecYT.clean();

      //**/ update store
      nOutputs_ = outCnt;
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
           VecSamOutputs_.getDVector(),
           VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      if (!strcmp(command, "odelete"))
        printf("odelete completed. Use 'write' to store.\n");
      else
        printf("okeep completed. Use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sdelete 
    //**/ remove selected sample points 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sdelete"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sdelete: delete a sample point \n");
        printf("Syntax: sdelete (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command deletes a subset of points ");
      printf("from the loaded sample.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nSamples_ == 1)
      {
        printf("You have only one sample point left -> no deletion.\n");
        continue;
      }
      printf("Options are: \n");
      printf("1. Delete a range of sample points (e.g. [11-20])\n");
      printf("2. Delete a noncontiguous subset (e.g. [8 14 17 22 ..]\n");
      sprintf(pString,"Which option ? (1 or 2 ) ");
      int option = getInt(1, 2, pString);
      
      int rl, ru;
      vecIT.setLength(nSamples_);
      if (option == 1)
      {
        printf("Select sample range (lower and upper bounds) :\n");
        sprintf(pString,"Lower bound ? (1 - %d) ", nSamples_-2);
        rl = getInt(1, nSamples_-2, pString);
        if (rl > 1)
        {
          sprintf(pString,"Upper bound ? (%d - %d) ",rl,nSamples_);
          ru = getInt(1, nSamples_, pString);
        }
        else
        {
          sprintf(pString,"Upper bound ? (%d - %d) ",rl,nSamples_-1);
          ru = getInt(1, nSamples_-1, pString);
        }
        for (ii = rl-1; ii < ru; ii++) vecIT[ii] = 1;
      }
      else
      {
        sprintf(pString,"Enter sample number to delete (0 when done) : ");
        while (vecIT.sum() < nSamples_-1)
        {
          kk = getInt(0, nSamples_, pString);
          if (kk == 0) break;
          vecIT[kk-1] = 1; 
        }
      }
      count = 0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (vecIT[ss] == 0)
        {
          for (ii = 0; ii < nInputs_; ii++)
            VecSamInputs_[count*nInputs_+ii] = 
                         VecSamInputs_[ss*nInputs_+ii];
          for (jj = 0; jj < nOutputs_; jj++)
            VecSamOutputs_[count*nOutputs_+jj] = 
                         VecSamOutputs_[ss*nOutputs_+jj];
          VecSamStates_[count] = VecSamStates_[ss]; 
          count++;
        }
      }
      nSamples_ = count;
      //**/ update the database and session
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                        NULL,NULL,VecSamInputs_.getDVector(),NULL, 
                        NULL, NULL, NULL, NULL);
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
            VecSamOutputs_.getDVector(),
            VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      psuadeIO_->updateMethodSection(-1,nSamples_,-1,-1,-1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("sdelete completed. Use 'write' to store.\n");
      printf("Use 'sinfo' to get updated sample information.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ skeep 
    //**/ keep selected sample points 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "skeep"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("skeep: delete sample points \n");
        printf("Syntax: skeep (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command keeps a subset of points in the ");
      printf("loaded sample and removes\n");
      printf("the others.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nSamples_ == 1)
      {
        printf("You have only one sample point left -> no deletion.\n");
        continue;
      }
      printf("Options are: \n");
      printf("1. Keep a range of sample points (e.g. [11-20])\n");
      printf("2. Keep a noncontiguous subset (e.g. [8 14 17 22 ..]\n");
      sprintf(pString,"Which option ? (1 or 2 ) ");
      int option = getInt(1, 2, pString);
      
      int rl, ru;
      vecIT.setLength(nSamples_);
      if (option == 1)
      {
        printf("Select sample range (lower and upper bounds) :\n");
        sprintf(pString,"Lower bound ? (1 - %d) ", nSamples_-2);
        rl = getInt(1, nSamples_-2, pString);
        if (rl > 1)
        {
          sprintf(pString,"Upper bound ? (%d - %d) ",rl,nSamples_);
          ru = getInt(1, nSamples_, pString);
        }
        else
        {
          sprintf(pString,"Upper bound ? (%d - %d) ",rl,nSamples_-1);
          ru = getInt(1, nSamples_-1, pString);
        }
        for (ii = rl-1; ii < ru; ii++) vecIT[ii] = 1;
      }
      else
      {
        sprintf(pString,"Enter sample number to keep (0 when done) : ");
        while (vecIT.sum() < nSamples_-1)
        {
          kk = getInt(0, nSamples_, pString);
          if (kk == 0) break;
          vecIT[kk-1] = 1; 
        }
      }
      count = 0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (vecIT[ss] == 1)
        {
          for (ii = 0; ii < nInputs_; ii++)
            VecSamInputs_[count*nInputs_+ii] = 
                         VecSamInputs_[ss*nInputs_+ii];
          for (jj = 0; jj < nOutputs_; jj++)
            VecSamOutputs_[count*nOutputs_+jj] = 
                         VecSamOutputs_[ss*nOutputs_+jj];
          VecSamStates_[count] = VecSamStates_[ss]; 
          count++;
        }
      }
      nSamples_ = count;
      //**/ update the database and session
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                        NULL,NULL,VecSamInputs_.getDVector(),NULL, 
                        NULL, NULL, NULL, NULL);
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
            VecSamOutputs_.getDVector(),
            VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      psuadeIO_->updateMethodSection(-1,nSamples_,-1,-1,-1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("skeep completed. Use 'write' to store.\n");
      printf("Use 'sinfo' to get updated sample information.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sfilter 
    //**/ remove duplicate points from another file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sfilter"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sfilter: delete sample points that exist ");
        printf("in another file \n");
        printf("Syntax: sfilter (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command removes sample points in the loaded ");
      printf("sample that already\n");
      printf("exist in another file.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Name of the other sample file for comparison : ");
      getString(pString, dataFile);
      kk = strlen(dataFile);
      dataFile[kk-1] = '\0';
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->setOutputLevel(0);
      status = ioPtr->readPsuadeFile(dataFile);
      if (status == 0)
      {
        ioPtr->getParameter("input_ninputs", pPtr);
        int nInps = pPtr.intData_;
        if (nInps != nInputs_)
        {
          printf("ERROR: nInputs are different.\n");
          printf("       incoming nInputs = %d\n", nInps);
          printf("       expected nInputs = %d\n", nInputs_);
          cmdStatus = 1;
          continue;
        }
        ioPtr->getParameter("output_noutputs", pPtr);
        int nOuts = pPtr.intData_;
        if (nOuts != nOutputs_)
        {
          printf("ERROR: nOutputs are different.\n");
          printf("       incoming nOutputs = %d\n", nOuts);
          printf("       expected nOutputs = %d\n", nOutputs_);
          cmdStatus = 1;
          continue;
        }
      }

      //**/ fetch the other sample
      ioPtr->getParameter("method_nsamples", pPtr);
      int nSam2 = pPtr.intData_;
      pData pInps, pOuts;
      ioPtr->getParameter("input_sample", pInps);
      double *dPtrX = pInps.dbleArray_;
      ioPtr->getParameter("output_sample", pOuts);
      double *dPtrY = pOuts.dbleArray_;
      psIVector vecTags;
      vecTags.setLength(nSamples_);
     
      //**/ compare
      for (ss = 0; ss < nSamples_; ss++)
      {
        for (kk = 0; kk < nSam2; kk++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            if (VecSamInputs_[ss*nInputs_+ii] != 
                dPtrX[kk*nInputs_+ii]) break;
          if (ii == nInputs_)
          {
            vecTags[ss] = 1;
            break;
          }
          for (ii = 0; ii < nOutputs_; ii++)
            if (VecSamOutputs_[ss*nOutputs_+ii] != 
                dPtrY[kk*nOutputs_+ii]) break;
          if (ii == nOutputs_)
          {
            vecTags[ss] = 1;
            break;
          }
        }
      }
      count = vecTags.sum();
      if (count == 0)
      {
        printf("INFO: no duplicates.\n");
        continue;
      }
      if (count == nSamples_)
      {
        printf("INFO: the resident and the given ");
        printf("samples are the same.\n"); 
        printf("      NO ACTION TAKEN.\n");
        continue;
      }
      psVector vecXT, vecYT;
      psIVector vecST;
      vecXT = VecSamInputs_;
      vecYT = VecSamOutputs_;
      vecST = VecSamStates_;
      VecSamInputs_.setLength((nSamples_-count)*nInputs_);
      VecSamOutputs_.setLength((nSamples_-count)*nOutputs_);
      VecSamStates_.setLength(nSamples_-count);
      kk = 0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (vecTags[ss] == 0)
        {
          for (ii = 0; ii < nInputs_; ii++)
            VecSamInputs_[kk*nInputs_+ii] = 
                    vecXT[ss*nInputs_+ii]; 
          for (ii = 0; ii < nOutputs_; ii++)
            VecSamOutputs_[kk*nOutputs_+ii] = 
                    vecYT[ss*nOutputs_+ii]; 
          VecSamStates_[kk] = vecST[kk]; 
          kk++;
        }
      }
      nSamples_ -= count;
      printf("New nSamples = %d\n", nSamples_);
      //**/ update the database and session
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                        NULL,NULL,VecSamInputs_.getDVector(),NULL, 
                        NULL, NULL, NULL, NULL);
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
            VecSamOutputs_.getDVector(),
            VecSamStates_.getIVector(),StrOutNames_.getStrings()); 
      psuadeIO_->updateMethodSection(-1,nSamples_,-1,-1,-1);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("sfilter completed. Use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ishuffle 
    //**/ change the order of the input parameters
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ishuffle"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ishuffle: re-arrange the orderings of input parameters.\n");
        printf("Syntax: ishuffle (no argument needed)\n");
        printf("NOTE: this command does not change the data file.\n");
        printf("      until a 'write' is issued.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command puts the sample inputs into a different order.\n");
      printf("You are to give the new orderings for all inputs.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ <= 1) 
      {
        printf("Number of inputs == %d ==> no reshuffle\n", nInputs_);
        continue;
      }
      vecIT.setLength(nInputs_);
      kk  = 0;
      while (kk < nInputs_)
      {
        sprintf(pString,
                "Enter the %d-th input (1 - %d) : ",kk+1,nInputs_);
        vecIT[kk] = getInt(1, nInputs_, pString);
        vecIT[kk] = vecIT[kk] - 1;
        kk++;
      }
      //**/ shuffle the sample matrix and input names
      vecXT.setLength(nSamples_*nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        kk = vecIT[ii];
        for (jj = 0; jj < nSamples_; jj++)
          vecXT[jj*nInputs_+ii] = VecSamInputs_[jj*nInputs_+kk];
      }
      VecSamInputs_.load(nSamples_*nInputs_, vecXT.getDVector());

      StrInpNames_.shuffleStrings(vecIT);

      //**/ shuffle the bounds
      psVector VecTLBs, VecTUBs;
      VecTLBs.setLength(nInputs_);
      VecTUBs.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        kk = vecIT[ii];
        VecTLBs[ii] = VecILowerBs_[kk];
        VecTUBs[ii] = VecIUpperBs_[kk];
      }
      VecILowerBs_ = VecTLBs;
      VecIUpperBs_ = VecTUBs;
      //**/ shuffle the input PDF info
      psIVector VecTPDFs;
      psVector  VecTMeans, VecTStds;
      VecTPDFs.setLength(nInputs_);
      VecTMeans.setLength(nInputs_);
      VecTStds.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        kk = vecIT[ii];
        VecTPDFs[ii]  = VecInpPDFs_[kk];
        VecTMeans[ii] = VecInpMeans_[kk];
        VecTStds[ii]  = VecInpStds_[kk];
      }
      VecInpPDFs_ = VecTPDFs;
      VecInpMeans_ = VecTMeans;
      VecInpStds_ = VecTStds;

      if (inputCMat_ != NULL)
      {
        psMatrix *tmpMat = new psMatrix();
        tmpMat->setDim(nInputs_, nInputs_);
        for (ii = 0; ii < nInputs_; ii++)
        {
          kk = vecIT[ii];
          for (jj = 0; jj < nInputs_; jj++)
          {
            ll = vecIT[jj];
            ddata = inputCMat_->getEntry(kk,ll);
            tmpMat->setEntry(ii,jj,ddata);
          }
        }
        delete inputCMat_;
        inputCMat_ = tmpMat;
      }

      //**/ update the data base and session
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
                   VecInpPDFs_.getIVector(), VecInpMeans_.getDVector(),
                   VecInpStds_.getDVector(),inputCMat_); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ishuffle completed. Use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oshuffle 
    //**/ change the order of the output parameters
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oshuffle"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ishuffle: re-arrange the orderings of sample outputs.\n");
        printf("Syntax: oshuffle (no argument needed)\n");
        printf("NOTE: this command does not change the sample file.\n");
        printf("      until a 'write' is issued.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command puts the sample outputs into a different order.\n");
      printf("You are to give the new orderings for all outputs.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ <= 1) 
      {
        printf("Number of outputs == %d ==> no reshuffle\n", nOutputs_);
        continue;
      }
      vecIT.setLength(nOutputs_);
      kk  = 0;
      while (kk < nOutputs_)
      {
        sprintf(pString,
                "Enter the %d-th output (1 - %d) : ",kk+1,nOutputs_);
        vecIT[kk] = getInt(1, nOutputs_, pString);
        vecIT[kk] = vecIT[kk] - 1;
        kk++;
      }
      //**/ shuffle the sample outputs and output names
      vecYT.setLength(nSamples_*nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        kk = vecIT[ii];
        for (jj = 0; jj < nSamples_; jj++)
          vecYT[jj*nOutputs_+ii] = VecSamOutputs_[jj*nOutputs_+kk];
      }
      VecSamOutputs_ = vecYT;

      StrOutNames_.shuffleStrings(vecIT);

      //**/ update the data base and session
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(),VecSamStates_.getIVector(),
                   StrOutNames_.getStrings());
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("oshuffle completed. Use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ itag 
    //**/ tag sample points based on input values 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "itag"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("itag: tag sample points based on input values\n");
        printf("Syntax: itag (no argument needed)\n");
        printf("This command is useful for extracting indices ");
        printf("of sample points\n");
        printf("with input values falling in a given range. This ");
        printf("command can be\n");
        printf("called multiple times to impose filters (the AND ");
        printf("operation) on the\n");
        printf("sample. This tagged list can be reset to the full list ");
        printf("by a 'load'.\n");
        printf("NOTE: This command does not change the loaded sample. ");
        printf("It only\n");
        printf("      displays sample indices information after tagging.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a list of sample points with ");
      printf("the selected input\n");
      printf("that falls in a given range. At the end of this command, ");
      printf("the list will\n");
      printf("be displayed. You can do 'itag/otag' multiple times to ");
      printf("trim the list.\n");
      printf("NOTE: A 'load' or 'read' will reset the full list.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (VecTagArray_.length() == 0 && nSamples_ > 0) 
      {
        VecTagArray_.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++) VecTagArray_[ii] = 1;
      } 
      sprintf(pString, "Enter input number (1 - %d) : ", nInputs_);
      inputID = getInt(1, nInputs_, pString);
      inputID--;
      Xmin = VecSamInputs_[inputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamInputs_[sInd*nInputs_+inputID] < Xmin)
          Xmin = VecSamInputs_[sInd*nInputs_+inputID];
      Xmax = VecSamInputs_[inputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamInputs_[sInd*nInputs_+inputID] > Xmax)
          Xmax = VecSamInputs_[sInd*nInputs_+inputID];
      printf("Xmin and Xmax found = %e %e.\n", Xmin, Xmax);
      sprintf(pString,"Enter the lower threshold (Xmin = %e) : ",Xmin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper threshold (Xmax = %e) : ",Xmax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: Lower bound should be < upper bound.\n");
        cmdStatus = 1;
        continue;
      }
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamInputs_[sInd*nInputs_+inputID] < threshL ||
            VecSamInputs_[sInd*nInputs_+inputID] > threshU)
        {
          VecTagArray_[sInd] = 0;
        }
      }
      printf("Sample points left after thresholding: \n[");
      for (sInd = 0; sInd < nSamples_; sInd++)
        if (VecTagArray_[sInd] == 1) printf("%d ", sInd+1);
      printf("]\n");
      printf("itag completed.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ otag 
    //**/ tag sample points based on output values 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "otag"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("otag: tag sample points based on output values\n");
        printf("Syntax: otag (no argument needed)\n");
        printf("This command is useful for extracting indices ");
        printf("of sample points\n");
        printf("with output values falling in a given range. This ");
        printf("command can be\n");
        printf("called multiple times to impose filters (the AND ");
        printf("operation) on the\n");
        printf("sample. This tagged list can be reset to the full list ");
        printf("by a 'load'.\n");
        printf("NOTE: This command does not change the loaded sample. ");
        printf("It only\n");
        printf("      displays sample indices information after tagging.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a list of sample points with ");
      printf("the selected output\n");
      printf("that falls in a given range. At the end of this command, ");
      printf("the list will\n");
      printf("be displayed. You can do 'itag/otag' multiple times to ");
      printf("trim the list.\n");
      printf("NOTE: A 'load' or 'read' will reset the list.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (VecTagArray_.length() == 0 && nSamples_ > 0) 
      {
        VecTagArray_.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++) VecTagArray_[ii] = 1;
      } 
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      Ymin = VecSamOutputs_[outputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
          Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
      Ymax = VecSamOutputs_[outputID];
      for (sInd = 1; sInd < nSamples_; sInd++)
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
          Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
      printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
      sprintf(pString,"Enter the lower threshold (Ymin = %e) : ",Ymin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper threshold (Ymax = %e) : ",Ymax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: Lower bound should be < upper bound.\n");
        cmdStatus = 1;
        continue;
      }
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < threshL ||
            VecSamOutputs_[sInd*nOutputs_+outputID] > threshU)
        {
          VecTagArray_[sInd] = 0;
        }
      }
      printf("Sample points left after thresholding: \n [");
      for (sInd = 0; sInd < nSamples_; sInd++)
        if (VecTagArray_[sInd] == 1) printf("%d ", sInd+1);
      printf("]\n");
      printf("otag completed.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oreset 
    //**/ reset an output to certain value
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oreset"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oreset: reset certain output values to another value.\n");
        printf("Syntax: oreset (no argument needed)\n");
        printf("This command is useful when some sample outputs are to\n");
        printf("be set to a different value (e.g. set undefined to 0).\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command resets a selected output that falls ");
      printf("in a given range.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      oInd = getInt(1, nOutputs_, pString);
      oInd--;
      Ymin = PSUADE_UNDEFINED;
      Ymax = -PSUADE_UNDEFINED;
      for (ii = 0; ii < nSamples_; ii++)
      {
        if (VecSamOutputs_[ii*nOutputs_+oInd] < Ymin)
          Ymin = VecSamOutputs_[ii*nOutputs_+oInd];
        if (VecSamOutputs_[ii*nOutputs_+oInd] > Ymax)
          Ymax = VecSamOutputs_[ii*nOutputs_+oInd];
      }
      printf("Lower and upper values for this output : %e %e\n",Ymin, Ymax);
      printf("Now specify the range of current values to be reset: \n");
      sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
      Ymin = getDouble(pString);
      sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
      Ymax = getDouble(pString);
      sprintf(pString,"Enter the desired output value to be set to : ");
      ddata = getDouble(pString);
      for (ii = 0; ii < nSamples_; ii++)
      {
        if (VecSamOutputs_[ii*nOutputs_+oInd] >= Ymin &&
            VecSamOutputs_[ii*nOutputs_+oInd] <= Ymax)
          VecSamOutputs_[ii*nOutputs_+oInd] = ddata;
      }
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_, 
                   VecSamOutputs_.getDVector(),
                   VecSamStates_.getIVector(), NULL); 
      printf("oreset completed: use 'write' to store the modified sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iscale
    //**/ re-scale an input in the loaded sample 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iscale"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iscale: scale a sample input to a different range.\n");
        printf("Syntax: iscale (no argument needed)\n");
        printf("NOTE: Input values will be re-scaled to the new range.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command scales a single sample input ");
      printf("to a different range.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      printf("Current lower bound for input %d = %16.8e\n",
             iInd+1,VecILowerBs_[iInd]);
      sprintf(pString,"Enter new lower bound for input %d : ",iInd+1);
      Xmin = getDouble(pString);
      printf("Current upper bound for input %d = %16.8e\n",
             iInd+1,VecIUpperBs_[iInd]);
      sprintf(pString,"Enter new upper bound for input %d : ",iInd+1);
      Xmax = getDouble(pString);
      if (Xmin >= Xmax)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        cmdStatus = 1;
        continue;
      }
      for (ii = 0; ii < nSamples_; ii++)
      {
        dtemp = VecSamInputs_[ii*nInputs_+iInd];
        dtemp = (dtemp-VecILowerBs_[iInd]) / 
                (VecIUpperBs_[iInd]-VecILowerBs_[iInd]);
        VecSamInputs_[ii*nInputs_+iInd] = dtemp * (Xmax - Xmin) + Xmin;
      }
      if (VecILowerBs_.length() > 0)
      { 
        VecILowerBs_[iInd] = Xmin;
        VecIUpperBs_[iInd] = Xmax;
        VecInpPDFs_[iInd] = 0;
        VecInpMeans_[iInd] = 0;
        VecInpStds_[iInd] = 0;
        if (inputCMat_ != NULL)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            inputCMat_->setEntry(jj,iInd,0.0);
            inputCMat_->setEntry(iInd,jj,0.0);
          }
          inputCMat_->setEntry(iInd,iInd,1.0);
        }
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
                   VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
                   VecSamInputs_.getDVector(),NULL,NULL,NULL,NULL,NULL); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("iscale completed: use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oscale
    //**/ re-scale an output in the loaded sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oscale"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oscale: scale a sample output to a different range.\n");
        printf("Syntax: oscale (no argument needed)\n");
        printf("NOTE: Output values will be re-scaled to the new ");
        printf("range.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command scales a single sample output ");
      printf("to a different range.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter output number (1 - %d) : ", nOutputs_);
      oInd = getInt(1, nOutputs_, pString);
      oInd--;
      double Ymax=-PSUADE_UNDEFINED, Ymin=PSUADE_UNDEFINED;
      for (ss = 0; ss < nSamples_; ss++)
      {
        if (VecSamOutputs_[ss*nOutputs_+oInd] > Ymax)
          Ymax = VecSamOutputs_[ss*nOutputs_+oInd];
        if (VecSamOutputs_[ss*nOutputs_+oInd] < Ymin)
          Ymin = VecSamOutputs_[ss*nOutputs_+oInd];
      }
      printf("Current lower bound for output %d = %16.8e\n",oInd+1,Ymin);
      sprintf(pString,"Enter new lower bound for output %d : ",oInd+1);
      double newYmin = getDouble(pString);
      printf("Current upper bound for output %d = %16.8e\n",oInd+1,Ymax);
      sprintf(pString,"Enter new upper bound for output %d : ",oInd+1);
      double newYmax = getDouble(pString);
      if (newYmin >= newYmax)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        cmdStatus = 1;
        continue;
      }
      for (ss = 0; ss < nSamples_; ss++)
      {
        dtemp = VecSamOutputs_[ss*nOutputs_+oInd];
        dtemp = (dtemp - Ymin) / (Ymax - Ymin);
        VecSamOutputs_[ss*nOutputs_+oInd] = 
             dtemp * (newYmax - newYmin) + newYmin;
      }
      psuadeIO_->updateOutputSection(nSamples_,nOutputs_,
                   VecSamOutputs_.getDVector(), 
                   VecSamStates_.getIVector(), NULL); 
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("oscale completed: use 'write' to store.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ireset 
    //**/ reset an input to certain values
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ireset"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ireset: map an input parameter to some distinct value\n");
        printf("Syntax: ireset (no argument needed)\n");
        printf("This command is useful when a sample input is to be\n");
        printf("re-mapped from intervals to distinct values.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command resets a selected input that falls ");
      printf("in a given range\n");
      printf("to a fixed value. It is useful to convert continuous ");
      printf("inputs to discrete.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      Xmin = PSUADE_UNDEFINED;
      Xmax = -PSUADE_UNDEFINED;
      for (ii = 0; ii < nSamples_; ii++)
      {
        if (VecSamInputs_[ii*nInputs_+iInd] < Xmin)
          Xmin = VecSamInputs_[ii*nInputs_+iInd];
        if (VecSamInputs_[ii*nInputs_+iInd] > Xmax)
          Xmax = VecSamInputs_[ii*nInputs_+iInd];
      }
      printf("Lower and upper values for this input : %e %e\n",Xmin, Xmax);
      printf("Now specify the range of current values to be reset.\n");
      sprintf(pString,"Enter the lower bound (inclusive) of this range : ");
      Xmin = getDouble(pString);
      sprintf(pString,"Enter the upper bound (inclusive) of this range : ");
      Xmax = getDouble(pString);
      sprintf(pString,"Enter the desired input value to be set to : ");
      ddata = getDouble(pString);
      for (jj = 0; jj < nSamples_; jj++)
      {
        if (VecSamInputs_[jj*nInputs_+iInd] >= Xmin &&
            VecSamInputs_[jj*nInputs_+iInd] <= Xmax)
          VecSamInputs_[jj*nInputs_+iInd] = ddata;
      }
      if (VecILowerBs_.length() > 0)
      { 
        VecILowerBs_[iInd] = Xmin;
        VecIUpperBs_[iInd] = Xmax;
        VecInpPDFs_[iInd] = 0;
        VecInpMeans_[iInd] = 0;
        VecInpStds_[iInd] = 0;
        if (inputCMat_ != NULL)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            inputCMat_->setEntry(jj,iInd,0.0);
            inputCMat_->setEntry(iInd,jj,0.0);
          }
          inputCMat_->setEntry(iInd,iInd,1.0);
        }
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ireset completed: use 'write' to store the modified sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ifloor 
    //**/ truncate an input to integer
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ifloor"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ifloor: map an input to nearest smaller integer\n");
        printf("Syntax: ifloor (no argument needed)\n");
        printf("This command is useful when some of the inputs are\n");
        printf("discrete. PSUADE creates samples based on continuous\n");
        printf("variables. This command helps to modify the samples\n");
        printf("to accommodate discrete variables.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command converts some of the sample inputs to ");
      printf("discrete values\n");
      printf("through truncation.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      for (jj = 0; jj < nSamples_; jj++)
      {
        kk = (int) VecSamInputs_[jj*nInputs_+iInd];
        VecSamInputs_[jj*nInputs_+iInd] = (double) kk;
        if (VecILowerBs_.length() > 0 && (kk < VecILowerBs_[iInd]))
          VecILowerBs_[iInd] = (double) kk;
      }
      if (VecILowerBs_.length() > 0)
      { 
        VecInpPDFs_[iInd] = 0;
        VecInpMeans_[iInd] = 0;
        VecInpStds_[iInd] = 0;
        if (inputCMat_ != NULL)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            inputCMat_->setEntry(jj,iInd,0.0);
            inputCMat_->setEntry(iInd,jj,0.0);
          }
          inputCMat_->setEntry(iInd,iInd,1.0);
        }
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("ifloor completed: use 'write' to store the modified sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iceil 
    //**/ round an input to integer
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iceil"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iceil: map an input to nearest larger integer\n");
        printf("Syntax: iceil (no argument needed)\n");
        printf("This command is useful when some of the inputs are\n");
        printf("discrete. PSUADE creates samples based on continuous\n");
        printf("variables. This command helps to modify the samples\n");
        printf("to accommodate discrete variables.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command converts some of the sample inputs to ");
      printf("discrete values\n");
      printf("by taking their ceiling integer values.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      for (jj = 0; jj < nSamples_; jj++)
      {
        ddata = VecSamInputs_[jj*nInputs_+iInd];
        kk = (int) ddata;
        if (kk != ddata) kk++;
        VecSamInputs_[jj*nInputs_+iInd] = (double) kk;
        if (VecIUpperBs_.length() > 0 && (kk > VecIUpperBs_[iInd]))
          VecIUpperBs_[iInd] = (double) kk;
      }
      if (VecIUpperBs_.length() > 0)
      { 
        VecInpPDFs_[iInd] = 0;
        VecInpMeans_[iInd] = 0;
        VecInpStds_[iInd] = 0;
        if (inputCMat_ != NULL)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            inputCMat_->setEntry(jj,iInd,0.0);
            inputCMat_->setEntry(iInd,jj,0.0);
          }
          inputCMat_->setEntry(iInd,iInd,1.0);
        }
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("iceil completed: use 'write' to store the modified sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iround 
    //**/ round an input to integer
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iround"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iround: map an input to nearest integer\n");
        printf("Syntax: iround (no argument needed)\n");
        printf("This command is useful when some of the inputs are\n");
        printf("discrete. PSUADE creates samples based on continuous\n");
        printf("variables. This command helps to modify the samples\n");
        printf("to accommodate discrete variables.\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command converts some of the sample inputs to ");
      printf("discrete values\n");
      printf("by rounding to the nearest integers.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString,"Enter input number (1 - %d) : ", nInputs_);
      iInd = getInt(1, nInputs_, pString);
      iInd--;
      for (jj = 0; jj < nSamples_; jj++)
      {
        ddata = VecSamInputs_[jj*nInputs_+iInd];
        kk = (int) ddata;
        if (kk != ddata) kk = (int) (ddata + 0.5);
        VecSamInputs_[jj*nInputs_+iInd] = (double) kk;
        if (VecILowerBs_.length() > 0 && (kk < VecILowerBs_[iInd]))
          VecILowerBs_[iInd] = (double) kk;
        if (VecIUpperBs_.length() > 0 && (kk > VecIUpperBs_[iInd]))
          VecIUpperBs_[iInd] = (double) kk;
      }
      if (VecIUpperBs_.length() > 0)
      { 
        VecInpPDFs_[iInd] = 0;
        VecInpMeans_[iInd] = 0;
        VecInpStds_[iInd] = 0;
        if (inputCMat_ != NULL)
        {
          for (jj = 0; jj < nInputs_; jj++)
          {
            inputCMat_->setEntry(jj,iInd,0.0);
            inputCMat_->setEntry(iInd,jj,0.0);
          }
          inputCMat_->setEntry(iInd,iInd,1.0);
        }
      }
      psuadeIO_->updateInputSection(nSamples_,nInputs_,NULL,
              VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
              VecSamInputs_.getDVector(),NULL,VecInpPDFs_.getIVector(), 
              VecInpMeans_.getDVector(), VecInpStds_.getDVector(), 
              inputCMat_);
      if (currSession != NULL) delete currSession;
      currSession = new PsuadeSession();
      psuadeIO_->getSession(currSession);
      printf("iround completed: use 'write' to store the modified sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iotrace 
    //**/ trace input and output values for each run
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iotrace"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iotrace: trace input/output values for each run\n");
        printf("Syntax: iotrace (no argument needed)\n");
        printf("This command is useful for tracing which values of\n");
        printf("inputs correspond to which values of the outputs.\n");
        printf("NOTE: one line for each sample point connecting all\n");
        printf("      inputs and outputs.\n");
        continue;
      }
      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command creates a plot that has one line ");
      printf("connecting all input and\n");
      printf("output values for all sample points.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs_);
      kk = getInt(1, nInputs_, pString);
      psIVector vecII;
      if (kk < nInputs_)
      {
        vecII.setLength(kk);
        for (ii = 0; ii < kk; ii++)
        {
          sprintf(pString,"Enter the %d-th input (1 - %d): ",ii+1,nInputs_);
          vecII[ii] = getInt(1,nInputs_,pString);
          vecII[ii]--;
        }
      }
      else
      {
        kk = nInputs_;
        vecII.setLength(kk);
        for (ii = 0; ii < nInputs_; ii++) vecII[ii] = ii;
      }
      sprintf(pString,"Number of outputs to plot (1 - %d) : ",nOutputs_);
      ll = getInt(1, nOutputs_, pString);
      vecIT.clean();
      if (ll < nOutputs_)
      {
        vecIT.setLength(ll);
        for (ii = 0; ii < ll; ii++)
        {
          sprintf(pString,"Enter the %d-th output (1 - %d): ",
                  ii+1, nOutputs_);
          vecIT[ii] = getInt(1,nOutputs_,pString);
          vecIT[ii]--;
        }
      }
      else
      {
        ll = nOutputs_;
        vecIT.setLength(nOutputs_);
        for (ii = 0; ii < nOutputs_; ii++) vecIT[ii] = ii;
      }
      if (plotScilab())
      {
        fp = fopen("scilabiotrace.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabiotrace.sci\n");
          cmdStatus = 1;
          continue;
        }
      }
      else
      {
        fp = fopen("matlabiotrace.m","w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabiotrace.m.\n");
          cmdStatus = 1;
          continue;
        }
      }
      vecVT.setLength(nOutputs_);
      vecWT.setLength(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        vecVT[ii] =  PSUADE_UNDEFINED;
        vecWT[ii] = -PSUADE_UNDEFINED;
        for (ss = 0; ss < nSamples_; ss++)
        {
          ddata = VecSamOutputs_[ss*nOutputs_+ii];
          vecVT[ii] = (ddata < vecVT[ii]) ? ddata : vecVT[ii];
          vecWT[ii] = (ddata > vecWT[ii]) ? ddata : vecWT[ii];
        }
      }
      fprintf(fp, "Y = 1:%d;\n", kk+ll); 
      fprintf(fp, "X = [\n"); 
      for (ss = 0; ss < nSamples_; ss++)
      {
        for (ii = 0; ii < kk; ii++)
        {
          ddata = (VecSamInputs_[ss*nInputs_+vecII[ii]] - 
                   VecILowerBs_[vecII[ii]]) /
                  (VecIUpperBs_[vecII[ii]]-VecILowerBs_[vecII[ii]]);
          fprintf(fp, "%e ", ddata);
        }
        for (ii = 0; ii < ll; ii++)
        {
          ddata = (VecSamOutputs_[ss*nOutputs_+vecIT[ii]] - 
                   vecVT[vecIT[ii]]) /
                  (vecWT[vecIT[ii]] - vecVT[vecIT[ii]]);
          fprintf(fp, "%e ", ddata);
        }
        fprintf(fp, "\n"); 
      }
      fprintf(fp, "];\n"); 
      fprintf(fp,"plot(X',Y)\n");
      sprintf(pString,"Input/Output Values");
      fwritePlotXLabel(fp, pString);
      sprintf(pString,"Inputs 1:%d, Outputs %d:%d",kk,kk+1,kk+ll);
      fwritePlotYLabel(fp, pString);
      sprintf(pString,"Trace Inputs/Outputs per sample point");
      fwritePlotTitle(fp, pString);
      fwritePlotAxes(fp);
      if (plotScilab())
           printf("Scilab iotrace is now in scilabiotrace.sci\n");
      else printf("Matlab iotrace is now in matlabiotrace.m\n");
      fclose(fp);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iplot2all 
    //**/ plot all pairs of inputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot2_all"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot2_all: plot all pairs of inputs\n");
        printf("Syntax: iplot2_all (no argument needed)\n");
        printf("This command is useful to find out which parameters\n");
        printf("are responsible for the failed runs.\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      printAsterisks(PL_INFO, 0);
      printf("This command creates a plot with scatter subplots ");
      printf("for all pairs of selected\n");
      printf("inputs.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Number of inputs to use (2 - %d) : ",
              nInputs_);
      kk = getInt(2, nInputs_, pString);
      if (kk < nInputs_)
      {
        vecIT.setLength(kk);
        for (ii = 0; ii < kk; ii++)
        {
          sprintf(pString,"Enter the %d-th input (1 - %d): ",
                  ii+1, nInputs_);
          vecIT[ii] = getInt(1,nInputs_,pString);
          vecIT[ii] = vecIT[ii] - 1;
        }
      }
      else
      {
        kk = nInputs_;
        vecIT.setLength(kk);
        for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = ii;
      }

      if (plotScilab())
      {
        fp = fopen("scilabiplt2all.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabiplt2all.sci\n");
          cmdStatus = 1;
          continue;
        }
      }
      else
      {
        fp = fopen("matlabiplt2all.m","w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabiplt2all.m.\n");
          cmdStatus = 1;
          continue;
        }
      }
      //**/ obsolete (Aug 2009)
      //**/for (ii = 0; ii < nInputs_; ii++)
      //**/{
      //**/   fprintf(fp, "X%d = [\n", ii+1); 
      //**/   for (jj = 0; jj < nSamples_; jj++)
      //**/      fprintf(fp, "%e\n", VecSamInputs_[jj*nInputs_+ii]);
      //**/   fprintf(fp, "];"); 
      //**/}
      sprintf(pString," plotMode  = 0 : plot all in a single plot");
      fwriteComment(fp, pString);
      sprintf(pString," plotMode ~= 0 : plot one at a time");
      fwriteComment(fp, pString);
      sprintf(pString," ranflag: to distinguish between identical points");
      fwriteComment(fp, pString);
      sprintf(pString,"         by adding a small perturbation (when on)");
      fwriteComment(fp, pString);
      fprintf(fp,"plotMode = 0;\n");
      fprintf(fp,"ranflag  = 0;\n");
      for (ii = 0; ii < kk; ii++)
      {
        fprintf(fp, "X%d = [\n", ii+1); 
        for (jj = 0; jj < nSamples_; jj++)
          fprintf(fp, "%e\n", VecSamInputs_[jj*nInputs_+vecIT[ii]]);
        fprintf(fp, "];\n"); 
      }
      fprintf(fp, "S = [\n"); 
      for (jj = 0; jj < nSamples_; jj++)
         fprintf(fp, "%d\n", VecSamStates_[jj]);
      fprintf(fp, "];\n"); 
      fprintf(fp, "X = [\n"); 
      for (ii = 0; ii < kk; ii++) fprintf(fp, "X%d ", ii+1); 
      fprintf(fp, "S];\n"); 
      fprintf(fp,"ms=9;\n");
      fprintf(fp,"fs=9;\n");
      for (ii = 0; ii < kk; ii++)
      {
        for (jj = ii; jj < kk; jj++)
        {
          fprintf(fp,"if plotMode == 0\n");
          fprintf(fp,"subplot(%d,%d,%d)\n",kk,kk,ii*kk+jj+1);
          fprintf(fp,"else\n");
          if (ii + jj > 0)
          {
            fprintf(fp,"pause\n");
            fprintf(fp,"disp('Press enter to continue')\n");
          }
          fwritePlotCLF(fp);
          fprintf(fp,"end;\n");
          if (plotScilab())
          {
            fprintf(fp,"iset = find(S==0);\n");
            fprintf(fp,
               "plot(X(iset,%d),X(iset,%d),'rX','MarkerSize',ms)\n",
               jj+1,ii+1);
            fprintf(fp,"iset = find(S~=0);\n");
            fprintf(fp,
               "plot(X(iset,%d),X(iset,%d),'bX','MarkerSize',ms)\n",
               jj+1,ii+1);
          }
          else
          {
            fprintf(fp, "iset = find(S==0);\n");
            fprintf(fp, "plot(X(iset,%d).*(1+ranflag*", jj+1);
            fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",ii+1);
            fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'rX',");
            fprintf(fp, "'MarkerSize',ms)\n");
            if (plotScilab())
                 fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            else fprintf(fp, "hold on\n");
            fprintf(fp, "iset = find(S~=0);\n");
            fprintf(fp, "plot(X(iset,%d).*(1+ranflag*",jj+1);
            fprintf(fp, "rand(size(iset,1),1)/100),X(iset,%d)",
                    ii+1);
            fprintf(fp, ".*(1+ranflag*rand(size(iset,1),1)/100),'b*',");
            fprintf(fp, "'MarkerSize',ms)\n");
          }
          fwritePlotXLabel(fp, StrInpNames_[vecIT[jj]]);
          fwritePlotYLabel(fp, StrInpNames_[vecIT[ii]]);
          fwritePlotAxes(fp);
          if (plotScilab())
          {
            fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                    VecILowerBs_[vecIT[jj]], VecILowerBs_[vecIT[ii]],
                    VecIUpperBs_[vecIT[jj]], VecIUpperBs_[vecIT[ii]]);
          }
          else
          {
            fprintf(fp, "hold off\n");
            fprintf(fp,"axis([%e %e %e %e])\n",
                    VecILowerBs_[vecIT[jj]], VecIUpperBs_[vecIT[jj]],
                    VecILowerBs_[vecIT[ii]], VecIUpperBs_[vecIT[ii]]);
          }
          fprintf(fp,"disp('red X: failed runs.')\n");
        }
      }
      if (plotScilab())
           printf("The Scilab file is in scilabiplt2all.sci\n");
      else printf("The Matlab file is in matlabiplt2all.m\n");
      fclose(fp);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ihist 
    //**/ plot histogram for a selected input
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ihist"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ihist: plot histogram for an input\n");
        printf("Syntax: ihist (no argument needed)\n");
        printf("This command is useful to examine input distribution\n");
        printf("from a posterior sample, e.g. MCMCPostSample (loaded\n");
        printf("using 'iread').\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a histogram for the values ");
      printf("of a selected input.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Which input to generate histogram (1 - %d) : ",
              nInputs_);
      kk = getInt(1, nInputs_, pString);
      kk--;
      if (plotScilab()) fp = fopen("scilabihist.sci", "w");
      else              fp = fopen("matlabihist.m", "w");
      fprintf(fp, "X = [\n");
      for (ss = 0; ss < nSamples_; ss++)
        fprintf(fp, "%e \n",VecSamInputs_[ss*nInputs_+kk]);
      fprintf(fp, "];\n");
      fwritePlotCLF(fp);
      fwritePlotFigure(fp, 1);
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, X, style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,xk] = hist(X,10);\n");
        fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      sprintf(pString,"Sample Histogram for %s",StrInpNames_[kk]);
      fwritePlotTitle(fp, pString);
      fwritePlotXLabel(fp, "Input Value");
      fwritePlotYLabel(fp, "Probabilities");
      fclose(fp);
      if (plotScilab()) 
           printf("Histogram is available in scilabihist.sci\n");
      else printf("Histogram is available in matlabihist.m\n");
      cmdStatus = 0;
    }
       
    //**/ -------------------------------------------------------------
    // +++ ihist2 
    //**/ plot 2D-histogram for two selected inputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ihist2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ihist2: plot histogram for two inputs\n");
        printf("Syntax: ihist2 (no argument needed)\n");
        printf("This command is useful to examine input distributions\n");
        printf("from a posterior sample, e.g. MCMCPostSample (loaded\n");
        printf("using 'iread').\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 2)
      {
        printf("ERROR: this command is for nInputs_ >= 2.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates a histogram for the values ");
      printf("of 2 selected inputs\n");
      printf("using heat map instead of 3D bars.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Select the first input (1 - %d) : ", nInputs_);
      jj = getInt(1, nInputs_, pString);
      jj--;
      sprintf(pString,"Select the second input (1 - %d) : ", nInputs_);
      kk = getInt(1, nInputs_, pString);
      kk--;
      if (plotScilab())
      {
        printf("ERROR: scilab plot not currently available.\n");
        cmdStatus = 1;
        continue;
      }
      fp = fopen("matlabihist2.m", "w");
      fprintf(fp, "X = [\n");
      for (ss = 0; ss < nSamples_; ss++)
        fprintf(fp, "%e %e\n",VecSamInputs_[ss*nInputs_+jj],
                VecSamInputs_[ss*nInputs_+kk]);
      fprintf(fp, "];\n");
      fprintf(fp, "X1 = X(:,1);\n");
      fprintf(fp, "X2 = X(:,2);\n");
      fprintf(fp, "nb = 20;\n");
      fprintf(fp, "bins = zeros(nb,nb);\n");
      fprintf(fp, "x1min = %e;\n",VecILowerBs_[jj]);
      fprintf(fp, "x1max = %e;\n",VecIUpperBs_[jj]);
      fprintf(fp, "if (x1min == x1max)\n");
      fprintf(fp, "   if (x1min == 0)\n");
      fprintf(fp, "      x1min = -0.1;\n");
      fprintf(fp, "      x1max = 0.1;\n");
      fprintf(fp, "   else\n");
      fprintf(fp, "      x1wid = x1max;\n");
      fprintf(fp, "      x1min = x1min - x1wid * 0.1;\n");
      fprintf(fp, "      x1max = x1max + x1wid * 0.1;\n");
      fprintf(fp, "   end\n");
      fprintf(fp, "end\n");
      fprintf(fp, "x1wid = (x1max - x1min) / (nb - 1);\n");
      fprintf(fp, "x2min = %e;\n",VecILowerBs_[kk]);
      fprintf(fp, "x2max = %e;\n",VecIUpperBs_[kk]);
      fprintf(fp, "if (x2min == x2max)\n");
      fprintf(fp, "   if (x2min == 0)\n");
      fprintf(fp, "      x2min = -0.1;\n");
      fprintf(fp, "      x2max = 0.1;\n");
      fprintf(fp, "   else\n");
      fprintf(fp, "      x2wid = x2max - x2min;\n");
      fprintf(fp, "      x2min = x2min - x2wid * 0.1;\n");
      fprintf(fp, "      x2max = x2max + x2wid * 0.1;\n");
      fprintf(fp, "   end\n");
      fprintf(fp, "end\n");
      fprintf(fp, "x2wid = (x2max - x2min) / (nb - 1);\n");
      fprintf(fp, "for ii = 1 : nb\n");
      fprintf(fp, "   x1l = x1min + (ii-1) * x1wid;\n");
      fprintf(fp, "   x1u = x1min + ii * x1wid;\n");
      fprintf(fp, "   if ii == nb\n");
      fprintf(fp, "      iset = find(X1 >= x1l & X1 <= x1u);\n");
      fprintf(fp, "   else\n");
      fprintf(fp, "      iset = find(X1 >= x1l & X1 < x1u);\n");
      fprintf(fp, "   end\n");
      fprintf(fp, "   if (length(iset) > 0)\n");
      fprintf(fp, "      XX = X2(iset);\n");
      fprintf(fp, "      for jj = 1 : nb\n");
      fprintf(fp, "         x2l = x2min + (jj-1) * x2wid;\n");
      fprintf(fp, "         x2u = x2min + jj * x2wid;\n");
      fprintf(fp, "         if jj == nb\n");
      fprintf(fp, "            jset = find(XX >= x2l & XX <= x2u);\n");
      fprintf(fp, "         else\n");
      fprintf(fp, "            jset = find(XX >= x2l & XX < x2u);\n");
      fprintf(fp, "         end\n");
      fprintf(fp, "         bins(ii,jj) = length(jset);\n");
      fprintf(fp, "      end\n");
      fprintf(fp, "   end\n");
      fprintf(fp, "end\n");
      fwritePlotCLF(fp);
      fwritePlotFigure(fp, 1);
      fprintf(fp,"imagesc(bins)\n");
      fprintf(fp,"%%bar3(bins)\n");
      fprintf(fp,"xtick = x1min:(x1max-x1min)/2:x1max;\n");
      fprintf(fp,"set(gca,'XTick',0:nb/2:nb);\n");
      fprintf(fp,"set(gca,'XTickLabel', xtick);\n");
      fprintf(fp,"ytick = x2min:(x2max-x2min)/2:x2max;\n");
      fprintf(fp,"set(gca,'YTick',0:nb/2:nb);\n");
      fprintf(fp,"set(gca,'YTickLabel', ytick);\n");
      fprintf(fp,"set(gca,'YDir', 'normal');\n");
      fprintf(fp,"xlabel('%s','FontWeight','bold','FontSize',12)\n",
              StrInpNames_[jj]);
      fprintf(fp,"ylabel('%s','FontWeight','bold','FontSize',12)\n",
              StrInpNames_[kk]);
      fwritePlotAxesNoGrid(fp);
      sprintf(pString,"Sample Input Histogram");
      fwritePlotTitle(fp, pString);
      fclose(fp);
      printf("Histogram is available in matlabihist2.m\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ oplot_pdf 
    //**/ plot histogram for a selected output
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oplot_pdf"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oplot_pdf: plot histogram for an output\n");
        printf("Syntax: oplot_pdf (no argument needed)\n");
        printf("This command is useful to examine output distribution\n");
        printf("without using the 'ua' command.\n");
        continue;
      }
      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates histograms for the values ");
      printf("of selected outputs.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Which output to generate histogram (1 - %d) : ",
              nOutputs_);
      kk = getInt(1, nOutputs_, pString);
      kk--;
      if (plotScilab()) fp = fopen("scilabopltpdf.sci", "w");
      else              fp = fopen("matlabopltpdf.m", "w");
      fprintf(fp, "Y = [\n");
      for (ss = 0; ss < nSamples_; ss++)
        fprintf(fp, "%e \n",VecSamOutputs_[ss*nOutputs_+kk]);
      fprintf(fp, "];\n");
      fwritePlotCLF(fp);
      fwritePlotFigure(fp, 1);
      if (plotScilab())
      {
        fprintf(fp, "histplot(10, Y, style=2);\n");
        fprintf(fp, "a = gce();\n");
        fprintf(fp, "a.children.fill_mode = \"on\";\n");
        fprintf(fp, "a.children.thickness = 2;\n");
        fprintf(fp, "a.children.foreground = 0;\n");
        fprintf(fp, "a.children.background = 2;\n");
      }
      else
      {
        fprintf(fp, "[nk,yk] = hist(Y,10);\n");
        fprintf(fp, "bar(yk,nk/%d,1.0)\n",nSamples_);
      }
      fwritePlotAxes(fp);
      sprintf(pString,"Sample Output Histogram for %s",
              StrOutNames_[kk]);
      fwritePlotTitle(fp, pString);
      fwritePlotXLabel(fp, "Output Value");
      fwritePlotYLabel(fp, "Probabilities");
      fclose(fp);
      if (plotScilab()) 
           printf("Histogram is available in scilabopltpdf.sci\n");
      else printf("Histogram is available in matlabopltpdf.m\n");
      cmdStatus = 0;
    }
       
    //**/ -------------------------------------------------------------
    // +++ oplot2_pdf 
    //**/ plot 2D-histogram for two selected outputs
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "oplot2_pdf"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("oplot2_pdf: plot histogram for two outputs\n");
        printf("Syntax: oplot2_pdf (no argument needed)\n");
        continue;
      }
      if (psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ < 2)
      {
        printf("ERROR: this command is for nOutputs >= 2.\n");
        cmdStatus = 1;
        continue;
      }
      int nbins=20;
      if (plotScilab())
           strcpy(pString, "scilaboplt2pdf.sci");
      else strcpy(pString, "matlaboplt2pdf.m");
      vecWT.setLength(nOutputs_);
      vecVT.setLength(nOutputs_);
      for (ii = 0; ii < nOutputs_; ii++)
      {
        vecWT[ii] = PSUADE_UNDEFINED;
        vecVT[ii] = -PSUADE_UNDEFINED;
        for (jj = 0; jj < nSamples_; jj++)
        {
          if (VecSamOutputs_[jj*nOutputs_+ii] < vecWT[ii])
            vecWT[ii] = VecSamOutputs_[jj*nOutputs_+ii];
          if (VecSamOutputs_[jj*nOutputs_+ii] > vecVT[ii])
            vecVT[ii] = VecSamOutputs_[jj*nOutputs_+ii];
        }
      }
      genMatlabPlotFile(nOutputs_,vecWT.getDVector(),
          vecVT.getDVector(),nSamples_,
          VecSamOutputs_.getDVector(), StrOutNames_.getStrings(),
          pString,nbins);
      if (plotScilab()) strcpy(dataFile, "scilaboplt2pdf.sci");
      else              strcpy(dataFile, "matlaboplt2pdf.m");
      printf("Distribution PDFs are available in %s.\n",dataFile);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iplot_pdf 
    //**/ plot input pdf from sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot_pdf"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("iplot_pdf: plot pdf of selected inputs\n");
        printf("Syntax: iplot_pdf (no argument needed)\n");
        continue;
      }
      if (psuadeIO_ == NULL || nInputs_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command creates histograms for the values ");
      printf("of selected inputs.\n");
      printf("This command is similar to ihist except you can select ");
      printf("more than 1 input.\n");
      printDashes(PL_INFO, 0);
      sprintf(pString,"Number of inputs to plot (1 - %d) : ", nInputs_);
      kk = getInt(1, nInputs_, pString);
      vecIT.setLength(nInputs_);
      if (kk < nInputs_)
      {
        for (ii = 0; ii < kk; ii++)
        {
          sprintf(pString,"Enter the %d-th input (1 - %d): ",
                  ii+1, nInputs_);
          vecIT[ii] = getInt(1,nInputs_,pString);
          vecIT[ii]--;
        }
      }
      else
      {
        kk = nInputs_;
        for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = ii;
      }

      if (plotScilab())
      {
        fp = fopen("scilabipltpdf.sci", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file scilabipltpdf.sci\n");
          cmdStatus = 1;
          continue;
        }
      }
      else
      {
        fp = fopen("matlabipltpdf.m","w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabipltpdf.m.\n");
          cmdStatus = 1;
          continue;
        }
      }
      fprintf(fp, "XX = [\n"); 
      for (ii = 0; ii < nSamples_; ii++)
      {
        for (jj = 0; jj < kk; jj++)
          fprintf(fp, "%24.16e ",VecSamInputs_[ii*nInputs_+vecIT[jj]]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "];\n"); 
      for (ii = 0; ii < kk; ii++)
        fprintf(fp, "X%d = XX(:,%d);\n", ii+1, ii+1); 
      ll = (int) sqrt(1.0*kk);
      if (ll * ll != kk) ll++; 
      fwritePlotCLF(fp);
      for (ii = 0; ii < kk; ii++)
      {
        if (kk > 1) fprintf(fp,"subplot(%d,%d,%d)\n",ll,ll,ii+1);
        if (plotScilab())
        {
          fprintf(fp, "X = X%d;\n",ii+1);
          fprintf(fp, "ymin = min(X);\n");
          fprintf(fp, "ymax = max(X);\n");
          fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
          fprintf(fp, "if (ywid < 1.0e-12)\n");
          fprintf(fp, "   disp('range too small.')\n");
          fprintf(fp, "   halt\n");
          fprintf(fp, "end;\n");
          fprintf(fp, "histplot(10, X, style=2);\n");
          fprintf(fp, "a = gce();\n");
          fprintf(fp, "a.children.fill_mode = \"on\";\n");
          fprintf(fp, "a.children.thickness = 2;\n");
          fprintf(fp, "a.children.foreground = 0;\n");
          fprintf(fp, "a.children.background = 2;\n");
          fwritePlotAxes(fp);
          fwritePlotTitle(fp, "Input Distribution");
          if (StrInpNames_[vecIT[ii]] != NULL)
               fwritePlotXLabel(fp, StrInpNames_[vecIT[ii]]);
          else fwritePlotXLabel(fp, "Input Value");
          fwritePlotYLabel(fp, "Probabilities");
        }
        else
        {
          fprintf(fp, "X = X%d;\n",ii+1);
          fprintf(fp, "[nk,xk]=hist(X,10);\n");
          fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples_);
          fwritePlotAxes(fp);
          fwritePlotTitle(fp, "Input Distribution");
          if (StrInpNames_[vecIT[ii]] != NULL)
            fwritePlotXLabel(fp, StrInpNames_[vecIT[ii]]);
          else fwritePlotXLabel(fp, "Input Value");
            fprintf(fp, "axis tight\n");
          fwritePlotYLabel(fp, "Probabilities");
        }
      }
      fclose(fp);
      if (plotScilab())
           printf("Plot file is scilabipltpdf.sci\n");
      else printf("Plot file is matlabipltpdf.m\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ iplot2_pdf 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "iplot2_pdf"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command uses the loaded sample and generates\n");
        printf("single-input and pairwise histograms similar to.\n");
        printf("the posterior plots created in MCMC.\n");
        printf("So first load the posterior sample using 'iread'\n");
        printf("and then use this command.\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ < 2)
      {
        printf("ERROR: this command is for nInputs_ >= 2.\n");
        cmdStatus = 1;
        continue;
      }
      int nbins=20;
      if (plotScilab())
           strcpy(pString, "scilabiplt2pdf.sci");
      else strcpy(pString, "matlabiplt2pdf.m");
      genMatlabPlotFile(nInputs_,VecILowerBs_.getDVector(),
               VecIUpperBs_.getDVector(),nSamples_,
               VecSamInputs_.getDVector(),StrInpNames_.getStrings(),
               pString,nbins);
      if (plotScilab()) strcpy(dataFile, "scilabiplt2pdf.sci");
      else              strcpy(dataFile, "matlabiplt2pdf.m");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ printlevel 
    //**/ set diagnostics output level
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "printlevel"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("printlevel: set printlevel (0 - 5)\n");
        continue;
      }
      sscanf(lineIn, "%s %d", command, &ii);
      if (ii >= PL_MIN && ii <= PL_MAX)
      {
        printf("printlevel set from %d to %d\n", outputLevel_, ii);
        outputLevel_ = ii;
      }
      else
      {
        printf("printlevel out of range; set to %d\n", PL_INTERACTIVE);
        outputLevel_ = PL_INTERACTIVE;
      }
      setPrintLevelTS(outputLevel_);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ turn on printlevel=4 and interactive
    //**/ set diagnostics output level
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "on"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("on: set printlevel=4 and turn on all interactive modes\n");
        continue;
      }
      outputLevel_ = 4;
      psConfig_.RSExpertModeOn();
      psConfig_.AnaExpertModeOn();
      psConfig_.SamExpertModeOn();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gp_sa2 
    //**/ screening of parameters with layered GP
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gp_sa2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gp_sa2: GP SA on sub-divisions of each input\n");
        printf("        (not verified yet)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      faType = PSUADE_RS_GP1;
      printf("This command may be very time-consuming to execute.\n");
      count = nSamples_ / 100;
      if (count > 4) count = 4;
      sprintf(pString, "How many partitions (1 - %d) = ", count);
      int nParts = getInt(1, count, pString);
      count = nSamples_ / nParts + 1;
      vecXT.setLength(count*nInputs_);
      vecYT.setLength(count);
      vecWT.setLength(nInputs_);
      double *dPtrX = vecXT.getDVector();
      double *dPtrY = vecYT.getDVector();
      double *dPtrW = vecWT.getDVector();
      for (ii = 0; ii < nInputs_; ii++)
      {
        dPtrW[ii] = 0.0;
        ddata = (VecIUpperBs_[ii] - VecILowerBs_[ii]) / (double) nParts;
        for (jj = 0; jj < nParts; jj++)
        {
          count = 0;
          for (kk = 0; kk < nSamples_; kk++)
          {
            dtemp = VecSamInputs_[kk*nInputs_+ii];
            if (dtemp >= ddata*jj+VecILowerBs_[ii] &&
                dtemp <= ddata*(jj+1)+VecILowerBs_[ii])
            {
              for (ind = 0; ind < nInputs_; ind++)
                dPtrX[count*nInputs_+ind] = 
                            VecSamInputs_[kk*nInputs_+ii];
              dPtrY[count++] = VecSamOutputs_[kk*nOutputs_+outputID];
            }
          }
          if (count == 0)
          {
            printf("gp_sa2: ERROR encountered.\n");
            exit(1);
          }
          if (outputLevel_ > 1)
            printf("Input %3d Part %3d has size %d\n",ii+1,jj+1,count);
          FuncApprox *faPtr = genFA(faType, nInputs_, iOne, count);
          if (faPtr != NULL)
          {
            faPtr->setBounds(VecILowerBs_.getDVector(), 
                             VecIUpperBs_.getDVector());
            faPtr->setOutputLevel(outputLevel_);
            faPtr->initialize(dPtrX,dPtrY);
            strcpy(pString, "rank");
            targv[0] = (char *) pString;
            targv[1] = (char *) &ii;
            dtemp = faPtr->setParams(2, targv);
            dPtrW[ii] += dtemp;
            delete faPtr;
            faPtr = NULL;
          }
        }
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sot_sa2 
    //**/ sum-of-trees screening of parameters using boosting
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sot_sa2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sot_sa2: sum-of-trees screening using boosting\n");
        printf("Syntax: sot_sa2 (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      faType = PSUADE_RS_SOTS;
      FuncApprox *faPtr  = genFA(faType, nInputs_, iOne, nSamples_);
      if (faPtr != NULL)
      {
        faPtr->setBounds(VecILowerBs_.getDVector(), 
                         VecIUpperBs_.getDVector());
        faPtr->setOutputLevel(outputLevel_);
        vecYT.setLength(nSamples_);
        for (ii = 0; ii < nSamples_; ii++)
          vecYT[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
        status = faPtr->initialize(VecSamInputs_.getDVector(),
                                   vecYT.getDVector());
        strcpy(pString, "mode1");
        targv[0] = (char *) pString;
        faPtr->setParams(1, targv);
        strcpy(pString, "rank");
        targv[0] = (char *) pString;
        faPtr->setParams(1, targv);
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ delta_test 
    //**/ Delta test for variable selection
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "delta_test"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("delta_test: perform Delta test\n");
        printf("Syntax: delta_test (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("*     Delta Test for Variable Selection\n");
      printDashes(PL_INFO, 0);
      printf("* This test works on random or quasi-random sample ");
      printf("to identify inputs\n");
      printf("* that are important in driving output variability.\n");
      printf("* It works best for large samples (>>1000) with number ");
      printf("of inputs up to\n");
      printf("* maybe about 50.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin);
      if (lineIn2[0] != 'y') continue;

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ == 1)
      {
        printf("ERROR: no point of using this test for nInputs=1\n");
        cmdStatus = 1;
        continue;
      }
      if (nSamples_ < 1000)
      {
        printf("WARNING: sample size too small to give useful results\n");
      }
      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_DTEST;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ eta_test 
    //**/ Eta test for variable selection
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "eta_test"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("eta_test: perform Eta test (a variant of Delta test)\n");
        printf("This test is a derivative of delta_test.\n");
        printf("The power of this test has not been verified.\n");
        printf("Syntax: eta_test (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("*     Eta Test for Variable Selection\n");
      printDashes(PL_INFO, 0);
      printf("* This test is a variant of the delta test for ");
      printf("dimension reduction. It\n");
      printf("* works on random or quasi-random sample to identify ");
      printf("inputs that are\n");
      printf("* important in driving output variability.\n");
      printf("* It works best for large samples (>>1000) with number ");
      printf("of inputs up to\n");
      printf("* maybe about 50.\n");
      printf("* Set printlevel = 2 to show progress.\n");
      printf("* Set printlevel = 4 to show more details.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin);
      if (lineIn2[0] != 'y') continue;

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ == 1)
      {
        printf("ERROR: no point of using this test for nInputs=1\n");
        cmdStatus = 1;
        continue;
      }

      sprintf(pString, "Enter output number (1 - %d) = ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_ETEST;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gd_test 
    //**/ Gower distance analysis
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gd_test"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gd_test: perform Gower distance analysis\n");
        printf("The power of this test has not been verified.\n");
        printf("This test is not ready yet.\n");
        printf("Syntax: gd_test (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      int analysisMethod = PSUADE_ANA_GOWER;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psuadeIO_->getParameter("ana_diagnostics",pPtr);
      ii = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,-1,-1);
      anaManager->analyze(psuadeIO_, 0, NULL, outputID);
      psuadeIO_->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
      delete anaManager;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rscreate 
    //**/ create response surface 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rscreate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rscreate: create response surface with loaded sample\n");
        printf("Syntax: rscreate (no argument needed)\n");
        printf("This command is useful for creating a response surface\n");
        printf("on the fly, and then use rseval new sample points.\n");
        printf("The sample to be used for rscreate should have already\n");
        printf("been loaded to local memory using 'load'.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      faFlag = 3;
      sprintf(pString, "Which output to use (1 - %d, 0 - all) ? ",nOutputs_);
      outputID = getInt(0, nOutputs_, pString);
      psuadeIO_->getParameter("ana_outputid", pPtr);
      kk = pPtr.intData_;
      faPtrsRsEval = new FuncApprox*[nOutputs_];
      if (outputID != 0)
      {
        psuadeIO_->updateAnalysisSection(-1,-1,-1,outputLevel_,
                                         outputID-1,-1);
        faPtrsRsEval[0] = genFAInteractive(psuadeIO_, faFlag);
        for (ii = 1; ii < nOutputs_; ii++) faPtrsRsEval[ii] = NULL;
      }
      else
      {
        for (ii = 0; ii < nOutputs_; ii++)
        {
          printf("Creating response surface for output %d\n", ii+1);
          psuadeIO_->updateAnalysisSection(-1, -1, -1, outputLevel_,ii,-1);
          faPtrsRsEval[ii] = genFAInteractive(psuadeIO_, faFlag);
        }
      }
      psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, kk, -1);
      printf("Now you can use the following to evaluate a new point: \n");
      printf("(1) use rseval with ivec_create/ivec_modify/ivec_show\n");
      printf("(2) use rseval and a file containing the new point\n");
      printf("(3) use rseval_m and a file containing new points\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rseval 
    //**/ evaluate the response surface at a given point
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rseval"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rseval: evaluate sample point after rscreate\n");
        printf("Syntax: rseval (no argument needed)\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command should be called after 'rscreate' to ");
      printf("evaluate one or\n");
      printf("more sample points from a register or from a file.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (faPtrsRsEval == NULL) 
      {
        printf("ERROR: response surface has not been created yet.\n");
        printf("       Please load sample and use rscreate first.\n");
        cmdStatus = 1;
        continue;
      }
      count = 0;
      sprintf(pString,
        "Sample taken from a file (n - sample from register)? (y or n) ");
      getString(pString, winput);
      psVector vecXVals;
      if (winput[0] == 'y')
      {
        printf("File format: \n");
        printf("PSUADE_BEGIN (optional)\n");
        printf("<nPts> <nInputs> \n");
        printf("1   <data> <data> ... \n");
        printf("2   <data> <data> ... \n");
        printf("... <data> <data> ... \n");
        printf("PSUADE_END \n");
        printf("Enter sample file for evaluation : ");
        scanf("%s", dataFile);
        fgets(lineIn,5000,stdin); 
        fp = fopen(dataFile, "r");
        if (fp == NULL)
        {
          printf("ERROR: sample file %s not found.\n", dataFile);
          cmdStatus = 1;
          continue;
        }
        fscanf(fp, "%s", winput);
        if (strcmp(winput, "PSUADE_BEGIN"))
        {
          printf("INFO: optional PSUADE_BEGIN keyword not found.\n");
          printf("      First line should have nSamples and nInputs\n");
          fclose(fp);
          fp = fopen(dataFile, "r");
        }
        fscanf(fp, "%d %d", &count, &kk);
        if (count <= 0)
        {
          printf("ERROR: invalid sample size\n");
          fclose(fp);
          cmdStatus = 1;
          continue;
        }
        if (kk != nInputs_)
        {
          printf("ERROR: input size does not match nInputs.\n");
          printf("       input size in local memory = %d.\n",nInputs_);
          printf("       input size from file       = %d.\n",kk);
          count = 0;
          fclose(fp);
          cmdStatus = 1;
          continue;
        }
        printf("Number of points to be evaluated = %d\n", count);
        //**/ detect comment line
        fgets(pString, 5000, fp);
        while (1)
        {
          kk = getc(fp);
          if (kk != '#')
          {
            ungetc(kk, fp);
            break;
          }
          else
          {
            fgets(pString, 5000, fp);
          }
        }
        vecXVals.setLength(count*nInputs_);
        double *XVals = vecXVals.getDVector(); 
        for (jj = 0; jj < count; jj++)
        {
          fscanf(fp, "%lg", &ddata);
          ind = (int) ddata;
          if (ind != (jj+1))
          {
            printf("ERROR: input index mismatch (%d,%d)\n",jj+1,ind);
            printf("       read     index = %d\n", ind);
            printf("       expected index = %d\n", jj+1);
            count = 0;
            break;
          }
          for (ii = 0; ii < nInputs_; ii++)
            fscanf(fp, "%lg", &(XVals[jj*nInputs_+ii]));
        }
        if (count > 0)
        {
          fscanf(fp, "%s", winput);
          if (strcmp(winput, "PSUADE_END"))
            printf("INFO: the optional PSUADE_END keyword is not found.\n");
        }
        fclose(fp);
      }
      else
      {
        vecXVals.setLength(nInputs_);
        if (VecDataReg_.length() == 0)
        {
          printf("ERROR: register has not been created yet.\n");
          cmdStatus = 1;
          continue;
        }
        count = 1;
        for (ii = 0; ii < nInputs_; ii++) 
          vecXVals[ii] = VecDataReg_[ii];
      }
      if (count > 0)
      {
        int flag = 0;
        printf("Fuzzy evaluation ? (y or n) ");
        fgets(winput,10,stdin); 
        if (winput[0] == 'y') flag = 1;
        printf("Evaluated sample written to a file ? (y or n) ");
        fgets(winput,10,stdin); 
        if (winput[0] == 'y')
        {
          fp = fopen("eval_sample.std","w");
          if (fp == NULL) 
          {
            printf("rseval ERROR: cannot open file eval_sample\n");
            continue;
          }
          //fprintf(fp, "# Inputs Out1 <Std1> Out2 <Std2> ...\n");
          if (nOutputs_ == 1 || (nOutputs_ > 1 && faPtrsRsEval[1] == NULL))
          {
            if (flag == 1) fprintf(fp, "%d %d 2\n", count, nInputs_);
            else           fprintf(fp, "%d %d 1\n", count, nInputs_);
          }   
          else
          {
            if (flag == 1)
                 fprintf(fp, "%d %d %d\n", count, nInputs_, nOutputs_);
            else fprintf(fp, "%d %d %d\n", count, nInputs_, 2*nOutputs_);
          }
        }
        else fp = NULL;
        vecYT.setLength(nOutputs_*count);
        vecWT.setLength(nOutputs_*count);
        double *dPtrY = vecYT.getDVector();
        double *dPtrW = vecWT.getDVector();
        for (ii = 0 ; ii < nOutputs_; ii++)
        {
          if (faPtrsRsEval[ii] != NULL)
          {
            if (flag == 1)
              dtemp = faPtrsRsEval[ii]->evaluatePointFuzzy(count,
                           vecXVals.getDVector(),&(dPtrY[ii*count]),
                           &(dPtrW[ii*count]));
            else
              dtemp = faPtrsRsEval[ii]->evaluatePoint(count,
                           vecXVals.getDVector(),&(dPtrY[ii*count]));
          }
        }
        for (kk = 0 ; kk < count; kk++)
        {
          if (fp != NULL)
            for (ii = 0 ; ii < nInputs_; ii++)
              fprintf(fp, "%e ", vecXVals[kk*nInputs_+ii]);
          printf("Interpolated Point %d: ", kk+1);
          for (ii = 0 ; ii < nOutputs_; ii++)
          {
            if (faPtrsRsEval[ii] != NULL)
            {
              if (flag == 1)
              {
                printf("output %d = %e (stdev = %e) ",
                       ii+1,dPtrY[ii*count+kk], dPtrW[ii*count+kk]);
                if (fp != NULL)
                  fprintf(fp,"%e %e ",dPtrY[ii*count+kk],
                          dPtrW[ii*count+kk]);
              }
              else
              {
                printf("output %d = %e ", ii+1,dPtrY[ii*count+kk]);
                if (fp != NULL) fprintf(fp,"%e ",dPtrY[ii*count+kk]);
              }
            }
          }
          printf("\n");
          if (fp != NULL) fprintf(fp, "\n");
        }
        if (fp != NULL) 
        {
          printf("Evaluated sample file in 'eval_sample.std'\n");
          fclose(fp);
          fp = NULL;
          cmdStatus = 0;
        }
      }
    }

    //**/ -------------------------------------------------------------
    // +++ rseval_m 
    //**/ evaluate the response surface at given points repeatedly
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rseval_m"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rseval_m: evaluate points repeatedly after rscreate.\n");
        printf("          Data taken from a file with format:\n");
        printf("PSUADE_BEGIN \n");
        printf("<number of points> <nInputs> \n");
        printf("1  <input data 1> <input data 2> .. \n");
        printf("2  <input data 1> <input data 2> .. \n");
        printf("PSUADE_END \n");
        printf("Syntax: rseval_m <filename>)\n");
      }
      if (!strcmp(winput, "-h")) continue;
      printAsterisks(PL_INFO, 0);
      printf("This command should be called after 'rscreate' to ");
      printf("evaluate one or more\n");
      printf("sample points from a file.\n");
      printf("NOTE: This command enables PSUADE to be a response ");
      printf("surface server. So\n");
      printf("      first run this command in another terminal. ");
      printf("Then run your client\n");
      printf("      program that needs response surface ");
      printf("evaluation from the server.\n");
      printf("      (NOTE: the server and client must be run ");
      printf("in the same directory.)\n");
      printf("      The client program must put new points to ");
      printf("be evaluated in a file\n");
      printf("      called rsevalDataIn.1. PSUADE, upon ");
      printf("detection of this file,\n");
      printf("      reads this file, evaluates the sample, and ");
      printf("puts the results in\n");
      printf("      the rsEvalDataOut.1 file. You can then ");
      printf("continue with the next\n");
      printf("      set of new points in rsEvalDataIn.2 and ");
      printf("so on. When you are done\n");
      printf("      with the response surface, simply create ");
      printf("an empty file called\n");
      printf("      psComplete in the same run directory for ");
      printf("graceful termination.\n");
      printf("The required format of the rsEvalDataIn.x file: \n");
      printf("Line 1: <nPoints> <nInputs>\n");
      printf("Line 2: 1  <input data for sample 1\n");
      printf("Line 3: 2  <input data for sample 2.\n");
      printf("...\n");

      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (faPtrsRsEval == NULL) 
      {
        printf("ERROR: response surface has not been created yet.\n");
        printf("       Please load sample and use rscreate first.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(dataFile, "psComplete");
      if ((fp=fopen(dataFile,"r")) != NULL)
      {
        printf("Please remove psComplete file and re-do.\n");
        cmdStatus = 1;
        continue;
      }
      int count2, terminate=0;
      vecWT.setLength(nInputs_);
      FILE   *fpOut=NULL, *fpTmp;
      count = 1;
      while ((fp=fopen(dataFile,"r")) == NULL)
      {
        sprintf(winput, "rsevalDataIn.%d", count);
        printf("Now expecting a sample file called %s (or psComplete).\n", 
               winput);
        while ((fp=fopen(winput,"r")) == NULL)
        {
#ifdef WINDOWS
          Sleep(1000);
#else
          sleep(1);
#endif
          if ((fpTmp=fopen("psComplete","r")) != NULL)
          {
            fclose(fpTmp);
            printf("INFO: psComplete file detected ==> terminate.\n");
            terminate = 1;
            break;
          }
        }
        //fscanf(fp, "%s", winput);
        //if (strcmp(winput, "PSUADE_BEGIN"))
        //{
        //  printf("INFO: Optional PSUADE_BEGIN keyword not found.\n");
        //  printf("      First line should have nSamples and nInputs\n");
        //  fclose(fp);
        //  fp = fopen(winput, "r");
        //}
        if (terminate == 1) 
        {
          cmdStatus = 0;
          continue;
        }
        fscanf(fp, "%d %d", &count2, &kk);
        printf("   numData = %d, nInputs = %d\n",count2,kk);
	if (count2 <= 0 || kk != nInputs_)
        {
          printf("ERROR: invalid sample data in file.\n");
          fclose(fp);
          cmdStatus = 1;
          continue;
        }
        sprintf(winput, "rsevalDataOut.%d", count);
        if ((fpOut=fopen(winput,"w")) == NULL)
        {
          printf("ERROR: cannot open output file.\n");
          fclose(fp);
          cmdStatus = 1;
          continue;
        }
        for (ss = 0; ss < count2; ss++)
        {
          fscanf(fp, "%d", &ind);
          if (ind != ss+1)
          {
            printf("ERROR: invalid data in file.\n");
            printf("       sample number = %d\n", ss+1);
            fclose(fp);
            break;
          }
          printf("Sample %4d: input = ",ss+1);
          for (ii = 0; ii < nInputs_; ii++)
          {
            fscanf(fp, "%lg", &ddata);
            vecWT[ii] = ddata;
            printf(" %14.6e ",ddata);
          }
          printf("\n");
          for (ii = 0; ii < nOutputs_; ii++)
          {
            if (faPtrsRsEval[ii] != NULL)
            {
              dtemp = faPtrsRsEval[ii]->evaluatePoint(vecWT.getDVector());
              fprintf(fpOut, "%e\n", dtemp);
              printf("        Evaluated output %3d = %14.6e\n",
                     ii+1,dtemp);
            } 
          } 
        } 
        fclose(fpOut);
        if (ss != count2)
        {
          cmdStatus = 1;
          continue;
        }
        fscanf(fp, "%s", winput);
        fscanf(fp, "%s", winput);
        //if (strcmp(winput, "PSUADE_END"))
        //  printf("INFO: the optional PSUADE_END keyword is not found.\n");
        fclose(fp);
        printf("==> Evaluated sample is now in rsevalDataOut.%d\n",count);
        printf("==> If done, create an empty psComplete file to signal");
        printf(" completion.\n");
        count++;
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsevaluate
    //**/ evaluate the response surface at a given point
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsevaluate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsevaluate: evaluate response surface\n");
        printf("    NOTE: An evaluation sample will be needed.\n");
        printf("    NOTE: It does not take discrepancy function.\n");
        printf("Syntax: rsevaluate (no argument needed)\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ get more information
      sprintf(pString,"Enter output number (1-based) : ");
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      //**/ create response surface ==> faPtr
      psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID, -1);
      faFlag = 3;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) {printf("ERROR detected in RS.\n"); continue;}
      faPtr->setOutputLevel(outputLevel_);

      //**/ read sample file to be evaluated ==> rsevalnSampes, vecTX
      printf("Now enter the name of your sample file (to be evaluated)\n");
      printf("Expected file format (nInputs should be %d): \n",nInputs_);
      printf("PSUADE_BEGIN (optional)\n");
      printf("<nPts> <nInputs> \n");
      printf("1   <data> <data> ... \n");
      printf("2   <data> <data> ... \n");
      printf("... <data> <data> ... \n");
      printf("PSUADE_END (optional)\n");
      sprintf(pString, "File name ? ");
      getString(pString, dataFile);
      kk = strlen(dataFile);
      dataFile[kk-1] = '\0';
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: sample file %s not found.\n", dataFile);
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 1;
        continue;
      }
      fscanf(fp, "%s", winput);
      int rsevalnSamps;
      if (strcmp(winput, "PSUADE_BEGIN"))
      {
        printf("INFO: optional PSUADE_BEGIN keyword not found.\n");
        printf("      First line should have nSamples and nInputs\n");
        fclose(fp);
        fp = fopen(dataFile, "r");
      }
      fscanf(fp, "%d %d", &rsevalnSamps, &kk);
      if (rsevalnSamps <= 0)
      {
        printf("ERROR: invalid sample size\n");
        fclose(fp);
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 1;
        continue;
      }
      if (kk != nInputs_)
      {
        printf("ERROR: input size does not match nInputs.\n");
        printf("       input size from file       = %d.\n",kk);
        printf("       input size expected        = %d.\n",nInputs_);
        fclose(fp);
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 1;
        continue;
      }
      printf("Number of points to be evaluated = %d\n",rsevalnSamps);
      //**/ detect comment line
      fgets(pString, 5000, fp);
      while (1)
      {
        kk = getc(fp);
        if (kk != '#')
        {
          ungetc(kk, fp);
          break;
        }
        else
        {
          fgets(pString, 5000, fp);
        }
      }
      psVector vecTX;
      vecTX.setLength(rsevalnSamps*nInputs_);
      for (ss = 0; ss < rsevalnSamps; ss++)
      {
        fscanf(fp, "%lg", &ddata);
        ind = (int) ddata;
        if (ind != (ss+1))
        {
          printf("ERROR: input index mismatch (%d,%d)\n",ss+1,ind);
          printf("       read     index = %d\n", ind);
          printf("       expected index = %d\n", ss+1);
          delete faPtr;
          faPtr = NULL;
          break;
        }
        for (ii = 0; ii < nInputs_; ii++) 
        {
          fscanf(fp,"%lg", &dtemp);
          vecTX[ss*nInputs_+ii] = dtemp;
        }
      }
      fclose(fp);
      if (faPtr == NULL) 
      {
        cmdStatus = 1;
        continue;
      }

      //**/ open an evaluation output file
      fp = fopen("rsevaluate.out", "w");
      if (fp == NULL)
      {
        printf("ERROR: rsevaluate - cannot open output file.\n");
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 1;
        continue;
      }

      //**/ now we have faPtr, rsevalnSamps, vecTX
      psVector vecTY, vecTW;
      vecTY.setLength(rsevalnSamps);
      vecTW.setLength(rsevalnSamps);
      faPtr->evaluatePointFuzzy(rsevalnSamps,vecTX.getDVector(),
                   vecTY.getDVector(),vecTW.getDVector());
      for (ss = 0; ss < rsevalnSamps; ss++)
        fprintf(fp, "%24.16e %24.16e\n",vecTY[ss],vecTW[ss]);
      printf("Evaluated sample file in 'rsevaluate.out'\n");
      fclose(fp);
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ mcmc_predict 
    //**/ evaluate the response surface at a given point in a special
    //**/ way
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc_predict") ||
             !strcmp(command, "rsevaluate2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc_predict: prediction when some input ");
        printf("parameters are to be used as\n");
        printf("    uncertain parameters (i.e. prediction ");
        printf("uncertainty is assumed to be\n");
        printf("    due primarily to uncertain parameters and ");
        printf("not to RS fitting error.\n");
        printf("    This command is especially useful for ");
        printf("prediction after MCMC (i.e.\n");
        printf("    MCMCPostSample and optionally psDiscrepancyModel ");
        printf("are available).\n");
        printf("Syntax: mcmc_predict (no argument needed)\n");
        continue;
      }
      //**/-----------------------------------------------------
      //**/ give instructions
      //**/-----------------------------------------------------
      printf("The steps in executing this command is:\n");
      printf("0. Load the base sample (should have been done)\n");
      printf("1. Select which inputs are uncertain inputs\n");
      printf("2. If discrepancy is to be added, enter the model file\n");
      printf("3. Load the evaluation sample (no uncertain inputs)\n");
      printf("4. Load sample for uncertain inputs (e.g. from MCMC)\n");
      printf("5. Select RS types for base and discrepancy models\n"); 
      printf("6. Wait for answers in the mcmcpredict.out file.\n");
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      //**/-----------------------------------------------------
      //**/ error checking
      //**/-----------------------------------------------------
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ == 1)
      {
        printf("ERROR: This command only works with nInputs > 1\n");
        printf("       when at least 1 parameter is uncertain.\n");
        cmdStatus = 1;
        continue;
      }
      if (nOutputs_ > 1)
      {
        printf("ERROR: This command only works with nOutputs == 1\n");
        printf("SUGGESTION: use odelete to remove all but 1 output.\n");
        cmdStatus = 1;
        continue;
      }

      //**/-----------------------------------------------------
      //**/ select parameter type ==> vecISet, nUParams
      //**/-----------------------------------------------------
      psIVector vecISet;
      vecISet.setLength(nInputs_);
      printf("\n1. Select between the following two types of inputs:\n");
      printf("(a) regular (or design) inputs, and\n");
      printf("(b) uncertain inputs (with an uncertain sample)\n");
      printf("That is, you need to specify inputs are uncertain.\n\n");
      sprintf(pString,"How many inputs are uncertain inputs ? ");
      int nUParams = getInt(1, nInputs_, pString);
      printf("Now specify which inputs are uncertain inputs.\n");
      for (ii = 0; ii < nUParams; ii++)
      {
        sprintf(pString,"Enter uncertain input number (1-based) : ");
        kk = getInt(1, nInputs_, pString);
        kk--; 
        if (vecISet[kk] != 0)
          printf("ERROR: input %d has been selected already.\n",kk+1);
        while (vecISet[kk] != 0)
        {
          kk = getInt(1, nInputs_, pString);
          kk--; 
        }
        vecISet[kk] = 1;
      }

      //**/-----------------------------------------------------
      //**/ ask for discrepancy model, if any; and check
      //**/ ==> discFile, ioPtr, nInps, nOuts
      //**/-----------------------------------------------------
      int discFile=1, nInps, nOuts;
      printf("\n2. Enter discrepancy sample file (PSUADE format)\n");
      printf("Discrepancy model file name? (enter 'none' if none): ");
      scanf("%s", winput);
      fgets(lineIn,5000,stdin);
      PsuadeData *ioPtr = NULL;
      if (!strcmp(winput, "none")) discFile = 0;
      else
      {
        ioPtr = new PsuadeData();
        status = ioPtr->readPsuadeFile(winput);
        if (status == 0)
        {
          ioPtr->getParameter("input_ninputs", pPtr);
          nInps = pPtr.intData_;
          if (nInps != nInputs_-nUParams)
          {
            printf("ERROR: discrepancy model should have %d inputs,\n",
                   nInputs_-nUParams);
            printf("       but it has %d inputs instead.\n",nInps);
            delete ioPtr;
            ioPtr = NULL;
            continue;
          }
          ioPtr->getParameter("output_noutputs", pPtr);
          nOuts = pPtr.intData_;
          if (nOuts > 1)
          {
            printf("ERROR: Discrepancy model has nOutputs > 1.\n");
            printf("       This is currently not supported.\n");
            printf("Suggestion: modify the discrepancy file.\n");
            delete ioPtr;
            ioPtr = NULL;
            cmdStatus = 1;
            continue;
          }
        }
        else
        {
          printf("ERROR: in reading the discrepancy model file %s.\n",
                 winput);
          printf("Suggestion: use load to check the format is fine.\n");
          discFile = 0;
          cmdStatus = 1;
          continue;
        }
      }

      //**/-----------------------------------------------------
      //**/ read sample for evaluation => mcmcpredictnSamps, vecDX
      //**/-----------------------------------------------------
      printf("\n3. Enter the name of the file that contains SAMPLES ");
      printf("TO BE EVALUATED\n");
      printf("   with the associated uncertainties induced by the ");
      printf("uncertain inputs.\n");
      printf("   As such, the sample points in this file should ");
      printf("contain only regular\n");
      printf("   (or design) inputs, NOT uncertain inputs (that is, ");
      printf("the nInputs in\n");
      printf("   this file should be = %d)\n",nInputs_-nUParams);
      printf("Expected file format: \n");
      printf("PSUADE_BEGIN (optional)\n");
      printf("<nPts> <nInputs> \n");
      printf("1   <data> <data> ... \n");
      printf("2   <data> <data> ... \n");
      printf("... <data> <data> ... \n");
      printf("PSUADE_END (optional)\n");
      sprintf(pString, "Enter sample file name ? ");
      getString(pString, dataFile);
      kk = strlen(dataFile);
      dataFile[kk-1] = '\0';
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: evaluation sample file %s not found.\n", dataFile);
        cmdStatus = 1;
        continue;
      }
      fscanf(fp, "%s", winput);
      if (strcmp(winput, "PSUADE_BEGIN"))
      {
        printf("INFO: optional PSUADE_BEGIN keyword not found.\n");
        printf("      First line should have nSamples and nInputs\n");
        fclose(fp);
        fp = fopen(dataFile, "r");
      }
      int mcmcpredictnSamps;
      fscanf(fp, "%d %d", &mcmcpredictnSamps, &kk);
      if (mcmcpredictnSamps <= 0)
      {
        printf("ERROR: invalid evaluation sample size <= 0\n");
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      if (kk != nInputs_-nUParams)
      {
        printf("ERROR: evaluation sample nInputs is not valid.\n");
        printf("          input size from file       = %d.\n",kk);
        printf("          input size expected        = %d.\n",
               nInputs_-nUParams);
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      printf("Evaluation sample size = %d\n",mcmcpredictnSamps);
      //**/ detect comment line
      fgets(pString, 5000, fp);
      while (1)
      {
        kk = getc(fp);
        if (kk != '#')
        {
          ungetc(kk, fp);
          break;
        }
        else fgets(pString, 5000, fp);
      }
      psVector vecDX;
      vecDX.setLength(mcmcpredictnSamps*nInputs_);
      for (ss = 0; ss < mcmcpredictnSamps; ss++)
      {
        fscanf(fp, "%d", &ind);
        if (ind != (ss+1))
        {
          printf("ERROR: evaluation sample input index mismatch.\n");
          printf("          index read from file = %d\n", ind);
          printf("          expected index       = %d\n", ss+1);
          printf("Check your input sample file format: \n");
          printf("PSUADE_BEGIN (optional)\n");
          printf("<nPts> <nInputs> \n");
          printf("1   <data> <data> ... \n");
          printf("2   <data> <data> ... \n");
          printf("... <data> <data> ... \n");
          printf("PSUADE_END (optional)\n");
          break;
        }
        for (ii = 0; ii < nInputs_-nUParams; ii++) 
        {
          fscanf(fp,"%lg", &ddata);
          vecDX[ss*(nInputs_-nUParams)+ii] = ddata;
        }
      }
      fclose(fp);
      if (ss != mcmcpredictnSamps) continue;

      //**/ read uncertain parameter sample ==> nUSamples, vecTU
      printf("\n4. Enter sample file name for the uncertain ");
      printf("inputs. This should\n");
      printf("   normally be the MCMCPostSample file created ");
      printf("by the rsmcmc command.\n");
      printf("Expected file format (nInputs should be %d): \n",nUParams);
      printf("PSUADE_BEGIN (optional)\n");
      printf("<nPts> <nInputs> \n");
      printf("1   <data> <data> ... \n");
      printf("2   <data> <data> ... \n");
      printf("... <data> <data> ... \n");
      printf("PSUADE_END (optional)\n");
      sprintf(pString, "File name ? ");
      getString(pString, dataFile);
      kk = strlen(dataFile);
      dataFile[kk-1] = '\0';
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: uncertain input sample file %s not found.\n",
               dataFile);
        continue;
      }
      fscanf(fp, "%s", winput);
      if (strcmp(winput, "PSUADE_BEGIN"))
      {
        printf("INFO: optional PSUADE_BEGIN keyword not found.\n");
        printf("      First line should have nSamples and nInputs\n");
        fclose(fp);
        fp = fopen(dataFile, "r");
      }
      int nUSamples, nUInpsChk;
      fscanf(fp, "%d %d", &nUSamples, &nUInpsChk);
      if (nUSamples <= 0)
      {
        printf("ERROR: invalid uncertain input sample size (<= 0)\n");
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      if (nUInpsChk < nUParams)
      {
        printf("ERROR: invalid uncertain input sample nInputs.\n");
        printf("       input size from file       = %d.\n",nUInpsChk);
        printf("       input size expected       >= %d.\n",nUParams);
        fclose(fp);
        cmdStatus = 1;
        continue;
      }
      printf("Uncertain input sample size = %d\n", nUSamples);
      psIVector vecITrack;
      vecITrack.setLength(nUInpsChk);
      if (nUInpsChk > nUParams)
      {
        printf("WARNING: nInputs in the uncertain sample file (%d)\n",
               nUInpsChk);
        printf("         is larger than the number of declared\n");
        printf("         uncertain inputs %d.\n", nUParams);
        printf("You must delete %d inputs for compatibility.\n",
               nUInpsChk-nUParams);
        sprintf(pString, "Which input to delete (1-%d)? ",nUInpsChk);
        ii = 0;
        while (ii < nUInpsChk-nUParams)
        {
          jj = getInt(1, nUInpsChk, pString);
          if (vecITrack[jj-1] == 0)
          {
            vecITrack[jj-1] = 1;
            ii++;
          }
          else printf("ERROR: input %d has already been selected.",jj);
        } 
      } 

      //**/ detect comment line
      fgets(pString, 5000, fp);
      while (1)
      {
        kk = getc(fp);
        if (kk != '#')
        {
          ungetc(kk, fp);
          break;
        }
        else
        {
          fgets(pString, 5000, fp);
        }
      }
      psVector vecTU;
      vecTU.setLength(nUSamples*nUParams);
      for (ss = 0; ss < nUSamples; ss++)
      {
        fscanf(fp, "%lg", &ddata);
        ind = (int) ddata;
        if (ind != (ss+1))
        {
          printf("ERROR: uncertain input sample index mismatch:\n");
          printf("       read     sample number = %d\n", ind);
          printf("       expected sample number = %d\n", ss+1);
          count = 0;
          break;
        }
        kk = 0;
        for (ii = 0; ii < nUInpsChk; ii++) 
        {
          fscanf(fp,"%lg", &ddata);
          if (vecITrack[ii] == 0)
          {
            vecTU[ss*nUParams+kk] = ddata;
            kk++;
          }
        }
      }
      fclose(fp);

      //**/-----------------------------------------------------
      //**/ create response surface ==> faPtr
      //**/-----------------------------------------------------
      printf("\n5a. Building RS for the base sample\n");
      faFlag = -1;
      FuncApprox *faPtr = genFA(faFlag, nInputs_, iOne, nSamples_);
      if (faPtr == NULL) 
      {
        printf("ERROR: when building RS for the base sample.\n"); 
        continue;
      }
      faPtr->setBounds(VecILowerBs_.getDVector(),
                       VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      faPtr->initialize(VecSamInputs_.getDVector(),
                        VecSamOutputs_.getDVector());

      //**/-----------------------------------------------------
      //**/ set up discrepancy model, if any
      //**/ and evaluate at mcmcpredictnSamps of vecDX => vecDY
      //**/-----------------------------------------------------
      psVector vecDY;
      vecDY.setLength(mcmcpredictnSamps);
      if (discFile == 1)
      {
        printf("\n5b. Building RS for the discrepancy sample\n");
        FuncApprox *discPtr = genFAInteractive(ioPtr, 3);
        if (discPtr == NULL) 
        {
          printf("ERROR: when building RS for the discrepancy sample.\n"); 
          continue;
        }
        discPtr->evaluatePoint(mcmcpredictnSamps,vecDX.getDVector(), 
                               vecDY.getDVector());
        delete ioPtr;
        ioPtr = NULL;
        delete discPtr;
        discPtr = NULL;
      }

      //**/-----------------------------------------------------
      //**/ open an evaluation output file
      //**/-----------------------------------------------------
      fp = fopen("mcmcpredict.out", "w");
      if (fp == NULL)
      {
        printf("ERROR: cannot open output file to write results.\n");
        if (faPtr != NULL) 
        {
          delete faPtr;
          faPtr = NULL;
        }
        cmdStatus = 1;
        continue;
      }

      //**/-----------------------------------------------------
      //**/ now we have vecISet, faPtr
      //**/ mcmcpredictnSamps, vecDY for discrepancy
      //**/ nUSamples, nUParams, vecTU for uncertain sample
      //**/ next create a evaluation vector V and result vector Y
      //**/-----------------------------------------------------
      printf("\n6. Sample evaluation proceeds ... \n");
      printf("     Number of points to evaluate   = %d\n",
             mcmcpredictnSamps);
      printf("     Number of evaluations for each = %d\n",nUSamples);
      printf("\n   (total sample size = %d)\n", 
             mcmcpredictnSamps*nUSamples);
      psVector vecTX, vecTY;
      vecTX.setLength(nUSamples*nInputs_);
      vecTY.setLength(nUSamples);
      int count2;
      printf("Sample      mean        std. dev.    ");
      printf("Min          Max\n");
      fprintf(fp,"mean            std. deviation   ");
      fprintf(fp,"Min              Max\n");
      for (ss = 0; ss < mcmcpredictnSamps; ss++)
      {
        //**/ inject uncertain sample and evaluate via faPtr
        //**/ then add discrepancy vecDY
        for (kk = 0; kk < nUSamples; kk++)
        {
          count = count2 = 0;
          for (ii = 0; ii < nInputs_; ii++)
          {
            if (vecISet[ii] == 0)
            {
              vecTX[kk*nInputs_+ii] = 
                       vecDX[ss*(nInputs_-nUParams)+count2];
              count2++;
            }
            else
            {
              vecTX[kk*nInputs_+ii] = vecTU[kk*nUParams+count];
              count++;
            }
          }
        }
        faPtr->evaluatePoint(nUSamples,vecTX.getDVector(),
                             vecTY.getDVector());
        ddata = 0.0;
        double dmax = -PSUADE_UNDEFINED;
        double dmin = +PSUADE_UNDEFINED;
        for (kk = 0; kk < nUSamples; kk++)
        {
          ddata += vecTY[kk];
          if (vecTY[kk] > dmax) dmax = vecTY[kk];
          if (vecTY[kk] < dmin) dmin = vecTY[kk];
        }
        ddata /= (double) nUSamples;

        dtemp = 0.0;
        for (kk = 0; kk < nUSamples; kk++)
          dtemp += (vecTY[kk] - ddata) * (vecTY[kk] - ddata);
        fprintf(fp,"%16.8e %16.8e %16.8e %16.8e\n",ddata+vecDY[ss],
                sqrt(dtemp/(double) nUSamples),dmin+vecDY[ss],
                dmax+vecDY[ss]);
        printf("%7d  %12.4e %12.4e %12.4e %12.4e\n",ss+1,
               ddata+vecDY[ss],sqrt(dtemp/(double) nUSamples),
               dmin+vecDY[ss],dmax+vecDY[ss]);
      }
      printf("Evaluated sample is now in 'mcmcpredict.out'\n");
      fclose(fp);
      if (faPtr != NULL)
      {
        delete faPtr;
        faPtr = NULL;
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rs_splot 
    //**/ probe response surface and create scatter plots
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rs_splot"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rs_splot: create a RS-based output scatter plot\n");
        printf("Syntax: rs_splot (no argument needed)\n");
        printf("This command is similar to the splot command, except\n");
        printf("that the scatter plot is created from response surface.\n");
        printf("This command is useful when the sample size is small\n");
        printf("and you would like to see the trend with scatter plots.\n");
        printf("The idea is to create a large sample for scatter plots\n");
        printf("using response surfaces from small data set.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;
      Ymax = - 1.0e35;
      Ymin =   1.0e35;
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (VecSamOutputs_[sInd*nOutputs_+outputID] > Ymax)
          Ymax = VecSamOutputs_[sInd*nOutputs_+outputID];
        if (VecSamOutputs_[sInd*nOutputs_+outputID] < Ymin)
          Ymin = VecSamOutputs_[sInd*nOutputs_+outputID];
      }
      printf("Option to clip the created sample by bounding Y.\n");
      sprintf(pString,"Enter the lower bound constraint (Ymin=%e) : ",
              Ymin);
      threshL = getDouble(pString);
      sprintf(pString,"Enter the upper bound constraint (Ymax=%e) : ",
              Ymax);
      threshU = getDouble(pString);
      if (threshL >= threshU)
      {
        printf("ERROR: lower bound >= upper bound.\n");
        cmdStatus = 1;
        continue;
      }
      //**/ create response surface
      faFlag = 3;
      psuadeIO_->updateAnalysisSection(-1, -1, -1, -1, outputID,-1); 
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      faPtr->setOutputLevel(outputLevel_);
      //**/ create a large sample
      sprintf(pString,
            "Enter Monte Carlo sample size for probing (max=100000): ");
      kk = getInt(100, 100000, pString);
      Sampling *sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MC);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nInputs_, VecILowerBs_.getDVector(), 
                              VecIUpperBs_.getDVector());
      sampPtr->setOutputParams(1);
      sampPtr->setSamplingParams(kk, 1, 0);
      sampPtr->initialize(0);
      vecXT.setLength(kk*nInputs_);
      vecYT.setLength(kk);
      vecST.setLength(kk);
      double *dPtrX = vecXT.getDVector();
      for (ii = 0; ii < kk; ii++) vecST[ii] = 1;
      sampPtr->getSamples(kk,nInputs_,1,vecXT.getDVector(),
                          vecYT.getDVector(), vecST.getIVector());
      count = 0;
      for (ii = 0; ii < kk; ii++)
      {
        dtemp = faPtr->evaluatePoint(&dPtrX[ii*nInputs_]);
        vecYT[ii] = dtemp;
        vecST[ii] = 1;
        if (dtemp >= threshL && dtemp <= threshU) count++; 
        else                                      vecST[ii] = 0;
      } 
      //**/sprintf(pString,"Generate 1D or 2D scatter plots (1 or 2): ");
      //**/getString(pString, winput);
      winput[0] = '1';
      if (winput[0] == '1')
      {
        //**/ write to matlab/scilab file
        if (plotScilab()) fp = fopen("scilabrssp.sci", "w");
        else              fp = fopen("matlabrssp.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabrssp.m.\n");
          cmdStatus = 1;
          continue;
        }
        fprintf(fp, "Y = [\n");
        for (sInd = 0; sInd < kk; sInd++)
          if (vecST[sInd] == 1) fprintf(fp, "%e\n",vecYT[sInd]);
        fprintf(fp, "];\n");
        for (iInd = 0; iInd < nInputs_; iInd++)
        {
          fprintf(fp, "X%d = [\n", iInd+1);
          for (sInd = 0; sInd < kk; sInd++)
            if (vecST[sInd] == 1)
              fprintf(fp, "%e\n",dPtrX[sInd*nInputs_+iInd]);
          fprintf(fp, "];\n");
        }
        for (iInd = 0; iInd < nInputs_; iInd++)
        {
          fwritePlotCLF(fp);
          fprintf(fp, "plot(X%d,Y,'X','MarkerSize',13)\n",iInd+1);
          sprintf(pString, "%s vs %s",
                  StrOutNames_[outputID],StrInpNames_[iInd]);
          fwritePlotTitle(fp, pString);
          fwritePlotXLabel(fp, StrInpNames_[iInd]);
          fwritePlotYLabel(fp, StrOutNames_[outputID]);
          fwritePlotAxes(fp);
          if (iInd < nInputs_-1) 
          {
            fprintf(fp, "disp(\'Press enter to go to the next plot\')\n");
            fprintf(fp, "pause\n");
          }
        }
        fclose(fp);    
        if (plotScilab())
             printf("scilabrssp.sci is now ready for scatter plots.\n");
        else printf("matlabrssp.m is now ready for scatter plots.\n");
      }
      delete faPtr;
      delete sampPtr;
      faPtr = NULL;
      sampPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ interactive 
    //**/ set interactive mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "interactive"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off interactive mode.\n");
        printf("Interactive mode, when turned on, will trigger many\n");
        printf("prompts for additonal information during the course\n");
        printf("of setup and analysis.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  psConfig_.InteractiveOff();
      if (!strcmp(winput, "off")) psConfig_.InteractiveOn();
      if (!psConfig_.InteractiveIsOn())
      {
        psConfig_.InteractiveOn();
        printf("interactive mode on\n");
      }
      else
      {
        psConfig_.InteractiveOn();
        printf("interactive mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ dontask 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "dontask"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off interactive mode.\n");
        printf("Interactive mode, when turned on, will trigger many\n");
        printf("prompts for additonal information during the course\n");
        printf("of setup and analysis.\n");
        continue;
      }
      if (!psConfig_.InteractiveIsOn()) psConfig_.InteractiveOn();
      else                              psConfig_.InteractiveOff();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ io_expert 
    //**/ set IO expert mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "io_expert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off I/O expert mode.\n");
        printf("This mode, when turned on, will trigger prompts for\n");
        printf("additonal information regarding data input/output.\n");
        printf("When this mode is off, default setting will be used.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.IOExpertModeOn();
        printf("IO expert mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.IOExpertModeOff();
        printf("IO expert mode off\n");
      }
      else if (!psConfig_.IOExpertModeIsOn())
      {
        psConfig_.IOExpertModeOn();
        printf("IO expert mode on\n");
      }
      else
      {
        psConfig_.IOExpertModeOff();
        printf("IO expert mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rs_expert 
    //**/ set response surface expert mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rs_expert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off response surface expert mode.\n");
        printf("This mode, when turned on, will trigger prompts for\n");
        printf("additional information regarding creation of response\n");
        printf("surfaces.\n");
        printf("When this mode is off, default setting will be used.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.RSExpertModeOn();
        printf("response surface expert mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.RSExpertModeOff();
        printf("response surface expert mode off\n");
      }
      else if (!psConfig_.RSExpertModeIsOn())
      {
        psConfig_.RSExpertModeOn();
        printf("response surface expert mode on\n");
      }
      else
      {
        psConfig_.RSExpertModeOff();
        printf("response surface expert mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rs_codegen 
    //**/ set response surface code generator mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rs_codegen"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off response surface code\n");
        printf("generation mode. This mode, when turned on, will\n");
        printf("generate a psuade_rs.info file whenever a response\n");
        printf("surface check ('rsfit') is performed. This file\n");
        printf("contains information for constructing the stand-alone\n");
        printf("response surface interpolators apart from PSUADE.\n");
        printf("In addition, a psuade_rs.py file will also be created\n");
        printf("as a stand-alone Python interpolator.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.RSCodeGenOn();
        printf("response surface code generation mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.RSCodeGenOff();
        printf("response surface code generation mode off\n");
      }
      else if (!psConfig_.RSCodeGenIsOn())
      {
        psConfig_.RSCodeGenOn();
        printf("response surface code generation mode on\n");
      }
      else
      {
        psConfig_.RSCodeGenOff();
        printf("response surface code generation mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ana_expert 
    //**/ set analysis expert mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "ana_expert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off data analysis expert mode.\n");
        printf("This mode, when turned on, will trigger prompts for\n");
        printf("additional information regarding data analysis.\n");
        printf("When this mode is off, default setting will be used.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.AnaExpertModeOn();
        printf("analysis expert mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.AnaExpertModeOff();
        printf("analysis expert mode off\n");
      }
      else if (!psConfig_.AnaExpertModeIsOn())
      {
        psConfig_.AnaExpertModeOn();
        printf("analysis expert mode on\n");
      }
      else
      {
        psConfig_.AnaExpertModeOff();
        printf("analysis expert mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sam_expert 
    //**/ set sampling expert mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sam_expert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off sample design expert mode.\n");
        printf("This mode, when turned on, will trigger prompts for\n");
        printf("additional information regarding creation of samples.\n");
        printf("When this mode is off, default setting will be used.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.SamExpertModeOn();
        printf("sampling expert mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.SamExpertModeOff();
        printf("sampling expert mode off\n");
      }
      else if (!psConfig_.SamExpertModeIsOn())
      {
        psConfig_.SamExpertModeOn();
        printf("sampling expert mode on\n");
      }
      else
      {
        psConfig_.SamExpertModeOff();
        printf("sampling expert mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ opt_expert 
    //**/ set optimization expert mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "opt_expert"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off optimization expert mode.\n");
        printf("This mode, when turned on, will trigger prompts for\n");
        printf("additonal information regarding numerical optimization.\n");
        printf("When this mode is off, default setting will be used.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  
      {
        psConfig_.OptExpertModeOn();
        printf("optimization expert mode on\n");
      }
      else if (!strcmp(winput, "off")) 
      {
        psConfig_.OptExpertModeOff();
        printf("optimization expert mode off\n");
      }
      else if (!psConfig_.OptExpertModeIsOn())
      {
        psConfig_.OptExpertModeOn();
        printf("optimization expert mode on\n");
      }
      else
      {
        psConfig_.OptExpertModeOff();
        printf("optimization expert mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ genhistogram 
    //**/ generate histogram from a psuade sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "genhistogram"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command uses the loaded sample and generates\n");
        printf("a sample of scenarios each associated with a\n");
        printf("probability. This command is useful for building\n");
        printf("scenarios for optimization under uncertainty.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      vecIT.setLength(nInputs_);
      printf("Please enter the number of bins per input dimension.\n");
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString,"Number of histogram bins for input %d (1-1000): ",
                ii+1);
        vecIT[ii] = getInt(1, 1000, pString);
      }
      PDFHistogram *pdfhist = new PDFHistogram(nSamples_, nInputs_, 
                   VecSamInputs_.getDVector(),vecIT.getIVector(),1);
      vecXT.setLength(nInputs_);
      pdfhist->getMeans(vecXT.getDVector());
      for (ii = 0; ii < nInputs_; ii++)
        printf("Estimated mean of input %7d = %e\n",ii+1,vecXT[ii]);
      count = 100000;
      vecXT.setLength(count*nInputs_);
      pdfhist->genSample(count,vecXT.getDVector(),
               VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->updateInputSection(count, nInputs_, NULL, 
               VecILowerBs_.getDVector(),VecIUpperBs_.getDVector(),
               vecXT.getDVector(),StrInpNames_.getStrings(),
               NULL,NULL,NULL,NULL);
      kk = 1;
      vecYT.setLength(count);
      vecYT.setLength(count);
      for (ii = 0; ii < count; ii++) vecYT[ii] = PSUADE_UNDEFINED;
      ioPtr->updateOutputSection(count,kk,vecYT.getDVector(),
                   vecST.getIVector(),StrOutNames_.getStrings());
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
      strcpy(dataFile, "psuade_pdfhist_checksample");
      ioPtr->writePsuadeFile(dataFile, 0);
      delete ioPtr;
      ioPtr = NULL;
      printf("A large sample for cross-checking PDFs is in %s.\n",
             dataFile);
      delete pdfhist;
      fp = fopen("psuade_pdfhist_sample", "r");
      if (fp == NULL)
      {
        printf("ERROR: psuade_pdfhist_sample file not found.\n");
        cmdStatus = 1;
        continue;
      }
      else
      {
        fscanf(fp, "%d", &kk);
        fclose(fp);
        printf("Final iteration: scenario size = %d\n", kk);
        printf("Histogram has been created in psuade_pdfhist_sample.\n");
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ genhistogram2 
    //**/ generate histogram from a psuade sample given a target size
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "genhistogram2"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command uses the loaded sample and generates\n");
        printf("a sample of scenarios each associated with a\n");
        printf("probability. The target number of probable scenarios\n");
        printf("is given by user. Internally, PSUADE iterates to find\n");
        printf("the number of bins to be used for each inputs.\n");
        printf("This command is useful for building scenarios for\n");
        printf("optimization under uncertainty.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      vecIT.setLength(nInputs_);
      sprintf(pString, "Approximate number of scenarios (2-%d): ",
              nSamples_);
      int nScen = getInt(2, nSamples_, pString);
      int nbins = 2;
      int nScenLast=-1;
      PDFHistogram *pdfhist=NULL;
      while (1)
      {
        for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = nbins;
        pdfhist = new PDFHistogram(nSamples_,nInputs_,
                         VecSamInputs_.getDVector(),vecIT.getIVector(),1);
        delete pdfhist;
        fp = fopen("psuade_pdfhist_sample", "r");
        if (fp == NULL)
        {
          printf("ERROR: psuade_pdfhist_sample file not found.\n");
          break;
        }
        else
        {
          fscanf(fp, "%d", &kk);
          fclose(fp);
          printf("Current iteration: nbins = %d, scenario size = %d\n",
                 nbins,kk);
          //**/ stop here, but check to see which one to use
          if (kk > nScen)
          {
            if (nScenLast == -1 || (PABS(nScen-kk) < PABS(kk-nScenLast)))
            {
              printf("Final iteration: nbins = %d, scenario size = %d\n",
                     nbins,kk);
              break;
            }
            nbins = nbins - 1;
            if (nbins >= 2)
            {
              for (ii = 0; ii < nInputs_; ii++) vecIT[ii] = nbins;
              pdfhist = new PDFHistogram(nSamples_,nInputs_,
                             VecSamInputs_.getDVector(),vecIT.getIVector(),1);
              delete pdfhist;
              printf("Final iteration: nbins = %d, scenario size = %d\n",
                     nbins,nScenLast);
              break; 
            }
          }
        }
        nbins++;
      }
      printf("Histogram has been created in psuade_pdfhist_sample.\n");
      printf("INFO: To estimate the goodness of this histogram, run \n");
      printf("  the genhistogram command with the nbins information\n");
      printf("  given above. Thereafter, use the large sample in\n");
      printf("  psuade_pdfhist_checksample with iplot2_pdf to compare\n");
      printf("  against the iplot2_pdf plots from the original sample.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ master 
    //**/ set master mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "master"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off master mode.\n");
        printf("This mode, when turned on, will trigger prompts\n");
        printf("for additional information in certain analysis.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  psConfig_.MasterModeOff();
      if (!strcmp(winput, "off")) psConfig_.MasterModeOn();
      if (!psConfig_.MasterModeIsOn())
      {
        psConfig_.MasterModeOn();
        printf("Master mode on\n");
      }
      else
      {
        psConfig_.MasterModeOff();
        printf("Master mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ gm 
    //**/ set gm mode
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gm"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on/off gm mode.\n");
        printf("This mode, when turned on, will trigger prompts\n");
        printf("for additional information in certain analysis.\n");
        continue;
      }
      if (!strcmp(winput, "on"))  psConfig_.GMModeOff();
      if (!strcmp(winput, "off")) psConfig_.GMModeOn();
      if (!psConfig_.GMModeIsOn())
      {
        psConfig_.GMModeOn();
        printf("gm mode on\n");
      }
      else
      {
        psConfig_.GMModeOff();
        printf("gm mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ useconfigfile 
    //**/ use configure file instead of interactive 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "useconfigfile"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command turns on the UseConfigFile mode.\n");
        printf("Syntax: useconfigfile <file>  or\n");
        printf("        useconfigfile off     to turn off this mode.\n");
        printf("This mode, when used, will cause some functions\n");
        printf("to read their parameters from the configure file \n");
        printf("instead of prompting users for parameters.\n");
        printf("NOTE: use genconfigfile to see what are available.\n");
        continue;
      }
      if (!strcmp(winput, "off"))
      {
        printf("Turn off config file usage.\n");
        psConfig_.reset();
      }
      fp = fopen(winput, "r");
      if (fp == NULL)
      {
        printf("ERROR: config file not found.\n");
        cmdStatus = 1;
        continue;
      }
      psConfig_.initialize(winput, 1);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ scilab 
    //**/ set scilab as plot tool
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "scilab"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("This command switches between matlab/scilab mode.\n");
        printf("When PSUADE creates output graphics, this mode will\n");
        printf("direct PSUADE to create matlab or scilab graphics.\n");
        printf("The default is matlab.\n");
        continue;
      }
      if (plotMatlab())
      {
        psConfig_.PlotTool_ = 1;
        printf("scilab mode on\n");
      }
       else
      {
        psConfig_.PlotTool_ = 0;
        printf("scilab mode off\n");
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsmax 
    //**/ set max data size 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsmax"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsmax: set maximum sample size for response surface\n"); 
        printf("       (used to guard against large sample sizes that\n");
        printf("       make response surface generation too expensive.)\n");
        printf("       The default max is 5000.\n");
        continue;
      }
      sscanf(lineIn, "%s %d", command, &kk);
      if (kk <= 1 || kk > 100000)
           psConfig_.RSMaxPts_ = 100000;
      else psConfig_.RSMaxPts_ = kk;
      printf("psuade max data size for response surface set to %d\n", 
             psConfig_.RSMaxPts_);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ rsca 
    //**/ RS-based brute force calibration 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "rsca"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("rsca: RS-based brute-force calibration (use mesh data)\n");
        printf("When the input dimension is low, an alternative to MCMC\n");
        printf("is to do a brute force search on an interpolated mesh.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ > 4)
      {
        printf("INFO: nInputs_ must be <= 4.\n");
        cmdStatus = 1;
        continue;
      }
      if (plotScilab())
      {
        printf("INFO: rsca is currently not available for scilab.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ ask users to specify the inputs and one output
      vecWT.setLength(nInputs_);
      iplot1 = 0;
      iplot2 = 1;
      iplot3 = 2;
      kk = nInputs_;
      if (nInputs_ > 3)
      {
        printf("Currently, you can only select at most 3 inputs.\n");
        sprintf(pString, "Enter the number of inputs to study (1 - 3) : ");
        kk = getInt(1, 3, pString);
      }
      if (kk == 1) iplot2 = iplot3 = -2;
      if (kk == 2) iplot3 = -2;
      if (nInputs_ > kk)
      {
        sprintf(pString, "Enter the first input (1 - %d) : ", nInputs_);
        iplot1 = getInt(1, nInputs_, pString);
        iplot1--;
        if (kk > 1)
        {
          iplot2 = iplot1;
          while (iplot1 == iplot2)
          {
            sprintf(pString, "Enter the second input (1 - %d), not %d : ",
                    nInputs_, iplot1+1);
            iplot2 = getInt(1, nInputs_, pString);
            iplot2--;
            if (iplot1 == iplot2)
              printf("ERROR: duplicate input number %d.\n",iplot2+1);
          }
        }
        if (kk > 2)
        {
          iplot3 = 3 - iplot1 - iplot2;
          while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
          {
            sprintf(pString,
              "Enter the input for t axis (1 - %d), not %d nor %d: ",
                    nInputs_, iplot1+1, iplot2+1);
            iplot3 = getInt(1, nInputs_, pString);
            iplot3--;
            if (iplot3 == iplot1 || iplot3 == iplot2)
              printf("ERROR: duplicate input number %d.\n",iplot3+1);
          }
        }
      }

      if      (kk == 1) nPtsPerDim = 1024;
      else if (kk == 2) nPtsPerDim = 64;
      else if (kk == 3) nPtsPerDim = 32;
      faFlag = 3;
      FuncApprox *faPtr = genFAInteractive(psuadeIO_, faFlag);
      if (faPtr == NULL) 
      {
        printf("Error detected.\n"); 
        cmdStatus = 1;
        continue;
      }
      faPtr->setOutputLevel(outputLevel_);
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(VecILowerBs_.getDVector(), 
                       VecIUpperBs_.getDVector());
      if (nInputs_ > kk)
        printf("The other inputs will be set to their nominal values.\n");
      for (iInd1 = 0; iInd1 < nInputs_; iInd1++)
      {
        if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
             vecWT[iInd1] = 0.5*(VecILowerBs_[iInd1]+VecIUpperBs_[iInd1]);
        else vecWT[iInd1] = 1.0;
      }
      jplot = 0;
      vecYT.setLength(nSamples_);
      for (sInd = 0; sInd < nSamples_; sInd++)
        vecYT[sInd] = VecSamOutputs_[sInd*nOutputs_+jplot];

      double *faXOut = NULL;
      double *faYOut = NULL;
      if (kk == 3)
      {
        printf("Please wait while generating data \n");
        fp = fopen("matlabrsbca3.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabrsbca3.m.\n");
          cmdStatus = 1;
          continue;
        }
        fprintf(fp, "twoPlots = 1;\n");
        fwritePlotCLF(fp);
        double GYmax = -1.0e35;
        double GYmin =  1.0e35;
        for (ii = 0; ii < nPtsPerDim; ii++)
        {
          printf(".");
          fflush(stdout);
          vecWT[iplot3] = 
              (VecIUpperBs_[iplot3] - VecILowerBs_[iplot3]) * ii / 
              (nPtsPerDim - 1.0) + VecILowerBs_[iplot3];
          faPtr->gen2DGridData(VecSamInputs_.getDVector(),vecYT.getDVector(), 
                iplot1, iplot2,vecWT.getDVector(), &faLeng, &faXOut,&faYOut);
          //**/ x and y data only needed to output once
          if (ii == 0)
          {
            fprintf(fp, "X = [\n");
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
              fprintf(fp, "%e ", faXOut[sInd*2]);
            fprintf(fp, "];\n");
            fprintf(fp, "Y(:,:) = [\n");
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
               fprintf(fp, "%e ", faXOut[sInd*2+1]);
            fprintf(fp, "];\n");
          }
          for (sInd = 0; sInd < faLeng; sInd++)
            faYOut[sInd] = exp(-faYOut[sInd]);

          //**/ output the response data data
          fprintf(fp, "A%d = [\n", ii + 1);
          for (sInd = 0; sInd < nPtsPerDim; sInd++)
          {
            for (jj = 0; jj < nPtsPerDim; jj++)
            {
              dtemp = faYOut[sInd*nPtsPerDim+jj];
              count = 1;
              if (sInd > 0)
              {
                dtemp += faYOut[(sInd-1)*nPtsPerDim+jj];
                count++;
              }
              if (sInd < nPtsPerDim-1)
              {
                dtemp += faYOut[(sInd+1)*nPtsPerDim+jj];
                count++;
              }
              if (jj > 0)
              {
                dtemp += faYOut[sInd*nPtsPerDim+jj-1];
                count++;
              }
              if (jj < nPtsPerDim-1)
              {
                dtemp += faYOut[sInd*nPtsPerDim+jj+1];
                count++;
              }
              if (sInd > 0 && jj > 0)
              {
                dtemp += faYOut[(sInd-1)*nPtsPerDim+jj-1];
                count++;
              }
              if (sInd > 0 && jj < nPtsPerDim-1)
              {
                dtemp += faYOut[(sInd-1)*nPtsPerDim+jj+1];
                count++;
              }
              if (sInd < (nPtsPerDim-1) && jj > 0)
              {
                dtemp += faYOut[(sInd+1)*nPtsPerDim+jj-1];
                count++;
              }
              if (sInd < (nPtsPerDim-1) && jj < (nPtsPerDim-1))
              {
                dtemp += faYOut[(sInd+1)*nPtsPerDim+jj+1];
                count++;
              }
              dtemp /= ((double) count);
              fprintf(fp, "%e\n", dtemp);
              if (dtemp > Ymax) Ymax = dtemp;
              if (dtemp < Ymin) Ymin = dtemp;
            }
          }
          fprintf(fp, "];\n");
          fprintf(fp, "A%d = reshape(A%d,%d,%d);\n", ii+1, ii+1,
                  nPtsPerDim, nPtsPerDim);
                                                                              
          //**/ search for lower and upper limits
          if (Ymax > GYmax) GYmax = Ymax;
          if (Ymin < GYmin) GYmin = Ymin;
          delete [] faXOut;
          delete [] faYOut;
        }

        //**/ create matlab movie
        fprintf(fp, "disp(\'Please wait while loading data.\')\n");
        fprintf(fp, "hold off\n");
        fwritePlotCLF(fp);
        for (ii = 0; ii < nPtsPerDim; ii++)
        {
          dtemp = (VecIUpperBs_[iplot3] - VecILowerBs_[iplot3]) *
                  ii / (nPtsPerDim - 1.0) + VecILowerBs_[iplot3];
          fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",
                  ii+1,nPtsPerDim);
          fprintf(fp, "if twoPlots == 1\n");
          fprintf(fp, "subplot(1,2,1), mesh(X,Y,A%d)\n", ii+1);
          fprintf(fp, "axis([%e %e %e %e %e %e])\n",VecILowerBs_[iplot1],
                  VecIUpperBs_[iplot1],VecILowerBs_[iplot2],
                  VecIUpperBs_[iplot2], GYmin, GYmax);
          fwritePlotXLabel(fp, StrInpNames_[iplot1]);
          fwritePlotYLabel(fp, StrInpNames_[iplot2]);
          fwritePlotZLabel(fp, StrOutNames_[jplot]);
          fwritePlotAxes(fp);
          fprintf(fp, "colorbar\n");
          fprintf(fp, "title(\'%s Mesh plot, val(3) = %14.7e\')\n",
                  StrOutNames_[jplot], dtemp);
          fprintf(fp, "subplot(1,2,2)\n");
          fprintf(fp, "end\n");
          fprintf(fp, "contourf(X,Y,A%d)\n",ii+1);
          fprintf(fp, "axis([%e %e %e %e])\n",VecILowerBs_[iplot1],
                  VecIUpperBs_[iplot1],VecILowerBs_[iplot2],
                  VecIUpperBs_[iplot2]);
          fwritePlotXLabel(fp, StrInpNames_[iplot1]);
          fwritePlotYLabel(fp, StrInpNames_[iplot2]);
          fwritePlotAxes(fp);
          fprintf(fp, "colorbar\n");
          fprintf(fp, "colormap(jet)\n");
          fprintf(fp, "title(\'%s contour plot, val(3) = %14.7e\')\n",
                  StrOutNames_[jplot], dtemp);
          fprintf(fp,"pause(1)\n");
        }
        fprintf(fp, "rotate3d on\n");
        fclose(fp);
        printf("\nmatlabrsbca3.m is now available.\n");
        printf("You can identify the max and min from the plots.\n");
        delete faPtr;
        faPtr = NULL;
      }

      if (kk == 2)
      {
        printf("Please wait while generating data ");
        fp = fopen("matlabrsbca2.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabrsbca2.m.\n");
          cmdStatus = 1;
          continue;
        }
        fprintf(fp, "twoPlots = 1;\n");
        fwritePlotCLF(fp);
        faPtr->gen2DGridData(VecSamInputs_.getDVector(),vecYT.getDVector(),
                    iplot1,iplot2,vecWT.getDVector(), &faLeng, &faXOut,&faYOut);
                                                                              
        for (jj = 0; jj < faLeng; jj++) faYOut[jj] = exp(-faYOut[jj]);
        Ymax = -1.0e35;
        Ymin =  1.0e35;
        fprintf(fp, "A = [\n");
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
        {
          for (jj = 0; jj < nPtsPerDim; jj++)
          {
            dtemp = faYOut[sInd*nPtsPerDim+jj];
            count = 1;
            if (sInd > 0)
            {
              dtemp += faYOut[(sInd-1)*nPtsPerDim+jj];
              count++;
            }
            if (sInd < nPtsPerDim-1)
            {
              dtemp += faYOut[(sInd+1)*nPtsPerDim+jj];
              count++;
            }
            if (jj > 0)
            {
              dtemp += faYOut[sInd*nPtsPerDim+jj-1];
              count++;
            }
            if (jj < nPtsPerDim-1)
            {
              dtemp += faYOut[sInd*nPtsPerDim+jj+1];
              count++;
            }
            if (sInd > 0 && jj > 0)
            {
              dtemp += faYOut[(sInd-1)*nPtsPerDim+jj-1];
              count++;
            }
            if (sInd > 0 && jj < nPtsPerDim-1)
            {
              dtemp += faYOut[(sInd-1)*nPtsPerDim+jj+1];
              count++;
            }
            if (sInd < (nPtsPerDim-1) && jj > 0)
            {
              dtemp += faYOut[(sInd+1)*nPtsPerDim+jj-1];
              count++;
            }
            if (sInd < (nPtsPerDim-1) && jj < (nPtsPerDim-1))
            {
              dtemp += faYOut[(sInd+1)*nPtsPerDim+jj+1];
              count++;
            }
            dtemp /= ((double) count);
            fprintf(fp, "%e\n", dtemp);
            if (dtemp > Ymax) Ymax = dtemp;
            if (dtemp < Ymin) Ymin = dtemp;
          }
        }
        fprintf(fp, "];\n");
        fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
          fprintf(fp, "%e\n", faXOut[sInd*2]);
        fprintf(fp, "];\n");
        fprintf(fp, "Y = [\n");
        for (sInd = 0; sInd < nPtsPerDim; sInd++)
          fprintf(fp, "%e\n", faXOut[sInd*2+1]);
        fprintf(fp, "];\n");
        fwritePlotCLF(fp);
        fprintf(fp, "if twoPlots == 1\n");
        fprintf(fp, "subplot(1,2,1), mesh(X,Y,A)\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotZLabel(fp, StrOutNames_[jplot]);
        fwritePlotAxes(fp);
        fprintf(fp, "title(\'Mesh Plot for %s\')\n",StrOutNames_[jplot]);
        fprintf(fp, "colorbar\n");
        fprintf(fp, "subplot(1,2,2)\n");
        fprintf(fp, "end\n");
        fprintf(fp, "contourf(X,Y,A)\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrInpNames_[iplot2]);
        fwritePlotAxes(fp);
        fprintf(fp, "colorbar\n");
        fprintf(fp, "colormap(jet)\n");
        fprintf(fp, "title(\'Contour Plot for %s\')\n",StrOutNames_[jplot]);
        fclose(fp);
        printf("matlabrsbca2.m is now available.\n");
        printf("You can identify the max and min from the plots.\n");
        delete [] faXOut;
        delete [] faYOut;
        delete faPtr;
        faPtr = NULL;
      }

      if (kk == 1)
      {
        printf("Please wait while generating data \n");
        fp = fopen("matlabrsbca1.m", "w");
        if (fp == NULL)
        {
          printf("ERROR: cannot open file matlabrsbca1.m.\n");
          cmdStatus = 1;
          continue;
        }
        fwritePlotCLF(fp);
        faPtr->gen1DGridData(VecSamInputs_.getDVector(),vecYT.getDVector(),
                   iplot1, vecWT.getDVector(), &faLeng, &faXOut,&faYOut);
        for (jj = 0; jj < faLeng; jj++)
          faYOut[jj] = exp(-faYOut[jj]);
        Ymax = -1.0e35;
        for (sInd = 0; sInd < faLeng; sInd++)
        {
          if (PABS(faYOut[sInd]) > Ymax) Ymax = PABS(faYOut[sInd]);
          if (PABS(faYOut[sInd]) < Ymin) Ymin = PABS(faYOut[sInd]);
        }
        if (Ymax == 0.0) Ymax = 1.0;
        printf("Ymax = %e\n", Ymax);
        for (jj = 0; jj < faLeng; jj++) faYOut[jj] /= Ymax;
        fprintf(fp, "A = [\n");
        vecVT.setLength(faLeng);
        sprintf(pString,"How many smoothing step : (0 - 100) : ");
        int nSmooth = getInt(0, 1000, pString);
        for (jj = 0; jj < nSmooth; jj++)
        {
          for (sInd = 0; sInd < faLeng; sInd++)
            vecVT[sInd] = faYOut[sInd];
   
          for (sInd = 0; sInd < faLeng; sInd++)
          {
            dtemp = vecVT[sInd];
            count = 1;
            if (sInd > 0)
            {
              dtemp += vecVT[sInd-1];
              count++;
            }
            if (sInd < nPtsPerDim-1)
            {
              dtemp += vecVT[sInd+1];
              count++;
            }
            dtemp /= ((double) count);
            faYOut[sInd] = dtemp;
          }
        }
        Ymax = -1.0e35;
        Ymin =  1.0e35;
        for (sInd = 0; sInd < faLeng; sInd++)
        {
          fprintf(fp, "%e\n", faYOut[sInd]);
          if (dtemp > Ymax) Ymax = dtemp;
          if (dtemp < Ymin) Ymin = dtemp;
        }
        fprintf(fp, "];\n");
        fprintf(fp, "X = [\n");
        for (sInd = 0; sInd < faLeng; sInd++)
          fprintf(fp, "%e\n", faXOut[sInd]);
        fprintf(fp, "];\n");
        fprintf(fp, "plot(X,A)\n");
        fwritePlotXLabel(fp, StrInpNames_[iplot1]);
        fwritePlotYLabel(fp, StrOutNames_[jplot]);
        fwritePlotAxes(fp);
        fprintf(fp, "title(\'Likelihood Plot for %s\')\n",
                StrOutNames_[jplot]);
        fclose(fp);
        printf("matlabrsbca1.m is now available.\n");
        printf("You can identify the max and min from the plots.\n");
        delete [] faXOut;
        delete [] faYOut;
        delete faPtr;
        faPtr = NULL;
        cmdStatus = 0;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ pdfcheck 
    //**/ 1-sample and 2-sample tests
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "pdfcheck"))
    {
      int nSam;
      psMatrix corMat;
      psVector vIn, vOut, vLBs, vUBs;

      vecIT.setLength(4);
      vecWT.setLength(4);
      vecVT.setLength(4);
      vLBs.setLength(4);
      vUBs.setLength(4);
      nSam = 200000;
      vecYT.setLength(nSam);

      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("pdfcheck: 1-sample and 2-sample tests for distributions.\n");
        printf("Syntax: pdfcheck (no argument needed)\n");
        printf("This is a demonstration command only. It creates and\n");
        printf("reconstruct different probability distributions.\n");
        continue;
      }
      //**/ 1-input test
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### One-input test (Normal(1,4)): \n");
      printf("Sample mean     should be = 1\n");
      printf("Sample std dev  should be = 2\n");
      printf("Sample skewness should be = 0\n");
      printf("Sample Kurtosis should be = 3\n");
      printDashes(PL_INFO, 0);
      vecWT[0] = 1.0;
      vecVT[0] = 2.0;
      vLBs[0] = -6.0;
      vUBs[0] = 8.0;
      PDFManager *pdfman = new PDFManager();
      vecIT[0] = PSUADE_PDF_NORMAL;
      ddata = 1.0;
      corMat.setDim(1,1);
      corMat.setEntry(0,0,ddata);
      pdfman->initialize(1, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(1, vLBs.getDVector());
      vecUBs.load(1, vUBs.getDVector());
      vIn.setLength(nSam); 
      vOut.setLength(nSam); 
      pdfman->genSample(nSam, vOut, vecLBs, vecUBs);
      int analysisMethod = PSUADE_ANA_MOMENT;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,vOut,
                          outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 1-input test
      pdfman = new PDFManager();
      vOut.setLength(nSam); 
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### One-input test (LogNormal(1,2)): \n");
      printf("Sample mean     should be = 4.5\n");
      printf("Sample std dev  should be = 5.8\n");
      printDashes(PL_INFO, 0);
      vecWT[0] = 1.0;
      vecVT[0] = 1.0;
      vLBs[0] = 0.0;
      vUBs[0] = 200.0;
      vecIT[0] = PSUADE_PDF_LOGNORMAL;
      ddata = 1.0;
      corMat.setDim(1,1);
      corMat.setEntry(0,0,ddata);
      pdfman->initialize(1, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(1, vLBs.getDVector());
      vecUBs.load(1, vUBs.getDVector());
      vOut.setLength(nSam); 
      pdfman->genSample(nSam, vOut, vecLBs, vecUBs);
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      vIn.setLength(nSam); 
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,
                          vIn,vOut,outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 1-input test
      pdfman = new PDFManager();
      vOut.setLength(nSam); 
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### One-input test (Weibull(1,1)): \n");
      printf("Sample mean     should be = 1\n");
      printf("Sample std dev  should be = 1\n");
      printf("Sample skewness should be = 2\n");
      printDashes(PL_INFO, 0);
      vecWT[0] = 1.0;
      vecVT[0] = 1.0;
      vLBs[0] = 0.0;
      vUBs[0] = 10.0;
      vecIT[0] = PSUADE_PDF_WEIBULL;
      ddata = 1.0;
      corMat.setDim(1,1);
      corMat.setEntry(0,0,ddata);
      pdfman->initialize(1, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(1, vLBs.getDVector());
      vecUBs.load(1, vUBs.getDVector());
      vOut.setLength(nSam); 
      pdfman->genSample(nSam, vOut, vecLBs, vecUBs);
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      vIn.setLength(nSam); 
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,
                          vOut,outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 2-input test
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### Two-input test (Normal(0,1)): \n");
      printf("Sample mean     should be = 0\n");
      printf("Sample std dev  should be = 1.414\n");
      printf("Sample skewness should be = 0\n");
      printf("Sample Kurtosis should be = 3\n");
      printDashes(PL_INFO, 0);
      pdfman = new PDFManager();
      vecIT[0] = PSUADE_PDF_NORMAL;
      vecIT[1] = PSUADE_PDF_NORMAL;
      vecWT[0] = vecWT[1] = 0.0;
      vecVT[0] = vecVT[1] = 1.0;
      vLBs[0] = vLBs[1] = -3.0;
      vUBs[0] = vUBs[1] =  3.0;
      ddata = 1.0;
      corMat.setDim(2,2);
      corMat.setEntry(0,0,ddata);
      corMat.setEntry(1,1,ddata);
      pdfman->initialize(2, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(2, vLBs.getDVector());
      vecUBs.load(2, vUBs.getDVector());
      vIn.setLength(2*nSam); 
      pdfman->genSample(nSam, vIn, vecLBs, vecUBs);
      vOut.setLength(nSam); 
      for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,vOut,
                          outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 2-input test
      pdfman = new PDFManager();
      vOut.setLength(nSam); 
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### Two-input test (LogNormal(log(mean)=0,1)): \n");
      printf("Sample mean     should be = 3.30\n");
      printf("Sample std dev  should be = 3.07\n");
      printDashes(PL_INFO, 0);
      vecWT[0] = vecWT[1] = 0.0;
      vecVT[0] = vecVT[1] = 1.0;
      vLBs[0] = vLBs[1] =  0.0;
      vUBs[0] = vUBs[1] =  200.0;
      vecIT[0] = PSUADE_PDF_LOGNORMAL;
      vecIT[1] = PSUADE_PDF_LOGNORMAL;
      ddata = 1.0;
      corMat.setDim(2,2);
      corMat.setEntry(0,0,ddata);
      corMat.setEntry(1,1,ddata);
      pdfman->initialize(2, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(2, vLBs.getDVector());
      vecUBs.load(2, vUBs.getDVector());
      vIn.setLength(2*nSam); 
      pdfman->genSample(nSam, vIn, vecLBs, vecUBs);
      vOut.setLength(nSam); 
      for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,vOut,
                          outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 2-input test with correlation
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### Two-input test (Normal(0,1)) with cor = 0.5: \n");
      printf("Sample mean     should be = 0\n");
      printf("Sample std dev  should be = 1.73\n");
      printf("Sample skewness should be = 0\n");
      printf("Sample Kurtosis should be = 3\n");
      printDashes(PL_INFO, 0);
      pdfman = new PDFManager();
      vecIT[0] = PSUADE_PDF_NORMAL;
      vecIT[1] = PSUADE_PDF_NORMAL;
      vecWT[0] = vecWT[1] = 0.0;
      vecVT[0] = vecVT[1] = 1.0;
      vLBs[0] = vLBs[1] = -5.0;
      vUBs[0] = vUBs[1] =  5.0;
      ddata = 1.0;
      corMat.setDim(2,2);
      corMat.setEntry(0,0,ddata);
      corMat.setEntry(1,1,ddata);
      ddata = 0.5;
      corMat.setEntry(0,1,ddata);
      corMat.setEntry(1,0,ddata);
      pdfman->initialize(2, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(2, vLBs.getDVector());
      vecUBs.load(2, vUBs.getDVector());
      vIn.setLength(2*nSam); 
      pdfman->genSample(nSam, vIn, vecLBs, vecUBs);
      vOut.setLength(nSam); 
      for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs, vIn, vOut,
                          outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 2-input test with correlation
      pdfman = new PDFManager();
      vOut.setLength(nSam); 
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("##### Two-input test (LogNormal(logmean=0),1) with cor=0.5\n");
      printf("Sample mean     should be = 3.3\n");
      printf("Sample std dev  should be = 3.5\n");
      printDashes(PL_INFO, 0);
      vecWT[0] = vecWT[1] = 0.0;
      vecVT[0] = vecVT[1] = 1.0;
      vLBs[0] = vLBs[1] = 0.0;
      vUBs[0] = vUBs[1] = 200.0;
      vecIT[0] = PSUADE_PDF_LOGNORMAL;
      vecIT[1] = PSUADE_PDF_LOGNORMAL;
      ddata = 1.0;
      corMat.setDim(2,2);
      corMat.setEntry(0,0,ddata);
      corMat.setEntry(1,1,ddata);
      ddata = 0.5;
      corMat.setEntry(0,1,ddata);
      corMat.setEntry(1,0,ddata);
      pdfman->initialize(2, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(2, vLBs.getDVector());
      vecUBs.load(2, vUBs.getDVector());
      vIn.setLength(2*nSam); 
      pdfman->genSample(nSam, vIn, vecLBs, vecUBs);
      vOut.setLength(nSam); 
      for (ii = 0; ii < nSam; ii++) vOut[ii] = vIn[ii*2] + vIn[ii*2+1];
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,vOut,
                          outputLevel_);
      delete anaManager;
      delete pdfman;

      //**/ 4-input test with correlation
      printDashes(PL_INFO, 0);
      printDashes(PL_INFO, 0);
      printf("4-input test (Normal(0,1), cor = 0.5, ");
      printf("LogNormal(0,1), no cor: \n");
      printf("Sample mean     should be = 3.3\n");
      printf("Sample std dev  should be = 3.5\n");
      printDashes(PL_INFO, 0);
      pdfman = new PDFManager();
      vecIT[0] = PSUADE_PDF_NORMAL;
      vecIT[1] = PSUADE_PDF_NORMAL;
      vecIT[2] = PSUADE_PDF_LOGNORMAL;
      vecIT[3] = PSUADE_PDF_LOGNORMAL;
      vecWT[0] = vecWT[1] = 0.0;
      vecWT[2] = vecWT[3] = 0.0;
      vecVT[0] = vecVT[1] = 1.0;
      vecVT[2] = vecVT[3] = 1.0;
      vLBs[0] = vLBs[1] = -10.0;
      vLBs[2] = vLBs[3] =  0.0;
      vUBs[0] = vUBs[1] =  10.0;
      vUBs[2] = vUBs[3] =  200.0;
      corMat.setDim(4,4);
      ddata = 1.0;
      corMat.setEntry(0,0,ddata);
      corMat.setEntry(1,1,ddata);
      corMat.setEntry(2,2,ddata);
      corMat.setEntry(3,3,ddata);
      ddata = 0.5;
      corMat.setEntry(0,1,ddata);
      corMat.setEntry(1,0,ddata);
      pdfman->initialize(4, vecIT.getIVector(), vecWT.getDVector(), 
                         vecVT.getDVector(), corMat, NULL, NULL);
      vecLBs.load(4, vLBs.getDVector());
      vecUBs.load(4, vUBs.getDVector());
      vIn.setLength(4*nSam); 
      pdfman->genSample(nSam, vIn, vecLBs, vecUBs);
      vOut.setLength(nSam); 
      for (ii = 0; ii < nSam; ii++)
        vOut[ii] = vIn[ii*4] + vIn[ii*4+1] + vIn[ii*4+2] + vIn[ii*4+3];
      analysisMethod = PSUADE_ANA_MOMENT;
      anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      anaManager->analyze(analysisMethod,nSam,vecLBs,vecUBs,vIn,
                          vOut,outputLevel_);
      delete anaManager;
      delete pdfman;
      cmdStatus = 0;
    }

    //**/ *************************************************************
    //**/ Commands for design of experiments
    //**/ *************************************************************
    else if (!strncmp(command, "odoe", 4))
    {
      cmdStatus = ODOEAnalysis(lineIn);
    }

    //**/ *************************************************************
    //**/ Commands for KPCA-based MCMC 
    //**/ *************************************************************
    //**/ -------------------------------------------------------------
    // +++ kernel PCA create
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kpca_create"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kpca_create: use given snapshot information to create\n");
        printf("            a KPCA model file.\n");
        printf("Syntax: kpca_create (no argument needed)\n");
        printf("This command asks for a collection of snapshots, the\n");
        printf("desired reduced dimension and kernel; and create a\n");
        printf("model file that can be used to convert a snapshot\n");
        printf("to its feature representation (use kpca_forward) or\n");
        printf("convert a feature vector back to its snapshot space\n");
        printf("(kpca_inverse).\n");
      }
      printf("Provide a file that contains snapshots (std format).\n");
      printf("   line 1: snapshot 1 data\n");
      printf("   line 2: snapshot 2 data\n");
      printf("   ...\n");
      sprintf(pString,"Enter your snapshot file name: ");
      getString(pString, dataFile);
      dataFile[strlen(dataFile)-1] = '\0';
      sprintf(pString,"How many snapshots are in the file? ");
      int nSnaps = getInt(1, 1000000000, pString);
      sprintf(pString,"Snapshot dimension (size of each snapshot)? ");
      int fieldSize = getInt(1, 1000000000, pString);
      psMatrix matSnapshots;
      matSnapshots.setDim(fieldSize,nSnaps);
      KPCA *kpcaObj = new KPCA();
      status = kpcaObj->readSnapshotsFile(dataFile, matSnapshots);
      if (status != 0) 
      {
        delete kpcaObj;
        continue;
      }
      if (matSnapshots.nrows() <= 1)
      {
        printf("ERROR: something is wrong (input dimension < 1)\n");
        delete kpcaObj;
        continue;
      }
      sprintf(pString, "Enter the desire reduced dimension: (1 - %d) ",
              matSnapshots.nrows());
      int rdim = getInt(1, matSnapshots.nrows(), pString);
      int kernel = kpcaObj->getKernel();
      sprintf(pString,"Enter your desired model file name: ");
      getString(pString, dataFile);
      dataFile[strlen(dataFile)-1] = '\0';
      status = kpcaObj->genKPCAModelFile(matSnapshots,kernel,rdim,
                                         dataFile);
      delete kpcaObj;
      if (status == 0) 
      {
        printf("Model file %s has been created.\n", winput);
        cmdStatus = 0;
      }
      else
      {
        printf("kpca_create unsuccessful.\n"); 
        cmdStatus = 1;
      }
    }
     
    //**/ -------------------------------------------------------------
    // +++ forward KPCA 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kpca_forward"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kpca_forward: use the already-created KPCA model to\n");
        printf("             transform a snapshot to its feature vector.\n");
        printf("Syntax: kpca_forward (no argument needed)\n");
        printf("This command asks for the already-created KPCA model,\n");
        printf("a snapshot, the output file where the feature vector is\n");
        printf("to be stored, and then perform the transformation.\n");
      }
      cmdStatus = KPCAForward(NULL,NULL,NULL);
    }
     
    //**/ -------------------------------------------------------------
    // +++ inverse KPCA 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kpca_inverse"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kpca_inverse: use the already-created KPCA model to\n");
        printf("        reconstruct a snapshot from a feature vector.\n");
        printf("Syntax: kpca_inverse (no argument needed)\n");
        printf("This command asks for the already-created KPCA model, ");
        printf("a feature\n");
        printf("vector, the output file where the reconstructed ");
        printf("snapshot is to be\n");
        printf("stored, and perform the transformation.\n");
      }
      cmdStatus = KPCAInverse(NULL,NULL,NULL,NULL);
    }
     
    //**/ -------------------------------------------------------------
    // +++ KPCA server
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "kpca_server"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("kpca_server: run kpca_server in a server-client mode\n");
        printf("Syntax: kpca_server (no argument needed)\n");
      }
      printf("This command waits for the presence of a file ");
      printf("called kpca_request.[1-4]\n");
      printf("This file should contain the names of 3 files:\n");
      printf("Line 1: name of the KPCA model file created by kpca_create\n");
      printf("Line 2: name of the feature vector file (low dimensional)\n");
      printf("Line 3: name of the result file to be written to.\n");
      printf("It then performs inverse KPCA on the incoming feature vector.\n");
      if (!strcmp(winput, "-h")) continue;
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin); 
      if (lineIn2[0] != 'y') continue; 

      int  kpcaFlag, reqID;
      char modelFile[1000], featureFile[1000], resultFile[1000];
      char simString[10000], oldModelFile[10000], kpcaReqFile[1000];;
      FILE *kpcaFp;
      KPCA *kpcaObj=NULL;
      strcpy(oldModelFile, "NOMODEL");
      printf("kpca_server has been launched.\n");
      printf("Create a file called 'kpca_quit' to quit this command.\n");
      while (1)
      {
        sprintf(kpcaReqFile, "kpca_request");
        fp = fopen(kpcaReqFile, "r");
        reqID = 0;
        if (fp == NULL) 
        {
          reqID = -1;
          for (ii = 1; ii < 9; ii++)
          {
            sprintf(kpcaReqFile, "kpca_request.%d", ii);
            fp = fopen(kpcaReqFile,"r");
            if (fp != NULL) 
            {
              reqID = ii;
              fclose(fp);
              fp = NULL;
              break;
            }
          }
        }
        if (reqID == -1) 
        {
#ifdef WINDOWS
          Sleep(1000);
#else
          sleep(1);
#endif
          if (kpcaFlag == 0) 
            printf("kpca_server: waiting for next kpca_request(.x) file\n");
          kpcaFlag++;
        }
        if (reqID >= 0) 
        {
          printAsterisks(PL_INFO, 0);
          printf("Request file found = %s\n",kpcaReqFile);
          fp = fopen(kpcaReqFile,"r");
          fscanf(fp,"%s",winput);
          if (!strcmp(winput,"inverse"))
          {
            status = 0;

            //**/ read model file name
            fscanf(fp,"%s",modelFile);
            modelFile[strlen(modelFile)] = '\0'; 
            printf("KPCA model file = %s\n", modelFile);
            kpcaFp = fopen(modelFile,"r");
            if (kpcaFp == NULL)
            {
              printf("kpca_server ERROR: model file %s not found\n",
                     modelFile);
              status = 1;
              fclose(fp);
              unlink(kpcaReqFile);
              printf("ABORT\n");
            }
            else fclose(kpcaFp);

            //**/ read feature file name
            if (status == 0)
            {
              fscanf(fp,"%s",featureFile);
              printf("KPCA feature file = %s\n", featureFile);
              featureFile[strlen(featureFile)] = '\0'; 
              kpcaFp = fopen(featureFile,"r");
              if (kpcaFp == NULL)
              {
                printf("kpca_server ERROR: feature file %s not found\n",
                       featureFile);
                status = 1;
                fclose(fp);
                unlink(kpcaReqFile);
                printf("ABORT\n");
              }
              else fclose(kpcaFp);
              if (status == 0 && (featureFile[strlen(featureFile)-1] != 
                  kpcaReqFile[strlen(kpcaReqFile)-1])) 
              {
                printf("Something wrong in kpca_server\n");
                printf("%s\n", kpcaReqFile);
                printf("%s\n", featureFile);
                printf("ABORT\n");
                fclose(fp);
                unlink(kpcaReqFile);
                status = 1;
              }
            }
            if (status == 0)
            {
              fscanf(fp,"%s",resultFile);
              printf("KPCA result file  = %s\n", resultFile);
              resultFile[strlen(resultFile)] = '\0'; 
              kpcaFp = fopen(resultFile,"w");
              if (kpcaFp == NULL)
              {
                printf("kpca_server ERROR: cannot write to file %s\n",
                       resultFile);
                fclose(fp);
                unlink(kpcaReqFile);
                printf("ABORT\n");
                status = 1;
              }
              else fclose(kpcaFp);
            }

            if (status == 0) 
            {
              fclose(fp);
              fp = fopen(resultFile, "r");
              if (fp != NULL)
              {
                fclose(fp);
                unlink(resultFile);
                fp = NULL;
              }
              if (!strcmp(modelFile, oldModelFile))
              {
                 printf("kpca_server INFO: re-use KPCA model\n");
                 KPCAInverse(modelFile,featureFile,resultFile, &kpcaObj);
              }
              else
              {
                 printf("INFO: new model -> re-create KPCA model\n");
                 strcpy(oldModelFile, modelFile);
                 if (kpcaObj != NULL) delete kpcaObj;
                 kpcaObj = NULL;
                 KPCAInverse(modelFile,featureFile,resultFile, &kpcaObj);
              }  
              unlink(kpcaReqFile);
              printAsterisks(PL_INFO, 0);
            }
          }
          else if (!strcmp(winput,"simulation"))
          {
            ii = getc(fp);
            fgets(simString, 10000, fp);
            printf("kpca_server: running simulation.\n");
            system(simString);
            printf("kpca_server: simulation completed.\n");
            unlink("kpca_request");
          }
          else
          {
            printf("kpca_server ERROR: request %s not recognized\n",
                   winput);
          }
          kpcaFlag = 0;
        }
        fp = fopen("kpca_quit","r");
        if (fp != NULL)
        {
          fclose(fp);
          unlink("kpca_quit");
          printf("File kpca_quit detected ==> quit\n");
          break; 
        }
      }
      if (kpcaObj != NULL) delete kpcaObj;
      kpcaObj = NULL;
    }
     
    //**/ -------------------------------------------------------------
    // +++ generate mcmc likelihood function 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gen_likelihood"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gen_likelihood: create a Python program to\n");
        printf("     compute likelihood for KPCA-based MCMC.\n");
        printf("     This function is used in mcmc_with_kpca.\n");
        printf("Syntax: gen_likelihood (no argument needed)\n");
      }
      //**/ explanation
      printf("Bayesian inference for high-dimensional correlated\n");
      printf("calibration parameters U with high-dimensional\n");
      printf("reduction can be expressed as: \n");
      printf("      posterior(U) ~ L(D|X) p(X|U) p(U)\n");
      printf("where\n");
      printf("      U - high-dimensional calibration parameters\n");
      printf("      X - KPCA-reduced feature parameters\n");
      printf("      D - measurement data\n");
      printf("   p(U) - prior distribution of U\n");
      printf(" p(X|U) - prior distribution of KPCA-reduced parameters\n\n");
      printf(" L(D|X) - likelihood of D given X\n\n");
      printf("This command helps create L(D|U) p(X|U) that can be\n");
      printf("inserted into mcmc_with_dr for DR-based inference.\n\n");
      printf("Please follow the instructions to create this likelihood\n");
      printf("function.\n");
      genKPCAMcmcWorkflow();
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ run MCMC with dimension reduction
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mcmc_with_dr"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mcmc_with_dr: perform inference on high-dimensional\n");
        printf("     inputs via dimension reduction methods such\n");
        printf("     as KPCA (with gen_likelihood, kpca_create, ..).\n");
        printf("Syntax: mcmc_with_dr (no argument needed)\n");
      }
      //**/ explanation
      printAsterisks(PL_INFO, 0);
      printf("This function solves: \n");
      printf("      posterior(U) ~ L(D|G(X)) p(X|U) p(U)\n");
      printf("where\n");
      printf("         U - high-dimensional calibration parameters\n");
      printf("         X - KPCA-reduced low-dimensional feature parameters\n");
      printf("         D - measurement data\n");
      printf("      p(U) - prior distribution of U (in ensemble snapshots)\n");
      printf("    p(X|U) - prior distribution of KPCA-reduced parameters\n");
      printf("      G(X) - Inverse operator from X to U space\n");
      printf(" L(D|G(X)) - likelihood of D given X\n\n");
      printf("p(U) is assumed to be a top-hat function.\n");
      printf("Thus, the only information needed from users are\n");
      printf("A likelihood function to compute L(D|X) p(X|U)\n");
      printf("- Use gen_likelihood to create this function\n");
      printf("An optional file that has the initial guess for U\n");
      printf("MCMC parameters such as maximum iterations\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput, 5000, stdin);
      if (lineIn2[0] != 'y') continue; 

      //**/ run command
      cmdStatus = MCMCWithKPCA();
    }

    //**/ *************************************************************
    //**/ advanced commands
    //**/ *************************************************************
    //**/ -------------------------------------------------------------
    // +++ gen_discrete 
    //**/ generate discrete value sample
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "gen_discrete"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("gen_discrete: generate a sample of discrete variables.\n");
        printf("Syntax: gen_discrete (no argument needed)\n");
        printf("You will be prompted for a parameter file describing\n");
        printf("the variables, their ranges, and their probabilities.\n");
        continue;
      }
      Sampling *sampPtr = 
              (Sampling *) SamplingCreateFromID(PSUADE_SAMP_DISCRETE);
      sampPtr->setPrintLevel(outputLevel_);
      sampPtr->setOutputParams(iOne);
      sampPtr->setSamplingParams(iOne, -1, 0);
      sampPtr->initialize(0);
      count = sampPtr->getNumSamples();
      kk    = sampPtr->getNumInputs();
      vecXT.setLength(count * kk);
      vecYT.setLength(count);
      vecST.setLength(count);
      sampPtr->getSamples(count,kk,iOne,vecXT.getDVector(),
                          vecYT.getDVector(),vecST.getIVector());
      delete sampPtr;
      sampPtr = NULL;
      vecVT.setLength(kk);
      vecWT.setLength(kk);
      for (ii = 0; ii < kk; ii++) vecVT[ii] = PSUADE_UNDEFINED;
      for (ii = 0; ii < kk; ii++) vecWT[ii] = -PSUADE_UNDEFINED;
      for (ii = 0; ii < kk; ii++)
      {
        for (jj = 0; jj < count; jj++)
        {
          if (vecXT[jj*kk+ii] < vecVT[ii]) vecVT[ii] = vecXT[jj*kk+ii];
          if (vecXT[jj*kk+ii] > vecWT[ii]) vecWT[ii] = vecXT[jj*kk+ii];
        }
      }
      psStrings Xnames;
      Xnames.setNumStrings(kk);
      for (ii = 0; ii < kk; ii++)
      {
        sprintf(pString, "X%d", ii+1);
        Xnames.loadOneString(ii, pString);
      }
      psStrings Yname;
      Yname.setNumStrings(iOne);
      strcpy(pString, "Y");
      Yname.loadOneString(0, pString);

      psIVector vecIPDFs;
      psVector  vecIMeans, vecIStds;
      vecIPDFs.setLength(nInputs_);
      vecIMeans.setLength(nInputs_); 
      vecIStds.setLength(nInputs_); 
      psMatrix *iCMat = new psMatrix();
      iCMat->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        vecIPDFs[ii] = 0;
        vecIMeans[ii] = 0;
        vecIStds[ii] = 0;
        iCMat->setEntry(ii,ii,1.0);
      }
      PsuadeData *ioPtr = new PsuadeData();
      ioPtr->updateInputSection(count,kk,NULL,vecVT.getDVector(),
          vecWT.getDVector(),vecXT.getDVector(),Xnames.getStrings(),
          vecIPDFs.getIVector(), vecIMeans.getDVector(), 
          vecIStds.getDVector(), iCMat); 
      delete iCMat;
      ioPtr->updateOutputSection(count, iOne, vecYT.getDVector(), 
                           vecST.getIVector(), Yname.getStrings()); 
      ioPtr->updateMethodSection(PSUADE_SAMP_MC, count, 1, -1, -1);
      ioPtr->writePsuadeFile("psuade_discrete_sample",0);
      printf("The test sample is in file psuade_discrete_sample.\n"); 
      delete ioPtr;
      ioPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ mo_opt
    //**/ response-surface based multi-objective optimization 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "mo_opt"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("mo_opt: multi-objective optimization\n");
        printf("Syntax: mo_opt\n");
        printf("NOTE: no need to load sample file first\n");
        continue;
      }
      printf("Step 1: MOO needs the name of a file in PSUADE ");
      printf("data format that has\n");
      printf("        its INPUT and OUTPUT sections containing ");
      printf("information about\n");
      printf("        model inputs and outputs.\n");
      printf("        Also, the driver and opt_driver fields in this ");
      printf("file should\n");
      printf("        point to the model simulator.\n");
      sprintf(pString,"Enter the name of this PSUADE input file: ");
      getString(pString, dataFile);
      kk = strlen(dataFile);
      dataFile[kk-1] = '\0';
      printf("PSUADE input file = %s\n", dataFile);
      fp = fopen(dataFile, "r");
      if (fp == NULL)
      {
        printf("ERROR: file %s not found.\n", dataFile);
        cmdStatus = 1;
        continue;
      }
      fclose(fp);
      cleanUp();
      PsuadeData *ioPtr = new PsuadeData();
      status = ioPtr->readPsuadeFile(dataFile);
      if (status != 0) 
      {
        printf("ERROR: cannot read file %s.\n", dataFile);
        continue;
      }

      //**/ get information from psuade input file
      ioPtr->getParameter("input_ninputs", pPtr);
      nInputs_ = pPtr.intData_;
      ioPtr->getParameter("output_noutputs", pPtr);
      nOutputs_ = pPtr.intData_;

      //**/ set MOO method in the ioPtr object
      int optMethod = 10; /* MOO */
      ioPtr->updateOptimizationSection(optMethod,1,1.0e-6,0,-1);

      //**/ set other MOO parameters and initial guess in matOptData
      ioPtr->getParameter("ana_opt_nlocalmin", pPtr);
      int nInitX=1;
      psMatrix matOptData;
      matOptData.setFormat(PS_MAT2D);
      matOptData.setDim(4, pPtr.intData_*nInputs_);
      ioPtr->getParameter("input_lbounds", pLower);
      ioPtr->getParameter("input_ubounds", pUpper);
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = 0.5*(pLower.dbleArray_[ii]+pUpper.dbleArray_[ii]);
        matOptData.setEntry(0,ii,ddata);
      }
      FunctionInterface *funcIO = createFunctionInterfaceSimplified(ioPtr);
      printf("Step 2: Follow more instructions inside MOO\n");
      status = OptimizerSearch(ioPtr, funcIO, matOptData, nInitX);

      //**/ clean up
      delete funcIO;
      delete ioPtr;
      pLower.clean();
      pUpper.clean();
      funcIO = NULL;
      ioPtr  = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ so_ua or soua  
    //**/ UQ for aleatoric and epistemic UQ using inner-outer iteration
    //**/ -------------------------------------------------------------
    else if ((!strcmp(command, "so_ua")) || (!strcmp(command, "soua")))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("soua: UQ for second order uncertainty analysis\n");
        printf("      (uncertainty in input probability distributions)\n");
        printf("Syntax: soua (no argument needed).\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs second-order uncertainty analysis ");
      printf("on the response\n");
      printf("surface constructed from the loaded sample, whereby ");
      printf("second-order means\n");
      printf("that the distributions of the uncertain parameters ");
      printf("are also uncertain.\n");
      printf("Users are prompted for uncertainty information of the ");
      printf("uncertain inputs\n");
      printf("and outer-inner iterations will be performed to ");
      printf("build ensemble CDFs\n");
      printf("similar to the p-box approach in 'aeua'.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin);
      if (lineIn2[0] != 'y') continue;

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        continue;
      }
      if (nInputs_ == 1)
      {
        printf("ERROR: nInputs must be 2 or more.\n");
        cmdStatus = 1;
        continue;
      }

      //**/ select output
      sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      //**/ set up the function approximator
      printf("Step 1: construct response surface\n");
      faFlag = 3;
      faType = -1;
      FuncApprox *faPtr=NULL;
      while (faType < 0)
      {
        faPtr = genFAInteractive(psuadeIO_, faFlag);
        if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
        faType = faPtr->getID();
      }
      if (faPtr == NULL) {printf("ERROR detected.\n"); continue;}
      nPtsPerDim = 64;
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);

      printf("Step 2: obtain uncertainty information of input distributions.\n");
      printf("   For example, if an input has normal distribution, ");
      printf("then uncertainties\n");
      printf("   (upper and lower bounds) of the mean and std dev are ");
      printf("to be specified\n");
      printf("   If no 2nd order uncertainty is needed for an input, ");
      printf("simply set the lower\n");
      printf("   and upper bounds to be the same.\n");
       
      int nParams = 2 * nInputs_;
      vecVT.setLength(nParams);
      vecWT.setLength(nParams);
      for (ii = 0; ii < nInputs_; ii++)
      {
        printf("Enter uncertainties for input %d:\n", ii+1); 
        if (VecInpPDFs_.length() == 0 || VecInpPDFs_[ii] == 0)
        {
          printf("Lower bounds = %e\n", VecILowerBs_[ii]);
          sprintf(pString, "Enter lower bound for input lower bound: ");
          vecVT[ii*2] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input lower bound: ");
          vecWT[ii*2] = getDouble(pString); 
          printf("Upper bounds = %e\n", VecIUpperBs_[ii]);
          sprintf(pString, "Enter lower bound for input upper bound: ");
          vecVT[ii*2+1] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input upper bound: ");
          vecWT[ii*2+1] = getDouble(pString); 
        }
        else if (VecInpPDFs_[ii] == PSUADE_PDF_NORMAL)
        {
          printf("Normal distribution mean = %e\n", VecInpMeans_[ii]);
          sprintf(pString, "Enter lower bound for input mean : ");
          vecVT[ii*2] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input mean : ");
          vecWT[ii*2] = getDouble(pString); 
          printf("Normal distribution std dev = %e\n", VecInpStds_[ii]);
          sprintf(pString, "Enter lower bound for input std dev : ");
          vecVT[ii*2+1] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input std dev : ");
          vecWT[ii*2+1] = getDouble(pString); 
        }
        else if (VecInpPDFs_[ii] == PSUADE_PDF_LOGNORMAL)
        {
          printf("LogNormal distribution mean = %e\n",VecInpMeans_[ii]);
          sprintf(pString, "Enter lower bound for input mean : ");
          vecVT[ii*2] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input mean : ");
          vecWT[ii*2] = getDouble(pString); 
          printf("LogNormal distribution std dev = %e\n",VecInpStds_[ii]);
          sprintf(pString, "Enter lower bound for input std dev : ");
          vecVT[ii*2+1] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input std dev : ");
          vecWT[ii*2+1] = getDouble(pString); 
        }
        else if (VecInpPDFs_[ii] == PSUADE_PDF_BETA)
        {
          printf("Beta distribution alpha = %e\n",VecInpMeans_[ii]);
          sprintf(pString, "Enter lower bound for input alpha : ");
          vecVT[ii*2] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input alpha : ");
          vecWT[ii*2] = getDouble(pString); 
          printf("Beta distribution beta = %e\n", VecInpStds_[ii]);
          sprintf(pString, "Enter lower bound for input beta : ");
          vecVT[ii*2+1] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input beta : ");
          vecWT[ii*2+1] = getDouble(pString); 
        }
        else if (VecInpPDFs_[ii] == PSUADE_PDF_TRIANGLE)
        {
          printf("Triangle distribution mean = %e\n",VecInpMeans_[ii]);
          sprintf(pString, "Enter lower bound for input mean : ");
          vecVT[ii*2] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input mean : ");
          vecWT[ii*2] = getDouble(pString); 
          printf("Triangle distribution half width = %e\n",
                 VecInpStds_[ii]);
          sprintf(pString, "Enter lower bound for input half width : ");
          vecVT[ii*2+1] = getDouble(pString); 
          sprintf(pString, "Enter upper bound for input half width : ");
          vecWT[ii*2+1] = getDouble(pString); 
        }
        else 
        {
          printf("ERROR: Other distributions currently not supported.\n");
          cmdStatus = 1;
          continue;
        }
      }

      //**/ generate a sample
      printf("Step 3: generate samples for inner and outer iterations.\n");
      int nValid=0;
      psIVector vecValidFlags;
      vecValidFlags.setLength(nInputs_*2);
      for (ii = 0; ii < nInputs_*2; ii++)
      {
        if (vecVT[ii] != vecWT[ii])
        {
          nValid++; 
          vecValidFlags[ii] = 1;
        }
        else
        {
          vecValidFlags[ii] = 0;
          for (jj = ii+1; jj < nInputs_*2; jj++)
          {
            vecVT[nValid+jj-ii-1] = vecVT[jj];
            vecWT[nValid+jj-ii-1] = vecWT[jj];
          }
        }
      }
      if (nValid == 0) 
      {
        printf("INFO: no perturbation has been prescribed.\n");
        cmdStatus = 1;
        continue;
      }
      int nSams=100, nSams2=1000;
      psIVector vecSamStas;
      psVector  vecSamInps, vecSamOuts;
      vecSamInps.setLength(nSams*nParams);
      vecSamOuts.setLength(nSams);
      vecSamStas.setLength(nSams);
      Sampling *sampPtr = 
            (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampPtr->setPrintLevel(PL_INTERACTIVE);
      sampPtr->setInputBounds(nValid, vecVT.getDVector(), 
                              vecWT.getDVector());
      sampPtr->setSamplingParams(nSams,-1,0);
      ii = 1;
      sampPtr->setOutputParams(ii);
      sampPtr->initialize(0);
      sampPtr->getSamples(nSams,nValid,ii,vecSamInps.getDVector(),
                 vecSamOuts.getDVector(), vecSamStas.getIVector());
      delete sampPtr;
      sampPtr = NULL;

      //**/ run the sample
      printf("Step 4: perform inner and outer iterations.\n");
      psIVector vecPDFTypes;
      psVector  vecMeans, vecStds, vecOneSample;
      psMatrix corMat;
      vecPDFTypes.setLength(nInputs_);
      vecMeans.setLength(nInputs_);
      vecStds.setLength(nInputs_);
      vecOneSample.setLength(nInputs_);

      for (ii = 0; ii < VecInpPDFs_.length(); ii++) 
        vecPDFTypes[ii] = VecInpPDFs_[ii];
      corMat.setDim(nInputs_,nInputs_);
      ddata = 1.0;
      for (ii = 0; ii < nInputs_; ii++) corMat.setEntry(ii,ii,ddata);
      if (plotScilab()) fp = fopen("scilabsoua.sci", "w");
      else              fp = fopen("matlabsoua.m", "w");
      fwritePlotCLF(fp);
      for (ss = 0; ss < nSams; ss++)
      { 
        kk = 0;
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (vecValidFlags[2*ii] == 1) 
          {
            vecMeans[ii] = vecSamInps[ss*nParams+kk];
            kk++;
          }
          else vecMeans[ii] = VecInpMeans_[ii];
          if (vecValidFlags[2*ii+1] == 1) 
          {
            vecStds[ii] = vecSamInps[ss*nParams+kk];
            kk++;
          }
          else vecMeans[ii] = VecInpStds_[ii];
        }
        PDFManager *pdfman = new PDFManager();
        pdfman->initialize(nInputs_,vecPDFTypes.getIVector(),
           vecMeans.getDVector(),vecStds.getDVector(),corMat,NULL,NULL);
        vecLBs.load(nInputs_, VecILowerBs_.getDVector());
        vecUBs.load(nInputs_, VecIUpperBs_.getDVector());
        vecInps.setLength(nSams2*nInputs_);
        pdfman->genSample(nSams2, vecInps, vecLBs, vecUBs);
        fprintf(fp, "Y = [\n");
        count = 0;
        for (kk = 0; kk < nSams2; kk++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            vecOneSample[ii] = vecInps[kk*nInputs_+ii];
          dtemp = faPtr->evaluatePoint(vecOneSample.getDVector());
          fprintf(fp, "%e\n", dtemp);
          count++;
        }
        fprintf(fp, "];\n");
        if (plotScilab()) fprintf(fp, "Y = gsort(Y,'g','i');\n");
        else              fprintf(fp, "Y = sort(Y);\n");
        fprintf(fp, "X = 1:%d;\n", count);
        fprintf(fp, "X = 0.001 * X;\n");
        fprintf(fp, "plot(Y,X)\n");
        sprintf(winput, "Cumulative Distributions");
        fwritePlotTitle(fp, winput);
        fwritePlotAxes(fp);
        if (StrOutNames_[outputID] != NULL) 
             sprintf(winput, "%s", StrOutNames_[outputID]);
        else sprintf(winput, "Output Values");
        fwritePlotXLabel(fp, winput);
        sprintf(winput, "Probabilities");
        fwritePlotYLabel(fp, winput);
        if (ss == 0)
        {
          if (plotScilab())
               fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
          else fprintf(fp, "hold on\n");
        }
        delete pdfman;
      }
      fclose(fp);
      printf("Plot file for 2nd order uncertainty analysis is now");
      if (plotScilab()) printf("in scilabsoua.sci.\n");
      else              printf("in matlabsoua.m.\n");
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ae_ua or aeua   
    //**/ UQ for aleatoric and epistemic UQ using inner-outer iteration
    //**/ -------------------------------------------------------------
    else if ((!strcmp(command, "ae_ua")) || (!strcmp(command, "aeua")))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("aeua: UQ for aleatoric-epistemic uncertainty analysis\n");
        printf("Syntax: aeua (no argument needed).\n");
        continue;
      }
      printAsterisks(PL_INFO, 0);
      printf("This command performs uncertainty analysis ");
      printf("on the response surface\n");
      printf("constructed from the loaded sample whereby ");
      printf("some of the sample inputs\n");
      printf("are aleatory (with prescribed distributions) ");
      printf("and some are epistemic\n");
      printf("(with prescribed ranges). These distributions ");
      printf("and ranges are taken\n");
      printf("from the loaded sample file. The result is a ");
      printf("bundle of CDFs each\n");
      printf("representing the output distribution based ");
      printf("on the aleatory inputs\n");
      printf("at some given epistemic input values. This bundle ");
      printf("of CDFs can be\n");
      printf("enveloped and the envelop is called a p-box.\n");
      printDashes(PL_INFO, 0);
      printf("Proceed ? (y or n to abort) ");
      scanf("%s", lineIn2);
      fgets(winput,5000,stdin);
      if (lineIn2[0] != 'y') continue;

      if (nInputs_ <= 0 || psuadeIO_ == NULL || nSamples_ <= 0)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      if (nInputs_ == 1)
      {
        printf("ERROR: nInputs must be 2 or more.\n");
        continue;
      }

      //**/ select output
      sprintf(pString,"Enter output number (1 - %d) : ",nOutputs_);
      outputID = getInt(1, nOutputs_, pString);
      outputID--;

      //**/ select epistemic uncertain parameters
      printf("Step 1: select aleatoric and epistemic parameters\n");
      int nEpistemic=0;
      psIVector vecUTypes;
      vecUTypes.setLength(nInputs_);
      //**/ set all parameters to be aleatoric
      for (ii = 0; ii < nInputs_; ii++) vecUTypes[ii] = 0;
      kk = 1;
      while (kk > 0)
      {
        sprintf(pString,
                "Select epistemic parameters (1 - %d, 0 if done) : ",
                nInputs_); 
        kk = getInt(0, nInputs_, pString);
        if (kk > 0)
        {
          vecUTypes[kk-1] = 1;
          nEpistemic++;
        }
      } 
      if (nEpistemic == 0 || nEpistemic == nInputs_)
      {
        printf("At least 1 and at most %d epistemic parameters are\n",
               nInputs_-1);
        printf("required for this command.\n");
        continue;
      }
      printf("You have specified %d epistemic parameters.\n",nEpistemic);

      //**/ set up the function approximator
      printf("Step 2: construct response surface\n");
      psuadeIO_->getParameter("ana_outputid", pPtr);
      int iSave = pPtr.intData_;
      psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,outputID,-1);
      faType = -1;
      FuncApprox *faPtr = genFA(faType,nInputs_,outputLevel_,nSamples_);
      if (faPtr == NULL) 
      {
        printf("ERROR : cannot construct response surface.\n"); 
        cmdStatus = 1;
        continue;
      }
      nPtsPerDim = 64;
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(VecILowerBs_.getDVector(),VecIUpperBs_.getDVector());
      faPtr->setOutputLevel(outputLevel_);
      vecYT.setLength(nSamples_);
      for (ii = 0; ii < nSamples_; ii++) 
        vecYT[ii] = VecSamOutputs_[ii*nOutputs_+outputID];
      status = faPtr->initialize(VecSamInputs_.getDVector(),vecYT.getDVector());

      //**/ generate a sample for outer iteration (uniform)
      printf("Step 3: construct CDFs via outer-inner iterations\n");
      printf("   Outer iteration: iterate on the epistemic samples\n");
      printf("   Inner iteration: iterate on the aleatoric samples\n");
      printf("   Objective: iterate 1000 times in the outer loop and ");
      printf("look for the\n");
      printf("              lower and upper edges of the p-box.\n");
      int nSams=20000;
      psIVector vecPDFs1;
      psVector  vecStds1, vecMeans1, vecSamInps, veclbs, vecubs;
      psMatrix  corMat1;

      vecPDFs1.setLength(nEpistemic);
      vecMeans1.setLength(nEpistemic);
      vecStds1.setLength(nEpistemic);
      veclbs.setLength(nEpistemic);
      vecubs.setLength(nEpistemic);
      ddata = 1.0;
      kk = 0;
      corMat1.setDim(nEpistemic,nEpistemic);
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecUTypes[ii] == 1)
        {
          corMat1.setEntry(kk,kk,ddata);
          vecPDFs1[kk] = VecInpPDFs_[ii]; 
          vecMeans1[kk] = VecInpMeans_[ii]; 
          vecStds1[kk] = VecInpStds_[ii]; 
          veclbs[kk] = VecILowerBs_[ii];
          vecubs[kk] = VecIUpperBs_[ii];
          kk++;
        }
      }
      PDFManager *pdfman = new PDFManager();
      pdfman->initialize(nEpistemic,vecPDFs1.getIVector(),
              vecMeans1.getDVector(),vecStds1.getDVector(),corMat1,
              SamPDFFiles_,VecSamPDFIndices_.getIVector());
      vecLBs.load(nEpistemic, veclbs.getDVector());
      vecUBs.load(nEpistemic, vecubs.getDVector());
      vecSamInps.setLength(nSams*nEpistemic);
      pdfman->genSample(nSams, vecSamInps, vecLBs, vecUBs);
      delete pdfman;

      //**/ generate a sample for the inner iterations
      int nAleatoric = nInputs_ - nEpistemic, nSams2=2000;
      psIVector vecPDFs2;
      psVector  vecStds2, vecMeans2;
      psMatrix  corMat2;

      vecPDFs2.setLength(nAleatoric);
      vecMeans2.setLength(nAleatoric);
      vecStds2.setLength(nAleatoric);
      veclbs.setLength(nAleatoric);
      vecubs.setLength(nAleatoric);
      ddata = 1.0;
      kk = 0;
      corMat2.setDim(nAleatoric,nAleatoric);
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecUTypes[ii] != 1)
        {
          corMat2.setEntry(kk,kk,ddata);
          vecPDFs2[kk] = VecInpPDFs_[ii]; 
          vecMeans2[kk] = VecInpMeans_[ii]; 
          vecStds2[kk] = VecInpStds_[ii]; 
          veclbs[kk] = VecILowerBs_[ii];
          vecubs[kk] = VecIUpperBs_[ii];
          kk++;
        }
      }
      pdfman = new PDFManager();
      pdfman->initialize(nAleatoric,vecPDFs2.getIVector(),
              vecMeans2.getDVector(),vecStds2.getDVector(),corMat2,
              SamPDFFiles_,VecSamPDFIndices_.getIVector());
      vecLBs.load(nAleatoric, veclbs.getDVector());
      vecUBs.load(nAleatoric, vecubs.getDVector());
      vecInps.setLength(nSams2*nAleatoric);
      pdfman->genSample(nSams2, vecInps, vecLBs, vecUBs);
      delete pdfman;

      if (plotScilab()) fp = fopen("scilabaeua.sci", "w");
      else              fp = fopen("matlabaeua.m", "w");
      if (fp == NULL)
      {
        printf("aeua ERROR: cannot open plot file.\n");
        cmdStatus = 1;
        continue;
      }
      fwritePlotCLF(fp);

      //**/ outer iteration
      psVector vecYmaxs, vecYmins, vec1Sample;
      vec1Sample.setLength(nInputs_);
      vecYmaxs.setLength(nSams2);
      vecYmins.setLength(nSams2);
      for (ss = 0; ss < nSams2; ss++)
      {
        vecYmaxs[ss] = -PSUADE_UNDEFINED;
        vecYmins[ss] =  PSUADE_UNDEFINED;
      }
      int converged = 0;
      int upperCnt, lowerCnt, convergedCnt=0;
      double upperAcc, lowerAcc;
      ss = 0;
      while (ss < nSams && (converged == 0 || ss < 1000))
      {
        if (outputLevel_ > 2) printf("Epistemic sample %d\n",ss+1);
        if (ss < 50) fprintf(fp, "Y%d = [\n",ss+1);
        if (ss > 0 && (ss % 100 == 0))
        {
          printf("%5.1f%% ", 0.1*ss);
          fflush(stdout);
        }
        int count2 = 0;
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (vecUTypes[ii] == 1)
          {
            vec1Sample[ii] = vecSamInps[ss*nEpistemic+count2];
            count2++;
          }
        }
        count = lowerCnt = upperCnt = 0;
        upperAcc = lowerAcc = 0;
        for (kk = 0; kk < nSams2; kk++)
        {
          int flag = 0;
          for (ii = 0; ii < nAleatoric; ii++)
          {
            if (vecInps[kk*nAleatoric+ii] < veclbs[ii] ||
              vecInps[kk*nAleatoric+ii] > vecubs[ii]) flag++;
          }
          if (flag == 0)
          {
            count2 = 0;
            for (ii = 0; ii < nInputs_; ii++)
            {
              if (vecUTypes[ii] == 0)
              {
                vec1Sample[ii] = vecInps[kk*nAleatoric+count2];
                count2++;
              }
            }
            dtemp = faPtr->evaluatePoint(vec1Sample.getDVector());
            if (ss < 50) fprintf(fp, "%e\n", dtemp);
            if (dtemp > vecYmaxs[kk]) 
            {
              if (vecYmaxs[kk] != 0) 
                upperAcc += PABS((vecYmaxs[kk]-dtemp)/vecYmaxs[kk]);
              vecYmaxs[kk] = dtemp;
              upperCnt++;
            }
            if (dtemp < vecYmins[kk])
            {
              if (vecYmins[kk] != 0) 
                lowerAcc += PABS((vecYmins[kk]-dtemp)/vecYmins[kk]);
              vecYmins[kk] = dtemp;
              lowerCnt++;
            }
            count++;
          }
        }
        ddata = 100.0 * upperAcc / nSams;
        if (outputLevel_ > 2 && ddata > 0.0 && ss > 0) 
          printf("  Upper lifted  %7.3f %% on average from last time\n",
                 ddata);
        ddata = 100.0 * lowerAcc / nSams;
        if (outputLevel_ > 2 && ddata > 0.0 && ss > 0) 
          printf("  Lower dropped %7.3f %% on average from last time\n",
                 ddata);
        if (upperCnt+lowerCnt == 0) convergedCnt++;
        else                        convergedCnt = converged = 0;
        if (outputLevel_ > 3) 
          printf("  Convergence indicator = %5d (20 times => converged)\n",
                 upperCnt+lowerCnt);
        if (convergedCnt >= 20) converged = 1;
        if (count < 0.5 * nSams2)
        {
          printf("WARNING: < half of the points are within bounds.\n");
          printf("         Input ranges may need to be widened.\n");
        }
        if (count == 0)
        {
          printf("ERROR: none of the sample points are within bounds.\n");
          continue;
        }
        if (ss < 50) fprintf(fp, "];\n");
        if (ss < 50) 
        {
          if (plotScilab()) 
               fprintf(fp,"Y%d = gsort(Y%d,'g','i');\n",ss+1,ss+1);
          else fprintf(fp,"Y%d = sort(Y%d);\n",ss+1,ss+1);
          fprintf(fp, "X = 1:%d;\n",count);
          fprintf(fp, "X = X' / %d;\n", count);
          fprintf(fp, "plot(Y%d,X)\n",ss+1);
          fprintf(fp, "drawnow\n");
          if (ss == 0)
          {
            sprintf(winput, "Cumulative Distributions");
            fwritePlotTitle(fp, winput);
            fwritePlotAxes(fp);
            if (StrOutNames_[outputID] != NULL) 
                 sprintf(winput,"%s",StrOutNames_[outputID]);
            else sprintf(winput,"Output Values");
            fwritePlotXLabel(fp, winput);
            sprintf(winput, "Probabilities");
            fwritePlotYLabel(fp, winput);
            if (plotScilab())
                 fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            else fprintf(fp, "hold on\n");
          }
        }
        ss++;
      }
      printf("\n");
      count = 0;
      fprintf(fp, "YU = [\n");
      for (ss = 0; ss < nSams2; ss++)
      {
        if (vecYmaxs[ss] > -PSUADE_UNDEFINED) 
        {
          fprintf(fp, "%e\n", vecYmaxs[ss]);
          count++;
        }
      }
      fprintf(fp, "];\n");
      if (plotScilab()) 
           fprintf(fp,"YU = gsort(YU,'g','i');\n");
      else fprintf(fp,"YU = sort(YU);\n");
      fprintf(fp, "X = 1:%d;\n",count);
      fprintf(fp, "X = X' / %d;\n", count);
      fprintf(fp, "plot(YU,X,'r-','lineWidth',3)\n");

      count = 0;
      fprintf(fp, "YL = [\n");
      for (ss = 0; ss < nSams2; ss++)
      {
        if (vecYmins[ss] < PSUADE_UNDEFINED) 
        {
          fprintf(fp, "%e\n", vecYmins[ss]);
          count++;
        }
      }
      fprintf(fp, "];\n");
      if (plotScilab()) 
           fprintf(fp,"YL = gsort(YL,'g','i');\n");
      else fprintf(fp,"YL = sort(YL);\n");
      fprintf(fp, "X = 1:%d;\n",count);
      fprintf(fp, "X = X' / %d;\n", count);
      fprintf(fp, "plot(YL,X,'r-','lineWidth',3)\n");
      fclose(fp);
      printf("Plot file for aleatoric-epistemic analysis is now ");
      if (plotScilab()) printf("in scilabaeua.sci.\n");
      else              printf("in matlabaeua.m.\n");
      psuadeIO_->updateAnalysisSection(-1,-1,-1,-1,iSave,-1);
      delete faPtr;
      faPtr = NULL;
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ sys
    //**/ system command
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sys"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sys: execute a system command \n");
        printf("Syntax: sys <command>\n");
        continue;
      }
      count = 0;
      for (ii = 0; ii < 498; ii++)
        if (lineIn[ii] == 's' && lineIn[ii+1] == 'y' &&
          lineIn[ii+2] == 's') break;
      if (ii != 498) system(&lineIn[ii+4]);
    }

    //**/ -------------------------------------------------------------
    // +++ ivec_create
    //**/ create random input vector
    //**/ -------------------------------------------------------------
    else if (!strcmp(command,"ivec_create") || !strcmp(command,"vcreate"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ivec_create: create an internal vector (for rseval)\n");
        printf("Syntax: ivec_create or vcreate (no argument needed)\n");
        printf("This command is used together with rscreate and rseval\n");
        printf("to create a response surface on the fly and evaluate\n");
        printf("a new sample point placed into the local register.\n");
        continue;
      }
      if (nInputs_ <= 0 || psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      VecDataReg_.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++) 
        VecDataReg_[ii] = 0.5 * (VecIUpperBs_[ii] + VecILowerBs_[ii]); 
      printf("Internal vector created and the values have been set");
      printf(" to be the mid points.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ivec_modify
    //**/ modify an entry in the input vector
    //**/ -------------------------------------------------------------
    else if (!strcmp(command,"ivec_modify") || !strcmp(command,"vmodify"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ivec_modify: modify an entry in the internal vector\n");
        printf("Syntax: ivec_modify or vmodify (no argument needed)\n");
        printf("This command is used after ivec_create to modify\n");
        printf("individual entries in the internal register.\n");
        continue;
      }
      if (VecDataReg_.length() == 0)
      {
        printf("ERROR: use ivec_create before this command.\n");
        cmdStatus = 1;
        continue;
      }
      for (ii = 0; ii < nInputs_; ii++)
      {
        sprintf(pString,"New value for input %d? (%e - %e) ",ii+1,
                VecILowerBs_[ii], VecIUpperBs_[ii]);
        ddata = getDouble(pString);
        if (ddata < VecILowerBs_[ii] || ddata > VecIUpperBs_[ii])
          printf("WARNING: data out of range (extrapolation).\n");
        VecDataReg_[ii] = ddata;
      }
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ ivec_show
    //**/ show the input vector
    //**/ -------------------------------------------------------------
    else if (!strcmp(command,"ivec_show") || !strcmp(command, "vshow"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("ivec_show: display the internal vector\n");
        printf("Syntax: ivec_show or vshow (no argument needed)\n");
        printf("This command is used together with ivec_create and\n");
        printf("ivec_modify.\n");
        continue;
      }
      if (VecDataReg_.length() == 0)
      {
        printf("ERROR: use ivec_create before this command.\n");
        cmdStatus = 1;
        continue;
      }
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %3d = %e\n", ii+1, VecDataReg_[ii]);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ showformat
    //**/ show auxillary file formats
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "showformat"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("showformat: display different PSUADE file formats \n");
        printf("Syntax: showformat (no argument needed)\n");
        printf("User-provided information are often needed in design\n");
        printf("and analysis. These information are provided by users\n");
        printf("to PSUADE at various stages. This command lists many\n");
        printf("such file formats.\n");
        continue;
      }
      printDashes(PL_INFO, 0);
      printf("PSUADE data file format (can be read by load): \n");
      printf("  run 'geninputfile' command to create an example.\n");
      printDashes(PL_INFO, 0);
      printf("Standard file format (can be read by read_std): \n");
      printf("  line 1: nSamples nInputs nOutputs\n");
      printf("  line 2: <sample 1 inputs> < sample 1 outputs>\n");
      printf("  line 3: <sample 2 inputs> < sample 2 outputs>\n");
      printf("  .......\n");
      printDashes(PL_INFO, 0);
      printf("Xls file format (can be read by read_xls): \n");
      printf("  line 1: nSamples nInputs nOutputs\n");
      printf("  line 2: (optional) # sample input and output names\n");
      printf("  line 3: 1 <sample 1 inputs> < sample 1 outputs>\n");
      printf("  line 4: 2 <sample 2 inputs> < sample 2 outputs>\n");
      printf("  .......\n");
      printDashes(PL_INFO, 0);
      printf("MCMC experimental data file format (O1 = output 1): \n");
      printf("  line 1: PSUADE_BEGIN\n");
      printf("  line 2: nExps(p) nOutputs(n) nDesignInps designInpList\n");
      printf("  line 3: 1 <designInp ...> <O1 mean> <O1 std dev> ... \n");
      printf("  line 4: 2 <designInp ...> <O1 mean> <O1 std dev> ... \n");
      printf("  ...\n");
      printf("  line  : p <designInp ...> <O1 mean> <O1 std dev> ... \n");
      printf("  line  : PSUADE_END\n");
      printDashes(PL_INFO, 0);
      printf("Sample Input Only format (used in PDF S type,iread):\n");
      printf("  line 1: PSUADE_BEGIN\n");
      printf("  line 2: <number of sample points> <number of inputs>\n");
      printf("  line 3: (optional) : '#' followed by input names\n");
      printf("  line 4: 1 sample point 1 inputs \n");
      printf("  line 5: 2 sample point 2 inputs \n");
      printf("  line 6: 3 sample point 3 inputs \n");
      printf("  ...\n");
      printf("  line n: PSUADE_END\n");
      printDashes(PL_INFO, 0);
      printf("RSConstraints file format (used in Analysis: rs_constraint):\n");
      printf("  line 1: nInputs\n ");
      printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n ");
      printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n ");
      printf("  ... \n");
      printDashes(PL_INFO, 0);
      printf("RS index file format (used in Analysis: rs_index_file): \n");
      printf("  line 1: nInputs in rs data (driver) file\n");
      printf("  line 2: 1 <num> <default if num == 0>\n");
      printf("  line 3: 2 <num> <0 if num != 0>\n");
      printf("  line 4: 3 <num> <default if num == 0>\n");
      printf("  line 5: 4 <num> <0 if num != 0>\n");
      printf("  ...\n");
      printDashes(PL_INFO, 0);
      printf("MOATConstraints file format ");
      printf("(used in Analysis: moat_constraint):\n");
      printf("  line 1: nInputs \n");
      printf("  line 2: <input (or 0)> <value (nominal val if 0)> \n");
      printf("  line 3: <input (or 0)> <value (nominal val if 0)> \n");
      printf("  ... \n");
      printDashes(PL_INFO, 0);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ setdriver 
    //**/ set the driver field in the data file
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "setdriver"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("setdriver: set the driver field in the PSUADE file\n");
        printf("Syntax: setdriver\n");
        continue;
      }
      if (psuadeIO_ == NULL)
      {
        printf("ERROR: sample data not loaded yet.\n");
        cmdStatus = 1;
        continue;
      }
      printf("Which driver do you want to set: \n");
      printf("1. simulation driver (driver) \n");
      printf("2. optimization driver (opt_driver) \n");
      printf("3. auxiliary optimization driver (aux_opt_driver) \n");
      printf("4. ensemble simulation driver (ensemble_driver) \n");
      printf("5. ensemble optimization driver (ensemble_opt_driver) \n");
      sprintf(pString, "Enter (1 - 5): ");
      kk = getInt(1,5,pString);
      printf("Enter name of the driver : ");
      fgets(dataFile, 5000, stdin);
      dataFile[strlen(dataFile)-1] = '\0';
      for (ii = 0; ii < strlen(dataFile)-1; ii++)
        if (dataFile[ii] == ' ') break;
      if (ii == strlen(dataFile)-1 && (fp=fopen(dataFile,"r")) == NULL)
      {
        printf("WARNING: file %s not found.\n", dataFile);
      }
      if (kk == 1)
        psuadeIO_->updateApplicationSection(dataFile,NULL,NULL,NULL,
                                            NULL,-1);
      else if (kk == 2)
        psuadeIO_->updateApplicationSection(NULL,dataFile,NULL,NULL,
                                            NULL,-1);
      else if (kk == 3)
        psuadeIO_->updateApplicationSection(NULL,NULL,dataFile,NULL,
                                            NULL,-1);
      else if (kk == 4)
        psuadeIO_->updateApplicationSection(NULL,NULL,NULL,dataFile,
                                            NULL,-1);
      else if (kk == 5)
        psuadeIO_->updateApplicationSection(NULL,NULL,NULL,NULL,
                                            dataFile,-1);
      printf("Use 'write' to update your PSUADE input file.\n");
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ start_matlab 
    //**/ start matlab
    //**/ -------------------------------------------------------------
    else if ((!strcmp(command, "start_matlab")))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("start_matlab: start matlab from within PSUADE\n");
        printf("              interactive mode.\n");
        printf("Syntax: start_matlab (no argument needed)\n");
        continue;
      }
      strcpy(command, "which matlab");
      status = system(command);
      if (status != 0)
        printf("Matlab not found (have you set the path?)\n");
      else
      {   
        //**/ Jan 2020 - it seems like -nodesktop does not work
        //**/printf("Two modes to start Matlab: \n");
        //**/printf("1. desktop mode (create a matlab window)\n");
        //**/printf("2. nodesktop mode (use current window for matlab)\n");
        //**/sprintf(pString,"Start Matlab in desktop mode (y or n)? ");
        //**/getString(pString, winput);
        winput[0] = 'y';
        if (winput[0] == 'y') strcpy(command, "matlab");
        else                  strcpy(command, "matlab -nodesktop");
        status = system(command);
      }
    }

    //**/ -------------------------------------------------------------
    // +++ checkformat 
    //**/ check PSUADE file formats
    //**/ -------------------------------------------------------------
    else if ((!strcmp(command, "checkformat")))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("checkformat: check format of files used by PSUADE.\n");
        printf("Syntax: checkformat (no argument needed)\n");
        continue;
      }
      printf("Please select from the following file formats to check: \n");
      printf("1. OUU Z3 sample file.\n");
      printf("2. OUU Z4 sample file.\n");
      printf("3. psuade input file\n");
      printf("4. MCMC specification file\n");
      printf("5. MCMC posterior file\n");
      printf("6. S-type PDF sample file\n");
      sprintf(pString, "Enter your selection : ");
      int option = getInt(1,6,pString);
      printf("Please enter your sample file: ");
      scanf("%s", cString);
      fgets(winput, 1000, stdin);
      if (option == 1 || option == 2)
      {
        sprintf(pString, "How many inputs are in this sample file? ");
        kk = getInt(1,10000000,pString);
        status = checkOUUFileFormat(cString, option, kk, outputLevel_);
      }
      else if (option == 3)
      {
        PsuadeData *ioPtr = new PsuadeData();
        status = ioPtr->readPsuadeFile(cString);
        if (status != 0) 
        {
          printf("PSUADE data file format is not valid.\n");
          status = -1;
        }
        ioPtr->getParameter("input_ninputs", pPtr);
        if (pPtr.intData_ <= 0) 
        {
          printf("PSUADE data file has no inputs.\n");
          status = -1;
        }
        ioPtr->getParameter("input_sample", pPtr);
        if (pPtr.dbleArray_ == NULL) 
          printf("PSUADE data file has no sample.\n");
        else delete [] pPtr.dbleArray_;
        pPtr.dbleArray_ = NULL;
      }
      else if (option == 4)
      {
        status = checkMCMCFileFormat(cString, 0, outputLevel_);
      }
      else if (option == 5)
      {
        status = checkMCMCFileFormat(cString, 1, outputLevel_);
      }
      else if (option == 6)
      {
        status = checkSPDFFileFormat(cString, outputLevel_);
      }
      if (status == 0) 
      {
        printf("PASSED: file format validated.\n");
        cmdStatus = 0;
      }
      else
      {
        printf("FAILED: invalid file format.\n");
        cmdStatus = 1;
      }
    }

    //**/ -------------------------------------------------------------
    // +++ sample_info or sinfo
    //**/ outputs current sample in memory
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "sample_info") || !strcmp(command, "sinfo"))
    {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
        printf("sample_info (or sinfo): display sample information.\n");
        printf("Syntax: sample_info or sinfo (no argument needed)\n");
        continue;
      }
      printf("Sample in memory: nSamples = %d\n", nSamples_);
      printf("                  nInputs  = %d\n", nInputs_);
      printf("                  nOutputs = %d\n", nOutputs_);
      printf("Input names: \n");
      for (ii = 0; ii < nInputs_; ii++)
        if (StrInpNames_[ii] != NULL)
          printf("  Input %4d %s = %12.5e %12.5e\n",ii+1,
                 StrInpNames_[ii], VecILowerBs_[ii], VecIUpperBs_[ii]);
      printf("Output names: \n");
      for (ii = 0; ii < nOutputs_; ii++)
        if (StrOutNames_[ii] != NULL)
          printf("  Output %4d %s\n",ii+1,StrOutNames_[ii]);
      count = 0;
      for (ss = 0; ss < nSamples_; ss++) 
        if (VecSamStates_[ss] == 1) count++;
      printf("Number of valid sample points (state = 1)  = %d\n",count);
      count = 0;
      for (ss = 0; ss < nSamples_; ss++)
      {
        kk = 1;
        for (ii = 0; ii < nOutputs_; ii++)
          if (VecSamOutputs_[ss*nOutputs_+ii] == PSUADE_UNDEFINED) 
            kk = 0;
        if (kk == 1) count++;
      }
      printf("Number of sample points with valid outputs = %d\n",count);
      cmdStatus = 0;
    }

    //**/ -------------------------------------------------------------
    // +++ get status 
    //**/ -------------------------------------------------------------
    else if (!strcmp(command, "getstatus"))
    {
      if (!strcmp(winput, "-h"))
      {
        printf("getstatus: get status of previous command.\n");
        printf("Syntax: getstatus (no argument needed)\n");
        printf("NOTE: if an error occurs previously, status = 1\n");
        continue;
      }
      if (cmdStatus > 0) 
        printf("An error occurs in the last command.\n");
      else
        printf("No error occurs in the last command.\n");
    }

    //**/ -------------------------------------------------------------
    // +++ quit 
    //**/ quit psuade interactive session
    //**/ -------------------------------------------------------------
    else if ((!strcmp(command, "quit")) || (!strcmp(command, "q")))
    {
      printf("psuade terminates ...\n");
      break;
    }
    else if ((!strcmp(command, "exit")))
    {
      break;
    }
    else if (command[0] == '#')
    {
      printf("%s", lineIn);
    }
    else if ((!strcmp(command, "script")))
    {
      char scriptName[5001];
      sscanf(lineIn,"%s %s",command,scriptName);
      fp = fopen(scriptName, "r");
      if (fp != NULL && scriptMode == 0)
      {
        scriptFp = fp;
        printf("Script file %s found.\n", scriptName);
        scriptMode = 1;
      }
      else if (fp != NULL && scriptMode == 1)
      {
        printf("ERROR: only one level of script interpretation allowed.\n");
        cmdStatus = 1;
        fclose(fp);
      }
    }
    else
    {
      printf("command %s not recognized\n", command);
      fflush(stdout);
      cmdStatus = 1;
    }
  }

  //**/ ----------------------------------------------------------------
  // quit psuade interactive session
  //**/ ----------------------------------------------------------------
  if (faPtrsRsEval != NULL)
  {
    for (ii = 0; ii < nOutputs_; ii++)
      if (faPtrsRsEval[ii] != NULL) delete faPtrsRsEval[ii];
    delete faPtrsRsEval;
  }
  delete currSession;
  faPtrsRsEval = NULL;
  currSession = NULL;
  return 0;
}

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
void PsuadeBase::displayHelpMenu(char *winput, char *pString)
{
  if (!strcmp(winput, "info"))
  {
    printf("Useful information for using PSUADE :\n");
    printf("\tI.    Uncertainty analysis: \n");
    printf("\t\t 1. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
    printf("\t\t    Analyzer: use 'ua' in command line mode\n");
    printf("\t\t 2. Sampling: first construct response surface (rs)\n");
    printf("\t\t    Analyzer: use rsua\n");
    printf("\t\t 3. Sampling: first construct response surface (rs)\n");
    printf("\t\t    Analyzer: use rsuab (rsua with bootstrapping)\n");
    printf("\tII.   Parameter Screening: \n");
    printf("\t\t 1. Sampling: MOAT, GMOAT\n");
    printf("\t\t    Analyzer: moat, moatmo\n");
    printf("\t\t 2. Sampling: LH, LPTAU\n");
    printf("\t\t    Analyzer: mars_sa/gp_sa (MARS/GP screening)\n");
    printf("\t\t 3. Sampling: LH, LPTAU\n");
    printf("\t\t    Analyzer: delta_test (large sample/low dimension)\n");
    printf("\t\t 4. Sampling: LH, MC\n");
    printf("\t\t    Analyzer: sot_sa (sum-of-trees screening)\n");
    printf("\t\t 5. Sampling: FF4, FF5\n");
    printf("\t\t    Analyzer: ff (fractional factorial analysis)\n");
    printf("\t\t 6. Sampling: LSA\n");
    printf("\t\t    Analyzer: lsa (local sensitivity analysis)\n");
    printf("\tIII.  Classical regression/sensitivity analysis: \n");
    printf("\t\t 1. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
    printf("\t\t    Regression-based correlation analysis (use ca)\n");
    printf("\t\t 2. Sampling: MC, LPTAU, METIS, LH, OA, or OALH\n");
    printf("\t\t    Regression-based (use rsfit/examine SRCs)\n");
    printf("\tIV.   Response surface analysis: \n");
    printf("\t\t 1. Sampling: LPTAU, METIS, LH, OA, or OALH\n");
    printf("\t\t    Response surface validation: a few options\n");
    printf("\t\t    (a) examine R-squared (may not be reliable)\n");
    printf("\t\t    (b) test it on the training set (use rstest_ts)\n");
    printf("\t\t    (c) test it on hold-out set (use rstest_ts)\n");
    printf("\t\t    (d) perform cross validation (use rstest_cv)\n");
    printf("\t\t    (e) perform generalization test (rstest_gt)\n");
    printf("\tV.   Global sensitivity analysis (first order): \n");
    printf("\t\t 1. Sampling: replicated Latin hypercube\n");
    printf("\t\t    Analyzer: use 'me'\n");
    printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
    printf("\t\t    Analyzer: use 'me'\n");
    printf("\t\t 3. Sampling: Sobol' sampling method (large)\n");
    printf("\t\t    Analyzer: use 'sobol'\n");
    printf("\t\t 4. Sampling: any space-filling sample\n");
    printf("\t\t    Analyzer: use 'rsfit' with Legendre (and scaling)\n");
    printf("\t\t 5. Sampling: first construct response surface\n");
    printf("\t\t    Analyzer: use 'rssobol1b, rsmeb' (faster)\n");
    printf("\tVI.  Global sensitivity analysis (second order): \n");
    printf("\t\t 1. Sampling: replicated OA\n");
    printf("\t\t    Analyzer: use 'ie' in command line mode\n");
    printf("\t\t 2. Sampling: any space-filling (large enough) sample\n");
    printf("\t\t    Analyzer: use 'ie'\n");
    printf("\t\t 3. Sampling: first construct response surface\n");
    printf("\t\t    Analyzer: use 'rssobol2, rssobol2b, rsieb' (less robust)\n");
    printf("\tVII.  Global sensitivity analysis (total order): \n");
    printf("\t\t 1. Sampling: use any space-filling (large) sample\n");
    printf("\t\t    Analyzer: use 'tsi' (coarse analysis for low dim)\n");
    printf("\t\t 2. Sampling: Sobol' sampling method (large)\n");
    printf("\t\t    Analyzer: use 'sobol'\n");
    printf("\t\t 3. Sampling: first construct response surface\n");
    printf("\t\t    Analyzer: use 'rssoboltsi, rssoboltsib'\n");
    printf("\t\t 4. Sampling: use FAST sampling\n");
    printf("\t\t    Analyzer: use 'fast'\n");
    printf("\tVIII. Global sensitivity analysis (group): \n");
    printf("\t\t 1. Sampling: first construct response surface\n");
    printf("\t\t    Analyzer: use 'rssobolg'\n");
    printf("\tIX.   Hypothesis testing: \n");
    printf("\t\t 1. Sampling: any of your choice \n");
    printf("\t\t    Analyzer: '1test' or '2test'\n");
    printf("\tX.    Principal component analysis: \n");
    printf("\t\t 1. Sampling: any of your choice \n");
    printf("\t\t    Analyzer: 'pca' in command line mode\n");
    printf("\tXI.   Bayesian inverse UQ: (response surface-based)\n");
    printf("\t\t 1. Sampling: LPTAU, LH, OA, METIS for RS\n");
    printf("\t\t    Analyzer: 'rsmcmc' (in command line mode), or\n");
    printf("\t\t 2. Sampling: LPTAU, LH, OA, METIS for RS\n");
    printf("\t\t    Analyzer: 'mcmc' (simulator-based)\n");
    printf("\tXII.  Advanced features: \n");
    printf("\t\t * Impose constraints in sampling/analysis\n");
    printf("\t\t * Plot PDF of response surface std dev. (rssd_ua)\n");
    printf("\t\t * Intersection or Bayes-like rules (e.g. rsi2)\n");
    printf("\t\t * Multi-objective optimization (mo_opt)\n");
    printf("\t\t * Aleatoric-epistemic uncertainty analysis (aeua)\n");
    printf("\t\t * 2nd order analysis (soua) - input PDF uncertainty\n");
    printf("\t\t * Tools for setting up user application with PSUADE\n");
    printf("\t\t * Sample refinement:\n");
    printf("\t\t       (uniform/adaptive: refine/arefine)\n");
    printf("\t\t * Optimal experimental design\n");
    printf("\t\t * Kernel PCA-based MCMC method\n");
    printf("\t\t * Methods for heteroskedastic models\n");
  }
  else if (!strcmp(winput, "io"))
  {
    printf("Commands for reading/write/updating data to/from files:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\tclear               (Clear PSUADE's workspace) \n");
    printf("\tload/read    <file> (Load data file in PSUADE format) \n");
    printf("\tloadmore     <file> (Add more data to current data set)\n");
    printf("\twrite        <file> (Write to file in PSUADE format)\n");
    printf("\tnwrite       <file> (Write unevaluated points only)\n");
    printf("\tiread        <file> (Read from file with only inputs)\n");
    printf("\tiwrite       <file> (Write to file with only inputs)\n");
    printf("\tiwrite2      <file> (same as iwrite but no header/footer)\n");
    printf("\towrite       <file> (Write to file with only outputs)\n");
    printf("\towrite2      <file> (same as owrite but no header/footer)\n");
    printf("\tread_std     <file> (Read data in standard format)\n");
    printf("\twrite_std    <file> (Write to file in standard format)\n");
    printf("\tread_xls     <file> (Read data in Excel format)\n");
    printf("\twrite_xls    <file> (Write to file in Excel format)\n");
    printf("\tread_csv     <file> (Read data in special CSV format)\n");
    printf("\twrite_csv    <file> (Write to file in special CSV format)\n");
    printf("\twrite_matlab <file> (Write to file in matlab format)\n");
    //printf("\twrite_ultra  <file> (Write to file in ultra format)\n");
    printf("\toupdate      <file> (Update OUTPUTS from a data file)\n");
    printf("\tiadd         <file> (Add more inputs from a data file)\n");
    printf("\toadd         <file> (Add more outputs from a data file)\n");
    printf("\tireplace     <file> (Replace one input from <file>)\n");
    printf("\toreplace     <file> (Replace all outputs from <file>)\n");
    printf("\tssplit              (Split sample into 2 files)\n");
  }
  else if (!strcmp(winput, "stat"))
  {
    printf("Commands for basic statistic on raw sample data:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\tua         (Uncertainty analysis on raw sample)\n");
    printf("\tsem        (Standard error of mean on raw sample)\n");
    printf("\tca         (Correlation analysis on raw sample)\n");
    printf("\tme         (Main effect analysis on raw sample)\n");
    printf("\tie         (pairwise effect analysis on raw sample)\n");
    printf("\ttsi        (Total sensitivity analysis on raw sample)\n");
    printf("\tsobol      (Sensitivity analysis on a Sobol' sample)\n");
    printf("\tfast       (Sensitivity analysis on a FAST sample)\n");
    printf("\tanova      (ANOVA using response surface)\n");
    printf("\t1stest     (1-sample test (Chi-squared, dist fit))\n");
    printf("\t2stest     (2-sample test (T-test,K-S,Mann-Whitney))\n");
    printf("\tgendist    (Create a 1-input sample given an PDF)\n");
    printf("\tpdfconvert (Convert a sample using selected PDFs)\n");
    printf("\trand_draw  (Draw a sample from the loaded sample)\n");
    printf("\trand_draw2 (Draw a sample from 2 files - 2 input sets)\n");
    printf("\trand_drawb (Draw a sample of blocks from the loaded sample)\n");
    printf("\tgensample  (Create a sample from the loaded PDFs)\n");
    printf("\tcdf_lookup (Look up CDF given the random variable value)\n");
    printf("\ticdf_lookup(Look up random variable value given CDF value)\n");
    printf("\tkde        (kernel density estimation: mean, s.d.)\n");
    printf("\tkde2       (kde but works on blocks of identical points)\n");
  }
  else if (!strcmp(winput, "screen"))
  {
    printf("Commands for parameter screening:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\tlsa        (Parameter screening with local SA method)\n");
    printf("\tmoat       (Parameter screening with the Morris method)\n");
    printf("\tmoatmo     (Morris screening for multiple outputs)\n");
    printf("\tff         (Parameter screening with fract. factorial)\n");
    printf("\tmars_sa    (Parameter screening with MARS)\n");
    printf("\tgp_sa      (Parameter screening with Gaussian process)\n");
    printf("\tdelta_test (Parameter screening with Delta test)\n");
    printf("\teta_test   (parameter screening with Eta test)\n");
    printf("\tsot_sa     (Parameter screening with sum-of-trees)\n");
  }
  else if (!strcmp(winput, "rs"))
  { 
    printf("Commands for response surface analysis:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\tset_rstype  (Set rstype in resident sample)\n");
    printf("\tsvmfind     (Search for good parameters for SVM)\n");
    printf("\trsfit       (or rscheck/rsvalidate - check RS quality)\n");
    printf("\trstest_hs   (Check RS quality with another test set)\n");
    printf("\trstest_ts   (Check RS quality with the training set)\n");
    printf("\trstest_cv   (same as rsfit)\n");
    printf("\trstest_rcv  (similar to rstest_cv but more rigorous)\n");
    printf("\trstest_gt   (Check RS quality by generalization test)\n");
    printf("\trstest_sd   (Compute prediction errors on a (large) sample)\n");
    printf("\trscreate    (Create RS to be used by rseval)\n");
    printf("\trseval      (Evaluate RS at given points)\n");
    printf("\trseval_m    (use PSUADE as response surface server)\n");
    printf("\trsevaluate  (Evaluate RS (simpler than rseval))\n");
    printf("\trs_splot    (Create scatter plots - splot - on RS)\n");
    printf("\trsvol       (Compute volume in constrained region)\n");
    printf("\trsint       (Compute volume under response surface)\n");
    printf("\trstgen      (Create a sample (FF/FACT) for rstest_hs)\n");
    printf("\tivec_create (Create input register - used with rseval)\n");
    printf("\tivec_modify (Modify input register - used with rseval)\n");
    printf("\tivec_show   (Display input register - used with rseval)\n");
  }
  else if (!strcmp(winput, "uqsa") && !strcmp(pString, "long"))
  {
    printf("Commands for RS-based uncertainty/sensitivity analysis:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\trsua       (RS-based UA on fuzzy response surface)\n");
  //printf("\trsua2      (rsua but PSUADE is to create evaluation sample)\n");
    printf("\trsuab      (RS-based UA on response surface+bootstrap)\n");
  //printf("\trsuap      (RS-based UA on response surface after rsmcmc)\n");
    printf("\trsmeb      (RS-based McKay main effect + bootstrap)\n");
    printf("\trsieb      (RS-based McKay pairwise effect + bootstrap)\n");
    printf("\trssobol1   (RS-based Sobol' main effect)\n");
    printf("\trssobol2   (RS-based Sobol' interaction effect)\n");
    printf("\trssobolg   (RS-based Sobol' group main effect)\n");
    printf("\trssoboltsi (RS-based Sobol' total effect)\n");
    printf("\trssobol1b  (RS-based Sobol' main effect + bootstrap)\n");
    printf("\trssobol2b  (RS-based Sobol' 2-way analysis + bootstrap)\n");
    printf("\trssoboltsib(RS-based Sobol' total effect + bootstrap)\n");
    printf("\taeua       (RS-based aleatoric-epistemic analysis)\n");
    printf("\tsoua       (RS-based 2nd order analysis: PDF variation)\n");
  }
  else if (!strcmp(winput, "uqsa"))
  {
    printf("Commands for RS-based uncertainty/sensitivity analysis:\n");
    printf("(to see more uqsa commands, use 'help uqsa long')\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\trsua       (RS-based UA)\n");
    printf("\trsuab      (RS-based UA with bootstrap)\n");
    printf("\trsmeb      (RS-based McKay main effect + bootstrapping)\n");
    printf("\trsieb      (RS-based McKay pairwise effect + bootstrap)\n");
    printf("\trssobol1b  (RS-based Sobol' main effect + bootstrap)\n");
    printf("\trssobol2b  (RS-based Sobol' 2-way analysis + bootstrap)\n");
    printf("\trssoboltsib(RS-based Sobol' total effect + bootstrap)\n");
    printf("\taeua       (RS-based aleatoric-epistemic analysis)\n");
    printf("\tsoua       (RS-based 2nd order analysis: PDF variation)\n");
    printf("\tNOTE: use 'help uqsa long' for more detailed information\n");
  }
  else if (!strcmp(winput, "calibration"))
  {
    printf("Commands for Statistical calibration:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\trsmcmc      (Perform RS-based Bayesian inference)\n");
    printf("\tmcmc_predict(Evaluate new points using MCMC posteriors)\n");
    printf("\tmcmc_genplot(Create Matlab plot for posterior samples)\n");
    printf("\tmcmc_dm     (Search good discrepancy model for inference)\n");
    printf("\tmcmc_set_option (set rsmcmc option)\n");
    printf("\tmcmc        (Perform simulator-based - not RS-based MCMC)\n");
    printf("\tmo_opt      (Perform RS-based multi-objective optimization)\n");
  }
  else if (!strcmp(winput, "optimization"))
  {
    printf("Commands for Numerical Optimization:\n");
    printf("There are several numerical optimization available in PSUADE.\n");
    printf("These are, however, not available in command mode. They are\n");
    printf("available in batch mode by using a standard PSUADE input file\n");
    printf("and set the appropriate 'optimization method' field.\n");
    printf("The available optimization methods are:\n");
    printf("- Bobyqa for bound-constrained continuous optimization\n");
    printf("- Lincoa for linearly-constrained continuous optimization\n");
    printf("- Cobya  for inequality-constrained continuous optimization\n");
    printf("- Newuoa for unconstrained continuous optimization\n");
    printf("- SCE    for continuous genetic algorithm-like optimization\n");
    printf("- LBFGS  for derivative-based continuous optimization\n");
    printf("- MOO    for multi-objective optimization\n");
    printf("- Nomad  for continuous and/or integer optimization\n");
    printf("- OUU    for optimization under uncertainty\n");
  }
  else if (!strcmp(winput, "plot") && !strcmp(pString, "long"))
  { 
    printf("Commands for plotting sample data:\n");
    printf("(to see details of each command, use '-h' option)\n");
    printf("\tsplot      (1-input/1 output scatter plot)\n");
    printf("\tsplot2     (2-input/1 output scatter plot)\n");
    printf("\tsplot3     (3-input/1 output scatter plot)\n");
    printf("\tsplot3m    (3-input/1 output scatter plot movie)\n");
    printf("\trs1        (1-input response surface)\n");
    printf("\trs1s       (1-input response surface with std dev)\n");
    printf("\trs2        (2-input response surface)\n");
    printf("\trs3        (3-input response surface)\n");
    printf("\trs3m       (3-input RS (movie), output in z-axis)\n");
    printf("\trs4        (4-input RS plot (movie))\n");
    printf("\trssd       (RS plot of response surface errors)\n");
    printf("\trssd_ua    (Histogram of response surface errors)\n");
    printf("\trsi2       (2-input RS intersection plot)\n");
    printf("\trsi3       (3-input RS intersection plot)\n");
    printf("\trsi3m      (3-input RS intersection plot, 2D movie)\n");
    printf("\trawi2      (2-input intersection plot on raw data)\n");
    printf("\trawi3      (3-input intersection plots on raw data)\n");
    printf("\trspairs    (2-input RS plots for all input pairs)\n");
    printf("\trsipairs   (2-input intersection RS plots)\n");
    printf("\tiplot1     (1-input (input only) scatter plot)\n");
    printf("\tiplot2     (2-input (input only) scatter plot)\n");
    printf("\tiplot3     (3-input (input only) scatter plot)\n");
    printf("\tiplot4m    (4-input (input only) scatter plot): movie\n");
    printf("\tiplot2_all (all pairs 2-input plots)\n");
    printf("\tiplot_pdf  (Sample PDF of selected inputs)\n");
    printf("\tiplot2_pdf (Sample PDFs of all input pairs)\n");
    printf("\toplot2     (2-output scatter plot)\n");
    printf("\toplot_pdf  (Sample PDF of selected outputs)\n");
    printf("\toplot2_pdf (Sample PDF of all output pairs)\n");
    printf("\tihist      (Histogram for a selected input)\n");
    printf("\tihist2     (Histogram for a selected input pair)\n");
    printf("\tiotrace    (Inputs/outputs plots for each sample)\n");
    printf("\tmetis_2dmap(create matlab plot of metis partition (2d only))\n");
    //printf("\tmeplot  (plot main effects with Pgplot)\n");
    //printf("\tmeplot2 (plot main effect/interaction with Pgplot)\n");
    //printf("\trsplot  (plot response surface with Pgplot)\n");
  }
  else if (!strcmp(winput, "plot"))
  { 
    printf("Commands for plotting sample data:\n");
    printf("(to see more plot commands, use 'help plot long')\n");
    printf("\tsplot      (1-input/1 output scatter plot)\n");
    printf("\tsplot2     (2-input/1 output scatter plot)\n");
    printf("\tsplot3     (3-input/1 output scatter plot)\n");
    printf("\tsplot3m    (3-input/1 output scatter plot movie)\n");
    printf("\trs1        (1-input response surface)\n");
    printf("\trs1s       (1-input response surface with std dev)\n");
    printf("\trs2        (2-input response surface)\n");
    printf("\trs3        (3-input response surface)\n");
    printf("\trs3m       (3-input RS (movie), output in z-axis)\n");
    printf("\trs4        (4-input RS plot (movie))\n");
    printf("\trssd       (RS plot of response surface errors)\n");
    printf("\trssd_ua    (Histogram of response surface errors)\n");
    printf("\trsi2       (2-input RS intersection plot)\n");
    printf("\trsi3       (3-input RS intersection plot)\n");
    printf("\trsi3m      (3-input RS intersection plot, 2D movie)\n");
    printf("\trspairs    (2-input RS plots for all input pairs)\n");
    printf("\tiplot1     (1-input (input only) scatter plot)\n");
    printf("\tiplot2     (2-input (input only) scatter plot)\n");
    printf("\tiplot3     (3-input (input only) scatter plot)\n");
    printf("\tiplot4m    (4-input (input only) scatter plot): movie\n");
    printf("\tiplot2_all (all pairs 2-input plots)\n");
    printf("\tiplot_pdf  (Sample PDF of selected inputs)\n");
    printf("\tiplot2_pdf (Sample PDFs of all input pairs)\n");
    printf("\toplot2     (2-output scatter plot)\n");
    printf("\toplot_pdf  (Sample PDF of selected outputs)\n");
    printf("\toplot2_pdf (Sample PDF of all output pairs)\n");
    printf("\tmetis_2dmap(create matlab plot of metis partition (2d only))\n");
    //printf("\tihist      (Histogram for a selected input)\n");
    //printf("\tihist2     (Histogram for a selected input pair)\n");
    printf("\tNOTE: use 'help plot long' for more detailed information\n");
  }
  else if (!strcmp(winput, "setup"))
  {
    printf("Commands for setting up/monitoring work flow:\n");
    printf("\tsetupguide     (Info on how to set up application)\n");
    printf("\tgeninputfile   (Create an input file for psuade)\n");
    printf("\tgenbatchfile   (Create a LLNL-specific batch file)\n");
    printf("\tgendriver      (Create an application driver)\n");
    printf("\tgenexample     (Create a demonstration example)\n");
    printf("\tchkjobs        (Check job status and create a report)\n");
  }
  else if (!strcmp(winput, "edit"))
  {
    printf("Commands for manipulating/displaying resident sample:\n");
    printf("\tsvalidate     (Set 'valid' flag of selected sample points)\n");
    printf("\tsinvalidate   (Reset 'valid' flag of selected sample points)\n");
    printf("\tsrandomize    (Shuffle sample points randomly)\n");
    printf("\tiadd1         (Add 1 input to the loaded sample)\n");
    printf("\toadd1         (Add 1 output to the loaded sample)\n");
    printf("\timodify       (Modify an input of a selected sample point)\n");
    printf("\tomodify       (Modify an output of a selected sample point)\n");
    printf("\tifilter       (Take out points outside some input bounds)\n");
    printf("\tofilter       (Take out points outside some output bounds)\n");
    printf("\tsfilter       (Take out points already exist in another file)\n");
    printf("\tidelete       (Delete a subset of inputs from the sample)\n");
    printf("\todelete       (Delete a subset of outputs from the sample)\n");
    printf("\tsdelete       (Delete selected sample points)\n");
    printf("\tikeep         (Keep a subset of inputs in the sample)\n");
    printf("\tokeep         (Keep a subset of outputs in the sample)\n");
    printf("\tskeep         (Keep a subset of sample points/remove the rest)\n");
    printf("\tspurge        (Take out invalid sample points)\n");
    printf("\tishuffle      (re-order the sample inputs)\n");
    printf("\toshuffle      (re-order the sample outputs)\n");
    printf("\tsshow         (Display one sample point)\n");
    printf("\tsinfo         (Display information on resident sample)\n");
    printf("\tlist1         (List 1-input/1-output data pair)\n");
    printf("\tlist2         (List 2-inputs/1-output data pair)\n");
    printf("\tlistall       (List all inputs/all outputs data)\n");
    printf("\tisort         (sort sample based on input values)\n");
    printf("\tosort         (sort sample based on output values)\n");
    printf("\timax          (Find sample point with max input value)\n");
    printf("\tomax          (Find sample point with max output value)\n");
    printf("\timin          (Find sample point with min input value)\n");
    printf("\tomin          (Find sample point with min output value)\n");
    printf("\tonorm         (Compute the 2-norm of a sample output)\n");
    printf("\tosum          (Compute the sum of a sample output)\n");
    printf("\tinormalize    (normalize a selected sample input)\n");
    printf("\tonormalize    (normalize a selected sample output)\n");
    printf("\tiscale        (rescale a selected input in the sample)\n");
    printf("\toscale        (rescale a selected output in the sample)\n");
    printf("\tireset        (Reset a selected input to some value)\n");
    printf("\toreset        (Reset a selected output to some value)\n");
    printf("\tifloor        (Truncate an input to integer)\n");
    printf("\ticeil         (Round an input to integer)\n");
    printf("\tiround        (Round an input to the nearest integer)\n");
    printf("\titag          (Tag sample points based on input values)\n");
    printf("\totag          (Tag sample points based on output value)\n");
    printf("\trm_dup        (Remove duplicate sample points)\n");
    printf("\tsopt_mmd      (Optimize the loaded sample based on MMD)\n");
    printf("\tiop           (Perform transformation on sample inputs)\n");
    printf("\toop           (Perform transformation on sample outputs)\n");
    printf("\tioop          (Replace Out1 = a * Inp1 + b * Inp2)\n");
    printf("\tsetdriver <s> (Set the driver field in the sample to <s>)\n");
  }
  else if (!strcmp(winput, "misc"))
  {
    printf("Miscellaneous commands:\n");
    printf("\trename <1> <2> (rename a file <1> to <2>) \n");
    printf("\trun <file>     (Run a psuade input script) \n");
    printf("\tquit (exit)    (Terminate command line session)\n");
    printf("\tgetstatus      (Show previous command status: 1 - error)\n");
    printf("\trsmax <d>      (Set maximum no. of data points for RS)\n");
    printf("\tsys <command>  (Execute a system command)\n");
    printf("\tprintlevel <d> (Set print level)\n");
    //printf("\tinteractive    (Turn on/off interactive mode: obsolete)\n");
    printf("\toutput_file    (Set the default output file name)\n");
    printf("==> These commands are for checking sample:\n");
    printf("\tsqc            (Sample quality check: distance metric)\n");
    printf("\tssc            (Sample smoothness check)\n");
    printf("\tnna            (Nearest neighbor analysis: for outliers\n");
    printf("==> Other miscellaneous commands:\n");
    printf("\tsetranseed <d> (set random number generator seed\n");
    printf("\tshowformat     (Formats for RS/MOAT constraint file)\n");
    printf("\tscilab         (Turn on/off scilab (limited support))\n");
    printf("\tstart_matlab   (Start matlab in PSUADE command mode)\n");
    printf("\tcheckformat    (Check format of various PSUADE files)\n");
    printf("\tset_sam_method (Set sampling method in loaded session)\n");
    printf("\tconvhull2d     (Find out if a point is inside a convex hull)\n");
  }
  else if (!strcmp(winput, "advanced"))
  {
    printf("Advanced analysis and control commands:\n");
    printf("\tscript <file>   (Interpret commands in <file>)\n");
    printf("==> These commands turns on/off different expert modes:\n");
    printf("\tio_expert       (Turn on/off IO expert mode)\n");
    printf("\trs_expert       (Turn on/off RS expert mode)\n");
    printf("\trs_codegen      (Turn on/off RS code generator)\n");
    printf("\tana_expert      (Turn on/off analysis expert mode)\n");
    printf("\tsam_expert      (Turn on/off sampling expert mode)\n");
    printf("\topt_expert      (Turn on/off optimization expert mode)\n");
    printf("==> These commands are for generating scenarios for OUU:\n");
    printf("\tgenhistogram    (Create histogram from loaded sample)\n");
    printf("\tgenhistogram2   (genhistogram with user-specified nbins)\n");
    printf("==> These commands are for managing the configuration table:\n");
    printf("\tgenconfigfile   (Create a template configuration file)\n");
    printf("\tsetconfigoption (Set option in configuration table)\n");
    printf("\tshowconfigtable (show content of configuration table)\n");
    printf("\tuseconfigfile <file> (Use a user config file)\n");
    printf("\tresetconfig     (remove all configuration table contents)\n");
    printf("==> These commands are for adaptive sample refinement:\n");
    printf("\turefine         (Refine a sample (uniform refinement))\n");
    printf("\tarefine         (Std.dev.-based adaptive sample refinement)\n");
    printf("\tarefine_metis   (arefine with initial sample=METIS)\n");
    printf("\tarefine_cv      (CV-based adaptive sample refinement)\n");
    printf("\tarel_fact       (Adaptive reliability using factorial design)\n");
    printf("\tarel_fact2      (Same as arel_fact but for nInputs=2)\n");
    printf("\tarel_metis      (Adaptive reliability using METIS design)\n");
    //printf("\tarefinem_thresh (Same as arel_fact except use METIS\n");
    printf("==> These commands are for modifying MOAT with constraints:\n");
    printf("\tmoat_adjust  <file> (Modify MOAT sample from <file>)\n");
    printf("\tgmoat_adjust <file> (Modify GMOAT sample from <file>)\n");
    printf("\tmoatgen             (Create MOAT adjust file: 1 constr)\n");
    printf("\tmoatgen2            (Moatgen with multiple constraints)\n");
    printf("\tmoat_concat  <file> (Combine 2 MOAT input samples)\n");
    printf("==> This command is for finding the interface plane in the space\n");
    printf("    (nInputs-1 dimensions) that separates outputs that are either\n");
    printf("    0 or 1. Then it uses polynomial regression to characterize\n");
    printf("    this interface.\n");
    printf("\tinterface_track (Track interface separated by a threshold)\n");
  }
  else if (!strcmp(winput, "future"))
  {
    printf("\tgen_discrete (Generate a sample of discrete variables)\n");
    printf("\tpdfcheck     (Internal self-check for accuracy of PDFs)\n");
    printf("\tgp_sa2       (Parameter screening via layered GP)\n");
    printf("\tsot_sa2      (Screening via boosted sum-of-trees)\n");
    printf("\tgd_test      (Gower/Mahalanobis extrapolation analysis)\n");
    printf("\trsevaluate2  (Special RS evaluation at given points)\n");
  }
  else if (!strcmp(winput, "odoe") && !strcmp(pString, "long"))
  {
    printAsterisks(PL_INFO, 0);
    displayHelpMenuODOE(1);
  }
  else if (!strcmp(winput, "odoe"))
  {
    displayHelpMenuODOE(0);
  }
  else if (!strcmp(winput, "kpca") && !strcmp(pString, "long"))
  {
    printEquals(PL_INFO, 0);
    printf("\tkpca_create   (Create a KPCA model file from snapshots.\n");
    printf("* This command takes an ensemble of snapshots from user ");
    printf("and creates\n");
    printf("* a KPCA model file for kpca_forward and kpca_inverse.\n");
    printDashes(PL_INFO, 0);
    printf("\tkpca_forward  (Transform from physical to KPCA space)\n");
    printf("\tkpca_inverse  (Transform from KPCA to physical space)\n");
    printf("\tkpca_server   (Run kpca_inverse continuously in server mode)\n");
    printf("* This command runs continuously waiting for a ");
    printf("user-generated reduced\n");
    printf("* parameter file, performing inverse KPCA, and ");
    printf("writing the result to\n");
    printf("a file.\n");
    printDashes(PL_INFO, 0);
    printf("\tgen_likelihood(Create Python script to compute likelihood\n");
    printf("* This command creates a Python script for ");
    printf("computing likelihood values\n");
    printf("* for MCMC. It is especially made to work with 'mcmc_with_dr'.\n");
    printDashes(PL_INFO, 0);
    printf("\tmcmc_with_dr  (Perform MCMC with dimension reduction)\n");
    printf("* This command calls the Metropolis-Hastings MCMC using likelihood\n");
    printf("* function created by gen_likelihood.\n");
    printDashes(PL_INFO, 0);
    printf("\topca          (Principal component analysis: on outputs)\n");
    printEquals(PL_INFO, 0);
  }
  else if (!strcmp(winput, "kpca"))
  {
    printf("\tkpca_create   (Create a KPCA model file from snapshots.\n");
    printf("\tkpca_forward  (Transform from physical to KPCA space)\n");
    printf("\tkpca_inverse  (Transform from KPCA to physical space)\n");
    printf("\tkpca_server   (Run kpca_inverse continuously in server mode)\n");
    printf("\tgen_likelihood(Create Python script to compute likelihood)\n");
    printf("\tmcmc_with_dr  (Perform MCMC with dimension reduction)\n");
    printf("\topca          (Principal component analysis: on outputs)\n");
    printf("\tNOTE: use 'help kpca long' for more detailed information\n");
  }
  else
  {
    printf("Help topics:\n");
    printf("\tinfo         (Information about the use of PSUADE)\n");
    printf("\tio           (File read/write commands)\n");
    printf("\tstat         (Basic statistics)\n");
    printf("\tscreen       (Parameter screening commands)\n");
    printf("\trs           (Response surface analysis commands)\n");
    printf("\tuqsa         (Quantitative UA/SA commands)\n");
    printf("\tcalibration  (Statistical Calibration commands)\n");
    //printf("\toptimization (Numerical optimization)\n");
    printf("\todoe         (Optimal experimental design)\n");
    printf("\tkpca         (Commands that support KPCA-based MCMC)\n");
    printf("\tplot         (Commands for creating plots)\n");
    printf("\tsetup        (Commands to set up PSUADE work flow)\n");
    printf("\tedit         (Commands to edit resident sample data)\n");
    printf("\tmisc         (Miscellaneous commands)\n");
    printf("\tadvanced     (Advanced analysis and control commands)\n");
    printf("\t<command -h> (Help for a specific command)\n");
  }
}

// ************************************************************************
// get sample data from psuadeIO object
// ------------------------------------------------------------------------
void PsuadeBase::getSampleFromPsuadeIO()
{
  pData  pPtr, pLower, pUpper, pINames, pONames, pMeans, pStds;
  pData  pPDFSIndices, pPDFs;
  
  psuadeIO_->getParameter("input_ninputs", pPtr);
  nInputs_ = pPtr.intData_;
  psuadeIO_->getParameter("output_noutputs", pPtr);
  nOutputs_ = pPtr.intData_;
  psuadeIO_->getParameter("input_lbounds", pLower);
  VecILowerBs_.load(nInputs_, pLower.dbleArray_);
  psuadeIO_->getParameter("input_ubounds", pUpper);
  VecIUpperBs_.load(nInputs_, pUpper.dbleArray_);
  psuadeIO_->getParameter("input_names", pINames);
  StrInpNames_.load(nInputs_, (const char **) pINames.strArray_);
  psuadeIO_->getParameter("output_names", pONames);
  StrOutNames_.load(nOutputs_, (const char **) pONames.strArray_);
  psuadeIO_->getParameter("input_pdfs", pPDFs);
  VecInpPDFs_.load(nInputs_, pPDFs.intArray_);
  psuadeIO_->getParameter("input_means", pMeans);
  VecInpMeans_.load(nInputs_, pMeans.dbleArray_);
  psuadeIO_->getParameter("input_stdevs", pStds);
  VecInpStds_.load(nInputs_, pStds.dbleArray_);
  psuadeIO_->getParameter("input_sample_files", pPDFs);
  SamPDFFiles_ = pPDFs.strArray_;
  pPDFs.strArray_ = NULL;
  psuadeIO_->getParameter("input_sample_indices",pPDFSIndices);
  VecSamPDFIndices_.load(nInputs_, pPDFSIndices.intArray_);
  psuadeIO_->getParameter("input_cor_matrix", pPtr);
  if (inputCMat_ != NULL) delete inputCMat_;
  inputCMat_ = new psMatrix();
  psMatrix *tmpMat = (psMatrix *) pPtr.psObject_;
  inputCMat_->load(*tmpMat);

  psuadeIO_->getParameter("method_sampling", pPtr);
  SamplingMethod_ = pPtr.intData_;
  psuadeIO_->getParameter("method_nsamples", pPtr);
  nSamples_ = pPtr.intData_;
  psuadeIO_->getParameter("method_nreplications",pPtr);
  nReplications_ = pPtr.intData_;
  psuadeIO_->getParameter("input_sample", pPtr);
  if (pPtr.dbleArray_ != NULL)
    VecSamInputs_.load(nSamples_*nInputs_, pPtr.dbleArray_);
  pPtr.clean();
  psuadeIO_->getParameter("output_sample", pPtr);
  if (pPtr.dbleArray_ != NULL)
    VecSamOutputs_.load(nSamples_*nOutputs_,pPtr.dbleArray_);
  pPtr.clean();
  psuadeIO_->getParameter("output_states", pPtr);
  if (pPtr.intArray_ != NULL)
    VecSamStates_.load(nSamples_,pPtr.intArray_);
  pPtr.clean();
}

// ************************************************************************
// collapse sample
// ------------------------------------------------------------------------
void PsuadeBase::collapseSample()
{
  int    sInd, ii, kk, jj, iInd, oInd, option, iOne=1;
  double ddata;
  char   charString[10000], pString[10000];
  FILE   *fp=NULL;
  psVector  vecXT, vecYT, vecWT;
  psIVector vecTags;

  printf("Select from the options below: \n");
  printf("(1) remove duplicate sample points (keep first instance)\n");
  printf("(2) combine duplicate sample points by taking average\n");
  printf("(3) combine duplicate sample points by taking median\n");
  printf("(4) combine duplicate sample points by taking std dev\n");
  printf("(5) combine duplicate sample points by taking max\n");
  printf("(6) partition duplicate sample points into bins (quantiles)\n");
  sprintf(charString,"Select option (1-6): ");
  option = getInt(1, 6, charString);
  vecTags.setLength(nSamples_);
  for (sInd = 0; sInd < nSamples_; sInd++) vecTags[sInd] = 1;
  if (option == 1)
  {
    kk = 0;
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,nSamples_,nInputs_,
              VecSamInputs_.getDVector(),vecTags.getIVector());
      if (jj < 0 || jj > sInd)
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
          VecSamInputs_[kk*nInputs_+iInd] = 
                 VecSamInputs_[sInd*nInputs_+iInd];
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[kk*nOutputs_+oInd] = 
                 VecSamOutputs_[sInd*nOutputs_+oInd];
        VecSamStates_[kk] = VecSamStates_[sInd]; 
        kk++;
      }
    }
    nSamples_ = kk;
  }
  else if (option == 2)
  {
    vecTags.setLength(nSamples_);
    for (sInd = 0; sInd < nSamples_; sInd++) vecTags[sInd] = 0;
    kk = 0;
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,kk,nInputs_,VecSamInputs_.getDVector(),
                          vecTags.getIVector());
      if (jj < 0 || jj > sInd)
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
          VecSamInputs_[kk*nInputs_+iInd] = 
                 VecSamInputs_[sInd*nInputs_+iInd];
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[kk*nOutputs_+oInd] = 
                 VecSamOutputs_[sInd*nOutputs_+oInd];
        VecSamStates_[kk] = VecSamStates_[sInd]; 
        vecTags[kk] = 1;
        kk++;
      }
      else
      {
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[jj*nOutputs_+oInd] += 
                 VecSamOutputs_[sInd*nOutputs_+oInd];
        vecTags[jj]++;
        if (outputLevel_ > 2)
          printf("Sample %d has the same input values as sample %d\n",
                 sInd+1, jj+1);
      }
    }
    for (sInd = 0; sInd < kk; sInd++)
      for (oInd = 0; oInd < nOutputs_; oInd++)
        VecSamOutputs_[sInd*nOutputs_+oInd] /= (double) vecTags[sInd];
    nSamples_ = kk;
  }
  //**/ find median
  else if (option == 3)
  {
    int nUniques=0;
    vecTags.setLength(nSamples_);
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,nSamples_,nInputs_,
              VecSamInputs_.getDVector(),vecTags.getIVector());
      //**/ first encounter
      if (jj < 0 || jj > sInd) vecTags[sInd] = nUniques++;
      else                     vecTags[sInd] = vecTags[jj];
    }
#if 0
    //**/ write out the groups 
    if (outputLevel_ > 1) 
    {
      fp = fopen("psuadeRmDupMedian.m", "w");
      for (kk = 0; kk < nUniques; kk++)
      {
        sprintf(charString, "A%d = [", kk+1); 
        if (fp != NULL) fprintf(fp, "%s\n", charString);
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (vecTags[sInd] == kk)
          {
            for (ii = 0; ii < nInputs_; ii++)
              fprintf(fp,"%16.8e ",VecSamInputs_[sInd*nInputs_+ii]); 
            for (ii = 0; ii < nOutputs_; ii++)
              fprintf(fp,"%16.8e ",VecSamOutputs_[sInd*nOutputs_+ii]); 
            fprintf(fp,"\n");
          }
        }
        fprintf(fp, "];\n");
      }
      fclose(fp);
      printf("psuadeRmDupMedians.m has compressed sample\n");
    }
#endif
    //**/ copy unique inputs
    for (kk = 0; kk < nUniques; kk++)
    {
      for (sInd = 0; sInd < nSamples_; sInd++)
      {
        if (vecTags[sInd] == kk)
        {
          for (ii = 0; ii < nInputs_; ii++)
            VecSamInputs_[kk*nInputs_+ii] = 
               VecSamInputs_[sInd*nInputs_+ii]; 
          break;
        }
      }
    }
    //**/ copy unique inputs
    vecYT.setLength(nSamples_*nOutputs_);
    for (ii = 0; ii < nOutputs_; ii++)
    {
      for (kk = 0; kk < nUniques; kk++)
      {
        jj = 0;
        for (sInd = 0; sInd < nSamples_; sInd++)
        {
          if (vecTags[sInd] == kk)
            vecYT[jj++] = VecSamOutputs_[sInd*nOutputs_+ii];
        }
        sortDbleList(jj, vecYT.getDVector());
        VecSamOutputs_[kk*nOutputs_+ii] = vecYT[jj/2];
      }
    }
    vecXT.setLength(nUniques*nInputs_);
    for (ii = 0; ii < nUniques*nInputs_; ii++)
      vecXT[ii] = VecSamInputs_[ii];
    VecSamInputs_ = vecXT; 
    vecYT.setLength(nUniques*nOutputs_);
    for (ii = 0; ii < nUniques*nOutputs_; ii++)
      vecYT[ii] = VecSamOutputs_[ii];
    VecSamOutputs_ = vecYT; 
    kk = nUniques;
    nSamples_ = kk;
  }
  //**/ find standard deviation 
  else if (option == 4)
  {
    vecTags.setLength(nSamples_);
    for (sInd = 0; sInd < nSamples_; sInd++) vecTags[sInd] = 0;
    vecXT.setLength(nSamples_*nInputs_);
    vecYT.setLength(nSamples_*nOutputs_);
    for (sInd = 0; sInd < nSamples_*nInputs_; sInd++)
      vecXT[sInd] = VecSamInputs_[sInd];
    for (sInd = 0; sInd < nSamples_*nOutputs_; sInd++)
      vecYT[sInd] = VecSamOutputs_[sInd];
    // compute average
    kk = 0;
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,kk,nInputs_,vecXT.getDVector(),
                          vecTags.getIVector());
      if (jj < 0 || jj > sInd)
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
          vecXT[kk*nInputs_+iInd] = vecXT[sInd*nInputs_+iInd];
        for (oInd = 0; oInd < nOutputs_; oInd++)
          vecYT[kk*nOutputs_+oInd] = vecYT[sInd*nOutputs_+oInd];
        vecTags[kk] = 1;
        kk++;
      }
      else
      {
        for (oInd = 0; oInd < nOutputs_; oInd++)
          vecYT[jj*nOutputs_+oInd] += vecYT[sInd*nOutputs_+oInd];
        vecTags[jj]++;
      }
    }
    for (sInd = 0; sInd < kk; sInd++)
      for (oInd = 0; oInd < nOutputs_; oInd++)
        vecYT[sInd*nOutputs_+oInd] /= (double) vecTags[sInd];
    // now vecYT has the average
    kk = 0;
    for (sInd = 0; sInd < nSamples_; sInd++) vecTags[sInd] = 0;
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,kk,nInputs_,VecSamInputs_.getDVector(),
                          vecTags.getIVector());
      if (jj < 0 || jj > sInd)
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
          VecSamInputs_[kk*nInputs_+iInd] = 
                 VecSamInputs_[sInd*nInputs_+iInd];
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[kk*nOutputs_+oInd] = 
             pow(VecSamOutputs_[sInd*nOutputs_+oInd]-
                 vecYT[kk*nOutputs_+oInd],2.0);
        VecSamStates_[kk] = VecSamStates_[sInd]; 
        vecTags[kk] = 1;
        kk++;
      }
      else
      {
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[jj*nOutputs_+oInd] += 
              pow(VecSamOutputs_[sInd*nOutputs_+oInd]-
                  vecYT[jj*nOutputs_+oInd],2.0);
        vecTags[jj]++;
      }
    }
    for (sInd = 0; sInd < kk; sInd++)
    {
      for (oInd = 0; oInd < nOutputs_; oInd++)
      {
        if (vecTags[sInd] > 1)
        {
          VecSamOutputs_[sInd*nOutputs_+oInd] = 
            sqrt(VecSamOutputs_[sInd*nOutputs_+oInd]/
            (double) vecTags[sInd]);
        }
        else
        {
          printf("New sample %d has no duplicate ==> s.d. = 0\n",sInd+1);
          VecSamOutputs_[sInd*nOutputs_+oInd] = 0;
        } 
      }
    }  
    nSamples_ = kk;
  }
  //**/ find max 
  else if (option == 5)
  {
    kk = 0;
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,kk,nInputs_,VecSamInputs_.getDVector(),
                          vecTags.getIVector());
      if (jj < 0 || jj > sInd)
      {
        for (iInd = 0; iInd < nInputs_; iInd++)
          VecSamInputs_[kk*nInputs_+iInd] = 
                 VecSamInputs_[sInd*nInputs_+iInd];
        for (oInd = 0; oInd < nOutputs_; oInd++)
          VecSamOutputs_[kk*nOutputs_+oInd] = 
                 VecSamOutputs_[sInd*nOutputs_+oInd];
        VecSamStates_[kk] = VecSamStates_[sInd]; 
        kk++;
      }
      else
      {
        for (oInd = 0; oInd < nOutputs_; oInd++)
          if (VecSamOutputs_[sInd*nOutputs_+oInd] > 
              VecSamOutputs_[jj*nOutputs_+oInd])
            VecSamOutputs_[jj*nOutputs_+oInd] = 
                    VecSamOutputs_[sInd*nOutputs_+oInd];
      }
    }
    nSamples_ = kk;
  }
  //**/ generate quantiles 
  else if (option == 6)
  {
    printf("NOTE: This option supports only \n");
    printf("(N = sample size, B = block size) :\n");
    printf("      (a) N/B blocks - each has unique sample point.\n");
    printf("      (b) B blocks - each block has the same set ");
    printf("of sample points.\n");
    //**/ check for blocks of constant size
    //**/ blkPattern = 0 : contiguous
    //**/ blkPattern = 1 : alternate 
    int blksize=1, blkPattern=0;
    //**/ first check contiguity
    for (sInd = 1; sInd < nSamples_; sInd++)
    {
      for (ii = 0; ii < nInputs_; ii++)
        if (VecSamInputs_[sInd*nInputs_+ii] != VecSamInputs_[ii])
          break;
      if (ii < nInputs_) break;
      blksize++;
    }
    //**/ fails contiguity, look for alternate repetition
    int nSam = 0, failFlag=0;
    if (blksize == 1)
    {
      blkPattern = 1;
      for (sInd = 1; sInd < nSamples_; sInd++)
      { 
        for (ii = 0; ii < nInputs_; ii++)
          if (VecSamInputs_[sInd*nInputs_+ii] != VecSamInputs_[ii])
            break;
        if (ii == nInputs_) break;
      }
      nSam = sInd;
      for (sInd = 0; sInd < nSam; sInd++)
      {
        for (kk = sInd+nSam; kk < nSamples_; kk+=nSam)
        {
          for (ii = 0; ii < nInputs_; ii++)
            if (VecSamInputs_[sInd*nInputs_+ii] != 
                VecSamInputs_[kk*nInputs_+ii])
              break;
          if (ii != nInputs_)
          {
            failFlag = 1; 
            break;
          }
        }
        if (failFlag == 1) break; 
      }
      blksize = nSamples_ / nSam;
    }
      
    //**/ if blkPattern == 0
    if (blkPattern == 0 && blksize == 1)
    {
      printf("ERROR: The sample does not have contiguous ");
      printf("blocks of identical sample\n");
      printf("       points with equal sizes. So this option ");
      printf("is not allowed.\n");
      return;
    }
    if (blkPattern == 1 && (failFlag == 1 || blksize == 1))
    {
      printf("ERROR: The sample does not have identical sample ");
      printf("points with equal\n");
      printf("       size. So this option is not allowed.\n");
      return;
    }
    if (nSamples_/blksize*blksize != nSamples_)
    {
      printf("ERROR: The sample does not have identical sample ");
      printf("points with equal\n");
      printf("       sizes. So this option is not allowed.\n");
      return;
    }
    if (blkPattern == 1)
    {
      psVector vecXT, vecYT;
      psIVector vecST;
      vecXT = VecSamInputs_;
      vecYT = VecSamOutputs_;
      vecST = VecSamStates_;
      for (sInd = 0; sInd < nSam; sInd++)
      {
        for (kk = 0; kk < blksize; kk++)
        {
          for (ii = 0; ii < nInputs_; ii++)
            VecSamInputs_[(sInd*blksize+kk)*nInputs_+ii] = 
                        vecXT[sInd*nInputs_+ii]; 
          for (ii = 0; ii < nOutputs_; ii++)
            VecSamOutputs_[(sInd*blksize+kk)*nOutputs_+ii] = 
                        vecYT[(kk*nSam+sInd)*nOutputs_+ii]; 
          VecSamStates_[sInd*blksize+kk] = vecST[kk*nSam+sInd]; 
        }
      }
    }

    int ss2;
    for (sInd = blksize; sInd < nSamples_; sInd+=blksize)
    {
      for (ss2 = 1; ss2 < blksize; ss2++)
      {
        for (ii = 0; ii < nInputs_; ii++)
          if (VecSamInputs_[(sInd+ss2)*nInputs_+ii] != 
              VecSamInputs_[sInd*nInputs_+ii])
            break;
        if (ii < nInputs_) break;
      }
      if (ss2 != blksize) break;
    }
    if (sInd != nSamples_)
    {
      printf("ERROR: The sample does not have continguous ");
      printf("blocks of identical sample\n");
      printf("       points with equal size. So this option ");
      printf("is not allowed.\n");
      return;
    } 

    int outputID, nlevels, nUniques=0, maxReps=0, entryValid;
    double PLower=0, PUpper=1;
    sprintf(charString, "Which output is to be used? (1 - %d) ",
            nOutputs_);
    outputID = getInt(1, nOutputs_, charString);
    outputID--;
    sprintf(charString,"How many levels to partition the output? (2-20) ");
    nlevels = getInt(2, 20, charString);
    sprintf(charString,"Lower probability bound for binning (>0, <=0.25) : ");
    while (PLower <= 0 || PLower > 0.5)
      PLower = getDouble(charString);
    sprintf(charString,"Upper probability bound for binning (>=0.75, <1) : ");
    while (PUpper >= 1 || PUpper <= 0.5)
      PUpper = getDouble(charString);
    sprintf(charString, "Bins are equally spaced? (y or n) ");
    getString(charString, pString);
    psVector vecBinLocs;
    vecBinLocs.setLength(nlevels);
    if (pString[0] == 'y')
    {
      for (ii = 0; ii < nlevels; ii++)
        vecBinLocs[ii] = PLower + ii * (PUpper - PLower) / (nlevels - 1); 
    }
    else
    {
      vecBinLocs[0] = PLower;
      vecBinLocs[nlevels-1] = PUpper;
      printf("Level 1 = %e\n", PLower);
      for (ii = 1; ii < nlevels-1; ii++)
      {
        sprintf(charString, "Enter level %d value = ", ii+1);
        entryValid = 0;
        while (entryValid == 0)
        {
          vecBinLocs[ii] = getDouble(charString);
          if (vecBinLocs[ii] <= vecBinLocs[ii-1] || vecBinLocs[ii] >= PUpper)
            printf("ERROR: level %d value has to be > %e and < %e\n",ii+1,
                   vecBinLocs[ii-1],PUpper);
          else entryValid = 1;
        }
      }
    }

    // counting size of each bin 
    vecTags.setLength(nSamples_);
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      jj = compareSamples(sInd,nSamples_,nInputs_,
              VecSamInputs_.getDVector(),vecTags.getIVector());
      if (jj < 0 || jj > sInd) vecTags[sInd]++;
      else
      {
        if (outputLevel_ > 2) 
          printf("Sample %d same as sample %d\n",sInd+1,jj+1);
        vecTags[jj]++;
        vecTags[sInd] = - jj - 1;
      }
    }
    for (sInd = 0; sInd < nSamples_; sInd++)
      if (vecTags[sInd] > maxReps) maxReps = vecTags[sInd];
    printf("Maximum degree of replications = %d\n", maxReps);
    for (sInd = 0; sInd < nSamples_; sInd++)
      if (vecTags[sInd] > 0) nUniques++;
    printf("Number of unique samples = %d\n", nUniques);

    //**/ process each unique points
    if (nlevels > 1) vecXT.setLength(nUniques*nlevels*(nInputs_+1));
    else             vecXT.setLength(nUniques*nInputs_);
    vecYT.setLength(nUniques*nlevels);
    KSDensity *ksd = new KSDensity();
    psVector vecX1, vecP1;
    nUniques = 0;
    fp = fopen("matlab_quantiles.m", "w");
    for (sInd = 0; sInd < nSamples_; sInd++)
    {
      if (vecTags[sInd] >= 0)
      {
        vecWT.setLength(maxReps);
        kk = 0;
        for (ii = sInd+1; ii < nSamples_; ii++)
        {
          if (vecTags[ii] < 0 && (-vecTags[ii]-1) == sInd)
          {
            vecWT[kk] = VecSamOutputs_[ii*nOutputs_+outputID];
            kk++;
          }
        }
        vecWT.subvector(0,kk-1);
        ksd->genDensity1D(vecWT, vecX1, vecP1);
        if (fp != NULL)
        {
          fprintf(fp, "A%d = [\n", nUniques+1);
          for (ii = 0; ii < vecP1.length(); ii++) 
            fprintf(fp, "%e %e\n", vecX1[ii], vecP1[ii]);
          fprintf(fp, "];\n");
          fprintf(fp, "hold off\n");
          fprintf(fp, "plot(A%d(:,1),A%d(:,2)/max(A%d(:,2)))\n", 
                  nUniques+1, nUniques+1, nUniques+1);
          fprintf(fp, "B%d = [\n", nUniques+1);
          for (ii = 0; ii < vecP1.length(); ii++) 
            fprintf(fp, "%e\n", vecWT[ii]);
          fprintf(fp, "];\n");
          fprintf(fp, "hold on\n");
          fprintf(fp, "[nk,xk] = hist(B%d,10);\n", nUniques+1);
          fprintf(fp, "bar(xk,nk/max(nk),1.0)\n");
          fprintf(fp, "pause\n");
        }
        for (ii = 1; ii < vecP1.length(); ii++) 
          vecP1[ii] += vecP1[ii-1];
        for (ii = 0; ii < nlevels; ii++)
        {
          for (jj = 0; jj < nInputs_; jj++)
            vecXT[nUniques*(nInputs_+1)*nlevels+ii*(nInputs_+1)+jj] =
                VecSamInputs_[sInd*nInputs_+jj];
          ddata = vecBinLocs[ii];
          vecXT[nUniques*(nInputs_+1)*nlevels+ii*(nInputs_+1)+nInputs_] = 
               ddata;
          for (jj = 0; jj < vecP1.length(); jj++)
          {
            if (ddata < vecP1[jj]) 
            {
              vecYT[nUniques*nlevels+ii] = vecX1[jj];
              break;
            }
          }
          if (jj == vecP1.length()) 
            printf("WARNING: probability not found.\n");
        }
        nUniques++;
      }
    }
    if (fp != NULL) fclose(fp);
    printf("NOTE: Pointwise distribution plots are in ");
    printf("matlab_quantiles.m for\n");
    printf("      visualization.\n");

    nSamples_ = nUniques * nlevels;

    VecSamInputs_ = vecXT;
    VecSamOutputs_ = vecYT;
    VecSamStates_.setLength(nSamples_);
    for (ii = 0; ii < nSamples_; ii++) VecSamStates_[ii] = 1;

    //**/ if more than one quantile, add one more inputs
    if (nlevels > 1) 
    {
      nInputs_++;
      StrInpNames_.addMoreStrings(iOne);
      strcpy(pString, "Quantile");
      StrInpNames_.loadOneString(nInputs_-1, pString);
    }

    //**/ set output names to be the one pointed to by outputID
    strcpy(pString, StrOutNames_[outputID]);
    StrOutNames_.setNumStrings(iOne);
    StrOutNames_.loadOneString(0, pString);
    
    vecXT.setLength(nInputs_);
    for (ii = 0; ii < nInputs_-1; ii++) vecXT[ii] = VecILowerBs_[ii];
    vecXT[nInputs_-1] = 0;
    VecILowerBs_ = vecXT;
    vecXT.setLength(nInputs_);
    for (ii = 0; ii < nInputs_-1; ii++) vecXT[ii] = VecIUpperBs_[ii];
    vecXT[nInputs_-1] = 1;
    VecIUpperBs_ = vecXT;

    VecInpPDFs_.setLength(nInputs_);
    VecInpMeans_.setLength(nInputs_);
    VecInpStds_.setLength(nInputs_);
    if (inputCMat_ != NULL) 
    {
      delete inputCMat_;
      inputCMat_ = new psMatrix();
      inputCMat_->setDim(nInputs_, nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
        inputCMat_->setEntry(ii,ii,1.0);
    }
  }
}

// ************************************************************************
// check resident sample
// ------------------------------------------------------------------------
int PsuadeBase::checkResidentSampleOutputs()
{
  int ii, jj;
  for (ii = 0; ii < nSamples_; ii++)
  {
    for (jj = 0; jj < nOutputs_; jj++)
      if (VecSamOutputs_[ii*nOutputs_+jj] >= PSUADE_UNDEFINED*0.5)
        break;
    if (jj != nOutputs_) break;
  }
  if (ii != nSamples_) 
  {
    printf("WARNING: some outputs are undefined.\n");
    return 1;
  }
  return 0;
}

// ************************************************************************
// update input bounds
// ------------------------------------------------------------------------
void PsuadeBase::updateSampleInputBounds()
{
  int ii, kk;
  for (kk = 0; kk < nInputs_; kk++)
  {
    for (ii = 0; ii < nSamples_; ii++)
    {
      if (VecSamInputs_[ii*nInputs_+kk] < VecILowerBs_[kk])
        VecILowerBs_[kk] = VecSamInputs_[ii*nInputs_+kk]; 
      if (VecSamInputs_[ii*nInputs_+kk] > VecIUpperBs_[kk])
        VecIUpperBs_[kk] = VecSamInputs_[ii*nInputs_+kk]; 
    }
  }
}
 

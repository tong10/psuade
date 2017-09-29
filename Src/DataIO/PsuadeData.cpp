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
// Functions for PsuadeData 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include "Util/dtype.h"
#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Base/Globals.h"
#include "PsuadeData.h"
#include "Psuade.h"
#include "PsuadeConfig.h"

// ************************************************************************
// global variables 
// ------------------------------------------------------------------------ 
extern char         *psConfigFileName_;
extern PsuadeConfig *psConfig_;

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
PsuadeData::PsuadeData()
{
   strcpy(psuadeFileName_, "psuadeData");
   outputLevel_ = 0;
   writeCnt_ = 0;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------ 
PsuadeData::~PsuadeData()
{ 
}

// ************************************************************************
// read from a Psuade input file
// ------------------------------------------------------------------------ 
int PsuadeData::readPsuadeFile(char *fname)
{
   int   lineLeng=300;
   char  lineIn[300], winput[200];
   const char* keywords[] = {"PSUADE", "INPUT", "OUTPUT", "METHOD",
                             "APPLICATION", "ANALYSIS", "END"};
   FILE* fIn;

   pInput_.reset();
   pOutput_.reset();
   pMethod_.reset();
   pApplication_.reset();
   pAnalysis_.reset();

   fIn = fopen(fname, "r");
   if (fIn == NULL) 
   {
      printf("readPsuadeFile ERROR: file %s not found.\n",
             fname);
      printf("You are currently in directory: ");
      fflush(stdout);
      system("pwd");
      printf("\n");
      fflush(stdout);
      return -1;
   }
   fclose(fIn);

   readPsuadeIO(fname);

   fIn = fopen(fname, "r");
   while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(lineIn,"%s", winput);
      if (strcmp(winput, keywords[0]) == 0) break;
   }
   if (feof(fIn) != 0)
   {
      printf("readPsuadeFile ERROR: keyword %s not found.\n",
             keywords[0]);
      fclose(fIn);
      exit(1);
   }
   while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(lineIn,"%s", winput);
      if      (strcmp(winput, keywords[1]) == 0) readInputSection(fIn);
      else if (strcmp(winput, keywords[2]) == 0) readOutputSection(fIn);
      else if (strcmp(winput, keywords[3]) == 0) readMethodSection(fIn);
      else if (strcmp(winput, keywords[4]) == 0) readApplicationSection(fIn);
      else if (strcmp(winput, keywords[5]) == 0) readAnalysisSection(fIn);
      else if (strcmp(winput, keywords[6]) == 0) break;
      else if (strcmp(winput,"#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else 
      {
         printf("readPsuadeFile ERROR: \n");
         printf("\t\t unrecognized line - %s\n", lineIn);
         exit(1);
      }
   }
   fclose(fIn);
   return 0;
}

// ************************************************************************
// writeToFile
// ------------------------------------------------------------------------ 
void PsuadeData::writePsuadeFile(const char *inFilename, int flag)
{
   int    length, cnt;
   char   command[100], fname[100], fname2a[100], fname2b[100];
   FILE   *fOut;

   if (pInput_.sampleInputs_ == NULL || pOutput_.sampleOutputs_ == NULL)
   {
      if (pInput_.sampleInputs_ == NULL) 
         printf("writePsuadeFile ERROR: no input\n");
      else if (pOutput_.sampleOutputs_ == NULL) 
         printf("writePsuadeFile ERROR: no outputs\n");
      exit(1);
   }

   if      (inFilename      != NULL) strcpy(fname, inFilename);
   else if (psuadeFileName_ != NULL) strcpy(fname, psuadeFileName_);
   else 
   {
      printf("writePsuadeFile ERROR: no file given.\n");
      return;
   }
   length = strlen(fname);
   fOut = fopen(fname, "r");
   if (fOut != NULL)
   {
      fclose(fOut);
      if (writeCnt_ == 0 && flag == 1)
      {
         strcpy(fname2a, fname);
         strcpy(fname2b, fname);
         strcpy(&fname2a[length], ".sav.5");
         strcpy(&fname2b[length], ".sav.6");
         cnt = 0;
         while (cnt <= 4)
         {
            if ((fOut=fopen(fname2a, "r")))
            {
               fclose(fOut);
               fOut = NULL;
               strcpy(command, "mv -f ");
               strncpy(&command[6], fname2a, length+6);
               strcpy(&command[length+12], " ");
               strncpy(&command[length+13], fname2b, length+6);
               strcpy(&command[2*length+19], "\n");
               system(command);
            }
            (fname2a[length+5])--;
            (fname2b[length+5])--;
            cnt++;
         }
         if (fOut != NULL) fclose(fOut);
         fOut = NULL;
         strcpy(command, "mv -f ");
         strcpy(&fname2b[length], ".sav.1");
         strncpy(&command[6], fname, length);
         strcpy(&command[length+6], " ");
         strncpy(&command[length+7], fname2b, length+6);
         strcpy(&command[2*length+13], "\n");
         system(command);
         writeCnt_++;
      }
   }

   fOut = fopen(fname, "w");
   if (fOut == NULL)
   {
      printf("Psuade::writePsuadeFile ERROR opening file %s\n",fname);
      exit(1);
   } 

   writePsuadeIO(fOut, 0);

   fprintf(fOut, "PSUADE\n");
   writeInputSection(fOut);
   writeOutputSection(fOut);
   writeMethodSection(fOut);
   writeApplicationSection(fOut);
   writeAnalysisSection(fOut);
   fprintf(fOut, "END\n");
   fclose(fOut);
   if (outputLevel_ > 3)
      printf("\nwritePsuadeFile: data written in file %s.\n", fname);
}

// ************************************************************************
// request data from this object 
// ------------------------------------------------------------------------ 
int PsuadeData::getParameter(const char *keyword, pData &pd)
{ 
   if (!strncmp(keyword,"input_",6))
        return(getInputParameter(keyword, pd));
   else if (!strncmp(keyword,"output_",7))
        return(getOutputParameter(keyword, pd));
   else if (!strncmp(keyword,"method_",7))
        return(getMethodParameter(keyword, pd));
   else if (!strncmp(keyword,"app_",4))
        return(getApplicationParameter(keyword,pd));
   else if (!strncmp(keyword,"ana_",4))
        return(getAnalysisParameter(keyword,pd));
   return 0;
}

// ************************************************************************
// request for INPUT section parameter 
// ------------------------------------------------------------------------ 
int PsuadeData::getInputParameter(const char *keyword, pData &pd)
{ 
   int  ii, jj, kk, retflag=0;

   if (!strcmp(keyword, "input_ninputs")) pd.intData_ = pInput_.nInputs_;
   else if (!strcmp(keyword, "input_lbounds"))
   {
      pd.dbleArray_ = new double[pInput_.nInputs_];
      for (ii = 0; ii < pInput_.nInputs_; ii++)
         pd.dbleArray_[ii] = pInput_.inputLBounds_[ii];
      pd.nDbles_ = pInput_.nInputs_;
   }
   else if (!strcmp(keyword, "input_ubounds"))
   {
      pd.dbleArray_ = new double[pInput_.nInputs_];
      for (ii = 0; ii < pInput_.nInputs_; ii++)
         pd.dbleArray_[ii] = pInput_.inputUBounds_[ii];
      pd.nDbles_ = pInput_.nInputs_;
   }
   else if (!strcmp(keyword, "input_names"))
   {
      pd.strArray_ = new char*[pInput_.nInputs_];
      for (ii = 0; ii < pInput_.nInputs_; ii++)
      {
         pd.strArray_[ii] = new char[100];
         strcpy(pd.strArray_[ii], pInput_.inputNames_[ii]);
      }
      pd.nStrings_ = pInput_.nInputs_;
   }
   else if (!strcmp(keyword, "input_pdfs"))
   {
      if (pInput_.inputPDFs_ != NULL)
      {
         pd.intArray_ = new int[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            pd.intArray_[ii] = pInput_.inputPDFs_[ii];
         pd.nInts_ = pInput_.nInputs_;
      }
      else pd.nInts_ = 0;
   }
   else if (!strcmp(keyword, "input_means"))
   {
      if (pInput_.inputMeans_ != NULL)
      {
         pd.dbleArray_ = new double[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            pd.dbleArray_[ii] = pInput_.inputMeans_[ii];
         pd.nDbles_ = pInput_.nInputs_;
      }
      else pd.nDbles_ = 0;
   }
   else if (!strcmp(keyword, "input_stdevs"))
   {
      if (pInput_.inputStdevs_ != NULL)
      {
         pd.dbleArray_ = new double[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            pd.dbleArray_[ii] = pInput_.inputStdevs_[ii];
         pd.nDbles_ = pInput_.nInputs_;
      } pd.nDbles_ = 0;
   }
   else if (!strcmp(keyword, "input_aux"))
   {
      if (pInput_.inputAuxs_ != NULL)
      {
         pd.dbleArray_ = new double[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            pd.dbleArray_[ii] = pInput_.inputAuxs_[ii];
         pd.nDbles_ = pInput_.nInputs_;
      } pd.nDbles_ = 0;
   }
   else if (!strcmp(keyword, "input_cor_matrix"))
   {
      pd.psObject_ = (void *) &(pInput_.corMatrix_);
   }
   else if (!strcmp(keyword, "input_symtable"))
   {
      if (pInput_.symbolTable_ != NULL)
      {
         pd.intArray_ = new int[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            pd.intArray_[ii] = pInput_.symbolTable_[ii];
         pd.nInts_ = pInput_.nInputs_;
      }
      else pd.nInts_ = 0;
   }
   else if (!strcmp(keyword, "input_settings"))
   {
      if (pInput_.inputNumSettings_ != NULL)
      {
         pd.nInts_ = pInput_.nInputs_;
         pd.intArray_ = new int[pInput_.nInputs_];
         pd.dbleArray2D_ = new double*[pInput_.nInputs_];
         for (ii = 0; ii < pInput_.nInputs_; ii++)
         {
            jj = pInput_.inputNumSettings_[ii];
            pd.intArray_[ii] = jj;
            if (jj > 0)
            {
               pd.dbleArray2D_[ii] = new double[jj];
               for (kk = 0; kk < jj; kk++)
                  pd.dbleArray2D_[ii][kk] = pInput_.inputSettings_[ii][kk];
            }
            else pd.dbleArray2D_[ii] = NULL;
         }
      }
      else pd.nInts_ = 0;
   }
   else if (!strcmp(keyword, "input_sample"))
   {
      if (pInput_.sampleInputs_ != NULL)
      {
         pd.dbleArray_ = new double[pInput_.nInputs_*pMethod_.nSamples_];
         for (ii = 0; ii < pInput_.nInputs_*pMethod_.nSamples_; ii++)
            pd.dbleArray_[ii] = pInput_.sampleInputs_[ii];
         pd.nDbles_ = pInput_.nInputs_ * pMethod_.nSamples_;
      }
      else pd.nDbles_ = 0;
   }
   else if (!strcmp(keyword, "input_use_input_pdfs"))
      pd.intData_ = pInput_.useInputPDFs_;
   else retflag = 1;
   return retflag;
}

// ************************************************************************
// request for OUTPUT section parameters 
// ------------------------------------------------------------------------ 
int PsuadeData::getOutputParameter(const char *keyword, pData &pd)
{ 
   int  ii, retflag=0;

   if (!strcmp(keyword, "output_noutputs")) pd.intData_ = pOutput_.nOutputs_;
   else if (!strcmp(keyword, "output_names"))
   {
      pd.strArray_ = new char*[pOutput_.nOutputs_];
      for (ii = 0; ii < pOutput_.nOutputs_; ii++)
      {
         pd.strArray_[ii] = new char[100];
         strcpy(pd.strArray_[ii], pOutput_.outputNames_[ii]);
      }
      pd.nStrings_ = pOutput_.nOutputs_;
   }
   else if (!strcmp(keyword, "output_sample"))
   {
      if (pOutput_.sampleOutputs_ != NULL)
      {
         pd.dbleArray_ = new double[pOutput_.nOutputs_*pMethod_.nSamples_];
         for (ii = 0; ii < pOutput_.nOutputs_*pMethod_.nSamples_; ii++)
            pd.dbleArray_[ii] = pOutput_.sampleOutputs_[ii];
         pd.nDbles_ = pOutput_.nOutputs_ * pMethod_.nSamples_;
      }
      else pd.nDbles_ = 0;
   }
   else if (!strcmp(keyword, "output_states"))
   {
      if (pOutput_.sampleStates_ != NULL)
      {
         pd.intArray_ = new int[pMethod_.nSamples_];
         for (ii = 0; ii < pMethod_.nSamples_; ii++)
            pd.intArray_[ii] = pOutput_.sampleStates_[ii];
         pd.nInts_ = pMethod_.nSamples_;
      }
      else pd.nInts_ = 0;
   }
   else retflag = 1;
   return retflag;
}

// ************************************************************************
// request for METHOD section parameters 
// ------------------------------------------------------------------------ 
int PsuadeData::getMethodParameter(const char *keyword, pData &pd)
{
   int  retflag=0;

   if (!strcmp(keyword, "method_sampling"))
      pd.intData_ = pMethod_.samplingMethod_;
   else if (!strcmp(keyword, "method_nsamples"))
      pd.intData_ = pMethod_.nSamples_;
   else if (!strcmp(keyword, "method_nreplications"))
      pd.intData_ = pMethod_.nReplications_;
   else if (!strcmp(keyword, "method_nrefinements"))
      pd.intData_ = pMethod_.nRefinements_;
   else if (!strcmp(keyword, "method_refnrefinements"))
      pd.intData_ = pMethod_.refNumRefinements_;
   else if (!strcmp(keyword, "method_refine_type"))
      pd.intData_ = pMethod_.refinementType_;
   else if (!strcmp(keyword, "method_refine_size"))
      pd.intData_ = pMethod_.refinementSize_;
   else if (!strcmp(keyword, "method_randomflag"))
      pd.intData_ = pMethod_.sampleRandomize_;
   else retflag = 1;
   return retflag;
}

// ************************************************************************
// request for APPLICATION section parameters 
// ------------------------------------------------------------------------ 
int PsuadeData::getApplicationParameter(const char *keyword, pData &pd)
{ 
   int  ii, retflag=0;

   if (!strcmp(keyword, "app_files"))
   {
      pd.strArray_ = new char*[5];
      for (ii = 0; ii < 5; ii++) pd.strArray_[ii] = new char[500];
      strcpy(pd.strArray_[0], pApplication_.appDriver_);
      strcpy(pd.strArray_[1], pApplication_.inputTemplate_);
      strcpy(pd.strArray_[2], pApplication_.outputTemplate_);
      strcpy(pd.strArray_[3], pApplication_.optDriver_);
      strcpy(pd.strArray_[4], pApplication_.auxOptDriver_);
      pd.nStrings_ = 5;
   }
   else if (!strcmp(keyword, "app_maxparalleljobs"))
      pd.intData_ = pApplication_.maxParallelJobs_;
   else if (!strcmp(keyword, "app_minjobwaittime"))
      pd.intData_ = pApplication_.minJobWaitTime_;
   else if (!strcmp(keyword, "app_maxjobwaittime"))
      pd.intData_ = pApplication_.maxJobWaitTime_;
   else if (!strcmp(keyword, "app_runtype"))
      pd.intData_ = pApplication_.runType_;
   else if (!strcmp(keyword, "app_launchinterval"))
      pd.intData_ = pApplication_.launchInterval_;
   else if (!strcmp(keyword, "app_savefrequency"))
      pd.intData_ = pApplication_.saveFrequency_;
   else if (!strcmp(keyword, "app_usefunction"))
      pd.intData_ = pApplication_.useFunction_;
   else retflag = 1;
   return retflag;
}

// ************************************************************************
// request for ANALYSIS section parameters 
// ------------------------------------------------------------------------ 
int PsuadeData::getAnalysisParameter(const char *keyword, pData &pd)
{ 
   int  ii, retflag=0;

   if (!strcmp(keyword, "ana_opt_switch"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[0];
   else if (!strcmp(keyword, "ana_opt_method"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[1];
   else if (!strcmp(keyword, "ana_opt_nlocalmin"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[2];
   else if (!strcmp(keyword, "ana_opt_rstype"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[3];
   else if (!strcmp(keyword, "ana_opt_printlevel"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[4];
   else if (!strcmp(keyword, "ana_opt_outputid"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[6];
   else if (!strcmp(keyword, "ana_opt_numfmin"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[5];
   else if (!strcmp(keyword, "ana_opt_fmin"))
      pd.dbleData_ = pAnalysis_.optimizeDbleOptions_[0];
   else if (!strcmp(keyword, "ana_opt_cutoff"))
      pd.dbleData_ = pAnalysis_.optimizeDbleOptions_[1];
   else if (!strcmp(keyword, "ana_opt_maxfeval"))
      pd.intData_ = pAnalysis_.optimizeIntOptions_[7];
   else if (!strcmp(keyword, "ana_opt_deltax"))
      pd.dbleData_ = pAnalysis_.optimizeDbleOptions_[3];
   else if (!strcmp(keyword, "ana_opt_tolerance"))
      pd.dbleData_ = pAnalysis_.optimizeDbleOptions_[2];
   else if (!strcmp(keyword, "ana_method"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[0];
   else if (!strcmp(keyword, "ana_outputid"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[1];
   else if (!strcmp(keyword, "ana_rstype"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[2];
   else if (!strcmp(keyword, "ana_poly_order"))
      pd.intData_ = pAnalysis_.legendreOrder_;
   else if (!strcmp(keyword, "ana_rsindexfile"))
   {
      pd.nStrings_ = 1;
      pd.strArray_ = new char*[1];
      pd.strArray_[0] = new char[200];
      strcpy(pd.strArray_[0], pAnalysis_.rsIndexFile_);
   }
   else if (!strcmp(keyword, "ana_diagnostics"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[3];
   else if (!strcmp(keyword, "ana_graphicsflag"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[4];
   else if (!strcmp(keyword, "ana_transform"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[5];
   else if (!strcmp(keyword, "ana_regressionwgtid"))
      pd.intData_ = pAnalysis_.analysisIntOptions_[6];
   else if (!strcmp(keyword, "ana_threshold"))
      pd.dbleData_ = pAnalysis_.analysisDbleOptions_[0];
   else if (!strcmp(keyword, "ana_opt_smtargetfile"))
   {
      pd.nStrings_ = 1;
      pd.strArray_ = new char*[1];
      pd.strArray_[0] = new char[200];
      strcpy(pd.strArray_[0], pAnalysis_.specsFile_);
   }
   else if (!strcmp(keyword, "ana_num_rsfilters"))
      pd.intData_ = pAnalysis_.numRSFilters_;
   else if (!strcmp(keyword, "ana_rsfilterlbounds"))
   {
      if (pAnalysis_.numRSFilters_ > 0)
      {
         pd.dbleArray_ = new double[pAnalysis_.numRSFilters_];
         for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
            pd.dbleArray_[ii] = pAnalysis_.RSFilters_[ii]->FilterLBound_;
      }
   }
   else if (!strcmp(keyword, "ana_rsfilterubounds"))
   {
      if (pAnalysis_.numRSFilters_ > 0)
      {
         pd.dbleArray_ = new double[pAnalysis_.numRSFilters_];
         for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
            pd.dbleArray_[ii] = pAnalysis_.RSFilters_[ii]->FilterUBound_;
      }
   }
   else if (!strcmp(keyword, "ana_rsfilterdatafile"))
   {
      if (pAnalysis_.numRSFilters_ > 0)
      {
         pd.nStrings_ = pAnalysis_.numRSFilters_;
         pd.strArray_ = new char*[pAnalysis_.numRSFilters_];
         for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
         {
            pd.strArray_[ii] = new char[500];
            strcpy(pd.strArray_[ii],
                   pAnalysis_.RSFilters_[ii]->FilterDataFile_);
         }
      }
   }
   else if (!strcmp(keyword, "ana_rsfilterindexfile"))
   {
      if (pAnalysis_.numRSFilters_ > 0)
      {
         pd.nStrings_ = pAnalysis_.numRSFilters_;
         pd.strArray_ = new char*[pAnalysis_.numRSFilters_];
         for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
         {
            pd.strArray_[ii] = new char[500];
            strcpy(pd.strArray_[ii],
                   pAnalysis_.RSFilters_[ii]->FilterIndexFile_);
         }
      }
   }
   else if (!strcmp(keyword, "ana_num_moatfilters"))
      pd.intData_ = pAnalysis_.numMOATFilters_;
   else if (!strcmp(keyword, "ana_moatfilterlbounds"))
   {
      if (pAnalysis_.numMOATFilters_ > 0)
      {
         pd.dbleArray_ = new double[pAnalysis_.numMOATFilters_];
         for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
            pd.dbleArray_[ii] = pAnalysis_.MOATFilters_[ii]->FilterLBound_;
      }
   }
   else if (!strcmp(keyword, "ana_moatfilterubounds"))
   {
      if (pAnalysis_.numMOATFilters_ > 0)
      {
         pd.dbleArray_ = new double[pAnalysis_.numMOATFilters_];
         for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
            pd.dbleArray_[ii] = pAnalysis_.MOATFilters_[ii]->FilterUBound_;
      }
   }
   else if (!strcmp(keyword, "ana_moatfilterdatafile"))
   {
      if (pAnalysis_.numMOATFilters_ > 0)
      {
         pd.nStrings_ = pAnalysis_.numMOATFilters_;
         pd.strArray_ = new char*[pAnalysis_.numMOATFilters_];
         for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
         {
            pd.strArray_[ii] = new char[500];
            strcpy(pd.strArray_[ii],
                   pAnalysis_.MOATFilters_[ii]->FilterDataFile_);
         }
      }
   }
   else if (!strcmp(keyword, "ana_moatfilterindexfile"))
   {
      if (pAnalysis_.numMOATFilters_ > 0)
      {
         pd.nStrings_ = pAnalysis_.numMOATFilters_;
         pd.strArray_ = new char*[pAnalysis_.numMOATFilters_];
         for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
         {
            pd.strArray_[ii] = new char[500];
            strcpy(pd.strArray_[ii],
                   pAnalysis_.MOATFilters_[ii]->FilterIndexFile_);
         }
      }
   }
   else if (!strcmp(keyword, "ana_use_input_pdfs"))
      pd.intData_ = pAnalysis_.useInputPDFs_;
   else retflag = 1;
   return retflag;
}

// ************************************************************************
// update the input section
// ------------------------------------------------------------------------ 
void PsuadeData::updateInputSection(int nSamples,int nInputs,int *symTable,
                                    double *lowerB, double *upperB, 
                                    double *sampleInputs, char **names)
{
   int ii;

   if (names != NULL) 
   {
      if (pInput_.inputNames_ != NULL)
      {
         for (ii = 0; ii < pInput_.nInputs_; ii++)
            if (pInput_.inputNames_[ii] != NULL)
               delete [] pInput_.inputNames_[ii]; 
         delete [] pInput_.inputNames_; 
      }
      pInput_.inputNames_ = new char*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         pInput_.inputNames_[ii] = new char[200];
         strcpy(pInput_.inputNames_[ii], names[ii]);
      }
   }
   pInput_.nInputs_ = nInputs;

   if (pInput_.sampleInputs_ != NULL) delete [] pInput_.sampleInputs_; 
   pMethod_.nSamples_ = nSamples;
   pInput_.sampleInputs_ = new double[nSamples * nInputs];
   for (ii = 0; ii < nInputs*nSamples; ii++) 
      pInput_.sampleInputs_[ii] = sampleInputs[ii]; 

   if (lowerB != NULL)
   {
      if (pInput_.inputLBounds_ != NULL) delete [] pInput_.inputLBounds_;
      pInput_.inputLBounds_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++) pInput_.inputLBounds_[ii] = lowerB[ii]; 
   }
   if (upperB != NULL)
   {
      if (pInput_.inputUBounds_ != NULL) delete [] pInput_.inputUBounds_;
      pInput_.inputUBounds_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++) pInput_.inputUBounds_[ii] = upperB[ii]; 
   }

   if (symTable != NULL)
   {
      if (pInput_.symbolTable_ != NULL) delete [] pInput_.symbolTable_;
      pInput_.symbolTable_ = new int[nInputs];
      for (ii = 0; ii < nInputs; ii++)
         pInput_.symbolTable_[ii] = symTable[ii];
   }
   if (pInput_.inputPDFs_ == NULL)
   {
      pInput_.inputPDFs_    = new int[nInputs];
      pInput_.inputMeans_   = new double[nInputs];
      pInput_.inputStdevs_  = new double[nInputs];
      pInput_.inputAuxs_    = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         pInput_.inputPDFs_[ii] = 0;
         pInput_.inputMeans_[ii] = 0.0;
         pInput_.inputStdevs_[ii] = 0.0;
         pInput_.inputAuxs_[ii] = 0.0;
      }
      pInput_.corMatrix_.setDim(nInputs, nInputs);
   }
   pInput_.corMatrix_.setDim(nInputs, nInputs);
} 

// ************************************************************************
// update the output section
// ------------------------------------------------------------------------ 
void PsuadeData::updateOutputSection(int nSamples, int nOutputs,  
                                     double *sampleOutputs, int *sampleStates, 
                                     char **names)
{
   int ii, ss;

   if (names != NULL) 
   {
      if (pOutput_.outputNames_ != NULL)
      {
         for (ii = 0; ii < pOutput_.nOutputs_; ii++)
            if (pOutput_.outputNames_[ii] != NULL)
               delete [] pOutput_.outputNames_[ii]; 
         delete [] pOutput_.outputNames_; 
      }
      pOutput_.outputNames_ = new char*[nOutputs];
      for (ii = 0; ii < nOutputs; ii++)
      {
         pOutput_.outputNames_[ii] = new char[200];
         strcpy(pOutput_.outputNames_[ii], names[ii]);
      }
   }
   pOutput_.nOutputs_ = nOutputs;
   pMethod_.nSamples_ = nSamples;
   if (pOutput_.sampleOutputs_ != NULL) delete [] pOutput_.sampleOutputs_; 
   if (pOutput_.sampleStates_  != NULL) delete [] pOutput_.sampleStates_; 
   pOutput_.sampleOutputs_ = new double[nSamples * nOutputs];
   pOutput_.sampleStates_  = new int[nSamples];
   for (ss = 0; ss < nOutputs*nSamples; ss++)
      pOutput_.sampleOutputs_[ss] = sampleOutputs[ss]; 
   for (ss = 0; ss < nSamples; ss++)
      pOutput_.sampleStates_[ss] = sampleStates[ss]; 
} 

// ************************************************************************
// update method information 
// ------------------------------------------------------------------------ 
void PsuadeData::updateMethodSection(int method, int nSamples, int nReps,
                                     int nRefine, int randomFlag)
{
   if (method     >= 0) pMethod_.samplingMethod_  = method;
   if (nSamples   >= 0) pMethod_.nSamples_        = nSamples;
   if (nReps      >= 0) pMethod_.nReplications_   = nReps;
   if (nRefine    >= 0) pMethod_.nRefinements_    = nRefine;
   if (randomFlag >= 0) pMethod_.sampleRandomize_ = randomFlag;
}

// ************************************************************************
// update application information 
// ------------------------------------------------------------------------ 
void PsuadeData::updateApplicationSection(char *appDriver, char *optDriver,
                                          int maxJobs) 
{
   if (appDriver != NULL) strcpy(pApplication_.appDriver_, appDriver);
   if (optDriver != NULL) strcpy(pApplication_.optDriver_, optDriver);
   if (maxJobs > 0) pApplication_.maxParallelJobs_ = maxJobs;
}

// ************************************************************************
// update analysis information 
// ------------------------------------------------------------------------ 
void PsuadeData::updateAnalysisSection(int method, int transform, int rstype,
                                       int diag, int outputID, int usePDFs)
{
   if (method    >=  0) pAnalysis_.analysisIntOptions_[0] = method;
   if (transform >=  0) pAnalysis_.analysisIntOptions_[5] = transform;
   if (rstype    >=  0) pAnalysis_.analysisIntOptions_[2] = rstype;
   if (diag      >= -1) pAnalysis_.analysisIntOptions_[3] = diag;
   if (outputID  >=  0) pAnalysis_.analysisIntOptions_[1] = outputID;
   if (usePDFs   >=  0) pAnalysis_.useInputPDFs_ = usePDFs;
}

// ************************************************************************
// update optimization information 
// ------------------------------------------------------------------------ 
void PsuadeData::updateOptimizationSection(int method, double tolerance)
{
   if (method >= 0 && method <= 10)
   {
      pAnalysis_.optimizeIntOptions_[1] = method;
      pAnalysis_.optimizeIntOptions_[0] = 1;
   }
   if (tolerance >= 0)
      pAnalysis_.optimizeDbleOptions_[2] = tolerance;
}

// ************************************************************************
// set output level for error output
// ------------------------------------------------------------------------ 
void PsuadeData::setOutputLevel(int level)
{
   outputLevel_ = level;
   if (outputLevel_ < 0) outputLevel_ = 1;
}

// ************************************************************************
// process output data - write to matlab file
// ------------------------------------------------------------------------ 
void PsuadeData::processOutputData()
{
   int    nInputs, nSamples, nOutputs, ss, ii;
   double *sampleInputs, *sampleOutputs;
   FILE   *fp;

   nInputs  = pInput_.nInputs_;
   nOutputs = pOutput_.nOutputs_;
   nSamples = pMethod_.nSamples_;
   sampleInputs  = pInput_.sampleInputs_;
   sampleOutputs = pOutput_.sampleOutputs_;

   if ((pAnalysis_.fileWriteFlag_ & 1) != 0)
   {
      fp = fopen("psuade_matlab.m", "w");
      if ( fp == NULL ) return;
      fprintf(fp, "XY = [\n");
      for (ss = 0; ss < nSamples; ss++) 
      {
         for (ii = 0; ii < nInputs; ii++) 
            fprintf(fp, "   %24.16e\n",
                    sampleInputs[ss*nInputs+ii]);
         for (ii = 0; ii < nOutputs; ii++) 
            fprintf(fp, "   %24.16e\n",
                    sampleOutputs[ss*nOutputs+ii]);
      }
      fprintf(fp, "];\n");
      fprintf(fp, "X = [\n");
      for (ii = 0; ii < nInputs*nSamples; ii++) 
         fprintf(fp, "   %24.16e\n", sampleInputs[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nOutputs*nSamples; ii++) 
         fprintf(fp, "   %24.16e\n", sampleOutputs[ii]);
      fprintf(fp, "];\n");
      for (ii = 0; ii < nInputs; ii++) 
         fprintf(fp, "X%d = X(%d:%d:%d);\n", ii+1, ii+1, nInputs,
                 nSamples*nInputs);
      for (ii = 0; ii < nOutputs; ii++) 
         fprintf(fp, "Y%d = Y(%d:%d:%d);\n", ii+1, ii+1,
                 nOutputs, nSamples*nOutputs);
      fclose(fp);
   }
} 

// ************************************************************************
// internal functions
// ************************************************************************

// ************************************************************************
// A function for reading PSUADE IO data
// ------------------------------------------------------------------------ 
void PsuadeData::readPsuadeIO(char *fname) 
{
   int    nInputs, nOutputs, *sampleStates, nSamples, ss, ii, idata;
   int    found=0;
   double *sampleInputs, *sampleOutputs;
   char   lineInput[500], keyword[500];
   FILE   *fIn;

   fIn = fopen(fname, "r");
   assert(fIn != NULL);
   fgets(lineInput, 500, fIn);
   sscanf(lineInput, "%s", keyword);
   while (keyword[0] == '#')
   {
      fgets(lineInput, 500, fIn);
      sscanf(lineInput, "%s", keyword);
   }
   if (!strcmp(keyword, "PSUADE_IO")) /* data is in this section */
   {
      found = 1;
      fscanf(fIn, "%d %d %d\n", &nInputs, &nOutputs, &nSamples);
      if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
      {
         printf("readPsuadeIO ERROR: first parameters <= 0.\n");
         exit(1);
      }
      sampleInputs  = new double[nInputs*nSamples];
      sampleOutputs = new double[nOutputs*nSamples];
      sampleStates  = new int[nSamples];
      for (ss = 0; ss < nSamples; ss++)
      {
         fscanf(fIn,"%d %d", &idata, &sampleStates[ss]);
         if (idata != (ss+1))
         {
            printf("readPsuadeIO ERROR: incorrect sample no.\n");
            printf("        Incoming/expected sample number = %d %d\n",
                   idata, ss+1);
            exit(1);
         }
         if (sampleStates[ss] != 1) sampleStates[ss] = 0;
         for (ii = 0; ii < nInputs; ii++) 
            fscanf(fIn,"%lg",&sampleInputs[ss*nInputs+ii]);
         for (ii = 0; ii < nOutputs; ii++) 
            fscanf(fIn,"%lg",&sampleOutputs[ss*nOutputs+ii]);
      }
      pInput_.nInputs_        = nInputs;
      pInput_.sampleInputs_   = sampleInputs;
      pOutput_.nOutputs_      = nOutputs;
      pOutput_.sampleOutputs_ = sampleOutputs;
      pOutput_.sampleStates_  = sampleStates;
      pMethod_.nSamples_      = nSamples;
      if (outputLevel_ > 1)
      {
         printf("readPsuadeIO: read sample data completed.\n");
         printf("   nInputs, nOutputs, nSamples = %d %d %d\n", nInputs, 
                nOutputs, nSamples);
      }
   }
   else
   {
      if (outputLevel_ > 0)
         printf("readPsuadeIO: PSUADE_IO section absent.\n");
   }
   fclose(fIn);
}

// ************************************************************************
// A function for reading PSUADE parameters from an input file.
// ------------------------------------------------------------------------ 
void PsuadeData::readInputSection(FILE *fp) 
{
   int    ii, idata, itmp, lineLeng=200, nInputs=0;
   double ddata;
   char   line[200], winput[200], winput2[200], winput3[200];
   const char *keywords[] = {"dimension", "variable", "PDF", "COR", "NAME", 
                             "END"};

   if (fp == NULL)
   {
      printf("readInputSection ERROR: file = NULL\n");
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      strcpy(winput, "#");
      sscanf(line,"%s", winput);
      if (strcmp(winput, keywords[0]) == 0) /* dimension */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2, &nInputs);
         else sscanf(line,"%s %d", winput, &nInputs);
         if (nInputs <= 0)
         {
            printf("readInputSection ERROR: nInputs <= 0.\n");
            exit(1);
         } 
         if (pInput_.nInputs_ != 0 && pInput_.nInputs_ != nInputs)
         {
            printf("readInputSection ERROR: nInputs mismatch.\n");
            printf("        Check nInputs in INPUT & PSUADE_IO sections.\n");
            exit(1);
         } 
         pInput_.nInputs_ = nInputs;
         pInput_.inputLBounds_ = new double[nInputs];
         pInput_.inputUBounds_ = new double[nInputs];
	 pInput_.inputNames_   = new char*[nInputs];
	 pInput_.inputPDFs_    = new int[nInputs];
	 pInput_.inputMeans_   = new double[nInputs];
	 pInput_.inputStdevs_  = new double[nInputs];
	 pInput_.inputAuxs_    = new double[nInputs];
         for (ii = 0; ii < nInputs; ii++)
         {
            pInput_.inputLBounds_[ii] = 0.0;
            pInput_.inputUBounds_[ii] = 1.0;
	    pInput_.inputNames_[ii] = NULL;
	    pInput_.inputPDFs_[ii] = 0;
	    pInput_.inputMeans_[ii] = 0.0;
	    pInput_.inputStdevs_[ii] = 0.0;
	    pInput_.inputAuxs_[ii] = 0.0;
         }
         pInput_.corMatrix_.setDim(nInputs, nInputs);
         ddata = 1.0;
         for (ii = 0; ii < nInputs; ii++)
            pInput_.corMatrix_.setEntry(ii, ii, ddata);
      }
      else if (strcmp(winput, keywords[1]) == 0) /* variable */
      {
         if (nInputs <= 0)
         {
            printf("readInputSection ERROR: nInputs not set.\n");
            exit(1);
         } 
         sscanf(line,"%s %d", winput, &idata);
         idata--;
         if (idata < 0 || idata >= pInput_.nInputs_)
         {
            printf("readInputSection ERROR- invalid input number %d\n",
                   idata+1);
            printf("     input number should be between 1 and %d\n", nInputs);
            exit(1);
         }
         sscanf(line,"%s %d %s", winput, &itmp, winput2);
         if (pInput_.inputNames_[idata] != NULL)
         {
            printf("readInputSection WARNING: \n");
            printf("\t\tinput number %d may have been re-defined.\n", idata+1);
            delete [] pInput_.inputNames_[idata];
         }
         pInput_.inputNames_[idata] = new char[strlen(winput2)+10];
         strncpy(pInput_.inputNames_[idata], winput2, strlen(winput2)+1);
         sscanf(line,"%s %d %s %s", winput, &itmp, winput2, winput3);
         if (strcmp(winput3, "=") != 0 )
         {
            printf("readInputSection ERROR: input format for input %d\n",itmp);
            printf("        syntax: variable 1 X1 = 0.0 1.0\n");
            exit(1);
         }
         sscanf(line,"%s %d %s %s %lg %lg", winput, &itmp, winput2,
                winput3, &(pInput_.inputLBounds_[idata]),
                &(pInput_.inputUBounds_[idata]));
         if (pInput_.inputLBounds_[idata] >= pInput_.inputUBounds_[idata])
         {
            printf("readInputSection ERROR: \n");
            printf("\t\tinput lbound >= ubound (%d %e %e).\n", idata+1,
                   pInput_.inputLBounds_[idata], pInput_.inputUBounds_[idata]);
            exit(1);
         }
      }
      else if (strcmp(winput, keywords[2]) == 0) /* PDF */
      {
         sscanf(line,"%s %d", winput, &idata);
         idata--;
         if (idata < 0 || idata >= pInput_.nInputs_)
         {
            printf("readInputSection ERROR: invalid input number %d.\n",
                   idata+1);
            printf("     input number should be between 1 and %d\n", nInputs);
            exit(1);
         } 
         sscanf(line,"%s %d %s", winput, &itmp, winput2);
         if      ( !strcmp(winput2, "U"))
         {
            pInput_.inputPDFs_[idata] = 0;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            if (pInput_.inputMeans_[idata] != pInput_.inputLBounds_[idata])
            {
               printf("readInputSection ERROR: for uniform distributions, the\n");
               printf("    first data should be the same as the lower bound\n");
               printf("    in the variable definition.\n");
               printf("    lower bound defined  = %e\n", pInput_.inputMeans_[idata]);
               printf("    lower bound expected = %e\n", pInput_.inputLBounds_[idata]);
               exit(1);
            }
            if (pInput_.inputStdevs_[idata] != pInput_.inputLBounds_[idata])
            {
               printf("readInputSection ERROR: for uniform distributions, the\n");
               printf("    second data should be the same as the upper bound\n");
               printf("    in the variable definition.\n");
               printf("    upper bound defined  = %e\n", pInput_.inputStdevs_[idata]);
               printf("    upper bound expected = %e\n", pInput_.inputLBounds_[idata]);
               exit(1);
            }
         }
         else if ( !strcmp(winput2, "N")) 
         {
            pInput_.inputPDFs_[idata] = 1;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "L")) 
         {
            pInput_.inputPDFs_[idata] = 2;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "T")) 
         {
            pInput_.inputPDFs_[idata] = 3;
            sscanf(line,"%s %d %s %lg %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]),
                   &(pInput_.inputAuxs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "B")) 
         {
            pInput_.inputPDFs_[idata] = 4;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "W")) 
         {
            pInput_.inputPDFs_[idata] = 5;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "G")) 
         {
            pInput_.inputPDFs_[idata] = 6;
            sscanf(line,"%s %d %s %lg %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]),
                   &(pInput_.inputStdevs_[idata]));
            pInput_.useInputPDFs_ = 1;
         }
         else if ( !strcmp(winput2, "E")) 
         {
            pInput_.inputPDFs_[idata] = 7;
            sscanf(line,"%s %d %s %lg",winput,&itmp,winput2,
                   &(pInput_.inputMeans_[idata]));
            pInput_.inputStdevs_[idata] = 0.0;
            pInput_.useInputPDFs_ = 1;
         }
         else
         {
            printf("readInputSection ERROR: input format (2).\n");
            exit(1);
         } 
      } 
      else if (strcmp(winput, keywords[3]) == 0) /* correlation */
      {
         if (nInputs <= 0)
         {
            printf("readInputSection ERROR: nInputs not set.\n");
            exit(1);
         } 
         sscanf(line,"%s %d %d %lg", winput, &idata, &itmp, &ddata);
         idata--;
         itmp--;
         if (idata < 0 || idata >= pInput_.nInputs_ ||
             itmp < 0 || itmp >= pInput_.nInputs_)
         {
            printf("readInputSection ERROR: invalid input numbers (%d,%d)\n",
                   idata+1, itmp+1);
            printf("     input numbers should be between 1 and %d\n", nInputs);
            exit(1);
         }
         if (idata == itmp && ddata != 1.0)
         { 
            printf("readInputSection ERROR: Cor(%d,%d) should be = 1.\n",
                   idata+1, idata+1);
            exit(1);
         }
         if (idata != itmp && (ddata <= -1.0 || ddata >= 1.0))
         { 
            printf("readInputSection ERROR: |Cor(%d,%d)| should be < 1.\n",
                   idata+1, itmp+1);
            exit(1);
         }
         pInput_.corMatrix_.setEntry(idata, itmp, ddata);
         pInput_.corMatrix_.setEntry(itmp, idata, ddata);
      }
      else if (strcmp(winput, keywords[4]) == 0) /* NAME */
      {
         /* display name for the variable - not used by psuade */
      }
      else if (strcmp(winput, keywords[5]) == 0) /* END */
      {
         break;
      }
      else if (strcmp(winput,"#") == 0) 
      {
         /* comment line */
      }
      else if (winput[0] == '#') /* comments */
      {
         /* comment line */
      }
      else 
      {
         printf("readInputSection ERROR: unrecognized line: %s\n",line);
         exit(1);
      }
   }
   for (ii = 0; ii < pInput_.nInputs_; ii++)
   {
      if (pInput_.inputNames_[ii] == NULL)
      {
         printf("readInputSection ERROR: input %d not defined\n",ii+1);
         pInput_.inputNames_[ii] = new char[100];
         sprintf(pInput_.inputNames_[ii], "X%d", ii+1);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      printf("readInputSection ERROR: END not found.\n");
      exit(1);
   }
}

// ************************************************************************
// A function for reading PSUADE parameters from an input file.
// ------------------------------------------------------------------------ 
void PsuadeData::readOutputSection(FILE *fp) 
{
   int  ii, idata, itmp, lineLeng=200, nOutputs;
   char line[200], winput[200], winput2[200];
   const char *keywords[] = {"dimension", "variable", "NAME", "END"}; 

   if (fp == NULL)
   {
      printf("readOutputSection ERROR: file = null\n");
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      strcpy(winput, "#");
      sscanf(line,"%s", winput);
      if (strcmp(winput, keywords[0]) == 0) /* dimension */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2, &nOutputs);
         else sscanf(line,"%s %d", winput, &nOutputs);
         if (nOutputs <= 0)
         {
            printf("readOutputSection ERROR: nOutputs %d <= 0.\n", nOutputs);
            exit(1);
         } 
         if (pOutput_.nOutputs_ != 0 && pOutput_.nOutputs_ != nOutputs)
         {
            printf("readOutputSection ERROR: nOutputs mismatch.\n");
            printf("        check nOutputs in OUTPUT & PSUADE_IO sections.\n");
            exit(1);
         } 
         pOutput_.nOutputs_ = nOutputs;
         pOutput_.outputNames_ = new char*[pOutput_.nOutputs_];
         for (ii = 0; ii < pOutput_.nOutputs_; ii++)
            pOutput_.outputNames_[ii] = NULL;
      } 
      else if ( strcmp(winput, keywords[1]) == 0 ) /* variable */
      {
         if (pOutput_.nOutputs_ <= 0)
         {
            printf("readOutputSection ERROR: nOutputs not set.\n");
            exit(1);
         } 
         sscanf(line,"%s %d", winput, &idata);
         idata--;
         if (idata < 0 || idata >= pOutput_.nOutputs_)
         {
            printf("readOutputSection ERROR: invalid output no. %d.\n",
                   idata+1);
            printf("        Example format: variable 1 Y\n");
            exit(1);
         } 
         sscanf(line,"%s %d %s", winput, &itmp, winput2);
         if (pOutput_.outputNames_[idata] != NULL)
            delete [] pOutput_.outputNames_[idata];
         pOutput_.outputNames_[idata] = new char[strlen(winput2)+10];
         strncpy(pOutput_.outputNames_[idata], winput2, strlen(winput2)+1);
      }
      else if (strcmp(winput, keywords[2]) == 0) /* NAME */
      {
         /* display name for the variable - not used by psuade */
      }
      else if (strcmp(winput, keywords[3]) == 0) /* END */
      {
         break;
      }
      else if (strcmp(winput,"#") == 0) /* comments */
      {
         /* comment line */
      }
      else if (winput[0] == '#') /* comments */
      {
         /* comment line */
      }
      else 
      {
         printf("readOutputSection ERROR: unrecognized line - %s\n",
                line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      printf("readOutputSection ERROR: END not found.\n");
      exit(1);
   }
   for (ii = 0; ii < pOutput_.nOutputs_; ii++)
   {
      if (pOutput_.outputNames_[ii] == NULL)
      {
         printf("readOutputSection WARNING: variable %d undeclared.\n",
                ii+1);
         pOutput_.outputNames_[ii] = new char[100];
         sprintf(pOutput_.outputNames_[ii], "Y%d", ii+1);
      } 
   } 
}

// ************************************************************************
// A function for reading PSUADE parameters from an input file.
// ------------------------------------------------------------------------ 
void PsuadeData::readMethodSection(FILE *fp) 
{
   int  methodID, ii, idata, idata2, lineLeng=200, nSamples;
   char line[200], winput[200], winput2[200], winput3[200];
   const char *keywords[] = {"sampling", "randomize", "randomize_more",
                             "num_samples", "num_replications", 
                             "num_refinements", "reference_num_refinements",
                             "refinement_type", "input_settings",
                             "random_seed", "refinement_size", "END"}; 
   int  nMethods=30;
   const char *methods[] = {"MC","FACT","LH","OA","OALH","MOAT","SOBOL",
                            "LPTAU", "METIS","FAST","BBD","PBD","FF4","FF5",
                            "CCI4","CCI5", "CCIF","CCF4","CCF5","CCFF",
                            "CCC4","CCC5","CCCF", "SFAST","UMETIS","GMOAT",
                            "GMETIS","SPARSEGRID","DISCRETE", "LSA"};

   if (fp == NULL)
   {
      printf("readMethodSection ERROR: file = null\n");
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      strcpy(winput, "#");
      sscanf(line,"%s", winput);
      if (strcmp(winput, keywords[0]) == 0) /* sampling */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") != 0)
         {
            printf("readMethodSection ERROR: sampling format.\n");
            exit(1);
         } 
         sscanf(line,"%s %s %s", winput, winput2, winput3);
         for (methodID = 0; methodID < nMethods; methodID++)
         {
            if (strcmp(winput3, methods[methodID]) == 0)
            {
               pMethod_.samplingMethod_ = PSUADE_SAMP_MC + methodID;
               break;
            }
         }
         if (methodID >= nMethods) 
         {
            printf("readMethodSection ERROR: invalid method %s.\n",
                   winput3);
            exit(1);
         }
      }   
      else if (strcmp(winput, keywords[1]) == 0) /* randomize */
      {
         pMethod_.sampleRandomize_ = 1;
      }
      else if (strcmp(winput, keywords[2]) == 0) /* random perturbation */
      {
         pMethod_.sampleRandomize_ = 3;
      }
      else if (strcmp(winput, keywords[3]) == 0) /* number of samples */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2, &nSamples);
         else sscanf(line,"%s %d", winput, &nSamples);
         if (nSamples <= 0)
         {
            printf("readMethodSection ERROR: nSamples %d <= 0.\n",
                   nSamples);
            exit(1);
         }
         if (pMethod_.nSamples_ != 0 && pMethod_.nSamples_ != nSamples) 
         {
            printf("readMethodSection ERROR: nSamples mismatch.\n");
            printf("        check nSamples in METHOD & PSUADE_IO sections.\n");
            exit(1);
         } 
         pMethod_.nSamples_ = nSamples;
      }
      else if (strcmp(winput, keywords[4]) == 0) /* number of replications */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2,&pMethod_.nReplications_);
         else sscanf(line,"%s %d", winput, &pMethod_.nReplications_);
         if (pMethod_.nReplications_ <= 0)
         {
            printf("readMethodSection ERROR: nReplications <= 0.\n");
            exit(1);
         } 
      }
      else if (strcmp(winput, keywords[5]) == 0) /* number of refinements */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2,&pMethod_.nRefinements_);
         else sscanf(line,"%s %d", winput, &pMethod_.nRefinements_);
         if (pMethod_.nRefinements_ < 0) 
         {
            printf("readMethodSection ERROR: nRefinements < 0.\n");
            exit(1);
         } 
      }
      else if (strcmp(winput, keywords[6]) == 0) /* reference no. refinements */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d",winput,winput2,&pMethod_.refNumRefinements_);
         else sscanf(line,"%s %d", winput, &pMethod_.refNumRefinements_);
         if (pMethod_.refNumRefinements_ < 0) 
         {
            printf("readMethodSection ERROR: ref. nRefinements < 0.\n");
            exit(1);
         } 
      }
      else if (strcmp(winput, keywords[7]) == 0) /* refinement type */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         if (!strcmp(winput3, "adaptive")) pMethod_.refinementType_ = 1;
      }
      else if (strcmp(winput, keywords[8]) == 0) /* input settings */
      {
         sscanf(line,"%s %d %d", winput, &idata, &idata2);
         if (idata < 1 || idata > pInput_.nInputs_ ||
             pInput_.nInputs_ == 0 || idata2 <= 0)
         {
            printf("readMethodSection ERROR: invalid input %d\n",idata);
            printf(" in input_settings.\n");
            exit(1);
         } 
         if (pInput_.inputSettings_ == NULL)
         {
            pInput_.inputSettings_ = new double*[pInput_.nInputs_];
            for (ii = 0; ii < pInput_.nInputs_; ii++)
               pInput_.inputSettings_[ii] = NULL;
         }
         if (pInput_.inputNumSettings_ == NULL)
         {
            pInput_.inputNumSettings_ = new int[pInput_.nInputs_];
            for (ii = 0; ii < pInput_.nInputs_; ii++)
               pInput_.inputNumSettings_[ii] = 0;
         }
         pInput_.inputNumSettings_[idata-1] = idata2;
         pInput_.inputSettings_[idata-1] = new double[idata2];
         for (ii = 0; ii < idata2; ii++) 
            fscanf(fp, "%lg", &(pInput_.inputSettings_[idata-1][ii])); 
         fgets(line, lineLeng, fp);
      }
      else if (strcmp(winput, keywords[9]) == 0) /* random seed */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d", winput, winput2, &psRandomSeed_);
         else sscanf(line,"%s %d", winput, &psRandomSeed_);
         if (psRandomSeed_ <= 0) psRandomSeed_ = -1; 
      }
      else if (strcmp(winput, keywords[10]) == 0) /* refine size */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %d",winput,winput2,&pMethod_.refinementSize_);
         else sscanf(line,"%s %d", winput, &pMethod_.refinementSize_);
         if (pMethod_.refinementSize_ <= 0) pMethod_.refinementSize_ = 100000; 
      }
      else if (strcmp(winput, keywords[11]) == 0) /* END */
      {
         break;
      }
      else if (strcmp(winput,"#") == 0) /* comments */
      {
         /* comment line */
      }
      else if (winput[0] == '#') /* comments */
      {
         /* comment line */
      }
      else 
      {
         printf("readMethodSection ERROR: unrecognized line - %s\n",
                line);
         exit(1);
      }
   }
   if (pMethod_.samplingMethod_ == 1 && pMethod_.nReplications_ > 1)
   {
      printf("readMethodSection ERROR: \n");
      printf("\t\tFor Factorial sampling, nReplications should be 1.\n");
      exit(1);
   } 
   idata = pMethod_.nSamples_ / pMethod_.nReplications_ * 
           pMethod_.nReplications_;
   if (idata != pMethod_.nSamples_)
   {
      printf("readMethodSection ERROR: \n");
      printf("\t\tnSamples should be multiples of nReplications.\n");
      exit(1);
   } 
   if (feof(fp) != 0)
   {
      printf("readMethodSection ERROR: END not found.\n");
      exit(1);
   }
}

// ************************************************************************
// A function for reading PSUADE parameters from an input file.
// ------------------------------------------------------------------------ 
void PsuadeData::readApplicationSection(FILE *fp) 
{
   int  lineLeng=200;
   char line[200], winput[200], winput2[200], winput3[200];
   const char *keywords[] = {"driver","opt_driver","aux_opt_driver",
                             "output_template","input_template",
                             "max_parallel_jobs", "max_job_wait_time", 
                             "min_job_wait_time", "nondeterministic", 
                             "launch_interval", "save_frequency", 
                             "launch_only", "limited_launch_only",
                             "gen_inputfile_only", "function",
                             "use_function", "module", "launch_function",
		             "use_launch_script", "rs_driver", "use_rs",
                             "END"};
   FILE *fp2;

   if (fp == NULL)
   {
      printf("readApplicationSection ERROR: file = null\n");
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      strcpy(winput, "#");
      sscanf(line,"%s", winput);
      if (strcmp(winput, keywords[0]) == 0) /* driver */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.appDriver_, winput3, strlen(winput3)+1); 
         if (strcmp(pApplication_.appDriver_, "NONE"))
         {
            if (outputLevel_ > 1)
            {
               fp2 = fopen(pApplication_.appDriver_, "r");
               if (fp2 == NULL && outputLevel_ > 1)
               {
                  printf("readApplicationSection WARNING: ");
                  printf("app driver %s not found.\n",pApplication_.appDriver_);
               }
               else fclose(fp2);
            }
         }
      }   
      else if (strcmp(winput, keywords[1]) == 0) /* optimization driver */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.optDriver_, winput3, strlen(winput3)+1); 
         if (strcmp(pApplication_.optDriver_, "NONE"))
         {
            if (outputLevel_ > 1)
            {
               fp2 = fopen(pApplication_.optDriver_, "r");
               if (fp2 == NULL && outputLevel_ > 1)
               {
                  printf("readApplicationSection WARNING: ");
                  printf("opt driver %s not found.\n",pApplication_.optDriver_);
               }
               else fclose(fp2);
            }
         }
      }   
      else if (strcmp(winput, keywords[2]) == 0) /* aux opt driver */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.auxOptDriver_, winput3, strlen(winput3)+1); 
         if (strcmp(pApplication_.auxOptDriver_, "NONE"))
         {
            if (outputLevel_ > 1)
            {
               fp2 = fopen(pApplication_.auxOptDriver_, "r");
               if (fp2 == NULL && outputLevel_ > 1)
               {
                  printf("readApplicationSection WARNING: ");
                  printf("opt driver %s not found.\n",
                          pApplication_.auxOptDriver_);
               }
               else fclose(fp2);
            }
         }
      }   
      else if (strcmp(winput, keywords[3]) == 0) /* output template */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.outputTemplate_, winput3, strlen(winput3)+1); 
      }
      else if (strcmp(winput, keywords[4]) == 0) /* input template */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.inputTemplate_, winput3, strlen(winput3)+1); 
         if (strcmp(pApplication_.inputTemplate_, "NONE"))
         {
            if (outputLevel_ > 1)
            {
               fp2 = fopen(pApplication_.inputTemplate_, "r");
               if (fp2 == NULL)
               {
                  printf("readApplicationSection WARNING: ");
                  printf("input template %s not found.\n",
                         pApplication_.inputTemplate_);
               }
               else fclose(fp2);
            }
         }
      }
      else if (strcmp(winput, keywords[5]) == 0) /* max parallel jobs */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
              sscanf(line,"%s %s %d", winput, winput2,
                   &pApplication_.maxParallelJobs_);
         else sscanf(line,"%s %d", winput, &pApplication_.maxParallelJobs_);

         if (pApplication_.maxParallelJobs_ < 1)
             pApplication_.maxParallelJobs_ = 1;
         if (pApplication_.maxParallelJobs_ > 200) 
         {
            if (outputLevel_ > 1)
               printf("readApplicationSection: parallel jobs = 200.\n");
            pApplication_.maxParallelJobs_ = 200;
         }
      }
      else if (strcmp(winput, keywords[6]) == 0) /* maximum job wait time */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
              sscanf(line,"%s %s %d", winput, winput2,
                     &pApplication_.maxJobWaitTime_);
         else sscanf(line,"%s %d", winput, &pApplication_.maxJobWaitTime_);
      }
      else if (strcmp(winput, keywords[7]) == 0) /* minimum job wait time */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
              sscanf(line,"%s %s %d", winput, winput2,
                     &pApplication_.minJobWaitTime_);
         else sscanf(line,"%s %d", winput, &pApplication_.minJobWaitTime_);
      }
      else if (strcmp(winput, keywords[8]) == 0) /* nondeterministic */
      {
         pApplication_.runType_ |= 1;
         break;
      }
      else if (strcmp(winput, keywords[9]) == 0) /* launch interval */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
              sscanf(line,"%s %s %d", winput, winput2,
                     &pApplication_.launchInterval_);
         else sscanf(line,"%s %d", winput, &pApplication_.launchInterval_);
         if (pApplication_.launchInterval_ < 0)
            pApplication_.launchInterval_ = 0;
      }
      else if (strcmp(winput, keywords[10]) == 0) /* save frequency */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
              sscanf(line,"%s %s %d", winput, winput2,
                     &pApplication_.saveFrequency_);
         else sscanf(line,"%s %d", winput, &pApplication_.saveFrequency_);
         if (pApplication_.saveFrequency_ < 1)
            pApplication_.saveFrequency_ = 10;
      }
      else if (strcmp(winput, keywords[11]) == 0) /* launch only */
      {
         pApplication_.runType_ |= 2;
      }
      else if (strcmp(winput, keywords[12]) == 0) /* limited launch only */
      {
         pApplication_.runType_ |= 4;
      }
      else if (strcmp(winput, keywords[13]) == 0) /* generate input file only */
      {
         pApplication_.runType_ |= 8;
      }
      else if (strcmp(winput, keywords[14]) == 0) /* function */
      {
	 /* python driver function, not used by psuade */
      }
      else if (strcmp(winput, keywords[15]) == 0) /* use_function */
      {
         pApplication_.useFunction_ = 1;
      }
      else if (strcmp(winput, keywords[16]) == 0) /* module */
      {
	 /* python launch script module, not used by psuade */
      }
      else if (strcmp(winput, keywords[17]) == 0) /* launch_function */
      {
	 /* python launch script function name, not used by psuade */
      }
      else if (strcmp(winput, keywords[18]) == 0) /* use_launch_script */
      {
         pApplication_.useFunction_ = 2;
      }
      else if (strcmp(winput, keywords[19]) == 0) /* rs_driver */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0 )
            sscanf(line,"%s %s %s", winput, winput2, winput3);
         else sscanf(line,"%s %s", winput, winput3);
         strncpy(pApplication_.rsDriver_, winput3, strlen(winput3)+1); 
         if (strcmp(pApplication_.rsDriver_, "NONE"))
         {
            if (outputLevel_ > 1)
            {
               fp2 = fopen(pApplication_.rsDriver_, "r");
               if (fp2 == NULL && outputLevel_ > 1)
               {
                  printf("readApplicationSection WARNING: ");
                  printf("rs driver %s not found.\n",pApplication_.rsDriver_);
               }
               else fclose(fp2);
            }
         }
      }   
      else if (strcmp(winput, keywords[20]) == 0) /* use_rs */
      {
         pApplication_.useFunction_ = 3;
      }
      else if (strcmp(winput, keywords[21]) == 0) /* END */
      {
         break;
      }
      else if (strcmp(winput,"#") == 0) 
      {
         /* comment line */
      }
      else if (winput[0] == '#') /* comments */
      {
         /* comment line */
      }
      else 
      {
         printf("readApplicationSection ERROR: \n");
         printf("\t\t unrecognized line - %s\n", line);
         exit(1);
      }
   }

   if ( pApplication_.useFunction_ == 1 )
      strcpy( pApplication_.appDriver_, "psuadeFunctionDriver" );
   else if ( pApplication_.useFunction_ == 2 )
      strcpy( pApplication_.appDriver_, "psuadeLaunchDriver" );
   else if ( pApplication_.useFunction_ == 3 )
      strcpy( pApplication_.appDriver_, pApplication_.rsDriver_ );

   if ( pApplication_.minJobWaitTime_ > 0 &&
        pApplication_.minJobWaitTime_ >= pApplication_.maxJobWaitTime_ )
   {
      pApplication_.maxJobWaitTime_ = pApplication_.minJobWaitTime_ * 5;
      if (outputLevel_ > 1)
      {
         printf("readApplicationSection INFO: max job wait time\n");
         printf("            has been reset to 5*(min job wait time)\n");
      }
   }
   if (feof(fp) != 0)
   {
      printf("readApplicationSection ERROR: END not found.\n");
      exit(1);
   }
}

// ************************************************************************
// A function for reading PSUADE analysis parameters from an input file.
// ------------------------------------------------------------------------ 
void PsuadeData::readAnalysisSection(FILE *fp) 
{
   int    lineLeng=200, outputID, rstype, transform, ii;
   double threshold, lbound, ubound;
   char   line[200], winput[200], winput2[200], winput3[200], winput4[200];
   char   winput5[200];
   const char *keywords[] = {
              "analyzer", "optimization", "graphics", "sample_graphics",
              "file_write", "diagnostics", "config_file", "interactive", 
              "ana_expert", "rs_expert", "opt_expert", "sam_expert", 
              "expert", "io_expert", "rs_max_pts", "scilab", "use_input_pdfs",
              "constraint_op_and", "END"};
   const char *analysisOptions[] = {
              "method", "threshold", "output_id", "rstype", "transform", 
              "regression_wgt_id", "rs_constraint", "moat_constraint", 
              "rs_index_file", "rs_legendre_order"};
   const char   *deprecateOptions[] = {
              "rsfilter", "moat_filter", "fileWrite", "rsMaxPts", 
              "sampleGraphics"};
   const char *analysisMethods[] = {
              "Moment", "MainEffect", "TwoParamEffect", "ANOVA", "GLSA", 
              "RSFA","MOAT", "Sobol", "Correlation", "Integration",
              "FAST", "FF", "PCA", "ARSM0", "FORM", "RSMSobol1",
              "RSMSobol2", "RSMSobolTSI", "Bootstrap", "RSMSobolG", "ARSMNN",
              "ARSM1", "REL", "AOPT", "GOWER", "DELTA", "ETA", "ARSMG", "LSA"};
   const char *resSurfTypes[] = {
              "MARS","linear","quadratic","cubic", "quartic", "ANN", 
              "selective_regression", "GP1", "GP2", "SVM", "PWL", "TGP",
              "MARSBag", "EARTH", "sum_of_trees", "Legendre", "user_regression",
              "sparse_grid_regression"};
   const char *transformTypes[] = {"logx","logy"};
   const char *optimizeOptions[] = {
              "method", "fmin", "num_local_minima", "use_response_surface", 
              "cutoff", "tolerance", "print_level", "num_fmin", "output_id",
              "max_feval", "deltax", "target_file"};
   const char *optimizeSchemes[] = {
              "crude", "txmath", "appspack", "minpack", "cobyla", "sm", "mm",
              "mm_adaptive", "bobyqa", "sce", "moo"};
   psuadeFilter **filters;
   FILE   *fconf;

   if (fp == NULL)
   {
      printf("readAnalysisSection ERROR: file = null\n");
      exit(1);
   }
   while ((fgets(line, lineLeng, fp) != NULL) && (feof(fp) == 0))
   {
      strcpy(winput, "#");
      sscanf(line,"%s", winput);

      if (strcmp(winput, keywords[0]) == 0)  /* analyzer */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (!strcmp(winput2,analysisOptions[0])) /* == analysis method == */
         {
            sscanf(line,"%s %s %s %s", winput, winput2, winput3, winput4);
            if (!strcmp(winput4,analysisMethods[0]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_MOMENT;
            else if (!strcmp(winput4,analysisMethods[1]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ME;
            else if (!strcmp(winput4,analysisMethods[2]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_IE;
            else if (!strcmp(winput4,analysisMethods[3]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ANOVA;
            else if (!strcmp(winput4,analysisMethods[4]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_GLSA;
            else if (!strcmp(winput4,analysisMethods[5]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_RSFA;
            else if (!strcmp(winput4,analysisMethods[6]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_MOAT;
            else if (!strcmp(winput4,analysisMethods[7]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_SOBOL;
            else if (!strcmp(winput4,analysisMethods[8]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_CORRELATION;
            else if (!strcmp(winput4,analysisMethods[9]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_INTEGRATION;
            else if (!strcmp(winput4,analysisMethods[10]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_FAST;
            else if (!strcmp(winput4,analysisMethods[11]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_FF;
            else if (!strcmp(winput4,analysisMethods[12]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_PCA;
            else if (!strcmp(winput4,analysisMethods[13]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ONESIGMA;
            else if (!strcmp(winput4,analysisMethods[14]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_FORM;
            else if (!strcmp(winput4,analysisMethods[15]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_RSSOBOL1;
            else if (!strcmp(winput4,analysisMethods[16]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_RSSOBOL2;
            else if (!strcmp(winput4,analysisMethods[17]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_RSSOBOLTSI;
            else if (!strcmp(winput4,analysisMethods[18]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_BSTRAP;
            else if (!strcmp(winput4,analysisMethods[19]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_RSSOBOLG;
            else if (!strcmp(winput4,analysisMethods[20]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ARSM;
            else if (!strcmp(winput4,analysisMethods[21]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ARSMMB;
            else if (!strcmp(winput4,analysisMethods[22]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_REL;
            else if (!strcmp(winput4,analysisMethods[23]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_AOPT;
            else if (!strcmp(winput4,analysisMethods[24]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_GOWER;
            else if (!strcmp(winput4,analysisMethods[25]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_DTEST;
            else if (!strcmp(winput4,analysisMethods[26]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ETEST;
            else if (!strcmp(winput4,analysisMethods[27]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_ARSMMBBS;
            else if (!strcmp(winput4,analysisMethods[28]))
               pAnalysis_.analysisIntOptions_[0] |= PSUADE_ANA_LSA;
            else
            {
               printf("readAnalysisSection ERROR: invalid method %s\n",
                      winput4);
               exit(1);
            }
         }
         else if (!strcmp(winput2,analysisOptions[1]))/* ==analysis thresh== */
         {
            sscanf(line,"%s %s %s %lg", winput, winput2, winput3, 
                   &threshold);
            if (threshold < 0.0 || threshold > 1.0)
            {
               printf("readAnalysisSection ERROR: invalid threshold %e.\n",
                      threshold);
               printf("           (threshold has to be > 0.0 and < 1.0.\n");
               exit(1);
            }
            pAnalysis_.analysisDbleOptions_[0] = threshold;
         }
         else if (!strcmp(winput2,analysisOptions[2]))/* == analysis outID == */
         {
            sscanf(line,"%s %s %s %d", winput, winput2, winput3, &outputID);
            if (outputID <= 0 || outputID > pOutput_.nOutputs_)
            {
               printf("readAnalysisSection ERROR: invalid outputID %d\n",
                      outputID);
               printf("           outputID has to be between 1 and %d\n",
                      pOutput_.nOutputs_);
               exit(1);
            }
            pAnalysis_.analysisIntOptions_[1] = outputID - 1;
         }
         else if (!strcmp(winput2,analysisOptions[3]))/* ==analysis rs type== */
         {
            sscanf(line,"%s %s %s %s", winput, winput2, winput3, winput4);
            rstype = 0;
            if      (!strcmp(winput4,resSurfTypes[0])) rstype = PSUADE_RS_MARS;
            else if (!strcmp(winput4,resSurfTypes[1])) rstype = PSUADE_RS_REGR1;
            else if (!strcmp(winput4,resSurfTypes[2])) rstype = PSUADE_RS_REGR2;
            else if (!strcmp(winput4,resSurfTypes[3])) rstype = PSUADE_RS_REGR3;
            else if (!strcmp(winput4,resSurfTypes[4])) rstype = PSUADE_RS_REGR4;
            else if (!strcmp(winput4,resSurfTypes[5])) rstype = PSUADE_RS_ANN;
            else if (!strcmp(winput4,resSurfTypes[6])) rstype = PSUADE_RS_REGRU;
            else if (!strcmp(winput4,resSurfTypes[7])) rstype = PSUADE_RS_GP1;
            else if (!strcmp(winput4,resSurfTypes[8])) rstype = PSUADE_RS_GP2;
            else if (!strcmp(winput4,resSurfTypes[9])) rstype = PSUADE_RS_SVM;
            else if (!strcmp(winput4,resSurfTypes[10])) rstype = PSUADE_RS_PWL;
            else if (!strcmp(winput4,resSurfTypes[11])) rstype = PSUADE_RS_TGP;
            else if (!strcmp(winput4,resSurfTypes[12])) rstype = PSUADE_RS_MARSB;
            else if (!strcmp(winput4,resSurfTypes[13])) rstype = PSUADE_RS_EARTH;
            else if (!strcmp(winput4,resSurfTypes[14])) rstype = PSUADE_RS_SOTS;
            else if (!strcmp(winput4,resSurfTypes[15])) rstype = PSUADE_RS_REGRL;
            else if (!strcmp(winput4,resSurfTypes[16])) rstype = PSUADE_RS_REGRU;
            else if (!strcmp(winput4,resSurfTypes[17])) rstype = PSUADE_RS_REGSG;
            else
            {
               printf("readAnalysisSection ERROR: invalid RS type %s\n",
                      winput4);
               exit(1);
            }
            pAnalysis_.analysisIntOptions_[2] = rstype;
         }
         else if (!strcmp(winput2,analysisOptions[4]))/* ==analysis xsform == */
         {
            sscanf(line,"%s %s %s %s", winput, winput2, winput3, winput4);
            transform = 0;
            if      (!strcmp(winput4,transformTypes[0])) transform |= 1;
            else if (!strcmp(winput4,transformTypes[1])) transform |= 2;
            else
            {
               printf("readAnalysisSection ERROR: invalid transform %s\n",
                      winput4);
               exit(1);
            }
            pAnalysis_.analysisIntOptions_[5] |= transform;
         }
         else if (!strcmp(winput2,analysisOptions[5]))/*=analysis regr wgt =*/
         {
            sscanf(line,"%s %s %s %d",winput,winput2,winput3,
                   &pAnalysis_.analysisIntOptions_[6]);
            pAnalysis_.analysisIntOptions_[6]--;
            printf("readAnalysisSection INFO: regression wgt outputID = %d.\n",
                   pAnalysis_.analysisIntOptions_[6]+1);
            printf("        Please make sure it is correct.\n");
         }
         else if (strcmp(winput2, analysisOptions[6]) == 0 ||
                  strcmp(winput2, deprecateOptions[0]) == 0) /* rs_constraint */
         {
            filters = pAnalysis_.RSFilters_;
            ii = pAnalysis_.numRSFilters_ + 1;
            pAnalysis_.RSFilters_ = new psuadeFilter*[ii];
            for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
               pAnalysis_.RSFilters_[ii] = filters[ii];
            if (filters != NULL) delete [] filters;
            pAnalysis_.numRSFilters_++;
            ii = pAnalysis_.numRSFilters_ - 1;
            pAnalysis_.RSFilters_[ii] = new psuadeFilter();
            sscanf(line,"%s %s %s %s %s %lg %lg", winput, winput2, winput3,
                   winput4, winput5, &lbound, &ubound);
            pAnalysis_.RSFilters_[ii]->FilterDataFile_ = new char[500];
            pAnalysis_.RSFilters_[ii]->FilterIndexFile_ = new char[500];
            strcpy(pAnalysis_.RSFilters_[ii]->FilterDataFile_, winput4);
            strcpy(pAnalysis_.RSFilters_[ii]->FilterIndexFile_, winput5);
            pAnalysis_.RSFilters_[ii]->FilterLBound_ = lbound;
            pAnalysis_.RSFilters_[ii]->FilterUBound_ = ubound;
            if (winput3[0] != '=')
            {
               printf("analyzer rs_contraint format ERROR: \n");
               printf("   format:  analyzer rs_contraint = ");
               printf("dataFile indexFile lbound ubound\n");
               exit(1);
            }
         }
         else if (strcmp(winput2, analysisOptions[7]) == 0 || 
                  strcmp(winput2, deprecateOptions[1]) == 0) /* moat_constraint */
         {
            filters = pAnalysis_.MOATFilters_;
            ii = pAnalysis_.numMOATFilters_ + 1;
            pAnalysis_.MOATFilters_ = new psuadeFilter*[ii];
            for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
               pAnalysis_.MOATFilters_[ii] = filters[ii];
            if (filters != NULL) delete [] filters;
            pAnalysis_.numMOATFilters_++;
            ii = pAnalysis_.numMOATFilters_ - 1;
            pAnalysis_.MOATFilters_[ii] = new psuadeFilter();
            sscanf(line,"%s %s %s %s %s %lg %lg", winput, winput2, winput3,
                   winput4, winput5, &lbound, &ubound);
            pAnalysis_.MOATFilters_[ii]->FilterDataFile_ = new char[500];
            pAnalysis_.MOATFilters_[ii]->FilterIndexFile_ = new char[500];
            strcpy(pAnalysis_.MOATFilters_[ii]->FilterDataFile_, winput4);
            strcpy(pAnalysis_.MOATFilters_[ii]->FilterIndexFile_, winput5);
            pAnalysis_.MOATFilters_[ii]->FilterLBound_ = lbound;
            pAnalysis_.MOATFilters_[ii]->FilterUBound_ = ubound;
            if (winput3[0] != '=')
            {
               printf("analyzer moat_constraint format ERROR: \n");
               printf("   format:  analyzer moat_constraint = ");
               printf("dataFile indexFile lbound ubound\n");
               exit(1);
            }
         }
         else if (strcmp(winput2, analysisOptions[8]) == 0) /* rs index file */
         {
            sscanf(line,"%s %s %s %s", winput, winput2, winput3, 
                   pAnalysis_.rsIndexFile_);
         }
         else if (strcmp(winput2, analysisOptions[9]) == 0) /* rs_legendre_order */
         {
            sscanf(line,"%s %s %s %d", winput, winput2, winput3, 
                   &pAnalysis_.legendreOrder_);
         }
         else 
         {
            printf("readAnalysisSection ERROR: \n");
            printf("unrecognized analyzer option - %s\n", winput2);
            exit(1);
         }
      }
      else if (strcmp(winput, keywords[1]) == 0)  /* optimization */
      {
         sscanf(line,"%s %s", winput, winput2);
        
         if (!strcmp(winput2, optimizeOptions[0])) /* opt method */
         {
            sscanf(line,"%s %s %s %s",winput,winput2,winput3,winput4);
            if (!strcmp(winput4, optimizeSchemes[0])) /* crude method */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 0;
            }
            else if (!strcmp(winput4, optimizeSchemes[1])) /* txmath */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 1;
            }
            else if (!strcmp(winput4, optimizeSchemes[2])) /* appspack */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 2;
            }
            else if (!strcmp(winput4, optimizeSchemes[3])) /* minpack */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 3;
            }
            else if (!strcmp(winput4, optimizeSchemes[4])) /* cobyla */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 4;
            }
            else if (!strcmp(winput4, optimizeSchemes[5])) /* sm */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 5;
            }
            else if (!strcmp(winput4, optimizeSchemes[6])) /* mm */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 6;
            }
            else if (!strcmp(winput4, optimizeSchemes[7])) /* adaptive mm */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 7;
            }
            else if (!strcmp(winput4, optimizeSchemes[8])) /* bobyqa */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 8;
            }
            else if (!strcmp(winput4, optimizeSchemes[9])) /* SCE */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 9;
            }
            else if (!strcmp(winput4, optimizeSchemes[10])) /* MOO */
            {
               pAnalysis_.optimizeIntOptions_[0] = 1;
               pAnalysis_.optimizeIntOptions_[1] = 10;
            }
            else
            {
               printf("readAnalysisSection ERROR: invalid opt scheme %s\n",
                      winput4);
               exit(1);
            }
         }
         else if (!strcmp(winput2, optimizeOptions[1])) /* opt fmin */
         {
            sscanf(line,"%s %s %s %lg",winput,winput2,winput3,
                   &(pAnalysis_.optimizeDbleOptions_[0]));
         }
         else if (!strcmp(winput2, optimizeOptions[2])) /* opt numpts */
         {
            sscanf(line,"%s %s %s %d",winput,winput2,winput3,
                   &(pAnalysis_.optimizeIntOptions_[2]));
            if (pAnalysis_.optimizeIntOptions_[2] < 1)
            {
               printf("readAnalysisSection ERROR: opt nPts < 1.\n");
               exit(1);
            }
         }
         else if (!strcmp(winput2, optimizeOptions[3])) /* opt RS */
         {
            pAnalysis_.optimizeIntOptions_[3] = 1;
         }
         else if (!strcmp(winput2, optimizeOptions[4])) /* opt cutoff */
         {
            sscanf(line,"%s %s %s %lg",winput,winput2,winput3,
                   &(pAnalysis_.optimizeDbleOptions_[1]));
         }
         else if (!strcmp(winput2, optimizeOptions[5])) /* opt tolerance */
         {
            sscanf(line,"%s %s %s %lg",winput,winput2,winput3,
                   &(pAnalysis_.optimizeDbleOptions_[2]));
         }
         else if (!strcmp(winput2, optimizeOptions[6])) /* opt outlevel */
         {
            sscanf(line,"%s %s %s %d",winput,winput2,winput3,
                   &(pAnalysis_.optimizeIntOptions_[4]));
         }
         else if (!strcmp(winput2, optimizeOptions[7])) /* opt num_fmin */
         {
            sscanf(line,"%s %s %s %d",winput,winput2,winput3,
                   &(pAnalysis_.optimizeIntOptions_[5]));
            if (pAnalysis_.optimizeIntOptions_[5] <= 0)
            {
               printf("readAnalysisSection ERROR: opt nfmin <= 0.\n");
               exit(1);
            }
         }
         else if (!strcmp(winput2, optimizeOptions[8])) /* opt output ID */
         {
            sscanf(line,"%s %s %s %d", winput, winput2, winput3, &outputID);
            if (outputID <= 0 || outputID > pOutput_.nOutputs_)
            {
               printf("readAnalysisSection ERROR: invalid outputID %d\n",
                      outputID);
               exit(1);
            }
            pAnalysis_.optimizeIntOptions_[6] = outputID - 1;
         }
         else if (!strcmp(winput2, optimizeOptions[9])) /* opt max feval */
         {
            sscanf(line,"%s %s %s %d", winput, winput2, winput3, 
                   &(pAnalysis_.optimizeIntOptions_[7]));
         }
         else if (!strcmp(winput2, optimizeOptions[10])) /* opt delta X */
         {
            sscanf(line,"%s %s %s %lg", winput, winput2, winput3, 
                   &(pAnalysis_.optimizeDbleOptions_[3]));
         }
         else if (!strcmp(winput2, optimizeOptions[11])) /* opt TARGET FILE */
         {
            sscanf(line,"%s %s %s %s", winput, winput2, winput3, 
                   pAnalysis_.specsFile_);
         }
         else
         {
            printf("readAnalysisSection ERROR: unrecognized optimize options %s.\n",
                   winput2);
            exit(1);
         }
      }
      else if (strcmp(winput, keywords[2]) == 0) /* graphics */
      {
         pAnalysis_.analysisIntOptions_[4] |= 1;
      }
      else if (strcmp(winput, keywords[3]) == 0 ||
               strcmp(winput, deprecateOptions[4]) == 0) /* sample Graphics */
      {
         pAnalysis_.analysisIntOptions_[4] |= 2;
      }
      else if (strcmp(winput, keywords[4]) == 0 ||
               strcmp(winput, deprecateOptions[2]) == 0) /* file write */
      {
         sscanf(line,"%s %s", winput, winput2);
         if (!strcmp(winput2, "matlab")) pAnalysis_.fileWriteFlag_ |= 1;
      }
      else if (strcmp(winput, keywords[5]) == 0) /* diagnostics */
      {
         sscanf(line,"%s %d", winput, &(pAnalysis_.analysisIntOptions_[3]));
         if (pAnalysis_.analysisIntOptions_[3] < 0)
            pAnalysis_.analysisIntOptions_[3] = 0;
         else if (pAnalysis_.analysisIntOptions_[3] > 5)
            pAnalysis_.analysisIntOptions_[3] = 5;
         outputLevel_ = pAnalysis_.analysisIntOptions_[3];
      }
      else if (!strcmp(winput, keywords[6])) /* configure file */
      {
         sscanf(line,"%s %s %s", winput, winput2, winput3); 
         if (strcmp(winput3, "NONE"))
         {
            fconf = fopen(winput3,"r");
            if (fconf == NULL)
            {
               printf("readAnalysisSection ERROR: config file %s not found.\n",
                      winput3);
            }
            else
            {
               fclose(fconf); 
               if (psConfigFileName_ == NULL) 
                  psConfigFileName_ = new char[500];
               strcpy(psConfigFileName_, winput3);
               psConfig_ = new PsuadeConfig(winput3,1);
            }
         }
      }
      else if (strcmp(winput, keywords[7]) == 0) /* interactive */
      {
         psAnalysisInteractive_ = 1;
      }
      else if (strcmp(winput, keywords[8]) == 0) /* analysis expert mode */
      {
         psAnaExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[9]) == 0)/*response surface expert mode*/
      {
         psRSExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[10]) == 0)/*optimization expert mode*/
      {
         psOptExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[11]) == 0)/*sampling expert mode*/
      {
         psSamExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[12]) == 0)/* expert mode*/
      {
         psExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[13]) == 0)/* IO expert mode*/
      {
         psIOExpertMode_ = 1;
      }
      else if (strcmp(winput, keywords[14]) == 0 ||
               strcmp(winput, deprecateOptions[3]) == 0) /* RS max points */
      {
         sscanf(line,"%s %s %d", winput, winput2, &psFAMaxDataPts_); 
         if (psFAMaxDataPts_ <= 1 || psFAMaxDataPts_ > 10000)
            psFAMaxDataPts_ = 4000;
         printf("Max number of data points for response surface is set to %d\n",
                psFAMaxDataPts_);
      }
      else if (strcmp(winput, keywords[15]) == 0) /* scilab mode */
      {
         psPlotTool_ = 1;
         break;
      }
      else if (strcmp(winput, keywords[16]) == 0) /* use_input_pdfs */
      {
         pAnalysis_.useInputPDFs_ = 1;
         printf("The use_input_pdfs option has been disabled.\n");
         printf("Instead, whenever PDFs are specified in the INPUT section,\n");
         printf("they will be used in analysis if relevant. To not use PDFs\n");
         printf("in analysis, you will have to comment out the PDF definitions.\n");
         break;
      }
      else if (strcmp(winput, keywords[17]) == 0) /* constraint_op_and */
      {
         psConstraintSetOp_ = 1;
         break;
      }
      else if (strcmp(winput, keywords[18]) == 0) /* END */
      {
         break;
      }
      else if (strcmp(winput,"#") == 0) 
      {
         /* comment line */
      }
      else if (winput[0] == '#') /* comments */
      {
         /* comment line */
      }
      else 
      {
         printf("readAnalysisSection ERROR: \n");
         printf("\t\t unrecognized line - %s\n", line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      printf("readAnalysisSection ERROR: END not found.\n");
      exit(1);
   }
}

// ************************************************************************
// A function for writing PSUADE data to an output file.
// ------------------------------------------------------------------------ 
void PsuadeData::writePsuadeIO(FILE *fOut, int flag) 
{
   int ss, ii;

   if (pMethod_.nSamples_ > 100000 && psIOExpertMode_ == 0)
   {
      printf("INFO: Data too large to be written (%d).\n", pMethod_.nSamples_);
      printf("      To enforce writing the large data set, use io_expert.\n");
      return;
   }
   fprintf(fOut,"PSUADE_IO (Note : inputs not true inputs if pdf ~=U)\n");
   fprintf(fOut,"%d %d %d\n", pInput_.nInputs_, pOutput_.nOutputs_, 
           pMethod_.nSamples_);
   for (ss = 0; ss < pMethod_.nSamples_; ss++)
   {
      if (flag == 0)
           fprintf(fOut,"%d %d\n", ss+1, pOutput_.sampleStates_[ss]);
      else fprintf(fOut,"%d 0\n", ss+1);
      for (ii = 0; ii < pInput_.nInputs_; ii++)
         fprintf(fOut,"%24.16e\n", 
              pInput_.sampleInputs_[ss*pInput_.nInputs_+ii]);
      for (ii = 0; ii < pOutput_.nOutputs_; ii++)
         fprintf(fOut,"%24.16e\n",
              pOutput_.sampleOutputs_[ss*pOutput_.nOutputs_+ii]);
   }
   fprintf(fOut,"PSUADE_IO\n");
}

// ************************************************************************
// A function for writing input section to a file
// ------------------------------------------------------------------------ 
void PsuadeData::writeInputSection(FILE *fOut) 
{
   int    ii, jj, testFlag;
   double ddata;

   fprintf(fOut,"INPUT\n");
   fprintf(fOut,"   dimension = %d\n", pInput_.nInputs_);
   fflush(fOut);
   for (ii = 0; ii < pInput_.nInputs_; ii++)
   {
      fprintf(fOut,"   variable %d %s ",ii+1, pInput_.inputNames_[ii]);
      fprintf(fOut," = %24.16e %24.16e\n",pInput_.inputLBounds_[ii],
              pInput_.inputUBounds_[ii]);
   }
   testFlag = 1;
   for (ii = 0; ii < pInput_.nInputs_; ii++)
   {
      if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 1)
      {
         fprintf(fOut,"   PDF %d N %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 2)
      {
         fprintf(fOut,"   PDF %d L %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 3)
      {
         fprintf(fOut,"   PDF %d T %12.5e %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii],
                 pInput_.inputAuxs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 4)
      {
         fprintf(fOut,"   PDF %d B %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 5)
      {
         fprintf(fOut,"   PDF %d W %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 6)
      {
         fprintf(fOut,"   PDF %d G %12.5e %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii], pInput_.inputStdevs_[ii]);
         testFlag = 0;
      }
      else if (pInput_.inputPDFs_ != NULL && pInput_.inputPDFs_[ii] == 7)
      {
         fprintf(fOut,"   PDF %d E %12.5e\n", ii+1,
                 pInput_.inputMeans_[ii]);
         testFlag = 0;
      }
      else
      {
      }
   }
   for (ii = 0; ii < pInput_.nInputs_; ii++)
   {
      for (jj = ii+1; jj < pInput_.nInputs_; jj++)
      {
         ddata = pInput_.corMatrix_.getEntry(ii, jj);
         if (ddata != 0.0)
            fprintf(fOut,"   COR %d %d %12.5e\n", ii+1, jj+1, ddata);
      }
   }
   fprintf(fOut,"END\n");
}

// ************************************************************************
// A function for writing output section to a file
// ------------------------------------------------------------------------ 
void PsuadeData::writeOutputSection(FILE *fOut) 
{
   int ii;

   fprintf(fOut,"OUTPUT\n");
   fprintf(fOut,"   dimension = %d\n", pOutput_.nOutputs_);
   for (ii = 0; ii < pOutput_.nOutputs_; ii++)
      fprintf(fOut,"   variable %d %s\n", ii+1, pOutput_.outputNames_[ii]);
   fprintf(fOut,"END\n");
}

// ************************************************************************
// A function for writing method section to a file
// ------------------------------------------------------------------------ 
void PsuadeData::writeMethodSection(FILE *fOut) 
{
   int ii, kk;

   fprintf(fOut,"METHOD\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_MC)
        fprintf(fOut,"   sampling = MC\n");
   else fprintf(fOut,"#  sampling = MC\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_FACT)
        fprintf(fOut,"   sampling = FACT\n");
   else fprintf(fOut,"#  sampling = FACT\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_LHS)
        fprintf(fOut,"   sampling = LH\n");
   else fprintf(fOut,"#  sampling = LH\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_OA)
        fprintf(fOut,"   sampling = OA\n");
   else fprintf(fOut,"#  sampling = OA\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_OALH)
        fprintf(fOut,"   sampling = OALH\n");
   else fprintf(fOut,"#  sampling = OALH\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_MOAT)
        fprintf(fOut,"   sampling = MOAT\n");
   else fprintf(fOut,"#  sampling = MOAT\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_SOBOL)
        fprintf(fOut,"   sampling = SOBOL\n");
   else fprintf(fOut,"#  sampling = SOBOL\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_LPTAU)
        fprintf(fOut,"   sampling = LPTAU\n");
   else fprintf(fOut,"#  sampling = LPTAU\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_METIS)
        fprintf(fOut,"   sampling = METIS\n");
   else fprintf(fOut,"#  sampling = METIS\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_FAST)
        fprintf(fOut,"   sampling = FAST\n");
   else fprintf(fOut,"#  sampling = FAST\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_BBD)
        fprintf(fOut,"   sampling = BBD\n");
   else fprintf(fOut,"#  sampling = BBD\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_PBD)
        fprintf(fOut,"   sampling = PBD\n");
   else fprintf(fOut,"#  sampling = PBD\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_FF4)
        fprintf(fOut,"   sampling = FF4\n");
   else fprintf(fOut,"#  sampling = FF4\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_FF5)
        fprintf(fOut,"   sampling = FF5\n");
   else fprintf(fOut,"#  sampling = FF5\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCI4)
        fprintf(fOut,"   sampling = CCI4\n");
   else fprintf(fOut,"#  sampling = CCI4\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCI5)
        fprintf(fOut,"   sampling = CCI5\n");
   else fprintf(fOut,"#  sampling = CCI5\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCIF)
        fprintf(fOut,"   sampling = CCIF\n");
   else fprintf(fOut,"#  sampling = CCIF\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCF4)
        fprintf(fOut,"   sampling = CCF4\n");
   else fprintf(fOut,"#  sampling = CCF4\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCF5)
        fprintf(fOut,"   sampling = CCF5\n");
   else fprintf(fOut,"#  sampling = CCF5\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCFF)
        fprintf(fOut,"   sampling = CCFF\n");
   else fprintf(fOut,"#  sampling = CCFF\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCC4)
        fprintf(fOut,"   sampling = CCC4\n");
   else fprintf(fOut,"#  sampling = CCC4\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCC5)
        fprintf(fOut,"   sampling = CCC5\n");
   else fprintf(fOut,"#  sampling = CCC5\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_CCCF)
        fprintf(fOut,"   sampling = CCCF\n");
   else fprintf(fOut,"#  sampling = CCCF\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_SFAST)
        fprintf(fOut,"   sampling = SFAST\n");
   else fprintf(fOut,"#  sampling = SFAST\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_UMETIS)
        fprintf(fOut,"   sampling = UMETIS\n");
   else fprintf(fOut,"#  sampling = UMETIS\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_GMOAT)
        fprintf(fOut,"   sampling = GMOAT\n");
   else fprintf(fOut,"#  sampling = GMOAT\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_GMETIS)
        fprintf(fOut,"   sampling = GMETIS\n");
   else fprintf(fOut,"#  sampling = GMETIS\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_SG)
        fprintf(fOut,"   sampling = SPARSEGRID\n");
   else fprintf(fOut,"#  sampling = SPARSEGRID\n");
   if (pMethod_.samplingMethod_ == PSUADE_SAMP_LSA)
        fprintf(fOut,"   sampling = LSA\n");
   else fprintf(fOut,"#  sampling = LSA\n");
   fprintf(fOut,"   num_samples = %d\n", pMethod_.nSamples_);
   fprintf(fOut,"   num_replications = %d\n", pMethod_.nReplications_);

   fprintf(fOut,"   num_refinements = %d\n", pMethod_.nRefinements_);
   fprintf(fOut,"   refinement_size = %d\n", pMethod_.refinementSize_);
   fprintf(fOut,"   reference_num_refinements = %d\n", 
           pMethod_.refNumRefinements_);
   if (pMethod_.refinementType_ == 1)
      fprintf(fOut,"   refinement_type = adaptive\n"); 
   else
      fprintf(fOut,"#  refinement_type = adaptive\n"); 
   if (pInput_.symbolTable_ != NULL) 
   {
      for (ii = 0; ii < pInput_.nInputs_; ii++)
         fprintf(fOut,"   num_symbols %d = %d\n",ii, 
                 pInput_.symbolTable_[ii]);
   }
   if (pMethod_.sampleRandomize_ == 0)
   {
      fprintf(fOut,"#  randomize\n");
      fprintf(fOut,"#  randomize_more\n");
   }
   else if (pMethod_.sampleRandomize_ == 1)
   {
      fprintf(fOut,"   randomize\n");
      fprintf(fOut,"#  randomize_more\n");
   }
   else if (pMethod_.sampleRandomize_ == 3)
   {
      fprintf(fOut,"   randomize\n");
      fprintf(fOut,"   randomize_more\n");
   }
   if (pInput_.inputNumSettings_ != NULL) 
   {
      for (ii = 0; ii < pInput_.nInputs_; ii++)
      {
         if (pInput_.inputNumSettings_[ii] > 0) 
         {
            fprintf(fOut,"   input_settings %4d %4d\n", ii+1,
                    pInput_.inputNumSettings_[ii]);
            for (kk = 0; kk < pInput_.inputNumSettings_[ii]; kk++)
               fprintf(fOut,"                  %12.4e\n", 
                       pInput_.inputSettings_[ii][kk]);
         }
      }
   }
   if (psRandomSeed_ != -1)
        fprintf(fOut,"   random_seed = %d\n", psRandomSeed_);
   else fprintf(fOut,"#  random_seed = 2147483647\n");
   fprintf(fOut,"END\n");
}

// ************************************************************************
// A function for writing application section to a file
// ------------------------------------------------------------------------ 
void PsuadeData::writeApplicationSection(FILE *fOut) 
{
   fprintf(fOut,"APPLICATION\n");
   if (strcmp(pApplication_.appDriver_, "NONE"))
        fprintf(fOut,"   driver = %s\n", pApplication_.appDriver_);
   else fprintf(fOut,"   driver = NONE\n");
   if (strcmp(pApplication_.optDriver_, "NONE"))
        fprintf(fOut,"   opt_driver = %s\n", pApplication_.optDriver_);
   else fprintf(fOut,"   opt_driver = NONE\n");
   if (strcmp(pApplication_.auxOptDriver_, "NONE"))
        fprintf(fOut,"   aux_opt_driver = %s\n",pApplication_.auxOptDriver_);
   else fprintf(fOut,"   aux_opt_driver = NONE\n");
   if (strcmp(pApplication_.inputTemplate_, "NONE"))
        fprintf(fOut,"   input_template = %s\n",pApplication_.inputTemplate_);
   else fprintf(fOut,"#  input_template = NONE\n");
   if (strcmp(pApplication_.outputTemplate_, "NONE"))
        fprintf(fOut,"   output_template = %s\n",pApplication_.outputTemplate_);
   else fprintf(fOut,"#  output_template = NONE\n");
   if (pApplication_.maxParallelJobs_ != 1)
        fprintf(fOut,"   max_parallel_jobs = %d\n",pApplication_.maxParallelJobs_);
   else fprintf(fOut,"#  max_parallel_jobs = 1\n");
   if (pApplication_.minJobWaitTime_ != 1)
        fprintf(fOut,"   min_job_wait_time = %d\n",pApplication_.minJobWaitTime_);
   else fprintf(fOut,"#  min_job_wait_time = 1\n");
   if (pApplication_.maxJobWaitTime_ != 1)
        fprintf(fOut,"   max_job_wait_time = %d\n",pApplication_.maxJobWaitTime_);
   else fprintf(fOut,"#  max_job_wait_time = 1000000\n");
   if (pApplication_.runType_ & 1) fprintf(fOut,"   nondeterministic\n");
   else                            fprintf(fOut,"#  nondeterministic\n");
   if (pApplication_.runType_ & 2) fprintf(fOut,"   launch_only\n");
   else                            fprintf(fOut,"#  launch_only\n");
   if (pApplication_.runType_ & 4) fprintf(fOut,"   limited_launch_only\n");
   else                            fprintf(fOut,"#  limited_launch_only\n");
   if (pApplication_.runType_ & 8) fprintf(fOut,"   gen_inputfile_only\n");
   else                            fprintf(fOut,"#  gen_inputfile_only\n");
   if (pApplication_.launchInterval_ != 1)
        fprintf(fOut,"   launch_interval = %d\n",pApplication_.launchInterval_);
   else fprintf(fOut,"#  launch_interval = 1\n");
   if (pApplication_.saveFrequency_ != 1000000)
        fprintf(fOut,"   save_frequency = %d\n",pApplication_.saveFrequency_);
   else fprintf(fOut,"#  save_frequency = 1000000\n");
   fprintf(fOut,"END\n");
}

// ************************************************************************
// A function for writing analysis section to a file
// ------------------------------------------------------------------------ 
void PsuadeData::writeAnalysisSection(FILE *fOut) 
{
   int ii;

   fprintf(fOut,"ANALYSIS\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_MOMENT) != 0) 
        fprintf(fOut,"   analyzer method = Moment\n"); 
   else fprintf(fOut,"#  analyzer method = Moment\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ME) != 0) 
        fprintf(fOut,"   analyzer method = MainEffect\n"); 
   else fprintf(fOut,"#  analyzer method = MainEffect\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_IE) != 0) 
        fprintf(fOut,"   analyzer method = TwoParamEffect\n"); 
   else fprintf(fOut,"#  analyzer method = TwoParamEffect\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ANOVA) != 0) 
        fprintf(fOut,"   analyzer method = ANOVA\n"); 
   else fprintf(fOut,"#  analyzer method = ANOVA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_GLSA)!= 0) 
        fprintf(fOut,"   analyzer method = GLSA\n");
   else fprintf(fOut,"#  analyzer method = GLSA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_RSFA) != 0)
        fprintf(fOut,"   analyzer method = RSFA\n");
   else fprintf(fOut,"#  analyzer method = RSFA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_MOAT) != 0)
        fprintf(fOut,"   analyzer method = MOAT\n");
   else fprintf(fOut,"#  analyzer method = MOAT\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_SOBOL) != 0)
        fprintf(fOut,"   analyzer method = Sobol\n");
   else fprintf(fOut,"#  analyzer method = Sobol\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_CORRELATION) != 0)
        fprintf(fOut,"   analyzer method = Correlation\n");
   else fprintf(fOut,"#  analyzer method = Correlation\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_INTEGRATION) != 0)
        fprintf(fOut,"   analyzer method = Integration\n");
   else fprintf(fOut,"#  analyzer method = Integration\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_FAST) != 0)
        fprintf(fOut,"   analyzer method = FAST\n");
   else fprintf(fOut,"#  analyzer method = FAST\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_FF) != 0)
        fprintf(fOut,"   analyzer method = FF\n");
   else fprintf(fOut,"#  analyzer method = FF\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_PCA) != 0)
        fprintf(fOut,"   analyzer method = PCA\n");
   else fprintf(fOut,"#  analyzer method = PCA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ONESIGMA) != 0)
        fprintf(fOut,"   analyzer method = ARSM0\n");
   else fprintf(fOut,"#  analyzer method = ARSM0\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_FORM) != 0)
        fprintf(fOut,"   analyzer method = FORM\n");
   else fprintf(fOut,"#  analyzer method = FORM\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_RSSOBOL1) != 0)
        fprintf(fOut,"   analyzer method = RSMSobol1\n");
   else fprintf(fOut,"#  analyzer method = RSMSobol1\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_RSSOBOL2) != 0)
        fprintf(fOut,"   analyzer method = RSMSobol2\n");
   else fprintf(fOut,"#  analyzer method = RSMSobol2\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_RSSOBOLTSI) != 0)
        fprintf(fOut,"   analyzer method = RSMSobolTSI\n");
   else fprintf(fOut,"#  analyzer method = RSMSobolTSI\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_BSTRAP) != 0)
        fprintf(fOut,"   analyzer method = Bootstrap\n");
   else fprintf(fOut,"#  analyzer method = Bootstrap\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_RSSOBOLG) != 0)
        fprintf(fOut,"   analyzer method = RSMSobolG\n");
   else fprintf(fOut,"#  analyzer method = RSMSobolG\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ARSM) != 0)
        fprintf(fOut,"   analyzer method = ARSMNN\n");
   else fprintf(fOut,"#  analyzer method = ARSMNN\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ARSMMB) != 0)
        fprintf(fOut,"   analyzer method = ARSM1\n");
   else fprintf(fOut,"#  analyzer method = ARSM1\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_REL) != 0)
        fprintf(fOut,"   analyzer method = REL\n");
   else fprintf(fOut,"#  analyzer method = REL\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_AOPT) != 0)
        fprintf(fOut,"   analyzer method = AOPT\n");
   else fprintf(fOut,"#  analyzer method = AOPT\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_GOWER) != 0)
        fprintf(fOut,"   analyzer method = GOWER\n");
   else fprintf(fOut,"#  analyzer method = GOWER\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_DTEST) != 0)
        fprintf(fOut,"   analyzer method = DELTA\n");
   else fprintf(fOut,"#  analyzer method = DELTA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ETEST) != 0)
        fprintf(fOut,"   analyzer method = ETA\n");
   else fprintf(fOut,"#  analyzer method = ETA\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_ARSMMBBS) != 0)
        fprintf(fOut,"   analyzer method = ARSMG\n");
   else fprintf(fOut,"#  analyzer method = ARSMG\n");
   if ((pAnalysis_.analysisIntOptions_[0] & PSUADE_ANA_LSA) != 0)
        fprintf(fOut,"   analyzer method = LSA\n");
   else fprintf(fOut,"#  analyzer method = LSA\n");
   if (pAnalysis_.analysisIntOptions_[1]+1 <= pOutput_.nOutputs_)
        fprintf(fOut,"   analyzer output_id  = %d\n", 
                pAnalysis_.analysisIntOptions_[1]+1);
   else fprintf(fOut,"   analyzer output_id  = 1\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_MARS)
        fprintf(fOut,"   analyzer rstype = MARS\n");
   else fprintf(fOut,"#  analyzer rstype = MARS\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGR1)
        fprintf(fOut,"   analyzer rstype = linear\n"); 
   else fprintf(fOut,"#  analyzer rstype = linear\n"); 
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGR2)
        fprintf(fOut,"   analyzer rstype = quadratic\n");
   else fprintf(fOut,"#  analyzer rstype = quadratic\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGR3)
        fprintf(fOut,"   analyzer rstype = cubic\n");
   else fprintf(fOut,"#  analyzer rstype = cubic\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGR4)
        fprintf(fOut,"   analyzer rstype = quartic\n");
   else fprintf(fOut,"#  analyzer rstype = quartic\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_ANN)
        fprintf(fOut,"   analyzer rstype = ANN\n");
   else fprintf(fOut,"#  analyzer rstype = ANN\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGRS)
        fprintf(fOut,"   analyzer rstype = selective_regression\n");
   else fprintf(fOut,"#  analyzer rstype = selective_regression\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_GP1)
        fprintf(fOut,"   analyzer rstype = GP1\n");
   else fprintf(fOut,"#  analyzer rstype = GP1\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_GP2)
        fprintf(fOut,"   analyzer rstype = GP2\n");
   else fprintf(fOut,"#  analyzer rstype = GP2\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_SVM)
        fprintf(fOut,"   analyzer rstype = SVM\n");
   else fprintf(fOut,"#  analyzer rstype = SVM\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_PWL)
        fprintf(fOut,"   analyzer rstype = PWL\n");
   else fprintf(fOut,"#  analyzer rstype = PWL\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_TGP)
        fprintf(fOut,"   analyzer rstype = TGP\n");
   else fprintf(fOut,"#  analyzer rstype = TGP\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_MARSB)
        fprintf(fOut,"   analyzer rstype = MARSBag\n");
   else fprintf(fOut,"#  analyzer rstype = MARSBag\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_EARTH)
        fprintf(fOut,"   analyzer rstype = EARTH\n");
   else fprintf(fOut,"#  analyzer rstype = EARTH\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_SOTS)
        fprintf(fOut,"   analyzer rstype = sum_of_trees\n");
   else fprintf(fOut,"#  analyzer rstype = sum_of_trees\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGRL)
        fprintf(fOut,"   analyzer rstype = Legendre\n");
   else fprintf(fOut,"#  analyzer rstype = Legendre\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGRU)
        fprintf(fOut,"   analyzer rstype = user_regression\n");
   else fprintf(fOut,"#  analyzer rstype = user_regression\n");
   if (pAnalysis_.analysisIntOptions_[2] == PSUADE_RS_REGSG)
        fprintf(fOut,"   analyzer rstype = sparse_grid_regression\n");
   else fprintf(fOut,"#  analyzer rstype = sparse_grid_regression\n");
   if (pAnalysis_.legendreOrder_ > 0)
        fprintf(fOut,"   analyzer rs_legendre_order = %d\n", pAnalysis_.legendreOrder_);
   else fprintf(fOut,"#  analyzer rs_legendre_order = -1\n");
   if (pAnalysis_.analysisIntOptions_[6] >= 0 && 
       pAnalysis_.analysisIntOptions_[6] < pOutput_.nOutputs_)
        fprintf(fOut,"   analyzer regression_wgt_id = %d\n",
           pAnalysis_.analysisIntOptions_[6]+1);
   else fprintf(fOut,"#  analyzer regression_wgt_id = -1\n");
   if ((pAnalysis_.analysisIntOptions_[4] & PSUADE_GRAPHICS) != 0) 
        fprintf(fOut,"   graphics\n");
   else fprintf(fOut,"#  graphics\n");
   if ((pAnalysis_.analysisIntOptions_[4] & PSUADE_SAMPLE_GRAPHICS) != 0)
        fprintf(fOut,"   sample_graphics\n");
   else fprintf(fOut,"#  sample_graphics\n");
   fprintf(fOut,"   analyzer threshold = %e\n", 
           pAnalysis_.analysisDbleOptions_[0]);
   if (psFAMaxDataPts_ != 4000)
      fprintf(fOut,"   rs_max_pts = %d\n", psFAMaxDataPts_);

   if (pAnalysis_.numRSFilters_ > 0)
   {
      for (ii = 0; ii < pAnalysis_.numRSFilters_; ii++)
      {
        fprintf(fOut,"   analyzer rs_constraint = %s %s %e %e\n",
                pAnalysis_.RSFilters_[ii]->FilterDataFile_,
                pAnalysis_.RSFilters_[ii]->FilterIndexFile_,
                pAnalysis_.RSFilters_[ii]->FilterLBound_,
                pAnalysis_.RSFilters_[ii]->FilterUBound_);
      }
   }
   else
   {
      fprintf(fOut,"#  analyzer rs_constraint = psData indexFile Lbound Ubound\n");
   }
   if (pAnalysis_.numMOATFilters_ > 0)
   {
      for (ii = 0; ii < pAnalysis_.numMOATFilters_; ii++)
      {
        fprintf(fOut,"   analyzer moat_constraint = %s %s %e %e\n",
                pAnalysis_.MOATFilters_[ii]->FilterDataFile_,
                pAnalysis_.MOATFilters_[ii]->FilterIndexFile_,
                pAnalysis_.MOATFilters_[ii]->FilterLBound_,
                pAnalysis_.MOATFilters_[ii]->FilterUBound_);
      }
   }
   else
   {
      fprintf(fOut,"#  analyzer moat_constraint = psData indexFile Lbound Ubound\n");
   }
   if (strcmp(pAnalysis_.rsIndexFile_, "NONE"))
      fprintf(fOut,"   analyzer rs_index_file = %s\n",
              pAnalysis_.rsIndexFile_);
   else
      fprintf(fOut,"#  analyzer rs_index_file = indexFile\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) && 
       (pAnalysis_.optimizeIntOptions_[1] == 0)) 
        fprintf(fOut, "   optimization method = crude\n");
   else fprintf(fOut, "#  optimization method = crude\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 1)) 
        fprintf(fOut, "   optimization method = txmath\n");
   else fprintf(fOut, "#  optimization method = txmath\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 2)) 
        fprintf(fOut, "   optimization method = appspack\n");
   else fprintf(fOut, "#  optimization method = appspack\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 3)) 
        fprintf(fOut, "   optimization method = minpack\n");
   else fprintf(fOut, "#  optimization method = minpack\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 4)) 
        fprintf(fOut, "   optimization method = cobyla\n");
   else fprintf(fOut, "#  optimization method = cobyla\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 5)) 
        fprintf(fOut, "   optimization method = sm\n");
   else fprintf(fOut, "#  optimization method = sm\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 6)) 
        fprintf(fOut, "   optimization method = mm\n");
   else fprintf(fOut, "#  optimization method = mm\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 7)) 
        fprintf(fOut, "   optimization method = mm_adaptive\n");
   else fprintf(fOut, "#  optimization method = mm_adaptive\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 8)) 
        fprintf(fOut, "   optimization method = bobyqa\n");
   else fprintf(fOut, "#  optimization method = bobyqa\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[1] == 9)) 
        fprintf(fOut, "   optimization method = sce\n");
   else fprintf(fOut, "#  optimization method = sce\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[2] > 0)) 
        fprintf(fOut, "   optimization num_local_minima = %d\n",
                pAnalysis_.optimizeIntOptions_[2]);
   else fprintf(fOut, "#  optimization num_local_minima = 0\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[3] == 1)) 
        fprintf(fOut, "   optimization use_response_surface\n");
   else fprintf(fOut, "#  optimization use_response_surface\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[4] > 0)) 
        fprintf(fOut, "   optimization print_level = %d\n", 
                pAnalysis_.optimizeIntOptions_[4]);
   else fprintf(fOut, "#  optimization print_level = 0\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[5] > 0)) 
        fprintf(fOut, "   optimization num_fmin = %d\n", 
                pAnalysis_.optimizeIntOptions_[5]);
   else fprintf(fOut, "#  optimization num_fmin = 0\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeIntOptions_[6] >= 0)) 
        fprintf(fOut, "   optimization output_id = %d\n", 
                pAnalysis_.optimizeIntOptions_[6]+1);
   else fprintf(fOut, "#  optimization output_id = 0\n");

   if (pAnalysis_.optimizeIntOptions_[0] > 0)
        fprintf(fOut, "   optimization max_feval = %d\n", 
                pAnalysis_.optimizeIntOptions_[7]);
   else fprintf(fOut, "#  optimization max_feval = 10000\n");

   if (pAnalysis_.optimizeIntOptions_[0] > 0)
        fprintf(fOut, "   optimization deltax = %e\n",
                pAnalysis_.optimizeDbleOptions_[3]);
   else fprintf(fOut, "#  optimization deltax = 1.0e-6\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeDbleOptions_[0] != -1.0e50)) 
        fprintf(fOut, "   optimization fmin = %e\n",
                pAnalysis_.optimizeDbleOptions_[0]);
   else fprintf(fOut, "#  optimization fmin = not defined\n");

   if ((pAnalysis_.optimizeIntOptions_[0] > 0) &&
       (pAnalysis_.optimizeDbleOptions_[1] != 1.0e50)) 
        fprintf(fOut, "   optimization cutoff = %e\n",
                pAnalysis_.optimizeDbleOptions_[1]);
   else fprintf(fOut, "#  optimization cutoff = not defined\n");

   if (pAnalysis_.optimizeIntOptions_[0] > 0)
        fprintf(fOut, "   optimization tolerance = %e\n",
                pAnalysis_.optimizeDbleOptions_[2]);
   else fprintf(fOut, "#  optimization tolerance = not defined\n");

   if (pAnalysis_.analysisIntOptions_[3] == 0)
        fprintf(fOut,"#  diagnostics\n"); 
   else fprintf(fOut,"   diagnostics %d\n",pAnalysis_.analysisIntOptions_[3]); 

   if ((pAnalysis_.fileWriteFlag_ & 1) != 0)
        fprintf(fOut,"   file_write matlab\n");
   else fprintf(fOut,"#  file_write matlab\n");
   if (psConfig_ != NULL)
        fprintf(fOut,"   config_file = %s\n",psConfigFileName_);
   else fprintf(fOut,"#  config_file = NONE\n");
   if (pAnalysis_.useInputPDFs_ == 1)
        fprintf(fOut,"   use_input_pdfs\n");
   else fprintf(fOut,"#  use_input_pdfs\n");
   if (psConstraintSetOp_ == 1)
        fprintf(fOut,"   constraint_op_and\n");
   else fprintf(fOut,"#  constraint_op_and\n");
   if (psAnalysisInteractive_ == 1)
        fprintf(fOut,"   interactive\n");
   else fprintf(fOut,"#  interactive\n");
   if (psExpertMode_ == 1)
      fprintf(fOut,"   expert\n");
   if (psIOExpertMode_ == 1)
      fprintf(fOut,"   io_expert\n");
   if (psRSExpertMode_ == 1)
      fprintf(fOut,"   rs_expert\n");
   if (psAnaExpertMode_ == 1)
      fprintf(fOut,"   ana_expert\n");
   if (psOptExpertMode_ == 1)
      fprintf(fOut,"   opt_expert\n");
   if (psSamExpertMode_ == 1)
      fprintf(fOut,"   sam_expert\n");
   if (psPlotTool_ == 1)
      fprintf(fOut,"   scilab\n");
   fprintf(fOut,"END\n");
}

// ************************************************************************
// Return the state of a single sample, from pOutput_.sampleStates_ array
// If given id = -1, return the number of samples
// ------------------------------------------------------------------------ 
int PsuadeData::getSampleState(int sampleid)
{
   if (sampleid == -1)
      return pMethod_.nSamples_;
   else
      return pOutput_.sampleStates_[sampleid];
}

// ************************************************************************
// Return the value of a single sample input, from pInput_.sampleInputs_ array
// If given id = -1, return the number of inputs
// ------------------------------------------------------------------------ 
double PsuadeData::getSampleInput(int inputNumber, int sampleid)
{
   if (sampleid == -1)
      return pInput_.nInputs_;
   else
      return pInput_.sampleInputs_[sampleid * pInput_.nInputs_ + inputNumber];
}

// ************************************************************************
// Return the value of a single sample output, from pInput_.sampleInputs_ array
// If given id = -1, return the number of outputs
// ------------------------------------------------------------------------ 
double PsuadeData::getSampleOutput(int outputNumber, int sampleid)
{
   if (sampleid == -1)
      return pOutput_.nOutputs_;
   else
      return pOutput_.sampleOutputs_[sampleid * pOutput_.nOutputs_ + outputNumber];
}

// ************************************************************************
// ************************************************************************
// Subclass function definition
// ************************************************************************
// ************************************************************************

// ************************************************************************
// functions for psuadeInputSection
// ------------------------------------------------------------------------ 
psuadeInputSection::psuadeInputSection(): nInputs_(0), inputLBounds_(NULL), 
                    inputUBounds_(NULL),     inputNames_(NULL),
                    inputPDFs_(NULL),        inputMeans_(NULL), 
                    inputStdevs_(NULL),      inputAuxs_(NULL),
                    symbolTable_(NULL),      inputSettings_(NULL), 
                    inputNumSettings_(NULL), sampleInputs_(NULL)
{
}

psuadeInputSection::~psuadeInputSection()
{ 
   reset();
}

void psuadeInputSection::reset()
{ 
   int ii;

   if (inputLBounds_ != NULL) delete [] inputLBounds_; 
   if (inputUBounds_ != NULL) delete [] inputUBounds_; 
   if (inputNames_ != NULL) 
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (inputNames_[ii] != NULL) delete [] inputNames_[ii]; 
      delete [] inputNames_; 
   }
   if (inputPDFs_   != NULL) delete [] inputPDFs_;
   if (inputMeans_  != NULL) delete [] inputMeans_;
   if (inputStdevs_ != NULL) delete [] inputStdevs_;
   if (inputAuxs_   != NULL) delete [] inputAuxs_;
   if (symbolTable_ != NULL) delete [] symbolTable_; 
   if (inputSettings_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (inputSettings_[ii] != NULL) delete [] inputSettings_[ii]; 
      delete [] inputSettings_;
   }
   if (inputNumSettings_ != NULL) delete [] inputNumSettings_;
   if (sampleInputs_ != NULL) delete [] sampleInputs_;

   nInputs_          = 0;
   inputLBounds_     = NULL;
   inputUBounds_     = NULL;
   inputNames_       = NULL;
   inputPDFs_        = NULL;
   inputMeans_       = NULL;
   inputStdevs_      = NULL;
   inputAuxs_        = NULL;
   symbolTable_      = NULL;
   inputSettings_    = NULL;
   inputNumSettings_ = NULL;
   sampleInputs_     = NULL;
   corMatrix_.setDim(0,0);
   useInputPDFs_     = 0;
}

// ************************************************************************
// functions for psuadeOutputSection
// ------------------------------------------------------------------------ 
psuadeOutputSection::psuadeOutputSection() : nOutputs_(0),
                     outputNames_(NULL), sampleOutputs_(NULL),
                     sampleStates_(NULL)
{
}

psuadeOutputSection::~psuadeOutputSection()
{ 
   reset();
}

void psuadeOutputSection::reset()
{ 
   if (outputNames_ != NULL) 
   {
      for (int ii = 0; ii < nOutputs_; ii++)
         if (outputNames_[ii] != NULL) delete [] outputNames_[ii]; 
      delete [] outputNames_; 
   }
   if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
   if (sampleStates_  != NULL) delete [] sampleStates_;
   nOutputs_      = 0;
   outputNames_   = NULL;
   sampleOutputs_ = NULL;
   sampleStates_  = NULL;
}

// ************************************************************************
// functions for psuadeMethodSection
// ------------------------------------------------------------------------ 
psuadeMethodSection::psuadeMethodSection()
{ 
   reset();
}

psuadeMethodSection::~psuadeMethodSection()
{ 
}

void psuadeMethodSection::reset()
{ 
   nSamples_          = 0;
   samplingMethod_    = PSUADE_SAMP_MC;
   nReplications_     = 1;
   nRefinements_      = 0;
   refinementType_    = 0;
   refinementSize_    = 10000000;
   refNumRefinements_ = 0;
   sampleRandomize_   = 0;
}

// ************************************************************************
// functions for psuadeApplicationSection
// ------------------------------------------------------------------------ 
psuadeApplicationSection::psuadeApplicationSection()
{ 
   reset();
}

psuadeApplicationSection::~psuadeApplicationSection()
{ 
}

void psuadeApplicationSection::reset()
{ 
   strcpy(appDriver_, "NONE");
   strcpy(rsDriver_, "NONE");
   strcpy(optDriver_, "NONE");
   strcpy(auxOptDriver_, "NONE");
   strcpy(inputTemplate_, "NONE");
   strcpy(outputTemplate_, "NONE");
   maxParallelJobs_ = 1;
   minJobWaitTime_  = 1;
   maxJobWaitTime_  = 1000000;
   launchInterval_  = 1;
   saveFrequency_   = 1000000;
   runType_         = 0;
   useFunction_     = 0;
}

// ************************************************************************
// functions for psuadeAnalysisSection
// ------------------------------------------------------------------------ 
psuadeAnalysisSection::psuadeAnalysisSection()
{ 
   numRSFilters_   = 0;
   numMOATFilters_ = 0;
   RSFilters_   = NULL;
   MOATFilters_ = NULL;
   useInputPDFs_ = 0;
   legendreOrder_ = -1;
   reset();
}

psuadeAnalysisSection::~psuadeAnalysisSection()
{ 
}

void psuadeAnalysisSection::reset()
{ 
   for (int ii = 0; ii < 10; ii++)
   {
      optimizeIntOptions_[ii]  = 0;
      optimizeDbleOptions_[ii] = 0.0;
      analysisIntOptions_[ii]  = 0;
      analysisDbleOptions_[ii] = 0.0;
   }
   /* analysis method (default is none) */
   analysisIntOptions_[0]  = 0;
   /* analysis output ID */
   analysisIntOptions_[1]  = 0;
   /* analysis RS type (default MARS) */
   analysisIntOptions_[2]  = 0;
   /* analysis output level */
   analysisIntOptions_[3]  = 0;
   /* analysis sample graphics flags */
   analysisIntOptions_[4]  = 0;
   /* analysis sample transformation flags */
   analysisIntOptions_[5]  = 0;
   /* analysis regression weight output id */
   analysisIntOptions_[6]  = -1;
   /* analysis threshold */
   analysisDbleOptions_[0] = 1.0;
   /* turn off optimization */
   optimizeIntOptions_[0]  = 0;
   /* optimization method - none */
   optimizeIntOptions_[1]  = 0;
   /* number of local minima to be stored */
   optimizeIntOptions_[2]  = 1;
   /* use RS for optimization (default : no) */
   optimizeIntOptions_[3]  = 0;
   /* optimization output level */
   optimizeIntOptions_[4]  = 0;
   /* num_fmin to be found */
   optimizeIntOptions_[5]  = 0;
   /* output ID to be analyzed */
   optimizeIntOptions_[6]  = 0;
   /* maximum number of function evaluation */
   optimizeIntOptions_[7]  = 10000;
   /* function mimumum */
   optimizeDbleOptions_[0] = 0.0;
   /* cut off point for optimization */
   optimizeDbleOptions_[1] = 1.0e50;
   /* optimization tolerance */
   optimizeDbleOptions_[2] = 1.0e-6;
   /* optimization deltaX */
   optimizeDbleOptions_[3] = 1.0e-6;
   fileWriteFlag_ = 0;
   strcpy(specsFile_, "NONE");
   strcpy(rsIndexFile_, "NONE");
   if (numRSFilters_ > 0)
   {
      for (int jj = 0; jj < numRSFilters_; jj++)
      {
         if (RSFilters_[jj]->FilterDataFile_ != NULL)
            delete [] RSFilters_[jj]->FilterDataFile_;
         if (RSFilters_[jj]->FilterIndexFile_ != NULL)
            delete [] RSFilters_[jj]->FilterIndexFile_;
         if (RSFilters_[jj] != NULL) delete RSFilters_[jj];
      }
   }
   if (numMOATFilters_ > 0)
   {
      for (int kk = 0; kk < numMOATFilters_; kk++)
      {
         if (MOATFilters_[kk]->FilterDataFile_ != NULL)
            delete [] MOATFilters_[kk]->FilterDataFile_;
         if (MOATFilters_[kk]->FilterIndexFile_ != NULL)
            delete [] MOATFilters_[kk]->FilterIndexFile_;
         if (MOATFilters_[kk] != NULL) delete MOATFilters_[kk];
      }
   }
}


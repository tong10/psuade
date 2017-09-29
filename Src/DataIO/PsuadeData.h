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
// Definition for the class PsuadeData
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#ifndef __PSUADEDATAH__
#define __PSUADEDATAH__

#include <math.h>
//#include <stream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pData.h"
#include "Util/Matrix.h"

// ************************************************************************
// subclasses
// ************************************************************************

class psuadeInputSection 
{
public:
   int    nInputs_;
   double *inputLBounds_;
   double *inputUBounds_;
   char   **inputNames_;
   int    *inputPDFs_;
   double *inputMeans_;
   double *inputStdevs_;
   double *inputAuxs_;
   int    *symbolTable_;
   double **inputSettings_;
   int    *inputNumSettings_;
   double *sampleInputs_;
   Matrix corMatrix_;
   int    useInputPDFs_;

   psuadeInputSection();
   ~psuadeInputSection();
   void reset();
};

class psuadeOutputSection 
{
public:
   int    nOutputs_;
   char   **outputNames_;
   double *sampleOutputs_;
   int    *sampleStates_;

   psuadeOutputSection();
   ~psuadeOutputSection();
   void reset();
};

class psuadeMethodSection 
{
public:
   int    samplingMethod_;
   int    nSamples_;
   int    nReplications_;
   int    nRefinements_;
   int    refNumRefinements_;
   int    refinementType_;
   int    refinementSize_;
   int    sampleRandomize_;

   psuadeMethodSection();
   ~psuadeMethodSection();
   void reset();
};

class psuadeApplicationSection
{
public:
   char   appDriver_[200];
   char   rsDriver_[200];
   char   optDriver_[200];
   char   auxOptDriver_[200];
   char   inputTemplate_[200];
   char   outputTemplate_[200];
   int    maxParallelJobs_;
   int    maxJobWaitTime_;
   int    minJobWaitTime_;
   int    runType_;
   int    useFunction_;
   int    launchInterval_;
   int    saveFrequency_;

   psuadeApplicationSection();
   ~psuadeApplicationSection();
   void reset();
};

typedef struct 
{
   double FilterLBound_;
   double FilterUBound_;
   char   *FilterDataFile_;  /* psuade data file name as RS filter    */
   char   *FilterIndexFile_; /* index file that gives input indices   */
} psuadeFilter;

// ------------------------------------------------------------------------

class psuadeAnalysisSection 
{
public:
   int    analysisIntOptions_[10];
   double analysisDbleOptions_[10];
   int    optimizeIntOptions_[10];
   double optimizeDbleOptions_[10];
   int    fileWriteFlag_;
   char   specsFile_[200];
   int    numRSFilters_;
   int    numMOATFilters_;
   char   rsIndexFile_[200];
   int    useInputPDFs_;
   int    legendreOrder_;
   psuadeFilter **RSFilters_;
   psuadeFilter **MOATFilters_;

   psuadeAnalysisSection();
   ~psuadeAnalysisSection();
   void reset();
};

// ************************************************************************
// ************************************************************************
// main class definition 
// ************************************************************************
// ************************************************************************

class PsuadeData 
{

   char                     psuadeFileName_[200];
   int                      outputLevel_;
   int                      writeCnt_;
   psuadeInputSection       pInput_;
   psuadeOutputSection      pOutput_;
   psuadeMethodSection      pMethod_;
   psuadeApplicationSection pApplication_;
   psuadeAnalysisSection    pAnalysis_;
   pData                    auxData_;

public:

   PsuadeData();

   ~PsuadeData();

   int readPsuadeFile(char *fname);

   void writePsuadeFile(const char *fname, int);

   int getParameter(const char *, pData &);

   void updateInputSection(int nSamples, int nInputs, int *nSymbols,
                           double *lowerB, double *upperB,
                           double *sampleInputs, char **names);

   void updateOutputSection(int nSamples, int nOutputs,
                            double *sampleOutputs, int *sampleStates,
                            char **names);

   void updateMethodSection(int method, int nSamples, int nReps,
                            int nRefines, int randomize);

   void updateApplicationSection(char *appDriver, char *optDriver,
                                 int maxJobs);

   void updateAnalysisSection(int methods, int transform, int rstype,
			      int diag, int outputID, int usePDFs);

   void updateOptimizationSection(int methods, double tolerance);

   void processOutputData();

   void setOutputLevel(int);

   void readPsuadeIO(char *);

   int    getSampleState(int sampleid);
   double getSampleInput(int inputNumber, int sampleid);
   double getSampleOutput(int outputNumber, int sampleid);
   pData  *getAuxData() {return &auxData_;}

private:
    void readInputSection(FILE *);
    void readOutputSection(FILE *);
    void readMethodSection(FILE *);
    void readApplicationSection(FILE *);
    void readAnalysisSection(FILE *);
    int  getInputParameter(const char *, pData &);
    int  getOutputParameter(const char *, pData &);
    int  getMethodParameter(const char *, pData &);
    int  getApplicationParameter(const char *, pData &);
    int  getAnalysisParameter(const char *, pData &);

public:
    void writePsuadeIO(FILE *, int);
    void writeInputSection(FILE *);
    void writeOutputSection(FILE *);
    void writeMethodSection(FILE *);
    void writeSamplingSection(FILE *);
    void writeApplicationSection(FILE *);
    void writeAnalysisSection(FILE *);
};

#endif // __PSUADEDATAH__


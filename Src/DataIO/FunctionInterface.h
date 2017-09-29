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
// Definition for the class FunctionInterface
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#ifndef __FUNCTIONINTERFACEH__
#define __FUNCTIONINTERFACEH__

#include <math.h>
//#include <stream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Util/sysdef.h"
#include "DataIO/PsuadeData.h"
#include "FuncApprox/FuncApprox.h"

// ************************************************************************
// Class Definition 
// ************************************************************************

class FunctionInterface 
{
   int    nInputs_;
   int    nOutputs_;
   char   **inputNames_;
   char   **outputNames_;
   char   appDriver_[200];
   char   optDriver_[200];
   char   auxOptDriver_[200];
   char   appInputTemplate_[200];
   char   appOutputTemplate_[200];
   int    executionMode_;
   int    launchInterval_;
   int    appOptFlag_;
   int    useRSModel_;
   FuncApprox **rsPtrs_;
   int    printLevel_;
   int    *rsIndices_;
   double *rsValues_;
   int    rsRanFlag_;

public:
   FunctionInterface();
   ~FunctionInterface();

   int  loadInputData(int, char **);
   int  loadOutputData(int, char **);
   int  loadFunctionData(int, char **);
   int  setAsynchronousMode();
   int  setSynchronousMode();
   int  setOutputLevel(int);
   int  setLaunchInterval(int);
   int  setDriver(int);
   int  evaluate(int,int,double *,int,double *,int);
   int  getNumInputs();
   int  getNumOutputs();
   int  getDriver();
   char **getInputNames();
   char **getOutputNames();
   char *getApplicationDriver();
   char *getOptimizationDriver();
   char *getAuxOptimizationDriver();

private :
   int createInputFile(int, double *dinputs);
   int removePattern(char *, char *, char *);
   int substitutePattern(char *, char *, char *, double);
   int getPatternData(char *, char *, double*);
   int psLocalFunction(int, double *, int, double *);
   void *genRSModel(char *);
};

// ************************************************************************
// function to instantiate this class
// ------------------------------------------------------------------------

FunctionInterface *createFunctionInterface(PsuadeData *psuadeIO);
FunctionInterface *createFunctionInterfaceSimplified(PsuadeData *psuadeIO);
FunctionInterface *createFunctionInterfaceGivenAppDriver(int, int, char *fname);

#endif //  __FUNCTIONINTERFACEH__


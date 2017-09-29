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
// Functions for the class FuncApprox  
// Definitions for the class ANN
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "Main/Psuade.h"
#include "FuncApprox/FuncApprox.h"
#include "FuncApprox/Mars.h"
#include "FuncApprox/Earth.h"
#include "FuncApprox/GP1.h"
#include "FuncApprox/GP2.h"
#include "FuncApprox/TBGP.h"
#include "FuncApprox/SVM.h"
#include "FuncApprox/SelectiveRegression.h"
#include "FuncApprox/UserRegression.h"
#include "FuncApprox/LegendreRegression.h"
#include "FuncApprox/Regression.h"
#include "FuncApprox/Ann.h"
#include "FuncApprox/PWLinear.h"
#include "FuncApprox/MarsBagg.h"
#include "FuncApprox/SumOfTrees.h"
#include "FuncApprox/SGRegression.h"
#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "PDFLib/PDFBase.h"
#include "DataIO/pData.h"
#include "DataIO/PsuadeData.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FuncApprox::FuncApprox(int nInputs, int nSamples)
{
   if (nSamples <= 0 || nInputs <= 0) 
   {
      printf("FuncApprox::FuncApprox ERROR - invalid inputs. \n");
      printf("            nSamples = %d\n", nSamples);
      printf("            nInputs  = %d\n", nInputs);
      exit(-1);
   }
   outputLevel_ = 0;
   nSamples_    = nSamples;
   nInputs_     = nInputs;
   nPtsPerDim_  = 10;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   for (int i = 0 ; i < nInputs_; i++) lowerBounds_[i] = upperBounds_[i] = 0.0;
   weights_     = new double[nSamples_];;
   for (int j = 0 ; j < nSamples_; j++) weights_[j] = 1.0;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FuncApprox::~FuncApprox()
{
   if (lowerBounds_ != NULL) delete [] lowerBounds_;
   if (upperBounds_ != NULL) delete [] upperBounds_;
   if (weights_     != NULL) delete [] weights_;
}

// ************************************************************************
// get function approximator identifier
// ------------------------------------------------------------------------
int FuncApprox::getID()
{
   return faID_;
}

// ************************************************************************
// Set print level
// ------------------------------------------------------------------------
int FuncApprox::setOutputLevel(int level)
{
   if (level >= -1) outputLevel_ = level;
   return 0;
}

// ************************************************************************
// Set bounds for object class FuncApprox
// ------------------------------------------------------------------------
int FuncApprox::setBounds( double *lower, double *upper )
{
   for (int ii=0 ; ii<nInputs_; ii++) 
   {
      lowerBounds_[ii] = lower[ii];
      upperBounds_[ii] = upper[ii];
   }
   return 0;
}

// ************************************************************************
// load output weights
// ------------------------------------------------------------------------
int FuncApprox::loadWeights(int n, double *wgts)
{
   if (n != nSamples_)
   {
      printf("FuncApprox::loadWeights ERROR : invalid length %d.\n",n);
      exit(1);
   }
   if (weights_ != NULL) delete [] weights_;
   weights_ = NULL;
   for (int ii = 0 ; ii < n; ii++) 
   {
      if (wgts[ii] < 0.0)
      {
         printf("FuncApprox::loadWeights WARNING : weight < 0 - not used.\n");
         return 0;
      }
   }
   weights_ = new double[nSamples_];
   for (int jj = 0 ; jj < n; jj++) weights_[jj] = wgts[jj];
   return 0;
}

// ************************************************************************
// set number of points to generate in each dimension
// ------------------------------------------------------------------------
void FuncApprox::setNPtsPerDim(int npoints)
{
   if (npoints > 0) nPtsPerDim_ = npoints;
}

// ************************************************************************
// get number of points to generate in each dimension
// ------------------------------------------------------------------------
int FuncApprox::getNPtsPerDim()
{
   return nPtsPerDim_;
}

// ************************************************************************
// generate 1 dimensional data
// ------------------------------------------------------------------------
int FuncApprox::gen1DGridData(double *dble1, double *dble2, int int1, 
                      double *dble3,int *int3,double **dble4,double**dble5) 
{
   (void) dble1;
   (void) dble2;
   (void) dble3;
   (void) dble4;
   (void) dble5;
   (void) int1;
   (void) int3;
   return -1;
}

// ************************************************************************
// generate 2 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen2DGridData(double *dble1, double *dble2, int int1, 
                     int int2, double *dble3, int *int3, double **dble4,
                     double**dble5) 
{
   (void) dble1;
   (void) dble2;
   (void) dble3;
   (void) dble4;
   (void) dble5;
   (void) int1;
   (void) int2;
   (void) int3;
   return -1;
}

// ************************************************************************
// generate 3 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen3DGridData(double *dble1, double *dble2, int int1, 
                     int int2, int int3, double *dble3, int *int4, 
                     double **dble4,double **dble5) 
{
   (void) dble1;
   (void) dble2;
   (void) dble3;
   (void) dble4;
   (void) dble5;
   (void) int1;
   (void) int2;
   (void) int3;
   (void) int4;
   return -1;
}

// ************************************************************************
// generate 4 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen4DGridData(double *dble1, double *dble2, int int1, 
                     int int2, int int3, int int4, double *dble3, 
                     int *int5, double **dble4, double **dble5) 
{
   (void) dble1;
   (void) dble2;
   (void) dble3;
   (void) dble4;
   (void) dble5;
   (void) int1;
   (void) int2;
   (void) int3;
   (void) int4;
   (void) int5;
   return -1;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(double *X, double &std)
{
   std = 0.0;
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(int npts, double *X, double *Y,
                                      double *Ystd)
{
   for (int ii = 0; ii < npts; ii++)
   {
      Y[ii]  = 0.0;
      Ystd[ii]  = 0.0;
   }
   return 0.0;
}

// ************************************************************************
// Set parameters
// ------------------------------------------------------------------------
double FuncApprox::setParams(int argc, char **argv)
{
   (void) argc;
   (void) argv;
   return -1.0;
}

// ************************************************************************
// friend function (print current function approximator)
// ------------------------------------------------------------------------
extern "C" 
int getFAType(char *pString)
{
   int faType;

   faType = getInt(0, PSUADE_NUM_RS-1, pString);
#ifndef HAVE_MARS
   if (faType == 0) faType = -1;
#endif
#ifndef HAVE_SNNS
   if (faType == 5) faType = -1;
#endif
#ifndef HAVE_TPROS
   if (faType == 7) faType = -1;
#endif
#ifndef HAVE_GPMC
   if (faType == 8) faType = -1;
#endif
#ifndef HAVE_SVM
   if (faType == 9) faType = -1;
#endif
#ifndef HAVE_EARTH
   if (faType == 13) faType = -1;
#endif
   return faType;
}

// ************************************************************************
// friend function (print current function approximator)
// ------------------------------------------------------------------------
extern "C" 
void printThisFA(int faType)
{
   switch( faType )
   {
#ifdef HAVE_MARS
      case 0: printf("MARS model\n"); break;
#else
      case 0: printf("MARS model (not installed)\n"); break;
#endif
      case 1: printf("Linear regression model\n"); break;
      case 2: printf("Quadratic regression model\n"); break;
      case 3: printf("Cubic regression model\n"); break;
      case 4: printf("Quartic regression model\n"); break;
#ifdef HAVE_SNNS
      case 5: printf("Artificial neural network model\n"); break;
#else
      case 5: printf("Artificial neural network model (not installed)\n"); 
              break;
#endif
      case 6: printf("User-defined regression model\n"); break;
#ifdef HAVE_TPROS
      case 7: printf("Gaussian process (MacKay) model\n"); break;
#else
      case 7: printf("Gaussian process (MacKay) model (not installed) \n");
              break;
#endif
#ifdef HAVE_GPMC
      case 8: printf("Gaussian process (Rasmussen) model\n"); break;
#else
      case 8: printf("Gaussian process (Rasmussen) model (not installed)\n");
              break;
#endif
#ifdef HAVE_SVM
      case 9: printf("SVM-light (Joachims) model\n"); break;
#else
      case 9: printf("SVM-light (Joachims) model (not installed)\n");
              break;
#endif
      case 10: printf("Piecewise-linear model\n"); break;
#ifdef HAVE_TGP
      case 11: printf("11. Tree-based Gaussian Process\n"); break;
#else
      case 11: printf("Tree-based Gaussian Process (not installed)\n");
               break;
#endif
#ifdef HAVE_MARS
      case 12: printf("MARS with bagging model\n"); break;
#else
      case 12: printf("MARS with bagging model (not installed)\n"); break;
#endif
#ifdef HAVE_EARTH
      case 13: printf("Earth model\n"); break;
#else
      case 13: printf("Earth model (not installed)\n"); break;
#endif
      case 14: printf("Sum-of-trees model\n"); break;
      case 15: printf("Legendre polynomial regression\n"); break;
      case 16: printf("User-defined (nonpolynomial) regression\n"); break;
      case 17: printf("Sparse Grid polynomial regression\n"); break;
   }
}

// ************************************************************************
// friend function (print function approximator information)
// ------------------------------------------------------------------------
extern "C" 
int writeFAInfo()
{
   printDashes(0);
   printf("Available response surface tools: \n");
   printDashes(0);
   printf("Expert advices: \n");
#ifdef HAVE_MARS
   printf(" MARS - may have accuracy problem near domain boundary. Use\n"); 
   printf("   this option if sample size is sufficiently large (>100).\n"); 
#endif
   printf(" LINEAR, QUADRATIC, CUBIC, QUARTIC - good for small sample\n");
   printf("   sizes; and when the function is sufficiently smooth.\n");
   printf("   For higher than fourth order, use LEGENDRE (option 15) with\n");
   printf("   response surface expert mode turned on to select order.\n");
#ifdef HAVE_SNNS
   printf(" ANN - supported but currently not maintained.\n");
#endif
   printf(" SELECTIVE POLYNOMIAL REGRESSION - for high order polynomials\n");
   printf("   but your sample size is too small. (So you select certain\n");
   printf("   terms, provided you know which ones.)\n");
#ifdef HAVE_TPROS
   printf(" GAUSSIAN PROCESS - may  encounter non-definite covariance\n");
   printf("   matrix problem. However, if no such problem appears, it\n");
   printf("   especially good for small sample sizes (a few to a few tens).\n");   printf("   GP is relatively slow, so for sample sizes of more than a\n");
   printf("   few hundred, be patiet. For nonsmooth functions, try TGP.\n");
#endif
#ifdef HAVE_SVM
   printf(" SVM - provides 3 options: turn on rs_expert to select. Also,\n");
   printf("   use svmfind to search for best settings.\n");
#endif
   printf(" BOOTSTRAPPED MARS - intended to be used with adaptive sample\n");
   printf("   refinement, which adds more sample points near the boundary\n");
   printf("   of the parameter space.\n");
   printf(" SUM-OF-TREES REGRESSION - usually gives non-smooth response\n");
   printf("   responses. It is provided here for completeness, but is not\n");
   printf("   generally recommended.\n");
   printf(" SPARSE GRID REGRESSION - has to use sparse grid designs. Also,\n");
   printf("   you cannot use cross validation on sparse grid regression.\n");
   printf("   Use a test set (rstest) to validate your response surface.\n");
#ifdef HAVE_MARS
   printDashes(0);
   printf("0. MARS \n");
#endif
   printf("1. Linear regression \n");
   printf("2. Quadratic regression \n");
   printf("3. Cubic regression \n");
   printf("4. Quartic regression \n");
#ifdef HAVE_SNNS
   printf("5. Artificial neural network \n");
#endif
   printf("6. Selective polynomial regression \n");
#ifdef HAVE_TPROS
   printf("7. Gaussian process (MacKay)\n");
#endif
#ifdef HAVE_GPMC
   printf("8. Gaussian process (Rasmussen)\n");
#endif
#ifdef HAVE_SVM
   printf("9. SVM-light (Joachims)\n");
#endif
   printf("10. Piecewise-linear model\n");
#ifdef HAVE_TGP
   printf("11. Tree-based Gaussian Process\n");
#endif
#ifdef HAVE_MARS
   printf("12. MARS with bootstrap aggregating (bagging)\n");
#endif
#ifdef HAVE_EARTH
   printf("13. Earth (another MARS)\n");
#endif
   printf("14. Sum-of-trees model\n");
   printf("15. Legendre polynomial regression\n");
   printf("16. User-defined (nonpolynomial) regression\n");
   printf("17. Sparse Grid polynomial regression\n"); 
   return PSUADE_NUM_RS;
}

// ************************************************************************
// friend function (create a function approximator from a few parameters)
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFA(int faType, int nInputs, int nSamples)
{
   int        rsType;
   FuncApprox *faPtr;
   char       *params[1], winput[100];

   if (faType >= 0) rsType = faType;
   else
   {
      rsType = -1;
      while (rsType < 0 || rsType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(winput, "Please enter your choice ? ");
         rsType = getInt(0, PSUADE_NUM_RS-1, winput);
      }
   }
   if      (rsType == PSUADE_RS_MARS) faPtr = new Mars(nInputs, nSamples);
   else if (rsType == PSUADE_RS_ANN)  faPtr = new Ann(nInputs, nSamples);
   else if (rsType == PSUADE_RS_REGRS)
           faPtr = new SelectiveRegression(nInputs, nSamples);
   else if (rsType == PSUADE_RS_GP1)   faPtr = new GP1(nInputs, nSamples);
   else if (rsType == PSUADE_RS_GP2)   faPtr = new GP2(nInputs, nSamples);
   else if (rsType == PSUADE_RS_SVM)   faPtr = new SVM(nInputs, nSamples);
   else if (rsType == PSUADE_RS_PWL)   faPtr = new PWLinear(nInputs, nSamples);
   else if (rsType == PSUADE_RS_TGP)   faPtr = new TGP(nInputs, nSamples);
   else if (rsType == PSUADE_RS_MARSB) faPtr = new MarsBagg(nInputs, nSamples);
   else if (rsType == PSUADE_RS_EARTH) faPtr = new Earth(nInputs, nSamples);
   else if (rsType == PSUADE_RS_SOTS)  faPtr = new SumOfTrees(nInputs,nSamples);
   else if (rsType == PSUADE_RS_REGRL)
           faPtr = new LegendreRegression(nInputs,nSamples);
   else if (rsType == PSUADE_RS_REGRU)
           faPtr = new UserRegression(nInputs, nSamples);
   else if (rsType == PSUADE_RS_REGSG)
           faPtr = new SparseGridRegression(nInputs, nSamples);
   else if (rsType < 0 || rsType >= PSUADE_NUM_RS)
   {
      printf("ERROR: unknown response surface type.\n");
      return NULL;
   }
   else
   {
      faPtr = new Regression(nInputs, nSamples);
      params[0] = (char *) &rsType;
      faPtr->setParams(1, params);
   }
   return faPtr;
}

// ************************************************************************
// friend function (create a function approximator from a data file)
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFAInteractive(PsuadeData *psuadeIO, int flag)
{
   int        faType, nInputs, nSamples, nOutputs, wgtID, ii, nPtsPerDim;
   int        totPts, outputID, length;
   double     *sampleOutputs, *wghts, *Y;
   FuncApprox *faPtr;
   char       *params[1], winput[501];
   pData      pPtr, pInputs, pOutputs, pStates, pLower, pUpper;

   if (psuadeIO == NULL)
   {
      printf("ERROR: PsuadeData does not exist.\n");
      return NULL;
   }

   if ((flag & 1) == 0)
   {
      assert(psuadeIO->getParameter("ana_rstype", pPtr) == 0);
      faType = pPtr.intData_;
      if (faType < 0 || faType >= PSUADE_NUM_RS)
      {
         printf("createFA : faType (%d) not valid.\n", faType);
         exit(1);
      }
   }
   else
   {
      faType = -1;
      while (faType < 0 || faType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(winput, "Please enter your choice ? ");
         faType = getInt(0, PSUADE_NUM_RS-1, winput);
      }
   }

   psuadeIO->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   if (nSamples > psFAMaxDataPts_)
   {
      printf("PSUADE WARNING: For nSamples > %d,\n", psFAMaxDataPts_);
      printf("  it can be extremely expensive to create\n");
      printf("  response surfaces.\n");
      printf("  Please consult PSUADE developers before moving on.\n");
      exit(1);
   }
   psuadeIO->getParameter("ana_outputid", pPtr);
   outputID = pPtr.intData_;
   psuadeIO->getParameter("ana_regressionwgtid", pPtr);
   wgtID = pPtr.intData_;
   psuadeIO->getParameter("output_sample", pOutputs);
   psuadeIO->getParameter("input_sample", pInputs);
   psuadeIO->getParameter("output_states", pStates);
   psuadeIO->getParameter("input_lbounds", pLower);
   psuadeIO->getParameter("input_ubounds", pUpper);

   if      (faType == PSUADE_RS_MARS) faPtr = new Mars(nInputs, nSamples);
   else if (faType == PSUADE_RS_ANN)  faPtr = new Ann(nInputs, nSamples);
   else if (faType == PSUADE_RS_REGRS)
        faPtr = new SelectiveRegression(nInputs, nSamples);
   else if (faType == PSUADE_RS_GP1)   faPtr = new GP1(nInputs, nSamples);
   else if (faType == PSUADE_RS_GP2)   faPtr = new GP2(nInputs, nSamples);
   else if (faType == PSUADE_RS_SVM)   faPtr = new SVM(nInputs, nSamples);
   else if (faType == PSUADE_RS_PWL)   faPtr = new PWLinear(nInputs, nSamples);
   else if (faType == PSUADE_RS_TGP)   faPtr = new TGP(nInputs, nSamples);
   else if (faType == PSUADE_RS_MARSB) faPtr = new MarsBagg(nInputs, nSamples);
   else if (faType == PSUADE_RS_EARTH) faPtr = new Earth(nInputs, nSamples);
   else if (faType == PSUADE_RS_SOTS)  faPtr = new SumOfTrees(nInputs,nSamples);
   else if (faType == PSUADE_RS_REGRL)
   {
      faPtr = new LegendreRegression(nInputs,nSamples);
      psuadeIO->getParameter("ana_poly_order", pPtr);
      ii = pPtr.intData_;
      if (ii > 0)
      {
         params[0] = (char *) &ii;
         faPtr->setParams(1, params);
      }
   }
   else if (faType == PSUADE_RS_REGRU)
           faPtr = new UserRegression(nInputs, nSamples);
   else if (faType == PSUADE_RS_REGSG)
           faPtr = new SparseGridRegression(nInputs, nSamples);
   else
   {
      faPtr = new Regression(nInputs, nSamples);
      params[0] = (char *) &faType;
      faPtr->setParams(1, params);
   }
   faPtr->setBounds(pLower.dbleArray_, pUpper.dbleArray_);

   nPtsPerDim = 256;
   totPts = 1000001;
   while (totPts > 1000000)
   {
      nPtsPerDim = nPtsPerDim / 2;
      totPts = 1;
      for (ii = 0; ii < nInputs; ii++)
      {
         totPts *= nPtsPerDim;
         if (totPts > 1000000) break;
      }
   }
   faPtr->setNPtsPerDim(nPtsPerDim);

   if (wgtID >= 0 && wgtID < nOutputs)
   {
      wghts = new double[nSamples];
      for (ii = 0; ii < nSamples; ii++)
         wghts[ii] = sampleOutputs[ii*nOutputs+wgtID];
      faPtr->loadWeights(nSamples, wghts);
      delete [] wghts;
   }

   if (flag & 2)
   {
      Y = new double[nSamples];
      for (ii = 0; ii < nSamples; ii++)
         Y[ii] = pOutputs.dbleArray_[nOutputs*ii+outputID];
      length = -999;
      faPtr->genNDGridData(pInputs.dbleArray_, Y, &length,
                           NULL, NULL);
      delete [] Y;
   }
   return faPtr;
}

// ************************************************************************
// friend function (create a function approximator given a file name)
// perform PDF transformation
// check invalid sample points
// RS type from file 
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFAFromFile(char *fname, int outputID)
{
   int        ii, nInputs, nOutputs, nSamples, status, *sampleStates;
   double     *sampleInputs, *sampleOutputs;
   PsuadeData *psuadeIO = new PsuadeData();
   FuncApprox *faPtr;
   pData      pPtr, pInpData, pOutData, pStates;

   for (ii = strlen(fname)-1; ii >= 0; ii--) if (fname[ii] == '/') break;
   status = psuadeIO->readPsuadeFile(fname);
   if (status != 0) return NULL;

   psuadeIO->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO->getParameter("input_sample", pInpData);
   sampleInputs = pInpData.dbleArray_;
   psuadeIO->getParameter("output_sample", pOutData);
   sampleOutputs = pOutData.dbleArray_;
   psuadeIO->getParameter("output_states", pStates);
   sampleStates = pStates.intArray_;

   for (ii = 0; ii < nSamples; ii++)
   {
      if (sampleStates[ii] != 1 || 
          sampleOutputs[nOutputs*ii+outputID] == PSUADE_UNDEFINED)
      {
         printf("FunctionInterface::genRSModel ERROR - invalid output.\n");
         exit(1);
      }
   }

   PDFTransform(psuadeIO, nSamples, nInputs, sampleInputs);
   psuadeIO->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                sampleInputs,NULL);
   psuadeIO->getParameter("ana_outputid", pPtr);
   ii = pPtr.intData_;
   psuadeIO->updateAnalysisSection(-1,-1,-1,-1,outputID,-1);
   faPtr = genFAInteractive(psuadeIO, 2);
   psuadeIO->updateAnalysisSection(-1,-1,-1,-1,ii,-1);
   delete psuadeIO;
   return faPtr;
}


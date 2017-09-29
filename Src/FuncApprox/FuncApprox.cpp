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

#include "Psuade.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "Earth.h"
#include "GP1.h"
#include "GP2.h"
#include "TBGP.h"
#include "SVM.h"
#include "SelectiveRegression.h"
#include "UserRegression.h"
#include "LegendreRegression.h"
#include "Regression.h"
#include "Ann.h"
#include "PWLinear.h"
#include "MarsBagg.h"
#include "SumOfTrees.h"
#include "SGRegression.h"
#include "GradLegendreRegression.h"
#include "Kriging.h"
#include "NPLearning.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFBase.h"
#include "pData.h"
#include "PsuadeData.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

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
   for (int ii = 0 ; ii < nInputs_; ii++)
      lowerBounds_[ii] = upperBounds_[ii] = 0.0;
   weights_ = new double[nSamples_];;
   for (int jj = 0 ; jj < nSamples_; jj++) weights_[jj] = 1.0;
   XMeans_ = NULL;
   XStds_  = NULL;
   YMean_ = 0.0;
   YStd_ = 0.0;
}

// ************************************************************************
// Copy constructor by Bill Oliver 
// ------------------------------------------------------------------------
FuncApprox::FuncApprox(const FuncApprox & fa)
{
   outputLevel_ = fa.outputLevel_;
   nSamples_ = fa.nSamples_;
   nInputs_ = fa.nInputs_;
   nPtsPerDim_ = fa.nPtsPerDim_;
   faID_ = fa.faID_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   weights_     = new double[nSamples_];
 
   for(int ii = 0; ii < nInputs_; ii++)
   {
      lowerBounds_[ii] = fa.lowerBounds_[ii];
      upperBounds_[ii] = fa.upperBounds_[ii];
   }
   for(int ii = 0; ii < nSamples_; ii++) weights_[ii] = fa.weights_[ii];
}

// ************************************************************************
// operator= by Bill Oliver 
// ------------------------------------------------------------------------
FuncApprox & FuncApprox::operator=(const FuncApprox & fa)
{
   if(this == &fa)  return *this;
   // free lowerBounds_, upperBounds_, and weights_
   delete [] lowerBounds_;
   delete [] upperBounds_;
   delete [] weights_;
  
   outputLevel_ = fa.outputLevel_;
   nSamples_ = fa.nSamples_;
   nInputs_ = fa.nInputs_;
   nPtsPerDim_ = fa.nPtsPerDim_;
   faID_ = fa.faID_;
   lowerBounds_ = new double[nInputs_];
   upperBounds_ = new double[nInputs_];
   weights_     = new double[nSamples_];

   for (int ii = 0; ii < nInputs_; ii++)
   {
      lowerBounds_[ii] = fa.lowerBounds_[ii];
      upperBounds_[ii] = fa.upperBounds_[ii];
   }
   for(int ii = 0; ii < nSamples_; ii++) weights_[ii] = fa.weights_[ii];
   return *this;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FuncApprox::~FuncApprox()
{
   if (lowerBounds_ != NULL) delete [] lowerBounds_;
   if (upperBounds_ != NULL) delete [] upperBounds_;
   if (weights_     != NULL) delete [] weights_;
   if (XMeans_      != NULL) delete [] XMeans_;
   if (XStds_       != NULL) delete [] XStds_;
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
int FuncApprox::gen1DGridData(double*,double*,int,double*,int*,double**,
                              double**) 
{
   return -1;
}

// ************************************************************************
// generate 2 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen2DGridData(double*,double*,int,int,double*,int*,
                              double**,double**) 
{
   return -1;
}

// ************************************************************************
// generate 3 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen3DGridData(double*,double*,int,int,int,double*,int*, 
                              double**,double**) 
{
   return -1;
}

// ************************************************************************
// generate 4 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen4DGridData(double*,double*,int,int,int,int,double*, 
                              int*,double**,double**) 
{
   return -1;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(double *, double &std)
{
   std = 0.0;
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(int npts, double *, double *Y,
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
double FuncApprox::setParams(int, char **)
{
   return -1.0;
}

// ************************************************************************
// initialize scaling for inputs
// ------------------------------------------------------------------------
int FuncApprox::initInputScaling(double *XIn, double *XOut, int flag)
{
   int    ii, jj;
   double ddata;
   char   pString[500], response[100];
                                                                                
   if (XMeans_ != NULL) delete [] XMeans_;
   if (XStds_  != NULL) delete [] XStds_;
   XMeans_ = new double[nInputs_];
   XStds_ = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
   {
      XMeans_[ii] = 0.0;
      XStds_[ii] = 1.0;
   }
   for (ii = 0; ii < nInputs_*nSamples_; ii++) XOut[ii] = XIn[ii];
   if (flag == 1)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nSamples_; jj++) ddata += XIn[jj*nInputs_+ii];
         XMeans_[ii] = ddata / (double) nSamples_;
         ddata = 0.0;
         for (jj = 0; jj < nSamples_; jj++)
            ddata += pow(XIn[jj*nInputs_+ii] - XMeans_[ii], 2.0);
         XStds_[ii] = sqrt(ddata / (double) (nSamples_ - 1));
         if (XStds_[ii] == 0.0) XStds_[ii] = 1.0;
         for (jj = 0; jj < nSamples_; jj++)
            XOut[jj*nInputs_+ii] = (XIn[jj*nInputs_+ii]-XMeans_[ii])/
                                   XStds_[ii];
         printf("Input %d scaling info : mean, std = %e %e\n",ii+1,
                XMeans_[ii], XStds_[ii]);
      }
   }
   else if (psRSExpertMode_ == 1)
   {
      sprintf(pString, "Scale the sample matrix ? (y or n) ");
      getString(pString, response);
      if (response[0] == 'y')
      {
         for (ii = 0; ii < nInputs_; ii++)
         {
            ddata = 0.0;
            for (jj = 0; jj < nSamples_; jj++) ddata += XIn[jj*nInputs_+ii];
            XMeans_[ii] = ddata / (double) nSamples_;
            ddata = 0.0;
            for (jj = 0; jj < nSamples_; jj++)
               ddata += pow(XIn[jj*nInputs_+ii] - XMeans_[ii], 2.0);
            XStds_[ii] = sqrt(ddata / (double) (nSamples_ - 1));
            if (XStds_[ii] == 0.0) XStds_[ii] = 1.0;
            for (jj = 0; jj < nSamples_; jj++)
               XOut[jj*nInputs_+ii] = (XIn[jj*nInputs_+ii]-XMeans_[ii])/
                                      XStds_[ii];
            printf("Input %d scaling info : mean, std = %e %e\n",ii+1,
                   XMeans_[ii], XStds_[ii]);
         }
      }
   }
   return 0;
}

// ************************************************************************
// initialize scaling for output
// ------------------------------------------------------------------------
int FuncApprox::initOutputScaling(double *YIn, double *YOut)
{
   int    ii;
   double ddata;
                                                                                
   ddata = 0.0;
   for (ii = 0; ii < nSamples_; ii++) ddata += YIn[ii];
   YMean_ = ddata / (double) nSamples_;
   ddata = 0.0;
   for (ii = 0; ii < nSamples_; ii++)
      ddata += pow(YIn[ii] - YMean_, 2.0);
   YStd_ = sqrt(ddata / (double) (nSamples_ - 1));
   if (YStd_ == 0.0) YStd_ = 1.0;
   for (ii = 0; ii < nSamples_; ii++)
      YOut[ii] = (YIn[ii] - YMean_) / YStd_;
   return 0;
}

// ************************************************************************
// generate m-dimensional grid data
// ------------------------------------------------------------------------
int FuncApprox::genNDGrid(int *nPts, double **XOut) 
{
   int    ii, mm, totPts;
   double *HX, *Xloc;

   if (nInputs_ > 21)
   {
      printf("FuncApprox genNDGrid INFO: nInputs > 21 not supported.\n");
      (*nPts) = 0;
      (*XOut) = NULL;
      return 0;
   }

   if (nInputs_ == 21 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 20 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 19 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 18 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 17 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 16 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 15 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 14 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
   if (nInputs_ == 13 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 12 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 11 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
   if (nInputs_ == 10 && nPtsPerDim_ >    4) nPtsPerDim_ =  4;
   if (nInputs_ ==  9 && nPtsPerDim_ >    5) nPtsPerDim_ =  5;
   if (nInputs_ ==  8 && nPtsPerDim_ >    6) nPtsPerDim_ =  6;
   if (nInputs_ ==  7 && nPtsPerDim_ >    8) nPtsPerDim_ =  8;
   if (nInputs_ ==  6 && nPtsPerDim_ >   10) nPtsPerDim_ = 10;
   if (nInputs_ ==  5 && nPtsPerDim_ >   16) nPtsPerDim_ = 16;
   if (nInputs_ ==  4 && nPtsPerDim_ >   32) nPtsPerDim_ = 32;
   if (nInputs_ ==  3 && nPtsPerDim_ >   64) nPtsPerDim_ = 64;
   if (nInputs_ ==  2 && nPtsPerDim_ > 1024) nPtsPerDim_ = 1024;
   if (nInputs_ ==  1 && nPtsPerDim_ > 8192) nPtsPerDim_ = 8192;
   totPts = nPtsPerDim_;
   for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
   HX = new double[nInputs_];
   for (ii = 0; ii < nInputs_; ii++)
      HX[ii] = (upperBounds_[ii] - lowerBounds_[ii]) /
               (double) (nPtsPerDim_ - 1);

   (*XOut) = new double[nInputs_ * totPts];
   (*nPts) = totPts;
   Xloc  = new double[nInputs_];

   for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];

   for (mm = 0; mm < totPts; mm++)
   {
      for (ii = 0; ii < nInputs_; ii++ ) (*XOut)[mm*nInputs_+ii] = Xloc[ii];
      for (ii = 0; ii < nInputs_; ii++ )
      {
         Xloc[ii] += HX[ii];
         if (Xloc[ii] < upperBounds_[ii] ||
             PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
         else Xloc[ii] = lowerBounds_[ii];
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// ************************************************************************
// ************************************************************************
// friend function (print current function approximator)
// ------------------------------------------------------------------------
extern "C" 
int getFAType(char *pString)
{
   int faType;

   faType = getInt(0, PSUADE_NUM_RS, pString);
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
      case 10: printf("Derivative-based Legendre polynomial regression\n"); break;
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
      case 18: printf("Kriging\n"); break;
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
   printf("   especially good for small sample sizes (a few to a few tens).\n");
   printf("   GP is relatively slow, so for sample sizes of more than a\n");
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
   printf(" KRIGING - This is another form of the Gaussian process which\n");
   printf("   uses deterministic optimization to compute the hyperparameters.\n");
   printf("   This method is good for up to about 2000 sample points;\n");
   printf("   otherwise it may be computationally expensive.\n");
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
   printf("10. Derivative-based Legendre polynomial regression\n");
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
   printf("18. Kriging\n"); 
   return PSUADE_NUM_RS;
}

// ************************************************************************
// friend function (create a function approximator from a few parameters)
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFA(int faType, int nInputs, int, int nSamples)
{
   int        rsType;
   FuncApprox *faPtr=NULL;
   char       *params[1], winput[100];

   if (faType >= 0) rsType = faType;
   else
   {
      rsType = -1;
      while (rsType < 0 || rsType >= PSUADE_NUM_RS)
      {
         writeFAInfo();
         sprintf(winput, "Please enter your choice ? ");
         rsType = getInt(0, PSUADE_NUM_RS, winput);
      }
   }
   if      (rsType == PSUADE_RS_MARS) faPtr = new Mars(nInputs, nSamples);
   else if (rsType == PSUADE_RS_ANN)  faPtr = new Ann(nInputs, nSamples);
   else if (rsType == PSUADE_RS_REGRS)
           faPtr = new SelectiveRegression(nInputs, nSamples);
   else if (rsType == PSUADE_RS_GP1)   faPtr = new GP1(nInputs, nSamples);
   else if (rsType == PSUADE_RS_GP2)   faPtr = new GP2(nInputs, nSamples);
   else if (rsType == PSUADE_RS_SVM)   faPtr = new SVM(nInputs, nSamples);
   else if (rsType == PSUADE_RS_REGRGL)
           faPtr = new GradLegendreRegression(nInputs, nSamples);
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
   else if (rsType == PSUADE_RS_KR)
           faPtr = new Kriging(nInputs, nSamples);
   else if (faType == PSUADE_RS_NPL)
           faPtr = new NPLearning(nInputs, nSamples);
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
   int        totPts, outputID, length, printLevel;
   double     *wghts, *Y;
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
   psuadeIO->getParameter("ana_diagnostics", pPtr);
   printLevel = pPtr.intData_;
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
   else if (faType == PSUADE_RS_REGRGL)
           faPtr = new GradLegendreRegression(nInputs, nSamples);
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
   else if (faType == PSUADE_RS_KR)
           faPtr = new Kriging(nInputs, nSamples);
   else if (faType == PSUADE_RS_NPL)
           faPtr = new NPLearning(nInputs, nSamples);
   else
   {
      faPtr = new Regression(nInputs, nSamples);
      params[0] = (char *) &faType;
      faPtr->setParams(1, params);
   }
   faPtr->setBounds(pLower.dbleArray_, pUpper.dbleArray_);
   faPtr->setOutputLevel(printLevel);

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
         wghts[ii] = pOutputs.dbleArray_[nOutputs*ii+wgtID];
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
   if (status != 0){
     delete psuadeIO;
     return NULL;
   }

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
         printf("FuncApprox::genRSModel ERROR - invalid output.\n");
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


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
// Functions for the class CobylaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "CobylaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

#ifdef HAVE_COBYLA
#include "cobyla.h"
#endif
#define psCobylaMaxSaved_ 10000
int     psCobylaNSaved_=0;
double  psCobylaSaveX_[psCobylaMaxSaved_*10];
double  psCobylaSaveY_[psCobylaMaxSaved_*10];
int     psCobylaSaveHistory_=0;
char    psCobylaExec_[5000];
int     psCobylaNConstr_=-1;
int     psCCurrDriver_=-1;
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
   int evaluateFunction(int nInputs, int nConstraints, double *XValues, 
                        double *YValue, double *constraints, void *data)
   {
      int    ii, jj, kk, funcID, nOuts, ignoreFlag, found, index;
      double *localY, ddata; 
      oData  *odata;
      FILE   *fp;

      odata    = (oData *) data;
      nOuts    = psCobylaNConstr_ + 1;
      localY   = (double *) malloc(nOuts * sizeof(double));

      found = 0;
      for (ii = 0; ii < psCobylaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psCobylaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            index = ii;
            printf("Cobyla: simulation results reuse (%d).\n",ii+1);
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
         odata->funcIO_->evaluate(funcID,nInputs,XValues,nOuts,localY,0);
         funcID = odata->numFuncEvals_++;
         (*YValue) = localY[0];
         for (ii = 0; ii < nOuts-1; ii++) constraints[ii] = localY[ii+1];
      }
      else
      {
         localY[0] = psCobylaSaveY_[index*nOuts];
         (*YValue) = localY[0];
         for (ii = 0; ii < nOuts-1; ii++) 
            constraints[ii] = psCobylaSaveY_[index*nOuts+ii+1];
      }

      if (psCobylaSaveHistory_ == 1 && found == 0 &&
          (psCobylaNSaved_+1)*nInputs < psCobylaMaxSaved_*10 &&
          (psCobylaNSaved_+1)*nOuts < psCobylaMaxSaved_*10 &&
          psCobylaNSaved_ < psCobylaMaxSaved_)
      {
         for (jj = 0; jj < nInputs; jj++)
            psCobylaSaveX_[psCobylaNSaved_*nInputs+jj] = XValues[jj];
         for (jj = 0; jj < nOuts; jj++)
            psCobylaSaveY_[psCobylaNSaved_*nOuts+jj] = localY[jj];
         psCobylaNSaved_++;
         fp = fopen("psuade_cobyla_history","w");
         if (fp != NULL)
         {
            for (ii = 0; ii < psCobylaNSaved_; ii++)
            {
               fprintf(fp, "999 %d %d ", nInputs, nOuts);
               for (kk = 0; kk < nInputs; kk++)
                  fprintf(fp, "%24.16e ", psCobylaSaveX_[ii*nInputs+kk]);
               for (kk = 0; kk < nOuts; kk++)
                  fprintf(fp, "%24.16e ", psCobylaSaveY_[ii*nOuts+kk]);
               fprintf(fp, "\n");
            }
            fclose(fp);
         }
      }
      if (odata->outputLevel_ > 4)
      {
         printf("CobylaOptimizer %6d : \n", odata->numFuncEvals_); 
         for (ii = 0; ii < nInputs; ii++)
            printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
         printf("    ObjY  = %16.8e\n", localY[0]);
      }

      ignoreFlag = 0;
      for (ii = 0; ii < nConstraints; ii++)
         if (constraints[ii] < 0.0) ignoreFlag = 1;
      if ((ignoreFlag == 0) && ((*YValue) < odata->optimalY_))
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
      }
      free(localY);
      return 0;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
CobylaOptimizer::CobylaOptimizer()
{
   printAsterisks(PL_INFO, 0);
   printf("*   COBYLA Optimizer Usage Information\n");
   printEquals(PL_INFO, 0);
   printf("* - To run this optimizer, first make sure opt_driver has\n");
   printf("*   been initialized to point to your optimization objective\n");
   printf("*   function evaluator\n");
   printf("* - Set optimization tolerance in your PSUADE input file\n");
   printf("* - Set maximum number of iterations in PSUADE input file\n");
   printf("* - Set num_local_minima to perform multistart optimization\n");
   printf("* - Set optimization print_level to give additonal outputs\n");
   printf("* - In Opt EXPERT mode, the optimization history log will be\n");
   printf("*   turned on automatically. Previous psuade_cobyla_history\n");
   printf("*   file will also be reused.\n");
   printf("* - If your opt_driver is a response surface which has more\n");
   printf("*   inputs than the number of optimization inputs, you can fix\n");
   printf("*   some driver inputs by creating a (analyzer) rs_index_file.\n");
   printf("* - The optimization objective function is: \n\n");
   printf("        F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)) \n\n");
   printf("    where C1, C2, ..., CM are the results of evaluating the\n");
   printf("    constraint functions. The constraint function evaluation\n");
   printf("    should be performed in the function evaluator so make \n");
   printf("    sure the number of outputs in your evaluator is M+1.\n");
   printAsterisks(PL_INFO, 0);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CobylaOptimizer::~CobylaOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void CobylaOptimizer::optimize(oData *odata)
{
   int    nInputs, nConstraints, ii, jj, nIns, nOuts, maxfun;
   int    ntimes=3, printLevel=0, nOutputs, cmpFlag, token;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp;
   char   filename[500], *cString;
   FILE   *fp=NULL;

   printLevel = odata->outputLevel_;
   nInputs = odata->nInputs_;
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   maxfun = odata->maxFEval_;
   odata->optimalY_ = 1.0e50;
   XValues = new double[nInputs+1];
   for (ii = 0; ii < nInputs; ii++) XValues[ii] = odata->initialX_[ii];
   rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
   for (ii = 1; ii < nInputs; ii++) 
   {
      dtemp = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
      if (dtemp < rhobeg) rhobeg = dtemp;
   }
   rhobeg *= 0.5;
   rhoend = rhobeg * odata->tolerance_;
   if (rhobeg < rhoend)
   {
      printf("CobylaOptimizer WARNING: tolerance too large.\n");
      printf("                         tolerance reset to 1.0e-4.\n");
      rhoend = rhobeg * 1.0e-4;
   }
   nOutputs = odata->nOutputs_;

   if ((odata->setOptDriver_ & 1))
   {
      printf("Cobyla: setting optimization simulation driver.\n");
      psCCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   printAsterisks(PL_INFO, 0);
   printf("Cobyla optimizer: max fevals   = %d\n", odata->maxFEval_);
   printf("Cobyla optimizer: tolerance    = %e\n", odata->tolerance_);
   printEquals(PL_INFO, 0);
   psCobylaNConstr_ = nOutputs - 1;
   printf("Cobyla optimizer: nConstraints = %d\n",psCobylaNConstr_); 

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psCobylaSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL)
      {
         printf("Cobyla: use history has been turned on.\n");
         fp = fopen("psuade_cobyla_history","r");
         if (fp != NULL)
         {
            psCobylaNSaved_ = 0;
            token = 999;
            while (token == 999)
            {
               fscanf(fp, "%d", &token);
               if (token != 999) 
               {
                  fclose(fp);
                  break;
               }
               fscanf(fp, "%d %d", &nIns, &nOuts);
               if (nIns != nInputs)
               {
                  printf("Cobyla: history file has invalid nInputs.\n");
                  fclose(fp);
                  psCobylaNSaved_ = 0;
                  break;
               }
               else if (nOuts != psCobylaNConstr_+1)
               {
                  printf("Cobyla: history file has invalid nOutputs.\n");
                  fclose(fp);
                  psCobylaNSaved_ = 0;
                  break;
               }
               ii = psCobylaNSaved_;
               for (jj = 0; jj < nInputs; jj++)
                  fscanf(fp, "%lg", &psCobylaSaveX_[ii*nInputs+jj]);
               for (jj = 0; jj < nOuts; jj++)
                  fscanf(fp, "%lg", &psCobylaSaveY_[ii*nOuts+jj]);
               psCobylaNSaved_++;
               if (psCobylaNSaved_ > psCobylaMaxSaved_*10/nInputs)
               {
                  printf("Cobyla: History file has too much data.\n");
                  printf("        Truncate after %d samples.\n",
                         psCobylaNSaved_);
                  fclose(fp);
                  break;
               }
               if (psCobylaNSaved_ > psCobylaMaxSaved_*10/nOuts)
               {
                  printf("PSUADE Cobyla: history file has too much data.\n");
                  printf("               Ask PSUADE developer for help.\n");
                  fclose(fp);
                  break;
               }
            }
         }
      }
   }

#ifdef HAVE_COBYLA
   nConstraints = psCobylaNConstr_; 
   odata->numFuncEvals_ = 0;
   cobyla(nInputs, nConstraints, XValues, rhobeg, rhoend, printLevel, 
          &maxfun, evaluateFunction, (void *) odata);
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
   odata->optimalY_ = XValues[nInputs];
   printf("Cobyla optimizer: number of function evaluation = %d\n",
          odata->numFuncEvals_);
#else
   printf("ERROR : Cobyla optimizer not installed.\n");
   exit(1);
#endif

   if ((odata->setOptDriver_ & 2))
   {
      printf("Cobyla: setting back to original simulation driver.\n");
      odata->funcIO_->setDriver(psCCurrDriver_);
   }
   delete [] XValues;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
CobylaOptimizer& CobylaOptimizer::operator=(const CobylaOptimizer &)
{
   printf("CobylaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


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
#include <iostream>
#include <fstream>
#include <string>

#include "CobylaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

#ifdef HAVE_COBYLA
#include "cobyla.h"
#endif
#define psCobylaMaxSaved_ 10000
int     psCobylaNSaved_=0;
double  psCobylaSaveX_[psCobylaMaxSaved_*10];
double  psCobylaSaveY_[psCobylaMaxSaved_];
int     psNumOVars_=0;
int     psNumLVars_=0;
int     psCCurrDriver_=-1;
int     *psOVars_=NULL;
int     *psLVars_=NULL;
double  *psOWghts_=NULL;
double  *psLWghts_=NULL;
double  *psLVals_=NULL;
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
      int    ii, jj, kk, funcID, nOutputs, outputID, ignoreFlag, found;
      int    index;
      double *localY, ddata; 
      oData  *odata;
      FILE   *infile;

      if (psOptExpertMode_ != 0 && psCobylaNSaved_ == 0)
      {
         infile = fopen("psuade_cobyla_data","r");
         if (infile != NULL)
         {
            fscanf(infile, "%d %d", &psCobylaNSaved_, &kk);
            if ((psCobylaNSaved_ <= 0) || 
                (psCobylaNSaved_+1 > psCobylaMaxSaved_*10/nInputs))
            {
               printf("PSUADE Cobyla: history file has too much data.\n");
               printf("               Ask PSUADE developer for help.\n");
               fclose(infile);
               psCobylaNSaved_ = 0;
            }
            else if (kk != nInputs)
            {
               printf("PSUADE Cobyla: history file has invalid nInputs.\n");
               fclose(infile);
               psCobylaNSaved_ = 0;
            }
            else
            {
               for (ii = 0; ii < psCobylaNSaved_; ii++)
               {
                  fscanf(infile, "%d", &kk);
                  if (kk != ii+1)
                  {
                     printf("PSUADE Cobyla: data index mismatch.\n");
                     psCobylaNSaved_ = 0;
                     break;
                  }
                  for (jj = 0; jj < nInputs; jj++)
                     fscanf(infile, "%lg", &psCobylaSaveX_[ii*nInputs+jj]);
                  fscanf(infile, "%lg", &psCobylaSaveY_[ii]);
               }
               fclose(infile);
            }
         } 
      } 

      odata    = (oData *) data;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));
      outputID = odata->outputID_;

      found = 0;
      for (ii = 0; ii < psCobylaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psCobylaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            index = ii;
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
         odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
         funcID = odata->numFuncEvals_++;
         if (psNumOVars_ + psNumLVars_ > 0)
         {
            ddata = 0.0;
            for (ii = 0; ii < psNumOVars_; ii++)
            {
               kk = psOVars_[ii];
               ddata += psOWghts_[ii] * localY[kk];
            }
            (*YValue) = ddata;
            ddata = 0.0;
            for (ii = 0; ii < psNumLVars_; ii++)
            {
               kk = psLVars_[ii];
               ddata += pow(psLVals_[ii] - localY[kk], 2.0) * psLWghts_[ii];
            }
            (*YValue) += ddata;
         }
         else (*YValue) = localY[outputID];
      }
      else
      {
         localY[outputID] = psCobylaSaveY_[index];
         (*YValue) = localY[outputID];
      }

      if (psOptExpertMode_ != 0 && found == 0)
      {
         for (jj = 0; jj < nInputs; jj++)
            psCobylaSaveX_[psCobylaNSaved_*nInputs+jj] = XValues[jj];
         psCobylaSaveY_[psCobylaNSaved_] = localY[outputID];
         psCobylaNSaved_++;
      }

      ignoreFlag = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         constraints[ii] = XValues[ii] - odata->lowerBounds_[ii];
         constraints[ii+nInputs] = odata->upperBounds_[ii] - XValues[ii];
         if (constraints[ii] < 0.0 || constraints[ii+nInputs] < 0.0)
            ignoreFlag = 1;
      }
      if ((ignoreFlag == 0) && ((*YValue) < odata->optimalY_))
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
         if (odata->outputLevel_ > 2)
         {
            printf("CobylaOptimizer %6d : \n", odata->numFuncEvals_); 
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            printf("    Ymin  = %16.8e\n", odata->optimalY_);
         }
      }
      free(localY);

      if (psOptExpertMode_ != 0 && found == 0)
      {
         infile = fopen("psuade_cobyla_data","w");
         if (infile != NULL)
         {
            fprintf(infile, "%d %d\n", psCobylaNSaved_, nInputs);
            for (ii = 0; ii < psCobylaNSaved_; ii++)
            {
               fprintf(infile, "%d ", ii+1);
               for (jj = 0; jj < nInputs; jj++)
                  fprintf(infile, "%24.16e ", psCobylaSaveX_[ii*nInputs+jj]);
               fprintf(infile, "%24.16e\n", psCobylaSaveY_[ii]);
            }
            fclose(infile);
         }
      }
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
   int    nInputs, nConstraints, ii, maxfun;
   int    ntimes=3, printLevel=0, nOutputs, outputID, cmpFlag;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp;
   char   filename[500];
   string   iLine;
   ifstream iFile;

   nInputs = odata->nInputs_;
   nConstraints = 2 * nInputs;
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
      printf("                         tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }

   nOutputs = odata->nOutputs_;
   outputID = odata->outputID_;
   if (nOutputs > 1)
   {
      printAsterisks(PL_INFO, 0);
      printf("* You have multiple outputs in your problem definition.\n");
      printf("* Minimization will be performed with respect to output %d.\n",
             outputID+1);
      printf("* NOTE: However, by turning on the Opt EXPERT mode, \n");
      printf("* You may specialize the objective function by creating a\n");
      printf("* file called psuadeCobylaSpecial in your current directory.\n");
      printf("*\n");
      printf("* The objective function is in the form: \n");
      printf("*\n");
      printf("*  sum_{i=1}^m w_i O_i + sum_{j=1}^n (O_j - C_j)^2\n");
      printf("*\n");
      printf("* where\n");
      printf("*   m - number of outputs to be used to form linear sum.\n");
      printf("*   n - number of outputs to form Lagrange multiplier.\n");
      printf("*   w_i - weight of output i.\n");
      printf("*   C_j - constraint for output j.\n");
      printf("* The psuadeCobylaSpecial should have the following format: \n");
      printf("*\n");
      printf("\tPSUADE_BEGIN\n");
      printf("\t<m>          /* m in the above formula */\n");
      printf("\t1  <double>  /* w_1 if output 1 is used. */\n");
      printf("\t3  <double>  /* w_3 if output 3 is used (and O2 not used)*/\n");
      printf("\t...\n");
      printf("\t<n>          /* n in the above formula */\n");
      printf("\t2  <double>  /* C_2 if output 2 is used in Lagrange */\n");
      printf("\t...\n");
      printf("\tPSUADE_END\n");
      printDashes(PL_INFO, 0);
      if (psOptExpertMode_ != 0)
      {
         strcpy(filename, "psuadeCobylaSpecial");
         iFile.open(filename);
         if (iFile.is_open())
         {
            printf("INFO: psuadeCobylaSpecial file found.\n");
            getline (iFile, iLine);
            cmpFlag = iLine.compare("PSUADE_BEGIN");
            if (cmpFlag != 0)
            {
               printf("Cobyla ERROR: PSUADE_BEGIN not found.\n");
               iFile.close();
               return;
            }
            iFile >> psNumOVars_;
            if (psNumOVars_ <= 0)
            {
               printf("Cobyla ERROR: numOVars <= 0\n");
               iFile.close();
               return;
            }
            psOVars_  = new int[psNumOVars_];
            psOWghts_ = new double[psNumOVars_];
            for (ii = 0; ii < psNumOVars_; ii++)
            {
               iFile >> psOVars_[ii];
               if (psOVars_[ii] <= 0 || psOVars_[ii] > nOutputs)
               {
                  printf("Cobyla ERROR: invalid variable index %d.\n",
                         psOVars_[ii]);
                  iFile.close();
                  return;
               }
               psOVars_[ii]--;
               iFile >> psOWghts_[ii];
            }
            printf("%d outputs selected: \n", psNumOVars_);
            for (ii = 0; ii < psNumOVars_; ii++)
               printf("%4d   weight = %16.8e\n", psOVars_[ii]+1, psOWghts_[ii]);
            iFile >> psNumLVars_;
            if (psNumLVars_ < 0)
            {
               printf("Cobyla ERROR: numLVars < 0\n");
               iFile.close();
               return;
            }
            if (psNumLVars_ > 0)
            {
               psLVars_  = new int[psNumLVars_];
               psLWghts_ = new double[psNumLVars_];
               psLVals_  = new double[psNumLVars_];
            }
            for (ii = 0; ii < psNumLVars_; ii++)
            {
               iFile >> psLVars_[ii];
               if (psLVars_[ii] <= 0 || psLVars_[ii] > nOutputs)
               {
                  printf("Cobyla ERROR: invalid variable index %d.\n",
                         psLVars_[ii]);
                  iFile.close();
                  return;
               }
               psLVars_[ii]--;
               iFile >> psLVals_[ii];
               iFile >> psLWghts_[ii];
            }
            printf("%d Lagrange outputs selected: \n", psNumLVars_);
            for (ii = 0; ii < psNumLVars_; ii++)
               printf("%4d   value = %16.8e, weight = %16.8e\n",
                      psLVars_[ii]+1, psLVals_[ii], psLWghts_[ii]);
            getline (iFile, iLine);
            cmpFlag = iLine.compare("PSUADE_END");
            if (cmpFlag != 0)
            {
               getline (iFile, iLine);
               cmpFlag = iLine.compare("PSUADE_END");
               if (cmpFlag != 0)
               {
                  printf("Bobyqa ERROR: PSUADE_END not found.\n");
                  iFile.close();
                  return;
               }
            }
         }
      }
   }

   if ((odata->setOptDriver_ & 1))
   {
      printf("Cobyla: setting optimization simulation driver.\n");
      psCCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   printAsterisks(PL_INFO, 0);
   printf("Cobyla optimizer: max fevals = %d\n", odata->maxFEval_*ntimes);
   printf("Cobyla optimizer: tolerance  = %e\n", odata->tolerance_);
   if (printLevel > 1)
      printf("Cobyla optimizer: rho1, rho2 = %e %e\n",rhobeg,rhoend);
   if (psNumOVars_ + psNumLVars_ == 0)
      printf("Cobyla optimizer: selected output = %d\n", outputID+1);
   printEquals(PL_INFO, 0);

#ifdef HAVE_COBYLA
   int tfevals = 0;
   for (int kk = 0; kk < ntimes; kk++)
   {
      cobyla(nInputs, nConstraints, XValues, rhobeg, rhoend, 
             printLevel, &maxfun, evaluateFunction, 
             (void *) odata);
      tfevals += odata->numFuncEvals_;
   }
   odata->numFuncEvals_ = tfevals;
   printf("Cobyla optimizer: number of function evaluation = %d\n",tfevals);
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
   delete [] psOVars_;
   delete [] psOWghts_;
   delete [] psLVars_;
   delete [] psLWghts_;
   delete [] psLVals_;
   psLVals_ = NULL;
   psLWghts_ = NULL;
   psOWghts_ = NULL;
   psOVars_ = NULL;
   psLVars_ = NULL;
   psNumOVars_ = psNumLVars_ = 0;
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


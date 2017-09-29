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
// Functions for the class BobyqaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "BobyqaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);

#define psBobyqaMaxSaved_ 10000
int     psBobyqaNSaved_=0;
double  psBobyqaSaveX_[psBobyqaMaxSaved_*10];
double  psBobyqaSaveY_[psBobyqaMaxSaved_];
void    *psBobyqaObj_=NULL;
int     psNumBOVars_=0;
int     psNumBLVars_=0;
int     *psBOVars_=NULL;
int     *psBLVars_=NULL;
int     psBCurrDriver_ = -1;
double  *psBOWghts_=NULL;
double  *psBLWghts_=NULL;
double  *psBLVals_=NULL;
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
   void *bobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
      double *localY, ddata;
      char   pString[1000], lineIn[1000];
      oData  *odata;
      FILE   *infile;

      nInputs = (*nInps);
#if 0
      if (psOptExpertMode_ != 0 && psBobyqaNSaved_ == 0)
      {
         infile = fopen("psuade_bobyqa_history","r");
         if (infile != NULL)
         {
            fscanf(infile, "%s", pString);
            if (!strcmp(pString, "PSUADE_BEGIN"))
            {
               fscanf(infile, "%d", &kk);
               if (kk != *nInps)
               {
                  printf("Bobyqa: history file found but with mismatched nInputs.\n");
               }
               else
               {
                  fgets(lineIn,500,stdin);
                  sscanf(infile, "%s", pString);
                  if (!strcmp(pString, "PSUADE_END"))

             
            fscanf(infile, "%d %d", &psBobyqaNSaved_, &kk);
            if ((psBobyqaNSaved_ <= 0) ||
                (psBobyqaNSaved_+1 > psBobyqaMaxSaved_*10/nInputs) ||
                (psBobyqaNSaved_ > psBobyqaMaxSaved_))
            {
               printf("PSUADE Bobyqa: history file has too much data.\n");
               printf("               Ask PSUADE developer for help.\n"); 
               fclose(infile);
               psBobyqaNSaved_ = 0;
            }
            else if (kk != nInputs)
            {
               printf("PSUADE Bobyqa: history file has invalid input count.\n");
               fclose(infile);
               psBobyqaNSaved_ = 0;
            }
            else
            {
               for (ii = 0; ii < psBobyqaNSaved_; ii++)
               {
                  fscanf(infile, "%d", &kk);
                  if (kk != ii+1)
                  {
                     printf("PSUADE Bobyqa save: data index mismatch (%d).\n",ii+1);
                     psBobyqaNSaved_ = 0;
                     break;
                  }
                  for (jj = 0; jj < nInputs; jj++)
                     fscanf(infile, "%lg", &psBobyqaSaveX_[ii*nInputs+jj]);
                  fscanf(infile, "%lg", &psBobyqaSaveY_[ii]);
               }
               fclose(infile);
            }
         }
      }
#endif

      odata    = (oData *) psBobyqaObj_;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));
      outputID = odata->outputID_;

      found = 0;
      for (ii = 0; ii < psBobyqaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psBobyqaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            printf("Bobyqa: simulation results reuse.\n");
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
         odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
         funcID = odata->numFuncEvals_++;
         if (psNumBOVars_ + psNumBLVars_ > 0)
         {
            ddata = 0.0;
            for (ii = 0; ii < psNumBOVars_; ii++)
            {
               kk = psBOVars_[ii];
               ddata += psBOWghts_[ii] * localY[kk];
            }
            (*YValue) = ddata;
            ddata = 0.0;
            for (ii = 0; ii < psNumBLVars_; ii++)
            {
               kk = psBLVars_[ii];
               ddata += pow(psBLVals_[ii] - localY[kk], 2.0) * psBLWghts_[ii];
            }
            (*YValue) += ddata;
         }
         else (*YValue) = localY[outputID];
      }
      else
      {
         localY[outputID] = psBobyqaSaveY_[ii];
         (*YValue) = localY[outputID];
      }

      if ((psOptExpertMode_ != 0 && found == 0 &&
          (psBobyqaNSaved_+1)*nInputs < psBobyqaMaxSaved_*10) && 
          psBobyqaNSaved_ < psBobyqaMaxSaved_)
      {
         for (jj = 0; jj < nInputs; jj++)
            psBobyqaSaveX_[psBobyqaNSaved_*nInputs+jj] = XValues[jj];
         psBobyqaSaveY_[psBobyqaNSaved_] = (*YValue);
         psBobyqaNSaved_++;
      }

      if ((*YValue) < odata->optimalY_)
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
         if (psOptExpertMode_ == 1 && odata->outputLevel_ > 1)
         {
            printf("BobyqaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            printf("    Ymin  = %16.8e\n", odata->optimalY_);
            for (ii = 0; ii < nOutputs; ii++) 
               printf("    Y %6d = %16.8e\n", ii+1, localY[ii]);
         }
      }
      free(localY);
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
BobyqaOptimizer::BobyqaOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BobyqaOptimizer::~BobyqaOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void BobyqaOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0;
   int    nOutputs, outputID, cmpFlag, bobyqaFlag=7777;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *workArray;
   char   filename[500], cinput[500];;
   FILE   *infile=NULL;
   string   iLine;
   ifstream iFile;

   printLevel = odata->outputLevel_;
   nInputs = odata->nInputs_;
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
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
      printf("BobyqaOptimizer WARNING: tolerance too large.\n");
      printf("                         tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   nOutputs = odata->nOutputs_;
   outputID = odata->outputID_;
   if (nOutputs > 1)
   {
      printAsterisks(PL_INFO, 0);
      if (psOptExpertMode_ == 1)
      {
         printf("* You have multiple outputs in your problem definition.\n");
         printf("* Minimization will be performed with respect to output %d.\n",
                outputID+1);
         printf("* NOTE: However, by turning on the Opt EXPERT mode, \n");
         printf("* You may specialize the objective function by creating a\n");
         printf("* file called psuadeBobyqaSpecial in your current directory.\n");
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
         printf("* The psuadeBobyqaSpecial should have the following format: \n");
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
         printAsterisks(PL_INFO, 0);
      }
      if (psOptExpertMode_ == 1)
      {
         strcpy(filename, "psuadeBobyqaSpecial");
         iFile.open(filename);
         if (iFile.is_open())
         {
            printf("INFO: psuadeBobyqaSpecial file found ***.\n");
            getline (iFile, iLine);
            cmpFlag = iLine.compare("PSUADE_BEGIN");
            if (cmpFlag != 0)
            {
               printf("Bobyqa ERROR: PSUADE_BEGIN not found.\n");
               iFile.close();
               return;
            }
            iFile >> psNumBOVars_;
            if (psNumBOVars_ <= 0)
            {
               printf("Bobyqa ERROR: numOVars <= 0\n");
               iFile.close();
               return;
            }
            psBOVars_  = new int[psNumBOVars_];
            psBOWghts_ = new double[psNumBOVars_];
            for (ii = 0; ii < psNumBOVars_; ii++)
            {
               iFile >> psBOVars_[ii];
               if (psBOVars_[ii] <= 0 || psBOVars_[ii] > nOutputs)
               {
                  printf("Bobyqa ERROR: invalid variable index %d.\n",
                         psBOVars_[ii]);
                  iFile.close();
                  return;
               }
               psBOVars_[ii]--;
               iFile >> psBOWghts_[ii];
            }
            printf("%d outputs selected: \n", psNumBOVars_);
            for (ii = 0; ii < psNumBOVars_; ii++)
               printf("%4d   weight = %16.8e\n",
                      psBOVars_[ii]+1,psBOWghts_[ii]);
            iFile >> psNumBLVars_;
            if (psNumBLVars_ < 0)
            {
               printf("Bobyqa ERROR: numLVars < 0\n");
               iFile.close();
               return;
            }
            if (psNumBLVars_ > 0)
            {
               psBLVars_  = new int[psNumBLVars_];
               psBLWghts_ = new double[psNumBLVars_];
               psBLVals_  = new double[psNumBLVars_];
            }
            for (ii = 0; ii < psNumBLVars_; ii++)
            {
               iFile >> psBLVars_[ii];
               if (psBLVars_[ii] <= 0 || psBLVars_[ii] > nOutputs)
               {
                  printf("Bobyqa ERROR: invalid variable index %d.\n",
                         psBLVars_[ii]);
                  iFile.close();
                  return;
               }
               psBLVars_[ii]--;
               iFile >> psBLVals_[ii];
               iFile >> psBLWghts_[ii];
            }
            printf("%d Lagrange outputs selected: \n", psNumBLVars_);
            for (ii = 0; ii < psNumBLVars_; ii++)
               printf("%4d   value = %16.8e, weight = %16.8e\n", 
                      psBLVars_[ii]+1, psBLVals_[ii], psBLWghts_[ii]);
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

   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1))
   {
      printf("Bobyla: setting optimization simulation driver.\n");
      psBCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psBobyqaObj_= (void *) odata;
   printAsterisks(PL_INFO, 0);
   printf("Bobyqa optimizer: max fevals = %d\n", odata->maxFEval_);
   printf("Bobyqa optimizer: tolerance  = %e\n", odata->tolerance_);
   if (printLevel > 1)
      printf("Bobyqa optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
   if (psNumBOVars_ + psNumBLVars_ == 0)
      printf("Bobyqa optimizer: selected output = %d\n", outputID+1);
   printEquals(PL_INFO, 0);

   nPts = (nInputs + 1) * (nInputs + 2) / 2;
   workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];

#ifdef HAVE_BOBYQA
   if (psOptExpertMode_ != 0)
   {
      infile = fopen("psuade_bobyqa_history","r");
      if (infile != NULL)
      {
         psBobyqaNSaved_ = 0;
         while (feof(infile) == 0)
         {
            fscanf(infile, "%d %d", &ii, &kk);
            if (ii != 999 || kk != nInputs)
            {
               break;
            }
            else
            {
               for (ii = 0; ii < nInputs; ii++) 
                  fscanf(infile, "%lg",&psBobyqaSaveX_[psBobyqaNSaved_*nInputs+ii]);
               fscanf(infile, "%lg",&psBobyqaSaveY_[psBobyqaNSaved_]);
               psBobyqaNSaved_++;
            }
            if (((psBobyqaNSaved_+1)*nInputs > psBobyqaMaxSaved_*10) ||
                psBobyqaNSaved_ > psBobyqaMaxSaved_) break;
         } 
         fclose(infile);
      }
   }
   for (ii = 0; ii < nInputs; ii++) 
      printf("Bobyqa initial X %3d = %e\n", ii+1, XValues[ii]);
   bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_,
           odata->upperBounds_, &rhobeg, &rhoend, &bobyqaFlag, &maxfun, 
           workArray);
   printf("Bobyqa optimizer: total number of evaluations = %d\n",
           odata->numFuncEvals_);

   if (psOptExpertMode_ != 0 && psBobyqaNSaved_ > 0)
   {
      infile = fopen("psuade_bobyqa_history","w");
      if (infile != NULL)
      {
         for (ii = 0; ii < psBobyqaNSaved_; ii++)
         {
            fprintf(infile, "999 %d ", nInputs);
            for (kk = 0; kk < nInputs; kk++)
               fprintf(infile, "%24.16e ", psBobyqaSaveX_[ii*nInputs+kk]);
            fprintf(infile, "%24.16e\n", psBobyqaSaveY_[ii]);
         }
         fclose(infile);
      }
      printf("Bobyqa: history saved in psuade_bobyqa_history\n");
   }
#else
   printf("ERROR : Bobyqa optimizer not installed.\n");
   exit(1);
#endif

   if (odata->setOptDriver_ & 2)
   {
      printf("Bobyla: setting back to original simulation driver.\n");
      odata->funcIO_->setDriver(psBCurrDriver_);
   }
   delete [] XValues;
   delete [] workArray;
   if (psBOVars_  != NULL) delete [] psBOVars_;
   if (psBOWghts_ != NULL) delete [] psBOWghts_;
   if (psBLVars_  != NULL) delete [] psBLVars_;
   if (psBLWghts_ != NULL) delete [] psBLWghts_;
   if (psBLVals_  != NULL) delete [] psBLVals_;
   psBLVals_ = NULL;
   psBLWghts_ = NULL;
   psBOWghts_ = NULL;
   psBOVars_ = NULL;
   psBLVars_ = NULL;
   psNumBOVars_ = psNumBLVars_ = 0;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
BobyqaOptimizer& BobyqaOptimizer::operator=(const BobyqaOptimizer &)
{
   printf("BobyqaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


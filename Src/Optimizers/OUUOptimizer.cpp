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
// Functions for the class OUUOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "PDFManager.h"
#include "OUUOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Sampling.h"
#include "FuncApprox.h"
#include "PrintingTS.h"

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);
extern "C" void obobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);

// ************************************************************************
// Internal 'global' variables (for passing parameters to Fortran and C)
// ------------------------------------------------------------------------
void    *psOUUObj_=NULL;
int     psOUUM1_=-1;
int     psOUUM2_=-1;
int     psOUUM3_=-1;
int     psOUUM4_=-1;
int     psOUUZ3nSamples_=-1;
int     psOUUZ4nSamples_=-1;
int     psOUUUserOpt_ = 0;
int     psOUUUseRS_ = 0;
int     psOUUZ4RSType_ = 0;
int     psOUUZ4RSAux_ = 0;
int     psOUUValidateRS_ = 0;
double  *psOUUZ3SamInputs_=NULL;
double  *psOUUZ4SamInputs_=NULL;
double  *psOUUZ4LBounds_=NULL;
double  *psOUUZ4UBounds_=NULL;
double  *psOUUSamOutputs_=NULL;
double  *psOUUSamProbs_=NULL;
double  *psOUUXValues_=NULL;
double  *psOUUWValues_=NULL;
double  *psOUUOptimalX_=NULL;
double  psOUUOptimalY_=0.0;
FuncApprox *psOUUfaPtr_=NULL;
int     psOUULargeSampleSize_=0;
double  *psOUULargeSamInputs_=NULL;
double  *psOUULargeSamOutputs_=NULL;
int     psOUUMode_=1;
double  psOUUPercentile_=0.5;
double  psOUUStdevMultiplier_=0;
int     psOUUCounter_=0;
int     psOUUParallel_=0;
#define psOUUMaxSaved_ 10000
int     psOUUNSaved_=0;
double  psOUUSaveX_[psOUUMaxSaved_*10];
double  psOUUSaveY_[psOUUMaxSaved_*10];
int     psOUUEnsembleEval_=0;
int     *psOUUInputTypes_=NULL;
int     psOUUMasterMode_=0;
int     psOUUSaveHistory_=0;

#define PABS(x)  ((x) > 0 ? x : -(x))
#define psOUUType1 1
#define psOUUType2 2
#define psOUUType3 3
#define psOUUType4 4

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
   /* ****************************************************************** */
   /* This function will evaluate the simulation given M2 variables by   */
   /* the optimizer appended with the M3 variables fixed in ouu2evalfunc */
   /* -------------------------------------------------------------------*/
   void *ouuevalfunc2_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, kk, M1, M2, M3, M4, M, index, iOne=1, found, funcID;
      double ddata;
      oData  *odata;
      FILE   *fp=NULL;

      odata = (oData *) psOUUObj_;
      M1    = psOUUM1_;
      M2    = psOUUM2_;
      M3    = psOUUM3_;
      M4    = psOUUM4_;
      M     = M1 + M2 + M3 + M4;
      funcID = odata->numFuncEvals_;

      index = 0;
      for (ii = 0; ii < M; ii++) 
      {
         if (psOUUInputTypes_[ii] == psOUUType2)
         {
            psOUUXValues_[ii] = XValues[index];
            index++;
         }
      }
      index = 0;
      for (ii = 0; ii < M; ii++) 
      {
         if (psOUUInputTypes_[ii] == psOUUType3)
         {
            psOUUXValues_[ii] = psOUUWValues_[index];
            index++;
         }
      }
      for (ii = 0; ii < M; ii++) 
      {
         if (psOUUInputTypes_[ii] == psOUUType4)
         {
            psOUUXValues_[ii] = psOUUWValues_[index];
            index++;
         }
      }

      found = 0;
      for (ii = 0; ii < psOUUNSaved_; ii++)
      {
         for (kk = 0; kk < M; kk++)
            if (PABS(psOUUSaveX_[ii*M+kk]-psOUUXValues_[kk])>1.0e-14) 
               break;
         if (kk == M)
         {
            found = 1;
            if (odata->outputLevel_ > 2)
               printf("OUUOptimizer: simulation results reuse.\n");
            ddata = psOUUSaveY_[ii];
            break;
         }
      }

      if (found == 0)
      {
         odata->funcIO_->evaluate(funcID,M,psOUUXValues_,iOne,&ddata,0);
         odata->numFuncEvals_++;
         if (psOUUSaveHistory_ == 1 && (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
         {
            for (kk = 0; kk < M; kk++)
               psOUUSaveX_[psOUUNSaved_*M+kk] = psOUUXValues_[kk];
            psOUUSaveY_[psOUUNSaved_] = ddata;
            psOUUNSaved_++;
            fp = fopen("psuade_ouu_history","w");
            if (fp != NULL)
            {
               for (ii = 0; ii < psOUUNSaved_; ii++)
               {
                  fprintf(fp, "999 %d ", M);
                  for (kk = 0; kk < M; kk++)
                     fprintf(fp, "%24.16e ", psOUUSaveX_[ii*M+kk]);
                  fprintf(fp, "%24.16e\n", psOUUSaveY_[ii]);
               }
               fclose(fp);
            }
         }
      }
      (*YValue) = ddata;

      if (ddata < psOUUOptimalY_)
      {
         psOUUOptimalY_ = ddata;
         for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = psOUUXValues_[ii];
         if (odata->outputLevel_ > 2)
         {
            printf("    OUUOptimizer inner loop New Ymin = %16.8e (%d)\n",
                   ddata, odata->numFuncEvals_);
            for (ii = 0; ii < M; ii++) 
               if (psOUUInputTypes_[ii] == psOUUType2)
                  printf("     Input %3d = %e\n", ii+1, psOUUOptimalX_[ii]);
         }
      }
      return NULL;
   }

   /* ****************************************************************** */
   /* This function will evaluate the simulation with variables          */
   /* corresponding only to the design variables (stage 1)               */
   /* It runs a number of optimizations to find the best M2 variable     */
   /* values for each of the sample points corresponding to M3 variables */
   /* -------------------------------------------------------------------*/
   void *ouuevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, kk, ss, funcID, M, M1, M2, M3, M4, bobyqaFlag=1112, nPts;
      int    maxfun, iOne=1, *readys, status, index, nSamp;
      double rhobeg, rhoend, ddata, *workArray, mean, stdev, *XLocal;
      double *lowers, *uppers;
      char   winput[1000];
      oData  *odata;
      FILE   *fp=NULL;

      odata = (oData *) psOUUObj_;
      M1    = psOUUM1_;
      M2    = psOUUM2_;
      M3    = psOUUM3_;
      M4    = psOUUM4_;
      M     = M1 + M2 + M3 + M4;

      fp = fopen("psuade_ouu_stop","r");
      if (fp != NULL)
      {
         unlink("psuade_ouu_stop");
         printf("OUUOptimizer: psuade_ouu_stop file found.\n");
         printf("              Abrupt termination.\n");
         if (psOUUNSaved_ > 0)
         {
            fp = fopen("psuade_ouu_history","w");
            if (fp != NULL)
            {
               for (ii = 0; ii < psOUUNSaved_; ii++)
               {
                  fprintf(fp, "999 %d ", M);
                  for (kk = 0; kk < M; kk++)
                     fprintf(fp, "%24.16e ", psOUUSaveX_[ii*M+kk]);
                  fprintf(fp, "%24.16e\n", psOUUSaveY_[ii]);
               }
               fclose(fp);
            }
            printf("OUUOptimizer: history saved in psuade_ouu_history\n");
         }
         exit(1);
      }

      nSamp = psOUUZ3nSamples_ * psOUUZ4nSamples_;
      XLocal = (double *) malloc(M*sizeof(double));
      readys = (int *) malloc(nSamp*sizeof(int));

      index = 0;
      for (ii = 0; ii < M; ii++)
      {
         if (psOUUInputTypes_[ii] == psOUUType1)
         {
            psOUUXValues_[ii] = XValues[index];
            index++;
         }
      }
      printf("OUUOptimizer: Outer optimization iteration = %d\n",
             odata->numFuncEvals_/nSamp+1);
      if (odata->outputLevel_ > 1) 
      {
         printf("OUUOptimizer: Outer optimization loop FuncEval %d.\n",
                odata->numFuncEvals_+1);
         index = 0;
         for (ii = 0; ii < M; ii++)
         {
            if (psOUUInputTypes_[ii] == psOUUType1)
            {
               printf("    Current Level 1 input %3d = %e\n", ii+1, 
                      XValues[index]);
               index++;
            }
         }
      }
      if (psOUUEnsembleEval_ == 0)
      {
         for (ss = 0; ss < nSamp; ss++)
         {
            if (odata->outputLevel_ > 3) 
               printf("OUUOptimizer sample %d (of %d)\n",ss+1,nSamp); 
            readys[ss] = -1;
            psOUUOptimalY_ = PSUADE_UNDEFINED;
            index = 0;
            for (ii = 0; ii < M; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType3)
               {
                  kk = ss / psOUUZ4nSamples_;
                  psOUUWValues_[index] = psOUUZ3SamInputs_[kk*M3+index];
                  index++;
               }
            }
            index = 0;
            for (ii = 0; ii < M; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                  kk = ss % psOUUZ4nSamples_;
                  psOUUWValues_[index] = psOUUZ4SamInputs_[kk*M4+index];
                  index++;
               }
            }
            if (psOUUUserOpt_ == 0 && psOUUM2_ > 0)
            {
               lowers = (double *) malloc(M*sizeof(double));
               uppers = (double *) malloc(M*sizeof(double));
               rhobeg = 1.0e35;
               index = 0;
               for (ii = 0; ii < M; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType2)
                  {
                     lowers[index] = odata->lowerBounds_[ii];
                     uppers[index] = odata->upperBounds_[ii];
                     ddata = 0.5 * (lowers[index] + uppers[index]);
                     XLocal[index] = ddata;
                     ddata = uppers[index] - lowers[index];
                     if (ddata < rhobeg) rhobeg = ddata;
                     index++;
                  }
               }
               rhobeg *= 0.5;
               rhoend = rhobeg * odata->tolerance_;
               if (rhobeg < rhoend) rhoend = rhobeg * 1.0e-6;
               maxfun = odata->maxFEval_;
               nPts = (M2 + 1) * (M2 + 2) / 2;
               kk = (nPts + 5) * (nPts + M2) + 3 * M2 * (M2 + 5) / 2 + 1;
               workArray = (double *) malloc(kk * sizeof(double));
               if (odata->outputLevel_ > 3) 
                  printf("OUU: inner optimization begins (sample %d of %d)\n",
                         ss+1,nSamp);
               bobyqaFlag = 1112;
               obobyqa_(&M2, &nPts, XLocal, lowers, uppers, &rhobeg, &rhoend, 
                        &bobyqaFlag, &maxfun, workArray);
               psOUUSamOutputs_[ss] = psOUUOptimalY_;
               if (odata->outputLevel_ > 3) 
                  printf("OUU: inner optimization sample %d ends, best Y = %e\n",
                         ss+1,psOUUSamOutputs_[ss]); 
               free(workArray);
               free(lowers);
               free(uppers);
            }
            else
            {
               index = 0;
               for (ii = 0; ii < M; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType1)
                  {
                     XLocal[ii] = psOUUXValues_[index];
                     index++;
                  }
               }
               index = 0;
               for (ii = 0; ii < M; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType2)
                  {
                     XLocal[ii] = 0.5 * (odata->lowerBounds_[ii] + 
                                         odata->upperBounds_[ii]);
                     index++;
                  }
               }
               index = 0;
               for (ii = 0; ii < M; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType3)
                  {
                     kk = ss / psOUUZ4nSamples_;
                     XLocal[ii] = psOUUZ3SamInputs_[kk*M3+index];
                     index++;
                  }
               }
               index = 0;
               for (ii = 0; ii < M; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType4)
                  {
                     kk = ss % psOUUZ4nSamples_;
                     XLocal[ii] = psOUUZ4SamInputs_[kk*M4+index];
                     index++;
                  }
               }
               int found = 0;
               for (kk = 0; kk < psOUUNSaved_; kk++)
               {
                  for (ii = 0; ii < M; ii++)
                  {
                     if (psOUUInputTypes_[ii] != psOUUType2 &&
                         (PABS(psOUUSaveX_[kk*M+ii]-XLocal[ii])>1.0e-14)) 
                        break;
                  }
                  if (ii == M)
                  {
                     found = 1;
                     if (odata->outputLevel_ > 2)
                        printf("OUUOptimizer: simulation results reuse.\n");
                     for (ii = 0; ii < M; ii++)
                     {
                        if (psOUUInputTypes_[ii] == psOUUType2)
                           XLocal[ii] = psOUUSaveX_[kk*M+ii];
                     }   
                     psOUUSamOutputs_[ss] = psOUUSaveY_[kk];
                     for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = XLocal[ii];
                     readys[ss] = 0;
                     break;
                  }
               }
               if (found == 0)
               {
                  funcID = psOUUCounter_ * nSamp + ss;
                  readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,iOne,
                                                        &ddata,0);
                  psOUUSamOutputs_[ss] = ddata;
                  for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = XLocal[ii];
                  if (odata->outputLevel_ > 3 && readys[ss] == 0) 
                     printf("OUUOptimizer (2) sample %d completed, best Y = %e\n",
                         ss+1,psOUUSamOutputs_[ss]); 
                  if (psOUUSaveHistory_ == 1 && (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
                  {
                     for (ii = 0; ii < M; ii++)
                        psOUUSaveX_[psOUUNSaved_*M+ii] = XLocal[ii];
                     psOUUSaveY_[psOUUNSaved_] = psOUUSamOutputs_[ss];
                     psOUUNSaved_++;
                  }
                  if (psOUUParallel_ == 0) odata->numFuncEvals_++;
               }
            }
         } 

         if (psOUUUserOpt_ == 1 && psOUUParallel_ == 1)
         {
            for (ss = 0; ss < nSamp; ss++)
            {
               funcID = psOUUCounter_ * nSamp + ss;
               if (readys[ss] != 0)
               {
                  while (readys[ss] != 0)
                  {
#ifdef WINDOWS
                     Sleep(1000);
#else
                     sleep(1);
#endif
                     readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,
                                                           iOne,&ddata,2);
                  }
                  psOUUSamOutputs_[ss] = ddata;
                  if (odata->outputLevel_ > 3)
                     printf("OUUOptimizer (3) sample %d completed, best Y = %e\n",
                            ss+1,psOUUSamOutputs_[ss]); 
                  odata->numFuncEvals_++;
               }
            }
         }
      }
      else
      {
         index = 0;
         for (ii = 0; ii < M; ii++)
         {
            if (psOUUInputTypes_[ii] == psOUUType1)
            {
               XLocal[index] = psOUUXValues_[ii];
               index++;
            }
         }
         for (ss = 0; ss < nSamp; ss++)
         {
            index = 0;
            for (ii = 0; ii < M; ii++) 
            {
               if (psOUUInputTypes_[ii] == psOUUType1)
               {
                  psOUUXValues_[ss*M+ii] = XLocal[index];
                  index++;
               }
            }
            for (ii = 0; ii < M; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType2)
               {
                  psOUUXValues_[ss*M+ii] = 0.5*(odata->lowerBounds_[ii] + 
                                                odata->upperBounds_[ii]);
               }
            }
            index = 0;
            for (ii = 0; ii < M; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType3)
               {
                  kk = ss / psOUUZ4nSamples_;
                  psOUUXValues_[ss*M+ii] = psOUUZ3SamInputs_[kk*M3+index];
                  index++;
               }
            }
            index = 0;
            for (ii = 0; ii < M; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                  kk = ss % psOUUZ4nSamples_;
                  psOUUXValues_[ss*M+ii] = psOUUZ4SamInputs_[kk*M4+index];
                  index++;
               }
            }
         }
         funcID = odata->numFuncEvals_;
         odata->funcIO_->ensembleEvaluate(nSamp,M,psOUUXValues_,
                                          iOne,psOUUSamOutputs_,funcID);
         odata->numFuncEvals_ += nSamp;
         for (ss = 0; ss < nSamp; ss++)
         {
            if (psOUUSaveHistory_ == 1 && (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
            {
               for (kk = 0; kk < M; kk++)
                  psOUUSaveX_[psOUUNSaved_*M+kk] = psOUUXValues_[ss*M+kk];
               psOUUSaveY_[psOUUNSaved_] = psOUUSamOutputs_[ss];
               psOUUNSaved_++;
            }
         }
         for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = psOUUXValues_[ii];
      }

      int failCnt=0;
      for (ss = 0; ss < nSamp; ss++) 
         if (psOUUSamOutputs_[ss] >= 0.98*PSUADE_UNDEFINED) failCnt++;
      if (failCnt != 0)
         printf("WARNING: there are %d failed runs out of %d\n",failCnt,nSamp);
      if (failCnt == nSamp || (failCnt > 0 && psOUUSamProbs_ != NULL))
      {   
         printf("ERROR: X3 sample cannot admit any failures ==> terminate.\n");
         exit(1);
      }
      if (failCnt > 0 && psOUUMode_ > 2)
      {   
         printf("ERROR: OUU mode %d cannot admit any failures ==> terminate.\n",
                psOUUMode_);
         exit(1);
      }
      if (psOUUUseRS_ == 0 || nSamp == 1)
      {
         if (odata->outputLevel_ > 2) 
            printf("OUUOptimizer: computing objective (no RS), nFuncEval = %d\n",
                   odata->numFuncEvals_);
         if (psOUUMode_ == 1 || psOUUMode_ == 2)
         {
            mean = 0.0;
            if (psOUUZ3nSamples_ == 1)
            {
               for (ss = 0; ss < nSamp; ss++) 
                  if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED) 
                     mean += psOUUSamOutputs_[ss] / (double) (nSamp - failCnt);
            }
            else
            {
               for (ss = 0; ss < nSamp; ss++) 
               {
                  index = ss / psOUUZ4nSamples_;
                  mean += psOUUSamOutputs_[ss] / psOUUZ4nSamples_ * 
                          psOUUSamProbs_[index];
               }
            }
            (*YValue) = mean;
         }
         if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
         {
            stdev = 0.0;
            for (ss = 0; ss < nSamp; ss++) 
            {
               index = ss / psOUUZ4nSamples_;
               if (psOUUZ3nSamples_ == 1)
                  if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED) 
                     stdev += pow(psOUUSamOutputs_[ss]-mean,2.0)/
                              (double) (nSamp - failCnt);
               else
                  stdev += pow(psOUUSamOutputs_[ss]-mean, 2.0)*
                           psOUUSamProbs_[index] / (double) psOUUZ4nSamples_;
            }
            (*YValue) = mean + psOUUStdevMultiplier_ * sqrt(stdev);
         }
         if (psOUUMode_ == 3)
         {
            sortDbleList(nSamp, psOUUSamOutputs_);
            kk = (int) (psOUUPercentile_ * nSamp);
            if (kk >= nSamp) kk = nSamp - 1;
            (*YValue) = psOUUSamOutputs_[kk];
         }
         if (odata->outputLevel_ > 2) 
         {
            printf("OUUOptimizer: computed  objective (no RS) = %e.\n",
                   (*YValue));
         }
      }
      else
      {
         if (odata->outputLevel_ > 2) 
            printf("OUUOptimizer: computing objective with RS, nFuncEval = %d\n",
                   odata->numFuncEvals_);
         
         double totalMean=0.0, totalStdv=0.0, *resultStore, cverrors[3];
         resultStore=(double *) malloc(psOUUZ3nSamples_*sizeof(double));
         for (ii = 0; ii < psOUUZ3nSamples_; ii++)
         {
            if (psGMMode_ == 1)
            {
               fp = fopen("ouuZ4Sample", "w");
               if (fp != NULL)
               {
                  fprintf(fp,"%d %d 1\n",psOUUZ4nSamples_,M4);
                  for (ss = 0; ss < psOUUZ4nSamples_; ss++) 
                  {
                     for (kk = 0; kk < M4; kk++) 
                        fprintf(fp,"%24.16e ",psOUUZ4SamInputs_[ss*M4+kk]);
                     fprintf(fp,"%24.16e\n",psOUUSamOutputs_[ii*psOUUZ4nSamples_+ss]);
                  }
               }
               fclose(fp);
               printf("OUU INFO: a Z4 sample is ready for viewing in ouuZ4Sample.\n");
               printf("To continue, enter 0 (or 1 if no more interruption) : ");
               scanf("%d", &kk);
               if (kk == 1) psGMMode_ = 0;
            }
            status = psOUUfaPtr_->initialize(psOUUZ4SamInputs_,
                                    &psOUUSamOutputs_[ii*psOUUZ4nSamples_]);
            if (psOUUValidateRS_ == 1)
            {
               validate(psOUUZ4nSamples_,psOUUZ4SamInputs_,
                        &psOUUSamOutputs_[ii*psOUUZ4nSamples_], cverrors);
               printf("OUUOptimizer: Z3 sample %d (of %d)\n",ii+1,psOUUZ3nSamples_);
               printf("RS CV avg error = %e (scaled) \n", cverrors[0]);
               printf("RS CV rms error = %e (scaled) \n", cverrors[1]);
               printf("RS CV max error = %e (scaled) \n", cverrors[2]);
            }
            psOUUfaPtr_->evaluatePoint(psOUULargeSampleSize_,
                               psOUULargeSamInputs_,psOUULargeSamOutputs_);
            if (psOUUMode_ == 1 || psOUUMode_ == 2)
            {
               mean = 0.0;
               for (ss = 0; ss < psOUULargeSampleSize_; ss++) 
                  mean += psOUULargeSamOutputs_[ss];
               mean /= (double) psOUULargeSampleSize_;
               resultStore[ii] = mean;
            }
            if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
            {
               stdev = 0.0;
               for (ss = 0; ss < psOUULargeSampleSize_; ss++) 
                  stdev += pow(psOUULargeSamOutputs_[ss] - mean, 2.0);
               stdev /= (double) psOUULargeSampleSize_;
               resultStore[ii] = mean + psOUUStdevMultiplier_ * stdev;
            }
            if (psOUUMode_ == 3)
            {
               sortDbleList(psOUULargeSampleSize_, psOUULargeSamOutputs_);
               kk = (int) (psOUUPercentile_ * psOUULargeSampleSize_);
               if (kk >= psOUULargeSampleSize_) 
                  kk = psOUULargeSampleSize_ - 1;
               resultStore[ii] = psOUULargeSamOutputs_[kk];
            }
         }
         if (psOUUMode_ == 1 || psOUUMode_ == 2)
         {
            mean = 0.0;
            for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
            {
               if (psOUUSamProbs_ == NULL)
                  mean += resultStore[ii] / (double) psOUUZ3nSamples_;
               else
                  mean += resultStore[ii] * psOUUSamProbs_[ii];
            }
            (*YValue) = mean;
         }
         if (psOUUMode_ == 3)
         {
            mean = 0.0;
            for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
               mean += resultStore[ii] / (double) psOUUZ3nSamples_;
            (*YValue) = mean;
         }
         if (odata->outputLevel_ > 2) 
         {
            printf("OUUOptimizer: computed  objective (with RS) = %e.\n",
                   (*YValue));
         }
         free(resultStore);
      }
 
      if ((*YValue) < odata->optimalY_)
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < M; ii++) 
            odata->optimalX_[ii] = psOUUOptimalX_[ii];
         if (odata->outputLevel_ > 0)
         {
            printf("    OUUOptimizer outer loop new Ymin = %16.8e (***)\n",
                   (*YValue));
            if (psOUUUserOpt_ == 0)
            {
               for (ii = 0; ii < M1+M2; ii++) 
                  printf("        Input %3d at min = %e\n", ii+1, 
                         odata->optimalX_[ii]);
            }
         }
      }
      psOUUCounter_++;

      free(readys);
      free(XLocal);
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
OUUOptimizer::OUUOptimizer()
{
   psOUUM1_ = -1;
   psOUUM2_ = -1;
   psOUUM3_ = -1;
   psOUUZ3nSamples_ = -1;
   psOUUZ4nSamples_ = -1;
   psOUUUserOpt_  = 0;
   psOUUObj_ = NULL;
   psOUUZ3SamInputs_ = NULL;
   psOUUZ4SamInputs_ = NULL;
   psOUUZ4LBounds_ = NULL;
   psOUUZ4UBounds_ = NULL;
   psOUUSamOutputs_ = NULL;
   psOUUXValues_ = NULL;
   psOUUWValues_ = NULL;
   psOUUOptimalX_ = NULL;
   psOUUOptimalY_ = 0.0;
   psOUUUseRS_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
OUUOptimizer::~OUUOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void OUUOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0, bobyqaFlag=1111;
   int    M1, M2, M3, M4, index, iOne=1, iZero=0, printHeader=1;
   int    currDriver, count, *inputPDFs;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, ddata, *workArray;
   char   pString[1000], *cString, lineIn[20001], filename[1000];
   FILE   *fp=NULL;
   pData  pdata;

   nInputs = odata->nInputs_;
   printLevel = odata->outputLevel_;
   if (((psOUUM1_+psOUUM2_+psOUUM3_+psOUUM4_) == nInputs) && 
       psOUUM1_ > 0) printHeader = 0;
   if (printLevel >= 0 && printHeader == 1)
   {
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printf("*     ONE- OR TWO-STAGE OPTIMIZATION UNDER UNCERTAINTY (OUU)\n");
      printEquals(PL_INFO, 0);
      printf("This optimization capability solves the following problem:\n");
      printf("\n   minimize_Z1 { Phi_{Z3,Z4} [ G(Z1,Z2,Z3,Z4) ] } \n\n");
      printf("   subject to bound constraints on Z1, Z2, and Z4; and\n\n");
      printf("   Z3 is a set of discrete parameters for which a sample\n\n");
      printf("   is to be provided by the user. \n\n");
      printf("   (0) How to perform regular optimization? \n");
      printf("       In this case \n"); 
      printf("       - Z1 will be the optimization variables\n");
      printf("       - Z2 should be an empty set\n");
      printf("       - Z3 should be an empty set\n");
      printf("       - Z4 should be an empty set\n");
      printf("       - G(Z1,Z3,Z4)=G(Z1) is the simulation function (opt_driver)\n");
      printf("   (1) How to perform 1-level OUU? \n");
      printf("       In this case \n"); 
      printf("       - Z1 will be the optimization variables\n");
      printf("       - Z2 should be an empty set\n");
      printf("       - Z3 are parameters that you will provide a sample for.\n");
      printf("       - Z4 are parameters that you do not provide a sample.\n");
      printf("         (optionally, you can choose sampling schemed and\n");
      printf("          whether you want response surface for Z4).\n");
      printf("       - G(Z1,Z3,Z4) is the simulation function (opt_driver)\n");
      printf("\n");
      printf("   (2) How to perform 2-level OUU? \n");
      printf("       In this case \n"); 
      printf("       - Z1 will be the optimization variables\n");
      printf("       - Z2 is the set of second level optimization parameters.\n");
      printf("         If Z2 is entirely embedded in your level 2 optimization\n");
      printf("         driver (that is, their initial and final values do not\n");
      printf("         need to be published outside of the G function), then it\n");
      printf("         can be specified as an empty set.\n");
      printf("       - Z3 are parameters that you will provide a sample for.\n");
      printf("       - Z4 are parameters that you do not provide a sample.\n");
      printf("       - There are two options to how to set up G(Z1,Z2,Z3,Z4): \n");
      printf("         (a) It is a user-provided level 2 optimizer that does\n\n");
      printf("               G(Z1,Z2,Z3,Z4) = minimize_Z2 { F(Z1,Z2,Z3,Z4) }\n\n");
      printf("             where F(Z1,Z2,Z3,Z4) is embedded in G(Z1,Z2,Z3,Z4)\n\n");
      printf("             In this case, the user will provide G(Z1,Z2,Z3,Z4)\n");
      printf("             via opt_driver (that is, opt_driver = gfunction)\n");
      printf("         (b) It is again solving\n\n");
      printf("               G(Z1,Z2,Z3,Z4) = minimize_Z2 { F(Z1,Z2,Z3,Z4) }\n\n");
      printf("             but in this case, user is expected to provide the\n");
      printf("             F(Z1,Z2,Z3,Z4) function via 'opt_driver = ffunction'\n");
      printf("             and OUU will provide the optimization solver (BOBYQA).\n");
      printf("\n   In case 1 and 2, Phi_{Z3,Z4} is a functional on G(Z1,Z2,Z3,Z4)\n");
      printf("   with respect to Z3 and Z4. For example, Phi_{Z3,Z4} can be:\n");
      printf("   1. mean of G(Z1,Z2,Z3,Z4) with respect to Z3,Z4 (default)\n"); 
      printf("   2. mean of G(Z1,Z2,Z3,Z4) + alpha * std dev of G(Z1,Z2,Z3,Z4)\n"); 
      printf("   3. G(Z1,Z2,Z3*,Z4*) such that \n");
      printf("           Prob(G(Z1,Z2,Z3,Z4)>G(Z1,Z2,Z3*,Z4*)) = epsilon\n");
      printf("   4. min_{Z3,Z4} G(Z1,Z2,Z3,Z4) given Z1 and Z2 (robust opt.)\n");
      printEquals(PL_INFO, 0);
      printf("In the above formulation, the total number of parameters M = %d\n", 
             nInputs);
      printf("These parameters are to be divided into four groups:\n");
      printf("(1) Stage 1 optimization parameters Z1 (M1 >= 1) \n");
      printf("(2) Stage 2 optimization parameters Z2\n");
      printf("(3) uncertain parameters Z3 (with a user-provided sample) \n");
      printf("(4) uncertain parameters Z4 \n");
      printf("    - be continuous parameters (PSUADE to generate sample), or\n");
      printf("    - a large sample (PSUADE to create RS using a small subset\n");
      printf("Thus, the first M1+M2 parameters are considered to be\n");
      printf("optimization parameters, and M3+M4 are uncertain parameters\n");
      printf("so that M = M1 + M2 + M3 + M4.\n");
      printEquals(PL_INFO, 0);
      printf("To reuse the simulations (e.g. restart due to abrupt termination),\n");
      printf("turn on save_history and use_history optimization options in the\n");
      printf("ANALYSIS section (e.g. optimization save_history). You will see\n");
      printf("a file created called 'psuade_bobyqa_history' afterward.\n");
      printAsterisks(PL_INFO, 0);
      printf("IF YOU ARE READY TO MOVE ON, ENTER 'y' AND RETURN : ");
      lineIn[0] = '0';
      while (lineIn[0] != 'y' && lineIn[0] != 'n')
      {
         scanf("%s", lineIn);
         if (lineIn[0] == 'n') 
         {
            odata->optimalY_ = 0;
            for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
            printf("OUUOptimizer INFO: abort.\n");
            exit(1);
         }
      }
      printEquals(PL_INFO, 0);
   }
   if (psMasterMode_ != 0)
   {
      printf("OUUOptimizer INFO: Master mode to be turned off.\n");
      psMasterMode_ = 0;
      psOUUMasterMode_ = 1;
   }
   if ((psOUUM1_+psOUUM2_+psOUUM3_+psOUUM4_) == nInputs && (psOUUM1_ > 0))
   {
      M1 = psOUUM1_;
      M2 = psOUUM2_;
      M3 = psOUUM3_ - psOUUM4_;
      M4 = psOUUM4_;
   }
   else
   {
      M1 = M2 = M3 = M4 = 0;
      printf("M1 = number of design (level 1 optimization) parameters\n");
      while (M1 <= 0 || M1 > nInputs)
      {
         printf("Enter M1 (between 1 and %d) : ", nInputs);
         scanf("%d", &M1);
      }
      if (M1 < nInputs)
      {
         printf("M2 = number of recourse (level 2 optimization) parameters\n");
         M2 = -1;
         while (M2 < 0 || M1+M2 > nInputs)
         {
            printf("Enter M2 (between 0 and %d) : ", nInputs-M1);
            scanf("%d", &M2);
         }
         if ((M1 + M2) < nInputs)
         {
            printf("M3 = number of discrete (scenario) parameters\n");
            M3 = -1;
            while (M3 < 0 || M1+M2+M3 > nInputs)
            {
               printf("Enter M3 (between 0 and %d) : ", nInputs-M1-M2);
               scanf("%d", &M3);
            }
            M4 = nInputs - M1 - M2 - M3;
         }
      }
      psOUUM1_ = M1;
      psOUUM2_ = M2;
      psOUUM3_ = M3;
      psOUUM4_ = M4;
   }
   if (printLevel >= 0)
   {
      printDashes(PL_INFO, 0);
      printf("Number of first  stage optimization parameters = %d\n", M1);
      printf("Number of second stage optimization parameters = %d\n", M2);
      printf("Number of discrete   uncertain parameters      = %d\n", M3);
      printf("Number of continuous uncertain parameters      = %d\n", M4);
      printDashes(PL_INFO, 0);
   }

   if (M1 == nInputs)
   {
      psOUUInputTypes_ = new int[nInputs];
      for (ii = 0; ii < nInputs; ii++) psOUUInputTypes_[ii] = 1;
   }
   else
   {
      printf("In the following, please select type for each variable:\n");
      printf("  1. design variable (level 1 optimization parameter)\n");
      printf("  2. operating variable (level 2 optimization parameter)\n");
      printf("  3. discrete uncertain variable (a sample will be given)\n");
      printf("  4. continuous uncertain variable\n");
      printf("NOTE: make sure your specification matches the data above.\n");
      fgets(lineIn, 5000, stdin);
      psOUUInputTypes_ = new int[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         sprintf(pString, "Type for variable %d ? ", ii+1);
         psOUUInputTypes_[ii] = getInt(1,4,pString); 
      }
   }
   printEquals(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++)
   {
      if (psOUUInputTypes_[ii] == psOUUType1) 
         printf("Input %4d is a design parameter.\n",ii+1);
      else if (psOUUInputTypes_[ii] == psOUUType2) 
         printf("Input %4d is an operating parameter.\n",ii+1);
      else if (psOUUInputTypes_[ii] == psOUUType3) 
         printf("Input %4d is a discrete uncertain parameter.\n",ii+1);
      else if (psOUUInputTypes_[ii] == psOUUType4) 
         printf("Input %4d is a continuous uncertain parameter.\n",ii+1);
   }
   int c1=0, c2=0, c3=0, c4=0;
   for (ii = 0; ii < nInputs; ii++)
   {
      if (psOUUInputTypes_[ii] == psOUUType1) c1++;
      if (psOUUInputTypes_[ii] == psOUUType2) c2++;
      if (psOUUInputTypes_[ii] == psOUUType3) c3++;
      if (psOUUInputTypes_[ii] == psOUUType4) c4++;
   }
   if (c1 != M1 || c2 != M2 || c3 != M3 || c4 != M4)
   {
      printf("OUUOptimizer ERROR: input type counts do not match.\n");
      printf("    Number of type 1 = %d (expected = %d)\n",c1,M1);
      printf("    Number of type 2 = %d (expected = %d)\n",c2,M2);
      printf("    Number of type 3 = %d (expected = %d)\n",c3,M3);
      printf("    Number of type 4 = %d (expected = %d)\n",c4,M4);
      exit(1);
   }
   if (M4 > 0)
   {
      odata->psIO_->getParameter("input_pdfs", pdata);
      inputPDFs = pdata.intArray_;
      pdata.intArray_ = NULL;
      kk = index = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (psOUUInputTypes_[ii] == psOUUType4)
         {
            if (inputPDFs[ii] == PSUADE_PDF_UNIFORM) 
               strcpy(pString, "Uniform");            
            if (inputPDFs[ii] == PSUADE_PDF_NORMAL) 
               strcpy(pString, "Normal");            
            if (inputPDFs[ii] == PSUADE_PDF_LOGNORMAL) 
               strcpy(pString, "Lognormal");            
            if (inputPDFs[ii] == PSUADE_PDF_TRIANGLE) 
               strcpy(pString, "Triangle");            
            if (inputPDFs[ii] == PSUADE_PDF_BETA) 
               strcpy(pString, "Beta");            
            if (inputPDFs[ii] == PSUADE_PDF_WEIBULL) 
               strcpy(pString, "Weibull");            
            if (inputPDFs[ii] == PSUADE_PDF_GAMMA) 
               strcpy(pString, "Gamma");            
            if (inputPDFs[ii] == PSUADE_PDF_EXPONENTIAL) 
               strcpy(pString, "Exponential");            
            if (inputPDFs[ii] == PSUADE_PDF_USER) 
               strcpy(pString, "User");            
            if (inputPDFs[ii] == PSUADE_PDF_F) 
               strcpy(pString, "F");            
            printf("PDF type for Input %5d = %s\n",ii+1,pString);
         }
      }
      delete [] inputPDFs;
   }
   printEquals(PL_INFO, 0);

   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   odata->optimalY_ = 1.0e50;
   XValues = new double[M1+1];
   index = 0;
   for (ii = 0; ii < nInputs; ii++)
   {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
         XValues[index] = odata->initialX_[ii];
         index++;
      }
   }
   rhobeg = 1e35;
   index = 0;
   for (ii = 0; ii < nInputs; ii++) 
   {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
         ddata = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
         if (ddata < rhobeg) rhobeg = ddata;
         index++;
      }
   }
   rhobeg *= 0.5;
   rhoend = rhobeg * odata->tolerance_;
   if (rhobeg < rhoend)
   {
      printf("OUUOptimizer WARNING: tolerance too large.\n");
      printf("                      tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }

   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1))
   {
      printf("OUUOptimizer: setting optimization simulation driver.\n");
      currDriver = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psOUUObj_= (void *) odata;
   printAsterisks(PL_INFO, 0);
   printf("OUUOptimizer: max fevals = %d\n", odata->maxFEval_);
   printf("OUUOptimizer: tolerance  = %e\n", odata->tolerance_);
   printEquals(PL_INFO, 0);

   psOUUUserOpt_ = 0;
   psOUUUseRS_ = 0;
   if ((psOptExpertMode_ == 1) && (nInputs != (M1+M2)))
   {
      printEquals(PL_INFO, 0);
      printf("Select which functional Phi_{Z3,Z4} to use: \n");
      printf("  1. mean of G(Z1,Z2,Z3,Z4) with respect to Z3,Z4 (default)\n"); 
      printf("  2. mean of G(Z1,Z2,Z3,Z4) + beta * std dev of G(Z1,Z2,Z3,Z4)\n"); 
      printf("  3. G(Z1,Z2,Z3*,Z4*) such that \n");
      printf("           Prob(G(Z1,Z2,Z3,Z4)>G(Z1,Z2,Z3*,Z4*)) = 1 - alpha\n");
      printf("     This is also called value-at-risk with confidence level alpha.\n");
      sprintf(pString,"Enter your choice of functional (1, 2 or 3) : ");
      psOUUMode_ = getInt(1, 3, pString);
      if (psOUUMode_ == 2)
      {
         sprintf(pString,"Enter beta (>= 0) : ");
         psOUUStdevMultiplier_ = -1;
         while (psOUUStdevMultiplier_ < 0)
         { 
            psOUUStdevMultiplier_ = getDouble(pString);
         } 
      } 
      if (psOUUMode_ == 3)
      {
         psOUUPercentile_ = 0.0;
         sprintf(pString,"Enter the confidence interval : [0.5 - 1.0] ");
         while ((psOUUPercentile_ < 0.5) || (psOUUPercentile_ > 1.0)) 
         { 
            psOUUPercentile_ = getDouble(pString);
         } 
      }
   }

   psOUUZ3nSamples_ = 1;
   if (M3 > 0)
   {
      printEquals(PL_INFO, 0);
      printf("A sample for Z3 is needed from you. Data format should be :\n");
      printf("line 1: <nSamples> <nInputs> \n");
      printf("line 2: <sample 1 input 1> <input 2> ... <probability>\n");
      printf("line 3: <sample 2 input 1> <input 2> ... <probability>\n");
      printf("...\n");
      printf("Enter user sample file name : ");
      scanf("%s", filename);
      fgets(lineIn, 5000, stdin);
      fp = fopen(filename, "r");
      if (fp == NULL)
      {
         printf("OUUOptimizer ERROR: user sample file %s not found.\n",
                 filename);
            exit(1);
      }
      fgets(lineIn, 10000, fp);
      sscanf(lineIn, "%d %d", &psOUUZ3nSamples_, &ii);
      if (ii != M3)
      {
         printf("OUUOptimizer ERROR: user sample nInputs %d != %d\n",ii, M3);
         fclose(fp);
         exit(1);
      } 
      if (psOUUZ3nSamples_ < M3+1)
      {
         printf("OUUOptimizer ERROR: user sample size should be >= %d\n",M3+1);
         fclose(fp);
         exit(1);
      } 
      psOUUZ3SamInputs_ = new double[psOUUZ3nSamples_ * M3];
      psOUUSamProbs_    = new double[psOUUZ3nSamples_];
      ddata = 0.0;
      if (odata->outputLevel_ > 4)
         printf("Reading sample for X3.\n");
      for (ii = 0; ii < psOUUZ3nSamples_; ii++)
      {
         fgets(lineIn, 20000, fp);
         index = kk = 0;
         while (kk < M3)
         {
            while (lineIn[index] == ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               printf("ERROR: reading sample file %s line %d.\n",filename,ii+2);
               fclose(fp);
               exit(1);
            }
            sscanf(&lineIn[index],"%lg",&psOUUZ3SamInputs_[ii*M3+kk]);
            while (lineIn[index] != ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               printf("ERROR: reading sample file %s line %d.\n",filename,ii+2);
               fclose(fp);
               exit(1);
            }
            if (odata->outputLevel_ > 4)
               printf("%12.4e ",psOUUZ3SamInputs_[ii*M3+kk]);
            kk++;
         }
         while (lineIn[index] == ' ') index++;
         if (lineIn[index] == '\0' || lineIn[index] == '\n')
         {
            printf("WARNING: reading probability at line %d of %s.\n",ii+2,
                   filename);
            printf("         Will set probabilities to be equiprobable.\n");
            psOUUSamProbs_[ii] = PSUADE_UNDEFINED;
            ddata += 1.0;
         }
         else sscanf(&lineIn[index],"%lg",&psOUUSamProbs_[ii]);
         if (odata->outputLevel_ > 4)
            printf("%12.4e\n",psOUUSamProbs_[ii]);
      }
      fclose(fp);
      if (ddata != 0)
      {
         for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
            psOUUSamProbs_[ii] = 1.0 / (double) psOUUZ3nSamples_;
      }
      else
      {
         ddata = 0.0;
         for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
            ddata += psOUUSamProbs_[ii];
      }
      printf("User sample for Z3 has %d points\n", psOUUZ3nSamples_);
      printf("User sample for Z3 CDF = %e (should be ~1)\n", ddata);
   }

   int    methodZ4=1,*samStates,samplingOption,*iPdfs,index2;
   double *samOuts,*lowers,*uppers,*inputMeans=NULL,*inputStdevs=NULL;
   double *iMeans, *iStdvs;
   char   **targv;
   Sampling   *sampler = NULL;
   PDFManager *pdfman=NULL;
   Vector     vecLB, vecUB, vecOut;
   Matrix     *corMat1, corMat2;
   psOUUZ4nSamples_ = 1;
   if (M4 > 0)
   {
      printf("OUUOptimizer: generating a sample for Z4. Two options:\n");
      printf("(1) Users can upload a sample to PSUADE\n");
      printf("(2) PSUADE can internally create a sample\n");
      sprintf(pString, "Select option 1 or 2 : ");
      kk = getInt(1,2,pString);
      if (kk == 1)
      {
         printEquals(PL_INFO, 0);
         printf("A Z4 sample is needed from you. The file format should be:\n");
         printf("line 1: <nSamples> <nInputs> \n");
         printf("line 2: <sample 1 input 1> <input 2> \n");
         printf("line 3: <sample 2 input 1> <input 2> \n");
         printf("...\n");
         printf("Enter user sample file name : ");
         scanf("%s", filename);
         fgets(lineIn, 5000, stdin);
         fp = fopen(filename, "r");
         if (fp == NULL)
         {
            printf("OUUOptimizer ERROR: user sample file %s not found.\n",
                    filename);
            exit(1);
         }
         fgets(lineIn, 5000, fp);
         sscanf(lineIn, "%d %d", &psOUUZ4nSamples_, &ii);
         if (ii != M4)
         {
            printf("OUUOptimizer ERROR: user sample nInputs %d != %d\n",ii, M4);
            fclose(fp);
            exit(1);
         }
         if (psOUUZ4nSamples_ < M4+1)
         {
            printf("OUUOptimizer ERROR: user sample size should be >= %d\n",M4+1);
            fclose(fp);
            exit(1);
         }
         psOUUZ4SamInputs_ = new double[psOUUZ4nSamples_ * M4];
         ddata = 0.0;
         for (ii = 0; ii < psOUUZ4nSamples_; ii++)
         {
            for (kk = 0; kk < M4; kk++)
               fscanf(fp,"%lg",&psOUUZ4SamInputs_[ii*M4+kk]);
         }
         fclose(fp);
         printf("The user sample for Z4 has %d points\n", psOUUZ4nSamples_);
         printf("You have the option to select a subset of Z4 for building\n");
         printf("response surfaces and use the original larger Z4 sample\n");
         printf("to estimate the statistics from the response surface.\n");
         printf("Use response surface for Z4 to compute statistics? (y or n) ");
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUUseRS_ = 1;
         if (psOUUUseRS_ == 1)
         {
            fgets(lineIn, 5000, stdin);
            psOUULargeSampleSize_ = psOUUZ4nSamples_;
            psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
            psOUULargeSamOutputs_ = new double[psOUULargeSampleSize_];
            for (ii = 0; ii < psOUULargeSampleSize_*M4; ii++)
               psOUULargeSamInputs_[ii] = psOUUZ4SamInputs_[ii]; 
            for (ii = 0; ii < psOUULargeSampleSize_; ii++)
               psOUULargeSamOutputs_[ii] = 0.0;
            printf("Your Z4 sample size is %d.\n", psOUUZ4nSamples_);
            printf("This sample size may be too large for building a RS.\n");
            sprintf(pString,"How many points to use for building RS? (%d - %d) ",
                    M4+1, psOUUZ4nSamples_);
            psOUUZ4nSamples_ = getInt(M4+1, psOUUZ4nSamples_, pString);
            printf("You have 2 options on how to generate this RS set:\n");
            printf("(1) You upload another Z4 sample of size %d\n",psOUUZ4nSamples_);
            printf("(2) PSUADE automatically selects another Z4 sample of size %d\n",
                   psOUUZ4nSamples_);
            sprintf(pString,"Select option 1 or 2 : ");
            kk = getInt(1, 2, pString);
            if (kk == 2)
            {
               printf("OUU will select a random subset of sample points from\n");
               printf("your original Z4 sample for building response surface.\n");
               count = 0;
               for (ii = 0; ii < psOUULargeSampleSize_; ii++)
               {
                  if (psOUUZ4nSamples_ == psOUULargeSampleSize_) index = ii;
                  else index = PSUADE_rand() % (psOUULargeSampleSize_-ii) + ii; 
                  for (kk = 0; kk < M4; kk++)
                     psOUUZ4SamInputs_[count*M4+kk] = 
                                  psOUULargeSamInputs_[index*M4+kk];
                  count++;
                  if (count == psOUUZ4nSamples_) break;
               }
            }
            else if (kk == 1)
            {
               printEquals(PL_INFO, 0);
               printf("A Z4 RS sample is needed from you. The file format should be:\n");
               printf("line 1: <nSamples> <nInputs> \n");
               printf("line 2: <sample 1 input 1> <input 2> \n");
               printf("line 3: <sample 2 input 1> <input 2> \n");
               printf("...\n");
               printf("Enter user sample file name : ");
               scanf("%s", filename);
               fgets(lineIn, 5000, stdin);
               fp = fopen(filename, "r");
               if (fp == NULL)
               {
                  printf("OUUOptimizer ERROR: user sample file %s not found.\n",
                         filename);
                  exit(1);
               }
               fgets(lineIn, 5000, fp);
               sscanf(lineIn, "%d %d", &kk, &ii);
               if (ii != M4 || kk != psOUUZ4nSamples_)
               {
                  printf("OUUOptimizer ERROR: parameter mismatch.\n");
                  fclose(fp);
                  exit(1);
               }
               for (ii = 0; ii < psOUUZ4nSamples_; ii++)
               {
                  for (kk = 0; kk < M4; kk++)
                     fscanf(fp,"%lg",&psOUUZ4SamInputs_[ii*M4+kk]);
               }
               fclose(fp);
            }
         }
      }
      else
      {
         psOUUZ4nSamples_ = 100;
         printEquals(PL_INFO, 0);
         printf("PSUADE will create a sample (size will be selected later)\n");
         printf("for Z4. Since the sample should be small for computational\n");
         printf("efficiency, and accuracy of statistics depends on sample\n");
         printf("size, you have the option to add one more step in OUU by\n");
         printf("estimating the statistics with a large sample evaluated on\n");
         printf("the response surfaces built from the small sample.\n");
         printf("Use response surface for Z4 to compute statistics? (y or n) ");
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUUseRS_ = 1;
         if (psOUUMasterMode_ == 1 && psOUUUseRS_ == 1)
         {
            printf("For diagnostics purposes, you may choose to validate\n");
            printf("the response surface constructed at each OUU iteration.\n");
            printf("The more accurate the cross validation result is, the\n");
            printf("more accurate is the statistics.\n"); 
            printf("Validate the response surfaces for Z4? (y or n) ");
            scanf("%s", lineIn);
            if (lineIn[0] == 'y') psOUUValidateRS_ = 1;
         }
         fgets(lineIn, 5000, stdin);
         samplingOption = 0;
         if (psOUUUseRS_ == 0)
         {
            odata->psIO_->getParameter("input_pdfs", pdata);
            inputPDFs = pdata.intArray_;
            pdata.intArray_ = NULL;
            kk = index = 0;
            if (inputPDFs != NULL)
            {
               for (ii = 0; ii < nInputs; ii++) 
               {
                  if (psOUUInputTypes_[ii] == psOUUType4)
                  {
                     kk += inputPDFs[ii];
                     index++;
                  }
               }
            }
            if (kk > 0) samplingOption = 1;
         }
         if (samplingOption == 0)
         {
            printEquals(PL_INFO, 0);
            printf("OUUOptimizer uses a Z4 sample to estimate the objective\n");
            printf("Default sampling method = Latin hypercube\n");
            printf("Default sample size     = %d\n",psOUUZ4nSamples_);
            printf("Available sampling method: \n");
            printf("   (1) LHS, \n");
            printf("   (2) factorial, or\n");
            printf("   (3) quasi-MC.\n");
            sprintf(pString,"Select sampling method (1, 2 or 3) : ");
            methodZ4 = getInt(1, 3, pString);
            if (methodZ4 == 1 || methodZ4 == 3)
            {
               kk = 2;
               if (psOUUUseRS_ == 1) kk = M4 + 1;
               sprintf(pString, "Enter sample size (>= %d, <= 1000) : ", kk);
               psOUUZ4nSamples_ = getInt(2, 10000, pString);
               printf("Latin hypercube/QMC has sample size = %d\n",
                      psOUUZ4nSamples_);  
            }
            else if (methodZ4 == 2)
            {
               sprintf(pString,
                       "Enter number of levels per variable (>=2, <=100) : ");
               psOUUZ4nSamples_ = getInt(2, 100, pString);
               kk = psOUUZ4nSamples_;
               for (ii = 1; ii < M4; ii++) psOUUZ4nSamples_ *= kk;
               printf("Factorial design has sample size = %d\n", 
                      psOUUZ4nSamples_);  
            }
            if (methodZ4 == 1)
                 sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
            else if (methodZ4 == 2)
                 sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
            else sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
            sampler->setPrintLevel(0);
            lowers = new double[M4];
            uppers = new double[M4];
            index = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                  lowers[index] = odata->lowerBounds_[ii];
                  uppers[index] = odata->upperBounds_[ii];
                  index++;
               }
            }
            sampler->setInputBounds(M4, lowers, uppers);
            sampler->setOutputParams(iOne);
            sampler->setSamplingParams(psOUUZ4nSamples_, iOne, iZero);
            sampler->initialize(0);
            psOUUZ4nSamples_ = sampler->getNumSamples();
            psOUUZ4SamInputs_ = new double[psOUUZ4nSamples_ * M4];
            samStates = new int[psOUUZ4nSamples_];
            samOuts = new double[psOUUZ4nSamples_];
            sampler->getSamples(psOUUZ4nSamples_, M4, iOne, psOUUZ4SamInputs_, 
                                samOuts, samStates);
            delete sampler;
            delete [] lowers;
            delete [] uppers;
            delete [] samStates;
            delete [] samOuts;
            printEquals(PL_INFO, 0);
         }
         else
         {
            printEquals(PL_INFO, 0);
            printf("OUUOptimizer uses a Z4 sample to estimate the objective\n");
            printf("Default sample size     = %d\n",psOUUZ4nSamples_);
            sprintf(pString,"Enter your desired sample size (>=10, <=1000) : ");
            psOUUZ4nSamples_ = getInt(10, 10000, pString);
            odata->psIO_->getParameter("input_means", pdata);
            inputMeans = pdata.dbleArray_;
            if (inputMeans == NULL)
            {
               inputMeans = new double[nInputs];
               for (ii = 0; ii < nInputs; ii++) inputMeans[ii] = 0;
            }
            pdata.dbleArray_ = NULL;
            odata->psIO_->getParameter("input_stdevs", pdata);
            inputStdevs = pdata.dbleArray_;
            if (inputStdevs == NULL)
            {
               inputStdevs = new double[nInputs];
               for (ii = 0; ii < nInputs; ii++) inputStdevs[ii] = 1;
            }
            pdata.dbleArray_ = NULL;
            odata->psIO_->getParameter("input_cor_matrix", pdata);
            corMat1 = (Matrix *) pdata.psObject_;
            pdata.psObject_ = NULL;

            corMat2.setDim(M4,M4);
            iPdfs  = new int[nInputs];
            iMeans = new double[nInputs];
            iStdvs = new double[nInputs];
            lowers = new double[nInputs];
            uppers = new double[nInputs];
            index = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               { 
                  iMeans[index] = inputMeans[ii];
                  iStdvs[index] = inputStdevs[ii];
                  iPdfs[index]  = inputPDFs[ii];
                  lowers[index] = odata->lowerBounds_[ii];
                  uppers[index] = odata->upperBounds_[ii];
                  index2 = 0;
                  for (kk = 0; kk < nInputs; kk++)
                  {
                     if (psOUUInputTypes_[kk] == psOUUType4)
                     {
                        ddata = corMat1->getEntry(ii,kk);
                        corMat2.setEntry(index,index2,ddata);
                        index2++;
                     }
                  }
                  index++;
               }
            }
            pdfman = new PDFManager();
            pdfman->initialize(M4,iPdfs,iMeans,iStdvs,corMat2,NULL,NULL);
            vecLB.load(M4, lowers);
            vecUB.load(M4, uppers);
            vecOut.setLength(psOUUZ4nSamples_*M4);
            pdfman->genSample(psOUUZ4nSamples_, vecOut, vecLB, vecUB);
            psOUUZ4SamInputs_= new double[psOUUZ4nSamples_*M4];
            psOUUSamOutputs_ = new double[psOUUZ4nSamples_];
            for (ii = 0; ii < psOUUZ4nSamples_*M4; ii++)
               psOUUZ4SamInputs_[ii] = vecOut[ii];
            delete [] inputPDFs;
            delete [] inputMeans;
            delete [] inputStdevs;
            delete [] iPdfs;
            delete [] iMeans;
            delete [] iStdvs;
            delete [] lowers;
            delete [] uppers;
            delete pdfman;
         }
         if (psOUUUseRS_ == 1)
         {
            odata->psIO_->getParameter("input_pdfs", pdata);
            inputPDFs = pdata.intArray_;
            pdata.intArray_ = NULL;
            if (inputPDFs == NULL)
            {
               inputPDFs = new int[nInputs];
               for (ii = 0; ii < nInputs; ii++) inputPDFs[ii] = 0;
            }
            kk = 0; index = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                  kk += inputPDFs[ii];
                  index++;
               }
            }
            lowers = new double[M4];
            uppers = new double[M4];
            index = 0;
            for (ii = 0; ii < nInputs; ii++)
            {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                  lowers[index] = odata->lowerBounds_[ii];
                  uppers[index] = odata->upperBounds_[ii];
                  index++;
               }
            }
            if (kk == 0)
            {
               sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
               sampler->setPrintLevel(0);
               sampler->setInputBounds(M4, lowers, uppers);
               sampler->setOutputParams(iOne);
               psOUULargeSampleSize_ = 100000;
               sampler->setSamplingParams(psOUULargeSampleSize_, iOne, iZero);
               sampler->initialize(0);
               psOUULargeSampleSize_ = sampler->getNumSamples();
               psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
               psOUULargeSamOutputs_ = new double[psOUULargeSampleSize_];
               samStates = new int[psOUULargeSampleSize_];
               sampler->getSamples(psOUULargeSampleSize_, M4, iOne,
                            psOUULargeSamInputs_,psOUULargeSamOutputs_,
                            samStates);
               delete [] samStates;
               delete sampler;
            }
            else
            {
               odata->psIO_->getParameter("input_means", pdata);
               inputMeans = pdata.dbleArray_;
               if (inputMeans == NULL)
               {
                  inputMeans = new double[nInputs];
                  for (ii = 0; ii < nInputs; ii++) inputMeans[ii] = 0;
               }
               pdata.dbleArray_ = NULL;
               odata->psIO_->getParameter("input_stdevs", pdata);
               inputStdevs = pdata.dbleArray_;
               if (inputStdevs == NULL)
               {
                  inputStdevs = new double[nInputs];
                  for (ii = 0; ii < nInputs; ii++) inputStdevs[ii] = 1;
               }
               pdata.dbleArray_ = NULL;
               odata->psIO_->getParameter("input_cor_matrix", pdata);
               corMat1 = (Matrix *) pdata.psObject_;
               pdata.psObject_ = NULL;

               corMat2.setDim(M4,M4);
               iPdfs  = new int[nInputs];
               iMeans = new double[nInputs];
               iStdvs = new double[nInputs];
               index = 0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  if (psOUUInputTypes_[ii] == psOUUType4)
                  {
                     iMeans[index] = inputMeans[ii];
                     iStdvs[index] = inputStdevs[ii];
                     iPdfs[index]  = inputPDFs[ii];
                     index2 = 0;
                     for (kk = 0; kk < nInputs; kk++)
                     {
                        if (psOUUInputTypes_[kk] == psOUUType4)
                        {
                           ddata = corMat1->getEntry(ii,kk);
                           corMat2.setEntry(index,index2,ddata);
                           index2++;
                        }
                     }
                     index++;
                  }
               }
               pdfman = new PDFManager();
               pdfman->initialize(M4,iPdfs,iMeans,iStdvs,corMat2,NULL,NULL);
               vecLB.load(M4, lowers);
               vecUB.load(M4, uppers);
               psOUULargeSampleSize_ = 100000;
               vecOut.setLength(psOUULargeSampleSize_*M4);
               pdfman->genSample(psOUULargeSampleSize_, vecOut, vecLB, vecUB);
               psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
               psOUULargeSamOutputs_ = new double[psOUULargeSampleSize_];
               for (ii = 0; ii < psOUULargeSampleSize_*M4; ii++)
                  psOUULargeSamInputs_[ii] = vecOut[ii];
               delete [] inputMeans;
               delete [] inputStdevs;
               delete [] iPdfs;
               delete [] iMeans;
               delete [] iStdvs;
               delete pdfman;
            }
            delete [] inputPDFs;
            delete [] lowers;
            delete [] uppers;
         }
      }
   }

   int nSamp=psOUUZ3nSamples_*psOUUZ4nSamples_;
   psOUUSamOutputs_ = new double[nSamp];
   if (M2 > 0)
   {
      printf("For 2-stage OUU, you have 2 options for inner optimization:\n");
      printf("(1) You can use BOBYQA available in PSUADE (in this case\n");
      printf("    opt_driver should point to your original function), or\n");
      printf("(2) You can provide your own optimizer (in opt_driver).\n");
      printf("Use your own optimizer instead of BOBYQA? (y or n) ");
      scanf("%s", lineIn);
      if (lineIn[0] == 'y') 
      {
         psOUUUserOpt_ = 1;
         printf("NOTE: Make sure your optimizer executable has been\n");
         printf("      assigned to 'opt_driver' and it optimizes with\n");
         printf("      respect to the %d-th to %d-th parameters.\n", 
                M1+1, M1+M2);
      }
   }
   else
   {
      printf("Since no recourse variable (for level 2 optimization) has\n");
      printf("been specified, PSUADE presumes that you are either doing\n");
      printf("1-level OUU, or you are providing the inner optimization\n");
      printf("solver in opt_driver.\n");
      psOUUUserOpt_ = 1;
   }
   printEquals(PL_INFO, 0);
 
   if (psOUUUserOpt_ == 1 && (M3+M4) > 0)
   {
      printf("Each simulation will call your opt_driver with one sample point.\n");
      printf("However, for higher efficiency (less I/O), you have the option\n");
      printf("to provide in 'ensemble_opt_driver' an executable that can run\n");
      printf("multiple sample points. In this case, OUU will call your ensemble\n");
      printf("executable with the following sequence: \n");
      printf("      <Your ensemble_opt_driver> <sampleFile> <outputFile>\n\n");
      printf("where <sampleFile> is in the following format:\n");
      printf("   line 1: <nSamples>\n");
      printf("   line 2: Sample point 1 input values\n");
      printf("   line 3: Sample point 2 input values\n");
      printf("   line n: ...\n\n");
      printf("and <outputFile> should have all sample output values.\n\n");
      printf("Use ensemble opt driver for ensemble runs ? (y or n) ");
      lineIn[0] = '1';
      while (lineIn[0] != 'n' && lineIn[0] != 'y')
      {
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUEnsembleEval_ = 1;
         printEquals(PL_INFO, 0);
      }
   }

   if ((psOUUEnsembleEval_ == 0) && ((M3+M4)>0))
   {
      printf("You can configure OUU to run the ensemble simulations in\n");
      printf("parallel/asynchronous using the Linux fork/join.\n");
      printf("If you say 'n', your simulator (opt_driver) will be evaluated\n");
      printf("sequentially (one sample point at a time). But if you say 'y',\n");
      printf("be careful; because PSUADE will launch %d jobs simultaneously.\n",
             psOUUZ3nSamples_*psOUUZ4nSamples_);
      printf("Turn on asynchronous mode ? (y or n) ");
      lineIn[0] = '1';
      while (lineIn[0] != 'n' && lineIn[0] != 'y')
      {
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUParallel_ = 1;
         printEquals(PL_INFO, 0);
      }
   }
   fgets(lineIn, 500, stdin);

   psOUUWValues_ = new double[nInputs];
   psOUUXValues_ = new double[nInputs*nSamp];
   psOUUOptimalX_ = new double[nInputs];

   int rstype=0, rsaux=0;
   psOUUfaPtr_ = NULL;
   if (psOUUUseRS_ == 1 && M4 > 0)
   {
      if (psOUUMasterMode_ == 0)
      {
         if (printLevel > 2) 
            printf("OUUOptimizer: setting up response surface\n");
         if (psOUUZ4nSamples_ > 600)
         {
            rstype = PSUADE_RS_MARS;
            printf("OUUOptimizer: use MARS since nSamples > 400\n");
         }
         else if (psOUUZ4nSamples_ > 300)
         {
            rstype = PSUADE_RS_RBF;
            printf("OUUOptimizer: use RBF response surface\n");
         }
         else if (psOUUZ4nSamples_ > 200)
         {
            rstype = PSUADE_RS_KR;
            printf("OUUOptimizer: use Kriging (fast) response surface\n");
         }
         else
         {
            rstype = PSUADE_RS_KR;
            rsaux = 1;
            printf("OUUOptimizer: use Kriging (slow) response surface\n");
         }
      }
      else
      {
         printf("Response surface available: \n");
         printf("1. MARS\n");
         printf("2. Kriging (fast)\n");
         printf("3. Kriging (slow)\n");
         printf("4. Radial basis function\n");
         sprintf(pString, "Which response surface? ");
         rstype = getInt(1,4,pString); 
         if (rstype == 1) rstype = PSUADE_RS_MARS;
         if (rstype == 2)
         {
            rstype = PSUADE_RS_KR;
            rsaux = 0;
         }
         if (rstype == 3)
         {
             rstype = PSUADE_RS_KR;
             rsaux = 1;
         }
         if (rstype == 4) rstype = PSUADE_RS_RBF;
      }
      kk = psInteractive_;
      psInteractive_ = 0;
      psOUUfaPtr_ = genFA(rstype, M4, -1, psOUUZ4nSamples_);
      lowers = new double[M4];
      uppers = new double[M4];
      psOUUZ4RSType_ = rstype;
      psOUUZ4RSAux_ = rsaux;
      psOUUZ4LBounds_ = new double[M4];
      psOUUZ4UBounds_ = new double[M4];
      index = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (psOUUInputTypes_[ii] == psOUUType4)
         {
            lowers[index] = odata->lowerBounds_[ii];
            uppers[index] = odata->upperBounds_[ii];
            psOUUZ4LBounds_[index] = lowers[index];
            psOUUZ4UBounds_[index] = uppers[index];
            index++;
         }
      }
      psOUUfaPtr_->setBounds(lowers,uppers);
      psOUUfaPtr_->setOutputLevel(0);
      psInteractive_ = kk;
      if (rstype == PSUADE_RS_KR) 
      {
         targv = new char*[1];
         targv[0] = new char[100];
         if (rsaux == 0) strcpy(targv[0], "setMode2");
         else            strcpy(targv[0], "setMode3");
         psOUUfaPtr_->setParams(1, targv);
         delete [] targv[0];
         delete [] targv;
      }
      delete [] lowers;
      delete [] uppers;
   }
   if (M3+M4 == 0)
   {
      psOUUEnsembleEval_ = 0;
      psOUUMode_ = 1;
   }

   nPts = (M1 + 1) * (M1 + 2) / 2;
   workArray = new double[(nPts+5)*(nPts+M1)+3*M1*(M1+5)/2+1];
   psOUUCounter_ = 0;
   if (psOUUParallel_ == 1) odata->funcIO_->setAsynchronousMode();

#ifdef HAVE_BOBYQA
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psOUUSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL) 
      {
         printf("OUU: use history has been turned on.\n");
         fp = fopen("psuade_ouu_history", "r");
         if (fp != NULL)
         {
            printf("OUUOptimizer history file found.\n");
            psOUUNSaved_ = 0;
            while (feof(fp) == 0)
            {
               fscanf(fp, "%d %d", &ii, &kk);
               if (ii != 999 || kk != nInputs)
               {
                  printf("OUU history file not used - nInputs mismatch.\n");
                  break;
               }
               else
               {
                  for (ii = 0; ii < nInputs; ii++)
                     fscanf(fp, "%lg",&psOUUSaveX_[psOUUNSaved_*nInputs+ii]);
                  fscanf(fp, "%lg",&psOUUSaveY_[psOUUNSaved_]);
                  psOUUNSaved_++;
               }
               if (((psOUUNSaved_+1)*nInputs > psOUUMaxSaved_*10) ||
                   psOUUNSaved_ > psOUUMaxSaved_) break;
            }
            fclose(fp);
         }
      }
   }

   for (ii = 0; ii < M1; ii++) 
      printf("OUUOptimizer initial X %3d = %e\n", ii+1, XValues[ii]);
   bobyqa_(&M1, &nPts, XValues, odata->lowerBounds_, odata->upperBounds_, 
           &rhobeg, &rhoend, &bobyqaFlag, &maxfun, workArray);
   printf("OUUOptimizer: total number of evaluations = %d\n",
           odata->numFuncEvals_);

   if (psOUUSaveHistory_ == 1 && psOUUNSaved_ > 0)
   {
      fp = fopen("psuade_ouu_history","w");
      if (fp != NULL)
      {
         for (ii = 0; ii < psOUUNSaved_; ii++)
         {
            fprintf(fp, "999 %d ", nInputs);
            for (kk = 0; kk < nInputs; kk++)
               fprintf(fp, "%24.16e ", psOUUSaveX_[ii*nInputs+kk]);
            fprintf(fp, "%24.16e\n", psOUUSaveY_[ii]);
         }
         fclose(fp);
      }
      printf("OUUOptimizer: history saved in psuade_ouu_history\n");
   }
#else
   printf("ERROR : Bobyqa optimizer not installed.\n");
   exit(1);
#endif
   if (psOUUParallel_ == 1) odata->funcIO_->setSynchronousMode();

   if (odata->setOptDriver_ & 2)
   {
      odata->funcIO_->setDriver(currDriver);
   }
   delete [] XValues;
   delete [] workArray;
   if (psOUUZ3SamInputs_ != NULL) delete [] psOUUZ3SamInputs_;
   if (psOUUZ4SamInputs_ != NULL) delete [] psOUUZ4SamInputs_;
   if (psOUUZ4LBounds_ != NULL) delete [] psOUUZ4LBounds_;
   if (psOUUZ4UBounds_ != NULL) delete [] psOUUZ4UBounds_;
   if (psOUUSamOutputs_  != NULL) delete [] psOUUSamOutputs_;
   if (psOUUSamProbs_    != NULL) delete [] psOUUSamProbs_;
   if (psOUULargeSamInputs_  != NULL) delete [] psOUULargeSamInputs_;
   if (psOUULargeSamOutputs_ != NULL) delete [] psOUULargeSamOutputs_;
   if (psOUUXValues_  != NULL) delete [] psOUUXValues_;
   if (psOUUWValues_  != NULL) delete [] psOUUWValues_;
   if (psOUUOptimalX_ != NULL) delete [] psOUUOptimalX_;
   if (psOUUInputTypes_ != NULL) delete psOUUInputTypes_;
   if (psOUUfaPtr_ != NULL) delete psOUUfaPtr_;
   psOUUZ3SamInputs_ = NULL;
   psOUUZ4SamInputs_ = NULL;
   psOUUZ4LBounds_ = NULL;
   psOUUZ4UBounds_ = NULL;
   psOUUSamOutputs_ = NULL;
   psOUUSamProbs_ = NULL;
   psOUULargeSamInputs_ = NULL;
   psOUULargeSamOutputs_ = NULL;
   psOUUXValues_ = NULL;
   psOUUWValues_ = NULL;
   psOUUOptimalX_ = NULL;
   psOUUInputTypes_ = NULL;
   psOUUfaPtr_ = NULL;
   odata->intData_ = M1;
   if (psOUUMasterMode_ != 0) psMasterMode_ = 1;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
OUUOptimizer& OUUOptimizer::operator=(const OUUOptimizer &)
{
   printf("OUUOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************ 
// validate
// ------------------------------------------------------------------------
extern "C" 
{
void validate(int nSamps,double *samInputs,double *samOutputs,double *errors)
{
   int        nInps, nSubSamples, ii, ss, *iArray;
   double     *X, *Y;
   char       **targv;
   FuncApprox *fa;

   for (ii = 0; ii < 3; ii++) errors[ii] = 1e12;
   nInps = psOUUM4_;
   if (nSamps <= 20) nSubSamples = 1;
   else
   {
      nSubSamples = nSamps / 10;
      if (nSubSamples*10 != nSamps) nSubSamples++;
   }
   int rsState = psRSExpertMode_;
   int saveInteractive = psInteractive_;
   psRSExpertMode_ = 0;
   psInteractive_ = 0;
   fa = genFA(psOUUZ4RSType_, nInps, -1, nSamps-nSubSamples);
   if (fa == NULL)
   {
      printOutTS(PL_ERROR,"ERROR: cannot create function approximator.\n");
      return;
   }
   fa->setBounds(psOUUZ4LBounds_,psOUUZ4UBounds_);
   fa->setOutputLevel(0);
   if (psOUUZ4RSType_ == PSUADE_RS_KR) 
   {
      targv = new char*[1];
      targv[0] = new char[100];
      if (psOUUZ4RSAux_ == 0) strcpy(targv[0], "setMode2");
      else                    strcpy(targv[0], "setMode3");
      fa->setParams(1, targv);
      delete [] targv[0];
      delete [] targv;
   }
   iArray = new int[nSamps];
   if (nSubSamples > 1) generateRandomIvector(nSamps, iArray);
   else
   {
      for (ss = 0; ss < nSamps; ss++) iArray[ss] = ss;
   }
   X = new double[nSamps*nInps];
   for (ii = 0; ii < nInps; ii++)
   {
      for (ss = 0; ss < nSamps; ss++)
         X[iArray[ss]*nInps+ii] = samInputs[ss*nInps+ii];
   }
   Y = new double[nSamps];
   for (ss = 0; ss < nSamps; ss++) Y[iArray[ss]] = samOutputs[ss];
   oData *odata = (oData *) psOUUObj_;
   if (odata != NULL && odata->outputLevel_ > 2)
      printOutTS(PL_INFO, "Cross validation begins ...\n");
   int    status, ss2, count;
   double cvErr1s,cvErr2s,cvMaxs,CVErr1s=0,CVErr2s=0,CVMaxs=0,ddata;
   double *YT = new double[nSamps];
   double *X2 = new double[nSamps*nInps];
   double *Y2 = new double[nSamps];
   for (ss = 0; ss < nSamps; ss+=nSubSamples)
   {
      if (odata != NULL && odata->outputLevel_ > 2)
         printOutTS(PL_INFO, "Cross validation group %d (of 10)\n",
                    ss/nSubSamples+1);
      count = 0;
      for (ss2 = 0; ss2 < nSamps; ss2++)
      {
         if (ss2 < ss || ss2 >= (ss+nSubSamples))
         {
            for (ii = 0; ii < nInps; ii++)
               X2[count*nInps+ii] = X[ss2*nInps+ii];
            Y2[count++] = Y[ss2];
         }
      }
      status = fa->initialize(X2, Y2);
      if (status == -1) break;
      count = nSubSamples;
      if ((ss + nSubSamples) > nSamps) count = nSamps - ss;
      fa->evaluatePoint(count, &(X[ss*nInps]), YT);
      cvErr1s = cvErr2s = cvMaxs = 0.0;
      for (ss2 = 0; ss2 < count; ss2++)
      {
         ddata = YT[ss2] - Y[ss+ss2];
         if (Y[ss+ss2] != 0.0) ddata = ddata / PABS(Y[ss+ss2]);
         cvErr1s += ddata;
         cvErr2s += (ddata * ddata);
         if (PABS(ddata) > cvMaxs) cvMaxs = PABS(ddata);
      }
      CVErr1s += cvErr1s;
      CVErr2s += cvErr2s;
      if (cvMaxs > CVMaxs ) CVMaxs = cvMaxs; 
   }
   if (status >= 0)
   {
      errors[0] = CVErr1s / (double) nSamps;
      errors[1] = sqrt(CVErr2s / nSamps);
      errors[2] = CVMaxs;
   }
   psRSExpertMode_ = rsState;
   psInteractive_ = saveInteractive; 
   delete fa;
   delete [] X;
   delete [] Y;
   delete [] YT;
   delete [] X2;
   delete [] Y2;
   delete [] iArray;
}
}


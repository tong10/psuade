// ************************************************************************
// Copyright (c) 2019   Lawrence Livermore National Security, LLC.
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
// Functions for the class IntOptimizer
// AUTHOR : Amer Abdulla 
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeUtil.h"
#include "IntOptimizer.h"
#include "PsuadeUtil.h"
#include "Sampling.h"
#include "sysdef.h"
#include "Psuade.h"
#include "PrintingTS.h"

// ------------------------------------------------------------------------
#include <math.h> // for standev and georange functions
#include <time.h> // for random number generator
// ------------------------------------------------------------------------

#define psINTMaxSaved_ 100000
int     psINTNSaved_=0;
double  psINTSaveX_[psINTMaxSaved_*10];
double  psINTSaveY_[psINTMaxSaved_];
void    *psINTObj_=NULL;
int     *psINTInputTypes_=NULL;
int     psINTnInputs_=0;
#define psINTnHist_ 10
double  psINTHistory_[psINTnHist_];
int     psINTHistCnt_=0;
int     psINTHistIndex_=0;
int     psINTCurrDriver_=-1;

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
  void *INTevalfunc_(int *nInps, double *XValues, double *YValue)
  {
    int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
    double *localY;
    oData  *odata;
    FILE   *infile;

    //**/ ------ get history ------
    nInputs = (*nInps);
    if (psConfig_.OptExpertModeIsOn() && psINTNSaved_ == 0)
    {
      infile = fopen("psuade_int_data","r");
      if (infile != NULL)
      {
        fscanf(infile, "%d %d", &psINTNSaved_, &kk);
        if ((psINTNSaved_ <= 0) ||
            (psINTNSaved_+1 > psINTMaxSaved_*10/nInputs))
        {
          printf("PSUADE INT: history file has too much data.\n");
          printf("            Ask PSUADE developer for help.\n"); 
          fclose(infile);
          psINTNSaved_ = 0;
        }
        else if (kk != nInputs)
        {
          printf("PSUADE INT: history file has invalid input count.\n");
          fclose(infile);
          psINTNSaved_ = 0;
        }
        else
        {
          for (ii = 0; ii < psINTNSaved_; ii++)
          {
            fscanf(infile, "%d", &kk);
            if (kk != ii+1)
            {
              printf("PSUADE INT: data index mismatch.\n");
              psINTNSaved_ = 0;
              break;
            }
            for (jj = 0; jj < nInputs; jj++)
              fscanf(infile, "%lg", &psINTSaveX_[ii*nInputs+jj]);
            fscanf(infile, "%lg", &psINTSaveY_[ii]);
          }
          fclose(infile);
        }
      }
    }

    //**/ ------ fetch data ------
    odata    = (oData *) psINTObj_;
    nOutputs = odata->nOutputs_;
    localY   = (double *) malloc(nOutputs * sizeof(double));
    outputID = odata->outputID_;

    //**/ ------ search to see if it has already been evaluated ------
    found = 0;
    for (ii = 0; ii < psINTNSaved_; ii++)
    {
      for (jj = 0; jj < nInputs; jj++)
        if (PABS(psINTSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
      if (jj == nInputs)
      {
        found = 1;
        break;
      }
    }

    //**/ ------ run simulation ------
    funcID = odata->numFuncEvals_;
    if (found == 0)
    {
      odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
      odata->numFuncEvals_++;
    }
    else localY[outputID] = psINTSaveY_[ii];
    (*YValue) = localY[outputID];

    //**/ ------ save the data ------
    if (psConfig_.OptExpertModeIsOn() && found == 0)
    {
      for (jj = 0; jj < nInputs; jj++)
        psINTSaveX_[psINTNSaved_*nInputs+jj] = XValues[jj];
      psINTSaveY_[psINTNSaved_] = localY[outputID];
      psINTNSaved_++;
    }
    // Use for Diagnostics
    //for (ii = 0; ii < nInputs; ii++)
    //printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
    //printf("    Y = %16.8e\n", localY[outputID]);

    //**/ ------ store optimal information ------
    if ((*YValue) < odata->optimalY_)
    {
      odata->optimalY_ = (*YValue);
      for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
      psINTHistIndex_ = funcID;
      if (odata->outputLevel_ > 0)
      {
        printf("INTOptimizer %6d : \n", odata->numFuncEvals_);
        if (odata->outputLevel_ > 1)
          for (ii = 0; ii < nInputs; ii++)
            printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
        printf("    Ymin  = %16.8e\n", odata->optimalY_);
      }
      if (psINTHistCnt_ < psINTnHist_) 
        psINTHistory_[psINTHistCnt_++] = (*YValue); 
      else
      {
        for (ii = 1; ii < psINTnHist_; ii++) 
          psINTHistory_[ii-1] = psINTHistory_[ii];
        psINTHistory_[psINTnHist_-1] = (*YValue);
        psINTHistCnt_++;
      }
    }
    free(localY);

    //**/ ------ save history ------
    if (psConfig_.OptExpertModeIsOn() && found == 0)
    {
      infile = fopen("psuade_int_data","w");
      if (infile != NULL)
      {
        fprintf(infile, "%d %d\n", psINTNSaved_, nInputs);
        for (ii = 0; ii < psINTNSaved_; ii++)
        {
          fprintf(infile, "%d ", ii+1);
          for (jj = 0; jj < nInputs; jj++)
            fprintf(infile, "%24.16e ", psINTSaveX_[ii*nInputs+jj]);
          fprintf(infile, "%24.16e\n", psINTSaveY_[ii]);
        }
        fclose(infile);
      }
    }
    return NULL;
  }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
INTOptimizer::INTOptimizer()
{
  psINTNSaved_ = 0;
  psINTObj_ = NULL;
  psINTInputTypes_ = NULL;
  psINTnInputs_ = 0;
  psINTHistCnt_ = 0;
  psINTHistIndex_ = 0;
  psINTCurrDriver_ = -1;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
INTOptimizer::~INTOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void INTOptimizer::optimize(oData *odata) 
{
  int ii, jj;

  //**/ ------ prepare object variables ------
  int printLevel = odata->outputLevel_;
  int nInputs  = odata->nInputs_;
  int nOutputs = odata->nOutputs_;
  if (nOutputs > 1)
  {
    printOutTS(PL_ERROR,"INT ERROR: nOutputs = %d.\n",nOutputs);
    printOutTS(PL_ERROR,"       Only nOutputs=1 is allowed.\n");
    printOutTS(PL_ERROR,"       INT cannot handle constraints.\n");
    printOutTS(PL_ERROR,"       Suggestion: use penalty term.\n");
    exit(1);
  }
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 1.0;
  int maxfun = odata->maxFEval_;
  odata->optimalY_ = 1.0e50;
  odata->numFuncEvals_ = 0;

  //**/ set optimization flag
  if ((odata->setOptDriver_ & 1))
  {
    psINTCurrDriver_ = odata->funcIO_->getDriver();
    odata->funcIO_->setDriver(1);
  }
  psINTObj_= (void *) odata;

  // Initialize INT parameters
  int nPts = (int) pow(3.0, nInputs) + 1;
  psVector vecSwarm;
  psINTHistCnt_ = 0;
  psINTHistIndex_ = -1;
  vecSwarm.setLength(nPts*nInputs);

  // For each parameter, determine range of points from which 
  // we can sample
  psIVector vecLBs, vecUBs, vecRanges;
  vecLBs.setLength(nInputs);
  vecUBs.setLength(nInputs);
  vecRanges.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++)
  {
    vecLBs[ii] = (int) odata->lowerBounds_[ii];
    vecUBs[ii] = (int) odata->upperBounds_[ii];
    vecRanges[ii] = vecUBs[ii] - vecLBs[ii];
  }

  //**/ determine whether to include initial point in 
  //**/ starting population  

  //**/ loop
  int icnt, iter = 0, nEvals;
  int stagnationCount=0, minInd;
  double minVal, ddata;
  psIVector vecCounters, vecMaxs, vecBest;
  vecCounters.setLength(nInputs);
  vecMaxs.setLength(nInputs);
  vecBest.setLength(nInputs);
  double *XSwarm = vecSwarm.getDVector(); 
  psVector vecYT;
  vecYT.setLength(nPts);
  double *YT = vecYT.getDVector(); 
  for (ii = 0; ii < nInputs; ++ii) 
    vecSwarm[ii] = (int) odata->initialX_[ii];

  minVal = PSUADE_UNDEFINED;
  int restart=1, repeatCnt=0;
  double localMinVal = PSUADE_UNDEFINED;
  int mutateInp=0, upDown=0, nRestart=0;
  while (iter < maxfun)
  {
    iter++;
    //**/ create simplex
    for (ii = 0; ii < nInputs; ii++)
    {
      if (vecSwarm[ii]-1 >= vecLBs[ii]) 
           vecCounters[ii] = (int) vecSwarm[ii] - 1;
      else vecCounters[ii] = (int) vecSwarm[ii];
      vecMaxs[ii] = (int) vecSwarm[ii] + 1;
      if (vecMaxs[ii] > vecUBs[ii]) vecMaxs[ii]--;
    }
    nEvals = 1;
    while (vecCounters[nInputs-1] < vecMaxs[nInputs-1])
    {
      for (ii = 0; ii < nInputs; ii++)
        if (vecCounters[ii] != vecSwarm[ii]) break;
      if (ii != nInputs)
      {
        for (ii = 0; ii < nInputs; ii++)
          vecSwarm[nEvals*nInputs+ii] = 1.0 * vecCounters[ii];
        nEvals++;
      }
      icnt = 0;
      vecCounters[icnt]++;
      while (icnt < nInputs-1)
      {
        if (vecCounters[icnt] <= vecMaxs[icnt])
          break;
        else 
        {
          vecCounters[icnt] = vecMaxs[icnt] - 2;
          if (vecCounters[icnt] < vecLBs[icnt])
            vecCounters[icnt] = vecMaxs[icnt] - 1;
          vecCounters[icnt+1]++;
        }
        icnt++;
      } 
    } 
    //**/ evaluate
    for (ii = 0; ii < nEvals; ii++)
    {
      INTevalfunc_(&nInputs, &(XSwarm[ii*nInputs]), &(YT[ii]));
    }

    //**/ see if this a a better one 
    minInd = -1;
    for (ii = 0; ii < nEvals; ++ii)
    {
      if (YT[ii] < localMinVal)
      {
        localMinVal = YT[ii];
        minInd = ii;
      }
    }
    if (minInd == -1)
    {
      for (ii = 0; ii < nEvals; ++ii)
      {
        if (YT[ii] <= localMinVal)
        {
          localMinVal = YT[ii];
          minInd = ii;
        }
      }
      repeatCnt++;
    }

    //**/ if this is a better one than before
    if (minInd > 0 && repeatCnt < 3)
    {
      if (localMinVal < minVal)
      {
        minVal = localMinVal;
        for (ii = 0; ii < nInputs; ii++)
          vecBest[ii] = XSwarm[minInd*nInputs+ii];
        printf("New optimal value at iteration %8d = %e\n",iter,minVal);
        for (jj = 0; jj < nInputs; jj++)
          printf("  Input %4d = %d\n", jj+1,(int) vecBest[jj]);
      }
      for (ii = 0; ii < nInputs; ii++)
        XSwarm[ii] = XSwarm[minInd*nInputs+ii];
      restart = 0;
      if (printLevel > 0)
        for (jj = 0; jj < nInputs; jj++)
          printf("Move Base Input %4d = %5d (%4d = %12.4e,%6d)\n", jj+1,
                (int) vecSwarm[jj],minInd,localMinVal,repeatCnt);
      mutateInp = 0;
      upDown = 0;
    }
    else
    {
      if (mutateInp < nInputs)
      {
        if (upDown == 0 && vecSwarm[mutateInp] > (vecLBs[mutateInp]+1))
        {
          vecSwarm[mutateInp] = vecSwarm[mutateInp] - 2;
          upDown = 1;
          if (printLevel > 0)
            printf("(1) Mutate Input %4d = %d (%d)\n", mutateInp+1,
                   (int) vecSwarm[mutateInp],repeatCnt);
        }
        else if (upDown == 1 && vecSwarm[mutateInp] < (vecUBs[mutateInp]-1))
        {
          vecSwarm[mutateInp] = vecSwarm[mutateInp] + 2;
          upDown = 0;
          if (printLevel > 0)
            printf("(2) Mutate Input %4d = %d (%d)\n", mutateInp+1,
                   (int) vecSwarm[mutateInp],repeatCnt);
          mutateInp++;
        }
        else
        {
          vecSwarm[mutateInp] = (double) (PSUADE_rand() % 
                        vecRanges[mutateInp] + vecLBs[mutateInp]);
          upDown = 0;
          if (printLevel > 0)
            printf("(3) Mutate Input %4d = %d (%d)\n", mutateInp+1,
                   (int) vecSwarm[mutateInp],repeatCnt);
          mutateInp++;
        }
      }
      else
      {
        mutateInp = 0;
        upDown = 0;
        for (ii = 0; ii < nInputs; ii++)
          vecSwarm[ii] = (double) (PSUADE_rand() % vecRanges[ii] + 
                                   vecLBs[ii]);
        restart = 1;
        localMinVal = PSUADE_UNDEFINED;
        repeatCnt = 0;
        nRestart++;
        printf("Restart %d\n", nRestart);
        if (printLevel > 0)
          for (jj = 0; jj < nInputs; jj++)
            printf("Restart Input %4d = %d\n", jj+1,(int) vecSwarm[jj]);
      }
    }
    if (printLevel > 0)
    {
      printf("INT current iter %d min value = %e\n",iter,minVal);
      for (jj = 0; jj < nInputs; jj++)
        printf("  Input %4d = %d\n", jj+1,(int) vecBest[jj]);
    }
  }

  psINTnInputs_ = nInputs;
  if (printLevel > 0)
  {
    printf("INTOptimizer: number of function evaluations = %d\n",
           odata->numFuncEvals_);
  }

  //**/ ------ reset things ------
  if ((odata->setOptDriver_ & 2) && psINTCurrDriver_ >= 0)
  {
    odata->funcIO_->setDriver(psINTCurrDriver_);
  }
}


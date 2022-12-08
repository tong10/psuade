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
// ************************************************************************
// Internal global variables 
// ------------------------------------------------------------------------
#define psCobylaMaxSaved_ 10000
int     psCobylaNSaved_=0;
double  psCobylaSaveX_[psCobylaMaxSaved_*10];
double  psCobylaSaveY_[psCobylaMaxSaved_*10];
psVector VecPsCobylaXLbounds_;
psVector VecPsCobylaXUbounds_;
double  *psCobylaXLbounds_=NULL;
double  *psCobylaXUbounds_=NULL;
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
    int    ii, jj, kk, funcID, nOuts, nOuts2, ignoreFlag, found, index;
    double *localY; 
    oData  *odata;
    FILE   *fp;

    //**/ ---- fetch object and allocate memory to store outputs ------
    odata  = (oData *) data;
    nOuts  = nConstraints + 1;
    nOuts2 = nOuts - 2 * nInputs;
    localY = (double *) malloc(nOuts * sizeof(double));

    //**/ ------ search to see if it has already been evaluated ------
    found = 0;
    for (ii = 0; ii < psCobylaNSaved_; ii++)
    {
      for (jj = 0; jj < nInputs; jj++)
        if (PABS(psCobylaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) 
          break;
      if (jj == nInputs)
      {
        found = 1;
        index = ii;
        printf("Cobyla: simulation results reuse (%d).\n",ii+1);
        break;
      }
    }

    //**/ ------ run simulation (if not already been evaluated) ------
    funcID = odata->numFuncEvals_;
    if (found == 0)
    {
      if (odata->optFunction_ != NULL)
        //**/ if users have loaded a simulation function, use it
        odata->optFunction_(nInputs, XValues, nOuts2, localY);
      else if (odata->funcIO_ != NULL)
        //**/ if not, use the simulation function set up in odata
        odata->funcIO_->evaluate(funcID,nInputs,XValues,nOuts2,localY,0);
      else
      {
        printf("CobylaOptimizer ERROR: no function evaluator.\n");
        exit(1);
      }
      funcID = odata->numFuncEvals_++;
      (*YValue) = localY[0];
      //**/ negation because Cobya considers g(X) >= 0 as feasible
      for (ii = 0; ii < nOuts2-1; ii++) constraints[ii] = - localY[ii+1];
      for (ii = 0; ii < nInputs; ii++) 
      {
        constraints[nOuts2-1+ii*2] = XValues[ii]-VecPsCobylaXLbounds_[ii];;
        constraints[nOuts2-1+ii*2+1] = VecPsCobylaXUbounds_[ii]-XValues[ii];
      }
    }
    else
    {
      localY[0] = psCobylaSaveY_[index*nOuts];
      (*YValue) = localY[0];
      for (ii = 0; ii < nOuts2-1; ii++) 
        constraints[ii] = psCobylaSaveY_[index*nOuts2+ii+1];
      for (ii = 0; ii < nInputs; ii++) 
      {
        constraints[nOuts2-1+ii*2] = XValues[ii] - VecPsCobylaXLbounds_[ii];;
        constraints[nOuts2-1+ii*2+1] = VecPsCobylaXUbounds_[ii] - XValues[ii];
      }
    }

    //**/ ------ save the simulation data ------
    if (psCobylaSaveHistory_ == 1 && found == 0 &&
      (psCobylaNSaved_+1)*nInputs < psCobylaMaxSaved_*10 &&
       (psCobylaNSaved_+1)*nOuts2 < psCobylaMaxSaved_*10 &&
       psCobylaNSaved_ < psCobylaMaxSaved_)
    {
      for (jj = 0; jj < nInputs; jj++)
        psCobylaSaveX_[psCobylaNSaved_*nInputs+jj] = XValues[jj];
      psCobylaSaveY_[psCobylaNSaved_*nOuts2] = localY[0];
      for (jj = 1; jj < nOuts2; jj++)
        psCobylaSaveY_[psCobylaNSaved_*nOuts2+jj] = - localY[jj];
      psCobylaNSaved_++;
      fp = fopen("psuade_cobyla_history","w");
      if (fp != NULL)
      {
        for (ii = 0; ii < psCobylaNSaved_; ii++)
        {
          fprintf(fp, "999 %d %d ", nInputs, nOuts2);
          for (jj = 0; jj < nInputs; jj++)
            fprintf(fp, "%24.16e ", psCobylaSaveX_[ii*nInputs+jj]);
          for (jj = 0; jj < nOuts2; jj++)
            fprintf(fp, "%24.16e ", psCobylaSaveY_[ii*nOuts2+jj]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }
    }
    if (psConfig_.InteractiveIsOn() && odata->outputLevel_ > 1)
    {
      printf("CobylaOptimizer %6d : \n", odata->numFuncEvals_+1); 
      for (ii = 0; ii < nInputs; ii++)
        printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
      printf("    ObjY  = %16.8e\n", localY[0]);
      if (odata->outputLevel_ > 4)
        for (ii = 0; ii < nConstraints; ii++)
          printf("    Constraint %3d = %16.8e\n",ii+1,constraints[ii]);
    }

    //**/ ------ update optimal (feasible) solution ------
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
  if (psConfig_.InteractiveIsOn())
  {
    printAsterisks(PL_INFO, 0);
    printAsterisks(PL_INFO, 0);
    printf("*   COBYLA Inequality-Constrained Optimization\n");
    printEquals(PL_INFO, 0);
    printf("* - To run this optimizer in batch mode, first ");
    printf("make sure opt_driver in\n");
    printf("*   your PSUADE input file has been set to point ");
    printf("to your objective\n");
    printf("*   function evaluator.\n");
    printf("* - Set optimization tolerance in your PSUADE ");
    printf("input file\n");
    printf("* - Set maximum number of iterations in PSUADE ");
    printf("input file\n");
    printf("* - Set num_starts to perform multi-start ");
    printf("optimization\n");
    printf("* - Set optimization print_level to give more ");
    printf("screen outputs\n");
    printf("* - In opt_expert mode, the optimization history ");
    printf("log will be turned on\n");
    printf("*   automatically. Previous psuade_cobyla_history ");
    printf("file will also be\n");
    printf("*   reused.\n");
    printf("* - If your opt_driver is a response surface ");
    printf("which has more inputs\n");
    printf("*   than the number of optimization inputs, you ");
    printf("can fix some driver\n");
    printf("*   driver inputs by creating an rs_index_file.\n");
    printf("* - If nOutput=1 and if only input bound ");
    printf("constraints are imposed, then\n");
    printf("*   nConstraints is expected to be 2m where ");
    printf("m is the number of inputs.\n");
    printf("*   NOTE: FEASIBLE CONSTRAINTS have non-positive values.\n");
    printf("*         (i.e. simulation outputs for any ");
    printf("constraint > 0 means\n");
    printf("*          the current search point is non-feasible).\n");
    printAsterisks(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CobylaOptimizer::~CobylaOptimizer()
{
}

// ************************************************************************
// optimize (this function should be used in the library mode)
// ------------------------------------------------------------------------
void CobylaOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                       double *ubds, int nOutputs, int maxfun, double tol)
{
  psVector vecOptX;
  if (nInputs <= 0)
  {
    printf("CobylaOptimizer ERROR: nInputs <= 0.\n");
    exit(1);
  }
  //**/ create and load up the oData object before calling optimize
  oData *odata = new oData();
  odata->outputLevel_ = 0;
  odata->nInputs_ = nInputs;
  vecOptX.setLength(nInputs);
  odata->optimalX_ = vecOptX.getDVector();
  odata->initialX_ = XValues;
  odata->lowerBounds_ = lbds;
  odata->upperBounds_ = ubds;
  odata->tolerance_ = tol;
  if (odata->tolerance_ <= 0) odata->tolerance_ = 1e-6;
  odata->nOutputs_ = nOutputs;
  odata->outputID_ = 0;
  odata->maxFEval_ = maxfun;
  odata->tolerance_ = tol;
  odata->setOptDriver_ = 0;
  //**/ for this function to work, setObjectiveFunction has to be called
  //**/ beforehand
  odata->optFunction_ = objFunction_;
  odata->funcIO_ = NULL;
  optimize(odata);
  odata->initialX_ = NULL;
  odata->lowerBounds_ = NULL;
  odata->upperBounds_ = NULL;
  odata->optFunction_ = NULL;
  for (int ii = 0; ii < nInputs; ii++)
    XValues[ii] = VecOptX_[ii] = odata->optimalX_[ii];
  odata->optimalX_ = NULL;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void CobylaOptimizer::optimize(oData *odata)
{
  int    nInputs, nConstraints, ii, jj, nIns, nOuts, maxfun;
  int    ntimes=3, printLevel=0, nOutputs, cmpFlag, token;
  double rhobeg=1.0, rhoend=1.0e-4, dtemp;
  char   filename[500], *cString;
  FILE   *fp=NULL;
  psVector vecXVals;

  //**/ ----------------------------------------------------------
  //**/ set up optimization parameters
  //**/ ----------------------------------------------------------
  printLevel = odata->outputLevel_;
  nInputs = odata->nInputs_;
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
  odata->optimalY_ = 1.0e50;
  maxfun = odata->maxFEval_;
  vecXVals.setLength(nInputs+1);
  for (ii = 0; ii < nInputs; ii++) vecXVals[ii] = odata->initialX_[ii];
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
  VecPsCobylaXLbounds_.setLength(nInputs);
  VecPsCobylaXUbounds_.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++) 
  {
    VecPsCobylaXLbounds_[ii] = odata->lowerBounds_[ii];
    VecPsCobylaXUbounds_[ii] = odata->upperBounds_[ii];
  }

  //**/ ----------------------------------------------------------
  //**/ output optimization information
  //**/ ----------------------------------------------------------
  if (psConfig_.InteractiveIsOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("Cobyla optimizer: max fevals   = %d\n",odata->maxFEval_);
    printf("Cobyla optimizer: tolerance    = %e\n",odata->tolerance_);
    printEquals(PL_INFO, 0);
  }
  //**/ Note: this is user-specified number of constraints
  nConstraints = nOutputs - 1;
  //**/ add bound constraints for the inputs
  psCobylaNConstr_ = nConstraints + 2 * nInputs;
  if (psConfig_.InteractiveIsOn())
    printf("Cobyla optimizer: nConstraints = %d\n",psCobylaNConstr_); 

  //**/ ----------------------------------------------------------
  //**/ retrieve history 
  //**/ ----------------------------------------------------------
  cString = psConfig_.getParameter("opt_save_history");
  if (cString != NULL) psCobylaSaveHistory_ = 1;
  cString = psConfig_.getParameter("opt_use_history");
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
        else if (nOuts != nConstraints+1)
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

  //**/ ----------------------------------------------------------
  //**/ ------ call optimizer ------
  //**/ ----------------------------------------------------------
#ifdef HAVE_COBYLA
  int numFuncEvals = odata->numFuncEvals_;
  printLevel = 0;
  cobyla(nInputs, psCobylaNConstr_, vecXVals.getDVector(), rhobeg, 
         rhoend, printLevel, &maxfun, evaluateFunction, (void *) odata);
  VecOptX_.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++) VecOptX_[ii] = odata->optimalX_[ii];
  optimalY_ = odata->optimalY_;
  if (psConfig_.InteractiveIsOn())
  {
    printf("Cobyla optimizer: number of function evaluation = %d\n",
           odata->numFuncEvals_-numFuncEvals);
    printf("                  optimal Y so far = %e\n",optimalY_);
  }
#else
  printf("ERROR : Cobyla optimizer not installed.\n");
  exit(1);
#endif

  //**/ ----------------------------------------------------------
  //**/ ------ set return values and clean up ------
  //**/ ----------------------------------------------------------
  if ((odata->setOptDriver_ & 2) && psCCurrDriver_ >= 0)
  {
    printf("Cobyla: setting back to original simulation driver.\n");
    odata->funcIO_->setDriver(psCCurrDriver_);
  }
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


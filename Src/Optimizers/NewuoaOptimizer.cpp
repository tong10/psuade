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
// Functions for the class NewuoaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>

#include "NewuoaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" 
{
  void newuoa_(int *,int *,double *,double *,double *,int *,int *,double*);
}

// ************************************************************************
// Internal global variables
// ------------------------------------------------------------------------
#define psNewuoaMaxSaved_ 10000
int     psNewuoaNSaved_=0;
double  psNewuoaSaveX_[psNewuoaMaxSaved_*10];
double  psNewuoaSaveY_[psNewuoaMaxSaved_*10];
void    *psNewuoaObj_=NULL;
int     psNCurrDriver_=-1;
int     psNewuoaSaveHistory_=0;
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
  void *newuoaevalfunc_(int *nInps, double *XValues, double *YValue)
  {
    int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
    double *localY, ddata;
    char   pString[1000], lineIn[1000];
    oData  *odata;
    FILE   *infile;

    //**/ ------ fetch data ------
    nInputs = (*nInps);
    odata    = (oData *) psNewuoaObj_;
    nOutputs = odata->nOutputs_;
    localY   = (double *) malloc(nOutputs * sizeof(double));
    outputID = odata->outputID_;

    //**/ ------ search to see if it has already been evaluated ------
    found = 0;
    for (ii = 0; ii < psNewuoaNSaved_; ii++)
    {
      for (jj = 0; jj < nInputs; jj++)
        if (PABS(psNewuoaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) 
            break;
      if (jj == nInputs)
      {
        found = 1;
        printf("Newuoa: simulation results reuse.\n");
        break;
      }
    }

    //**/ ------ run simulation ------
    funcID = odata->numFuncEvals_;
    if (found == 0)
    {
      if (odata->optFunction_ != NULL)
        //**/ if users have loaded a simulation function, use it
        odata->optFunction_(nInputs, XValues, nOutputs, localY);
      else if (odata->funcIO_ != NULL)
        //**/ if not, use the simulation function set up in odata
        odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
      else
      {
        printf("NewuoaOptimizer ERROR: no function evaluator.\n");
        exit(1);
      }
      if (psConfig_.InteractiveIsOn() && odata->outputLevel_ > 1)
      {
        printf("NewuoaOptimizer %6d : \n", odata->numFuncEvals_+1);
        for (ii = 0; ii < nInputs; ii++)
          printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
        printf("    Y = %16.8e\n", localY[outputID]);
      }
      funcID = odata->numFuncEvals_++;
      (*YValue) = localY[outputID];
      if (odata->outputLevel_ > 4)
        printf("    ==> Objective function value = %16.8e\n", (*YValue));
    }
    else
    {
      localY[outputID] = psNewuoaSaveY_[ii];
      (*YValue) = localY[outputID];
    }

    //**/ ------ save the data ------
    if ((psNewuoaSaveHistory_ == 1 && found == 0 &&
        (psNewuoaNSaved_+1)*nInputs < psNewuoaMaxSaved_*10) && 
        psNewuoaNSaved_ < psNewuoaMaxSaved_)
    {
      for (jj = 0; jj < nInputs; jj++)
        psNewuoaSaveX_[psNewuoaNSaved_*nInputs+jj] = XValues[jj];
      psNewuoaSaveY_[psNewuoaNSaved_] = (*YValue);
      psNewuoaNSaved_++;
      infile = fopen("psuade_newuoa_history","w");
      if (infile != NULL)
      {
        for (ii = 0; ii < psNewuoaNSaved_; ii++)
        {
          fprintf(infile, "999 %d ", nInputs);
          for (kk = 0; kk < nInputs; kk++)
            fprintf(infile, "%24.16e ", psNewuoaSaveX_[ii*nInputs+kk]);
          fprintf(infile, "%24.16e\n", psNewuoaSaveY_[ii]);
        }
        fclose(infile);
      }
    }

    //**/ ------ store optimal information ------
    if ((*YValue) < odata->optimalY_)
    {
      odata->optimalY_ = (*YValue);
      for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
      if (psConfig_.OptExpertModeIsOn() && odata->outputLevel_ > 1)
      {
         printf("NewuoaOptimizer %6d : \n", odata->numFuncEvals_);
         for (ii = 0; ii < nInputs; ii++)
            printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
         if (found == 0) printf("    Y = %16.8e\n", localY[outputID]);
         printf("    *** Current best objective function value = %16.8e\n", 
                odata->optimalY_);
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
NewuoaOptimizer::NewuoaOptimizer()
{
  if (psConfig_.InteractiveIsOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("*   NEWUOA Unconstrained Optimization \n");
    printEquals(PL_INFO, 0);
    printf("* NEWUOA is an optimization software ");
    printf("developed by Professor Michael\n");
    printf("* Powell to solve unconstrained optimization ");
    printf("problems without using\n");
    printf("* suitable for problems with hundreds ");
    printf("of variables. The formulation\n");
    printf("* of the problem is as follow:\n");
    printf("\n*                min_X F(X)\n");
    printf("*\n");
    printf("* To run this optimizer, do the following: \n");
    printf("* (1) Prepare a PSUADE input file with ");
    printf("all variable defined.\n");
    printf("* (2) In the PSUADE input file, make sure ");
    printf("opt_driver has been\n");
    printf("*     initialized to point your ");
    printf("optimization objective function\n");
    printf("*     evaluator\n");
    printf("* (3) Since it is unconstrained, you must ");
    printf("provide an initial guess in\n");
    printf("      the PSUADE_IO section, or you can ");
    printf("ask PSUADE to generate initial\n");
    printf("      guesses using one of the available ");
    printf("sampling methods.\n");
    printf("* (4) Set optimization tolerance in PSUADE input file\n");
    printf("* (5) Set maximum number of iterations ");
    printf("in PSUADE input file\n");
    printf("* (6) Set optimization print_level to ");
    printf("give additonal outputs\n");
    printf("* (7) In Opt EXPERT mode, the optimization ");
    printf("history log will be turned\n");
    printf("*     on automatically. Previous psuade_newuoa_history ");
    printf("file will also\n");
    printf("*     be reused to save evaluations.\n");
    printf("* (8) If your opt_driver is a response surface ");
    printf("which has more inputs\n");
    printf("*     than the number of optimization inputs, ");
    printf("you can fix some driver\n");
    printf("*     inputs by using a rs_index_file (ANALYSIS).\n");
    printEquals(PL_INFO, 0);
    printf("To reuse the results (e.g. restart from ");
    printf("abrupt termination), turn on\n");
    printf("save_history and use_history optimization ");
    printf("options in the ANALYSIS\n");
    printf("section (e.g. optimization save_history). ");
    printf("You will see a file called\n");
    printf("'psuade_newuoa_history' afterward.\n");
    printAsterisks(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
NewuoaOptimizer::~NewuoaOptimizer()
{
}

// ************************************************************************
// optimize (this function should be used in the library call mode)
// ------------------------------------------------------------------------
void NewuoaOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                               double *ubds, int maxfun, double tol)
{
  psVector vecOptX;
  if (nInputs <= 0)
  {
    printf("NewuoaOptimizer ERROR: nInputs <= 0.\n");
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
  odata->nOutputs_ = 1;
  odata->outputID_ = 0;
  odata->maxFEval_ = maxfun;
  odata->tolerance_ = tol;
  odata->setOptDriver_ = 0;
  //**/ for this function to work, setObjectiveFunction has to be 
  //**/ called beforehand
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
void NewuoaOptimizer::optimize(oData *odata)
{
  int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0;
  int    nOutputs, outputID, newuoaFlag=0;
  double rhobeg=1.0, rhoend=1.0e-4, dtemp;
  char   cinput[5000], *cString, winput[5001];;
  FILE   *infile=NULL;
  psVector vecXVals, vecWT;

  //**/ ------ prepare object variables ------
  printLevel = odata->outputLevel_;
  nInputs = odata->nInputs_;
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
  odata->optimalY_ = 1.0e50;
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
    printf("NewuoaOptimizer WARNING: tolerance too large.\n");
    printf("                         tolerance reset to 1.0e-6.\n");
    rhoend = rhobeg * 1.0e-6;
  }
  nOutputs = odata->nOutputs_;
  outputID = odata->outputID_;
  maxfun = odata->maxFEval_;
  if ((odata->setOptDriver_ & 1))
  {
    printf("Newuoa: setting optimization simulation driver.\n");
    psNCurrDriver_ = odata->funcIO_->getDriver();
    odata->funcIO_->setDriver(1);
  }
  psNewuoaObj_= (void *) odata;
  if (psConfig_.InteractiveIsOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("Newuoa optimizer: max fevals = %d\n", odata->maxFEval_);
    printf("Newuoa optimizer: tolerance  = %e\n", odata->tolerance_);
    if (printLevel > 1)
      printf("Newuoa optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
    printEquals(PL_INFO, 0);
  }

  //**/ ------ call optimizer ------
  nPts = (nInputs + 1) * (nInputs + 2) / 2;
  kk   = (nPts + 13) * (nPts + nInputs) + 3*nInputs*(nInputs+3)/2;
  vecWT.setLength(kk);

#ifdef HAVE_NEWUOA
  //**/ ----------------------------------------------------------
  //**/ ------ retrieve history ------
  //**/ ----------------------------------------------------------
  cString = psConfig_.getParameter("opt_save_history");
  if (cString != NULL) psNewuoaSaveHistory_ = 1;
  cString = psConfig_.getParameter("opt_use_history");
  if (cString != NULL)
  {
    printf("Newuoa: use history has been turned on.\n");
    infile = fopen("psuade_newuoa_history","r");
    if (infile != NULL)
    {
      psNewuoaNSaved_ = 0;
      while (feof(infile) == 0)
      {
        fscanf(infile, "%d %d", &ii, &kk);
        if (ii != 999 || kk != nInputs) break;
        else
        {
          for (ii = 0; ii < nInputs; ii++) 
            fscanf(infile, "%lg",
                   &psNewuoaSaveX_[psNewuoaNSaved_*nInputs+ii]);
          fscanf(infile, "%lg",&psNewuoaSaveY_[psNewuoaNSaved_]);
          psNewuoaNSaved_++;
        }
        if (((psNewuoaNSaved_+1)*nInputs > psNewuoaMaxSaved_*10) ||
              psNewuoaNSaved_ > psNewuoaMaxSaved_) break;
      } 
      fclose(infile);
    }
  }
  if (psConfig_.InteractiveIsOn())
  {
    for (ii = 0; ii < nInputs; ii++) 
      printf("Newuoa initial X %3d = %e\n", ii+1, vecXVals[ii]);
  }
  int numFuncEvals = odata->numFuncEvals_;
  newuoaFlag = 7777;
#ifdef HAVE_NEWUOA
  newuoa_(&nInputs, &nPts, vecXVals.getDVector(), &rhobeg, &rhoend, 
          &newuoaFlag, &maxfun, vecWT.getDVector());
#else
  printf("Newuoa optimizer ERROR: newuoa not installed.\n");
  exit(1);
#endif

  if (psConfig_.InteractiveIsOn())
  {
    printf("Newuoa optimizer: total number of evaluations = %d\n",
            odata->numFuncEvals_-numFuncEvals);
  }
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = vecXVals[ii];
  VecOptX_.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++) VecOptX_[ii] = odata->optimalX_[ii];
  optimalY_ = odata->optimalY_;

  //**/ ------ save history ------
  if (psNewuoaSaveHistory_ == 1 && psNewuoaNSaved_ > 0)
  {
    infile = fopen("psuade_newuoa_history","w");
    if (infile != NULL)
    {
      for (ii = 0; ii < psNewuoaNSaved_; ii++)
      {
        fprintf(infile, "999 %d ", nInputs);
        for (kk = 0; kk < nInputs; kk++)
          fprintf(infile, "%24.16e ", psNewuoaSaveX_[ii*nInputs+kk]);
        fprintf(infile, "%24.16e\n", psNewuoaSaveY_[ii]);
      }
      fclose(infile);
    }
    printf("Newuoa: history saved in psuade_newuoa_history\n");
  }
#else
  printf("ERROR : Newuoa optimizer not installed.\n");
  exit(1);
#endif

  //**/ ------ set return values and clean up ------
  if ((odata->setOptDriver_ & 2) && psNCurrDriver_ >= 0)
  {
    printf("Newuoa: setting back to original simulation driver.\n");
    odata->funcIO_->setDriver(psNCurrDriver_);
  }
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
NewuoaOptimizer& NewuoaOptimizer::operator=(const NewuoaOptimizer &)
{
   printf("NewuoaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


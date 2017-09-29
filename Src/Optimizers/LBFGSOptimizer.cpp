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
// Functions for the class LBFGSOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "LBFGSOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

#ifdef HAVE_LBFGS
extern "C" {
#include "../../External/L-BFGS-B-C/src/lbfgsb.h"
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
LBFGSOptimizer::LBFGSOptimizer()
{
   printAsterisks(PL_INFO, 0);
   printf("*   LBFGS Optimizer Usage Information\n");
   printEquals(PL_INFO, 0);
   printf("* - To run this optimizer, first make sure opt_driver has\n");
   printf("*   been initialized to point to your optimization objective\n");
   printf("*   function evaluator\n");
   printf("* - Set optimization tolerance in your PSUADE input file\n");
   printf("* - Set maximum number of iterations in PSUADE input file\n");
   printf("* - Set num_local_minima to perform multistart optimization\n");
   printf("* - Set optimization print_level to give additonal outputs\n");
   printf("* - In Opt EXPERT mode, the optimization history log will be\n");
   printf("*   turned on automatically. Previous psuade_lbfgs_history\n");
   printf("*   file will also be reused.\n");
   printf("* - Your opt_driver must have the number of outputs equal to\n");
   printf("*   the number of inputs plus 1 whereby the last nInputs\n");
   printf("*   outputs are the derivatives of the outputs with respect\n");
   printf("*   to each input.\n");
   printAsterisks(PL_INFO, 0);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LBFGSOptimizer::~LBFGSOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void LBFGSOptimizer::optimize(oData *odata)
{
#if HAVE_LBFGS
   int     funcID, ii, kk, nOuts, its, printLevel;
   double  *XValues, *lbounds, *ubounds, *GValues, FValue, *SValues;
   integer nInps, iprint=1, *task=&iprint, lsave[4], isave[44], csave[60];
   integer *iwork, nCorr=5, *nbds;
   double  factr, pgtol, *work, dsave[29]; 

   printLevel = odata->outputLevel_;
   nInps = odata->nInputs_;
   nOuts = odata->nOutputs_;
   if (nInps+1 != nOuts)
   {
      printf("ERROR: LBFGS optimizer needs nOuts = nInps+1.\n");
      exit(1);
   }
   lbounds = odata->lowerBounds_;
   ubounds = odata->upperBounds_;
   XValues = new double[nInps];
   for (ii = 0; ii < nInps; ii++) XValues[ii] = odata->initialX_[ii];
   funcID = odata->numFuncEvals_;

   odata->optimalY_ = 1.0e50;
   odata->funcIO_->setDriver(1);
   GValues = new double[nInps];
   SValues = new double[nInps+1];
   for (ii = 0; ii < nInps; ii++) GValues[ii] = 0.0;
   for (ii = 0; ii < nInps+1; ii++) SValues[ii] = 0.0;
   nbds    = new integer[nInps];
   for (ii = 0; ii < nInps; ii++) nbds[ii] = 2;
   factr = 1e7;
   pgtol = 1e-6;
   kk = (2 * nCorr + 5) * nInps + 11 * nCorr * nCorr + 8 * nCorr;
   work  = new double[kk];
   iwork = new integer[3*nInps];
   *task = (integer) START;
   its   = 0;
 
   while (1)
   {
      its++;
      setulb(&nInps, &nCorr, XValues, lbounds, ubounds, nbds, &FValue,
          GValues, &factr, &pgtol, work, iwork, task, &iprint, csave,
          lsave, isave, dsave);
      if (IS_FG(*task))
      {
         if (printLevel > 0)
            for (ii = 0; ii < nInps; ii++)
               printf("Current Input X %4d = %24.16e\n",ii+1,XValues[ii]);
         odata->funcIO_->evaluate(funcID,nInps,XValues,nOuts,SValues,0);
         funcID = odata->numFuncEvals_++;
         FValue = SValues[0];
         for (ii = 0; ii < nInps; ii++) GValues[ii] = SValues[ii+1];
      }
      else if (*task != NEW_X)
      {
         if (printLevel > 0)
            for (ii = 0; ii < nInps; ii++)
               printf("Final Input X %4d = %24.16e\n", ii+1, XValues[ii]);
         printf("Final objective function value       = %16.8e\n",FValue);
         printf("Total number of function evaluations = %d\n", its);
         break;
      }
   }

   for (ii = 0; ii < nInps; ii++) odata->optimalX_[ii] = XValues[ii];
   odata->optimalY_ = FValue;
   delete [] work;
   delete [] iwork;
   delete [] nbds;
   delete [] XValues;
   delete [] GValues;
#else
   printf("ERROR : LBFGS optimizer not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
LBFGSOptimizer& LBFGSOptimizer::operator=(const LBFGSOptimizer &)
{
   printf("LBFGSOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


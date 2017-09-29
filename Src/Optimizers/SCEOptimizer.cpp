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
// Functions for the class SCEOptimizer
// AUTHOR : Amer Abdulla 
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeUtil.h"
#include "SCEOptimizer.h"
#include "PsuadeUtil.h"
#include "Sampling.h"
#include "sysdef.h"
#include "Psuade.h"

// ------------------------------------------------------------------------
#include <math.h> // for standev and georange functions
#include <time.h> // for random number generator
// ------------------------------------------------------------------------

#define psSCEMaxSaved_ 100000
int     psSCENSaved_=0;
double  psSCESaveX_[psSCEMaxSaved_*10];
double  psSCESaveY_[psSCEMaxSaved_];
void    *psSCEObj_=NULL;

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
   void *SCEevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
      double *localY;
      oData  *odata;
      FILE   *infile;

      nInputs = (*nInps);
      if (psOptExpertMode_ == 1 && psSCENSaved_ == 0)
      {
         infile = fopen("psuade_sce_data","r");
         if (infile != NULL)
         {
            fscanf(infile, "%d %d", &psSCENSaved_, &kk);
            if ((psSCENSaved_ <= 0) ||
                (psSCENSaved_+1 > psSCEMaxSaved_*10/nInputs))
            {
               printf("PSUADE SCE: history file has too much data.\n");
               printf("            Ask PSUADE developer for help.\n"); 
               fclose(infile);
               psSCENSaved_ = 0;
            }
            else if (kk != nInputs)
            {
               printf("PSUADE SCE: history file has invalid input count.\n");
               fclose(infile);
               psSCENSaved_ = 0;
            }
            else
            {
               for (ii = 0; ii < psSCENSaved_; ii++)
               {
                  fscanf(infile, "%d", &kk);
                  if (kk != ii+1)
                  {
                     printf("PSUADE SCE: data index mismatch.\n");
                     psSCENSaved_ = 0;
                     break;
                  }
                  for (jj = 0; jj < nInputs; jj++)
                     fscanf(infile, "%lg", &psSCESaveX_[ii*nInputs+jj]);
                  fscanf(infile, "%lg", &psSCESaveY_[ii]);
               }
               fclose(infile);
            }
         }
      }

      odata    = (oData *) psSCEObj_;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));
      outputID = odata->outputID_;

      found = 0;
      for (ii = 0; ii < psSCENSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psSCESaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
           odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
           funcID = odata->numFuncEvals_++;
      }
      else localY[outputID] = psSCESaveY_[ii];
      (*YValue) = localY[outputID];

      if (psOptExpertMode_ != 0 && found == 0)
      {
         for (jj = 0; jj < nInputs; jj++)
            psSCESaveX_[psSCENSaved_*nInputs+jj] = XValues[jj];
         psSCESaveY_[psSCENSaved_] = localY[outputID];
         psSCENSaved_++;
      }
      // Use for Diagnostics
      //for (ii = 0; ii < nInputs; ii++)
      //printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
      //printf("    Y = %16.8e\n", localY[outputID]);

      if ((*YValue) < odata->optimalY_)
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
         if (odata->outputLevel_ > 2)
         {
            printf("SCEOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            printf("    Ymin  = %16.8e\n", odata->optimalY_);
         }
      }
      free(localY);

      if (psOptExpertMode_ != 0 && found == 0)
      {
         infile = fopen("psuade_bobyqa_data","w");
         if (infile != NULL)
         {
            fprintf(infile, "%d %d\n", psSCENSaved_, nInputs);
            for (ii = 0; ii < psSCENSaved_; ii++)
            {
               fprintf(infile, "%d ", ii+1);
               for (jj = 0; jj < nInputs; jj++)
                  fprintf(infile, "%24.16e ", psSCESaveX_[ii*nInputs+jj]);
               fprintf(infile, "%24.16e\n", psSCESaveY_[ii]);
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
SCEOptimizer::SCEOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SCEOptimizer::~SCEOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void SCEOptimizer::optimize(oData *odata) 
{
   int    nInputs, printLevel=0, ii, jj,zz, maxfun, currDriver, nPts;
   int    nComplex=2, includeInitialPoint=0, maxEvoLoop = 10;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, sum, Mean;
   double paramSpaceConvergence = 0.001, percentChange = .1;

   nInputs = odata->nInputs_;
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   maxfun = odata->maxFEval_;
   odata->optimalY_ = 1.0e50;
   odata->numFuncEvals_ = 0;
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
      printf("SCEOptimizer WARNING: tolerance too large.\n");
      printf("                      tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   currDriver = odata->funcIO_->getDriver();
   odata->funcIO_->setDriver(1);
   psSCEObj_= (void *) odata;

   // Initialize SCE parameters
   int nMemPerComplex = 2 * nInputs + 1; 
   int nMemPerSimplex = nInputs + 1;
   int nEvoStep = nMemPerComplex;
   nPts = nMemPerComplex * nComplex;

   // For each parameter, determine range of points from which we can sample
   double* bound = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
     bound[ii] = odata->upperBounds_[ii] - odata->lowerBounds_[ii];

   // Create an initial population to fill array x1[npt*nInputs]
   // Use PSUADE sampling method of choice
   Sampling *sampler;
   sampler = SamplingCreateFromID(PSUADE_SAMP_MC);
   sampler->setInputBounds(nInputs,odata->lowerBounds_,odata->upperBounds_);
   sampler->setSamplingParams (nPts, 1, 0);
   sampler->setOutputParams (1);
   sampler->initialize(0);
   nPts = sampler -> getNumSamples();
   double *x1 = new double[nPts*nInputs];
   double *xf = new double[nPts];
   int    *S = new int[nPts];
   sampler -> getSamples(nPts, nInputs, 1, x1, xf, S);
   delete [] S;
   S = NULL;

   if ( includeInitialPoint == 1 )
     for (ii = 0; ii < nInputs; ++ii) x1[ii] = odata->initialX_[ii];

   int iCall = 0;
   int nLoop = 0;
   double *xCoor = new double[nInputs];
   for (ii = 0; ii < nPts; ++ii)
   {
      for (jj = 0; jj < nInputs; ++jj) xCoor[jj]=x1[ii*nInputs+jj];
      SCEevalfunc_(&nInputs, xCoor, &xf[ii]);
      iCall += 1;
   }

   int *sortArray = new int[nPts];
   for (ii = 0; ii < nPts; ++ii) sortArray[ii] = ii;
   sortDbleList2a(nPts, xf, sortArray); 

   double *x = new double[nPts*nInputs];
   for (ii = 0; ii != nPts; ++ii)
      for (jj = 0; jj != nInputs; ++jj)
         x[ii*nInputs + jj] = x1[sortArray[ii]*nInputs + jj];
   delete [] x1;
   x1 = NULL;
   delete [] sortArray;
   sortArray = NULL;

   double bestx[nInputs]; 
   double worstx[nInputs]; 
   for (ii = 0; ii < nInputs; ++ii)
   {
      bestx[ii] = x[ii];
      worstx[ii] = x[(nPts-1)*nInputs+ii];
   }

   double bestf = xf[0];
   double worstf = xf[nPts -1];

   double geoRange = georange(nInputs, nPts, x, bound);


   // Check for convergency 
   if (odata->numFuncEvals_ >= maxfun)
   {
      printf("Optimization search terminated because the limit on the\n");
      printf("maximum number of trials %d has been exceeded. Search\n",maxfun);
      printf("was stopped at trial number %d of the initial loop!\n", iCall);
   }

   if (geoRange < paramSpaceConvergence)
      printf("Population has converged to a small parameter space.\n");

   if (printLevel > 0)
   {
      printf("SCE optimizer: max fevals = %d\n", odata->maxFEval_);
      printf("SCE optimizer: tolerance  = %e\n", odata->tolerance_);
      printf("SCE optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
   }

   odata->funcIO_->setDriver(currDriver);
   delete [] XValues;

   double criter[1000];
   double criter_change = 100000;
   int    ic, k2[nMemPerComplex], flag, iter,iz;
   double *cx = new double[nMemPerComplex*nInputs];
   double *cf = new double[nMemPerComplex];
   int    *simIndex = new int[nMemPerSimplex];
   double simPosition;
   double *s = new double[nMemPerSimplex*nInputs];
   double *sf = new double[nMemPerSimplex];
   double *sNew = new double[nInputs];
   double sfNew;
   double* cx2 = new double[nMemPerComplex*nInputs];
   double* x2 = new double[nInputs*nPts];
   sortArray = new int[nPts];

   while ((odata->numFuncEvals_ < maxfun) && 
          (geoRange > paramSpaceConvergence) && 
          (criter_change > percentChange)) 
   {
      nLoop += 1;

      for (ic = 0; ic < nComplex; ++ic)
      {
         for (ii = 0; ii < nMemPerComplex; ++ii)
            k2[ii] = (ii * nComplex) + ic;
       
         for (ii = 0; ii < nMemPerComplex; ++ii)
         {
            cf[ii] = xf[k2[ii]];
            for (jj = 0; jj < nInputs; ++jj)
              cx[ii*nInputs + jj] = x[k2[ii]*nInputs+jj];
         }

         for (ii = 0; ii < nEvoStep; ++ii)
         {
            simIndex[0] = 0; 
            flag = 0;
            for (jj = 1; jj < nMemPerSimplex; ++jj)
            {
               for (iter = 0; iter < 1000; ++iter)
               {
                  simPosition = floor(nMemPerComplex+0.5-sqrt((nMemPerComplex + 
                                0.5)*(nMemPerComplex+0.5)-nMemPerComplex * 
                                (nMemPerComplex + 1)*PSUADE_drand()));
 	          for (iz = 0; iz <= jj-1; ++iz)
                  {
 	             if (simIndex[iz] == simPosition)
                     {
                         flag = 1;
                         break;
                     }
 	          }
 	          if (flag != 1) break;
 	       }
 	       simIndex[jj] = (int) simPosition;
            }
            sortIntList(nMemPerSimplex, simIndex);
	
            for (jj = 0; jj != nMemPerSimplex; ++jj)
            {
	       sf[jj] = cf[simIndex[jj]];
	       for (zz = 0; zz != nInputs; ++zz)
	          s[jj*nInputs + zz] = cx[simIndex[jj]*nInputs+zz];
            }
            newPoint(&sfNew,&iCall,nInputs,nMemPerSimplex,maxfun, 
                     s, sf, sNew, odata);
      
            for (jj = 0; jj < nInputs; ++jj)
               s[nInputs*(nMemPerSimplex-1) + jj] = sNew[jj];
            sf[nMemPerSimplex-1] = sfNew;  

            for (jj = 0; jj != nMemPerSimplex; ++jj)
            {
	       cf[simIndex[jj]] = sf[jj];
	       for (zz = 0; zz != nInputs; ++zz)
	          cx[simIndex[jj]*nInputs+zz] = s[jj*nInputs + zz];
            }

            for (jj = 0; jj < nPts; ++jj) sortArray[jj] = jj;

            sortDbleList2a(nMemPerComplex, cf, sortArray);

            for (jj = 0; jj != nMemPerComplex*nInputs; ++jj) cx2[jj] = cx[jj];

            for (jj = 0; jj != nMemPerComplex; ++jj)
            {
	       for (zz = 0; zz != nInputs; ++zz)
	          cx[jj*nInputs + zz] = cx2[sortArray[jj]*nInputs + zz];
            }
         }
	  
      
         for (ii = 0; ii != nMemPerComplex; ++ii)
         {
            xf[k2[ii]] = cf[ii];
	    for (jj = 0; jj != nInputs; ++jj)
	       x[k2[ii]*nInputs+jj] = cx[ii*nInputs + jj];
         }
      }

      for (ii = 0; ii < nPts; ++ii) sortArray[ii] = ii;
      sortDbleList2a(nPts, xf, sortArray); 

      for (ii = 0; ii != nPts*nInputs; ++ii) x2[ii] = x[ii];

      for (ii = 0; ii != nPts; ++ii)
      {
         for (jj = 0; jj != nInputs; ++jj)
            x[ii*nInputs + jj] = x2[sortArray[ii]*nInputs + jj];
      }

      // Record the best and worst points
      for (ii = 0; ii != nInputs; ++ii)
      {
          bestx[ii] = x[ii];
          worstx[ii] = x[(nPts-1)*nInputs+ii];
      }

      bestf = xf[0];
      worstf = xf[nPts -1];

      geoRange = georange(nInputs, nPts, x,  bound);


      if (odata->numFuncEvals_ >= maxfun)
      {
         printf("*** Optimization search terminated because the limit on\n");
	 printf("maximum number of trials %d has been exceeded.\n", maxfun);
      
      }
      if (geoRange < paramSpaceConvergence)
         printf("The population has converged to a small parameter space.\n");
    
      criter[nLoop-1] = bestf;
      if (nLoop >= maxEvoLoop)
      {
         if (bestf != criter[nLoop - maxEvoLoop])
         {
	    criter_change = fabs(criter[nLoop-1]-criter[nLoop-maxEvoLoop])*100;
            sum = 0.0; 
	    for (ii = nLoop-maxEvoLoop; ii != nLoop-1; ++ii)
	       sum += fabs(criter[ii]);
	    Mean = sum/maxEvoLoop;
	    criter_change = criter_change/Mean;
	    if (criter_change < percentChange)
            {
	       printf("The best point has improved in last %d loops by less\n",
                      maxEvoLoop);
	       printf("than the threshold %e => convergence achieved.\n",
                      percentChange);
	    }
         }
      }
   } 
   delete [] sortArray;
   delete [] cx2;
   delete [] cx;
   delete [] cf;
   delete [] x2;
 
   for (ii = 0; ii < 60; ii++) printf("*");
   printf("\n");
   printf("Search was stopped at trial number %d\n", iCall);
   printf("Normalized geometric range = %e\n", geoRange);
   printf("Calculated global optimum = %e\n", bestf);

}

// ***********************************************************************
// georange function: calculates geometric range of each parameter
// -----------------------------------------------------------------------
double SCEOptimizer::georange(int nInputs,int nPts,double* x,double* bound)
{
   int    ii, jj;
   double y[nInputs], maxval, minval, sum, mean;
   for (ii = 0; ii != nInputs; ++ii)
   {
      maxval = x[ii];
      minval = x[ii];
      for (jj = ii; jj < (nPts*nInputs); jj += nInputs)
      {
         if      (x[jj] > maxval) maxval = x[jj];
         else if (x[jj] < minval) minval = x[jj];
      }
      y[ii] = maxval - minval;
   }

   for (ii = 0; ii != nInputs; ++ii) y[ii] = y[ii] / bound[ii];

   for (ii = 0; ii != nInputs; ++ii) y[ii] = log(y[ii]);

   sum = 0.0;
   for (ii = 0; ii != nInputs; ++ii) sum += y[ii];
   mean = sum / nInputs;
   return exp(mean);
}

// **********************************************************************
// newPoint function: generates a new point in a simplex
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
void SCEOptimizer::newPoint(double *sfNew, int* iCall, int nInputs,
                            int nMemPerSimplex, int maxfun, double* s, 
                            double* sf, double* sNew, oData *odata)
{
   int    ii, jj, outOfBound = 0;
   double alpha = 1.0, sum;
   double beta = 0.5;
   double sBest[nInputs]; 
   double sWorst[nInputs]; 
   double sfNew1;
   double centroid[nInputs];

   for (ii = 0; ii !=nInputs; ++ii)
   {
      sBest[ii] = s[ii];
      sWorst[ii] = s[nInputs*(nMemPerSimplex-1)+ii];
   }
   double sfWorst = sf[nMemPerSimplex-1];

   for (ii = 0; ii != nInputs; ++ii)
   { 
      sum = 0.0;
      for (jj = ii; jj < ((nMemPerSimplex-1)*nInputs); jj += nInputs)
         sum += s[jj];
      centroid[ii]= sum/(nMemPerSimplex-1);
   }

   for (ii = 0; ii != nInputs; ++ii)
      sNew[ii] = centroid[ii] + alpha *(centroid[ii] - sWorst[ii]);

   double s1[nInputs]; // double* s1 = new double[nInputs];
   double s2[nInputs]; // double* s2 = new double[nInputs]; 
   for (ii = 0; ii != nInputs; ++ii)
   {
      s1[ii] = sNew[ii] - odata->lowerBounds_[ii];
      s2[ii] = odata->upperBounds_[ii] - sNew[ii];
   }

   for (ii = 0; ii != nInputs; ++ii)
   {
      if ((s1[ii]<0) || (s2[ii]<0))
      {
         outOfBound = 1;
         break;
      }
   }
   if (outOfBound == 1)
   {
      for (ii = 0; ii != nInputs ; ++ii)
         sNew[ii] = odata->lowerBounds_[ii] + PSUADE_drand() *
                    (odata->upperBounds_[ii] - odata->lowerBounds_[ii]);
   }
   SCEevalfunc_(&nInputs, sNew, &sfNew1);
   *iCall += 1;

   if (sfNew1 > sfWorst)
   {
      for (ii = 0; ii != nInputs; ++ii)
         sNew[ii] = sWorst[ii] + beta *(centroid[ii] - sWorst[ii]);
      SCEevalfunc_(&nInputs, sNew, &sfNew1);
      *iCall += 1;

      if (sfNew1 > sfWorst)
      {
         for (ii = 0; ii != nInputs ; ++ii)
            sNew[ii] = odata->lowerBounds_[ii] + PSUADE_drand() *
                       (odata->upperBounds_[ii] - odata->lowerBounds_[ii]);
         SCEevalfunc_(&nInputs, sNew, &sfNew1);
         *iCall += 1;
      }
   }
   *sfNew = sfNew1;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
SCEOptimizer& SCEOptimizer::operator=(const SCEOptimizer &)
{
   printf("SCEOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


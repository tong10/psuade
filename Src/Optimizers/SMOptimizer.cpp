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
// Functions for the class SMOptimizer
// AUTHOR : David Echeverria Ciaurri - Charles Tong
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "SMOptimizer.h"
#ifdef HAVE_COBYLA
#include "../../External/COBYLA/cobyla.h"
#endif

double *px0, *B0, *x0, *fx, *y, normy;
int coarseoptimtype = 1;
//void savesmaux(int, int, double *, double *, double *);

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
int evaluateFunctionSM(int nInputs, int nConstraints, double *XValues,
                       double *YValue, double *constraints, void *data)
{
   int    i, j, funcID, nOutputs, outputID, currDriver;
   double *localY, *aux; /*, *localX2, *localY2, deltaX;*/
   oData  *odata;
   //extern double *px0, *B0, *x0, *fx, y[4], normy;
   //extern int coarseoptimtype;
   double *xaux;


   odata    = (oData *) data;
   nOutputs = odata->nOutputs_;
   localY   = (double *) malloc(nOutputs * sizeof(double));
   aux      = (double *) malloc( nInputs *sizeof(double));
   xaux     = (double *) malloc( nInputs *sizeof(double));
   funcID   = odata->numFuncEvals_++;
   outputID = odata->outputID_;

   // check
   /*for (i = 0; i < nInputs; i++)
         printf("%lg ", px0[i]);
         printf("\n");
     for (i = 0; i < nInputs; i++){
        for(j = 0; j < nInputs; j++)
           printf("%lg ", B0[j + i*nInputs]);
           printf("\n");
      }
      for (i = 0; i < nInputs; i++)
         printf("%lg ", x0[i]);
      printf("\n");*/



   for (i=0; i < nInputs; i++) xaux[i] = XValues[i];

   if (coarseoptimtype == 1)      // modify values
   {
      /* Computing pxk + Bk * (S - xk) */
      for (i = 0; i < nInputs; i++)
      {
         aux[i] = 0.0;
         for (j = 0; j < nInputs; j++)
            aux[i] = aux[i] + B0[i + j*nInputs]*(xaux[j]-x0[j]);
      };
      for (i = 0; i < nInputs; i++) xaux[i] = px0[i] + aux[i];
   }


   currDriver = odata->funcIO_->getDriver();
   odata->funcIO_->setDriver(2);      /* we always optimize the coarse model */
   odata->funcIO_->evaluate(funcID,nInputs,xaux,nOutputs,localY,0);
   odata->funcIO_->setDriver(currDriver);

   (*YValue) = 0.0;
   if (coarseoptimtype == 1){       // optim
      for (i=0; i<nOutputs; i++)
         (*YValue) += (localY[i] - y[i])*(localY[i] - y[i]);
      (*YValue) = sqrt((*YValue))/normy*100.00;
   }
   else
   {
      for (i=0; i<nOutputs; i++)
         (*YValue) += (localY[i] - fx[i])*(localY[i] - fx[i]);
      (*YValue) = sqrt((*YValue));
   }



   if ((*YValue) < odata->optimalY_)
   {
      odata->optimalY_ = (*YValue);
      for (i = 0; i < nInputs; i++) odata->optimalX_[i] = XValues[i];
      if (odata->outputLevel_ > 0)
         printf("CobylaOptimizer %6d : Ymin = %16.8e\n",
                odata->numFuncEvals_, odata->optimalY_);
   }
   for (i = 0; i < nInputs; i++)
   {
      constraints[i] = XValues[i] - odata->lowerBounds_[i];
      constraints[i+nInputs] = odata->upperBounds_[i] - XValues[i];
   }
   free(localY);
   free(aux);
   return 0;
}
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
SMOptimizer::SMOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SMOptimizer::~SMOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void SMOptimizer::optimize(oData *odata)
{
   int    nInputs, nOutputs, nConstraints, i,j, maxfun=1000;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp;
   double *zstar, *px1, *B1, *x1, *B0h;
   int    coboutput = 0, nIter = 1, nIterMAX = 10, ncevals = 0, nfevals = 0;
   double Fx, h = 1e50, tolX = 1e-4;
   extern double *px0, *B0, *x0, *fx, *y, normy;
   extern int coarseoptimtype;

/* --- VOY AQUI ------------------------------------- *
   odata->targetFile_ = "./specs";
   FILE *fspecs = fopen(odata->targetFile_, "r");
 * -------------------------------------------------- */

   /*FILE  *ffineaux    = fopen(  "fineaux", "rw"); */
   /*FILE  *fcoarseaux  = fopen("coarseaux", "rw"); */

   nInputs = odata->nInputs_;
   nOutputs = odata->nOutputs_;
   nConstraints = 2 * nInputs;
   for (i = 0; i < nInputs; i++) odata->optimalX_[i] = 0.0;
   odata->numFuncEvals_ = 0;
   odata->optimalY_     = 1.0e50;
   tolX                 = odata->tolerance_;

   if (!strcmp(odata->targetFile_,"NONE"))
   {
      printf("SMOptimizer ERROR : no target file.\n");
      exit(1);
   }
   FILE *fspecs = fopen(odata->targetFile_, "r");

   y      = (double *) malloc(nOutputs * sizeof(double));
   for (i = 0; i < nOutputs; i++) fscanf(fspecs, "%lg", &y[i]);
   fclose(fspecs);

   normy = 0.0;
   for (i=0; i < nOutputs; i++) normy = normy+ y[i]*y[i];
   normy = sqrt(normy);

   printf("\n\n");
   printf(" ****************************************************\n");
   printf(" SPACE MAPPING: Specifications: \n");
   for(i=0; i<nOutputs; i++)
   printf("                                      %20.14f\n", y[i]);
   printf(" ****************************************************\n");

   XValues = (double *) malloc((nInputs+1) * sizeof(double));

   /* sm auxiliary variables: initialization + storing in file smaux */
   fx      = (double *) malloc(nOutputs * sizeof(double));
   x0      = (double *) malloc((nInputs+1) * sizeof(double));
   px0     = (double *) malloc((nInputs+1) * sizeof(double));
   B0      = (double *) malloc(nInputs * nInputs * sizeof(double));

   x1      = (double *) malloc((nInputs+1) * sizeof(double));
   px1     = (double *) malloc((nInputs+1) * sizeof(double));
   B1      = (double *) malloc(nInputs * nInputs * sizeof(double));

   B0h     = (double *) malloc(nInputs * sizeof(double));

   zstar   = (double *) malloc(nInputs * sizeof(double));

   /* Cobyla STOP Criteria   */
   rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
   for (i = 0; i < nInputs; i++)
   {
      dtemp = odata->upperBounds_[i] - odata->lowerBounds_[i];
      if (dtemp > rhobeg) rhobeg = dtemp;
   }

   rhobeg *= 0.05;                    /* orig: 0.5    */
   rhoend = rhobeg * 0.00001/2;         /* orig: 0.0001 */

   /* Initializing state variables */

   for (i = 0; i < nInputs; i++){
      x0[i]  = 0.0;
      px0[i] = 0.0;
      for (j = 0; j < nInputs; j++)
         if (i==j) B0[i + j*nInputs] = 1.0;
	 else      B0[i + j*nInputs] = 0.0;
   };

   /* next is optimization => 1 */
   //savesmaux(1, nInputs, px0, B0, x0);

////--------------------------------------------------------------------
//   odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
//   odata->funcIO_->evaluate(odata->numFuncEvals_++,2,x0,nOutputs,Fx,0);
////--------------------------------------------------------------------

   /* I. INITIAL GUESS */
   for (i = 0; i < nInputs; i++) XValues[i] = odata->initialX_[i];

   printf(" ****************************************************\n");
   printf(" SPACE MAPPING: Initial Guess: \n");
   for(i=0; i<nInputs; i++)
   printf("                                      %20.14f\n", XValues[i]);
   printf(" ****************************************************\n");

   /* 0. COARSE MODEL OPTIMIZATION */

   // ------ call optimizer ------

#ifdef HAVE_COBYLA
   cobyla(nInputs, nConstraints, XValues, rhobeg, rhoend,
          coboutput, &maxfun, evaluateFunctionSM,
          (void *) odata);
#else
   printf("ERROR : Cobyla optimizer not installed.\n");
   exit(1);
#endif
   ncevals += maxfun;
   coarseoptimtype = 0;

/*   zstar = XValues; */
   for (i=0;i<nInputs;i++){
      x0[i]       = XValues[i];
      zstar[i]    = XValues[i];
   }

   printf("\n\n");
   printf(" ****************************************************\n");
   printf(" SPACE MAPPING: Coarse model optimum: \n");
   for(i=0; i<nInputs; i++)
   printf("                                      %20.14f\n", x0[i]);
   printf(" ****************************************************\n");
   printf("\n\n");

   /* 1. EVALUATING THE FINE MODEL AT THE COARSE OPTIMUM */

   /*currDriver = odata->funcIO_->getDriver(); */
   odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
   odata->funcIO_->evaluate(odata->numFuncEvals_++,nInputs,x0,nOutputs,fx,0);

   Fx = 0.0;
   for (i=0; i<nOutputs; i++)
      Fx += (fx[i] - y[i])*(fx[i] - y[i]);
   Fx = sqrt(Fx)/normy*100.00;

   nfevals += 1;
   /* odata->funcIO_->setDriver(currDriver); */
   /* (*YValue) = localY[outputID]; */
   /*printf("                                      %20.14f\n", Fx); */
   printf("   Iteration      # f evals   # c eval   Fine Cost Function    ||x_{k+1} - x_k||\n");
   printf(" ---------------------------------------------------------------------------------\n");
   printf("   %5d           %5d       %5d         %8.3f\n",nIter,nfevals,ncevals,Fx);


   /* Now f(x0) is stored in the file fileaux (format: m f(x0)) */

   /* 2. EVALUATING THE SPACE-MAPPING FUNCTION AT x0            */
   //savesmaux(0, nInputs, px0, B0, x0);       /* parameter extraction */
   for (i=0;i<nInputs;i++)
      px0[i]       = x0[i];         /* makes sense for parameter extraction */
#ifdef HAVE_COBYLA
   maxfun=1000;
   cobyla(nInputs, nConstraints, px0, rhobeg, rhoend,
          coboutput, &maxfun, evaluateFunctionSM, (void *) odata);
#else
   printf("ERROR : Cobyla optimizer not installed.\n");
   exit(1);
#endif
   ncevals += maxfun;
   coarseoptimtype = 1;
   //savesmaux(1, nInputs, px0, B0, x0);

while (1){

   /* 2. OPTIMIZING THE MAPPED COARSE MODEL          */

   for (i=0;i<nInputs;i++)
      x1[i]       = x0[i];         /* makes sense for parameter extraction */


#ifdef HAVE_COBYLA
   maxfun=1000;
   cobyla(nInputs, nConstraints, x1, rhobeg, rhoend,
          coboutput, &maxfun, evaluateFunctionSM, (void *) odata);
#else
   printf("ERROR : Cobyla optimizer not installed.\n");
   exit(1);
#endif
   ncevals += maxfun;
   coarseoptimtype = 0;

   //printf("x = %10.6f   %10.6f\n",  x1[0], x1[1]);

   /* 3. EVALUATING THE FINE MODEL AT x1 */

   odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
   odata->funcIO_->evaluate(odata->numFuncEvals_++,nInputs,x1,nOutputs,fx,0);
   Fx   = 0.0;
   for (i=0; i<nOutputs; i++)
      Fx   += (fx[i] - y[i])*(fx[i] - y[i]);
   Fx   = sqrt(Fx)/normy*100.00;
   nfevals += 1;
   /*printf("                                      %20.14f\n", Fx);*/

   /* Now f(x1) is stored in the file fileaux (format: m f(x1)) */

   /* S. STOPPING CRITERIA              */
   h = 0.0;
   for (i=0;i<nInputs;i++)
      h += (x1[i]-x0[i])*(x1[i]-x0[i]);
   h = sqrt(h);
   /*printf(" h = %10.6f\n", h); */
   printf("   %5d           %5d       %5d         %8.3f              %4.2e \n", nIter+1,nfevals,ncevals,Fx,h);
   if (nIter == nIterMAX || h < tolX)
      break;


   /* 4. EVALUATING THE SPACE-MAPPING FUNCTION AT x1            */
   //savesmaux(0, nInputs, px0, B0, x0);       /* parameter extraction */
   for (i=0;i<nInputs;i++)
      px1[i]       = px0[i];         /* makes sense for parameter extraction */
#ifdef HAVE_COBYLA
   maxfun=1000;
   cobyla(nInputs, nConstraints, px1, rhobeg, rhoend,
          coboutput, &maxfun, evaluateFunctionSM, (void *) odata);
#else
   printf("ERROR : Cobyla optimizer not installed.\n");
   exit(1);
#endif
   ncevals += maxfun;
   coarseoptimtype = 1;

   /*savesmaux(1, nInputs, px0, B0, x0);*/

   nIter = nIter + 1;

   /* 5. UPDATING THE APPROX. OF p */
   /*    B1 = B0 + (px1-px0-B0*h)/(h'*h)*h' */

   for (i=0;i<nInputs;i++){
      B0h[i] = 0.0;
      for (j=0;j<nInputs;j++)
         B0h[i] += B0[i+j*nInputs]*(x1[j]-x0[j]);
   }

   for (i=0;i<nInputs;i++)
      for (j=0;j<nInputs;j++)
         B1[i+j*nInputs] = B0[i+j*nInputs] + 
                           (px1[i]-px0[i]-B0h[i])/(h*h)*(x1[j]-x0[j]);

/* - CHECKINGS ------------------------------------------------------ */
/*

   for (i=0;i<nInputs;i++){
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", B0[i+j*nInputs]);
      printf("\n");
   }
      printf("\n");


   for (i=0;i<nInputs;i++){
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", B1[i+j*nInputs]);
      printf("\n");
   }
      printf("\n");


      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", x0[j]);
      printf("\n");
      printf("\n");
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", x1[j]);
      printf("\n");
      printf("\n");
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", px0[j]);
      printf("\n");
      printf("\n");
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", px1[j]);
      printf("\n");
      printf("\n");
      for (j=0;j<nInputs;j++)
         printf("%10.6f   ", B0h[j]);
      printf("\n");
      printf("\n"); 
*/

   /* updating variables */

   for (i=0;i<nInputs;i++){
        x0[i]  = x1[i];
	px0[i] = px1[i];
	for (j=0;j<nInputs;j++)
           B0[i+j*nInputs] = B1[i+j*nInputs];
   }

   //savesmaux(1, nInputs, px0, B0, x0);
/*   x1  = x0;
   B1  = B0;
   px1 = px0;*/
}

/* DIAGNOSTICS */
   if (nIter == nIterMAX)
      printf("\n Maximum number of %4d SM iterations reached.\n\n", nIterMAX);
   else
      printf("\n Normal termination.\n\n");

   printf("\n\n");
   printf(" ****************************************************\n");
   printf(" SPACE MAPPING: Space-mapping solution:\n");
   for(i=0; i<nInputs; i++)
   printf("                                      %20.14f\n", x1[i]);
   printf(" ****************************************************\n");
   printf("\n\n");


   odata->optimalY_ = Fx;
   for (i = 0; i < nInputs; i++) odata->optimalX_[i] = x1[i];

   free(x0);
   free(px0);
   free(B0);
   free(x1);
   /*free(px1);
   free(B1);*/
   free(zstar);
   free(XValues);
   free(fx);
   free(y);

   /*fclose(ffineaux);*/
   /*fclose(fcoarseaux);*/
   // ------ set return values and clean up ------
}

//void savesmaux(int coarsetype,int nInputs,double *px0,double *B0,double *x0)
//{
//   int i,j;
//   FILE  *fsmaux      = fopen(    "smaux", "w");
//
//
//if (coarsetype == 0){      /* parameter extraction */
//
//   fprintf(fsmaux, "%d ", 0);
//   fprintf(fsmaux, "\n");
//   for (i = 0; i < nInputs; i++)
//      fprintf(fsmaux, "%lg ", 0.0);
//   fprintf(fsmaux, "\n");
//   for (i = 0; i < nInputs; i++){
//      for(j = 0; j < nInputs; j++)
//         if (i==j)
//            fprintf(fsmaux, "%lg ", 1.0);
//         else
//            fprintf(fsmaux, "%lg ", 0.0);
//      fprintf(fsmaux, "\n");
//   }
//   for (i = 0; i < nInputs; i++)
//      fprintf(fsmaux, "%lg ", 0.0);
//   fprintf(fsmaux, "\n");
//}
//else{
//
//   fprintf(fsmaux, "%d ", 1);
//   fprintf(fsmaux, "\n");
//
//   for (i = 0; i < nInputs; i++)
//      fprintf(fsmaux, "%lg ", px0[i]);
//   fprintf(fsmaux, "\n");
//   for (i = 0; i < nInputs; i++){
//      for(j = 0; j < nInputs; j++)
//         fprintf(fsmaux, "%lg ", B0[j + i*nInputs]);
//         fprintf(fsmaux, "\n");
//   }
//   for (i = 0; i < nInputs; i++)
//      fprintf(fsmaux, "%lg ", x0[i]);
//   fprintf(fsmaux, "\n");
//
//}
//
//   fclose(fsmaux);
//}
//

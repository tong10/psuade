// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class CorrelationAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "CorrelationAnalyzer.h"
#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
CorrelationAnalyzer::CorrelationAnalyzer() : Analyzer()
{
   setName("CORRELATION");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CorrelationAnalyzer::~CorrelationAnalyzer()
{
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double CorrelationAnalyzer::analyze(aData &adata)
{
   int     nInputs, nOutputs, nSamples, outputID, ii, jj, ss, info, idata;
   int     printLevel, count, iOne=1;
   double  xmean, xvar, *ymeans, *yvars, *Ylocal, *X, *XX, *Xvals, *Wlocal;
   double  ymean, yvar, *xmeans, *xvars, *Xlocal, *Y, *YY, *Yvals, ddata;
   FILE    *fp;
   FuncApprox *faPtr=NULL;
   PsuadeData *ioPtr=NULL;
   pData   pData;

   // extract data
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   outputID = adata.outputID_;
   printLevel = adata.printLevel_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printf("Correlation INFO: non-uniform probability distributions\n");
         printf("           have been defined in the data file, but they\n");
         printf("           will not be used in this analysis.\n");
      }
   }
   ioPtr = adata.ioPtr_;
   if (ioPtr != NULL) ioPtr->getParameter("input_names", pData);

   // error checking
   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("Correlation ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID < 0 || outputID >= nOutputs)
   {
      printf("Correlation ERROR: invalid outputID.\n");
      printf("    nOutputs = %d\n", nOutputs);
      printf("    outputID = %d\n", outputID+1);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples == 1)
   {
      printf("Correlation INFO: analysis not meaningful for nSamples=1\n");
      return PSUADE_UNDEFINED;
   } 
   info = 0; 
   for (ss = 0; ss < nSamples; ss++)
      if (Y[nOutputs*ss+outputID] == PSUADE_UNDEFINED) info = 1;
   if (info == 1)
   {
      printf("Correlation ERROR: Some outputs are undefined.\n");
      printf("                   Prune the undefined's first.\n");
      return PSUADE_UNDEFINED;
   } 
   
   // first find the mean of the current set of samples
   printAsterisks(0);
   printf("*                   Correlation Analysis\n");
   printEquals(0);
   printf("*  Basic Statistics\n");
   printDashes(0);
   printf("* Output of interest = %d\n", outputID+1);
   printDashes(0);
   computeMeanVariance(nSamples,nOutputs,Y,&ymean,&yvar,outputID,1);

   // compute the Pearson product moment correlation coefficient (PEAR)
   printEquals(0);
   printf("*  Pearson correlation coefficients (PEAR) - linear -\n");
   printf("*  which gives a measure of relationship between X_i's & Y.\n");
   printDashes(0);
   xmeans = new double[nInputs];
   xvars  = new double[nInputs];
   Xvals  = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      computeMeanVariance(nSamples, nInputs, X, &(xmeans[ii]),
                          &(xvars[ii]), ii, 0);
   computeCovariance(nSamples,nInputs,X,nOutputs,Y,xmeans,xvars,
                     ymean,yvar,outputID,Xvals);
   for (ii = 0; ii < nInputs; ii++)
      printf("* Pearson Correlation coeff.  (Input %3d) = %e\n", ii+1, 
             Xvals[ii]);

   // now write these information to a plot file
   if (psPlotTool_ == 1)
   {
      fp = fopen("scilabca.sci","w");
      if (fp == NULL)
         printOutTS(0, "CorrelationAnalysis: cannot write to scilab file.\n");
      else
         fprintf(fp,"// This file contains correlation coefficients.\n");
   }
   else
   {
      fp = fopen("matlabca.m","w");
      if (fp == NULL)
         printOutTS(0, "CorrelationAnalysis: cannot write to matlab file.\n");
      else
         fprintf(fp,"%% This file contains correlation coefficients.\n");
   }
   if (fp != NULL)
   {
      fprintf(fp, "sortFlag = 0;\n");
      fprintf(fp, "nn = %d;\n", nInputs);
      fprintf(fp, "PCC = [\n");
      for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", Xvals[ii]);
      fprintf(fp, "];\n");
      if (pData.strArray_ != NULL)
      {
         fprintf(fp, "  Str = {");
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
         fprintf(fp,"'X%d'};\n",nInputs);
      }
      else
      {
         fprintf(fp, "  Str = {");
         for (ii = 0; ii < nInputs-1; ii++)
         {
            if (pData.strArray_[ii] != NULL) fprintf(fp,"'%s',",pData.strArray_[ii]);
            else                             fprintf(fp,"'X%d',",ii+1);
         }
         if (pData.strArray_[nInputs-1] != NULL) 
              fprintf(fp,"'%s'};\n",pData.strArray_[nInputs-1]);
         else fprintf(fp,"'X%d'};\n",nInputs);
      }
      fwritePlotCLF(fp);
      fprintf(fp, "if (sortFlag == 1)\n");
      if (psPlotTool_ == 1)
           fprintf(fp, "  [PCC, II] = gsort(PCC,'g','d');\n");
      else fprintf(fp, "  [PCC, II] = sort(PCC,'descend');\n");
      fprintf(fp, "  II   = II(1:nn);\n");
      fprintf(fp, "  PCC  = PCC(1:nn);\n");
      fprintf(fp, "  Str1 = Str(II);\n");
      fprintf(fp, "else\n");
      fprintf(fp, "  Str1 = Str;\n");
      fprintf(fp, "end\n");
      fprintf(fp, "ymin = min(PCC);\n");
      fprintf(fp, "ymax = max(PCC);\n");
      fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
      fprintf(fp, "subplot(1,2,1)\n");
      fprintf(fp, "bar(PCC,0.8);\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp,"Pearson Correlation Coefficients");
      fwritePlotYLabel(fp, "Correlation Coefficient");
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "a=gca();\n");
         fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
         fprintf(fp, "a.x_ticks(2) = [1:nn]';\n");
         fprintf(fp, "a.x_ticks(3) = Str1';\n");
         fprintf(fp, "a.x_label.font_size = 3;\n");
         fprintf(fp, "a.x_label.font_style = 4;\n");
      }
      else
      {
         fprintf(fp, "axis([0 nn+1 ymin ymax])\n");
         fprintf(fp, "set(gca,'XTickLabel',[]);\n");
         fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str1,");
         fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp, "set(th, 'fontsize', 12)\n");
         fprintf(fp, "set(th, 'fontweight', 'bold')\n");
      }
   }

   // compute write these information to a plot file
   ymeans = new double[nOutputs];
   yvars  = new double[nOutputs];
   Yvals  = new double[nOutputs];
   if (nOutputs > 1)
   {
      printEquals(0);
      printf("*  PEAR (linear) for Y_i's versus Y.\n");
      printDashes(0);
      for (ii = 0; ii < nOutputs; ii++)
         computeMeanVariance(nSamples, nOutputs, Y, &(ymeans[ii]),
                             &(yvars[ii]), ii, 0);
      computeCovariance(nSamples,nOutputs,Y,nOutputs,Y,ymeans,yvars,
                        ymean,yvar,outputID,Yvals);
      for (ii = 0; ii < nOutputs; ii++)
         if (ii != outputID)
            printf("* Pearson Correlation coeff. (Output %2d) = %e\n", ii+1,
                   Yvals[ii]);
   }
   printEquals(0);


   printf("*  Spearman coefficients (SPEA) - nonlinear relationship -  *\n");
   printf("*  which gives a measure of relationship between X_i's & Y. *\n");
   printDashes(0);
   YY = new double[nSamples];
   Xlocal = new double[nSamples];
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ss = 0; ss < nSamples; ss++)
      {
         Xlocal[ss] = X[ss*nInputs+ii];
         YY[ss] = Y[ss*nOutputs+outputID];
      }
      sortDbleList2(nSamples, Xlocal, YY);
      for (ss = 0; ss < nSamples; ss++) Xlocal[ss] = (double) ss;
      sortDbleList2(nSamples, YY, Xlocal);
      for (ss = 0; ss < nSamples; ss++) YY[ss] = (double) ss;
      computeMeanVariance(nSamples,1,Xlocal,&xmean,&xvar,0,0);
      computeMeanVariance(nSamples,1,YY,&ymean,&yvar,0,0);
      computeCovariance(nSamples,1,Xlocal,1,YY,&xmean,&xvar,
                        ymean,yvar,0,&(Xvals[ii]));
      printf("* Spearman coefficient         (Input %3d ) = %e\n", ii+1, 
             Xvals[ii]);
   }
   if (fp != NULL)
   {
      fprintf(fp, "SPEA = [\n");
      for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n", Xvals[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "if (sortFlag == 1)\n");
      if (psPlotTool_ == 1)
           fprintf(fp, "  [SPEA, II] = gsort(SPEA,'g','d');\n");
      else fprintf(fp, "  [SPEA, II] = sort(SPEA,'descend');\n");
      fprintf(fp, "  II   = II(1:nn);\n");
      fprintf(fp, "  SPEA = SPEA(1:nn);\n");
      fprintf(fp, "  Str2 = Str(II);\n");
      fprintf(fp, "else\n");
      fprintf(fp, "  Str2 = Str;\n");
      fprintf(fp, "end\n");
      fprintf(fp, "ymin = min(SPEA);\n");
      fprintf(fp, "ymax = max(SPEA);\n");
      fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
      fprintf(fp, "subplot(1,2,2)\n");
      fprintf(fp, "bar(SPEA,0.8);\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp,"Spearman Correlation Coefficients");
      fwritePlotYLabel(fp, "Correlation Coefficient");
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "a=gca();\n");
         fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
         fprintf(fp, "a.x_ticks(2) = [1:nn]';\n");
         fprintf(fp, "a.x_ticks(3) = Str2';\n");
         fprintf(fp, "a.x_label.font_size = 3;\n");
         fprintf(fp, "a.x_label.font_style = 4;\n");
         printOutTS(2, " Correlation analysis plot file = scilab.sci.\n");
      }
      else
      {
         fprintf(fp, "axis([0 nn+1 ymin ymax])\n");
         fprintf(fp, "set(gca,'XTickLabel',[]);\n");
         fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str2,");
         fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp, "set(th, 'fontsize', 12)\n");
         fprintf(fp, "set(th, 'fontweight', 'bold')\n");
         printOutTS(2, " Correlation analysis plot file = matlabca.m.\n");
      }
      fclose(fp);
      fp = NULL;
   }
   printDashes(0);
   for (ii = 0; ii < nInputs; ii++) xmeans[ii] = (double) ii;
   for (ii = 0; ii < nInputs; ii++) Xvals[ii] = PABS(Xvals[ii]);
   sortDbleList2(nInputs, Xvals, xmeans);
   for (ii = nInputs-1; ii >= 0; ii--)
      printf("* Spearman coefficient(ordered) (Input %3d ) = %e\n", 
             (int) (xmeans[ii]+1), Xvals[ii]);

   Ylocal = new double[nSamples];
   if (nOutputs > 1)
   {
      printEquals(0);
      printf("*  SPEA (nonlinear) for Y_i's versus Y.                     *\n");
      printDashes(0);
      for (ii = 0; ii < nOutputs; ii++)
      {
         if (ii != outputID)
         {
            for (ss = 0; ss < nSamples; ss++)
            {
               Ylocal[ss] = Y[ss*nOutputs+ii];
               YY[ss] = Y[ss*nOutputs+outputID];
            }
            sortDbleList2(nSamples, Ylocal, YY);
            for (ss = 0; ss < nSamples; ss++) Ylocal[ss] = (double) ss;
            sortDbleList2(nSamples, YY, Ylocal);
            for (ss = 0; ss < nSamples; ss++) YY[ss] = (double) ss;
            computeMeanVariance(nSamples,1,Ylocal,&xmean,&xvar,0,0);
            computeMeanVariance(nSamples,1,YY,&ymean,&yvar,0,0);
            computeCovariance(nSamples,1,Ylocal,1,YY,&xmean,&xvar,
                              ymean,yvar,0,&(Yvals[ii]));
            printf("* Spearman coefficient        (Input %3d ) = %e\n", ii+1, 
                   Yvals[ii]);
         }
      }
   }

   if (printLevel > 1)
   {
      int nc=0, nd=0;
      printEquals(0);
      for (ii = 0; ii < nInputs; ii++)
      {
         nc = nd = 0;
         for (ss = 0; ss < nSamples; ss++)
         {
            Xlocal[ss] = X[ss*nInputs+ii];
            YY[ss] = Y[ss*nOutputs+outputID];
         }
         sortDbleList2(nSamples, YY, Xlocal);
         for (ss = 0; ss < nSamples; ss++) YY[ss] = (double) ss;
         sortDbleList2(nSamples, Xlocal, YY);
         for (ss = 0; ss < nSamples; ss++) Xlocal[ss] = (double) ss;
         for (ss = 0; ss < nSamples; ss++)
         {
            for (jj = ss+1; jj < nSamples; jj++)
            {
               ddata = Xlocal[ss] - Xlocal[jj];
               if (ddata != 0.0)
               {
                  ddata = (YY[ss] - YY[jj]) / ddata;
                  if (ddata >= 0.0) nc++;
                  else              nd++;
               }
               else
               {
                  ddata = YY[ss] - YY[jj];
                  if (ddata >= 0.0) nc++;
                  else              nd++;
               }
            }
         }
         printf("* Kendall coefficient         (Input %3d ) = %10.2e \n",
                ii+1, 2.0*(nc-nd)/(nSamples*(nSamples-1)));
         Xvals[ii] = 2.0 * (nc - nd) / (nSamples * (nSamples - 1));
      }
      printDashes(0);
      for (ii = 0; ii < nInputs; ii++) xmeans[ii] = (double) ii;
      for (ii = 0; ii < nInputs; ii++) Xvals[ii] = PABS(Xvals[ii]);
      sortDbleList2(nInputs, Xvals, xmeans);
      for (ii = nInputs-1; ii >= 0; ii--)
         printf("* Kendall coefficient(ordered) (Input %3d ) = %e\n", 
                (int) (xmeans[ii]+1), Xvals[ii]);
   }

   if (printLevel > 2)
   {
      printEquals(0);
      printf("*  Regression analysis on rank-ordered inputs/outptus       *\n");
      XX = new double[nSamples*nInputs];
      Wlocal = new double[nSamples];
      for (ss = 0; ss < nSamples*nInputs; ss++) XX[ss] = X[ss];
      for (ss = 0; ss < nSamples; ss++)
      {
         YY[ss] = Y[ss*nOutputs+outputID];
         Wlocal[ss] = (double) ss;
      }
      sortDbleList2(nSamples, YY, Wlocal);
      for (ss = 0; ss < nSamples; ss++)
      {
         idata = (int) Wlocal[ss];
         YY[idata] = (double) (ss + 1);
      }
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ss = 0; ss < nSamples; ss++)
         {
            Wlocal[ss] = (double) ss;
            Xlocal[ss] = XX[nSamples*ii+ss];
         }
         sortDbleList2(nSamples, Xlocal, Wlocal);
         for (ss = 0; ss < nSamples; ss++)
         {
            idata = (int) Wlocal[ss];
            XX[ii*nSamples+idata] = (double) (ss + 1);
         }
      }
      faPtr = genFA(PSUADE_RS_REGR1, nInputs, iOne, nSamples);
      faPtr->setOutputLevel(0);
      info = -999;
      faPtr->genNDGridData(XX, YY, &info, NULL, NULL);
      delete faPtr;
      delete [] XX;
      delete [] Wlocal;
   }
   printAsterisks(0);

   delete [] YY;
   delete [] Xlocal;
   delete [] Ylocal;
   delete [] xmeans;
   delete [] xvars;
   delete [] Xvals;
   delete [] ymeans;
   delete [] yvars;
   delete [] Yvals;
   return 0.0;
}

// *************************************************************************
// Compute mean and variance
// -------------------------------------------------------------------------
int CorrelationAnalyzer::computeMeanVariance(int nSamples, int xDim, 
              double *X, double *xmean, double *xvar, int xID, int flag)
{
   int    ss;
   double mean, variance;

   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++) mean += X[xDim*ss+xID];
   mean /= (double) nSamples;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
      variance += ((X[xDim*ss+xID] - mean) * (X[xDim*ss+xID] - mean));
   variance /= (double) (nSamples - 1);
   (*xmean) = mean;
   (*xvar)  = variance;
   if (flag == 1)
   {
      printf("Correlation: mean     = %e\n", mean);
      printf("Correlation: variance = %e\n", variance);
   }
   return 0;
}

// *************************************************************************
// Compute agglomerated covariances
// -------------------------------------------------------------------------
int CorrelationAnalyzer::computeCovariance(int nSamples,int nX,double *X,
             int nY, double *Y, double *xmeans, double *xvars, double ymean,
             double yvar, int yID, double *Rvalues)
{
   int    ii, ss;
   double denom, numer;

   for (ii = 0; ii < nX; ii++)
   {
      numer = 0.0;
      for (ss = 0; ss < nSamples; ss++)
         numer += ((X[ss*nX+ii] - xmeans[ii]) * (Y[ss*nY+yID] - ymean));
      numer /= (double) (nSamples - 1);
      denom = sqrt(xvars[ii] * yvar);
      if (denom == 0.0)
      {
         printf("Correlation ERROR: denom=0 for input %d\n", ii+1);
         printf("denom = xvar * yvar : xvar = %e, yvar = %e\n",xvars[ii],yvar);
         Rvalues[ii] = 0.0;
      }
      else Rvalues[ii] = numer / denom;
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
CorrelationAnalyzer& CorrelationAnalyzer::operator=(const CorrelationAnalyzer &)
{
   printf("Correlation operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}


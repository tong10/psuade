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
// Functions for the class PCAnalyzer (principal component analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "PCAnalyzer.h"
#include "Util/PsuadeUtil.h"
#include "Util/sysdef.h"
#include "Main/Psuade.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C" {
   void dgels_(char *, int *, int *, int *, double *A, int *LDA,
               double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}
                                                                                
// ************************************************************************
// constructor
// ------------------------------------------------------------------------
PCAnalyzer::PCAnalyzer() : Analyzer()
{
   setName("PCA");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PCAnalyzer::~PCAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double PCAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, jj, ss, ii, M, N, info;
   int    wlen, pcCnt, pcaFlag=0;
   double *Y, *YY, mean, stdev, *WW, *SS, *VV, *UU, dtemp;
   char   jobu, jobvt, winput1[500], winput2[500], pcaFile[500], *cString;
   char   pString[500];
   FILE   *fp;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   Y        = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      printf("PCAAnalyzer INFO: some inputs have non-uniform PDFs, but\n");
      printf("                  they will not be relevant in this analysis.n");
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printf("PCAnalyzer ERROR: invalid arguments.\n");
      printf("    nInputs  = %d\n", nInputs);
      printf("    nOutputs = %d\n", nOutputs);
      printf("    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nOutputs == 1)
   {
      printf("PCAnalyzer ERROR: analysis not done since nOutputs=1.\n");
      return PSUADE_UNDEFINED;
   }
   
   if (psAnaExpertMode_ == 1 && psPlotTool_ == 0)
   {
      printf("Compute contribution of principal ");
      printf("components to each output vector (y or n) ");
      sprintf(pString, "? ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the matlab file name : ");
         getString(pString, pcaFile);
         pcaFile[strlen(pcaFile)-1] = '\0';
         fp = fopen(pcaFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            pcaFlag = 1;
         }
         else printf("PCAnalyzer: cannot open matlab file %s\n", pcaFile);
      }
   }
   else
   {
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("PCA_matlab_file");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %s", winput1, winput2, pcaFile);
            fp = fopen(pcaFile, "w");
            if (fp != NULL)
            {
               fclose(fp);
               pcaFlag = 1;
            }
            else printf("PCAnalyzer ERROR: cannot open matlab file %s\n", 
                        pcaFile);
         }
      }
   }

   YY = new double[nSamples*nOutputs];
   for (jj = 0; jj < nOutputs; jj++)
      for (ss = 0; ss < nSamples; ss++)
         YY[nSamples*jj+ss] = Y[ss*nOutputs+jj];

   for (jj = 0; jj < nOutputs; jj++)
   {
      mean = 0.0;
      for (ss = 0; ss < nSamples; ss++) mean += YY[nSamples*jj+ss];
      mean /= (double) nSamples;
      for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] -= mean;
      stdev = 0.0;
      for (ss = 0; ss < nSamples; ss++)
         stdev += YY[nSamples*jj+ss] * YY[nSamples*jj+ss];
      stdev /= (double) (nSamples - 1);
      stdev = sqrt(stdev);
      if (PABS(stdev) > 1.0e-14)
         for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] /= stdev;
      printf("PCA: output %5d has mean = %e, std. dev = %e\n", jj+1, 
             mean, stdev);
   }

   M = nSamples;
   N = nOutputs;
   jobu = 'S';
   jobvt = 'A';
   wlen = 5 * M;
   SS = new double[N];
   UU = new double[M*N];
   VV = new double[N*N];
   WW = new double[wlen];
   dgesvd_(&jobu, &jobvt, &M, &N, YY, &M, SS, UU, &M, VV, &N, WW,
           &wlen, &info);
   if (info != 0)
      printf("* PCAnalyzer INFO: dgesvd returns a nonzero (%d).\n",info);
   for (ii = 0; ii < N; ii++) SS[ii] = SS[ii] * SS[ii];
   for (ii = 0; ii < N; ii++)
      printf("principal component %3d has variance = %16.8e\n",ii+1,
             SS[ii]);
   sprintf(pString,"Enter how many principal components to keep : ");
   pcCnt = getInt(1, N, pString);

   for (jj = 0; jj < nOutputs; jj++)
      for (ss = 0; ss < nSamples; ss++)
         YY[nSamples*jj+ss] = Y[ss*nOutputs+jj];
   for (jj = 0; jj < nOutputs; jj++)
   {
      mean = 0.0;
      for (ss = 0; ss < nSamples; ss++) mean += YY[nSamples*jj+ss];
      mean /= (double) nSamples;
      for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] -= mean;
      stdev = 0.0;
      for (ss = 0; ss < nSamples; ss++)
         stdev += YY[nSamples*jj+ss] * YY[nSamples*jj+ss];
      stdev /= (double) (nSamples - 1);
      stdev = sqrt(stdev);
      if (PABS(stdev) > 1.0e-14)
         for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] /= stdev;
   }


   if (pcaFlag == 1)
   {
      fp = fopen(pcaFile, "w");
      if (fp == NULL)
      {
         fprintf(fp, "%% Contribution of principal component to each output\n");
         fprintf(fp, "%% First display the Scree diagram\n");
         fprintf(fp, "%% A - original sample outputs (normalized)\n");
         fprintf(fp, "%% V - eigenvector matrix\n");
         fprintf(fp, "%% U - left singular vector matrix\n");
         fprintf(fp, "%% S - singular values (column vector)\n");

         fprintf(fp, "A = [\n");
         for (ss = 0; ss < M; ss++) 
         {
            for (jj = 0; jj < N; jj++)
               fprintf(fp, "%e ", YY[M*jj+ss]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "];\n");
         fprintf(fp, "S = [\n");
         for (ii = 0; ii < N; ii++) fprintf(fp, " %e\n", SS[ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "hold off\n");
         fprintf(fp, "plot(S)\n");
         fprintf(fp, "grid\n");
         fprintf(fp, "xlabel('Principal component number')\n");
         fprintf(fp, "ylabel('Variances')\n");
         fprintf(fp, "title('Scree Diagram')\n");
         fprintf(fp, "disp('Press enter to continue')\n");
         fprintf(fp, "pause\n");

         fprintf(fp, "V = [\n");
         for (ii = 0; ii < N; ii++)
         {
            for (jj = 0; jj < N; jj++)
               fprintf(fp, " %e", VV[jj*N+ii]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "];\n");
         fprintf(fp, "U = [\n");
         for (ii = 0; ii < nSamples; ii++)
         {
            for (jj = 0; jj < pcCnt; jj++)
               fprintf(fp, " %e", UU[jj*M+ii]);
            fprintf(fp, "\n");
         }
         fprintf(fp, "];\n");
         fprintf(fp, "YA = diag(S(1:%d)) * U';\n",pcCnt);
         fprintf(fp, "XA = 1:%d;\n",pcCnt);
         fprintf(fp, "hold off\n");
         fprintf(fp, "disp('Each output projected on the PC space')\n");
         fprintf(fp, "for ii = 1 : %d\n", N);
         fprintf(fp, "   plot(XA,YA(:,ii),XA,YA(:,ii),'x')\n");
         fprintf(fp, "   disp(['output = ' int2str(ii) ' projected on PCs'])\n");
         fprintf(fp, "   disp('Press enter to continue')\n");
         fprintf(fp, "   pause\n");
         fprintf(fp, "   if ii == 1\n");
         fprintf(fp, "      hold on\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "hold off\n");
         fprintf(fp, "disp('Next display each principal component')\n");
         fprintf(fp, "disp('Press enter to continue')\n");
         fprintf(fp, "pause\n");
         fprintf(fp, "XV = 1:%d;\n",N);
         fprintf(fp, "for ii = 1 : %d\n", pcCnt);
         fprintf(fp, "   plot(XV,V(:,ii),XV,V(:,ii),'x')\n");
         fprintf(fp, "   disp(['principal component = ' int2str(ii)])\n");
         fprintf(fp, "   disp('Press enter to continue')\n");
         fprintf(fp, "   pause\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "hold off\n");
         fprintf(fp, "disp('Now all principal components')\n");
         fprintf(fp, "plot(XV,V,XV,V,'x')\n");
         fprintf(fp, "%%barh(abs(A)','stacked')\n");
         fprintf(fp, "%%area(XA,abs(A))\n");
         fprintf(fp, "hold off\n");
         fclose(fp);
      }
      else printf("WARNING: cannot open file for plotting.\n");
   }

   for (ss = 0; ss < nSamples; ss++)
   {
      for (jj = 0; jj < pcCnt; jj++)
      {
         dtemp = 0.0;
         for (ii = 0; ii < N; ii++)
            dtemp += VV[jj*N+ii] * YY[ii*nSamples+ss];
         Y[pcCnt*ss+jj] = dtemp;
      }
   }
   fp = fopen("psPCA.out", "w");
   if (fp == NULL)
   {
      printf("PCAnalyzer ERROR: failed to store principal component.\n");
   }
   else
   {
      fprintf(fp, "%d %d\n", nSamples, pcCnt);
      for (ss = 0; ss < nSamples; ss++)
      {
         for (jj = 0; jj < pcCnt; jj++)
            fprintf(fp, "%24.16e ", Y[pcCnt*ss+jj]);
         fprintf(fp, "\n");
      }
      fclose(fp);
      printf("PCAnalyzer: principal components stored in psPCA.out file.\n");
   }

   adata.nOutputs_ = pcCnt;
   delete [] SS;
   delete [] WW;
   delete [] VV;
   delete [] UU;
   delete [] YY;
   printf("PCAnalysis completed: you can use write to update file.\n");
   return 0.0;
}


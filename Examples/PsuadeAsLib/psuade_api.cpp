#include <memory>
#include <iostream>
#include <fstream>
#include "Sampling.h"
#include "LPtauSampling.h"
#include "MOATSampling.h"
#include "MomentAnalyzer.h"
#include "MOATAnalyzer.h"
#include "RSMSobol1Analyzer.h"
#include "RSMSobol2Analyzer.h"
#include "RSMSobolTSIAnalyzer.h"
#include "PsuadeData.h"
#include "aData.h"
#include "sysdef.h"

using namespace std;

void runSimulation (double *inputs, double *results);

int main(int argc, char *argv[])
{
   // test problem has 3 inputs and 1 output and some bounds
   int    printLevel=1, rstype = PSUADE_RS_GP3;
   int    nInputs = 3, nOutputs = 1;
   double *lowerB = new double[3];
   double *upperB = new double[3];
   lowerB[0] = 40.0; // parameter 1
   upperB[0] = 60.0;
   lowerB[1] = 67.0; // parameter 2
   upperB[1] = 74.0;
   lowerB[2] = 20.0; // parameter 3
   upperB[2] = 40.0;

   // generate quasi-Monte Carlo sample
   int nSamples = 200;
   printf("Creating a LPTAU sample of size %d ...\n", nSamples);
   Sampling *sampler = (Sampling *) new LPtauSampling();
   sampler->setPrintLevel(printLevel);
   sampler->setInputBounds(nInputs, lowerB, upperB);
   sampler->setOutputParams(nOutputs);
   sampler->setSamplingParams(nSamples, 1, 0);
   sampler->initialize(0);
   nSamples = sampler->getNumSamples();
   double *sampleInputs  = new double[nSamples * nInputs];
   double *sampleOutputs = new double[nSamples * nOutputs];
   int    *sampleStates  = new int[nSamples];
   sampleOutputs = new double[nSamples * nOutputs];
   sampleStates  = new int[nSamples];
   sampler->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                       sampleOutputs, sampleStates);

   // run samples
   printf("Running %d simulations ...\n", nSamples);
   for (int ii = 0; ii < nSamples; ii++)
   {
      runSimulation(&sampleInputs[ii*nInputs],&sampleOutputs[ii*nOutputs]);
      printf("Simulation %4d : Inputs = ", ii+1);
      printf("%4d ", ii+1);
      for (int jj = 0; jj < nInputs; jj++)
         printf("%12.4e ", sampleInputs[ii*nInputs+jj]); 
      printf(", Outputs = ");
      for (int jj = 0; jj < nOutputs; jj++)
         printf("%12.4e ", sampleOutputs[ii*nOutputs+jj]); 
      printf("\n");
      sampleStates[ii] = 1;
   }

   // Global sensitivity analysis: first, second, and total order 
   RSMSobol1Analyzer   *analyzer1 = new RSMSobol1Analyzer();
   RSMSobol2Analyzer   *analyzer2 = new RSMSobol2Analyzer();
   RSMSobolTSIAnalyzer *analyzer3 = new RSMSobolTSIAnalyzer();
   int *pdfFlags = new int[nInputs];
   for (int kk = 0; kk < nInputs; kk++) pdfFlags[kk] = 0;
   double *means  = new double[nInputs];
   double *stdevs = new double[nInputs];
   for (int kk = 0; kk < nInputs; kk++) means[kk] = stdevs[kk] = 0;
   aData adata;
   adata.printLevel_ = printLevel;
   adata.nInputs_  = nInputs;
   adata.nOutputs_ = nOutputs;
   adata.nSamples_ = nSamples;
   adata.iLowerB_  = lowerB;
   adata.iUpperB_  = upperB;
   adata.sampleInputs_  = sampleInputs;
   adata.sampleOutputs_ = sampleOutputs;
   adata.sampleStates_  = sampleStates;
   adata.outputID_      = 0;
   adata.inputPDFs_     = pdfFlags;
   adata.inputMeans_    = means;
   adata.inputStdevs_   = stdevs;
   PsuadeData *pIO = new PsuadeData();
   pIO->updateInputSection(nSamples,nInputs,NULL,lowerB,upperB,
              sampleInputs, NULL, NULL,NULL,NULL,NULL);
   pIO->updateOutputSection(nSamples,1,sampleOutputs,sampleStates,
              NULL);
   pIO->updateMethodSection(PSUADE_SAMP_MC,nSamples,-1,-1,-1);
   pIO->updateAnalysisSection(-1,-1,rstype,-1,-1,-1);
   adata.ioPtr_ = pIO;

   printf("Running first-order Sobol Analysis ...\n");
   analyzer1->analyze(adata);
   for (int ii = 0; ii < nInputs; ii++)
     printf("Sobol first-order index for input %d = %10.2e\n",
            ii+1,analyzer1->get_vce(ii));
   printf("Running second-order Sobol Analysis ...\n");
   analyzer2->analyze(adata);
   for (int ii = 0; ii < nInputs; ii++)
     for (int jj = ii+1; jj < nInputs; jj++)
       printf("Sobol second-order index for inputs %d %d = %10.2e\n",
              ii+1,jj+1,analyzer2->get_vce(ii,jj));
   printf("Running total-order Sobol Analysis ...\n");
   analyzer3->analyze(adata);
   for (int ii = 0; ii < nInputs; ii++)
     printf("Sobol total-order index for input %d = %10.2e\n",
            ii+1,analyzer3->get_tsi(ii));
   printf("Running uncertainty Analysis ...\n");
   MomentAnalyzer *analyzer4 = new MomentAnalyzer();
   analyzer4->analyze(adata);

   // clean up
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   delete sampler;
   delete analyzer1;
   delete analyzer2;
   delete analyzer3;
   delete analyzer4;
   delete [] pdfFlags;
   delete [] means;
   delete [] stdevs;

   // MOAT analysis
   sampler = (Sampling *) new MOATSampling();
   nSamples = 40;
   printf("Creating a MOAT sample of size %d ...\n", nSamples);
   sampler->setPrintLevel(printLevel);
   sampler->setInputBounds(nInputs, lowerB, upperB);
   sampler->setOutputParams(nOutputs);
   sampler->setSamplingParams(nSamples, 1, 0);
   sampler->initialize(0);
   nSamples = sampler->getNumSamples();
   sampleInputs  = new double[nSamples * nInputs];
   sampleOutputs = new double[nSamples * nOutputs];
   sampleStates  = new int[nSamples];
   sampleOutputs = new double[nSamples * nOutputs];
   sampleStates  = new int[nSamples];
   sampler->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                       sampleOutputs, sampleStates);
   printf("Running %d simulations ...\n", nSamples);
   for (int ii = 0; ii < nSamples; ii++)
   {
      runSimulation(&sampleInputs[ii*nInputs],&sampleOutputs[ii*nOutputs]);
      printf("Simulation %4d : Inputs = ", ii+1);
      for (int jj = 0; jj < nInputs; jj++)
         printf("%12.4e ", sampleInputs[ii*nInputs+jj]); 
      printf(", Outputs = ");
      for (int jj = 0; jj < nOutputs; jj++)
         printf("%12.4e ", sampleOutputs[ii*nOutputs+jj]); 
      printf("\n");
      sampleStates[ii] = 1;
   }
   printf("Running MOAT Analysis ...\n");
   MOATAnalyzer *analyzer5 = new MOATAnalyzer();
   adata.nSamples_ = nSamples;
   adata.sampleInputs_  = sampleInputs;
   adata.sampleOutputs_ = sampleOutputs;
   adata.sampleOutputs_ = sampleOutputs;
   adata.sampleStates_  = sampleStates;
   adata.inputPDFs_     = NULL;
   adata.inputMeans_    = NULL;
   adata.inputStdevs_   = NULL;
   analyzer5->analyze(adata);
   for (int ii = 0; ii < nInputs; ii++)
     printf("MOAT modified mean/std for input %d = %10.2e %10.2e\n",
            ii+1,analyzer5->get_modifiedMean(ii),
            analyzer5->get_modifiedStdev(ii));

   // clean up
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   delete sampler;
   delete analyzer5;
   delete pIO;

   delete [] lowerB;
   delete [] upperB;
   return 0;
}

//**/ user example: Bungee
void runSimulation (double *inputs, double *results)
{
   double Y, kel=1.5, g=9.8;
   Y = inputs[0] - 2.0 * inputs[1] * g / (kel * inputs[2]);
   (*results) = Y;
}


#include <memory>
#include <iostream>
#include <fstream>
#include "aData.h"
#include "Sampling.h"
#include "LPtauSampling.h"
#include "MOATSampling.h"
#include "MomentAnalyzer.h"
#include "MOATAnalyzer.h"
#include "RSMSobol1Analyzer.h"
#include "RSMSobol2Analyzer.h"
#include "RSMSobolTSIAnalyzer.h"
#include "Globals.h"

using namespace std;

void runSimulation (double *inputs, double *results);

int main(int argc, char *argv[])
{
   // test problem has 3 inputs and 1 output and some bounds
   int    printLevel=0, targc, ii;
   int    nInputs = 3, nOutputs = 1;
   double *lowerB = new double[3];
   double *upperB = new double[3];
   char   **targv;

   lowerB[0] = 40.0; // parameter 1
   upperB[0] = 60.0;
   lowerB[1] = 67.0; // parameter 2
   upperB[1] = 74.0;
   lowerB[2] = 20.0; // parameter 3
   upperB[2] = 40.0;

   // generate quasi-Monte Carlo sample
   int nSamples = 300;
   Sampling *sampler = (Sampling *) new LPtauSampling();
   sampler->setPrintLevel(printLevel);
   sampler->setInputBounds(nInputs, lowerB, upperB);
   sampler->setOutputParams(nOutputs);
   sampler->setSamplingParams(nSamples, 1, 0);
   sampler->initialize(0);
   nSamples = sampler->getNumSamples();
   double *sampleInputs  = new double[nSamples * nInputs];
   double *sampleOutputs = new double[nSamples * nOutputs];
   sampler->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                       sampleOutputs);
   delete sampler;

   // run samples
   printf("Running samples: \n");
   for (int ii = 0; ii < nSamples; ii++)
   {
      runSimulation(&sampleInputs[ii*nInputs],
                    &sampleOutputs[ii*nOutputs]);
      //printf("%4d ", ii+1);
      //for (int jj = 0; jj < nInputs; jj++)
      //   printf("%12.4e ", sampleInputs[ii*nInputs+jj]); 
      //for (int jj = 0; jj < nOutputs; jj++)
      //   printf("%12.4e ", sampleOutputs[ii*nOutputs+jj]); 
      //printf("\n");
   }

   // Global sensitivity analysis: first order
   printf("====> Begin analysis: RSMSobol1\n");
   RSMSobol1Analyzer *analyzer1 = new RSMSobol1Analyzer();
   targc = 2;
   targv = new char*[2];;
   targv[0] = new char[100];
   strcpy(targv[0], "rstype");
   targv[1] = new char[100];
   strcpy(targv[1], "mars");
   analyzer1->setParams(targc, targv);
   psConfig_.InteractiveOn();
   analyzer1->analyze(nInputs,nSamples,lowerB,upperB,sampleInputs,
                      sampleOutputs);
   for (ii = 0; ii < nInputs; ii++)
      printf("RSM VCE %2d = %e\n", ii+1, analyzer1->get_vce(ii));

   // Global sensitivity analysis: second order 
   printf("====> Begin analysis: RSMSobol2\n");
   RSMSobol2Analyzer *analyzer2 = new RSMSobol2Analyzer();
   analyzer2->setParams(targc, targv);
   analyzer2->analyze(nInputs,nSamples,lowerB,upperB,sampleInputs,
                      sampleOutputs);
   for (int ii = 0; ii < nInputs; ii++)
      for (int ii2 = ii+1; ii2 < nInputs; ii2++)
         printf("RSM2 VCE %2d %2d = %e\n", ii+1, ii2+1,
                analyzer2->get_vce(ii,ii2));

   // Global sensitivity analysis: total order 
   printf("====> Begin analysis: RSMSobolTSI\n");
   RSMSobolTSIAnalyzer *analyzer3 = new RSMSobolTSIAnalyzer();
   analyzer3->setParams(targc, targv);
   analyzer3->analyze(nInputs,nSamples,lowerB,upperB,sampleInputs,
                      sampleOutputs);
   for (ii = 0; ii < nInputs; ii++)
      printf("RSM TSI %2d = %e\n", ii+1, analyzer3->get_tsi(ii));

   // Moment analysis
   printf("====> Begin analysis: Moment\n");
   MomentAnalyzer *analyzer4 = new MomentAnalyzer();
   analyzer4->analyze(nInputs,nSamples,lowerB,upperB,sampleInputs,
                      sampleOutputs);
   printf("Sample mean     = %e\n", analyzer4->get_mean());
   printf("Sample stdev    = %e\n", analyzer4->get_stdev());
   printf("Sample skewness = %e\n", analyzer4->get_skewness());
   printf("Sample kurtosis = %e\n", analyzer4->get_kurtosis());

   // clean up
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete analyzer1;
   delete analyzer2;
   delete analyzer3;
   delete analyzer4;

   // MOAT analysis
   sampler = (Sampling *) new MOATSampling();
   nSamples = 40;
   sampler->setPrintLevel(printLevel);
   sampler->setInputBounds(nInputs, lowerB, upperB);
   sampler->setOutputParams(nOutputs);
   sampler->setSamplingParams(nSamples, 1, 0);
   sampler->initialize(0);
   nSamples = sampler->getNumSamples();
   sampleInputs  = new double[nSamples * nInputs];
   sampleOutputs = new double[nSamples * nOutputs];
   sampler->getSamples(nSamples, nInputs, nOutputs, sampleInputs,
                       sampleOutputs);
   delete sampler;
   printf("Running simulations for MOAT:\n");
   for (int ii = 0; ii < nSamples; ii++)
   {
      runSimulation(&sampleInputs[ii*nInputs],
                    &sampleOutputs[ii*nOutputs]);
      //printf("%4d ", ii+1);
      //for (int jj = 0; jj < nInputs; jj++)
      //   printf("%12.4e ", sampleInputs[ii*nInputs+jj]); 
      //for (int jj = 0; jj < nOutputs; jj++)
      //   printf("%12.4e ", sampleOutputs[ii*nOutputs+jj]); 
      //printf("\n");
   }
   printf("====> Begin analysis: MOAT\n");
   MOATAnalyzer *analyzer5 = new MOATAnalyzer();
   analyzer5->analyze(nInputs,nSamples,lowerB,upperB,sampleInputs,
                      sampleOutputs);
   for (ii = 0; ii < nInputs; ii++)
      printf("MOAT %2d = %e\n", ii+1, analyzer5->get_modifiedMean(ii));

   // clean up
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete analyzer5;
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


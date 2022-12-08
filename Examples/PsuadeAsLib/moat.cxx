#include <sysdef.h>  //Where the ENUMS are hidden in PSUADE
#include <Sampling.h>
#include <MOATSampling.h>
#include <PsuadeData.h>
#include <FunctionInterface.h>
#include <stdio.h>
#include <PrintingTS.h>
#include <AnalysisManager.h>
#include <Globals.h>
#include <PsuadeConfig.h>
#include <fstream>

#define CPPMOATOUT (char *)"C++MOAT.out"

int main() 
{

  PsuadeData psuadeIO;

  //this has no effect on analysis output
  setPrintLevelTS(3); //Full Interactive Output!
  
  psConfig_.setRandomSeed(1);
  psConfig_.InteractiveOn();

  int numSamples = 180;

  //// Input section
  { 
    int nInputs = 8;
    double lowerB[nInputs];
    double upperB[nInputs];
    char*  names[nInputs];

    for(int ii = 0; ii < nInputs; ++ii) {
      lowerB[ii] = 0;
      upperB[ii] = 1.0;
      char* test = new char[10];
      snprintf(test, 10, "X%d", ii+1);
      names[ii] = test;
    }
    //symbols can be null
    psuadeIO.createInputSection(nInputs, NULL, lowerB, upperB, names);

    //Update also calls create, but it turns out to be unecessary for 
    //what we're doing.
    //psuadeIO.updateInputSection(numSamples, nInputs, NULL, lowerB, 
    //upperB, NULL, names, NULL, NULL, NULL, NULL);
  }  
  //End input section

  //// Output section
  { 
    int nOutputs = 1;
    char*  names[nOutputs];

    for(int ii = 0; ii < nOutputs; ++ii) {
      char* test = new char[10];
      snprintf(test, 10, "Y%d", ii+1);
      names[ii] = test;
    }
    
    psuadeIO.createOutputSection(nOutputs, names);

  }
  //End output section

  //Method Section
  {
    //MOAT needs 1 replication and no random
    psuadeIO.updateMethodSection(PSUADE_SAMP_MOAT, numSamples, 0, 0, 0);
  }
  //End method section
  Sampling* sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_MOAT);

  //  MOATSampling moat;
  sampler->doSampling(&psuadeIO);
  
  //Application run section
  {
    psuadeIO.updateApplicationSection(const_cast<char*>("simulator"), 
                                      const_cast<char*>(""), 
                                      const_cast<char*>(""), 
                                      const_cast<char*>(""), 
                                      const_cast<char*>(""), 1);
  }


  psuadeIO.writePsuadeFile("MOATSample.in", 0);

  //Run section  Obviously this could be done by my effectively by 
  //something other than PSUADE
  {
    //The Function interface runs jobs... kinda  We should improve this
    FunctionInterface* runner = createFunctionInterface(&psuadeIO);
    
    int numSamples = psuadeIO.getMethod_nSamples();
    int numInputs = psuadeIO.getInput_nInputs();
    double* sampleInputs = psuadeIO.getInput_sampleInputs();
    int numOutputs = psuadeIO.getOutput_nOutputs();
    double* sampleOutputs = psuadeIO.getInput_sampleInputs();
    int* sampleStates = psuadeIO.getOutput_sampleStates();

    for (int sampleID = 0; sampleID < numSamples; sampleID++) {
      //launches the run (blocking in this case)
      int status = runner->evaluate(sampleID, numInputs,
                                    &sampleInputs[sampleID*numInputs], 
                                    numOutputs, 
                                    &sampleOutputs[sampleID*numOutputs], 0);   
      if(status) {
        printf("ERROR, sample %d failed with %d\n", sampleID+1, status);
      }

      sampleStates[sampleID] = 1;

    }
    
    psuadeIO.updateOutputSection(numSamples, numOutputs, sampleOutputs, 
                sampleStates, psuadeIO.getOutput_outputNames());
    
  }

  //the way analysis output is controlled is thru an outputLevel
  int outputLevel = 0;
  psuadeIO.updateAnalysisSection(PSUADE_ANA_MOAT, 0,0,outputLevel,0,0);

  psuadeIO.writePsuadeFile("MOATPostRun.in", 0);
  

  //Analysis

  {
    AnalysisManager analysisManager;
    analysisManager.setup(&psuadeIO);
    analysisManager.analyze(&psuadeIO, 0, NULL, 0);

    //check MOAT analysis results
    //MOATAnalyzer  *MOAT_analyzer = analysisManager.getMOATAnalyzer();
    //MOAT_analyzer->createOutputString();
    //MOAT_analyzer->saveOutputString(CPPMOATOUT);
    //if (MOAT_analyzer->checkOutputString(CPPMOATOUT))
    //	std::cout << "MOAT Anaylsis agrees with " << CPPMOATOUT << endl;
  }
  return 0;
}


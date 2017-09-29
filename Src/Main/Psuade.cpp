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
// This is the main program for experimenting with different samplings
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

// -------------------------------------------------------------------------
// internal and external library declarations
// -------------------------------------------------------------------------
#ifdef HAVE_MPICH
#include <mpi.h>
#endif
#include <math.h>
#include <stdio.h>
#include "PsuadeBase.h"
#include "Comm/CommManager.h"
#include "DataIO/PsuadeConfig.h"
#include "DataIO/PsuadeData.h"
#include "dtype.h"
#include "Util/Exceptions.h"
#include "Util/PsuadeUtil.h"
#include "Base/Globals.h"

// should be changed if fixed
#ifdef SUN4
  extern "C" {
    void MAIN_(int argc, void** argv) { }
  }
#endif
extern "C" {
int RunParallel(char *);
int RunParallelLocal(char *);
}

#ifdef  HAVE_MPICH
#ifndef HAVE_MPI_USR_FCN
extern "C" {
int UserFunction(MPI_Comm mpiComm, int index, char *workdir)
{
   (void) mpiComm;
   (void) index;
   (void) workdir;
   return -1;
}
}
#endif
#endif

// ************************************************************************
// Main driver
// ------------------------------------------------------------------------
main(int argc, char** argv)
{
   int        ind, status=0, mypid=0, nprocs=1, continueFlag=1, mode=0;
   int        one=1, root=0;
   char       *inputString, psuadeFile[200], *inputFileName;
   PsuadeBase *psuade;
   FILE       *fp;

   psCommMgr_ = new CommManager(argc, (void **) argv);
   mypid      = psCommMgr_->getPID();
   nprocs     = psCommMgr_->getNumProcs();

   for (ind = 0; ind < argc; ind++)
   {
      inputString = (char *) argv[ind]; 
      if (!strcmp(inputString, "-help"))
      {
         printf("How to run psuade:\n");
         printf("    Normal      mode: psuade <file>\n");
         printf("    Interactive mode: psuade\n");
         printf("    MP          mode: mpirun/psub/srun ... -mp \n");
         printf("    MPP         mode: psub/srun psfile\n");
         continueFlag = 0;
         break;
      }
   }

   if (nprocs > 1 && argc >= 2)
   {
      if (mypid == 0)
      {
         for (ind = 0; ind < argc; ind++)
         {
            inputString = (char *) argv[ind]; 
            if (!strcmp(inputString, "-mp"))
            {
               break;
            }
            else if (!strcmp(inputString, "-mpp"))
            {
               break;
            }
         }
      }
      if (mypid == 0 && mode == 1 && argc > (ind+1))
      {
         inputFileName = (char *) argv[ind+1];
         mode = 2;
      }
      if (mypid == 0 && mode == 3 && argc > (ind+1))
      {
         inputFileName = (char *) argv[ind+1];
         mode = 4;
      }
      psCommMgr_->bcast((void *) &mode, one, INT, root);
      if (mode == 2)
      {
         RunParallel(inputFileName);
         continueFlag = 0;
      }
      else if (mode == 4)
      {
         RunParallelLocal(inputFileName);
         continueFlag = 0;
      }
      else if (mode == 1)
      {
         if (mypid == 0)
            printf("PSUADE syntax: %s -mp <psuadeInFile>\n", (char *) argv[0]);
         continueFlag = 0;
      }
   }

   if (continueFlag == 1 && mypid == 0)
   {
      strcpy(psuadeFile, "NOFILE");
      printAsterisks(0);
      printf("*      Welcome to PSUADE (version 1.4e)\n");
      printAsterisks(0);
      fflush(stdout);

      if      (argc == 1) inputFileName = NULL;
      else if (argc == 2) inputFileName = (char *) argv[1]; 
      else
      {
         printf("PSUADE syntax : %s <psuadeInFile>\n", (char *) argv[0]);
         exit(1);
      }
      if (inputFileName != NULL) strcpy(psuadeFile, inputFileName);


      psuade = new PsuadeBase();
      if (strcmp(psuadeFile, "NOFILE"))
      {
         status = psuade->getInputFromFile(psuadeFile);
         if (status != 0) exit(1);
      }

      try 
      {
         if (inputFileName == NULL) psuade->analyzeInteractive();
         else                       psuade->run();
      }
      catch ( Psuade_Stop_Exception ) 
      {
         exit(1);
      }

      delete psuade;
   }

   if (psConfigFileName_ != NULL) delete psConfigFileName_;
   if (psConfig_  != NULL) delete psConfig_;
   if (psCommMgr_ != NULL) delete psCommMgr_;
}

// ************************************************************************
// run jobs in MPI parallel mode
// ------------------------------------------------------------------------
int RunParallel(char *inFileName)
{
   int  mypid=0, nprocs=1, root=0, one=1, sampleID, lineLeng=200;
   int  nSamples=0, strLeng, *statusArray=NULL, nActive, index;
   char *keywords[] = {"workdir", "executable", "num_samples",
                       "PSUADE_PARALLEL"}; 
   char lineIn[501], inString[500], inString2[500];
   char workdir[500], executable[500], runLine[500];
   FILE *inFile;
                                                                                
   mypid  = psCommMgr_->getPID();
   nprocs = psCommMgr_->getNumProcs();

   if (mypid == 0)
   {
      inFile = fopen(inFileName, "r");
      if (inFile == NULL)
      {
         printf("RunParallel read ERROR - file %s = NULL\n",
                inFileName);
         exit(1);
      }
      fgets(lineIn, lineLeng, inFile);
      sscanf(lineIn, "%s", inString);
      if (strcmp(inString, "PSUADE_PARALLEL")) 
      {
         printf("RunParallel ERROR - keyword missing.\n");
         exit(1);
      }
      while ((fgets(lineIn, lineLeng, inFile) != NULL) && (feof(inFile) == 0))
      {
         strcpy(inString, "#");
         sscanf(lineIn,"%s", inString);
         if (strcmp(inString, keywords[0]) == 0) // work directory 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %s", inString, inString2, workdir);
            else sscanf(lineIn, "%s %s", inString, workdir);
         }
         else if (strcmp(inString, keywords[1]) == 0) // executable 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %s", inString, inString2, executable);
            else sscanf(lineIn, "%s %s", inString, executable);
         }
         else if (strcmp(inString, keywords[2]) == 0) // nSamples 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %d", inString, inString2, &nSamples);
            else sscanf(lineIn, "%s %d", inString, &nSamples);
            if (nSamples <= 0)
            {
               printf("RunParallel ERROR - nSamples %d <= 0.\n",
                      nSamples);
               exit(1);
            }
         }
         else if (strcmp(inString, keywords[3]) == 0) // end 
         {
            break;
         }
         else if (inString[0] == '#') /* comments */
         {
            /* comment line */
         }
         else
         {
            printf("RunParallel ERROR - unrecognized line - %s\n",lineIn);
            exit(1);
         }
      }
      fclose(inFile);
   }

   if (mypid == 0)
   {
      statusArray = new int[nSamples];
      for (sampleID = 0; sampleID < nSamples; sampleID++)
         statusArray[sampleID] = -1;
      nActive = 0;
      for (sampleID = 0; sampleID < nSamples; sampleID++)
      {
         sprintf(runLine, "%s.%d", workdir, sampleID+1);
         inFile = fopen(runLine, "r");
         if (inFile == NULL)
         {
            printf("RunParallel ERROR: %s does not exist.\n",runLine);
            exit(1);
         }
         fclose(inFile);
         sprintf(runLine, "%s.%d/completed", workdir, sampleID+1);
         inFile = fopen(runLine, "r");
         if (inFile == NULL) statusArray[nActive++] = sampleID;
         if (inFile != NULL) fclose(inFile);
      }
   }

   if (mypid == 0) strLeng = strlen(workdir) + 1;
   else            strLeng = 10;
   workdir[strLeng] = '\0';
   psCommMgr_->bcast((void *) &strLeng, one, INT, root);
   psCommMgr_->bcast((void *) workdir, strLeng, CHAR, root);
   if (mypid == 0) strLeng = strlen(executable) + 1;
   else            strLeng = 10;
   workdir[strLeng] = '\0';
   psCommMgr_->bcast((void *) &strLeng, one, INT, root);
   psCommMgr_->bcast((void *) executable, strLeng, CHAR, root);
   psCommMgr_->bcast((void *) &nActive, one, INT, root);
   if (nActive > 0 && mypid != 0) statusArray = new int[nActive]; 
   psCommMgr_->bcast((void *) statusArray, nActive, INT, root);

   for (sampleID = 0; sampleID < nActive; sampleID++)
   {
      index = statusArray[sampleID];
      if (sampleID % nprocs == mypid)
      {
         printf("Processor %d running job %d\n", mypid, index+1);
         sprintf(runLine, "(cd %s.%d; %s)", workdir, index+1, executable);
         system(runLine);
      }
   }
   if (statusArray != NULL) delete [] statusArray;
   return 0;
}

// ************************************************************************
// run jobs in MPI parallel mode
// ------------------------------------------------------------------------
int RunParallelLocal(char *inFileName)
{
   int  mypid=0, nprocs=1, root=0, one=1, sampleID, ip, ii, lineLeng=200;
   int  nSamples=0, strLeng, *statusArray, nActive, index, npPerJob=0;
   char *keywords[] = {"workdir", "num_samples", "nprocs_per_job", 
                       "PSAUDE_PARALLEL"}; 
   char lineIn[501], inString[500], inString2[500], outString[500];
   char workdir[500];
   FILE *inFile;
                                                                                
   mypid  = psCommMgr_->getPID();
   nprocs = psCommMgr_->getNumProcs();

   if (mypid == 0)
   {
      inFile = fopen(inFileName, "r");
      if (inFile == NULL)
      {
         printf("RunParallelLocal read ERROR - file %s = NULL\n",
                inFileName);
         exit(1);
      }
      fgets(lineIn, lineLeng, inFile);
      sscanf(lineIn, "%s", inString);
      if (strcmp(inString, "PSUADE_PARALLEL")) 
      {
         printf("RunParallelLocal ERROR - keyword missing.\n");
         exit(1);
      }
      while ((fgets(lineIn, lineLeng, inFile) != NULL) && (feof(inFile) == 0))
      {
         strcpy(inString, "#");
         sscanf(lineIn,"%s", inString);
         if (strcmp(inString, keywords[0]) == 0) // work directory 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %s", inString, inString2, workdir);
            else sscanf(lineIn, "%s %s", inString, workdir);
         }
         else if (strcmp(inString, keywords[1]) == 0) // nSamples 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %d", inString, inString2, &nSamples);
            else sscanf(lineIn, "%s %d", inString, &nSamples);
            if (nSamples <= 0)
            {
               printf("RunParallelLocal ERROR - nSamples %d <= 0.\n",
                      nSamples);
               exit(1);
            }
         }
         else if (strcmp(inString, keywords[2]) == 0) // number of procs 
         {
            sscanf(lineIn, "%s %s", inString, inString2);
            if (strcmp(inString2, "=") == 0 )
               sscanf(lineIn, "%s %s %d", inString, inString2, &npPerJob);
            else sscanf(lineIn, "%s %d", inString, &npPerJob);
            if (npPerJob <= 0)
            {
               printf("RunParallelLocal ERROR - npPerJob %d <= 0.\n",
                      npPerJob);
               exit(1);
            }
         }
         else if (strcmp(inString, keywords[3]) == 0) // end 
         {
            break;
         }
         else if (inString[0] == '#') /* comments */
         {
            /* comment line */
         }
         else
         {
            printf("RunParallelLocal ERROR - unrecognized line - %s\n",lineIn);
            exit(1);
         }
      }
      fclose(inFile);
   }
   if (nSamples <= 0 || npPerJob <= 0)
   {
      printf("RunParallelLocal ERROR: nSamples or npPerJob not specified.\n");
      exit(1);
   }

   if (mypid == 0)
   {
      statusArray = new int[nSamples];
      for (sampleID = 0; sampleID < nSamples; sampleID++)
         statusArray[sampleID] = -1;
      nActive = 0;
      for (sampleID = 0; sampleID < nSamples; sampleID++)
      {
         sprintf(outString, "%s.%d", workdir, sampleID+1);
         inFile = fopen(outString, "r");
         if (inFile == NULL)
         {
            printf("RunParallelLocal ERROR: %s does not exist.\n",outString);
            exit(1);
         }
         fclose(inFile);
         sprintf(outString, "%s.%d/completed", workdir, sampleID+1);
         inFile = fopen(outString, "r");
         if (inFile == NULL) statusArray[nActive++] = sampleID;
         if (inFile != NULL) fclose(inFile);
      }
   }

   if (mypid == 0) strLeng = strlen(workdir) + 1;
   else            strLeng = 10;
   workdir[strLeng] = '\0';
   psCommMgr_->bcast((void *) &strLeng, one, INT, root);
   psCommMgr_->bcast((void *) workdir, strLeng, CHAR, root);
   psCommMgr_->bcast((void *) &nActive, one, INT, root);
   if (nActive > 0 && mypid != 0) statusArray = new int[nActive]; 
   psCommMgr_->bcast((void *) statusArray, nActive, INT, root);
   psCommMgr_->bcast((void *) &npPerJob, one, INT, root);

#ifdef HAVE_MPICH
   MPI_Comm  newComm;
   MPI_Group baseGroup, newGroup;
   int  myGroup=-1, nGroups, *pgroup;

   nGroups = nprocs / npPerJob;
   if (nGroups <= 0)
   {
      printf("RunParallelLocal ERROR: not enough processors.\n");
      exit(1);
   }
   pgroup  = new int[npPerJob];
   for (ip = 0; ip < nGroups; ip++)
   {
      for (ii = 0; ii < npPerJob; ii++)
         pgroup[ii] = ip * npPerJob + ii;
      for (ii = 0; ii < npPerJob; ii++)
      {
         if (mypid == pgroup[ii])
         {
            myGroup = ip;
            MPI_Comm_group(MPI_COMM_WORLD, &baseGroup);
            MPI_Group_incl(baseGroup, npPerJob, pgroup, &newGroup);
            MPI_Comm_create(MPI_COMM_WORLD,newGroup,&newComm);
            break;
         }
      }
   }
   for (sampleID = 0; sampleID < nActive; sampleID++)
   {
      index = statusArray[sampleID];
      if (sampleID % nGroups == myGroup)
      {
         printf("Processor %d running job %d\n", mypid, index+1);
         UserFunction(newComm, index+1, workdir);
      }
   }
   delete [] pgroup;
   if (statusArray != NULL) delete [] statusArray;
#else
   printf("RunParallelLocal ERROR: MPI build not enabled.\n");
   exit(1);
#endif
   return 0;
}


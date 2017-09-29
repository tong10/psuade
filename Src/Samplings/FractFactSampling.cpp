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
// Functions for the fractional factorial sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <sstream>
using namespace std;

#include "Util/sysdef.h"
#include "Util/PsuadeUtil.h"
#include "FractFactSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FractFactSampling::FractFactSampling() : Sampling()
{
   samplingID_ = PSUADE_SAMP_FF4;
   resolution_ = 4;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
FractFactSampling::~FractFactSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int FractFactSampling::initialize(int initLevel)
{
   int maxInputs, nInputCnt, ii, increment, count, sampleID, inputID, sampleID2;

   if (nInputs_ == 0 || lowerBounds_ == NULL || upperBounds_ == NULL)
   {
      printf("FractFactSampling::initialize ERROR - input not set up.\n");
      exit(1);
   }
   if (nInputs_ > 32)
   {
      printf("FractFactSampling::initialize ERROR - nInputs > 32 not");
      printf(" supported.\n");
      exit(1);
   }

   deleteSampleData();
   //if (initLevel != 0) return 0;

resolution_ = 5;
printf("enter nInputs: ");
scanf("%d", &nInputs_);
   if (resolution_ == 4)
   {
      if (nInputs_ <= 3) nInputCnt = nInputs_;
      else
      {
         nInputCnt = 2;
         maxInputs = 0;
         while (maxInputs < nInputs_)
         {
            nInputCnt++;
            increment = 2;
            maxInputs = nInputCnt;
            while (nInputCnt-increment > 0)
            {
               increment++;
               count = 1;
               for (ii = nInputCnt; ii > nInputCnt-increment; ii--)
                  count = count * ii / (nInputCnt - ii + 1); 
               maxInputs += count;
            }
         }
      }
   }
   else if (resolution_ == 5)
   {
      if (nInputs_ <= 4) nInputCnt = nInputs_;
      else
      {
         nInputCnt = 3;
         maxInputs = 0;
         while (maxInputs < nInputs_)
         {
            nInputCnt++;
            increment = 3;
            maxInputs = nInputCnt;
            while (nInputCnt-increment > 0)
            {
               increment++;
               count = 1;
               for (ii = nInputCnt; ii > nInputCnt-increment; ii--)
                  count = count * ii / (nInputCnt - ii + 1); 
               maxInputs += count;
            }
         }
      }
   }
   nSamples_ = (int) pow(2.0001, (double) nInputCnt);
   printf("nSamples = %d\n", nSamples_);
   if (nInputs_ == 1) nSamples_ = 2;
   if (nInputs_ == 2) nSamples_ = 4;
   if (nInputs_ == 3) nSamples_ = 8;
   if (nInputs_ == 4 && resolution_ == 4) nSamples_ = 8;
   if (nInputs_ == 4 && resolution_ == 5) nSamples_ = 16;
   if (nInputs_ == 5) nSamples_ = 16;
   if (nInputs_ == 6 && resolution_ == 4) nSamples_ = 16;
   if (nInputs_ == 6 && resolution_ != 4) nSamples_ = 32;
   if (nInputs_ == 7) nSamples_ = 16;
   if (nInputs_ == 8) nSamples_ = 16;
   if (nInputs_ == 9 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 9 && resolution_ == 5) nSamples_ = 128;
   if (nInputs_ == 10 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 10 && resolution_ == 5) nSamples_ = 128;
   if (nInputs_ == 11 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 11 && resolution_ == 5) nSamples_ = 128;
   if (nInputs_ == 12 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 12 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = 12 Resolution V\n");
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   if (nInputs_ == 13 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 13 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = 13 Resolution V\n");
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   if (nInputs_ == 14 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 14 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = 14 Resolution V\n");
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   if (nInputs_ == 15 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 15 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = 15 Resolution V\n");
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   if (nInputs_ == 16 && resolution_ == 4) nSamples_ = 32;
   if (nInputs_ == 16 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = 16 Resolution V\n");
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   if (nInputs_ >= 17 && nInputs_ <= 32 && resolution_ == 4) nSamples_ = 64;
   if (nInputs_ >= 17 && nInputs_ <= 32 && resolution_ == 5)
   {
      printf("FractFactSampling::initialize ERROR - nInputs = %d Resolution V\n",
             nInputs_);
      printf("                              not supported. Try Resolution IV.\n");
      exit(1);
   }
   allocSampleData();

   if (printLevel_ > 4)
   {
      printf("FractFactSampling::initialize: nSamples = %d\n", nSamples_);
      printf("FractFactSampling::initialize: nInputs  = %d\n", nInputs_);
      printf("FractFactSampling::initialize: nOutputs = %d\n", nOutputs_);
      for (inputID = 0; inputID < nInputs_; inputID++)
         printf("    FractFactSampling input %3d = [%e %e]\n", inputID+1,
                lowerBounds_[inputID], upperBounds_[inputID]);
   }

   if (nInputs_ == 1)
   {
      sampleMatrix_[0][0] = -1; sampleMatrix_[1][0] = 1;
   }
   if (nInputs_ == 2)
   {
      sampleMatrix_[0][0] = -1; sampleMatrix_[0][1] = -1;
      sampleMatrix_[1][0] =  1; sampleMatrix_[1][1] = -1;
      sampleMatrix_[2][0] = -1; sampleMatrix_[2][1] =  1;
      sampleMatrix_[3][0] =  1; sampleMatrix_[3][1] =  1;
   }
   if (nSamples_ == 8)
   {
      for (sampleID = 0; sampleID < 8; sampleID+=2)
      {
         sampleMatrix_[sampleID][2] = -1; 
         sampleMatrix_[sampleID+1][2] = 1;
      }
      for (sampleID = 0; sampleID < 8; sampleID+=4)
      {
         sampleMatrix_[sampleID][1]   = -1;
         sampleMatrix_[sampleID+1][1] = -1;
         sampleMatrix_[sampleID+2][1] = 1;
         sampleMatrix_[sampleID+3][1] = 1;
      }
      for (sampleID = 0; sampleID < 4; sampleID++)
      {
         sampleMatrix_[sampleID][0] = -1;
         sampleMatrix_[sampleID+4][0] = 1;
      }
      if (nInputs_ == 4) // generator : 4 = 123
      {
         for (sampleID = 0; sampleID < 8; sampleID++)
            setPattern3(sampleMatrix_[sampleID],3,0,1,2);
      }
   }
   if (nSamples_ == 16)
   {
      for (sampleID = 0; sampleID < 16; sampleID+=2)
      {
         sampleMatrix_[sampleID][3] = -1;
         sampleMatrix_[sampleID+1][3] = 1;
      }
      for (sampleID = 0; sampleID < 16; sampleID+=4)
      {
         for (sampleID2 = 0; sampleID2 < 2; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][2] = -1;
            sampleMatrix_[sampleID+2+sampleID2][2] = 1;
         }
      }
      for (sampleID = 0; sampleID < 16; sampleID+=8)
      {
         for (sampleID2 = 0; sampleID2 < 4; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][1] = -1;
            sampleMatrix_[sampleID+4+sampleID2][1] = 1;
         }
      }
      for (sampleID = 0; sampleID < 8; sampleID++)
      {
         sampleMatrix_[sampleID][0] = -1;
         sampleMatrix_[sampleID+8][0] = 1;
      }
      if (nInputs_ == 5) // generator : 5 = 1234
      {
         for (sampleID = 0; sampleID < 16; sampleID++)
            setPattern4(sampleMatrix_[sampleID],4,0,1,2,3);
      }
      if (nInputs_ == 6) // generator : 5 = 123, 6 = 124
      {
         for (sampleID = 0; sampleID < 16; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],4,0,1,2);
            setPattern3(sampleMatrix_[sampleID],5,0,1,3);
         }
      }
      if (nInputs_ == 7) // generator : 5 = 123, 6 = 124, 7 = 134
      {
         for (sampleID = 0; sampleID < 16; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],4,0,1,2);
            setPattern3(sampleMatrix_[sampleID],5,0,1,3);
            setPattern3(sampleMatrix_[sampleID],6,0,2,3);
         }
      }
      if (nInputs_ == 8) // generator : 5 = 123, 6 = 124, 7 = 134, 8 = 234
      {
         for (sampleID = 0; sampleID < 16; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],4,0,1,2);
            setPattern3(sampleMatrix_[sampleID],5,0,1,3);
            setPattern3(sampleMatrix_[sampleID],6,0,2,3);
            setPattern3(sampleMatrix_[sampleID],7,1,2,3);
         }
      }
   }
   if (nSamples_ == 32)
   {
      for (sampleID = 0; sampleID < 32; sampleID+=2)
      {
         sampleMatrix_[sampleID][4] = -1;
         sampleMatrix_[sampleID+1][4] = 1;
      }
      for (sampleID = 0; sampleID < 32; sampleID+=4)
      {
         for (sampleID2 = 0; sampleID2 < 2; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][3] = -1;
            sampleMatrix_[sampleID+2+sampleID2][3] = 1;
         }
      }
      for (sampleID = 0; sampleID < 32; sampleID+=8)
      {
         for (sampleID2 = 0; sampleID2 < 4; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][2] = -1;
            sampleMatrix_[sampleID+4+sampleID2][2] = 1;
         }
      }
      for (sampleID = 0; sampleID < 32; sampleID+=16)
      {
         for (sampleID2 = 0; sampleID2 < 8; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][1] = -1;
            sampleMatrix_[sampleID+8+sampleID2][1] = 1;
         }
      }
      for (sampleID = 0; sampleID < 16; sampleID++)
      {
         sampleMatrix_[sampleID][0] = -1;
         sampleMatrix_[sampleID+16][0] = 1;
      }
      if (nInputs_ == 6) // generator : 6 = 12345
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
            setPattern5(sampleMatrix_[sampleID],5,0,1,2,3,4);
      }
      if (nInputs_ == 9) // generator : 6=123,7=124,8=125, 9=1345
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,1,4);
            setPattern4(sampleMatrix_[sampleID],8,0,2,3,4);
         }
      }
      if (nInputs_ == 10) // generator : 6=123,7=124,8=125,9=1345,A=2345
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,1,4);
            setPattern4(sampleMatrix_[sampleID],8,0,2,3,4);
            setPattern4(sampleMatrix_[sampleID],9,1,2,3,4);
         }
      }
      if (nInputs_ == 11) // 6=123,7=124,8=134,9=125,A=135,B=145
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,0,1,4);
            setPattern3(sampleMatrix_[sampleID],9,0,2,4);
            setPattern3(sampleMatrix_[sampleID],10,0,3,4);
         }
      }
      if (nInputs_ == 12) // 6=123,7=124,8=134,9=234,A=125,B=135,C=145
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,1,2,3);
            setPattern3(sampleMatrix_[sampleID],9,0,1,4);
            setPattern3(sampleMatrix_[sampleID],10,0,2,4);
            setPattern3(sampleMatrix_[sampleID],11,0,3,4);
         }
      }
      if (nInputs_ == 13) // 6=123,7=124,8=134,9=234,A=125,B=135,C=235,B=145
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,1,2,3);
            setPattern3(sampleMatrix_[sampleID],9,0,1,4);
            setPattern3(sampleMatrix_[sampleID],10,0,2,4);
            setPattern3(sampleMatrix_[sampleID],11,1,2,4);
            setPattern3(sampleMatrix_[sampleID],12,0,3,4);
         }
      }
      if (nInputs_ == 14) // 6=123,7=124,8=134,9=234,A=125,B=135,C=235,D=145,
                          // E=245
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,1,2,3);
            setPattern3(sampleMatrix_[sampleID],9,0,1,4);
            setPattern3(sampleMatrix_[sampleID],10,0,2,4);
            setPattern3(sampleMatrix_[sampleID],11,1,2,4);
            setPattern3(sampleMatrix_[sampleID],12,0,3,4);
            setPattern3(sampleMatrix_[sampleID],13,1,3,4);
         }
      }
      if (nInputs_ == 15) // 6=123,7=124,8=134,9=234,A=125,B=135,C=235,D=145,
                          // E=245,F=345
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,1,2,3);
            setPattern3(sampleMatrix_[sampleID],9,0,1,4);
            setPattern3(sampleMatrix_[sampleID],10,0,2,4);
            setPattern3(sampleMatrix_[sampleID],11,1,2,4);
            setPattern3(sampleMatrix_[sampleID],12,0,3,4);
            setPattern3(sampleMatrix_[sampleID],13,1,3,4);
            setPattern3(sampleMatrix_[sampleID],14,2,3,4);
         }
      }
      if (nInputs_ == 16) // 6=123,7=124,8=134,9=234,A=125,B=135,C=235,D=145,
                          // E=245,F=345,G=12345
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],5,0,1,2);
            setPattern3(sampleMatrix_[sampleID],6,0,1,3);
            setPattern3(sampleMatrix_[sampleID],7,0,2,3);
            setPattern3(sampleMatrix_[sampleID],8,1,2,3);
            setPattern3(sampleMatrix_[sampleID],9,0,1,4);
            setPattern3(sampleMatrix_[sampleID],10,0,2,4);
            setPattern3(sampleMatrix_[sampleID],11,1,2,4);
            setPattern3(sampleMatrix_[sampleID],12,0,3,4);
            setPattern3(sampleMatrix_[sampleID],13,1,3,4);
            setPattern3(sampleMatrix_[sampleID],14,2,3,4);
            setPattern5(sampleMatrix_[sampleID],15,0,1,2,3,4);
         }
      }
   }
   if (nSamples_ == 64)
   {
      for (sampleID = 0; sampleID < nSamples_; sampleID+=2)
      {
         sampleMatrix_[sampleID][5] = -1;
         sampleMatrix_[sampleID+1][5] = 1;
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=4)
      {
         for (sampleID2 = 0; sampleID2 < 2; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][4] = -1;
            sampleMatrix_[sampleID+2+sampleID2][4] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=8)
      {
         for (sampleID2 = 0; sampleID2 < 4; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][3] = -1;
            sampleMatrix_[sampleID+4+sampleID2][3] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=16)
      {
         for (sampleID2 = 0; sampleID2 < 8; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][2] = -1;
            sampleMatrix_[sampleID+8+sampleID2][2] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=32)
      {
         for (sampleID2 = 0; sampleID2 < 16; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][1] = -1;
            sampleMatrix_[sampleID+16+sampleID2][1] = 1;
         }
      }
      for (sampleID = 0; sampleID < 32; sampleID++)
      {
         sampleMatrix_[sampleID][0] = -1;
         sampleMatrix_[sampleID+32][0] = 1;
      }
      if (nInputs_ == 17) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=145,
                          // F=245,G=345,H=123456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            sampleMatrix_[sampleID][16] = sampleMatrix_[sampleID][0] * 
                sampleMatrix_[sampleID][1] * sampleMatrix_[sampleID][2] *
                sampleMatrix_[sampleID][3] * sampleMatrix_[sampleID][4] *
                sampleMatrix_[sampleID][5];
         }
      }
      if (nInputs_ == 18) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=126,
                          // F=136,G=1456,H=2456,I=3456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,1,5);
            setPattern3(sampleMatrix_[sampleID],14,0,2,5);
            setPattern4(sampleMatrix_[sampleID],15,0,3,4,5);
            setPattern4(sampleMatrix_[sampleID],16,1,3,4,5);
            setPattern4(sampleMatrix_[sampleID],17,2,3,4,5);
         }
      }
      if (nInputs_ == 19) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=126,
                          // F=136,G=236,H=1456,I=2456,J=3456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,1,5);
            setPattern3(sampleMatrix_[sampleID],14,0,2,5);
            setPattern3(sampleMatrix_[sampleID],15,1,2,5);
            setPattern4(sampleMatrix_[sampleID],16,0,3,4,5);
            setPattern4(sampleMatrix_[sampleID],17,1,3,4,5);
            setPattern4(sampleMatrix_[sampleID],18,2,3,4,5);
         }
      }
      if (nInputs_ == 20) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=126,
                          // F=136,G=236,H=1456,I=2456,J=3456,K=123456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,1,5);
            setPattern3(sampleMatrix_[sampleID],14,0,2,5);
            setPattern3(sampleMatrix_[sampleID],15,1,2,5);
            setPattern4(sampleMatrix_[sampleID],16,0,3,4,5);
            setPattern4(sampleMatrix_[sampleID],17,1,3,4,5);
            setPattern4(sampleMatrix_[sampleID],18,2,3,4,5);
            sampleMatrix_[sampleID][19] = sampleMatrix_[sampleID][0] * 
                sampleMatrix_[sampleID][1] * sampleMatrix_[sampleID][2] *
                sampleMatrix_[sampleID][3] * sampleMatrix_[sampleID][4] *
                sampleMatrix_[sampleID][5];
         }
      }
      if (nInputs_ == 21) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=145,
                          // F=126,G=146,H=246,I=156,J=356,K=456,L=23456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,0,1,5);
            setPattern3(sampleMatrix_[sampleID],15,0,3,5);
            setPattern3(sampleMatrix_[sampleID],16,1,3,5);
            setPattern3(sampleMatrix_[sampleID],17,0,4,5);
            setPattern3(sampleMatrix_[sampleID],18,2,4,5);
            setPattern3(sampleMatrix_[sampleID],19,3,4,5);
            setPattern5(sampleMatrix_[sampleID],20,1,2,3,4,5);
         }
      }
      if (nInputs_ == 22) // 7=123,8=124,9=134,A=234,B=125,C=135,D=235,E=145,
                          // F=126,G=136,H=146,I=246,J=156,K=356,L=456,M=23456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,0,1,5);
            setPattern3(sampleMatrix_[sampleID],15,0,2,5);
            setPattern3(sampleMatrix_[sampleID],16,0,3,5);
            setPattern3(sampleMatrix_[sampleID],17,1,3,5);
            setPattern3(sampleMatrix_[sampleID],18,0,4,5);
            setPattern3(sampleMatrix_[sampleID],19,2,4,5);
            setPattern3(sampleMatrix_[sampleID],20,3,4,5);
            setPattern5(sampleMatrix_[sampleID],21,1,2,3,4,5);
         }
      }
      if (nInputs_ == 23)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,0,1,5);
            setPattern3(sampleMatrix_[sampleID],16,0,2,5);
            setPattern3(sampleMatrix_[sampleID],17,0,3,5);
            setPattern3(sampleMatrix_[sampleID],18,2,3,5);
            setPattern3(sampleMatrix_[sampleID],19,0,4,5);
            setPattern3(sampleMatrix_[sampleID],20,2,4,5);
            setPattern3(sampleMatrix_[sampleID],21,3,4,5);
            setPattern5(sampleMatrix_[sampleID],22,1,2,3,4,5);
         }
      }
      if (nInputs_ == 24)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,0,1,5);
            setPattern3(sampleMatrix_[sampleID],16,0,2,5);
            setPattern3(sampleMatrix_[sampleID],17,1,2,5);
            setPattern3(sampleMatrix_[sampleID],18,0,3,5);
            setPattern3(sampleMatrix_[sampleID],19,1,3,5);
            setPattern3(sampleMatrix_[sampleID],20,0,4,5);
            setPattern3(sampleMatrix_[sampleID],21,2,4,5);
            setPattern3(sampleMatrix_[sampleID],22,3,4,5);
            setPattern5(sampleMatrix_[sampleID],23,1,2,3,4,5);
         }
      }
      if (nInputs_ == 25)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern3(sampleMatrix_[sampleID],16,0,1,5);
            setPattern3(sampleMatrix_[sampleID],17,0,2,5);
            setPattern3(sampleMatrix_[sampleID],18,1,2,5);
            setPattern3(sampleMatrix_[sampleID],19,0,3,5);
            setPattern3(sampleMatrix_[sampleID],20,1,3,5);
            setPattern3(sampleMatrix_[sampleID],21,0,4,5);
            setPattern3(sampleMatrix_[sampleID],22,2,4,5);
            setPattern3(sampleMatrix_[sampleID],23,3,4,5);
            setPattern5(sampleMatrix_[sampleID],24,1,2,3,4,5);
         }
      }
      if (nInputs_ == 26)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern3(sampleMatrix_[sampleID],16,0,1,5);
            setPattern3(sampleMatrix_[sampleID],17,0,2,5);
            setPattern3(sampleMatrix_[sampleID],18,1,2,5);
            setPattern3(sampleMatrix_[sampleID],19,0,3,5);
            setPattern3(sampleMatrix_[sampleID],20,1,3,5);
            setPattern3(sampleMatrix_[sampleID],21,2,3,5);
            setPattern3(sampleMatrix_[sampleID],22,0,4,5);
            setPattern3(sampleMatrix_[sampleID],23,1,4,5);
            setPattern3(sampleMatrix_[sampleID],24,2,4,5);
            setPattern3(sampleMatrix_[sampleID],25,3,4,5);
         }
      }
      if (nInputs_ == 27)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern3(sampleMatrix_[sampleID],23,0,4,5);
            setPattern3(sampleMatrix_[sampleID],24,1,4,5);
            setPattern3(sampleMatrix_[sampleID],25,2,4,5);
            setPattern3(sampleMatrix_[sampleID],26,3,4,5);
         }
      }
      if (nInputs_ == 28)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern5(sampleMatrix_[sampleID],23,0,1,2,3,5);
            setPattern3(sampleMatrix_[sampleID],24,0,4,5);
            setPattern3(sampleMatrix_[sampleID],25,1,4,5);
            setPattern3(sampleMatrix_[sampleID],26,2,4,5);
            setPattern3(sampleMatrix_[sampleID],27,3,4,5);
         }
      }
      if (nInputs_ == 29)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern5(sampleMatrix_[sampleID],23,0,1,2,3,5);
            setPattern3(sampleMatrix_[sampleID],24,0,4,5);
            setPattern3(sampleMatrix_[sampleID],25,1,4,5);
            setPattern3(sampleMatrix_[sampleID],26,2,4,5);
            setPattern5(sampleMatrix_[sampleID],27,0,1,2,4,5);
            setPattern3(sampleMatrix_[sampleID],28,3,4,5);
         }
      }
      if (nInputs_ == 30)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern5(sampleMatrix_[sampleID],23,0,1,2,3,5);
            setPattern3(sampleMatrix_[sampleID],24,0,4,5);
            setPattern3(sampleMatrix_[sampleID],25,1,4,5);
            setPattern3(sampleMatrix_[sampleID],26,2,4,5);
            setPattern5(sampleMatrix_[sampleID],27,0,1,2,4,5);
            setPattern3(sampleMatrix_[sampleID],28,3,4,5);
            setPattern5(sampleMatrix_[sampleID],29,0,1,3,4,5);
         }
      }
      if (nInputs_ == 31)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern5(sampleMatrix_[sampleID],23,0,1,2,3,5);
            setPattern3(sampleMatrix_[sampleID],24,0,4,5);
            setPattern3(sampleMatrix_[sampleID],25,1,4,5);
            setPattern3(sampleMatrix_[sampleID],26,2,4,5);
            setPattern5(sampleMatrix_[sampleID],27,0,1,2,4,5);
            setPattern3(sampleMatrix_[sampleID],28,3,4,5);
            setPattern5(sampleMatrix_[sampleID],29,0,1,3,4,5);
            setPattern5(sampleMatrix_[sampleID],30,0,2,3,4,5);
         }
      }
      if (nInputs_ == 32)
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern3(sampleMatrix_[sampleID],6,0,1,2);
            setPattern3(sampleMatrix_[sampleID],7,0,1,3);
            setPattern3(sampleMatrix_[sampleID],8,0,2,3);
            setPattern3(sampleMatrix_[sampleID],9,1,2,3);
            setPattern3(sampleMatrix_[sampleID],10,0,1,4);
            setPattern3(sampleMatrix_[sampleID],11,0,2,4);
            setPattern3(sampleMatrix_[sampleID],12,1,2,4);
            setPattern3(sampleMatrix_[sampleID],13,0,3,4);
            setPattern3(sampleMatrix_[sampleID],14,1,3,4);
            setPattern3(sampleMatrix_[sampleID],15,2,3,4);
            setPattern5(sampleMatrix_[sampleID],16,0,1,2,3,4);
            setPattern3(sampleMatrix_[sampleID],17,0,1,5);
            setPattern3(sampleMatrix_[sampleID],18,0,2,5);
            setPattern3(sampleMatrix_[sampleID],19,1,2,5);
            setPattern3(sampleMatrix_[sampleID],20,0,3,5);
            setPattern3(sampleMatrix_[sampleID],21,1,3,5);
            setPattern3(sampleMatrix_[sampleID],22,2,3,5);
            setPattern5(sampleMatrix_[sampleID],23,0,1,2,3,5);
            setPattern3(sampleMatrix_[sampleID],24,0,4,5);
            setPattern3(sampleMatrix_[sampleID],25,1,4,5);
            setPattern3(sampleMatrix_[sampleID],26,2,4,5);
            setPattern5(sampleMatrix_[sampleID],27,0,1,2,4,5);
            setPattern3(sampleMatrix_[sampleID],28,3,4,5);
            setPattern5(sampleMatrix_[sampleID],29,0,1,3,4,5);
            setPattern5(sampleMatrix_[sampleID],30,0,2,3,4,5);
            setPattern5(sampleMatrix_[sampleID],31,1,2,3,4,5);
         }
      }
   }
   if (nSamples_ == 128)
   {
      for (sampleID = 0; sampleID < nSamples_; sampleID+=2)
      {
         sampleMatrix_[sampleID][6] = -1;
         sampleMatrix_[sampleID+1][6] = 1;
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=4)
      {
         for (sampleID2 = 0; sampleID2 < 2; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][5] = -1;
            sampleMatrix_[sampleID+2+sampleID2][5] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=8)
      {
         for (sampleID2 = 0; sampleID2 < 4; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][4] = -1;
            sampleMatrix_[sampleID+4+sampleID2][4] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=16)
      {
         for (sampleID2 = 0; sampleID2 < 8; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][3] = -1;
            sampleMatrix_[sampleID+8+sampleID2][3] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=32)
      {
         for (sampleID2 = 0; sampleID2 < 16; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][2] = -1;
            sampleMatrix_[sampleID+16+sampleID2][2] = 1;
         }
      }
      for (sampleID = 0; sampleID < nSamples_; sampleID+=64)
      {
         for (sampleID2 = 0; sampleID2 < 32; sampleID2++)
         {
            sampleMatrix_[sampleID+sampleID2][1] = -1;
            sampleMatrix_[sampleID+32+sampleID2][1] = 1;
         }
      }
      for (sampleID = 0; sampleID < 64; sampleID++)
      {
         sampleMatrix_[sampleID][0] = -1;
         sampleMatrix_[sampleID+64][0] = 1;
      }
      if (nInputs_ == 9) // 8=13457, 9=12467
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern5(sampleMatrix_[sampleID],7,0,2,3,4,6);
            setPattern5(sampleMatrix_[sampleID],8,0,1,3,5,6);
         }
      }
      if (nInputs_ == 10) // 8=13457, 9=12467, A=3456
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern5(sampleMatrix_[sampleID],7,0,2,3,4,6);
            setPattern5(sampleMatrix_[sampleID],8,0,1,3,5,6);
            setPattern4(sampleMatrix_[sampleID],9,2,3,4,5);
         }
      }
      if (nInputs_ == 11) // 8=13457, 9=12467, A=3456, B=2567
      {
         for (sampleID = 0; sampleID < nSamples_; sampleID++)
         {
            setPattern5(sampleMatrix_[sampleID],7,0,2,3,4,6);
            setPattern5(sampleMatrix_[sampleID],8,0,1,3,5,6);
            setPattern4(sampleMatrix_[sampleID],9,2,3,4,5);
            setPattern4(sampleMatrix_[sampleID],10,1,4,5,6);
         }
      }
   }

   for (sampleID = 0; sampleID < nSamples_; sampleID++)
   {
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
         if (sampleMatrix_[sampleID][inputID] == -1)
            sampleMatrix_[sampleID][inputID] = lowerBounds_[inputID];
         else if (sampleMatrix_[sampleID][inputID] == 1)
            sampleMatrix_[sampleID][inputID] = upperBounds_[inputID];
      }
   }
   return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int FractFactSampling::refine(int refineRatio, int randomize, double thresh,
                              int nSamples, double *sampleErrors)
{
   (void) refineRatio;
   (void) randomize;
   (void) thresh;
   (void) nSamples;
   (void) sampleErrors;

   printf("FractFactSampling::refine ERROR - not available.\n");
   exit(1);
   return 0;
}

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int FractFactSampling::setParam(string sparam)
{
   istringstream buffer;
   int           pos = sparam.find("setResolution");
   string        substr;
                                                                                
   if (pos >= 0)
   {
      substr = sparam.substr(14);
      buffer.str(substr);
      buffer >> resolution_;
      if (resolution_ != 4 && resolution_ != 5) resolution_ = 4;
      if (resolution_ == 4) samplingID_ = PSUADE_SAMP_FF4;
      else                  samplingID_ = PSUADE_SAMP_FF5;
   }
   return 0;
}

// ************************************************************************
// set sample pattern 
// ------------------------------------------------------------------------
void FractFactSampling::setPattern3(double *sample, int ind, int ind1,
                                    int ind2, int ind3)
{
   sample[ind] = sample[ind1] * sample[ind2] * sample[ind3];
}

// ************************************************************************
// set sample pattern 
// ------------------------------------------------------------------------
void FractFactSampling::setPattern4(double *sample, int ind, int ind1,
                                    int ind2, int ind3, int ind4)
{
   sample[ind] = sample[ind1] * sample[ind2] * sample[ind3] * sample[ind4];
}

// ************************************************************************
// set sample pattern 
// ------------------------------------------------------------------------
void FractFactSampling::setPattern5(double *sample, int ind, int ind1,
                                    int ind2, int ind3, int ind4, int ind5)
{
   sample[ind] = sample[ind1]*sample[ind2]*sample[ind3]*sample[ind4]*
                 sample[ind5];
}


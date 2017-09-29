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
// Functions for PsuadeConfig 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PsuadeConfig.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
PsuadeConfig::PsuadeConfig(char *fname, int printLevel)
{
   int  ii, lineLeng=500;
   char lineIn[500], winput[500];
   FILE *fIn;

   fIn = fopen(fname, "r");
   if (fIn == NULL) 
   {
      printf("PsuadeConfig ERROR:: configure file %s not found.\n",fname);
      fclose(fIn);
      exit(1);
   }

   while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(lineIn, "%s", winput);
      if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
   }
   if (feof(fIn) != 0)
   {
      printf("PsuadeConfig ERROR:: keyword PSUADE_CONFIG not found.\n");
      fclose(fIn);
      exit(1);
   }

   nLines_ = 0;
   while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(lineIn,"%s", winput);
      if   (strcmp(winput, "PSUADE_END") == 0) break;
      else nLines_++;
   }
   if (strcmp(winput, "PSUADE_END") != 0)
   {
      printf("PsuadeConfig ERROR:: keyword PSUADE_END not found.\n");
      fclose(fIn);
      exit(1);
   }
   fclose(fIn);

   if (nLines_ == 0) fileData_ = NULL;
   else
   {
      fIn = fopen(fname, "r");
      while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
      {
         sscanf(lineIn,"%s", winput);
         if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
      }
      fileData_ = new char*[nLines_];
      for (ii = 0; ii < nLines_; ii++)
      {
         fileData_[ii] = new char[lineLeng];
         fgets(fileData_[ii], lineLeng, fIn);
      }
      fclose(fIn);
   }
   printLevel_ = printLevel;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------ 
PsuadeConfig::~PsuadeConfig()
{ 
   int ii;

   if (nLines_ > 0 && fileData_ != NULL)
   {
      for (ii = 0; ii < nLines_; ii++)
         if (fileData_[ii] != NULL) delete [] fileData_[ii];
      delete [] fileData_;
   }
}

// ************************************************************************
// request data from this object 
// ------------------------------------------------------------------------ 
char *PsuadeConfig::getParameter(const char *keyword)
{ 
   int  ii;
   char firstWord[100];

   for (ii = 0; ii < nLines_; ii++)
   {
      sscanf(fileData_[ii], "%s", firstWord);
      if (strcmp(keyword,firstWord) == 0)
         return fileData_[ii];
   }
   return NULL;
}

// ************************************************************************
// write to file 
// ------------------------------------------------------------------------ 
void PsuadeConfig::writeToFile(char *fname)
{
   int  ii;
   FILE *fOut;

   fOut = fopen(fname, "w");
   if (fOut == NULL) 
   {
      printf("PsuadeConfig ERROR:: cannot write to configure file %s.\n",
             fname);
      fclose(fOut);
      exit(1);
   }

   fprintf(fOut, "PSUADE_CONFIG\n");
   for (ii = 0; ii < nLines_; ii++)
      fprintf(fOut, "%s\n", fileData_[ii]);
   fprintf(fOut, "PSUADE_END\n");
   fclose(fOut);
}


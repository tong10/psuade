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
// Utility functions 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include <string>

// ------------------------------------------------------------------------
// local defines
// ------------------------------------------------------------------------
static double  elapsedTime_;
static randctx rctx_; 
static int     randFlag_=0;
static int     randMask_;
static double  pgplotXmin_, pgplotXmax_, pgplotYmin_, pgplotYmax_;
static long    pgplotFlag_; 

// ************************************************************************
// random number generator initializer 
// ------------------------------------------------------------------------
void PSUADE_randInit()
{
   int    ntime, ii, nbits;
   double dtime;

   if (psRandomSeed_ == -1)
   {
      for (ii = 0; ii < RANDSIZ; ii++) 
      {
         dtime = getClock();
         ntime = (int) (dtime * 1.0E6);
         ntime = ntime % 10000;
         dtime = getClock();
         ntime = ntime * 10000 + (int) (dtime * 1.0E6);
         rctx_.randrsl[ii] = ntime;
         dtime = getClock();
         ntime = (int) (dtime * 1.0E6);
         ntime = ntime % 10000;
         dtime = getClock();
         ntime = ntime * 10000 + (int) (dtime * 1.0E6);
         rctx_.randmem[ii] = ntime;
      }
      randinit(&rctx_, TRUE);
   }
   else
   {
      randinitBySeed(&rctx_, psRandomSeed_);
   }
   nbits = sizeof(int) * 8;  
   randMask_ = (1 << (nbits-2)) + ((1 << (nbits-2)) - 1);
   randFlag_ = 1;
}

// ************************************************************************
// random number generator 
// ------------------------------------------------------------------------
int PSUADE_rand()
{
   int  irand;

   if (randFlag_ == 0) PSUADE_randInit(); 

   irand = ISAAC_RAND_RAND(&rctx_);
   irand = irand & randMask_;
      
   return irand;
}

// ************************************************************************
// random number generator (double) 
// ------------------------------------------------------------------------
double PSUADE_drand()
{
   ub4 irand;
   double drand;

   irand = PSUADE_rand();
   drand = ((double) irand / (double) randMask_);
   if (drand > 1 - 1.2e-7) drand = 1.0 - 1.2e-7;
   return drand;
}

// ************************************************************************
// copies a string into a new malloc'd pointer
// ************************************************************************
char* PSUADE_strdup(const char *s)
{
  size_t len = strlen (s) + 1;
  char *result = (char*) malloc (len);
  if (result == (char*) 0)
    return (char*) 0;
  return (char*) memcpy (result, s, len);
}

// ************************************************************************
// copies a n characters of a string into a new malloc'd pointer
// ************************************************************************
char* PSUADE_strndup(const char *s, size_t n)
{
  char *result;
  size_t len = strlen (s);

  if (n < len)
    len = n;

  result = (char *) malloc (len + 1);
  if (!result)
    return 0;

  result[len] = '\0';
  return (char *) memcpy (result, s, len);
}

// ************************************************************************
// generate an array of random number using standard normal distribution
// ------------------------------------------------------------------------
void PSUADE_drandn(int length, double *rArray)
{    
   int    ii, newLength;
   double R, *rData, pi, cTheta;

   newLength = length;
   pi = 4.0 * atan(1.0e0);
   if (length % 2 == 1) newLength++;
   rData = new double[newLength];
   for (ii = 0; ii < newLength; ii++) rArray[ii] = PSUADE_drand();
   for (ii = 0; ii < newLength; ii=+2)
   {
      R = -2.0 * log(rArray[ii]);
      if (R < 0.0) R = 0.0;
      R = sqrt(R);
      cTheta = 2.0 * pi * rArray[ii+1];
      rArray[ii] = R * cos(cTheta);
      if ((ii+1) < length) rArray[ii+1] = R * sin(cTheta);
   }
   delete [] rData;
}

// ************************************************************************
// sortIntList : sort an integer list
// ------------------------------------------------------------------------
void sortIntList(int length, int *intList)
{
   int i, last, mid, idata;

   if ( length <= 1 ) return;
   mid = ( length - 1 ) / 2;
   idata = intList[0];
   intList[0] = intList[mid];
   intList[mid] = idata;
   last = 0;
   for (i = 1; i < length; i++)
   {
      if (intList[i] < intList[0])
      {
         last++;
         idata = intList[last];
         intList[last] = intList[i];
         intList[i] = idata;
      }
   }
   idata = intList[0];
   intList[0] = intList[last];
   intList[last] = idata;
   sortIntList(last,intList);
   sortIntList(length-last-1,&(intList[last+1]));
}

// ************************************************************************
// sort integer list with an integer and a double array 
// ------------------------------------------------------------------------
void sortIntList2a(int length, int *intList, double *valueList)
{
   int    mid, i, idata, last, isort, *intList2, *intList3;
   double ddata, *valueList2, *valueList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   idata = intList[0];
   intList[0] = intList[mid];
   intList[mid] = idata;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   last = 0;
   isort = intList[0];
   intList2 = &(intList[last]);
   valueList2 = &(valueList[last]);
   intList3 = &(intList[1]);
   valueList3 = &(valueList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *intList3 < isort )
      {
         last++;
         intList2++; valueList2++;
         idata = *intList2;
         *intList2 = *intList3;
         *intList3 = idata;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
      }
      intList3++; valueList3++;
   }
   idata = intList[0];
   intList[0] = intList[last];
   intList[last] = idata;
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   sortIntList2a(last, intList, valueList);
   sortIntList2a(length-last-1, &(intList[last+1]), &(valueList[last+1]));
}

// ************************************************************************
// sort double list 
// ------------------------------------------------------------------------
void sortDbleList(int length, double *valueList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
      }
      valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   sortDbleList(last, valueList);
   sortDbleList(length-last-1,&(valueList[last+1]));
}

// ************************************************************************
// sort double list with 2 double arrays 
// ------------------------------------------------------------------------
void sortDbleList2(int length, double *valueList, double *dataList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   ddata = dataList[0];
   dataList[0] = dataList[mid];
   dataList[mid] = ddata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   dataList2 = &(dataList[last]);
   valueList3 = &(valueList[1]);
   dataList3 = &(dataList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         dataList2++; valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         ddata = *dataList2;
         *dataList2 = *dataList3;
         *dataList3 = ddata;
      }
      dataList3++; valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   ddata = dataList[0];
   dataList[0] = dataList[last];
   dataList[last] = ddata;
   sortDbleList2(last, valueList, dataList);
   sortDbleList2(length-last-1,&(valueList[last+1]),&(dataList[last+1]));
}

// ************************************************************************
// sort double list with 1 double array and 1 int array
// ------------------------------------------------------------------------
void sortDbleList2a(int length, double *valueList, int *dataList)
{
   int    mid, i, last, idata, *dataList2, *dataList3;
   double ddata, *valueList2, *valueList3, dsort;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   idata = dataList[0];
   dataList[0] = dataList[mid];
   dataList[mid] = idata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   dataList2 = &(dataList[last]);
   valueList3 = &(valueList[1]);
   dataList3 = &(dataList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         dataList2++; valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         idata = *dataList2;
         *dataList2 = *dataList3;
         *dataList3 = idata;
      }
      dataList3++; valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   idata = dataList[0];
   dataList[0] = dataList[last];
   dataList[last] = idata;
   sortDbleList2a(last, valueList, dataList);
   sortDbleList2a(length-last-1,&(valueList[last+1]),&(dataList[last+1]));
}

// ************************************************************************
// sort double list with 3 double arrays 
// ------------------------------------------------------------------------
void sortDbleList3(int length, double *valueList, double *dataList,
                   double *auxList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;
   double *auxList2, *auxList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[mid];
      dataList[mid] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[mid];
      auxList[mid] = ddata;
   }
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   if (dataList != NULL)
   {
      dataList2 = &(dataList[last]);
      dataList3 = &(dataList[1]);
   }
   if (auxList != NULL)
   {
      auxList2 = &(auxList[last]);
      auxList3 = &(auxList[1]);
   }
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         if (dataList != NULL)
         {
            dataList2++;
            ddata = *dataList2;
            *dataList2 = *dataList3;
            *dataList3 = ddata;
         }
         if (auxList != NULL)
         {
            auxList2++;
            ddata = *auxList2;
            *auxList2 = *auxList3;
            *auxList3 = ddata;
         }
      }
      valueList3++; 
      if (dataList != NULL) dataList3++; 
      if (auxList  != NULL) auxList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[last];
      dataList[last] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[last];
      auxList[last] = ddata;
   }
   sortDbleList3(last, valueList, dataList, auxList);
   if (dataList != NULL && auxList != NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]));
   else if (dataList != NULL && auxList == NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL);
   else if (dataList == NULL && auxList != NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]));
   else if (dataList == NULL && auxList == NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),NULL,NULL);
}

// ************************************************************************
// sort double list with 4 double arrays 
// ------------------------------------------------------------------------
void sortDbleList4(int length, double *valueList, double *dataList,
                   double *auxList, double *moreList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;
   double *auxList2, *auxList3, *moreList2, *moreList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[mid];
      dataList[mid] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[mid];
      auxList[mid] = ddata;
   }
   if (moreList != NULL)
   {
      ddata = moreList[0];
      moreList[0] = moreList[mid];
      moreList[mid] = ddata;
   }
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   if (dataList != NULL)
   {
      dataList2 = &(dataList[last]);
      dataList3 = &(dataList[1]);
   }
   if (auxList != NULL)
   {
      auxList2 = &(auxList[last]);
      auxList3 = &(auxList[1]);
   }
   if (moreList != NULL)
   {
      moreList2 = &(moreList[last]);
      moreList3 = &(moreList[1]);
   }
   for (i = 1; i < length; i++)
   {
      if (*valueList3 < dsort)
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         if (dataList != NULL)
         {
            dataList2++;
            ddata = *dataList2;
            *dataList2 = *dataList3;
            *dataList3 = ddata;
         }
         if (auxList != NULL)
         {
            auxList2++;
            ddata = *auxList2;
            *auxList2 = *auxList3;
            *auxList3 = ddata;
         }
         if (moreList != NULL)
         {
            moreList2++;
            ddata = *moreList2;
            *moreList2 = *moreList3;
            *moreList3 = ddata;
         }
      }
      valueList3++;
      if (dataList != NULL) dataList3++;
      if (auxList  != NULL) auxList3++;
      if (moreList != NULL) moreList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[last];
      dataList[last] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[last];
      auxList[last] = ddata;
   }
   if (moreList != NULL)
   {
      ddata = moreList[0];
      moreList[0] = moreList[last];
      moreList[last] = ddata;
   }
   sortDbleList4(last, valueList, dataList, auxList, moreList);
   if (dataList != NULL && auxList != NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]), &(moreList[last+1]));
   else if (dataList != NULL && auxList != NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]), NULL);
   else if (dataList != NULL && auxList == NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL, &(moreList[last+1]));
   else if (dataList != NULL && auxList == NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL, NULL);
   else if (dataList == NULL && auxList != NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]), &(moreList[last+1]));
   else if (dataList == NULL && auxList != NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]), NULL);
   else if (dataList == NULL && auxList == NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL, NULL,
                    &(moreList[last+1]));
   else if (dataList == NULL && auxList == NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,NULL, NULL);
}

// ************************************************************************
// special sort and delet (very specialized for sparse grid)
// sort double list with 4 double arrays 
// ------------------------------------------------------------------------
int sortAndDelete(int nrows, int ncols, double **ddata)
{
   int    *isortList, ii, jj, kk, ll, ibegin, iclose, nrows2;
   double *dsortList;
   
   // range checking by Bill Oliver
   if(nrows <= 0){ 
     printf("First parameter to function sortAndDelete is <= 0 in file %s \
              and must be >0 returning -1", __FILE__);
     return -1;
   }

   if(ncols <= 0){ 
     printf("Second parameter to function sortAndDelete is <= 0 in file %s \
              and must be >0 returning -1", __FILE__);
     return -1;
   }
             
   dsortList = (double *) malloc(nrows * sizeof(double));
   isortList = (int *) malloc(nrows * sizeof(int));

   // Bill Oliver check for out of memory
   if(dsortList == NULL || isortList == NULL){
     printf("out of memory in file %s line %d aborting\n", __FILE__, __LINE__);
     abort();
   }

   for (ii = 0; ii < ncols-1; ii++)
   {
      if (ii == 0)
      {
         for (jj = 0; jj < nrows; jj++)
         {
            dsortList[jj] = ddata[jj][ii];
            isortList[jj] = jj;
         }
         sortDbleList2a(nrows, dsortList, isortList);
         for (jj = 0; jj < nrows; jj++)
            ddata[jj][ii] = dsortList[jj] ;
         for (kk = ii+1; kk < ncols; kk++)
         {
            for (jj = 0; jj < nrows; jj++)
               dsortList[jj] = ddata[isortList[jj]][kk];
            for (jj = 0; jj < nrows; jj++) ddata[jj][kk] = dsortList[jj];
         }
      }
      else
      {
         ibegin = 0;
         while (ibegin < nrows)
         {
           iclose = ibegin + 1;
           for (jj = iclose; jj < nrows; jj++)
           {
              if (ddata[jj][ii-1] != ddata[jj-1][ii-1]) break;
              else
              {
                 for (kk = 0; kk < ii-1; kk++)
                    if (ddata[jj][kk] != ddata[jj-1][kk]) break;
                 if (kk < ii-1) break;
              }
           }
           iclose = jj;
           if (iclose - ibegin > 1)
           {
              for (jj = ibegin; jj < iclose; jj++)
              {
                 dsortList[jj-ibegin] = ddata[jj][ii];
                 isortList[jj-ibegin] = jj;
              }
              sortDbleList2a(iclose-ibegin, dsortList, isortList);
              for (jj = ibegin; jj < iclose; jj++)
                 ddata[jj][ii] = dsortList[jj-ibegin];
              for (kk = ii+1; kk < ncols; kk++)
              {
                 for (jj = 0; jj < iclose-ibegin; jj++)
                    dsortList[jj] = ddata[isortList[jj]][kk];
                 for (jj = 0; jj < iclose-ibegin; jj++)
                     ddata[jj+ibegin][kk] = dsortList[jj];
              }
            }
            ibegin = iclose;
         }
      }
   }
   free(dsortList);
   free(isortList);
   /* delete */
   kk = 1;
   nrows2 = nrows;
   while (kk < nrows2)
   {
      for (jj = 0; jj < ncols-1; jj++)
         if (ddata[kk-1][jj] != ddata[kk][jj]) break;
      if (jj == ncols-1)
      {
         ddata[kk-1][ncols-1] += ddata[kk][ncols-1];
         for (jj = kk+1; jj < nrows2; jj++)
            for (ll = 0; ll < ncols; ll++) ddata[jj-1][ll] = ddata[jj][ll];
         nrows2--;
      }
      else kk++;
   }
   return nrows2;
}
// ************************************************************************
// search
// ------------------------------------------------------------------------
int binarySearchInt(int sKey, int *sData, int length)
{
   int  left, right, mid, found, index;

   if (length <= 0) return -1;
   left  = 0;
   right = length - 1;
   if (sKey > sData[right]) return -(right+1);
   if (sKey < sData[left])  return -(left+1);
   found = 0;
   while ((found == 0) && ((right-left)>1))
   {
      mid = (left + right) / 2;
      if      (sKey == sData[mid]) {index  = mid; found = 1;}
      else if (sKey >  sData[mid])  left  = mid;
      else                          right = mid;
   }
   if (found == 1)                return index;
   else if (sKey == sData[left])  return left;
   else if (sKey == sData[right]) return right;
   else                           return -(left+1);
}

// ************************************************************************
// output a line of the provided character
// ------------------------------------------------------------------------
void printCharLine(int printLevel, int length, char printchar)
{
  if(length <= 0) {
    length = 70;
  }
  std::string outString = std::string(length, printchar);
  printOutTS(printLevel, "%s\n", outString.c_str());
}


// ************************************************************************
// output a line of asterisks
// ------------------------------------------------------------------------
void printAsterisks(int length)
{
  printCharLine(2, length, '*');
}

// ************************************************************************
// output a line of dashes
// ------------------------------------------------------------------------
void printDashes(int length)
{
  printCharLine(2, length, '-');
}

// ************************************************************************
// output a line of equal signs
// ------------------------------------------------------------------------
void printEquals(int length)
{
  printCharLine(2, length, '=');
}

// ************************************************************************
// search
// ------------------------------------------------------------------------
int binarySearchDble(double sKey, double *sData, int length)
{
   int  left, right, mid, found, index;

   if (length <= 0) return -1;
   left  = 0;
   right = length - 1;
   if (sKey > sData[right]) return -(right+1);
   if (sKey < sData[left])  return -(left+1);
   found = 0;
   while ((found == 0) && ((right-left)>1))
   {
      mid = (left + right) / 2;
      if      (sKey == sData[mid]) {index  = mid; found = 1;}
      else if (sKey >  sData[mid])  left  = mid;
      else                          right = mid;
   }
   if (found == 1)                return index;
   else if (sKey == sData[left])  return left;
   else if (sKey == sData[right]) return right;
   else                           return -(left+1);
}

// ************************************************************************
// timer functions
// ************************************************************************

// ************************************************************************
// get the time of the day in seconds
// ------------------------------------------------------------------------
double getClock()
{
   double time_i;
   double time_d;
   struct timeval tp;
   struct timezone tzp;
   Gettimeofday(&tp,&tzp);
   time_i = tp.tv_sec % 10;
   time_d = (double) time_i;
   time_i = tp.tv_usec;
   time_d = time_d + (double) time_i / 1000000.0; 
   return(time_d);
}

// ************************************************************************
// startTimer 
// ------------------------------------------------------------------------
void startTimer()
{
   elapsedTime_ = getClock();
}

// ************************************************************************
// StopTimer 
// ------------------------------------------------------------------------
void StopTimer()
{
   elapsedTime_ = getClock() - elapsedTime_;
}

// ************************************************************************
// Get the elapsed time (from startTimer to stopTimer)
// ------------------------------------------------------------------------
double getElapsedTime()
{
   return elapsedTime_;
}

// ************************************************************************
// plot facilities : call this before calling anything else
// ------------------------------------------------------------------------
void Plotbegin(double xmin, double xmax, double ymin, double ymax)
{
   pgplotFlag_ = 0;
   pgplotXmin_ = xmin;
   pgplotXmax_ = xmax;
   pgplotYmin_ = ymin;
   pgplotYmax_ = ymax;
}

// ************************************************************************
// plot 2D data
// ------------------------------------------------------------------------
void PlotSamples2D(long n, double *x, double *y, long *nstat)
{
#ifdef IRIX
   plot2dirix_(&n, x, y, nstat, &pgplotXmin_, &pgplotXmax_, 
               &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
#else
   plotsamples2d_(&n, x, y, nstat, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
#endif
}

// ************************************************************************
// 2D scatter plots
// ------------------------------------------------------------------------
void PlotScatter2D(long n, double *x, double *y)
{
   plotscatter2d_(&n, x, y, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
}

// ************************************************************************
// 2D multiple scatter plots
// ------------------------------------------------------------------------
void PlotScatterM2D(long n, double *x, double *x2, double *y)
{
   plotscatterm2d_(&n, x, x2, y, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
}

// ************************************************************************
// plot 3D data
// ------------------------------------------------------------------------
void Plot3d(long n, double *, double *, double *z)
{
   plot3d_(&n, z, &pgplotFlag_);
}

void Plotend()
{
   plotend_();
}

// ************************************************************************
// generate random integer vector
// ------------------------------------------------------------------------
void generateRandomIvector(int leng, int *vec)
{
   int    k, j, irand, *iArray;

   iArray = new int[leng];
   for ( j = 0; j < leng; j++ ) iArray[j] = j; 
   k = leng;
   while (k > 1)
   {
      irand = PSUADE_rand() % k;
      vec[leng-k] = iArray[irand];
      for ( j = irand; j < k-1; j++ ) iArray[j] = iArray[j+1];
      k--; 
   }
   if (k == 1) vec[leng-1] = iArray[0];
   delete [] iArray;
}

// ************************************************************************
// check a number is prime
// ------------------------------------------------------------------------
int checkPrime(int idata)
{
  int  i, upper;

  upper = (int) (sqrt((double) idata) + 1.0);
  if (idata <= 5) return 1;
  else 
  {
     for ( i = 2; i < upper; i++ )
       if (((i & 1) != 0) && (idata / i * i == idata)) return 0;
  }
  return 1;
}

// ************************************************************************
// factorize a number
// ------------------------------------------------------------------------
int *factorize(int iValue)
{
   int  ind, itemp, sqrtN, *factors, nFactors;

   factors = new int[201];
   sqrtN    = (int)(sqrt((double) iValue));
   itemp    = iValue;
   ind      = 2;
   nFactors = 0;

   while (ind <= sqrtN)
   {
      if (itemp % ind == 0)
      {
         if (nFactors > 200)
         {
            printf("Factorize ERROR : factors > 200.\n");
            exit(1);
         } 
         factors[nFactors+1] = ind;
         nFactors++;
         itemp /= ind;
      }
      else ind++;
   }
   factors[0] = nFactors;
   return factors;
}

// ************************************************************************
// compare samples to see if there are duplicates
// ------------------------------------------------------------------------
int compareSamples(int sampleID, int nSamples, int nInputs,
                   double *sampleData, int *sampleStates)
{
   int  sampleBase, indexBase, found=0, ss, ii;

   sampleBase = sampleID * nInputs;
   for (ss = 0; ss < nSamples; ss++)
   {
      if ((ss != sampleID) && (sampleStates[ss] > 0))
      {
         indexBase = ss * nInputs;
         for (ii = 0; ii < nInputs; ii++)
            if (sampleData[sampleBase+ii] != sampleData[indexBase+ii]) break;
         if (ii == nInputs)
         {
            found = 1;
            break;
         }
      }
   } 
   if (found == 0) return -1;
   else            return ss;
}

// ************************************************************************
// compute the number of terms in the permutations
// ------------------------------------------------------------------------
int computeNumPCEPermutations(int nRVs, int pOrder)
{
   int ii, nTerms;

   nTerms = 1;
   for (ii = nRVs; ii < nRVs+pOrder; ii++) nTerms *= (ii + 1);
   for (ii = nRVs; ii < nRVs+pOrder; ii++) nTerms /= (ii - nRVs + 1);
   return nTerms;
}


// ************************************************************************
// plot stuff for scilab or matlab
// ------------------------------------------------------------------------
int fwritePlotAxes(FILE *fp)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.thickness = 2;\n");
      fprintf(fp, "a.font_size = 3;\n");
      fprintf(fp, "a.font_style = 4;\n");
      fprintf(fp, "a.box = \"on\";\n");
      fprintf(fp, "a.grid = [1 1];\n");
   }
   else
   {
      fprintf(fp, "set(gca,'linewidth',2)\n");
      fprintf(fp, "set(gca,'fontweight','bold')\n");
      fprintf(fp, "set(gca,'fontsize',12)\n");
      fprintf(fp, "grid on\n");
      fprintf(fp, "box on\n");
   }
   return 0;
}
// ************************************************************************
int fwritePlotAxesNoGrid(FILE *fp)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.thickness = 2;\n");
      fprintf(fp, "a.font_size = 3;\n");
      fprintf(fp, "a.font_style = 4;\n");
      fprintf(fp, "a.box = \"on\";\n");
   }
   else
   {
      fprintf(fp, "set(gca,'linewidth',2)\n");
      fprintf(fp, "set(gca,'fontweight','bold')\n");
      fprintf(fp, "set(gca,'fontsize',12)\n");
      fprintf(fp, "box on\n");
   }
   return 0;
}
// ************************************************************************
int fwritePlotXLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.x_label.text = \"%s\";\n", label);
      fprintf(fp, "a.x_label.font_size = 3;\n");
      fprintf(fp, "a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "xlabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotYLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.y_label.text = \"%s\";\n", label);
      fprintf(fp, "a.y_label.font_size = 3;\n");
      fprintf(fp, "a.y_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "ylabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotZLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.z_label.text = \"%s\";\n", label);
      fprintf(fp, "a.z_label.font_size = 3;\n");
      fprintf(fp, "a.z_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "zlabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotTitle(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.title.text = \"%s\";\n", label);
      fprintf(fp, "a.title.font_size = 4;\n");
      fprintf(fp, "a.title.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "title('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotCLF(FILE *fp)
{
   if (psPlotTool_ == 1) fprintf(fp, "clf();\n");
   else                  fprintf(fp, "clf\n");
   return 0;
}
// ************************************************************************
int fwritePlotFigure(FILE *fp, int num)
{
   if (psPlotTool_ == 1) fprintf(fp, "scf(%d);\n", num);
   else                  fprintf(fp, "figure(%d)\n", num);
   return 0;
}


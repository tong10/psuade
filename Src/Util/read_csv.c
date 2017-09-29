#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LBUFSIZ 1000

int read_csv(char *filename, int *nSamples, int *nInOuts, double **XX)
{
   FILE *myfile, *outfile;
   char line1[LBUFSIZ];
   char line2[LBUFSIZ];
   char line3[LBUFSIZ]; /* Each field in the line */
   char *stptr;
   int flag = 0, counter=0, first=1;
   int idx = 0, xsize, kk;
   int lcount = 0; /* Cell Seperator */
   double ddata;
   double *X, *Xt;

   if (!(myfile = fopen(filename, "r")))
   {
      printf("read_csv ERROR: Could not open file for reading\n");
      return -1;
   }
   fgets(line1,sizeof line1,myfile);

   xsize = 1000;
   X = (double *) malloc(xsize * sizeof(double));
   /* Get a line from file */
   while (fgets(line1,sizeof line1,myfile) != NULL)
   { 
      lcount = 0;
      strcpy(line2,line1);
      stptr = line2;

      /* start going character by character thro the line */
      while (*stptr != '\0')
      { 
         lcount++;
         /* If field begins with " */
         if (*stptr == '"')
         {
            int flag = 0;
            idx = 0;
            while (flag == 0)
            {

               stptr++;
               /* Find corresponding closing " */
               while (*stptr != '"')
               { line3[idx] = *stptr;
                  idx++;
                  stptr++;
               }
               stptr++;
               if (*stptr != '\0' && *stptr == ',')
               {
                  line3[idx] = '\0';
                  flag = 1;
               }
               else if (*stptr != '\0' && *stptr == '"')
               { line3[idx] = *stptr;
                  idx++;
               }
               else
               {
                  line3[idx] = '\0';
                  flag = 1;
               }
            }
         }
         else
         { 
            idx = 0;
            while (*stptr != '\0' && *stptr != ',')
            {  
               line3[idx] = *stptr;
               idx++;
               stptr++;
            }
            line3[idx] = '\0';
            if (lcount >= 2 && lcount <=9)
            {
               sscanf(line3, "%lg", &ddata);
               X[counter++] = ddata;
               if (counter >= xsize)
               {
                  Xt = X;
                  xsize += 1000;
                  X = (double *) malloc(xsize*sizeof(double));
                  for (kk = 0; kk < counter; kk++) X[kk] = Xt[kk];
                  free(Xt);
               }
            }
            if (lcount >= 12 && lcount <=17)
            {
               sscanf(line3, "%lg", &ddata);
               X[counter++] = ddata;
               if (counter >= xsize)
               {
                  Xt = X;
                  xsize += 1000;
                  X = (double *) malloc(xsize*sizeof(double));
                  for (kk = 0; kk < counter; kk++) X[kk] = Xt[kk];
                  free(Xt);
               }
            }
         }
         if (*stptr != '\0' && *stptr == ',')
            stptr++;
         strcpy(line2,stptr);
         stptr = line2;
      }
      if (first == 1)
      {
         (*nInOuts) = counter;
         first = 0;
      }
   }
   fclose(myfile);
   (*XX) = X;
   (*nSamples) = counter / (*nInOuts);
   return 0;
}


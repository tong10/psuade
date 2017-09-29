/*
------------------------------------------------------------------------------
Support for R functions
------------------------------------------------------------------------------
*/
#include <math.h>

static double R_NaN=1.0/0.0;
static double R_NegInf=-HUGE;
static double R_PostInf=HUGE;

void Rf_warning(int id, char *str)
{
   printf("R warning: code = %d, message = %s\n",id,str);
}


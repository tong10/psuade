/* ******************************************************************
 * This subroutine tests to see if a given point (x,y) is inside the
 * convex hull formed by the n points given in xarray and yarray    *
 * ******************************************************************/
int checkvertical0(double x1, double y1, double x2, double y2, 
                   double x3, double y3, double x4, double y4);

int TestPtConvexHull(double x,double y,int n,double *xarray,
                     double *yarray)
{
  int    i, j, k, found=0, flag=0, good;
  double x1, y1, x2, y2, x3, y3, x4, y4, xint, yint, m1, m2, mi1, mi2;
  double xtmp, ytmp, yt1, yt2, yt3, yt4;

  /* form a line between (x,y) and a vertex point on the convex hull */

  for (i=0; i<n; i++) 
  {
    x1 = x; y1 = y;
    x2 = xarray[i]; y2 = yarray[i];
    found = 1;

    /* if it is a point, then (x,y) is on the boundary */
    /* if it is indeed a line, then proceed to test whether the line
       intersects with any line formed by any other two vertices in the
       convex hull */

    if (x1 != x2 || y1 != y2) 
    {
      /* re-orient the 2 points of line #1 so that the 2 points are 
         from left to right */
      if (x2 < x1) 
      {
        xtmp = x1; ytmp = y1;
        x1 = x2; y1 = y2;
        x2 = xtmp; y2 = ytmp;
      }

      /* now loop over all posible lines formed by any two other 
         vertices of the convex hull */

      found = 0;
      for (j=0; j<n; j++) 
      {
        for (k=j+1; k<n; k++) 
        {
          if ((j != i) && (k != i)) 
          {
            x3 = xarray[j]; y3 = yarray[j]; 
            x4 = xarray[k]; y4 = yarray[k]; 

            /* re-orient the 2 points of line #1 so that the 2 
               points are from left to right */
            if (x4 < x3) 
            {
              xtmp = x3; ytmp = y3;
              x3 = x4; y3 = y4;
              x4 = xtmp; y4 = ytmp;
            }

            /* now begin the search */
         
            if (x2 == x1) 
            {
              if (x3 != x4) 
              {
                /* if only line 1 is a vertical line */
                if (((x1>x3) && (x1<x4)) || ((x1>x4) && (x1<x3)))
                  found = 1;
              }
            } 
            else 
            {
              if (x3 == x4) 
              {
                /* if only line 2 is a vertical line */
                if (y1 == y2) 
                {
                  /* and if line 1 is a horizontal line */
                  ytmp = y1;
                } 
                else 
                { 
                  /* and if line 1 is not a horizontal line */
                  ytmp = (y2 - y1) / (x2 - x1) * (x3 - x1) + y1;
                }
                if (((ytmp>y3) && (ytmp<y4)) || ((ytmp>y4) && (ytmp<y3)))
                  found = 1;
              }
              else 
              {
                /* if both lines are not vertical */
                found = checkvertical0(x1,y1,x2,y2,x3,y3,x4,y4);
              }
            }
          }
          if (found == 1) break;
        } /* for k */
        if (found == 1) break;
      } /* for j */
    } /* if */ else return 2;
    if (found != 1) break;
  } /* for i */

  /* now if (x,y) is inside the convex hull, return 1 */
  if (found == 1) return 1;

  /*     if (x,y) is on a vertex, return 2            */
  if (found == 2) return 2;

  /* but if found is equal to 0, (x,y) may still lie on the */
  /* boundary, so further test is needed.                   */
  if (found == 0) 
  {
    for (j=0; j<n; j++) 
    {
      for (k=j+1; k<n; k++) 
      {
        x1 = xarray[j]; y1 = yarray[j]; 
        x2 = xarray[k]; y2 = yarray[k]; 
        if (x2 < x1) 
        {
          xtmp = x1; ytmp = y1;
          x1 = x2; y1 = y2;
          x2 = xtmp; y2 = ytmp;
        }
        if (x1 == x2) 
        {
          if (x == x1) 
          {
            if (((y >= y1) && (y <= y2)) || ((y >= y2) && (y <= y1))) 
            {
              found = 3; break;
            }
          }
        }
        else 
        {
          if (x != x1) 
          {
            m1 = (y - y1) / (x - x1);
            m2 = (y2 - y1) / (x2 - x1);
            if (m1 == m2) 
            {
              if (x >= x1 && x <= x2) 
              {
                if (((y >= y1) && (y <= y2)) || ((y >= y2) && (y <= y1))) 
                {
                  found = 3; break;
                }
              }
            }
          }
        } 
      } 
    } 
  }
  return found;
}

/* ******************************************************************
 * if there are no vertical lines, then use this subroutine to check 
   for intersection.  Assume x2 > x1, and x4 > x3.
 * *****************************************************************/
int checkvertical0(double x1, double y1, double x2, double y2, 
                   double x3, double y3, double x4, double y4) 
{
  int found=0;
  double m1, m2, mi1, mi2, yt1, yt2, yt3, yt4, ytmp, xtmp;

  /* compute the slope of line 1 and 2 */

  m1 = (y2 - y1) / (x2 - x1);
  m2 = (y4 - y3) / (x4 - x3);

  if (m1 != m2) 
  {
    /* if the two lines are not parallel */
    if (y1 == y2) 
    {
      /* if line 1 is horizontal */
      if (((y1>y3) && (y1<y4)) || ((y1>y4) && (y1<y3)))
        found = 1;
    } 
    else if (y3 == y4) 
    {
      /* if line 2 is horizontal */
      ytmp = (y3 - y1) * (x2 - x1) / (y2 - y1) + x1;
      if (((ytmp>y3) && (ytmp<y4)) || ((ytmp>y4) && (ytmp<y3)))
        found = 1;
    } 
    else 
    {
      xtmp = (y3 - y1 + m1 * x1 - m2 * x3) / (m1 - m2);
      mi1 = 1.0 / m1; mi2 = 1.0 / m2;
      ytmp = (x3 - x1 + mi1 * y1 - mi2 * y3) / (mi1 - mi2);
      if (xtmp > x3 && xtmp < x4)
        found = 1;
    }
  }
  return found;
}


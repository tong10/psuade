/* mars36_fort.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#ifndef NO_F2C

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static real c_b196 = (float)1.;
static integer c__0 = 0;
static doublereal c_b214 = -1.;
static doublereal c_b395 = 0.;
static doublereal c_b1099 = 1.;
static doublereal c_b1305 = 2147483647.;


/* Multivariate Adaptive Regression Splines (MARS modeling, version 3.6). */


/* Coded and copywrite (c) by Jerome H. Friedman (3/25/93). */



/*                         A Micro User's Guide */
/*                                  to */
/*                               MARS 3.6 */

/*                          Jerome H. Friedman */
/*                          Stanford University */

/*   MARS 3.6 is a collection of subroutines that implement the multivariate*/
/* adaptive regression spline strategy for data fitting and function */
/*approximation described in Friedman (1991a, 1991b, 1993). It is a general-*/
/*ization of MARS 1.0 described in Friedman (1988) and MARS 2.5 described in*/
/* Friedman (1991a). It implements as a special case a palindromically */
/*invariant version of the TURBO fitting technique for smoothing and additive
*/
/* modeling described in Friedman and Silverman (1989). */

/*   These subroutines represent a set of tools that can be invoked from a use
r*/
/*coded program to perform various analyses. The user routine is responsible f
or*/
/*reading the data into memory and passing it to the MARS subroutines, along*/
/*with the various parameter settings, as arguments. This set of subroutines*/
/* can also form the basis for incorporating this methodology into a */
/* statistical language or package. */

/*    The user interface subroutines are: */

/* MARS: computes the mars model from the data and provides summary */
/*     information for interpreting it. */

/* MISS: sets up mars to handle missing values. */

/* PLOT: constructs graphical output, useful for interpreting the */
/*     continuous part of the mars model, in a format that can be plotted */
/*     with a local graphics package. */

/*CATPRT: prints tables useful for interpreting the purely categorical parts*/
/*     of the mars model. */

/* SLICE: produces user selected lower dimensional representations of the */
/*    mars model to aid in interpretation. */

/* FMOD: computes response estimates for given predictor variable vectors */
/*    given the mars model. */

/*  It should be noted that this is a new methodology with which there is,*/
/* as of this time, little collective experience. */


/* References: */

/* [1] Friedman, J. H. (1988). Fitting functions to noisy data in high */
/*    dimensions. Proc., Twentyth Symposium on the Interface, Wegman, Gantz,*/
/*    and Miller, eds. American Statistical Association, Alexandria, VA. 3-43.
*/

/* [2] Friedman, J. H. (1991a). Multivariate adaptive regression splines */
/*     (with discussion).  Annals of Statistics, 19, 1-141 (March). */

/* [3] Friedman, J. H. (1991b). Estimating functions of mixed ordinal and */
/*    categorical variables using adaptive splines. Department of Statistics,
*/
/*     Stanford University, Tech. Report LCS108. */

/* [4] Friedman, J. H. (1993). Fast MARS. Department of Statistics, */
/*     Stanford University, Tech. Report LCS110. */

/* [5] Friedman, J. H. and Silverman, B. W. (1989). Flexible parsimonious */
/*    smoothing and additive modeling (with discussion). TECHNOMETRICS, 31,*/
/*     3-39 (Feburary). */



/* User interface subroutines: */

/*       All arguments in the calling sequence of user called MARS */
/*       subroutines must be of the same type in the calling program */
/*       as indicated in the subroutine's documentation below. All */
/*       calling sequence arrays must be dimensioned  in the calling */
/*       program as indicated in the documentation below. This includes */
/*       all workspace arrays. The leading dimensions of multidimensional */
/*       arrays must match, whereas the last dimension and singley */
/*       dimensioned arrays can be as large or larger than that indicated. */

/*       More detailed explanations for some of the quanities below can */
/*       be found in the indicated sections of references [2], [3], and/or */
/*       [4] above. */



/* call mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm): */

/* input: */
/* n = number of observations. */
/* p = number of predictor variables per observation (integer). */
/* x(n,p) = predictor variable data matrix. */
/* y(n) = response value for each observation. */
/* w(n) = weight (mass) for each observation. */
/*nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)*/
/*mi = maximum number of variables per basis function (interaction level).*/
/*      mi=1 => additive modeling (main effects only); */
/*      mi>1 => up to mi-variable interactions allowed. */
/*lx(p) = predictor variable flags; lx(i) corresponds to the ith variable:*/
/*    lx(i) =  0 : exclude variable from model. */
/*             1 : ordinal variable - no restriction. */
/*             2 : ordinal variable that can only enter additively; */
/*                 no interactions with other variables. */
/*             3 : ordinal variable that can enter only linearly. */
/*            -1 : categorical variable - no restriction. */
/*            -2 : categorical variable that can only enter additively; */
/*                 no interactions with other variables. */

/* output: */
/* fm(3+nk*(5*mi+nmcv+6)+2*p+ntcv), im(21+nk*(3*mi+8)) = mars model. */
/*   (nmcv = maximum number of distinct values for any categorical variable;*/
/*    ntcv = total number of distinct values over all categorical variables.)
*/

/*   note: upon return im(1) and im(2) contain the lengths of the fm and im*/
/*          arrays (respectively) actually used by the program. */

/* workspace: */
/* sp(n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk) : real. */
/*dp(max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)) : double precision.*/
/* mm(n*p+2*max(mi,nmcv)) : integer. */


/* defaults: */
/*   the following quanities are set to default values in mars. their values*/
/*   can be changed by executing any of the following statements before the*/
/*    call to mars, with the new value as the argument. */

/* call speed(is): */
/* is = speed acceleration factor (1-5). */
/*    larger values progressively sacrifice optimization thuroughness for */
/*   computational speed advantage. this usually results in marked decrease*/
/*   in computing time with little or no effect on resulting approximation*/
/*    accuracy (especially useful for exploratory work). */
/*    is = 1 => no acceleration. */
/*    is = 5 => maximum speed advantage. */
/*    (default: is=4) (Ref [4] Secs. 3.0 - 4.0) */

/* call logit(il): */
/* il = ordinary/logistic regression flag. */
/*    il=0 => ordinary least-squares regression. */
/*    il=1 => logistic regression. */
/*    (default: il=0). If logistic regression is selected (il=1) then */
/*    the response variable is assumed to take on only the values 0/1 and */
/*   the mars model is for the log-odds: f(x) = log (Pr(Y=1:x)/Pr(Y=0:x)).*/
/*    (Ref[2] Sec. 4.5) */

/* call setdf(df): */
/* df = number of degrees-of-freedom charged for (unrestricted) */
/*      knot optimization. (default: df=3.0) */
/*      (Ref[2] Sec. 3.6, Ref[3] Sec. 2.3) */

/* call xvalid(ix): */
/* ix = control parameter for sample reuse technique used to automatically */
/*    estimate smoothing parameter df (see above) from the data. */
/*    ix = 0 => no effect (default). value used for df is set by user if */
/*            setdf(df) is called (see above), otherwise default value */
/*            (df=3.0) is used. */
/*    ix > 0 => ix - fold cross-validation. */
/*   ix < 0 => single validation pass using every (-ix)th (randomly selected)
*/
/*            observation as an independent test set. */
/*    if ix.ne.0 then call setdf(df) (see above) has no effect. if ix > 0, */
/*   computation increases roughly by a factor of ix over that for ix = 0.*/
/*    for ix < 0 computation increases approximately by a factor of two. */
/*    (Ref[3] Sec. 2.3) */

/* call stseed(is): */
/*is = seed for internal random number generator used to group observation*/
/*      subsets for validation (ix.ne.0). (default: is=987654321). */

/* call print(it): */
/*it = fortran file number for printed output. (it.le.0 => no printed output.)
*/
/*      note that this controls printed output for all user called mars */
/*      routines. (default: it=6). */

/* call setfv(fv): */
/*fv = (fractional) incremental penalty for increasing the number of variables
*/
/*     in the mars model. sometimes useful with highly collinear designs as it
*/
/*     may produce nearly equivalent models with fewer predictor variables,*/
/*      aiding in interpretation. (fv .ge. 0) */
/*    fv=0.0  => no penalty (default). */
/*    fv=0.05 => moderate penalty. */
/*    fv=0.1  => heavy penality. */
/*    the best value depends on the specific situation and some user */
/*   experimentation using different values is usually required. this option*/
/*    should be used with some care. (Ref[2] Sec. 5.3) */

/* call setic(ic): */
/* ic = flag restricting categorical - ordinal interactions. */
/*    ic=0 => no effect (default). */
/*   ic=1 => interactions between categorical and ordinal variables prohibited
.*/
/*    ic=2 => maximum number of ordinal variables participating in any */
/*           interaction is restricted to two. categorical interactions are*/
/*            unrestricted. */
/*   the restrictions associated with a value of ic are imposed in addition*/
/*    to those that are controlled by the mi and lx flags (see above). */

/* call setint(i,j,k): */
/* i,j = predictor variable numbers. */
/*   k = interaction flag: */
/*      k.eq.0 => interactions between variables i and j are prohibited. */
/*      k.ne.0 => interactions between variables i and j are permitted */
/*               if allowed by the mi, lx, and ic parameter values (see above)
.*/
/*                (default) */
/*   i.eq.0 .or. j.eq.0 => reset to defaults for all predictor variables. */

/*   this call can be executed repeatedly to invoke (or remove) multiple */
/*   constraints. */

/* call nest (n,i,j,nv,vals); */
/* nests predictor variable i within categorical predictor variable j. */
/*links the existance of a value for var(i) to a subset of values for var(j).
*/
/*variable j must be categorical (lx(j) < 0, see above). (Ref[3] Sec. 3.3)*/

/* n = same as in mars (see above). */
/* i = variable to be nested (ordinal or categorical). */
/* j = variable to which i is nested (categorical). */
/* nv = size of corresponding subset of values for variable j. */
/* vals(nv) = specific values of variable j for which values of variable i */
/*            are defined to exist. */

/*setting nv = 0 removes the nesting associated with i and j previously set.*/
/* setting i = 0 and/or j = 0 removes all nesting (default). */

/* this call can be executed repeatedly to invoke (or remove) nestings. */
/* (recursive nesting is not implemented. j may not itself be nested to */
/* another variable.) */

/*note: variable nesting does NOT override the interaction constraints set by
*/
/*      calling setic or setint (see above). the mi and lx flags (see above)*/
/*      do not limit interactions on variables to which others are nested.*/
/*      they control nested and other nonnested variables as described above*/
/*      except that interactions with variables to which others are nested*/
/*       are ignored in applying the constraints. */


/* call setms(ms): */
/* ms = minimum span (minimum number of observations between each knot). */
/*      ms .le. 0 => default value (depending on n and p) is used. */
/*      (default: ms=0). (Ref[2] Sec. 3.8) */


/* the following three routines (miss, mkmiss, xmiss) are used to enable */
/* mars to deal with various aspects of missing predictor values in the */
/* training data and/or future data to be predicted. for problems with */
/* no such missing values, these three routines are of no concern. */


/* call miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp); */

/*called  (optionally) before mars (see above) to indicate the presence of*/
/* missing values in the predictor variable data matrix. */

/* sets up mars to accomodate missing values. */
/* produces as output transformations of some of the original mars input */
/* quanities defining the problem. the new transformed quanities define */
/*the same problem - but with missing values - and replace the corresponding*/
/* original quanities in the call to mars (see above) and plot and slice */
/* (see below). in particular, a new predictor data matrix is created */
/* with extra (dummy) variables - one for each original variable with */
/* missing values - to indicate observations for which the corresponding */
/* original variable is not missing. each such original variable is */
/*automatically nested (see above) to this (corresponding) dummy variable.*/
/* also produces as output a slicing vector to be used as input to */
/* slice (see below) to produce the mars model corresponding to no missing */
/* values on any of the original predictor variables, for interpretation. */
/* (Ref[3] Sec. 3.4) */

/* input: */
/* n,p,x,lx = original intended input to mars (see above). */
/* xm(p) = vector giving missing value flag for each original variable. */
/*         x(i,j) = xm(j) => x(i,j) is missing. */
/* flg = input to slice (see below). */

/* output: */
/*pn,xn,lxn = corresponding transformed quanities to be used as input to mars.
*/
/* pn = number of predictor variables in transformed data matrix (integer, */
/*      .le. 2*p) */
/* xn(n,pn) = transformed data matrix. */
/* lxn(pn) = predictor variable flags for transformed data matrix. */
/*xs(pn) = input for slice (see below). slicing vector that produces a slice*/
/*         of the augmented predictor variable space to produce the nonmissing
*/
/*          mars model. */
/* xp(2*p+1) = input for xmiss (see below). */

/* notes: */
/*    (1) the value of the output quanity pn is less than or equal to 2*p. */
/*        the other output quanities (arrays) should be dimensioned */
/*        xn(n,2*p), lxn(2*p), xs(2*p) in the calling program to be safe. */
/*    (2) if there are no missing values then the output quanities */
/*        pn,xn,lxn will be identical to p,x,lx respectively, and */
/*        xs(j)=flg, j=1,p. */
/*   (3) the corresponding quanities for input and output can be the same in*/
/*        the calling program. however, they should have the corresponding */
/*        dimension equal to 2*p, and they will be altered. */
/*    (4) dimensions of all relevant workspace arrays in mars and other */
/*        user callable routines must be large enough to accomodate the */
/*        the increased number of variables pn (.le. 2*p) produced. */


/* call mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc); */

/* called (optionally) before miss (see above) to generate additional data */
/* with missing values by resampling from original (input) data. */

/* used to train mars for missing predictor values when they are under */
/* represented in the training data. takes as input original data and */
/* produces as output a new (larger) data set containing both the original */
/*data and additional data resampled from it with predictor variable values*/
/* flagged as missing. (Ref[3] Sec. 3.4) */

/* input: */
/* n,p,x,y,w = original data set (see mars above). */
/* xm(p) = same as in miss (see above). */
/* pm(p) = vector of fractions of missing values in output sample: */
/*   pm(j) = fraction of the total output data with jth predictor variable*/
/*            value missing. */
/* nnx = maximum sample size of new (output) data. */

/* output: */
/* nn = sample size of generated output data set (input to miss, mars, and */
/*      other user called routines, in place of n). */
/* xnn(nn,p) = new generated predictor data matrix to be used as input to */
/*             miss (see above) in place of x. */
/*yn(nn),wn(nn) = new generated data to be used as input to mars and other*/
/*                 user called routines in place of y and w. */

/* workspace: */
/* sc(p*nnx): real. */

/* notes: */
/*   (1) the output arrays should be dimensioned xnn(nnx,p),yn(nnx),wn(nnx)*/
/*        in the calling program to be safe. */
/*    (2) if much larger fractions of missing values are requested than */
/*        exist in the orginal training data then the size of the output */
/*        (resampled) data set can become very large (nn = big). */
/*    (3) if the size of the output (resampled) data set reaches the */
/*        specified maximum (nnx) then the actual fractions of missing */
/*       values for each variable will be less than those specified in the*/
/*        input vector pm(p). */
/*    (4) dimensions of all relevant workspace arrays in mars and other */
/*        user callable routines must be large enough to accomodate the */
/*        the increased number of observations nn (.le. nnx) produced. */


/* call xmiss (n,x,xm,xp,xn); */

/* must be called before fmod (see below) if miss (see above) was called */
/* before mars (see above) to handle missing values. produces a new set */
/*of covariate vectors, from the original set intended for fmod, to replace*/
/* the original set in the call to fmod. (Ref[3] Sec. 3.4) */

/* input: */
/* n = number of covariate vectors (see fmod below). */
/* x(n,p) = original covariate vectors intended for fmod. */
/* xm, xp = same as in miss (see above). */

/* output: */
/* xn(n,pn) = new set of covariate vectors to be used as input to fmod. */

/* notes: */
/*   (1) the value of the quanity pn is less than or equal to 2*p (p = the*/
/*       number of original predictor variables). the output quanity (array)*/
/*       should be dimensioned xn(n,2*p) in the calling program to be safe.*/
/*    (2) the corresponding quanities x and xn can be the same in */
/*        the calling program. however, it should have the second */
/*        dimension equal to 2*p, and it will be altered. */




/* The following subroutines can be called only after mars. */



/* call plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm): */

/* computes plots of all purely additive and all purely bivariate ordinal */
/* contributions to the mars model,in a form suitable for displaying */
/* with a computer graphics package. If there are no interactions between */
/* categorical and ordinal variables present in the mars model, this */
/*subroutine returns plots of the (compound) anova functions (Ref[2] Sec. 3.5,
*/
/* eqns (25) and (27)) for those ordinal variables involved in at most two */
/*(ordinal) variable interactions. If categorical-ordinal interactions are*/
/*present, it returns the corresponding plots for each ordinal contribution*/
/* in the categorical-ordinal decomposition (Ref[3] Sec. 3.2).  since the */
/*locations of the plotted functions are arbitrary, they are all translated*/
/* to have zero minimum value. */

/* input: */
/*   m = model flag: */
/*     = 1 => plot piecewise-linear mars model. */
/*     = 2 => plot piecewise-cubic mars model. (Ref[2] Sec. 3.7) */
/*   x,fm,im = same as in mars (see above). */
/*   ngc = number of raster points for computing curve estimates. */
/*  ngs = number of raster points on each axis for computing surface estimates
.*/
/*   icx = convex hull flag: */
/*      = 0 => plot surface estimates over entire range of argument limits.*/
/*       > 0 => plot surface estimates only inside the convex hull of the */
/*             (bivariate) point set. */

/* output: */
/*    nc = number of curves (purely additive ordinal anova functions). */
/*    crv(ngc,2,nc) = additive ordinal anova functions in anova */
/*                    decomposition order. */
/*       crv(.,1,m) = ordered abscissa values for mth anova function. */
/*       crv(.,2,m) = corresponding ordinate values. */
/*    ns = number of surfaces (purely bivariate ordinal anova functions). */
/*    srf(ngs,ngs,ns) = bivariate (plus associated univariate) ordinal */
/*                      anova functions in anova decomposition order. */
/*      srf(.,.,m) = total contribution (bivariate + univariate) of ordinal*/
/*                 variables associated with mth bivariate anova function of t
he*/
/*                 mars model, in a square raster format (ngs x ngs) over the
*/
/*                 ranges of the two variables. the first and second indicies
*/
/*                 correspond to the first and second variables respectively.
*/

/* workspace: */
/*    sp(max(4*ngs*ngs,ngc,2*n)) : real. */
/*    mm(max(2*(mi+1),nmcv)) : integer. */



/* call catprt (m,fm,im,sp,mm): */

/*prints all univariate and bivariate contributions to the purely categorical
*/
/*part of the mars model that are not involved in higher order (categorical)*/
/*interactions. These are the (compound) anova functions (Ref[2] Sec. 3.5,*/
/* eqns (25) and (27)) for the purely categorical part of the categorical- */
/* ordinal decomposition (Ref[3] Sec. 3.2, eqns (32b) and (38b)). since */
/* the locations of the printed functions are arbitrary they are all */
/*translated to have zero minimum value. the functions are printed as tables*/
/* of integers in the range [0,99] or [0,9]. the function value for each */
/* entry is the product of the corresponding integer and the scale factor */
/* printed above the table. */

/* input: */
/* m = model flag: */
/*   = 1 => piecewise-linear model. */
/*   = 2 => piecewise-cubic  model. (Ref[2] Sec. 3.7) */
/* fm,im = same as in mars (see above). */

/* workspace: */
/* sp(nmcv*nmcv) : real. */
/* mm(2*nk+nmcv) : integer. */



/* call slice (flg,xs,x,fm,im,fmn,imn,sp,mm): */

/*computes the mars model within a lower dimensional subspace defined by an*/
/*axis oriented slice of the predictor variable space. the slice is selected*/
/*by assigning specific values to a subset of the predictor variables. the*/
/*returned model is a function of the variables complement to the selected*/
/*subset, and represents the mars model conditioned on the specified values*/
/* for the selected variables. this new lower dimensional sliced model can */
/*be input to plot and/or catprt for interpretation in the same manner as the
*/
/* original mars model. (Ref[2] Sec. 4.7, Ref[3] Sec. 2.4) */

/* input: */
/*flg = flag for indicating that a predictor variable is not in the subset*/
/*   defining the slice. its value should be outside the range of values for*/
/*    all predictor variables. */
/* xs(p) = vector defining the slice: */
/*    xs(i).eq.flg => do not condition on ith variable. */
/*    xs(i).ne.flg => condition on ith variable at value stored in xs(i). */
/* x = same as in mars (see above). */
/* fm,im = arrays defining mars model (output from mars - see above). */

/* output: */
/* fmn,imn = corresponding arrays defining the sliced model */
/*             (dimensioned same as in mars - see above). */

/* workspace: */
/* sp(2*nk+2*p+max(nk,3*p)) : real. */
/* mm(2*p) : integer. */



/* call fmod (m,n,x,fm,im,f,sp): */

/* calculates mars model response estimates for sets of covariate vectors. */

/* input: */
/* m = model flag: */
/*   = 1 => piecewise-linear mars model. */
/*   = 2 => piecewise-cubic mars model. (Ref[2] Sec. 3.7) */
/* n = number of covariate vectors. */
/* x(n,p) = covariate vectors. */
/* fm,im = same as in mars (see above). */

/* output: */
/* f(n) = value of the mars model estimate for each covariate vector. */

/* workspace: */
/* sp(n,2) : real. */



/* call cvinfo (dfs,pse,nbf): */

/*returns results of sample reuse procedure for estimating optimal smoothing*/
/* parameter df (see above). can only be called if xvalid(ix) was called */
/* before mars with ix.ne.0 (see above). (Ref[2] Sec 3.6, Ref[3] Sec. 2.3) */

/* output: */
/* dfs = optimal smoothing parameter estimate. */
/* pse = estimate of corresponding predictive-squared-error. */
/* nbf = estimate of associated number of (nonconstant) basis functions. */



/*<       subroutine mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm)                >*/
/* Subroutine */ int mars_(n, p, x, y, w, nk, mi, lx, fm, im, sp, dp, mm)
integer *n, *p;
real *x, *y, *w;
integer *nk, *mi, *lx;
real *fm;
integer *im;
real *sp;
doublereal *dp;
integer *mm;
{
    extern /* Subroutine */ int mars1_();
    extern integer lcm_();

/*<       integer p,lx(*),im(*),mm(*)                                        >*/
/*<       real x(*),y(*),w(*),fm(*),sp(*)                                    >*/
/*<       double precision dp(*)                                             >*/
/*<       im(3)=n                                                            >*/
    /* Parameter adjustments */
    --mm;
    --dp;
    --sp;
    --im;
    --fm;
    --lx;
    --w;
    --y;
    --x;

    /* Function Body */
    im[3] = *n;
/*<       im(4)=p                                                            >*/
    im[4] = *p;
/*<       im(5)=nk                                                           >*/
    im[5] = *nk;
/*<       im(6)=mi                                                           >*/
    im[6] = *mi;
/*<       im(7)=16                                                           >*/
    im[7] = 16;
/*<       im(8)=im(7)+5*nk                                                   >*/
    im[8] = im[7] + *nk * 5;
/*<       im(9)=im(8)+2*nk*mi                                                >*/
    im[9] = im[8] + (*nk << 1) * *mi;
/*<       im(10)=im(9)+3*(nk+2)                                              >*/
    im[10] = im[9] + (*nk + 2) * 3;
/*<       im(2)=im(10)+nk*mi-1                                               >*/
    im[2] = im[10] + *nk * *mi - 1;
/*<       im(11)=1                                                           >*/
    im[11] = 1;
/*<       im(12)=2                                                           >*/
    im[12] = 2;
/*<       im(13)=im(12)+5*nk                                                 >*/
    im[13] = im[12] + *nk * 5;
/*<       im(14)=im(13)+1                                                    >*/
    im[14] = im[13] + 1;
/*<       im(15)=im(14)+nk*(5*mi+1)                                          >*/
    im[15] = im[14] + *nk * (*mi * 5 + 1);
/*<    >*/
    mars1_(n, p, &x[1], &y[1], &w[1], nk, mi, &lx[1], &fm[im[11]], &fm[im[12]]
        , &fm[im[15]], &im[im[7]], &im[im[8]], &im[im[9]], &im[im[10]], &
        fm[im[13]], &fm[im[14]], &sp[1], &dp[1], &mm[1]);
/*<       im(1)=im(15)+lcm(p,nk,fm(im(12)),fm(im(15)))-1                     >*/
    im[1] = im[15] + lcm_(p, nk, &fm[im[12]], &fm[im[15]]) - 1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* mars_ */

/*<       subroutine plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm)        >*/
/* Subroutine */ int plot_(m, x, fm, im, ngc, ngs, icx, nc, crv, ns, srf, sp, 
    mm)
integer *m;
real *x, *fm;
integer *im, *ngc, *ngs, *icx, *nc;
real *crv;
integer *ns;
real *srf, *sp;
integer *mm;
{
    extern /* Subroutine */ int plotc_(), plotl_();

/*<       integer im(*),mm(*)                                                >*/
/*<       real x(*),fm(*),crv(*),srf(*),sp(*)                                >*/
/*<       if(m .ne. 1) go to 1                                               >*/
    /* Parameter adjustments */
    --mm;
    --sp;
    --srf;
    --crv;
    --im;
    --fm;
    --x;

    /* Function Body */
    if (*m != 1) {
    goto L1;
    }
/*<    >*/
    plotl_(&im[3], &im[4], &x[1], &im[5], &im[im[7]], &im[im[8]], &im[im[9]], 
        &im[im[10]], &fm[im[12]], &fm[im[15]], ngc, ngs, icx, nc, &crv[1],
         ns, &srf[1], &sp[1], &mm[1]);
/*<       return                                                             >*/
    return 0;
/*<    >*/
L1:
    plotc_(&im[3], &im[4], &x[1], &im[5], &im[im[7]], &im[im[8]], &im[im[9]], 
        &im[im[10]], &fm[im[14]], &fm[im[15]], ngc, ngs, icx, nc, &crv[1],
         ns, &srf[1], &sp[1], &mm[1]);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* plot_ */

/*<       subroutine catprt (m,fm,im,sp,mm)                                  >*/
/* Subroutine */ int catprt_(m, fm, im, sp, mm)
integer *m;
real *fm;
integer *im;
real *sp;
integer *mm;
{
    extern /* Subroutine */ int ctprt1_();

/*<       integer im(*),mm(*)                                                >*/
/*<       real fm(*),sp(*)                                                   >*/
/*<    >*/
    /* Parameter adjustments */
    --mm;
    --sp;
    --im;
    --fm;

    /* Function Body */
    ctprt1_(m, &im[5], &im[im[7]], &im[im[8]], &fm[im[12]], &fm[im[15]], &fm[
        im[14]], &sp[1], &mm[1]);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* catprt_ */

/*<       subroutine slice (flg,xs,x,fm,im,fmn,imn,sp,mm)                    >*/
/* Subroutine */ int slice_(flg, xs, x, fm, im, fmn, imn, sp, mm)
real *flg, *xs, *x, *fm;
integer *im;
real *fmn;
integer *imn;
real *sp;
integer *mm;
{
    static integer i__;
    extern /* Subroutine */ int slice1_();

/*<       integer im(*),imn(*),mm(*)                                         >*/
/*<       real xs(*),x(*),fm(*),fmn(*),sp(*)                                 >*/
/*<       do 1 i=1,15                                                        >*/
    /* Parameter adjustments */
    --mm;
    --sp;
    --imn;
    --fmn;
    --im;
    --fm;
    --x;
    --xs;

    /* Function Body */
    for (i__ = 1; i__ <= 15; ++i__) {
/*<       imn(i)=im(i)                                                       >*/
    imn[i__] = im[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       i=im(15)                                                           >*/
    i__ = im[15];
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 i=i+1                                                              >*/
L2:
    ++i__;
/*<     3 if((i).gt.(im(1))) go to 4                                         >*/
L3:
    if (i__ > im[1]) {
    goto L4;
    }
/*<       fmn(i)=fm(i)                                                       >*/
    fmn[i__] = fm[i__];
/*<       go to 2                                                            >*/
    goto L2;
/*<    >*/
L4:
    slice1_(flg, &xs[1], &im[3], &im[4], &x[1], &im[5], &fm[im[11]], &fm[im[
        12]], &fm[im[15]], &im[im[7]], &im[im[8]], &im[im[9]], &im[im[10]]
        , &fm[im[13]], &fm[im[14]], &fmn[im[11]], &fmn[im[12]], &imn[im[7]
        ], &imn[im[8]], &imn[im[9]], &imn[im[10]], &fmn[im[13]], &fmn[im[
        14]], &sp[1], &mm[1]);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* slice_ */

/*<       subroutine fmod (m,n,x,fm,im,f,sp)                                 >*/
/* Subroutine */ int fmod_(m, n, x, fm, im, f, sp)
integer *m, *n;
real *x, *fm;
integer *im;
real *f, *sp;
{
    extern /* Subroutine */ int cmrs_(), fmrs_();

/*<       integer im(*)                                                      >*/
/*<       real x(*),fm(*),f(*),sp(*)                                         >*/
/*<       if(m .ne. 1) go to 1                                               >*/
    /* Parameter adjustments */
    --sp;
    --f;
    --im;
    --fm;
    --x;

    /* Function Body */
    if (*m != 1) {
    goto L1;
    }
/*<       call fmrs(n,x,im(5),fm(im(11)),fm(im(12)),fm(im(15)),f)            >*/
    fmrs_(n, &x[1], &im[5], &fm[im[11]], &fm[im[12]], &fm[im[15]], &f[1]);
/*<       return                                                             >*/
    return 0;
/*<    >*/
L1:
    cmrs_(n, &x[1], &fm[im[15]], &im[im[7]], &im[im[8]], &im[im[9]], &im[im[
        10]], &fm[im[13]], &fm[im[14]], &f[1], &sp[1]);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* fmod_ */

/*<       subroutine print(it)                                               >*/
/* Subroutine */ int print_(it)
integer *it;
{
    extern /* Subroutine */ int printc_(), printg_(), prtslc_(), printm_();

/*<       call printm(it)                                                    >*/
    printm_(it);
/*<       call printg(it)                                                    >*/
    printg_(it);
/*<       call printc(it)                                                    >*/
    printc_(it);
/*<       call prtslc(it)                                                    >*/
    prtslc_(it);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* print_ */

/*<       subroutine setint(i,j,k)                                           >*/
/* Subroutine */ int setint_0_(n__, i__, j, k, it)
int n__;
integer *i__, *j, *k, *it;
{
    /* Initialized data */

    static integer il = 0;

    /* Builtin functions */
    /* Subroutine */ int s_stop();
    integer s_wsfe(), e_wsfe();

    /* Local variables */
    static integer l, m[2000]   /* was [2][1000] */, m1, m2, ig, ll;

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 0, 0, "(/,' interactions prohibited between:\
')", 0 };


/*<       parameter(mlist=1000)                                              >*/
/*<       integer m(2,mlist)                                                 >*/
/*<       save m                                                             >*/
/*<       data il /0/                                                        >*/
    switch(n__) {
    case 1: goto L_intlst;
    case 2: goto L_intalw;
    }

/*<       if((i .ne. 0) .and. (j .ne. 0)) go to 1                            >*/
    if (*i__ != 0 && *j != 0) {
    goto L1;
    }
/*<       il=0                                                               >*/
    il = 0;
/*<       return                                                             >*/
    return 0;
/*<     1 if(i.eq.j) return                                                  >*/
L1:
    if (*i__ == *j) {
    return 0;
    }
/*<       m1=min0(i,j)                                                       >*/
    m1 = min(*i__,*j);
/*<       m2=max0(i,j)                                                       >*/
    m2 = max(*i__,*j);
/*<       if(k .ne. 0) go to 6                                               >*/
    if (*k != 0) {
    goto L6;
    }
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 l=l+1                                                              >*/
L2:
    ++l;
/*<     3 if((l).gt.(il)) go to 4                                            >*/
L3:
    if (l > il) {
    goto L4;
    }
/*<       if(m1.eq.m(1,l).and.m2.eq.m(2,l)) return                           >*/
    if (m1 == m[(l << 1) - 2] && m2 == m[(l << 1) - 1]) {
    return 0;
    }
/*<       go to 2                                                            >*/
    goto L2;
/*<     4 il=il+1                                                            >*/
L4:
    ++il;
/*<       if(il .le. mlist) go to 5                                          >*/
    if (il <= 1000) {
    goto L5;
    }
/*    write(6,  '('' increase parameter mlist in subroutine setint to gr  
 98*/
/*   1eater than'',           i5,/,'' and recompile.'')') il              
 99*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     5 m(1,il)=m1                                                         >*/
L5:
    m[(il << 1) - 2] = m1;
/*<       m(2,il)=m2                                                         >*/
    m[(il << 1) - 1] = m2;
/*<       return                                                             >*/
    return 0;
/*<     6 ig=0                                                               >*/
L6:
    ig = 0;
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 l=l+1                                                              >*/
L7:
    ++l;
/*<     8 if((l).gt.(il)) go to 10                                           >*/
L8:
    if (l > il) {
    goto L10;
    }
/*<       if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 7                     >*/
    if (m1 != m[(l << 1) - 2] || m2 != m[(l << 1) - 1]) {
    goto L7;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<    10 if(ig.eq.0) return                                                 >*/
L10:
    if (ig == 0) {
    return 0;
    }
/*<       il=il-1                                                            >*/
    --il;
/*<       ll=l                                                               >*/
    ll = l;
/*<       go to 12                                                           >*/
    goto L12;
/*<    11 ll=ll+1                                                            >*/
L11:
    ++ll;
/*<    12 if((ll).gt.(il)) go to 13                                          >*/
L12:
    if (ll > il) {
    goto L13;
    }
/*<       m(1,ll)=m(1,ll+1)                                                  >*/
    m[(ll << 1) - 2] = m[(ll + 1 << 1) - 2];
/*<       m(2,ll)=m(2,ll+1)                                                  >*/
    m[(ll << 1) - 1] = m[(ll + 1 << 1) - 1];
/*<       go to 11                                                           >*/
    goto L11;
/*<    13 return                                                             >*/
L13:
    return 0;
/*<       entry intlst(it)                                                   >*/

L_intlst:
/*<       if(it.le.0) return                                                 >*/
    if (*it <= 0) {
    return 0;
    }
/*<       if(il.eq.0) return                                                 >*/
    if (il == 0) {
    return 0;
    }
/*<       write(it,'(/,'' interactions prohibited between:'')')              >*/
    io___9.ciunit = *it;
    s_wsfe(&io___9);
    e_wsfe();
/*    do 14 l=1,il                                                        
125*/
/*    write(it,'(''    var('',i3,'')  and  var('',i3,'')'')') m(1,l),m(2  
126*/
/*   1,l)                                                                 
127*/
/* 14 continue                                                            
128*/
/*<       return                                                             >*/
    return 0;
/*<       entry intalw(i,j,k)                                                >*/

L_intalw:
/*<       k=1                                                                >*/
    *k = 1;
/*<       m1=min0(i,j)                                                       >*/
    m1 = min(*i__,*j);
/*<       m2=max0(i,j)                                                       >*/
    m2 = max(*i__,*j);
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 l=l+1                                                              >*/
L15:
    ++l;
/*<    16 if((l).gt.(il)) go to 18                                           >*/
L16:
    if (l > il) {
    goto L18;
    }
/*<       if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 15                    >*/
    if (m1 != m[(l << 1) - 2] || m2 != m[(l << 1) - 1]) {
    goto L15;
    }
/*<       k=0                                                                >*/
    *k = 0;
/*<    18 return                                                             >*/
L18:
    return 0;
/*<       end                                                                >*/
} /* setint_ */

/* Subroutine */ int setint_(i__, j, k)
integer *i__, *j, *k;
{
    return setint_0_(0, i__, j, k, (integer *)0);
    }

/* Subroutine */ int intlst_(it)
integer *it;
{
    return setint_0_(1, (integer *)0, (integer *)0, (integer *)0, it);
    }

/* Subroutine */ int intalw_(i__, j, k)
integer *i__, *j, *k;
{
    return setint_0_(2, i__, j, k, (integer *)0);
    }

/*<    >*/
/* Subroutine */ int mars1_0_(n__, n, p, x, y, w, nk, mi, lx, az, tb, cm, kp, 
    kv, lp, lv, bz, tc, sp, dp, mm, mal, val)
int n__;
integer *n, *p;
real *x, *y, *w;
integer *nk, *mi, *lx;
real *az, *tb, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *sp;
doublereal *dp;
integer *mm, *mal;
real *val;
{
    /* Initialized data */

    static integer ms = 0;
    static real df = (float)3.;
    static integer il = 0;
    static real fv = (float)0.;
    static integer it = 6;
    static integer ic = 0;
    static integer ix = 0;

    /* System generated locals */
    integer mm_dim1, mm_offset, x_dim1, x_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    extern /* Subroutine */ int cmrs_(), fmrs_();
    static integer i__, j, k;
    static real s, t;
    extern /* Subroutine */ int cubic_(), ccoll_(), anova_(), catpr_(), 
        orgpc_();
    extern doublereal stelg_();
    extern /* Subroutine */ int orgpl_(), ordpr_();
    static integer i1, i2;
    extern /* Subroutine */ int psort_();
    static real ef;
    static integer im, is;
    static real wn, sw, yh;
    extern /* Subroutine */ int anoval_(), logitc_(), atoscl_(), sclato_(), 
        marsgo_(), logitl_(), cvmars_(), varimp_(), oknest_(), intlst_(), 
        xvmrgo_(), rspnpr_(), nstlst_();
    static real gcv, z00001;

/*<       integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),mm(n,*),lx(p)              >*/
/*<       real x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),sp(*)                   >*/
/*<       double precision dp(*)                                             >*/
/*<       data ms,df,il,fv,it,ic,ix /0,3.0,0,0.0,6,0,0/                      >*/
    /* Parameter adjustments */
    if (y) {
    --y;
    }
    if (w) {
    --w;
    }
    if (mm) {
    mm_dim1 = *n;
    mm_offset = mm_dim1 + 1;
    mm -= mm_offset;
    }
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (lx) {
    --lx;
    }
    if (tb) {
    tb -= 6;
    }
    if (cm) {
    --cm;
    }
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (lp) {
    lp -= 4;
    }
    if (lv) {
    --lv;
    }
    if (tc) {
    --tc;
    }
    if (sp) {
    --sp;
    }
    if (dp) {
    --dp;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_setms;
    case 2: goto L_setdf;
    case 3: goto L_printm;
    case 4: goto L_logit;
    case 5: goto L_setfv;
    case 6: goto L_setic;
    case 7: goto L_xvalid;
    }

/*    if(it.gt.0) write(it,11)                                            
148*/
/*    if(it.gt.0) write(it,10) n,p,nk,ms,mi,df,il,fv,ic                   
149*/
/*    if(it.gt.0) write(it,12)                                            
150*/
/*    if(it.gt.0) write(it,'('' var: '',5('' '',20i3,/))') (i,i=1,p)      
151*/
/*    if(it.gt.0) write(it,'('' flag:'',5('' '',20i3,/))') (lx(i),i=1,p)  
152*/
/*     print *, ' ' */
/*     do 321 i = 1, n */
/*        print *,'M1 ',x(i,1),' ',x(i,2),' ',x(i,3),' ',x(i,4), */
/*    1         ' ',y(i) */
/* 321  continue */
/*<       call intlst(it)                                                    >*/
    intlst_(&it);
/*<       call nstlst(it)                                                    >*/
    nstlst_(&it);
/*<       i1=max0(n*(nk+1),2*n)+1                                            >*/
/* Computing MAX */
    i__1 = *n * (*nk + 1), i__2 = *n << 1;
    i1 = max(i__1,i__2) + 1;
/*<       im=i1+n+max0(3*n+5*nk,2*p,4*n,2*n+5*nk+p)                          >*/
/* Computing MAX */
    i__1 = *n * 3 + *nk * 5, i__2 = *p << 1, i__1 = max(i__1,i__2), i__2 = *n 
        << 2, i__1 = max(i__1,i__2), i__2 = (*n << 1) + *nk * 5 + *p;
    im = i1 + *n + max(i__1,i__2);
/*<       is=im+p                                                            >*/
    is = im + *p;
/*<       i2=max0(n*nk,(nk+1)*(nk+1))+1                                      >*/
/* Computing MAX */
    i__1 = *n * *nk, i__2 = (*nk + 1) * (*nk + 1);
    i2 = max(i__1,i__2) + 1;
/*<       call rspnpr(it,il,n,y,w,mm)                                        >*/
    rspnpr_(&it, &il, n, &y[1], &w[1], &mm[mm_offset]);
/*<       do 2 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       do 1 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       mm(i,j)=i                                                          >*/
        mm[i__ + j * mm_dim1] = i__;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       call psort(x(1,j),mm(1,j),1,n)                                     >*/
    psort_(&x[j * x_dim1 + 1], &mm[j * mm_dim1 + 1], &c__1, n);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       call ordpr(it,n,p,x,lx,mm)                                         >*/
    ordpr_(&it, n, p, &x[x_offset], &lx[1], &mm[mm_offset]);
/*<       call atoscl (n,p,w,x,lx,mm,sp(im),sp(is),cm,x)                     >*/
    atoscl_(n, p, &w[1], &x[x_offset], &lx[1], &mm[mm_offset], &sp[im], &sp[
        is], &cm[1], &x[x_offset]);
/*<       call catpr(it,n,p,x,cm,mm(1,p+1))                                  >*/
    catpr_(&it, n, p, &x[x_offset], &cm[1], &mm[(*p + 1) * mm_dim1 + 1]);
/*<       call oknest(it,p,lx,cm)                                            >*/
    oknest_(&it, p, &lx[1], &cm[1]);
/*<    >*/
    if (ix != 0) {
    cvmars_(&ix, n, p, &x[x_offset], &y[1], &w[1], nk, &ms, &df, &fv, mi, 
        &lx[1], &it, &sp[im], &sp[is], &tb[6], &cm[1], &sp[1], &dp[1],
         &dp[i2], &mm[mm_offset], &sp[is + *p], &sp[is + *p + (*n << 
        1)]);
    }
/*<    >*/
    marsgo_(n, p, &x[x_offset], &y[1], &w[1], nk, &ms, &df, &fv, mi, &lx[1], &
        it, &sp[im], &sp[is], az, &tb[6], &cm[1], &sp[1], &dp[1], &dp[i2],
         &mm[mm_offset]);
/*<       if(il .le. 0) go to 6                                              >*/
    if (il <= 0) {
    goto L6;
    }
/*<       call logitl(n,x,y,w,nk,il,az,tb,cm,sp,dp)                          >*/
    logitl_(n, &x[x_offset], &y[1], &w[1], nk, &il, az, &tb[6], &cm[1], &sp[1]
        , &dp[1]);
/*<       if(it .le. 0) go to 6                                              >*/
    if (it <= 0) {
    goto L6;
    }
/*<       sw=0.0                                                             >*/
    sw = (float)0.;
/*<       wn=sw                                                              >*/
    wn = sw;
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    r__1 = sw;
    wn = r__1 * r__1 / wn;
/*<       ef=1.0                                                             >*/
    ef = (float)1.;
/*<       do 4 k=1,nk                                                        >*/
    i__1 = *nk;
    for (k = 1; k <= i__1; ++k) {
/*<       if(tb(1,k).ne.0.0) ef=ef+tb(5,k)                                   >*/
    if (tb[k * 5 + 1] != (float)0.) {
        ef += tb[k * 5 + 5];
    }
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       ef=1.0/(1.0-ef/wn)**2                                              >*/
/* Computing 2nd power */
    r__1 = (float)1. - ef / wn;
    ef = (float)1. / (r__1 * r__1);
/*<       s=0.0                                                              >*/
    s = (float)0.;
/*<       t=s                                                                >*/
    t = s;
/*<       call fmrs(n,x,nk,az,tb,cm,sp)                                      >*/
    fmrs_(n, &x[x_offset], nk, az, &tb[6], &cm[1], &sp[1]);
/*<       do 5 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       yh=1.0/(1.0+exp(-sp(i)))                                           >*/
    yh = (float)1. / (exp(-sp[i__]) + (float)1.);
/*<       gcv=ef*(y(i)-yh)**2                                                >*/
/* Computing 2nd power */
    r__1 = y[i__] - yh;
    gcv = ef * (r__1 * r__1);
/*<       s=s+w(i)*gcv                                                       >*/
    s += w[i__] * gcv;
/*<       t=t+w(i)*yh*(1.0-yh)                                               >*/
    t += w[i__] * yh * ((float)1. - yh);
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       t=t/sw                                                             >*/
    t /= sw;
/*    write(it,13) s,t                                                    
200*/
/*<     6 if(it .le. 0) go to 7                                              >*/
L6:
    if (it <= 0) {
    goto L7;
    }
/*<       if(il.eq.0) call anova (n,x,y,w,nk,it,tb,cm,lp,lv,sp,dp)           >*/
    if (il == 0) {
    anova_(n, &x[x_offset], &y[1], &w[1], nk, &it, &tb[6], &cm[1], &lp[4],
         &lv[1], &sp[1], &dp[1]);
    }
/*<       if(il.gt.0) call anoval(n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sp,dp)     >*/
    if (il > 0) {
    anoval_(n, &x[x_offset], &y[1], &w[1], nk, &il, &it, az, &tb[6], &cm[
        1], &lp[4], &lv[1], &sp[1], &dp[1]);
    }
/*<     7 call ccoll (nk,tb,cm,kp,kv,lp,lv,mm)                               >*/
L7:
    ccoll_(nk, &tb[6], &cm[1], &kp[6], &kv[3], &lp[4], &lv[1], &mm[mm_offset])
        ;
/*<    >*/
    cubic_(n, p, &x[x_offset], &y[1], &w[1], nk, &it, &tb[6], &cm[1], &kp[6], 
        &kv[3], &lp[4], &lv[1], bz, &tc[1], &sp[1], &sp[i1], &sp[i1 + (*p 
        << 1)], &mm[mm_offset], &dp[1]);
/*<       if(il .le. 0) go to 9                                              >*/
    if (il <= 0) {
    goto L9;
    }
/*<       call logitc(n,x,y,w,nk,il,cm,kp,kv,lp,lv,bz,tc,sp,sp(i1+4*n),dp)   >*/
    logitc_(n, &x[x_offset], &y[1], &w[1], nk, &il, &cm[1], &kp[6], &kv[3], &
        lp[4], &lv[1], bz, &tc[1], &sp[1], &sp[i1 + (*n << 2)], &dp[1]);
/*<       if(it .le. 0) go to 9                                              >*/
    if (it <= 0) {
    goto L9;
    }
/*<       call cmrs(n,x,cm,kp,kv,lp,lv,bz,tc,sp,sp(n+1))                     >*/
    cmrs_(n, &x[x_offset], &cm[1], &kp[6], &kv[3], &lp[4], &lv[1], bz, &tc[1],
         &sp[1], &sp[*n + 1]);
/*<       s=0.0                                                              >*/
    s = (float)0.;
/*<       t=s                                                                >*/
    t = s;
/*<       do 8 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       yh=1.0/(1.0+exp(-sp(i)))                                           >*/
    yh = (float)1. / (exp(-sp[i__]) + (float)1.);
/*<       gcv=ef*(y(i)-yh)**2                                                >*/
/* Computing 2nd power */
    r__1 = y[i__] - yh;
    gcv = ef * (r__1 * r__1);
/*<       s=s+w(i)*gcv                                                       >*/
    s += w[i__] * gcv;
/*<       t=t+w(i)*yh*(1.0-yh)                                               >*/
    t += w[i__] * yh * ((float)1. - yh);
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       t=t/sw                                                             >*/
    t /= sw;
/*    write(it,14) s,t                                                    
221*/
/*<    >*/
L9:
    if (it > 0) {
    varimp_(n, p, &x[x_offset], &y[1], &w[1], nk, &il, &it, az, &tb[6], &
        cm[1], &sp[1], &sp[*p + 1], &dp[1]);
    }
/*<       call orgpl(sp(im),sp(is),nk,tb,cm)                                 >*/
    orgpl_(&sp[im], &sp[is], nk, &tb[6], &cm[1]);
/*<       call orgpc(sp(im),sp(is),lp,lv,tc)                                 >*/
    orgpc_(&sp[im], &sp[is], &lp[4], &lv[1], &tc[1]);
/*<       call sclato(n,p,x,sp(im),sp(is),cm,x)                              >*/
    sclato_(n, p, &x[x_offset], &sp[im], &sp[is], &cm[1], &x[x_offset]);
/*<       return                                                             >*/
    return 0;
/*<       entry setms(mal)                                                   >*/

L_setms:
/*<       ms=mal                                                             >*/
    ms = *mal;
/*<       return                                                             >*/
    return 0;
/*<       entry setdf(val)                                                   >*/

L_setdf:
/*<       df=val                                                             >*/
    df = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry printm(mal)                                                  >*/

L_printm:
/*<       it=mal                                                             >*/
    it = *mal;
/*<       return                                                             >*/
    return 0;
/*<       entry logit(mal)                                                   >*/

L_logit:
/*<       il=mal                                                             >*/
    il = *mal;
/*<       return                                                             >*/
    return 0;
/*<       entry setfv(val)                                                   >*/

L_setfv:
/*<       fv=val                                                             >*/
    fv = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry setic(mal)                                                   >*/

L_setic:
/*<       ic=mal                                                             >*/
    ic = *mal;
/*<       z00001=stelg(ic)                                                   >*/
    z00001 = stelg_(&ic);
/*<       return                                                             >*/
    return 0;
/*<       entry xvalid(mal)                                                  >*/

L_xvalid:
/*<       ix=mal                                                             >*/
    ix = *mal;
/*<       call xvmrgo(ix)                                                    >*/
    xvmrgo_(&ix);
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L10: */
/*<    11 format(//,' MARS modeling, version 3.6 (3/25/93)',/)               >*/
/* L11: */
/*<    12 format(/' predictor variable flags:')                              >*/
/* L12: */
/*<    >*/
/* L13: */
/*<    >*/
/* L14: */
/*<       end                                                                >*/
} /* mars1_ */

/* Subroutine */ int mars1_(n, p, x, y, w, nk, mi, lx, az, tb, cm, kp, kv, lp,
     lv, bz, tc, sp, dp, mm)
integer *n, *p;
real *x, *y, *w;
integer *nk, *mi, *lx;
real *az, *tb, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *sp;
doublereal *dp;
integer *mm;
{
    return mars1_0_(0, n, p, x, y, w, nk, mi, lx, az, tb, cm, kp, kv, lp, lv, 
        bz, tc, sp, dp, mm, (integer *)0, (real *)0);
    }

/* Subroutine */ int setms_(mal)
integer *mal;
{
    return mars1_0_(1, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, mal, (real *)0);
    }

/* Subroutine */ int setdf_(val)
real *val;
{
    return mars1_0_(2, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, (integer *)0, val);
    }

/* Subroutine */ int printm_(mal)
integer *mal;
{
    return mars1_0_(3, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, mal, (real *)0);
    }

/* Subroutine */ int logit_(mal)
integer *mal;
{
    return mars1_0_(4, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, mal, (real *)0);
    }

/* Subroutine */ int setfv_(val)
real *val;
{
    return mars1_0_(5, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, (integer *)0, val);
    }

/* Subroutine */ int setic_(mal)
integer *mal;
{
    return mars1_0_(6, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, mal, (real *)0);
    }

/* Subroutine */ int xvalid_(mal)
integer *mal;
{
    return mars1_0_(7, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        integer *)0, mal, (real *)0);
    }

/*<    >*/
/* Subroutine */ int plotc_0_(n__, n, p, x, nk, kp, kv, lp, lv, tc, cm, ngc, 
    ngs, icx, nc, crv, ns, srf, sp, mm, tb, nal)
int n__;
integer *n, *p;
real *x;
integer *nk, *kp, *kv, *lp, *lv;
real *tc, *cm;
integer *ngc, *ngs, *icx, *nc;
real *crv;
integer *ns;
real *srf, *sp;
integer *mm;
real *tb;
integer *nal;
{
    /* Initialized data */

    static real big = (float)1e30;
    static integer it = 6;

    /* System generated locals */
    integer x_dim1, x_offset, crv_dim1, crv_offset, srf_dim1, srf_dim2, 
        srf_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int cfun_(), pair_();
    static integer ngsq;
    static real d__;
    static integer i__, j, k, l, m;
    static real r__;
    extern /* Subroutine */ int cpair_();
    static real d1, d2;
    static integer k1, k2, l1, k4, l2;
    static real dc, dl;
    static integer ne, nf, jj, ll, ko;
    static real fx;
    static integer jv, iz;
    static real zl[2];
    static integer nh;
    static real zu[2];
    extern /* Subroutine */ int hulset_(), cvxhul_();
    static integer ncx;
    extern /* Subroutine */ int fun_();
    static integer jnt, nxs;

/*<       integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),mm(*)                      >*/
/*<    >*/
/*<       data big,it /1.e30,6/                                              >*/
    /* Parameter adjustments */
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (tb) {
    tb -= 6;
    }
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (lp) {
    lp -= 4;
    }
    if (lv) {
    --lv;
    }
    if (tc) {
    --tc;
    }
    if (cm) {
    --cm;
    }
    if (crv) {
    crv_dim1 = *ngc;
    crv_offset = crv_dim1 * 3 + 1;
    crv -= crv_offset;
    }
    if (srf) {
    srf_dim1 = *ngs;
    srf_dim2 = *ngs;
    srf_offset = srf_dim1 * (srf_dim2 + 1) + 1;
    srf -= srf_offset;
    }
    if (sp) {
    --sp;
    }
    if (mm) {
    --mm;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_plotl;
    case 2: goto L_printg;
    }

/*    if(it.gt.0) write(it,'(/'' mars graphics (piecewise-cubic):'',/)')  
267*/
/*<       jnt=2                                                              >*/
    jnt = 2;
/*<       go to 1                                                            >*/
    goto L1;
/*<    >*/

L_plotl:
/*    if(it.gt.0) write(it,'(/'' mars graphics (piecewise-linear):'',/)'  
272*/
/*   1)                                                                   
273*/
/*<       jnt=1                                                              >*/
    jnt = 1;
/*<     1 ngsq=ngs**2                                                        >*/
L1:
/* Computing 2nd power */
    i__1 = *ngs;
    ngsq = i__1 * i__1;
/*<       iz=2*ngsq                                                          >*/
    iz = ngsq << 1;
/*<       d=1.0/(ngs-1)                                                      >*/
    d__ = (float)1. / (*ngs - 1);
/*<       dc=1.0/(ngc-1)                                                     >*/
    dc = (float)1. / (*ngc - 1);
/*<       ll=1                                                               >*/
    ll = 1;
/*<       nc=0                                                               >*/
    *nc = 0;
/*<       ns=nc                                                              >*/
    *ns = *nc;
/*<     2 if(kp(1,ll).lt.0) go to 36                                         >*/
L2:
    if (kp[ll * 5 + 1] < 0) {
    goto L36;
    }
/*<       if(kp(3,ll) .gt. 0) go to 3                                        >*/
    if (kp[ll * 5 + 3] > 0) {
    goto L3;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     3 nf=kp(3,ll)                                                        >*/
L3:
    nf = kp[ll * 5 + 3];
/*<       k4=kp(4,ll)-1                                                      >*/
    k4 = kp[ll * 5 + 4] - 1;
/*<       k1=kp(1,ll)                                                        >*/
    k1 = kp[ll * 5 + 1];
/*<       k2=kp(2,ll)                                                        >*/
    k2 = kp[ll * 5 + 2];
/*<       if(it .le. 0) go to 7                                              >*/
    if (it <= 0) {
    goto L7;
    }
/*<       if(k1 .ne. 0) go to 4                                              >*/
    if (k1 != 0) {
    goto L4;
    }
/*    write(it,'('' pure ordinal contribution:'')')                       
292*/
/*<       go to 7                                                            >*/
    goto L7;
/*<     4 continue >*/
L4:
/*    write(it,'('' categorical - ordinal interaction:'')')               
294*/
/*<       do 6 i=1,k1                                                        >*/
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       jj=kv(1,k2+i-1)                                                    >*/
    jj = kv[(k2 + i__ - 1 << 1) + 1];
/*<       j=iabs(jj)                                                         >*/
    j = abs(jj);
/*<       k=kv(2,k2+i-1)                                                     >*/
    k = kv[(k2 + i__ - 1 << 1) + 2];
/*<       ncx=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                            >*/
    ncx = (integer) (cm[(j << 1) + 1] + (float).1) - (integer) (cm[j * 2] 
        + (float).1) + 1;
/*<       do 5 l=1,ncx                                                       >*/
    i__2 = ncx;
    for (l = 1; l <= i__2; ++l) {
/*<       mm(l)=cm(k+l)+.1                                                   >*/
        mm[l] = cm[k + l] + (float).1;
/*<       if(jj.lt.0) mm(l)=mod(mm(l)+1,2)                                   >*/
        if (jj < 0) {
        mm[l] = (mm[l] + 1) % 2;
        }
/*<     5 continue                                                           >*/
/* L5: */
    }
/*    write(it,'('' x('',i3,'') ='',70i1/80i1)') j,(mm(l),l=1,ncx)    
    304*/
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<     7 do 35 k=1,nf                                                       >*/
L7:
    i__1 = nf;
    for (k = 1; k <= i__1; ++k) {
/*<       l=lp(1,k+k4)                                                       >*/
    l = lp[(k + k4) * 3 + 1];
/*<       if(l.gt.2) go to 35                                                >*/
    if (l > 2) {
        goto L35;
    }
/*<       ko=lp(2,k+k4)                                                      >*/
    ko = lp[(k + k4) * 3 + 2];
/*<       if(l .ne. 1) go to 17                                              >*/
    if (l != 1) {
        goto L17;
    }
/*<       j=0                                                                >*/
    j = 0;
/*<       jv=lv(ko)                                                          >*/
    jv = lv[ko];
/*<       do 9 m=k,nf                                                        >*/
    i__2 = nf;
    for (m = k; m <= i__2; ++m) {
/*<       l1=lp(1,m+k4)                                                      >*/
        l1 = lp[(m + k4) * 3 + 1];
/*<       if(l1.eq.1) go to 9                                                >*/
        if (l1 == 1) {
        goto L9;
        }
/*<       l2=lp(2,m+k4)-1                                                    >*/
        l2 = lp[(m + k4) * 3 + 2] - 1;
/*<       do 8 i=1,l1                                                        >*/
        i__3 = l1;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       if(jv.eq.lv(l2+i)) j=1                                             >*/
        if (jv == lv[l2 + i__]) {
            j = 1;
        }
/*<     8 continue                                                           >*/
/* L8: */
        }
/*<       if(j.eq.1) go to 10                                                >*/
        if (j == 1) {
        goto L10;
        }
/*<     9 continue                                                           >*/
L9:
        ;
    }
/*<    10 if(j.eq.1) go to 35                                                >*/
L10:
    if (j == 1) {
        goto L35;
    }
/*<       nc=nc+1                                                            >*/
    ++(*nc);
/*<       zl(1)=big                                                          >*/
    zl[0] = big;
/*<       zu(1)=-big                                                         >*/
    zu[0] = -big;
/*<       do 11 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       r=x(i,jv)                                                          >*/
        r__ = x[i__ + jv * x_dim1];
/*<       zl(1)=amin1(zl(1),r)                                               >*/
        zl[0] = dmin(zl[0],r__);
/*<       zu(1)=amax1(zu(1),r)                                               >*/
        zu[0] = dmax(zu[0],r__);
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       dl=(zu(1)-zl(1))*dc                                                >*/
    dl = (zu[0] - zl[0]) * dc;
/*<       do 12 i=1,ngc                                                      >*/
    i__2 = *ngc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       crv(i,1,nc)=zl(1)+dl*(i-1)                                         >*/
        crv[i__ + ((*nc << 1) + 1) * crv_dim1] = zl[0] + dl * (i__ - 1);
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<       if(jnt .ne. 1) go to 13                                            >*/
    if (jnt != 1) {
        goto L13;
    }
/*<       call fun(l,jv,ngc,crv(1,1,nc),nk,tb,cm,k1,kv(1,k2),crv(1,2,nc),mm) >*/
    fun_(&l, &jv, ngc, &crv[((*nc << 1) + 1) * crv_dim1 + 1], nk, &tb[6], 
        &cm[1], &k1, &kv[(k2 << 1) + 1], &crv[((*nc << 1) + 2) * 
        crv_dim1 + 1], &mm[1]);
/*<       go to 14                                                           >*/
    goto L14;
/*<    >*/
L13:
    cfun_(&l, &jv, ngc, &crv[((*nc << 1) + 1) * crv_dim1 + 1], &nf, &lp[(
        k4 + 1) * 3 + 1], &lv[1], &tc[kp[ll * 5 + 5]], &crv[((*nc << 
        1) + 2) * crv_dim1 + 1], &sp[1], &mm[1]);
/*<    14 dl=big                                                             >*/
L14:
    dl = big;
/*<       do 15 i=1,ngc                                                      >*/
    i__2 = *ngc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       dl=amin1(dl,crv(i,2,nc))                                           >*/
/* Computing MIN */
        r__1 = dl, r__2 = crv[i__ + ((*nc << 1) + 2) * crv_dim1];
        dl = dmin(r__1,r__2);
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<       fx=0.0                                                             >*/
    fx = (float)0.;
/*<       do 16 i=1,ngc                                                      >*/
    i__2 = *ngc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       crv(i,2,nc)=crv(i,2,nc)-dl                                         >*/
        crv[i__ + ((*nc << 1) + 2) * crv_dim1] -= dl;
/*<       fx=amax1(fx,crv(i,2,nc))                                           >*/
/* Computing MAX */
        r__1 = fx, r__2 = crv[i__ + ((*nc << 1) + 2) * crv_dim1];
        fx = dmax(r__1,r__2);
/*<    16 continue                                                           >*/
/* L16: */
    }
/*    if(it.gt.0) write(it,39) nc,jv,fx                               
    349*/
/*<       go to 35                                                           >*/
    goto L35;
/*<    17 j=0                                                                >*/
L17:
    j = 0;
/*<       mm(1)=lv(ko)                                                       >*/
    mm[1] = lv[ko];
/*<       mm(2)=lv(ko+1)                                                     >*/
    mm[2] = lv[ko + 1];
/*<       do 19 m=k,nf                                                       >*/
    i__2 = nf;
    for (m = k; m <= i__2; ++m) {
/*<       l1=lp(1,m+k4)                                                      >*/
        l1 = lp[(m + k4) * 3 + 1];
/*<       if(l1.le.2) go to 19                                               >*/
        if (l1 <= 2) {
        goto L19;
        }
/*<       l2=lp(2,m+k4)-1                                                    >*/
        l2 = lp[(m + k4) * 3 + 2] - 1;
/*<       do 18 i=1,l1                                                       >*/
        i__3 = l1;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       if(mm(1).eq.lv(l2+i).or.mm(2).eq.lv(l2+i)) j=1                     >*/
        if (mm[1] == lv[l2 + i__] || mm[2] == lv[l2 + i__]) {
            j = 1;
        }
/*<    18 continue                                                           >*/
/* L18: */
        }
/*<       if(j.eq.1) go to 20                                                >*/
        if (j == 1) {
        goto L20;
        }
/*<    19 continue                                                           >*/
L19:
        ;
    }
/*<    20 if(j.eq.1) go to 35                                                >*/
L20:
    if (j == 1) {
        goto L35;
    }
/*<       ns=ns+1                                                            >*/
    ++(*ns);
/*<       zl(1)=big                                                          >*/
    zl[0] = big;
/*<       zl(2)=zl(1)                                                        >*/
    zl[1] = zl[0];
/*<       zu(1)=-big                                                         >*/
    zu[0] = -big;
/*<       zu(2)=zu(1)                                                        >*/
    zu[1] = zu[0];
/*<       do 22 j=1,2                                                        >*/
    for (j = 1; j <= 2; ++j) {
/*<       do 21 i=1,n                                                        >*/
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
/*<       r=x(i,mm(j))                                                       >*/
        r__ = x[i__ + mm[j] * x_dim1];
/*<       zl(j)=amin1(zl(j),r)                                               >*/
/* Computing MIN */
        r__1 = zl[j - 1];
        zl[j - 1] = dmin(r__1,r__);
/*<       zu(j)=amax1(zu(j),r)                                               >*/
/* Computing MAX */
        r__1 = zu[j - 1];
        zu[j - 1] = dmax(r__1,r__);
/*<    21 continue                                                           >*/
/* L21: */
        }
/*<    22 continue                                                           >*/
/* L22: */
    }
/*<       do 23 j=1,2                                                        >*/
    for (j = 1; j <= 2; ++j) {
/*<       dl=(zu(j)-zl(j))/(ngs-3)                                           >*/
        dl = (zu[j - 1] - zl[j - 1]) / (*ngs - 3);
/*<       zu(j)=zu(j)+dl                                                     >*/
        zu[j - 1] += dl;
/*<       zl(j)=zl(j)-dl                                                     >*/
        zl[j - 1] -= dl;
/*<    23 continue                                                           >*/
/* L23: */
    }
/*<       ne=0                                                               >*/
    ne = 0;
/*<       d1=d*(zu(1)-zl(1))                                                 >*/
    d1 = d__ * (zu[0] - zl[0]);
/*<       d2=d*(zu(2)-zl(2))                                                 >*/
    d2 = d__ * (zu[1] - zl[1]);
/*<       do 25 j=1,ngs                                                      >*/
    i__2 = *ngs;
    for (j = 1; j <= i__2; ++j) {
/*<       do 24 i=1,ngs                                                      >*/
        i__3 = *ngs;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       ne=ne+1                                                            >*/
        ++ne;
/*<       sp(iz+ne)=zl(1)+d1*(i-1)                                           >*/
        sp[iz + ne] = zl[0] + d1 * (i__ - 1);
/*<       sp(iz+ngsq+ne)=zl(2)+d2*(j-1)                                      >*/
        sp[iz + ngsq + ne] = zl[1] + d2 * (j - 1);
/*<    24 continue                                                           >*/
/* L24: */
        }
/*<    25 continue                                                           >*/
/* L25: */
    }
/*<       dl=big                                                             >*/
    dl = big;
/*<       if(jnt .ne. 1) go to 26                                            >*/
    if (jnt != 1) {
        goto L26;
    }
/*<    >*/
    pair_(&mm[1], &ngsq, &sp[iz + 1], nk, &tb[6], &cm[1], &k1, &kv[(k2 << 
        1) + 1], &srf[(*ns * srf_dim2 + 1) * srf_dim1 + 1], &sp[1], &
        mm[3]);
/*<       go to 27                                                           >*/
    goto L27;
/*<    >*/
L26:
    cpair_(&mm[1], &ngsq, &sp[iz + 1], &nf, &lp[(k4 + 1) * 3 + 1], &lv[1],
         &tc[kp[ll * 5 + 5]], &srf[(*ns * srf_dim2 + 1) * srf_dim1 + 
        1], &sp[1]);
/*<    27 if(icx .le. 0) go to 29                                            >*/
L27:
    if (*icx <= 0) {
        goto L29;
    }
/*<       call cvxhul(n,x(1,mm(1)),x(1,mm(2)),big,nh,sp)                     >*/
    cvxhul_(n, &x[mm[1] * x_dim1 + 1], &x[mm[2] * x_dim1 + 1], &big, &nh, 
        &sp[1]);
/*<       if(it .le. 0 .or. 3*nh .lt. iz) go to 28                           >*/
    if (it <= 0 || nh * 3 < iz) {
        goto L28;
    }
/*<       nxs=sqrt(float(3*nh)*0.5)+1.1                                      >*/
    nxs = sqrt((real) (nh * 3) * (float).5) + (float)1.1;
/*    write(it,38) nxs                                                
    402*/
/*<    28 call hulset(ngsq,sp(iz+1),big,nh,sp,srf(1,1,ns))                   >*/
L28:
    hulset_(&ngsq, &sp[iz + 1], &big, &nh, &sp[1], &srf[(*ns * srf_dim2 + 
        1) * srf_dim1 + 1]);
/*<    29 do 31 j=1,ngs                                                      >*/
L29:
    i__2 = *ngs;
    for (j = 1; j <= i__2; ++j) {
/*<       do 30 i=1,ngs                                                      >*/
        i__3 = *ngs;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<    >*/
        if (i__ == 1 || j == 1 || i__ == *ngs || j == *ngs || srf[i__ 
            + (j + *ns * srf_dim2) * srf_dim1] >= big) {
            goto L30;
        }
/*<       dl=amin1(dl,srf(i,j,ns))                                           >*/
/* Computing MIN */
        r__1 = dl, r__2 = srf[i__ + (j + *ns * srf_dim2) * srf_dim1];
        dl = dmin(r__1,r__2);
/*<    30 continue                                                           >*/
L30:
        ;
        }
/*<    31 continue                                                           >*/
/* L31: */
    }
/*<       fx=0.0                                                             >*/
    fx = (float)0.;
/*<       do 34 j=1,ngs                                                      >*/
    i__2 = *ngs;
    for (j = 1; j <= i__2; ++j) {
/*<       do 33 i=1,ngs                                                      >*/
        i__3 = *ngs;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<    >*/
        if (i__ != 1 && (j != 1 && (i__ != *ngs && (j != *ngs && srf[
            i__ + (j + *ns * srf_dim2) * srf_dim1] < big)))) {
            goto L32;
        }
/*<       srf(i,j,ns)=0.0                                                    >*/
        srf[i__ + (j + *ns * srf_dim2) * srf_dim1] = (float)0.;
/*<       go to 33                                                           >*/
        goto L33;
/*<    32 srf(i,j,ns)=srf(i,j,ns)-dl                                         >*/
L32:
        srf[i__ + (j + *ns * srf_dim2) * srf_dim1] -= dl;
/*<       fx=amax1(fx,srf(i,j,ns))                                           >*/
/* Computing MAX */
        r__1 = fx, r__2 = srf[i__ + (j + *ns * srf_dim2) * srf_dim1];
        fx = dmax(r__1,r__2);
/*<    33 continue                                                           >*/
L33:
        ;
        }
/*<    34 continue                                                           >*/
/* L34: */
    }
/*    if(it.gt.0) write(it,40) ns,mm(1),mm(2),fx                      
    422*/
/*<    35 continue                                                           >*/
L35:
    ;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<    36 continue >*/
L36:
/*    if(it.gt.0) write(it,37) nc,ns                                      
426*/
/*<       return                                                             >*/
    return 0;
/*<       entry printg(nal)                                                  >*/

L_printg:
/*<       it=nal                                                             >*/
    it = *nal;
/*<       return                                                             >*/
    return 0;
/*<    37 format(/,' ',i3,' curves and',i3,' surfaces.'/)                    >*/
/* L37: */
/*<    38 format(' plot: convex hull too large. increase ngs to',i6)         >*/
/* L38: */
/*<    39 format('   crv',i3,':  x(',i2,').  max =',g12.4)                   >*/
/* L39: */
/*<    40 format('   srf',i3,':  x(',i2,'), x(',i2,').  max =',g12.4)        >*/
/* L40: */
/*<       end                                                                >*/
} /* plotc_ */

/* Subroutine */ int plotc_(n, p, x, nk, kp, kv, lp, lv, tc, cm, ngc, ngs, 
    icx, nc, crv, ns, srf, sp, mm)
integer *n, *p;
real *x;
integer *nk, *kp, *kv, *lp, *lv;
real *tc, *cm;
integer *ngc, *ngs, *icx, *nc;
real *crv;
integer *ns;
real *srf, *sp;
integer *mm;
{
    return plotc_0_(0, n, p, x, nk, kp, kv, lp, lv, tc, cm, ngc, ngs, icx, nc,
         crv, ns, srf, sp, mm, (real *)0, (integer *)0);
    }

/* Subroutine */ int plotl_(n, p, x, nk, kp, kv, lp, lv, tb, cm, ngc, ngs, 
    icx, nc, crv, ns, srf, sp, mm)
integer *n, *p;
real *x;
integer *nk, *kp, *kv, *lp, *lv;
real *tb, *cm;
integer *ngc, *ngs, *icx, *nc;
real *crv;
integer *ns;
real *srf, *sp;
integer *mm;
{
    return plotc_0_(1, n, p, x, nk, kp, kv, lp, lv, (real *)0, cm, ngc, ngs, 
        icx, nc, crv, ns, srf, sp, mm, tb, (integer *)0);
    }

/* Subroutine */ int printg_(nal)
integer *nal;
{
    return plotc_0_(2, (integer *)0, (integer *)0, (real *)0, (integer *)0, (
        integer *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, 
        (real *)0, (integer *)0, (integer *)0, (integer *)0, (integer *)0,
         (real *)0, (integer *)0, (real *)0, (real *)0, (integer *)0, (
        real *)0, nal);
    }

/*<       subroutine ctprt1 (m,nk,kp,kv,tb,cm,tc,sc,js)                      >*/
/* Subroutine */ int ctprt1_0_(n__, m, nk, kp, kv, tb, cm, tc, sc, js, nal)
int n__;
integer *m, *nk, *kp, *kv;
real *tb, *cm, *tc, *sc;
integer *js, *nal;
{
    /* Initialized data */

    static real big = (float)9.9e30;
    static integer it = 6;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    extern integer ncat_();
    extern /* Subroutine */ int catv_();
    extern doublereal cvlv_();
    static integer i__, j, k;
    static real x;
    static integer j1, j2, n1, n2;
    static real s1;
    static integer ja, jb, na, nc, nb, jj, nl, nv;
    static real xm, px, rx, xx, rxp;

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 0, 0, "(/,' there are',i3,' purely categori\
cal basis functions.')", 0 };


/*<       integer kp(5,*),kv(2,*),js(*)                                      >*/
/*<       real cm(*),tc(*),sc(*)                                             >*/
/*<       data big,it /9.9e30,6/                                             >*/
    /* Parameter adjustments */
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (cm) {
    --cm;
    }
    if (tc) {
    --tc;
    }
    if (sc) {
    --sc;
    }
    if (js) {
    --js;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_printc;
    }

/*<       if(it.le.0) return                                                 >*/
    if (it <= 0) {
    return 0;
    }
/*<       nc=ncat(kp)                                                        >*/
    nc = ncat_(&kp[6]);
/*<       if(nc.eq.0) return                                                 >*/
    if (nc == 0) {
    return 0;
    }
/*<    >*/
    io___68.ciunit = it;
    s_wsfe(&io___68);
    do_fio(&c__1, (char *)&nc, (ftnlen)sizeof(integer));
    e_wsfe();
/*    write(it,'('' purely additive and bivariate contributions follow''  
445*/
/*   1)')                                                                 
446*/
/*<       if(m .ne. 1) go to 1                                               >*/
    if (*m != 1) {
    goto L1;
    }
/*    write(it,'('' (piecewise-linear fit):'')')                          
448*/
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 continue >*/
L1:
/*    write(it,'('' (piecewise-cubic fit):'')')                           
450*/
/*<     2 call catv(1,kp,kv,nv,js)                                           >*/
L2:
    catv_(&c__1, &kp[6], &kv[3], &nv, &js[1]);
/*<       do 8 jj=1,nv                                                       >*/
    i__1 = nv;
    for (jj = 1; jj <= i__1; ++jj) {
/*<       j=js(jj)                                                           >*/
    j = js[jj];
/*<       xm=big                                                             >*/
    xm = big;
/*<       xx=-big                                                            >*/
    xx = -big;
/*<       nl=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                             >*/
    nl = (integer) (cm[(j << 1) + 1] + (float).1) - (integer) (cm[j * 2] 
        + (float).1) + 1;
/*<       do 3 i=1,nl                                                        >*/
    i__2 = nl;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       sc(i)=cvlv(m,1,j,i,nk,kp,kv,tb,cm,tc)                              >*/
        sc[i__] = cvlv_(m, &c__1, &j, &i__, nk, &kp[6], &kv[3], tb, &cm[1]
            , &tc[1]);
/*<       xm=amin1(xm,sc(i))                                                 >*/
/* Computing MIN */
        r__1 = xm, r__2 = sc[i__];
        xm = dmin(r__1,r__2);
/*<       xx=amax1(xx,sc(i))                                                 >*/
/* Computing MAX */
        r__1 = xx, r__2 = sc[i__];
        xx = dmax(r__1,r__2);
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       px=99.0                                                            >*/
    px = (float)99.;
/*<       if(nl.gt.26) px=9.0                                                >*/
    if (nl > 26) {
        px = (float)9.;
    }
/*<       rx=xx-xm                                                           >*/
    rx = xx - xm;
/*<       rxp=rx/px                                                          >*/
    rxp = rx / px;
/*    write(it,'(/,'' f( x('',i3,'') ) : scale ='',g12.4)') j,rxp     
    466*/
/*<       if(rxp.le.0.0) go to 8                                             >*/
    if (rxp <= (float)0.) {
        goto L8;
    }
/*<       do 4 i=1,nl                                                        >*/
    i__2 = nl;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       js(i+nv)=(sc(i)-xm)/rxp+.5                                         >*/
        js[i__ + nv] = (sc[i__] - xm) / rxp + (float).5;
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       if(nl .gt. 26) go to 5                                             >*/
    if (nl > 26) {
        goto L5;
    }
/*    write(it,28) (i,i=1,nl)                                         
    472*/
/*    write(it,28) (js(i+nv),i=1,nl)                                  
    473*/
/*<       go to 8                                                            >*/
    goto L8;
/*<     5 if(nl .gt. 38) go to 6                                             >*/
L5:
    if (nl > 38) {
        goto L6;
    }
/*    write(it,29) (i,i=1,nl)                                         
    476*/
/*    write(it,29) (js(i+nv),i=1,nl)                                  
    477*/
/*<       go to 8                                                            >*/
    goto L8;
/*<     6 if(nl .gt. 78) go to 7                                             >*/
L6:
    if (nl > 78) {
        goto L7;
    }
/*    write(it,30) (mod(i,10),i=1,nl)                                 
    480*/
/*    write(it,30) (js(i+nv),i=1,nl)                                  
    481*/
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 continue >*/
L7:
/*    write(it,37) 78                                                 
    483*/
/*<     8 continue                                                           >*/
L8:
    ;
    }
/*<       call catv(2,kp,kv,nv,js)                                           >*/
    catv_(&c__2, &kp[6], &kv[3], &nv, &js[1]);
/*<       do 27 jj=1,nv                                                      >*/
    i__1 = nv;
    for (jj = 1; jj <= i__1; ++jj) {
/*<       j1=js(2*jj-1)                                                      >*/
    j1 = js[(jj << 1) - 1];
/*<       j2=js(2*jj)                                                        >*/
    j2 = js[jj * 2];
/*<       xm=big                                                             >*/
    xm = big;
/*<       xx=-big                                                            >*/
    xx = -big;
/*<       n1=int(cm(2*j1+1)+.1)-int(cm(2*j1)+.1)+1                           >*/
    n1 = (integer) (cm[(j1 << 1) + 1] + (float).1) - (integer) (cm[j1 * 2]
         + (float).1) + 1;
/*<       n2=int(cm(2*j2+1)+.1)-int(cm(2*j2)+.1)+1                           >*/
    n2 = (integer) (cm[(j2 << 1) + 1] + (float).1) - (integer) (cm[j2 * 2]
         + (float).1) + 1;
/*<       k=0                                                                >*/
    k = 0;
/*<       do 10 i=1,n1                                                       >*/
    i__2 = n1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s1=cvlv(m,1,j1,i,nk,kp,kv,tb,cm,tc)                                >*/
        s1 = cvlv_(m, &c__1, &j1, &i__, nk, &kp[6], &kv[3], tb, &cm[1], &
            tc[1]);
/*<       js(2*nv+1)=i                                                       >*/
        js[(nv << 1) + 1] = i__;
/*<       do 9 j=1,n2                                                        >*/
        i__3 = n2;
        for (j = 1; j <= i__3; ++j) {
/*<       js(2*nv+2)=j                                                       >*/
        js[(nv << 1) + 2] = j;
/*<       k=k+1                                                              >*/
        ++k;
/*<       sc(k)=s1+cvlv(m,2,js(2*jj-1),js(2*nv+1),nk,kp,kv,tb,cm,tc)         >*/
        sc[k] = s1 + cvlv_(m, &c__2, &js[(jj << 1) - 1], &js[(nv << 1)
             + 1], nk, &kp[6], &kv[3], tb, &cm[1], &tc[1]);
/*<     9 continue                                                           >*/
/* L9: */
        }
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       do 12 j=1,n2                                                       >*/
    i__2 = n2;
    for (j = 1; j <= i__2; ++j) {
/*<       s1=cvlv(m,1,j2,j,nk,kp,kv,tb,cm,tc)                                >*/
        s1 = cvlv_(m, &c__1, &j2, &j, nk, &kp[6], &kv[3], tb, &cm[1], &tc[
            1]);
/*<       do 11 i=1,n1                                                       >*/
        i__3 = n1;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       k=j+n2*(i-1)                                                       >*/
        k = j + n2 * (i__ - 1);
/*<       sc(k)=sc(k)+s1                                                     >*/
        sc[k] += s1;
/*<    11 continue                                                           >*/
/* L11: */
        }
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<       k=0                                                                >*/
    k = 0;
/*<       do 14 i=1,n1                                                       >*/
    i__2 = n1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       do 13 j=1,n2                                                       >*/
        i__3 = n2;
        for (j = 1; j <= i__3; ++j) {
/*<       k=k+1                                                              >*/
        ++k;
/*<       x=sc(k)                                                            >*/
        x = sc[k];
/*<       xx=amax1(xx,x)                                                     >*/
        xx = dmax(xx,x);
/*<       xm=amin1(xm,x)                                                     >*/
        xm = dmin(xm,x);
/*<    13 continue                                                           >*/
/* L13: */
        }
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       na=min0(n1,n2)                                                     >*/
    na = min(n1,n2);
/*<       nb=max0(n1,n2)                                                     >*/
    nb = max(n1,n2);
/*<       if(na .ne. n1) go to 15                                            >*/
    if (na != n1) {
        goto L15;
    }
/*<       ja=j1                                                              >*/
    ja = j1;
/*<       jb=j2                                                              >*/
    jb = j2;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 ja=j2                                                              >*/
L15:
    ja = j2;
/*<       jb=j1                                                              >*/
    jb = j1;
/*<    16 px=99.0                                                            >*/
L16:
    px = (float)99.;
/*<       if(na.gt.25) px=9.0                                                >*/
    if (na > 25) {
        px = (float)9.;
    }
/*<       rx=xx-xm                                                           >*/
    rx = xx - xm;
/*<       rxp=rx/px                                                          >*/
    rxp = rx / px;
/*    write(it,'(/,'' f( x('',i3,''), x('',i3,'') ) : scale ='',g12.4)
')  531*/
/*   1ja,jb,rxp                                                       
    532*/
/*<       if(rxp.le.0.0) go to 27                                            >*/
    if (rxp <= (float)0.) {
        goto L27;
    }
/*<       if(na .le. 75) go to 17                                            >*/
    if (na <= 75) {
        goto L17;
    }
/*    write(it,37) 75                                                 
    535*/
/*<       go to 27                                                           >*/
    goto L27;
/*<    17 if(na .gt. 25) go to 18                                            >*/
L17:
    if (na > 25) {
        goto L18;
    }
/*    write(it,34) (i,i=1,na)                                         
    538*/
/*<       go to 20                                                           >*/
    goto L20;
/*<    18 if(na .gt. 37) go to 19                                            >*/
L18:
    if (na > 37) {
        goto L19;
    }
/*    write(it,35) (i,i=1,na)                                         
    541*/
/*<       go to 20                                                           >*/
    goto L20;
/*<    19 continue >*/
L19:
/*    write(it,36) (mod(i,10),i=1,na)                                 
    543*/
/*<    20 do 26 j=1,nb                                                       >*/
L20:
    i__2 = nb;
    for (j = 1; j <= i__2; ++j) {
/*<       do 23 i=1,na                                                       >*/
        i__3 = na;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       if(na .ne. n1) go to 21                                            >*/
        if (na != n1) {
            goto L21;
        }
/*<       k=j+n2*(i-1)                                                       >*/
        k = j + n2 * (i__ - 1);
/*<       go to 22                                                           >*/
        goto L22;
/*<    21 k=i+n2*(j-1)                                                       >*/
L21:
        k = i__ + n2 * (j - 1);
/*<    22 js(i+2*nv)=(sc(k)-xm)/rxp+.5                                       >*/
L22:
        js[i__ + (nv << 1)] = (sc[k] - xm) / rxp + (float).5;
/*<    23 continue                                                           >*/
/* L23: */
        }
/*<       if(na .gt. 25) go to 24                                            >*/
        if (na > 25) {
        goto L24;
        }
/*    write(it,31) j,(js(i+2*nv),i=1,na)                          
        553*/
/*<       go to 26                                                           >*/
        goto L26;
/*<    24 if(na .gt. 37) go to 25                                            >*/
L24:
        if (na > 37) {
        goto L25;
        }
/*    write(it,32) j,(js(i+2*nv),i=1,na)                          
        556*/
/*<       go to 26                                                           >*/
        goto L26;
/*<    25 continue >*/
L25:
/*    write(it,33) j,(js(i+2*nv),i=1,na)                          
        558*/
/*<    26 continue                                                           >*/
L26:
        ;
    }
/*<    27 continue                                                           >*/
L27:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       entry printc(nal)                                                  >*/

L_printc:
/*<       it=nal                                                             >*/
    it = *nal;
/*<       return                                                             >*/
    return 0;
/*<    28 format(' ',26i3)                                                   >*/
/* L28: */
/*<    29 format(' ',38i2)                                                   >*/
/* L29: */
/*<    30 format(' ',78i1)                                                   >*/
/* L30: */
/*<    31 format(' ',i3,' ',25i3)                                            >*/
/* L31: */
/*<    32 format(' ',i3,' ',37i2)                                            >*/
/* L32: */
/*<    33 format(' ',i3,' ',75i1)                                            >*/
/* L33: */
/*<    34 format('     ',25i3)                                               >*/
/* L34: */
/*<    35 format('     ',37i2)                                               >*/
/* L35: */
/*<    36 format('     ',75i1)                                               >*/
/* L36: */
/*<    >*/
/* L37: */
/*<       end                                                                >*/
} /* ctprt1_ */

/* Subroutine */ int ctprt1_(m, nk, kp, kv, tb, cm, tc, sc, js)
integer *m, *nk, *kp, *kv;
real *tb, *cm, *tc, *sc;
integer *js;
{
    return ctprt1_0_(0, m, nk, kp, kv, tb, cm, tc, sc, js, (integer *)0);
    }

/* Subroutine */ int printc_(nal)
integer *nal;
{
    return ctprt1_0_(1, (integer *)0, (integer *)0, (integer *)0, (integer *)
        0, (real *)0, (real *)0, (real *)0, (real *)0, (integer *)0, nal);
    }

/*<    >*/
/* Subroutine */ int slice1_0_(n__, flg, xs, n, p, x, nk, az, tb, cm, kp, kv, 
    lp, lv, bz, tc, azn, tbn, kpn, kvn, lpn, lvn, bzn, tcn, sp, mm, ig)
int n__;
real *flg, *xs;
integer *n, *p;
real *x;
integer *nk;
real *az, *tb, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *azn, *tbn;
integer *kpn, *kvn, *lpn, *lvn;
real *bzn, *tcn, *sp;
integer *mm, *ig;
{
    /* Initialized data */

    static integer it = 6;
    static real big = (float)9.9e30;

    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, m;
    extern /* Subroutine */ int ccoll_(), slova_();
    static integer i1, i2, i3, ni;
    static real xl, xr;
    extern /* Subroutine */ int reducl_(), qslice_(), reducq_();

/*<    >*/
/*<       real xs(p),x(n,p),tb(5,nk),cm(*),tc(*),tbn(5,nk),tcn(*),sp(*)      >*/
/*<       data it,big /6,9.9e30/                                             >*/
    /* Parameter adjustments */
    if (xs) {
    --xs;
    }
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (tb) {
    tb -= 6;
    }
    if (tbn) {
    tbn -= 6;
    }
    if (cm) {
    --cm;
    }
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (lp) {
    lp -= 4;
    }
    if (lv) {
    --lv;
    }
    if (tc) {
    --tc;
    }
    if (kpn) {
    kpn -= 6;
    }
    if (kvn) {
    kvn -= 3;
    }
    if (lpn) {
    lpn -= 4;
    }
    if (lvn) {
    --lvn;
    }
    if (tcn) {
    --tcn;
    }
    if (sp) {
    --sp;
    }
    if (mm) {
    --mm;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_prtslc;
    }

/*<       ni=0                                                               >*/
    ni = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).ne.0.0) ni=ni+1                                         >*/
    if (tb[m * 5 + 1] != (float)0.) {
        ++ni;
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       if(ni .ne. 0) go to 3                                              >*/
    if (ni != 0) {
    goto L3;
    }
/*<       kpn(1,1)=-1                                                        >*/
    kpn[6] = -1;
/*<       lpn(1,1)=0                                                         >*/
    lpn[4] = 0;
/*<       do 2 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       tbn(1,m)=0.0                                                       >*/
    tbn[m * 5 + 1] = (float)0.;
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       azn=0.0                                                            >*/
    *azn = (float)0.;
/*<       bzn=azn                                                            >*/
    *bzn = *azn;
/*    if(it.gt.0) write(it,'('' slice: original mars model = constant.''  
595*/
/*   1)')                                                                 
596*/
/*<       return                                                             >*/
    return 0;
/*<     3 if(it .le. 0) go to 5                                              >*/
L3:
    if (it <= 0) {
    goto L5;
    }
/*    write(it,'(/,'' sliced mars model: flag ='',g12.4)') flg            
599*/
/*    write(it,'(/,'' slice:'')')                                         
600*/
/*<       do 4 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       if(xs(j).eq.flg) go to 4                                           >*/
    if (xs[j] == *flg) {
        goto L4;
    }
/*    write(it,'('' x('',i3,'') ='',g12.4)') j,xs(j)                  
    603*/
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<     5 i1=2*nk+1                                                          >*/
L5:
    i1 = (*nk << 1) + 1;
/*<       i2=i1+2*p                                                          >*/
    i2 = i1 + (*p << 1);
/*<       i3=max0(i2+p,i1+nk)                                                >*/
/* Computing MAX */
    i__1 = i2 + *p, i__2 = i1 + *nk;
    i3 = max(i__1,i__2);
/*<       do 7 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       xl=big                                                             >*/
    xl = big;
/*<       xr=-xl                                                             >*/
    xr = -xl;
/*<       do 6 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       xl=amin1(xl,x(i,j))                                                >*/
/* Computing MIN */
        r__1 = xl, r__2 = x[i__ + j * x_dim1];
        xl = dmin(r__1,r__2);
/*<       xr=amax1(xr,x(i,j))                                                >*/
/* Computing MAX */
        r__1 = xr, r__2 = x[i__ + j * x_dim1];
        xr = dmax(r__1,r__2);
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       sp(j+i3-1)=xr-xl                                                   >*/
    sp[j + i3 - 1] = xr - xl;
/*<       sp(j+i3-1+p)=xl                                                    >*/
    sp[j + i3 - 1 + *p] = xl;
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<    >*/
    reducq_(flg, &xs[1], nk, &tb[6], &cm[1], &tc[1], &kp[6], &kv[3], &lp[4], &
        lv[1], &sp[i3], &sp[1], &sp[i1], &sp[i2]);
/*<       call reducl(flg,xs,nk,az,tb,cm,bz,sp,sp(i3),azn,tbn,bzn,sp(i1))    >*/
    reducl_(flg, &xs[1], nk, az, &tb[6], &cm[1], bz, &sp[1], &sp[i3], azn, &
        tbn[6], bzn, &sp[i1]);
/*<       ni=0                                                               >*/
    ni = 0;
/*<       do 8 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tbn(1,m).ne.0.0) ni=ni+1                                        >*/
    if (tbn[m * 5 + 1] != (float)0.) {
        ++ni;
    }
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       if(ni .ne. 0) go to 10                                             >*/
    if (ni != 0) {
    goto L10;
    }
/*<       kpn(1,1)=-1                                                        >*/
    kpn[6] = -1;
/*<       lpn(1,1)=0                                                         >*/
    lpn[4] = 0;
/*<       do 9 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       tbn(1,m)=0.0                                                       >*/
    tbn[m * 5 + 1] = (float)0.;
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       azn=0.0                                                            >*/
    *azn = (float)0.;
/*<       bzn=azn                                                            >*/
    *bzn = *azn;
/*    if(it.gt.0) write(it,'('' sliced mars model = constant.'')')        
633*/
/*<       return                                                             >*/
    return 0;
/*<    10 if(it.gt.0) call slova(nk,it,tbn,ni,lpn,lvn)                       >*/
L10:
    if (it > 0) {
    slova_(nk, &it, &tbn[6], &ni, &lpn[4], &lvn[1]);
    }
/*<       call ccoll(nk,tbn,cm,kpn,kvn,lpn,lvn,mm)                           >*/
    ccoll_(nk, &tbn[6], &cm[1], &kpn[6], &kvn[3], &lpn[4], &lvn[1], &mm[1]);
/*<       call qslice(p,nk,tbn,cm,sp,kpn,kvn,lpn,lvn,tcn,sp(i3),sp(i1),mm)   >*/
    qslice_(p, nk, &tbn[6], &cm[1], &sp[1], &kpn[6], &kvn[3], &lpn[4], &lvn[1]
        , &tcn[1], &sp[i3], &sp[i1], &mm[1]);
/*<       return                                                             >*/
    return 0;
/*<       entry prtslc(ig)                                                   >*/

L_prtslc:
/*<       it=ig                                                              >*/
    it = *ig;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* slice1_ */

/* Subroutine */ int slice1_(flg, xs, n, p, x, nk, az, tb, cm, kp, kv, lp, lv,
     bz, tc, azn, tbn, kpn, kvn, lpn, lvn, bzn, tcn, sp, mm)
real *flg, *xs;
integer *n, *p;
real *x;
integer *nk;
real *az, *tb, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *azn, *tbn;
integer *kpn, *kvn, *lpn, *lvn;
real *bzn, *tcn, *sp;
integer *mm;
{
    return slice1_0_(0, flg, xs, n, p, x, nk, az, tb, cm, kp, kv, lp, lv, bz, 
        tc, azn, tbn, kpn, kvn, lpn, lvn, bzn, tcn, sp, mm, (integer *)0);
    }

/* Subroutine */ int prtslc_(ig)
integer *ig;
{
    return slice1_0_(1, (real *)0, (real *)0, (integer *)0, (integer *)0, (
        real *)0, (integer *)0, (real *)0, (real *)0, (real *)0, (integer 
        *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)
        0, (real *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0,
         (integer *)0, (real *)0, (real *)0, (real *)0, (integer *)0, ig);
    }

/*<       subroutine cmrs (n,x,cm,kp,kv,lp,lv,bz,tc,y,sc)                    >*/
/* Subroutine */ int cmrs_0_(n__, n, x, cm, kp, kv, lp, lv, bz, tc, y, sc, i1)
int n__;
integer *n;
real *x, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *y, *sc;
integer *i1;
{
    /* Initialized data */

    static integer ifg = 0;

    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, i__1, i__2, i__3;

    /* Local variables */
    extern integer icat_();
    static integer i__, j, k, l, m, l1, ic, la, lb, jj, il, ll, jl, kk, nt, 
        kp3;
    extern /* Subroutine */ int que_();

/*<       integer kp(5,*),kv(2,*),lp(3,*),lv(*)                              >*/
/*<       real x(n,*),tc(*),cm(*),y(n),sc(n,2)                               >*/
/*<       data ifg /0/                                                       >*/
    /* Parameter adjustments */
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (y) {
    --y;
    }
    if (sc) {
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    }
    if (cm) {
    --cm;
    }
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (lp) {
    lp -= 4;
    }
    if (lv) {
    --lv;
    }
    if (tc) {
    --tc;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_stcmrs;
    }

/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       y(i)=bz                                                            >*/
    y[i__] = *bz;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       ll=1                                                               >*/
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<     2 if(kp(1,ll).lt.0) go to 19                                         >*/
L2:
    if (kp[ll * 5 + 1] < 0) {
    goto L19;
    }
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sc(i,1)=1.0                                                        >*/
    sc[i__ + sc_dim1] = (float)1.;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       if(kp(1,ll) .le. 0) go to 11                                       >*/
    if (kp[ll * 5 + 1] <= 0) {
    goto L11;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       do 10 il=1,jl                                                      >*/
    i__1 = jl;
    for (il = 1; il <= i__1; ++il) {
/*<       k=kp(2,ll)+il-1                                                    >*/
    k = kp[ll * 5 + 2] + il - 1;
/*<       jj=kv(1,k)                                                         >*/
    jj = kv[(k << 1) + 1];
/*<       j=iabs(jj)                                                         >*/
    j = abs(jj);
/*<       kk=kv(2,k)                                                         >*/
    kk = kv[(k << 1) + 2];
/*<       do 9 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(sc(i,1).eq.0.0) go to 9                                         >*/
        if (sc[i__ + sc_dim1] == (float)0.) {
        goto L9;
        }
/*<       if(ifg .ne. 0) go to 4                                             >*/
        if (ifg != 0) {
        goto L4;
        }
/*<       ic=icat(x(i,j),j,cm)                                               >*/
        ic = icat_(&x[i__ + j * x_dim1], &j, &cm[1]);
/*<       go to 5                                                            >*/
        goto L5;
/*<     4 ic=x(i,j)+.1                                                       >*/
L4:
        ic = x[i__ + j * x_dim1] + (float).1;
/*<     5 if(ic .ne. 0) go to 6                                              >*/
L5:
        if (ic != 0) {
        goto L6;
        }
/*<       sc(i,1)=0.0                                                        >*/
        sc[i__ + sc_dim1] = (float)0.;
/*<       go to 7                                                            >*/
        goto L7;
/*<     6 sc(i,1)=cm(ic+kk)                                                  >*/
L6:
        sc[i__ + sc_dim1] = cm[ic + kk];
/*<     7 if(jj .ge. 0) go to 9                                              >*/
L7:
        if (jj >= 0) {
        goto L9;
        }
/*<       if(sc(i,1) .ne. 0.0) go to 8                                       >*/
        if (sc[i__ + sc_dim1] != (float)0.) {
        goto L8;
        }
/*<       sc(i,1)=1.0                                                        >*/
        sc[i__ + sc_dim1] = (float)1.;
/*<       go to 9                                                            >*/
        goto L9;
/*<     8 sc(i,1)=0.0                                                        >*/
L8:
        sc[i__ + sc_dim1] = (float)0.;
/*<     9 continue                                                           >*/
L9:
        ;
    }
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       go to 12                                                           >*/
    goto L12;
/*<    11 if(kp(3,ll) .gt. 0) go to 12                                       >*/
L11:
    if (kp[ll * 5 + 3] > 0) {
    goto L12;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<    12 if(kp(3,ll) .ge. 0) go to 14                                       >*/
L12:
    if (kp[ll * 5 + 3] >= 0) {
    goto L14;
    }
/*<       k=-kp(3,ll)                                                        >*/
    k = -kp[ll * 5 + 3];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       do 13 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(sc(i,1).eq.0.0) go to 13                                        >*/
    if (sc[i__ + sc_dim1] == (float)0.) {
        goto L13;
    }
/*<       y(i)=y(i)+tc(k)                                                    >*/
    y[i__] += tc[k];
/*<    13 continue                                                           >*/
L13:
    ;
    }
/*<       go to 2                                                            >*/
    goto L2;
/*<    14 kp3=kp(3,ll)                                                       >*/
L14:
    kp3 = kp[ll * 5 + 3];
/*<       do 18 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       lb=la+5*l*nt-1                                                     >*/
    lb = la + l * 5 * nt - 1;
/*<       do 17 j=1,nt                                                       >*/
    i__2 = nt;
    for (j = 1; j <= i__2; ++j) {
/*<       do 15 i=1,n                                                        >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       sc(i,2)=sc(i,1)                                                    >*/
        sc[i__ + (sc_dim1 << 1)] = sc[i__ + sc_dim1];
/*<    15 continue                                                           >*/
/* L15: */
        }
/*<       call que(j,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,2))                   >*/
        que_(&j, &l, &nt, &lv[lp[l1 * 3 + 2]], n, &x[x_offset], &tc[la], &
            sc[(sc_dim1 << 1) + 1]);
/*<       do 16 i=1,n                                                        >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       y(i)=y(i)+tc(lb+j)*sc(i,2)                                         >*/
        y[i__] += tc[lb + j] * sc[i__ + (sc_dim1 << 1)];
/*<    16 continue                                                           >*/
/* L16: */
        }
/*<    17 continue                                                           >*/
/* L17: */
    }
/*<       la=lb+nt+1                                                         >*/
    la = lb + nt + 1;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<    18 continue                                                           >*/
/* L18: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<    19 return                                                             >*/
L19:
    return 0;
/*<       entry stcmrs(i1)                                                   >*/

L_stcmrs:
/*<       ifg=i1                                                             >*/
    ifg = *i1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* cmrs_ */

/* Subroutine */ int cmrs_(n, x, cm, kp, kv, lp, lv, bz, tc, y, sc)
integer *n;
real *x, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *y, *sc;
{
    return cmrs_0_(0, n, x, cm, kp, kv, lp, lv, bz, tc, y, sc, (integer *)0);
    }

/* Subroutine */ int stcmrs_(i1)
integer *i1;
{
    return cmrs_0_(1, (integer *)0, (real *)0, (real *)0, (integer *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, i1);
    }

/*<       subroutine fmrs (n,x,nk,az,tb,cm,y)                                >*/
/* Subroutine */ int fmrs_0_(n__, n, x, nk, az, tb, cm, y, i1)
int n__;
integer *n;
real *x;
integer *nk;
real *az, *tb, *cm, *y;
integer *i1;
{
    /* Initialized data */

    static integer ifg = 0;

    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    extern integer icat_();
    static integer i__, j, k, m;
    static doublereal s;
    static real t, u;
    static integer ip;
    static real phi;

/*<       real x(n,*),tb(5,nk),cm(*),y(n)                                    >*/
/*<       double precision s                                                 >*/
/*<       data ifg /0/                                                       >*/
    /* Parameter adjustments */
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (y) {
    --y;
    }
    if (tb) {
    tb -= 6;
    }
    if (cm) {
    --cm;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_stfmrs;
    }

/*<       do 13 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=az                                                               >*/
    s = *az;
/*<       do 12 m=1,nk                                                       >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 12                                        >*/
        if (tb[m * 5 + 1] == (float)0.) {
        goto L12;
        }
/*<       phi=1.0                                                            >*/
        phi = (float)1.;
/*<       ip=m                                                               >*/
        ip = m;
/*<     1 if(ip.le.0) go to 11                                               >*/
L1:
        if (ip <= 0) {
        goto L11;
        }
/*<       t=tb(2,ip)                                                         >*/
        t = tb[ip * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
        j = dabs(t) + (float).1;
/*<       if(cm(2*j) .le. 0.0) go to 8                                       >*/
        if (cm[j * 2] <= (float)0.) {
        goto L8;
        }
/*<       if(ifg .ne. 0) go to 2                                             >*/
        if (ifg != 0) {
        goto L2;
        }
/*<       k=icat(x(i,j),j,cm)                                                >*/
        k = icat_(&x[i__ + j * x_dim1], &j, &cm[1]);
/*<       go to 3                                                            >*/
        goto L3;
/*<     2 k=x(i,j)+.1                                                        >*/
L2:
        k = x[i__ + j * x_dim1] + (float).1;
/*<     3 if(k .ne. 0) go to 4                                               >*/
L3:
        if (k != 0) {
        goto L4;
        }
/*<       u=0.0                                                              >*/
        u = (float)0.;
/*<       go to 5                                                            >*/
        goto L5;
/*<     4 u=cm(k+int(tb(3,ip)+.1))                                           >*/
L4:
        u = cm[k + (integer) (tb[ip * 5 + 3] + (float).1)];
/*<     5 if(t .ge. 0.0) go to 9                                             >*/
L5:
        if (t >= (float)0.) {
        goto L9;
        }
/*<       if(u .ne. 0.0) go to 6                                             >*/
        if (u != (float)0.) {
        goto L6;
        }
/*<       u=1.0                                                              >*/
        u = (float)1.;
/*<       go to 9                                                            >*/
        goto L9;
/*<     6 u=0.0                                                              >*/
L6:
        u = (float)0.;
/*<       go to 9                                                            >*/
        goto L9;
/*<     8 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,ip)))                         >*/
L8:
/* Computing MAX */
        r__1 = (float)0., r__2 = r_sign(&c_b196, &t) * (x[i__ + j * 
            x_dim1] - tb[ip * 5 + 3]);
        u = dmax(r__1,r__2);
/*<     9 if(u .ne. 0.0) go to 10                                            >*/
L9:
        if (u != (float)0.) {
        goto L10;
        }
/*<       phi=0.0                                                            >*/
        phi = (float)0.;
/*<       go to 11                                                           >*/
        goto L11;
/*<    10 phi=phi*u                                                          >*/
L10:
        phi *= u;
/*<       ip=tb(4,ip)+.1                                                     >*/
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
        goto L1;
/*<    11 s=s+tb(1,m)*phi                                                    >*/
L11:
        s += tb[m * 5 + 1] * phi;
/*<    12 continue                                                           >*/
L12:
        ;
    }
/*<       y(i)=s                                                             >*/
    y[i__] = s;
/*<    13 continue                                                           >*/
/* L13: */
    }
/*<       return                                                             >*/
    return 0;
/*<       entry stfmrs(i1)                                                   >*/

L_stfmrs:
/*<       ifg=i1                                                             >*/
    ifg = *i1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* fmrs_ */

/* Subroutine */ int fmrs_(n, x, nk, az, tb, cm, y)
integer *n;
real *x;
integer *nk;
real *az, *tb, *cm, *y;
{
    return fmrs_0_(0, n, x, nk, az, tb, cm, y, (integer *)0);
    }

/* Subroutine */ int stfmrs_(i1)
integer *i1;
{
    return fmrs_0_(1, (integer *)0, (real *)0, (integer *)0, (real *)0, (real 
        *)0, (real *)0, (real *)0, i1);
    }

/*<    >*/
/* Subroutine */ int marsgo_0_(n__, n, p, x, y, w, nk, ms, df, fv, mi, lx, it,
     xm, xs, az, tb, cm, sc, db, d__, mm, val, nal)
int n__;
integer *n, *p;
real *x, *y, *w;
integer *nk, *ms;
real *df, *fv;
integer *mi, *lx, *it;
real *xm, *xs, *az, *tb, *cm, *sc;
doublereal *db, *d__;
integer *mm;
real *val;
integer *nal;
{
    /* Initialized data */

    static integer ix = 0;
    static doublereal alr = 1e-7;
    static doublereal eps = 1e-4;
    static real big = (float)9.9e30;
    static real fln = (float)-1.;
    static integer nmin = 5;
    static real alf = (float).05;
    static real vcst[3] = { (float)1.,(float).666667,(float).333333 };

    /* System generated locals */
    integer mm_dim1, mm_offset, x_dim1, x_offset, sc_dim1, sc_offset, db_dim1,
         db_offset, d_dim1, d_offset, i__1, i__2;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern integer newb_();
    extern /* Subroutine */ int holl_(), sscp_();
    static real tcmx, tcst;
    static integer nopt, mtot;
    static doublereal a, b;
    static real h__;
    static integer i__, j, k, l, m;
    static doublereal s, t, u, v;
    extern /* Subroutine */ int isfac_();
    static logical newbf;
    extern /* Subroutine */ int array_();
    extern integer nnord_();
    extern /* Subroutine */ int bkstp_();
    static integer j0, k1;
    static real tcsts;
    static integer ja;
    extern integer jf_();
    static integer me, nc;
    static doublereal yb, we, xb, xd, dx, wn, se, tt, st, sw;
    static real tx[5];
    static doublereal sy, xt, yv, xx, su, yc, dy, dv;
    static integer mk, kr;
    extern /* Subroutine */ int addpar_();
    static integer jq;
    static real df1;
    extern /* Subroutine */ int updpar_();
    static integer jp;
    extern /* Subroutine */ int mnspan_();
    static integer mn, jd1;
    extern /* Subroutine */ int itrpar_();
    static integer jd2;
    extern /* Subroutine */ int update_();
    static real xa;
    static integer mj;
    static real sj;
    static integer mp;
    extern /* Subroutine */ int selpar_(), getnst_();
    static integer jn, mm1;
    extern /* Subroutine */ int nxtpar_();
    extern integer ibfext_();
    static real xk;
    static integer kl, ll;
    extern /* Subroutine */ int coefpr_();
    static real tx1;
    static integer lbf;
    extern /* Subroutine */ int blf_();
    extern logical elg_();
    static integer kcp, mel, ict;
    static char hol[28];
    static integer nep;
    extern integer jft_();
    extern /* Subroutine */ int csp_();
    static integer jas, nnl;
    static real cst;
    extern doublereal phi_();
    static integer nop;
    static real fvr;
    static integer nnr, nnt;
    static real fjn, fkr, gcv;
    static integer nli;
    static real txi;
    static doublereal rsq, ssq;
    static real txl, txm;
    static integer nst;
    static real asm__;
    static doublereal txt;
    extern /* Subroutine */ int blf0_();
    static integer kcp0;
    static doublereal asq0;
    extern /* Subroutine */ int lsf1_();
    static integer mkp1, mkp2;
    static real cfac;

/*<       integer p,mm(n,p),lx(p)                                            >*/
/*<       logical elg,newbf                                                  >*/
/*<    >*/
/*<       double precision db(n,*),d(nk,*)                                   >*/
/*<       double precision yb,yv,sw,s,t,u,v,we,sy,a,b,xb,xx,xd,ssq,alr       >*/
/*<       double precision dx,wn,se,tt,txt,xt,st,su,yc,eps,rsq,dy,dv,asq0    >*/
/*<       character*28 hol                                                   >*/
/*<    >*/
    /* Parameter adjustments */
    if (y) {
    --y;
    }
    if (w) {
    --w;
    }
    if (sc) {
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    }
    if (db) {
    db_dim1 = *n;
    db_offset = db_dim1 + 1;
    db -= db_offset;
    }
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (lx) {
    --lx;
    }
    if (xm) {
    --xm;
    }
    if (xs) {
    --xs;
    }
    if (mm) {
    mm_dim1 = *n;
    mm_offset = mm_dim1 + 1;
    mm -= mm_offset;
    }
    if (tb) {
    tb -= 6;
    }
    if (d__) {
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    }
    if (cm) {
    --cm;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_setfln;
    case 2: goto L_setalf;
    case 3: goto L_setmin;
    case 4: goto L_setcta;
    case 5: goto L_setctl;
    case 6: goto L_xvmrgo;
    case 7: goto L_setalr;
    }

/*    if(it.gt.0) write(it,97)                                            
773*/
/*<       mk=nk                                                              >*/
    mk = *nk;
/*<       df1=0.0                                                            >*/
    df1 = (float)0.;
/*<       nep=0                                                              >*/
    nep = 0;
/*<       do 1 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       if(lx(j).eq.0) go to 1                                             >*/
    if (lx[j] == 0) {
        goto L1;
    }
/*<       if(x(mm(1,j),j).ge.x(mm(n,j),j)) go to 1                           >*/
    if (x[mm[j * mm_dim1 + 1] + j * x_dim1] >= x[mm[*n + j * mm_dim1] + j 
        * x_dim1]) {
        goto L1;
    }
/*<       nep=nep+1                                                          >*/
    ++nep;
/*<       cst=vcst(iabs(lx(j)))                                              >*/
    cst = vcst[(i__2 = lx[j], abs(i__2)) - 1];
/*<       if(mi.eq.1) cst=amin1(cst,vcst(2))                                 >*/
    if (*mi == 1) {
        cst = dmin(cst,vcst[1]);
    }
/*<       df1=df1+cst                                                        >*/
    df1 += cst;
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<       if(nep .ne. 0) go to 2                                             >*/
    if (nep != 0) {
    goto L2;
    }
/*    if(it.gt.0) write(it,'('' no predictor variables.'')')              
786*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     2 if(nep.eq.1) df1=vcst(3)                                           >*/
L2:
    if (nep == 1) {
    df1 = vcst[2];
    }
/*<       cfac=df1/nep                                                       >*/
    cfac = df1 / nep;
/*<       df1=df*cfac                                                        >*/
    df1 = *df * cfac;
/*<       mkp1=mk+1                                                          >*/
    mkp1 = mk + 1;
/*<       mkp2=mk+2                                                          >*/
    mkp2 = mk + 2;
/*<       sw=0.d0                                                            >*/
    sw = 0.;
/*<       wn=sw                                                              >*/
    wn = sw;
/*<       yb=wn                                                              >*/
    yb = wn;
/*<       s=yb                                                               >*/
    s = yb;
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<       yb=yb+w(i)*y(i)                                                    >*/
    yb += w[i__] * y[i__];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       yb=yb/sw                                                           >*/
    yb /= sw;
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    d__1 = sw;
    wn = d__1 * d__1 / wn;
/*<       do 4 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=s+w(i)*(y(i)-yb)**2                                              >*/
/* Computing 2nd power */
    d__1 = y[i__] - yb;
    s += w[i__] * (d__1 * d__1);
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       yv=s/sw                                                            >*/
    yv = s / sw;
/*<       tcst=1.0                                                           >*/
    tcst = (float)1.;
/*<       tcmx=wn-df1*vcst(1)-2.0                                            >*/
    tcmx = wn - df1 * vcst[0] - (float)2.;
/*<       if(cm(1) .le. 0.0) go to 7                                         >*/
    if (cm[1] <= (float)0.) {
    goto L7;
    }
/*<       i=2                                                                >*/
    i__ = 2;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 i=i+(2)                                                            >*/
L5:
    i__ += 2;
/*<     6 if((2)*((i)-(2*p)).gt.0) go to 7                                   >*/
L6:
    if (i__ - (*p << 1) << 1 > 0) {
    goto L7;
    }
/*<       if(cm(i).gt.0.0) kcp0=cm(i+1)+.1                                   >*/
    if (cm[i__] > (float)0.) {
    kcp0 = cm[i__ + 1] + (float).1;
    }
/*<       go to 5                                                            >*/
    goto L5;
/*<     7 m=0                                                                >*/
L7:
    m = 0;
/*<       mtot=m                                                             >*/
    mtot = m;
/*<       txm=yv/(1.d0-1.d0/wn)**2                                           >*/
/* Computing 2nd power */
    d__1 = 1. - 1. / wn;
    txm = yv / (d__1 * d__1);
/*<       rsq=yv*sw                                                          >*/
    rsq = yv * sw;
/*<       kr=0                                                               >*/
    kr = 0;
/*<       nopt=0                                                             >*/
    nopt = 0;
/*    if(it.gt.0) write(it,98) m,txm,0.0,1.0                              
823*/
/*<       if(fln.lt.0.0) fln=1.0+4.0/wn                                      >*/
    if (fln < (float)0.) {
    fln = (float)4. / wn + (float)1.;
    }
/*<       call addpar(0)                                                     >*/
    addpar_(&c__0);
/*<     8 if(m.ge.mk.or.tcst.ge.tcmx) go to 69                               >*/
L8:
    if (m >= mk || tcst >= tcmx) {
    goto L69;
    }
/*<       nopt=nopt+1                                                        >*/
    ++nopt;
/*<       call itrpar(nopt)                                                  >*/
    itrpar_(&nopt);
/*<       mm1=m                                                              >*/
    mm1 = m;
/*<       m=m+1                                                              >*/
    ++m;
/*<       txi=big                                                            >*/
    txi = big;
/*<       kcp=kcp0                                                           >*/
    kcp = kcp0;
/*<       asq0=rsq/sw                                                        >*/
    asq0 = rsq / sw;
/*<     9 call nxtpar(l,jq)                                                  >*/
L9:
    nxtpar_(&l, &jq);
/*<       if(l.lt.0) go to 53                                                >*/
    if (l < 0) {
    goto L53;
    }
/*<       txl=big                                                            >*/
    txl = big;
/*<       if(nnord(l,tb) .lt. mi) go to 10                                   >*/
    if (nnord_(&l, &tb[6]) < *mi) {
    goto L10;
    }
/*<       call updpar(0,-1.d0)                                               >*/
    updpar_(&c__0, &c_b214);
/*<       go to 9                                                            >*/
    goto L9;
/*<    10 call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))                          >*/
L10:
    blf0_(&l, &c__0, n, &x[x_offset], &w[1], &cm[1], &sc[sc_offset], &nnt, &
        sc[mkp1 * sc_dim1 + 1]);
/*<       lbf=0                                                              >*/
    lbf = 0;
/*<       if(nnt .gt. nmin) go to 11                                         >*/
    if (nnt > nmin) {
    goto L11;
    }
/*<       call updpar(0,-1.d0)                                               >*/
    updpar_(&c__0, &c_b214);
/*<       go to 9                                                            >*/
    goto L9;
/*<    11 nep=0                                                              >*/
L11:
    nep = 0;
/*<       do 12 jp=1,p                                                       >*/
    i__1 = *p;
    for (jp = 1; jp <= i__1; ++jp) {
/*<       if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 12                      >*/
    if (x[mm[jp * mm_dim1 + 1] + jp * x_dim1] >= x[mm[*n + jp * mm_dim1] 
        + jp * x_dim1]) {
        goto L12;
    }
/*<       if(jf(l,jp,tb).ne.0) go to 12                                      >*/
    if (jf_(&l, &jp, &tb[6]) != 0) {
        goto L12;
    }
/*<       call isfac(l,jp,mm1,tb,cm,ja)                                      >*/
    isfac_(&l, &jp, &mm1, &tb[6], &cm[1], &ja);
/*<       if(ja.lt.0) go to 12                                               >*/
    if (ja < 0) {
        goto L12;
    }
/*<       if(.not.elg(jp,l,lx,tb,cm)) go to 12                               >*/
    if (! elg_(&jp, &l, &lx[1], &tb[6], &cm[1])) {
        goto L12;
    }
/*<       nep=nep+1                                                          >*/
    ++nep;
/*<    12 continue                                                           >*/
L12:
    ;
    }
/*<       if(nep .ne. 0) go to 13                                            >*/
    if (nep != 0) {
    goto L13;
    }
/*<       call updpar(0,-1.d0)                                               >*/
    updpar_(&c__0, &c_b214);
/*<       go to 9                                                            >*/
    goto L9;
/*<    13 call mnspan(ms,alf,nep,nnt,mn,me,mel)                              >*/
L13:
    mnspan_(ms, &alf, &nep, &nnt, &mn, &me, &mel);
/*<       if(nnt .gt. max0(me,mel)) go to 14                                 >*/
    if (nnt > max(me,mel)) {
    goto L14;
    }
/*<       call updpar(0,-1.d0)                                               >*/
    updpar_(&c__0, &c_b214);
/*<       go to 9                                                            >*/
    goto L9;
/*<    14 if(jq .ne. 0) go to 15                                             >*/
L14:
    if (jq != 0) {
    goto L15;
    }
/*<       jd1=1                                                              >*/
    jd1 = 1;
/*<       jd2=p                                                              >*/
    jd2 = *p;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 jd1=jq                                                             >*/
L15:
    jd1 = jq;
/*<       jd2=jd1                                                            >*/
    jd2 = jd1;
/*<    16 do 52 jp=jd1,jd2                                                   >*/
L16:
    i__1 = jd2;
    for (jp = jd1; jp <= i__1; ++jp) {
/*<       if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 52                      >*/
    if (x[mm[jp * mm_dim1 + 1] + jp * x_dim1] >= x[mm[*n + jp * mm_dim1] 
        + jp * x_dim1]) {
        goto L52;
    }
/*<       if(jf(l,jp,tb).ne.0) go to 52                                      >*/
    if (jf_(&l, &jp, &tb[6]) != 0) {
        goto L52;
    }
/*<       call isfac(l,jp,mm1,tb,cm,ja)                                      >*/
    isfac_(&l, &jp, &mm1, &tb[6], &cm[1], &ja);
/*<       if(ja.lt.0) go to 52                                               >*/
    if (ja < 0) {
        goto L52;
    }
/*<       if(.not.elg(jp,l,lx,tb,cm)) go to 52                               >*/
    if (! elg_(&jp, &l, &lx[1], &tb[6], &cm[1])) {
        goto L52;
    }
/*<       if(ja .ne. 0) go to 18                                             >*/
    if (ja != 0) {
        goto L18;
    }
/*<       if(lbf .eq. 0) go to 19                                            >*/
    if (lbf == 0) {
        goto L19;
    }
/*<       call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))                          >*/
    blf0_(&l, &c__0, n, &x[x_offset], &w[1], &cm[1], &sc[sc_offset], &nnt,
         &sc[mkp1 * sc_dim1 + 1]);
/*<       lbf=0                                                              >*/
    lbf = 0;
/*<       call mnspan(ms,alf,nep,nnt,mn,me,mel)                              >*/
    mnspan_(ms, &alf, &nep, &nnt, &mn, &me, &mel);
/*<       go to 19                                                           >*/
    goto L19;
/*<    18 call blf0(l,ja,n,x,w,cm,sc,nnt,sc(1,mkp1))                         >*/
L18:
    blf0_(&l, &ja, n, &x[x_offset], &w[1], &cm[1], &sc[sc_offset], &nnt, &
        sc[mkp1 * sc_dim1 + 1]);
/*<       lbf=1                                                              >*/
    lbf = 1;
/*<       if(nnt.le.nmin) go to 52                                           >*/
    if (nnt <= nmin) {
        goto L52;
    }
/*<       call mnspan(ms,alf,nep,nnt,mn,me,mel)                              >*/
    mnspan_(ms, &alf, &nep, &nnt, &mn, &me, &mel);
/*<       if(nnt.le.max0(me,mel)) go to 52                                   >*/
    if (nnt <= max(me,mel)) {
        goto L52;
    }
/*<    19 fvr=1.0                                                            >*/
L19:
    fvr = (float)1.;
/*<       if(jft(mm1,jp,tb).eq.0) fvr=1.0+fv                                 >*/
    if (jft_(&mm1, &jp, &tb[6]) == 0) {
        fvr = *fv + (float)1.;
    }
/*<       ict=0                                                              >*/
    ict = 0;
/*<       if(lx(jp) .ge. 0) go to 20                                         >*/
    if (lx[jp] >= 0) {
        goto L20;
    }
/*<       ict=1                                                              >*/
    ict = 1;
/*<       nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1                           >*/
    nc = (integer) (cm[(jp << 1) + 1] + (float).1) - (integer) (cm[jp * 2]
         + (float).1) + 1;
/*<    >*/
    csp_(&jp, &nc, &m, n, &x[x_offset], &y[1], &w[1], nk, &tb[6], &cm[1], 
        &kcp, &yb, &d__[d_offset], &kr, &nnt, &sw, &me, &mkp2, &nop, &
        sc[mkp1 * sc_dim1 + 1], &db[db_offset], &d__[d_dim1 * 3 + 1], 
        &mm[(*p + 1) * mm_dim1 + 1]);
/*<       if(nop.eq.0) go to 52                                              >*/
    if (nop == 0) {
        goto L52;
    }
/*<       go to 45                                                           >*/
    goto L45;
/*<    20 tb(2,m)=jp                                                         >*/
L20:
    tb[m * 5 + 2] = (real) jp;
/*<       tb(3,m)=x(mm(1,jp),jp)                                             >*/
    tb[m * 5 + 3] = x[mm[jp * mm_dim1 + 1] + jp * x_dim1];
/*<       tb(4,m)=l                                                          >*/
    tb[m * 5 + 4] = (real) l;
/*<       k1=kr                                                              >*/
    k1 = kr;
/*<       ssq=rsq                                                            >*/
    ssq = rsq;
/*<       call update(1,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  >*/
    update_(&c__1, n, &m, &kr, &x[x_offset], &y[1], &w[1], &sw, &yb, &tb[
        6], &cm[1], &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1], &db[
        db_offset], &d__[d_offset], &d__[d_dim1 * 3 + 1]);
/*<       if(kr .le. k1) go to 21                                            >*/
    if (kr <= k1) {
        goto L21;
    }
/*<       rsq=rsq-d(kr,1)**2                                                 >*/
/* Computing 2nd power */
    d__1 = d__[kr + d_dim1];
    rsq -= d__1 * d__1;
/*<       tb(1,m)=rsq/sw                                                     >*/
    tb[m * 5 + 1] = rsq / sw;
/*<       go to 22                                                           >*/
    goto L22;
/*<    21 tb(1,m)=big                                                        >*/
L21:
    tb[m * 5 + 1] = big;
/*<    >*/
L22:
    if (lx[jp] != 3 && (m < mk && nnt > me + mel)) {
        goto L26;
    }
/*<       tb(1,m)=rsq/sw                                                     >*/
    tb[m * 5 + 1] = rsq / sw;
/*<       newbf=newb(m,tb).eq.0                                              >*/
    newbf = newb_(&m, &tb[6]) == 0;
/*<       if(fvr*tb(1,m) .gt. txl .or. .not.(newbf)) go to 23                >*/
    if (fvr * tb[m * 5 + 1] > txl || ! newbf) {
        goto L23;
    }
/*<       txl=fvr*tb(1,m)                                                    >*/
    txl = fvr * tb[m * 5 + 1];
/*<       tx1=tb(1,m)                                                        >*/
    tx1 = tb[m * 5 + 1];
/*<       jq=jp                                                              >*/
    jq = jp;
/*<    23 if(fvr*tb(1,m) .gt. txi .or. .not.(newbf)) go to 25                >*/
L23:
    if (fvr * tb[m * 5 + 1] > txi || ! newbf) {
        goto L25;
    }
/*<       txi=fvr*tb(1,m)                                                    >*/
    txi = fvr * tb[m * 5 + 1];
/*<       tx(1)=tb(1,m)                                                      >*/
    tx[0] = tb[m * 5 + 1];
/*<       do 24 i=2,4                                                        >*/
    for (i__ = 2; i__ <= 4; ++i__) {
/*<       tx(i)=tb(i,m)                                                      >*/
        tx[i__ - 1] = tb[i__ + m * 5];
/*<    24 continue                                                           >*/
/* L24: */
    }
/*<       jas=ja                                                             >*/
    jas = ja;
/*<    25 kr=k1                                                              >*/
L25:
    kr = k1;
/*<       rsq=ssq                                                            >*/
    rsq = ssq;
/*<       go to 52                                                           >*/
    goto L52;
/*<    26 mm1=m                                                              >*/
L26:
    mm1 = m;
/*<       m=m+1                                                              >*/
    ++m;
/*<       tb(1,m)=big                                                        >*/
    tb[m * 5 + 1] = big;
/*<       xa=0.0                                                             >*/
    xa = (float)0.;
/*<       j=n                                                                >*/
    j = *n;
/*<       nnl=nnt                                                            >*/
    nnl = nnt;
/*<       nst=0                                                              >*/
    nst = 0;
/*<       nnr=-1                                                             >*/
    nnr = -1;
/*<    27 j0=j                                                               >*/
L27:
    j0 = j;
/*<    28 mj=mm(j,jp)                                                        >*/
L28:
    mj = mm[j + jp * mm_dim1];
/*<       h=sc(mj,mkp1)                                                      >*/
    h__ = sc[mj + mkp1 * sc_dim1];
/*<       if(w(mj) .le. 0.0 .or. h .le. 0.0) go to 29                        >*/
    if (w[mj] <= (float)0. || h__ <= (float)0.) {
        goto L29;
    }
/*<       nst=nst+1                                                          >*/
    ++nst;
/*<       nnl=nnl-1                                                          >*/
    --nnl;
/*<       nnr=nnr+1                                                          >*/
    ++nnr;
/*<    >*/
L29:
    if (x[mm[j - 1 + jp * mm_dim1] + jp * x_dim1] < x[mm[j + jp * mm_dim1]
         + jp * x_dim1] && nst >= mn && nnl >= mel && nnr >= me) {
        goto L30;
    }
/*<       j=j-1                                                              >*/
    --j;
/*<       if(j.le.1) go to 30                                                >*/
    if (j <= 1) {
        goto L30;
    }
/*<       go to 28                                                           >*/
    goto L28;
/*<    30 if(j.le.1) go to 45                                                >*/
L30:
    if (j <= 1) {
        goto L45;
    }
/*<       nst=0                                                              >*/
    nst = 0;
/*<       xb=xa                                                              >*/
    xb = xa;
/*<       xa=x(mm(j,jp),jp)                                                  >*/
    xa = x[mm[j + jp * mm_dim1] + jp * x_dim1];
/*<       if(j0 .ne. n) go to 34                                             >*/
    if (j0 != *n) {
        goto L34;
    }
/*<       v=0.d0                                                             >*/
    v = 0.;
/*<       u=v                                                                >*/
    u = v;
/*<       t=u                                                                >*/
    t = u;
/*<       we=t                                                               >*/
    we = t;
/*<       se=we                                                              >*/
    se = we;
/*<       sy=se                                                              >*/
    sy = se;
/*<       dy=sy                                                              >*/
    dy = sy;
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 32                                                           >*/
    goto L32;
/*<    31 i=i+1                                                              >*/
L31:
    ++i__;
/*<    32 if((i).gt.(kr)) go to 33                                           >*/
L32:
    if (i__ > kr) {
        goto L33;
    }
/*<       d(i,2)=0.d0                                                        >*/
    d__[i__ + (d_dim1 << 1)] = 0.;
/*<       d(i,3)=d(i,2)                                                      >*/
    d__[i__ + d_dim1 * 3] = d__[i__ + (d_dim1 << 1)];
/*<       go to 31                                                           >*/
    goto L31;
/*<    33 txt=x(mm(1,jp),jp)+x(mm(n,jp),jp)                                  >*/
L33:
    txt = x[mm[jp * mm_dim1 + 1] + jp * x_dim1] + x[mm[*n + jp * mm_dim1] 
        + jp * x_dim1];
/*<       xt=0.5*txt                                                         >*/
    xt = txt * (float).5;
/*<       go to 37                                                           >*/
    goto L37;
/*<    34 dx=xb-xa                                                           >*/
L34:
    dx = xb - xa;
/*<       dy=dy+dx*sy                                                        >*/
    dy += dx * sy;
/*<       we=we+dx*se                                                        >*/
    we += dx * se;
/*<       v=v+dx*(2.d0*u-(xb+xa-txt)*t)                                      >*/
    v += dx * (u * 2. - (xb + xa - txt) * t);
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 36                                                           >*/
    goto L36;
/*<    35 i=i+1                                                              >*/
L35:
    ++i__;
/*<    36 if((i).gt.(kr)) go to 37                                           >*/
L36:
    if (i__ > kr) {
        goto L37;
    }
/*<       d(i,2)=d(i,2)+dx*d(i,3)                                            >*/
    d__[i__ + (d_dim1 << 1)] += dx * d__[i__ + d_dim1 * 3];
/*<       go to 35                                                           >*/
    goto L35;
/*<    37 do 40 k=j,j0                                                       >*/
L37:
    i__2 = j0;
    for (k = j; k <= i__2; ++k) {
/*<       mj=mm(k,jp)                                                        >*/
        mj = mm[k + jp * mm_dim1];
/*<       h=sc(mj,mkp1)                                                      >*/
        h__ = sc[mj + mkp1 * sc_dim1];
/*<       if(w(mj).le.0.0.or.h.le.0.0) go to 40                              >*/
        if (w[mj] <= (float)0. || h__ <= (float)0.) {
        goto L40;
        }
/*<       xx=x(mj,jp)                                                        >*/
        xx = x[mj + jp * x_dim1];
/*<       xd=xx-xa                                                           >*/
        xd = xx - xa;
/*<       su=w(mj)*h                                                         >*/
        su = w[mj] * h__;
/*<       st=su*xd                                                           >*/
        st = su * xd;
/*<       yc=y(mj)-yb                                                        >*/
        yc = y[mj] - yb;
/*<       dy=dy+st*yc                                                        >*/
        dy += st * yc;
/*<       sy=sy+su*yc                                                        >*/
        sy += su * yc;
/*<       we=we+st                                                           >*/
        we += st;
/*<       se=se+su                                                           >*/
        se += su;
/*<       sj=w(mj)*h**2                                                      >*/
/* Computing 2nd power */
        r__1 = h__;
        sj = w[mj] * (r__1 * r__1);
/*<       v=v+sj*xd**2                                                       >*/
/* Computing 2nd power */
        d__1 = xd;
        v += sj * (d__1 * d__1);
/*<       t=t+sj                                                             >*/
        t += sj;
/*<       u=u+sj*(xx-xt)                                                     >*/
        u += sj * (xx - xt);
/*<       i=1                                                                >*/
        i__ = 1;
/*<       go to 39                                                           >*/
        goto L39;
/*<    38 i=i+1                                                              >*/
L38:
        ++i__;
/*<    39 if((i).gt.(kr)) go to 40                                           >*/
L39:
        if (i__ > kr) {
        goto L40;
        }
/*<       tt=db(mj,i)                                                        >*/
        tt = db[mj + i__ * db_dim1];
/*<       d(i,2)=d(i,2)+st*tt                                                >*/
        d__[i__ + (d_dim1 << 1)] += st * tt;
/*<       d(i,3)=d(i,3)+su*tt                                                >*/
        d__[i__ + d_dim1 * 3] += su * tt;
/*<       go to 38                                                           >*/
        goto L38;
/*<    40 continue                                                           >*/
L40:
        ;
    }
/*<       dv=v-we**2/sw                                                      >*/
/* Computing 2nd power */
    d__1 = we;
    dv = v - d__1 * d__1 / sw;
/*<       if(dv .le. 0.d0) go to 44                                          >*/
    if (dv <= 0.) {
        goto L44;
    }
/*<       a=0.d0                                                             >*/
    a = 0.;
/*<       b=a                                                                >*/
    b = a;
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 42                                                           >*/
    goto L42;
/*<    41 i=i+1                                                              >*/
L41:
    ++i__;
/*<    42 if((i).gt.(kr)) go to 43                                           >*/
L42:
    if (i__ > kr) {
        goto L43;
    }
/*<       s=d(i,2)                                                           >*/
    s = d__[i__ + (d_dim1 << 1)];
/*<       a=a+s*d(i,1)                                                       >*/
    a += s * d__[i__ + d_dim1];
/*<       b=b+s**2                                                           >*/
/* Computing 2nd power */
    d__1 = s;
    b += d__1 * d__1;
/*<       go to 41                                                           >*/
    goto L41;
/*<    43 b=dv-b                                                             >*/
L43:
    b = dv - b;
/*<       if(b .le. eps*dv) go to 44                                         >*/
    if (b <= eps * dv) {
        goto L44;
    }
/*<       b=-(dy-a)**2/b                                                     >*/
/* Computing 2nd power */
    d__1 = dy - a;
    b = -(d__1 * d__1) / b;
/*<       if(b .ge. tb(1,m)) go to 44                                        >*/
    if (b >= tb[m * 5 + 1]) {
        goto L44;
    }
/*<       tb(1,m)=b                                                          >*/
    tb[m * 5 + 1] = b;
/*<       tb(3,m)=xa                                                         >*/
    tb[m * 5 + 3] = xa;
/*<    44 j=j-1                                                              >*/
L44:
    --j;
/*<       if(j.le.1) go to 45                                                >*/
    if (j <= 1) {
        goto L45;
    }
/*<       go to 27                                                           >*/
    goto L27;
/*<    45 tb(2,m)=jp                                                         >*/
L45:
    tb[m * 5 + 2] = (real) jp;
/*<       tb(4,m)=l                                                          >*/
    tb[m * 5 + 4] = (real) l;
/*<       tb(1,m)=(rsq+tb(1,m))/sw                                           >*/
    tb[m * 5 + 1] = (rsq + tb[m * 5 + 1]) / sw;
/*<       if(ict .ne. 0 .or. tb(1,mm1) .gt. fln*tb(1,m)) go to 46            >*/
    if (ict != 0 || tb[mm1 * 5 + 1] > fln * tb[m * 5 + 1]) {
        goto L46;
    }
/*<       mp=mm1                                                             >*/
    mp = mm1;
/*<       go to 47                                                           >*/
    goto L47;
/*<    46 mp=m                                                               >*/
L46:
    mp = m;
/*<    47 newbf=newb(mp,tb).eq.0                                             >*/
L47:
    newbf = newb_(&mp, &tb[6]) == 0;
/*<       if(fvr*tb(1,mp) .ge. txl .or. .not.(newbf)) go to 48               >*/
    if (fvr * tb[mp * 5 + 1] >= txl || ! newbf) {
        goto L48;
    }
/*<       txl=fvr*tb(1,mp)                                                   >*/
    txl = fvr * tb[mp * 5 + 1];
/*<       tx1=tb(1,mp)                                                       >*/
    tx1 = tb[mp * 5 + 1];
/*<       jq=jp                                                              >*/
    jq = jp;
/*<    48 if(fvr*tb(1,mp) .ge. txi .or. .not.(newbf)) go to 51               >*/
L48:
    if (fvr * tb[mp * 5 + 1] >= txi || ! newbf) {
        goto L51;
    }
/*<       txi=fvr*tb(1,mp)                                                   >*/
    txi = fvr * tb[mp * 5 + 1];
/*<       tx(1)=tb(1,mp)                                                     >*/
    tx[0] = tb[mp * 5 + 1];
/*<       do 49 i=2,4                                                        >*/
    for (i__ = 2; i__ <= 4; ++i__) {
/*<       tx(i)=tb(i,mp)                                                     >*/
        tx[i__ - 1] = tb[i__ + mp * 5];
/*<    49 continue                                                           >*/
/* L49: */
    }
/*<       jas=ja                                                             >*/
    jas = ja;
/*<       if(ict .eq. 0) go to 51                                            >*/
    if (ict == 0) {
        goto L51;
    }
/*<       do 50 i=1,nc                                                       >*/
    i__2 = nc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       cm(kcp0+i)=cm(kcp+i)                                               >*/
        cm[kcp0 + i__] = cm[kcp + i__];
/*<    50 continue                                                           >*/
/* L50: */
    }
/*<       kcp=kcp0+nc                                                        >*/
    kcp = kcp0 + nc;
/*<       tx(3)=kcp0                                                         >*/
    tx[2] = (real) kcp0;
/*<    51 if(ict .ne. 0) go to 52                                            >*/
L51:
    if (ict != 0) {
        goto L52;
    }
/*<       m=mm1                                                              >*/
    m = mm1;
/*<       mm1=m-1                                                            >*/
    mm1 = m - 1;
/*<       kr=k1                                                              >*/
    kr = k1;
/*<       rsq=ssq                                                            >*/
    rsq = ssq;
/*<    52 continue                                                           >*/
L52:
    ;
    }
/*<       call updpar(jq,asq0-tx1)                                           >*/
    d__1 = asq0 - tx1;
    updpar_(&jq, &d__1);
/*<       go to 9                                                            >*/
    goto L9;
/*<    53 jp=tx(2)+.1                                                        >*/
L53:
    jp = tx[1] + (float).1;
/*<       call selpar(int(tx(4)+.1))                                         >*/
    i__1 = (integer) (tx[3] + (float).1);
    selpar_(&i__1);
/*<       if(cm(2*jp) .le. 0.) go to 54                                      >*/
    if (cm[jp * 2] <= (float)0.) {
    goto L54;
    }
/*<       nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1                           >*/
    nc = (integer) (cm[(jp << 1) + 1] + (float).1) - (integer) (cm[jp * 2] + (
        float).1) + 1;
/*<       kcp0=kcp0+nc                                                       >*/
    kcp0 += nc;
/*<    54 if(jas .le. 0) go to 60                                            >*/
L54:
    if (jas <= 0) {
    goto L60;
    }
/*<       call getnst(jas,cm,jn,kcp,cm(kcp0+1))                              >*/
    getnst_(&jas, &cm[1], &jn, &kcp, &cm[kcp0 + 1]);
/*<       tb(2,m)=jn                                                         >*/
    tb[m * 5 + 2] = (real) jn;
/*<       tb(3,m)=kcp0                                                       >*/
    tb[m * 5 + 3] = (real) kcp0;
/*<       kcp0=kcp0+kcp                                                      >*/
    kcp0 += kcp;
/*<       tb(4,m)=tx(4)                                                      >*/
    tb[m * 5 + 4] = tx[3];
/*<       k1=kr                                                              >*/
    k1 = kr;
/*<       call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))                            >*/
    i__1 = (integer) (tx[3] + (float).1);
    blf_(&i__1, n, &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1]);
/*<       tx(4)=m                                                            >*/
    tx[3] = (real) m;
/*<       call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  >*/
    update_(&c__2, n, &m, &kr, &x[x_offset], &y[1], &w[1], &sw, &yb, &tb[6], &
        cm[1], &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1], &db[db_offset], &
        d__[d_offset], &d__[d_dim1 * 3 + 1]);
/*<       if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    >*/
    if (kr > k1) {
/* Computing 2nd power */
    d__1 = d__[kr + d_dim1];
    rsq -= d__1 * d__1;
    }
/*<       call addpar(m)                                                     >*/
    addpar_(&m);
/*<       if(m .ge. mk) go to 58                                             >*/
    if (m >= mk) {
    goto L58;
    }
/*<       m=m+1                                                              >*/
    ++m;
/*<       tb(2,m)=-tb(2,m-1)                                                 >*/
    tb[m * 5 + 2] = -tb[(m - 1) * 5 + 2];
/*<       do 55 i=3,4                                                        >*/
    for (i__ = 3; i__ <= 4; ++i__) {
/*<       tb(i,m)=tb(i,m-1)                                                  >*/
    tb[i__ + m * 5] = tb[i__ + (m - 1) * 5];
/*<    55 continue                                                           >*/
/* L55: */
    }
/*<       if(ibfext(m,tb,cm) .eq. 0) go to 56                                >*/
    if (ibfext_(&m, &tb[6], &cm[1]) == 0) {
    goto L56;
    }
/*<       m=m-1                                                              >*/
    --m;
/*<       go to 58                                                           >*/
    goto L58;
/*<    56 do 57 i=1,n                                                        >*/
L56:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sc(i,m)=phi(m,i,n,x,tb,cm)                                         >*/
    sc[i__ + m * sc_dim1] = phi_(&m, &i__, n, &x[x_offset], &tb[6], &cm[1]
        );
/*<    57 continue                                                           >*/
/* L57: */
    }
/*<       call addpar(m)                                                     >*/
    addpar_(&m);
/*<    58 if(it .le. 0) go to 59                                             >*/
L58:
    if (*it <= 0) {
    goto L59;
    }
/*<       mp=m-1                                                             >*/
    mp = m - 1;
/*<       tcst=(nopt-1)*df1+kr+1.0                                           >*/
    tcst = (nopt - 1) * df1 + kr + (float)1.;
/*<       fjn=jn                                                             >*/
    fjn = (real) jn;
/*<       fkr=kr                                                             >*/
    fkr = (real) kr;
/*<       gcv=(rsq/sw)/(1.d0-tcst/wn)**2                                     >*/
/* Computing 2nd power */
    d__1 = 1. - tcst / wn;
    gcv = rsq / sw / (d__1 * d__1);
/*<       call holl(jn,cm,tb(3,m),hol)                                       >*/
    holl_(&jn, &cm[1], &tb[m * 5 + 3], hol, 28L);
/*    if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,fjn,hol,tb(4,m)       1
092*/
/*    if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,fjn,hol,tb(4,m)     1
093*/
/*<    59 mtot=m                                                             >*/
L59:
    mtot = m;
/*<       m=m+1                                                              >*/
    ++m;
/*<       if(m.gt.mk) go to 69                                               >*/
    if (m > mk) {
    goto L69;
    }
/*<    60 do 61 i=1,5                                                        >*/
L60:
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       tb(i,m)=tx(i)                                                      >*/
    tb[i__ + m * 5] = tx[i__ - 1];
/*<    61 continue                                                           >*/
/* L61: */
    }
/*<       k1=kr                                                              >*/
    k1 = kr;
/*<       call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))                            >*/
    i__1 = (integer) (tx[3] + (float).1);
    blf_(&i__1, n, &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1]);
/*<       call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  >*/
    update_(&c__2, n, &m, &kr, &x[x_offset], &y[1], &w[1], &sw, &yb, &tb[6], &
        cm[1], &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1], &db[db_offset], &
        d__[d_offset], &d__[d_dim1 * 3 + 1]);
/*<       if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    >*/
    if (kr > k1) {
/* Computing 2nd power */
    d__1 = d__[kr + d_dim1];
    rsq -= d__1 * d__1;
    }
/*<       call addpar(m)                                                     >*/
    addpar_(&m);
/*<    >*/
    if (m >= mk || cm[jp * 2] <= (float)0. && tx[2] <= x[mm[jp * mm_dim1 + 1] 
        + jp * x_dim1]) {
    goto L66;
    }
/*<       m=m+1                                                              >*/
    ++m;
/*<       do 62 i=1,4                                                        >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<       tb(i,m)=tx(i)                                                      >*/
    tb[i__ + m * 5] = tx[i__ - 1];
/*<    62 continue                                                           >*/
/* L62: */
    }
/*<       tb(2,m)=-tb(2,m)                                                   >*/
    tb[m * 5 + 2] = -tb[m * 5 + 2];
/*<       if(cm(2*jp) .le. 0.0) go to 64                                     >*/
    if (cm[jp * 2] <= (float)0.) {
    goto L64;
    }
/*<       do 63 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sc(i,m)=phi(m,i,n,x,tb,cm)                                         >*/
    sc[i__ + m * sc_dim1] = phi_(&m, &i__, n, &x[x_offset], &tb[6], &cm[1]
        );
/*<    63 continue                                                           >*/
/* L63: */
    }
/*<       go to 65                                                           >*/
    goto L65;
/*<    64 k1=kr                                                              >*/
L64:
    k1 = kr;
/*<       call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  >*/
    update_(&c__2, n, &m, &kr, &x[x_offset], &y[1], &w[1], &sw, &yb, &tb[6], &
        cm[1], &sc[sc_offset], &sc[mkp1 * sc_dim1 + 1], &db[db_offset], &
        d__[d_offset], &d__[d_dim1 * 3 + 1]);
/*<       if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    >*/
    if (kr > k1) {
/* Computing 2nd power */
    d__1 = d__[kr + d_dim1];
    rsq -= d__1 * d__1;
    }
/*<    65 call addpar(m)                                                     >*/
L65:
    addpar_(&m);
/*<    66 tcst=nopt*df1+kr+1.0                                               >*/
L66:
    tcst = nopt * df1 + kr + (float)1.;
/*<       if(it .le. 0) go to 68                                             >*/
    if (*it <= 0) {
    goto L68;
    }
/*<       mp=m-1                                                             >*/
    mp = m - 1;
/*<       jp=abs(tx(2))+.1                                                   >*/
    jp = dabs(tx[1]) + (float).1;
/*<       fkr=kr                                                             >*/
    fkr = (real) kr;
/*<       gcv=(rsq/sw)/(1.d0-tcst/wn)**2                                     >*/
/* Computing 2nd power */
    d__1 = 1. - tcst / wn;
    gcv = rsq / sw / (d__1 * d__1);
/*<       if(cm(2*jp) .le. 0.0) go to 67                                     >*/
    if (cm[jp * 2] <= (float)0.) {
    goto L67;
    }
/*<       call holl(jp,cm,tx(3),hol)                                         >*/
    holl_(&jp, &cm[1], &tx[2], hol, 28L);
/*    if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,tx(2),hol,tx(4)       1
129*/
/*    if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,tx(2),hol,tx(4)     1
130*/
/*<       go to 68                                                           >*/
    goto L68;
/*<    67 xk=xm(jp)+xs(jp)*tx(3)                                             >*/
L67:
    xk = xm[jp] + xs[jp] * tx[2];
/*    if(m.eq.mtot+1) write(it,93) m,gcv,fkr,tcst,tx(2),xk,tx(4)         1
133*/
/*    if(m.eq.mtot+2) write(it,94) m,mp,gcv,fkr,tcst,tx(2),xk,tx(4)      1
134*/
/*<    68 mtot=m                                                             >*/
L68:
    mtot = m;
/*<       go to 8                                                            >*/
    goto L8;
/*<    69 mk=min0(m,mk)                                                      >*/
L69:
    mk = min(m,mk);
/*<       m=mk+1                                                             >*/
    m = mk + 1;
/*<       k=m                                                                >*/
    k = m;
/*<       go to 71                                                           >*/
    goto L71;
/*<    70 k=k+1                                                              >*/
L70:
    ++k;
/*<    71 if((k).gt.(nk)) go to 72                                           >*/
L71:
    if (k > *nk) {
    goto L72;
    }
/*<       tb(1,k)=0.0                                                        >*/
    tb[k * 5 + 1] = (float)0.;
/*<       go to 70                                                           >*/
    goto L70;
/*<    72 call sscp(n,m,sc,y,w,yb,yv,sw,db,d)                                >*/
L72:
    sscp_(n, &m, &sc[sc_offset], &y[1], &w[1], &yb, &yv, &sw, &db[db_offset], 
        &d__[d_offset]);
/*<       call lsf1(db,m,d,yb,alr,b,d(1,2),a,d(1,3))                         >*/
    lsf1_(&db[db_offset], &m, &d__[d_offset], &yb, &alr, &b, &d__[(d_dim1 << 
        1) + 1], &a, &d__[d_dim1 * 3 + 1]);
/*<       nli=0                                                              >*/
    nli = 0;
/*<       do 73 k=1,mk                                                       >*/
    i__1 = mk;
    for (k = 1; k <= i__1; ++k) {
/*<       if(d(k,2).ne.0.d0) nli=nli+1                                       >*/
    if (d__[k + (d_dim1 << 1)] != 0.) {
        ++nli;
    }
/*<    73 continue                                                           >*/
/* L73: */
    }
/*<       df1=df1*nopt+nli                                                   >*/
    df1 = df1 * nopt + nli;
/*<       tcst=df1+1.0                                                       >*/
    tcst = df1 + (float)1.;
/*<       df1=df1/nli                                                        >*/
    df1 /= nli;
/*<       do 74 k=1,nk                                                       >*/
    i__1 = *nk;
    for (k = 1; k <= i__1; ++k) {
/*<       tb(5,k)=df1                                                        >*/
    tb[k * 5 + 5] = df1;
/*<    74 continue                                                           >*/
/* L74: */
    }
/*<       asm=(b/sw)/(1.d0-tcst/wn)**2                                       >*/
/* Computing 2nd power */
    d__1 = 1. - tcst / wn;
    asm__ = b / sw / (d__1 * d__1);
/*<       tcsts=tcst                                                         >*/
    tcsts = tcst;
/*<       az=a                                                               >*/
    *az = a;
/*<       do 75 k=1,mk                                                       >*/
    i__1 = mk;
    for (k = 1; k <= i__1; ++k) {
/*<       tb(1,k)=0.0                                                        >*/
    tb[k * 5 + 1] = (float)0.;
/*<       if(d(k,2).ne.0.d0) tb(1,k)=d(k,2)                                  >*/
    if (d__[k + (d_dim1 << 1)] != 0.) {
        tb[k * 5 + 1] = d__[k + (d_dim1 << 1)];
    }
/*<    75 continue                                                           >*/
/* L75: */
    }
/*<       if(ix .eq. 0) go to 81                                             >*/
    if (ix == 0) {
    goto L81;
    }
/*<       sc(1,1)=(cfac*nopt)/nli                                            >*/
    sc[sc_dim1 + 1] = cfac * nopt / nli;
/*<       sc(2,1)=wn                                                         >*/
    sc[sc_dim1 + 2] = wn;
/*<       sc(3,1)=yv                                                         >*/
    sc[sc_dim1 + 3] = yv;
/*<       sc(4,1)=yb                                                         >*/
    sc[sc_dim1 + 4] = yb;
/*<       do 80 k=nli,nk                                                     >*/
    i__1 = *nk;
    for (k = nli; k <= i__1; ++k) {
/*<       call array(k+4,n,i,j)                                              >*/
    i__2 = k + 4;
    array_(&i__2, n, &i__, &j);
/*<       sc(i,j)=b/sw                                                       >*/
    sc[i__ + j * sc_dim1] = b / sw;
/*<       k1=k*(nk+1)+3                                                      >*/
    k1 = k * (*nk + 1) + 3;
/*<       l=0                                                                >*/
    l = 0;
/*<       go to 77                                                           >*/
    goto L77;
/*<    76 l=l+1                                                              >*/
L76:
    ++l;
/*<    77 if((l).gt.(nk)) go to 80                                           >*/
L77:
    if (l > *nk) {
        goto L80;
    }
/*<       k1=k1+1                                                            >*/
    ++k1;
/*<       call array(k1,n,i,j)                                               >*/
    array_(&k1, n, &i__, &j);
/*<       if(l .ne. 0) go to 78                                              >*/
    if (l != 0) {
        goto L78;
    }
/*<       sc(i,j)=a                                                          >*/
    sc[i__ + j * sc_dim1] = a;
/*<       go to 76                                                           >*/
    goto L76;
/*<    78 if(l .le. mk) go to 79                                             >*/
L78:
    if (l <= mk) {
        goto L79;
    }
/*<       sc(i,j)=0.0                                                        >*/
    sc[i__ + j * sc_dim1] = (float)0.;
/*<       go to 76                                                           >*/
    goto L76;
/*<    79 sc(i,j)=d(l,2)                                                     >*/
L79:
    sc[i__ + j * sc_dim1] = d__[l + (d_dim1 << 1)];
/*<       go to 76                                                           >*/
    goto L76;
/*<    80 continue                                                           >*/
L80:
    ;
    }
/*<       call array((nk+1)**2+4,n,i,j)                                      >*/
/* Computing 2nd power */
    i__2 = *nk + 1;
    i__1 = i__2 * i__2 + 4;
    array_(&i__1, n, &i__, &j);
/*<       sc(i,j)=mk                                                         >*/
    sc[i__ + j * sc_dim1] = (real) mk;
/*<       kl=nli                                                             >*/
    kl = nli;
/*<    81 do 88 ll=2,nli                                                     >*/
L81:
    i__1 = nli;
    for (ll = 2; ll <= i__1; ++ll) {
/*<       call bkstp(db,m,d,yb,alr,b,d(1,2),a,k,d(1,3))                      >*/
    bkstp_(&db[db_offset], &m, &d__[d_offset], &yb, &alr, &b, &d__[(
        d_dim1 << 1) + 1], &a, &k, &d__[d_dim1 * 3 + 1]);
/*<       if(k.eq.0) go to 89                                                >*/
    if (k == 0) {
        goto L89;
    }
/*<       if(ix .eq. 0) go to 86                                             >*/
    if (ix == 0) {
        goto L86;
    }
/*<       call array(kl+3,n,i,j)                                             >*/
    i__2 = kl + 3;
    array_(&i__2, n, &i__, &j);
/*<       sc(i,j)=b/sw                                                       >*/
    sc[i__ + j * sc_dim1] = b / sw;
/*<       kl=kl-1                                                            >*/
    --kl;
/*<       k1=kl*(nk+1)+3                                                     >*/
    k1 = kl * (*nk + 1) + 3;
/*<       l=0                                                                >*/
    l = 0;
/*<       go to 83                                                           >*/
    goto L83;
/*<    82 l=l+1                                                              >*/
L82:
    ++l;
/*<    83 if((l).gt.(nk)) go to 86                                           >*/
L83:
    if (l > *nk) {
        goto L86;
    }
/*<       k1=k1+1                                                            >*/
    ++k1;
/*<       call array(k1,n,i,j)                                               >*/
    array_(&k1, n, &i__, &j);
/*<       if(l .ne. 0) go to 84                                              >*/
    if (l != 0) {
        goto L84;
    }
/*<       sc(i,j)=a                                                          >*/
    sc[i__ + j * sc_dim1] = a;
/*<       go to 82                                                           >*/
    goto L82;
/*<    84 if(l .le. mk) go to 85                                             >*/
L84:
    if (l <= mk) {
        goto L85;
    }
/*<       sc(i,j)=0.0                                                        >*/
    sc[i__ + j * sc_dim1] = (float)0.;
/*<       go to 82                                                           >*/
    goto L82;
/*<    85 sc(i,j)=d(l,2)                                                     >*/
L85:
    sc[i__ + j * sc_dim1] = d__[l + (d_dim1 << 1)];
/*<       go to 82                                                           >*/
    goto L82;
/*<    86 tcst=tcst-df1                                                      >*/
L86:
    tcst -= df1;
/*<       b=(b/sw)/(1.d0-tcst/wn)**2                                         >*/
/* Computing 2nd power */
    d__1 = 1. - tcst / wn;
    b = b / sw / (d__1 * d__1);
/*<       if(b.ge.asm) go to 88                                              >*/
    if (b >= asm__) {
        goto L88;
    }
/*<       asm=b                                                              >*/
    asm__ = b;
/*<       tcsts=tcst                                                         >*/
    tcsts = tcst;
/*<       az=a                                                               >*/
    *az = a;
/*<       do 87 i=1,mk                                                       >*/
    i__2 = mk;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       tb(1,i)=0.0                                                        >*/
        tb[i__ * 5 + 1] = (float)0.;
/*<       if(d(i,2).ne.0.d0) tb(1,i)=d(i,2)                                  >*/
        if (d__[i__ + (d_dim1 << 1)] != 0.) {
        tb[i__ * 5 + 1] = d__[i__ + (d_dim1 << 1)];
        }
/*<    87 continue                                                           >*/
/* L87: */
    }
/*<    88 continue                                                           >*/
L88:
    ;
    }
/*<    89 if(txm .gt. asm) go to 91                                          >*/
L89:
    if (txm > asm__) {
    goto L91;
    }
/*<       asm=txm                                                            >*/
    asm__ = txm;
/*<       tcsts=1.0                                                          >*/
    tcsts = (float)1.;
/*<       az=yb                                                              >*/
    *az = yb;
/*<       do 90 i=1,nk                                                       >*/
    i__1 = *nk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       tb(1,i)=0.0                                                        >*/
    tb[i__ * 5 + 1] = (float)0.;
/*<    90 continue                                                           >*/
/* L90: */
    }
/*<    91 if(it .le. 0) go to 92                                             >*/
L91:
    if (*it <= 0) {
    goto L92;
    }
/*    write(it,95)                                                       1
232*/
/*<       call coefpr(it,mk,az,tb,cm,xs)                                     >*/
    coefpr_(it, &mk, az, &tb[6], &cm[1], &xs[1]);
/*    write(it,96) asm,tcsts                                             1
234*/
/*<    92 return                                                             >*/
L92:
    return 0;
/*<       entry setfln(val)                                                  >*/

L_setfln:
/*<       fln=val                                                            >*/
    fln = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry setalf(val)                                                  >*/

L_setalf:
/*<       alf=val                                                            >*/
    alf = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry setmin(nal)                                                  >*/

L_setmin:
/*<       nmin=nal                                                           >*/
    nmin = *nal;
/*<       return                                                             >*/
    return 0;
/*<       entry setcta(val)                                                  >*/

L_setcta:
/*<       vcst(2)=val                                                        >*/
    vcst[1] = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry setctl(val)                                                  >*/

L_setctl:
/*<       vcst(3)=val                                                        >*/
    vcst[2] = *val;
/*<       return                                                             >*/
    return 0;
/*<       entry xvmrgo(nal)                                                  >*/

L_xvmrgo:
/*<       ix=nal                                                             >*/
    ix = *nal;
/*<       return                                                             >*/
    return 0;
/*<       entry setalr(val)                                                  >*/

L_setalr:
/*<       alr=val                                                            >*/
    alr = *val;
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L93: */
/*<    >*/
/* L94: */
/*<    95 format(/,' final model after backward stepwise elimination:')      >*/
/* L95: */
/*<    >*/
/* L96: */
/*<    >*/
/* L97: */
/*<    98 format('   ',i3,'    ',g12.4,2('   ',f5.1))                        >*/
/* L98: */
/*<    >*/
/* L99: */
/*<    >*/
/* L100: */
/*<       end                                                                >*/
} /* marsgo_ */

/* Subroutine */ int marsgo_(n, p, x, y, w, nk, ms, df, fv, mi, lx, it, xm, 
    xs, az, tb, cm, sc, db, d__, mm)
integer *n, *p;
real *x, *y, *w;
integer *nk, *ms;
real *df, *fv;
integer *mi, *lx, *it;
real *xm, *xs, *az, *tb, *cm, *sc;
doublereal *db, *d__;
integer *mm;
{
    return marsgo_0_(0, n, p, x, y, w, nk, ms, df, fv, mi, lx, it, xm, xs, az,
         tb, cm, sc, db, d__, mm, (real *)0, (integer *)0);
    }

/* Subroutine */ int setfln_(val)
real *val;
{
    return marsgo_0_(1, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, val, (integer *)0);
    }

/* Subroutine */ int setalf_(val)
real *val;
{
    return marsgo_0_(2, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, val, (integer *)0);
    }

/* Subroutine */ int setmin_(nal)
integer *nal;
{
    return marsgo_0_(3, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, (real *)0, nal);
    }

/* Subroutine */ int setcta_(val)
real *val;
{
    return marsgo_0_(4, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, val, (integer *)0);
    }

/* Subroutine */ int setctl_(val)
real *val;
{
    return marsgo_0_(5, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, val, (integer *)0);
    }

/* Subroutine */ int xvmrgo_(nal)
integer *nal;
{
    return marsgo_0_(6, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, (real *)0, nal);
    }

/* Subroutine */ int setalr_(val)
real *val;
{
    return marsgo_0_(7, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (
        doublereal *)0, (integer *)0, val, (integer *)0);
    }

/*<       subroutine addpar (ib)                                             >*/
/* Subroutine */ int addpar_0_(n__, ib, l, jq, val, iarg, arg)
int n__;
integer *ib, *l, *jq;
doublereal *val;
integer *iarg;
real *arg;
{
    /* Initialized data */

    static real big = (float)9.9e30;
    static integer mpr = 10;
    static integer mtr = 5;
    static real beta = (float)1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    static integer i__, j, k, m[1000], n[1000];
    extern /* Subroutine */ int psort_();
    static integer jp[2000] /* was [2][1000] */, kp, lq;
    static real sp[1000], que[2000] /* was [2][1000] */;
    static integer itr, ktr;

/*<       parameter(maxdph=1000)                                             >*/
/*<       real que(2,maxdph),sp(maxdph)                                      >*/
/*<       integer m(maxdph),n(maxdph),jp(2,maxdph)                           >*/
/*<       double precision val                                               >*/
/*<       save mpr,mtr,ktr,big,beta,lq,kp,itr,jp,que,n,m                     >*/
/*<       data big,mpr,mtr,beta /9.9e30,10,5,1.0/                            >*/
    switch(n__) {
    case 1: goto L_nxtpar;
    case 2: goto L_updpar;
    case 3: goto L_selpar;
    case 4: goto L_itrpar;
    case 5: goto L_setmpr;
    case 6: goto L_setbta;
    case 7: goto L_setfrq;
    }

/*<       if(ib .ne. 0) go to 1                                              >*/
    if (*ib != 0) {
    goto L1;
    }
/*<       lq=1                                                               >*/
    lq = 1;
/*<       que(1,1)=big                                                       >*/
    que[0] = big;
/*<       que(2,1)=0.0                                                       >*/
    que[1] = (float)0.;
/*<       m(1)=1                                                             >*/
    m[0] = 1;
/*<       kp=0                                                               >*/
    kp = 0;
/*<       itr=kp                                                             >*/
    itr = kp;
/*<       n(1)=itr                                                           >*/
    n[0] = itr;
/*<       jp(1,1)=n(1)                                                       >*/
    jp[0] = n[0];
/*<       jp(2,1)=jp(1,1)                                                    >*/
    jp[1] = jp[0];
/*<       ktr=0.5*(mpr-1)+.1                                                 >*/
    ktr = (mpr - 1) * (float).5 + (float).1;
/*<       return                                                             >*/
    return 0;
/*<     1 if(que(1,lq).ge.-0.5) go to 2                                      >*/
L1:
    if (que[(lq << 1) - 2] >= (float)-.5) {
    goto L2;
    }
/*<       lq=lq-1                                                            >*/
    --lq;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 i=1                                                                >*/
L2:
    i__ = 1;
/*<       go to 4                                                            >*/
    goto L4;
/*<     3 i=i+1                                                              >*/
L3:
    ++i__;
/*<     4 if((i).gt.(lq)) go to 7                                            >*/
L4:
    if (i__ > lq) {
    goto L7;
    }
/*<       if(que(1,i).ge.-0.5) go to 3                                       >*/
    if (que[(i__ << 1) - 2] >= (float)-.5) {
    goto L3;
    }
/*<       lq=lq-1                                                            >*/
    --lq;
/*<       do 6 j=i,lq                                                        >*/
    i__1 = lq;
    for (j = i__; j <= i__1; ++j) {
/*<       n(j)=n(j+1)                                                        >*/
    n[j - 1] = n[j];
/*<       do 5 k=1,2                                                         >*/
    for (k = 1; k <= 2; ++k) {
/*<       jp(k,j)=jp(k,j+1)                                                  >*/
        jp[k + (j << 1) - 3] = jp[k + (j + 1 << 1) - 3];
/*<       que(k,j)=que(k,j+1)                                                >*/
        que[k + (j << 1) - 3] = que[k + (j + 1 << 1) - 3];
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       i=i-1                                                              >*/
    --i__;
/*<       go to 3                                                            >*/
    goto L3;
/*<     7 lq=lq+1                                                            >*/
L7:
    ++lq;
/*<       if(lq .le. maxdph) go to 8                                         >*/
    if (lq <= 1000) {
    goto L8;
    }
/*    write(6, '('' increase parameter maxdph in subroutine addpar to '' 1
312*/
/*   1)')                                                                1
313*/
/*    write(6,'('' '',i10,''  or larger, and recompile mars.'')') lq     1
314*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     8 que(1,lq)=big                                                      >*/
L8:
    que[(lq << 1) - 2] = big;
/*<       que(2,lq)=0.0                                                      >*/
    que[(lq << 1) - 1] = (float)0.;
/*<       n(lq)=ib                                                           >*/
    n[lq - 1] = *ib;
/*<       jp(1,lq)=0                                                         >*/
    jp[(lq << 1) - 2] = 0;
/*<       jp(2,lq)=jp(1,lq)                                                  >*/
    jp[(lq << 1) - 1] = jp[(lq << 1) - 2];
/*<       do 9 i=1,lq                                                        >*/
    i__1 = lq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       m(i)=i                                                             >*/
    m[i__ - 1] = i__;
/*<       sp(i)=que(1,i)                                                     >*/
    sp[i__ - 1] = que[(i__ << 1) - 2];
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       call psort(sp,m,1,lq)                                              >*/
    psort_(sp, m, &c__1, &lq);
/*<       do 10 i=1,lq                                                       >*/
    i__1 = lq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       j=m(i)                                                             >*/
    j = m[i__ - 1];
/*<       sp(j)=i+beta*(itr-que(2,j))                                        >*/
    sp[j - 1] = i__ + beta * (itr - que[(j << 1) - 1]);
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       call psort(sp,m,1,lq)                                              >*/
    psort_(sp, m, &c__1, &lq);
/*<       kp=max0(0,lq-mpr)                                                  >*/
/* Computing MAX */
    i__1 = 0, i__2 = lq - mpr;
    kp = max(i__1,i__2);
/*<       return                                                             >*/
    return 0;
/*<       entry nxtpar (l,jq)                                                >*/

L_nxtpar:
/*<       kp=kp+1                                                            >*/
    ++kp;
/*<       if(kp .le. lq) go to 11                                            >*/
    if (kp <= lq) {
    goto L11;
    }
/*<       l=-1                                                               >*/
    *l = -1;
/*<       return                                                             >*/
    return 0;
/*<    11 l=n(m(kp))                                                         >*/
L11:
    *l = n[m[kp - 1] - 1];
/*<       if(itr-jp(2,m(kp)).gt.mtr.or.itr.le.ktr) jp(1,m(kp))=0             >*/
    if (itr - jp[(m[kp - 1] << 1) - 1] > mtr || itr <= ktr) {
    jp[(m[kp - 1] << 1) - 2] = 0;
    }
/*<       jq=jp(1,m(kp))                                                     >*/
    *jq = jp[(m[kp - 1] << 1) - 2];
/*<       return                                                             >*/
    return 0;
/*<       entry updpar (jq,val)                                              >*/

L_updpar:
/*<       que(1,m(kp))=val                                                   >*/
    que[(m[kp - 1] << 1) - 2] = *val;
/*<       que(2,m(kp))=itr                                                   >*/
    que[(m[kp - 1] << 1) - 1] = (real) itr;
/*<       if(jp(1,m(kp)) .ne. 0) go to 12                                    >*/
    if (jp[(m[kp - 1] << 1) - 2] != 0) {
    goto L12;
    }
/*<       jp(1,m(kp))=jq                                                     >*/
    jp[(m[kp - 1] << 1) - 2] = *jq;
/*<       jp(2,m(kp))=itr                                                    >*/
    jp[(m[kp - 1] << 1) - 1] = itr;
/*<    12 return                                                             >*/
L12:
    return 0;
/*<       entry selpar(ib)                                                   >*/

L_selpar:
/*<       do 13 i=lq,1,-1                                                    >*/
    for (i__ = lq; i__ >= 1; --i__) {
/*<       if(n(i).ne.ib) go to 13                                            >*/
    if (n[i__ - 1] != *ib) {
        goto L13;
    }
/*<       jp(1,i)=0                                                          >*/
    jp[(i__ << 1) - 2] = 0;
/*<       go to 14                                                           >*/
    goto L14;
/*<    13 continue                                                           >*/
L13:
    ;
    }
/*<    14 return                                                             >*/
L14:
    return 0;
/*<       entry itrpar(iarg)                                                 >*/

L_itrpar:
/*<       itr=iarg                                                           >*/
    itr = *iarg;
/*<       return                                                             >*/
    return 0;
/*<       entry setmpr(iarg)                                                 >*/

L_setmpr:
/*<       mpr=iarg                                                           >*/
    mpr = *iarg;
/*<       return                                                             >*/
    return 0;
/*<       entry setbta(arg)                                                  >*/

L_setbta:
/*<       beta=arg                                                           >*/
    beta = *arg;
/*<       return                                                             >*/
    return 0;
/*<       entry setfrq(arg)                                                  >*/

L_setfrq:
/*<       mtr=1.0/amax1(arg,0.01)+.1                                         >*/
    mtr = (float)1. / dmax(*arg,(float).01) + (float).1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* addpar_ */

/* Subroutine */ int addpar_(ib)
integer *ib;
{
    return addpar_0_(0, ib, (integer *)0, (integer *)0, (doublereal *)0, (
        integer *)0, (real *)0);
    }

/* Subroutine */ int nxtpar_(l, jq)
integer *l, *jq;
{
    return addpar_0_(1, (integer *)0, l, jq, (doublereal *)0, (integer *)0, (
        real *)0);
    }

/* Subroutine */ int updpar_(jq, val)
integer *jq;
doublereal *val;
{
    return addpar_0_(2, (integer *)0, (integer *)0, jq, val, (integer *)0, (
        real *)0);
    }

/* Subroutine */ int selpar_(ib)
integer *ib;
{
    return addpar_0_(3, ib, (integer *)0, (integer *)0, (doublereal *)0, (
        integer *)0, (real *)0);
    }

/* Subroutine */ int itrpar_(iarg)
integer *iarg;
{
    return addpar_0_(4, (integer *)0, (integer *)0, (integer *)0, (doublereal 
        *)0, iarg, (real *)0);
    }

/* Subroutine */ int setmpr_(iarg)
integer *iarg;
{
    return addpar_0_(5, (integer *)0, (integer *)0, (integer *)0, (doublereal 
        *)0, iarg, (real *)0);
    }

/* Subroutine */ int setbta_(arg)
real *arg;
{
    return addpar_0_(6, (integer *)0, (integer *)0, (integer *)0, (doublereal 
        *)0, (integer *)0, arg);
    }

/* Subroutine */ int setfrq_(arg)
real *arg;
{
    return addpar_0_(7, (integer *)0, (integer *)0, (integer *)0, (doublereal 
        *)0, (integer *)0, arg);
    }

/*<       subroutine speed(is)                                               >*/
/* Subroutine */ int speed_(is)
integer *is;
{
    /* Initialized data */

    static integer lque[5] = { 9999,20,20,10,5 };
    static real freq[5] = { (float)9e30,(float)9e30,(float).2,(float).2,(
        float).2 };

    static integer j;
    extern /* Subroutine */ int setfrq_(), setmpr_();

/*<       integer lque(5)                                                    >*/
/*<       real freq(5)                                                       >*/
/*<       save lque,freq                                                     >*/
/*<       data lque /9999,20,20,10,5/                                        >*/
/*<       data freq /9.e30,9.e30,0.2,0.2,0.2/                                >*/
/*<       j=is                                                               >*/
    j = *is;
/*<       if(is.lt.1) j=1                                                    >*/
    if (*is < 1) {
    j = 1;
    }
/*<       if(is.gt.5) j=5                                                    >*/
    if (*is > 5) {
    j = 5;
    }
/*<       call setmpr(lque(j))                                               >*/
    setmpr_(&lque[j - 1]);
/*<       call setfrq(freq(j))                                               >*/
    setfrq_(&freq[j - 1]);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* speed_ */

/*<       subroutine atoscl(n,p,w,x,lx,mm,xm,xs,cm,z)                        >*/
/* Subroutine */ int atoscl_(n, p, w, x, lx, mm, xm, xs, cm, z__)
integer *n, *p;
real *w, *x;
integer *lx, *mm;
real *xm, *xs, *cm, *z__;
{
    /* System generated locals */
    integer mm_dim1, mm_offset, x_dim1, x_offset, z_dim1, z_offset, i__1, 
        i__2;
    real r__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer j0, n2, nc, ip;
    static doublereal sw;
    extern /* Subroutine */ int stcmrs_(), stfmrs_();
    static integer n2p1, nct;

/*<       integer p,lx(p),mm(n,p)                                            >*/
/*<       real w(n),x(n,p),z(n,p),xm(p),xs(p),cm(*)                          >*/
/*<       double precision s,t,sw                                            >*/
/*<       sw=0.d0                                                            >*/
    /* Parameter adjustments */
    --w;
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --xs;
    --xm;
    mm_dim1 = *n;
    mm_offset = mm_dim1 + 1;
    mm -= mm_offset;
    --lx;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --cm;

    /* Function Body */
    sw = 0.;
/*<       do 1 j=1,n                                                         >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<       sw=sw+w(j)                                                         >*/
    sw += w[j];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       ip=0                                                               >*/
    ip = 0;
/*<       nct=ip                                                             >*/
    nct = ip;
/*<       do 12 i=1,p                                                        >*/
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(lx(i) .ne. 0) go to 2                                           >*/
    if (lx[i__] != 0) {
        goto L2;
    }
/*<       xm(i)=0.0                                                          >*/
    xm[i__] = (float)0.;
/*<       xs(i)=xm(i)                                                        >*/
    xs[i__] = xm[i__];
/*<       go to 12                                                           >*/
    goto L12;
/*<     2 if(lx(i) .ge. 0) go to 8                                           >*/
L2:
    if (lx[i__] >= 0) {
        goto L8;
    }
/*<       nc=0                                                               >*/
    nc = 0;
/*<       xm(i)=ip                                                           >*/
    xm[i__] = (real) ip;
/*<       j=1                                                                >*/
    j = 1;
/*<       nct=nct+1                                                          >*/
    ++nct;
/*<     3 j0=j                                                               >*/
L3:
    j0 = j;
/*<       if(j .ge. n) go to 5                                               >*/
    if (j >= *n) {
        goto L5;
    }
/*<     4 if(x(mm(j+1,i),i).gt.x(mm(j,i),i)) go to 5                         >*/
L4:
    if (x[mm[j + 1 + i__ * mm_dim1] + i__ * x_dim1] > x[mm[j + i__ * 
        mm_dim1] + i__ * x_dim1]) {
        goto L5;
    }
/*<       j=j+1                                                              >*/
    ++j;
/*<       if(j.ge.n) go to 5                                                 >*/
    if (j >= *n) {
        goto L5;
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<     5 ip=ip+1                                                            >*/
L5:
    ++ip;
/*<       cm(ip)=x(mm(j,i),i)                                                >*/
    cm[ip] = x[mm[j + i__ * mm_dim1] + i__ * x_dim1];
/*<       nc=nc+1                                                            >*/
    ++nc;
/*<       do 6 k=j0,j                                                        >*/
    i__2 = j;
    for (k = j0; k <= i__2; ++k) {
/*<       z(mm(k,i),i)=nc                                                    >*/
        z__[mm[k + i__ * mm_dim1] + i__ * z_dim1] = (real) nc;
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       j=j+1                                                              >*/
    ++j;
/*<       if(j.gt.n) go to 7                                                 >*/
    if (j > *n) {
        goto L7;
    }
/*<       go to 3                                                            >*/
    goto L3;
/*<     7 xs(i)=nc                                                           >*/
L7:
    xs[i__] = (real) nc;
/*<       go to 12                                                           >*/
    goto L12;
/*<     8 s=0.d0                                                             >*/
L8:
    s = 0.;
/*<       t=s                                                                >*/
    t = s;
/*<       do 9 j=1,n                                                         >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       s=s+w(j)*x(j,i)                                                    >*/
        s += w[j] * x[j + i__ * x_dim1];
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       xm(i)=s                                                            >*/
    xm[i__] = s;
/*<       do 10 j=1,n                                                        >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       z(j,i)=x(j,i)-s                                                    >*/
        z__[j + i__ * z_dim1] = x[j + i__ * x_dim1] - s;
/*<       t=t+w(j)*z(j,i)**2                                                 >*/
/* Computing 2nd power */
        r__1 = z__[j + i__ * z_dim1];
        t += w[j] * (r__1 * r__1);
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       xs(i)=1.0                                                          >*/
    xs[i__] = (float)1.;
/*<       if(t.le.0.d0) go to 12                                             >*/
    if (t <= 0.) {
        goto L12;
    }
/*<       t=dsqrt(t/sw)                                                      >*/
    t = sqrt(t / sw);
/*<       xs(i)=t                                                            >*/
    xs[i__] = t;
/*<       t=1.d0/t                                                           >*/
    t = 1. / t;
/*<       do 11 j=1,n                                                        >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       z(j,i)=t*z(j,i)                                                    >*/
        z__[j + i__ * z_dim1] = t * z__[j + i__ * z_dim1];
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<    12 continue                                                           >*/
L12:
    ;
    }
/*<       n2=2*p+1                                                           >*/
    n2 = (*p << 1) + 1;
/*<       if(nct .ne. 0) go to 14                                            >*/
    if (nct != 0) {
    goto L14;
    }
/*<       do 13 i=1,n2                                                       >*/
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       cm(i)=0.0                                                          >*/
    cm[i__] = (float)0.;
/*<    13 continue                                                           >*/
/* L13: */
    }
/*<       return                                                             >*/
    return 0;
/*<    14 n2p1=n2+1                                                          >*/
L14:
    n2p1 = n2 + 1;
/*<       i=ip                                                               >*/
    i__ = ip;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 i=i+(-1)                                                           >*/
L15:
    --i__;
/*<    16 if((-1)*((i)-(1)).gt.0) go to 17                                   >*/
L16:
    if (-(i__ - 1) > 0) {
    goto L17;
    }
/*<       cm(i+n2)=cm(i)                                                     >*/
    cm[i__ + n2] = cm[i__];
/*<       go to 15                                                           >*/
    goto L15;
/*<    17 j=0                                                                >*/
L17:
    j = 0;
/*<       i=2                                                                >*/
    i__ = 2;
/*<       go to 19                                                           >*/
    goto L19;
/*<    18 i=i+(2)                                                            >*/
L18:
    i__ += 2;
/*<    19 if((2)*((i)-(n2)).gt.0) go to 22                                   >*/
L19:
    if (i__ - n2 << 1 > 0) {
    goto L22;
    }
/*<       j=j+1                                                              >*/
    ++j;
/*<       if(lx(j) .ge. 0) go to 20                                          >*/
    if (lx[j] >= 0) {
    goto L20;
    }
/*<       cm(i)=xm(j)+n2p1                                                   >*/
    cm[i__] = xm[j] + n2p1;
/*<       cm(i+1)=cm(i)+xs(j)-1.0                                            >*/
    cm[i__ + 1] = cm[i__] + xs[j] - (float)1.;
/*<       go to 18                                                           >*/
    goto L18;
/*<    20 cm(i)=0.0                                                          >*/
L20:
    cm[i__] = (float)0.;
/*<       cm(i+1)=cm(i)                                                      >*/
    cm[i__ + 1] = cm[i__];
/*<       go to 18                                                           >*/
    goto L18;
/*<    22 cm(1)=nct                                                          >*/
L22:
    cm[1] = (real) nct;
/*<       call stfmrs(1)                                                     >*/
    stfmrs_(&c__1);
/*<       call stcmrs(1)                                                     >*/
    stcmrs_(&c__1);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* atoscl_ */

/*<       subroutine orgpl (xm,xs,nk,tb,cm)                                  >*/
/* Subroutine */ int orgpl_(xm, xs, nk, tb, cm)
real *xm, *xs;
integer *nk;
real *tb, *cm;
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer j, m, ip;
    static real scl;

/*<       real xm(*),xs(*),tb(5,nk),cm(*)                                    >*/
/*<       do 1 m=1,nk                                                        >*/
    /* Parameter adjustments */
    --xm;
    --xs;
    tb -= 6;
    --cm;

    /* Function Body */
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       j=abs(tb(2,m))+.1                                                  >*/
    j = (r__1 = tb[m * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(cm(2*j).gt.0.0) go to 1                                         >*/
    if (cm[j * 2] > (float)0.) {
        goto L1;
    }
/*<       tb(3,m)=xm(j)+xs(j)*tb(3,m)                                        >*/
    tb[m * 5 + 3] = xm[j] + xs[j] * tb[m * 5 + 3];
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<       do 4 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 4                                         >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L4;
    }
/*<       scl=1.0                                                            >*/
    scl = (float)1.;
/*<       ip=m                                                               >*/
    ip = m;
/*<     2 if(ip.le.0) go to 3                                                >*/
L2:
    if (ip <= 0) {
        goto L3;
    }
/*<       j=abs(tb(2,ip))+.1                                                 >*/
    j = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(cm(2*j).eq.0.0) scl=scl*xs(j)                                   >*/
    if (cm[j * 2] == (float)0.) {
        scl *= xs[j];
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 2                                                            >*/
    goto L2;
/*<     3 tb(1,m)=tb(1,m)/scl                                                >*/
L3:
    tb[m * 5 + 1] /= scl;
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* orgpl_ */

/*<       subroutine anova (n,x,y,w,nk,it,tb,cm,lp,lv,t,d)                   >*/
/* Subroutine */ int anova_(n, x, y, w, nk, it, tb, cm, lp, lv, t, d__)
integer *n;
real *x, *y, *w;
integer *nk, *it;
real *tb, *cm;
integer *lp, *lv;
real *t;
doublereal *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, t_dim1, t_offset, d_dim1, d_offset, i__1, i__2, 
        i__3;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int exch_(), coll_();
    extern doublereal varf_();
    extern integer nord_();
    static integer i__, j, k, l, m;
    static doublereal s, u;
    static integer i2, k2, na;
    extern integer jf_();
    static integer ni, ll, lm;
    static doublereal yb;
    static integer np, im;
    static doublereal wn, sw, yv;
    static real efm;
    extern doublereal efp_();
    static real eft;
    extern doublereal phi_();
    extern /* Subroutine */ int lsf_();
    static integer nim1, nkp1, nkp2, nkp3;

/*<       integer lp(3,*),lv(*)                                              >*/
/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),t(n,nk)                       >*/
/*<       double precision d(nk,*),s,u,sw,yv,wn,yb                           >*/
/*<       if(it.le.0) return                                                 >*/
    /* Parameter adjustments */
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    t_dim1 = *n;
    t_offset = t_dim1 + 1;
    t -= t_offset;
    tb -= 6;
    --cm;
    lp -= 4;
    --lv;

    /* Function Body */
    if (*it <= 0) {
    return 0;
    }
/*<       nkp1=nk+1                                                          >*/
    nkp1 = *nk + 1;
/*<       nkp2=nkp1+1                                                        >*/
    nkp2 = nkp1 + 1;
/*<       nkp3=nkp2+1                                                        >*/
    nkp3 = nkp2 + 1;
/*<       lm=nkp3+nk                                                         >*/
    lm = nkp3 + *nk;
/*<       sw=0.d0                                                            >*/
    sw = 0.;
/*<       wn=sw                                                              >*/
    wn = sw;
/*<       s=wn                                                               >*/
    s = wn;
/*<       u=s                                                                >*/
    u = s;
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<       s=s+w(i)*y(i)                                                      >*/
    s += w[i__] * y[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       yb=s                                                               >*/
    yb = s;
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    d__1 = sw;
    wn = d__1 * d__1 / wn;
/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       u=u+w(i)*(y(i)-s)**2                                               >*/
/* Computing 2nd power */
    d__1 = y[i__] - s;
    u += w[i__] * (d__1 * d__1);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       yv=u/sw                                                            >*/
    yv = u / sw;
/*<       eft=1.0                                                            >*/
    eft = (float)1.;
/*<       do 3 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).ne.0.0) eft=eft+tb(5,m)                                 >*/
    if (tb[m * 5 + 1] != (float)0.) {
        eft += tb[m * 5 + 5];
    }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       ni=0                                                               >*/
    ni = 0;
/*<       do 9 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 9                                         >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L9;
    }
/*<       ni=ni+1                                                            >*/
    ++ni;
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 4 j=1,n                                                         >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       t(j,ni)=phi(m,j,n,x,tb,cm)                                         >*/
        t[j + ni * t_dim1] = phi_(&m, &j, n, &x[x_offset], &tb[6], &cm[1])
            ;
/*<       s=s+w(j)*t(j,ni)                                                   >*/
        s += w[j] * t[j + ni * t_dim1];
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       do 5 j=1,n                                                         >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       t(j,ni)=t(j,ni)-s                                                  >*/
        t[j + ni * t_dim1] -= s;
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       do 7 i=1,ni                                                        >*/
    i__2 = ni;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=0.d0                                                             >*/
        s = 0.;
/*<       do 6 j=1,n                                                         >*/
        i__3 = *n;
        for (j = 1; j <= i__3; ++j) {
/*<       s=s+w(j)*t(j,i)*t(j,ni)                                            >*/
        s += w[j] * t[j + i__ * t_dim1] * t[j + ni * t_dim1];
/*<     6 continue                                                           >*/
/* L6: */
        }
/*<       d(i,ni)=s                                                          >*/
        d__[i__ + ni * d_dim1] = s;
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 8 j=1,n                                                         >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       s=s+w(j)*t(j,ni)*(y(j)-yb)                                         >*/
        s += w[j] * t[j + ni * t_dim1] * (y[j] - yb);
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       d(ni,nkp1)=s                                                       >*/
    d__[ni + nkp1 * d_dim1] = s;
/*<       d(ni,nkp2)=tb(1,m)                                                 >*/
    d__[ni + nkp2 * d_dim1] = tb[m * 5 + 1];
/*<       lv(ni)=m                                                           >*/
    lv[ni] = m;
/*<     9 continue                                                           >*/
L9:
    ;
    }
/*<       if(ni .ne. 0) go to 10                                             >*/
    if (ni != 0) {
    goto L10;
    }
/*    write(it,26)                                                       1
548*/
/*<       return                                                             >*/
    return 0;
/*<    10 do 11 m=1,ni                                                       >*/
L10:
    i__1 = ni;
    for (m = 1; m <= i__1; ++m) {
/*<       t(m,1)=lv(m)                                                       >*/
    t[m + t_dim1] = (real) lv[m];
/*<    11 continue                                                           >*/
/* L11: */
    }
/*    write(it,24) ni                                                    1
553*/
/*<       call coll(nk,tb,lp,lv,lp(1,nkp1))                                  >*/
    coll_(nk, &tb[6], &lp[4], &lv[1], &lp[nkp1 * 3 + 1]);
/*<       m=1                                                                >*/
    m = 1;
/*<    12 if(lp(1,m).eq.0) go to 13                                          >*/
L12:
    if (lp[m * 3 + 1] == 0) {
    goto L13;
    }
/*<       m=m+1                                                              >*/
    ++m;
/*<       go to 12                                                           >*/
    goto L12;
/*<    13 na=m-1                                                             >*/
L13:
    na = m - 1;
/*<       m=1                                                                >*/
    m = 1;
/*<       nim1=ni-1                                                          >*/
    nim1 = ni - 1;
/*<       if(na .ne. 1) go to 14                                             >*/
    if (na != 1) {
    goto L14;
    }
/*<       k2=lp(2,m)                                                         >*/
    k2 = lp[m * 3 + 2];
/*<       i2=lp(1,m)+k2-1                                                    >*/
    i2 = lp[m * 3 + 1] + k2 - 1;
/*<       efm=eft-1.0                                                        >*/
    efm = eft - (float)1.;
/*<       u=yv/(1.d0-1.d0/wn)**2                                             >*/
/* Computing 2nd power */
    d__1 = 1. - 1. / wn;
    u = yv / (d__1 * d__1);
/*<       s=sqrt(varf(nk,d,d(1,nkp2),sw,1,ni))                               >*/
    s = sqrt(varf_(nk, &d__[d_offset], &d__[nkp2 * d_dim1 + 1], &sw, &c__1, &
        ni));
/*    write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)                     1
568*/
/*<       return                                                             >*/
    return 0;
/*<    14 do 23 m=1,na                                                       >*/
L14:
    i__1 = na;
    for (m = 1; m <= i__1; ++m) {
/*<       k2=lp(2,m)                                                         >*/
    k2 = lp[m * 3 + 2];
/*<       l=lp(1,m)                                                          >*/
    l = lp[m * 3 + 1];
/*<       i2=l+k2-1                                                          >*/
    i2 = l + k2 - 1;
/*<       ll=k2-1                                                            >*/
    ll = k2 - 1;
/*<       np=ni                                                              >*/
    np = ni;
/*<       do 19 im=1,ni                                                      >*/
    i__2 = ni;
    for (im = 1; im <= i__2; ++im) {
/*<       i=t(im,1)+.1                                                       >*/
        i__ = t[im + t_dim1] + (float).1;
/*<       if(nord(i,tb) .eq. l) go to 15                                     >*/
        if (nord_(&i__, &tb[6]) == l) {
        goto L15;
        }
/*<       t(im,2)=0.0                                                        >*/
        t[im + (t_dim1 << 1)] = (float)0.;
/*<       go to 19                                                           >*/
        goto L19;
/*<    15 k=0                                                                >*/
L15:
        k = 0;
/*<       do 16 j=1,l                                                        >*/
        i__3 = l;
        for (j = 1; j <= i__3; ++j) {
/*<       if(jf(i,lv(ll+j),tb).eq.1) go to 16                                >*/
        if (jf_(&i__, &lv[ll + j], &tb[6]) == 1) {
            goto L16;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 17                                                           >*/
        goto L17;
/*<    16 continue                                                           >*/
L16:
        ;
        }
/*<    17 if(k .ne. 1) go to 18                                              >*/
L17:
        if (k != 1) {
        goto L18;
        }
/*<       t(im,2)=0.0                                                        >*/
        t[im + (t_dim1 << 1)] = (float)0.;
/*<       go to 19                                                           >*/
        goto L19;
/*<    18 t(im,2)=1.0                                                        >*/
L18:
        t[im + (t_dim1 << 1)] = (float)1.;
/*<       np=np-1                                                            >*/
        --np;
/*<    19 continue                                                           >*/
L19:
        ;
    }
/*<    20 k=0                                                                >*/
L20:
    k = 0;
/*<       do 21 i=1,nim1                                                     >*/
    i__2 = nim1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(t(i,2) .le. t(i+1,2)) go to 21                                  >*/
        if (t[i__ + (t_dim1 << 1)] <= t[i__ + 1 + (t_dim1 << 1)]) {
        goto L21;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       call exch(nk,ni,i,d,t,t(1,2))                                      >*/
        exch_(nk, &ni, &i__, &d__[d_offset], &t[t_offset], &t[(t_dim1 << 
            1) + 1]);
/*<    21 continue                                                           >*/
L21:
        ;
    }
/*<       if(k.eq.0) go to 22                                                >*/
    if (k == 0) {
        goto L22;
    }
/*<       go to 20                                                           >*/
    goto L20;
/*<    22 call lsf(nk,np,nkp1,0.d0,d,d(1,lm),s,u,d(1,nkp3),1)                >*/
L22:
    lsf_(nk, &np, &nkp1, &c_b395, &d__[d_offset], &d__[lm * d_dim1 + 1], &
        s, &u, &d__[nkp3 * d_dim1 + 1], &c__1);
/*<       efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)                                 >*/
    efm = efp_(&lp[m * 3 + 1], &lv[lp[m * 3 + 2]], nk, &tb[6]);
/*<       u=(u/sw+yv)/(1.d0-(eft-efm)/wn)**2                                 >*/
/* Computing 2nd power */
    d__1 = 1. - (eft - efm) / wn;
    u = (u / sw + yv) / (d__1 * d__1);
/*<       s=sqrt(varf(nk,d,d(1,nkp2),sw,np+1,ni))                            >*/
    i__2 = np + 1;
    s = sqrt(varf_(nk, &d__[d_offset], &d__[nkp2 * d_dim1 + 1], &sw, &
        i__2, &ni));
/*    write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)                  
   1605*/
/*<    23 continue                                                           >*/
/* L23: */
    }
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L24: */
/*<    25 format(' ',i3,' ',2g12.4,'  ',i2,'      ',f4.1,'  ',20i4)          >*/
/* L25: */
/*<    26 format(/' estimated optimal model = response mean.')               >*/
/* L26: */
/*<       end                                                                >*/
} /* anova_ */

/*<       subroutine anoval (n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sc,d)           >*/
/* Subroutine */ int anoval_(n, x, y, w, nk, il, it, az, tb, cm, lp, lv, sc, 
    d__)
integer *n;
real *x, *y, *w;
integer *nk, *il, *it;
real *az, *tb, *cm;
integer *lp, *lv;
real *sc;
doublereal *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1, 
        i__2, i__3;
    real r__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int cptb_(), coll_();
    extern integer nord_();
    extern /* Subroutine */ int setz_();
    static integer i__, j, k, l, m;
    static real u, a0;
    static integer i2, k2, na;
    extern integer jf_();
    static integer ni, ip;
    static doublereal yb;
    static integer ll;
    static doublereal wn;
    extern /* Subroutine */ int vp_();
    static doublereal sw, yv;
    static real efm;
    extern doublereal efp_();
    static real eft;

/*<       integer lp(3,*),lv(*)                                              >*/
/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*)                       >*/
/*<       double precision d(nk,*),sw,yv,wn,yb                               >*/
/*<       if(it.le.0) return                                                 >*/
    /* Parameter adjustments */
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    tb -= 6;
    --cm;
    lp -= 4;
    --lv;

    /* Function Body */
    if (*it <= 0) {
    return 0;
    }
/*<       sw=0.d0                                                            >*/
    sw = 0.;
/*<       yb=sw                                                              >*/
    yb = sw;
/*<       yv=yb                                                              >*/
    yv = yb;
/*<       wn=yv                                                              >*/
    wn = yv;
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<       yb=yb+w(i)*y(i)                                                    >*/
    yb += w[i__] * y[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       yb=yb/sw                                                           >*/
    yb /= sw;
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    d__1 = sw;
    wn = d__1 * d__1 / wn;
/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       yv=yv+w(i)*(y(i)-yb)**2                                            >*/
/* Computing 2nd power */
    d__1 = y[i__] - yb;
    yv += w[i__] * (d__1 * d__1);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       yv=yv/sw                                                           >*/
    yv /= sw;
/*<       eft=1.0                                                            >*/
    eft = (float)1.;
/*<       ni=0                                                               >*/
    ni = 0;
/*<       do 3 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 3                                         >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L3;
    }
/*<       ni=ni+1                                                            >*/
    ++ni;
/*<       eft=eft+tb(5,m)                                                    >*/
    eft += tb[m * 5 + 5];
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<       if(ni .ne. 0) go to 4                                              >*/
    if (ni != 0) {
    goto L4;
    }
/*    write(it,14)                                                       1
641*/
/*<       return                                                             >*/
    return 0;
/*<     4 continue >*/
L4:
/*    write(it,12) ni                                                    1
643*/
/*<       call coll(nk,tb,lp,lv,lp(1,nk+1))                                  >*/
    coll_(nk, &tb[6], &lp[4], &lv[1], &lp[(*nk + 1) * 3 + 1]);
/*<       m=1                                                                >*/
    m = 1;
/*<     5 if(lp(1,m).eq.0) go to 6                                           >*/
L5:
    if (lp[m * 3 + 1] == 0) {
    goto L6;
    }
/*<       m=m+1                                                              >*/
    ++m;
/*<       go to 5                                                            >*/
    goto L5;
/*<     6 na=m-1                                                             >*/
L6:
    na = m - 1;
/*<       m=1                                                                >*/
    m = 1;
/*<       if(na .ne. 1) go to 7                                              >*/
    if (na != 1) {
    goto L7;
    }
/*<       k2=lp(2,m)                                                         >*/
    k2 = lp[m * 3 + 2];
/*<       i2=lp(1,m)+k2-1                                                    >*/
    i2 = lp[m * 3 + 1] + k2 - 1;
/*<       efm=eft-1.0                                                        >*/
    efm = eft - (float)1.;
/*<       u=yv/(1.d0-1.d0/wn)**2                                             >*/
/* Computing 2nd power */
    d__1 = 1. - 1. / wn;
    u = yv / (d__1 * d__1);
/*    write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)                       1
656*/
/*<       return                                                             >*/
    return 0;
/*<     7 ip=nk+4                                                            >*/
L7:
    ip = *nk + 4;
/*<       do 11 m=1,na                                                       >*/
    i__1 = na;
    for (m = 1; m <= i__1; ++m) {
/*<       k2=lp(2,m)                                                         >*/
    k2 = lp[m * 3 + 2];
/*<       l=lp(1,m)                                                          >*/
    l = lp[m * 3 + 1];
/*<       i2=l+k2-1                                                          >*/
    i2 = l + k2 - 1;
/*<       ll=k2-1                                                            >*/
    ll = k2 - 1;
/*<       call cptb(nk,tb,sc(1,ip))                                          >*/
    cptb_(nk, &tb[6], &sc[ip * sc_dim1 + 1]);
/*<       do 10 i=1,nk                                                       >*/
    i__2 = *nk;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(tb(1,i).eq.0.0) go to 10                                        >*/
        if (tb[i__ * 5 + 1] == (float)0.) {
        goto L10;
        }
/*<       if(nord(i,tb).ne.l) go to 10                                       >*/
        if (nord_(&i__, &tb[6]) != l) {
        goto L10;
        }
/*<       k=0                                                                >*/
        k = 0;
/*<       do 8 j=1,l                                                         >*/
        i__3 = l;
        for (j = 1; j <= i__3; ++j) {
/*<       if(jf(i,lv(ll+j),tb).eq.1) go to 8                                 >*/
        if (jf_(&i__, &lv[ll + j], &tb[6]) == 1) {
            goto L8;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 9                                                            >*/
        goto L9;
/*<     8 continue                                                           >*/
L8:
        ;
        }
/*<     9 if(k.eq.1) go to 10                                                >*/
L9:
        if (k == 1) {
        goto L10;
        }
/*<       call setz(i,sc(1,ip))                                              >*/
        setz_(&i__, &sc[ip * sc_dim1 + 1]);
/*<    10 continue                                                           >*/
L10:
        ;
    }
/*<       a0=az                                                              >*/
    a0 = *az;
/*<       call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,u,sc,d)                 >*/
    vp_(n, &x[x_offset], &y[1], &w[1], nk, il, &yb, &sw, &a0, &sc[ip * 
        sc_dim1 + 1], &cm[1], &u, &sc[sc_offset], &d__[d_offset]);
/*<       efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)                                 >*/
    efm = efp_(&lp[m * 3 + 1], &lv[lp[m * 3 + 2]], nk, &tb[6]);
/*<       u=u/(1.d0-(eft-efm)/wn)**2                                         >*/
/* Computing 2nd power */
    d__1 = 1. - (eft - efm) / wn;
    u /= d__1 * d__1;
/*    write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)                    
   1681*/
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L12: */
/*<    13 format(' ',i3,' ',g12.4,'   ',i2,'     ',f4.1,'    ',20i4)         >*/
/* L13: */
/*<    14 format(/' estimated optimal model = response mean.')               >*/
/* L14: */
/*<       end                                                                >*/
} /* anoval_ */

/*<       subroutine cptb(nk,tb,ub)                                          >*/
/* Subroutine */ int cptb_0_(n__, nk, tb, ub, l)
int n__;
integer *nk;
real *tb, *ub;
integer *l;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, m;

/*<       real tb(5,nk),ub(5,*)                                              >*/
/*<       do 2 m=1,nk                                                        >*/
    /* Parameter adjustments */
    if (tb) {
    tb -= 6;
    }
    ub -= 6;

    /* Function Body */
    switch(n__) {
    case 1: goto L_setz;
    }

    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       do 1 k=1,5                                                         >*/
    for (k = 1; k <= 5; ++k) {
/*<       ub(k,m)=tb(k,m)                                                    >*/
        ub[k + m * 5] = tb[k + m * 5];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       return                                                             >*/
    return 0;
/*<       entry setz(l,ub)                                                   >*/

L_setz:
/*<       ub(1,l)=0.0                                                        >*/
    ub[*l * 5 + 1] = (float)0.;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* cptb_ */

/* Subroutine */ int cptb_(nk, tb, ub)
integer *nk;
real *tb, *ub;
{
    return cptb_0_(0, nk, tb, ub, (integer *)0);
    }

/* Subroutine */ int setz_(l, ub)
integer *l;
real *ub;
{
    return cptb_0_(1, (integer *)0, (real *)0, ub, l);
    }

/*<       subroutine fun (l,jv,n,x,nk,tb,cm,jl,kv,t,js)                      >*/
/* Subroutine */ int fun_(l, jv, n, x, nk, tb, cm, jl, kv, t, js)
integer *l, *jv, *n;
real *x;
integer *nk;
real *tb, *cm;
integer *jl, *kv;
real *t;
integer *js;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal s;
    static real u;
    extern integer nordc_(), jf_();
    static integer ip;
    extern integer icf_();
    static real phi;

/*<       integer jv(l),kv(2,jl),js(*)                                       >*/
/*<       real x(n,l),tb(5,nk),cm(*),t(n)                                    >*/
/*<       double precision s                                                 >*/
/*<       do 8 i=1,n                                                         >*/
    /* Parameter adjustments */
    --jv;
    --t;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    tb -= 6;
    --cm;
    kv -= 3;
    --js;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 7 m=1,nk                                                        >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7                             >*/
        if (icf_(&m, &tb[6], &cm[1], jl, &kv[3], &js[1]) == 0) {
        goto L7;
        }
/*<       if(nordc(1,m,tb,cm).ne.l) go to 7                                  >*/
        if (nordc_(&c__1, &m, &tb[6], &cm[1]) != *l) {
        goto L7;
        }
/*<       k=0                                                                >*/
        k = 0;
/*<       do 1 j=1,l                                                         >*/
        i__3 = *l;
        for (j = 1; j <= i__3; ++j) {
/*<       if(jf(m,jv(j),tb).eq.1) go to 1                                    >*/
        if (jf_(&m, &jv[j], &tb[6]) == 1) {
            goto L1;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 continue                                                           >*/
L1:
        ;
        }
/*<     2 if(k.eq.1) go to 7                                                 >*/
L2:
        if (k == 1) {
        goto L7;
        }
/*<       phi=1.0                                                            >*/
        phi = (float)1.;
/*<       ip=m                                                               >*/
        ip = m;
/*<     3 if(ip.le.0) go to 6                                                >*/
L3:
        if (ip <= 0) {
        goto L6;
        }
/*<       u=tb(2,ip)                                                         >*/
        u = tb[ip * 5 + 2];
/*<       j=abs(u)+.1                                                        >*/
        j = dabs(u) + (float).1;
/*<       if(cm(2*j) .eq. 0.0) go to 4                                       >*/
        if (cm[j * 2] == (float)0.) {
        goto L4;
        }
/*<       ip=tb(4,ip)+.1                                                     >*/
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
        goto L3;
/*<     4 do 5 k=1,l                                                         >*/
L4:
        i__3 = *l;
        for (k = 1; k <= i__3; ++k) {
/*<       if(j.eq.jv(k)) j=k                                                 >*/
        if (j == jv[k]) {
            j = k;
        }
/*<     5 continue                                                           >*/
/* L5: */
        }
/*<       phi=phi*amax1(0.0,sign(1.0,u)*(x(i,j)-tb(3,ip)))                   >*/
/* Computing MAX */
        r__1 = (float)0., r__2 = r_sign(&c_b196, &u) * (x[i__ + j * 
            x_dim1] - tb[ip * 5 + 3]);
        phi *= dmax(r__1,r__2);
/*<       ip=tb(4,ip)+.1                                                     >*/
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
        goto L3;
/*<     6 s=s+tb(1,m)*phi                                                    >*/
L6:
        s += tb[m * 5 + 1] * phi;
/*<     7 continue                                                           >*/
L7:
        ;
    }
/*<       t(i)=s                                                             >*/
    t[i__] = s;
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* fun_ */

/*<    >*/
/* Subroutine */ int cubic_(n, p, x, y, w, nk, it, tb, cm, kp, kv, lp, lv, bz,
     tc, t, z__, sc, js, d__)
integer *n, *p;
real *x, *y, *w;
integer *nk, *it;
real *tb, *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *t, *z__, *sc;
integer *js;
doublereal *d__;
{
    /* Initialized data */

    static real big = (float)9.9e30;

    /* System generated locals */
    integer x_dim1, x_offset, t_dim1, t_offset, d_dim1, d_offset, i__1, i__2, 
        i__3;
    real r__1, r__2;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int side_(), knts_();
    static integer i__, j, k, l, m;
    static doublereal s, u;
    static integer l1, ic, la, le, jj, il, ni;
    static doublereal yb;
    static integer lm;
    static real xl;
    static doublereal wn;
    static integer ll, lt, jl, kk;
    static doublereal sw;
    static real xr;
    static integer nt, jp;
    static doublereal yv;
    static integer kp3;
    static real eft;
    extern /* Subroutine */ int lsf_(), que_();
    static integer nkp1, nkp2, nkp3;

/*<       integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*)                      >*/
/*<       real x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),t(n,nk),z(2,p),sc(n)    >*/
/*<       double precision d(nk,*),s,u,sw,yb,wn,yv                           >*/
/*<       data big /9.9e30/                                                  >*/
    /* Parameter adjustments */
    --sc;
    --w;
    --y;
    z__ -= 3;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    t_dim1 = *n;
    t_offset = t_dim1 + 1;
    t -= t_offset;
    tb -= 6;
    --cm;
    kp -= 6;
    kv -= 3;
    lp -= 4;
    --lv;
    --tc;
    --js;

    /* Function Body */
/*<       yb=0.d0                                                            >*/
    yb = 0.;
/*<       sw=yb                                                              >*/
    sw = yb;
/*<       wn=sw                                                              >*/
    wn = sw;
/*<       yv=wn                                                              >*/
    yv = wn;
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<       yb=yb+w(i)*y(i)                                                    >*/
    yb += w[i__] * y[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       yb=yb/sw                                                           >*/
    yb /= sw;
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    d__1 = sw;
    wn = d__1 * d__1 / wn;
/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       yv=yv+w(i)*(y(i)-yb)**2                                            >*/
/* Computing 2nd power */
    d__1 = y[i__] - yb;
    yv += w[i__] * (d__1 * d__1);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       yv=yv/sw                                                           >*/
    yv /= sw;
/*<       ni=0                                                               >*/
    ni = 0;
/*<       do 3 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).ne.0.0) ni=ni+1                                         >*/
    if (tb[m * 5 + 1] != (float)0.) {
        ++ni;
    }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       if(ni .ne. 0) go to 4                                              >*/
    if (ni != 0) {
    goto L4;
    }
/*<       bz=yb                                                              >*/
    *bz = yb;
/*<       u=yv/(1.0-1.0/wn)**2                                               >*/
/* Computing 2nd power */
    d__1 = (float)1. - (float)1. / wn;
    u = yv / (d__1 * d__1);
/*    if(it.gt.0) write(it,34) ni,u                                      1
765*/
/*<       return                                                             >*/
    return 0;
/*<     4 nkp1=nk+1                                                          >*/
L4:
    nkp1 = *nk + 1;
/*<       nkp2=nk+2                                                          >*/
    nkp2 = *nk + 2;
/*<       nkp3=nk+3                                                          >*/
    nkp3 = *nk + 3;
/*<       lm=nkp3+nk                                                         >*/
    lm = nkp3 + *nk;
/*<       do 6 i=1,p                                                         >*/
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       xl=big                                                             >*/
    xl = big;
/*<       xr=-xl                                                             >*/
    xr = -xl;
/*<       do 5 j=1,n                                                         >*/
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
/*<       xl=amin1(xl,x(j,i))                                                >*/
/* Computing MIN */
        r__1 = xl, r__2 = x[j + i__ * x_dim1];
        xl = dmin(r__1,r__2);
/*<       xr=amax1(xr,x(j,i))                                                >*/
/* Computing MAX */
        r__1 = xr, r__2 = x[j + i__ * x_dim1];
        xr = dmax(r__1,r__2);
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       z(1,i)=xl                                                          >*/
    z__[(i__ << 1) + 1] = xl;
/*<       z(2,i)=xr                                                          >*/
    z__[(i__ << 1) + 2] = xr;
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       ll=1                                                               >*/
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<       lt=0                                                               >*/
    lt = 0;
/*<     7 if(kp(1,ll).lt.0) go to 20                                         >*/
L7:
    if (kp[ll * 5 + 1] < 0) {
    goto L20;
    }
/*<       do 8 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sc(i)=1.0                                                          >*/
    sc[i__] = (float)1.;
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       if(kp(1,ll) .le. 0) go to 12                                       >*/
    if (kp[ll * 5 + 1] <= 0) {
    goto L12;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       do 11 il=1,jl                                                      >*/
    i__1 = jl;
    for (il = 1; il <= i__1; ++il) {
/*<       k=kp(2,ll)+il-1                                                    >*/
    k = kp[ll * 5 + 2] + il - 1;
/*<       jj=kv(1,k)                                                         >*/
    jj = kv[(k << 1) + 1];
/*<       j=iabs(jj)                                                         >*/
    j = abs(jj);
/*<       kk=kv(2,k)                                                         >*/
    kk = kv[(k << 1) + 2];
/*<       do 10 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(sc(i).eq.0.0) go to 10                                          >*/
        if (sc[i__] == (float)0.) {
        goto L10;
        }
/*<       ic=x(i,j)+.1                                                       >*/
        ic = x[i__ + j * x_dim1] + (float).1;
/*<       sc(i)=cm(ic+kk)                                                    >*/
        sc[i__] = cm[ic + kk];
/*<       if(jj .ge. 0) go to 10                                             >*/
        if (jj >= 0) {
        goto L10;
        }
/*<       if(sc(i) .ne. 0.0) go to 9                                         >*/
        if (sc[i__] != (float)0.) {
        goto L9;
        }
/*<       sc(i)=1.0                                                          >*/
        sc[i__] = (float)1.;
/*<       go to 10                                                           >*/
        goto L10;
/*<     9 sc(i)=0.0                                                          >*/
L9:
        sc[i__] = (float)0.;
/*<    10 continue                                                           >*/
L10:
        ;
    }
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       go to 13                                                           >*/
    goto L13;
/*<    12 if(kp(3,ll) .gt. 0) go to 13                                       >*/
L12:
    if (kp[ll * 5 + 3] > 0) {
    goto L13;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 7                                                            >*/
    goto L7;
/*<    13 if(kp(3,ll) .gt. 0) go to 15                                       >*/
L13:
    if (kp[ll * 5 + 3] > 0) {
    goto L15;
    }
/*<       lt=lt+1                                                            >*/
    ++lt;
/*<       kp(5,ll)=0                                                         >*/
    kp[ll * 5 + 5] = 0;
/*<       do 14 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       t(i,lt)=sc(i)                                                      >*/
    t[i__ + lt * t_dim1] = sc[i__];
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       go to 19                                                           >*/
    goto L19;
/*<    15 kp3=kp(3,ll)                                                       >*/
L15:
    kp3 = kp[ll * 5 + 3];
/*<       kp(5,ll)=la                                                        >*/
    kp[ll * 5 + 5] = la;
/*<       do 18 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<    >*/
    knts_(&l, &nt, &lv[lp[l1 * 3 + 2]], &kp[ll * 5 + 1], &kv[(kp[ll * 5 + 
        2] << 1) + 1], nk, &tb[6], &cm[1], &tc[la], &js[1]);
/*<       call side(l,nt,lv(lp(2,l1)),z,tc(la))                              >*/
    side_(&l, &nt, &lv[lp[l1 * 3 + 2]], &z__[3], &tc[la]);
/*<       do 17 jp=1,nt                                                      >*/
    i__2 = nt;
    for (jp = 1; jp <= i__2; ++jp) {
/*<       lt=lt+1                                                            >*/
        ++lt;
/*<       do 16 i=1,n                                                        >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       t(i,lt)=sc(i)                                                      >*/
        t[i__ + lt * t_dim1] = sc[i__];
/*<    16 continue                                                           >*/
/* L16: */
        }
/*<       call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),t(1,lt))                  >*/
        que_(&jp, &l, &nt, &lv[lp[l1 * 3 + 2]], n, &x[x_offset], &tc[la], 
            &t[lt * t_dim1 + 1]);
/*<    17 continue                                                           >*/
/* L17: */
    }
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       la=la+nt*(5*l+1)                                                   >*/
    la += nt * (l * 5 + 1);
/*<    18 continue                                                           >*/
/* L18: */
    }
/*<    19 ll=ll+1                                                            >*/
L19:
    ++ll;
/*<       go to 7                                                            >*/
    goto L7;
/*<    20 do 26 j=1,lt                                                       >*/
L20:
    i__1 = lt;
    for (j = 1; j <= i__1; ++j) {
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       u=s                                                                >*/
    u = s;
/*<       do 21 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+w(i)*t(i,j)                                                    >*/
        s += w[i__] * t[i__ + j * t_dim1];
/*<    21 continue                                                           >*/
/* L21: */
    }
/*<       s=s/sw                                                             >*/
    s /= sw;
/*<       d(j,nkp2)=s                                                        >*/
    d__[j + nkp2 * d_dim1] = s;
/*<       do 22 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       t(i,j)=t(i,j)-s                                                    >*/
        t[i__ + j * t_dim1] -= s;
/*<    22 continue                                                           >*/
/* L22: */
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 23 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+w(i)*(y(i)-yb)*t(i,j)                                          >*/
        s += w[i__] * (y[i__] - yb) * t[i__ + j * t_dim1];
/*<    23 continue                                                           >*/
/* L23: */
    }
/*<       d(j,nkp1)=s                                                        >*/
    d__[j + nkp1 * d_dim1] = s;
/*<       do 25 k=1,j                                                        >*/
    i__2 = j;
    for (k = 1; k <= i__2; ++k) {
/*<       s=0.d0                                                             >*/
        s = 0.;
/*<       do 24 i=1,n                                                        >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       s=s+w(i)*t(i,k)*t(i,j)                                             >*/
        s += w[i__] * t[i__ + k * t_dim1] * t[i__ + j * t_dim1];
/*<    24 continue                                                           >*/
/* L24: */
        }
/*<       d(k,j)=s                                                           >*/
        d__[k + j * d_dim1] = s;
/*<    25 continue                                                           >*/
/* L25: */
    }
/*<    26 continue                                                           >*/
/* L26: */
    }
/*<       call lsf(nk,lt,nkp1,yb,d,d(1,lm),s,u,d(1,nkp3),1)                  >*/
    lsf_(nk, &lt, &nkp1, &yb, &d__[d_offset], &d__[lm * d_dim1 + 1], &s, &u, &
        d__[nkp3 * d_dim1 + 1], &c__1);
/*<       eft=1.0                                                            >*/
    eft = (float)1.;
/*<       do 27 i=1,nk                                                       >*/
    i__1 = *nk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(tb(1,i).ne.0.0) eft=eft+tb(5,i)                                 >*/
    if (tb[i__ * 5 + 1] != (float)0.) {
        eft += tb[i__ * 5 + 5];
    }
/*<    27 continue                                                           >*/
/* L27: */
    }
/*<       u=(u/sw+yv)/(1.0-eft/wn)**2                                        >*/
/* Computing 2nd power */
    d__1 = (float)1. - eft / wn;
    u = (u / sw + yv) / (d__1 * d__1);
/*<       bz=s                                                               >*/
    *bz = s;
/*<       ll=1                                                               >*/
    ll = 1;
/*<       l1=ll                                                              >*/
    l1 = ll;
/*<       le=la-1                                                            >*/
    le = la - 1;
/*<       la=0                                                               >*/
    la = 0;
/*<       lt=la                                                              >*/
    lt = la;
/*<    28 if(kp(1,ll).lt.0) go to 33                                         >*/
L28:
    if (kp[ll * 5 + 1] < 0) {
    goto L33;
    }
/*<       if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 29                  >*/
    if (kp[ll * 5 + 1] != 0 || kp[ll * 5 + 3] > 0) {
    goto L29;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 28                                                           >*/
    goto L28;
/*<    29 if(kp(3,ll) .gt. 0) go to 30                                       >*/
L29:
    if (kp[ll * 5 + 3] > 0) {
    goto L30;
    }
/*<       le=le+1                                                            >*/
    ++le;
/*<       kp(3,ll)=-le                                                       >*/
    kp[ll * 5 + 3] = -le;
/*<       lt=lt+1                                                            >*/
    ++lt;
/*<       tc(le)=d(lt,lm)                                                    >*/
    tc[le] = d__[lt + lm * d_dim1];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 28                                                           >*/
    goto L28;
/*<    30 kp3=kp(3,ll)                                                       >*/
L30:
    kp3 = kp[ll * 5 + 3];
/*<       do 32 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       la=la+5*lp(1,l1)*nt                                                >*/
    la += lp[l1 * 3 + 1] * 5 * nt;
/*<       do 31 i=1,nt                                                       >*/
    i__2 = nt;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       lt=lt+1                                                            >*/
        ++lt;
/*<       tc(i+la)=d(lt,lm)                                                  >*/
        tc[i__ + la] = d__[lt + lm * d_dim1];
/*<    31 continue                                                           >*/
/* L31: */
    }
/*<       la=la+nt                                                           >*/
    la += nt;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<    32 continue                                                           >*/
/* L32: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 28                                                           >*/
    goto L28;
/*<    33 continue >*/
L33:
/*    if(it.gt.0) write(it,34) lt,u                                      1
898*/
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L34: */
/*<       end                                                                >*/
} /* cubic_ */

/*<       subroutine cfun (l,jv,n,x,nf,lp,lv,tc,t,sc,jw)                     >*/
/* Subroutine */ int cfun_(l, jv, n, x, nf, lp, lv, tc, t, sc, jw)
integer *l, *jv, *n;
real *x;
integer *nf, *lp, *lv;
real *tc, *t, *sc;
integer *jw;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m, l1, l2, la, lb, nt;
    extern /* Subroutine */ int que_();

/*<       integer jv(l),lp(3,*),lv(*),jw(l)                                  >*/
/*<       real x(n,l),tc(*),t(n),sc(n)                                       >*/
/*<       do 1 i=1,n                                                         >*/
    /* Parameter adjustments */
    --jw;
    --jv;
    --sc;
    --t;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    lp -= 4;
    --lv;
    --tc;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       t(i)=0.0                                                           >*/
    t[i__] = (float)0.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       la=1                                                               >*/
    la = 1;
/*<       do 10 l1=1,nf                                                      >*/
    i__1 = *nf;
    for (l1 = 1; l1 <= i__1; ++l1) {
/*<       if(lp(1,l1).ne.l) go to 9                                          >*/
    if (lp[l1 * 3 + 1] != *l) {
        goto L9;
    }
/*<       l2=lp(2,l1)-1                                                      >*/
    l2 = lp[l1 * 3 + 2] - 1;
/*<       do 3 j=1,l                                                         >*/
    i__2 = *l;
    for (j = 1; j <= i__2; ++j) {
/*<       m=0                                                                >*/
        m = 0;
/*<       do 2 k=1,l                                                         >*/
        i__3 = *l;
        for (k = 1; k <= i__3; ++k) {
/*<       if(jv(j).eq.lv(k+l2)) m=1                                          >*/
        if (jv[j] == lv[k + l2]) {
            m = 1;
        }
/*<     2 continue                                                           >*/
/* L2: */
        }
/*<       if(m.eq.0) go to 9                                                 >*/
        if (m == 0) {
        goto L9;
        }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       lb=la+5*l*nt-1                                                     >*/
    lb = la + *l * 5 * nt - 1;
/*<       do 8 j=1,nt                                                        >*/
    i__2 = nt;
    for (j = 1; j <= i__2; ++j) {
/*<       do 5 k=1,l                                                         >*/
        i__3 = *l;
        for (k = 1; k <= i__3; ++k) {
/*<       do 4 i=1,l                                                         >*/
        i__4 = *l;
        for (i__ = 1; i__ <= i__4; ++i__) {
/*<       if(lv(k+l2).eq.jv(i)) jw(k)=i                                      >*/
            if (lv[k + l2] == jv[i__]) {
            jw[k] = i__;
            }
/*<     4 continue                                                           >*/
/* L4: */
        }
/*<     5 continue                                                           >*/
/* L5: */
        }
/*<       do 6 i=1,n                                                         >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       sc(i)=1.0                                                          >*/
        sc[i__] = (float)1.;
/*<     6 continue                                                           >*/
/* L6: */
        }
/*<       call que(j,l,nt,jw,n,x,tc(la),sc)                                  >*/
        que_(&j, l, &nt, &jw[1], n, &x[x_offset], &tc[la], &sc[1]);
/*<       do 7 i=1,n                                                         >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       t(i)=t(i)+tc(lb+j)*sc(i)                                           >*/
        t[i__] += tc[lb + j] * sc[i__];
/*<     7 continue                                                           >*/
/* L7: */
        }
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       go to 11                                                           >*/
    goto L11;
/*<     9 la=la+lp(3,l1)*(5*lp(1,l1)+1)                                      >*/
L9:
    la += lp[l1 * 3 + 3] * (lp[l1 * 3 + 1] * 5 + 1);
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<    11 return                                                             >*/
L11:
    return 0;
/*<       end                                                                >*/
} /* cfun_ */

/*<       subroutine orgpc (xm,xs,lp,lv,tc)                                  >*/
/* Subroutine */ int orgpc_(xm, xs, lp, lv, tc)
real *xm, *xs;
integer *lp, *lv;
real *tc;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int scpc_();
    static integer j, l, l1, la, lb, nt;

/*<       integer lp(3,*),lv(*)                                              >*/
/*<       real xm(*),xs(*),tc(*)                                             >*/
/*<       la=1                                                               >*/
    /* Parameter adjustments */
    --tc;
    --lv;
    lp -= 4;
    --xs;
    --xm;

    /* Function Body */
    la = 1;
/*<       l1=la                                                              >*/
    l1 = la;
/*<     1 if(lp(1,l1).eq.0) go to 3                                          >*/
L1:
    if (lp[l1 * 3 + 1] == 0) {
    goto L3;
    }
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       lb=la+5*l*nt-1                                                     >*/
    lb = la + l * 5 * nt - 1;
/*<       do 2 j=1,nt                                                        >*/
    i__1 = nt;
    for (j = 1; j <= i__1; ++j) {
/*<       call scpc(xm,xs,j,l,nt,lv(lp(2,l1)),tc(la),tc(lb+j))               >*/
    scpc_(&xm[1], &xs[1], &j, &l, &nt, &lv[lp[l1 * 3 + 2]], &tc[la], &tc[
        lb + j]);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       la=lb+nt+1                                                         >*/
    la = lb + nt + 1;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     3 return                                                             >*/
L3:
    return 0;
/*<       end                                                                >*/
} /* orgpc_ */

/*<       subroutine pair (jv,n,x,nk,tb,cm,jl,kv,f,sc,js)                    >*/
/* Subroutine */ int pair_(jv, n, x, nk, tb, cm, jl, kv, f, sc, js)
integer *jv, *n;
real *x;
integer *nk;
real *tb, *cm;
integer *jl, *kv;
real *f, *sc;
integer *js;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int fun_();

/*<       integer jv(2),kv(2,jl),js(*)                                       >*/
/*<       real x(n,*),tb(5,nk),cm(*),f(n),sc(n)                              >*/
/*<       call fun(2,jv,n,x,nk,tb,cm,jl,kv,f,js)                             >*/
    /* Parameter adjustments */
    --jv;
    --sc;
    --f;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    tb -= 6;
    --cm;
    kv -= 3;
    --js;

    /* Function Body */
    fun_(&c__2, &jv[1], n, &x[x_offset], nk, &tb[6], &cm[1], jl, &kv[3], &f[1]
        , &js[1]);
/*<       do 2 k=1,2                                                         >*/
    for (k = 1; k <= 2; ++k) {
/*<       call fun(1,jv(k),n,x(1,k),nk,tb,cm,jl,kv,sc,js)                    >*/
    fun_(&c__1, &jv[k], n, &x[k * x_dim1 + 1], nk, &tb[6], &cm[1], jl, &
        kv[3], &sc[1], &js[1]);
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       f(i)=f(i)+sc(i)                                                    >*/
        f[i__] += sc[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* pair_ */

/*<       subroutine cpair (jv,n,x,nf,lp,lv,tc,f,sc)                         >*/
/* Subroutine */ int cpair_(jv, n, x, nf, lp, lv, tc, f, sc)
integer *jv, *n;
real *x;
integer *nf, *lp, *lv;
real *tc, *f, *sc;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int cfun_();
    static integer i__, k, jw[2];

/*<       integer jv(2),lp(3,*),lv(*),jw(2)                                  >*/
/*<       real x(n,*),tc(*),f(n),sc(n,2)                                     >*/
/*<       call cfun(2,jv,n,x,nf,lp,lv,tc,f,sc,jw)                            >*/
    /* Parameter adjustments */
    --jv;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --f;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    lp -= 4;
    --lv;
    --tc;

    /* Function Body */
    cfun_(&c__2, &jv[1], n, &x[x_offset], nf, &lp[4], &lv[1], &tc[1], &f[1], &
        sc[sc_offset], jw);
/*<       do 2 k=1,2                                                         >*/
    for (k = 1; k <= 2; ++k) {
/*<       call cfun(1,jv(k),n,x(1,k),nf,lp,lv,tc,sc,sc(1,2),jw)              >*/
    cfun_(&c__1, &jv[k], n, &x[k * x_dim1 + 1], nf, &lp[4], &lv[1], &tc[1]
        , &sc[sc_offset], &sc[(sc_dim1 << 1) + 1], jw);
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       f(i)=f(i)+sc(i,1)                                                  >*/
        f[i__] += sc[i__ + sc_dim1];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* cpair_ */

/*<       subroutine logitl (n,x,y,w,nk,il,az,tb,cm,sc,d)                    >*/
/* Subroutine */ int logitl_0_(n__, n, x, y, w, nk, il, az, tb, cm, sc, d__, 
    kp, kv, lp, lv, bz, tc, ss)
int n__;
integer *n;
real *x, *y, *w;
integer *nk, *il;
real *az, *tb, *cm, *sc;
doublereal *d__;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *ss;
{
    /* Initialized data */

    static integer niter = 30;
    static real wm = (float)1e-4;
    static real thr = (float)1e-4;

    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1, 
        i__2, i__3;
    real r__1;

    /* Builtin functions */
    double log(), exp();

    /* Local variables */
    static integer iter;
    static doublereal a, b;
    static integer i__, j, k, l, m;
    static doublereal s;
    static integer l1, ic, la, jj, mk;
    static doublereal yb;
    static integer ll, lt, jl, kk, nt, jp;
    static doublereal sw;
    static real pp, ww;
    static integer mm1, kp3;
    extern doublereal phi_();
    extern /* Subroutine */ int lsf_(), que_();
    static integer jnt, mkp1, mkp2, mkp3, mkp4;

/*<       integer kp(5,*),kv(2,*),lp(3,*),lv(*)                              >*/
/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*),tc(*),ss(n)           >*/
/*<       double precision d(nk,*),a,b,s,sw,yb                               >*/
/*<       data niter,wm,thr /30,0.0001,0.0001/                               >*/
    /* Parameter adjustments */
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --y;
    --w;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    if (ss) {
    --ss;
    }
    if (tb) {
    tb -= 6;
    }
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    --cm;
    if (kp) {
    kp -= 6;
    }
    if (kv) {
    kv -= 3;
    }
    if (lp) {
    lp -= 4;
    }
    if (lv) {
    --lv;
    }
    if (tc) {
    --tc;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_logitc;
    }

/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       k=0                                                                >*/
    k = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 1                                         >*/
        if (tb[m * 5 + 1] == (float)0.) {
        goto L1;
        }
/*<       k=k+1                                                              >*/
        ++k;
/*<       sc(i,k)=phi(m,i,n,x,tb,cm)                                         >*/
        sc[i__ + k * sc_dim1] = phi_(&m, &i__, n, &x[x_offset], &tb[6], &
            cm[1]);
/*<     1 continue                                                           >*/
L1:
        ;
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       if(k .ne. 0) go to 3                                               >*/
    if (k != 0) {
    goto L3;
    }
/*<       az=alog(az/(1.0-az))                                               >*/
    *az = log(*az / ((float)1. - *az));
/*<       return                                                             >*/
    return 0;
/*<     3 mk=k                                                               >*/
L3:
    mk = k;
/*<       a=az                                                               >*/
    a = *az;
/*<       jnt=1                                                              >*/
    jnt = 1;
/*<       go to 19                                                           >*/
    goto L19;
/*<       entry logitc (n,x,y,w,nk,il,cm,kp,kv,lp,lv,bz,tc,sc,ss,d)          >*/

L_logitc:
/*<       ll=1                                                               >*/
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<       lt=0                                                               >*/
    lt = 0;
/*<     4 if(kp(1,ll).lt.0) go to 17                                         >*/
L4:
    if (kp[ll * 5 + 1] < 0) {
    goto L17;
    }
/*<       do 5 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       ss(i)=1.0                                                          >*/
    ss[i__] = (float)1.;
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       if(kp(1,ll) .le. 0) go to 9                                        >*/
    if (kp[ll * 5 + 1] <= 0) {
    goto L9;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       do 8 il=1,jl                                                       >*/
    i__1 = jl;
    for (*il = 1; *il <= i__1; ++(*il)) {
/*<       k=kp(2,ll)+il-1                                                    >*/
    k = kp[ll * 5 + 2] + *il - 1;
/*<       jj=kv(1,k)                                                         >*/
    jj = kv[(k << 1) + 1];
/*<       j=iabs(jj)                                                         >*/
    j = abs(jj);
/*<       kk=kv(2,k)                                                         >*/
    kk = kv[(k << 1) + 2];
/*<       do 7 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(ss(i).eq.0.0) go to 7                                           >*/
        if (ss[i__] == (float)0.) {
        goto L7;
        }
/*<       ic=x(i,j)+.1                                                       >*/
        ic = x[i__ + j * x_dim1] + (float).1;
/*<       ss(i)=cm(ic+kk)                                                    >*/
        ss[i__] = cm[ic + kk];
/*<       if(jj .ge. 0) go to 7                                              >*/
        if (jj >= 0) {
        goto L7;
        }
/*<       if(ss(i) .ne. 0.0) go to 6                                         >*/
        if (ss[i__] != (float)0.) {
        goto L6;
        }
/*<       ss(i)=1.0                                                          >*/
        ss[i__] = (float)1.;
/*<       go to 7                                                            >*/
        goto L7;
/*<     6 ss(i)=0.0                                                          >*/
L6:
        ss[i__] = (float)0.;
/*<     7 continue                                                           >*/
L7:
        ;
    }
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 if(kp(3,ll) .gt. 0) go to 10                                       >*/
L9:
    if (kp[ll * 5 + 3] > 0) {
    goto L10;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 4                                                            >*/
    goto L4;
/*<    10 if(kp(3,ll) .gt. 0) go to 12                                       >*/
L10:
    if (kp[ll * 5 + 3] > 0) {
    goto L12;
    }
/*<       lt=lt+1                                                            >*/
    ++lt;
/*<       do 11 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sc(i,lt)=ss(i)                                                     >*/
    sc[i__ + lt * sc_dim1] = ss[i__];
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       go to 16                                                           >*/
    goto L16;
/*<    12 kp3=kp(3,ll)                                                       >*/
L12:
    kp3 = kp[ll * 5 + 3];
/*<       do 15 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       do 14 jp=1,nt                                                      >*/
    i__2 = nt;
    for (jp = 1; jp <= i__2; ++jp) {
/*<       lt=lt+1                                                            >*/
        ++lt;
/*<       do 13 i=1,n                                                        >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       sc(i,lt)=ss(i)                                                     >*/
        sc[i__ + lt * sc_dim1] = ss[i__];
/*<    13 continue                                                           >*/
/* L13: */
        }
/*<       call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,lt))                 >*/
        que_(&jp, &l, &nt, &lv[lp[l1 * 3 + 2]], n, &x[x_offset], &tc[la], 
            &sc[lt * sc_dim1 + 1]);
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       la=la+nt*(5*l+1)                                                   >*/
    la += nt * (l * 5 + 1);
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<    16 ll=ll+1                                                            >*/
L16:
    ++ll;
/*<       go to 4                                                            >*/
    goto L4;
/*<    17 if(lt .ne. 0) go to 18                                             >*/
L17:
    if (lt != 0) {
    goto L18;
    }
/*<       bz=alog(bz/(1.0-bz))                                               >*/
    *bz = log(*bz / ((float)1. - *bz));
/*<       return                                                             >*/
    return 0;
/*<    18 mk=lt                                                              >*/
L18:
    mk = lt;
/*<       a=bz                                                               >*/
    a = *bz;
/*<       jnt=2                                                              >*/
    jnt = 2;
/*<    19 mkp1=mk+1                                                          >*/
L19:
    mkp1 = mk + 1;
/*<       mkp2=mk+2                                                          >*/
    mkp2 = mk + 2;
/*<       mkp3=mk+3                                                          >*/
    mkp3 = mk + 3;
/*<       mkp4=mk+4                                                          >*/
    mkp4 = mk + 4;
/*<       iter=0                                                             >*/
    iter = 0;
/*<       if(jnt .ne. 1) go to 21                                            >*/
    if (jnt != 1) {
    goto L21;
    }
/*<       k=0                                                                >*/
    k = 0;
/*<       do 20 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 20                                        >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L20;
    }
/*<       k=k+1                                                              >*/
    ++k;
/*<       d(k,mkp3)=tb(1,m)                                                  >*/
    d__[k + mkp3 * d_dim1] = tb[m * 5 + 1];
/*<    20 continue                                                           >*/
L20:
    ;
    }
/*<       go to 27                                                           >*/
    goto L27;
/*<    21 ll=1                                                               >*/
L21:
    ll = 1;
/*<       l1=ll                                                              >*/
    l1 = ll;
/*<       la=0                                                               >*/
    la = 0;
/*<       lt=la                                                              >*/
    lt = la;
/*<    22 if(kp(1,ll).lt.0) go to 27                                         >*/
L22:
    if (kp[ll * 5 + 1] < 0) {
    goto L27;
    }
/*<       if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 23                  >*/
    if (kp[ll * 5 + 1] != 0 || kp[ll * 5 + 3] > 0) {
    goto L23;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 22                                                           >*/
    goto L22;
/*<    23 if(kp(3,ll) .gt. 0) go to 24                                       >*/
L23:
    if (kp[ll * 5 + 3] > 0) {
    goto L24;
    }
/*<       lt=lt+1                                                            >*/
    ++lt;
/*<       d(lt,mkp3)=tc(-kp(3,ll))                                           >*/
    d__[lt + mkp3 * d_dim1] = tc[-kp[ll * 5 + 3]];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 22                                                           >*/
    goto L22;
/*<    24 kp3=kp(3,ll)                                                       >*/
L24:
    kp3 = kp[ll * 5 + 3];
/*<       do 26 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       la=la+5*lp(1,l1)*nt                                                >*/
    la += lp[l1 * 3 + 1] * 5 * nt;
/*<       do 25 i=1,nt                                                       >*/
    i__2 = nt;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       lt=lt+1                                                            >*/
        ++lt;
/*<       d(lt,mkp3)=tc(i+la)                                                >*/
        d__[lt + mkp3 * d_dim1] = tc[i__ + la];
/*<    25 continue                                                           >*/
/* L25: */
    }
/*<       la=la+nt                                                           >*/
    la += nt;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<    26 continue                                                           >*/
/* L26: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 22                                                           >*/
    goto L22;
/*<    27 iter=iter+1                                                        >*/
L27:
    ++iter;
/*<       b=0.d0                                                             >*/
    b = 0.;
/*<       sw=b                                                               >*/
    sw = b;
/*<       yb=sw                                                              >*/
    yb = sw;
/*<       do 29 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=a                                                                >*/
    s = a;
/*<       do 28 m=1,mk                                                       >*/
    i__2 = mk;
    for (m = 1; m <= i__2; ++m) {
/*<       s=s+d(m,mkp3)*sc(i,m)                                              >*/
        s += d__[m + mkp3 * d_dim1] * sc[i__ + m * sc_dim1];
/*<    28 continue                                                           >*/
/* L28: */
    }
/*<       sc(i,mkp3)=s                                                       >*/
    sc[i__ + mkp3 * sc_dim1] = s;
/*<       pp=1.0/(1.0+exp(-sc(i,mkp3)))                                      >*/
    pp = (float)1. / (exp(-sc[i__ + mkp3 * sc_dim1]) + (float)1.);
/*<       ww=amax1(pp*(1.0-pp),wm)                                           >*/
/* Computing MAX */
    r__1 = pp * ((float)1. - pp);
    ww = dmax(r__1,wm);
/*<       sc(i,mkp3)=sc(i,mkp3)+(y(i)-pp)/ww                                 >*/
    sc[i__ + mkp3 * sc_dim1] += (y[i__] - pp) / ww;
/*<       if(il.eq.2) ww=ww**2                                               >*/
    if (*il == 2) {
/* Computing 2nd power */
        r__1 = ww;
        ww = r__1 * r__1;
    }
/*<       ww=ww*w(i)                                                         >*/
    ww *= w[i__];
/*<       sc(i,mkp2)=ww                                                      >*/
    sc[i__ + mkp2 * sc_dim1] = ww;
/*<       sw=sw+ww                                                           >*/
    sw += ww;
/*<       yb=yb+ww*sc(i,mkp3)                                                >*/
    yb += ww * sc[i__ + mkp3 * sc_dim1];
/*<       if(iter.gt.1) b=b+abs(pp-sc(i,mkp1))                               >*/
    if (iter > 1) {
        b += (r__1 = pp - sc[i__ + mkp1 * sc_dim1], dabs(r__1));
    }
/*<       sc(i,mkp1)=pp                                                      >*/
    sc[i__ + mkp1 * sc_dim1] = pp;
/*<    29 continue                                                           >*/
/* L29: */
    }
/*<       if(iter.gt.niter.or.(iter.gt.1.and.b/n.lt.thr)) go to 37           >*/
    if (iter > niter || iter > 1 && b / *n < thr) {
    goto L37;
    }
/*<       yb=yb/sw                                                           >*/
    yb /= sw;
/*<       do 36 m=1,mk                                                       >*/
    i__1 = mk;
    for (m = 1; m <= i__1; ++m) {
/*<       b=0.d0                                                             >*/
    b = 0.;
/*<       do 30 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       b=b+sc(i,mkp2)*sc(i,m)                                             >*/
        b += sc[i__ + mkp2 * sc_dim1] * sc[i__ + m * sc_dim1];
/*<    30 continue                                                           >*/
/* L30: */
    }
/*<       b=b/sw                                                             >*/
    b /= sw;
/*<       mm1=m-1                                                            >*/
    mm1 = m - 1;
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 32                                                           >*/
    goto L32;
/*<    31 l=l+1                                                              >*/
L31:
    ++l;
/*<    32 if((l).gt.(mm1)) go to 34                                          >*/
L32:
    if (l > mm1) {
        goto L34;
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 33 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+sc(i,mkp2)*(sc(i,m)-b)*sc(i,l)                                 >*/
        s += sc[i__ + mkp2 * sc_dim1] * (sc[i__ + m * sc_dim1] - b) * sc[
            i__ + l * sc_dim1];
/*<    33 continue                                                           >*/
/* L33: */
    }
/*<       d(l,m)=s                                                           >*/
    d__[l + m * d_dim1] = s;
/*<       go to 31                                                           >*/
    goto L31;
/*<    34 a=0.d0                                                             >*/
L34:
    a = 0.;
/*<       s=a                                                                >*/
    s = a;
/*<       do 35 i=1,n                                                        >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       ww=sc(i,mkp2)                                                      >*/
        ww = sc[i__ + mkp2 * sc_dim1];
/*<       pp=sc(i,m)-b                                                       >*/
        pp = sc[i__ + m * sc_dim1] - b;
/*<       s=s+ww*pp**2                                                       >*/
/* Computing 2nd power */
        r__1 = pp;
        s += ww * (r__1 * r__1);
/*<       a=a+ww*pp*sc(i,mkp3)                                               >*/
        a += ww * pp * sc[i__ + mkp3 * sc_dim1];
/*<    35 continue                                                           >*/
/* L35: */
    }
/*<       d(m,m)=s                                                           >*/
    d__[m + m * d_dim1] = s;
/*<       d(m,mkp1)=a                                                        >*/
    d__[m + mkp1 * d_dim1] = a;
/*<       d(m,mkp2)=b                                                        >*/
    d__[m + mkp2 * d_dim1] = b;
/*<    36 continue                                                           >*/
/* L36: */
    }
/*<       call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)                >*/
    lsf_(nk, &mk, &mkp1, &yb, &d__[d_offset], &d__[mkp3 * d_dim1 + 1], &a, &s,
         &d__[mkp4 * d_dim1 + 1], &c__1);
/*<       go to 27                                                           >*/
    goto L27;
/*<    37 if(jnt .ne. 1) go to 39                                            >*/
L37:
    if (jnt != 1) {
    goto L39;
    }
/*<       az=a                                                               >*/
    *az = a;
/*<       k=0                                                                >*/
    k = 0;
/*<       do 38 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 38                                        >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L38;
    }
/*<       k=k+1                                                              >*/
    ++k;
/*<       tb(1,m)=d(k,mkp3)                                                  >*/
    tb[m * 5 + 1] = d__[k + mkp3 * d_dim1];
/*<    38 continue                                                           >*/
L38:
    ;
    }
/*<       go to 45                                                           >*/
    goto L45;
/*<    39 bz=a                                                               >*/
L39:
    *bz = a;
/*<       ll=1                                                               >*/
    ll = 1;
/*<       l1=ll                                                              >*/
    l1 = ll;
/*<       la=0                                                               >*/
    la = 0;
/*<       lt=la                                                              >*/
    lt = la;
/*<    40 if(kp(1,ll).lt.0) go to 45                                         >*/
L40:
    if (kp[ll * 5 + 1] < 0) {
    goto L45;
    }
/*<       if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 41                  >*/
    if (kp[ll * 5 + 1] != 0 || kp[ll * 5 + 3] > 0) {
    goto L41;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 40                                                           >*/
    goto L40;
/*<    41 if(kp(3,ll) .gt. 0) go to 42                                       >*/
L41:
    if (kp[ll * 5 + 3] > 0) {
    goto L42;
    }
/*<       lt=lt+1                                                            >*/
    ++lt;
/*<       tc(-kp(3,ll))=d(lt,mkp3)                                           >*/
    tc[-kp[ll * 5 + 3]] = d__[lt + mkp3 * d_dim1];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 40                                                           >*/
    goto L40;
/*<    42 kp3=kp(3,ll)                                                       >*/
L42:
    kp3 = kp[ll * 5 + 3];
/*<       do 44 m=1,kp3                                                      >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       la=la+5*lp(1,l1)*nt                                                >*/
    la += lp[l1 * 3 + 1] * 5 * nt;
/*<       do 43 i=1,nt                                                       >*/
    i__2 = nt;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       lt=lt+1                                                            >*/
        ++lt;
/*<       tc(i+la)=d(lt,mkp3)                                                >*/
        tc[i__ + la] = d__[lt + mkp3 * d_dim1];
/*<    43 continue                                                           >*/
/* L43: */
    }
/*<       la=la+nt                                                           >*/
    la += nt;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<    44 continue                                                           >*/
/* L44: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 40                                                           >*/
    goto L40;
/*<    45 return                                                             >*/
L45:
    return 0;
/*<       end                                                                >*/
} /* logitl_ */

/* Subroutine */ int logitl_(n, x, y, w, nk, il, az, tb, cm, sc, d__)
integer *n;
real *x, *y, *w;
integer *nk, *il;
real *az, *tb, *cm, *sc;
doublereal *d__;
{
    return logitl_0_(0, n, x, y, w, nk, il, az, tb, cm, sc, d__, (integer *)0,
         (integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0);
    }

/* Subroutine */ int logitc_(n, x, y, w, nk, il, cm, kp, kv, lp, lv, bz, tc, 
    sc, ss, d__)
integer *n;
real *x, *y, *w;
integer *nk, *il;
real *cm;
integer *kp, *kv, *lp, *lv;
real *bz, *tc, *sc, *ss;
doublereal *d__;
{
    return logitl_0_(1, n, x, y, w, nk, il, (real *)0, (real *)0, cm, sc, d__,
         kp, kv, lp, lv, bz, tc, ss);
    }

/*<       subroutine varimp (n,p,x,y,w,nk,il,it,az,tb,cm,vip,sc,d)           >*/
/* Subroutine */ int varimp_(n, p, x, y, w, nk, il, it, az, tb, cm, vip, sc, 
    d__)
integer *n, *p;
real *x, *y, *w;
integer *nk, *il, *it;
real *az, *tb, *cm, *vip, *sc;
doublereal *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int varz_();
    static real g;
    static integer i__, j;
    static real a0, g0;
    static integer nd, ip;
    static doublereal yb, wn;
    extern /* Subroutine */ int vp_();
    static doublereal sw, yv;
    extern /* Subroutine */ int numprt_();
    static real cst;

/*<       integer p                                                          >*/
/*<       real x(n,p),y(n),w(n),tb(5,nk),cm(*),vip(p),sc(n,*)                >*/
/*<       double precision d(nk,*),sw,yb,yv,wn                               >*/
/*<       sw=0.d0                                                            >*/
    /* Parameter adjustments */
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --w;
    --y;
    --vip;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    tb -= 6;
    --cm;

    /* Function Body */
    sw = 0.;
/*<       yb=sw                                                              >*/
    yb = sw;
/*<       yv=yb                                                              >*/
    yv = yb;
/*<       wn=yv                                                              >*/
    wn = yv;
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       sw=sw+w(i)                                                         >*/
    sw += w[i__];
/*<       wn=wn+w(i)**2                                                      >*/
/* Computing 2nd power */
    r__1 = w[i__];
    wn += r__1 * r__1;
/*<       yb=yb+w(i)*y(i)                                                    >*/
    yb += w[i__] * y[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       yb=yb/sw                                                           >*/
    yb /= sw;
/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       yv=yv+w(i)*(y(i)-yb)**2                                            >*/
/* Computing 2nd power */
    d__1 = y[i__] - yb;
    yv += w[i__] * (d__1 * d__1);
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       yv=yv/sw                                                           >*/
    yv /= sw;
/*<       wn=sw**2/wn                                                        >*/
/* Computing 2nd power */
    d__1 = sw;
    wn = d__1 * d__1 / wn;
/*<       ip=nk+4                                                            >*/
    ip = *nk + 4;
/*<       call varz(0,nk,tb,sc(1,ip),cst,nd)                                 >*/
    varz_(&c__0, nk, &tb[6], &sc[ip * sc_dim1 + 1], &cst, &nd);
/*<       if(cst .ne. 1.0) go to 3                                           >*/
    if (cst != (float)1.) {
    goto L3;
    }
/*<       g0=0.0                                                             >*/
    g0 = (float)0.;
/*<       if(il.gt.0) g0=yv                                                  >*/
    if (*il > 0) {
    g0 = yv;
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<     3 a0=az                                                              >*/
L3:
    a0 = *az;
/*<       call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g0,sc,d)                >*/
    vp_(n, &x[x_offset], &y[1], &w[1], nk, il, &yb, &sw, &a0, &sc[ip * 
        sc_dim1 + 1], &cm[1], &g0, &sc[sc_offset], &d__[d_offset]);
/*<     4 cst=1.d0/(1.d0-cst/wn)**2                                          >*/
L4:
/* Computing 2nd power */
    d__1 = 1. - cst / wn;
    cst = 1. / (d__1 * d__1);
/*<       if(il .ne. 0) go to 5                                              >*/
    if (*il != 0) {
    goto L5;
    }
/*<       g0=(g0+yv)*cst                                                     >*/
    g0 = (g0 + yv) * cst;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 g0=g0*cst                                                          >*/
L5:
    g0 *= cst;
/*<     6 do 12 j=1,p                                                        >*/
L6:
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       call varz(j,nk,tb,sc(1,ip),cst,nd)                                 >*/
    varz_(&j, nk, &tb[6], &sc[ip * sc_dim1 + 1], &cst, &nd);
/*<       if(nd .ne. 0) go to 7                                              >*/
    if (nd != 0) {
        goto L7;
    }
/*<       vip(j)=g0                                                          >*/
    vip[j] = g0;
/*<       go to 12                                                           >*/
    goto L12;
/*<     7 if(cst .ne. 1.0) go to 8                                           >*/
L7:
    if (cst != (float)1.) {
        goto L8;
    }
/*<       g=0.0                                                              >*/
    g = (float)0.;
/*<       if(il.gt.0) g=yv                                                   >*/
    if (*il > 0) {
        g = yv;
    }
/*<       go to 9                                                            >*/
    goto L9;
/*<     8 a0=az                                                              >*/
L8:
    a0 = *az;
/*<       call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g,sc,d)                 >*/
    vp_(n, &x[x_offset], &y[1], &w[1], nk, il, &yb, &sw, &a0, &sc[ip * 
        sc_dim1 + 1], &cm[1], &g, &sc[sc_offset], &d__[d_offset]);
/*<     9 cst=1.d0/(1.d0-cst/wn)**2                                          >*/
L9:
/* Computing 2nd power */
    d__1 = 1. - cst / wn;
    cst = 1. / (d__1 * d__1);
/*<       if(il .ne. 0) go to 10                                             >*/
    if (*il != 0) {
        goto L10;
    }
/*<       g=(g+yv)*cst                                                       >*/
    g = (g + yv) * cst;
/*<       go to 11                                                           >*/
    goto L11;
/*<    10 g=g*cst                                                            >*/
L10:
    g *= cst;
/*<    11 vip(j)=g                                                           >*/
L11:
    vip[j] = g;
/*<    12 continue                                                           >*/
L12:
    ;
    }
/*<       if(it .le. 0) go to 13                                             >*/
    if (*it <= 0) {
    goto L13;
    }
/*    write(it,17)                                                       2
243*/
/*<       call numprt(it,p,vip)                                              >*/
    numprt_(it, p, &vip[1]);
/*<    13 a0=0.0                                                             >*/
L13:
    a0 = (float)0.;
/*<       do 14 j=1,p                                                        >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       vip(j)=sqrt(amax1(0.0,vip(j)-g0))                                  >*/
/* Computing MAX */
    r__1 = (float)0., r__2 = vip[j] - g0;
    vip[j] = sqrt((dmax(r__1,r__2)));
/*<       a0=amax1(a0,vip(j))                                                >*/
/* Computing MAX */
    r__1 = a0, r__2 = vip[j];
    a0 = dmax(r__1,r__2);
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       if(a0.le.0.0) return                                               >*/
    if (a0 <= (float)0.) {
    return 0;
    }
/*<       do 15 j=1,p                                                        >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       vip(j)=100.0*vip(j)/a0                                             >*/
    vip[j] = vip[j] * (float)100. / a0;
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<       if(it .le. 0) go to 16                                             >*/
    if (*it <= 0) {
    goto L16;
    }
/*    write(it,18)                                                       2
255*/
/*<       call numprt(it,p,vip)                                              >*/
    numprt_(it, p, &vip[1]);
/*<    16 return                                                             >*/
L16:
    return 0;
/*<    17 format(/,' -gcv removing each variable:')                          >*/
/* L17: */
/*<    18 format(/,' relative variable importance:')                         >*/
/* L18: */
/*<       end                                                                >*/
} /* varimp_ */

/*<       subroutine numprt(it,n,a)                                          >*/
/* Subroutine */ int numprt_(it, n, a)
integer *it, *n;
real *a;
{
    static integer i1, i2;

/*<       real a(*)                                                          >*/
/*<       i2=0                                                               >*/
    /* Parameter adjustments */
    --a;

    /* Function Body */
    i2 = 0;
/*<     1 if(i2.ge.n) go to 2                                                >*/
L1:
    if (i2 >= *n) {
    goto L2;
    }
/*<       i1=i2+1                                                            >*/
    i1 = i2 + 1;
/*<       i2=i2+6                                                            >*/
    i2 += 6;
/*<       if(i2.gt.n) i2=n                                                   >*/
    if (i2 > *n) {
    i2 = *n;
    }
/*    write(it,'(/,'' '',6(''    '',i4,''    ''))') (i,i=i1,i2)          2
268*/
/*    write(it,'('' '',6g12.4)') (a(i),i=i1,i2)                          2
269*/
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 return                                                             >*/
L2:
    return 0;
/*<       end                                                                >*/
} /* numprt_ */

/*<       subroutine varz(j,nk,tb,ub,cst,nd)                                 >*/
/* Subroutine */ int varz_(j, nk, tb, ub, cst, nd)
integer *j, *nk;
real *tb, *ub, *cst;
integer *nd;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, m;
    extern integer jf_();

/*<       real tb(5,nk),ub(5,nk)                                             >*/
/*<       do 2 m=1,nk                                                        >*/
    /* Parameter adjustments */
    ub -= 6;
    tb -= 6;

    /* Function Body */
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       do 1 k=1,5                                                         >*/
    for (k = 1; k <= 5; ++k) {
/*<       ub(k,m)=tb(k,m)                                                    >*/
        ub[k + m * 5] = tb[k + m * 5];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       nd=0                                                               >*/
    *nd = 0;
/*<       if(j .le. 0) go to 4                                               >*/
    if (*j <= 0) {
    goto L4;
    }
/*<       do 3 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(ub(1,m).eq.0.0) go to 3                                         >*/
    if (ub[m * 5 + 1] == (float)0.) {
        goto L3;
    }
/*<       if(jf(m,j,ub) .eq. 0) go to 3                                      >*/
    if (jf_(&m, j, &ub[6]) == 0) {
        goto L3;
    }
/*<       ub(1,m)=0.0                                                        >*/
    ub[m * 5 + 1] = (float)0.;
/*<       nd=nd+1                                                            >*/
    ++(*nd);
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<     4 cst=1.0                                                            >*/
L4:
    *cst = (float)1.;
/*<       do 5 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(ub(1,m).ne.0.0) cst=cst+ub(5,m)                                 >*/
    if (ub[m * 5 + 1] != (float)0.) {
        *cst += ub[m * 5 + 5];
    }
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* varz_ */

/*<       subroutine vp (n,x,y,w,nk,il,yb,sw,az,tb,cm,gof,sc,d)              >*/
/* Subroutine */ int vp_(n, x, y, w, nk, il, yb, sw, az, tb, cm, gof, sc, d__)
integer *n;
real *x, *y, *w;
integer *nk, *il;
doublereal *yb, *sw;
real *az, *tb, *cm, *gof, *sc;
doublereal *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1, 
        i__2;
    real r__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real a;
    static integer i__, k, m;
    static doublereal s, t;
    static real pp;
    extern /* Subroutine */ int logitl_(), lstsqr_();

/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)                      >*/
/*<       double precision d(nk,*),s,t,yb,sw                                 >*/
/*<       if(il .ne. 0) go to 1                                              >*/
    /* Parameter adjustments */
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    tb -= 6;
    --cm;

    /* Function Body */
    if (*il != 0) {
    goto L1;
    }
/*<       call lstsqr(n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)                       >*/
    lstsqr_(n, &x[x_offset], &y[1], &w[1], nk, yb, sw, &tb[6], &cm[1], gof, &
        sc[sc_offset], &d__[d_offset]);
/*<       return                                                             >*/
    return 0;
/*<     1 call logitl(n,x,y,w,nk,il,az,tb,cm,sc,d)                           >*/
L1:
    logitl_(n, &x[x_offset], &y[1], &w[1], nk, il, az, &tb[6], &cm[1], &sc[
        sc_offset], &d__[d_offset]);
/*<       t=0.d0                                                             >*/
    t = 0.;
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=az                                                               >*/
    s = *az;
/*<       k=0                                                                >*/
    k = 0;
/*<       do 2 m=1,nk                                                        >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 2                                         >*/
        if (tb[m * 5 + 1] == (float)0.) {
        goto L2;
        }
/*<       k=k+1                                                              >*/
        ++k;
/*<       s=s+tb(1,m)*sc(i,k)                                                >*/
        s += tb[m * 5 + 1] * sc[i__ + k * sc_dim1];
/*<     2 continue                                                           >*/
L2:
        ;
    }
/*<       a=s                                                                >*/
    a = s;
/*<       pp=1.0/(1.0+exp(-a))                                               >*/
    pp = (float)1. / (exp(-a) + (float)1.);
/*<       t=t+w(i)*(y(i)-pp)**2                                              >*/
/* Computing 2nd power */
    r__1 = y[i__] - pp;
    t += w[i__] * (r__1 * r__1);
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       gof=t/sw                                                           >*/
    *gof = t / *sw;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* vp_ */

/*<       subroutine lstsqr (n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)                >*/
/* Subroutine */ int lstsqr_(n, x, y, w, nk, yb, sw, tb, cm, gof, sc, d__)
integer *n;
real *x, *y, *w;
integer *nk;
doublereal *yb, *sw;
real *tb, *cm, *gof, *sc;
doublereal *d__;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1, 
        i__2;
    real r__1;

    /* Local variables */
    static doublereal a, b;
    static integer i__, k, l, m;
    static doublereal s;
    static integer mk;
    static real pp, ww;
    static integer mm1;
    extern doublereal phi_();
    extern /* Subroutine */ int lsf_();
    static integer mkp1, mkp2, mkp3, mkp4;

/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)                      >*/
/*<       double precision d(nk,*),a,b,s,yb,sw                               >*/
/*<       do 2 i=1,n                                                         >*/
    /* Parameter adjustments */
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    tb -= 6;
    --cm;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       k=0                                                                >*/
    k = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 1                                         >*/
        if (tb[m * 5 + 1] == (float)0.) {
        goto L1;
        }
/*<       k=k+1                                                              >*/
        ++k;
/*<       sc(i,k)=phi(m,i,n,x,tb,cm)                                         >*/
        sc[i__ + k * sc_dim1] = phi_(&m, &i__, n, &x[x_offset], &tb[6], &
            cm[1]);
/*<     1 continue                                                           >*/
L1:
        ;
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       mk=k                                                               >*/
    mk = k;
/*<       mkp1=mk+1                                                          >*/
    mkp1 = mk + 1;
/*<       mkp2=mk+2                                                          >*/
    mkp2 = mk + 2;
/*<       mkp3=mk+3                                                          >*/
    mkp3 = mk + 3;
/*<       mkp4=mk+4                                                          >*/
    mkp4 = mk + 4;
/*<       do 9 m=1,mk                                                        >*/
    i__1 = mk;
    for (m = 1; m <= i__1; ++m) {
/*<       b=0.d0                                                             >*/
    b = 0.;
/*<       do 3 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       b=b+w(i)*sc(i,m)                                                   >*/
        b += w[i__] * sc[i__ + m * sc_dim1];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       b=b/sw                                                             >*/
    b /= *sw;
/*<       mm1=m-1                                                            >*/
    mm1 = m - 1;
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 l=l+1                                                              >*/
L4:
    ++l;
/*<     5 if((l).gt.(mm1)) go to 7                                           >*/
L5:
    if (l > mm1) {
        goto L7;
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 6 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+w(i)*(sc(i,m)-b)*sc(i,l)                                       >*/
        s += w[i__] * (sc[i__ + m * sc_dim1] - b) * sc[i__ + l * sc_dim1];
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       d(l,m)=s                                                           >*/
    d__[l + m * d_dim1] = s;
/*<       go to 4                                                            >*/
    goto L4;
/*<     7 a=0.d0                                                             >*/
L7:
    a = 0.;
/*<       s=a                                                                >*/
    s = a;
/*<       do 8 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       ww=w(i)                                                            >*/
        ww = w[i__];
/*<       pp=sc(i,m)-b                                                       >*/
        pp = sc[i__ + m * sc_dim1] - b;
/*<       s=s+ww*pp**2                                                       >*/
/* Computing 2nd power */
        r__1 = pp;
        s += ww * (r__1 * r__1);
/*<       a=a+ww*pp*y(i)                                                     >*/
        a += ww * pp * y[i__];
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       d(m,m)=s                                                           >*/
    d__[m + m * d_dim1] = s;
/*<       d(m,mkp1)=a                                                        >*/
    d__[m + mkp1 * d_dim1] = a;
/*<       d(m,mkp2)=b                                                        >*/
    d__[m + mkp2 * d_dim1] = b;
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)                >*/
    lsf_(nk, &mk, &mkp1, yb, &d__[d_offset], &d__[mkp3 * d_dim1 + 1], &a, &s, 
        &d__[mkp4 * d_dim1 + 1], &c__1);
/*<       gof=s/sw                                                           >*/
    *gof = s / *sw;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* lstsqr_ */

/*<       function efp(l,jv,nk,tb)                                           >*/
doublereal efp_(l, jv, nk, tb)
integer *l, *jv, *nk;
real *tb;
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    extern integer nord_();
    static integer j, k, m;
    extern integer jf_();

/*<       integer jv(l)                                                      >*/
/*<       real tb(5,nk)                                                      >*/
/*<       efp=0.0                                                            >*/
    /* Parameter adjustments */
    --jv;
    tb -= 6;

    /* Function Body */
    ret_val = (float)0.;
/*<       do 3 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 3                                         >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L3;
    }
/*<       if(nord(m,tb).ne.l) go to 3                                        >*/
    if (nord_(&m, &tb[6]) != *l) {
        goto L3;
    }
/*<       k=0                                                                >*/
    k = 0;
/*<       do 1 j=1,l                                                         >*/
    i__2 = *l;
    for (j = 1; j <= i__2; ++j) {
/*<       if(jf(m,jv(j),tb).eq.1) go to 1                                    >*/
        if (jf_(&m, &jv[j], &tb[6]) == 1) {
        goto L1;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 continue                                                           >*/
L1:
        ;
    }
/*<     2 if(k.eq.1) go to 3                                                 >*/
L2:
    if (k == 1) {
        goto L3;
    }
/*<       efp=efp+tb(5,m)                                                    >*/
    ret_val += tb[m * 5 + 5];
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* efp_ */

/*<       function elg(jv,l,lx,tb,cm)                                        >*/
VOID elg_0_(n__, ret_val, jv, l, lx, tb, cm, i1)
int n__;
Multitype *ret_val;
integer *jv, *l, *lx;
real *tb, *cm;
integer *i1;
{
    /* Initialized data */

    static integer ic = 0;

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer k;
    extern integer nordc_(), nnord_();
    static integer jb, jl, ip, kx;
    extern /* Subroutine */ int intalw_(), isnstr_();

/*<       real tb(5,*),cm(*)                                                 >*/
/*<       integer lx(*)                                                      >*/
/*<       logical elg                                                        >*/
/*<       data ic /0/                                                        >*/
    /* Parameter adjustments */
    if (lx) {
    --lx;
    }
    if (tb) {
    tb -= 6;
    }
    if (cm) {
    --cm;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_stelg;
    }

/*<       elg=.false.                                                        >*/
    (*ret_val).i = FALSE_;
/*<       kx=iabs(lx(jv))                                                    >*/
    kx = (i__1 = lx[*jv], abs(i__1));
/*<       if(kx.eq.0) return                                                 >*/
    if (kx == 0) {
    return ;
    }
/*<       if(l .ne. 0) go to 1                                               >*/
    if (*l != 0) {
    goto L1;
    }
/*<       elg=.true.                                                         >*/
    (*ret_val).i = TRUE_;
/*<       return                                                             >*/
    return ;
/*<     1 if((kx .ne. 2) .and. (kx .ne. 3)) go to 2                          >*/
L1:
    if (kx != 2 && kx != 3) {
    goto L2;
    }
/*<       if(nnord(l,tb).gt.0) return                                        >*/
    if (nnord_(l, &tb[6]) > 0) {
    return ;
    }
/*<     2 ip=l                                                               >*/
L2:
    ip = *l;
/*<     3 if(ip.le.0) go to 4                                                >*/
L3:
    if (ip <= 0) {
    goto L4;
    }
/*<       jl=abs(tb(2,ip))+.1                                                >*/
    jl = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     4 k=iabs(lx(jl))                                                     >*/
L4:
    k = (i__1 = lx[jl], abs(i__1));
/*<       call isnstr(jl,jb)                                                 >*/
    isnstr_(&jl, &jb);
/*<       if((k.eq.2.or.k.eq.3).and.jb.eq.0) return                          >*/
    if ((k == 2 || k == 3) && jb == 0) {
    return ;
    }
/*<       if(ic .ne. 1) go to 5                                              >*/
    if (ic != 1) {
    goto L5;
    }
/*<       if(lx(jv).lt.0.and.nordc(1,l,tb,cm).gt.0) return                   >*/
    if (lx[*jv] < 0 && nordc_(&c__1, l, &tb[6], &cm[1]) > 0) {
    return ;
    }
/*<       if(lx(jv).gt.0.and.nordc(2,l,tb,cm).gt.0) return                   >*/
    if (lx[*jv] > 0 && nordc_(&c__2, l, &tb[6], &cm[1]) > 0) {
    return ;
    }
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 if(ic .ne. 2) go to 6                                              >*/
L5:
    if (ic != 2) {
    goto L6;
    }
/*<       if(lx(jv).gt.0.and.nordc(1,l,tb,cm).ge.2) return                   >*/
    if (lx[*jv] > 0 && nordc_(&c__1, l, &tb[6], &cm[1]) >= 2) {
    return ;
    }
/*<     6 ip=l                                                               >*/
L6:
    ip = *l;
/*<     7 if(ip.le.0) go to 8                                                >*/
L7:
    if (ip <= 0) {
    goto L8;
    }
/*<       jl=abs(tb(2,ip))+.1                                                >*/
    jl = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       call intalw(jv,jl,k)                                               >*/
    intalw_(jv, &jl, &k);
/*<       if(k.eq.0) return                                                  >*/
    if (k == 0) {
    return ;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 7                                                            >*/
    goto L7;
/*<     8 elg=.true.                                                         >*/
L8:
    (*ret_val).i = TRUE_;
/*<       return                                                             >*/
    return ;
/*<       entry stelg(i1)                                                    >*/

L_stelg:
/*<       ic=i1                                                              >*/
    ic = *i1;
/*<       return                                                             >*/
    return ;
/*<       end                                                                >*/
} /* elg_ */

logical elg_(jv, l, lx, tb, cm)
integer *jv, *l, *lx;
real *tb, *cm;
{
    Multitype ret_val;
    elg_0_(0, &ret_val, jv, l, lx, tb, cm, (integer *)0);
    return ret_val.i;
    }

doublereal stelg_(i1)
integer *i1;
{
    Multitype ret_val;
    elg_0_(1, &ret_val, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, i1);
    return ret_val.r;
    }

/*<       function phi(m,i,n,x,tb,cm)                                        >*/
doublereal phi_(m, i__, n, x, tb, cm)
integer *m, *i__, *n;
real *x, *tb, *cm;
{
    /* System generated locals */
    integer x_dim1, x_offset;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    static integer j;
    static real t, u;
    static integer ip;

/*<       real tb(5,*),x(n,*),cm(*)                                          >*/
/*<       phi=1.0                                                            >*/
    /* Parameter adjustments */
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    tb -= 6;
    --cm;

    /* Function Body */
    ret_val = (float)1.;
/*<       ip=m                                                               >*/
    ip = *m;
/*<     1 if(ip.le.0) go to 7                                                >*/
L1:
    if (ip <= 0) {
    goto L7;
    }
/*<       t=tb(2,ip)                                                         >*/
    t = tb[ip * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(cm(2*j) .le. 0.0) go to 4                                       >*/
    if (cm[j * 2] <= (float)0.) {
    goto L4;
    }
/*<       u=cm(int(x(i,j)+.1)+int(tb(3,ip)+.1))                              >*/
    u = cm[(integer) (x[*i__ + j * x_dim1] + (float).1) + (integer) (tb[ip * 
        5 + 3] + (float).1)];
/*<       if(t .ge. 0.0) go to 5                                             >*/
    if (t >= (float)0.) {
    goto L5;
    }
/*<       if(u .ne. 0.0) go to 2                                             >*/
    if (u != (float)0.) {
    goto L2;
    }
/*<       u=1.0                                                              >*/
    u = (float)1.;
/*<       go to 5                                                            >*/
    goto L5;
/*<     2 u=0.0                                                              >*/
L2:
    u = (float)0.;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,ip)))                         >*/
L4:
/* Computing MAX */
    r__1 = (float)0., r__2 = r_sign(&c_b196, &t) * (x[*i__ + j * x_dim1] - tb[
        ip * 5 + 3]);
    u = dmax(r__1,r__2);
/*<     5 if(u .gt. 0.0) go to 6                                             >*/
L5:
    if (u > (float)0.) {
    goto L6;
    }
/*<       phi=0.0                                                            >*/
    ret_val = (float)0.;
/*<       return                                                             >*/
    return ret_val;
/*<     6 phi=phi*u                                                          >*/
L6:
    ret_val *= u;
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     7 return                                                             >*/
L7:
    return ret_val;
/*<       end                                                                >*/
} /* phi_ */

/*<       function nord(m,tb)                                                >*/
integer nord_(m, tb)
integer *m;
real *tb;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer ip;

/*<       real tb(5,*)                                                       >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       nord=0                                                             >*/
    ret_val = 0;
/*<     1 if(ip.le.0) go to 2                                                >*/
L1:
    if (ip <= 0) {
    goto L2;
    }
/*<       nord=nord+1                                                        >*/
    ++ret_val;
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* nord_ */

/*<       function jf(m,j,tb)                                                >*/
integer jf_(m, j, tb)
integer *m, *j;
real *tb;
{
    /* System generated locals */
    integer ret_val;
    real r__1;

    /* Local variables */
    static integer ip, jp;

/*<       real tb(5,*)                                                       >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       jf=0                                                               >*/
    ret_val = 0;
/*<     1 if(ip.le.0) go to 2                                                >*/
L1:
    if (ip <= 0) {
    goto L2;
    }
/*<       jp=abs(tb(2,ip))+.1                                                >*/
    jp = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(jp.eq.j) jf=1                                                   >*/
    if (jp == *j) {
    ret_val = 1;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* jf_ */

/*<       subroutine lsf(nk,m,mkp1,yb,d,a,a0,gf,dp,k1)                       >*/
/* Subroutine */ int lsf_0_(n__, nk, m, mkp1, yb, d__, a, a0, gf, dp, k1, dps)
int n__;
integer *nk, *m, *mkp1;
doublereal *yb, *d__, *a, *a0, *gf, *dp;
integer *k1;
doublereal *dps;
{
    /* Initialized data */

    static real big = (float)9.9e30;
    static doublereal eps = 1e-5;

    /* System generated locals */
    integer d_dim1, d_offset, dp_dim1, dp_offset, i__1, i__2;

    /* Local variables */
    static integer info, i__, k;
    static doublereal s, t;
    extern /* Subroutine */ int spofa_(), sposl_();
    static doublereal sl;
    static integer mkp2;

/*<       double precision a(m),d(nk,*),dp(nk,*),eps,yb,a0,gf,s,t,sl,dps     >*/
/*<       data big,eps /9.9e30,1.d-05/                                       >*/
    /* Parameter adjustments */
    if (d__) {
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    }
    if (dp) {
    dp_dim1 = *nk;
    dp_offset = dp_dim1 + 1;
    dp -= dp_offset;
    }
    if (a) {
    --a;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_seteps;
    }

/*<       mkp2=mkp1+1                                                        >*/
    mkp2 = *mkp1 + 1;
/*<       gf=big                                                             >*/
    *gf = big;
/*<       if(d(m,m).le.0.d0) return                                          >*/
    if (d__[*m + *m * d_dim1] <= 0.) {
    return 0;
    }
/*<       sl=1.d0+eps                                                        >*/
    sl = eps + 1.;
/*<       do 2 k=k1,m                                                        >*/
    i__1 = *m;
    for (k = *k1; k <= i__1; ++k) {
/*<       do 1 i=1,k                                                         >*/
    i__2 = k;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       dp(i,k)=d(i,k)                                                     >*/
        dp[i__ + k * dp_dim1] = d__[i__ + k * d_dim1];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       dp(k,k)=dp(k,k)*sl                                                 >*/
    dp[k + k * dp_dim1] *= sl;
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       do 3 k=1,m                                                         >*/
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
/*<       a(k)=d(k,mkp1)                                                     >*/
    a[k] = d__[k + *mkp1 * d_dim1];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       info=k1                                                            >*/
    info = *k1;
/*<       call spofa(dp,nk,m,info)                                           >*/
    spofa_(&dp[dp_offset], nk, m, &info);
/*<       if(info.ne.0) return                                               >*/
    if (info != 0) {
    return 0;
    }
/*<       call sposl(dp,nk,m,a)                                              >*/
    sposl_(&dp[dp_offset], nk, m, &a[1]);
/*<       s=yb                                                               >*/
    s = *yb;
/*<       t=0.d0                                                             >*/
    t = 0.;
/*<       do 4 i=1,m                                                         >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       s=s-a(i)*d(i,mkp2)                                                 >*/
    s -= a[i__] * d__[i__ + mkp2 * d_dim1];
/*<       t=t-a(i)*(d(i,mkp1)+eps*d(i,i)*a(i))                               >*/
    t -= a[i__] * (d__[i__ + *mkp1 * d_dim1] + eps * d__[i__ + i__ * 
        d_dim1] * a[i__]);
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       a0=s                                                               >*/
    *a0 = s;
/*<       gf=t                                                               >*/
    *gf = t;
/*<       return                                                             >*/
    return 0;
/*<       entry seteps(dps)                                                  >*/

L_seteps:
/*<       eps=dps                                                            >*/
    eps = *dps;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* lsf_ */

/* Subroutine */ int lsf_(nk, m, mkp1, yb, d__, a, a0, gf, dp, k1)
integer *nk, *m, *mkp1;
doublereal *yb, *d__, *a, *a0, *gf, *dp;
integer *k1;
{
    return lsf_0_(0, nk, m, mkp1, yb, d__, a, a0, gf, dp, k1, (doublereal *)0)
        ;
    }

/* Subroutine */ int seteps_(dps)
doublereal *dps;
{
    return lsf_0_(1, (integer *)0, (integer *)0, (integer *)0, (doublereal *)
        0, (doublereal *)0, (doublereal *)0, (doublereal *)0, (doublereal 
        *)0, (doublereal *)0, (integer *)0, dps);
    }

/*<       subroutine coll(nk,tb,lp,lv,jv)                                    >*/
/* Subroutine */ int coll_(nk, tb, lp, lv, jv)
integer *nk;
real *tb;
integer *lp, *lv, *jv;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern integer nord_();
    static integer i__, j, k, m, l1, l2, l10, ig, jg, mo, mt, l1m1;
    extern /* Subroutine */ int jfv_();

/*<       integer lp(3,*),lv(*),jv(*)                                        >*/
/*<       real tb(5,nk)                                                      >*/
/*<       mo=0                                                               >*/
    /* Parameter adjustments */
    tb -= 6;
    lp -= 4;
    --lv;
    --jv;

    /* Function Body */
    mo = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).ne.0.0) mo=max0(mo,nord(m,tb))                          >*/
    if (tb[m * 5 + 1] != (float)0.) {
/* Computing MAX */
        i__2 = mo, i__3 = nord_(&m, &tb[6]);
        mo = max(i__2,i__3);
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       if(mo .ne. 0) go to 2                                              >*/
    if (mo != 0) {
    goto L2;
    }
/*<       lp(1,1)=0                                                          >*/
    lp[4] = 0;
/*<       return                                                             >*/
    return 0;
/*<     2 l1=1                                                               >*/
L2:
    l1 = 1;
/*<       l2=l1                                                              >*/
    l2 = l1;
/*<       do 11 mt=1,mo                                                      >*/
    i__1 = mo;
    for (mt = 1; mt <= i__1; ++mt) {
/*<       l10=l1                                                             >*/
    l10 = l1;
/*<       do 10 m=1,nk                                                       >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(tb(1,m).eq.0.0.or.nord(m,tb).ne.mt) go to 10                    >*/
        if (tb[m * 5 + 1] == (float)0. || nord_(&m, &tb[6]) != mt) {
        goto L10;
        }
/*<       call jfv(m,tb,jv)                                                  >*/
        jfv_(&m, &tb[6], &jv[1]);
/*<       jg=0                                                               >*/
        jg = 0;
/*<       l1m1=l1-1                                                          >*/
        l1m1 = l1 - 1;
/*<       i=l10                                                              >*/
        i__ = l10;
/*<       go to 4                                                            >*/
        goto L4;
/*<     3 i=i+1                                                              >*/
L3:
        ++i__;
/*<     4 if((i).gt.(l1m1)) go to 8                                          >*/
L4:
        if (i__ > l1m1) {
        goto L8;
        }
/*<       k=lp(2,i)-1                                                        >*/
        k = lp[i__ * 3 + 2] - 1;
/*<       ig=0                                                               >*/
        ig = 0;
/*<       do 5 j=1,mt                                                        >*/
        i__3 = mt;
        for (j = 1; j <= i__3; ++j) {
/*<       if(jv(j).eq.lv(k+j)) go to 5                                       >*/
        if (jv[j] == lv[k + j]) {
            goto L5;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 6                                                            >*/
        goto L6;
/*<     5 continue                                                           >*/
L5:
        ;
        }
/*<     6 if(ig .ne. 0) go to 3                                              >*/
L6:
        if (ig != 0) {
        goto L3;
        }
/*<       jg=1                                                               >*/
        jg = 1;
/*<       lp(3,i)=lp(3,i)+1                                                  >*/
        ++lp[i__ * 3 + 3];
/*<     8 if(jg .ne. 0) go to 10                                             >*/
L8:
        if (jg != 0) {
        goto L10;
        }
/*<       lp(1,l1)=mt                                                        >*/
        lp[l1 * 3 + 1] = mt;
/*<       lp(2,l1)=l2                                                        >*/
        lp[l1 * 3 + 2] = l2;
/*<       lp(3,l1)=1                                                         >*/
        lp[l1 * 3 + 3] = 1;
/*<       k=l2-1                                                             >*/
        k = l2 - 1;
/*<       do 9 i=1,mt                                                        >*/
        i__3 = mt;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       lv(i+k)=jv(i)                                                      >*/
        lv[i__ + k] = jv[i__];
/*<     9 continue                                                           >*/
/* L9: */
        }
/*<       l1=l1+1                                                            >*/
        ++l1;
/*<       l2=l2+mt                                                           >*/
        l2 += mt;
/*<    10 continue                                                           >*/
L10:
        ;
    }
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       lp(1,l1)=0                                                         >*/
    lp[l1 * 3 + 1] = 0;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* coll_ */

/*<       subroutine jfv(m,tb,jv)                                            >*/
/* Subroutine */ int jfv_(m, tb, jv)
integer *m;
real *tb;
integer *jv;
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, j, k, l, ip;

/*<       integer jv(*)                                                      >*/
/*<       real tb(5,*)                                                       >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    --jv;
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       j=0                                                                >*/
    j = 0;
/*<     1 if(ip.le.0) go to 2                                                >*/
L1:
    if (ip <= 0) {
    goto L2;
    }
/*<       j=j+1                                                              >*/
    ++j;
/*<       jv(j)=abs(tb(2,ip))+.1                                             >*/
    jv[j] = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 if(j.eq.1) return                                                  >*/
L2:
    if (j == 1) {
    return 0;
    }
/*<       j=j-1                                                              >*/
    --j;
/*<     3 l=0                                                                >*/
L3:
    l = 0;
/*<       do 4 i=1,j                                                         >*/
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(jv(i) .le. jv(i+1)) go to 4                                     >*/
    if (jv[i__] <= jv[i__ + 1]) {
        goto L4;
    }
/*<       k=jv(i)                                                            >*/
    k = jv[i__];
/*<       jv(i)=jv(i+1)                                                      >*/
    jv[i__] = jv[i__ + 1];
/*<       jv(i+1)=k                                                          >*/
    jv[i__ + 1] = k;
/*<       l=1                                                                >*/
    l = 1;
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<       if(l.eq.0) go to 5                                                 >*/
    if (l == 0) {
    goto L5;
    }
/*<       go to 3                                                            >*/
    goto L3;
/*<     5 return                                                             >*/
L5:
    return 0;
/*<       end                                                                >*/
} /* jfv_ */

/*<       subroutine side(l,nt,jv,xe,x)                                      >*/
/* Subroutine */ int side_(l, nt, jv, xe, x)
integer *l, *nt, *jv;
real *xe, *x;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;

    /* Local variables */
    static real a;
    static integer j, k, m;
    static real z__;
    static integer l2, l3, l4;
    static real x1, x2, dl, dr, dx;
    extern /* Subroutine */ int pr_();
    static real xl, xr;

/*<       integer jv(l)                                                      >*/
/*<       real xe(2,*),x(nt,*)                                               >*/
/*<       l2=l+l                                                             >*/
    /* Parameter adjustments */
    --jv;
    x_dim1 = *nt;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    xe -= 3;

    /* Function Body */
    l2 = *l + *l;
/*<       l3=l2+l                                                            >*/
    l3 = l2 + *l;
/*<       l4=l3+l                                                            >*/
    l4 = l3 + *l;
/*<       do 7 k=1,l                                                         >*/
    i__1 = *l;
    for (k = 1; k <= i__1; ++k) {
/*<       xl=xe(1,jv(k))                                                     >*/
    xl = xe[(jv[k] << 1) + 1];
/*<       xr=xe(2,jv(k))                                                     >*/
    xr = xe[(jv[k] << 1) + 2];
/*<       do 6 j=1,nt                                                        >*/
    i__2 = *nt;
    for (j = 1; j <= i__2; ++j) {
/*<       z=x(j,k)                                                           >*/
        z__ = x[j + k * x_dim1];
/*<       if(z .gt. xl) go to 1                                              >*/
        if (z__ > xl) {
        goto L1;
        }
/*<       x(j,k+l)=xl                                                        >*/
        x[j + (k + *l) * x_dim1] = xl;
/*<       x(j,k+l2)=x(j,k+l)                                                 >*/
        x[j + (k + l2) * x_dim1] = x[j + (k + *l) * x_dim1];
/*<       x(j,k+l3)=0.0                                                      >*/
        x[j + (k + l3) * x_dim1] = (float)0.;
/*<       x(j,k+l4)=x(j,k+l3)                                                >*/
        x[j + (k + l4) * x_dim1] = x[j + (k + l3) * x_dim1];
/*<       go to 6                                                            >*/
        goto L6;
/*<     1 dl=z-xl                                                            >*/
L1:
        dl = z__ - xl;
/*<       dr=xr-z                                                            >*/
        dr = xr - z__;
/*<       x1=xl                                                              >*/
        x1 = xl;
/*<       x2=xr                                                              >*/
        x2 = xr;
/*<       do 3 m=1,nt                                                        >*/
        i__3 = *nt;
        for (m = 1; m <= i__3; ++m) {
/*<       a=x(m,k)                                                           >*/
        a = x[m + k * x_dim1];
/*<       if(a.eq.z) go to 3                                                 >*/
        if (a == z__) {
            goto L3;
        }
/*<       dx=a-z                                                             >*/
        dx = a - z__;
/*<       if(dx .ge. 0.0 .or. -dx .ge. dl) go to 2                           >*/
        if (dx >= (float)0. || -dx >= dl) {
            goto L2;
        }
/*<       dl=-dx                                                             >*/
        dl = -dx;
/*<       x1=a                                                               >*/
        x1 = a;
/*<     2 if(dx .le. 0.0 .or. dx .ge. dr) go to 3                            >*/
L2:
        if (dx <= (float)0. || dx >= dr) {
            goto L3;
        }
/*<       dr=dx                                                              >*/
        dr = dx;
/*<       x2=a                                                               >*/
        x2 = a;
/*<     3 continue                                                           >*/
L3:
        ;
        }
/*<       x1=0.5*(x1+z)                                                      >*/
        x1 = (x1 + z__) * (float).5;
/*<       x2=0.5*(x2+z)                                                      >*/
        x2 = (x2 + z__) * (float).5;
/*<       if(x(j,k+l) .le. 0.0) go to 4                                      >*/
        if (x[j + (k + *l) * x_dim1] <= (float)0.) {
        goto L4;
        }
/*<       x(j,k+l)=x1                                                        >*/
        x[j + (k + *l) * x_dim1] = x1;
/*<       x(j,k+l2)=x2                                                       >*/
        x[j + (k + l2) * x_dim1] = x2;
/*<       go to 5                                                            >*/
        goto L5;
/*<     4 x(j,k+l)=x2                                                        >*/
L4:
        x[j + (k + *l) * x_dim1] = x2;
/*<       x(j,k+l2)=x1                                                       >*/
        x[j + (k + l2) * x_dim1] = x1;
/*<     5 call pr(x(j,k+l),x(j,k),x(j,k+l2),x(j,k+l3),x(j,k+l4))             >*/
L5:
        pr_(&x[j + (k + *l) * x_dim1], &x[j + k * x_dim1], &x[j + (k + l2)
             * x_dim1], &x[j + (k + l3) * x_dim1], &x[j + (k + l4) * 
            x_dim1]);
/*<     6 continue                                                           >*/
L6:
        ;
    }
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* side_ */

/*<       subroutine que(jp,l,nt,jv,n,x,tc,t)                                >*/
/* Subroutine */ int que_(jp, l, nt, jv, n, x, tc, t)
integer *jp, *l, *nt, *jv, *n;
real *x, *tc, *t;
{
    /* System generated locals */
    integer x_dim1, x_offset, tc_dim1, tc_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static real q;
    static integer l2, l3, l4;
    extern doublereal cue_();

/*<       integer jv(l)                                                      >*/
/*<       real x(n,*),tc(nt,*),t(n)                                          >*/
/*<       l2=l+l                                                             >*/
    /* Parameter adjustments */
    --jv;
    tc_dim1 = *nt;
    tc_offset = tc_dim1 + 1;
    tc -= tc_offset;
    --t;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    l2 = *l + *l;
/*<       l3=l2+l                                                            >*/
    l3 = l2 + *l;
/*<       l4=l3+l                                                            >*/
    l4 = l3 + *l;
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(t(i).eq.0.0) go to 3                                            >*/
    if (t[i__] == (float)0.) {
        goto L3;
    }
/*<       q=1.0                                                              >*/
    q = (float)1.;
/*<       do 1 k=1,l                                                         >*/
    i__2 = *l;
    for (k = 1; k <= i__2; ++k) {
/*<       j=jv(k)                                                            >*/
        j = jv[k];
/*<    >*/
        q *= cue_(&x[i__ + j * x_dim1], &tc[*jp + (k + *l) * tc_dim1], &
            tc[*jp + k * tc_dim1], &tc[*jp + (k + l2) * tc_dim1], &tc[
            *jp + (k + l3) * tc_dim1], &tc[*jp + (k + l4) * tc_dim1]);
/*<       if(q.eq.0.0) go to 2                                               >*/
        if (q == (float)0.) {
        goto L2;
        }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 t(i)=q                                                             >*/
L2:
    t[i__] = q;
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* que_ */

/*<       subroutine scpc(xm,xs,jp,l,nt,jv,tc,b)                             >*/
/* Subroutine */ int scpc_(xm, xs, jp, l, nt, jv, tc, b)
real *xm, *xs;
integer *jp, *l, *nt, *jv;
real *tc, *b;
{
    /* System generated locals */
    integer tc_dim1, tc_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal g, h__;
    static integer j, k;
    static doublereal q;
    static integer l2, l3, l4;

/*<       integer jv(l)                                                      >*/
/*<       real xm(*),xs(*),tc(nt,*)                                          >*/
/*<       double precision g,h,q                                             >*/
/*<       l2=l+l                                                             >*/
    /* Parameter adjustments */
    --xm;
    --xs;
    --jv;
    tc_dim1 = *nt;
    tc_offset = tc_dim1 + 1;
    tc -= tc_offset;

    /* Function Body */
    l2 = *l + *l;
/*<       l3=l2+l                                                            >*/
    l3 = l2 + *l;
/*<       l4=l3+l                                                            >*/
    l4 = l3 + *l;
/*<       q=1.d0                                                             >*/
    q = 1.;
/*<       do 1 k=1,l                                                         >*/
    i__1 = *l;
    for (k = 1; k <= i__1; ++k) {
/*<       j=jv(k)                                                            >*/
    j = jv[k];
/*<       g=xm(j)                                                            >*/
    g = xm[j];
/*<       h=xs(j)                                                            >*/
    h__ = xs[j];
/*<       q=q*h                                                              >*/
    q *= h__;
/*<       tc(jp,k+l)=g+h*tc(jp,k+l)                                          >*/
    tc[*jp + (k + *l) * tc_dim1] = g + h__ * tc[*jp + (k + *l) * tc_dim1];
/*<       tc(jp,k)=g+h*tc(jp,k)                                              >*/
    tc[*jp + k * tc_dim1] = g + h__ * tc[*jp + k * tc_dim1];
/*<       tc(jp,k+l2)=g+h*tc(jp,k+l2)                                        >*/
    tc[*jp + (k + l2) * tc_dim1] = g + h__ * tc[*jp + (k + l2) * tc_dim1];
/*<       tc(jp,k+l3)=tc(jp,k+l3)/h                                          >*/
    tc[*jp + (k + l3) * tc_dim1] /= h__;
/*<       tc(jp,k+l4)=tc(jp,k+l4)/h**2                                       >*/
/* Computing 2nd power */
    d__1 = h__;
    tc[*jp + (k + l4) * tc_dim1] /= d__1 * d__1;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       b=b/q                                                              >*/
    *b /= q;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* scpc_ */

/*<       subroutine update(il,n,m,kr,x,y,w,sw,yb,tb,cm,sc,bl,d,dy,db)       >*/
/* Subroutine */ int update_(il, n, m, kr, x, y, w, sw, yb, tb, cm, sc, bl, 
    d__, dy, db)
integer *il, *n, *m, *kr;
real *x, *y, *w;
doublereal *sw, *yb;
real *tb, *cm, *sc, *bl;
doublereal *d__, *dy, *db;
{
    /* Initialized data */

    static doublereal eps = 1e-4;

    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, d_dim1, d_offset, i__1;
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double r_sign(), sqrt();

    /* Local variables */
    static doublereal b;
    static real h__;
    static integer i__, j, k;
    static doublereal q, s;
    static real t, u;
    static doublereal v;
    static integer n0;
    static doublereal dv;
    static integer kp;
    static real tk, sg;
    static integer nw;

/*<       real x(n,*),y(n),w(n),tb(5,*),cm(*),sc(n,*),bl(n)                  >*/
/*<       double precision d(n,*),dy(*),db(*),b,s,yb,sw,dv,eps,v,q           >*/
/*<       data eps /1.d-4/                                                   >*/
    /* Parameter adjustments */
    d_dim1 = *n;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    --bl;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    tb -= 6;
    --cm;
    --dy;
    --db;

    /* Function Body */
/*<       kp=kr+1                                                            >*/
    kp = *kr + 1;
/*<       b=0.d0                                                             >*/
    b = 0.;
/*<       t=tb(2,m)                                                          >*/
    t = tb[*m * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(il .ne. 1) go to 3                                              >*/
    if (*il != 1) {
    goto L3;
    }
/*<       tk=tb(3,m)                                                         >*/
    tk = tb[*m * 5 + 3];
/*<       do 2 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       h=bl(i)                                                            >*/
    h__ = bl[i__];
/*<       if(h .gt. 0.0) go to 1                                             >*/
    if (h__ > (float)0.) {
        goto L1;
    }
/*<       sc(i,m)=0.0                                                        >*/
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 sc(i,m)=h*(x(i,j)-tk)                                              >*/
L1:
    sc[i__ + *m * sc_dim1] = h__ * (x[i__ + j * x_dim1] - tk);
/*<       b=b+w(i)*sc(i,m)                                                   >*/
    b += w[i__] * sc[i__ + *m * sc_dim1];
/*<     2 continue                                                           >*/
L2:
    ;
    }
/*<       go to 17                                                           >*/
    goto L17;
/*<     3 if(cm(2*j) .le. 0.0) go to 12                                      >*/
L3:
    if (cm[j * 2] <= (float)0.) {
    goto L12;
    }
/*<       k=tb(3,m)+.1                                                       >*/
    k = tb[*m * 5 + 3] + (float).1;
/*<       nw=0                                                               >*/
    nw = 0;
/*<       n0=nw                                                              >*/
    n0 = nw;
/*<       do 11 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       h=bl(i)                                                            >*/
    h__ = bl[i__];
/*<       if(h .gt. 0.0) go to 4                                             >*/
    if (h__ > (float)0.) {
        goto L4;
    }
/*<       sc(i,m)=0.0                                                        >*/
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 11                                                           >*/
    goto L11;
/*<     4 u=cm(int(x(i,j)+.1)+k)                                             >*/
L4:
    u = cm[(integer) (x[i__ + j * x_dim1] + (float).1) + k];
/*<       if(w(i) .le. 0.0) go to 5                                          >*/
    if (w[i__] <= (float)0.) {
        goto L5;
    }
/*<       nw=nw+1                                                            >*/
    ++nw;
/*<       if(u.eq.0.0) n0=n0+1                                               >*/
    if (u == (float)0.) {
        ++n0;
    }
/*<     5 if(t .ge. 0.0) go to 8                                             >*/
L5:
    if (t >= (float)0.) {
        goto L8;
    }
/*<       if(u .ne. 0.0) go to 6                                             >*/
    if (u != (float)0.) {
        goto L6;
    }
/*<       sc(i,m)=h                                                          >*/
    sc[i__ + *m * sc_dim1] = h__;
/*<       go to 10                                                           >*/
    goto L10;
/*<     6 sc(i,m)=0.0                                                        >*/
L6:
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 10                                                           >*/
    goto L10;
/*<     8 if(u .gt. 0.0) go to 9                                             >*/
L8:
    if (u > (float)0.) {
        goto L9;
    }
/*<       sc(i,m)=0.0                                                        >*/
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 sc(i,m)=h                                                          >*/
L9:
    sc[i__ + *m * sc_dim1] = h__;
/*<    10 b=b+w(i)*sc(i,m)                                                   >*/
L10:
    b += w[i__] * sc[i__ + *m * sc_dim1];
/*<    11 continue                                                           >*/
L11:
    ;
    }
/*<       if(n0.eq.0.or.n0.eq.nw) return                                     >*/
    if (n0 == 0 || n0 == nw) {
    return 0;
    }
/*<       go to 17                                                           >*/
    goto L17;
/*<    12 tk=tb(3,m)                                                         >*/
L12:
    tk = tb[*m * 5 + 3];
/*<       sg=sign(1.0,t)                                                     >*/
    sg = r_sign(&c_b196, &t);
/*<       do 16 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       h=bl(i)                                                            >*/
    h__ = bl[i__];
/*<       if(h .gt. 0.0) go to 13                                            >*/
    if (h__ > (float)0.) {
        goto L13;
    }
/*<       sc(i,m)=0.0                                                        >*/
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 16                                                           >*/
    goto L16;
/*<    13 u=amax1(0.0,sg*(x(i,j)-tk))                                        >*/
L13:
/* Computing MAX */
    r__1 = (float)0., r__2 = sg * (x[i__ + j * x_dim1] - tk);
    u = dmax(r__1,r__2);
/*<       if(u .gt. 0.0) go to 14                                            >*/
    if (u > (float)0.) {
        goto L14;
    }
/*<       sc(i,m)=0.0                                                        >*/
    sc[i__ + *m * sc_dim1] = (float)0.;
/*<       go to 15                                                           >*/
    goto L15;
/*<    14 sc(i,m)=h*u                                                        >*/
L14:
    sc[i__ + *m * sc_dim1] = h__ * u;
/*<    15 b=b+w(i)*sc(i,m)                                                   >*/
L15:
    b += w[i__] * sc[i__ + *m * sc_dim1];
/*<    16 continue                                                           >*/
L16:
    ;
    }
/*<    17 b=b/sw                                                             >*/
L17:
    b /= *sw;
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       v=s                                                                >*/
    v = s;
/*<       do 18 j=1,kr                                                       >*/
    i__1 = *kr;
    for (j = 1; j <= i__1; ++j) {
/*<       db(j)=0.d0                                                         >*/
    db[j] = 0.;
/*<    18 continue                                                           >*/
/* L18: */
    }
/*<       do 21 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       d(i,kp)=sc(i,m)-b                                                  >*/
    d__[i__ + kp * d_dim1] = sc[i__ + *m * sc_dim1] - b;
/*<       if(sc(i,m).le.0.0) go to 21                                        >*/
    if (sc[i__ + *m * sc_dim1] <= (float)0.) {
        goto L21;
    }
/*<       q=w(i)*sc(i,m)                                                     >*/
    q = w[i__] * sc[i__ + *m * sc_dim1];
/*<       s=s+q*(sc(i,m)-b)                                                  >*/
    s += q * (sc[i__ + *m * sc_dim1] - b);
/*<       v=v+q*(y(i)-yb)                                                    >*/
    v += q * (y[i__] - *yb);
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 20                                                           >*/
    goto L20;
/*<    19 j=j+1                                                              >*/
L19:
    ++j;
/*<    20 if((j).gt.(kr)) go to 21                                           >*/
L20:
    if (j > *kr) {
        goto L21;
    }
/*<       db(j)=db(j)+q*d(i,j)                                               >*/
    db[j] += q * d__[i__ + j * d_dim1];
/*<       go to 19                                                           >*/
    goto L19;
/*<    21 continue                                                           >*/
L21:
    ;
    }
/*<       if(s.le.0.d0) return                                               >*/
    if (s <= 0.) {
    return 0;
    }
/*<       dv=s                                                               >*/
    dv = s;
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 23                                                           >*/
    goto L23;
/*<    22 j=j+1                                                              >*/
L22:
    ++j;
/*<    23 if((j).gt.(kr)) go to 24                                           >*/
L23:
    if (j > *kr) {
    goto L24;
    }
/*<       s=s-db(j)**2                                                       >*/
/* Computing 2nd power */
    d__1 = db[j];
    s -= d__1 * d__1;
/*<       go to 22                                                           >*/
    goto L22;
/*<    24 if(s.lt.eps*dv) return                                             >*/
L24:
    if (s < eps * dv) {
    return 0;
    }
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 26                                                           >*/
    goto L26;
/*<    25 j=j+1                                                              >*/
L25:
    ++j;
/*<    26 if((j).gt.(kr)) go to 28                                           >*/
L26:
    if (j > *kr) {
    goto L28;
    }
/*<       do 27 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       d(i,kp)=d(i,kp)-db(j)*d(i,j)                                       >*/
    d__[i__ + kp * d_dim1] -= db[j] * d__[i__ + j * d_dim1];
/*<    27 continue                                                           >*/
/* L27: */
    }
/*<       go to 25                                                           >*/
    goto L25;
/*<    28 s=1.d0/dsqrt(s)                                                    >*/
L28:
    s = 1. / sqrt(s);
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 30                                                           >*/
    goto L30;
/*<    29 j=j+1                                                              >*/
L29:
    ++j;
/*<    30 if((j).gt.(kr)) go to 31                                           >*/
L30:
    if (j > *kr) {
    goto L31;
    }
/*<       v=v-db(j)*dy(j)                                                    >*/
    v -= db[j] * dy[j];
/*<       go to 29                                                           >*/
    goto L29;
/*<    31 dy(kp)=v*s                                                         >*/
L31:
    dy[kp] = v * s;
/*<       do 32 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       d(i,kp)=d(i,kp)*s                                                  >*/
    d__[i__ + kp * d_dim1] *= s;
/*<    32 continue                                                           >*/
/* L32: */
    }
/*<       kr=kp                                                              >*/
    *kr = kp;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* update_ */

/*<       subroutine pr(um,u,up,p,r)                                         >*/
/* Subroutine */ int pr_(um, u, up, p, r__)
real *um, *u, *up, *p, *r__;
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real s;

/*<       s=1.0                                                              >*/
    s = (float)1.;
/*<       if(um.gt.up) s=-1.0                                                >*/
    if (*um > *up) {
    s = (float)-1.;
    }
/*<       p=s*(2.0*up+um-3.0*u)/(up-um)**2                                   >*/
/* Computing 2nd power */
    r__1 = *up - *um;
    *p = s * (*up * (float)2. + *um - *u * (float)3.) / (r__1 * r__1);
/*<       r=s*(2.0*u-up-um)/(up-um)**3                                       >*/
/* Computing 3rd power */
    r__1 = *up - *um, r__2 = r__1;
    *r__ = s * (*u * (float)2. - *up - *um) / (r__2 * (r__1 * r__1));
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* pr_ */

/*<       function cue(x,um,u,up,p,r)                                        >*/
doublereal cue_(x, um, u, up, p, r__)
real *x, *um, *u, *up, *p, *r__;
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static real s, y;

/*<       s=1.0                                                              >*/
    s = (float)1.;
/*<       if(um.gt.up) s=-1.0                                                >*/
    if (*um > *up) {
    s = (float)-1.;
    }
/*<       y=s*x                                                              >*/
    y = s * *x;
/*<       if(y .gt. s*um) go to 1                                            >*/
    if (y > s * *um) {
    goto L1;
    }
/*<       cue=0.0                                                            >*/
    ret_val = (float)0.;
/*<       return                                                             >*/
    return ret_val;
/*<     1 if(y .lt. s*up) go to 2                                            >*/
L1:
    if (y < s * *up) {
    goto L2;
    }
/*<       cue=y-s*u                                                          >*/
    ret_val = y - s * *u;
/*<       return                                                             >*/
    return ret_val;
/*<     2 cue=p*(x-um)**2+r*(x-um)**3                                        >*/
L2:
/* Computing 2nd power */
    r__1 = *x - *um;
/* Computing 3rd power */
    r__2 = *x - *um, r__3 = r__2;
    ret_val = *p * (r__1 * r__1) + *r__ * (r__3 * (r__2 * r__2));
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* cue_ */

/*<       function varf(nk,d,a,sw,k1,k2)                                     >*/
doublereal varf_(nk, d__, a, sw, k1, k2)
integer *nk;
doublereal *d__, *a, *sw;
integer *k1, *k2;
{
    /* System generated locals */
    integer d_dim1, d_offset, i__1, i__2;
    real ret_val;

    /* Local variables */
    static integer i__, j;
    static doublereal s, t, u;

/*<       double precision d(nk,*),a(nk),sw,s,t,u                            >*/
/*<       s=0.d0                                                             >*/
    /* Parameter adjustments */
    --a;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;

    /* Function Body */
    s = 0.;
/*<       do 4 i=k1,k2                                                       >*/
    i__1 = *k2;
    for (i__ = *k1; i__ <= i__1; ++i__) {
/*<       t=0.d0                                                             >*/
    t = 0.;
/*<       do 3 j=k1,k2                                                       >*/
    i__2 = *k2;
    for (j = *k1; j <= i__2; ++j) {
/*<       if(j .gt. i) go to 1                                               >*/
        if (j > i__) {
        goto L1;
        }
/*<       u=d(j,i)                                                           >*/
        u = d__[j + i__ * d_dim1];
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 u=d(i,j)                                                           >*/
L1:
        u = d__[i__ + j * d_dim1];
/*<     2 t=t+a(j)*u                                                         >*/
L2:
        t += a[j] * u;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       s=s+a(i)*t                                                         >*/
    s += a[i__] * t;
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       varf=s/sw                                                          >*/
    ret_val = s / *sw;
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* varf_ */

/*<       subroutine exch(nk,m,k,d,a,b)                                      >*/
/* Subroutine */ int exch_(nk, m, k, d__, a, b)
integer *nk, *m, *k;
doublereal *d__;
real *a, *b;
{
    /* System generated locals */
    integer d_dim1, d_offset;

    /* Local variables */
    static integer i__, j, l;
    static real r__;
    static doublereal t;
    static integer km1, kp2;

/*<       real a(m),b(m)                                                     >*/
/*<       double precision d(nk,m),t                                         >*/
/*<       l=k+1                                                              >*/
    /* Parameter adjustments */
    --b;
    --a;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;

    /* Function Body */
    l = *k + 1;
/*<       km1=k-1                                                            >*/
    km1 = *k - 1;
/*<       kp2=k+2                                                            >*/
    kp2 = *k + 2;
/*<       r=a(k)                                                             >*/
    r__ = a[*k];
/*<       a(k)=a(l)                                                          >*/
    a[*k] = a[l];
/*<       a(l)=r                                                             >*/
    a[l] = r__;
/*<       r=b(k)                                                             >*/
    r__ = b[*k];
/*<       b(k)=b(l)                                                          >*/
    b[*k] = b[l];
/*<       b(l)=r                                                             >*/
    b[l] = r__;
/*<       do 1 j=1,2                                                         >*/
    for (j = 1; j <= 2; ++j) {
/*<       i=nk+j                                                             >*/
    i__ = *nk + j;
/*<       t=d(k,i)                                                           >*/
    t = d__[*k + i__ * d_dim1];
/*<       d(k,i)=d(l,i)                                                      >*/
    d__[*k + i__ * d_dim1] = d__[l + i__ * d_dim1];
/*<       d(l,i)=t                                                           >*/
    d__[l + i__ * d_dim1] = t;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       t=d(k,k)                                                           >*/
    t = d__[*k + *k * d_dim1];
/*<       d(k,k)=d(l,l)                                                      >*/
    d__[*k + *k * d_dim1] = d__[l + l * d_dim1];
/*<       d(l,l)=t                                                           >*/
    d__[l + l * d_dim1] = t;
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 j=j+1                                                              >*/
L2:
    ++j;
/*<     3 if((j).gt.(km1)) go to 4                                           >*/
L3:
    if (j > km1) {
    goto L4;
    }
/*<       t=d(j,k)                                                           >*/
    t = d__[j + *k * d_dim1];
/*<       d(j,k)=d(j,l)                                                      >*/
    d__[j + *k * d_dim1] = d__[j + l * d_dim1];
/*<       d(j,l)=t                                                           >*/
    d__[j + l * d_dim1] = t;
/*<       go to 2                                                            >*/
    goto L2;
/*<     4 j=kp2                                                              >*/
L4:
    j = kp2;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 j=j+1                                                              >*/
L5:
    ++j;
/*<     6 if((j).gt.(m)) go to 7                                             >*/
L6:
    if (j > *m) {
    goto L7;
    }
/*<       t=d(k,j)                                                           >*/
    t = d__[*k + j * d_dim1];
/*<       d(k,j)=d(l,j)                                                      >*/
    d__[*k + j * d_dim1] = d__[l + j * d_dim1];
/*<       d(l,j)=t                                                           >*/
    d__[l + j * d_dim1] = t;
/*<       go to 5                                                            >*/
    goto L5;
/*<     7 return                                                             >*/
L7:
    return 0;
/*<       end                                                                >*/
} /* exch_ */

/*<       function jft(m,j,tb)                                               >*/
integer jft_(m, j, tb)
integer *m, *j;
real *tb;
{
    /* System generated locals */
    integer ret_val;
    real r__1;

    /* Local variables */
    static integer k;

/*<       real tb(5,*)                                                       >*/
/*<       k=1                                                                >*/
    /* Parameter adjustments */
    tb -= 6;

    /* Function Body */
    k = 1;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 k=k+1                                                              >*/
L1:
    ++k;
/*<     2 if((k).gt.(m)) go to 4                                             >*/
L2:
    if (k > *m) {
    goto L4;
    }
/*<       if(int(abs(tb(2,k))+.1) .ne. j) go to 1                            >*/
    if ((integer) ((r__1 = tb[k * 5 + 2], dabs(r__1)) + (float).1) != *j) {
    goto L1;
    }
/*<       jft=1                                                              >*/
    ret_val = 1;
/*<       return                                                             >*/
    return ret_val;
/*<     4 jft=0                                                              >*/
L4:
    ret_val = 0;
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* jft_ */

/*<       subroutine coefpr (it,nk,az,tb,cm,xs)                              >*/
/* Subroutine */ int coefpr_(it, nk, az, tb, cm, xs)
integer *it, *nk;
real *az, *tb, *cm, *xs;
{
    static real a[6];
    static integer i1, i2, l2;
    extern /* Subroutine */ int org_();

/*<       real tb(5,*),cm(*),xs(*),a(6)                                      >*/
/*<       i2=0                                                               >*/
    /* Parameter adjustments */
    --xs;
    --cm;
    tb -= 6;

    /* Function Body */
    i2 = 0;
/*<     1 if(i2.ge.nk) go to 4                                               >*/
L1:
    if (i2 >= *nk) {
    goto L4;
    }
/*<       if(i2 .ne. 0) go to 2                                              >*/
    if (i2 != 0) {
    goto L2;
    }
/*<       i1=0                                                               >*/
    i1 = 0;
/*<       i2=min0(5,nk)                                                      >*/
    i2 = min(5,*nk);
/*<       l2=i2+1                                                            >*/
    l2 = i2 + 1;
/*<       a(1)=az                                                            >*/
    a[0] = *az;
/*<       call org(1,i2,tb,cm,xs,a(2))                                       >*/
    org_(&c__1, &i2, &tb[6], &cm[1], &xs[1], &a[1]);
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 i1=i2+1                                                            >*/
L2:
    i1 = i2 + 1;
/*<       i2=i2+6                                                            >*/
    i2 += 6;
/*<       if(i2.gt.nk) i2=nk                                                 >*/
    if (i2 > *nk) {
    i2 = *nk;
    }
/*<       l2=i2-i1+1                                                         >*/
    l2 = i2 - i1 + 1;
/*<       call org(i1,i2,tb,cm,xs,a)                                         >*/
    org_(&i1, &i2, &tb[6], &cm[1], &xs[1], a);
/*<     3 continue >*/
L3:
/*    write(it,'(/,'' bsfn:'',6(''    '',i4,''    ''))') (i,i=i1,i2)     2
874*/
/*    write(it,'('' coef:'',6g12.4)') (a(i),i=1,l2)                      2
875*/
/*<       go to 1                                                            >*/
    goto L1;
/*<     4 return                                                             >*/
L4:
    return 0;
/*<       end                                                                >*/
} /* coefpr_ */

/*<       subroutine org(m1,m2,tb,cm,xs,a)                                   >*/
/* Subroutine */ int org_(m1, m2, tb, cm, xs, a)
integer *m1, *m2;
real *tb, *cm, *xs, *a;
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer j, k, m;
    static real s;
    static integer ip;

/*<       real xs(*),tb(5,*),cm(*),a(*)                                      >*/
/*<       k=0                                                                >*/
    /* Parameter adjustments */
    --a;
    --xs;
    --cm;
    tb -= 6;

    /* Function Body */
    k = 0;
/*<       do 4 m=m1,m2                                                       >*/
    i__1 = *m2;
    for (m = *m1; m <= i__1; ++m) {
/*<       k=k+1                                                              >*/
    ++k;
/*<       if(tb(1,m) .ne. 0.0) go to 1                                       >*/
    if (tb[m * 5 + 1] != (float)0.) {
        goto L1;
    }
/*<       a(k)=0.0                                                           >*/
    a[k] = (float)0.;
/*<       go to 4                                                            >*/
    goto L4;
/*<     1 s=1.0                                                              >*/
L1:
    s = (float)1.;
/*<       ip=m                                                               >*/
    ip = m;
/*<     2 if(ip.le.0) go to 3                                                >*/
L2:
    if (ip <= 0) {
        goto L3;
    }
/*<       j=abs(tb(2,ip))+.1                                                 >*/
    j = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(cm(2*j).eq.0.0) s=s*xs(j)                                       >*/
    if (cm[j * 2] == (float)0.) {
        s *= xs[j];
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 2                                                            >*/
    goto L2;
/*<     3 a(k)=tb(1,m)/s                                                     >*/
L3:
    a[k] = tb[m * 5 + 1] / s;
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* org_ */

/*<       subroutine hulset (n,x,big,nh,xh,y)                                >*/
/* Subroutine */ int hulset_(n, x, big, nh, xh, y)
integer *n;
real *x, *big;
integer *nh;
real *xh, *y;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static real a, b;
    static integer i__, j, k;
    static real s, x1, x2, sg;

/*<       real x(n,*),y(n),xh(3,nh)                                          >*/
/*<       do 5 j=1,n                                                         >*/
    /* Parameter adjustments */
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    xh -= 4;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<       k=0                                                                >*/
    k = 0;
/*<       x1=x(j,1)                                                          >*/
    x1 = x[j + x_dim1];
/*<       x2=x(j,2)                                                          >*/
    x2 = x[j + (x_dim1 << 1)];
/*<       do 3 i=1,nh                                                        >*/
    i__2 = *nh;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       a=xh(1,i)                                                          >*/
        a = xh[i__ * 3 + 1];
/*<       b=xh(2,i)                                                          >*/
        b = xh[i__ * 3 + 2];
/*<       sg=xh(3,i)                                                         >*/
        sg = xh[i__ * 3 + 3];
/*<       if(a .lt. big) go to 1                                             >*/
        if (a < *big) {
        goto L1;
        }
/*<       s=x1-b                                                             >*/
        s = x1 - b;
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 s=x2-a*x1-b                                                        >*/
L1:
        s = x2 - a * x1 - b;
/*<     2 if(s*sg .ge. 0.0) go to 3                                          >*/
L2:
        if (s * sg >= (float)0.) {
        goto L3;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 4                                                            >*/
        goto L4;
/*<     3 continue                                                           >*/
L3:
        ;
    }
/*<     4 if(k.eq.1) y(j)=big                                                >*/
L4:
    if (k == 1) {
        y[j] = *big;
    }
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* hulset_ */

/*<       subroutine cvxhul (n,x1,x2,big,nh,xh)                              >*/
/* Subroutine */ int cvxhul_(n, x1, x2, big, nh, xh)
integer *n;
real *x1, *x2, *big;
integer *nh;
real *xh;
{
    /* Initialized data */

    static real eps = (float).001;

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    static real a, b;
    static integer i__, k;
    static real s, x, y, a0;
    static integer k0;
    static real x0, y0, am, xb, yb;
    static integer kq, lq, iq;
    static real xm, ym;

/*<       real x1(n),x2(n),xh(3,*)                                           >*/
/*<       data eps /1.e-3/                                                   >*/
    /* Parameter adjustments */
    --x2;
    --x1;
    xh -= 4;

    /* Function Body */
/*<       x0=big                                                             >*/
    x0 = *big;
/*<       y0=x0                                                              >*/
    y0 = x0;
/*<       xm=-big                                                            >*/
    xm = -(*big);
/*<       ym=xm                                                              >*/
    ym = xm;
/*<       xb=0.0                                                             >*/
    xb = (float)0.;
/*<       yb=xb                                                              >*/
    yb = xb;
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       xb=xb+x1(i)                                                        >*/
    xb += x1[i__];
/*<       yb=yb+x2(i)                                                        >*/
    yb += x2[i__];
/*<       x0=amin1(x0,x1(i))                                                 >*/
/* Computing MIN */
    r__1 = x0, r__2 = x1[i__];
    x0 = dmin(r__1,r__2);
/*<       y0=amin1(y0,x2(i))                                                 >*/
/* Computing MIN */
    r__1 = y0, r__2 = x2[i__];
    y0 = dmin(r__1,r__2);
/*<       xm=amax1(xm,x1(i))                                                 >*/
/* Computing MAX */
    r__1 = xm, r__2 = x1[i__];
    xm = dmax(r__1,r__2);
/*<       ym=amax1(ym,x2(i))                                                 >*/
/* Computing MAX */
    r__1 = ym, r__2 = x2[i__];
    ym = dmax(r__1,r__2);
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       x0=x0-eps*(xm-x0)                                                  >*/
    x0 -= eps * (xm - x0);
/*<       y0=y0-eps*(ym-y0)                                                  >*/
    y0 -= eps * (ym - y0);
/*<       xb=xb/n                                                            >*/
    xb /= *n;
/*<       yb=yb/n                                                            >*/
    yb /= *n;
/*<       nh=0                                                               >*/
    *nh = 0;
/*<       a0=0.0                                                             >*/
    a0 = (float)0.;
/*<       lq=1                                                               >*/
    lq = 1;
/*<     2 am=big                                                             >*/
L2:
    am = *big;
/*<       kq=4                                                               >*/
    kq = 4;
/*<       k=0                                                                >*/
    k = 0;
/*<       do 15 i=1,n                                                        >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       x=x1(i)                                                            >*/
    x = x1[i__];
/*<       y=x2(i)                                                            >*/
    y = x2[i__];
/*<       if(x .ne. x0) go to 5                                              >*/
    if (x != x0) {
        goto L5;
    }
/*<       if(y.eq.y0) go to 15                                               >*/
    if (y == y0) {
        goto L15;
    }
/*<       if(y .le. y0) go to 3                                              >*/
    if (y <= y0) {
        goto L3;
    }
/*<       iq=2                                                               >*/
    iq = 2;
/*<       go to 10                                                           >*/
    goto L10;
/*<     3 iq=4                                                               >*/
L3:
    iq = 4;
/*<       go to 10                                                           >*/
    goto L10;
/*<     5 if(x .le. x0) go to 8                                              >*/
L5:
    if (x <= x0) {
        goto L8;
    }
/*<       if(y .lt. y0) go to 6                                              >*/
    if (y < y0) {
        goto L6;
    }
/*<       iq=1                                                               >*/
    iq = 1;
/*<       go to 10                                                           >*/
    goto L10;
/*<     6 iq=4                                                               >*/
L6:
    iq = 4;
/*<       go to 10                                                           >*/
    goto L10;
/*<     8 if(y .le. y0) go to 9                                              >*/
L8:
    if (y <= y0) {
        goto L9;
    }
/*<       iq=2                                                               >*/
    iq = 2;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 iq=3                                                               >*/
L9:
    iq = 3;
/*<    10 if(iq.gt.kq) go to 15                                              >*/
L10:
    if (iq > kq) {
        goto L15;
    }
/*<       if(iq.lt.lq) go to 15                                              >*/
    if (iq < lq) {
        goto L15;
    }
/*<       if((iq .ne. 1) .and. (iq .ne. 3)) go to 11                         >*/
    if (iq != 1 && iq != 3) {
        goto L11;
    }
/*<       a=abs((y-y0)/(x-x0))                                               >*/
    a = (r__1 = (y - y0) / (x - x0), dabs(r__1));
/*<       go to 12                                                           >*/
    goto L12;
/*<    11 a=abs((x-x0)/(y-y0))                                               >*/
L11:
    a = (r__1 = (x - x0) / (y - y0), dabs(r__1));
/*<    12 if(iq.eq.lq.and.a.lt.a0) go to 15                                  >*/
L12:
    if (iq == lq && a < a0) {
        goto L15;
    }
/*<       if(iq .ge. kq) go to 13                                            >*/
    if (iq >= kq) {
        goto L13;
    }
/*<       kq=iq                                                              >*/
    kq = iq;
/*<       am=a                                                               >*/
    am = a;
/*<       k=i                                                                >*/
    k = i__;
/*<       go to 15                                                           >*/
    goto L15;
/*<    13 if(a .ge. am) go to 14                                             >*/
L13:
    if (a >= am) {
        goto L14;
    }
/*<       am=a                                                               >*/
    am = a;
/*<       k=i                                                                >*/
    k = i__;
/*<       go to 15                                                           >*/
    goto L15;
/*<    14 if(a .ne. am) go to 15                                             >*/
L14:
    if (a != am) {
        goto L15;
    }
/*<       if((x-x0)**2+(y-y0)**2.gt.(x1(k)-x0)**2+(x2(k)-y0)**2) k=i         >*/
/* Computing 2nd power */
    r__1 = x - x0;
/* Computing 2nd power */
    r__2 = y - y0;
/* Computing 2nd power */
    r__3 = x1[k] - x0;
/* Computing 2nd power */
    r__4 = x2[k] - y0;
    if (r__1 * r__1 + r__2 * r__2 > r__3 * r__3 + r__4 * r__4) {
        k = i__;
    }
/*<    15 continue                                                           >*/
L15:
    ;
    }
/*<       if(k .ne. 0) go to 16                                              >*/
    if (k != 0) {
    goto L16;
    }
/*<       a0=0.0                                                             >*/
    a0 = (float)0.;
/*<       lq=1                                                               >*/
    lq = 1;
/*<       go to 2                                                            >*/
    goto L2;
/*<    16 if(nh .ne. 0) go to 17                                             >*/
L16:
    if (*nh != 0) {
    goto L17;
    }
/*<       k0=k                                                               >*/
    k0 = k;
/*<       go to 20                                                           >*/
    goto L20;
/*<    17 if(x1(k) .ne. x0) go to 18                                         >*/
L17:
    if (x1[k] != x0) {
    goto L18;
    }
/*<       a=big                                                              >*/
    a = *big;
/*<       b=x0                                                               >*/
    b = x0;
/*<       s=xb-b                                                             >*/
    s = xb - b;
/*<       go to 19                                                           >*/
    goto L19;
/*<    18 a=(x2(k)-y0)/(x1(k)-x0)                                            >*/
L18:
    a = (x2[k] - y0) / (x1[k] - x0);
/*<       b=y0-a*x0                                                          >*/
    b = y0 - a * x0;
/*<       s=yb-a*xb-b                                                        >*/
    s = yb - a * xb - b;
/*<    19 xh(1,nh)=a                                                         >*/
L19:
    xh[*nh * 3 + 1] = a;
/*<       xh(2,nh)=b                                                         >*/
    xh[*nh * 3 + 2] = b;
/*<       xh(3,nh)=sign(1.0,s)                                               >*/
    xh[*nh * 3 + 3] = r_sign(&c_b196, &s);
/*<       if(k.eq.k0) go to 21                                               >*/
    if (k == k0) {
    goto L21;
    }
/*<    20 nh=nh+1                                                            >*/
L20:
    ++(*nh);
/*<       x0=x1(k)                                                           >*/
    x0 = x1[k];
/*<       y0=x2(k)                                                           >*/
    y0 = x2[k];
/*<       lq=kq                                                              >*/
    lq = kq;
/*<       a0=am                                                              >*/
    a0 = am;
/*<       go to 2                                                            >*/
    goto L2;
/*<    21 return                                                             >*/
L21:
    return 0;
/*<       end                                                                >*/
} /* cvxhul_ */

/*<       subroutine knts (l,nt,jv,jl,kv,nk,tb,cm,x,js)                      >*/
/* Subroutine */ int knts_(l, nt, jv, jl, kv, nk, tb, cm, x, js)
integer *l, *nt, *jv, *jl, *kv, *nk;
real *tb, *cm, *x;
integer *js;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    static integer j, k, m;
    static real t;
    extern integer nordc_();
    static integer l1;
    extern integer jf_();
    static integer ip;
    extern integer icf_();

/*<       integer jv(l),kv(2,jl),js(*)                                       >*/
/*<       real tb(5,nk),cm(*),x(nt,l)                                        >*/
/*<       l1=0                                                               >*/
    /* Parameter adjustments */
    --jv;
    x_dim1 = *nt;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    kv -= 3;
    tb -= 6;
    --cm;
    --js;

    /* Function Body */
    l1 = 0;
/*<       do 7 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7                             >*/
    if (icf_(&m, &tb[6], &cm[1], jl, &kv[3], &js[1]) == 0) {
        goto L7;
    }
/*<       if(nordc(1,m,tb,cm).ne.l) go to 7                                  >*/
    if (nordc_(&c__1, &m, &tb[6], &cm[1]) != *l) {
        goto L7;
    }
/*<       k=0                                                                >*/
    k = 0;
/*<       do 1 j=1,l                                                         >*/
    i__2 = *l;
    for (j = 1; j <= i__2; ++j) {
/*<       if(jf(m,jv(j),tb).eq.1) go to 1                                    >*/
        if (jf_(&m, &jv[j], &tb[6]) == 1) {
        goto L1;
        }
/*<       k=1                                                                >*/
        k = 1;
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 continue                                                           >*/
L1:
        ;
    }
/*<     2 if(k.eq.1) go to 7                                                 >*/
L2:
    if (k == 1) {
        goto L7;
    }
/*<       ip=m                                                               >*/
    ip = m;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<     3 if(ip.le.0) go to 7                                                >*/
L3:
    if (ip <= 0) {
        goto L7;
    }
/*<       t=tb(2,ip)                                                         >*/
    t = tb[ip * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(cm(2*j) .eq. 0.0) go to 4                                       >*/
    if (cm[j * 2] == (float)0.) {
        goto L4;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     4 k=1                                                                >*/
L4:
    k = 1;
/*<     5 if(jv(k).eq.j) go to 6                                             >*/
L5:
    if (jv[k] == j) {
        goto L6;
    }
/*<       k=k+1                                                              >*/
    ++k;
/*<       go to 5                                                            >*/
    goto L5;
/*<     6 x(l1,k)=tb(3,ip)                                                   >*/
L6:
    x[l1 + k * x_dim1] = tb[ip * 5 + 3];
/*<       x(l1,l+k)=sign(1.0,t)                                              >*/
    x[l1 + (*l + k) * x_dim1] = r_sign(&c_b196, &t);
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     7 continue                                                           >*/
L7:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* knts_ */

/*<       subroutine sclato (n,p,x,xm,xs,cm,z)                               >*/
/* Subroutine */ int sclato_(n, p, x, xm, xs, cm, z__)
integer *n, *p;
real *x, *xm, *xs, *cm, *z__;
{
    /* System generated locals */
    integer x_dim1, x_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, j1;
    extern /* Subroutine */ int stcmrs_(), stfmrs_();

/*<       integer p                                                          >*/
/*<       real x(n,p),xm(p),xs(p),cm(*),z(n,p)                               >*/
/*<       do 4 j=1,p                                                         >*/
    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --xs;
    --xm;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --cm;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       j1=cm(2*j)+.1                                                      >*/
    j1 = cm[j * 2] + (float).1;
/*<       if(j1 .ne. 0) go to 2                                              >*/
    if (j1 != 0) {
        goto L2;
    }
/*<       if(xs(j).le.0.0) go to 4                                           >*/
    if (xs[j] <= (float)0.) {
        goto L4;
    }
/*<       do 1 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       z(i,j)=xs(j)*x(i,j)+xm(j)                                          >*/
        z__[i__ + j * z_dim1] = xs[j] * x[i__ + j * x_dim1] + xm[j];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<     2 j1=j1-1                                                            >*/
L2:
    --j1;
/*<       do 3 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       l=x(i,j)+.1                                                        >*/
        l = x[i__ + j * x_dim1] + (float).1;
/*<       z(i,j)=cm(l+j1)                                                    >*/
        z__[i__ + j * z_dim1] = cm[l + j1];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<       call stfmrs(0)                                                     >*/
    stfmrs_(&c__0);
/*<       call stcmrs(0)                                                     >*/
    stcmrs_(&c__0);
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* sclato_ */

/*<       function ncat(kp)                                                  >*/
integer ncat_(kp)
integer *kp;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer ll;

/*<       integer kp(5,*)                                                    >*/
/*<       ncat=0                                                             >*/
    /* Parameter adjustments */
    kp -= 6;

    /* Function Body */
    ret_val = 0;
/*<       ll=1                                                               >*/
    ll = 1;
/*<     1 if(kp(1,ll).lt.0) go to 2                                          >*/
L1:
    if (kp[ll * 5 + 1] < 0) {
    goto L2;
    }
/*<       if(kp(1,ll).gt.0.and.kp(3,ll).le.0) ncat=ncat+1                    >*/
    if (kp[ll * 5 + 1] > 0 && kp[ll * 5 + 3] <= 0) {
    ++ret_val;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* ncat_ */

/*<       subroutine catv (jl,kp,kv,nv,jv)                                   >*/
/* Subroutine */ int catv_(jl, kp, kv, nv, jv)
integer *jl, *kp, *kv, *nv, *jv;
{
    /* System generated locals */
    integer jv_dim1, jv_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, k1, l1, ig, jg, ll;

/*<       integer kp(5,*),kv(2,*),jv(jl,*)                                   >*/
/*<       nv=0                                                               >*/
    /* Parameter adjustments */
    jv_dim1 = *jl;
    jv_offset = jv_dim1 + 1;
    jv -= jv_offset;
    kp -= 6;
    kv -= 3;

    /* Function Body */
    *nv = 0;
/*<       ll=1                                                               >*/
    ll = 1;
/*<     1 if(kp(1,ll).lt.0) go to 20                                         >*/
L1:
    if (kp[ll * 5 + 1] < 0) {
    goto L20;
    }
/*<       if(kp(3,ll) .le. 0) go to 2                                        >*/
    if (kp[ll * 5 + 3] <= 0) {
    goto L2;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 if(kp(1,ll) .eq. jl) go to 3                                       >*/
L2:
    if (kp[ll * 5 + 1] == *jl) {
    goto L3;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     3 jg=0                                                               >*/
L3:
    jg = 0;
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 j=j+1                                                              >*/
L4:
    ++j;
/*<     5 if((j).gt.(nv)) go to 9                                            >*/
L5:
    if (j > *nv) {
    goto L9;
    }
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 6 i=1,jl                                                        >*/
    i__1 = *jl;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(jv(i,j).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 6                    >*/
    if (jv[i__ + j * jv_dim1] == (i__2 = kv[(kp[ll * 5 + 2] + i__ - 1 << 
        1) + 1], abs(i__2))) {
        goto L6;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<       go to 7                                                            >*/
    goto L7;
/*<     6 continue                                                           >*/
L6:
    ;
    }
/*<     7 if(ig .ne. 0) go to 4                                              >*/
L7:
    if (ig != 0) {
    goto L4;
    }
/*<       jg=1                                                               >*/
    jg = 1;
/*<     9 if(jg .ne. 1) go to 10                                             >*/
L9:
    if (jg != 1) {
    goto L10;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<    10 l1=ll+1                                                            >*/
L10:
    l1 = ll + 1;
/*<       ig=0                                                               >*/
    ig = 0;
/*<    11 if(kp(1,l1).lt.0) go to 17                                         >*/
L11:
    if (kp[l1 * 5 + 1] < 0) {
    goto L17;
    }
/*<       k1=kp(1,l1)                                                        >*/
    k1 = kp[l1 * 5 + 1];
/*<       if((k1 .ne. jl) .and. (kp(3,l1) .le. 0)) go to 12                  >*/
    if (k1 != *jl && kp[l1 * 5 + 3] <= 0) {
    goto L12;
    }
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       go to 11                                                           >*/
    goto L11;
/*<    12 do 15 i=1,k1                                                       >*/
L12:
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       k=iabs(kv(1,kp(2,l1)+i-1))                                         >*/
    k = (i__2 = kv[(kp[l1 * 5 + 2] + i__ - 1 << 1) + 1], abs(i__2));
/*<       do 13 j=1,jl                                                       >*/
    i__2 = *jl;
    for (j = 1; j <= i__2; ++j) {
/*<       l=iabs(kv(1,kp(2,ll)+j-1))                                         >*/
        l = (i__3 = kv[(kp[ll * 5 + 2] + j - 1 << 1) + 1], abs(i__3));
/*<       if(l .ne. k) go to 13                                              >*/
        if (l != k) {
        goto L13;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 14                                                           >*/
        goto L14;
/*<    13 continue                                                           >*/
L13:
        ;
    }
/*<    14 if(ig.eq.1) go to 16                                               >*/
L14:
    if (ig == 1) {
        goto L16;
    }
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<    16 if(ig.eq.1) go to 17                                               >*/
L16:
    if (ig == 1) {
    goto L17;
    }
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       go to 11                                                           >*/
    goto L11;
/*<    17 if(ig .ne. 1) go to 18                                             >*/
L17:
    if (ig != 1) {
    goto L18;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<    18 nv=nv+1                                                            >*/
L18:
    ++(*nv);
/*<       do 19 i=1,jl                                                       >*/
    i__1 = *jl;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       jv(i,nv)=iabs(kv(1,kp(2,ll)+i-1))                                  >*/
    jv[i__ + *nv * jv_dim1] = (i__2 = kv[(kp[ll * 5 + 2] + i__ - 1 << 1) 
        + 1], abs(i__2));
/*<    19 continue                                                           >*/
/* L19: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<    20 return                                                             >*/
L20:
    return 0;
/*<       end                                                                >*/
} /* catv_ */

/*<       function cvlv (m,jl,jv,lv,nk,kp,kv,tb,cm,tc)                       >*/
doublereal cvlv_(m, jl, jv, lv, nk, kp, kv, tb, cm, tc)
integer *m, *jl, *jv, *lv, *nk, *kp, *kv;
real *tb, *cm, *tc;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern doublereal cvll_(), cvlq_();

/*<       integer jv(jl),lv(jl),kp(5,*),kv(2,*)                              >*/
/*<       real tb(5,nk),cm(*),tc(*)                                          >*/
/*<       if(m .ne. 1) go to 1                                               >*/
    /* Parameter adjustments */
    --lv;
    --jv;
    tb -= 6;
    kp -= 6;
    kv -= 3;
    --cm;
    --tc;

    /* Function Body */
    if (*m != 1) {
    goto L1;
    }
/*<       cvlv=cvll(jl,jv,lv,nk,tb,cm)                                       >*/
    ret_val = cvll_(jl, &jv[1], &lv[1], nk, &tb[6], &cm[1]);
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 cvlv=cvlq(jl,jv,lv,kp,kv,cm,tc)                                    >*/
L1:
    ret_val = cvlq_(jl, &jv[1], &lv[1], &kp[6], &kv[3], &cm[1], &tc[1]);
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* cvlv_ */

/*<       function cvll (jl,jv,lv,nk,tb,cm)                                  >*/
doublereal cvll_(jl, jv, lv, nk, tb, cm)
integer *jl, *jv, *lv, *nk;
real *tb, *cm;
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    extern integer nord_();
    static integer i__, j, m;
    static real t, u;
    extern integer nordc_();
    static integer ig, ip, iv[2];
    static real phi;
    extern /* Subroutine */ int jfv_();

/*<       integer jv(jl),lv(jl),iv(2)                                        >*/
/*<       real tb(5,nk),cm(*)                                                >*/
/*<       cvll=0.0                                                           >*/
    /* Parameter adjustments */
    --lv;
    --jv;
    tb -= 6;
    --cm;

    /* Function Body */
    ret_val = (float)0.;
/*<       if(jl.gt.2) return                                                 >*/
    if (*jl > 2) {
    return ret_val;
    }
/*<       do 11 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 11                                        >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L11;
    }
/*<       if(nordc(1,m,tb,cm).gt.0) go to 11                                 >*/
    if (nordc_(&c__1, &m, &tb[6], &cm[1]) > 0) {
        goto L11;
    }
/*<       if(nord(m,tb).ne.jl) go to 11                                      >*/
    if (nord_(&m, &tb[6]) != *jl) {
        goto L11;
    }
/*<       call jfv(m,tb,iv)                                                  >*/
    jfv_(&m, &tb[6], iv);
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 1 i=1,jl                                                        >*/
    i__2 = *jl;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(jv(i).eq.iv(i)) go to 1                                         >*/
        if (jv[i__] == iv[i__ - 1]) {
        goto L1;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 continue                                                           >*/
L1:
        ;
    }
/*<     2 if(ig.eq.1) go to 11                                               >*/
L2:
    if (ig == 1) {
        goto L11;
    }
/*<       phi=1.0                                                            >*/
    phi = (float)1.;
/*<       ip=m                                                               >*/
    ip = m;
/*<     3 if(ip.le.0) go to 10                                               >*/
L3:
    if (ip <= 0) {
        goto L10;
    }
/*<       t=tb(2,ip)                                                         >*/
    t = tb[ip * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 i=i+1                                                              >*/
L4:
    ++i__;
/*<     5 if((i).gt.(jl)) go to 6                                            >*/
L5:
    if (i__ > *jl) {
        goto L6;
    }
/*<       if(jv(i).eq.j) go to 6                                             >*/
    if (jv[i__] == j) {
        goto L6;
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<     6 u=cm(lv(i)+int(tb(3,ip)+.1))                                       >*/
L6:
    u = cm[lv[i__] + (integer) (tb[ip * 5 + 3] + (float).1)];
/*<       if(t .ge. 0.0) go to 8                                             >*/
    if (t >= (float)0.) {
        goto L8;
    }
/*<       if(u .ne. 0.0) go to 7                                             >*/
    if (u != (float)0.) {
        goto L7;
    }
/*<       u=1.0                                                              >*/
    u = (float)1.;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 u=0.0                                                              >*/
L7:
    u = (float)0.;
/*<     8 if(u .ne. 0.0) go to 9                                             >*/
L8:
    if (u != (float)0.) {
        goto L9;
    }
/*<       phi=0.0                                                            >*/
    phi = (float)0.;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 ip=tb(4,ip)+.1                                                     >*/
L9:
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 3                                                            >*/
    goto L3;
/*<    10 if(phi.gt.0.0) cvll=cvll+tb(1,m)                                   >*/
L10:
    if (phi > (float)0.) {
        ret_val += tb[m * 5 + 1];
    }
/*<    11 continue                                                           >*/
L11:
    ;
    }
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* cvll_ */

/*<       function cvlq (jl,jv,lv,kp,kv,cm,tc)                               >*/
doublereal cvlq_(jl, jv, lv, kp, kv, cm, tc)
integer *jl, *jv, *lv, *kp, *kv;
real *cm, *tc;
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static integer i__, j, k, ig, ll, jt, kt;

/*<       integer jv(jl),lv(jl),kp(5,*),kv(2,*)                              >*/
/*<       real cm(*),tc(*)                                                   >*/
/*<       ll=1                                                               >*/
    /* Parameter adjustments */
    --lv;
    --jv;
    kp -= 6;
    kv -= 3;
    --cm;
    --tc;

    /* Function Body */
    ll = 1;
/*<       cvlq=0.0                                                           >*/
    ret_val = (float)0.;
/*<     1 if(kp(1,ll).lt.0) go to 12                                         >*/
L1:
    if (kp[ll * 5 + 1] < 0) {
    goto L12;
    }
/*<       if(kp(3,ll) .le. 0) go to 2                                        >*/
    if (kp[ll * 5 + 3] <= 0) {
    goto L2;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 if(kp(1,ll) .eq. jl) go to 3                                       >*/
L2:
    if (kp[ll * 5 + 1] == *jl) {
    goto L3;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     3 ig=0                                                               >*/
L3:
    ig = 0;
/*<       do 4 i=1,jl                                                        >*/
    i__1 = *jl;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(jv(i).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 4                      >*/
    if (jv[i__] == (i__2 = kv[(kp[ll * 5 + 2] + i__ - 1 << 1) + 1], abs(
        i__2))) {
        goto L4;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<     5 if(ig .ne. 1) go to 6                                              >*/
L5:
    if (ig != 1) {
    goto L6;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     6 kt=1                                                               >*/
L6:
    kt = 1;
/*<       do 9 j=1,jl                                                        >*/
    i__1 = *jl;
    for (j = 1; j <= i__1; ++j) {
/*<       k=kp(2,ll)+j-1                                                     >*/
    k = kp[ll * 5 + 2] + j - 1;
/*<       jt=cm(lv(j)+kv(2,k))+.1                                            >*/
    jt = cm[lv[j] + kv[(k << 1) + 2]] + (float).1;
/*<       if(kv(1,k) .ge. 0) go to 8                                         >*/
    if (kv[(k << 1) + 1] >= 0) {
        goto L8;
    }
/*<       if(jt .ne. 0) go to 7                                              >*/
    if (jt != 0) {
        goto L7;
    }
/*<       jt=1                                                               >*/
    jt = 1;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 jt=0                                                               >*/
L7:
    jt = 0;
/*<     8 if(jt .ne. 0) go to 9                                              >*/
L8:
    if (jt != 0) {
        goto L9;
    }
/*<       kt=0                                                               >*/
    kt = 0;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 continue                                                           >*/
L9:
    ;
    }
/*<    10 if(kt .ne. 1) go to 11                                             >*/
L10:
    if (kt != 1) {
    goto L11;
    }
/*<       cvlq=cvlq+tc(-kp(3,ll))                                            >*/
    ret_val += tc[-kp[ll * 5 + 3]];
/*<    11 ll=ll+1                                                            >*/
L11:
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<    12 return                                                             >*/
L12:
    return ret_val;
/*<       end                                                                >*/
} /* cvlq_ */

/*<       subroutine ccoll (nk,tb,cm,kp,kv,lp,lv,jv)                         >*/
/* Subroutine */ int ccoll_(nk, tb, cm, kp, kv, lp, lv, jv)
integer *nk;
real *tb, *cm;
integer *kp, *kv, *lp, *lv, *jv;
{
    extern /* Subroutine */ int collc_(), collf_();
    static integer l1, l2, li, ll;
    extern /* Subroutine */ int purcat_();

/*<       integer kp(5,*),kv(2,*),lp(3,*),lv(*),jv(*)                        >*/
/*<       real tb(5,*),cm(*)                                                 >*/
/*<       call collc(nk,tb,cm,kp,kv,jv)                                      >*/
    /* Parameter adjustments */
    --jv;
    --lv;
    lp -= 4;
    kv -= 3;
    kp -= 6;
    --cm;
    tb -= 6;

    /* Function Body */
    collc_(nk, &tb[6], &cm[1], &kp[6], &kv[3], &jv[1]);
/*<       call purcat(nk,tb,cm,kp,kv,li,jv)                                  >*/
    purcat_(nk, &tb[6], &cm[1], &kp[6], &kv[3], &li, &jv[1]);
/*<       ll=li+1                                                            >*/
    ll = li + 1;
/*<       l1=1                                                               >*/
    l1 = 1;
/*<       l2=l1                                                              >*/
    l2 = l1;
/*<     1 if(kp(1,ll).lt.0) go to 2                                          >*/
L1:
    if (kp[ll * 5 + 1] < 0) {
    goto L2;
    }
/*<       kp(4,ll)=l1                                                        >*/
    kp[ll * 5 + 4] = l1;
/*<       call collf(nk,tb,cm,kp(1,ll),kv(1,kp(2,ll)),l1,l2,lp,lv,jv)        >*/
    collf_(nk, &tb[6], &cm[1], &kp[ll * 5 + 1], &kv[(kp[ll * 5 + 2] << 1) + 1]
        , &l1, &l2, &lp[4], &lv[1], &jv[1]);
/*<       kp(3,ll)=l1-kp(4,ll)                                               >*/
    kp[ll * 5 + 3] = l1 - kp[ll * 5 + 4];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 lp(1,l1)=0                                                         >*/
L2:
    lp[l1 * 3 + 1] = 0;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* ccoll_ */

/*<       subroutine collc (nk,tb,cm,kp,kv,jv)                               >*/
/* Subroutine */ int collc_(nk, tb, cm, kp, kv, jv)
integer *nk;
real *tb, *cm;
integer *kp, *kv, *jv;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    extern /* Subroutine */ int jfvc_();
    static integer i__, j, k, m;
    static real z__;
    extern integer nordc_();
    static integer l1, l2, m1, m2, l10, mc, jg, ig, jj, nc, kk, jk, mt, nv, 
        l1m1;

/*<       integer kp(5,*),kv(2,*),jv(*)                                      >*/
/*<       real tb(5,nk),cm(*)                                                >*/
/*<       kp(1,1)=0                                                          >*/
    /* Parameter adjustments */
    tb -= 6;
    --cm;
    kp -= 6;
    kv -= 3;
    --jv;

    /* Function Body */
    kp[6] = 0;
/*<       kp(2,1)=1                                                          >*/
    kp[7] = 1;
/*<       l1=2                                                               >*/
    l1 = 2;
/*<       l2=1                                                               >*/
    l2 = 1;
/*<       mc=0                                                               >*/
    mc = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).ne.0.0) mc=max0(mc,nordc(2,m,tb,cm))                    >*/
    if (tb[m * 5 + 1] != (float)0.) {
/* Computing MAX */
        i__2 = mc, i__3 = nordc_(&c__2, &m, &tb[6], &cm[1]);
        mc = max(i__2,i__3);
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       mt=1                                                               >*/
    mt = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 mt=mt+1                                                            >*/
L2:
    ++mt;
/*<     3 if((mt).gt.(mc)) go to 18                                          >*/
L3:
    if (mt > mc) {
    goto L18;
    }
/*<       l10=l1                                                             >*/
    l10 = l1;
/*<       do 17 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.mt) go to 17              >*/
    if (tb[m * 5 + 1] == (float)0. || nordc_(&c__2, &m, &tb[6], &cm[1]) !=
         mt) {
        goto L17;
    }
/*<       call jfvc(2,m,tb,cm,nv,jv,jv(mt+1))                                >*/
    jfvc_(&c__2, &m, &tb[6], &cm[1], &nv, &jv[1], &jv[mt + 1]);
/*<       jg=0                                                               >*/
    jg = 0;
/*<       l1m1=l1-1                                                          >*/
    l1m1 = l1 - 1;
/*<       i=l10                                                              >*/
    i__ = l10;
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 i=i+1                                                              >*/
L4:
    ++i__;
/*<     5 if((i).gt.(l1m1)) go to 15                                         >*/
L5:
    if (i__ > l1m1) {
        goto L15;
    }
/*<       k=kp(2,i)-1                                                        >*/
    k = kp[i__ * 5 + 2] - 1;
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 6 j=1,mt                                                        >*/
    i__2 = mt;
    for (j = 1; j <= i__2; ++j) {
/*<       if(iabs(jv(j)).eq.iabs(kv(1,k+j))) go to 6                         >*/
        if ((i__3 = jv[j], abs(i__3)) == (i__4 = kv[(k + j << 1) + 1], 
            abs(i__4))) {
        goto L6;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 7                                                            >*/
        goto L7;
/*<     6 continue                                                           >*/
L6:
        ;
    }
/*<     7 if(ig .ne. 0) go to 13                                             >*/
L7:
    if (ig != 0) {
        goto L13;
    }
/*<       do 12 j=1,mt                                                       >*/
    i__2 = mt;
    for (j = 1; j <= i__2; ++j) {
/*<       m1=kv(2,k+j)                                                       >*/
        m1 = kv[(k + j << 1) + 2];
/*<       m2=jv(mt+j)                                                        >*/
        m2 = jv[mt + j];
/*<       jj=iabs(jv(j))                                                     >*/
        jj = (i__3 = jv[j], abs(i__3));
/*<       nc=int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1                           >*/
        nc = (integer) (cm[(jj << 1) + 1] + (float).1) - (integer) (cm[jj 
            * 2] + (float).1) + 1;
/*<       kk=jv(j)*kv(1,k+j)                                                 >*/
        kk = jv[j] * kv[(k + j << 1) + 1];
/*<       do 10 jk=1,nc                                                      >*/
        i__3 = nc;
        for (jk = 1; jk <= i__3; ++jk) {
/*<       z=cm(jk+m2)                                                        >*/
        z__ = cm[jk + m2];
/*<       if(kk .ge. 0) go to 9                                              >*/
        if (kk >= 0) {
            goto L9;
        }
/*<       if(z .ne. 0.0) go to 8                                             >*/
        if (z__ != (float)0.) {
            goto L8;
        }
/*<       z=1.0                                                              >*/
        z__ = (float)1.;
/*<       go to 9                                                            >*/
        goto L9;
/*<     8 z=0.0                                                              >*/
L8:
        z__ = (float)0.;
/*<     9 if(cm(jk+m1).eq.z) go to 10                                        >*/
L9:
        if (cm[jk + m1] == z__) {
            goto L10;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 11                                                           >*/
        goto L11;
/*<    10 continue                                                           >*/
L10:
        ;
        }
/*<    11 if(ig.eq.1) go to 13                                               >*/
L11:
        if (ig == 1) {
        goto L13;
        }
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<    13 if(ig .ne. 0) go to 4                                              >*/
L13:
    if (ig != 0) {
        goto L4;
    }
/*<       jg=1                                                               >*/
    jg = 1;
/*<    15 if(jg .ne. 0) go to 17                                             >*/
L15:
    if (jg != 0) {
        goto L17;
    }
/*<       kp(1,l1)=mt                                                        >*/
    kp[l1 * 5 + 1] = mt;
/*<       kp(2,l1)=l2                                                        >*/
    kp[l1 * 5 + 2] = l2;
/*<       k=l2-1                                                             >*/
    k = l2 - 1;
/*<       do 16 i=1,mt                                                       >*/
    i__2 = mt;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       kv(1,i+k)=jv(i)                                                    >*/
        kv[(i__ + k << 1) + 1] = jv[i__];
/*<       kv(2,i+k)=jv(i+mt)                                                 >*/
        kv[(i__ + k << 1) + 2] = jv[i__ + mt];
/*<    16 continue                                                           >*/
/* L16: */
    }
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       l2=l2+mt                                                           >*/
    l2 += mt;
/*<    17 continue                                                           >*/
L17:
    ;
    }
/*<       go to 2                                                            >*/
    goto L2;
/*<    18 kp(1,l1)=-1                                                        >*/
L18:
    kp[l1 * 5 + 1] = -1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* collc_ */

/*<       function nordc (l,m,tb,cm)                                         >*/
integer nordc_(l, m, tb, cm)
integer *l, *m;
real *tb, *cm;
{
    /* System generated locals */
    integer ret_val;
    real r__1;

    /* Local variables */
    static integer j, ip;

/*<       real tb(5,*),cm(*)                                                 >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    --cm;
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       nordc=0                                                            >*/
    ret_val = 0;
/*<     1 if(ip.le.0) go to 4                                                >*/
L1:
    if (ip <= 0) {
    goto L4;
    }
/*<       j=abs(tb(2,ip))+.1                                                 >*/
    j = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(l .ne. 1) go to 2                                               >*/
    if (*l != 1) {
    goto L2;
    }
/*<       if(cm(2*j).eq.0.0) nordc=nordc+1                                   >*/
    if (cm[j * 2] == (float)0.) {
    ++ret_val;
    }
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 if(cm(2*j).gt.0.0) nordc=nordc+1                                   >*/
L2:
    if (cm[j * 2] > (float)0.) {
    ++ret_val;
    }
/*<     3 ip=tb(4,ip)+.1                                                     >*/
L3:
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     4 return                                                             >*/
L4:
    return ret_val;
/*<       end                                                                >*/
} /* nordc_ */

/*<       subroutine jfvc (l,m,tb,cm,nv,jv,jp)                               >*/
/* Subroutine */ int jfvc_(l, m, tb, cm, nv, jv, jp)
integer *l, *m;
real *tb, *cm;
integer *nv, *jv, *jp;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k, ll, ip;

/*<       integer jv(*),jp(*)                                                >*/
/*<       real tb(5,*),cm(*)                                                 >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    --jp;
    --jv;
    --cm;
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       nv=0                                                               >*/
    *nv = 0;
/*<     1 if(ip.le.0) go to 5                                                >*/
L1:
    if (ip <= 0) {
    goto L5;
    }
/*<       j=abs(tb(2,ip))+.1                                                 >*/
    j = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(l .ne. 1) go to 3                                               >*/
    if (*l != 1) {
    goto L3;
    }
/*<       if(cm(2*j) .le. 0.0) go to 4                                       >*/
    if (cm[j * 2] <= (float)0.) {
    goto L4;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     3 if(cm(2*j) .ne. 0.0) go to 4                                       >*/
L3:
    if (cm[j * 2] != (float)0.) {
    goto L4;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     4 nv=nv+1                                                            >*/
L4:
    ++(*nv);
/*<       jv(nv)=j                                                           >*/
    jv[*nv] = j;
/*<       if(l.ne.1.and.tb(2,ip).lt.0.0) jv(nv)=-j                           >*/
    if (*l != 1 && tb[ip * 5 + 2] < (float)0.) {
    jv[*nv] = -j;
    }
/*<       if(l.ne.1) jp(nv)=tb(3,ip)+.1                                      >*/
    if (*l != 1) {
    jp[*nv] = tb[ip * 5 + 3] + (float).1;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     5 if(nv.le.1) return                                                 >*/
L5:
    if (*nv <= 1) {
    return 0;
    }
/*<       j=nv-1                                                             >*/
    j = *nv - 1;
/*<     6 ll=0                                                               >*/
L6:
    ll = 0;
/*<       do 7 i=1,j                                                         >*/
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(iabs(jv(i)) .le. iabs(jv(i+1))) go to 7                         >*/
    if ((i__2 = jv[i__], abs(i__2)) <= (i__3 = jv[i__ + 1], abs(i__3))) {
        goto L7;
    }
/*<       ll=1                                                               >*/
    ll = 1;
/*<       k=jv(i)                                                            >*/
    k = jv[i__];
/*<       jv(i)=jv(i+1)                                                      >*/
    jv[i__] = jv[i__ + 1];
/*<       jv(i+1)=k                                                          >*/
    jv[i__ + 1] = k;
/*<       if(l .eq. 1) go to 7                                               >*/
    if (*l == 1) {
        goto L7;
    }
/*<       k=jp(i)                                                            >*/
    k = jp[i__];
/*<       jp(i)=jp(i+1)                                                      >*/
    jp[i__] = jp[i__ + 1];
/*<       jp(i+1)=k                                                          >*/
    jp[i__ + 1] = k;
/*<     7 continue                                                           >*/
L7:
    ;
    }
/*<       if(ll.eq.0) go to 8                                                >*/
    if (ll == 0) {
    goto L8;
    }
/*<       go to 6                                                            >*/
    goto L6;
/*<     8 return                                                             >*/
L8:
    return 0;
/*<       end                                                                >*/
} /* jfvc_ */

/*<       subroutine purcat (nk,tb,cm,kp,kv,li,jv)                           >*/
/* Subroutine */ int purcat_(nk, tb, cm, kp, kv, li, jv)
integer *nk;
real *tb, *cm;
integer *kp, *kv, *li, *jv;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern integer nord_();
    static integer i__, j, m, jl, ll, lm;
    extern integer icf_();
    static integer ifg, jfg;

/*<       integer kp(5,*),kv(2,*),jv(*)                                      >*/
/*<       real tb(5,nk),cm(*)                                                >*/
/*<       lm=1                                                               >*/
    /* Parameter adjustments */
    tb -= 6;
    --cm;
    kp -= 6;
    kv -= 3;
    --jv;

    /* Function Body */
    lm = 1;
/*<     1 if(kp(1,lm).lt.0) go to 2                                          >*/
L1:
    if (kp[lm * 5 + 1] < 0) {
    goto L2;
    }
/*<       lm=lm+1                                                            >*/
    ++lm;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 ll=1                                                               >*/
L2:
    ll = 1;
/*<       li=0                                                               >*/
    *li = 0;
/*<     3 if(kp(1,ll).lt.0) go to 20                                         >*/
L3:
    if (kp[ll * 5 + 1] < 0) {
    goto L20;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       if(jl .gt. 0) go to 4                                              >*/
    if (jl > 0) {
    goto L4;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 3                                                            >*/
    goto L3;
/*<     4 ifg=0                                                              >*/
L4:
    ifg = 0;
/*<       jfg=ifg                                                            >*/
    jfg = ifg;
/*<       do 6 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(icf(m,tb,cm,jl,kv(1,kp(2,ll)),jv).eq.0) go to 6                 >*/
    if (icf_(&m, &tb[6], &cm[1], &jl, &kv[(kp[ll * 5 + 2] << 1) + 1], &jv[
        1]) == 0) {
        goto L6;
    }
/*<       if(nord(m,tb) .ne. jl) go to 5                                     >*/
    if (nord_(&m, &tb[6]) != jl) {
        goto L5;
    }
/*<       ifg=1                                                              >*/
    ifg = 1;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 jfg=1                                                              >*/
L5:
    jfg = 1;
/*<     6 continue                                                           >*/
L6:
    ;
    }
/*<       if(ifg .ne. 0) go to 9                                             >*/
    if (ifg != 0) {
    goto L9;
    }
/*<       if(jfg .ne. 0) go to 8                                             >*/
    if (jfg != 0) {
    goto L8;
    }
/*    write(6,7)                                                         3
390*/
/*<     7 format (' bug in purcat - term not found.')                        >*/
/* L7: */
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     8 ll=ll+1                                                            >*/
L8:
    ++ll;
/*<       go to 3                                                            >*/
    goto L3;
/*<     9 li=li+1                                                            >*/
L9:
    ++(*li);
/*<       j=lm                                                               >*/
    j = lm;
/*<       go to 11                                                           >*/
    goto L11;
/*<    10 j=j+(-1)                                                           >*/
L10:
    --j;
/*<    11 if((-1)*((j)-(li)).gt.0) go to 13                                  >*/
L11:
    if (-(j - *li) > 0) {
    goto L13;
    }
/*<       do 12 i=1,5                                                        >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       kp(i,j+1)=kp(i,j)                                                  >*/
    kp[i__ + (j + 1) * 5] = kp[i__ + j * 5];
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<       go to 10                                                           >*/
    goto L10;
/*<    13 lm=lm+1                                                            >*/
L13:
    ++lm;
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       do 14 i=1,5                                                        >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       kp(i,li)=kp(i,ll)                                                  >*/
    kp[i__ + *li * 5] = kp[i__ + ll * 5];
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       kp(3,li)=0                                                         >*/
    kp[*li * 5 + 3] = 0;
/*<       kp(4,li)=1                                                         >*/
    kp[*li * 5 + 4] = 1;
/*<       kp(5,li)=0                                                         >*/
    kp[*li * 5 + 5] = 0;
/*<       if(jfg .ne. 1) go to 15                                            >*/
    if (jfg != 1) {
    goto L15;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 3                                                            >*/
    goto L3;
/*<    15 j=ll+1                                                             >*/
L15:
    j = ll + 1;
/*<       go to 17                                                           >*/
    goto L17;
/*<    16 j=j+1                                                              >*/
L16:
    ++j;
/*<    17 if((j).gt.(lm)) go to 19                                           >*/
L17:
    if (j > lm) {
    goto L19;
    }
/*<       do 18 i=1,5                                                        >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       kp(i,j-1)=kp(i,j)                                                  >*/
    kp[i__ + (j - 1) * 5] = kp[i__ + j * 5];
/*<    18 continue                                                           >*/
/* L18: */
    }
/*<       go to 16                                                           >*/
    goto L16;
/*<    19 lm=lm-1                                                            >*/
L19:
    --lm;
/*<       go to 3                                                            >*/
    goto L3;
/*<    20 return                                                             >*/
L20:
    return 0;
/*<       end                                                                >*/
} /* purcat_ */

/*<       subroutine collf (nk,tb,cm,jl,kv,l1,l2,lp,lv,jv)                   >*/
/* Subroutine */ int collf_(nk, tb, cm, jl, kv, l1, l2, lp, lv, jv)
integer *nk;
real *tb, *cm;
integer *jl, *kv, *l1, *l2, *lp, *lv, *jv;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int jfvc_();
    static integer i__, j, k, m;
    extern integer nordc_();
    static integer l10, jg, ig, mo, mt, nv;
    extern integer icf_();
    static integer l1m1;

/*<       integer kv(2,*),lp(3,*),lv(*),jv(*)                                >*/
/*<       real tb(5,*),cm(*)                                                 >*/
/*<       mo=0                                                               >*/
    /* Parameter adjustments */
    --jv;
    --lv;
    lp -= 4;
    kv -= 3;
    --cm;
    tb -= 6;

    /* Function Body */
    mo = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(icf(m,tb,cm,jl,kv,jv).ne.0) mo=max0(mo,nordc(1,m,tb,cm))        >*/
    if (icf_(&m, &tb[6], &cm[1], jl, &kv[3], &jv[1]) != 0) {
/* Computing MAX */
        i__2 = mo, i__3 = nordc_(&c__1, &m, &tb[6], &cm[1]);
        mo = max(i__2,i__3);
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       if(mo.eq.0) return                                                 >*/
    if (mo == 0) {
    return 0;
    }
/*<       do 10 mt=1,mo                                                      >*/
    i__1 = mo;
    for (mt = 1; mt <= i__1; ++mt) {
/*<       l10=l1                                                             >*/
    l10 = *l1;
/*<       do 9 m=1,nk                                                        >*/
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       if(icf(m,tb,cm,jl,kv,jv).eq.0) go to 9                             >*/
        if (icf_(&m, &tb[6], &cm[1], jl, &kv[3], &jv[1]) == 0) {
        goto L9;
        }
/*<       if(nordc(1,m,tb,cm).ne.mt) go to 9                                 >*/
        if (nordc_(&c__1, &m, &tb[6], &cm[1]) != mt) {
        goto L9;
        }
/*<       call jfvc(1,m,tb,cm,nv,jv,jv)                                      >*/
        jfvc_(&c__1, &m, &tb[6], &cm[1], &nv, &jv[1], &jv[1]);
/*<       jg=0                                                               >*/
        jg = 0;
/*<       l1m1=l1-1                                                          >*/
        l1m1 = *l1 - 1;
/*<       i=l10                                                              >*/
        i__ = l10;
/*<       go to 3                                                            >*/
        goto L3;
/*<     2 i=i+1                                                              >*/
L2:
        ++i__;
/*<     3 if((i).gt.(l1m1)) go to 7                                          >*/
L3:
        if (i__ > l1m1) {
        goto L7;
        }
/*<       k=lp(2,i)-1                                                        >*/
        k = lp[i__ * 3 + 2] - 1;
/*<       ig=0                                                               >*/
        ig = 0;
/*<       do 4 j=1,mt                                                        >*/
        i__3 = mt;
        for (j = 1; j <= i__3; ++j) {
/*<       if(jv(j).eq.lv(k+j)) go to 4                                       >*/
        if (jv[j] == lv[k + j]) {
            goto L4;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 5                                                            >*/
        goto L5;
/*<     4 continue                                                           >*/
L4:
        ;
        }
/*<     5 if(ig .ne. 0) go to 2                                              >*/
L5:
        if (ig != 0) {
        goto L2;
        }
/*<       jg=1                                                               >*/
        jg = 1;
/*<       lp(3,i)=lp(3,i)+1                                                  >*/
        ++lp[i__ * 3 + 3];
/*<     7 if(jg .ne. 0) go to 9                                              >*/
L7:
        if (jg != 0) {
        goto L9;
        }
/*<       lp(1,l1)=mt                                                        >*/
        lp[*l1 * 3 + 1] = mt;
/*<       lp(2,l1)=l2                                                        >*/
        lp[*l1 * 3 + 2] = *l2;
/*<       lp(3,l1)=1                                                         >*/
        lp[*l1 * 3 + 3] = 1;
/*<       k=l2-1                                                             >*/
        k = *l2 - 1;
/*<       do 8 i=1,mt                                                        >*/
        i__3 = mt;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       lv(i+k)=jv(i)                                                      >*/
        lv[i__ + k] = jv[i__];
/*<     8 continue                                                           >*/
/* L8: */
        }
/*<       l1=l1+1                                                            >*/
        ++(*l1);
/*<       l2=l2+mt                                                           >*/
        *l2 += mt;
/*<     9 continue                                                           >*/
L9:
        ;
    }
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* collf_ */

/*<       function icf (m,tb,cm,jl,kv,jv)                                    >*/
integer icf_(m, tb, cm, jl, kv, jv)
integer *m;
real *tb, *cm;
integer *jl, *kv, *jv;
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int jfvc_();
    static integer i__, j, k;
    static real z__;
    extern integer nordc_();
    static integer l1, l2, nc, kk, nv;

/*<       integer kv(2,jl),jv(*)                                             >*/
/*<       real tb(5,*),cm(*)                                                 >*/
/*<       icf=0                                                              >*/
    /* Parameter adjustments */
    tb -= 6;
    --cm;
    kv -= 3;
    --jv;

    /* Function Body */
    ret_val = 0;
/*<       if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.jl) return                >*/
    if (tb[*m * 5 + 1] == (float)0. || nordc_(&c__2, m, &tb[6], &cm[1]) != *
        jl) {
    return ret_val;
    }
/*<       if(jl .ne. 0) go to 1                                              >*/
    if (*jl != 0) {
    goto L1;
    }
/*<       icf=1                                                              >*/
    ret_val = 1;
/*<       return                                                             >*/
    return ret_val;
/*<     1 call jfvc(2,m,tb,cm,nv,jv,jv(jl+1))                                >*/
L1:
    jfvc_(&c__2, m, &tb[6], &cm[1], &nv, &jv[1], &jv[*jl + 1]);
/*<       do 2 j=1,jl                                                        >*/
    i__1 = *jl;
    for (j = 1; j <= i__1; ++j) {
/*<       if(iabs(jv(j)).ne.iabs(kv(1,j))) return                            >*/
    if ((i__2 = jv[j], abs(i__2)) != (i__3 = kv[(j << 1) + 1], abs(i__3)))
         {
        return ret_val;
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       do 6 j=1,jl                                                        >*/
    i__1 = *jl;
    for (j = 1; j <= i__1; ++j) {
/*<       l1=kv(2,j)                                                         >*/
    l1 = kv[(j << 1) + 2];
/*<       l2=jv(jl+j)                                                        >*/
    l2 = jv[*jl + j];
/*<       k=2*iabs(jv(j))                                                    >*/
    k = (i__2 = jv[j], abs(i__2)) << 1;
/*<       kk=jv(j)*kv(1,j)                                                   >*/
    kk = jv[j] * kv[(j << 1) + 1];
/*<       nc=int(cm(k+1)+.1)-int(cm(k)+.1)+1                                 >*/
    nc = (integer) (cm[k + 1] + (float).1) - (integer) (cm[k] + (float).1)
         + 1;
/*<       do 5 i=1,nc                                                        >*/
    i__2 = nc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       z=cm(i+l2)                                                         >*/
        z__ = cm[i__ + l2];
/*<       if(kk .ge. 0) go to 4                                              >*/
        if (kk >= 0) {
        goto L4;
        }
/*<       if(z .ne. 0.0) go to 3                                             >*/
        if (z__ != (float)0.) {
        goto L3;
        }
/*<       z=1.0                                                              >*/
        z__ = (float)1.;
/*<       go to 4                                                            >*/
        goto L4;
/*<     3 z=0.0                                                              >*/
L3:
        z__ = (float)0.;
/*<     4 if(cm(i+l1).ne.z) return                                           >*/
L4:
        if (cm[i__ + l1] != z__) {
        return ret_val;
        }
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       icf=1                                                              >*/
    ret_val = 1;
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* icf_ */

/*<       function icat (x,j,cm)                                             >*/
integer icat_(x, j, cm)
real *x;
integer *j;
real *cm;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer k, j0, j1, j2;

/*<       real cm(*)                                                         >*/
/*<       j0=cm(2*j)+.1                                                      >*/
    /* Parameter adjustments */
    --cm;

    /* Function Body */
    j0 = cm[*j * 2] + (float).1;
/*<       j1=j0                                                              >*/
    j1 = j0;
/*<       j2=cm(2*j+1)+.1                                                    >*/
    j2 = cm[(*j << 1) + 1] + (float).1;
/*<     1 if(j2.eq.j1+1) go to 5                                             >*/
L1:
    if (j2 == j1 + 1) {
    goto L5;
    }
/*<       k=(j1+j2)/2                                                        >*/
    k = (j1 + j2) / 2;
/*<       if(cm(k) .ne. x) go to 2                                           >*/
    if (cm[k] != *x) {
    goto L2;
    }
/*<       icat=k-j0+1                                                        >*/
    ret_val = k - j0 + 1;
/*<       return                                                             >*/
    return ret_val;
/*<     2 if(cm(k) .ge. x) go to 3                                           >*/
L2:
    if (cm[k] >= *x) {
    goto L3;
    }
/*<       j1=k                                                               >*/
    j1 = k;
/*<       go to 1                                                            >*/
    goto L1;
/*<     3 j2=k                                                               >*/
L3:
    j2 = k;
/*<       go to 1                                                            >*/
    goto L1;
/*<     5 if(x .ne. cm(j1)) go to 6                                          >*/
L5:
    if (*x != cm[j1]) {
    goto L6;
    }
/*<       icat=j1-j0+1                                                       >*/
    ret_val = j1 - j0 + 1;
/*<       go to 8                                                            >*/
    goto L8;
/*<     6 if(x .ne. cm(j2)) go to 7                                          >*/
L6:
    if (*x != cm[j2]) {
    goto L7;
    }
/*<       icat=j2-j0+1                                                       >*/
    ret_val = j2 - j0 + 1;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 icat=0                                                             >*/
L7:
    ret_val = 0;
/*<     8 return                                                             >*/
L8:
    return ret_val;
/*<       end                                                                >*/
} /* icat_ */

/*<    >*/
/* Subroutine */ int csp_(jp, nc, m, n, x, y, w, nk, tb, cm, kcp, yb, d__, kr,
     ntt, sw, me, mkp2, nop, sc, db, sp, mm)
integer *jp, *nc, *m, *n;
real *x, *y, *w;
integer *nk;
real *tb, *cm;
integer *kcp;
doublereal *yb, *d__;
integer *kr, *ntt;
doublereal *sw;
integer *me, *mkp2, *nop;
real *sc;
doublereal *db, *sp;
integer *mm;
{
    /* Initialized data */

    static doublereal eps = 1e-4;
    static real big = (float)9.9e30;

    /* System generated locals */
    integer mm_dim1, mm_offset, x_dim1, x_offset, d_dim1, d_offset, db_dim1, 
        db_offset, sp_dim1, sp_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal a, b;
    static real h__;
    static integer i__, j, k;
    static doublereal s;
    static integer k1, n1, jj, mk;
    static doublereal dv, dy;
    static real wh;
    static integer ns, js, nr, nrt;
    static real bof0, bof1;
    static integer mkp1;

/*<       integer mm(nc,2)                                                   >*/
/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n)                         >*/
/*<       double precision yb,sw,d(nk,*),db(n,*),sp(mkp2,*),a,b,s,eps,dv,dy  >*/
/*<       data eps,big /1.d-4,9.9e30/                                        >*/
    /* Parameter adjustments */
    mm_dim1 = *nc;
    mm_offset = mm_dim1 + 1;
    mm -= mm_offset;
    db_dim1 = *n;
    db_offset = db_dim1 + 1;
    db -= db_offset;
    --sc;
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    tb -= 6;
    --cm;
    sp_dim1 = *mkp2;
    sp_offset = sp_dim1 + 1;
    sp -= sp_offset;

    /* Function Body */
/*<       nop=0                                                              >*/
    *nop = 0;
/*<       if(nc .gt. 1) go to 1                                              >*/
    if (*nc > 1) {
    goto L1;
    }
/*<       tb(1,m)=big                                                        >*/
    tb[*m * 5 + 1] = big;
/*<       return                                                             >*/
    return 0;
/*<     1 mk=mkp2-2                                                          >*/
L1:
    mk = *mkp2 - 2;
/*<       mkp1=mk+1                                                          >*/
    mkp1 = mk + 1;
/*<       n1=nc+1                                                            >*/
    n1 = *nc + 1;
/*<       do 3 j=1,n1                                                        >*/
    i__1 = n1;
    for (j = 1; j <= i__1; ++j) {
/*<       do 2 i=1,mkp2                                                      >*/
    i__2 = *mkp2;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       sp(i,j)=0.d0                                                       >*/
        sp[i__ + j * sp_dim1] = 0.;
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       do 4 j=1,nc                                                        >*/
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
/*<       mm(j,2)=0                                                          >*/
    mm[j + (mm_dim1 << 1)] = 0;
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       do 7 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       h=sc(i)                                                            >*/
    h__ = sc[i__];
/*<       if(h.le.0.0.or.w(i).le.0.0) go to 7                                >*/
    if (h__ <= (float)0. || w[i__] <= (float)0.) {
        goto L7;
    }
/*<       wh=w(i)*h                                                          >*/
    wh = w[i__] * h__;
/*<       k=x(i,jp)+.1                                                       >*/
    k = x[i__ + *jp * x_dim1] + (float).1;
/*<       mm(k,2)=mm(k,2)+1                                                  >*/
    ++mm[k + (mm_dim1 << 1)];
/*<       sp(mkp2,k)=sp(mkp2,k)+wh                                           >*/
    sp[*mkp2 + k * sp_dim1] += wh;
/*<       sp(mkp1,k)=sp(mkp1,k)+wh*(y(i)-yb)                                 >*/
    sp[mkp1 + k * sp_dim1] += wh * (y[i__] - *yb);
/*<       sp(m,k)=sp(m,k)+wh*h                                               >*/
    sp[*m + k * sp_dim1] += wh * h__;
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 j=j+1                                                              >*/
L5:
    ++j;
/*<     6 if((j).gt.(kr)) go to 7                                            >*/
L6:
    if (j > *kr) {
        goto L7;
    }
/*<       sp(j,k)=sp(j,k)+wh*db(i,j)                                         >*/
    sp[j + k * sp_dim1] += wh * db[i__ + j * db_dim1];
/*<       go to 5                                                            >*/
    goto L5;
/*<     7 continue                                                           >*/
L7:
    ;
    }
/*<       do 8 j=1,nc                                                        >*/
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
/*<       mm(j,1)=j                                                          >*/
    mm[j + mm_dim1] = j;
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       bof0=big                                                           >*/
    bof0 = big;
/*<       ns=0                                                               >*/
    ns = 0;
/*<       jj=nc                                                              >*/
    jj = *nc;
/*<       nrt=0                                                              >*/
    nrt = 0;
/*<       k1=1                                                               >*/
    k1 = 1;
/*<     9 bof1=big                                                           >*/
L9:
    bof1 = big;
/*<       js=0                                                               >*/
    js = 0;
/*<       do 18 j=1,jj                                                       >*/
    i__1 = jj;
    for (j = 1; j <= i__1; ++j) {
/*<       k=mm(j,1)                                                          >*/
    k = mm[j + mm_dim1];
/*<       if(mm(k,2).eq.0) go to 18                                          >*/
    if (mm[k + (mm_dim1 << 1)] == 0) {
        goto L18;
    }
/*<       nr=nrt+mm(k,2)                                                     >*/
    nr = nrt + mm[k + (mm_dim1 << 1)];
/*<       if(nr.le.me.or.ntt-nr.le.me) go to 18                              >*/
    if (nr <= *me || *ntt - nr <= *me) {
        goto L18;
    }
/*<       dy=sp(mkp1,n1)+sp(mkp1,k)                                          >*/
    dy = sp[mkp1 + n1 * sp_dim1] + sp[mkp1 + k * sp_dim1];
/*<       a=sp(mkp2,n1)+sp(mkp2,k)                                           >*/
    a = sp[*mkp2 + n1 * sp_dim1] + sp[*mkp2 + k * sp_dim1];
/*<       dv=sp(m,n1)+sp(m,k)-a**2/sw                                        >*/
/* Computing 2nd power */
    d__1 = a;
    dv = sp[*m + n1 * sp_dim1] + sp[*m + k * sp_dim1] - d__1 * d__1 / *sw;
/*<       if(dv .le. 0.d0) go to 17                                          >*/
    if (dv <= 0.) {
        goto L17;
    }
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 11                                                           >*/
    goto L11;
/*<    10 i=i+1                                                              >*/
L10:
    ++i__;
/*<    11 if((i).gt.(kr)) go to 12                                           >*/
L11:
    if (i__ > *kr) {
        goto L12;
    }
/*<       d(i,2)=sp(i,n1)+sp(i,k)                                            >*/
    d__[i__ + (d_dim1 << 1)] = sp[i__ + n1 * sp_dim1] + sp[i__ + k * 
        sp_dim1];
/*<       go to 10                                                           >*/
    goto L10;
/*<    12 a=0.d0                                                             >*/
L12:
    a = 0.;
/*<       b=a                                                                >*/
    b = a;
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 14                                                           >*/
    goto L14;
/*<    13 i=i+1                                                              >*/
L13:
    ++i__;
/*<    14 if((i).gt.(kr)) go to 15                                           >*/
L14:
    if (i__ > *kr) {
        goto L15;
    }
/*<       s=d(i,2)                                                           >*/
    s = d__[i__ + (d_dim1 << 1)];
/*<       a=a+s*d(i,1)                                                       >*/
    a += s * d__[i__ + d_dim1];
/*<       b=b+s**2                                                           >*/
/* Computing 2nd power */
    d__1 = s;
    b += d__1 * d__1;
/*<       go to 13                                                           >*/
    goto L13;
/*<    15 b=dv-b                                                             >*/
L15:
    b = dv - b;
/*<       if(b .le. eps*dv) go to 17                                         >*/
    if (b <= eps * dv) {
        goto L17;
    }
/*<       nop=nop+1                                                          >*/
    ++(*nop);
/*<       b=-(dy-a)**2/b                                                     >*/
/* Computing 2nd power */
    d__1 = dy - a;
    b = -(d__1 * d__1) / b;
/*<       if(b .ge. bof1) go to 16                                           >*/
    if (b >= bof1) {
        goto L16;
    }
/*<       bof1=b                                                             >*/
    bof1 = b;
/*<       js=j                                                               >*/
    js = j;
/*<    16 if(b .ge. bof0) go to 17                                           >*/
L16:
    if (b >= bof0) {
        goto L17;
    }
/*<       bof0=b                                                             >*/
    bof0 = b;
/*<       ns=jj                                                              >*/
    ns = jj;
/*<    17 if(nc.eq.2) go to 19                                               >*/
L17:
    if (*nc == 2) {
        goto L19;
    }
/*<    18 continue                                                           >*/
L18:
    ;
    }
/*<    19 if(js.eq.0) go to 23                                               >*/
L19:
    if (js == 0) {
    goto L23;
    }
/*<       k=mm(js,1)                                                         >*/
    k = mm[js + mm_dim1];
/*<       mm(js,1)=mm(jj,1)                                                  >*/
    mm[js + mm_dim1] = mm[jj + mm_dim1];
/*<       mm(jj,1)=k                                                         >*/
    mm[jj + mm_dim1] = k;
/*<       sp(mkp1,n1)=sp(mkp1,n1)+sp(mkp1,k)                                 >*/
    sp[mkp1 + n1 * sp_dim1] += sp[mkp1 + k * sp_dim1];
/*<       sp(mkp2,n1)=sp(mkp2,n1)+sp(mkp2,k)                                 >*/
    sp[*mkp2 + n1 * sp_dim1] += sp[*mkp2 + k * sp_dim1];
/*<       nrt=nrt+mm(k,2)                                                    >*/
    nrt += mm[k + (mm_dim1 << 1)];
/*<       sp(m,n1)=sp(m,n1)+sp(m,k)                                          >*/
    sp[*m + n1 * sp_dim1] += sp[*m + k * sp_dim1];
/*<       i=1                                                                >*/
    i__ = 1;
/*<       go to 21                                                           >*/
    goto L21;
/*<    20 i=i+1                                                              >*/
L20:
    ++i__;
/*<    21 if((i).gt.(kr)) go to 22                                           >*/
L21:
    if (i__ > *kr) {
    goto L22;
    }
/*<       sp(i,n1)=sp(i,n1)+sp(i,k)                                          >*/
    sp[i__ + n1 * sp_dim1] += sp[i__ + k * sp_dim1];
/*<       go to 20                                                           >*/
    goto L20;
/*<    22 jj=jj-1                                                            >*/
L22:
    --jj;
/*<       if(jj.le.2) go to 23                                               >*/
    if (jj <= 2) {
    goto L23;
    }
/*<       go to 9                                                            >*/
    goto L9;
/*<    23 tb(1,m)=bof0                                                       >*/
L23:
    tb[*m * 5 + 1] = bof0;
/*<       tb(3,m)=kcp                                                        >*/
    tb[*m * 5 + 3] = (real) (*kcp);
/*<       do 24 j=1,nc                                                       >*/
    i__1 = *nc;
    for (j = 1; j <= i__1; ++j) {
/*<       cm(j+kcp)=0.0                                                      >*/
    cm[j + *kcp] = (float)0.;
/*<    24 continue                                                           >*/
/* L24: */
    }
/*<       if(ns.eq.0) return                                                 >*/
    if (ns == 0) {
    return 0;
    }
/*<       do 25 j=ns,nc                                                      >*/
    i__1 = *nc;
    for (j = ns; j <= i__1; ++j) {
/*<       cm(mm(j,1)+kcp)=1.0                                                >*/
    cm[mm[j + mm_dim1] + *kcp] = (float)1.;
/*<    25 continue                                                           >*/
/* L25: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* csp_ */

/*<       subroutine rspnpr (it,il,n,y,w,m)                                  >*/
/* Subroutine */ int rspnpr_(it, il, n, y, w, m)
integer *it, *il, *n;
real *y, *w;
integer *m;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int psort_();
    static real wm[2], wt;

/*<       integer m(n)                                                       >*/
/*<       real y(n),w(n),wm(2)                                               >*/
/*<       if(it.le.0) return                                                 >*/
    /* Parameter adjustments */
    --m;
    --w;
    --y;

    /* Function Body */
    if (*it <= 0) {
    return 0;
    }
/*<       if(il .ne. 1) go to 2                                              >*/
    if (*il != 1) {
    goto L2;
    }
/*<       wm(1)=0.0                                                          >*/
    wm[0] = (float)0.;
/*<       wm(2)=wm(1)                                                        >*/
    wm[1] = wm[0];
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       k=y(i)+1.1                                                         >*/
    k = y[i__] + (float)1.1;
/*<       wm(k)=wm(k)+w(i)                                                   >*/
    wm[k - 1] += w[i__];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       wt=wm(1)+wm(2)                                                     >*/
    wt = wm[0] + wm[1];
/*<       wm(1)=wm(1)/wt                                                     >*/
    wm[0] /= wt;
/*<       wm(2)=wm(2)/wt                                                     >*/
    wm[1] /= wt;
/*    write(it,'(/,'' binary (0/1) response:  mass(0) ='',g12.4,         3
652*/
/*   1                           ''   mass(1) ='',g12.4)') wm(1),wm(2)   3
653*/
/*<       return                                                             >*/
    return 0;
/*<     2 continue >*/
L2:
/*    write(it,'(/,'' ordinal response:'')')                             3
655*/
/*    write(it,'(''      min         n/4         n/2        3n/4         3
656*/
/*   1 max'')')                                                          3
657*/
/*<       do 3 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       m(i)=i                                                             >*/
    m[i__] = i__;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       call psort(y,m,1,n)                                                >*/
    psort_(&y[1], &m[1], &c__1, n);
/*    write(it,'('' '',5g12.4)') y(m(1)),y(m(n/4)),y(m(n/2)),y(m(n-n/4)) 3
662*/
/*   1,y(m(n))                                                           3
663*/
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* rspnpr_ */

/*<       subroutine ordpr (it,n,p,x,lx,m)                                   >*/
/* Subroutine */ int ordpr_(it, n, p, x, lx, m)
integer *it, *n, *p;
real *x;
integer *lx, *m;
{
    /* System generated locals */
    integer m_dim1, m_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    static integer j, n1, n2, n3, no;

/*<       integer p,lx(p),m(n,p)                                             >*/
/*<       real x(n,p)                                                        >*/
/*<       if(it.le.0) return                                                 >*/
    /* Parameter adjustments */
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    --lx;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    if (*it <= 0) {
    return 0;
    }
/*<       no=0                                                               >*/
    no = 0;
/*<       do 1 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       if(lx(j).gt.0) no=no+1                                             >*/
    if (lx[j] > 0) {
        ++no;
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       if(no.eq.0) return                                                 >*/
    if (no == 0) {
    return 0;
    }
/*    write(it,'(/,'' there are'',i3,'' ordinal predictor variables.'',/ 3
675*/
/*   1)') no                                                             3
676*/
/*    write(it,'(''  var     min         n/4         n/2        3n/4     3
677*/
/*   1     max'')')                                                      3
678*/
/*<       n1=n/4                                                             >*/
    n1 = *n / 4;
/*<       n2=n/2                                                             >*/
    n2 = *n / 2;
/*<       n3=n-n1                                                            >*/
    n3 = *n - n1;
/*<       do 2 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       if(lx(j).le.0) go to 2                                             >*/
    if (lx[j] <= 0) {
        goto L2;
    }
/*    write(it,'('' '',i3,'' '',5g12.4)') j,  x(m(1,j),j),x(m(n1,j),j)
,x 3684*/
/*   1(m(n2,j),j),x(m(n3,j),j),x(m(n,j),j)                            
   3685*/
/*<     2 continue                                                           >*/
L2:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* ordpr_ */

/*<       subroutine catpr(it,n,p,x,cm,mm)                                   >*/
/* Subroutine */ int catpr_(it, n, p, x, cm, mm)
integer *it, *n, *p;
real *x, *cm;
integer *mm;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1;

    /* Local variables */
    static integer i__, j, k, j1, j2, n2, ic, np, nv, nct;

/*<       integer p,mm(*)                                                    >*/
/*<       real x(n,p),cm(*)                                                  >*/
/*<       if(it.le.0) return                                                 >*/
    /* Parameter adjustments */
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --cm;
    --mm;

    /* Function Body */
    if (*it <= 0) {
    return 0;
    }
/*<       nct=cm(1)+.1                                                       >*/
    nct = cm[1] + (float).1;
/*<       if(nct.eq.0) return                                                >*/
    if (nct == 0) {
    return 0;
    }
/*<       n2=2*p+1                                                           >*/
    n2 = (*p << 1) + 1;
/*<       np=0                                                               >*/
    np = 0;
/*    write(it,'(/,'' there are'',i3,'' categorical predictor variables. 3
697*/
/*   1'')') nct                                                          3
698*/
/*<       i=2                                                                >*/
    i__ = 2;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 i=i+(2)                                                            >*/
L1:
    i__ += 2;
/*<     2 if((2)*((i)-(n2)).gt.0) go to 6                                    >*/
L2:
    if (i__ - n2 << 1 > 0) {
    goto L6;
    }
/*<       np=np+1                                                            >*/
    ++np;
/*<       j1=cm(i)+.1                                                        >*/
    j1 = cm[i__] + (float).1;
/*<       if(j1.eq.0) go to 1                                                >*/
    if (j1 == 0) {
    goto L1;
    }
/*<       j2=cm(i+1)+.1                                                      >*/
    j2 = cm[i__ + 1] + (float).1;
/*<       nv=j2-j1+1                                                         >*/
    nv = j2 - j1 + 1;
/*<       do 3 j=1,nv                                                        >*/
    i__1 = nv;
    for (j = 1; j <= i__1; ++j) {
/*<       mm(j)=0                                                            >*/
    mm[j] = 0;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       do 4 j=1,n                                                         >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<       ic=x(j,np)+.1                                                      >*/
    ic = x[j + np * x_dim1] + (float).1;
/*<       mm(ic)=mm(ic)+1                                                    >*/
    ++mm[ic];
/*<     4 continue                                                           >*/
/* L4: */
    }
/*    write(it,'(/,'' categorical variable'',i3,'' has'',i3,'' values.'' 3
715*/
/*   1)') np,nv                                                          3
716*/
/*    write(it,'(''  value     internal code     counts'')')             3
717*/
/*<       k=0                                                                >*/
    k = 0;
/*<       do 5 j=j1,j2                                                       >*/
    i__1 = j2;
    for (j = j1; j <= i__1; ++j) {
/*<       k=k+1                                                              >*/
    ++k;
/*    write(it,'(f6.0,i13,i15)') cm(j),k,mm(k)                        
   3721*/
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       go to 1                                                            >*/
    goto L1;
/*<     6 return                                                             >*/
L6:
    return 0;
/*<       end                                                                >*/
} /* catpr_ */

/*<       subroutine holl (jp,cm,t,h)                                        >*/
/* Subroutine */ int holl_(jp, cm, t, h__, h_length)
integer *jp;
real *cm, *t;
char *h__;
ftnlen h_length;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer j, k, j1, j2;

/*<       real cm(*)                                                         >*/
/*<       character*28 h                                                     >*/
/*<       j1=cm(2*jp)+.1                                                     >*/
    /* Parameter adjustments */
    --cm;

    /* Function Body */
    j1 = cm[*jp * 2] + (float).1;
/*<       j2=cm(2*jp+1)+.1                                                   >*/
    j2 = cm[(*jp << 1) + 1] + (float).1;
/*<       j2=j2-j1+1                                                         >*/
    j2 = j2 - j1 + 1;
/*<       if(j2 .le. 28) go to 1                                             >*/
    if (j2 <= 28) {
    goto L1;
    }
/*<       h='   cat. factor > 28 values  '                                   >*/
    s_copy(h__, "   cat. factor > 28 values  ", 28L, 28L);
/*<       return                                                             >*/
    return 0;
/*<     1 h='                            '                                   >*/
L1:
    s_copy(h__, "                            ", 28L, 28L);
/*<       j1=(28-j2)/2                                                       >*/
    j1 = (28 - j2) / 2;
/*<       j2=j1+j2-1                                                         >*/
    j2 = j1 + j2 - 1;
/*<       k=t+.1                                                             >*/
    k = *t + (float).1;
/*<       do 3 j=j1,j2                                                       >*/
    i__1 = j2;
    for (j = j1; j <= i__1; ++j) {
/*<       if(cm(k+j-j1+1) .le. 0.0) go to 2                                  >*/
    if (cm[k + j - j1 + 1] <= (float)0.) {
        goto L2;
    }
/*<       h(j:j)='1'                                                         >*/
    *(unsigned char *)&h__[j - 1] = '1';
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 h(j:j)='0'                                                         >*/
L2:
    *(unsigned char *)&h__[j - 1] = '0';
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* holl_ */

/*<       subroutine slova (nk,it,tb,ni,lp,lv)                               >*/
/* Subroutine */ int slova_(nk, it, tb, ni, lp, lv)
integer *nk, *it;
real *tb;
integer *ni, *lp, *lv;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int coll_();
    static integer m, i2, k2, na;

/*<       integer lp(3,*),lv(*)                                              >*/
/*<       real tb(5,nk)                                                      >*/
/*    write(it,4) ni                                                     3
750*/
/*<       call coll(nk,tb,lp,lv,lp(1,nk+1))                                  >*/
    /* Parameter adjustments */
    tb -= 6;
    lp -= 4;
    --lv;

    /* Function Body */
    coll_(nk, &tb[6], &lp[4], &lv[1], &lp[(*nk + 1) * 3 + 1]);
/*<       m=1                                                                >*/
    m = 1;
/*<     1 if(lp(1,m).eq.0) go to 2                                           >*/
L1:
    if (lp[m * 3 + 1] == 0) {
    goto L2;
    }
/*<       m=m+1                                                              >*/
    ++m;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 na=m-1                                                             >*/
L2:
    na = m - 1;
/*<       do 3 m=1,na                                                        >*/
    i__1 = na;
    for (m = 1; m <= i__1; ++m) {
/*<       k2=lp(2,m)                                                         >*/
    k2 = lp[m * 3 + 2];
/*<       i2=lp(1,m)+k2-1                                                    >*/
    i2 = lp[m * 3 + 1] + k2 - 1;
/*    write(it,5) m,lp(3,m),(lv(i),i=k2,i2)                           
   3760*/
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       return                                                             >*/
    return 0;
/*<    >*/
/* L4: */
/*<     5 format('  ',i3,'         ',i2,'       ',20i4)                      >*/
/* L5: */
/*<       end                                                                >*/
} /* slova_ */

/*<       subroutine reducq (flg,x,nk,tb,cm,tc,kp,kv,lp,lv,r,td,sc,fc)       >*/
/* Subroutine */ int reducq_(flg, x, nk, tb, cm, tc, kp, kv, lp, lv, r__, td, 
    sc, fc)
real *flg, *x;
integer *nk;
real *tb, *cm, *tc;
integer *kp, *kv, *lp, *lv;
real *r__, *td, *sc, *fc;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int gtrm_();
    static integer k, l, m;
    extern integer match_();
    static integer l1, la, il, ll, jl, jp, nt, nv, kp3, laa;
    extern /* Subroutine */ int std_();

/*<       integer kp(5,*),kv(2,*),lp(3,*),lv(*)                              >*/
/*<       real x(*),tb(5,nk),cm(*),tc(*),r(*),td(2,nk),sc(2,*),fc(*)         >*/
/*<       ll=1                                                               >*/
    /* Parameter adjustments */
    --x;
    td -= 3;
    tb -= 6;
    --cm;
    --tc;
    kp -= 6;
    kv -= 3;
    lp -= 4;
    --lv;
    --r__;
    sc -= 3;
    --fc;

    /* Function Body */
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<       laa=0                                                              >*/
    laa = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       td(1,m)=0.0                                                        >*/
    td[(m << 1) + 1] = (float)0.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 if(kp(1,ll).lt.0) go to 9                                          >*/
L2:
    if (kp[ll * 5 + 1] < 0) {
    goto L9;
    }
/*<       nv=0                                                               >*/
    nv = 0;
/*<       if(kp(1,ll) .le. 0) go to 4                                        >*/
    if (kp[ll * 5 + 1] <= 0) {
    goto L4;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       do 3 il=1,jl                                                       >*/
    i__1 = jl;
    for (il = 1; il <= i__1; ++il) {
/*<       k=kp(2,ll)+il-1                                                    >*/
    k = kp[ll * 5 + 2] + il - 1;
/*<       nv=nv+1                                                            >*/
    ++nv;
/*<       sc(1,nv)=kv(1,k)                                                   >*/
    sc[(nv << 1) + 1] = (real) kv[(k << 1) + 1];
/*<       sc(2,nv)=kv(2,k)                                                   >*/
    sc[(nv << 1) + 2] = (real) kv[(k << 1) + 2];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       go to 5                                                            >*/
    goto L5;
/*<     4 if(kp(3,ll) .gt. 0) go to 5                                        >*/
L4:
    if (kp[ll * 5 + 3] > 0) {
    goto L5;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     5 if(kp(3,ll) .gt. 0) go to 6                                        >*/
L5:
    if (kp[ll * 5 + 3] > 0) {
    goto L6;
    }
/*<       m=match(nv,sc,nk,tb,cm,r,0)                                        >*/
    m = match_(&nv, &sc[3], nk, &tb[6], &cm[1], &r__[1], &c__0);
/*<       td(1,m)=tc(-kp(3,ll))                                              >*/
    td[(m << 1) + 1] = tc[-kp[ll * 5 + 3]];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     6 kp3=kp(3,ll)                                                       >*/
L6:
    kp3 = kp[ll * 5 + 3];
/*<       do 8 k=1,kp3                                                       >*/
    i__1 = kp3;
    for (k = 1; k <= i__1; ++k) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       laa=laa+5*l*nt                                                     >*/
    laa += l * 5 * nt;
/*<       do 7 jp=1,nt                                                       >*/
    i__2 = nt;
    for (jp = 1; jp <= i__2; ++jp) {
/*<       call gtrm(1,jp,l,nt,lv(lp(2,l1)),flg,x,nk,tb,tc(la),sc(1,nv+1),fc) >*/
        gtrm_(&c__1, &jp, &l, &nt, &lv[lp[l1 * 3 + 2]], flg, &x[1], nk, &
            tb[6], &tc[la], &sc[(nv + 1 << 1) + 1], &fc[1]);
/*<       m=match(nv+l,sc,nk,tb,cm,r,0)                                      >*/
        i__3 = nv + l;
        m = match_(&i__3, &sc[3], nk, &tb[6], &cm[1], &r__[1], &c__0);
/*<       td(1,m)=tc(jp+laa)                                                 >*/
        td[(m << 1) + 1] = tc[jp + laa];
/*<       call std(m,flg,x,l,sc(1,nv+1),fc,nk,tb,r,td)                       >*/
        std_(&m, flg, &x[1], &l, &sc[(nv + 1 << 1) + 1], &fc[1], nk, &tb[
            6], &r__[1], &td[3]);
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<       laa=laa+nt                                                         >*/
    laa += nt;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       la=la+nt*(5*l+1)                                                   >*/
    la += nt * (l * 5 + 1);
/*<     8 continue                                                           >*/
/* L8: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     9 return                                                             >*/
L9:
    return 0;
/*<       end                                                                >*/
} /* reducq_ */

/*<       subroutine gtrm (il,jp,l,nt,jv,flg,x,nk,tb,tc,te,fc)               >*/
/* Subroutine */ int gtrm_(il, jp, l, nt, jv, flg, x, nk, tb, tc, te, fc)
integer *il, *jp, *l, *nt, *jv;
real *flg, *x;
integer *nk;
real *tb, *tc, *te, *fc;
{
    /* System generated locals */
    integer tc_dim1, tc_offset, i__1;

    /* Local variables */
    static integer j, k, l2, l3, l4, nf, jj;
    extern doublereal cue_();

/*<       integer jv(l)                                                      >*/
/*<       real x(*),tb(5,nk),tc(nt,*),te(2,*),fc(*)                          >*/
/*<       l2=l+l                                                             >*/
    /* Parameter adjustments */
    --jv;
    tc_dim1 = *nt;
    tc_offset = tc_dim1 + 1;
    tc -= tc_offset;
    --x;
    tb -= 6;
    te -= 3;
    --fc;

    /* Function Body */
    l2 = *l + *l;
/*<       nf=0                                                               >*/
    nf = 0;
/*<       l3=l2+l                                                            >*/
    l3 = l2 + *l;
/*<       l4=l3+l                                                            >*/
    l4 = l3 + *l;
/*<       do 1 k=1,l                                                         >*/
    i__1 = *l;
    for (k = 1; k <= i__1; ++k) {
/*<       j=jv(k)                                                            >*/
    j = jv[k];
/*<       jj=j                                                               >*/
    jj = j;
/*<       if(tc(jp,k+l).gt.tc(jp,k+l2)) jj=-jj                               >*/
    if (tc[*jp + (k + *l) * tc_dim1] > tc[*jp + (k + l2) * tc_dim1]) {
        jj = -jj;
    }
/*<       te(1,k)=jj                                                         >*/
    te[(k << 1) + 1] = (real) jj;
/*<       te(2,k)=tc(jp,k)                                                   >*/
    te[(k << 1) + 2] = tc[*jp + k * tc_dim1];
/*<       if(il.eq.2) go to 1                                                >*/
    if (*il == 2) {
        goto L1;
    }
/*<       if(x(j).eq.flg) go to 1                                            >*/
    if (x[j] == *flg) {
        goto L1;
    }
/*<       nf=nf+1                                                            >*/
    ++nf;
/*<    >*/
    fc[nf] = cue_(&x[j], &tc[*jp + (k + *l) * tc_dim1], &tc[*jp + k * 
        tc_dim1], &tc[*jp + (k + l2) * tc_dim1], &tc[*jp + (k + l3) * 
        tc_dim1], &tc[*jp + (k + l4) * tc_dim1]);
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* gtrm_ */

/*<       function match (nv,te,nk,tb,cm,r,iz)                               >*/
integer match_(nv, te, nk, tb, cm, r__, iz)
integer *nv;
real *te;
integer *nk;
real *tb, *cm, *r__;
integer *iz;
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern integer nord_();
    static integer i__, j, m;
    static real t, u;
    static integer i1, i2, j1, j2, jg, ig, nc, kg, jp, ip, jq, jp2, jp21;
    extern integer ieq_();
    static integer jpp, jqq;

/*<       real te(2,nv),tb(5,nk),cm(*),r(*)                                  >*/
/*<       match=0                                                            >*/
    /* Parameter adjustments */
    te -= 3;
    tb -= 6;
    --cm;
    --r__;

    /* Function Body */
    ret_val = 0;
/*<       do 15 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       if(tb(1,m).eq.0.0) go to 15                                        >*/
    if (tb[m * 5 + 1] == (float)0.) {
        goto L15;
    }
/*<       if(nord(m,tb).ne.nv) go to 15                                      >*/
    if (nord_(&m, &tb[6]) != *nv) {
        goto L15;
    }
/*<       jg=0                                                               >*/
    jg = 0;
/*<       do 13 j=1,nv                                                       >*/
    i__2 = *nv;
    for (j = 1; j <= i__2; ++j) {
/*<       t=te(1,j)                                                          >*/
        t = te[(j << 1) + 1];
/*<       u=te(2,j)                                                          >*/
        u = te[(j << 1) + 2];
/*<       jp=abs(t)+.1                                                       >*/
        jp = dabs(t) + (float).1;
/*<       jp2=2*jp                                                           >*/
        jp2 = jp << 1;
/*<       jp21=jp2+1                                                         >*/
        jp21 = jp2 + 1;
/*<       jpp=jp                                                             >*/
        jpp = jp;
/*<       if(t.lt.0.0) jpp=-jpp                                              >*/
        if (t < (float)0.) {
        jpp = -jpp;
        }
/*<       ig=0                                                               >*/
        ig = 0;
/*<       ip=m                                                               >*/
        ip = m;
/*<     1 if(ip.le.0) go to 12                                               >*/
L1:
        if (ip <= 0) {
        goto L12;
        }
/*<       t=tb(2,ip)                                                         >*/
        t = tb[ip * 5 + 2];
/*<       jq=abs(t)+.1                                                       >*/
        jq = dabs(t) + (float).1;
/*<       jqq=jq                                                             >*/
        jqq = jq;
/*<       if(t.lt.0.0) jqq=-jqq                                              >*/
        if (t < (float)0.) {
        jqq = -jqq;
        }
/*<       if(jp .eq. jq) go to 2                                             >*/
        if (jp == jq) {
        goto L2;
        }
/*<       ip=tb(4,ip)+.1                                                     >*/
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
        goto L1;
/*<     2 if(cm(jp2) .ne. 0.0) go to 4                                       >*/
L2:
        if (cm[jp2] != (float)0.) {
        goto L4;
        }
/*<       if(jpp .ne. jqq .or. ieq(tb(3,ip),u,r(jp)) .ne. 1) go to 3         >*/
        if (jpp != jqq || ieq_(&tb[ip * 5 + 3], &u, &r__[jp]) != 1) {
        goto L3;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 12                                                           >*/
        goto L12;
/*<     3 ip=tb(4,ip)+.1                                                     >*/
L3:
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
        goto L1;
/*<     4 nc=int(cm(jp21)+.1)-int(cm(jp2)+.1)+1                              >*/
L4:
        nc = (integer) (cm[jp21] + (float).1) - (integer) (cm[jp2] + (
            float).1) + 1;
/*<       i1=u+.1                                                            >*/
        i1 = u + (float).1;
/*<       i2=tb(3,ip)+.1                                                     >*/
        i2 = tb[ip * 5 + 3] + (float).1;
/*<       kg=0                                                               >*/
        kg = 0;
/*<       do 9 i=1,nc                                                        >*/
        i__3 = nc;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       j1=cm(i1+i)                                                        >*/
        j1 = cm[i1 + i__];
/*<       j2=cm(i2+i)                                                        >*/
        j2 = cm[i2 + i__];
/*<       if(jpp .ge. 0) go to 6                                             >*/
        if (jpp >= 0) {
            goto L6;
        }
/*<       if(j1 .ne. 0) go to 5                                              >*/
        if (j1 != 0) {
            goto L5;
        }
/*<       j1=1                                                               >*/
        j1 = 1;
/*<       go to 6                                                            >*/
        goto L6;
/*<     5 j1=0                                                               >*/
L5:
        j1 = 0;
/*<     6 if(jqq .ge. 0) go to 8                                             >*/
L6:
        if (jqq >= 0) {
            goto L8;
        }
/*<       if(j2 .ne. 0) go to 7                                              >*/
        if (j2 != 0) {
            goto L7;
        }
/*<       j2=1                                                               >*/
        j2 = 1;
/*<       go to 8                                                            >*/
        goto L8;
/*<     7 j2=0                                                               >*/
L7:
        j2 = 0;
/*<     8 if(j1 .eq. j2) go to 9                                             >*/
L8:
        if (j1 == j2) {
            goto L9;
        }
/*<       kg=1                                                               >*/
        kg = 1;
/*<       go to 10                                                           >*/
        goto L10;
/*<     9 continue                                                           >*/
L9:
        ;
        }
/*<    10 if(kg .ne. 0) go to 11                                             >*/
L10:
        if (kg != 0) {
        goto L11;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 12                                                           >*/
        goto L12;
/*<    11 ip=tb(4,ip)+.1                                                     >*/
L11:
        ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
        goto L1;
/*<    12 if(ig .ne. 0) go to 13                                             >*/
L12:
        if (ig != 0) {
        goto L13;
        }
/*<       jg=1                                                               >*/
        jg = 1;
/*<       go to 14                                                           >*/
        goto L14;
/*<    13 continue                                                           >*/
L13:
        ;
    }
/*<    14 if(jg .ne. 0) go to 15                                             >*/
L14:
    if (jg != 0) {
        goto L15;
    }
/*<       match=m                                                            >*/
    ret_val = m;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 continue                                                           >*/
L15:
    ;
    }
/*<    16 if(match.gt.0.or.iz.ne.0) return                                   >*/
L16:
    if (ret_val > 0 || *iz != 0) {
    return ret_val;
    }
/*    write(6,17)                                                        3
902*/
/*<    17 format (' bug in match - term not found.')                         >*/
/* L17: */
/*<       do 19 j=1,nv                                                       >*/
    i__1 = *nv;
    for (j = 1; j <= i__1; ++j) {
/*    write(6,18)j,te(1,j),te(2,j)                                    
   3905*/
/*<    18 format (' te(',i2,')=',2g12.4)                                     >*/
/* L18: */
/*<    19 continue                                                           >*/
/* L19: */
    }
/*<       do 21 j=1,nk                                                       >*/
    i__1 = *nk;
    for (j = 1; j <= i__1; ++j) {
/*    write(6,20)j,(tb(i,j),i=1,4)                                    
   3909*/
/*<    20 format (' tb(',i2,')=',4g12.4)                                     >*/
/* L20: */
/*<    21 continue                                                           >*/
/* L21: */
    }
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<       end                                                                >*/
    return ret_val;
} /* match_ */

/*<       subroutine std (m,flg,x,l,te,fc,nk,tb,r,td)                        >*/
/* Subroutine */ int std_(m, flg, x, l, te, fc, nk, tb, r__, td)
integer *m;
real *flg, *x;
integer *l;
real *te, *fc;
integer *nk;
real *tb, *r__, *td;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static real t, u;
    static integer ig, jj, ip, jp;
    extern integer ieq_();

/*<       real x(*),te(2,*),fc(*),tb(5,nk),r(*),td(2,nk)                     >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    --x;
    te -= 3;
    --fc;
    td -= 3;
    tb -= 6;
    --r__;

    /* Function Body */
    ip = *m;
/*<     1 if(ip.le.0) go to 4                                                >*/
L1:
    if (ip <= 0) {
    goto L4;
    }
/*<       t=tb(2,ip)                                                         >*/
    t = tb[ip * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       jj=j                                                               >*/
    jj = j;
/*<       if(t.lt.0.0) jj=-jj                                                >*/
    if (t < (float)0.) {
    jj = -jj;
    }
/*<       u=tb(3,ip)                                                         >*/
    u = tb[ip * 5 + 3];
/*<       k=0                                                                >*/
    k = 0;
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 2 i=1,l                                                         >*/
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       t=te(1,i)                                                          >*/
    t = te[(i__ << 1) + 1];
/*<       jp=abs(t)+.1                                                       >*/
    jp = dabs(t) + (float).1;
/*<       if(x(jp).ne.flg) k=k+1                                             >*/
    if (x[jp] != *flg) {
        ++k;
    }
/*<       if(t.lt.0.0) jp=-jp                                                >*/
    if (t < (float)0.) {
        jp = -jp;
    }
/*<       if(jj .ne. jp .or. ieq(te(2,i),u,r(j)) .ne. 1) go to 2             >*/
    if (jj != jp || ieq_(&te[(i__ << 1) + 2], &u, &r__[j]) != 1) {
        goto L2;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 continue                                                           >*/
L2:
    ;
    }
/*<     3 if(ig.eq.1.and.x(j).ne.flg) td(2,ip)=fc(k)                         >*/
L3:
    if (ig == 1 && x[j] != *flg) {
    td[(ip << 1) + 2] = fc[k];
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     4 return                                                             >*/
L4:
    return 0;
/*<       end                                                                >*/
} /* std_ */

/*<       subroutine reducl (flg,x,nk,az,tb,cm,bz,td,r,azn,tbn,bzn,sc)       >*/
/* Subroutine */ int reducl_(flg, x, nk, az, tb, cm, bz, td, r__, azn, tbn, 
    bzn, sc)
real *flg, *x;
integer *nk;
real *az, *tb, *cm, *bz, *td, *r__, *azn, *tbn, *bzn, *sc;
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    extern integer icat_();
    static integer i__, j, k, m;
    static real t, u;
    extern integer match_();
    static integer ip, iq, no, nv;

/*<       real x(*),tb(5,nk),cm(*),td(2,nk),r(*),tbn(5,nk),sc(*)             >*/
/*<       azn=az                                                             >*/
    /* Parameter adjustments */
    --x;
    tbn -= 6;
    td -= 3;
    tb -= 6;
    --cm;
    --r__;
    --sc;

    /* Function Body */
    *azn = *az;
/*<       do 2 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       do 1 i=1,5                                                         >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       tbn(i,m)=tb(i,m)                                                   >*/
        tbn[i__ + m * 5] = tb[i__ + m * 5];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       bzn=bz                                                             >*/
    *bzn = *bz;
/*<       do 9 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       t=tb(2,m)                                                          >*/
    t = tb[m * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(x(j).eq.flg) go to 9                                            >*/
    if (x[j] == *flg) {
        goto L9;
    }
/*<       if(cm(2*j) .le. 0.0) go to 7                                       >*/
    if (cm[j * 2] <= (float)0.) {
        goto L7;
    }
/*<       k=icat(x(j),j,cm)                                                  >*/
    k = icat_(&x[j], &j, &cm[1]);
/*<       if(k .ne. 0) go to 3                                               >*/
    if (k != 0) {
        goto L3;
    }
/*<       u=0.0                                                              >*/
    u = (float)0.;
/*<       go to 4                                                            >*/
    goto L4;
/*<     3 u=cm(k+int(tb(3,m)+.1))                                            >*/
L3:
    u = cm[k + (integer) (tb[m * 5 + 3] + (float).1)];
/*<     4 if(t .ge. 0.0) go to 6                                             >*/
L4:
    if (t >= (float)0.) {
        goto L6;
    }
/*<       if(u .ne. 0.0) go to 5                                             >*/
    if (u != (float)0.) {
        goto L5;
    }
/*<       u=1.0                                                              >*/
    u = (float)1.;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 u=0.0                                                              >*/
L5:
    u = (float)0.;
/*<     6 td(2,m)=u                                                          >*/
L6:
    td[(m << 1) + 2] = u;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 u=amax1(0.0,sign(1.0,t)*(x(j)-tb(3,m)))                            >*/
L7:
/* Computing MAX */
    r__1 = (float)0., r__2 = r_sign(&c_b196, &t) * (x[j] - tb[m * 5 + 3]);
    u = dmax(r__1,r__2);
/*<     8 sc(m)=u                                                            >*/
L8:
    sc[m] = u;
/*<     9 continue                                                           >*/
L9:
    ;
    }
/*<       m=nk                                                               >*/
    m = *nk;
/*<       go to 11                                                           >*/
    goto L11;
/*<    10 m=m+(-1)                                                           >*/
L10:
    --m;
/*<    11 if((-1)*((m)-(1)).gt.0) go to 21                                   >*/
L11:
    if (-(m - 1) > 0) {
    goto L21;
    }
/*<       ip=tbn(4,m)+.1                                                     >*/
    ip = tbn[m * 5 + 4] + (float).1;
/*<       t=tbn(2,m)                                                         >*/
    t = tbn[m * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(x(j) .ne. flg) go to 15                                         >*/
    if (x[j] != *flg) {
    goto L15;
    }
/*<       if(tbn(1,m) .eq. 0.0) go to 10                                     >*/
    if (tbn[m * 5 + 1] == (float)0.) {
    goto L10;
    }
/*<       iq=ip                                                              >*/
    iq = ip;
/*<    12 if(iq.le.0) go to 10                                               >*/
L12:
    if (iq <= 0) {
    goto L10;
    }
/*<       t=tbn(2,iq)                                                        >*/
    t = tbn[iq * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(x(j) .eq. flg) go to 13                                         >*/
    if (x[j] == *flg) {
    goto L13;
    }
/*<       tbn(1,m)=tbn(1,m)*sc(iq)                                           >*/
    tbn[m * 5 + 1] *= sc[iq];
/*<       td(1,m)=td(1,m)*td(2,iq)                                           >*/
    td[(m << 1) + 1] *= td[(iq << 1) + 2];
/*<    13 iq=tbn(4,iq)+.1                                                    >*/
L13:
    iq = tbn[iq * 5 + 4] + (float).1;
/*<       go to 12                                                           >*/
    goto L12;
/*<    15 k=m+1                                                              >*/
L15:
    k = m + 1;
/*<       go to 17                                                           >*/
    goto L17;
/*<    16 k=k+1                                                              >*/
L16:
    ++k;
/*<    17 if((k).gt.(nk)) go to 18                                           >*/
L17:
    if (k > *nk) {
    goto L18;
    }
/*<       if(int(tbn(4,k)+.1).eq.m) tbn(4,k)=tbn(4,m)                        >*/
    if ((integer) (tbn[k * 5 + 4] + (float).1) == m) {
    tbn[k * 5 + 4] = tbn[m * 5 + 4];
    }
/*<       go to 16                                                           >*/
    goto L16;
/*<    18 if(tbn(1,m).eq.0.0) go to 10                                       >*/
L18:
    if (tbn[m * 5 + 1] == (float)0.) {
    goto L10;
    }
/*<       if(ip .ne. 0) go to 19                                             >*/
    if (ip != 0) {
    goto L19;
    }
/*<       azn=azn+tbn(1,m)*sc(m)                                             >*/
    *azn += tbn[m * 5 + 1] * sc[m];
/*<       bzn=bzn+td(1,m)*td(2,m)                                            >*/
    *bzn += td[(m << 1) + 1] * td[(m << 1) + 2];
/*<       go to 10                                                           >*/
    goto L10;
/*<    19 tbn(1,ip)=tbn(1,ip)+tbn(1,m)*sc(m)                                 >*/
L19:
    tbn[ip * 5 + 1] += tbn[m * 5 + 1] * sc[m];
/*<       td(1,ip)=td(1,ip)+td(1,m)*td(2,m)                                  >*/
    td[(ip << 1) + 1] += td[(m << 1) + 1] * td[(m << 1) + 2];
/*<       go to 10                                                           >*/
    goto L10;
/*<    21 no=nk                                                              >*/
L21:
    no = *nk;
/*<       m=nk                                                               >*/
    m = *nk;
/*<       go to 23                                                           >*/
    goto L23;
/*<    22 m=m+(-1)                                                           >*/
L22:
    --m;
/*<    23 if((-1)*((m)-(1)).gt.0) go to 31                                   >*/
L23:
    if (-(m - 1) > 0) {
    goto L31;
    }
/*<       t=tb(2,m)                                                          >*/
    t = tb[m * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(x(j).eq.flg) go to 22                                           >*/
    if (x[j] == *flg) {
    goto L22;
    }
/*<       k=m+1                                                              >*/
    k = m + 1;
/*<       go to 25                                                           >*/
    goto L25;
/*<    24 k=k+1                                                              >*/
L24:
    ++k;
/*<    25 if((k).gt.(no)) go to 30                                           >*/
L25:
    if (k > no) {
    goto L30;
    }
/*<       td(1,k-1)=td(1,k)                                                  >*/
    td[(k - 1 << 1) + 1] = td[(k << 1) + 1];
/*<       do 26 i=1,5                                                        >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<       tbn(i,k-1)=tbn(i,k)                                                >*/
    tbn[i__ + (k - 1) * 5] = tbn[i__ + k * 5];
/*<    26 continue                                                           >*/
/* L26: */
    }
/*<       i=k+1                                                              >*/
    i__ = k + 1;
/*<       go to 28                                                           >*/
    goto L28;
/*<    27 i=i+1                                                              >*/
L27:
    ++i__;
/*<    28 if((i).gt.(no)) go to 24                                           >*/
L28:
    if (i__ > no) {
    goto L24;
    }
/*<       if(int(tbn(4,i)+.1).eq.k) tbn(4,i)=k-1                             >*/
    if ((integer) (tbn[i__ * 5 + 4] + (float).1) == k) {
    tbn[i__ * 5 + 4] = (real) (k - 1);
    }
/*<       go to 27                                                           >*/
    goto L27;
/*<    30 no=no-1                                                            >*/
L30:
    --no;
/*<       go to 22                                                           >*/
    goto L22;
/*<    31 m=no+1                                                             >*/
L31:
    m = no + 1;
/*<       go to 33                                                           >*/
    goto L33;
/*<    32 m=m+1                                                              >*/
L32:
    ++m;
/*<    33 if((m).gt.(nk)) go to 34                                           >*/
L33:
    if (m > *nk) {
    goto L34;
    }
/*<       tbn(1,m)=0.0                                                       >*/
    tbn[m * 5 + 1] = (float)0.;
/*<       go to 32                                                           >*/
    goto L32;
/*<    34 m=no                                                               >*/
L34:
    m = no;
/*<       go to 36                                                           >*/
    goto L36;
/*<    35 m=m+(-1)                                                           >*/
L35:
    --m;
/*<    36 if((-1)*((m)-(2)).gt.0) go to 39                                   >*/
L36:
    if (-(m - 2) > 0) {
    goto L39;
    }
/*<       if(tbn(1,m).eq.0.0) go to 35                                       >*/
    if (tbn[m * 5 + 1] == (float)0.) {
    goto L35;
    }
/*<       nv=0                                                               >*/
    nv = 0;
/*<       ip=m                                                               >*/
    ip = m;
/*<    37 if(ip.le.0) go to 38                                               >*/
L37:
    if (ip <= 0) {
    goto L38;
    }
/*<       nv=nv+1                                                            >*/
    ++nv;
/*<       sc(2*nv-1)=tbn(2,ip)                                               >*/
    sc[(nv << 1) - 1] = tbn[ip * 5 + 2];
/*<       sc(2*nv)=tbn(3,ip)                                                 >*/
    sc[nv * 2] = tbn[ip * 5 + 3];
/*<       ip=tbn(4,ip)+.1                                                    >*/
    ip = tbn[ip * 5 + 4] + (float).1;
/*<       go to 37                                                           >*/
    goto L37;
/*<    38 k=match(nv,sc,m-1,tbn,cm,r,1)                                      >*/
L38:
    i__1 = m - 1;
    k = match_(&nv, &sc[1], &i__1, &tbn[6], &cm[1], &r__[1], &c__1);
/*<       if(k.eq.0) go to 35                                                >*/
    if (k == 0) {
    goto L35;
    }
/*<       tbn(1,k)=tbn(1,k)+tbn(1,m)                                         >*/
    tbn[k * 5 + 1] += tbn[m * 5 + 1];
/*<       td(1,k)=td(1,k)+td(1,m)                                            >*/
    td[(k << 1) + 1] += td[(m << 1) + 1];
/*<       tbn(1,m)=0.0                                                       >*/
    tbn[m * 5 + 1] = (float)0.;
/*<       go to 35                                                           >*/
    goto L35;
/*<    39 return                                                             >*/
L39:
    return 0;
/*<       end                                                                >*/
} /* reducl_ */

/*<       subroutine qslice (p,nk,tb,cm,td,kp,kv,lp,lv,tc,r,sc,js)           >*/
/* Subroutine */ int qslice_(p, nk, tb, cm, td, kp, kv, lp, lv, tc, r__, sc, 
    js)
integer *p, *nk;
real *tb, *cm, *td;
integer *kp, *kv, *lp, *lv;
real *tc, *r__, *sc;
integer *js;
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int side_(), gtrm_(), knts_();
    static integer j, k, l, m;
    extern integer match_();
    static integer l1, la, le, il, ll, jl, jp, nt, nv, kp3, laa;
    static real dum;

/*<       integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*)                      >*/
/*<       real tb(5,nk),cm(*),td(2,*),tc(*),r(p,2),sc(2,p)                   >*/
/*<       do 1 j=1,p                                                         >*/
    /* Parameter adjustments */
    sc -= 3;
    r_dim1 = *p;
    r_offset = r_dim1 + 1;
    r__ -= r_offset;
    tb -= 6;
    --cm;
    td -= 3;
    kp -= 6;
    kv -= 3;
    lp -= 4;
    --lv;
    --tc;
    --js;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       sc(1,j)=r(j,2)                                                     >*/
    sc[(j << 1) + 1] = r__[j + (r_dim1 << 1)];
/*<       sc(2,j)=sc(1,j)+r(j,1)                                             >*/
    sc[(j << 1) + 2] = sc[(j << 1) + 1] + r__[j + r_dim1];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       ll=1                                                               >*/
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<     2 if(kp(1,ll).lt.0) go to 5                                          >*/
L2:
    if (kp[ll * 5 + 1] < 0) {
    goto L5;
    }
/*<       if(kp(3,ll) .gt. 0) go to 3                                        >*/
    if (kp[ll * 5 + 3] > 0) {
    goto L3;
    }
/*<       kp(5,ll)=0                                                         >*/
    kp[ll * 5 + 5] = 0;
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     3 kp3=kp(3,ll)                                                       >*/
L3:
    kp3 = kp[ll * 5 + 3];
/*<       kp(5,ll)=la                                                        >*/
    kp[ll * 5 + 5] = la;
/*<       do 4 m=1,kp3                                                       >*/
    i__1 = kp3;
    for (m = 1; m <= i__1; ++m) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<    >*/
    knts_(&l, &nt, &lv[lp[l1 * 3 + 2]], &kp[ll * 5 + 1], &kv[(kp[ll * 5 + 
        2] << 1) + 1], nk, &tb[6], &cm[1], &tc[la], &js[1]);
/*<       call side(l,nt,lv(lp(2,l1)),sc,tc(la))                             >*/
    side_(&l, &nt, &lv[lp[l1 * 3 + 2]], &sc[3], &tc[la]);
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       la=la+nt*(5*l+1)                                                   >*/
    la += nt * (l * 5 + 1);
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 2                                                            >*/
    goto L2;
/*<     5 le=la-1                                                            >*/
L5:
    le = la - 1;
/*<       ll=1                                                               >*/
    ll = 1;
/*<       la=ll                                                              >*/
    la = ll;
/*<       l1=la                                                              >*/
    l1 = la;
/*<       laa=0                                                              >*/
    laa = 0;
/*<     6 if(kp(1,ll).lt.0) go to 13                                         >*/
L6:
    if (kp[ll * 5 + 1] < 0) {
    goto L13;
    }
/*<       nv=0                                                               >*/
    nv = 0;
/*<       if(kp(1,ll) .le. 0) go to 8                                        >*/
    if (kp[ll * 5 + 1] <= 0) {
    goto L8;
    }
/*<       jl=kp(1,ll)                                                        >*/
    jl = kp[ll * 5 + 1];
/*<       do 7 il=1,jl                                                       >*/
    i__1 = jl;
    for (il = 1; il <= i__1; ++il) {
/*<       k=kp(2,ll)+il-1                                                    >*/
    k = kp[ll * 5 + 2] + il - 1;
/*<       nv=nv+1                                                            >*/
    ++nv;
/*<       sc(1,nv)=kv(1,k)                                                   >*/
    sc[(nv << 1) + 1] = (real) kv[(k << 1) + 1];
/*<       sc(2,nv)=kv(2,k)                                                   >*/
    sc[(nv << 1) + 2] = (real) kv[(k << 1) + 2];
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<       go to 9                                                            >*/
    goto L9;
/*<     8 if(kp(3,ll) .gt. 0) go to 9                                        >*/
L8:
    if (kp[ll * 5 + 3] > 0) {
    goto L9;
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 6                                                            >*/
    goto L6;
/*<     9 if(kp(3,ll) .gt. 0) go to 10                                       >*/
L9:
    if (kp[ll * 5 + 3] > 0) {
    goto L10;
    }
/*<       m=match(nv,sc,nk,tb,cm,r,0)                                        >*/
    m = match_(&nv, &sc[3], nk, &tb[6], &cm[1], &r__[r_offset], &c__0);
/*<       le=le+1                                                            >*/
    ++le;
/*<       kp(3,ll)=-le                                                       >*/
    kp[ll * 5 + 3] = -le;
/*<       tc(le)=td(1,m)                                                     >*/
    tc[le] = td[(m << 1) + 1];
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 6                                                            >*/
    goto L6;
/*<    10 kp3=kp(3,ll)                                                       >*/
L10:
    kp3 = kp[ll * 5 + 3];
/*<       do 12 k=1,kp3                                                      >*/
    i__1 = kp3;
    for (k = 1; k <= i__1; ++k) {
/*<       l=lp(1,l1)                                                         >*/
    l = lp[l1 * 3 + 1];
/*<       nt=lp(3,l1)                                                        >*/
    nt = lp[l1 * 3 + 3];
/*<       laa=laa+5*l*nt                                                     >*/
    laa += l * 5 * nt;
/*<       do 11 jp=1,nt                                                      >*/
    i__2 = nt;
    for (jp = 1; jp <= i__2; ++jp) {
/*<    >*/
        gtrm_(&c__2, &jp, &l, &nt, &lv[lp[l1 * 3 + 2]], &dum, &dum, nk, &
            tb[6], &tc[la], &sc[(nv + 1 << 1) + 1], &dum);
/*<       m=match(nv+l,sc,nk,tb,cm,r,0)                                      >*/
        i__3 = nv + l;
        m = match_(&i__3, &sc[3], nk, &tb[6], &cm[1], &r__[r_offset], &
            c__0);
/*<       tc(jp+laa)=td(1,m)                                                 >*/
        tc[jp + laa] = td[(m << 1) + 1];
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       laa=laa+nt                                                         >*/
    laa += nt;
/*<       l1=l1+1                                                            >*/
    ++l1;
/*<       la=la+nt*(5*l+1)                                                   >*/
    la += nt * (l * 5 + 1);
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<       ll=ll+1                                                            >*/
    ++ll;
/*<       go to 6                                                            >*/
    goto L6;
/*<    13 return                                                             >*/
L13:
    return 0;
/*<       end                                                                >*/
} /* qslice_ */

/*<       function ieq(a,b,r)                                                >*/
integer ieq_(a, b, r__)
real *a, *b, *r__;
{
    /* System generated locals */
    integer ret_val;
    real r__1;

/*<       ieq=0                                                              >*/
    ret_val = 0;
/*<       if(abs((a-b)/r).lt.1.e-5) ieq=1                                    >*/
    if ((r__1 = (*a - *b) / *r__, dabs(r__1)) < (float)1e-5) {
    ret_val = 1;
    }
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* ieq_ */

/*<       function lcm (p,nk,tb,cm)                                          >*/
integer lcm_(p, nk, tb, cm)
integer *p, *nk;
real *tb, *cm;
{
    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static integer j, m, jj, ix;

/*<       integer p                                                          >*/
/*<       real tb(5,nk),cm(*)                                                >*/
/*<       ix=0                                                               >*/
    /* Parameter adjustments */
    tb -= 6;
    --cm;

    /* Function Body */
    ix = 0;
/*<       do 1 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       j=abs(tb(2,m))+.1                                                  >*/
    j = (r__1 = tb[m * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(cm(2*j).eq.0.0) go to 1                                         >*/
    if (cm[j * 2] == (float)0.) {
        goto L1;
    }
/*<       if(int(tb(3,m)+.1) .le. ix) go to 1                                >*/
    if ((integer) (tb[m * 5 + 3] + (float).1) <= ix) {
        goto L1;
    }
/*<       ix=tb(3,m)+.1                                                      >*/
    ix = tb[m * 5 + 3] + (float).1;
/*<       jj=j                                                               >*/
    jj = j;
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<       if(ix .le. 0) go to 2                                              >*/
    if (ix <= 0) {
    goto L2;
    }
/*<       lcm=ix+int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1                       >*/
    ret_val = ix + (integer) (cm[(jj << 1) + 1] + (float).1) - (integer) (cm[
        jj * 2] + (float).1) + 1;
/*<       return                                                             >*/
    return ret_val;
/*<     2 lcm=2*p+1                                                          >*/
L2:
    ret_val = (*p << 1) + 1;
/*<       do 3 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       if(cm(2*j).eq.0.0) go to 3                                         >*/
    if (cm[j * 2] == (float)0.) {
        goto L3;
    }
/*<       lcm=lcm+int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                        >*/
    ret_val = ret_val + (integer) (cm[(j << 1) + 1] + (float).1) - (
        integer) (cm[j * 2] + (float).1) + 1;
/*<     3 continue                                                           >*/
L3:
    ;
    }
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* lcm_ */

/*<       function newb (m,tb)                                               >*/
integer newb_(m, tb)
integer *m;
real *tb;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer k, mm1;
    extern integer ieq_();

/*<       real tb(5,m)                                                       >*/
/*<       newb=0                                                             >*/
    /* Parameter adjustments */
    tb -= 6;

    /* Function Body */
    ret_val = 0;
/*<       mm1=m-1                                                            >*/
    mm1 = *m - 1;
/*<       do 1 k=1,mm1                                                       >*/
    i__1 = mm1;
    for (k = 1; k <= i__1; ++k) {
/*<       if(ieq(tb(2,k),tb(2,m),1.0).eq.0) go to 1                          >*/
    if (ieq_(&tb[k * 5 + 2], &tb[*m * 5 + 2], &c_b196) == 0) {
        goto L1;
    }
/*<       if(ieq(tb(3,k),tb(3,m),1.0).eq.0) go to 1                          >*/
    if (ieq_(&tb[k * 5 + 3], &tb[*m * 5 + 3], &c_b196) == 0) {
        goto L1;
    }
/*<       if(ieq(tb(4,k),tb(4,m),1.0).eq.0) go to 1                          >*/
    if (ieq_(&tb[k * 5 + 4], &tb[*m * 5 + 4], &c_b196) == 0) {
        goto L1;
    }
/*<       newb=1                                                             >*/
    ret_val = 1;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* newb_ */

/*<       subroutine sscp (n,m,sc,y,w,yb,yv,sw,d,da)                         >*/
/* Subroutine */ int sscp_(n, m, sc, y, w, yb, yv, sw, d__, da)
integer *n, *m;
real *sc, *y, *w;
doublereal *yb, *yv, *sw, *d__, *da;
{
    /* System generated locals */
    integer sc_dim1, sc_offset, d_dim1, d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static integer mm1;

/*<       real sc(n,*),y(n),w(n)                                             >*/
/*<       double precision d(m,m),da(*),yb,yv,sw,s                           >*/
/*<       mm1=m-1                                                            >*/
    /* Parameter adjustments */
    --w;
    --y;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    d_dim1 = *m;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    --da;

    /* Function Body */
    mm1 = *m - 1;
/*<       do 6 k=1,mm1                                                       >*/
    i__1 = mm1;
    for (k = 1; k <= i__1; ++k) {
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 1 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+w(i)*sc(i,k)                                                   >*/
        s += w[i__] * sc[i__ + k * sc_dim1];
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       s=s/sw                                                             >*/
    s /= *sw;
/*<       da(k)=s                                                            >*/
    da[k] = s;
/*<       do 2 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       sc(i,k)=sc(i,k)-s                                                  >*/
        sc[i__ + k * sc_dim1] -= s;
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       do 4 j=1,k                                                         >*/
    i__2 = k;
    for (j = 1; j <= i__2; ++j) {
/*<       s=0.d0                                                             >*/
        s = 0.;
/*<       do 3 i=1,n                                                         >*/
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<       s=s+w(i)*sc(i,j)*sc(i,k)                                           >*/
        s += w[i__] * sc[i__ + j * sc_dim1] * sc[i__ + k * sc_dim1];
/*<     3 continue                                                           >*/
/* L3: */
        }
/*<       d(j,k)=s                                                           >*/
        d__[j + k * d_dim1] = s;
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 5 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       s=s+w(i)*sc(i,k)*(y(i)-yb)                                         >*/
        s += w[i__] * sc[i__ + k * sc_dim1] * (y[i__] - *yb);
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<       d(k,m)=s                                                           >*/
    d__[k + *m * d_dim1] = s;
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       d(m,m)=sw*yv                                                       >*/
    d__[*m + *m * d_dim1] = *sw * *yv;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* sscp_ */

/*<       subroutine lsf1 (d,m,xb,yb,al,rss,a,a0,dp)                         >*/
/* Subroutine */ int lsf1_(d__, m, xb, yb, al, rss, a, a0, dp)
doublereal *d__;
integer *m;
doublereal *xb, *yb, *al, *rss, *a, *a0, *dp;
{
    /* Initialized data */

    static doublereal eps = 1e-4;

    /* System generated locals */
    integer d_dim1, d_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    extern /* Subroutine */ int sweep_();
    static integer im1, mm1;

/*<       double precision d(m,m),xb(*),yb,al,rss,a(*),a0,dp(*),eps,s        >*/
/*<       data eps /1.d-4/                                                   >*/
    /* Parameter adjustments */
    d_dim1 = *m;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    --xb;
    --a;
    --dp;

    /* Function Body */
/*<       mm1=m-1                                                            >*/
    mm1 = *m - 1;
/*<       do 1 i=1,mm1                                                       >*/
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       dp(i)=d(i,i)                                                       >*/
    dp[i__] = d__[i__ + i__ * d_dim1];
/*<       d(i,i)=d(i,i)*(1.d0+al)                                            >*/
    d__[i__ + i__ * d_dim1] *= *al + 1.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       do 5 i=1,mm1                                                       >*/
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(dp(i).le.0.d0) go to 5                                          >*/
    if (dp[i__] <= 0.) {
        goto L5;
    }
/*<       im1=i-1                                                            >*/
    im1 = i__ - 1;
/*<       s=dp(i)                                                            >*/
    s = dp[i__];
/*<       j=1                                                                >*/
    j = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 j=j+1                                                              >*/
L2:
    ++j;
/*<     3 if((j).gt.(im1)) go to 4                                           >*/
L3:
    if (j > im1) {
        goto L4;
    }
/*<       if(d(j,j).lt.0.d0) s=s+dp(j)*d(j,i)**2                             >*/
    if (d__[j + j * d_dim1] < 0.) {
/* Computing 2nd power */
        d__1 = d__[j + i__ * d_dim1];
        s += dp[j] * (d__1 * d__1);
    }
/*<       go to 2                                                            >*/
    goto L2;
/*<     4 if((d(i,i)-al*s)/dp(i).lt.eps) go to 5                             >*/
L4:
    if ((d__[i__ + i__ * d_dim1] - *al * s) / dp[i__] < eps) {
        goto L5;
    }
/*<       call sweep(d,m,i,-1.d0,dp(m))                                      >*/
    sweep_(&d__[d_offset], m, &i__, &c_b214, &dp[*m]);
/*<     5 continue                                                           >*/
L5:
    ;
    }
/*<       rss=0.d0                                                           >*/
    *rss = 0.;
/*<       a0=yb                                                              >*/
    *a0 = *yb;
/*<       do 6 i=1,mm1                                                       >*/
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       a(i)=0.d0                                                          >*/
    a[i__] = 0.;
/*<       if(d(i,i).ge.0.d0) go to 6                                         >*/
    if (d__[i__ + i__ * d_dim1] >= 0.) {
        goto L6;
    }
/*<       a(i)=d(i,m)                                                        >*/
    a[i__] = d__[i__ + *m * d_dim1];
/*<       a0=a0-a(i)*xb(i)                                                   >*/
    *a0 -= a[i__] * xb[i__];
/*<       rss=rss+dp(i)*a(i)**2                                              >*/
/* Computing 2nd power */
    d__1 = a[i__];
    *rss += dp[i__] * (d__1 * d__1);
/*<     6 continue                                                           >*/
L6:
    ;
    }
/*<       rss=d(m,m)-al*rss                                                  >*/
    *rss = d__[*m + *m * d_dim1] - *al * *rss;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* lsf1_ */

/*<       subroutine bkstp (d,m,xb,yb,al,rss,a,a0,k,dp)                      >*/
/* Subroutine */ int bkstp_(d__, m, xb, yb, al, rss, a, a0, k, dp)
doublereal *d__;
integer *m;
doublereal *xb, *yb, *al, *rss, *a, *a0;
integer *k;
doublereal *dp;
{
    /* Initialized data */

    static real big = (float)9.9e30;

    /* System generated locals */
    integer d_dim1, d_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    extern /* Subroutine */ int sweep_();
    static integer mm1;

/*<       double precision d(m,m),xb(*),yb,al                                >*/
/*<       double precision a(*),a0,rss,dp(*),s                               >*/
/*<       data big /9.9e30/                                                  >*/
    /* Parameter adjustments */
    d_dim1 = *m;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    --xb;
    --a;
    --dp;

    /* Function Body */
/*<       mm1=m-1                                                            >*/
    mm1 = *m - 1;
/*<       rss=big                                                            >*/
    *rss = big;
/*<       k=0                                                                >*/
    *k = 0;
/*<       do 4 i=1,mm1                                                       >*/
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(d(i,i).ge.0.d0) go to 4                                         >*/
    if (d__[i__ + i__ * d_dim1] >= 0.) {
        goto L4;
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       do 3 j=1,mm1                                                       >*/
    i__2 = mm1;
    for (j = 1; j <= i__2; ++j) {
/*<       if(d(j,j).ge.0.d0) go to 3                                         >*/
        if (d__[j + j * d_dim1] >= 0.) {
        goto L3;
        }
/*<       if(j.eq.i) go to 3                                                 >*/
        if (j == i__) {
        goto L3;
        }
/*<       if(j .ge. i) go to 1                                               >*/
        if (j >= i__) {
        goto L1;
        }
/*<       a0=d(j,i)                                                          >*/
        *a0 = d__[j + i__ * d_dim1];
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 a0=d(i,j)                                                          >*/
L1:
        *a0 = d__[i__ + j * d_dim1];
/*<     2 s=s+dp(j)*(d(j,m)-a0*d(i,m)/d(i,i))**2                             >*/
L2:
/* Computing 2nd power */
        d__1 = d__[j + *m * d_dim1] - *a0 * d__[i__ + *m * d_dim1] / d__[
            i__ + i__ * d_dim1];
        s += dp[j] * (d__1 * d__1);
/*<     3 continue                                                           >*/
L3:
        ;
    }
/*<       s=d(m,m)-d(i,m)**2/d(i,i)-al*s                                     >*/
/* Computing 2nd power */
    d__1 = d__[i__ + *m * d_dim1];
    s = d__[*m + *m * d_dim1] - d__1 * d__1 / d__[i__ + i__ * d_dim1] - *
        al * s;
/*<       if(s .gt. rss) go to 4                                             >*/
    if (s > *rss) {
        goto L4;
    }
/*<       rss=s                                                              >*/
    *rss = s;
/*<       k=i                                                                >*/
    *k = i__;
/*<     4 continue                                                           >*/
L4:
    ;
    }
/*<       if(k.gt.0) call sweep(d,m,k,1.d0,dp(m))                            >*/
    if (*k > 0) {
    sweep_(&d__[d_offset], m, k, &c_b1099, &dp[*m]);
    }
/*<       a0=yb                                                              >*/
    *a0 = *yb;
/*<       rss=0.d0                                                           >*/
    *rss = 0.;
/*<       do 5 i=1,mm1                                                       >*/
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       a(i)=0.d0                                                          >*/
    a[i__] = 0.;
/*<       if(d(i,i).ge.0.d0) go to 5                                         >*/
    if (d__[i__ + i__ * d_dim1] >= 0.) {
        goto L5;
    }
/*<       a(i)=d(i,m)                                                        >*/
    a[i__] = d__[i__ + *m * d_dim1];
/*<       a0=a0-a(i)*xb(i)                                                   >*/
    *a0 -= a[i__] * xb[i__];
/*<       rss=rss+dp(i)*a(i)**2                                              >*/
/* Computing 2nd power */
    d__1 = a[i__];
    *rss += dp[i__] * (d__1 * d__1);
/*<     5 continue                                                           >*/
L5:
    ;
    }
/*<       rss=d(m,m)-al*rss                                                  >*/
    *rss = d__[*m + *m * d_dim1] - *al * *rss;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* bkstp_ */

/*<       subroutine sweep (a,m,k,fl,u)                                      >*/
/* Subroutine */ int sweep_(a, m, k, fl, u)
doublereal *a;
integer *m, *k;
doublereal *fl, *u;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;

/*<       double precision a(m,m),u(m),fl,c                                  >*/
/*<       c=a(k,k)                                                           >*/
    /* Parameter adjustments */
    --u;
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    c__ = a[*k + *k * a_dim1];
/*<       do 1 i=1,k                                                         >*/
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       u(i)=a(i,k)                                                        >*/
    u[i__] = a[i__ + *k * a_dim1];
/*<       a(i,k)=0.d0                                                        >*/
    a[i__ + *k * a_dim1] = 0.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       do 2 i=k,m                                                         >*/
    i__1 = *m;
    for (i__ = *k; i__ <= i__1; ++i__) {
/*<       u(i)=a(k,i)                                                        >*/
    u[i__] = a[*k + i__ * a_dim1];
/*<       a(k,i)=0.d0                                                        >*/
    a[*k + i__ * a_dim1] = 0.;
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<       u(k)=fl                                                            >*/
    u[*k] = *fl;
/*<       do 4 i=1,m                                                         >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       do 3 j=i,m                                                         >*/
    i__2 = *m;
    for (j = i__; j <= i__2; ++j) {
/*<       a(i,j)=a(i,j)-u(i)*u(j)/c                                          >*/
        a[i__ + j * a_dim1] -= u[i__] * u[j] / c__;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* sweep_ */

/*<       subroutine array(p,n,i,j)                                          >*/
/* Subroutine */ int array_(p, n, i__, j)
integer *p, *n, *i__, *j;
{
/*<       integer p                                                          >*/
/*<       i=mod(p,n)                                                         >*/
    *i__ = *p % *n;
/*<       if(i.eq.0) i=n                                                     >*/
    if (*i__ == 0) {
    *i__ = *n;
    }
/*<       j=(p-i)/n+1                                                        >*/
    *j = (*p - *i__) / *n + 1;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* array_ */

/*<    >*/
/* Subroutine */ int cvmars_0_(n__, ix, n, p, x, y, w, nk, ms, df, fv, mi, lx,
     it, xm, xs, tb, cm, sc, db, d__, mm, wt, cv, a1, a2, ia1)
int n__;
integer *ix, *n, *p;
real *x, *y, *w;
integer *nk, *ms;
real *df, *fv;
integer *mi, *lx, *it;
real *xm, *xs, *tb, *cm, *sc;
doublereal *db, *d__;
integer *mm;
real *wt, *cv, *a1, *a2;
integer *ia1;
{
    /* Initialized data */

    static real eps = (float)1e-6;
    static real big = (float)9.9e30;
    static real dfs = (float)0.;
    static real cvm = (float)0.;
    static integer im = 0;

    /* System generated locals */
    integer mm_dim1, mm_offset, x_dim1, x_offset, wt_dim1, wt_offset, cv_dim1,
         cv_offset, db_dim1, db_offset, d_dim1, d_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int rnms_();
    static integer i__, k, m;
    static real r__, t;
    extern /* Subroutine */ int cvmod_();
    static real fc, am;
    static integer nd, ir, nr;
    static real az, wn;
    static integer mk;
    static real sw, yv;
    extern /* Subroutine */ int marsgo_();
    static real am1, cv0, wn1, yv1, dfu, gcv, cvl, dmx;

/*<       integer p,mm(n,*),lx(p)                                            >*/
/*<    >*/
/*<       double precision db(n,*),d(nk,*)                                   >*/
/*<       data eps,big,dfs,cvm,im /1.e-6,9.9e30,2*0.0,0/                     >*/
    /* Parameter adjustments */
    if (y) {
    --y;
    }
    if (w) {
    --w;
    }
    if (db) {
    db_dim1 = *n;
    db_offset = db_dim1 + 1;
    db -= db_offset;
    }
    if (mm) {
    mm_dim1 = *n;
    mm_offset = mm_dim1 + 1;
    mm -= mm_offset;
    }
    if (wt) {
    wt_dim1 = *n;
    wt_offset = wt_dim1 + 1;
    wt -= wt_offset;
    }
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (lx) {
    --lx;
    }
    if (xm) {
    --xm;
    }
    if (xs) {
    --xs;
    }
    if (tb) {
    tb -= 6;
    }
    if (d__) {
    d_dim1 = *nk;
    d_offset = d_dim1 + 1;
    d__ -= d_offset;
    }
    if (cv) {
    cv_dim1 = *nk;
    cv_offset = cv_dim1 + 1;
    cv -= cv_offset;
    }
    if (cm) {
    --cm;
    }
    if (sc) {
    --sc;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_cvinfo;
    }

/*    if(it.gt.0) write(it,'(/,'' sample reuse to estimate df:'')')      4
296*/
/*<       if(ix .le. 0) go to 1                                              >*/
    if (*ix <= 0) {
    goto L1;
    }
/*<       nr=ix                                                              >*/
    nr = *ix;
/*<       nd=nr                                                              >*/
    nd = nr;
/*    if(it.gt.0) write(it,'('' '',i3,'' - fold cross-validation.'',/)') 4
300*/
/*   1 ix                                                                4
301*/
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 nr=1                                                               >*/
L1:
    nr = 1;
/*<       nd=-ix                                                             >*/
    nd = -(*ix);
/*    if(it.gt.0) write(it,'('' independent test set - every'',i4,'' obs 4
305*/
/*   1ervations.'',/)') nd                                               4
306*/
/*<     2 do 3 i=1,n                                                         >*/
L2:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       wt(i,1)=w(i)                                                       >*/
    wt[i__ + wt_dim1] = w[i__];
/*<       wt(i,2)=i                                                          >*/
    wt[i__ + (wt_dim1 << 1)] = (real) i__;
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       do 4 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       call rnms(r,1)                                                     >*/
    rnms_(&r__, &c__1);
/*<       k=(n-i+1)*r+i                                                      >*/
    k = (*n - i__ + 1) * r__ + i__;
/*<       t=wt(i,2)                                                          >*/
    t = wt[i__ + (wt_dim1 << 1)];
/*<       wt(i,2)=wt(k,2)                                                    >*/
    wt[i__ + (wt_dim1 << 1)] = wt[k + (wt_dim1 << 1)];
/*<       wt(k,2)=t                                                          >*/
    wt[k + (wt_dim1 << 1)] = t;
/*<     4 continue                                                           >*/
/* L4: */
    }
/*<       do 6 i=1,3                                                         >*/
    for (i__ = 1; i__ <= 3; ++i__) {
/*<       do 5 m=1,nk                                                        >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       cv(m,i)=0.0                                                        >*/
        cv[m + i__ * cv_dim1] = (float)0.;
/*<     5 continue                                                           >*/
/* L5: */
    }
/*<     6 continue                                                           >*/
/* L6: */
    }
/*<       cv0=0.0                                                            >*/
    cv0 = (float)0.;
/*<       sw=cv0                                                             >*/
    sw = cv0;
/*<       wn=sw                                                              >*/
    wn = sw;
/*<       yv=wn                                                              >*/
    yv = wn;
/*<       fc=yv                                                              >*/
    fc = yv;
/*<       do 14 ir=1,nr                                                      >*/
    i__1 = nr;
    for (ir = 1; ir <= i__1; ++ir) {
/*<       i=ir                                                               >*/
    i__ = ir;
/*<     7 if(i.gt.n) go to 8                                                 >*/
L7:
    if (i__ > *n) {
        goto L8;
    }
/*<       wt(int(wt(i,2)+.1),1)=0.0                                          >*/
    wt[(integer) (wt[i__ + (wt_dim1 << 1)] + (float).1) + wt_dim1] = (
        float)0.;
/*<       i=i+nd                                                             >*/
    i__ += nd;
/*<       go to 7                                                            >*/
    goto L7;
/*<    >*/
L8:
    marsgo_(n, p, &x[x_offset], &y[1], &wt[wt_offset], nk, ms, df, fv, mi,
         &lx[1], &c__0, &xm[1], &xs[1], &az, &tb[6], &cm[1], &sc[1], &
        db[db_offset], &d__[d_offset], &mm[mm_offset]);
/*<       yv1=sc(3)                                                          >*/
    yv1 = sc[3];
/*<       yv=yv+yv1                                                          >*/
    yv += yv1;
/*<       wn1=sc(2)                                                          >*/
    wn1 = sc[2];
/*<       wn=wn+wn1                                                          >*/
    wn += wn1;
/*<       fc=fc+sc(1)                                                        >*/
    fc += sc[1];
/*<       mk=sc((nk+1)**2+4)+.1                                              >*/
/* Computing 2nd power */
    i__2 = *nk + 1;
    mk = sc[i__2 * i__2 + 4] + (float).1;
/*<       i=ir                                                               >*/
    i__ = ir;
/*<     9 if(i.gt.n) go to 10                                                >*/
L9:
    if (i__ > *n) {
        goto L10;
    }
/*<       k=wt(i,2)+.1                                                       >*/
    k = wt[i__ + (wt_dim1 << 1)] + (float).1;
/*<       wt(k,1)=w(k)                                                       >*/
    wt[k + wt_dim1] = w[k];
/*<       sw=sw+w(k)                                                         >*/
    sw += w[k];
/*<       call cvmod(k,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv(1,3))                   >*/
    cvmod_(&k, n, &x[x_offset], &y[1], &w[1], nk, &mk, &tb[6], &cm[1], &
        sc[1], &cv0, &cv[cv_dim1 * 3 + 1]);
/*<       i=i+nd                                                             >*/
    i__ += nd;
/*<       go to 9                                                            >*/
    goto L9;
/*<    10 do 13 m=1,nk                                                       >*/
L10:
    i__2 = *nk;
    for (m = 1; m <= i__2; ++m) {
/*<       am=sc(m+4)                                                         >*/
        am = sc[m + 4];
/*<       cv(m,2)=cv(m,2)+am                                                 >*/
        cv[m + (cv_dim1 << 1)] += am;
/*<       am1=yv1                                                            >*/
        am1 = yv1;
/*<       if(m.gt.1) am1=sc(m+3)                                             >*/
        if (m > 1) {
        am1 = sc[m + 3];
        }
/*<       if(am1/yv1 .le. eps) go to 11                                      >*/
        if (am1 / yv1 <= eps) {
        goto L11;
        }
/*<       r=sqrt(am/am1)                                                     >*/
        r__ = sqrt(am / am1);
/*<       go to 12                                                           >*/
        goto L12;
/*<    11 r=1.0                                                              >*/
L11:
        r__ = (float)1.;
/*<    12 cv(m,1)=cv(m,1)+((wn1-1.0)*(1.0-r)/(m-r*(m-1))-1.0)/sc(1)          >*/
L12:
        cv[m + cv_dim1] += ((wn1 - (float)1.) * ((float)1. - r__) / (m - 
            r__ * (m - 1)) - (float)1.) / sc[1];
/*<    13 continue                                                           >*/
/* L13: */
    }
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<       do 15 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       cv(m,1)=cv(m,1)/nr                                                 >*/
    cv[m + cv_dim1] /= nr;
/*<       cv(m,2)=cv(m,2)/nr                                                 >*/
    cv[m + (cv_dim1 << 1)] /= nr;
/*<       cv(m,3)=cv(m,3)/sw                                                 >*/
    cv[m + cv_dim1 * 3] /= sw;
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<       fc=fc/nr                                                           >*/
    fc /= nr;
/*<       yv=yv/nr                                                           >*/
    yv /= nr;
/*<       wn=wn/nr                                                           >*/
    wn /= nr;
/*<       cv0=cv0/sw                                                         >*/
    cv0 /= sw;
/*    if(it.gt.0) write(it,21)                                           4
371*/
/*<       im=0                                                               >*/
    im = 0;
/*<       cvm=cv0                                                            >*/
    cvm = cv0;
/*<       dmx=-big                                                           >*/
    dmx = -big;
/*<       cvl=cv(nk,1)                                                       >*/
    cvl = cv[*nk + cv_dim1];
/*<       m=nk                                                               >*/
    m = *nk;
/*<       go to 17                                                           >*/
    goto L17;
/*<    16 m=m+(-1)                                                           >*/
L16:
    --m;
/*<    17 if((-1)*((m)-(1)).gt.0) go to 19                                   >*/
L17:
    if (-(m - 1) > 0) {
    goto L19;
    }
/*<       if(cv(m,1).le.dmx) go to 16                                        >*/
    if (cv[m + cv_dim1] <= dmx) {
    goto L16;
    }
/*<       dmx=cv(m,1)                                                        >*/
    dmx = cv[m + cv_dim1];
/*<       dfu=0.5*(cvl+cv(m,1))                                              >*/
    dfu = (cvl + cv[m + cv_dim1]) * (float).5;
/*<       cvl=cv(m,1)                                                        >*/
    cvl = cv[m + cv_dim1];
/*<       if(cv(m,3) .gt. cvm) go to 18                                      >*/
    if (cv[m + cv_dim1 * 3] > cvm) {
    goto L18;
    }
/*<       cvm=cv(m,3)                                                        >*/
    cvm = cv[m + cv_dim1 * 3];
/*<       df=dfu                                                             >*/
    *df = dfu;
/*<       im=m                                                               >*/
    im = m;
/*<    18 gcv=cv(m,2)/(1.0-((dfu*fc+1.0)*m+1.0)/wn)**2                       >*/
L18:
/* Computing 2nd power */
    r__1 = (float)1. - ((dfu * fc + (float)1.) * m + (float)1.) / wn;
    gcv = cv[m + (cv_dim1 << 1)] / (r__1 * r__1);
/*    if(it.gt.0) write(it,22) m,dfu,cv(m,2),gcv,cv(m,3)                 4
389*/
/*<       go to 16                                                           >*/
    goto L16;
/*<    19 if(cv0 .gt. cvm) go to 20                                          >*/
L19:
    if (cv0 > cvm) {
    goto L20;
    }
/*<       cvm=cv0                                                            >*/
    cvm = cv0;
/*<       df=dmx                                                             >*/
    *df = dmx;
/*<       im=0                                                               >*/
    im = 0;
/*<    20 dfs=df                                                             >*/
L20:
    dfs = *df;
/*<       gcv=yv/(1.0-1.0/wn)**2                                             >*/
/* Computing 2nd power */
    r__1 = (float)1. - (float)1. / wn;
    gcv = yv / (r__1 * r__1);
/*    if(it.gt.0) write(it,22) 0,dmx,yv,gcv,cv0                          4
397*/
/*    if(it.gt.0) write(it,'(/,'' estimated optimal df('',i3,'') ='',f7. 4
398*/
/*   12,               '' with (estimated) pse ='',g12.4)') im,df,cvm    4
399*/
/*<       return                                                             >*/
    return 0;
/*<       entry cvinfo(a1,a2,ia1)                                            >*/

L_cvinfo:
/*<       a1=dfs                                                             >*/
    *a1 = dfs;
/*<       a2=cvm                                                             >*/
    *a2 = cvm;
/*<       ia1=im                                                             >*/
    *ia1 = im;
/*<       return                                                             >*/
    return 0;
/*<    21 format('  #bsfns     df        asr           gcv           cv')    >*/
/* L21: */
/*<    22 format(' ',i5,f10.2,3g14.4)                                        >*/
/* L22: */
/*<       end                                                                >*/
} /* cvmars_ */

/* Subroutine */ int cvmars_(ix, n, p, x, y, w, nk, ms, df, fv, mi, lx, it, 
    xm, xs, tb, cm, sc, db, d__, mm, wt, cv)
integer *ix, *n, *p;
real *x, *y, *w;
integer *nk, *ms;
real *df, *fv;
integer *mi, *lx, *it;
real *xm, *xs, *tb, *cm, *sc;
doublereal *db, *d__;
integer *mm;
real *wt, *cv;
{
    return cvmars_0_(0, ix, n, p, x, y, w, nk, ms, df, fv, mi, lx, it, xm, xs,
         tb, cm, sc, db, d__, mm, wt, cv, (real *)0, (real *)0, (integer *
        )0);
    }

/* Subroutine */ int cvinfo_(a1, a2, ia1)
real *a1, *a2;
integer *ia1;
{
    return cvmars_0_(1, (integer *)0, (integer *)0, (integer *)0, (real *)0, (
        real *)0, (real *)0, (integer *)0, (integer *)0, (real *)0, (real 
        *)0, (integer *)0, (integer *)0, (integer *)0, (real *)0, (real *)
        0, (real *)0, (real *)0, (real *)0, (doublereal *)0, (doublereal *
        )0, (integer *)0, (real *)0, (real *)0, a1, a2, ia1);
    }

/*<       subroutine cvmod (i,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv)                 >*/
/* Subroutine */ int cvmod_(i__, n, x, y, w, nk, mk, tb, cm, sc, cv0, cv)
integer *i__, *n;
real *x, *y, *w;
integer *nk, *mk;
real *tb, *cm, *sc, *cv0, *cv;
{
    /* System generated locals */
    integer x_dim1, x_offset, cv_dim1, cv_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign();

    /* Local variables */
    static integer j, k, l, m;
    static real s, t, u;
    static integer kp;

/*<       real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(*),cv(nk,2)                >*/
/*<       do 8 m=1,mk                                                        >*/
    /* Parameter adjustments */
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    cv_dim1 = *nk;
    cv_offset = cv_dim1 + 1;
    cv -= cv_offset;
    tb -= 6;
    --cm;
    --sc;

    /* Function Body */
    i__1 = *mk;
    for (m = 1; m <= i__1; ++m) {
/*<       t=tb(2,m)                                                          >*/
    t = tb[m * 5 + 2];
/*<       j=abs(t)+.1                                                        >*/
    j = dabs(t) + (float).1;
/*<       if(cm(2*j) .le. 0.0) go to 5                                       >*/
    if (cm[j * 2] <= (float)0.) {
        goto L5;
    }
/*<       k=x(i,j)+.1                                                        >*/
    k = x[*i__ + j * x_dim1] + (float).1;
/*<       if(k .ne. 0) go to 1                                               >*/
    if (k != 0) {
        goto L1;
    }
/*<       u=0.0                                                              >*/
    u = (float)0.;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 u=cm(k+int(tb(3,m)+.1))                                            >*/
L1:
    u = cm[k + (integer) (tb[m * 5 + 3] + (float).1)];
/*<     2 if(t .ge. 0.0) go to 6                                             >*/
L2:
    if (t >= (float)0.) {
        goto L6;
    }
/*<       if(u .ne. 0.0) go to 3                                             >*/
    if (u != (float)0.) {
        goto L3;
    }
/*<       u=1.0                                                              >*/
    u = (float)1.;
/*<       go to 6                                                            >*/
    goto L6;
/*<     3 u=0.0                                                              >*/
L3:
    u = (float)0.;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,m)))                          >*/
L5:
/* Computing MAX */
    r__1 = (float)0., r__2 = r_sign(&c_b196, &t) * (x[*i__ + j * x_dim1] 
        - tb[m * 5 + 3]);
    u = dmax(r__1,r__2);
/*<     6 l=tb(4,m)+.1                                                       >*/
L6:
    l = tb[m * 5 + 4] + (float).1;
/*<       if(l .le. 0) go to 7                                               >*/
    if (l <= 0) {
        goto L7;
    }
/*<       cv(m,2)=u*cv(l,2)                                                  >*/
    cv[m + (cv_dim1 << 1)] = u * cv[l + (cv_dim1 << 1)];
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 cv(m,2)=u                                                          >*/
L7:
    cv[m + (cv_dim1 << 1)] = u;
/*<     8 continue                                                           >*/
L8:
    ;
    }
/*<       kp=nk+4                                                            >*/
    kp = *nk + 4;
/*<       cv0=cv0+w(i)*(y(i)-sc(4))**2                                       >*/
/* Computing 2nd power */
    r__1 = y[*i__] - sc[4];
    *cv0 += w[*i__] * (r__1 * r__1);
/*<       do 10 m=1,nk                                                       >*/
    i__1 = *nk;
    for (m = 1; m <= i__1; ++m) {
/*<       kp=kp+1                                                            >*/
    ++kp;
/*<       s=sc(kp)                                                           >*/
    s = sc[kp];
/*<       do 9 l=1,nk                                                        >*/
    i__2 = *nk;
    for (l = 1; l <= i__2; ++l) {
/*<       kp=kp+1                                                            >*/
        ++kp;
/*<       if(l.le.mk) s=s+sc(kp)*cv(l,2)                                     >*/
        if (l <= *mk) {
        s += sc[kp] * cv[l + (cv_dim1 << 1)];
        }
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       cv(m,1)=cv(m,1)+w(i)*(y(i)-s)**2                                   >*/
/* Computing 2nd power */
    r__1 = y[*i__] - s;
    cv[m + cv_dim1] += w[*i__] * (r__1 * r__1);
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* cvmod_ */

/*<       subroutine nest (n,i,j,nv,vals)                                    >*/
/* Subroutine */ int nest_0_(n__, n, i__, j, nv, vals, it, p, lx, cm, jb, lm, 
    mk, tb, ja, x, bl)
int n__;
integer *n, *i__, *j, *nv;
real *vals;
integer *it, *p, *lx;
real *cm;
integer *jb, *lm, *mk;
real *tb;
integer *ja;
real *x, *bl;
{
    /* Initialized data */

    static integer il = 0;
    static integer jl = 0;

    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern integer icat_(), nord_();
    static integer norm, k, l, m[800]   /* was [4][200] */;
    static real t;
    extern integer ieqbf_();
    static integer j1, k1, k2;
    static real t1;
    static integer ig, nc, jg, kg, ll, jn, jp, kk, ip, jv, kp, lp;
    static real vm[2000], ex;
    static integer lk, kx, lon;

/*<       parameter (mlist=200, nlist=2000)                                  >*/
/*<       real vals(*),vm(nlist),tb(5,*),cm(*),x(n,*),bl(*)                  >*/
/*<       integer m(4,mlist),p,lx(*)                                         >*/
/*<       save m,vm                                                          >*/
/*<       data il,jl /2*0/                                                   >*/
    /* Parameter adjustments */
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (vals) {
    --vals;
    }
    if (lx) {
    --lx;
    }
    if (cm) {
    --cm;
    }
    if (tb) {
    tb -= 6;
    }
    if (bl) {
    --bl;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_nstlst;
    case 2: goto L_oknest;
    case 3: goto L_isnstr;
    case 4: goto L_isfac;
    case 5: goto L_cmpnst;
    case 6: goto L_getnst;
    }

/*<       if((i .ne. 0) .and. (j .ne. 0)) go to 1                            >*/
    if (*i__ != 0 && *j != 0) {
    goto L1;
    }
/*<       il=0                                                               >*/
    il = 0;
/*<       jl=il                                                              >*/
    jl = il;
/*<       return                                                             >*/
    return 0;
/*<     1 if(i.eq.j) return                                                  >*/
L1:
    if (*i__ == *j) {
    return 0;
    }
/*<       ig=0                                                               >*/
    ig = 0;
/*<       if(nv .le. 0) go to 8                                              >*/
    if (*nv <= 0) {
    goto L8;
    }
/*<       k=1                                                                >*/
    k = 1;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 k=k+1                                                              >*/
L2:
    ++k;
/*<     3 if((k).gt.(il)) go to 4                                            >*/
L3:
    if (k > il) {
    goto L4;
    }
/*<       if(m(1,k).eq.i.or.m(1,k).eq.j) return                              >*/
    if (m[(k << 2) - 4] == *i__ || m[(k << 2) - 4] == *j) {
    return 0;
    }
/*<       go to 2                                                            >*/
    goto L2;
/*<     4 il=il+1                                                            >*/
L4:
    ++il;
/*<       if(il .le. mlist) go to 5                                          >*/
    if (il <= 200) {
    goto L5;
    }
/*    write(6,  '('' increase parameter mlist in subroutine nest to grea 4
467*/
/*   1ter than'',               i5,'' and recompile.'')') il             4
468*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     5 m(1,il)=i                                                          >*/
L5:
    m[(il << 2) - 4] = *i__;
/*<       m(2,il)=j                                                          >*/
    m[(il << 2) - 3] = *j;
/*<       m(3,il)=nv                                                         >*/
    m[(il << 2) - 2] = *nv;
/*<       m(4,il)=jl                                                         >*/
    m[(il << 2) - 1] = jl;
/*<       if(jl+nv .le. nlist) go to 6                                       >*/
    if (jl + *nv <= 2000) {
    goto L6;
    }
/*    write(6,  '('' increase parameter nlist in subroutine nest to grea 4
475*/
/*   1ter than'',               i5,'' and recompile.'')') jl+nv          4
476*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     6 do 7 k=1,nv                                                        >*/
L6:
    i__1 = *nv;
    for (k = 1; k <= i__1; ++k) {
/*<       jl=jl+1                                                            >*/
    ++jl;
/*<       vm(jl)=vals(k)                                                     >*/
    vm[jl - 1] = vals[k];
/*<     7 continue                                                           >*/
/* L7: */
    }
/*<       return                                                             >*/
    return 0;
/*<     8 k=1                                                                >*/
L8:
    k = 1;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 k=k+1                                                              >*/
L9:
    ++k;
/*<    10 if((k).gt.(il)) go to 12                                           >*/
L10:
    if (k > il) {
    goto L12;
    }
/*<       if(m(1,k) .ne. i .or. m(2,k) .ne. j) go to 9                       >*/
    if (m[(k << 2) - 4] != *i__ || m[(k << 2) - 3] != *j) {
    goto L9;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<    12 if(ig.eq.0) return                                                 >*/
L12:
    if (ig == 0) {
    return 0;
    }
/*<       il=il-1                                                            >*/
    --il;
/*<       ll=k                                                               >*/
    ll = k;
/*<       go to 14                                                           >*/
    goto L14;
/*<    13 ll=ll+1                                                            >*/
L13:
    ++ll;
/*<    14 if((ll).gt.(il)) go to 16                                          >*/
L14:
    if (ll > il) {
    goto L16;
    }
/*<       do 15 l=1,4                                                        >*/
    for (l = 1; l <= 4; ++l) {
/*<       m(l,ll)=m(l,ll+1)                                                  >*/
    m[l + (ll << 2) - 5] = m[l + (ll + 1 << 2) - 5];
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<       go to 13                                                           >*/
    goto L13;
/*<    16 return                                                             >*/
L16:
    return 0;
/*<       entry nstlst(it)                                                   >*/

L_nstlst:
/*<       if(it.le.0) return                                                 >*/
    if (*it <= 0) {
    return 0;
    }
/*<       if(il.eq.0) return                                                 >*/
    if (il == 0) {
    return 0;
    }
/*    write(it,'(/,'' variable nesting:'',/)')                           4
503*/
/*<       do 18 k=1,il                                                       >*/
    i__1 = il;
    for (k = 1; k <= i__1; ++k) {
/*<       if(m(3,k) .le. 5) go to 17                                         >*/
    if (m[(k << 2) - 2] <= 5) {
        goto L17;
    }
/*    write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') =
'' 4506*/
/*   1)')  k,m(1,k),m(2,k)                                            
   4507*/
/*    write(it,'(100('' '',10f7.1))') (vm(l),l=m(4,k)+1,m(4,k)+m(3,k))
   4508*/
/*<       go to 18                                                           >*/
    goto L18;
/*<    17 continue >*/
L17:
/*    write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') =
'' 4510*/
/*   1,5f7.1)')  k,m(1,k),m(2,k),(vm(l),l=m(4,k)+1,m(4,k)+m(3,k))     
   4511*/
/*<    18 continue                                                           >*/
L18:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       entry oknest(it,p,lx,cm)                                           >*/

L_oknest:
/*<       if(it.le.0) return                                                 >*/
    if (*it <= 0) {
    return 0;
    }
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 20                                                           >*/
    goto L20;
/*<    19 l=l+1                                                              >*/
L19:
    ++l;
/*<    20 if((l).gt.(il)) go to 24                                           >*/
L20:
    if (l > il) {
    goto L24;
    }
/*<       j1=m(1,l)                                                          >*/
    j1 = m[(l << 2) - 4];
/*<       jn=m(2,l)                                                          >*/
    jn = m[(l << 2) - 3];
/*<       jv=m(3,l)                                                          >*/
    jv = m[(l << 2) - 2];
/*<       jp=m(4,l)                                                          >*/
    jp = m[(l << 2) - 1];
/*    if(j1.lt.1.or.j1.gt.p) write(it,25) l,j1                           4
524*/
/*    if(jn.lt.1.or.jn.gt.p) write(it,26) l,jn                           4
525*/
/*    if(lx(jn).ge.0) write(it,27) l,jn,lx(jn)                           4
526*/
/*<       k1=cm(2*jn)+.1                                                     >*/
    k1 = cm[jn * 2] + (float).1;
/*<       k2=cm(2*jn+1)+.1                                                   >*/
    k2 = cm[(jn << 1) + 1] + (float).1;
/*<       do 23 k=jp+1,jp+jv                                                 >*/
    i__1 = jp + jv;
    for (k = jp + 1; k <= i__1; ++k) {
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 21 kk=k1,k2                                                     >*/
    i__2 = k2;
    for (kk = k1; kk <= i__2; ++kk) {
/*<       if(vm(k) .ne. cm(kk)) go to 21                                     >*/
        if (vm[k - 1] != cm[kk]) {
        goto L21;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 22                                                           >*/
        goto L22;
/*<    21 continue                                                           >*/
L21:
        ;
    }
/*<    22 continue >*/
L22:
/*    if(ig.eq.0) write(it,28) l,vm(k),jn                             
   4536*/
/*<    23 continue                                                           >*/
/* L23: */
    ;
    }
/*<       go to 19                                                           >*/
    goto L19;
/*<    24 return                                                             >*/
L24:
    return 0;
/*<    >*/
/* L25: */
/*<    26 format(' nesting entry',i3,', invalid nesting variable',i3,'.')    >*/
/* L26: */
/*<    27 format(' nesting entry',i3,', lx(',i3,') =',i2,'. must be < 0.')   >*/
/* L27: */
/*<    >*/
/* L28: */
/*<       entry isnstr(j,jb)                                                 >*/

L_isnstr:
/*<       jb=0                                                               >*/
    *jb = 0;
/*<       k=1                                                                >*/
    k = 1;
/*<       go to 30                                                           >*/
    goto L30;
/*<    29 k=k+1                                                              >*/
L29:
    ++k;
/*<    30 if((k).gt.(il)) go to 32                                           >*/
L30:
    if (k > il) {
    goto L32;
    }
/*<       if(m(2,k) .ne. j) go to 29                                         >*/
    if (m[(k << 2) - 3] != *j) {
    goto L29;
    }
/*<       jb=m(1,k)                                                          >*/
    *jb = m[(k << 2) - 4];
/*<    32 return                                                             >*/
L32:
    return 0;
/*<       entry isfac (lm,j,mk,tb,cm,ja)                                     >*/

L_isfac:
/*<       ja=0                                                               >*/
    *ja = 0;
/*<       ig=ja                                                              >*/
    ig = *ja;
/*<       l=1                                                                >*/
    l = 1;
/*<       go to 34                                                           >*/
    goto L34;
/*<    33 l=l+1                                                              >*/
L33:
    ++l;
/*<    34 if((l).gt.(il)) go to 36                                           >*/
L34:
    if (l > il) {
    goto L36;
    }
/*<       if(j .ne. m(1,l)) go to 33                                         >*/
    if (*j != m[(l << 2) - 4]) {
    goto L33;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<    36 if(ig.eq.0) return                                                 >*/
L36:
    if (ig == 0) {
    return 0;
    }
/*<       jn=m(2,l)                                                          >*/
    jn = m[(l << 2) - 3];
/*<       if(cm(2*jn).eq.0.0) return                                         >*/
    if (cm[jn * 2] == (float)0.) {
    return 0;
    }
/*<       jv=m(3,l)                                                          >*/
    jv = m[(l << 2) - 2];
/*<       jp=m(4,l)                                                          >*/
    jp = m[(l << 2) - 1];
/*<       ig=0                                                               >*/
    ig = 0;
/*<       ip=lm                                                              >*/
    ip = *lm;
/*<    37 if(ip.le.0) go to 39                                               >*/
L37:
    if (ip <= 0) {
    goto L39;
    }
/*<       j1=abs(tb(2,ip))+.1                                                >*/
    j1 = (r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1;
/*<       if(j1 .ne. jn) go to 38                                            >*/
    if (j1 != jn) {
    goto L38;
    }
/*<       ig=1                                                               >*/
    ig = 1;
/*<       go to 39                                                           >*/
    goto L39;
/*<    38 ip=tb(4,ip)+.1                                                     >*/
L38:
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 37                                                           >*/
    goto L37;
/*<    39 if(ig .eq. 0) go to 45                                             >*/
L39:
    if (ig == 0) {
    goto L45;
    }
/*<       nc=cm(2*jn+1)-cm(2*jn)+1.1                                         >*/
    nc = cm[(jn << 1) + 1] - cm[jn * 2] + (float)1.1;
/*<       t=tb(2,ip)                                                         >*/
    t = tb[ip * 5 + 2];
/*<       kp=tb(3,ip)+.1                                                     >*/
    kp = tb[ip * 5 + 3] + (float).1;
/*<       do 44 l=1,nc                                                       >*/
    i__1 = nc;
    for (l = 1; l <= i__1; ++l) {
/*<       lp=l+kp                                                            >*/
    lp = l + kp;
/*<       if(t .le. 0) go to 40                                              >*/
    if (t <= (float)0.) {
        goto L40;
    }
/*<       if(cm(lp).eq.0.0) go to 44                                         >*/
    if (cm[lp] == (float)0.) {
        goto L44;
    }
/*<       go to 41                                                           >*/
    goto L41;
/*<    40 if(cm(lp).ne.0.0) go to 44                                         >*/
L40:
    if (cm[lp] != (float)0.) {
        goto L44;
    }
/*<    41 ex=cm(int(cm(2*jn)+.1)+l-1)                                        >*/
L41:
    ex = cm[(integer) (cm[jn * 2] + (float).1) + l - 1];
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 42 k=jp+1,jp+jv                                                 >*/
    i__2 = jp + jv;
    for (k = jp + 1; k <= i__2; ++k) {
/*<       if(ex .ne. vm(k)) go to 42                                         >*/
        if (ex != vm[k - 1]) {
        goto L42;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 43                                                           >*/
        goto L43;
/*<    42 continue                                                           >*/
L42:
        ;
    }
/*<    43 if(ig .ne. 0) go to 44                                             >*/
L43:
    if (ig != 0) {
        goto L44;
    }
/*<       ja=-1                                                              >*/
    *ja = -1;
/*<       return                                                             >*/
    return 0;
/*<    44 continue                                                           >*/
L44:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<    45 ja=l                                                               >*/
L45:
    *ja = l;
/*<       norm=nord(lm,tb)+1                                                 >*/
    norm = nord_(lm, &tb[6]) + 1;
/*<       nc=cm(2*jn+1)-cm(2*jn)+1.1                                         >*/
    nc = cm[(jn << 1) + 1] - cm[jn * 2] + (float)1.1;
/*<       do 56 lk=1,mk                                                      >*/
    i__1 = *mk;
    for (lk = 1; lk <= i__1; ++lk) {
/*<       if(nord(lk,tb).ne.norm) go to 56                                   >*/
    if (nord_(&lk, &tb[6]) != norm) {
        goto L56;
    }
/*<       jg=0                                                               >*/
    jg = 0;
/*<       ip=lk                                                              >*/
    ip = lk;
/*<    46 if(ip.le.0) go to 55                                               >*/
L46:
    if (ip <= 0) {
        goto L55;
    }
/*<       t1=tb(2,ip)                                                        >*/
    t1 = tb[ip * 5 + 2];
/*<       j1=abs(t1)+.1                                                      >*/
    j1 = dabs(t1) + (float).1;
/*<       if(j1 .ne. jn) go to 54                                            >*/
    if (j1 != jn) {
        goto L54;
    }
/*<       kp=tb(3,ip)+.1                                                     >*/
    kp = tb[ip * 5 + 3] + (float).1;
/*<       kg=0                                                               >*/
    kg = 0;
/*<       do 52 l=1,nc                                                       >*/
    i__2 = nc;
    for (l = 1; l <= i__2; ++l) {
/*<       lp=l+kp                                                            >*/
        lp = l + kp;
/*<       lon=cm(lp)+.1                                                      >*/
        lon = cm[lp] + (float).1;
/*<       if(t1 .ge. 0.0) go to 48                                           >*/
        if (t1 >= (float)0.) {
        goto L48;
        }
/*<       if(lon .ne. 0) go to 47                                            >*/
        if (lon != 0) {
        goto L47;
        }
/*<       lon=1                                                              >*/
        lon = 1;
/*<       go to 48                                                           >*/
        goto L48;
/*<    47 lon=0                                                              >*/
L47:
        lon = 0;
/*<    48 ex=cm(int(cm(2*jn)+.1)+l-1)                                        >*/
L48:
        ex = cm[(integer) (cm[jn * 2] + (float).1) + l - 1];
/*<       ig=0                                                               >*/
        ig = 0;
/*<       do 49 k=jp+1,jp+jv                                                 >*/
        i__3 = jp + jv;
        for (k = jp + 1; k <= i__3; ++k) {
/*<       if(ex .ne. vm(k)) go to 49                                         >*/
        if (ex != vm[k - 1]) {
            goto L49;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 50                                                           >*/
        goto L50;
/*<    49 continue                                                           >*/
L49:
        ;
        }
/*<    50 if(lon .ne. 1 .or. ig .ne. 0) go to 51                             >*/
L50:
        if (lon != 1 || ig != 0) {
        goto L51;
        }
/*<       kg=1                                                               >*/
        kg = 1;
/*<       go to 53                                                           >*/
        goto L53;
/*<    51 if(lon .ne. 0 .or. ig .ne. 1) go to 52                             >*/
L51:
        if (lon != 0 || ig != 1) {
        goto L52;
        }
/*<       kg=1                                                               >*/
        kg = 1;
/*<       go to 53                                                           >*/
        goto L53;
/*<    52 continue                                                           >*/
L52:
        ;
    }
/*<    53 if(kg .ne. 0) go to 54                                             >*/
L53:
    if (kg != 0) {
        goto L54;
    }
/*<       jg=1                                                               >*/
    jg = 1;
/*<       go to 55                                                           >*/
    goto L55;
/*<    54 ip=tb(4,ip)+.1                                                     >*/
L54:
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 46                                                           >*/
    goto L46;
/*<    55 if(jg.eq.0) go to 56                                               >*/
L55:
    if (jg == 0) {
        goto L56;
    }
/*<       if(ieqbf(lk,lm,tb,cm) .ne. 1) go to 56                             >*/
    if (ieqbf_(&lk, lm, &tb[6], &cm[1]) != 1) {
        goto L56;
    }
/*<       ja=-1                                                              >*/
    *ja = -1;
/*<       return                                                             >*/
    return 0;
/*<    56 continue                                                           >*/
L56:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       entry cmpnst(ja,n,x,cm,bl)                                         >*/

L_cmpnst:
/*<       jn=m(2,ja)                                                         >*/
    jn = m[(*ja << 2) - 3];
/*<       jv=m(3,ja)                                                         >*/
    jv = m[(*ja << 2) - 2];
/*<       jp=m(4,ja)                                                         >*/
    jp = m[(*ja << 2) - 1];
/*<       do 59 l=1,n                                                        >*/
    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
/*<       kx=x(l,jn)+.1                                                      >*/
    kx = x[l + jn * x_dim1] + (float).1;
/*<       ex=cm(int(cm(2*jn)+.1)+kx-1)                                       >*/
    ex = cm[(integer) (cm[jn * 2] + (float).1) + kx - 1];
/*<       ig=0                                                               >*/
    ig = 0;
/*<       do 57 k=jp+1,jp+jv                                                 >*/
    i__2 = jp + jv;
    for (k = jp + 1; k <= i__2; ++k) {
/*<       if(ex .ne. vm(k)) go to 57                                         >*/
        if (ex != vm[k - 1]) {
        goto L57;
        }
/*<       ig=1                                                               >*/
        ig = 1;
/*<       go to 58                                                           >*/
        goto L58;
/*<    57 continue                                                           >*/
L57:
        ;
    }
/*<    58 if(ig.eq.1) go to 59                                               >*/
L58:
    if (ig == 1) {
        goto L59;
    }
/*<       bl(l)=0.0                                                          >*/
    bl[l] = (float)0.;
/*<    59 continue                                                           >*/
L59:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       entry getnst(ja,cm,j,nv,vals)                                      >*/

L_getnst:
/*<       j=m(2,ja)                                                          >*/
    *j = m[(*ja << 2) - 3];
/*<       jv=m(3,ja)                                                         >*/
    jv = m[(*ja << 2) - 2];
/*<       jp=m(4,ja)                                                         >*/
    jp = m[(*ja << 2) - 1];
/*<       nv=cm(2*j+1)-cm(2*j)+1.1                                           >*/
    *nv = cm[(*j << 1) + 1] - cm[*j * 2] + (float)1.1;
/*<       do 60 k=1,nv                                                       >*/
    i__1 = *nv;
    for (k = 1; k <= i__1; ++k) {
/*<       vals(k)=0.0                                                        >*/
    vals[k] = (float)0.;
/*<    60 continue                                                           >*/
/* L60: */
    }
/*<       do 61 l=jp+1,jp+jv                                                 >*/
    i__1 = jp + jv;
    for (l = jp + 1; l <= i__1; ++l) {
/*<       k=icat(vm(l),j,cm)                                                 >*/
    k = icat_(&vm[l - 1], j, &cm[1]);
/*<       if(k.gt.0) vals(k)=1.0                                             >*/
    if (k > 0) {
        vals[k] = (float)1.;
    }
/*<    61 continue                                                           >*/
/* L61: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* nest_ */

/* Subroutine */ int nest_(n, i__, j, nv, vals)
integer *n, *i__, *j, *nv;
real *vals;
{
    return nest_0_(0, n, i__, j, nv, vals, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, (integer *)0, (integer *)0, (integer *)0, 
        (real *)0, (integer *)0, (real *)0, (real *)0);
    }

/* Subroutine */ int nstlst_(it)
integer *it;
{
    return nest_0_(1, (integer *)0, (integer *)0, (integer *)0, (integer *)0, 
        (real *)0, it, (integer *)0, (integer *)0, (real *)0, (integer *)
        0, (integer *)0, (integer *)0, (real *)0, (integer *)0, (real *)0,
         (real *)0);
    }

/* Subroutine */ int oknest_(it, p, lx, cm)
integer *it, *p, *lx;
real *cm;
{
    return nest_0_(2, (integer *)0, (integer *)0, (integer *)0, (integer *)0, 
        (real *)0, it, p, lx, cm, (integer *)0, (integer *)0, (integer *)
        0, (real *)0, (integer *)0, (real *)0, (real *)0);
    }

/* Subroutine */ int isnstr_(j, jb)
integer *j, *jb;
{
    return nest_0_(3, (integer *)0, (integer *)0, j, (integer *)0, (real *)0, 
        (integer *)0, (integer *)0, (integer *)0, (real *)0, jb, (integer 
        *)0, (integer *)0, (real *)0, (integer *)0, (real *)0, (real *)0);
    }

/* Subroutine */ int isfac_(lm, j, mk, tb, cm, ja)
integer *lm, *j, *mk;
real *tb, *cm;
integer *ja;
{
    return nest_0_(4, (integer *)0, (integer *)0, j, (integer *)0, (real *)0, 
        (integer *)0, (integer *)0, (integer *)0, cm, (integer *)0, lm, 
        mk, tb, ja, (real *)0, (real *)0);
    }

/* Subroutine */ int cmpnst_(ja, n, x, cm, bl)
integer *ja, *n;
real *x, *cm, *bl;
{
    return nest_0_(5, n, (integer *)0, (integer *)0, (integer *)0, (real *)0, 
        (integer *)0, (integer *)0, (integer *)0, cm, (integer *)0, (
        integer *)0, (integer *)0, (real *)0, ja, x, bl);
    }

/* Subroutine */ int getnst_(ja, cm, j, nv, vals)
integer *ja;
real *cm;
integer *j, *nv;
real *vals;
{
    return nest_0_(6, (integer *)0, (integer *)0, j, nv, vals, (integer *)0, (
        integer *)0, (integer *)0, cm, (integer *)0, (integer *)0, (
        integer *)0, (real *)0, ja, (real *)0, (real *)0);
    }

/*<       subroutine blf0(l,ja,n,x,w,cm,sc,nnt,bl)                           >*/
/* Subroutine */ int blf0_(l, ja, n, x, w, cm, sc, nnt, bl)
integer *l, *ja, *n;
real *x, *w, *cm, *sc;
integer *nnt;
real *bl;
{
    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int cmpnst_(), blf_();

/*<       real x(n,*),w(n),cm(*),sc(n,*),bl(n)                               >*/
/*<       nnt=0                                                              >*/
    /* Parameter adjustments */
    --bl;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --w;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --cm;

    /* Function Body */
    *nnt = 0;
/*<       call blf(l,n,sc,bl)                                                >*/
    blf_(l, n, &sc[sc_offset], &bl[1]);
/*<       if(ja.gt.0) call cmpnst(ja,n,x,cm,bl)                              >*/
    if (*ja > 0) {
    cmpnst_(ja, n, &x[x_offset], &cm[1], &bl[1]);
    }
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       if(bl(i).gt.0.0.and.w(i).gt.0.0) nnt=nnt+1                         >*/
    if (bl[i__] > (float)0. && w[i__] > (float)0.) {
        ++(*nnt);
    }
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* blf0_ */

/*<       subroutine blf(l,n,sc,bl)                                          >*/
/* Subroutine */ int blf_(l, n, sc, bl)
integer *l, *n;
real *sc, *bl;
{
    /* System generated locals */
    integer sc_dim1, sc_offset, i__1;

    /* Local variables */
    static integer i__;

/*<       real sc(n,*),bl(n)                                                 >*/
/*<       if(l .gt. 0) go to 2                                               >*/
    /* Parameter adjustments */
    --bl;
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;

    /* Function Body */
    if (*l > 0) {
    goto L2;
    }
/*<       do 1 i=1,n                                                         >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       bl(i)=1.0                                                          >*/
    bl[i__] = (float)1.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<     2 do 3 i=1,n                                                         >*/
L2:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       bl(i)=sc(i,l)                                                      >*/
    bl[i__] = sc[i__ + *l * sc_dim1];
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<     4 return                                                             >*/
L4:
    return 0;
/*<       end                                                                >*/
} /* blf_ */

/*<       subroutine mnspan(ms,alf,nep,nnt,mn,me,mel)                        >*/
/* Subroutine */ int mnspan_(ms, alf, nep, nnt, mn, me, mel)
integer *ms;
real *alf;
integer *nep, *nnt, *mn, *me, *mel;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log();

    /* Local variables */
    static real allf, fme, fmn;
    static integer nnl, nnr, nst;

/*<       parameter(al2=0.693147,al25=1.732868)                              >*/
/*<       allf=-alog(1.0-alf)                                                >*/
    allf = -log((float)1. - *alf);
/*<       fmn=-alog(allf/(nep*nnt))/al25                                     >*/
    fmn = -log(allf / (*nep * *nnt)) / (float)1.732868;
/*<       fme=-alog(alf*0.125/nep)/al2                                       >*/
    fme = -log(*alf * (float).125 / *nep) / (float).693147;
/*<       if(ms .le. 0) go to 1                                              >*/
    if (*ms <= 0) {
    goto L1;
    }
/*<       me=ms*fme/fmn+0.5                                                  >*/
    *me = *ms * fme / fmn + (float).5;
/*<       mn=ms                                                              >*/
    *mn = *ms;
/*<       go to 2                                                            >*/
    goto L2;
/*<     1 me=fme+0.5                                                         >*/
L1:
    *me = fme + (float).5;
/*<       mn=fmn+0.5                                                         >*/
    *mn = fmn + (float).5;
/*<     2 me=max0(me,mn,2)                                                   >*/
L2:
/* Computing MAX */
    i__1 = max(*me,*mn);
    *me = max(i__1,2);
/*<       nst=nnt-2*me-1                                                     >*/
    nst = *nnt - (*me << 1) - 1;
/*<       nnr=nst/mn                                                         >*/
    nnr = nst / *mn;
/*<       nnl=nst-nnr*mn                                                     >*/
    nnl = nst - nnr * *mn;
/*<       nnr=(nnr+1)*mn-nst                                                 >*/
    nnr = (nnr + 1) * *mn - nst;
/*<       nst=min0(nnl,nnr)                                                  >*/
    nst = min(nnl,nnr);
/*<       if(nnl .gt. nnr) go to 3                                           >*/
    if (nnl > nnr) {
    goto L3;
    }
/*<       nnl=1                                                              >*/
    nnl = 1;
/*<       go to 4                                                            >*/
    goto L4;
/*<     3 nnl=-1                                                             >*/
L3:
    nnl = -1;
/*<     4 nnr=nst/2                                                          >*/
L4:
    nnr = nst / 2;
/*<       mel=me                                                             >*/
    *mel = *me;
/*<       me=me+nnl*nnr                                                      >*/
    *me += nnl * nnr;
/*<       mel=mel+nnl*nnr                                                    >*/
    *mel += nnl * nnr;
/*<       if(mod(nst,2).ne.0) mel=mel+nnl                                    >*/
    if (nst % 2 != 0) {
    *mel += nnl;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* mnspan_ */

/*<       function ieqbf(lk,lm,tb,cm)                                        >*/
integer ieqbf_(lk, lm, tb, cm)
integer *lk, *lm;
real *tb, *cm;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer l;
    static real t;
    static integer j1;
    static real t1;
    static integer ic, jg, nc, lg, kg, jo, ko, ip, kp, lo, lp;
    static real to;
    extern integer ieq_();
    static integer ipo, lon, lop;

/*<       real tb(5,*),cm(*)                                                 >*/
/*<       ipo=lm                                                             >*/
    /* Parameter adjustments */
    --cm;
    tb -= 6;

    /* Function Body */
    ipo = *lm;
/*<       lg=0                                                               >*/
    lg = 0;
/*<     1 if(ipo.le.0) go to 16                                              >*/
L1:
    if (ipo <= 0) {
    goto L16;
    }
/*<       to=tb(2,ipo)                                                       >*/
    to = tb[ipo * 5 + 2];
/*<       jo=abs(to)+.1                                                      >*/
    jo = dabs(to) + (float).1;
/*<       jg=0                                                               >*/
    jg = 0;
/*<       if(cm(2*jo) .ne. 0.0) go to 2                                      >*/
    if (cm[jo * 2] != (float)0.) {
    goto L2;
    }
/*<       t=tb(3,ipo)                                                        >*/
    t = tb[ipo * 5 + 3];
/*<       ic=0                                                               >*/
    ic = 0;
/*<       go to 3                                                            >*/
    goto L3;
/*<     2 ko=tb(3,ipo)+.1                                                    >*/
L2:
    ko = tb[ipo * 5 + 3] + (float).1;
/*<       nc=cm(2*jo+1)-cm(2*jo)+1.1                                         >*/
    nc = cm[(jo << 1) + 1] - cm[jo * 2] + (float)1.1;
/*<       ic=1                                                               >*/
    ic = 1;
/*<     3 ip=lk                                                              >*/
L3:
    ip = *lk;
/*<     4 if(ip.le.0) go to 14                                               >*/
L4:
    if (ip <= 0) {
    goto L14;
    }
/*<       t1=tb(2,ip)                                                        >*/
    t1 = tb[ip * 5 + 2];
/*<       j1=abs(t1)+.1                                                      >*/
    j1 = dabs(t1) + (float).1;
/*<       if(j1 .ne. jo) go to 13                                            >*/
    if (j1 != jo) {
    goto L13;
    }
/*<       if(ic .ne. 0) go to 6                                              >*/
    if (ic != 0) {
    goto L6;
    }
/*<       if(to*t1 .le. 0.0) go to 13                                        >*/
    if (to * t1 <= (float)0.) {
    goto L13;
    }
/*<       if(ieq(t,tb(3,ip),1.0) .ne. 1) go to 13                            >*/
    if (ieq_(&t, &tb[ip * 5 + 3], &c_b196) != 1) {
    goto L13;
    }
/*<       jg=1                                                               >*/
    jg = 1;
/*<       go to 14                                                           >*/
    goto L14;
/*<     6 kp=tb(3,ip)+.1                                                     >*/
L6:
    kp = tb[ip * 5 + 3] + (float).1;
/*<       kg=0                                                               >*/
    kg = 0;
/*<       do 11 l=1,nc                                                       >*/
    i__1 = nc;
    for (l = 1; l <= i__1; ++l) {
/*<       lo=l+ko                                                            >*/
    lo = l + ko;
/*<       lp=l+kp                                                            >*/
    lp = l + kp;
/*<       lon=cm(lo)+.1                                                      >*/
    lon = cm[lo] + (float).1;
/*<       lop=cm(lp)+.1                                                      >*/
    lop = cm[lp] + (float).1;
/*<       if(to .ge. 0.0) go to 8                                            >*/
    if (to >= (float)0.) {
        goto L8;
    }
/*<       if(lon .ne. 0) go to 7                                             >*/
    if (lon != 0) {
        goto L7;
    }
/*<       lon=1                                                              >*/
    lon = 1;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 lon=0                                                              >*/
L7:
    lon = 0;
/*<     8 if(t1 .ge. 0.0) go to 10                                           >*/
L8:
    if (t1 >= (float)0.) {
        goto L10;
    }
/*<       if(lop .ne. 0) go to 9                                             >*/
    if (lop != 0) {
        goto L9;
    }
/*<       lop=1                                                              >*/
    lop = 1;
/*<       go to 10                                                           >*/
    goto L10;
/*<     9 lop=0                                                              >*/
L9:
    lop = 0;
/*<    10 if(lon .eq. lop) go to 11                                          >*/
L10:
    if (lon == lop) {
        goto L11;
    }
/*<       kg=1                                                               >*/
    kg = 1;
/*<       go to 12                                                           >*/
    goto L12;
/*<    11 continue                                                           >*/
L11:
    ;
    }
/*<    12 if(kg .ne. 0) go to 13                                             >*/
L12:
    if (kg != 0) {
    goto L13;
    }
/*<       jg=1                                                               >*/
    jg = 1;
/*<       go to 14                                                           >*/
    goto L14;
/*<    13 ip=tb(4,ip)+.1                                                     >*/
L13:
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 4                                                            >*/
    goto L4;
/*<    14 if(jg .ne. 0) go to 15                                             >*/
L14:
    if (jg != 0) {
    goto L15;
    }
/*<       lg=1                                                               >*/
    lg = 1;
/*<       go to 16                                                           >*/
    goto L16;
/*<    15 ipo=tb(4,ipo)+.1                                                   >*/
L15:
    ipo = tb[ipo * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<    16 if(lg .ne. 0) go to 17                                             >*/
L16:
    if (lg != 0) {
    goto L17;
    }
/*<       ieqbf=1                                                            >*/
    ret_val = 1;
/*<       go to 18                                                           >*/
    goto L18;
/*<    17 ieqbf=0                                                            >*/
L17:
    ret_val = 0;
/*<    18 return                                                             >*/
L18:
    return ret_val;
/*<       end                                                                >*/
} /* ieqbf_ */

/*<       function ibfext(m,tb,cm)                                           >*/
integer ibfext_(m, tb, cm)
integer *m;
real *tb, *cm;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    extern integer nord_();
    static integer norm, l;
    extern integer ieqbf_();
    static integer mm1;

/*<       real tb(5,*),cm(*)                                                 >*/
/*<       mm1=m-1                                                            >*/
    /* Parameter adjustments */
    --cm;
    tb -= 6;

    /* Function Body */
    mm1 = *m - 1;
/*<       ibfext=0                                                           >*/
    ret_val = 0;
/*<       norm=nord(m,tb)                                                    >*/
    norm = nord_(m, &tb[6]);
/*<       do 1 l=1,mm1                                                       >*/
    i__1 = mm1;
    for (l = 1; l <= i__1; ++l) {
/*<       if(nord(l,tb).ne.norm) go to 1                                     >*/
    if (nord_(&l, &tb[6]) != norm) {
        goto L1;
    }
/*<       if(ieqbf(l,m,tb,cm) .eq. 0) go to 1                                >*/
    if (ieqbf_(&l, m, &tb[6], &cm[1]) == 0) {
        goto L1;
    }
/*<       ibfext=1                                                           >*/
    ret_val = 1;
/*<       return                                                             >*/
    return ret_val;
/*<     1 continue                                                           >*/
L1:
    ;
    }
/*<       return                                                             >*/
    return ret_val;
/*<       end                                                                >*/
} /* ibfext_ */

/*<       subroutine miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp)                  >*/
/* Subroutine */ int miss_(n, p, x, lx, xm, flg, pn, xn, lxn, xs, xp)
integer *n, *p;
real *x;
integer *lx;
real *xm, *flg;
integer *pn;
real *xn;
integer *lxn;
real *xs, *xp;
{
    /* System generated locals */
    integer x_dim1, x_offset, xn_dim1, xn_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int nest_();
    static integer i__, j;
    static doublereal s;
    static integer mf;
    static real ss;

/*<       integer p,pn,lx(*),lxn(*)                                          >*/
/*<       real x(n,*),xm(*),xn(n,*),xs(*),xp(*)                              >*/
/*<       double precision s                                                 >*/
/*<       pn=p                                                               >*/
    /* Parameter adjustments */
    xn_dim1 = *n;
    xn_offset = xn_dim1 + 1;
    xn -= xn_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --lx;
    --xm;
    --lxn;
    --xs;
    --xp;

    /* Function Body */
    *pn = *p;
/*<       xp(1)=p                                                            >*/
    xp[1] = (real) (*p);
/*<       do 1 j=2,2*p+1                                                     >*/
    i__1 = (*p << 1) + 1;
    for (j = 2; j <= i__1; ++j) {
/*<       xp(j)=0.0                                                          >*/
    xp[j] = (float)0.;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       ss=0.0                                                             >*/
    ss = (float)0.;
/*<       do 7 j=1,p                                                         >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       lxn(j)=lx(j)                                                       >*/
    lxn[j] = lx[j];
/*<       xs(j)=flg                                                          >*/
    xs[j] = *flg;
/*<       if(lx(j).eq.0) go to 7                                             >*/
    if (lx[j] == 0) {
        goto L7;
    }
/*<       s=0.d0                                                             >*/
    s = 0.;
/*<       mf=0                                                               >*/
    mf = 0;
/*<       do 4 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       xn(i,j)=x(i,j)                                                     >*/
        xn[i__ + j * xn_dim1] = x[i__ + j * x_dim1];
/*<       if(x(i,j) .ne. xm(j)) go to 2                                      >*/
        if (x[i__ + j * x_dim1] != xm[j]) {
        goto L2;
        }
/*<       mf=mf+1                                                            >*/
        ++mf;
/*<       go to 4                                                            >*/
        goto L4;
/*<     2 if(lx(j) .ge. 0) go to 3                                           >*/
L2:
        if (lx[j] >= 0) {
        goto L3;
        }
/*<       ss=x(i,j)                                                          >*/
        ss = x[i__ + j * x_dim1];
/*<       go to 4                                                            >*/
        goto L4;
/*<     3 s=s+x(i,j)                                                         >*/
L3:
        s += x[i__ + j * x_dim1];
/*<     4 continue                                                           >*/
L4:
        ;
    }
/*<       if(mf.eq.0) go to 7                                                >*/
    if (mf == 0) {
        goto L7;
    }
/*<       if(mf .ne. n) go to 5                                              >*/
    if (mf != *n) {
        goto L5;
    }
/*<       lxn(j)=0                                                           >*/
    lxn[j] = 0;
/*<       go to 7                                                            >*/
    goto L7;
/*<     5 s=s/(n-mf)                                                         >*/
L5:
    s /= *n - mf;
/*<       pn=pn+1                                                            >*/
    ++(*pn);
/*<       lxn(pn)=-1                                                         >*/
    lxn[*pn] = -1;
/*<       xs(pn)=1.0                                                         >*/
    xs[*pn] = (float)1.;
/*<       xp(j+1)=pn                                                         >*/
    xp[j + 1] = (real) (*pn);
/*<       call nest(n,j,pn,1,1.0)                                            >*/
    nest_(n, &j, pn, &c__1, &c_b196);
/*<       if(lx(j).gt.0) ss=s                                                >*/
    if (lx[j] > 0) {
        ss = s;
    }
/*<       xp(j+p+1)=ss                                                       >*/
    xp[j + *p + 1] = ss;
/*<       do 6 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       xn(i,pn)=1.0                                                       >*/
        xn[i__ + *pn * xn_dim1] = (float)1.;
/*<       if(x(i,j) .ne. xm(j)) go to 6                                      >*/
        if (x[i__ + j * x_dim1] != xm[j]) {
        goto L6;
        }
/*<       xn(i,j)=ss                                                         >*/
        xn[i__ + j * xn_dim1] = ss;
/*<       xn(i,pn)=0.0                                                       >*/
        xn[i__ + *pn * xn_dim1] = (float)0.;
/*<     6 continue                                                           >*/
L6:
        ;
    }
/*<     7 continue                                                           >*/
L7:
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* miss_ */

/*<       subroutine mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc)            >*/
/* Subroutine */ int mkmiss_0_(n__, n, p, x, y, w, xm, pm, nnx, nn, xnn, yn, 
    wn, sc, val)
int n__;
integer *n, *p;
real *x, *y, *w, *xm, *pm;
integer *nnx, *nn;
real *xnn, *yn, *wn, *sc, *val;
{
    /* Initialized data */

    static real tol = (float).001;

    /* System generated locals */
    integer x_dim1, x_offset, sc_dim1, sc_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern /* Subroutine */ int rnms_();
    static integer i__, j, k, m[500];
    static real r__;
    static integer in, km, jp;
    static real fin;
    static integer nnk;
    static real cvx;

/*<       parameter(nlist=500)                                               >*/
/*<       integer p,m(nlist)                                                 >*/
/*<       real pm(p),xm(p),x(n,p),y(n),w(n),xnn(*),yn(*),wn(*),sc(p,*)       >*/
/*<       data tol /0.001/                                                   >*/
    /* Parameter adjustments */
    if (y) {
    --y;
    }
    if (w) {
    --w;
    }
    if (x) {
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    }
    if (xm) {
    --xm;
    }
    if (pm) {
    --pm;
    }
    if (sc) {
    sc_dim1 = *p;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    }
    if (xnn) {
    --xnn;
    }
    if (yn) {
    --yn;
    }
    if (wn) {
    --wn;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_smktol;
    }

/*<       if(p .le. nlist) go to 1                                           >*/
    if (*p <= 500) {
    goto L1;
    }
/*    write(6,'('' increase parameter nlist in subroutine mkmiss to '',i 4
855*/
/*   15,                      '' and recompile.'')') p                   4
856*/
/*<       stop                                                               >*/
    s_stop("", 0L);
/*<     1 do 3 j=1,p                                                         >*/
L1:
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       m(j)=0                                                             >*/
    m[j - 1] = 0;
/*<       do 2 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(x(i,j).eq.xm(j)) m(j)=m(j)+1                                    >*/
        if (x[i__ + j * x_dim1] == xm[j]) {
        ++m[j - 1];
        }
/*<       sc(j,i)=x(i,j)                                                     >*/
        sc[j + i__ * sc_dim1] = x[i__ + j * x_dim1];
/*<       yn(i)=y(i)                                                         >*/
        yn[i__] = y[i__];
/*<       wn(i)=w(i)                                                         >*/
        wn[i__] = w[i__];
/*<     2 continue                                                           >*/
/* L2: */
    }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       nn=n                                                               >*/
    *nn = *n;
/*<       jp=0                                                               >*/
    jp = 0;
/*<       km=jp                                                              >*/
    km = jp;
/*<     4 if(nn.ge.nnx.or.km.gt.p) go to 13                                  >*/
L4:
    if (*nn >= *nnx || km > *p) {
    goto L13;
    }
/*<       jp=jp+1                                                            >*/
    ++jp;
/*<       if(jp.gt.p) jp=1                                                   >*/
    if (jp > *p) {
    jp = 1;
    }
/*<       fin=nn*pm(jp)-m(jp)                                                >*/
    fin = *nn * pm[jp] - m[jp - 1];
/*<       if(fin .le. 0.0) go to 5                                           >*/
    if (fin <= (float)0.) {
    goto L5;
    }
/*<       in=fin+0.5                                                         >*/
    in = fin + (float).5;
/*<       go to 6                                                            >*/
    goto L6;
/*<     5 in=0                                                               >*/
L5:
    in = 0;
/*<     6 in=min0(in,nnx-nn)                                                 >*/
L6:
/* Computing MIN */
    i__1 = in, i__2 = *nnx - *nn;
    in = min(i__1,i__2);
/*<       if(in .le. 0) go to 7                                              >*/
    if (in <= 0) {
    goto L7;
    }
/*<       km=0                                                               >*/
    km = 0;
/*<       go to 8                                                            >*/
    goto L8;
/*<     7 km=km+1                                                            >*/
L7:
    ++km;
/*<       go to 4                                                            >*/
    goto L4;
/*<     8 do 11 k=1,in                                                       >*/
L8:
    i__1 = in;
    for (k = 1; k <= i__1; ++k) {
/*<       call rnms(r,1)                                                     >*/
    rnms_(&r__, &c__1);
/*<       i=nn*r+1.0                                                         >*/
    i__ = *nn * r__ + (float)1.;
/*<       nnk=nn+k                                                           >*/
    nnk = *nn + k;
/*<       do 9 j=1,p                                                         >*/
    i__2 = *p;
    for (j = 1; j <= i__2; ++j) {
/*<       sc(j,nnk)=sc(j,i)                                                  >*/
        sc[j + nnk * sc_dim1] = sc[j + i__ * sc_dim1];
/*<     9 continue                                                           >*/
/* L9: */
    }
/*<       sc(jp,nnk)=xm(jp)                                                  >*/
    sc[jp + nnk * sc_dim1] = xm[jp];
/*<       yn(nnk)=yn(i)                                                      >*/
    yn[nnk] = yn[i__];
/*<       wn(nnk)=wn(i)                                                      >*/
    wn[nnk] = wn[i__];
/*<       do 10 j=1,p                                                        >*/
    i__2 = *p;
    for (j = 1; j <= i__2; ++j) {
/*<       if(sc(j,nnk).eq.xm(j)) m(j)=m(j)+1                                 >*/
        if (sc[j + nnk * sc_dim1] == xm[j]) {
        ++m[j - 1];
        }
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<    11 continue                                                           >*/
/* L11: */
    }
/*<       nn=nn+in                                                           >*/
    *nn += in;
/*<       cvx=-9.9e30                                                        >*/
    cvx = (float)-9.9e30;
/*<       do 12 j=1,p                                                        >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       cvx=amax1(cvx,(nn*pm(j)-m(j))/float(nn))                           >*/
/* Computing MAX */
    r__1 = cvx, r__2 = (*nn * pm[j] - m[j - 1]) / (real) (*nn);
    cvx = dmax(r__1,r__2);
/*<    12 continue                                                           >*/
/* L12: */
    }
/*<       if(cvx.lt.tol) go to 13                                            >*/
    if (cvx < tol) {
    goto L13;
    }
/*<       go to 4                                                            >*/
    goto L4;
/*<    13 k=0                                                                >*/
L13:
    k = 0;
/*<       do 15 j=1,p                                                        >*/
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
/*<       do 14 i=1,nn                                                       >*/
    i__2 = *nn;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       k=k+1                                                              >*/
        ++k;
/*<       xnn(k)=sc(j,i)                                                     >*/
        xnn[k] = sc[j + i__ * sc_dim1];
/*<    14 continue                                                           >*/
/* L14: */
    }
/*<    15 continue                                                           >*/
/* L15: */
    }
/*<       return                                                             >*/
    return 0;
/*<       entry smktol(val)                                                  >*/

L_smktol:
/*<       tol=val                                                            >*/
    tol = *val;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* mkmiss_ */

/* Subroutine */ int mkmiss_(n, p, x, y, w, xm, pm, nnx, nn, xnn, yn, wn, sc)
integer *n, *p;
real *x, *y, *w, *xm, *pm;
integer *nnx, *nn;
real *xnn, *yn, *wn, *sc;
{
    return mkmiss_0_(0, n, p, x, y, w, xm, pm, nnx, nn, xnn, yn, wn, sc, (
        real *)0);
    }

/* Subroutine */ int smktol_(val)
real *val;
{
    return mkmiss_0_(1, (integer *)0, (integer *)0, (real *)0, (real *)0, (
        real *)0, (real *)0, (real *)0, (integer *)0, (integer *)0, (real 
        *)0, (real *)0, (real *)0, (real *)0, val);
    }

/*<       subroutine xmiss (n,x,xm,xp,xn)                                    >*/
/* Subroutine */ int xmiss_(n, x, xm, xp, xn)
integer *n;
real *x, *xm, *xp, *xn;
{
    /* System generated locals */
    integer x_dim1, x_offset, xn_dim1, xn_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, p;

/*<       real x(n,*),xm(*),xp(*),xn(n,*)                                    >*/
/*<       integer p                                                          >*/
/*<       p=xp(1)+.1                                                         >*/
    /* Parameter adjustments */
    xn_dim1 = *n;
    xn_offset = xn_dim1 + 1;
    xn -= xn_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --xm;
    --xp;

    /* Function Body */
    p = xp[1] + (float).1;
/*<       do 3 j=1,p                                                         >*/
    i__1 = p;
    for (j = 1; j <= i__1; ++j) {
/*<       k=xp(j+1)+.1                                                       >*/
    k = xp[j + 1] + (float).1;
/*<       do 2 i=1,n                                                         >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<       if(x(i,j) .eq. xm(j)) go to 1                                      >*/
        if (x[i__ + j * x_dim1] == xm[j]) {
        goto L1;
        }
/*<       xn(i,j)=x(i,j)                                                     >*/
        xn[i__ + j * xn_dim1] = x[i__ + j * x_dim1];
/*<       if(k.gt.0) xn(i,k)=1.0                                             >*/
        if (k > 0) {
        xn[i__ + k * xn_dim1] = (float)1.;
        }
/*<       go to 2                                                            >*/
        goto L2;
/*<     1 xn(i,j)=xp(j+p+1)                                                  >*/
L1:
        xn[i__ + j * xn_dim1] = xp[j + p + 1];
/*<       if(k.gt.0) xn(i,k)=0.0                                             >*/
        if (k > 0) {
        xn[i__ + k * xn_dim1] = (float)0.;
        }
/*<     2 continue                                                           >*/
L2:
        ;
    }
/*<     3 continue                                                           >*/
/* L3: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* xmiss_ */

/*<       function nnord(m,tb)                                               >*/
integer nnord_(m, tb)
integer *m;
real *tb;
{
    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static integer jb, ip;
    extern /* Subroutine */ int isnstr_();

/*<       real tb(5,*)                                                       >*/
/*<       ip=m                                                               >*/
    /* Parameter adjustments */
    tb -= 6;

    /* Function Body */
    ip = *m;
/*<       nnord=0                                                            >*/
    ret_val = 0;
/*<     1 if(ip.le.0) go to 2                                                >*/
L1:
    if (ip <= 0) {
    goto L2;
    }
/*<       call isnstr(int(abs(tb(2,ip))+.1),jb)                              >*/
    i__1 = (integer) ((r__1 = tb[ip * 5 + 2], dabs(r__1)) + (float).1);
    isnstr_(&i__1, &jb);
/*<       if(jb.eq.0) nnord=nnord+1                                          >*/
    if (jb == 0) {
    ++ret_val;
    }
/*<       ip=tb(4,ip)+.1                                                     >*/
    ip = tb[ip * 5 + 4] + (float).1;
/*<       go to 1                                                            >*/
    goto L1;
/*<     2 return                                                             >*/
L2:
    return ret_val;
/*<       end                                                                >*/
} /* nnord_ */

/*<       subroutine spofa(a,m,n,info)                                       >*/
/* Subroutine */ int spofa_(a, m, n, info)
doublereal *a;
integer *m, *n, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t, u;
    static integer j1, jm1, km1;

/*<       double precision a(m,*),s,t,u                                      >*/
/*<          j1 = info                                                       >*/
    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    j1 = *info;
/*<          do 30 j = j1, n                                                 >*/
    i__1 = *n;
    for (j = j1; j <= i__1; ++j) {
/*<             info=j                                                       >*/
    *info = j;
/*<             s = 0.0d0                                                    >*/
    s = 0.;
/*<             jm1 = j - 1                                                  >*/
    jm1 = j - 1;
/*<             if (jm1 .lt. 1) go to 20                                     >*/
    if (jm1 < 1) {
        goto L20;
    }
/*<             do 10 k = 1, jm1                                             >*/
    i__2 = jm1;
    for (k = 1; k <= i__2; ++k) {
/*<                u=0.0                                                     >*/
        u = (float)0.;
/*<                km1=k-1                                                   >*/
        km1 = k - 1;
/*<                if(km1.le.0) go to 40                                     >*/
        if (km1 <= 0) {
        goto L40;
        }
/*<                do 50 i=1,km1                                             >*/
        i__3 = km1;
        for (i__ = 1; i__ <= i__3; ++i__) {
/*<                   u=u+a(i,k)*a(i,j)                                      >*/
        u += a[i__ + k * a_dim1] * a[i__ + j * a_dim1];
/*<    50          continue                                                  >*/
/* L50: */
        }
/*<    40          continue                                                  >*/
L40:
/*<                t = a(k,j) - u                                            >*/
        t = a[k + j * a_dim1] - u;
/*<                t = t/a(k,k)                                              >*/
        t /= a[k + k * a_dim1];
/*<                a(k,j) = t                                                >*/
        a[k + j * a_dim1] = t;
/*<                s = s + t*t                                               >*/
        s += t * t;
/*<    10       continue                                                     >*/
/* L10: */
    }
/*<    20       continue                                                     >*/
L20:
/*<             s = a(j,j) - s                                               >*/
    s = a[j + j * a_dim1] - s;
/*<             if (s .le. 0.0d0)  return                                    >*/
    if (s <= 0.) {
        return 0;
    }
/*<             a(j,j) = dsqrt(s)                                            >*/
    a[j + j * a_dim1] = sqrt(s);
/*<    30    continue                                                        >*/
/* L30: */
    }
/*<       info=0                                                             >*/
    *info = 0;
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* spofa_ */

/*<       subroutine sposl(a,m,n,b)                                          >*/
/* Subroutine */ int sposl_(a, m, n, b)
doublereal *a;
integer *m, *n;
doublereal *b;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal t;
    static integer kb, km1;

/*<       double precision a(m,*),b(*),t                                     >*/
/*<       do 10 k = 1, n                                                     >*/
    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*<          t = 0.0                                                         >*/
    t = (float)0.;
/*<          km1=k-1                                                         >*/
    km1 = k - 1;
/*<          if(km1.le.0) go to 30                                           >*/
    if (km1 <= 0) {
        goto L30;
    }
/*<          do 40 i=1,km1                                                   >*/
    i__2 = km1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             t=t+a(i,k)*b(i)                                              >*/
        t += a[i__ + k * a_dim1] * b[i__];
/*<    40    continue                                                        >*/
/* L40: */
    }
/*<    30    continue                                                        >*/
L30:
/*<          b(k) = (b(k) - t)/a(k,k)                                        >*/
    b[k] = (b[k] - t) / a[k + k * a_dim1];
/*<    10 continue                                                           >*/
/* L10: */
    }
/*<       do 20 kb=1,n                                                       >*/
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
/*<       k=n+1-kb                                                           >*/
    k = *n + 1 - kb;
/*<       b(k)=b(k)/a(k,k)                                                   >*/
    b[k] /= a[k + k * a_dim1];
/*<       t=-b(k)                                                            >*/
    t = -b[k];
/*<       km1=k-1                                                            >*/
    km1 = k - 1;
/*<       if(km1.le.0) go to 50                                              >*/
    if (km1 <= 0) {
        goto L50;
    }
/*<       if(t.eq.0.0) go to 50                                              >*/
    if (t == (float)0.) {
        goto L50;
    }
/*<       do 60 i=1,km1                                                      >*/
    i__2 = km1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<          b(i)=b(i)+t*a(i,k)                                              >*/
        b[i__] += t * a[i__ + k * a_dim1];
/*<    60 continue                                                           >*/
/* L60: */
    }
/*<    50 continue                                                           >*/
L50:
/*<    20 continue                                                           >*/
/* L20: */
    ;
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* sposl_ */

/*<       subroutine psort (v,a,ii,jj)                                       >*/
/* Subroutine */ int psort_(v, a, ii, jj)
real *v;
integer *a, *ii, *jj;
{
    static integer i__, j, k, l, m, t, ij, il[20], iu[20], tt;
    static real vt, vtt;


/*     puts into a the permutation vector which sorts v into */
/*     increasing order. the array v is not modified. */
/*     only elements from ii to jj are considered. */
/*     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements */

/*     this is a modification of cacm algorithm #347 by r. c. singleton, 
*/
/*     which is a modified hoare quicksort. */

/*<       dimension a(jj),v(jj),iu(20),il(20)                                >*/
/*<       integer t,tt                                                       >*/
/*<       integer a                                                          >*/
/*<       real v                                                             >*/
/*<       m=1                                                                >*/
    /* Parameter adjustments */
    --a;
    --v;

    /* Function Body */
    m = 1;
/*<       i=ii                                                               >*/
    i__ = *ii;
/*<       j=jj                                                               >*/
    j = *jj;
/*<  10   if (i.ge.j) go to 80                                               >*/
L10:
    if (i__ >= j) {
    goto L80;
    }
/*<  20   k=i                                                                >*/
L20:
    k = i__;
/*<       ij=(j+i)/2                                                         >*/
    ij = (j + i__) / 2;
/*<       t=a(ij)                                                            >*/
    t = a[ij];
/*<       vt=v(t)                                                            >*/
    vt = v[t];
/*<       if (v(a(i)).le.vt) go to 30                                        >*/
    if (v[a[i__]] <= vt) {
    goto L30;
    }
/*<       a(ij)=a(i)                                                         >*/
    a[ij] = a[i__];
/*<       a(i)=t                                                             >*/
    a[i__] = t;
/*<       t=a(ij)                                                            >*/
    t = a[ij];
/*<       vt=v(t)                                                            >*/
    vt = v[t];
/*<  30   l=j                                                                >*/
L30:
    l = j;
/*<       if (v(a(j)).ge.vt) go to 50                                        >*/
    if (v[a[j]] >= vt) {
    goto L50;
    }
/*<       a(ij)=a(j)                                                         >*/
    a[ij] = a[j];
/*<       a(j)=t                                                             >*/
    a[j] = t;
/*<       t=a(ij)                                                            >*/
    t = a[ij];
/*<       vt=v(t)                                                            >*/
    vt = v[t];
/*<       if (v(a(i)).le.vt) go to 50                                        >*/
    if (v[a[i__]] <= vt) {
    goto L50;
    }
/*<       a(ij)=a(i)                                                         >*/
    a[ij] = a[i__];
/*<       a(i)=t                                                             >*/
    a[i__] = t;
/*<       t=a(ij)                                                            >*/
    t = a[ij];
/*<       vt=v(t)                                                            >*/
    vt = v[t];
/*<       go to 50                                                           >*/
    goto L50;
/*<  40   a(l)=a(k)                                                          >*/
L40:
    a[l] = a[k];
/*<       a(k)=tt                                                            >*/
    a[k] = tt;
/*<  50   l=l-1                                                              >*/
L50:
    --l;
/*<       if (v(a(l)).gt.vt) go to 50                                        >*/
    if (v[a[l]] > vt) {
    goto L50;
    }
/*<       tt=a(l)                                                            >*/
    tt = a[l];
/*<       vtt=v(tt)                                                          >*/
    vtt = v[tt];
/*<  60   k=k+1                                                              >*/
L60:
    ++k;
/*<       if (v(a(k)).lt.vt) go to 60                                        >*/
    if (v[a[k]] < vt) {
    goto L60;
    }
/*<       if (k.le.l) go to 40                                               >*/
    if (k <= l) {
    goto L40;
    }
/*<       if (l-i.le.j-k) go to 70                                           >*/
    if (l - i__ <= j - k) {
    goto L70;
    }
/*<       il(m)=i                                                            >*/
    il[m - 1] = i__;
/*<       iu(m)=l                                                            >*/
    iu[m - 1] = l;
/*<       i=k                                                                >*/
    i__ = k;
/*<       m=m+1                                                              >*/
    ++m;
/*<       go to 90                                                           >*/
    goto L90;
/*<  70   il(m)=k                                                            >*/
L70:
    il[m - 1] = k;
/*<       iu(m)=j                                                            >*/
    iu[m - 1] = j;
/*<       j=l                                                                >*/
    j = l;
/*<       m=m+1                                                              >*/
    ++m;
/*<       go to 90                                                           >*/
    goto L90;
/*<  80   m=m-1                                                              >*/
L80:
    --m;
/*<       if (m.eq.0) return                                                 >*/
    if (m == 0) {
    return 0;
    }
/*<       i=il(m)                                                            >*/
    i__ = il[m - 1];
/*<       j=iu(m)                                                            >*/
    j = iu[m - 1];
/*<  90   if (j-i.gt.10) go to 20                                            >*/
L90:
    if (j - i__ > 10) {
    goto L20;
    }
/*<       if (i.eq.ii) go to 10                                              >*/
    if (i__ == *ii) {
    goto L10;
    }
/*<       i=i-1                                                              >*/
    --i__;
/*<  100  i=i+1                                                              >*/
L100:
    ++i__;
/*<       if (i.eq.j) go to 80                                               >*/
    if (i__ == j) {
    goto L80;
    }
/*<       t=a(i+1)                                                           >*/
    t = a[i__ + 1];
/*<       vt=v(t)                                                            >*/
    vt = v[t];
/*<       if (v(a(i)).le.vt) go to 100                                       >*/
    if (v[a[i__]] <= vt) {
    goto L100;
    }
/*<       k=i                                                                >*/
    k = i__;
/*<  110  a(k+1)=a(k)                                                        >*/
L110:
    a[k + 1] = a[k];
/*<       k=k-1                                                              >*/
    --k;
/*<       if (vt.lt.v(a(k))) go to 110                                       >*/
    if (vt < v[a[k]]) {
    goto L110;
    }
/*<       a(k+1)=t                                                           >*/
    a[k + 1] = t;
/*<       go to 100                                                          >*/
    goto L100;
/*<       end                                                                >*/
} /* psort_ */

/*<       subroutine stseed (iseed)                                          >*/
/* Subroutine */ int stseed_0_(n__, iseed, x, n)
int n__;
integer *iseed;
real *x;
integer *n;
{
    /* Initialized data */

    static integer i__ = 987654321;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_mod();

    /* Local variables */
    static integer j;
    static doublereal u;

/*<       real x(1)                                                          >*/
/*<       double precision u                                                 >*/
/*<       data i /987654321/                                                 >*/
    /* Parameter adjustments */
    if (x) {
    --x;
    }

    /* Function Body */
    switch(n__) {
    case 1: goto L_rnms;
    }

/*<       i=iseed                                                            >*/
    i__ = *iseed;
/*<       return                                                             >*/
    return 0;
/*<       entry rnms (x,n)                                                   >*/

L_rnms:
/*<       do 1 j=1,n                                                         >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<       i=dmod(i*16807.d0,2147483647.d0)                                   >*/
    d__1 = i__ * 16807.;
    i__ = (integer) d_mod(&d__1, &c_b1305);
/*<       u=i                                                                >*/
    u = (doublereal) i__;
/*<       u=u*.465661287d-9                                                  >*/
    u *= 4.65661287e-10;
/*<       x(j)=u                                                             >*/
    x[j] = u;
/*<     1 continue                                                           >*/
/* L1: */
    }
/*<       return                                                             >*/
    return 0;
/*<       end                                                                >*/
} /* stseed_ */

/* Subroutine */ int stseed_(iseed)
integer *iseed;
{
    return stseed_0_(0, iseed, (real *)0, (integer *)0);
    }

/* Subroutine */ int rnms_(x, n)
real *x;
integer *n;
{
    return stseed_0_(1, (integer *)0, x, n);
    }

#endif /* NO_F2C */


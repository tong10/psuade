c
c Sample user coded driver program for running MARS 3.6
c
c This is a simple program to run MARS on a small data set for testing
c purposes. The input data (provided below) is assumed to be on a file
c named 'mars.data' in the local directory. Printed output (also provided
c below) is sent to the standard output (Fortran file unit 6). Note that
c this output might be slightly different than that obtained by running
c this program on a computer with different floating point arithmetic.
c The resulting model ought to be (nearly) equivalent.
c
c This program provides only the minimal input necessary to run MARS. It
c thereby makes maximal use of internally supplied defaults, and therefore
c does not take advantage of the many user options available for guiding
c the analysis (see MARS documentation).
c
c set up problem dependent parameters for this (small) example:
c  50 observations (cases).
c   5 predictor variables (inputs).
c   2 variable interactions maximum.
c  15 basis functions maximum.
c
      parameter (n=500, np=8, mi=3, nk=15)
c
c set up (generic) graphics parameters:
c 100 raster points for curves.
c  40 raster points (on each axis) for surfaces.
c  plot piecewise-cubic (continuous derivative) model.
c  show surfaces inside respective bivariate convex hulls.
c 
      parameter (ngc=100, ngs=40, m=2, icx=1)
c
c  set up working storage:
c  constant dimensioned arrays may need to be set to larger values
c  for much larger problems (see MARS documentation).
c

      real x(n,np),y(n),w(n),fm(50000),sp(50000)
      real crv(ngc,2,nk),srf(ngs,ngs,nk)
      integer lx(np),im(50000),mm(50000),nsamples,ninputs
      double precision dp(500000)
c
c set predictor variable flags to indicate all are unrestricted
c ordinal variables, and set observation weights:
c
      data lx,w /np*1,n*1.0/
c
c write output header:
c
      write(6,'(/,''  simple driver to test MARS 3.6. '')')
c
c read in data:
c
      open(10,file='psData',status='unknown')
      read(10,*) nsamples, ninputs
      if (nsamples .gt. n) go to 90
      if (ninputs .gt. np) go to 90
      do 1 i=1,nsamples
      read(10,*) (x(i,j),j=1,ninputs),y(i)
    1 continue
c
c invoke MARS:
c
      print *, 'calling mars'
      call mars (nsamples,ninputs,x,y,w,nk,mi,lx,fm,im,sp,dp,mm)
c
c construct plots for interpreting resulting model:
c

c     call plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm)
c
c write plots to output files for plotting with local graphics package:
c
c     if(nc.le.0) go to 2
c     open(11,file='mars.curves',status='unknown')
c     write(11,*) (((crv(i,j,k),i=1,ngc),j=1,2),k=1,nc)
c   2 continue
c     if(ns.le.0) go to 3
c     open(12,file='mars.surfs',status='unknown')
c     write(12,*) (((srf(i,j,k),i=1,ngs),j=1,ngs),k=1,ns)
c   3 continue
   90 stop
      end
c


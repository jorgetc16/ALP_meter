************************************************************************
*   Authors: Jorge Martin Camalcih, Jorge Terol Calvo, Andrés Castillo *
*                                                                      *
*   Program that calculates the bootstrap of a time series, calculates *
*   the generalized Lomb-Scargle periodogram and identifies the        *
*   significant peaks. The data is characterized by a vector of times  *
*   tobs(ic), a vector of cv's of observations xobs(ic), a vector of   *
*   errors of the observations sxobs(ic) and an integer nobs with the  *
*   number of total obervations and length of the data set.            *
*                                                                      *
************************************************************************


      program Peak Finder
      implicit none
      double precision tobs,yobs,syobs ! Input data vectors, time series, pol. angle and 1sigma error of the pol. angle
      integer npuls1, nobs ! Number of data points
      integer ntot,n0 ! Total number of iterations and sampling density
      integer nmc !Total number MCs
      integer ic,jc ! Loop counters
      integer n68,n95,n99 ! Indices for 68% and 95% upper limit
      double precision dumb,pi,t0 ! Ancillary variable
      real hist,shist,h0,h1,h2,sh0,sh1,dh,dsh ! Histogram variables (must be real!)
      real histaux,shistaux
      integer ihist,nhist ! Histogram variables
      real ysint,sysint
      double precision W,wi,nu,dnu,numin,numax,ysint2
      double precision pLS,pLS68,pLS95,pLS99,pmax
!      double precision pLST,pLSaux,pLS68T,pLS95T,pLS99T ! OPTIONAL variables to study the nu distribution of the peaks

      parameter (pi=3.141592654d0)
      parameter (n0=5,nhist=20,nmc=1000)
      parameter (npuls1=100) !Number of data points of each pulsar
      dimension tobs(10000),yobs(10000),syobs(10000)
      dimension ysint(10000),sysint(10000),ysint2(10000)
      dimension wi(10000)
      dimension hist(1000),shist(1000),histaux(1000),shistaux(1000)
      dimension pmax(100000)
c      dimension pLST(160000,1000),pLSaux(160000),pLS68T(100000), ! OPTIONAL !
c     &pLS95T(100000),pLS99T(100000)


*************************************************************************
*                                                                       *
*     pulsar desired: list above                                        *
*                                                                       *
*************************************************************************
      nobs=npuls1
      open(1,file='InputData/XXXXX.dat',
     &status='old') ! data time series
      open(2,file='Output/Files/pLS_XXXXX.dat') !We create the file where we will record the values for the LS periodogram
      open(3,file='Output/Files/pLS_peaks_XXXXX.dat') !We create the file where we will record the peaks that are over the FAP above 95% C.L.
      open(4,file='Output/Files/pLS_FAP_XXXXX.dat') !We create the file where we will record the FAP values at 68%, 95% and 99% C.L. 
c      open(5,file='OutputData/pLS_boots_nu0.dat') !OPTIONAL!



*************************************************************************

      numax=20.d0 ! days^-1 maximum frequency
      n68=0.68*nmc
      n95=0.95*nmc
      n99=0.99*nmc

      do ic=1,nmc
        pmax(ic)=0.d0
      enddo
*
*   > Calculate the bootstrap
*
*
*        >> Histogram
*
      do ic=1,nhist
        hist(ic)=0.e0
      enddo


      h0=1.e6
      h1=-1.e6
      h2=1.e6
      sh0=1.e6
      sh1=-1.e6

      
      do ic=1,nobs
        read (1,*) tobs(ic),yobs(ic),syobs(ic)
        if (ic.eq.1) t0=tobs(ic)
        tobs(ic)=tobs(ic)-t0
        if (yobs(ic).le.h0) h0=yobs(ic)
        if (yobs(ic).ge.h1) h1=yobs(ic)
        if (syobs(ic).le.sh0) sh0=syobs(ic)
        if (syobs(ic).ge.sh1) sh1=syobs(ic)
      enddo

      !We calculate the difference between all the t_obs and the next one and take the minimum
      do ic=1,nobs-1
            dumb=tobs(ic+1)-tobs(ic)
            if (dumb.le.h2) h2=dumb
      enddo
      numax = 1/h2!maximum frequency
      numin=1.d0/tobs(nobs) !Minimum frequency
      dnu=1.d0/dfloat(n0)/tobs(nobs) !Frequency resolution
      ntot=int(numax/dnu) !Total number of frequencies

      !We see wether sho and sh1 have the same value to determine if all the errors are the same
      if (sh1.eq.sh0) then
            dh=(h1-h0)/dfloat(nhist)
            
            do ic=1,nobs
              ihist=(yobs(ic)-h0)/dh
              hist(ihist)=hist(ihist)+1.e0
            enddo
      
*
*       >> Calculate the cumulative distribution
*
            call RNHPRE(hist,nhist)
            do ic=1,nhist
              histaux(ic)=hist(ic)
            enddo
*
*   > Bootstraps
*
      
            do jc=1,nmc
              do ic=1,nobs
                  call RNHRAN(histaux,nhist,h0,dh,ysint(ic))

              enddo
      
              do ic=1,nobs
                  wi(ic)=1.d0
              enddo
              do ic=1,ntot+1
                  nu=numin+dnu*dfloat(ic-1)
                  dumb=pLS(nu,nobs,tobs,ysint2,wi)
                  if (dumb.ge.pmax(jc)) pmax(jc)=dumb
      !            pLST(ic,jc)=dumb !OPTIONAL!
              enddo
            enddo            
      else
            do ic=1,nhist
                  shist(ic)=0.e0
            enddo
            dh=(h1-h0)/dfloat(nhist)
            dsh=(sh1-sh0)/dfloat(nhist)
            
            do ic=1,nobs
              ihist=(yobs(ic)-h0)/dh
              hist(ihist)=hist(ihist)+1.e0
              ihist=(syobs(ic)-sh0)/dsh
              shist(ihist)=shist(ihist)+1.e0
            enddo
      
*
*       >> Calculate the cumulative distribution
*
            call RNHPRE(hist,nhist)
            call RNHPRE(shist,nhist)
            do ic=1,nhist
              histaux(ic)=hist(ic)
              shistaux(ic)=shist(ic)
            enddo
*
*   > Bootstraps
*
      
            do jc=1,nmc
              do ic=1,nobs
                  call RNHRAN(histaux,nhist,h0,dh,ysint(ic))
                  call RNHRAN(shistaux,nhist,sh0,dsh,sysint(ic))
                  if (sysint(ic).lt.0.d0) then
                      print *, "stop!"
                      stop
                  endif
              enddo
      
              W=0.d0
              do ic=1,nobs
                  wi(ic)=1.d0/sysint(ic)**2
                  W=W+wi(ic)
                  ysint2(ic)=ysint(ic)
              enddo
              do ic=1,nobs
                  wi(ic)=wi(ic)/W
              enddo
              do ic=1,ntot+1
                  nu=numin+dnu*dfloat(ic-1)
                  dumb=pLS(nu,nobs,tobs,ysint2,wi)
                  if (dumb.ge.pmax(jc)) pmax(jc)=dumb
      !            pLST(ic,jc)=dumb !OPTIONAL!
              enddo
            enddo
      end if
*
*   > Sort the pmax distribution and pick FPAs
*
      call SORTD(pmax,1,nmc,1)

      pLS68=pmax(n68)
      pLS95=pmax(n95)
      pLS99=pmax(n99)

      write(4,*) pLS68,pLS95,pLS99

*
*   > Next is OPTIONAL and allows to see nu distribution of maximum peaks
*
c      do ic=1,ntot+1
c        do jc=1,nmc
c            pLSaux(jc)=pLST(ic,jc)
c        enddo
c        call SORTD(pLSaux,1,nmc,1)
c        nu=numin+dnu*dfloat(ic-1)
c        pLS68T(ic)=pLSaux(n68)
c        pLS95T(ic)=pLSaux(n95)
c        pLS99T(ic)=pLSaux(n99)
c        write(5,*) nu,pLS68T(ic),pLS95T(ic),pLS95T(ic),pLSaux(nmc)
c      enddo

*
*   > Calculate data periodogram and identify peaks above with 5% FPA
*
      W=0.d0
      do ic=1,nobs
        wi(ic)=1.d0/syobs(ic)**2
        W=W+wi(ic)
      enddo
      do ic=1,nobs
        wi(ic)=wi(ic)/W
      enddo

      do ic=1,ntot+1
        nu=numin+dfloat(ic-1)*dnu
        dumb=pLS(nu,nobs,tobs,yobs,wi)
        write(2,*) nu,dumb
        if (dumb.ge.pLS95) write (3,*) nu,dumb
      enddo

      close(2)
      close(3)
      close(4)
c      close(5) !OPTIONAL!

      stop
      end

************************************************************************
*                                                                      *
*   > Calculating Generalized Lomb-Scargle periodogram                 *
*                                                                      *
************************************************************************


      function pLS(nu,nobs,tobs,yobs,wi)
      implicit none
      double precision nu,pLS,pi
      double precision tobs,yobs,wi,cosi,sini
      double precision Yi,Ci,Si,YYdi,YCdi,YSdi,CCdi,SSdi,CSdi
      double precision YY,YC,YS,CC,SS,CS,Dn
      integer ic,nobs
      parameter (pi=3.141592654d0)

      dimension tobs(10000),yobs(10000),wi(10000)

      Ci=0.d0
      Si=0.d0
      Yi=0.d0
      YYdi=0.0d0
      YCdi=0.d0
      YSdi=0.d0
      CCdi=0.d0
      SSdi=0.d0
      CSdi=0.d0

      do ic=1,nobs
        Yi=Yi+wi(ic)*yobs(ic)
        cosi=dcos(2.d0*pi*nu*tobs(ic))
        sini=dsin(2.d0*pi*nu*tobs(ic))
        Ci=Ci+wi(ic)*cosi
        Si=Si+wi(ic)*sini
        YYdi=YYdi+wi(ic)*yobs(ic)**2
        YCdi=YCdi+wi(ic)*yobs(ic)*cosi
        YSdi=YSdi+wi(ic)*yobs(ic)*sini
        CCdi=CCdi+wi(ic)*cosi**2
        CSdi=CSdi+wi(ic)*cosi*sini
c        CSdi=CSdi+wi(ic)*dsin(2.d0*pi*nu*tobs(ic))*
c     &  dcos(2.d0*pi*nu*tobs(ic))
      enddo
      SSdi=1.d0-CCdi
      YY=YYdi-Yi**2
      YC=YCdi-Yi*Ci
      YS=YSdi-Yi*Si
      CC=CCdi-Ci**2
      SS=SSdi-Si**2
      CS=CSdi-Ci*Si
      Dn=CC*SS-CS**2

      pLS=(SS*YC**2+CC*YS**2-2.d0*CS*YC*YS)/(Dn*YY)

      end
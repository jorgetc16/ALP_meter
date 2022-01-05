************************************************************************
*   Authors: Jorge Martin Camalcih, Andrés Castillo, Jorge Terol Calvo *
*                                                                      *
*   Program that calculates the bootstrap of a time series, calculates *
*   the generalized Lomb-Scargle periodogram and identifies the        *
*   significant peaks. The data is characterized by a vector of times  *
*   tobs(ic), a vector of cv's of observations xobs(ic), a vector of   *
*   errors of the observations sxobs(ic) and an integer nobs with the  *
*   number of total obervations and length of the data set.            *
*                                                                      *
************************************************************************


      program phi_limits
      implicit none
      double precision tobs,yobs,syobs, pi ! Input data vectors, time series, pol. angle and 1sigma error of the pol. angle
      integer npuls1,npuls2,npuls3,npuls4,npuls5,npuls6,npuls7,nobs ! Input number of observations/length of vectors (list above)
      integer npuls8,npuls9,npuls10,npuls11,npuls12,npuls13,npuls14
      integer npuls15,npuls16,npuls17,npuls18,npuls19,npuls20,npuls21
      integer ntot,n0 ! Total number of iterations and sampling density
      integer nmc ! Total number MCs
      integer nsc ! Numbers of points in the scan of nu
c      double precision Tmax,Tmin,logTmax,logTmin,dlogT,logT ! Scan in Period
      double precision lognumax,lognumin,dlognu,lognu ! Scan in frequency
      integer ic,jc,kc ! Loop counters
      double precision dumb,dumbb,t0 ! Ancillary variable
      real hist,shist,h0,h1,sh0,sh1,dh,dsh ! Histogram variables (must be real!)
      real histaux,shistaux
      integer ihist,nhist ! Histogram variables
      double precision W,wi,nu,dnu,numin,numax
      double precision pLS,pLS95,pLS99
      double precision pLSdata,phi95f,phi95
      double precision dzerox
!      double precision pLST,pLSaux,pLS68T,pLS95T,pLS99T ! OPTIONAL variables to study the nu distribution of the peaks
      double precision phi0,phi1,eps
      integer maxf

      common /hist/histaux,shistaux,nhist,h0,dh,sh0,dsh
      common /datas/tobs,nobs
      common /mc/nu,pLSdata,nmc

      parameter (n0=1,nsc=1)
      parameter (npuls1=393,npuls2=125,npuls3=47,npuls4=84,npuls5=78,
     &npuls6=110,npuls7=103,npuls8=101,npuls9=98,npuls10=112,npuls11=70,
     &npuls12=22,npuls13=104,npuls14=74,npuls15=59,npuls16=173,
     &npuls17=74,npuls18=100,npuls19=62,npuls20=67,npuls21=1099)


      dimension tobs(10000),yobs(10000),syobs(10000)
      dimension wi(10000)
      dimension hist(1000),shist(1000),histaux(1000),shistaux(1000)

      parameter (pi=3.141592654d0)

      external phi95f

*************************************************************************
*                                                                       *
*     pulsar desired: list above                                        *
*                                                                       *
*************************************************************************
      nobs=npuls1
      open(1,file='InputData/XXXXX.dat'
     & ,status='old') ! Data time series
      open(2,file='Output/Files/phi95_XXXXX.dat')

*************************************************************************
*************************************************************************
   
      nhist=20
      
      nmc=3000
      numax=3.0081532001495361 ! days^-1 maximum frequency
      phi0=0.d0 ! minimum possible value of phi
      phi1=10.d0 ! maximum possible value of phi
      eps=1.d-4
      maxf=30


*68
*
*        >> Histogram
*
      do ic=1,nhist
        hist(ic)=0.e0
        shist(ic)=0.e0
      enddo


      h0=1.e6
      h1=-1.e6
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


      numin=3.0081532001495361!Minimum frequency
      dnu=1.d0!Frequency resolution
      ntot=1!Total number of frequencies

      dh=(h1-h0)/dfloat(nhist)
      dsh=(sh1-sh0)/dfloat(nhist)
      !print *,sh0,sh1,dsh
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
      call RNHPRE(shist,nhist)      ! CERNLIB: Subroutines to initialize RNHRAN in function phi95f(x)
      do ic=1,nhist
        histaux(ic)=hist(ic)
        shistaux(ic)=shist(ic)
      enddo
*
*   > Signal + given MC
*

c      Tmin=1.d0/numax              ! Scan in period
c      Tmax=1.d0/numin
c      logTmax=dlog10(Tmax)
c      logTmin=dlog10(Tmin)
c      dlogT=(logTmax-logTmin)/dfloat(nsc)
      lognumax=dlog10(numax)
      lognumin=dlog10(numin)
      dlognu=(lognumax-lognumin)/dfloat(nsc)

c      nu=3.80077d0                 ! Block to test a particular value of the periodogram data
c      W=0.d0
c      do ic=1,nobs
c        wi(ic)=1.d0/syobs(ic)**2
c        W=W+wi(ic)
c      enddo
c      do ic=1,nobs
c        wi(ic)=wi(ic)/W
c      enddo
c
c      pLSdata=pLS(nu,nobs,tobs,yobs,wi)
c      print *,pLSdata
c      stop
      do kc=1,nsc+1
        lognu=lognumin+dlognu*dfloat(kc-1)
        nu=10.d0**lognu
        W=0.d0
        do ic=1,nobs
            wi(ic)=1.d0/syobs(ic)**2
            W=W+wi(ic)
        enddo
        do ic=1,nobs
            wi(ic)=wi(ic)/W
        enddo
        
        pLSdata=pLS(nu,nobs,tobs,yobs,wi)

        dumb=phi95f(phi0)
        dumbb=phi95f(phi1)
       if ((dumb*dumbb).ge.0.d0) then
            print *,"Crap=",nu,pLSdata
            write(2,*) nu,''
            cycle
        endif
        phi95=dzerox(phi0,phi1,eps,maxf,phi95f,1)  ! CERNLIB: Subroutine to find a zero of an external function
        print *, nu,phi95
        write(2,*) nu,phi95
      enddo
      close(2)
      stop
      end


************************************************************************
*                                                                      *
*   Function calculating 95% range of phi for a given nu               *
*                                                                      *
************************************************************************

      function phi95f(phi)
      implicit none
      double precision phi95f,phi,pi
      integer ic,jc,kc,nmc,nobs,n95, n68, n50
      real histaux,shistaux,h0,sh0,dh,dsh ! Histogram variables (must be real!)
      integer nhist ! Histogram variables
      real ysint,sysint
      double precision W,wi,nu,ysint2,tobs
      double precision pLS,pLSphi,pLSdata
      real delta

      common /hist/histaux,shistaux,nhist,h0,dh,sh0,dsh
      common /datas/tobs,nobs
      common /mc/nu,pLSdata,nmc
      parameter (pi=3.141592654d0)

      dimension tobs(10000),ysint(10000),sysint(10000),ysint2(10000)
      dimension wi(10000)
      dimension histaux(1000),shistaux(1000)
      dimension pLSphi(100000)


      n95=0.05*nmc
      n68=0.32*nmc
      n50=0.50*nmc

      call RANLUX(delta,1) ! CERNLIB subroutine: Calls random number [0,1]
      delta=2.d0*pi*delta ! the phase delta is a random number between 0 and 2pi

      do jc=1,nmc
        do ic=1,nobs
            call RNHRAN(histaux,nhist,h0,dh,ysint(ic))
            call RNHRAN(shistaux,nhist,sh0,dsh,sysint(ic)) ! CERNLIB subroutines: Calls random number with specified distro
            if (sysint(ic).lt.0.d0) then
                print *, "stop!"
                stop
            endif
        enddo

        W=0.d0
        do ic=1,nobs
            ysint(ic)=ysint(ic)+phi*dcos(2.d0*pi*nu*tobs(ic)+delta) !We add the periodic signal to the sintetic data
            wi(ic)=1.d0/sysint(ic)**2
            W=W+wi(ic)
            ysint2(ic)=ysint(ic)
        enddo
        do ic=1,nobs
            wi(ic)=wi(ic)/W
        enddo
         pLSphi(jc)=pLS(nu,nobs,tobs,ysint2,wi)
c        print *,pLSphi(jc)
      enddo
      call SORTD(pLSphi,1,nmc,1) ! Subroutine to sort lists
      
      phi95f=pLSphi(n95)-pLSdata

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

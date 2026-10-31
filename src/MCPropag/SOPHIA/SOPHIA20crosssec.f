c*****************************************************************************
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c**!!              IF YOU USE THIS PROGRAM, PLEASE CITE:                 !!***
c**!! A.M"ucke, Ralph Engel, J.P.Rachen, R.J.Protheroe and Todor Stanev, !!***
c**!!  1999, astro-ph/9903478, to appear in Comp.Phys.Commun.            !!***
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c*****************************************************************************
c** Further SOPHIA related papers:                                         ***
c** (1) M"ucke A., et al 1999, astro-ph/9808279, to appear in PASA.        ***
c** (2) M"ucke A., et al 1999, to appear in: Proc. of the                  ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (3) M"ucke A., et al 1999, astro-ph/9905153, to appear in: Proc. of    ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (4) M"ucke A., et al 1999, to appear in: Proc. of 26th Int.Cosmic Ray  ***
c**      Conf. (Salt Lake City, Utah)                                      ***
c*****************************************************************************



       program SOPHIA20

c************************************************************************
c*** Simulation Of PhotoHadronic Interactions in Astrophysics ***********
c***           VERSION 2.0                                    ***********
c************************************************************************

c************************************************************************
c** Main program for photopion production of relativistic nucleons     **
c**  in a radiation field (blackbody or power law)                     **
c************************************************************************
c** Date: 20/01/98          **
c** correct.:19/02/98       **
c** first release to        ** 
c**  collab.:      25/05/98 **
c** Version 1.3:   31/08/98 **
c** Version 1.4:   12/10/98 **
c** change to DP:  16/11/98 **
c** last corr.:       07/99 **
c** authors: A.G.F. Muecke  **
c**          R.R. Engel     **
c**   in collaboration with **
c**         R.J. Protheroe  **
c**         J.P. Rachen     **
c**         T.S. Stanev     **
c*****************************

c****** INPUT ***********************************************************************
c E0 = energy of incident proton (in lab frame) [in GeV]
c L0 = code number of the incident nucleon (L0=13: proton, L0=14: neutron)
c soft photon field:   for thermal spectrum: temperature tbb [in K]
c                      for PL: (set tbb = 0) alpha = PL index (S \sim nu^{-alpha}), 
c                                            epsmin = low energy cut off  [in eV]
c                                            epsmax = high energy cut off [in eV]
c****** OUTPUT **********************************************************************
c** energy distribution P(x) (logarithmic scale, x=E_particle/E0) for ***
c** photons, protons, neutrons, e-neutrinos, nu-neutrinos             ***
c************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       COMMON/input/ tbb,E0,alpha1,alpha2,
     &           epsm1,epsm2,epsb,L0
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_MASS1/ AM(49), AM2(49)
       COMMON /S_CHP/  S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)

      CHARACTER*6 NAMPRES
      COMMON /RES_PROP/ AMRES(9), SIG0(9),WIDTH(9), 
     +                    NAMPRES(0:9)

      CHARACTER*6 NAMPRESp
      COMMON /RES_PROPp/ AMRESp(9), BGAMMAp(9),WIDTHp(9),  
     +                    RATIOJp(9),NAMPRESp(0:9)

      CHARACTER*6 NAMPRESn
      COMMON /RES_PROPn/ AMRESn(9), BGAMMAn(9),WIDTHn(9),  
     +                    RATIOJn(9),NAMPRESn(0:9)


       DIMENSION DNg(201),DNnum(201),DNnue(201),DNp(201),DNn(201)
       DIMENSION DNnuma(201),DNnuea(201),E0_arr(101)
       DIMENSION Dg(101,201),Dnum(101,201),Dnue(101,201)
       DIMENSION Dp(101,201),Dn(101,201),Dnuma(101,201),Dnuea(101,201)
       DIMENSION Dem(101,201),Dep(101,201),DNem(201),DNep(201)

       dimension c_feps(10000),eps_arr(10000)

       character*6 nameinc
       character*1 ans

       external sample_eps,sample_s,eventgen,
     &  listdistr,output,initial,prob_epskt

       DATA pi /3.141593D0/ 

        do i=1,201
          DNg(i) = 0.D0
          DNnum(i) = 0.D0
          DNnue(i) = 0.D0
          DNp(i) = 0.D0
          DNn(i) = 0.D0
          DNem(i) = 0.D0
          DNep(i) = 0.D0
          DNnuma(i) = 0.D0
          DNnuea(i) = 0.D0

          do jm=1,101
            Dg(jm,i) = 0.D0
            Dnum(jm,i) = 0.D0
            Dnue(jm,i) = 0.D0
            Dp(jm,i) = 0.D0
            Dn(jm,i) = 0.D0
            Dem(jm,i) = 0.D0
            Dep(jm,i) = 0.D0
            Dnuma(jm,i) = 0.D0
            Dnuea(jm,i) = 0.D0
         enddo

       enddo


c****** INPUT **************************************************
       print*
       print*,'Give in code number of incident nucleon: '
       print*,' (13 = proton, 14 = neutron)'
       read(*,*) L0

       delE = 0.D0
c       print*
c       print*,'Incident nucleon spectrum:'
c       print*,'energy grid [s] or single nucleon [n] ? '
c       read(*,'(A1)') ans

c       if (ans.eq.'n') then
        jm = 1
        ninc = 1
c 4     print*,'Give energy [in GeV] of incident nucleon: '
c       read(*,*) E0

c       if (E0.lt.AM(L0)) then
c       print*,'No valid input !'
c       goto 4
c       endif
       E0=1.D11
       Emin = log10(E0)

       Emax=Emin

c 5     print*,'Give s [in GeV^2]: '
c       read(*,*) s

c       if (s.lt.1.1646D0) then
c        print*,'No valid input (s>1.1646)!'
c        goto 5
c       endif

c 6     print*,'Give cos(theta): '
c       read(*,*) theta

c       if (theta.lt.-1) then
c        print*,'No valid input !'
c        goto 6
c       endif

c       if (theta.gt.1) then
c        print*,'No valid input !'
c        goto 6
c       endif

c       endif
c       if (ans.eq.'s') then
c  5    print*,'Give in low-energy cut off (in GeV) of nucleon energy grid:'
c       read(*,*) Emin

c       if (Emin.lt.AM(L0)) then
c        print*,'No valid input !'
c        goto 5
c       endif
c       Emin = log10(Emin)
c       print*,'Give in high-energy cut off (in GeV) of nucleon energy grid:'
c       read(*,*) Emax
c
c       if (Emax.lt.Emin) then
c        print*,'No valid input !'
c        goto 5
c       endif
c       Emax = log10(Emax)
c       print*,'Give number of bins (< 101):'
c       read(*,*) ninc
c       delE = (Emax-Emin)/ninc
c       ninc = ninc+1
c       endif

c       print*
c  2    print*,'Give soft photon spectrum: '
c       print*,' blackbody spectrum ? (y/n): '
c       read(*,'(A1)') ans
c
c       if (ans.eq.'y') then
c         rad_den = 1.D0
c         print*,'Give temperature [in K]: '
c         read(*,*) tbb
c
c       else if (ans.eq.'n') then
c       tbb = 0.D0
c       print*,' for power law spectrum n(eps) ~ eps^-alpha:'
c       print*,'  (alpha = 0 ... 3)'
c  8    print*,'straight power law ? (y/n): '

c       read(*,'(A1)') ans
c       if (ans.eq.'y') then
c  3    print*,'   spectral index alpha = '
c       read(*,*) alpha2

c        if (alpha2.lt.0.D0.or.alpha2.gt.3.D0) then
c         print*,'no valid input !'
c         goto 3
c        endif
c       alpha1 = alpha2
c       else if (ans.eq.'n') then
c        print*,'broken power law:'
c        print*,' spectral index alpha1, alpha2 = '
c        read(*,*) alpha1,alpha2
c        if (alpha1.lt.0.D0.or.alpha1.gt.3.D0.or.
c     &  alpha2.lt.0.D0.or.alpha2.gt.3.D0) then
c         print*,'no valid input !'
c         goto 8
c        endif
c        print*,'break energy [eV] = '
c        read(*,*) epsb        
c       else
c         print*,'no valid input !'
c         goto 8
c       endif 
c       print*,'low-energy cut off [eV] = '
c       read(*,*) epsmin
c       print*,'high-energy cut off [eV] = '
c       read(*,*) epsmax

c       if (ans.eq.'y') epsb = epsmin/100.D0
c       else
c        print*,'no valid input !'
c        goto 2
c       endif
      
c       print*
c       print*,'Give number of trials: '
c       read(*,*) ntrial
       ntrial=2
       nbins=10
c       print*
c       print*,
c     &'Give number of bins (<201) for output particle spectra: '
c       print*,
c     &'(spectra in logarithmic equal bins, with stepsize Delta x'
c       print*,'  with x=E_particle/E_initial nucleon )'
c       read(*,*) nbins
c
      delx=0.1
c       print*
c       print*,
c     &'Give stepsize Delta x for output particle spectra: '
c       read(*,*) delx

      nameinc='aaa'
c       print*
c       print*,'Give filename (< 7 CH) for output: '
c       read(*,'(A6)') nameinc
c       print*

       call initial(L0)
c*******************************************************
c calculate epsilon = energy of photon in lab frame
c*******************************************************
       pm = AM(L0)
       gammap = E0/pm
       betap = sqrt(1.D0-1.D0/gammap/gammap)
       eps = (s-pm*pm)*0.5D0/E0/(1.D0-betap*theta)
       theta = acos(theta)*180.D0/pi 

c... nucleon energy loop:
c       do jm=1,ninc
         Elog = Emin+(jm-1)*delE
         E0 = 10.D0**Elog
         E0_arr(jm) = E0

c... trial loop:
       ntrial=280
       do nt=0,ntrial 

c  6    call sample_eps(epseV,epsmin,epsmax)
c       if (epseV.le.0.) goto 7
c *** eps in GeV:
c       eps = epseV/1.D9
c       Etot = E0+eps
       s=1.1646D0*10.D0**(0.05D0*nt)
c       if (s.gt.1.1646D0) then
        eps_prime=(s-pm*pm)/pm*0.5D0
        sig = crossection(eps_prime,3,L0)
c        write (6,922) '   vot: ',nt,eps_prime,s,sig
        WRITE(*,*) eps_prime,s,sig
c       endif
       
       enddo
  7    continue

c... nucleon energy loop:
c       enddo

c       call output(Dg,Dnum,Dnuma,Dnue,Dnuea,Dp,Dn,Dem,Dep,
c     &  nbins,ninc,nameinc,delx,Emin,Emax,E0_arr,epsmin,epsmax)

       STOP
       

 922     FORMAT(7X,A8,20X,2(F9.3,2X),2(E9.3,2X))
       END

       subroutine sample_photopion(L0,E0,eps,theta,iSec,energies,nTypes)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE
       DIMENSION energies(2000)
       DIMENSION nTypes(2000)
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb

       call eventgen(L0,E0,eps,theta,Imode)
       iSec=NP
       do i=1,NP
        nTypes(i) = abs(LLIST(i))
        energies(i) = abs(P(i,4))
       enddo

       RETURN
       END

       subroutine crossec(L0,eps_prime,sig)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE

       sig = crossection(eps_prime,3,L0)

       RETURN
       END

       subroutine sample_photopion_rel(L0,eps_prime,iSec,relEnergies,nTypes)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE
       DIMENSION relEnergies(2000)
       DIMENSION nTypes(2000)
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_MASS1/ AM(49), AM2(49)

c      Using proton rest frame
       E0=AM(L0)
c      we assume that lab frame speed is directed along z-axes
c      than with good accuracy initial photon momentum is directed
c      against z-axis
       theta=0

       call eventgen(L0,E0,eps_prime,theta,Imode)
       iSec=NP
       do i=1,NP
        nTypes(i) = abs(LLIST(i))
c relEnergies(i) - energy of final particle in units of
c initial nucleon energy in lab frame. Here we use the fact that nucleon
c rest frame speed in lab frame is directed along z-axis (against photon
c momentum in rest frame) with good accuracy and assume
c ultrarelativistic case: beta=1 and so the fraction doesn't depend on
c the initial nucleon gamma factor.
c The precise expression would be
c       relEnergies(i) = abs(P(i,4)+beta*P(i,3))/AM(L0)
c where beta is speed of initial nucleon in the rest frame
        relEnergies(i) = abs(P(i,4)+P(i,3))/AM(L0)
       enddo

       RETURN
       END

       subroutine get_mass(L0,pm)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       SAVE
       COMMON /S_MASS1/ AM(49), AM2(49)
       pm = AM(L0)
       RETURN
       END

       subroutine set_random_seed(iseed)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       COMMON/LUDATR/MRLU(6),RRLU(100)
       SAVE
       MRLU(1) = iseed
       RETURN
       END
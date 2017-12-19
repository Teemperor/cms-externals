      REAL*8 FUNCTION DISTJKWK(ID,S,T,H1,H2,KEY)
C
C     INFO: DISTJKWK is used through  nonSM_adopt which is located  
C           in tau-reweight-test.cxx this is the only reference.
C
C***************************************************************************
C* CALCULATES DSIGMA/DT FOR QQBAR -> TAU+TAU- VIA S-CHANNEL GAMMA+Z+SPIN2
C        ID=1 FOR UP-TYPE QQBAR, ID=2 FOR DOWN-TYPE QQBAR
C        S -- CM QQBAR ENERGY^2
C        T -- cosTheta (scattering angle of tau- with respect to quark
C        H1 - HELICITY OF TAU-:  1=R,   -1=L
C        H2 - HELICITY OF TAU+   is always opposite to TAU-
C        KEY=0  SM,   KEY=1 SM+SPIN2
C***************************************************************************

C***************************************************************************
C*    For explanations see   S. Banerjee et al  CERN-PH-TH-2012-347        *
C***************************************************************************

      IMPLICIT REAL*8(A-H,M,O-Z)
      INTEGER H1,H2,K
      SAVE K


      IF (K.lT.10) THEN 
       K=K+1
       IF(KEY.EQ.0) THEN  
            WRITE(*,*) 'SM' 
       ELSE IF(KEY.EQ.1) THEN 
            WRITE(*,*) 'SM+nonSM SPIN-2 state'
       ELSE
       WRITE(*,*) 'WRONG KEY'
       STOP
       ENDIF
      ENDIF

      cost=T
     
C*      ==========  Standard Model Parameters  ==========
      alphaEM=1.d0/128.d0      ! at mZ
      GF=1.1667D-5             ! Fermi constant
      S2W = .2315D0            ! sin^2(theta_W)
      C2W = 1.D0-S2W 
      CW  = DSQRT(1.D0-S2W)
      SW  = DSQRT(S2W)
      F=1000.d0                ! scale of spin-2 physics ????

      MZ = 91.187D0     ! Z mass
      wz = 2.d0         ! Z width
      mh=125.d0         ! spin-2 mass            choose your own       
      wh=1.5d0          ! spin-2 width           choose your own
      MW = MZ*CW
      Pi=4.d0*datan(1.d0)
      e=dsqrt(4.d0*Pi*alphaEM) 

c     couplings to tau lepton
c*****************************
c  photon couplings to tau
      gttL=-1.d0
      gttR=-1.d0
c  Z couplings to tau
      ZttL= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(-0.5D0+S2W)
      ZttR= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(       S2W)
c  spin-2 couplings to tau
      httL=1.0d0           !    RANGE????
      httR=1.0d0           !    RANGE????  


c  couplings to quarks
C **************************
      IF (ID.EQ.1) THEN  !  couplings to UP-TYPE QUARKS
c  z couplings
      ZqqL= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(0.5D0-2.D0*S2W/3.D0)
      ZqqR= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(     -2.D0*S2W/3.D0)
c  photon couplings   
      gqqL=2.d0/3.d0
      gqqR=2.d0/3.d0
c   spin-2 couplings
      hqqL=1.0d0         !     RANGE????
      hqqR=1.0d0         !     RANGE????
      ELSE IF (ID.EQ.2) THEN ! couplings to DOWN-TYPE QUARKS
c  z couplings
      ZqqL= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(-0.5D0+S2W/3.D0)
      ZqqR= DSQRT(8*GF*MZ**2/DSQRT(2.D0))*(       S2W/3.D0)
c  photon couplings
      gqqL=-1.d0/3.d0
      gqqR=-1.d0/3.d0
c  spin-2 couplings
      hqqL=1.0d0         !     RANGE?????
      hqqR=1.0d0         !     RANGE?????
      ELSE

      WRITE(*,*) 'WRONG ID'
      STOP   
      ENDIF

      IF(H1.EQ.H2) THEN
      distJKWK=0.d0
      return 
      ENDIF

C     SELECT DEMANDED HELICITY OF TAU
      IF(H1.EQ.1) THEN                             ! R-HANDED TAU-
      gttL=0.d0
      ZttL=0.d0
      httL=0.d0
      ELSE IF(H1.EQ.-1) THEN                       ! L-HANDED TAU-
      gttR=0.d0
      ZttR=0.d0
      httR=0.d0
      ELSE
      WRITE(*,*) 'WRONG H1'
      STOP
      ENDIF 



      distJKWK=
C PHOTON^2
     -    (alphaEM**2*(gqqL**2*                     
     -        ((1 + cost)**2*gttL**2 +
     -          (-1 + cost)**2*gttR**2) +
     -       gqqR**2*((-1 + cost)**2*gttL**2 +
     -          (1 + cost)**2*gttR**2))*Pi)/(24.*s) +
C PHOTON-Z interference
     -  (e**2*(-mZ**2 + s)*
     -     (gqqL*ZqqL*((1 + cost)**2*gttL*ZttL +
     -          (-1 + cost)**2*gttR*ZttR) +
     -       gqqR*ZqqR*
     -        ((-1 + cost)**2*gttL*ZttL +
     -          (1 + cost)**2*gttR*ZttR)))/
     -   (192.*Pi*((mZ**2 - s)**2 + mZ**2*wZ**2)) +
C Z^2
     -  (s*(ZqqL**2*((1 + cost)**2*ZttL**2 +
     -          (-1 + cost)**2*ZttR**2) +
     -       ZqqR**2*((-1 + cost)**2*ZttL**2 +
     -          (1 + cost)**2*ZttR**2)))/
     -   (384.*Pi*((mZ**2 - s)**2 + mZ**2*wZ**2))
c     the spin-2 contribution below
     - + KEY*(
c (spin-2)^2
     -  ((hqqL**2*((-1 + cost + 2*cost**2)**2*httL**2 +
     -         (1 + cost - 2*cost**2)**2*httR**2) +
     -      hqqR**2*((1 + cost - 2*cost**2)**2*
     -          httL**2 +
     -         (-1 + cost + 2*cost**2)**2*httR**2))*s**3
     -    )/(24576.*Pi*((mh**2 - s)**2 + mh**2*wh**2))/F**4
c spin-2 - photon interference
     -  - (e**2*(gqqR*hqqR*
     -       ((-1 + cost)**2*(1 + 2*cost)*gttL*httL +
     -         (1 + cost)**2*(-1 + 2*cost)*gttR*httR)
     -       + gqqL*hqqL*
     -       ((1 + cost)**2*(-1 + 2*cost)*gttL*httL +
     -         (-1 + cost)**2*(1 + 2*cost)*gttR*httR))
     -     *(mh**2 - s)*s)/
     -  (1536.*Pi*((mh**2 - s)**2 + mh**2*wh**2))/F**2
c spin-2 - Z interference
     -   + (s**2*((mh**2 - s)*(mZ**2 - s) + mh*mZ*wh*wZ)*
     -     (hqqR*ZqqR*((-1 + cost)**2*(1 + 2*cost)*httL*
     -          ZttL +
     -         (1 + cost)**2*(-1 + 2*cost)*httR*ZttR)
     -       + hqqL*ZqqL*
     -       ((1 + cost)**2*(-1 + 2*cost)*httL*ZttL +
     -         (-1 + cost)**2*(1 + 2*cost)*httR*ZttR)))
     -   /(1536.*Pi*((mh**2 - s)**2 + mh**2*wh**2)*
     -    ((mZ**2 - s)**2 + mZ**2*wZ**2))/F**2     )  

      return
      END



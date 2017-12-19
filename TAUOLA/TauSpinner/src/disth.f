C this is variant SM-only. We need it inside TauSpinner `always' but may be
C it is not an appropriate place for it. 
    
      REAL*8 FUNCTION DISTH(S,T,H1,H2)
C***************************************************************************
C* CALCULATES DSIGMA/DT FOR GluonGluon -> TAU+TAU- VIA:
C            s-channel SM Higgs if KEY=0
C            s-channel spin2 object if KEY=1
C        S -- CM gluon-gluon ENERGY^2
C        T -- cosTheta (scattering angle of tau- with respect to gluon
C        H1 - HELICITY OF TAU-:  1=R,   -1=L
C        H2 - HELICITY OF TAU+   
C        KEY=0  SM,   KEY=1 SPIN2
C***************************************************************************

C***************************************************************************
C*                           *
C***************************************************************************

      IMPLICIT REAL*8(A-H,M,O-Z)
      INTEGER H1,H2

 
!      IF(KEY.EQ.0) THEN  
!           WRITE(*,*) 'SM Higgs' 
!      ELSE IF(KEY.EQ.1) THEN 
!           WRITE(*,*) 'nonSM not ready'
!           STOP
!      ELSE
!      WRITE(*,*) 'WRONG KEY'
!      STOP
!      ENDIF



      cost=T
     
C*      ==========  Standard Model Parameters  ==========
      alphaEM=1.d0/128.d0      ! at mZ
      alphaS=0.112d0           ! at mZ
      GF=1.1667D-5             ! Fermi constant
      S2W = .2315D0            ! sin^2(theta_W)
      C2W = 1.D0-S2W 
      CW  = DSQRT(1.D0-S2W)
      SW  = DSQRT(S2W)
      F=1000.d0                ! scale of spin-2 physics ????

      MZ = 91.187D0     ! Z mass
      wz = 2.d0         ! Z width
      mH = 125.d0       ! SM Higgs mass
      mT=173.d0         ! top quark mass
      mX=125.d0         ! spin-2 mass            choose your own       
      wX=10.d0          ! spin-2 width           choose your own
      MW = MZ*CW
      Pi=4.d0*datan(1.d0)
      e=dsqrt(4.d0*Pi*alphaEM) 
      Br=  0.065d0                    ! Br(H->tautau) 




      ! SM Higgs via top quark loop

        IF(H1.EQ.-H2) THEN
           DISTH=0.d0           ! for spin 0 Higgs
           return 
        ELSE
           tauT=mH**2/4/mT**2
           tauQ=tauT
           stauQ=DSQRT(tauQ)


           DISTH=GF*alphaS**2/288.d0/DSQRT(2.d0)/PI*DABS(3.d0/2d0*
     .     (tauQ+(tauQ-1.0)*(DASIN(stauQ))**2)/tauQ**2)**2
     .     *Br/2.d0

C  Br/2 because of FB symmetry
           return
        ENDIF
      




      return
      END

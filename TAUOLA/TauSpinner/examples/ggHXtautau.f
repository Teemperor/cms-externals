 
     
      REAL*8 FUNCTION DISTHJKWK(Shat,T,H1,H2,KEY)
C***************************************************************************
C* CALCULATES DSIGMA/DT FOR GluonGluon -> TAU+TAU- VIA:
C            s-channel SM Higgs if KEY=0
C            s-channel spin2 object if KEY=1
c without zero-width approximation
C        Shat -- CM gluon-gluon ENERGY^2
C        T -- cosTheta (scattering angle of tau- with respect to gluon
C        H1 - HELICITY OF TAU-:  1=R,   -1=L
C        H2 - HELICITY OF TAU+   
C        KEY=0  SM,   KEY=1 SPIN2
C***************************************************************************

C***************************************************************************
C*                           *
C***************************************************************************

      IMPLICIT REAL*8(A-H,M,O-Z)
      complex*16 ftau,aftau
      INTEGER H1,H2
      DATA K /0/
      SAVE K


C*      ==========  Standard Model Parameters  ==========
       alphaEM=1.d0/128.d0      ! at mZ
       alphaS=0.118d0           ! at mZ
       GF=1.1667D-5             ! Fermi constant
       S2W = .2315D0            ! sin^2(theta_W)
       C2W = 1.D0-S2W 
       CW  = DSQRT(1.D0-S2W)
       SW  = DSQRT(S2W)
       F=1000.d0                ! scale of spin-2 physics ????

       MZ = 91.187D0     ! Z mass
       wz = 2.d0         ! Z width
       mH = 125.d0       ! SM Higgs mass
       wH=  0.0041d0     ! total decay width of SM Higgs
       mtau=1.777d0      ! tau mass
       mT=173.d0         ! top quark mass
       mX=125.d0         ! spin-2 mass            choose your own       
       wX=wH             ! spin-2 width           choose your own
       ce= 0.1d0         ! spin-2 coupling        choose your own
       MW = MZ*CW
       Pi=4.d0*datan(1.d0)
       e=dsqrt(4.d0*Pi*alphaEM) 
       Br=  0.065d0                    ! Br(H->tautau) 

       IF (K.lT.10) THEN 
       K=K+1
        IF(KEY.EQ.0) THEN  
           WRITE(*,*) 'SM Higgs' 
        ELSE IF(KEY.EQ.1) THEN 
            WRITE(*,*) 'SM+nonSM SPIN-2 state'
        ELSE
        WRITE(*,*) 'WRONG KEY'
        STOP
        ENDIF
       ENDIF

  

      cost=T

      if(Shat.lt.4d0*mtau**2) then
       write(*,*) 'Shat too low'
       DISTHJKWK=0.0
       RETURN
      endif


c  spin-2 couplings to tau
      httL=0.1d0           !    RANGE????
      httR=0.1d0           !    RANGE????  


C     SELECT DEMANDED HELICITY OF TAU
c here the tau mass is neglected, i.e. helicity=chirality
      IF(H1.EQ.1) THEN                             ! R-HANDED TAU-
      httL=0.d0
      ELSE IF(H1.EQ.-1) THEN                       ! L-HANDED TAU-
      httR=0.d0
      ELSE
      WRITE(*,*) 'WRONG H1'
      STOP
      ENDIF 




      IF(KEY.EQ.0) THEN    ! SM Higgs via top quark loop

        IF(H1.EQ.-H2) THEN
           distHJKWK=0.d0           ! for spin 0 Higgs
           return 
        ELSE
           tauT=Shat/4./mT**2
           tauQ=tauT
c   only top in the loop included in the production gg->H
           stauQ=DSQRT(tauQ)
           if(tauQ.lt.1.d0) then
             ftau=(DASIN(stauQ))**2*DCMPLX(1d0,0d0)
           else
             ftau=-(dlog((1d0+dsqrt(1d0-1d0/tauQ))/
     .                    (1d0-dsqrt(1d0-1d0/tauQ)))*DCMPLX(1d0,0d0)
     .                          -Pi*DCMPLX(0d0,1d0))**2/4d0
           endif
           sigma0=GF*alphaS**2/288.d0/DSQRT(2d0)/PI*CDABS(3.d0/2d0*
     .     (tauQ+(tauQ-1.d0)*ftau)/tauQ**2)**2

c  higgs width via Breit-Wigner
c  mH**2*delta(shat-mH**2) -> 
      prop=mH**3/((shat-mH**2)**2+wH**2*mH**2)/pi

c  partial decay width H-> tautau  
      gamtau=GF*mtau**2*dsqrt(shat)/4d0/dsqrt(2d0)/pi
     .   *dsqrt(1d0-4d0*mtau**2/shat)**3

c  dsigma(gg->H->tautau)/dcosTheta
      distHJKWK=sigma0*prop*gamtau/2d0
C  1/2 because of FB symmetry


           return
        ENDIF
      
      ELSE               !   SPIN-2 EXCHANGE - coupling to CP-even GG
        IF(H1.EQ.H2) THEN
           distHJKWK=0.d0           ! for spin-2
           return 
         ELSE

c dsigma(gg->X->tautau)/dcosTheta
         distHJKWK=ce**2*(1-t**4)*shat**3*(httL**2+httR**2)/(
     .   8192*pi*f**4*((shat-mX**2)**2+mX**2*wX**2))



      return 
      ENDIF


      endif

      return
      END



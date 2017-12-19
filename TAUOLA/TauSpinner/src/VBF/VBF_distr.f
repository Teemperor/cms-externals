      REAL*8 FUNCTION VBFDISTR(ID1,ID2,ID3,ID4,HH1,HH2,PP,KEYIN)
C***************************************************************************
C* CALCULATES (matrix element)^2 for 
C        I1 I2 -> I3 I4 TAU+TAU- with given tau polarizations H1 H2
C        for a given set of momenta P
C        I1,...,I4 are PDG particle codes
C        if I3 or I4=0, then summed over final jet flavors 
C        H1 and H2 are tau+ and tau- helicities R: 1, L: -1, R+L: 0 
C        KeyIn=0: no Higgs
C        KeyIn=1: SM Higgs
C        KeyIn=2: no Higgs
C        KeyIn=3: SM Higgs
C        KeyIn=4,5: non-standard state 
C***************************************************************************
      IMPLICIT NONE

      INTEGER I1,I2,I3,I4,ID1,ID2,ID3,ID4,H1,H2,HH1,HH2,BUF_H,KEY,KEYIN
      REAL*8 P(0:3,6), PP(0:3,6),ANS
      INTEGER ICP
      COMMON /CPSTATUS/ ICP
      REAL*8  BUF(0:3)
      INTEGER BUF_I,IGLU,I,J,K
      
      LOGICAL INITIALIZED
      DATA    INITIALIZED/.FALSE./
      SAVE    INITIALIZED
      LOGICAL FLIPER,TESTUJEMY
      SAVE KEYSTORED
      INTEGER KEYSTORED
      TESTUJEMY=ID1.eq.-222.and.ID2.eq.1.and.ID3.eq.-1.and.ID4.eq.1
      ICP=0

      IF(.NOT.INITIALIZED) THEN
          INITIALIZED = .TRUE.
          KEYSTORED=KEYIN
      ENDIF

      IF (KEYIN.NE.KEYSTORED) THEN
       CALL VBF_REINIT(KEYIN) 
       KEYSTORED=KEYIN
      ENDIF

      KEY=MOD(KEYIN,2)
      IF(KEYIN.GT.3) THEN
          WRITE(*,*) 'non-standard state -- implementation not finished'
          STOP
      ELSE IF(KEY.NE.0.AND.KEY.NE.1) THEN
          WRITE(*,*) 'WRONG KEY'
          STOP
      ENDIF
      I1=ID1
      I2=ID2
      I3=ID3
      I4=ID4
      IF (TESTUJEMY) WRITE(*,*) 'idsy=',id1,id2,id3,id4
C     ------------
C     Fast track for cases where results must be zero. 
C     ------------

C     bayron number conservation:
C     --------------------------

      IF(I1+I2.EQ.42.AND.I3+I4.EQ.0) THEN
C       do nothing
      ELSEIF(I3+I4.EQ.42.AND.I1+I2.EQ.0) THEN
C       do nothing
      ELSEIF(I1*I2*I3*I4.LT.0) THEN
         VBFDISTR=0.0
         RETURN
      ENDIF

      IF(MOD(I1+I2+I3+I4,2).EQ.1) THEN
         VBFDISTR=0.0
         RETURN
      ENDIF

      IF(I1.LT.0.and.I2.LT.0.AND.(I3.GT.0.OR.I4.GT.0)) THEN
         VBFDISTR=0.0
         RETURN
      ENDIF

      IF(I3.LT.0.and.I4.LT.0.AND.(I1.GT.0.OR.I2.GT.0)) THEN
         VBFDISTR=0.0
         RETURN
      ENDIF

C charge conservation
C -------------------

      IF(SIGN(MOD(I1,2),I1)+SIGN(MOD(I2,2),I2).NE.SIGN(MOD(I3,2),I3)+SIGN(MOD(I4,2),I4)) THEN
        IF(I1+I2+I3+I4.LT.20) THEN 
         VBFDISTR=0.0
         RETURN
        ENDIF
      ENDIF

C zero or two gluons
C ------------------

      IGLU=0
      IF(I1.EQ.21) IGLU=IGLU+1 
      IF(I2.EQ.21) IGLU=IGLU+1 
      IF(I3.EQ.21) IGLU=IGLU+1 
      IF(I4.EQ.21) IGLU=IGLU+1


      IF(IGLU.EQ.1.OR.IGLU.GT.2) THEN
         VBFDISTR=0.0
         RETURN
      ENDIF

C if gluons are present extra consitency check.
C important, that this is before 5 to 1 replacement for configs with gluons.
C --------------------------------------------------------------------------
      IF(IGLU.EQ.2.AND.I1+I2.NE.I3+I4) THEN

         IF(.NOT.(I1+I2.EQ.0.OR.I3+I4.EQ.0)) THEN
            VBFDISTR=0.0
            RETURN
         ENDIF
      ENDIF

C =============================================
C now we have a chance that result is non zero.
C we copy configuration to local variables
C =============================================
      H1=HH1
      H2=HH2
      P(0:3,1) = PP(0:3,1)
      P(0:3,2) = PP(0:3,2)
      P(0:3,3) = PP(0:3,3)
      P(0:3,4) = PP(0:3,4)
      P(0:3,5) = PP(0:3,5)
      P(0:3,6) = PP(0:3,6)


C ===========================================
C REARRANING order of 4-vectors and ID's
C Parity reflection
C ===========================================


C treat  quarks 5 as 1 if gluons present. It is possible
C --------------------------------------
      ! IF (IGLU.EQ.2) THEN
      !  IF (ABS(I1).EQ.5) I1=I1/ABS(I1)
      !  IF (ABS(I2).EQ.5) I2=I2/ABS(I2)
      !  IF (ABS(I3).EQ.5) I3=I3/ABS(I3)
      !  IF (ABS(I4).EQ.5) I4=I4/ABS(I4)
      !ELSEIF((I1*I1.EQ.25.OR.I2*I2.EQ.25.OR.I3*I3.EQ.25.OR.I4*I4.EQ.25)) THEN
      IF((I1*I1.EQ.25.OR.I2*I2.EQ.25.OR.I3*I3.EQ.25.OR.I4*I4.EQ.25)) THEN
C in other cases processes with b-quarks are not yet installed.
C -------------------------------------------------------------
         VBFDISTR=0.0
         RETURN
      ENDIF

 
      if(testujemy) write(*,*) 'doszlimy do stepX',i1,i2, i3, i4

C If  I1 I2 particle/antiparticle (no gluon), |I1|>|I2| enforce
C -------------------------------------------------------------
      IF(I1*I2.LT.0.AND.I1+I2.LT.11.AND.I1**2.LT.I2**2) THEN
          ! flip IDs
          BUF_I = I1
          I1    = I2
          I2    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,1)
          P(0:3,1) = P(0:3,2)
          P(0:3,2) = BUF(0:3)
      ENDIF
       
 
      if(testujemy) write(*,*) 'doszlimy do step2',i1,i2, i3, i4


C C-reflection if all quarks antipartices
C ---------------------------------------    
      IF(
     $    (I1.LT.0.OR.I1.EQ.21).AND.
     $    (I2.LT.0.OR.I2.EQ.21).AND.
     $    (I3.LT.0.OR.I3.EQ.21).AND.
     $    (I4.LT.0.OR.I4.EQ.21)     ) THEN
         IF(I1.NE.21) I1=-I1
         IF(I2.NE.21) I2=-I2
         IF(I3.NE.21) I3=-I3
         IF(I4.NE.21) I4=-I4

!        C-parity  on taus
         BUF_H=H2
         H2=-H1
         H1=-BUF_H

         BUF(0:3) = P(0:3,5)
         P(0:3,5) = P(0:3,6)
         P(0:3,6) = BUF(0:3)
C  and P-Parity
         DO K=1,3
          DO J=1,6
            P(K,J)=-P(K,J)
          enddo
         enddo
         ICP=ICP+1
      ENDIF

C for incoming particle antiparticle larger |id| must be first,
C EXCEPTION: ID1 ID2 = (1,-2) or (-2,1)  There is UDX.f but no DUX.f file
C -----------------------------------------------------------------------


      FLIPER=(I1*I2.LT.0.AND.I1+I2.LT.11)
      IF(FLIPER) THEN
        FLIPER=I2*I2.GT.I1*I1
        IF(ID1*ID2.EQ.-2.AND.(ID1.EQ.1.OR.ID1.EQ.-2)) FLIPER=.NOT.FLIPER
      ENDIF
      IF(FLIPER) THEN
         IF(I1.NE.21) I1=-I1
         IF(I2.NE.21) I2=-I2
         IF(I3.NE.21) I3=-I3
         IF(I4.NE.21) I4=-I4

!        C-parity  on taus
         BUF_H=H2
         H2=-H1
         H1=-BUF_H
         BUF(0:3) = P(0:3,5)
         P(0:3,5) = P(0:3,6)
         P(0:3,6) = BUF(0:3)

C  and P-Parity
         DO K=1,3
          DO J=1,6
            P(K,J)=-P(K,J)
          enddo
         enddo

         ICP=ICP+1
      ENDIF
     
      if(testujemy) write(*,*) 'doszlimy do step3',i1,i2,i3,i4

C First incoming can not be antiparticle, second can not be lone gluon 
C ---------------------------------------------------------------------
      IF(I1.LT.0.OR.(I2.EQ.21.AND.I1.NE.21)) THEN  

          ! flip IDs
          BUF_I = I1
          I1    = I2
          I2    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,1)
          P(0:3,1) = P(0:3,2)
          P(0:3,2) = BUF(0:3)
      ENDIF

C If both I1 I2 positive (no gluon) enforce that I1>I2
      IF(I1.GT.0.AND.I2.GT.0.AND.I1+I2.LT.11.AND.I1.LT.I2) THEN
          ! flip IDs
          BUF_I = I1
          I1    = I2
          I2    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,1)
          P(0:3,1) = P(0:3,2)
          P(0:3,2) = BUF(0:3)
      ENDIF

C ================
C NOW FINAL STATES
C ================

C I3 never negative and I4 never alone 21
C --------------------------------------- 
      IF(I3.LT.0.OR.(I4.EQ.21.AND.I3.NE.21)) THEN
          ! flip IDs
          BUF_I = I3
          I3    = I4
          I4    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,3)
          P(0:3,3) = P(0:3,4)
          P(0:3,4) = BUF(0:3)
      ENDIF

C sole posibility <I3 even I4 odd> if both positive and non gluon
C ---------------------------------------------------------------
      IF(MOD(I3,2).EQ.1.AND.MOD(I4,2).EQ.0.AND.I3*I4.GT.0.AND.I3.NE.21) THEN
          ! flip IDs
          BUF_I = I3
          I3    = I4
          I4    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,3)
          P(0:3,3) = P(0:3,4)
          P(0:3,4) = BUF(0:3)
      ENDIF


C if I3,I4 simultaneously odd I4 must be larger/equal  I3
C -------------------------------------------------------
      IF(MOD(I3,2).EQ.1.AND.MOD(I4,2).EQ.1.AND.I3*I4.GT.0.AND.I3.GT.I4.AND.I3.NE.21) THEN
          ! flip IDs
          BUF_I = I3
          I3    = I4
          I4    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,3)
          P(0:3,3) = P(0:3,4)
          P(0:3,4) = BUF(0:3)
      ENDIF

C 
C if I3,I4 simultaneously even I4 must be larger/equal I3
C -------------------------------------------------------
      IF(MOD(I3,2).EQ.0.AND.MOD(I4,2).EQ.0.AND.I3*I4.GT.0.AND.I3.GT.I4) THEN
          ! flip IDs
          BUF_I = I3
          I3    = I4
          I4    = BUF_I

          ! flip 4-vectors
          BUF(0:3) = P(0:3,3)
          P(0:3,3) = P(0:3,4)
          P(0:3,4) = BUF(0:3)
      ENDIF


      if(testujemy) write(*,*) 'doszlimy do case-a ',i1,i2,i3,i4

      !
      ! FINALLY select appropriate function

C the ANS=0.0 is not needed unless there is something wrong in list for CASE below
      ANS=0.0 

      SELECT CASE(I1+I2)
         CASE(0)     !   UUX  DDX CCX SSX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.1) CALL DDX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.2) CALL UUX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SSX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.4) CALL CCX(P,I3,I4,H1,H2,KEY,ANS)
      if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4,ans
         CASE(1)     !   UDX INITIAL STATE
           IF(ABS(I1).EQ.2) CALL UDX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.4) CALL CSX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SUX(P,I3,I4,H1,H2,KEY,ANS)
         CASE(2)     !   DD INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do 2 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.1) CALL DD(P,I3,I4,H1,H2,KEY,ANS)
            IF(ABS(I1).EQ.4) CALL CUX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SDX(P,I3,I4,H1,H2,KEY,ANS)
         CASE(3)     !   UD INITIAL STATE
           IF(ABS(I1).EQ.2) CALL UD(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.4) CALL CDX(P,I3,I4,H1,H2,KEY,ANS)
         CASE(4)     !   UU INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do 4 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.2) CALL UU(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.1) CALL DS(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SD(P,I3,I4,H1,H2,KEY,ANS)
        CASE(22)    !   GLUON D INITIAL STATE
           CALL GD(P,I3,I4,H1,H2,KEY,ANS)
         CASE(23)    !   GLUON U INITIAL STATE
           CALL GU(P,I3,I4,H1,H2,KEY,ANS)
         CASE(24)    !   GLUON S INITIAL STATE
           CALL GD(P,I3,I4,H1,H2,KEY,ANS)
         CASE(25)    !   GLUON C INITIAL STATE
           CALL GU(P,I3,I4,H1,H2,KEY,ANS)
         CASE(26)    !   GLUON B INITIAL STATE
           CALL GD(P,I3,I4,H1,H2,KEY,ANS)
         CASE(42)    !   GLUON GLUON INITIAL STATE
           CALL GG(P,I3,I4,H1,H2,KEY,ANS)
         CASE(8)      ! CC INITIAL STATE
           CALL CC(P,I3,I4,H1,H2,KEY,ANS)
         CASE(7)      ! CS INITIAL STATE
           CALL CS(P,I3,I4,H1,H2,KEY,ANS)
         CASE(5)      ! DC INITIAL STATE
           IF(ABS(I1).EQ.1) CALL DC(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.4) CALL CD(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SU(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.2) CALL US(P,I3,I4,H1,H2,KEY,ANS)
         CASE(6)      ! SS INITIAL STATE
           if(testujemy) write(*,*) 'doszlismy do cu i1,i2=',i1,i2,i3,i4
           IF(ABS(I1).EQ.3) CALL SS(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.4) CALL CU(P,I3,I4,H1,H2,KEY,ANS)
         CASE(-2)      ! UCX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -2  case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.2) CALL UCX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.1) CALL DSX(P,I3,I4,H1,H2,KEY,ANS)
         CASE(-1)        ! USX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -1  case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.2) CALL USX(P,I3,I4,H1,H2,KEY,ANS)
           IF(ABS(I1).EQ.3) CALL SCX(P,I3,I4,H1,H2,KEY,ANS)
         CASE(-3)        ! DCX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -3  case-a ',i1,i2,i3,i4
           CALL DCX(P,I3,I4,H1,H2,KEY,ANS)
         CASE DEFAULT
 !          WRITE(*,*) "VBFDISTR: UNSUPPORTED PROCESS:",I1,I2
           ANS=0.0
      END SELECT
      IF(I3.NE.I4) ANS=ANS/2.D0   ! we divide by 2 because we do not order I3,I4
                                  ! but we will take both I3,I4 and I4,I3
      VBFDISTR=ANS

      END FUNCTION VBFDISTR

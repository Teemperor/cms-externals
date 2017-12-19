      REAL*8 FUNCTION SPIN2DISTR(ID1,ID2,ID3,ID4,HH1,HH2,PP,KEYIN)
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
C  
C HOWEVER: in this prototype calculation
C          parameter  KeyIn even though provided, is not used at all.
C***************************************************************************
      IMPLICIT NONE

      INTEGER I1,I2,I3,I4,ID1,ID2,ID3,ID4,H1,H2,BUF_H,KEY,KEYIN,HEL,HH1,HH2
      REAL*8 P(0:3,6),ANS,PP(0:3,6)
      INTEGER ICP
      COMMON /CPSTATUS2/ ICP
    
      REAL*8  BUF(0:3)
      INTEGER BUF_I,IGLU,I,J,K
      
      LOGICAL INITIALIZED
      DATA    INITIALIZED/.FALSE./
      SAVE    INITIALIZED
      LOGICAL FLIPER,TESTUJEMY

      TESTUJEMY=ID1.eq.-222.and.ID2.eq.1.and.ID3.eq.-1.and.ID4.eq.1
      ICP=0
    
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
         SPIN2DISTR=0.0
         RETURN
      ENDIF

      IF(MOD(I1+I2+I3+I4,2).EQ.1) THEN
         SPIN2DISTR=0.0
         RETURN
      ENDIF

      IF(I1.LT.0.and.I2.LT.0.AND.(I3.GT.0.OR.I4.GT.0)) THEN
         SPIN2DISTR=0.0
         RETURN
      ENDIF

      IF(I3.LT.0.and.I4.LT.0.AND.(I1.GT.0.OR.I2.GT.0)) THEN
         SPIN2DISTR=0.0
         RETURN
      ENDIF

C charge conservation
C -------------------

      IF(SIGN(MOD(I1,2),I1)+SIGN(MOD(I2,2),I2).NE.SIGN(MOD(I3,2),I3)+SIGN(MOD(I4,2),I4)) THEN
        IF(I1+I2+I3+I4.LT.20) THEN 
         SPIN2DISTR=0.0
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


      IF(IGLU.EQ.1.OR.IGLU.EQ.3.OR.IGLU.GT.4) THEN
         SPIN2DISTR=0.0
         RETURN
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
         SPIN2DISTR=0.0
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



C the ANS=0.0 is not needed unless there is something wrong in list for CASE below
      ANS=0.0 
C SPIN2 MATRIX ELEMENT USED IN CASE OF GG PROCESS
      
     
      SELECT CASE(I1+I2)

         CASE(0)     !   UUX  DDX CCX SSX INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.1) CALL ddx_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.2) CALL uux_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.3) CALL ssx_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.4) CALL ccx_s2(P,I3,I4,H1,H2,ANS)

           if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4,ans

         CASE(1)     !   UDX INITIAL STATE
           IF(ABS(I1).EQ.2) CALL udx_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.4) CALL csx_s2(P,I3,I4,H1,H2,ANS)

           IF(ABS(I1).EQ.3) CALL sux_s2(P,I3,I4,H1,H2,ANS)
           if(testujemy) write(*,*) 'doszlimy do 1 case-a ',i1,i2,i3,i4

         CASE(2)     !   DD INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do 2 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.1) CALL dd_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.4) CALL cux_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.3) CALL sdx_s2(P,I3,I4,H1,H2,ANS)

         CASE(3)     !   UD INITIAL STATE
           IF(ABS(I1).EQ.2) CALL ud_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.4) CALL cdx_s2(P,I3,I4,H1,H2,ANS)

         CASE(4)     !   UU INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do 4 case-a ',i1,i2,i3,i4
           IF(ABS(I1).EQ.2) CALL uu_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.3) CALL sd_s2(P,I3,I4,H1,H2,ANS)
         CASE(23)    !   GLUON U INITIAL STATE
           CALL gu_s2(P,I3,I4,H1,H2,ANS)

         CASE(22)
           CALL gu_s2(P,I3,I4,H1,H2,ANS)

         CASE(25)    !   GLUON C INITIAL STATE
           CALL gc_s2(P,I3,I4,H1,H2,ANS)

         CASE(24)
           CALL gc_s2(P,I3,I4,H1,H2,ANS) 

         CASE(42)    !   GLUON GLUON INITIAL STATE
           CALL gg_s2(P,I3,I4,H1,H2,ANS)

         CASE(19)
           CALL gux_s2(P,I3,I4,H1,H2,ANS)

         CASE(20)
           CALL gux_s2(P,I3,I4,H1,H2,ANS)

         CASE(17)    !   GLUON cx INITIAL STATE
           CALL gcx_s2(P,I3,I4,H1,H2,ANS)

         CASE(18)
           CALL gcx_s2(P,I3,I4,H1,H2,ANS)

         CASE(8)     
           IF(ABS(I1).EQ.4)  CALL cc_s2(P,I3,I4,H1,H2,ANS)

         CASE(7)      ! CS INITIAL STATE
           IF(ABS(I1).EQ.4)  CALL cs_s2(P,I3,I4,H1,H2,ANS)
         CASE(5)      ! DC INITIAL STATE
           IF(ABS(I1).EQ.4) CALL cd_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.3) CALL su_s2(P,I3,I4,H1,H2,ANS)
         CASE(6)      ! SS INITIAL STATE
           if(testujemy) write(*,*) 'doszlismy do cu i1,i2=',i1,i2,i3,i4
           IF(ABS(I1).EQ.3) CALL ss_s2(P,I3,I4,H1,H2,ANS)
           IF(ABS(I1).EQ.4) CALL cu_s2(P,I3,I4,H1,H2,ANS)
         CASE(-2)      ! UCX INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do -2  case-a ',i1,i2,i3,i4
           IF(I1.EQ.2)  CALL ucx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.1)  CALL dsx_s2(P,I3,I4,H1,H2,ANS)

         CASE(-1)        ! USX INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do -1  case-a ',i1,i2,i3,i4
           IF(I1.EQ.2) CALL usx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.3) CALL scx_s2(P,I3,I4,H1,H2,ANS)
   
         CASE(-3)        ! DCX INITIAL STATE
           if(testujemy) write(*,*) 'doszlimy do -3  case-a ',i1,i2,i3,i4
           IF(I1.EQ.1 .AND. I2.EQ.-4) CALL dcx_s2(P,I3,I4,H1,H2,ANS)
     

         CASE DEFAULT
         
           ANS=0.0
      END SELECT
      IF(I3.NE.I4) ANS=ANS/2.D0   ! we divide by 2 because we do not order I3,I4
                                  ! but we will take both I3,I4 and I4,I3
      SPIN2DISTR=ANS

      END FUNCTION SPIN2DISTR

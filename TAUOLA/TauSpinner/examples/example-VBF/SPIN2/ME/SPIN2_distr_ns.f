      REAL*8 FUNCTION SPIN2DISTR(ID1,ID2,ID3,ID4,H1,H2,P,KEYIN)
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

      INTEGER I1,I2,I3,I4,ID1,ID2,ID3,ID4,H1,H2,BUF_H,KEY,KEYIN,HEL
      REAL*8 P(0:3,6),ANS
      INTEGER ICP
      COMMON /CPSTATUS/ ICP
    
      REAL*8  BUF(0:3)
      INTEGER BUF_I,IGLU,I,J,K
      
      LOGICAL INITIALIZED
      DATA    INITIALIZED/.FALSE./
      SAVE    INITIALIZED
      LOGICAL FLIPER,TESTUJEMY

      TESTUJEMY=ID1.eq.-222.and.ID2.eq.1.and.ID3.eq.-1.and.ID4.eq.1


      
      I1=ID1
      I2=ID2
      I3=ID3
      I4=ID4
      IF (TESTUJEMY) WRITE(*,*) 'idsy=',id1,id2,id3,id4
C     ------------
C     Fast track for cases where results must be zero. 
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


      IF(IGLU.EQ.1.OR.IGLU.GT.2) THEN
         SPIN2DISTR=0.0
         RETURN
      ENDIF

C if gluons are present extra consitency check.
C important, that this is before 5 to 1 replacement for configs with gluons.
C --------------------------------------------------------------------------
      IF(IGLU.EQ.2.AND.I1+I2.NE.I3+I4) THEN

         IF(.NOT.(I1+I2.EQ.0.OR.I3+I4.EQ.0)) THEN
            SPIN2DISTR=0.0
            RETURN
         ENDIF
      ENDIF



C ===========================================
C REARRANING order of 4-vectors and ID's
C Parity reflection
C ===========================================

      IF((I1*I1.EQ.25.OR.I2*I2.EQ.25.OR.I3*I3.EQ.25.OR.I4*I4.EQ.25)) THEN
C in other cases processes with b-quarks are not yet installed.
C -------------------------------------------------------------
         SPIN2DISTR=0.0
         RETURN
       ENDIF


 
      if(testujemy) write(*,*) 'doszlimy do stepX',i1,i2, i3, i4


C the ANS=0.0 is not needed unless there is something wrong in list for CASE below
      ANS=0.0 
C SPIN2 MATRIX ELEMENT USED IN CASE OF GG PROCESS
      
      
      SELECT CASE(I1+I2)
         CASE(0)     !   UUX  DDX CCX SSX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4
           IF(I1.EQ.1) CALL ddx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.2) CALL uux_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.3) CALL ssx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.4) CALL ccx_s2(P,I3,I4,H1,H2,ANS)

      if(testujemy) write(*,*) 'doszlimy do 0 case-a ',i1,i2,i3,i4,ans
         CASE(1)     !   UDX INITIAL STATE
           IF(I1.EQ.2 .AND. I2.EQ.-1) CALL udx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.4 .AND. I2.EQ.-3) CALL csx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.3 .AND. I2.EQ.-2) CALL sux_s2(P,I3,I4,H1,H2,ANS)
         CASE(2)     !   DD INITIAL STATE

      if(testujemy) write(*,*) 'doszlimy do 2 case-a ',i1,i2,i3,i4
           IF(I1.EQ.1 .AND. I2.EQ.1) CALL dd_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.4 .AND. I2.EQ.-2) CALL cux_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.3 .AND. I2.eq.-1) CALL sdx_s2(P,I3,I4,H1,H2,ANS)
         CASE(3)     !   UD INITIAL STATE
           IF(I1.EQ.2 .AND. I2.EQ.1) CALL ud_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.4 .AND. I2.EQ.-1) CALL cdx_s2(P,I3,I4,H1,H2,ANS)
         CASE(4)     !   UU INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do 4 case-a ',i1,i2,i3,i4
           IF(I1.EQ.2 .AND. I2.EQ.2) CALL uu_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.1 .AND. I2.EQ.3) CALL ds_s2(P,I3,I4,H1,H2,ANS)
c           IF(ABS(I1).EQ.3) CALL SD(P,I3,I4,H1,H2,KEY,ANS)
        CASE(22)    !   GLUON D INITIAL STATE
           CALL gd_s2(P,I3,I4,H1,H2,ANS)
         CASE(23)    !   GLUON U INITIAL STATE
           CALL gu_s2(P,I3,I4,H1,H2,ANS)
         CASE(25)    !   GLUON C INITIAL STATE
           CALL gc_s2(P,I3,I4,H1,H2,ANS)
         CASE(42)    !   GLUON GLUON INITIAL STATE
           CALL gg_s2(P,I3,I4,H1,H2,ANS)
         CASE(20)
           CALL gdx_s2(P,I3,I4,H1,H2,ANS)
         CASE(19)
           CALL gux_s2(P,I3,I4,H1,H2,ANS)
         CASE(17)    !   GLUON S INITIAL STATE
           CALL gcx_s2(P,I3,I4,H1,H2,ANS)
         CASE(8)     
         IF(I1.EQ.4.AND.I2.EQ.4)  CALL cc_s2(P,I3,I4,H1,H2,ANS)
         CASE(7)      ! CS INITIAL STATE
         IF(I1.EQ.4.AND.I2.EQ.3)  CALL cs_s2(P,I3,I4,H1,H2,ANS)
         CASE(5)      ! DC INITIAL STATE

           IF(I1.EQ.4 .AND. I2.EQ.1) CALL cd_s2(P,I3,I4,H1,H2,ANS)

           IF(I1.EQ.2 .AND. I2.EQ.3) CALL us_s2(P,I3,I4,H1,H2,ANS)
         CASE(6)      ! SS INITIAL STATE
           if(testujemy) write(*,*) 'doszlismy do cu i1,i2=',i1,i2,i3,i4
           IF(I1.EQ.3 .AND. I3.EQ.3) CALL ss_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.2 .AND. I2.EQ.4) CALL uc_s2(P,I3,I4,H1,H2,ANS)
         
         CASE(-2)      ! UCX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -2  case-a ',i1,i2,i3,i4
           IF(I1.EQ.2 .AND. I2.EQ.-4) CALL ucx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.-1 .AND. I2.EQ.-1) CALL dxdx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.1 .AND. I2.EQ.-3) CALL dsx_s2(P,I3,I4,H1,H2,ANS)
         CASE(-1)        ! USX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -1  case-a ',i1,i2,i3,i4
           IF(I1.EQ.2 .AND. I2.EQ.-3) CALL usx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.3 .AND. I2.EQ.-4) CALL scx_s2(P,I3,I4,H1,H2,ANS)
           IF(I1.EQ.1 .AND. I2.EQ.-2) CALL dux_s2(P,I3,I4,H1,H2,ANS)
         CASE(-3)        ! DCX INITIAL STATE
      if(testujemy) write(*,*) 'doszlimy do -3  case-a ',i1,i2,i3,i4
       IF(I1.EQ.1 .AND. I2.EQ.-4) CALL dcx_s2(P,I3,I4,H1,H2,ANS)
       IF(I1.EQ.-2 .AND. I2.EQ.-1) CALL uxdx_s2(P,I3,I4,H1,H2,ANS)
        CASE(-4)
       IF(I1.EQ.-2 .AND. I2.EQ.-2) CALL uxux_s2(P,I3,I4,H1,H2,ANS)
       IF(I1.EQ.-1 .AND. I2.EQ.-3) CALL dxsx_s2(P,I3,I4,H1,H2,ANS)
        CASE(-5)
       IF(I1.EQ.-4 .AND. I2.EQ.-1) CALL cxdx_s2(P,I3,I4,H1,H2,ANS)
       IF(I1.EQ.-2 .AND. I2.EQ.-3) CALL uxsx_s2(P,I3,I4,H1,H2,ANS)
        CASE(-6)
       IF(I1.EQ.-2 .AND. I2.EQ.-4) CALL uxcx_s2(P,I3,I4,H1,H2,ANS)
       IF(I1.EQ.-3 .AND. I2.EQ.-3) CALL sxsx_s2(P,I3,I4,H1,H2,ANS)
        CASE(-7)
       IF(I1.EQ.-4 .AND. I2.EQ.-3) CALL cxsx_s2(P,I3,I4,H1,H2,ANS)
        CASE(-8)
       IF(I1.EQ.-4 .AND. I2.EQ.-4) CALL cxcx_s2(P,I3,I4,H1,H2,ANS)       
         CASE DEFAULT
         
           ANS=0.0
      END SELECT
c      IF(I3.NE.I4) ANS=ANS/2.D0   ! we divide by 2 because we do not order I3,I4
                                  ! but we will take both I3,I4 and I4,I3
      SPIN2DISTR=ANS

      END FUNCTION SPIN2DISTR

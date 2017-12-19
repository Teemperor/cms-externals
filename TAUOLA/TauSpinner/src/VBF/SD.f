C  JK. done 27.10.2015
c------------------begin subprocess initiated by SD----------
      SUBROUTINE  SD(P,I3,I4,H1,H2,KEY,ANS)
      IMPLICIT NONE
 
      INTEGER I3,I4,   H1,H2,   KEY
      REAL*8 P(0:3,6), ANS

C      I3=0,1             I4=0,3    ! only DS->DS
      
      ANS=0.D0

      IF(KEY.EQ.1) THEN
         IF(I3.EQ.1 .AND. I4.EQ.3) CALL SD_DS_H(P,H1,H2,ANS)
         IF(I3.EQ.0 .AND. I4.EQ.0) CALL SD_DS_H(P,H1,H2,ANS)
      ELSE IF(KEY.EQ.0) THEN 
         IF(I3.EQ.1 .AND. I4.EQ.3) CALL SD_DS_NOH(P,H1,H2,ANS)
         IF(I3.EQ.0 .AND. I4.EQ.0) CALL SD_DS_NOH(P,H1,H2,ANS)
      ELSE 
         WRITE(*,*) 'spin=2  NOT FINISHED'
         STOP
      ENDIF
      END ! SUBROUTINE SD

C ---------begin subprocess SD->jjH with H-> tautau
c   ---------------------jj=DS only -------------

      SUBROUTINE SD_DS_H(P,H1,H2,ANS)
C     
C     Generated by MadGraph 5 v. 1.5.15, 2013-12-11
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: s d > d s h WEIGHTED=6
C     *   Decay: h > ta+ ta- WEIGHTED=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
	  INTEGER H1,H2
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX_SD_DS_H
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
          T=MATRIX_SD_DS_H(P ,H1,H2,NHEL(1,IHEL),JC(1))
          ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_SD_DS_H(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph 5 v. 1.5.15, 2013-12-11
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: s d > d s h WEIGHTED=6
C     *   Decay: h > ta+ ta- WEIGHTED=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=7, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    9/
C     1 T(3,2) T(4,1)
	  
      INTEGER H1,H2
      REAL*8 MATRIX
      MATRIX_SD_DS_H=0.D0
      IF(H1.EQ.0. OR .H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0. OR .H2.EQ.NHEL(6)) THEN
	  
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFS4_3(W(1,5),W(1,6),GC_99,MH,WH,W(1,7))
      CALL FFV2_3_3(W(1,1),W(1,4),GC_50,GC_58,MZ,WZ,W(1,6))
      CALL FFV2_3_3(W(1,2),W(1,3),GC_50,GC_58,MZ,WZ,W(1,4))
C     Amplitude(s) for diagram number 1
      CALL VVS1_0(W(1,6),W(1,4),W(1,7),GC_81,AMP(1))
      JAMP(1)=+AMP(1)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_SD_DS_H=MATRIX
      ENDIF             ! CLOSES H1 IF
      ENDIF             ! CLOSES H2 IF
      END
C-----------------------end subprocess SD->ddH-------------




c-----------------begin subprocess SD->jj_noH--------------
c                           jj=DS only

      SUBROUTINE SD_DS_NOH(P,H1,H2,ANS)
C     
C     Generated by MadGraph 5 v. 1.5.15, 2013-12-11
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: s d > d s ta+ ta- / h QED=4
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
	INTEGER H1,H2
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX_SD_DS_NOH
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
          T=MATRIX_SD_DS_NOH(P ,H1,H2,NHEL(1,IHEL),JC(1))
          ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_SD_DS_NOH(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph 5 v. 1.5.15, 2013-12-11
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: s d > d s ta+ ta- / h QED=4
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=32)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=13, NCOLOR=2)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C       
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  2) /    9,    3/
C     1 T(3,1) T(4,2)
      DATA DENOM(2)/1/
      DATA (CF(I,  2),I=  1,  2) /    3,    9/
C     1 T(3,2) T(4,1)
	
      INTEGER H1,H2
      REAL*8 MATRIX
      MATRIX_SD_DS_NOH=0.D0
      IF(H1.EQ.0. OR .H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0. OR .H2.EQ.NHEL(6)) THEN
	
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFV1P0_3(W(1,1),W(1,4),GC_1,ZERO,ZERO,W(1,7))
      CALL FFV1P0_3(W(1,2),W(1,3),GC_1,ZERO,ZERO,W(1,8))
      CALL FFV1_2(W(1,5),W(1,7),GC_3,MTA,ZERO,W(1,9))
C     Amplitude(s) for diagram number 1
      CALL FFV1_0(W(1,9),W(1,6),W(1,8),GC_3,AMP(1))
      CALL FFV1_1(W(1,6),W(1,7),GC_3,MTA,ZERO,W(1,10))
C     Amplitude(s) for diagram number 2
      CALL FFV1_0(W(1,5),W(1,10),W(1,8),GC_3,AMP(2))
      CALL FFV2_3_3(W(1,2),W(1,3),GC_50,GC_58,MZ,WZ,W(1,11))
C     Amplitude(s) for diagram number 3
      CALL FFV2_4_0(W(1,9),W(1,6),W(1,11),GC_50,GC_59,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL FFV2_4_0(W(1,5),W(1,10),W(1,11),GC_50,GC_59,AMP(4))
      CALL FFV2_3_3(W(1,1),W(1,4),GC_50,GC_58,MZ,WZ,W(1,10))
      CALL FFV2_4_2(W(1,5),W(1,10),GC_50,GC_59,MTA,ZERO,W(1,9))
C     Amplitude(s) for diagram number 5
      CALL FFV1_0(W(1,9),W(1,6),W(1,8),GC_3,AMP(5))
      CALL FFV2_4_1(W(1,6),W(1,10),GC_50,GC_59,MTA,ZERO,W(1,12))
C     Amplitude(s) for diagram number 6
      CALL FFV1_0(W(1,5),W(1,12),W(1,8),GC_3,AMP(6))
C     Amplitude(s) for diagram number 7
      CALL FFV2_4_0(W(1,9),W(1,6),W(1,11),GC_50,GC_59,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL FFV2_4_0(W(1,5),W(1,12),W(1,11),GC_50,GC_59,AMP(8))
      CALL FFV1P0_3(W(1,5),W(1,6),GC_3,ZERO,ZERO,W(1,12))
      CALL FFV1_2(W(1,2),W(1,7),GC_1,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 9
      CALL FFV1_0(W(1,9),W(1,3),W(1,12),GC_1,AMP(9))
      CALL FFV1_1(W(1,3),W(1,7),GC_1,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 10
      CALL FFV1_0(W(1,2),W(1,13),W(1,12),GC_1,AMP(10))
      CALL FFV2_4_3(W(1,5),W(1,6),GC_50,GC_59,MZ,WZ,W(1,7))
C     Amplitude(s) for diagram number 11
      CALL FFV2_3_0(W(1,9),W(1,3),W(1,7),GC_50,GC_58,AMP(11))
C     Amplitude(s) for diagram number 12
      CALL FFV2_3_0(W(1,2),W(1,13),W(1,7),GC_50,GC_58,AMP(12))
      CALL FFV1P0_3(W(1,1),W(1,4),GC_11,ZERO,ZERO,W(1,13))
      CALL FFV1_2(W(1,2),W(1,13),GC_11,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 13
      CALL FFV1_0(W(1,9),W(1,3),W(1,12),GC_1,AMP(13))
      CALL FFV1_1(W(1,3),W(1,13),GC_11,ZERO,ZERO,W(1,6))
C     Amplitude(s) for diagram number 14
      CALL FFV1_0(W(1,2),W(1,6),W(1,12),GC_1,AMP(14))
C     Amplitude(s) for diagram number 15
      CALL FFV2_3_0(W(1,9),W(1,3),W(1,7),GC_50,GC_58,AMP(15))
C     Amplitude(s) for diagram number 16
      CALL FFV2_3_0(W(1,2),W(1,6),W(1,7),GC_50,GC_58,AMP(16))
      CALL FFV2_3_2(W(1,2),W(1,10),GC_50,GC_58,ZERO,ZERO,W(1,6))
C     Amplitude(s) for diagram number 17
      CALL FFV1_0(W(1,6),W(1,3),W(1,12),GC_1,AMP(17))
      CALL FFV2_3_1(W(1,3),W(1,10),GC_50,GC_58,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 18
      CALL FFV1_0(W(1,2),W(1,9),W(1,12),GC_1,AMP(18))
C     Amplitude(s) for diagram number 19
      CALL FFV2_3_0(W(1,6),W(1,3),W(1,7),GC_50,GC_58,AMP(19))
C     Amplitude(s) for diagram number 20
      CALL FFV2_3_0(W(1,2),W(1,9),W(1,7),GC_50,GC_58,AMP(20))
      CALL FFV1_2(W(1,1),W(1,8),GC_1,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 21
      CALL FFV1_0(W(1,9),W(1,4),W(1,12),GC_1,AMP(21))
      CALL FFV1_2(W(1,1),W(1,12),GC_1,ZERO,ZERO,W(1,6))
C     Amplitude(s) for diagram number 22
      CALL FFV1_0(W(1,6),W(1,4),W(1,8),GC_1,AMP(22))
C     Amplitude(s) for diagram number 23
      CALL FFV2_3_0(W(1,9),W(1,4),W(1,7),GC_50,GC_58,AMP(23))
      CALL FFV2_3_2(W(1,1),W(1,7),GC_50,GC_58,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 24
      CALL FFV1_0(W(1,9),W(1,4),W(1,8),GC_1,AMP(24))
      CALL FFV1P0_3(W(1,2),W(1,3),GC_11,ZERO,ZERO,W(1,8))
      CALL FFV1_2(W(1,1),W(1,8),GC_11,ZERO,ZERO,W(1,3))
C     Amplitude(s) for diagram number 25
      CALL FFV1_0(W(1,3),W(1,4),W(1,12),GC_1,AMP(25))
C     Amplitude(s) for diagram number 26
      CALL FFV1_0(W(1,6),W(1,4),W(1,8),GC_11,AMP(26))
C     Amplitude(s) for diagram number 27
      CALL FFV2_3_0(W(1,3),W(1,4),W(1,7),GC_50,GC_58,AMP(27))
C     Amplitude(s) for diagram number 28
      CALL FFV1_0(W(1,9),W(1,4),W(1,8),GC_11,AMP(28))
      CALL FFV2_3_2(W(1,1),W(1,11),GC_50,GC_58,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 29
      CALL FFV1_0(W(1,8),W(1,4),W(1,12),GC_1,AMP(29))
C     Amplitude(s) for diagram number 30
      CALL FFV2_3_0(W(1,6),W(1,4),W(1,11),GC_50,GC_58,AMP(30))
C     Amplitude(s) for diagram number 31
      CALL FFV2_3_0(W(1,8),W(1,4),W(1,7),GC_50,GC_58,AMP(31))
C     Amplitude(s) for diagram number 32
      CALL FFV2_3_0(W(1,9),W(1,4),W(1,11),GC_50,GC_58,AMP(32))
      JAMP(1)=+1./2.*(+AMP(13)+AMP(14)+AMP(15)+AMP(16)+AMP(25)+AMP(26)
     $ +AMP(27)+AMP(28))
      JAMP(2)=+AMP(1)+AMP(2)+AMP(3)+AMP(4)+AMP(5)+AMP(6)+AMP(7)+AMP(8)
     $ +AMP(9)+AMP(10)+AMP(11)+AMP(12)-1./6.*AMP(13)-1./6.*AMP(14)
     $ -1./6.*AMP(15)-1./6.*AMP(16)+AMP(17)+AMP(18)+AMP(19)+AMP(20)
     $ +AMP(21)+AMP(22)+AMP(23)+AMP(24)-1./6.*AMP(25)-1./6.*AMP(26)
     $ -1./6.*AMP(27)-1./6.*AMP(28)+AMP(29)+AMP(30)+AMP(31)+AMP(32)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_SD_DS_NOH=MATRIX
      ENDIF
      ENDIF
      END
C---------------------end subprocess initiated by SD-------
C DONE 27.10.2015 

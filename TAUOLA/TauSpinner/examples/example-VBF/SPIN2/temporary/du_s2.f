      SUBROUTINE du_s2(P,I3,I4,H1,H2,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)

C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER I3,I4,H1,H2
C     
C     GLOBAL VARIABLES
C     
C     ----------
C     BEGIN CODE
C     ----------
      IF(I3.EQ.4 .AND. I4.EQ.1) CALL du_cd_s2(P,H1,H2,ANS)
      IF(I3.EQ.4 .AND. I4.EQ.3) CALL du_cs_s2(P,H1,H2,ANS)
      IF(I3.EQ.2 .AND. I4.EQ.1) CALL du_ud_s2(P,H1,H2,ANS)
      IF(I3.EQ.2 .AND. I4.EQ.3) CALL du_us_s2(P,H1,H2,ANS)
      END

      SUBROUTINE du_cd_s2(P,H1,H2,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d u > c d x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
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
      REAL*8 MATRIX_du_cd_s2
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     

      DATA (NHEL(I,   1),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
            T=MATRIX_du_cd_s2(P,H1,H2 ,NHEL(1,IHEL),JC(1))
            ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_du_cd_s2(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d u > c d x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=4)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=8, NCOLOR=1)
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
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl2.inc'
      INTEGER H1, H2
      REAL*8 MATRIX
      MATRIX_du_cd_s2=0.0D0

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    9/
C     1 T(3,1) T(4,2)
C     ----------
      IF(H1.EQ.0 .OR. H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0 .OR. H2.EQ.NHEL(6)) THEN
C     BEGIN CODE
C     ----------
      CALL IXXXXX_S(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX_S(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX_S(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX_S(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX_S(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX_S(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFT4_5_3(W(1,5),W(1,6),GC_14,GC_15,MX2,WX2,W(1,7))
      CALL FFV5_3_S(W(1,1),W(1,3),GC_29,MW,WW,W(1,6))
      CALL FFV5_3_S(W(1,2),W(1,4),GC_56,MW,WW,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVT2_0(W(1,5),W(1,6),W(1,7),GC_16,AMP(1))
      CALL FFT4_5_2(W(1,2),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 2
      CALL FFV5_0_S(W(1,8),W(1,4),W(1,6),GC_56,AMP(2))
      CALL FFT4_5_1(W(1,4),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 3
      CALL FFV5_0_S(W(1,2),W(1,8),W(1,6),GC_56,AMP(3))
      CALL FFT4_5_2(W(1,1),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 4
      CALL FFV5_0_S(W(1,8),W(1,3),W(1,5),GC_29,AMP(4))
      JAMP(1)=-AMP(1)-AMP(2)-AMP(3)-AMP(4)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_du_cd_s2=MATRIX
      ENDIF
      ENDIF
      END
     
C     ----------


   

      SUBROUTINE du_cs_s2(P,H1,H2,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d u > c s x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
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
      REAL*8 MATRIX_du_cs_s2
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     

      DATA (NHEL(I,   1),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
            T=MATRIX_du_cs_s2(P,H1,H2 ,NHEL(1,IHEL),JC(1))
            ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_du_cs_s2(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d u > c s x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=3)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=8, NCOLOR=1)
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
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl2.inc'
      INTEGER H1, H2
      REAL*8 MATRIX
      MATRIX_du_cs_s2=0.0D0

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    9/
C     1 T(3,1) T(4,2)
C     ----------
      IF(H1.EQ.0 .OR. H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0 .OR. H2.EQ.NHEL(6)) THEN
C     BEGIN CODE
C     ----------
      CALL IXXXXX_S(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX_S(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX_S(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX_S(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX_S(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX_S(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFT4_5_3(W(1,5),W(1,6),GC_14,GC_15,MX2,WX2,W(1,7))
      CALL FFV5_3_S(W(1,1),W(1,3),GC_29,MW,WW,W(1,6))
      CALL FFV5_3_S(W(1,2),W(1,4),GC_57,MW,WW,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVT2_0(W(1,5),W(1,6),W(1,7),GC_16,AMP(1))
      CALL FFT4_5_2(W(1,2),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 2
      CALL FFV5_0_S(W(1,8),W(1,4),W(1,6),GC_57,AMP(2))
      CALL FFT4_5_2(W(1,1),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 3
      CALL FFV5_0_S(W(1,8),W(1,3),W(1,5),GC_29,AMP(3))
      JAMP(1)=-AMP(1)-AMP(2)-AMP(3)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_du_cs_s2=MATRIX
      ENDIF
      ENDIF
      END
      


   

      SUBROUTINE du_ud_s2(P,H1,H2,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d u > u d x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
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
      REAL*8 MATRIX_du_ud_s2
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     

      DATA (NHEL(I,   1),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
            T=MATRIX_du_ud_s2(P,H1,H2 ,NHEL(1,IHEL),JC(1))
            ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_du_ud_s2(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d u > u d x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=26)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=15, NCOLOR=2)
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
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl2.inc'
      INTEGER H1, H2
      REAL*8 MATRIX
      MATRIX_du_ud_s2=0.0D0

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  2) /    9,    3/
C     1 T(3,1) T(4,2)
      DATA DENOM(2)/1/
      DATA (CF(I,  2),I=  1,  2) /    3,    9/
C     1 T(3,2) T(4,1)
C     ----------
      IF(H1.EQ.0 .OR. H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0 .OR. H2.EQ.NHEL(6)) THEN
C     BEGIN CODE
C     ----------
      CALL IXXXXX_S(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX_S(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX_S(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX_S(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX_S(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX_S(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFT4_5_3(W(1,5),W(1,6),GC_14,GC_15,MX2,WX2,W(1,7))
      CALL FFV5_3_S(W(1,1),W(1,3),GC_26,MW,WW,W(1,6))
      CALL FFV5_3_S(W(1,2),W(1,4),GC_56,MW,WW,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVT2_0(W(1,5),W(1,6),W(1,7),GC_16,AMP(1))
      CALL FFT4_5_2(W(1,2),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 2
      CALL FFV5_0_S(W(1,8),W(1,4),W(1,6),GC_56,AMP(2))
      CALL FFT4_5_1(W(1,4),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 3
      CALL FFV5_0_S(W(1,2),W(1,9),W(1,6),GC_56,AMP(3))
      CALL FFV4P0_3(W(1,1),W(1,4),GC_1,ZERO,ZERO,W(1,6))
      CALL FFV4P0_3(W(1,2),W(1,3),GC_2,ZERO,ZERO,W(1,10))
C     Amplitude(s) for diagram number 4
      CALL VVT2_0(W(1,6),W(1,10),W(1,7),GC_49,AMP(4))
      CALL FFV5_6_3(W(1,2),W(1,3),GC_43,GC_40,MZ,WZ,W(1,11))
C     Amplitude(s) for diagram number 5
      CALL VVT2_0(W(1,6),W(1,11),W(1,7),GC_46,AMP(5))
      CALL FFV4P0_3(W(1,1),W(1,4),GC_7,ZERO,ZERO,W(1,12))
      CALL FFV4P0_3(W(1,2),W(1,3),GC_7,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 6
      CALL VVT2_0(W(1,12),W(1,13),W(1,7),GC_9,AMP(6))
      CALL FFV5_6_3(W(1,1),W(1,4),GC_42,GC_39,MZ,WZ,W(1,14))
C     Amplitude(s) for diagram number 7
      CALL VVT2_0(W(1,10),W(1,14),W(1,7),GC_46,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL VVT2_0(W(1,14),W(1,11),W(1,7),GC_48,AMP(8))
      CALL FFT4_5_3(W(1,1),W(1,4),GC_12,GC_13,MX2,WX2,W(1,15))
C     Amplitude(s) for diagram number 9
      CALL FFT4_5_0(W(1,8),W(1,3),W(1,15),GC_12,GC_13,AMP(9))
C     Amplitude(s) for diagram number 10
      CALL FFV4_0_S(W(1,8),W(1,3),W(1,6),GC_2,AMP(10))
C     Amplitude(s) for diagram number 11
      CALL FFV4_0_S(W(1,8),W(1,3),W(1,12),GC_7,AMP(11))
C     Amplitude(s) for diagram number 12
      CALL FFV5_6_0(W(1,8),W(1,3),W(1,14),GC_43,GC_40,AMP(12))
      CALL FFT4_5_1(W(1,3),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 13
      CALL FFT4_5_0(W(1,2),W(1,8),W(1,15),GC_12,GC_13,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL FFV4_0_S(W(1,2),W(1,8),W(1,6),GC_2,AMP(14))
C     Amplitude(s) for diagram number 15
      CALL FFV4_0_S(W(1,2),W(1,8),W(1,12),GC_7,AMP(15))
C     Amplitude(s) for diagram number 16
      CALL FFV5_6_0(W(1,2),W(1,8),W(1,14),GC_43,GC_40,AMP(16))
      CALL FFT4_5_2(W(1,1),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,14))
      CALL FFT4_5_3(W(1,2),W(1,3),GC_12,GC_13,MX2,WX2,W(1,7))
C     Amplitude(s) for diagram number 17
      CALL FFT4_5_0(W(1,14),W(1,4),W(1,7),GC_12,GC_13,AMP(17))
C     Amplitude(s) for diagram number 18
      CALL FFV4_0_S(W(1,14),W(1,4),W(1,10),GC_1,AMP(18))
C     Amplitude(s) for diagram number 19
      CALL FFV4_0_S(W(1,14),W(1,4),W(1,13),GC_7,AMP(19))
C     Amplitude(s) for diagram number 20
      CALL FFV5_6_0(W(1,14),W(1,4),W(1,11),GC_42,GC_39,AMP(20))
C     Amplitude(s) for diagram number 21
      CALL FFV5_0_S(W(1,14),W(1,3),W(1,5),GC_26,AMP(21))
C     Amplitude(s) for diagram number 22
      CALL FFT4_5_0(W(1,1),W(1,9),W(1,7),GC_12,GC_13,AMP(22))
C     Amplitude(s) for diagram number 23
      CALL FFV4_0_S(W(1,1),W(1,9),W(1,10),GC_1,AMP(23))
C     Amplitude(s) for diagram number 24
      CALL FFV4_0_S(W(1,1),W(1,9),W(1,13),GC_7,AMP(24))
C     Amplitude(s) for diagram number 25
      CALL FFV5_6_0(W(1,1),W(1,9),W(1,11),GC_42,GC_39,AMP(25))
C     Amplitude(s) for diagram number 26
      CALL FFV5_0_S(W(1,1),W(1,8),W(1,5),GC_26,AMP(26))
      JAMP(1)=-AMP(1)-AMP(2)-AMP(3)+1D0/2D0*AMP(6)+1D0/2D0*AMP(11)+1D0
     $ /2D0*AMP(15)+1D0/2D0*AMP(19)-AMP(21)+1D0/2D0*AMP(24)-AMP(26)
      JAMP(2)=+AMP(4)+AMP(5)-1D0/6D0*AMP(6)+AMP(7)+AMP(8)+AMP(9)
     $ +AMP(10)-1D0/6D0*AMP(11)+AMP(12)+AMP(13)+AMP(14)-1D0/6D0*AMP(15)
     $ +AMP(16)+AMP(17)+AMP(18)-1D0/6D0*AMP(19)+AMP(20)+AMP(22)+AMP(23)
     $ -1D0/6D0*AMP(24)+AMP(25)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_du_ud_s2=MATRIX
      ENDIF
      ENDIF
      END
     



      SUBROUTINE du_us_s2(P,H1,H2,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: d u > u s x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
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
      REAL*8 MATRIX_du_us_s2
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     

      DATA (NHEL(I,   1),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA IDEN/36/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
            T=MATRIX_du_us_s2(P,H1,H2 ,NHEL(1,IHEL),JC(1))
            ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_du_us_s2(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: d u > u s x QED<=99 NPqq<=99 NPll<=99 QCD<=99 NPgg<=99
C      NPvv<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=4)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=8, NCOLOR=1)
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
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl2.inc'
      INTEGER H1, H2
      REAL*8 MATRIX
      MATRIX_du_us_s2=0.0D0

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    9/
C     1 T(3,1) T(4,2)
C     ----------
      IF(H1.EQ.0 .OR. H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0 .OR. H2.EQ.NHEL(6)) THEN
C     BEGIN CODE
C     ----------
      CALL IXXXXX_S(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL IXXXXX_S(P(0,2),ZERO,NHEL(2),+1*IC(2),W(1,2))
      CALL OXXXXX_S(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL OXXXXX_S(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX_S(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX_S(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFT4_5_3(W(1,5),W(1,6),GC_14,GC_15,MX2,WX2,W(1,7))
      CALL FFV5_3_S(W(1,1),W(1,3),GC_26,MW,WW,W(1,6))
      CALL FFV5_3_S(W(1,2),W(1,4),GC_57,MW,WW,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVT2_0(W(1,5),W(1,6),W(1,7),GC_16,AMP(1))
      CALL FFT4_5_2(W(1,2),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 2
      CALL FFV5_0_S(W(1,8),W(1,4),W(1,6),GC_57,AMP(2))
      CALL FFT4_5_2(W(1,1),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 3
      CALL FFV5_0_S(W(1,8),W(1,3),W(1,5),GC_26,AMP(3))
      CALL FFT4_5_1(W(1,3),W(1,7),GC_12,GC_13,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 4
      CALL FFV5_0_S(W(1,1),W(1,8),W(1,5),GC_26,AMP(4))
      JAMP(1)=-AMP(1)-AMP(2)-AMP(3)-AMP(4)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_du_us_s2=MATRIX
      ENDIF
      ENDIF
      END

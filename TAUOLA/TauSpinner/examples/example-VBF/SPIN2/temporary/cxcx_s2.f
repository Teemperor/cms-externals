      SUBROUTINE cxcx_s2(P,I3,I4,H1,H2,ANS)
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
      IF(I3.EQ.-4 .AND. I4.EQ.-4) CALL cxcx_cxcx_s2(P,H1,H2,ANS)

      END

      SUBROUTINE cxcx_cxcx_s2(P,H1,H2,ANS)
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
C     Process: c~ c~ > c~ c~ x QED<=99 NPqq<=99 NPll<=99 QCD<=99
C      NPgg<=99 NPVV<=99 @1
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
      REAL*8 MATRIX_cxcx_cxcx_s2
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     

      DATA (NHEL(I,   1),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA IDEN/72/
C     ----------
C     BEGIN CODE
C     ----------
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
            T=MATRIX_cxcx_cxcx_s2(P,H1,H2 ,NHEL(1,IHEL),JC(1))
            ANS=ANS+T
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_cxcx_cxcx_s2(P,H1,H2,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.4.3, 2016-08-01
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: c~ c~ > c~ c~ x QED<=99 NPqq<=99 NPll<=99 QCD<=99
C      NPgg<=99 NPVV<=99 @1
C     *   Decay: x > ta+ ta- WEIGHTED<=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=10)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=9, NCOLOR=2)
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
      MATRIX_cxcx_cxcx_s2=0.0D0

C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  2) /    9,    3/
C     1 T(1,3) T(2,4)
      DATA DENOM(2)/1/
      DATA (CF(I,  2),I=  1,  2) /    3,    9/
C     1 T(1,4) T(2,3)
C     ----------
      IF(H1.EQ.0 .OR. H1.EQ.NHEL(5)) THEN
      IF(H2.EQ.0 .OR. H2.EQ.NHEL(6)) THEN
C     BEGIN CODE
C     ----------
      CALL OXXXXX_S(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL OXXXXX_S(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX_S(P(0,3),ZERO,NHEL(3),-1*IC(3),W(1,3))
      CALL IXXXXX_S(P(0,4),ZERO,NHEL(4),-1*IC(4),W(1,4))
      CALL IXXXXX_S(P(0,5),MTA,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX_S(P(0,6),MTA,NHEL(6),+1*IC(6),W(1,6))
      CALL FFT4_5_3(W(1,5),W(1,6),GC_14,GC_15,MX2,WX2,W(1,7))
      CALL FFV4P0_3(W(1,3),W(1,1),GC_2,ZERO,ZERO,W(1,6))
      CALL FFV4P0_3(W(1,4),W(1,2),GC_2,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL VVT2_0(W(1,6),W(1,5),W(1,7),GC_49,AMP(1))
      CALL FFV5_6_3(W(1,4),W(1,2),GC_43,GC_40,MZ,WZ,W(1,8))
C     Amplitude(s) for diagram number 2
      CALL VVT2_0(W(1,6),W(1,8),W(1,7),GC_46,AMP(2))
      CALL FFV4P0_3(W(1,3),W(1,1),GC_7,ZERO,ZERO,W(1,6))
      CALL FFV4P0_3(W(1,4),W(1,2),GC_7,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 3
      CALL VVT2_0(W(1,6),W(1,9),W(1,7),GC_9,AMP(3))
      CALL FFV5_6_3(W(1,3),W(1,1),GC_43,GC_40,MZ,WZ,W(1,9))
C     Amplitude(s) for diagram number 4
      CALL VVT2_0(W(1,5),W(1,9),W(1,7),GC_46,AMP(4))
C     Amplitude(s) for diagram number 5
      CALL VVT2_0(W(1,9),W(1,8),W(1,7),GC_48,AMP(5))
      CALL FFV4P0_3(W(1,3),W(1,2),GC_2,ZERO,ZERO,W(1,9))
      CALL FFV4P0_3(W(1,4),W(1,1),GC_2,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 6
      CALL VVT2_0(W(1,9),W(1,8),W(1,7),GC_49,AMP(6))
      CALL FFV5_6_3(W(1,4),W(1,1),GC_43,GC_40,MZ,WZ,W(1,5))
C     Amplitude(s) for diagram number 7
      CALL VVT2_0(W(1,9),W(1,5),W(1,7),GC_46,AMP(7))
      CALL FFV4P0_3(W(1,3),W(1,2),GC_7,ZERO,ZERO,W(1,9))
      CALL FFV4P0_3(W(1,4),W(1,1),GC_7,ZERO,ZERO,W(1,6))
C     Amplitude(s) for diagram number 8
      CALL VVT2_0(W(1,9),W(1,6),W(1,7),GC_9,AMP(8))
      CALL FFV5_6_3(W(1,3),W(1,2),GC_43,GC_40,MZ,WZ,W(1,6))
C     Amplitude(s) for diagram number 9
      CALL VVT2_0(W(1,8),W(1,6),W(1,7),GC_46,AMP(9))
C     Amplitude(s) for diagram number 10
      CALL VVT2_0(W(1,6),W(1,5),W(1,7),GC_48,AMP(10))
      JAMP(1)=-AMP(1)-AMP(2)+1D0/6D0*AMP(3)-AMP(4)-AMP(5)+1D0/2D0
     $ *AMP(8)
      JAMP(2)=-1D0/2D0*AMP(3)+AMP(6)+AMP(7)-1D0/6D0*AMP(8)+AMP(9)
     $ +AMP(10)

      MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX = MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      MATRIX_cxcx_cxcx_s2=MATRIX
      ENDIF
      ENDIF
      END

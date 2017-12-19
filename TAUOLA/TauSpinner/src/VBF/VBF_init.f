ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE VBF_REINIT(KEYIN) 
      include 'input.inc'
      include 'coupl.inc'
      INTEGER KEYIN

C flipping defined earlier initialization variants for EW scheme used. 
C `new physics model' weight: VBFINIT(0,1)  
C  default:                   VBFINIT(0,0) 
      IF (KEYIN.GT.1) THEN
         CALL VBFINIT(0,1)
      ELSE
         CALL VBFINIT(0,0)
      ENDIF

      END

      SUBROUTINE VBFINIT(MODE,VARIANT)
C MODE=0 means use, otherwise actual initialization.
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INTEGER MODE
      INTEGER MODESAVED
      DATA    MODESAVED /0/
      INTEGER VARIANT
      INTEGER VARIANTSAVED
      DATA    VARIANTSAVED /0/
      include 'input.inc'
      include 'coupl.inc'

C primary initialization:
      IF(MODE.NE.0) THEN
         MODESAVED=MODE
         VARIANTSAVED=VARIANT
         CALL VBF_INI_OPT(MODE)
         RETURN
      ENDIF
      IF(MODESAVED.EQ.0) THEN
       WRITE(*,*) 'initialization for  VBF_INIT not chosen'
       STOP
      ENDIF

C re-initialization for `sample' EW-scheme
      IF(MODE.EQ.0.AND.VARIANT.EQ.0) THEN
        CALL VBF_INI_OPT(MODESAVED)
        RETURN
C re-initialization for `reweighted' EW-scheme
      ELSEIF(MODE.EQ.0.AND.VARIANT.EQ.1) THEN
        CALL VBF_INI_OPT(VARIANTSAVED)
        RETURN
      ENDIF

       WRITE(*,*) 'WARNING: VBF_INIT should not reach this point.'
      call testuja(1)

      END

      subroutine testuja(i)
C     This testing routine prints content of common blocks
C     used by Madgraph 
C     To activate set ifuse non-zero

      integer i,ifuse
      INCLUDE 'coupl.inc'
      INCLUDE 'input.inc'
      data ifuse /0/
      if (ifuse.eq.0) return

      write(*,*) 'a kuku w punkcie i=',i
      write(*,*) '======================'
      write(*,*) 'file couplings.inc'
      write(*,*) ' '
      write(*,*) 'COMMON/MASSES/'
      write(*,*) 'MB,MH,MT,MW,MTA,MZ'
      write(*,*)  MB,MH,MT,MW,MTA,MZ
      write(*,*) ' '
      write(*,*) 'COMMON/WIDTHS/'
      write(*,*) 'WW,WT,WZ,WH'
      write(*,*)  WW,WT,WZ,WH
      write(*,*) ' '
      write(*,*) 'COMMON/COUPLINGS/' 
      write(*,*) 'GC_1  =',GC_1 
      write(*,*) 'GC_2  =',GC_2 
      write(*,*) 'GC_3  =',GC_3 
      write(*,*) 'GC_4  =',GC_4 
      write(*,*) 'GC_10 =',GC_10
      write(*,*) 'GC_11 =',GC_11
      write(*,*) 'GC_44 =',GC_44
      write(*,*) 'GC_50 =',GC_50
      write(*,*) 'GC_59 =',GC_59
      write(*,*) 'GC_100=',GC_100
      write(*,*) 'GC_101=',GC_101
      write(*,*) 'GC_108=',GC_108
      write(*,*) 'GC_72 =',GC_72
      write(*,*) 'GC_81 =',GC_81
      write(*,*) 'GC_99 =',GC_99
      write(*,*) ' '
      write(*,*) '======================'
      write(*,*) 'file input.inc'
      write(*,*) 'COMMON/PARAMS_R/'
      write(*,*) 'SQRT__AS      =',SQRT__AS
      write(*,*) 'G__EXP__2     =',G__EXP__2
      write(*,*) 'CONJG__CKM3X3 =',CONJG__CKM3X3
      write(*,*) 'CKM3X3        =',CKM3X3
      write(*,*) 'LAMWS__EXP__2 =',LAMWS__EXP__2
      write(*,*) 'LAMWS__EXP__3 =',LAMWS__EXP__3
      write(*,*) 'MZ__EXP__2    =',MZ__EXP__2
      write(*,*) 'MZ__EXP__4    =',MZ__EXP__4
      write(*,*) 'SQRT__2       =',SQRT__2
      write(*,*) 'MH__EXP__2    =',MH__EXP__2
      write(*,*) 'AEW           =',AEW
      write(*,*) 'SQRT__AEW     =',SQRT__AEW
      write(*,*) 'EE            =',EE
      write(*,*) 'MW__EXP__2    =',MW__EXP__2
      write(*,*) 'SW2           =',SW2
      write(*,*) 'CW            =',CW
      write(*,*) 'SQRT__SW2     =',SQRT__SW2
      write(*,*) 'SW            =',SW
      write(*,*) 'G1            =',G1
      write(*,*) 'GW            =',GW
      write(*,*) 'VEV           =',VEV
      write(*,*) 'VEV__EXP__2   =',VEV__EXP__2
      write(*,*) 'LAM           =',LAM
      write(*,*) 'YB            =',YB
      write(*,*) 'YT            =',YT
      write(*,*) 'YTAU          =',YTAU
      write(*,*) 'MUH           =',MUH
      write(*,*) 'EE__EXP__2    =',EE__EXP__2
      write(*,*) 'SW__EXP__2    =',SW__EXP__2
      write(*,*) 'CW__EXP__2    =',CW__EXP__2
      write(*,*) 'AEWM1         =',AEWM1
      write(*,*) 'GF            =',GF
      write(*,*) 'AS            =',AS
      write(*,*) 'LAMWS         =',LAMWS
      write(*,*) 'AWS           =',AWS
      write(*,*) 'RHOWS         =',RHOWS
      write(*,*) 'ETAWS         =',ETAWS
      write(*,*) 'YMB           =',YMB
      write(*,*) 'YMT           =',YMT
      write(*,*) 'YMTAU         =',YMTAU
      write(*,*) ' '
      write(*,*) '======================'
      write(*,*) 'file input.inc'
      write(*,*) 'COMMON/PARAMS_C/'
      write(*,*) 'CKM1X1          =', CKM1X1
      write(*,*) 'CKM1X2          =', CKM1X2
      write(*,*) 'COMPLEXI        =', COMPLEXI
      write(*,*) 'CKM1X1          =', CKM1X1
      write(*,*) 'CKM2X1          =', CKM2X1
      write(*,*) 'CKM2X2          =', CKM2X2
      write(*,*) 'CKM2X3          =', CKM2X3
      write(*,*) 'CKM3X1          =', CKM3X1
      write(*,*) 'CKM3X2          =', CKM3X2
      write(*,*) 'CONJG__CKM1X3   =', CONJG__CKM1X3
      write(*,*) 'CONJG__CKM2X3   =', CONJG__CKM2X3
      write(*,*) 'CONJG__CKM2X1   =', CONJG__CKM2X1
      write(*,*) 'CONJG__CKM3X1   =', CONJG__CKM3X1
      write(*,*) 'CONJG__CKM2X2   =', CONJG__CKM2X2
      write(*,*) 'CONJG__CKM3X2   =', CONJG__CKM3X2
      write(*,*) 'CONJG__CKM1X1   =', CONJG__CKM1X1
      write(*,*) 'CONJG__CKM1X2   =', CONJG__CKM1X2
      write(*,*) 'I1X31           =', I1X31
      write(*,*) 'I1X32           =', I1X32
      write(*,*) 'I1X33           =', I1X33
      write(*,*) 'I1X31           =', I1X31
      write(*,*) 'I1X32           =', I1X32
      write(*,*) 'I1X33           =', I1X33
      write(*,*) 'I2X13           =', I2X13
      write(*,*) 'I2X23           =', I2X23
      write(*,*) 'I2X33           =', I2X33
      write(*,*) 'I3X31           =', I3X31
      write(*,*) 'I3X32           =', I3X32
      write(*,*) 'I3X33           =', I3X33
      write(*,*) 'I4X13           =', I4X13
      write(*,*) 'I4X23           =', I4X23
      write(*,*) 'I4X33           =', I4X33

      write(*,*) '======================'
      write(*,*) '======================'
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE VBF_INI_OPT(I)

      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INTEGER I
      include 'input.inc'
      include 'coupl.inc'

c sminputs
      aEWM1 = 1.325070D+02
      Gf    = 1.166390D-05
!      aS    = 1.180000D-01 now defined elsewhere in vbfdistr.cxx
c wolfenstein
      lamWS = 2.277360e-01 !2.253000D-01
      AWS   = 8.080000D-01
      rhoWS = 1.320000D-01
      etaWS = 3.410000D-01
C yukawa
      ymb   = 4.700000D+00
      ymt   = 1.730000D+02
      ymtau = 1.777000D+00
C widths
      WT  = 1.508336e+00 !1.491500e+00
      WZ  = 2.495200e+00 !2.441404e+00
      WW  = 2.085000e+00 !2.047600e+00
      WH  =  4.070000e-03 !6.382339e-03  ! 5.753088e-03
C masses
      MB  = 4.700000D+00
      MT  = 1.730000D+02
      MTA = 1.777000D+00
      MZ  = 9.118800D+01
      MH  = 1.250000D+02
      MW  = 80.419002    ! this value is recalculated below in 'input' section

c input
      G = 2 * DSQRT(AS*PI)  ! for the first init
      CONJG__CKM3X3 = 1.000000D+00
      CKM3X3 = 1.000000D+00
      LAMWS__EXP__2 = LAMWS**2

      CKM1X1 = 1.000000D+00-LAMWS__EXP__2/2.000000D+00
      CKM1X2 = LAMWS
      COMPLEXI = (0.000000D+00,1.000000D+00)
      LAMWS__EXP__3 = LAMWS**3

      CKM1X3 = AWS*LAMWS__EXP__3*(-(ETAWS*COMPLEXI)+RHOWS)
      CKM2X1 = -LAMWS
      CKM2X2 = 1.000000D+00-LAMWS__EXP__2/2.000000D+00
      CKM2X3 = AWS*LAMWS__EXP__2
      CKM3X1 = AWS*LAMWS__EXP__3*(1.000000D+00-ETAWS*COMPLEXI-RHOWS)
      CKM3X2 = -(AWS*LAMWS__EXP__2)

      MZ__EXP__2 = MZ**2
      MZ__EXP__4 = MZ**4
      SQRT__2 = SQRT(2.000000D+00)
      MH__EXP__2 = MH**2
 
      CONJG__CKM1X3 = CONJG(CKM1X3)
      CONJG__CKM2X3 = CONJG(CKM2X3)
      CONJG__CKM2X1 = CONJG(CKM2X1)

      CONJG__CKM3X1 = CONJG(CKM3X1)
      CONJG__CKM2X2 = CONJG(CKM2X2)
      CONJG__CKM3X2 = CONJG(CKM3X2)

      CONJG__CKM1X1 = CONJG(CKM1X1)
      CONJG__CKM1X2 = CONJG(CKM1X2)


      IF (I.EQ.1) THEN
c     IN: Gf, aEWM1, MZ; OUT: MW, SW2                             
         AEW = 1.000000D+00/AEWM1                                        
         MW    = SQRT(MZ__EXP__2/2.000000D+00+SQRT(MZ__EXP__4/4.000000D
     $           +00-(AEW*PI*MZ__EXP__2)/(GF*SQRT__2)))
         SW2   = 1.000000D+00-MW**2/MZ__EXP__2
      ELSEIF (I.EQ.2) THEN
c     IN: Gf, SW2, MZ; OUT: MW, aEWM1                                                                                                       
	 SW2   = 0.23147
         MW    = MZ*SQRT(1-SW2)
         aEWM1 = 4D0*PI/(4D0*SQRT(2D0)*Gf*MZ**2*(1-SW2)*SW2)
      ELSEIF (I.EQ.3) THEN
c     IN: Gf, MW, MZ; OUT: SW2, aEWM1                                                                                                     
	 SW2   = 1.000000D+00-MW**2/MZ__EXP__2
         aEWM1 = 4D0*PI/(4D0*SQRT(2D0)*Gf*MW**2*SW2)
      ELSEIF (i.EQ.4) THEN
c     IN: Gf, SW2, MW, MZ; OUT: aEWM1
         SW2   = 0.23147
         MW    = 80.4189
         aEWM1 = 4D0*PI/(4D0*SQRT(2D0)*Gf*MW**2*SW2)
      ELSEIF (i.EQ.5) THEN
c     IN: Gf, SW2, MW, MZ; OUT: aEWM1
         SW2   = 0.23147
         MW    = 80.4189
         aEWM1 = 4D0*PI/(4D0*SQRT(2D0)*Gf*MW**2*SW2)
C         also see later for GC_53:
C         arbitrary modification of the coupling for triple boson vertex WWZ
      ELSE
       WRITE(*,*)  'ERROR: VBF_INIT_VER: UNEXPECTED VERSION REQUIRED =',I 
       STOP
      ENDIF

      AEW = 1.000000D+00/AEWM1
      SQRT__AEW = SQRT(AEW)
      EE = 2.000000D+00*SQRT__AEW*SQRT(PI)
      MW__EXP__2 = MW**2
      CW = SQRT(1.000000D+00-SW2)
      SQRT__SW2 = SQRT(SW2)
      SW = SQRT__SW2
      G1 = EE/CW
      GW = EE/SW
      VEV = (2.000000D+00*MW*SW)/EE
      VEV__EXP__2 = VEV**2
      LAM = MH__EXP__2/(2.000000D+00*VEV__EXP__2)
      YB = (YMB*SQRT__2)/VEV
      YT = (YMT*SQRT__2)/VEV
      YTAU = (YMTAU*SQRT__2)/VEV
      MUH = SQRT(LAM*VEV__EXP__2)
      I1X31 = YB*CONJG__CKM1X3
      I1X32 = YB*CONJG__CKM2X3
      I1X33 = YB*CONJG__CKM3X3
      I2X13 = YT*CONJG__CKM3X1
      I2X23 = YT*CONJG__CKM3X2
      I2X33 = YT*CONJG__CKM3X3
      I3X31 = CKM3X1*YT
      I3X32 = CKM3X2*YT
      I3X33 = CKM3X3*YT
      I4X13 = CKM1X3*YB
      I4X23 = CKM2X3*YB
      I4X33 = CKM3X3*YB
      EE__EXP__2 = EE**2
      SW__EXP__2 = SW**2
      CW__EXP__2 = CW**2
c couplings
      GC_1 = -(EE*COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*EE*COMPLEXI)/3.000000D+00
      GC_3 = -(EE*COMPLEXI)
      GC_4 = EE*COMPLEXI
      GC_10 = -G
      GC_11 = COMPLEXI*G
      GC_44 = (CKM2X1*EE*COMPLEXI)/(SW*SQRT__2)
      GC_50 = -(CW*EE*COMPLEXI)/(2.000000D+00*SW)
      GC_51 = (CW*EE*COMPLEXI)/(2.000000D+00*SW)
      GC_53 = (CW*EE*COMPLEXI)/SW
C      arbitrary modification of the coupling for triple boson vertex WWZ
C      used for  `new physics model' weights:
      IF (i.EQ.5) THEN
       GC_53 = (CW*EE*COMPLEXI)/SW   *1.05
C       GC_53 = (CW*EE*COMPLEXI)/SW   *0.95
      ENDIF

      GC_58 = -(EE*COMPLEXI*SW)/(6.000000D+00*CW)
      GC_59 = (EE*COMPLEXI*SW)/(2.000000D+00*CW)
      GC_100 = (EE*COMPLEXI*CONJG__CKM1X1)/(SW*SQRT__2)
      GC_101 = (EE*COMPLEXI*CONJG__CKM1X2)/(SW*SQRT__2)
      GC_108 = (EE*COMPLEXI*CONJG__CKM3X3)/(SW*SQRT__2)
      GC_72 = (EE__EXP__2*COMPLEXI*VEV)/(2.000000D+00*SW__EXP__2)
      GC_81 = EE__EXP__2*COMPLEXI*VEV+(CW__EXP__2*EE__EXP__2*COMPLEXI
     $ *VEV)/(2.000000D+00*SW__EXP__2)+(EE__EXP__2*COMPLEXI*SW__EXP__2
     $ *VEV)/(2.000000D+00*CW__EXP__2)
      GC_99 = -((COMPLEXI*YTAU)/SQRT__2)

c after init
      AS = G**2/4/PI
      SQRT__AS = SQRT(AS)
      G__EXP__2 = G**2

      call testuja(2)

      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      SUBROUTINE SPIN2_REINIT(KEYIN) 
C      include 'input2.inc'
C      include 'coupl2.inc'
C      INTEGER KEYIN

C flipping defined earlier initialization variants for EW scheme used. 
C `new physics model' weight: SPIN2INIT(0,1)  
C  default:                   SPIN2INIT(0,0) 
c     IF (KEYIN.GT.3) THEN
c         CALL SPIN2INIT(0,1)
c      ELSE
       !  CALL SPIN2INIT(0,0)
c      ENDIF

C      END

      SUBROUTINE SPIN2INIT(MODE,VARIANT)
C MODE=0 means use, otherwise actual initialization.
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INTEGER MODE
      INTEGER MODESAVED
      DATA    MODESAVED /0/
      INTEGER VARIANT
      INTEGER VARIANTSAVED
      DOUBLE PRECISION consi
      DATA    VARIANTSAVED /0/
      include 'input2.inc'
      include 'coupl2.inc'

C primary initialization:
      IF(MODE.NE.0) THEN
         MODESAVED=MODE
         VARIANTSAVED=VARIANT
         CALL SPIN2_INI_OPT(MODE)
         RETURN
      ENDIF
      IF(MODESAVED.EQ.0) THEN
       WRITE(*,*) 'initialization for  SPIN2_INIT not chosen'
       STOP
      ENDIF



       WRITE(*,*) 'WARNING: SPIN2_INIT should not reach this point.'
      call testujaj(1)

      END

      subroutine testujaj(i)
C     This testing routine prints content of common blocks
C     used by Madgraph 
C     To activate set ifuse non-zero

      integer i,ifuse
      INCLUDE 'coupl2.inc'
      INCLUDE 'input2.inc'
      data ifuse /0/
      if (ifuse.eq.0) return
      WRITE(*,*)  ' External Params'
      WRITE(*,*)  ' ---------------------------------'
      WRITE(*,*)  ' '
      WRITE(*,*) 'lamws =      ', LAMWS
      WRITE(*,*) 'aws =        ', AWS
      WRITE(*,*) 'rhows =      ', RHOWS
      WRITE(*,*) 'etaws =      ', ETAWS
      WRITE(*,*) 'aEWM1 =      ', AEWM1
      WRITE(*,*) 'Gf =         ', GF
      WRITE(*,*) 'aS =         ', AS
      WRITE(*,*) 'gXtautauM =  ', GXTAUTAUM
      WRITE(*,*) 'gXtautauP =  ', GXTAUTAUP
      WRITE(*,*) 'gXqqM =      ', GXQQM
      WRITE(*,*) 'gXqqP =      ', GXQQP
      WRITE(*,*) 'gXggEven =   ', GXGGEVEN
      WRITE(*,*) 'gXggOdd =    ', GXGGODD
      WRITE(*,*) 'gXWWEven =   ', GXWWEVEN
      WRITE(*,*) 'gXWWOdd =    ', GXWWODD
      WRITE(*,*) 'gXBBEven =   ', GXBBEVEN
      WRITE(*,*) 'gXBBOdd =    ', GXBBODD
      WRITE(*,*) 'ymb =        ', YMB
      WRITE(*,*) 'ymt =        ', YMT
      WRITE(*,*) 'ymtau =      ', YMTAU
      WRITE(*,*) 'MZ =         ', MZ
      WRITE(*,*) 'MTA =        ', MTA
      WRITE(*,*) 'MT =         ', MT
      WRITE(*,*) 'MB =         ', MB
      WRITE(*,*) 'MH =         ', MH
      WRITE(*,*) 'MX2 =        ', MX2
      WRITE(*,*) 'WZ =         ', WZ
      WRITE(*,*) 'WW =         ', WW
      WRITE(*,*) 'WT =         ', WT
      WRITE(*,*) 'WH =         ', WH
      WRITE(*,*) 'WX2 =        ', WX2
      WRITE(*,*)  ' Internal Params'
      WRITE(*,*)  ' ---------------------------------'
      WRITE(*,*)  ' '
      WRITE(*,*) 'lamws__exp__2 =  ', LAMWS__EXP__2
      WRITE(*,*) 'CKM1x1 =         ', CKM1X1
      WRITE(*,*) 'CKM1x2 =         ', CKM1X2
      WRITE(*,*) 'complexi =       ', COMPLEXI
      WRITE(*,*) 'lamws__exp__3 =  ', LAMWS__EXP__3
      WRITE(*,*) 'CKM1x3 =         ', CKM1X3
      WRITE(*,*) 'CKM2x1 =         ', CKM2X1
      WRITE(*,*) 'CKM2x2 =         ', CKM2X2
      WRITE(*,*) 'CKM2x3 =         ', CKM2X3
      WRITE(*,*) 'CKM3x1 =         ', CKM3X1
      WRITE(*,*) 'CKM3x2 =         ', CKM3X2
      WRITE(*,*) 'CKM3x3 =         ', CKM3X3
      WRITE(*,*) 'MZ__exp__2 =     ', MZ__EXP__2
      WRITE(*,*) 'MZ__exp__4 =     ', MZ__EXP__4
      WRITE(*,*) 'sqrt__2 =        ', SQRT__2
      WRITE(*,*) 'MH__exp__2 =     ', MH__EXP__2
      WRITE(*,*) 'conjg__CKM1x1 =  ', CONJG__CKM1X1
      WRITE(*,*) 'conjg__CKM1x2 =  ', CONJG__CKM1X2
      WRITE(*,*) 'conjg__CKM1x3 =  ', CONJG__CKM1X3
      WRITE(*,*) 'conjg__CKM2x1 =  ', CONJG__CKM2X1
      WRITE(*,*) 'conjg__CKM2x2 =  ', CONJG__CKM2X2
      WRITE(*,*) 'conjg__CKM2x3 =  ', CONJG__CKM2X3
      WRITE(*,*) 'conjg__CKM3x1 =  ', CONJG__CKM3X1
      WRITE(*,*) 'conjg__CKM3x2 =  ', CONJG__CKM3X2
      WRITE(*,*) 'conjg__CKM3x3 =  ', CONJG__CKM3X3
      WRITE(*,*) 'aEW =            ', AEW
      WRITE(*,*) 'MW =             ', MW
      WRITE(*,*) 'sqrt__aEW =      ', SQRT__AEW
      WRITE(*,*) 'ee =             ', EE
      WRITE(*,*) 'MW__exp__2 =     ', MW__EXP__2
      WRITE(*,*) 'sw2 =            ', SW2
      WRITE(*,*) 'cw =             ', CW
      WRITE(*,*) 'sqrt__sw2 =      ', SQRT__SW2
      WRITE(*,*) 'sw =             ', SW
      WRITE(*,*) 'g1 =             ', G1
      WRITE(*,*) 'gw =             ', GW
      WRITE(*,*) 'vev =            ', VEV
      WRITE(*,*) 'vev__exp__2 =    ', VEV__EXP__2
      WRITE(*,*) 'lam =            ', LAM
      WRITE(*,*) 'yb =             ', YB
      WRITE(*,*) 'yt =             ', YT
      WRITE(*,*) 'ytau =           ', YTAU
      WRITE(*,*) 'muH =            ', MUH
      WRITE(*,*) 'ee__exp__2 =     ', EE__EXP__2
      WRITE(*,*) 'sw__exp__2 =     ', SW__EXP__2
      WRITE(*,*) 'cw__exp__2 =     ', CW__EXP__2
      WRITE(*,*)  ' Internal Params evaluated point by point'
      WRITE(*,*)  ' ----------------------------------------'
      WRITE(*,*)  ' '
      WRITE(*,*) 'sqrt__aS =       ', SQRT__AS
      WRITE(*,*) 'G__exp__2 =      ', G__EXP__2

      WRITE(*,*)  ' Couplings of spin2_w_CKM_UFO'
      WRITE(*,*)  ' ---------------------------------'
      WRITE(*,*)  ' '
      WRITE(*,*) 'GC_6 =    ', GC_6
      WRITE(*,*) 'GC_7 =    ', GC_7
      WRITE(*,*) 'GC_8 =    ', GC_8
      WRITE(*,*) 'GC_10 =   ', GC_10
      WRITE(*,*) 'GC_11 =   ', GC_11
      WRITE(*,*) 'GC_1 =    ', GC_1
      WRITE(*,*) 'GC_2 =    ', GC_2
      WRITE(*,*) 'GC_26 =   ', GC_26
      WRITE(*,*) 'GC_27 =   ', GC_27
      WRITE(*,*) 'GC_29 =   ', GC_29
      WRITE(*,*) 'GC_30 =   ', GC_30
      WRITE(*,*) 'GC_39 =   ', GC_39
      WRITE(*,*) 'GC_40 =   ', GC_40
      WRITE(*,*) 'GC_42 =   ', GC_42
      WRITE(*,*) 'GC_43 =   ', GC_43
      WRITE(*,*) 'GC_46 =   ', GC_46
      WRITE(*,*) 'GC_48 =   ', GC_48
      WRITE(*,*) 'GC_49 =   ', GC_49
      WRITE(*,*) 'GC_56 =   ', GC_56
      WRITE(*,*) 'GC_57 =   ', GC_57
      WRITE(*,*) 'GC_59 =   ', GC_59
      WRITE(*,*) 'GC_60 =   ', GC_60
      WRITE(*,*) 'GC_9 =    ', GC_9
      WRITE(*,*) 'GC_12 =   ', GC_12
      WRITE(*,*) 'GC_13 =   ', GC_13
      WRITE(*,*) 'GC_14 =   ', GC_14
      WRITE(*,*) 'GC_15 =   ', GC_15
      WRITE(*,*) 'GC_16 =   ', GC_16
      write(*,*) '======================'
      write(*,*) '======================'
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SPIN2_INI_OPT(I)

      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INTEGER I
      DOUBLE PRECISION consi
      include 'input2.inc'
      include 'coupl2.inc'

c sminputs
      aEWM1 = 1.325070e+02  ! 1.279000D+02
      Gf    =  1.166390e-05 ! 1.166370D-05
      aS    = 1.184000D-01 !now defined elsewhere in SPIN2distr.cxx
c wolfenstein
      lamWS =  2.277360e-01 !2.277359999D-01
      AWS   = 8.080000e-01! 0.0D0
      rhoWS = 1.320000e-01 !0.0D0
      etaWS = 3.410000e-01 !0.0D0
C yukawa
      ymb   = 4.700000D+00
      ymt   = 1.720000D+02
      ymtau = 1.777000D+00
C widths
      WT  = 1.508336D+00
      WZ  = 2.495200D+00
      WW  = 2.085000D+00
      WH  = 4.070000D-03 ! 5.753088e-03
      WX2 = 4.070000D-03
C masses
      MB  = 4.700000D+00
      MT  = 1.730000e+02 ! 1.720000D+02
      MTA = 1.777000D+00
      MZ  = 9.118760000D+01
      MH  = 1.250000D+02
      MX2 = 1.250000D+02
      MW  = 79.824360D+0   ! this value is recalculated below in 'input' section


c input
      G = 2 * DSQRT(AS*PI) 

      LAMWS__EXP__2 = LAMWS**2

      CKM1X1 = 1.000000D+00-LAMWS__EXP__2/2.000000D+00

      CKM1X2 = LAMWS

      COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

      LAMWS__EXP__3 = LAMWS**3

      CKM1X3 = AWS*LAMWS__EXP__3*(-(ETAWS
     $   *COMPLEXI)+RHOWS)

      CKM2X1 = -LAMWS

      CKM2X2 = 1.000000D+00-LAMWS__EXP__2/2.000000D+00

      CKM2X3 = AWS*LAMWS__EXP__2

      CKM3X1 = AWS*LAMWS__EXP__3*(1.000000D+00-ETAWS
     $   *COMPLEXI-RHOWS)

      CKM3X2 = -(AWS*LAMWS__EXP__2)

      CKM3X3 = 1.000000D+00

      MZ__EXP__2 = MZ**2

      MZ__EXP__4 = MZ**4

      SQRT__2 = SQRT(DCMPLX(2.000000D+00))

      MH__EXP__2 = MH**2

      CONJG__CKM1X1 = CONJG(DCMPLX(CKM1X1))

      CONJG__CKM1X2 = CONJG(DCMPLX(CKM1X2))

      CONJG__CKM1X3 = CONJG(DCMPLX(CKM1X3))

      CONJG__CKM2X1 = CONJG(DCMPLX(CKM2X1))

      CONJG__CKM2X2 = CONJG(DCMPLX(CKM2X2))

      CONJG__CKM2X3 = CONJG(DCMPLX(CKM2X3))

      CONJG__CKM3X1 = CONJG(DCMPLX(CKM3X1))

      CONJG__CKM3X2 = CONJG(DCMPLX(CKM3X2))

      CONJG__CKM3X3 = CONJG(DCMPLX(CKM3X3))

    
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

      MW = SQRT(DCMPLX(MZ__EXP__2/2.000000D+00
     $   +SQRT(DCMPLX(MZ__EXP__4/4.000000D+00-(AEW*PI
     $   *MZ__EXP__2)/(GF*SQRT__2)))))

      SQRT__AEW = SQRT(DCMPLX(AEW))

      EE = 2.000000D+00*SQRT__AEW*SQRT(DCMPLX(PI))

      MW__EXP__2 = MW**2

      SW2 = 1.000000D+00-MW__EXP__2/MZ__EXP__2

      CW = SQRT(DCMPLX(1.000000D+00-SW2))
 
      SQRT__SW2 = SQRT(DCMPLX(SW2))

      SW = SQRT__SW2

      G1 = EE/CW

      GW = EE/SW

      VEV = (2.000000D+00*MW*SW)/EE

      VEV__EXP__2 = VEV**2

      LAM = MH__EXP__2/(2.000000D+00*VEV__EXP__2)

      YB = (YMB*SQRT__2)/VEV

      YT = (YMT*SQRT__2)/VEV

      YTAU = (YMTAU*SQRT__2)/VEV

      MUH = SQRT(DCMPLX(LAM*VEV__EXP__2))

      EE__EXP__2 = EE**2

      SW__EXP__2 = SW**2

      CW__EXP__2 = CW**2

    
c       after init
      AS = G**2/4/PI
      SQRT__AS = SQRT(AS)
      G__EXP__2 = G**2

   
c marzieh 
      gXtautauM= 1.0D0
      gXtautauP= 1.0D0 
      gXqqM = 0.0D0 
      gXqqP = 0.0D0
      gXggEven = 0.0D0 
      gXggOdd = 0.0D0
      gXWWEven =1.0D0
      gXWWOdd =0.0D0 
      gXBBEven = 0.0D0 
      gXBBOdd = 0.0D0
      
c couplings
      GC_1 = -(EE*COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*EE*COMPLEXI)/3.000000D+00
      GC_26 = (CKM1X1*EE*COMPLEXI)/(SW*SQRT__2)
      GC_27 = (CKM1X2*EE*COMPLEXI)/(SW*SQRT__2)
      GC_29 = (CKM2X1*EE*COMPLEXI)/(SW*SQRT__2)
      GC_30 = (CKM2X2*EE*COMPLEXI)/(SW*SQRT__2)
      GC_39 = (EE*COMPLEXI*SW)/(3.000000D+00*CW)
      GC_40 = (-2.000000D+00*EE*COMPLEXI*SW)/(3.000000D+00
     $ *CW)
      GC_42 = -(CW*EE*COMPLEXI)/(2.000000D+00*SW)
     $ -(EE*COMPLEXI*SW)/(6.000000D+00*CW)
      GC_43 = (CW*EE*COMPLEXI)/(2.000000D+00*SW)
     $ -(EE*COMPLEXI*SW)/(6.000000D+00*CW)
      GC_46 = -(CW*COMPLEXI*GXBBEVEN*SW)/1.000000D+03
     $ +(CW*COMPLEXI*GXWWEVEN*SW)/1.000000D+03
      GC_48 = (CW__EXP__2*COMPLEXI*GXWWEVEN)/1.000000D+03
     $ +(COMPLEXI*GXBBEVEN*SW__EXP__2)/1.000000D+03
      GC_49 = (CW__EXP__2*COMPLEXI*GXBBEVEN)/1.000000D+03
     $ +(COMPLEXI*GXWWEVEN*SW__EXP__2)/1.000000D+03
      GC_56 = (EE*COMPLEXI*CONJG__CKM1X1)/(SW
     $ *SQRT__2)
      GC_57 = (EE*COMPLEXI*CONJG__CKM1X2)/(SW
     $ *SQRT__2)
      GC_59 = (EE*COMPLEXI*CONJG__CKM2X1)/(SW
     $ *SQRT__2)
      GC_60 = (EE*COMPLEXI*CONJG__CKM2X2)/(SW
     $ *SQRT__2)
      GC_9 = -(COMPLEXI*GXGGEVEN)/1.000000D+03
      GC_12 = -(COMPLEXI*GXQQM)/4.000000D+03
      GC_13 = -(COMPLEXI*GXQQP)/4.000000D+03
      GC_14 = -(COMPLEXI*GXTAUTAUM)/4.000000D+03
      GC_15 = -(COMPLEXI*GXTAUTAUP)/4.000000D+03
      GC_16 = (COMPLEXI*GXWWEVEN)/1.000000D+03      

      GC_6 = -G
      GC_7 = COMPLEXI*G
      GC_8 = COMPLEXI*G__EXP__2
      GC_10 = (G*GXGGEVEN)/1.000000D+03
      GC_11 = (COMPLEXI*G__EXP__2*GXGGEVEN)/1.000000D+03
      consi=1.0
      GC_1  = GC_1*consi
      GC_2  = GC_2*consi
      GC_26 = GC_26*consi
      GC_27 = GC_27*consi
      GC_29 = GC_29*consi
      GC_30 = GC_30*consi
      GC_39 = GC_39*consi
      GC_40 = GC_40*consi
      GC_42 = GC_42 *consi
      GC_43 = GC_43 *consi
      GC_46 =  GC_46*consi
      GC_48 = GC_48*consi
      GC_49 = GC_49*consi
      GC_56 =GC_56 *consi
      GC_57 = GC_57 *consi
      GC_59 = GC_59*consi
      GC_60 = GC_60*consi
      GC_9 = GC_9*consi
      GC_12 = GC_12 *consi
      GC_13 = GC_13*consi
      GC_14 = GC_14*consi
      GC_15 = GC_15*consi
      GC_16 = GC_16*consi
      GC_6  =GC_6 *consi
      GC_7  =GC_7 *consi
      GC_8  = GC_8*consi
      GC_10 =GC_10 *consi
      GC_11 = GC_11*consi
      call testujaj(2)

      END

      subroutine list_pythia
C
C     list jetset particle definitions
C
#include "pydat1.inc"

      INTEGER LNHOUT
C
      lnhout=22
      MSTU(11)=LNHOUT
      OPEN(UNIT=LNHOUT,FILE='examListPythia.lpt',STATUS='NEW')
C
C... 
      call pylist(11)
C...enable pylist(12) if you need a new pythia decay table
C      call pylist(12)
C
C
c...close output file
      CLOSE(UNIT=LNHOUT)
      return
      END

      subroutine getkf(kc,kf)
#include "pydat2.inc"
      integer kf,kc
      kf = KCHG(KC,4)
      return
      END

      subroutine getpyname(kf,idef,chap)
#include "pydat2.inc"
      integer kf
      character*16 chap
      integer idef
      idef = 1
      CALL PYNAME(KF,CHAP)
      if(CHAP.EQ.' ') idef=0
      return
      END
     
      

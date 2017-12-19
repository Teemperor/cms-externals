C###############################################################################
C
C Copyright (c) 2010 The ALOHA Development team and Contributors
C
C This file is a part of the MadGraph5_aMC@NLO project, an application which
C automatically generates Feynman diagrams and matrix elements for arbitrary
C high-energy processes in the Standard Model and beyond.
C
C It is subject to the ALOHA license which should accompany this
C distribution.
C
C###############################################################################
      subroutine ixxxxx_s(p, fmass, nhel, nsf ,fi)
c
c This subroutine computes a fermion wavefunction with the flowing-IN
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fi(6)          : fermion wavefunction               |fi>
c
      implicit none
      double complex fi(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,ip,im,nh

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )

c#ifdef HELAS_CHECK
c      double precision p2
c      double precision epsi
c      parameter( epsi = 2.0d-5 )
c      integer stdo
c      parameter( stdo = 6 )
c#endif
c
c#ifdef HELAS_CHECK
c      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
c      if ( abs(p(0))+pp.eq.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in ixxxxx is zero momentum'
c      endif
c      if ( p(0).le.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in ixxxxx has non-positive energy'
c         write(stdo,*)
c     &        '             : p(0) = ',p(0)
c      endif
c      p2 = (p(0)-pp)*(p(0)+pp)
c      if ( abs(p2-fmass**2).gt.p(0)**2*epsi ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in ixxxxx has inappropriate mass'
c         write(stdo,*)
c     &        '             : p**2 = ',p2,' : fmass**2 = ',fmass**2
c      endif
c      if (abs(nhel).ne.1) then
c         write(stdo,*) ' helas-error : nhel in ixxxxx is not -1,1'
c         write(stdo,*) '             : nhel = ',nhel
c      endif
c      if (abs(nsf).ne.1) then
c         write(stdo,*) ' helas-error : nsf in ixxxxx is not -1,1'
c         write(stdo,*) '             : nsf = ',nsf
c      endif
c#endif

      fi(1) = dcmplx(p(0),p(3))*nsf*-1
      fi(2) = dcmplx(p(1),p(2))*nsf*-1

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then

            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            ip = (1+nh)/2
            im = (1-nh)/2

            fi(3) = ip     * sqm(ip)
            fi(4) = im*nsf * sqm(ip)
            fi(5) = ip*nsf * sqm(im)
            fi(6) = im     * sqm(im)

         else

            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , p(2) )/dsqrt(rTwo*pp*pp3)
            endif

            fi(3) = sfomeg(1)*chi(im)
            fi(4) = sfomeg(1)*chi(ip)
            fi(5) = sfomeg(2)*chi(im)
            fi(6) = sfomeg(2)*chi(ip)

         endif

      else

         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fi(3) = dcmplx( rZero )
            fi(4) = dcmplx( rZero )
            fi(5) = chi(1)
            fi(6) = chi(2)
         else
            fi(3) = chi(2)
            fi(4) = chi(1)
            fi(5) = dcmplx( rZero )
            fi(6) = dcmplx( rZero )
         endif
      endif
c
      return
      end


    
      subroutine oxxxxx_s(p,fmass,nhel,nsf , fo)
c
c This subroutine computes a fermion wavefunction with the flowing-OUT
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fo(6)          : fermion wavefunction               <fo|
c
      implicit none
      double complex fo(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,nh,ip,im

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )

c#ifdef HELAS_CHECK
c      double precision p2
c      double precision epsi
c      parameter( epsi = 2.0d-5 )
c      integer stdo
c      parameter( stdo = 6 )
c#endif
c
c#ifdef HELAS_CHECK
c      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
c      if ( abs(p(0))+pp.eq.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in oxxxxx is zero momentum'
c      endif
c      if ( p(0).le.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in oxxxxx has non-positive energy'
c         write(stdo,*)
c     &        '         : p(0) = ',p(0)
c      endif
c      p2 = (p(0)-pp)*(p(0)+pp)
c      if ( abs(p2-fmass**2).gt.p(0)**2*epsi ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in oxxxxx has inappropriate mass'
c         write(stdo,*)
c     &        '             : p**2 = ',p2,' : fmass**2 = ',fmass**2
c      endif
c      if ( abs(nhel).ne.1 ) then
c         write(stdo,*) ' helas-error : nhel in oxxxxx is not -1,1'
c         write(stdo,*) '             : nhel = ',nhel
c      endif
c      if ( abs(nsf).ne.1 ) then
c         write(stdo,*) ' helas-error : nsf in oxxxxx is not -1,1'
c         write(stdo,*) '             : nsf = ',nsf
c      endif
c#endif

      fo(1) = dcmplx(p(0),p(3))*nsf
      fo(2) = dcmplx(p(1),p(2))*nsf

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then

            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            im = nhel * (1+nh)/2
            ip = nhel * -1 * ((1-nh)/2)
            fo(3) = im     * sqm(abs(ip))
            fo(4) = ip*nsf * sqm(abs(ip))
            fo(5) = im*nsf * sqm(abs(im))
            fo(6) = ip     * sqm(abs(im))
         else

            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif

            fo(3) = sfomeg(2)*chi(im)
            fo(4) = sfomeg(2)*chi(ip)
            fo(5) = sfomeg(1)*chi(im)
            fo(6) = sfomeg(1)*chi(ip)

         endif

      else

         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fo(3) = chi(1)
            fo(4) = chi(2)
            fo(5) = dcmplx( rZero )
            fo(6) = dcmplx( rZero )
         else
            fo(3) = dcmplx( rZero )
            fo(4) = dcmplx( rZero )
            fo(5) = chi(2)
            fo(6) = chi(1)
         endif

      endif
c
      return
      end

     

     

     
      subroutine vxxxxx_s(p,vmass,nhel,nsv , vc)
c
c This subroutine computes a VECTOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer nhel = -1, 0, 1: helicity      of vector boson
c                                (0 is forbidden if vmass=0.0)
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c
      implicit none
      double complex vc(6)
      double precision p(0:3),vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
      integer nhel,nsv,nsvahl

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )

c#ifdef HELAS_CHECK
c      double precision p2
c      double precision epsi
c      parameter( epsi = 2.0d-5 )
c      integer stdo
c      parameter( stdo = 6 )
c#endif
c
c#ifdef HELAS_CHECK
c      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
c      if ( abs(p(0))+pp.eq.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in vxxxxx is zero momentum'
c      endif
c      if ( p(0).le.rZero ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in vxxxxx has non-positive energy'
c         write(stdo,*)
c     &        '             : p(0) = ',p(0)
c      endif
c      p2 = (p(0)+pp)*(p(0)-pp)
c      if ( abs(p2-vmass**2).gt.p(0)**2*2.e-5 ) then
c         write(stdo,*)
c     &        ' helas-error : p(0:3) in vxxxxx has inappropriate mass'
c         write(stdo,*)
c     &        '             : p**2 = ',p2,' : vmass**2 = ',vmass**2
c      endif
c      if ( vmass.ne.rZero ) then
c         if ( abs(nhel).gt.1 ) then
c            write(stdo,*) ' helas-error : nhel in vxxxxx is not -1,0,1'
c            write(stdo,*) '             : nhel = ',nhel
c         endif
c      else
c         if ( abs(nhel).ne.1 ) then
c            write(stdo,*) ' helas-error : nhel in vxxxxx is not -1,1'
c            write(stdo,*) '             : nhel = ',nhel
c         endif
c      endif
c      if ( abs(nsv).ne.1 ) then
c         write(stdo,*) ' helas-error : nsv in vmxxxx is not -1,1'
c         write(stdo,*) '             : nsv = ',nsv
c      endif
c#endif

      sqh = dsqrt(rHalf)
      hel = dble(nhel)
      nsvahl = nsv*dabs(hel)
      pt2 = p(1)**2+p(2)**2
      pp = min(p(0),dsqrt(pt2+p(3)**2))
      pt = min(pp,dsqrt(pt2))

      vc(1) = dcmplx(p(0),p(3))*nsv
      vc(2) = dcmplx(p(1),p(2))*nsv

c#ifdef HELAS_CHECK
c nhel=4 option for scalar polarization
c      if( nhel.eq.4 ) then
c         if( vmass.eq.rZero ) then
c            vc(1) = rOne
c            vc(2) = p(1)/p(0)
c            vc(3) = p(2)/p(0)
c            vc(4) = p(3)/p(0)
c         else
c            vc(1) = p(0)/vmass
c            vc(2) = p(1)/vmass
c            vc(3) = p(2)/vmass
c            vc(4) = p(3)/vmass
c         endif
c         return
c      endif
c#endif

      if ( vmass.ne.rZero ) then

         hel0 = rOne-dabs(hel)

         if ( pp.eq.rZero ) then

            vc(3) = dcmplx( rZero )
            vc(4) = dcmplx(-hel*sqh )
            vc(5) = dcmplx( rZero , nsvahl*sqh )
            vc(6) = dcmplx( hel0 )

         else

            emp = p(0)/(vmass*pp)
            vc(3) = dcmplx( hel0*pp/vmass )
            vc(6) = dcmplx( hel0*p(3)*emp+hel*pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = p(3)/(pp*pt)*sqh*hel
               vc(4) = dcmplx( hel0*p(1)*emp-p(1)*pzpt ,
     &                         -nsvahl*p(2)/pt*sqh       )
               vc(5) = dcmplx( hel0*p(2)*emp-p(2)*pzpt ,
     &                          nsvahl*p(1)/pt*sqh       )
            else
               vc(4) = dcmplx( -hel*sqh )
               vc(5) = dcmplx( rZero , nsvahl*sign(sqh,p(3)) )
            endif

         endif

      else

         pp = p(0)
         pt = sqrt(p(1)**2+p(2)**2)
         vc(3) = dcmplx( rZero )
         vc(6) = dcmplx( hel*pt/pp*sqh )
         if ( pt.ne.rZero ) then
            pzpt = p(3)/(pp*pt)*sqh*hel
            vc(4) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
            vc(5) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
         else
            vc(4) = dcmplx( -hel*sqh )
            vc(5) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
         endif

      endif
c
      return
      end

    
   

      complex*16 function THETA_FUNCTION(cond, out_true, out_false)

      double precision cond
      complex*16 out_true, out_false

      if (cond.ge.0d0) then
        THETA_FUNCTION = out_true
      else
        THETA_FUNCTION = out_false
      endif

      return
      end

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_0(F1, F2, T3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP34
      COMPLEX*16 TMP37
      REAL*8 P1(0:3)
      COMPLEX*16 TMP36
      REAL*8 P2(0:3)
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      COMPLEX*16 TMP35
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      TMP37 = (F1(3)*(F2(5)*(P2(0)*(T3(3)+T3(6))+(P2(1)*(-1D0)*(T3(7)
     $ +T3(10))+(P2(2)*(-1D0)*(T3(11)+T3(14))-P2(3)*(T3(15)+T3(18)))))
     $ +F2(6)*(P2(0)*(T3(4)+CI*(T3(5)))+(P2(1)*(-1D0)*(T3(8)+CI*(T3(9))
     $ )+(P2(2)*(-1D0)*(T3(12)+CI*(T3(13)))-P2(3)*(T3(16)+CI*(T3(17))))
     $ )))+F1(4)*(F2(5)*(P2(0)*(T3(4)-CI*(T3(5)))+(P2(1)*(+CI*(T3(9))
     $ -T3(8))+(P2(2)*(+CI*(T3(13))-T3(12))+P2(3)*(+CI*(T3(17))-T3(16))
     $ )))+F2(6)*(P2(0)*(T3(3)-T3(6))+(P2(1)*(T3(10)-T3(7))+(P2(2)
     $ *(T3(14)-T3(11))+P2(3)*(T3(18)-T3(15)))))))
      TMP36 = (F1(3)*(F2(5)*(P1(0)*(T3(3)+T3(6))+(P1(1)*(-1D0)*(T3(7)
     $ +T3(10))+(P1(2)*(-1D0)*(T3(11)+T3(14))-P1(3)*(T3(15)+T3(18)))))
     $ +F2(6)*(P1(0)*(T3(4)+CI*(T3(5)))+(P1(1)*(-1D0)*(T3(8)+CI*(T3(9))
     $ )+(P1(2)*(-1D0)*(T3(12)+CI*(T3(13)))-P1(3)*(T3(16)+CI*(T3(17))))
     $ )))+F1(4)*(F2(5)*(P1(0)*(T3(4)-CI*(T3(5)))+(P1(1)*(+CI*(T3(9))
     $ -T3(8))+(P1(2)*(+CI*(T3(13))-T3(12))+P1(3)*(+CI*(T3(17))-T3(16))
     $ )))+F2(6)*(P1(0)*(T3(3)-T3(6))+(P1(1)*(T3(10)-T3(7))+(P1(2)
     $ *(T3(14)-T3(11))+P1(3)*(T3(18)-T3(15)))))))
      TMP35 = (F1(3)*(F2(5)*(P2(0)*(T3(3)+T3(15))+(P2(1)*(-1D0)*(T3(4)
     $ +T3(16))+(P2(2)*(-1D0)*(T3(5)+T3(17))-P2(3)*(T3(6)+T3(18)))))
     $ +F2(6)*(P2(0)*(T3(7)+CI*(T3(11)))+(P2(1)*(-1D0)*(T3(8)+CI
     $ *(T3(12)))+(P2(2)*(-1D0)*(T3(9)+CI*(T3(13)))-P2(3)*(T3(10)+CI
     $ *(T3(14)))))))+F1(4)*(F2(5)*(P2(0)*(T3(7)-CI*(T3(11)))+(P2(1)*(
     $ +CI*(T3(12))-T3(8))+(P2(2)*(+CI*(T3(13))-T3(9))+P2(3)*(+CI
     $ *(T3(14))-T3(10)))))+F2(6)*(P2(0)*(T3(3)-T3(15))+(P2(1)*(T3(16)
     $ -T3(4))+(P2(2)*(T3(17)-T3(5))+P2(3)*(T3(18)-T3(6)))))))
      TMP34 = (F1(3)*(F2(5)*(P1(0)*(T3(3)+T3(15))+(P1(1)*(-1D0)*(T3(4)
     $ +T3(16))+(P1(2)*(-1D0)*(T3(5)+T3(17))-P1(3)*(T3(6)+T3(18)))))
     $ +F2(6)*(P1(0)*(T3(7)+CI*(T3(11)))+(P1(1)*(-1D0)*(T3(8)+CI
     $ *(T3(12)))+(P1(2)*(-1D0)*(T3(9)+CI*(T3(13)))-P1(3)*(T3(10)+CI
     $ *(T3(14)))))))+F1(4)*(F2(5)*(P1(0)*(T3(7)-CI*(T3(11)))+(P1(1)*(
     $ +CI*(T3(12))-T3(8))+(P1(2)*(+CI*(T3(13))-T3(9))+P1(3)*(+CI
     $ *(T3(14))-T3(10)))))+F2(6)*(P1(0)*(T3(3)-T3(15))+(P1(1)*(T3(16)
     $ -T3(4))+(P1(2)*(T3(17)-T3(5))+P1(3)*(T3(18)-T3(6)))))))
      VERTEX = COUP*(-CI*(TMP34+TMP36)+CI*(TMP35+TMP37))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_5_0(F1, F2, T3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      REAL*8 P1(0:3)
      COMPLEX*16 COUP2
      REAL*8 P2(0:3)
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP
      COMPLEX*16 T3(*)
      CALL FFT4_0(F1,F2,T3,COUP1,VERTEX)
      CALL FFT5_0(F1,F2,T3,COUP2,TMP)
      VERTEX = VERTEX + TMP
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_1(F2, T3, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      REAL*8 W1
      COMPLEX*16 F1(6)
      COMPLEX*16 DENOM
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      F1(1) = +F2(1)+T3(1)
      F1(2) = +F2(2)+T3(2)
      P1(0) = -DBLE(F1(1))
      P1(1) = -DBLE(F1(2))
      P1(2) = -DIMAG(F1(2))
      P1(3) = -DIMAG(F1(1))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(3)= DENOM*CI * M1*(F2(5)*(P1(1)*(-1D0)*(T3(4)+T3(16)+T3(7)
     $ +T3(10))+(P1(2)*(-1D0)*(T3(5)+T3(17)+T3(11)+T3(14))+(P2(1)
     $ *(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(T3(5)+T3(17)+T3(11)+T3(14))
     $ +(T3(15)*(P1(0)+P2(3)-P2(0)-P1(3))+(T3(6)*(P2(3)+P1(0)-P1(3)
     $ -P2(0))+(T3(3)*2D0*(P1(0)-P2(0))+2D0*(T3(18)*(P2(3)-P1(3))))))))
     $ ))+F2(6)*(P1(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)*(-1D0)
     $ *(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(0)*(-1D0)*(T3(7)+T3(4)
     $ +CI*(T3(11)+T3(5)))+(P2(3)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))
     $ +(T3(12)*(P2(2)-CI*(P1(1))+CI*(P2(1))-P1(2))+(T3(9)*(P2(2)-CI
     $ *(P1(1))+CI*(P2(1))-P1(2))+(T3(13)*2D0*(-CI*(P1(2))+CI*(P2(2)))
     $ +2D0*(T3(8)*(P2(1)-P1(1)))))))))))
      F1(4)= DENOM*CI * M1*(F2(5)*(P1(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5))
     $ )+(P1(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P2(0)*(+CI*(T3(11)
     $ +T3(5))-T3(7)-T3(4))+(P2(3)*(T3(10)+T3(16)-CI*(T3(14)+T3(17)))
     $ +(T3(12)*(P2(2)-CI*(P2(1))+CI*(P1(1))-P1(2))+(T3(9)*(P2(2)-CI
     $ *(P2(1))+CI*(P1(1))-P1(2))+(T3(13)*2D0*(-CI*(P2(2))+CI*(P1(2)))
     $ +2D0*(T3(8)*(P2(1)-P1(1))))))))))+F2(6)*(P1(1)*(T3(16)+T3(10)
     $ -T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)-T3(11))+(P2(1)*(T3(4)
     $ +T3(7)-T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))
     $ +(T3(15)*(P2(0)+P2(3)-P1(0)-P1(3))+(T3(6)*(P2(3)+P2(0)-P1(3)
     $ -P1(0))+(T3(3)*2D0*(P1(0)-P2(0))+2D0*(T3(18)*(P1(3)-P2(3))))))))
     $ )))
      F1(5)= DENOM*(-CI)*(F2(6)*(P1(0)*(P1(3)*(-1D0)*(T3(10)+T3(7)
     $ +T3(16)+T3(4)+CI*(T3(14)+T3(11)+T3(17)+T3(5)))+(P1(1)*(-1D0)*(
     $ +2D0*(T3(8)+T3(3))+CI*(T3(12)+T3(9))-T3(15)-T3(6))+(P1(2)*(-1D0)
     $ *(T3(9)+T3(12)-CI*(T3(15)+T3(6))+2D0 * CI*(T3(13)+T3(3)))+(P1(0)
     $ *(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(0)*(-1D0)*(T3(7)+T3(4)+CI
     $ *(T3(11)+T3(5)))+(P2(3)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))
     $ +(P2(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))+P2(2)*(T3(9)+T3(12)
     $ +2D0 * CI*(T3(13))))))))))+(P1(1)*(P1(2)*(T3(5)+T3(11)-CI
     $ *(T3(16)+T3(10))+CI*(T3(4)+T3(7))-T3(17)-T3(14))+(P1(3)*(T3(6)
     $ +T3(15)+2D0*(T3(8))+CI*(T3(12)+T3(9))-2D0*(T3(18)))+(P1(1)
     $ *(T3(4)+T3(7)-T3(16)-T3(10))+(P2(1)*(T3(16)+T3(10)-T3(4)-T3(7))
     $ +(P2(2)*(T3(17)+T3(14)-T3(5)-T3(11))+(P2(0)*(-1D0)*(T3(15)+T3(6)
     $ -2D0*(T3(3)))-P2(3)*(T3(6)+T3(15)-2D0*(T3(18)))))))))+(P1(2)
     $ *(P1(3)*(T3(9)+T3(12)-2D0 * CI*(T3(18))+CI*(T3(6)+T3(15))+2D0 *
     $  CI*(T3(13)))+(P1(2)*(-CI*(T3(17)+T3(14))+CI*(T3(5)+T3(11)))
     $ +(P2(1)*(-CI*(T3(4)+T3(7))+CI*(T3(16)+T3(10)))+(P2(2)*(-CI
     $ *(T3(5)+T3(11))+CI*(T3(17)+T3(14)))+(P2(0)*(-1D0)*(-2D0 * CI
     $ *(T3(3))+CI*(T3(15)+T3(6)))-P2(3)*(-2D0 * CI*(T3(18))+CI*(T3(6)
     $ +T3(15))))))))+P1(3)*(P1(3)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))
     $ +(P2(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(-1D0)*(T3(10)
     $ +T3(16)+CI*(T3(14)+T3(17)))+(P2(1)*(-1D0)*(+2D0*(T3(8))+CI
     $ *(T3(12)+T3(9)))-P2(2)*(T3(9)+T3(12)+2D0 * CI*(T3(13))))))))))
     $ +F2(5)*(P1(0)*(P1(1)*(-1D0)*(T3(16)+T3(10)+2D0*(T3(4)+T3(7))-CI
     $ *(T3(11)+T3(5)))+(P1(2)*(-1D0)*(T3(17)+T3(14)+2D0*(T3(5)+T3(11))
     $ +CI*(T3(7)+T3(4)))+(P1(3)*(-2D0)*(T3(6)+T3(18)+T3(3)+T3(15))
     $ +(P2(1)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(T3(5)+T3(17)+T3(11)
     $ +T3(14))+(P1(0)*(T3(15)+T3(6)+2D0*(T3(3)))+(P2(0)*(-1D0)*(T3(15)
     $ +T3(6)+2D0*(T3(3)))+P2(3)*(T3(6)+T3(15)+2D0*(T3(18))))))))))
     $ +(P1(3)*(P1(1)*(T3(4)+T3(7)+2D0*(T3(10)+T3(16))-CI*(T3(14)
     $ +T3(17)))+(P1(2)*(T3(5)+T3(11)+2D0*(T3(14)+T3(17))+CI*(T3(10)
     $ +T3(16)))+(P2(1)*(-1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(
     $ -1D0)*(T3(5)+T3(17)+T3(11)+T3(14))+(P1(3)*(T3(6)+T3(15)+2D0
     $ *(T3(18)))+(P2(0)*(T3(15)+T3(6)+2D0*(T3(3)))-P2(3)*(T3(6)+T3(15)
     $ +2D0*(T3(18)))))))))+(P1(1)*(P1(2)*2D0*(T3(9)+T3(12)-CI*(T3(13))
     $ +CI*(T3(8)))+(P2(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P2(3)*(+CI
     $ *(T3(14)+T3(17))-T3(10)-T3(16))+(P1(1)*(-1D0)*(+CI*(T3(12)+T3(9)
     $ )-2D0*(T3(8)))+(P2(1)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))-P2(2)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(13))))))))+P1(2)*(P2(0)*(T3(11)
     $ +T3(5)+CI*(T3(7)+T3(4)))+(P2(3)*(-1D0)*(T3(14)+T3(17)+CI*(T3(10)
     $ +T3(16)))+(P1(2)*(+2D0*(T3(13))+CI*(T3(9)+T3(12)))+(P2(1)*(-1D0)
     $ *(T3(12)+T3(9)+2D0 * CI*(T3(8)))-P2(2)*(+2D0*(T3(13))+CI*(T3(9)
     $ +T3(12)))))))))))
      F1(6)= DENOM*(-CI)*(F2(5)*(P1(0)*(P1(3)*(T3(7)+T3(4)-CI*(T3(11)
     $ +T3(5))+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P1(1)*(-1D0)*(T3(15)
     $ +T3(6)+2D0*(T3(8)+T3(3))-CI*(T3(12)+T3(9)))+(P1(2)*(+CI*(T3(15)
     $ +T3(6))+2D0 * CI*(T3(13)+T3(3))-T3(9)-T3(12))+(P1(0)*(T3(7)
     $ +T3(4)-CI*(T3(11)+T3(5)))+(P2(0)*(+CI*(T3(11)+T3(5))-T3(7)-T3(4)
     $ )+(P2(3)*(T3(10)+T3(16)-CI*(T3(14)+T3(17)))+(P2(1)*(-1D0)*(+CI
     $ *(T3(12)+T3(9))-2D0*(T3(8)))+P2(2)*(T3(9)+T3(12)-2D0 * CI
     $ *(T3(13))))))))))+(P1(1)*(P1(2)*(T3(5)+T3(17)+T3(11)+T3(14)-CI
     $ *(T3(4)+T3(16)+T3(7)+T3(10)))+(P1(3)*(T3(6)+T3(15)+2D0*(T3(18))
     $ +CI*(T3(12)+T3(9))-2D0*(T3(8)))+(P1(1)*(T3(4)+T3(16)+T3(7)
     $ +T3(10))+(P2(1)*(-1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(-1D0)
     $ *(T3(5)+T3(17)+T3(11)+T3(14))+(P2(0)*(T3(15)+T3(6)+2D0*(T3(3)))
     $ -P2(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P1(2)*(P1(3)*(-1D0)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(13))+CI*(T3(6)+T3(15))+2D0 * CI
     $ *(T3(18)))+(P1(2)*(-1D0)*(+CI*(T3(5)+T3(17)+T3(11)+T3(14)))
     $ +(P2(1)*(+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P2(2)*(+CI*(T3(5)
     $ +T3(17)+T3(11)+T3(14)))+(P2(0)*(-1D0)*(+CI*(T3(15)+T3(6))+2D0 *
     $  CI*(T3(3)))+P2(3)*(+CI*(T3(6)+T3(15))+2D0 * CI*(T3(18))))))))
     $ +P1(3)*(P1(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P2(0)*(+CI
     $ *(T3(11)+T3(5))-T3(7)-T3(4))+(P2(3)*(T3(10)+T3(16)-CI*(T3(14)
     $ +T3(17)))+(P2(1)*(-1D0)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))+P2(2)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(13))))))))))+F2(6)*(P1(0)*(P1(1)*(
     $ -1D0)*(+2D0*(T3(4)+T3(7))+CI*(T3(11)+T3(5))-T3(16)-T3(10))
     $ +(P1(2)*(T3(17)+T3(14)+CI*(T3(7)+T3(4))-2D0*(T3(5)+T3(11)))
     $ +(P1(3)*2D0*(T3(18)+T3(3)-T3(6)-T3(15))+(P2(1)*(T3(4)+T3(7)
     $ -T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(P1(0)*(
     $ -1D0)*(T3(15)+T3(6)-2D0*(T3(3)))+(P2(0)*(T3(15)+T3(6)-2D0*(T3(3)
     $ ))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18))))))))))+(P1(3)*(P1(1)*(+2D0
     $ *(T3(10)+T3(16))+CI*(T3(14)+T3(17))-T3(4)-T3(7))+(P1(2)*(-1D0)
     $ *(T3(5)+T3(11)+CI*(T3(10)+T3(16))-2D0*(T3(14)+T3(17)))+(P2(1)
     $ *(T3(4)+T3(7)-T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))
     $ +(P1(3)*(-1D0)*(T3(6)+T3(15)-2D0*(T3(18)))+(P2(0)*(T3(15)+T3(6)
     $ -2D0*(T3(3)))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18)))))))))+(P1(1)
     $ *(P1(2)*2D0*(T3(9)+T3(12)-CI*(T3(8))+CI*(T3(13)))+(P2(0)*(T3(7)
     $ +T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(-1D0)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P1(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))+(P2(1)
     $ *(-1D0)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))-P2(2)*(T3(9)+T3(12)
     $ +2D0 * CI*(T3(13))))))))+P1(2)*(P2(0)*(T3(11)+T3(5)-CI*(T3(7)
     $ +T3(4)))+(P2(3)*(+CI*(T3(10)+T3(16))-T3(14)-T3(17))+(P1(2)*(
     $ -1D0)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))+(P2(1)*(-1D0)*(T3(12)
     $ +T3(9)-2D0 * CI*(T3(8)))+P2(2)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))
     $ ))))))))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_5_1(F2, T3, COUP1, COUP2, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      REAL*8 W1
      COMPLEX*16 F1(6)
      COMPLEX*16 COUP1
      COMPLEX*16 DENOM
      COMPLEX*16 COUP2
      INTEGER*4 I
      COMPLEX*16 T3(*)
      COMPLEX*16 FTMP(6)
      CALL FFT4_1(F2,T3,COUP1,M1,W1,F1)
      CALL FFT5_1(F2,T3,COUP2,M1,W1,FTMP)
      DO I = 3, 6
        F1(I) = F1(I) + FTMP(I)
      ENDDO
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_2(F1, T3, COUP, M2, W2,F2)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(6)
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      REAL*8 W2
      COMPLEX*16 F1(*)
      REAL*8 M2
      COMPLEX*16 DENOM
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      F2(1) = +F1(1)+T3(1)
      F2(2) = +F1(2)+T3(2)
      P2(0) = -DBLE(F2(1))
      P2(1) = -DBLE(F2(2))
      P2(2) = -DIMAG(F2(2))
      P2(3) = -DIMAG(F2(1))
      DENOM = COUP/(P2(0)**2-P2(1)**2-P2(2)**2-P2(3)**2 - M2 * (M2 -CI
     $ * W2))
      F2(3)= DENOM*CI*(F1(4)*(P2(0)*(P2(3)*(T3(10)+T3(7)+T3(16)+T3(4)
     $ -CI*(T3(14)+T3(11)+T3(17)+T3(5)))+(P2(1)*(-1D0)*(T3(15)+T3(6)
     $ +CI*(T3(12)+T3(9))-2D0*(T3(8)+T3(3)))+(P2(2)*(T3(9)+T3(12)-2D0 
     $ * CI*(T3(13)+T3(3))+CI*(T3(15)+T3(6)))+(P1(0)*(T3(7)+T3(4)-CI
     $ *(T3(11)+T3(5)))+(P1(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))
     $ +(P2(0)*(+CI*(T3(11)+T3(5))-T3(7)-T3(4))+(P1(1)*(+CI*(T3(12)
     $ +T3(9))-2D0*(T3(8)))-P1(2)*(T3(9)+T3(12)-2D0 * CI*(T3(13))))))))
     $ ))+(P2(1)*(P2(2)*(T3(17)+T3(14)-CI*(T3(16)+T3(10))+CI*(T3(4)
     $ +T3(7))-T3(5)-T3(11))+(P2(3)*(+2D0*(T3(18))+CI*(T3(12)+T3(9))
     $ -2D0*(T3(8))-T3(6)-T3(15))+(P1(1)*(T3(4)+T3(7)-T3(16)-T3(10))
     $ +(P1(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(P2(1)*(T3(16)+T3(10)-T3(4)
     $ -T3(7))+(P1(0)*(T3(15)+T3(6)-2D0*(T3(3)))+P1(3)*(T3(6)+T3(15)
     $ -2D0*(T3(18)))))))))+(P2(2)*(P2(3)*(-2D0 * CI*(T3(18))+CI*(T3(6)
     $ +T3(15))+2D0 * CI*(T3(13))-T3(9)-T3(12))+(P1(1)*(-CI*(T3(4)
     $ +T3(7))+CI*(T3(16)+T3(10)))+(P1(2)*(-CI*(T3(5)+T3(11))+CI
     $ *(T3(17)+T3(14)))+(P2(2)*(-CI*(T3(17)+T3(14))+CI*(T3(5)+T3(11)))
     $ +(P1(0)*(-1D0)*(-2D0 * CI*(T3(3))+CI*(T3(15)+T3(6)))-P1(3)*(
     $ -2D0 * CI*(T3(18))+CI*(T3(6)+T3(15))))))))+P2(3)*(P1(0)*(+CI
     $ *(T3(11)+T3(5))-T3(7)-T3(4))+(P1(3)*(T3(10)+T3(16)-CI*(T3(14)
     $ +T3(17)))+(P2(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P1(1)*(
     $ -1D0)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))+P1(2)*(T3(9)+T3(12)-2D0 
     $ * CI*(T3(13))))))))))+F1(3)*(P2(0)*(P2(1)*(T3(16)+T3(10)+2D0
     $ *(T3(4)+T3(7))+CI*(T3(11)+T3(5)))+(P2(2)*(T3(17)+T3(14)+2D0
     $ *(T3(5)+T3(11))-CI*(T3(7)+T3(4)))+(P1(1)*(-1D0)*(T3(4)+T3(16)
     $ +T3(7)+T3(10))+(P1(2)*(-1D0)*(T3(5)+T3(17)+T3(11)+T3(14))+(P2(3)
     $ *2D0*(T3(6)+T3(18)+T3(3)+T3(15))+(P1(0)*(T3(15)+T3(6)+2D0*(T3(3)
     $ ))+(P1(3)*(-1D0)*(T3(6)+T3(15)+2D0*(T3(18)))-P2(0)*(T3(15)+T3(6)
     $ +2D0*(T3(3))))))))))+(P2(3)*(P2(1)*(-1D0)*(T3(4)+T3(7)+2D0
     $ *(T3(10)+T3(16))+CI*(T3(14)+T3(17)))+(P2(2)*(-1D0)*(T3(5)+T3(11)
     $ -CI*(T3(10)+T3(16))+2D0*(T3(14)+T3(17)))+(P1(1)*(T3(4)+T3(16)
     $ +T3(7)+T3(10))+(P1(2)*(T3(5)+T3(17)+T3(11)+T3(14))+(P1(0)*(-1D0)
     $ *(T3(15)+T3(6)+2D0*(T3(3)))+(P1(3)*(T3(6)+T3(15)+2D0*(T3(18)))
     $ -P2(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P2(1)*(P1(0)*(-1D0)
     $ *(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P2(2)*(-2D0)*(T3(9)+T3(12)-CI*(T3(8))+CI
     $ *(T3(13)))+(P1(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))+(P1(2)*(T3(9)
     $ +T3(12)+2D0 * CI*(T3(13)))-P2(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9))
     $ ))))))+P2(2)*(P1(0)*(+CI*(T3(7)+T3(4))-T3(11)-T3(5))+(P1(3)
     $ *(T3(14)+T3(17)-CI*(T3(10)+T3(16)))+(P1(1)*(T3(12)+T3(9)-2D0 *
     $  CI*(T3(8)))+(P1(2)*(-1D0)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))
     $ +P2(2)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))))))))))
      F2(4)= DENOM*CI*(F1(3)*(P2(0)*(P2(3)*(T3(10)+T3(16)-CI*(T3(11)
     $ +T3(5))+CI*(T3(14)+T3(17))-T3(7)-T3(4))+(P2(1)*(T3(15)+T3(6)
     $ +2D0*(T3(8)+T3(3))+CI*(T3(12)+T3(9)))+(P2(2)*(T3(9)+T3(12)+CI
     $ *(T3(15)+T3(6))+2D0 * CI*(T3(13)+T3(3)))+(P1(0)*(T3(7)+T3(4)+CI
     $ *(T3(11)+T3(5)))+(P1(3)*(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17))
     $ )+(P2(0)*(-1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(1)*(-1D0)*(
     $ +2D0*(T3(8))+CI*(T3(12)+T3(9)))-P1(2)*(T3(9)+T3(12)+2D0 * CI
     $ *(T3(13))))))))))+(P2(1)*(P2(2)*(-1D0)*(T3(5)+T3(17)+T3(11)
     $ +T3(14)+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P2(3)*(+2D0*(T3(8))+CI
     $ *(T3(12)+T3(9))-2D0*(T3(18))-T3(6)-T3(15))+(P1(1)*(T3(4)+T3(16)
     $ +T3(7)+T3(10))+(P1(2)*(T3(5)+T3(17)+T3(11)+T3(14))+(P2(1)*(-1D0)
     $ *(T3(4)+T3(16)+T3(7)+T3(10))+(P1(0)*(-1D0)*(T3(15)+T3(6)+2D0
     $ *(T3(3)))+P1(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P2(2)*(P2(3)
     $ *(T3(9)+T3(12)-CI*(T3(6)+T3(15))-2D0 * CI*(T3(18))+2D0 * CI
     $ *(T3(13)))+(P1(1)*(+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P1(2)*(+CI
     $ *(T3(5)+T3(17)+T3(11)+T3(14)))+(P2(2)*(-1D0)*(+CI*(T3(5)+T3(17)
     $ +T3(11)+T3(14)))+(P1(0)*(-1D0)*(+CI*(T3(15)+T3(6))+2D0 * CI
     $ *(T3(3)))+P1(3)*(+CI*(T3(6)+T3(15))+2D0 * CI*(T3(18))))))))
     $ +P2(3)*(P1(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)*(-1D0)
     $ *(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(3)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P1(1)*(-1D0)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))
     $ -P1(2)*(T3(9)+T3(12)+2D0 * CI*(T3(13))))))))))+F1(4)*(P2(0)
     $ *(P2(1)*(-1D0)*(T3(16)+T3(10)+CI*(T3(11)+T3(5))-2D0*(T3(4)+T3(7)
     $ ))+(P2(2)*(+2D0*(T3(5)+T3(11))+CI*(T3(7)+T3(4))-T3(17)-T3(14))
     $ +(P1(1)*(T3(16)+T3(10)-T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)
     $ -T3(11))+(P2(3)*2D0*(T3(6)+T3(15)-T3(18)-T3(3))+(P1(0)*(-1D0)
     $ *(T3(15)+T3(6)-2D0*(T3(3)))+(P1(3)*(-1D0)*(T3(6)+T3(15)-2D0
     $ *(T3(18)))+P2(0)*(T3(15)+T3(6)-2D0*(T3(3))))))))))+(P2(3)*(P2(1)
     $ *(T3(4)+T3(7)+CI*(T3(14)+T3(17))-2D0*(T3(10)+T3(16)))+(P2(2)*(
     $ -1D0)*(+2D0*(T3(14)+T3(17))+CI*(T3(10)+T3(16))-T3(5)-T3(11))
     $ +(P1(1)*(T3(16)+T3(10)-T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)
     $ -T3(11))+(P1(0)*(-1D0)*(T3(15)+T3(6)-2D0*(T3(3)))+(P1(3)*(-1D0)
     $ *(T3(6)+T3(15)-2D0*(T3(18)))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18)))))
     $ ))))+(P2(1)*(P1(0)*(+CI*(T3(11)+T3(5))-T3(7)-T3(4))+(P1(3)
     $ *(T3(10)+T3(16)-CI*(T3(14)+T3(17)))+(P2(2)*(-2D0)*(T3(9)+T3(12)
     $ -CI*(T3(13))+CI*(T3(8)))+(P1(1)*(-1D0)*(+CI*(T3(12)+T3(9))-2D0
     $ *(T3(8)))+(P1(2)*(T3(9)+T3(12)-2D0 * CI*(T3(13)))+P2(1)*(+CI
     $ *(T3(12)+T3(9))-2D0*(T3(8))))))))+P2(2)*(P1(0)*(-1D0)*(T3(11)
     $ +T3(5)+CI*(T3(7)+T3(4)))+(P1(3)*(T3(14)+T3(17)+CI*(T3(10)+T3(16)
     $ ))+(P1(1)*(T3(12)+T3(9)+2D0 * CI*(T3(8)))+(P1(2)*(+2D0*(T3(13))
     $ +CI*(T3(9)+T3(12)))-P2(2)*(+2D0*(T3(13))+CI*(T3(9)+T3(12))))))))
     $ )))
      F2(5)= DENOM*CI * M2*(F1(3)*(P1(1)*(-1D0)*(T3(4)+T3(16)+T3(7)
     $ +T3(10))+(P1(2)*(-1D0)*(T3(5)+T3(17)+T3(11)+T3(14))+(P2(1)
     $ *(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(T3(5)+T3(17)+T3(11)+T3(14))
     $ +(T3(15)*(P1(0)+P2(3)-P2(0)-P1(3))+(T3(6)*(P2(3)+P1(0)-P1(3)
     $ -P2(0))+(T3(3)*2D0*(P1(0)-P2(0))+2D0*(T3(18)*(P2(3)-P1(3))))))))
     $ ))+F1(4)*(P1(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P1(3)*(+CI
     $ *(T3(14)+T3(17))-T3(10)-T3(16))+(P2(0)*(+CI*(T3(11)+T3(5))-T3(7)
     $ -T3(4))+(P2(3)*(T3(10)+T3(16)-CI*(T3(14)+T3(17)))+(T3(12)*(P2(2)
     $ -CI*(P2(1))+CI*(P1(1))-P1(2))+(T3(9)*(P2(2)-CI*(P2(1))+CI*(P1(1)
     $ )-P1(2))+(T3(13)*2D0*(-CI*(P2(2))+CI*(P1(2)))+2D0*(T3(8)*(P2(1)
     $ -P1(1)))))))))))
      F2(6)= DENOM*CI * M2*(F1(3)*(P1(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5))
     $ )+(P1(3)*(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(0)*(-1D0)
     $ *(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(T3(12)*(P2(2)-CI*(P1(1))+CI*(P2(1))-P1(2))
     $ +(T3(9)*(P2(2)-CI*(P1(1))+CI*(P2(1))-P1(2))+(T3(13)*2D0*(-CI
     $ *(P1(2))+CI*(P2(2)))+2D0*(T3(8)*(P2(1)-P1(1))))))))))+F1(4)
     $ *(P1(1)*(T3(16)+T3(10)-T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)
     $ -T3(11))+(P2(1)*(T3(4)+T3(7)-T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)
     $ -T3(17)-T3(14))+(T3(15)*(P2(0)+P2(3)-P1(0)-P1(3))+(T3(6)*(P2(3)
     $ +P2(0)-P1(3)-P1(0))+(T3(3)*2D0*(P1(0)-P2(0))+2D0*(T3(18)*(P1(3)
     $ -P2(3)))))))))))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_5_2(F1, T3, COUP1, COUP2, M2, W2,F2)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(6)
      REAL*8 P1(0:3)
      COMPLEX*16 COUP2
      REAL*8 P2(0:3)
      REAL*8 W2
      COMPLEX*16 F1(*)
      REAL*8 M2
      COMPLEX*16 DENOM
      COMPLEX*16 COUP1
      INTEGER*4 I
      COMPLEX*16 T3(*)
      COMPLEX*16 FTMP(6)
      CALL FFT4_2(F1,T3,COUP1,M2,W2,F2)
      CALL FFT5_2(F1,T3,COUP2,M2,W2,FTMP)
      DO I = 3, 6
        F2(I) = F2(I) + FTMP(I)
      ENDDO
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_3(F1, F2, COUP, M3, W3,T3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP33
      COMPLEX*16 TMP39
      COMPLEX*16 COUP
      REAL*8 P1(0:3)
      REAL*8 W3
      REAL*8 P2(0:3)
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 TMP21
      COMPLEX*16 F2(*)
      REAL*8 OM3
      COMPLEX*16 DENOM
      COMPLEX*16 T3(18)
      COMPLEX*16 TMP18
      COMPLEX*16 TMP38
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      T3(1) = +F1(1)+F2(1)
      T3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(T3(1))
      P3(1) = -DBLE(T3(2))
      P3(2) = -DIMAG(T3(2))
      P3(3) = -DIMAG(T3(1))
      TMP33 = (F1(3)*(F2(5)*(P3(0)+P3(3))+F2(6)*(P3(1)+CI*(P3(2))))
     $ +F1(4)*(F2(5)*(P3(1)-CI*(P3(2)))+F2(6)*(P3(0)-P3(3))))
      TMP39 = (F1(3)*(F2(5)*(P2(0)+P2(3))+F2(6)*(P2(1)+CI*(P2(2))))
     $ +F1(4)*(F2(5)*(P2(1)-CI*(P2(2)))+F2(6)*(P2(0)-P2(3))))
      TMP38 = (F1(3)*(F2(5)*(P1(0)+P1(3))+F2(6)*(P1(1)+CI*(P1(2))))
     $ +F1(4)*(F2(5)*(P1(1)-CI*(P1(2)))+F2(6)*(P1(0)-P1(3))))
      TMP18 = (P1(0)*P3(0)-P1(1)*P3(1)-P1(2)*P3(2)-P1(3)*P3(3))
      TMP21 = (P2(0)*P3(0)-P2(1)*P3(1)-P2(2)*P3(2)-P2(3)*P3(3))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      T3(3)= DENOM*2D0 * CI*(OM3*(P3(0)*(P3(0)*(OM3*2D0/3D0 * TMP33
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38)))+(F1(3)*F2(5)
     $ *(TMP21-TMP18)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(0)-P1(0)))))
     $ +1D0/3D0*(TMP33*(TMP18-TMP21)))+(F1(3)*F2(5)*(P1(0)-P2(0))
     $ +(F1(4)*F2(6)*(P1(0)-P2(0))+(-1D0/3D0*(TMP38)+1D0/3D0*(TMP39))))
     $ )
      T3(4)= DENOM*CI*(OM3*(P3(0)*(P3(1)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(TMP18
     $ -TMP21)+(F1(4)*F2(5)*(TMP18-TMP21)+TMP33*(P2(1)-P1(1)))))+P3(1)
     $ *(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33
     $ *(P2(0)-P1(0)))))+(F1(3)*(F2(5)*(P1(1)-P2(1))+F2(6)*(P2(0)-P1(0)
     $ ))+F1(4)*(F2(5)*(P2(0)-P1(0))+F2(6)*(P1(1)-P2(1)))))
      T3(5)= DENOM*CI*(OM3*(P3(0)*(P3(2)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(-CI
     $ *(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)+CI*(TMP21))
     $ +TMP33*(P2(2)-P1(2)))))+P3(2)*(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)
     $ *F2(6)*(TMP21-TMP18)+TMP33*(P2(0)-P1(0)))))+(F1(3)*(F2(5)*(P1(2)
     $ -P2(2))+F2(6)*(-CI*(P1(0))+CI*(P2(0))))+F1(4)*(F2(5)*(-CI*(P2(0)
     $ )+CI*(P1(0)))+F2(6)*(P1(2)-P2(2)))))
      T3(6)= DENOM*CI*(OM3*(P3(0)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33
     $ *(P2(0)-P1(0)))))+(F1(3)*F2(5)*(P1(3)+P2(0)-P1(0)-P2(3))+F1(4)
     $ *F2(6)*(P1(3)+P1(0)-P2(3)-P2(0))))
      T3(7)= DENOM*CI*(OM3*(P3(0)*(P3(1)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(TMP18
     $ -TMP21)+(F1(4)*F2(5)*(TMP18-TMP21)+TMP33*(P2(1)-P1(1)))))+P3(1)
     $ *(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33
     $ *(P2(0)-P1(0)))))+(F1(3)*(F2(5)*(P1(1)-P2(1))+F2(6)*(P2(0)-P1(0)
     $ ))+F1(4)*(F2(5)*(P2(0)-P1(0))+F2(6)*(P1(1)-P2(1)))))
      T3(8)= DENOM*2D0 * CI*(OM3*(P3(1)*(P3(1)*(OM3*2D0/3D0 * TMP33
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38)))+(F1(3)*F2(6)
     $ *(TMP18-TMP21)+(F1(4)*F2(5)*(TMP18-TMP21)+TMP33*(P2(1)-P1(1)))))
     $ +1D0/3D0*(TMP33*(TMP21-TMP18)))+(F1(3)*F2(6)*(P2(1)-P1(1))
     $ +(F1(4)*F2(5)*(P2(1)-P1(1))+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38))))
     $ )
      T3(9)= DENOM*CI*(OM3*(P3(1)*(P3(2)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(-CI
     $ *(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)+CI*(TMP21))
     $ +TMP33*(P2(2)-P1(2)))))+P3(2)*(F1(3)*F2(6)*(TMP18-TMP21)+(F1(4)
     $ *F2(5)*(TMP18-TMP21)+TMP33*(P2(1)-P1(1)))))+(F1(3)*F2(6)*(P2(2)
     $ -CI*(P1(1))+CI*(P2(1))-P1(2))+F1(4)*F2(5)*(P2(2)-CI*(P2(1))+CI
     $ *(P1(1))-P1(2))))
      T3(10)= DENOM*CI*(OM3*(P3(1)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(6)*(TMP18-TMP21)+(F1(4)*F2(5)*(TMP18-TMP21)+TMP33
     $ *(P2(1)-P1(1)))))+(F1(3)*(F2(5)*(P2(1)-P1(1))+F2(6)*(P2(3)-P1(3)
     $ ))+F1(4)*(F2(5)*(P2(3)-P1(3))+F2(6)*(P1(1)-P2(1)))))
      T3(11)= DENOM*CI*(OM3*(P3(0)*(P3(2)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(-CI
     $ *(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)+CI*(TMP21))
     $ +TMP33*(P2(2)-P1(2)))))+P3(2)*(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)
     $ *F2(6)*(TMP21-TMP18)+TMP33*(P2(0)-P1(0)))))+(F1(3)*(F2(5)*(P1(2)
     $ -P2(2))+F2(6)*(-CI*(P1(0))+CI*(P2(0))))+F1(4)*(F2(5)*(-CI*(P2(0)
     $ )+CI*(P1(0)))+F2(6)*(P1(2)-P2(2)))))
      T3(12)= DENOM*CI*(OM3*(P3(1)*(P3(2)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(6)*(-CI
     $ *(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)+CI*(TMP21))
     $ +TMP33*(P2(2)-P1(2)))))+P3(2)*(F1(3)*F2(6)*(TMP18-TMP21)+(F1(4)
     $ *F2(5)*(TMP18-TMP21)+TMP33*(P2(1)-P1(1)))))+(F1(3)*F2(6)*(P2(2)
     $ -CI*(P1(1))+CI*(P2(1))-P1(2))+F1(4)*F2(5)*(P2(2)-CI*(P2(1))+CI
     $ *(P1(1))-P1(2))))
      T3(13)= DENOM*2D0 * CI*(OM3*(P3(2)*(P3(2)*(OM3*2D0/3D0 * TMP33
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38)))+(F1(3)*F2(6)
     $ *(-CI*(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)+CI*(TMP21))
     $ +TMP33*(P2(2)-P1(2)))))+1D0/3D0*(TMP33*(TMP21-TMP18)))+(F1(3)
     $ *F2(6)*(-CI*(P1(2))+CI*(P2(2)))+(F1(4)*F2(5)*(-CI*(P2(2))+CI
     $ *(P1(2)))+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38)))))
      T3(14)= DENOM*CI*(OM3*(P3(2)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(6)*(-CI*(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)
     $ +CI*(TMP21))+TMP33*(P2(2)-P1(2)))))+(F1(3)*(F2(5)*(P2(2)-P1(2))
     $ +F2(6)*(-CI*(P1(3))+CI*(P2(3))))+F1(4)*(F2(5)*(-CI*(P2(3))+CI
     $ *(P1(3)))+F2(6)*(P1(2)-P2(2)))))
      T3(15)= DENOM*CI*(OM3*(P3(0)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(5)*(TMP21-TMP18)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33
     $ *(P2(0)-P1(0)))))+(F1(3)*F2(5)*(P1(3)+P2(0)-P1(0)-P2(3))+F1(4)
     $ *F2(6)*(P1(0)+P1(3)-P2(0)-P2(3))))
      T3(16)= DENOM*CI*(OM3*(P3(1)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(6)*(TMP18-TMP21)+(F1(4)*F2(5)*(TMP18-TMP21)+TMP33
     $ *(P2(1)-P1(1)))))+(F1(3)*(F2(5)*(P2(1)-P1(1))+F2(6)*(P2(3)-P1(3)
     $ ))+F1(4)*(F2(5)*(P2(3)-P1(3))+F2(6)*(P1(1)-P2(1)))))
      T3(17)= DENOM*CI*(OM3*(P3(2)*(P3(3)*(OM3*4D0/3D0 * TMP33*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP39)+2D0/3D0*(TMP38)))+(F1(3)*F2(5)*(TMP18
     $ -TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(3)*F2(6)*(-CI*(TMP21)+CI*(TMP18))+(F1(4)*F2(5)*(-CI*(TMP18)
     $ +CI*(TMP21))+TMP33*(P2(2)-P1(2)))))+(F1(3)*(F2(5)*(P2(2)-P1(2))
     $ +F2(6)*(-CI*(P1(3))+CI*(P2(3))))+F1(4)*(F2(5)*(-CI*(P2(3))+CI
     $ *(P1(3)))+F2(6)*(P1(2)-P2(2)))))
      T3(18)= DENOM*2D0 * CI*(OM3*(P3(3)*(P3(3)*(OM3*2D0/3D0 * TMP33
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38)))+(F1(3)*F2(5)
     $ *(TMP18-TMP21)+(F1(4)*F2(6)*(TMP21-TMP18)+TMP33*(P2(3)-P1(3)))))
     $ +1D0/3D0*(TMP33*(TMP21-TMP18)))+(F1(3)*F2(5)*(P2(3)-P1(3))
     $ +(F1(4)*F2(6)*(P1(3)-P2(3))+(-1D0/3D0*(TMP39)+1D0/3D0*(TMP38))))
     $ )
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjM(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjM(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjM(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFT4_5_3(F1, F2, COUP1, COUP2, M3, W3,T3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      REAL*8 P1(0:3)
      REAL*8 W3
      REAL*8 P2(0:3)
      COMPLEX*16 TTMP(18)
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 COUP1
      REAL*8 OM3
      INTEGER*4 I
      COMPLEX*16 DENOM
      COMPLEX*16 COUP2
      COMPLEX*16 T3(18)
      CALL FFT4_3(F1,F2,COUP1,M3,W3,T3)
      CALL FFT5_3(F1,F2,COUP2,M3,W3,TTMP)
      DO I = 3, 18
        T3(I) = T3(I) + TTMP(I)
      ENDDO
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjP(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjP(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjP(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFT5_1(F2, T3, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      REAL*8 W1
      COMPLEX*16 F1(6)
      COMPLEX*16 DENOM
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      F1(1) = +F2(1)+T3(1)
      F1(2) = +F2(2)+T3(2)
      P1(0) = -DBLE(F1(1))
      P1(1) = -DBLE(F1(2))
      P1(2) = -DIMAG(F1(2))
      P1(3) = -DIMAG(F1(1))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(3)= DENOM*(-CI)*(F2(4)*(P1(0)*(P1(3)*(T3(10)+T3(16)-CI*(T3(11)
     $ +T3(5))+CI*(T3(14)+T3(17))-T3(7)-T3(4))+(P1(1)*(T3(15)+T3(6)
     $ +2D0*(T3(8)+T3(3))+CI*(T3(12)+T3(9)))+(P1(2)*(T3(9)+T3(12)+CI
     $ *(T3(15)+T3(6))+2D0 * CI*(T3(13)+T3(3)))+(P1(0)*(-1D0)*(T3(7)
     $ +T3(4)+CI*(T3(11)+T3(5)))+(P2(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))
     $ +(P2(3)*(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(1)*(-1D0)
     $ *(+2D0*(T3(8))+CI*(T3(12)+T3(9)))-P2(2)*(T3(9)+T3(12)+2D0 * CI
     $ *(T3(13))))))))))+(P1(1)*(P1(2)*(-1D0)*(T3(5)+T3(17)+T3(11)
     $ +T3(14)+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P1(3)*(+2D0*(T3(8))+CI
     $ *(T3(12)+T3(9))-2D0*(T3(18))-T3(6)-T3(15))+(P1(1)*(-1D0)*(T3(4)
     $ +T3(16)+T3(7)+T3(10))+(P2(1)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)
     $ *(T3(5)+T3(17)+T3(11)+T3(14))+(P2(0)*(-1D0)*(T3(15)+T3(6)+2D0
     $ *(T3(3)))+P2(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P1(2)*(P1(3)
     $ *(T3(9)+T3(12)-CI*(T3(6)+T3(15))-2D0 * CI*(T3(18))+2D0 * CI
     $ *(T3(13)))+(P1(2)*(-1D0)*(+CI*(T3(5)+T3(17)+T3(11)+T3(14)))
     $ +(P2(1)*(+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P2(2)*(+CI*(T3(5)
     $ +T3(17)+T3(11)+T3(14)))+(P2(0)*(-1D0)*(+CI*(T3(15)+T3(6))+2D0 *
     $  CI*(T3(3)))+P2(3)*(+CI*(T3(6)+T3(15))+2D0 * CI*(T3(18))))))))
     $ +P1(3)*(P1(3)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(0)*(T3(7)
     $ +T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(-1D0)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P2(1)*(-1D0)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))
     $ -P2(2)*(T3(9)+T3(12)+2D0 * CI*(T3(13))))))))))+F2(3)*(P1(0)
     $ *(P1(1)*(T3(16)+T3(10)+CI*(T3(11)+T3(5))-2D0*(T3(4)+T3(7)))
     $ +(P1(2)*(-1D0)*(+2D0*(T3(5)+T3(11))+CI*(T3(7)+T3(4))-T3(17)
     $ -T3(14))+(P1(3)*2D0*(T3(18)+T3(3)-T3(6)-T3(15))+(P2(1)*(T3(4)
     $ +T3(7)-T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(P1(0)
     $ *(-1D0)*(T3(15)+T3(6)-2D0*(T3(3)))+(P2(0)*(T3(15)+T3(6)-2D0
     $ *(T3(3)))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18))))))))))+(P1(3)*(P1(1)
     $ *(-1D0)*(T3(4)+T3(7)+CI*(T3(14)+T3(17))-2D0*(T3(10)+T3(16)))
     $ +(P1(2)*(+2D0*(T3(14)+T3(17))+CI*(T3(10)+T3(16))-T3(5)-T3(11))
     $ +(P2(1)*(T3(4)+T3(7)-T3(16)-T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)
     $ -T3(14))+(P1(3)*(-1D0)*(T3(6)+T3(15)-2D0*(T3(18)))+(P2(0)
     $ *(T3(15)+T3(6)-2D0*(T3(3)))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18))))))
     $ )))+(P1(1)*(P1(2)*2D0*(T3(9)+T3(12)-CI*(T3(13))+CI*(T3(8)))
     $ +(P2(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P2(3)*(+CI*(T3(14)
     $ +T3(17))-T3(10)-T3(16))+(P1(1)*(-1D0)*(+CI*(T3(12)+T3(9))-2D0
     $ *(T3(8)))+(P2(1)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))-P2(2)*(T3(9)
     $ +T3(12)-2D0 * CI*(T3(13))))))))+P1(2)*(P2(0)*(T3(11)+T3(5)+CI
     $ *(T3(7)+T3(4)))+(P2(3)*(-1D0)*(T3(14)+T3(17)+CI*(T3(10)+T3(16)))
     $ +(P1(2)*(+2D0*(T3(13))+CI*(T3(9)+T3(12)))+(P2(1)*(-1D0)*(T3(12)
     $ +T3(9)+2D0 * CI*(T3(8)))-P2(2)*(+2D0*(T3(13))+CI*(T3(9)+T3(12)))
     $ ))))))))
      F1(4)= DENOM*CI*(F2(3)*(P1(0)*(P1(3)*(+CI*(T3(14)+T3(11)+T3(17)
     $ +T3(5))-T3(10)-T3(7)-T3(16)-T3(4))+(P1(1)*(T3(15)+T3(6)+CI
     $ *(T3(12)+T3(9))-2D0*(T3(8)+T3(3)))+(P1(2)*(-1D0)*(T3(9)+T3(12)
     $ -2D0 * CI*(T3(13)+T3(3))+CI*(T3(15)+T3(6)))+(P1(0)*(T3(7)+T3(4)
     $ -CI*(T3(11)+T3(5)))+(P2(0)*(+CI*(T3(11)+T3(5))-T3(7)-T3(4))
     $ +(P2(3)*(T3(10)+T3(16)-CI*(T3(14)+T3(17)))+(P2(1)*(-1D0)*(+CI
     $ *(T3(12)+T3(9))-2D0*(T3(8)))+P2(2)*(T3(9)+T3(12)-2D0 * CI
     $ *(T3(13))))))))))+(P1(1)*(P1(2)*(T3(5)+T3(11)-CI*(T3(4)+T3(7))
     $ +CI*(T3(16)+T3(10))-T3(17)-T3(14))+(P1(3)*(T3(6)+T3(15)+2D0
     $ *(T3(8))-CI*(T3(12)+T3(9))-2D0*(T3(18)))+(P1(1)*(T3(4)+T3(7)
     $ -T3(16)-T3(10))+(P2(1)*(T3(16)+T3(10)-T3(4)-T3(7))+(P2(2)
     $ *(T3(17)+T3(14)-T3(5)-T3(11))+(P2(0)*(-1D0)*(T3(15)+T3(6)-2D0
     $ *(T3(3)))-P2(3)*(T3(6)+T3(15)-2D0*(T3(18)))))))))+(P1(2)*(P1(3)
     $ *(T3(9)+T3(12)-CI*(T3(6)+T3(15))-2D0 * CI*(T3(13))+2D0 * CI
     $ *(T3(18)))+(P1(2)*(-CI*(T3(5)+T3(11))+CI*(T3(17)+T3(14)))+(P2(1)
     $ *(-CI*(T3(16)+T3(10))+CI*(T3(4)+T3(7)))+(P2(2)*(-CI*(T3(17)
     $ +T3(14))+CI*(T3(5)+T3(11)))+(P2(0)*(-2D0 * CI*(T3(3))+CI*(T3(15)
     $ +T3(6)))+P2(3)*(-2D0 * CI*(T3(18))+CI*(T3(6)+T3(15))))))))+P1(3)
     $ *(P1(3)*(T3(10)+T3(16)-CI*(T3(14)+T3(17)))+(P2(0)*(T3(7)+T3(4)
     $ -CI*(T3(11)+T3(5)))+(P2(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))
     $ +(P2(1)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))-P2(2)*(T3(9)+T3(12)
     $ -2D0 * CI*(T3(13))))))))))+F2(4)*(P1(0)*(P1(1)*(T3(16)+T3(10)
     $ +2D0*(T3(4)+T3(7))+CI*(T3(11)+T3(5)))+(P1(2)*(T3(17)+T3(14)+2D0
     $ *(T3(5)+T3(11))-CI*(T3(7)+T3(4)))+(P1(3)*2D0*(T3(6)+T3(18)+T3(3)
     $ +T3(15))+(P2(1)*(-1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(-1D0)
     $ *(T3(5)+T3(17)+T3(11)+T3(14))+(P1(0)*(-1D0)*(T3(15)+T3(6)+2D0
     $ *(T3(3)))+(P2(0)*(T3(15)+T3(6)+2D0*(T3(3)))-P2(3)*(T3(6)+T3(15)
     $ +2D0*(T3(18))))))))))+(P1(3)*(P1(1)*(-1D0)*(T3(4)+T3(7)+2D0
     $ *(T3(10)+T3(16))+CI*(T3(14)+T3(17)))+(P1(2)*(-1D0)*(T3(5)+T3(11)
     $ -CI*(T3(10)+T3(16))+2D0*(T3(14)+T3(17)))+(P2(1)*(T3(4)+T3(16)
     $ +T3(7)+T3(10))+(P2(2)*(T3(5)+T3(17)+T3(11)+T3(14))+(P1(3)*(-1D0)
     $ *(T3(6)+T3(15)+2D0*(T3(18)))+(P2(0)*(-1D0)*(T3(15)+T3(6)+2D0
     $ *(T3(3)))+P2(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P1(1)*(P1(2)
     $ *(-2D0)*(T3(9)+T3(12)-CI*(T3(8))+CI*(T3(13)))+(P2(0)*(-1D0)
     $ *(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P1(1)*(-1D0)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))
     $ +(P2(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))+P2(2)*(T3(9)+T3(12)
     $ +2D0 * CI*(T3(13))))))))+P1(2)*(P2(0)*(+CI*(T3(7)+T3(4))-T3(11)
     $ -T3(5))+(P2(3)*(T3(14)+T3(17)-CI*(T3(10)+T3(16)))+(P1(2)*(+CI
     $ *(T3(9)+T3(12))-2D0*(T3(13)))+(P2(1)*(T3(12)+T3(9)-2D0 * CI
     $ *(T3(8)))-P2(2)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))))))))))
      F1(5)= DENOM*CI * M1*(F2(3)*(P1(1)*(T3(16)+T3(10)-T3(4)-T3(7))
     $ +(P1(2)*(T3(17)+T3(14)-T3(5)-T3(11))+(P2(1)*(T3(4)+T3(7)-T3(16)
     $ -T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(T3(15)*(P2(0)
     $ +P2(3)-P1(0)-P1(3))+(T3(6)*(P2(3)+P2(0)-P1(3)-P1(0))+(T3(3)*2D0
     $ *(P1(0)-P2(0))+2D0*(T3(18)*(P1(3)-P2(3))))))))))+F2(4)*(P1(0)*(
     $ -1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)*(T3(10)+T3(16)+CI
     $ *(T3(14)+T3(17)))+(P2(0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(3)
     $ *(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(T3(12)*(P1(2)-CI
     $ *(P2(1))+CI*(P1(1))-P2(2))+(T3(9)*(P1(2)-CI*(P2(1))+CI*(P1(1))
     $ -P2(2))+(T3(13)*2D0*(-CI*(P2(2))+CI*(P1(2)))+2D0*(T3(8)*(P1(1)
     $ -P2(1)))))))))))
      F1(6)= DENOM*(-CI )* M1*(F2(3)*(P1(0)*(T3(7)+T3(4)-CI*(T3(11)
     $ +T3(5)))+(P1(3)*(+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P2(0)*(+CI
     $ *(T3(11)+T3(5))-T3(7)-T3(4))+(P2(3)*(T3(10)+T3(16)-CI*(T3(14)
     $ +T3(17)))+(T3(12)*(P2(2)-CI*(P2(1))+CI*(P1(1))-P1(2))+(T3(9)
     $ *(P2(2)-CI*(P2(1))+CI*(P1(1))-P1(2))+(T3(13)*2D0*(-CI*(P2(2))
     $ +CI*(P1(2)))+2D0*(T3(8)*(P2(1)-P1(1))))))))))+F2(4)*(P1(1)
     $ *(T3(4)+T3(16)+T3(7)+T3(10))+(P1(2)*(T3(5)+T3(17)+T3(11)+T3(14))
     $ +(P2(1)*(-1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(-1D0)*(T3(5)
     $ +T3(17)+T3(11)+T3(14))+(T3(15)*(P2(0)+P1(3)-P1(0)-P2(3))+(T3(6)
     $ *(P1(3)+P2(0)-P2(3)-P1(0))+(T3(3)*2D0*(P2(0)-P1(0))+2D0*(T3(18)
     $ *(P1(3)-P2(3)))))))))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjP(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjP(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjP(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFT5_2(F1, T3, COUP, M2, W2,F2)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(6)
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      REAL*8 W2
      COMPLEX*16 F1(*)
      REAL*8 M2
      COMPLEX*16 DENOM
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      F2(1) = +F1(1)+T3(1)
      F2(2) = +F1(2)+T3(2)
      P2(0) = -DBLE(F2(1))
      P2(1) = -DBLE(F2(2))
      P2(2) = -DIMAG(F2(2))
      P2(3) = -DIMAG(F2(1))
      DENOM = COUP/(P2(0)**2-P2(1)**2-P2(2)**2-P2(3)**2 - M2 * (M2 -CI
     $ * W2))
      F2(3)= DENOM*CI * M2*(F1(5)*(P1(1)*(T3(16)+T3(10)-T3(4)-T3(7))
     $ +(P1(2)*(T3(17)+T3(14)-T3(5)-T3(11))+(P2(1)*(T3(4)+T3(7)-T3(16)
     $ -T3(10))+(P2(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(T3(15)*(P2(0)
     $ +P2(3)-P1(0)-P1(3))+(T3(6)*(P2(3)+P2(0)-P1(3)-P1(0))+(T3(3)*2D0
     $ *(P1(0)-P2(0))+2D0*(T3(18)*(P1(3)-P2(3))))))))))+F1(6)*(P1(0)*(
     $ +CI*(T3(11)+T3(5))-T3(7)-T3(4))+(P1(3)*(T3(10)+T3(16)-CI*(T3(14)
     $ +T3(17)))+(P2(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P2(3)*(+CI
     $ *(T3(14)+T3(17))-T3(10)-T3(16))+(T3(12)*(P1(2)-CI*(P1(1))+CI
     $ *(P2(1))-P2(2))+(T3(9)*(P1(2)-CI*(P1(1))+CI*(P2(1))-P2(2))
     $ +(T3(13)*2D0*(-CI*(P1(2))+CI*(P2(2)))+2D0*(T3(8)*(P1(1)-P2(1))))
     $ )))))))
      F2(4)= DENOM*(-CI )* M2*(F1(5)*(P1(0)*(T3(7)+T3(4)+CI*(T3(11)
     $ +T3(5)))+(P1(3)*(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(0)
     $ *(-1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P2(3)*(T3(10)+T3(16)
     $ +CI*(T3(14)+T3(17)))+(T3(12)*(P2(2)-CI*(P1(1))+CI*(P2(1))-P1(2))
     $ +(T3(9)*(P2(2)-CI*(P1(1))+CI*(P2(1))-P1(2))+(T3(13)*2D0*(-CI
     $ *(P1(2))+CI*(P2(2)))+2D0*(T3(8)*(P2(1)-P1(1))))))))))+F1(6)
     $ *(P1(1)*(T3(4)+T3(16)+T3(7)+T3(10))+(P1(2)*(T3(5)+T3(17)+T3(11)
     $ +T3(14))+(P2(1)*(-1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P2(2)*(-1D0)
     $ *(T3(5)+T3(17)+T3(11)+T3(14))+(T3(15)*(P2(0)+P1(3)-P1(0)-P2(3))
     $ +(T3(6)*(P1(3)+P2(0)-P2(3)-P1(0))+(T3(3)*2D0*(P2(0)-P1(0))+2D0
     $ *(T3(18)*(P1(3)-P2(3)))))))))))
      F2(5)= DENOM*CI*(F1(6)*(P2(0)*(P2(3)*(T3(7)+T3(4)-CI*(T3(11)
     $ +T3(5))+CI*(T3(14)+T3(17))-T3(10)-T3(16))+(P2(1)*(-1D0)*(T3(15)
     $ +T3(6)+2D0*(T3(8)+T3(3))-CI*(T3(12)+T3(9)))+(P2(2)*(+CI*(T3(15)
     $ +T3(6))+2D0 * CI*(T3(13)+T3(3))-T3(9)-T3(12))+(P1(0)*(+CI
     $ *(T3(11)+T3(5))-T3(7)-T3(4))+(P1(3)*(T3(10)+T3(16)-CI*(T3(14)
     $ +T3(17)))+(P2(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P1(1)*(-1D0)
     $ *(+CI*(T3(12)+T3(9))-2D0*(T3(8)))+P1(2)*(T3(9)+T3(12)-2D0 * CI
     $ *(T3(13))))))))))+(P2(1)*(P2(2)*(T3(5)+T3(17)+T3(11)+T3(14)-CI
     $ *(T3(4)+T3(16)+T3(7)+T3(10)))+(P2(3)*(T3(6)+T3(15)+2D0*(T3(18))
     $ +CI*(T3(12)+T3(9))-2D0*(T3(8)))+(P1(1)*(-1D0)*(T3(4)+T3(16)
     $ +T3(7)+T3(10))+(P1(2)*(-1D0)*(T3(5)+T3(17)+T3(11)+T3(14))+(P2(1)
     $ *(T3(4)+T3(16)+T3(7)+T3(10))+(P1(0)*(T3(15)+T3(6)+2D0*(T3(3)))
     $ -P1(3)*(T3(6)+T3(15)+2D0*(T3(18)))))))))+(P2(2)*(P2(3)*(-1D0)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(13))+CI*(T3(6)+T3(15))+2D0 * CI
     $ *(T3(18)))+(P1(1)*(+CI*(T3(4)+T3(16)+T3(7)+T3(10)))+(P1(2)*(+CI
     $ *(T3(5)+T3(17)+T3(11)+T3(14)))+(P2(2)*(-1D0)*(+CI*(T3(5)+T3(17)
     $ +T3(11)+T3(14)))+(P1(0)*(-1D0)*(+CI*(T3(15)+T3(6))+2D0 * CI
     $ *(T3(3)))+P1(3)*(+CI*(T3(6)+T3(15))+2D0 * CI*(T3(18))))))))
     $ +P2(3)*(P1(0)*(+CI*(T3(11)+T3(5))-T3(7)-T3(4))+(P1(3)*(T3(10)
     $ +T3(16)-CI*(T3(14)+T3(17)))+(P2(3)*(+CI*(T3(14)+T3(17))-T3(10)
     $ -T3(16))+(P1(1)*(-1D0)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))+P1(2)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(13))))))))))+F1(5)*(P2(0)*(P2(1)*(
     $ +2D0*(T3(4)+T3(7))+CI*(T3(11)+T3(5))-T3(16)-T3(10))+(P2(2)*(
     $ -1D0)*(T3(17)+T3(14)+CI*(T3(7)+T3(4))-2D0*(T3(5)+T3(11)))+(P1(1)
     $ *(T3(16)+T3(10)-T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)-T3(11))
     $ +(P2(3)*2D0*(T3(6)+T3(15)-T3(18)-T3(3))+(P1(0)*(-1D0)*(T3(15)
     $ +T3(6)-2D0*(T3(3)))+(P1(3)*(-1D0)*(T3(6)+T3(15)-2D0*(T3(18)))
     $ +P2(0)*(T3(15)+T3(6)-2D0*(T3(3))))))))))+(P2(3)*(P2(1)*(-1D0)*(
     $ +2D0*(T3(10)+T3(16))+CI*(T3(14)+T3(17))-T3(4)-T3(7))+(P2(2)
     $ *(T3(5)+T3(11)+CI*(T3(10)+T3(16))-2D0*(T3(14)+T3(17)))+(P1(1)
     $ *(T3(16)+T3(10)-T3(4)-T3(7))+(P1(2)*(T3(17)+T3(14)-T3(5)-T3(11))
     $ +(P1(0)*(-1D0)*(T3(15)+T3(6)-2D0*(T3(3)))+(P1(3)*(-1D0)*(T3(6)
     $ +T3(15)-2D0*(T3(18)))+P2(3)*(T3(6)+T3(15)-2D0*(T3(18)))))))))
     $ +(P2(1)*(P1(0)*(-1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)
     $ *(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(2)*(-2D0)*(T3(9)+T3(12)
     $ -CI*(T3(8))+CI*(T3(13)))+(P1(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))
     $ +(P1(2)*(T3(9)+T3(12)+2D0 * CI*(T3(13)))-P2(1)*(+2D0*(T3(8))+CI
     $ *(T3(12)+T3(9))))))))+P2(2)*(P1(0)*(+CI*(T3(7)+T3(4))-T3(11)
     $ -T3(5))+(P1(3)*(T3(14)+T3(17)-CI*(T3(10)+T3(16)))+(P1(1)*(T3(12)
     $ +T3(9)-2D0 * CI*(T3(8)))+(P1(2)*(-1D0)*(+CI*(T3(9)+T3(12))-2D0
     $ *(T3(13)))+P2(2)*(+CI*(T3(9)+T3(12))-2D0*(T3(13)))))))))))
      F2(6)= DENOM*(-CI)*(F1(5)*(P2(0)*(P2(3)*(T3(10)+T3(7)+T3(16)
     $ +T3(4)+CI*(T3(14)+T3(11)+T3(17)+T3(5)))+(P2(1)*(+2D0*(T3(8)
     $ +T3(3))+CI*(T3(12)+T3(9))-T3(15)-T3(6))+(P2(2)*(T3(9)+T3(12)-CI
     $ *(T3(15)+T3(6))+2D0 * CI*(T3(13)+T3(3)))+(P1(0)*(T3(7)+T3(4)+CI
     $ *(T3(11)+T3(5)))+(P1(3)*(-1D0)*(T3(10)+T3(16)+CI*(T3(14)+T3(17))
     $ )+(P2(0)*(-1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(1)*(-1D0)*(
     $ +2D0*(T3(8))+CI*(T3(12)+T3(9)))-P1(2)*(T3(9)+T3(12)+2D0 * CI
     $ *(T3(13))))))))))+(P2(1)*(P2(2)*(T3(17)+T3(14)-CI*(T3(4)+T3(7))
     $ +CI*(T3(16)+T3(10))-T3(5)-T3(11))+(P2(3)*(-1D0)*(T3(6)+T3(15)
     $ +2D0*(T3(8))+CI*(T3(12)+T3(9))-2D0*(T3(18)))+(P1(1)*(T3(4)+T3(7)
     $ -T3(16)-T3(10))+(P1(2)*(T3(5)+T3(11)-T3(17)-T3(14))+(P2(1)
     $ *(T3(16)+T3(10)-T3(4)-T3(7))+(P1(0)*(T3(15)+T3(6)-2D0*(T3(3)))
     $ +P1(3)*(T3(6)+T3(15)-2D0*(T3(18)))))))))+(P2(2)*(P2(3)*(-1D0)
     $ *(T3(9)+T3(12)-2D0 * CI*(T3(18))+CI*(T3(6)+T3(15))+2D0 * CI
     $ *(T3(13)))+(P1(1)*(-CI*(T3(16)+T3(10))+CI*(T3(4)+T3(7)))+(P1(2)
     $ *(-CI*(T3(17)+T3(14))+CI*(T3(5)+T3(11)))+(P2(2)*(-CI*(T3(5)
     $ +T3(11))+CI*(T3(17)+T3(14)))+(P1(0)*(-2D0 * CI*(T3(3))+CI
     $ *(T3(15)+T3(6)))+P1(3)*(-2D0 * CI*(T3(18))+CI*(T3(6)+T3(15))))))
     $ ))+P2(3)*(P1(0)*(-1D0)*(T3(7)+T3(4)+CI*(T3(11)+T3(5)))+(P1(3)
     $ *(T3(10)+T3(16)+CI*(T3(14)+T3(17)))+(P2(3)*(-1D0)*(T3(10)+T3(16)
     $ +CI*(T3(14)+T3(17)))+(P1(1)*(+2D0*(T3(8))+CI*(T3(12)+T3(9)))
     $ +P1(2)*(T3(9)+T3(12)+2D0 * CI*(T3(13))))))))))+F1(6)*(P2(0)
     $ *(P2(1)*(-1D0)*(T3(16)+T3(10)+2D0*(T3(4)+T3(7))-CI*(T3(11)+T3(5)
     $ ))+(P2(2)*(-1D0)*(T3(17)+T3(14)+2D0*(T3(5)+T3(11))+CI*(T3(7)
     $ +T3(4)))+(P1(1)*(T3(4)+T3(16)+T3(7)+T3(10))+(P1(2)*(T3(5)+T3(17)
     $ +T3(11)+T3(14))+(P2(3)*(-2D0)*(T3(6)+T3(18)+T3(3)+T3(15))+(P1(0)
     $ *(-1D0)*(T3(15)+T3(6)+2D0*(T3(3)))+(P1(3)*(T3(6)+T3(15)+2D0
     $ *(T3(18)))+P2(0)*(T3(15)+T3(6)+2D0*(T3(3))))))))))+(P2(3)*(P2(1)
     $ *(T3(4)+T3(7)+2D0*(T3(10)+T3(16))-CI*(T3(14)+T3(17)))+(P2(2)
     $ *(T3(5)+T3(11)+2D0*(T3(14)+T3(17))+CI*(T3(10)+T3(16)))+(P1(1)*(
     $ -1D0)*(T3(4)+T3(16)+T3(7)+T3(10))+(P1(2)*(-1D0)*(T3(5)+T3(17)
     $ +T3(11)+T3(14))+(P1(0)*(T3(15)+T3(6)+2D0*(T3(3)))+(P1(3)*(-1D0)
     $ *(T3(6)+T3(15)+2D0*(T3(18)))+P2(3)*(T3(6)+T3(15)+2D0*(T3(18)))))
     $ ))))+(P2(1)*(P1(0)*(T3(7)+T3(4)-CI*(T3(11)+T3(5)))+(P1(3)*(+CI
     $ *(T3(14)+T3(17))-T3(10)-T3(16))+(P2(2)*2D0*(T3(9)+T3(12)-CI
     $ *(T3(13))+CI*(T3(8)))+(P1(1)*(+CI*(T3(12)+T3(9))-2D0*(T3(8)))
     $ +(P1(2)*(-1D0)*(T3(9)+T3(12)-2D0 * CI*(T3(13)))-P2(1)*(+CI
     $ *(T3(12)+T3(9))-2D0*(T3(8))))))))+P2(2)*(P1(0)*(T3(11)+T3(5)+CI
     $ *(T3(7)+T3(4)))+(P1(3)*(-1D0)*(T3(14)+T3(17)+CI*(T3(10)+T3(16)))
     $ +(P1(1)*(-1D0)*(T3(12)+T3(9)+2D0 * CI*(T3(8)))+(P1(2)*(-1D0)*(
     $ +2D0*(T3(13))+CI*(T3(9)+T3(12)))+P2(2)*(+2D0*(T3(13))+CI*(T3(9)
     $ +T3(12)))))))))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjP(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjP(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjP(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFT5_3(F1, F2, COUP, M3, W3,T3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 DENOM
      REAL*8 P1(0:3)
      REAL*8 W3
      REAL*8 P2(0:3)
      REAL*8 M3
      COMPLEX*16 TMP17
      REAL*8 P3(0:3)
      COMPLEX*16 TMP20
      COMPLEX*16 F1(*)
      COMPLEX*16 TMP21
      COMPLEX*16 F2(*)
      REAL*8 OM3
      COMPLEX*16 T3(18)
      COMPLEX*16 COUP
      COMPLEX*16 TMP19
      COMPLEX*16 TMP18
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      T3(1) = +F1(1)+F2(1)
      T3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(T3(1))
      P3(1) = -DBLE(T3(2))
      P3(2) = -DIMAG(T3(2))
      P3(3) = -DIMAG(T3(1))
      TMP19 = (F1(5)*(F2(3)*(P3(0)-P3(3))-F2(4)*(P3(1)+CI*(P3(2))))
     $ +F1(6)*(F2(3)*(+CI*(P3(2))-P3(1))+F2(4)*(P3(0)+P3(3))))
      TMP17 = (F1(5)*(F2(3)*(P1(0)-P1(3))-F2(4)*(P1(1)+CI*(P1(2))))
     $ +F1(6)*(F2(3)*(+CI*(P1(2))-P1(1))+F2(4)*(P1(0)+P1(3))))
      TMP20 = (F1(5)*(F2(3)*(P2(0)-P2(3))-F2(4)*(P2(1)+CI*(P2(2))))
     $ +F1(6)*(F2(3)*(+CI*(P2(2))-P2(1))+F2(4)*(P2(0)+P2(3))))
      TMP18 = (P1(0)*P3(0)-P1(1)*P3(1)-P1(2)*P3(2)-P1(3)*P3(3))
      TMP21 = (P2(0)*P3(0)-P2(1)*P3(1)-P2(2)*P3(2)-P2(3)*P3(3))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      T3(3)= DENOM*2D0 * CI*(OM3*(P3(0)*(P3(0)*(OM3*2D0/3D0 * TMP19
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17)))+(F1(5)*F2(3)
     $ *(TMP21-TMP18)+(F1(6)*F2(4)*(TMP21-TMP18)+TMP19*(P2(0)-P1(0)))))
     $ +1D0/3D0*(TMP19*(TMP18-TMP21)))+(F1(5)*F2(3)*(P1(0)-P2(0))
     $ +(F1(6)*F2(4)*(P1(0)-P2(0))+(-1D0/3D0*(TMP17)+1D0/3D0*(TMP20))))
     $ )
      T3(4)= DENOM*CI*(OM3*(P3(0)*(P3(1)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(TMP21
     $ -TMP18)+(F1(6)*F2(3)*(TMP21-TMP18)+TMP19*(P2(1)-P1(1)))))+P3(1)
     $ *(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)*F2(4)*(TMP21-TMP18)+TMP19
     $ *(P2(0)-P1(0)))))+(F1(5)*(F2(3)*(P1(1)-P2(1))+F2(4)*(P1(0)-P2(0)
     $ ))+F1(6)*(F2(3)*(P1(0)-P2(0))+F2(4)*(P1(1)-P2(1)))))
      T3(5)= DENOM*CI*(OM3*(P3(0)*(P3(2)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(-CI
     $ *(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)+CI*(TMP18))
     $ +TMP19*(P2(2)-P1(2)))))+P3(2)*(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)
     $ *F2(4)*(TMP21-TMP18)+TMP19*(P2(0)-P1(0)))))+(F1(5)*(F2(3)*(P1(2)
     $ -P2(2))+F2(4)*(-CI*(P2(0))+CI*(P1(0))))+F1(6)*(F2(3)*(-CI*(P1(0)
     $ )+CI*(P2(0)))+F2(4)*(P1(2)-P2(2)))))
      T3(6)= DENOM*CI*(OM3*(P3(0)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)*F2(4)*(TMP21-TMP18)+TMP19
     $ *(P2(0)-P1(0)))))+(F1(5)*F2(3)*(P1(3)+P1(0)-P2(3)-P2(0))+F1(6)
     $ *F2(4)*(P1(3)+P2(0)-P1(0)-P2(3))))
      T3(7)= DENOM*CI*(OM3*(P3(0)*(P3(1)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(TMP21
     $ -TMP18)+(F1(6)*F2(3)*(TMP21-TMP18)+TMP19*(P2(1)-P1(1)))))+P3(1)
     $ *(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)*F2(4)*(TMP21-TMP18)+TMP19
     $ *(P2(0)-P1(0)))))+(F1(5)*(F2(3)*(P1(1)-P2(1))+F2(4)*(P1(0)-P2(0)
     $ ))+F1(6)*(F2(3)*(P1(0)-P2(0))+F2(4)*(P1(1)-P2(1)))))
      T3(8)= DENOM*2D0 * CI*(OM3*(P3(1)*(P3(1)*(OM3*2D0/3D0 * TMP19
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17)))+(F1(5)*F2(4)
     $ *(TMP21-TMP18)+(F1(6)*F2(3)*(TMP21-TMP18)+TMP19*(P2(1)-P1(1)))))
     $ +1D0/3D0*(TMP19*(TMP21-TMP18)))+(F1(5)*F2(4)*(P1(1)-P2(1))
     $ +(F1(6)*F2(3)*(P1(1)-P2(1))+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17))))
     $ )
      T3(9)= DENOM*CI*(OM3*(P3(1)*(P3(2)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(-CI
     $ *(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)+CI*(TMP18))
     $ +TMP19*(P2(2)-P1(2)))))+P3(2)*(F1(5)*F2(4)*(TMP21-TMP18)+(F1(6)
     $ *F2(3)*(TMP21-TMP18)+TMP19*(P2(1)-P1(1)))))+(F1(5)*F2(4)*(P1(2)
     $ -CI*(P2(1))+CI*(P1(1))-P2(2))+F1(6)*F2(3)*(P1(2)-CI*(P1(1))+CI
     $ *(P2(1))-P2(2))))
      T3(10)= DENOM*CI*(OM3*(P3(1)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(4)*(TMP21-TMP18)+(F1(6)*F2(3)*(TMP21-TMP18)+TMP19
     $ *(P2(1)-P1(1)))))+(F1(5)*(F2(3)*(P1(1)-P2(1))+F2(4)*(P1(3)-P2(3)
     $ ))+F1(6)*(F2(3)*(P1(3)-P2(3))+F2(4)*(P2(1)-P1(1)))))
      T3(11)= DENOM*CI*(OM3*(P3(0)*(P3(2)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(-CI
     $ *(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)+CI*(TMP18))
     $ +TMP19*(P2(2)-P1(2)))))+P3(2)*(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)
     $ *F2(4)*(TMP21-TMP18)+TMP19*(P2(0)-P1(0)))))+(F1(5)*(F2(3)*(P1(2)
     $ -P2(2))+F2(4)*(-CI*(P2(0))+CI*(P1(0))))+F1(6)*(F2(3)*(-CI*(P1(0)
     $ )+CI*(P2(0)))+F2(4)*(P1(2)-P2(2)))))
      T3(12)= DENOM*CI*(OM3*(P3(1)*(P3(2)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(4)*(-CI
     $ *(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)+CI*(TMP18))
     $ +TMP19*(P2(2)-P1(2)))))+P3(2)*(F1(5)*F2(4)*(TMP21-TMP18)+(F1(6)
     $ *F2(3)*(TMP21-TMP18)+TMP19*(P2(1)-P1(1)))))+(F1(5)*F2(4)*(P1(2)
     $ -CI*(P2(1))+CI*(P1(1))-P2(2))+F1(6)*F2(3)*(P1(2)-CI*(P1(1))+CI
     $ *(P2(1))-P2(2))))
      T3(13)= DENOM*2D0 * CI*(OM3*(P3(2)*(P3(2)*(OM3*2D0/3D0 * TMP19
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17)))+(F1(5)*F2(4)
     $ *(-CI*(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)+CI*(TMP18))
     $ +TMP19*(P2(2)-P1(2)))))+1D0/3D0*(TMP19*(TMP21-TMP18)))+(F1(5)
     $ *F2(4)*(-CI*(P2(2))+CI*(P1(2)))+(F1(6)*F2(3)*(-CI*(P1(2))+CI
     $ *(P2(2)))+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17)))))
      T3(14)= DENOM*CI*(OM3*(P3(2)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(4)*(-CI*(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)
     $ +CI*(TMP18))+TMP19*(P2(2)-P1(2)))))+(F1(5)*(F2(3)*(P1(2)-P2(2))
     $ +F2(4)*(-CI*(P2(3))+CI*(P1(3))))+F1(6)*(F2(3)*(-CI*(P1(3))+CI
     $ *(P2(3)))+F2(4)*(P2(2)-P1(2)))))
      T3(15)= DENOM*CI*(OM3*(P3(0)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(3)*(TMP21-TMP18)+(F1(6)*F2(4)*(TMP21-TMP18)+TMP19
     $ *(P2(0)-P1(0)))))+(F1(5)*F2(3)*(P1(0)+P1(3)-P2(0)-P2(3))+F1(6)
     $ *F2(4)*(P1(3)+P2(0)-P1(0)-P2(3))))
      T3(16)= DENOM*CI*(OM3*(P3(1)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(4)*(TMP21-TMP18)+(F1(6)*F2(3)*(TMP21-TMP18)+TMP19
     $ *(P2(1)-P1(1)))))+(F1(5)*(F2(3)*(P1(1)-P2(1))+F2(4)*(P1(3)-P2(3)
     $ ))+F1(6)*(F2(3)*(P1(3)-P2(3))+F2(4)*(P2(1)-P1(1)))))
      T3(17)= DENOM*CI*(OM3*(P3(2)*(P3(3)*(OM3*4D0/3D0 * TMP19*(TMP18
     $ -TMP21)+(-2D0/3D0*(TMP20)+2D0/3D0*(TMP17)))+(F1(5)*F2(3)*(TMP21
     $ -TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))+P3(3)
     $ *(F1(5)*F2(4)*(-CI*(TMP18)+CI*(TMP21))+(F1(6)*F2(3)*(-CI*(TMP21)
     $ +CI*(TMP18))+TMP19*(P2(2)-P1(2)))))+(F1(5)*(F2(3)*(P1(2)-P2(2))
     $ +F2(4)*(-CI*(P2(3))+CI*(P1(3))))+F1(6)*(F2(3)*(-CI*(P1(3))+CI
     $ *(P2(3)))+F2(4)*(P2(2)-P1(2)))))
      T3(18)= DENOM*2D0 * CI*(OM3*(P3(3)*(P3(3)*(OM3*2D0/3D0 * TMP19
     $ *(TMP18-TMP21)+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17)))+(F1(5)*F2(3)
     $ *(TMP21-TMP18)+(F1(6)*F2(4)*(TMP18-TMP21)+TMP19*(P2(3)-P1(3)))))
     $ +1D0/3D0*(TMP19*(TMP21-TMP18)))+(F1(5)*F2(3)*(P1(3)-P2(3))
     $ +(F1(6)*F2(4)*(P2(3)-P1(3))+(-1D0/3D0*(TMP20)+1D0/3D0*(TMP17))))
     $ )
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,1)
C     
      SUBROUTINE FFV4_0_S(F1, F2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 TMP6
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      TMP6 = (F1(3)*(F2(5)*(V3(3)+V3(6))+F2(6)*(V3(4)+CI*(V3(5))))
     $ +(F1(4)*(F2(5)*(V3(4)-CI*(V3(5)))+F2(6)*(V3(3)-V3(6)))+(F1(5)
     $ *(F2(3)*(V3(3)-V3(6))-F2(4)*(V3(4)+CI*(V3(5))))+F1(6)*(F2(3)*(
     $ +CI*(V3(5))-V3(4))+F2(4)*(V3(3)+V3(6))))))
      VERTEX = COUP*(-CI * TMP6)
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,1)
C     
      SUBROUTINE FFV4_1_S(F2, V3, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      COMPLEX*16 V3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 W1
      COMPLEX*16 F1(6)
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      F1(1) = +F2(1)+V3(1)
      F1(2) = +F2(2)+V3(2)
      P1(0) = -DBLE(F1(1))
      P1(1) = -DBLE(F1(2))
      P1(2) = -DIMAG(F1(2))
      P1(3) = -DIMAG(F1(1))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(3)= DENOM*CI*(F2(3)*(P1(0)*(V3(6)-V3(3))+(P1(1)*(V3(4)-CI
     $ *(V3(5)))+(P1(2)*(V3(5)+CI*(V3(4)))+P1(3)*(V3(6)-V3(3)))))
     $ +(F2(4)*(P1(0)*(V3(4)+CI*(V3(5)))+(P1(1)*(-1D0)*(V3(3)+V3(6))
     $ +(P1(2)*(-1D0)*(+CI*(V3(3)+V3(6)))+P1(3)*(V3(4)+CI*(V3(5))))))
     $ +M1*(F2(5)*(V3(3)+V3(6))+F2(6)*(V3(4)+CI*(V3(5))))))
      F1(4)= DENOM*(-CI)*(F2(3)*(P1(0)*(+CI*(V3(5))-V3(4))+(P1(1)
     $ *(V3(3)-V3(6))+(P1(2)*(-CI*(V3(3))+CI*(V3(6)))+P1(3)*(V3(4)-CI
     $ *(V3(5))))))+(F2(4)*(P1(0)*(V3(3)+V3(6))+(P1(1)*(-1D0)*(V3(4)
     $ +CI*(V3(5)))+(P1(2)*(+CI*(V3(4))-V3(5))-P1(3)*(V3(3)+V3(6)))))
     $ +M1*(F2(5)*(+CI*(V3(5))-V3(4))+F2(6)*(V3(6)-V3(3)))))
      F1(5)= DENOM*(-CI)*(F2(5)*(P1(0)*(V3(3)+V3(6))+(P1(1)*(+CI*(V3(5)
     $ )-V3(4))+(P1(2)*(-1D0)*(V3(5)+CI*(V3(4)))-P1(3)*(V3(3)+V3(6)))))
     $ +(F2(6)*(P1(0)*(V3(4)+CI*(V3(5)))+(P1(1)*(V3(6)-V3(3))+(P1(2)*(
     $ -CI*(V3(3))+CI*(V3(6)))-P1(3)*(V3(4)+CI*(V3(5))))))+M1*(F2(3)
     $ *(V3(6)-V3(3))+F2(4)*(V3(4)+CI*(V3(5))))))
      F1(6)= DENOM*CI*(F2(5)*(P1(0)*(+CI*(V3(5))-V3(4))+(P1(1)*(V3(3)
     $ +V3(6))+(P1(2)*(-1D0)*(+CI*(V3(3)+V3(6)))+P1(3)*(+CI*(V3(5))
     $ -V3(4)))))+(F2(6)*(P1(0)*(V3(6)-V3(3))+(P1(1)*(V3(4)+CI*(V3(5)))
     $ +(P1(2)*(V3(5)-CI*(V3(4)))+P1(3)*(V3(6)-V3(3)))))+M1*(F2(3)*(
     $ +CI*(V3(5))-V3(4))+F2(4)*(V3(3)+V3(6)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,1)
C     
      SUBROUTINE FFV4_2_S(F1, V3, COUP, M2, W2,F2)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(6)
      COMPLEX*16 V3(*)
      REAL*8 P2(0:3)
      REAL*8 W2
      COMPLEX*16 F1(*)
      REAL*8 M2
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      F2(1) = +F1(1)+V3(1)
      F2(2) = +F1(2)+V3(2)
      P2(0) = -DBLE(F2(1))
      P2(1) = -DBLE(F2(2))
      P2(2) = -DIMAG(F2(2))
      P2(3) = -DIMAG(F2(1))
      DENOM = COUP/(P2(0)**2-P2(1)**2-P2(2)**2-P2(3)**2 - M2 * (M2 -CI
     $ * W2))
      F2(3)= DENOM*CI*(F1(3)*(P2(0)*(V3(3)+V3(6))+(P2(1)*(-1D0)*(V3(4)
     $ +CI*(V3(5)))+(P2(2)*(+CI*(V3(4))-V3(5))-P2(3)*(V3(3)+V3(6)))))
     $ +(F1(4)*(P2(0)*(V3(4)-CI*(V3(5)))+(P2(1)*(V3(6)-V3(3))+(P2(2)*(
     $ -CI*(V3(6))+CI*(V3(3)))+P2(3)*(+CI*(V3(5))-V3(4)))))+M2*(F1(5)
     $ *(V3(3)-V3(6))+F1(6)*(+CI*(V3(5))-V3(4)))))
      F2(4)= DENOM*(-CI)*(F1(3)*(P2(0)*(-1D0)*(V3(4)+CI*(V3(5)))+(P2(1)
     $ *(V3(3)+V3(6))+(P2(2)*(+CI*(V3(3)+V3(6)))-P2(3)*(V3(4)+CI*(V3(5)
     $ )))))+(F1(4)*(P2(0)*(V3(6)-V3(3))+(P2(1)*(V3(4)-CI*(V3(5)))
     $ +(P2(2)*(V3(5)+CI*(V3(4)))+P2(3)*(V3(6)-V3(3)))))+M2*(F1(5)
     $ *(V3(4)+CI*(V3(5)))-F1(6)*(V3(3)+V3(6)))))
      F2(5)= DENOM*(-CI)*(F1(5)*(P2(0)*(V3(6)-V3(3))+(P2(1)*(V3(4)+CI
     $ *(V3(5)))+(P2(2)*(V3(5)-CI*(V3(4)))+P2(3)*(V3(6)-V3(3)))))
     $ +(F1(6)*(P2(0)*(V3(4)-CI*(V3(5)))+(P2(1)*(-1D0)*(V3(3)+V3(6))
     $ +(P2(2)*(+CI*(V3(3)+V3(6)))+P2(3)*(V3(4)-CI*(V3(5))))))+M2
     $ *(F1(3)*(-1D0)*(V3(3)+V3(6))+F1(4)*(+CI*(V3(5))-V3(4)))))
      F2(6)= DENOM*CI*(F1(5)*(P2(0)*(-1D0)*(V3(4)+CI*(V3(5)))+(P2(1)
     $ *(V3(3)-V3(6))+(P2(2)*(-CI*(V3(6))+CI*(V3(3)))+P2(3)*(V3(4)+CI
     $ *(V3(5))))))+(F1(6)*(P2(0)*(V3(3)+V3(6))+(P2(1)*(+CI*(V3(5))
     $ -V3(4))+(P2(2)*(-1D0)*(V3(5)+CI*(V3(4)))-P2(3)*(V3(3)+V3(6)))))
     $ +M2*(F1(3)*(V3(4)+CI*(V3(5)))+F1(4)*(V3(3)-V3(6)))))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,1)
C     
      SUBROUTINE FFV4P0_3(F1, F2, COUP, M3, W3,V3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      COMPLEX*16 V3(6)
      REAL*8 W3
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      V3(1) = +F1(1)+F2(1)
      V3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(V3(1))
      P3(1) = -DBLE(V3(2))
      P3(2) = -DIMAG(V3(2))
      P3(3) = -DIMAG(V3(1))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      V3(3)= DENOM*(-CI)*(F1(3)*F2(5)+F1(4)*F2(6)+F1(5)*F2(3)+F1(6)
     $ *F2(4))
      V3(4)= DENOM*(-CI)*(F1(5)*F2(4)+F1(6)*F2(3)-F1(3)*F2(6)-F1(4)
     $ *F2(5))
      V3(5)= DENOM*(-CI)*(-CI*(F1(3)*F2(6)+F1(6)*F2(3))+CI*(F1(4)*F2(5)
     $ +F1(5)*F2(4)))
      V3(6)= DENOM*(-CI)*(F1(4)*F2(6)+F1(5)*F2(3)-F1(3)*F2(5)-F1(6)
     $ *F2(4))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFV5_0_S(F1, F2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 TMP32
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      TMP32 = (F1(3)*(F2(5)*(V3(3)+V3(6))+F2(6)*(V3(4)+CI*(V3(5))))
     $ +F1(4)*(F2(5)*(V3(4)-CI*(V3(5)))+F2(6)*(V3(3)-V3(6))))
      VERTEX = COUP*(-CI * TMP32)
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFV5_6_0(F1, F2, V3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP2
      COMPLEX*16 V3(*)
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP
      CALL FFV5_0_S(F1,F2,V3,COUP1,VERTEX)
      CALL FFV6_0(F1,F2,V3,COUP2,TMP)
      VERTEX = VERTEX + TMP
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFV5_3_S(F1, F2, COUP, M3, W3,V3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP33
      COMPLEX*16 V3(6)
      REAL*8 W3
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      REAL*8 OM3
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      V3(1) = +F1(1)+F2(1)
      V3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(V3(1))
      P3(1) = -DBLE(V3(2))
      P3(2) = -DIMAG(V3(2))
      P3(3) = -DIMAG(V3(1))
      TMP33 = (F1(3)*(F2(5)*(P3(0)+P3(3))+F2(6)*(P3(1)+CI*(P3(2))))
     $ +F1(4)*(F2(5)*(P3(1)-CI*(P3(2)))+F2(6)*(P3(0)-P3(3))))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      V3(3)= DENOM*(-CI)*(F1(3)*F2(5)+F1(4)*F2(6)-P3(0)*OM3*TMP33)
      V3(4)= DENOM*(-CI)*(-F1(3)*F2(6)-F1(4)*F2(5)-P3(1)*OM3*TMP33)
      V3(5)= DENOM*(-CI)*(-CI*(F1(3)*F2(6))+CI*(F1(4)*F2(5))-P3(2)*OM3
     $ *TMP33)
      V3(6)= DENOM*(-CI)*(F1(4)*F2(6)-F1(3)*F2(5)-P3(3)*OM3*TMP33)
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjM(-1,1)
C     
      SUBROUTINE FFV5_6_3(F1, F2, COUP1, COUP2, M3, W3,V3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V3(6)
      REAL*8 W3
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      COMPLEX*16 COUP2
      REAL*8 OM3
      INTEGER*4 I
      COMPLEX*16 DENOM
      COMPLEX*16 VTMP(6)
      CALL FFV5_3_S(F1,F2,COUP1,M3,W3,V3)
      CALL FFV6_3(F1,F2,COUP2,M3,W3,VTMP)
      DO I = 3, 6
        V3(I) = V3(I) + VTMP(I)
      ENDDO
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFV6_0(F1, F2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP22
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      TMP22 = (F1(5)*(F2(3)*(V3(3)-V3(6))-F2(4)*(V3(4)+CI*(V3(5))))
     $ +F1(6)*(F2(3)*(+CI*(V3(5))-V3(4))+F2(4)*(V3(3)+V3(6))))
      VERTEX = COUP*(-CI * TMP22)
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFV6_3(F1, F2, COUP, M3, W3,V3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 DENOM
      COMPLEX*16 V3(6)
      REAL*8 W3
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      REAL*8 OM3
      COMPLEX*16 COUP
      COMPLEX*16 TMP19
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      V3(1) = +F1(1)+F2(1)
      V3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(V3(1))
      P3(1) = -DBLE(V3(2))
      P3(2) = -DIMAG(V3(2))
      P3(3) = -DIMAG(V3(1))
      TMP19 = (F1(5)*(F2(3)*(P3(0)-P3(3))-F2(4)*(P3(1)+CI*(P3(2))))
     $ +F1(6)*(F2(3)*(+CI*(P3(2))-P3(1))+F2(4)*(P3(0)+P3(3))))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      V3(3)= DENOM*(-CI)*(F1(5)*F2(3)+F1(6)*F2(4)-P3(0)*OM3*TMP19)
      V3(4)= DENOM*(-CI)*(F1(5)*F2(4)+F1(6)*F2(3)-P3(1)*OM3*TMP19)
      V3(5)= DENOM*(-CI)*(-CI*(F1(6)*F2(3))+CI*(F1(5)*F2(4))-P3(2)*OM3
     $ *TMP19)
      V3(6)= DENOM*(-CI)*(F1(5)*F2(3)-F1(6)*F2(4)-P3(3)*OM3*TMP19)
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1003,2)*P(2003,1)*Metric(1,2) + P(1003,1)*P(2003,2)*Metric(1,2)
C      - P(2,1)*P(2003,2)*Metric(1,1003) - P(2,1)*P(1003,2)*Metric(1,20
C     03) - P(1,2)*P(2003,1)*Metric(2,1003) + P(-1,1)*P(-1,2)*Metric(1,
C     2003)*Metric(2,1003) - P(1,2)*P(1003,1)*Metric(2,2003) +
C      P(-1,1)*P(-1,2)*Metric(1,1003)*Metric(2,2003)
C     
      SUBROUTINE VVT2_0(V1, V2, T3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 COUP
      COMPLEX*16 TMP11
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 TMP23
      COMPLEX*16 TMP31
      COMPLEX*16 TMP30
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP26
      COMPLEX*16 TMP28
      COMPLEX*16 TMP27
      COMPLEX*16 TMP29
      COMPLEX*16 TMP24
      COMPLEX*16 T3(*)
      COMPLEX*16 TMP25
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      TMP24 = (P1(0)*(-1D0)*(P2(1)*T3(4)+P2(2)*T3(5)+P2(3)*T3(6)-P2(0)
     $ *T3(3))+(P1(1)*(P2(1)*T3(8)+P2(2)*T3(9)+P2(3)*T3(10)-P2(0)*T3(7)
     $ )+(P1(2)*(P2(1)*T3(12)+P2(2)*T3(13)+P2(3)*T3(14)-P2(0)*T3(11))
     $ +P1(3)*(P2(1)*T3(16)+P2(2)*T3(17)+P2(3)*T3(18)-P2(0)*T3(15)))))
      TMP25 = (P2(0)*(-1D0)*(V1(4)*T3(7)+V1(5)*T3(11)+V1(6)*T3(15)
     $ -V1(3)*T3(3))+(P2(1)*(V1(4)*T3(8)+V1(5)*T3(12)+V1(6)*T3(16)
     $ -V1(3)*T3(4))+(P2(2)*(V1(4)*T3(9)+V1(5)*T3(13)+V1(6)*T3(17)
     $ -V1(3)*T3(5))+P2(3)*(V1(4)*T3(10)+V1(5)*T3(14)+V1(6)*T3(18)
     $ -V1(3)*T3(6)))))
      TMP26 = (P2(0)*(-1D0)*(V1(4)*T3(4)+V1(5)*T3(5)+V1(6)*T3(6)-V1(3)
     $ *T3(3))+(P2(1)*(V1(4)*T3(8)+V1(5)*T3(9)+V1(6)*T3(10)-V1(3)*T3(7)
     $ )+(P2(2)*(V1(4)*T3(12)+V1(5)*T3(13)+V1(6)*T3(14)-V1(3)*T3(11))
     $ +P2(3)*(V1(4)*T3(16)+V1(5)*T3(17)+V1(6)*T3(18)-V1(3)*T3(15)))))
      TMP27 = (P1(0)*(-1D0)*(V2(4)*T3(7)+V2(5)*T3(11)+V2(6)*T3(15)
     $ -V2(3)*T3(3))+(P1(1)*(V2(4)*T3(8)+V2(5)*T3(12)+V2(6)*T3(16)
     $ -V2(3)*T3(4))+(P1(2)*(V2(4)*T3(9)+V2(5)*T3(13)+V2(6)*T3(17)
     $ -V2(3)*T3(5))+P1(3)*(V2(4)*T3(10)+V2(5)*T3(14)+V2(6)*T3(18)
     $ -V2(3)*T3(6)))))
      TMP23 = (P1(0)*(-1D0)*(P2(1)*T3(7)+P2(2)*T3(11)+P2(3)*T3(15)
     $ -P2(0)*T3(3))+(P1(1)*(P2(1)*T3(8)+P2(2)*T3(12)+P2(3)*T3(16)
     $ -P2(0)*T3(4))+(P1(2)*(P2(1)*T3(9)+P2(2)*T3(13)+P2(3)*T3(17)
     $ -P2(0)*T3(5))+P1(3)*(P2(1)*T3(10)+P2(2)*T3(14)+P2(3)*T3(18)
     $ -P2(0)*T3(6)))))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP29 = (V1(3)*(-1D0)*(V2(4)*T3(7)+V2(5)*T3(11)+V2(6)*T3(15)
     $ -V2(3)*T3(3))+(V1(4)*(V2(4)*T3(8)+V2(5)*T3(12)+V2(6)*T3(16)
     $ -V2(3)*T3(4))+(V1(5)*(V2(4)*T3(9)+V2(5)*T3(13)+V2(6)*T3(17)
     $ -V2(3)*T3(5))+V1(6)*(V2(4)*T3(10)+V2(5)*T3(14)+V2(6)*T3(18)
     $ -V2(3)*T3(6)))))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP31 = (V1(3)*(-1D0)*(V2(4)*T3(4)+V2(5)*T3(5)+V2(6)*T3(6)-V2(3)
     $ *T3(3))+(V1(4)*(V2(4)*T3(8)+V2(5)*T3(9)+V2(6)*T3(10)-V2(3)*T3(7)
     $ )+(V1(5)*(V2(4)*T3(12)+V2(5)*T3(13)+V2(6)*T3(14)-V2(3)*T3(11))
     $ +V1(6)*(V2(4)*T3(16)+V2(5)*T3(17)+V2(6)*T3(18)-V2(3)*T3(15)))))
      TMP30 = (P1(0)*(-1D0)*(V2(4)*T3(4)+V2(5)*T3(5)+V2(6)*T3(6)-V2(3)
     $ *T3(3))+(P1(1)*(V2(4)*T3(8)+V2(5)*T3(9)+V2(6)*T3(10)-V2(3)*T3(7)
     $ )+(P1(2)*(V2(4)*T3(12)+V2(5)*T3(13)+V2(6)*T3(14)-V2(3)*T3(11))
     $ +P1(3)*(V2(4)*T3(16)+V2(5)*T3(17)+V2(6)*T3(18)-V2(3)*T3(15)))))
      TMP11 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      VERTEX = COUP*(TMP11*(+CI*(TMP27+TMP30))+(TMP28*(-1D0)*(+CI
     $ *(TMP29+TMP31))+(TMP3*(-1D0)*(+CI*(TMP23+TMP24))+TMP9*(+CI
     $ *(TMP25+TMP26)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1003,2)*P(2003,1)*Metric(1,2) + P(1003,1)*P(2003,2)*Metric(1,2)
C      - P(2,1)*P(2003,2)*Metric(1,1003) - P(2,1)*P(1003,2)*Metric(1,20
C     03) - P(1,2)*P(2003,1)*Metric(2,1003) + P(-1,1)*P(-1,2)*Metric(1,
C     2003)*Metric(2,1003) - P(1,2)*P(1003,1)*Metric(2,2003) +
C      P(-1,1)*P(-1,2)*Metric(1,1003)*Metric(2,2003)
C     
      SUBROUTINE VVT2_3(V1, V2, COUP, M3, W3,T3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP12
      COMPLEX*16 TMP11
      REAL*8 P1(0:3)
      REAL*8 W3
      COMPLEX*16 TMP10
      REAL*8 P2(0:3)
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 TMP21
      COMPLEX*16 DENOM
      REAL*8 OM3
      COMPLEX*16 TMP28
      COMPLEX*16 T3(18)
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      COMPLEX*16 TMP18
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      T3(1) = +V1(1)+V2(1)
      T3(2) = +V1(2)+V2(2)
      P3(0) = -DBLE(T3(1))
      P3(1) = -DBLE(T3(2))
      P3(2) = -DIMAG(T3(2))
      P3(3) = -DIMAG(T3(1))
      TMP21 = (P2(0)*P3(0)-P2(1)*P3(1)-P2(2)*P3(2)-P2(3)*P3(3))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP18 = (P1(0)*P3(0)-P1(1)*P3(1)-P1(2)*P3(2)-P1(3)*P3(3))
      TMP11 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      TMP10 = (V2(3)*P3(0)-V2(4)*P3(1)-V2(5)*P3(2)-V2(6)*P3(3))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP12 = (V1(3)*P3(0)-V1(4)*P3(1)-V1(5)*P3(2)-V1(6)*P3(3))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      T3(3)= DENOM*2D0*(OM3*(P3(0)*(P3(0)*(OM3*(TMP10*2D0/3D0*(-CI
     $ *(TMP11*TMP18)+CI*(TMP12*TMP28))+2D0/3D0*(TMP21*(-CI*(TMP9
     $ *TMP12)+CI*(TMP3*TMP18))))+(-2D0/3D0 * CI*(TMP9*TMP11)+2D0/3D0 
     $ * CI*(TMP3*TMP28)))+(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))
     $ +(P2(0)*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI
     $ *(V2(3)*TMP12+V1(3)*TMP10))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11
     $ *TMP18))))))+(TMP10*1D0/3D0*(-CI*(TMP11*TMP18)+CI*(TMP12*TMP28))
     $ +1D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI*(TMP3*TMP18)))))+(TMP11*(
     $ -CI*(V2(3)*P1(0))+2D0/3D0 * CI*(TMP9))+(TMP28*(-2D0/3D0 * CI
     $ *(TMP3)+CI*(V2(3)*V1(3)))+P2(0)*(-CI*(V1(3)*TMP9)+CI*(TMP3*P1(0)
     $ )))))
      T3(4)= DENOM*(OM3*(P3(0)*(P3(1)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(4)
     $ *TMP10+V2(4)*TMP12))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+P3(1)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(3)
     $ *TMP12+V1(3)*TMP10))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(P1(1)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(3)*V1(4)+V2(4)
     $ *V1(3)))-TMP9*(+CI*(V1(3)*P2(1)+V1(4)*P2(0)))))))
      T3(5)= DENOM*(OM3*(P3(0)*(P3(2)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(5)
     $ *TMP10+V2(5)*TMP12))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+P3(2)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(3)
     $ *TMP12+V1(3)*TMP10))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(P1(2)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(3)*V1(5)+V2(5)
     $ *V1(3)))-TMP9*(+CI*(V1(3)*P2(2)+V1(5)*P2(0)))))))
      T3(6)= DENOM*(OM3*(P3(0)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(6)
     $ *TMP10+V2(6)*TMP12))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(3)
     $ *TMP12+V1(3)*TMP10))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(3)*V1(6)+V2(6)
     $ *V1(3)))-TMP9*(+CI*(V1(3)*P2(3)+V1(6)*P2(0)))))))
      T3(7)= DENOM*(OM3*(P3(0)*(P3(1)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(4)
     $ *TMP12+V1(4)*TMP10))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+P3(1)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(3)
     $ *TMP10+V2(3)*TMP12))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(P1(1)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(4)*V1(3)+V2(3)
     $ *V1(4)))-TMP9*(+CI*(V1(4)*P2(0)+V1(3)*P2(1)))))))
      T3(8)= DENOM*2D0*(OM3*(P3(1)*(P3(1)*(OM3*(TMP10*2D0/3D0*(-CI
     $ *(TMP11*TMP18)+CI*(TMP12*TMP28))+2D0/3D0*(TMP21*(-CI*(TMP9
     $ *TMP12)+CI*(TMP3*TMP18))))+(-2D0/3D0 * CI*(TMP9*TMP11)+2D0/3D0 
     $ * CI*(TMP3*TMP28)))+(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))
     $ +(P2(1)*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI
     $ *(V2(4)*TMP12+V1(4)*TMP10))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11
     $ *TMP18))))))+(TMP10*1D0/3D0*(-CI*(TMP12*TMP28)+CI*(TMP11*TMP18))
     $ +1D0/3D0*(TMP21*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12)))))+(TMP11*(
     $ -1D0)*(+CI*(V2(4)*P1(1))+2D0/3D0 * CI*(TMP9))+(TMP28*(+CI*(V2(4)
     $ *V1(4))+2D0/3D0 * CI*(TMP3))+P2(1)*(-CI*(V1(4)*TMP9)+CI*(TMP3
     $ *P1(1))))))
      T3(9)= DENOM*(OM3*(P3(1)*(P3(2)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(5)
     $ *TMP10+V2(5)*TMP12))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+P3(2)*(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(4)
     $ *TMP12+V1(4)*TMP10))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+(P1(1)*(-CI*(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(P1(2)*(-CI
     $ *(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(TMP28*(+CI*(V2(4)*V1(5)+V2(5)
     $ *V1(4)))-TMP9*(+CI*(V1(4)*P2(2)+V1(5)*P2(1)))))))
      T3(10)= DENOM*(OM3*(P3(1)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(6)
     $ *TMP10+V2(6)*TMP12))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(4)
     $ *TMP12+V1(4)*TMP10))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+(P1(1)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(TMP28*(+CI*(V2(4)*V1(6)+V2(6)
     $ *V1(4)))-TMP9*(+CI*(V1(4)*P2(3)+V1(6)*P2(1)))))))
      T3(11)= DENOM*(OM3*(P3(0)*(P3(2)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(5)
     $ *TMP12+V1(5)*TMP10))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+P3(2)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(3)
     $ *TMP10+V2(3)*TMP12))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(P1(2)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(5)*V1(3)+V2(3)
     $ *V1(5)))-TMP9*(+CI*(V1(5)*P2(0)+V1(3)*P2(2)))))))
      T3(12)= DENOM*(OM3*(P3(1)*(P3(2)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(5)
     $ *TMP12+V1(5)*TMP10))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+P3(2)*(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(4)
     $ *TMP10+V2(4)*TMP12))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+(P1(1)*(-CI*(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(P1(2)*(-CI
     $ *(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(TMP28*(+CI*(V2(5)*V1(4)+V2(4)
     $ *V1(5)))-TMP9*(+CI*(V1(5)*P2(1)+V1(4)*P2(2)))))))
      T3(13)= DENOM*2D0*(OM3*(P3(2)*(P3(2)*(OM3*(TMP10*2D0/3D0*(-CI
     $ *(TMP11*TMP18)+CI*(TMP12*TMP28))+2D0/3D0*(TMP21*(-CI*(TMP9
     $ *TMP12)+CI*(TMP3*TMP18))))+(-2D0/3D0 * CI*(TMP9*TMP11)+2D0/3D0 
     $ * CI*(TMP3*TMP28)))+(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))
     $ +(P2(2)*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI
     $ *(V2(5)*TMP12+V1(5)*TMP10))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11
     $ *TMP18))))))+(TMP10*1D0/3D0*(-CI*(TMP12*TMP28)+CI*(TMP11*TMP18))
     $ +1D0/3D0*(TMP21*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12)))))+(TMP11*(
     $ -1D0)*(+CI*(V2(5)*P1(2))+2D0/3D0 * CI*(TMP9))+(TMP28*(+CI*(V2(5)
     $ *V1(5))+2D0/3D0 * CI*(TMP3))+P2(2)*(-CI*(V1(5)*TMP9)+CI*(TMP3
     $ *P1(2))))))
      T3(14)= DENOM*(OM3*(P3(2)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(6)
     $ *TMP10+V2(6)*TMP12))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(5)
     $ *TMP12+V1(5)*TMP10))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+(P1(2)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(TMP28*(+CI*(V2(5)*V1(6)+V2(6)
     $ *V1(5)))-TMP9*(+CI*(V1(5)*P2(3)+V1(6)*P2(2)))))))
      T3(15)= DENOM*(OM3*(P3(0)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(6)
     $ *TMP12+V1(6)*TMP10))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(0)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(0)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(3)
     $ *TMP10+V2(3)*TMP12))+(+CI*(V1(3)*TMP9*TMP21+V2(3)*TMP11*TMP18)))
     $ )))+(P1(0)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(3)*TMP11)+CI*(TMP3*P2(0)))+(TMP28*(+CI*(V2(6)*V1(3)+V2(3)
     $ *V1(6)))-TMP9*(+CI*(V1(6)*P2(0)+V1(3)*P2(3)))))))
      T3(16)= DENOM*(OM3*(P3(1)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(6)
     $ *TMP12+V1(6)*TMP10))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(1)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(1)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(4)
     $ *TMP10+V2(4)*TMP12))+(+CI*(V1(4)*TMP9*TMP21+V2(4)*TMP11*TMP18)))
     $ )))+(P1(1)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(4)*TMP11)+CI*(TMP3*P2(1)))+(TMP28*(+CI*(V2(6)*V1(4)+V2(4)
     $ *V1(6)))-TMP9*(+CI*(V1(6)*P2(1)+V1(4)*P2(3)))))))
      T3(17)= DENOM*(OM3*(P3(2)*(P3(3)*(OM3*(TMP10*4D0/3D0*(-CI*(TMP11
     $ *TMP18)+CI*(TMP12*TMP28))+4D0/3D0*(TMP21*(-CI*(TMP9*TMP12)+CI
     $ *(TMP3*TMP18))))+(-4D0/3D0 * CI*(TMP9*TMP11)+4D0/3D0 * CI*(TMP3
     $ *TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(3)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V2(6)
     $ *TMP12+V1(6)*TMP10))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11*TMP18)))
     $ )))+P3(3)*(P1(2)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))+(P2(2)*(
     $ -CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI*(V1(5)
     $ *TMP10+V2(5)*TMP12))+(+CI*(V1(5)*TMP9*TMP21+V2(5)*TMP11*TMP18)))
     $ )))+(P1(2)*(-CI*(V2(6)*TMP11)+CI*(TMP3*P2(3)))+(P1(3)*(-CI
     $ *(V2(5)*TMP11)+CI*(TMP3*P2(2)))+(TMP28*(+CI*(V2(6)*V1(5)+V2(5)
     $ *V1(6)))-TMP9*(+CI*(V1(6)*P2(2)+V1(5)*P2(3)))))))
      T3(18)= DENOM*2D0*(OM3*(P3(3)*(P3(3)*(OM3*(TMP10*2D0/3D0*(-CI
     $ *(TMP11*TMP18)+CI*(TMP12*TMP28))+2D0/3D0*(TMP21*(-CI*(TMP9
     $ *TMP12)+CI*(TMP3*TMP18))))+(-2D0/3D0 * CI*(TMP9*TMP11)+2D0/3D0 
     $ * CI*(TMP3*TMP28)))+(P1(3)*(-CI*(TMP3*TMP21)+CI*(TMP10*TMP11))
     $ +(P2(3)*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12))+(TMP28*(-1D0)*(+CI
     $ *(V2(6)*TMP12+V1(6)*TMP10))+(+CI*(V1(6)*TMP9*TMP21+V2(6)*TMP11
     $ *TMP18))))))+(TMP10*1D0/3D0*(-CI*(TMP12*TMP28)+CI*(TMP11*TMP18))
     $ +1D0/3D0*(TMP21*(-CI*(TMP3*TMP18)+CI*(TMP9*TMP12)))))+(TMP11*(
     $ -1D0)*(+CI*(V2(6)*P1(3))+2D0/3D0 * CI*(TMP9))+(TMP28*(+CI*(V2(6)
     $ *V1(6))+2D0/3D0 * CI*(TMP3))+P2(3)*(-CI*(V1(6)*TMP9)+CI*(TMP3
     $ *P1(3))))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1003,2)*P(2003,1)*Metric(1,2) + P(1003,1)*P(2003,2)*Metric(1,2)
C      - P(2,1)*P(2003,2)*Metric(1,1003) - P(2,1)*P(1003,2)*Metric(1,20
C     03) - P(1,2)*P(2003,1)*Metric(2,1003) + P(-1,1)*P(-1,2)*Metric(1,
C     2003)*Metric(2,1003) - P(1,2)*P(1003,1)*Metric(2,2003) +
C      P(-1,1)*P(-1,2)*Metric(1,1003)*Metric(2,2003)
C     
      SUBROUTINE VVT2P0_1(V2, T3, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      COMPLEX*16 TMP23
      REAL*8 W1
      COMPLEX*16 TMP30
      COMPLEX*16 DENOM
      COMPLEX*16 TMP28
      COMPLEX*16 TMP27
      COMPLEX*16 TMP24
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      COMPLEX*16 V1(6)
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      V1(1) = +V2(1)+T3(1)
      V1(2) = +V2(2)+T3(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP24 = (P1(0)*(-1D0)*(P2(1)*T3(4)+P2(2)*T3(5)+P2(3)*T3(6)-P2(0)
     $ *T3(3))+(P1(1)*(P2(1)*T3(8)+P2(2)*T3(9)+P2(3)*T3(10)-P2(0)*T3(7)
     $ )+(P1(2)*(P2(1)*T3(12)+P2(2)*T3(13)+P2(3)*T3(14)-P2(0)*T3(11))
     $ +P1(3)*(P2(1)*T3(16)+P2(2)*T3(17)+P2(3)*T3(18)-P2(0)*T3(15)))))
      TMP27 = (P1(0)*(-1D0)*(V2(4)*T3(7)+V2(5)*T3(11)+V2(6)*T3(15)
     $ -V2(3)*T3(3))+(P1(1)*(V2(4)*T3(8)+V2(5)*T3(12)+V2(6)*T3(16)
     $ -V2(3)*T3(4))+(P1(2)*(V2(4)*T3(9)+V2(5)*T3(13)+V2(6)*T3(17)
     $ -V2(3)*T3(5))+P1(3)*(V2(4)*T3(10)+V2(5)*T3(14)+V2(6)*T3(18)
     $ -V2(3)*T3(6)))))
      TMP23 = (P1(0)*(-1D0)*(P2(1)*T3(7)+P2(2)*T3(11)+P2(3)*T3(15)
     $ -P2(0)*T3(3))+(P1(1)*(P2(1)*T3(8)+P2(2)*T3(12)+P2(3)*T3(16)
     $ -P2(0)*T3(4))+(P1(2)*(P2(1)*T3(9)+P2(2)*T3(13)+P2(3)*T3(17)
     $ -P2(0)*T3(5))+P1(3)*(P2(1)*T3(10)+P2(2)*T3(14)+P2(3)*T3(18)
     $ -P2(0)*T3(6)))))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP30 = (P1(0)*(-1D0)*(V2(4)*T3(4)+V2(5)*T3(5)+V2(6)*T3(6)-V2(3)
     $ *T3(3))+(P1(1)*(V2(4)*T3(8)+V2(5)*T3(9)+V2(6)*T3(10)-V2(3)*T3(7)
     $ )+(P1(2)*(V2(4)*T3(12)+V2(5)*T3(13)+V2(6)*T3(14)-V2(3)*T3(11))
     $ +P1(3)*(V2(4)*T3(16)+V2(5)*T3(17)+V2(6)*T3(18)-V2(3)*T3(15)))))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(TMP28*(V2(4)*(+CI*(T3(7)+T3(4)))+(V2(5)*(+CI
     $ *(T3(11)+T3(5)))+(V2(6)*(+CI*(T3(15)+T3(6)))-2D0 * CI*(V2(3)
     $ *T3(3)))))+(TMP9*(P2(1)*(-1D0)*(+CI*(T3(4)+T3(7)))+(P2(2)*(-1D0)
     $ *(+CI*(T3(5)+T3(11)))+(P2(3)*(-1D0)*(+CI*(T3(6)+T3(15)))+2D0 *
     $  CI*(P2(0)*T3(3)))))+(P2(0)*(+CI*(TMP27+TMP30))-V2(3)*(+CI
     $ *(TMP23+TMP24)))))
      V1(4)= DENOM*(TMP28*(V2(3)*(-1D0)*(+CI*(T3(4)+T3(7)))+(V2(5)*(
     $ +CI*(T3(12)+T3(9)))+(V2(6)*(+CI*(T3(16)+T3(10)))+2D0 * CI*(V2(4)
     $ *T3(8)))))+(TMP9*(P2(0)*(+CI*(T3(7)+T3(4)))+(P2(2)*(-1D0)*(+CI
     $ *(T3(9)+T3(12)))+(P2(3)*(-1D0)*(+CI*(T3(10)+T3(16)))-2D0 * CI
     $ *(P2(1)*T3(8)))))+(P2(1)*(+CI*(TMP27+TMP30))-V2(4)*(+CI*(TMP23
     $ +TMP24)))))
      V1(5)= DENOM*(TMP28*(V2(3)*(-1D0)*(+CI*(T3(5)+T3(11)))+(V2(4)*(
     $ +CI*(T3(9)+T3(12)))+(V2(6)*(+CI*(T3(17)+T3(14)))+2D0 * CI*(V2(5)
     $ *T3(13)))))+(TMP9*(P2(0)*(+CI*(T3(11)+T3(5)))+(P2(1)*(-1D0)*(
     $ +CI*(T3(12)+T3(9)))+(P2(3)*(-1D0)*(+CI*(T3(14)+T3(17)))-2D0 *
     $  CI*(P2(2)*T3(13)))))+(P2(2)*(+CI*(TMP27+TMP30))-V2(5)*(+CI
     $ *(TMP23+TMP24)))))
      V1(6)= DENOM*(TMP28*(V2(3)*(-1D0)*(+CI*(T3(6)+T3(15)))+(V2(4)*(
     $ +CI*(T3(10)+T3(16)))+(V2(5)*(+CI*(T3(14)+T3(17)))+2D0 * CI
     $ *(V2(6)*T3(18)))))+(TMP9*(P2(0)*(+CI*(T3(15)+T3(6)))+(P2(1)*(
     $ -1D0)*(+CI*(T3(16)+T3(10)))+(P2(2)*(-1D0)*(+CI*(T3(17)+T3(14)))
     $ -2D0 * CI*(P2(3)*T3(18)))))+(P2(3)*(+CI*(TMP27+TMP30))-V2(6)*(
     $ +CI*(TMP23+TMP24)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
C      P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
C     
      SUBROUTINE VVV2_0(V1, V2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP12
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP1
      REAL*8 P1(0:3)
      COMPLEX*16 TMP10
      REAL*8 P2(0:3)
      COMPLEX*16 TMP7
      REAL*8 P3(0:3)
      COMPLEX*16 TMP5
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP11
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      COMPLEX*16 TMP8
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP8 = (V3(3)*P2(0)-V3(4)*P2(1)-V3(5)*P2(2)-V3(6)*P2(3))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP7 = (V3(3)*P1(0)-V3(4)*P1(1)-V3(5)*P1(2)-V3(6)*P1(3))
      TMP11 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      TMP10 = (V2(3)*P3(0)-V2(4)*P3(1)-V2(5)*P3(2)-V2(6)*P3(3))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP12 = (V1(3)*P3(0)-V1(4)*P3(1)-V1(5)*P3(2)-V1(6)*P3(3))
      VERTEX = COUP*(TMP1*(-CI*(TMP10)+CI*(TMP9))+(TMP3*(-CI*(TMP7)+CI
     $ *(TMP8))+TMP5*(-CI*(TMP11)+CI*(TMP12))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
C      P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
C     
      SUBROUTINE VVV2P0_1(V2, V3, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 DENOM
      COMPLEX*16 V3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      COMPLEX*16 TMP10
      REAL*8 P2(0:3)
      COMPLEX*16 TMP7
      REAL*8 P3(0:3)
      REAL*8 W1
      COMPLEX*16 TMP5
      COMPLEX*16 COUP
      COMPLEX*16 TMP9
      COMPLEX*16 V1(6)
      COMPLEX*16 TMP8
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      V1(1) = +V2(1)+V3(1)
      V1(2) = +V2(2)+V3(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP7 = (V3(3)*P1(0)-V3(4)*P1(1)-V3(5)*P1(2)-V3(6)*P1(3))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP8 = (V3(3)*P2(0)-V3(4)*P2(1)-V3(5)*P2(2)-V3(6)*P2(3))
      TMP10 = (V2(3)*P3(0)-V2(4)*P3(1)-V2(5)*P3(2)-V2(6)*P3(3))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(TMP5*(-CI*(P2(0))+CI*(P3(0)))+(V2(3)*(-CI*(TMP7)
     $ +CI*(TMP8))+V3(3)*(-CI*(TMP10)+CI*(TMP9))))
      V1(4)= DENOM*(TMP5*(-CI*(P2(1))+CI*(P3(1)))+(V2(4)*(-CI*(TMP7)
     $ +CI*(TMP8))+V3(4)*(-CI*(TMP10)+CI*(TMP9))))
      V1(5)= DENOM*(TMP5*(-CI*(P2(2))+CI*(P3(2)))+(V2(5)*(-CI*(TMP7)
     $ +CI*(TMP8))+V3(5)*(-CI*(TMP10)+CI*(TMP9))))
      V1(6)= DENOM*(TMP5*(-CI*(P2(3))+CI*(P3(3)))+(V2(6)*(-CI*(TMP7)
     $ +CI*(TMP8))+V3(6)*(-CI*(TMP10)+CI*(TMP9))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2004,2)*Metric(1,1004)*Metric(2,3) - P(2004,3)*Metric(1,1004)*M
C     etric(2,3) + P(1004,2)*Metric(1,2004)*Metric(2,3) -
C      P(1004,3)*Metric(1,2004)*Metric(2,3) - P(2004,1)*Metric(1,3)*Met
C     ric(2,1004) + P(2004,3)*Metric(1,3)*Metric(2,1004) +
C      P(3,1)*Metric(1,2004)*Metric(2,1004) - P(3,2)*Metric(1,2004)*Met
C     ric(2,1004) - P(1004,1)*Metric(1,3)*Metric(2,2004) +
C      P(1004,3)*Metric(1,3)*Metric(2,2004) + P(3,1)*Metric(1,1004)*Met
C     ric(2,2004) - P(3,2)*Metric(1,1004)*Metric(2,2004) +
C      P(2004,1)*Metric(1,2)*Metric(3,1004) - P(2004,2)*Metric(1,2)*Met
C     ric(3,1004) - P(2,1)*Metric(1,2004)*Metric(3,1004) +
C      P(2,3)*Metric(1,2004)*Metric(3,1004) + P(1,2)*Metric(2,2004)*Met
C     ric(3,1004) - P(1,3)*Metric(2,2004)*Metric(3,1004) +
C      P(1004,1)*Metric(1,2)*Metric(3,2004) - P(1004,2)*Metric(1,2)*Met
C     ric(3,2004) - P(2,1)*Metric(1,1004)*Metric(3,2004) +
C      P(2,3)*Metric(1,1004)*Metric(3,2004) + P(1,2)*Metric(2,1004)*Met
C     ric(3,2004) - P(1,3)*Metric(2,1004)*Metric(3,2004)
C     
      SUBROUTINE VVVT2_0(V1, V2, V3, T4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP61
      COMPLEX*16 TMP1
      COMPLEX*16 TMP10
      COMPLEX*16 T4(*)
      COMPLEX*16 TMP56
      REAL*8 P3(0:3)
      COMPLEX*16 TMP5
      COMPLEX*16 TMP64
      COMPLEX*16 TMP52
      COMPLEX*16 TMP69
      COMPLEX*16 TMP9
      COMPLEX*16 TMP60
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP11
      COMPLEX*16 TMP57
      REAL*8 P2(0:3)
      COMPLEX*16 TMP67
      COMPLEX*16 TMP53
      COMPLEX*16 TMP68
      COMPLEX*16 TMP63
      COMPLEX*16 TMP3
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP12
      COMPLEX*16 TMP58
      REAL*8 P1(0:3)
      COMPLEX*16 TMP7
      COMPLEX*16 TMP66
      COMPLEX*16 TMP54
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP62
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP59
      COMPLEX*16 TMP55
      COMPLEX*16 TMP65
      COMPLEX*16 COUP
      COMPLEX*16 TMP8
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      TMP68 = (V1(3)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(V1(4)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(V1(5)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +V1(6)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP69 = (V2(3)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(V2(4)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(V2(5)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +V2(6)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP60 = (P3(0)*(-1D0)*(V2(4)*T4(4)+V2(5)*T4(5)+V2(6)*T4(6)-V2(3)
     $ *T4(3))+(P3(1)*(V2(4)*T4(8)+V2(5)*T4(9)+V2(6)*T4(10)-V2(3)*T4(7)
     $ )+(P3(2)*(V2(4)*T4(12)+V2(5)*T4(13)+V2(6)*T4(14)-V2(3)*T4(11))
     $ +P3(3)*(V2(4)*T4(16)+V2(5)*T4(17)+V2(6)*T4(18)-V2(3)*T4(15)))))
      TMP61 = (V1(3)*(-1D0)*(V2(4)*T4(4)+V2(5)*T4(5)+V2(6)*T4(6)-V2(3)
     $ *T4(3))+(V1(4)*(V2(4)*T4(8)+V2(5)*T4(9)+V2(6)*T4(10)-V2(3)*T4(7)
     $ )+(V1(5)*(V2(4)*T4(12)+V2(5)*T4(13)+V2(6)*T4(14)-V2(3)*T4(11))
     $ +V1(6)*(V2(4)*T4(16)+V2(5)*T4(17)+V2(6)*T4(18)-V2(3)*T4(15)))))
      TMP62 = (P1(0)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(P1(1)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(P1(2)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+P1(3)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP63 = (P2(0)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(P2(1)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(P2(2)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+P2(3)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP64 = (V1(3)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(V1(4)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(V1(5)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+V1(6)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP65 = (V2(3)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(V2(4)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(V2(5)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+V2(6)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP66 = (P1(0)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(P1(1)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(P1(2)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +P1(3)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP67 = (P2(0)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(P2(1)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(P2(2)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +P2(3)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP53 = (P3(0)*(-1D0)*(V1(4)*T4(7)+V1(5)*T4(11)+V1(6)*T4(15)
     $ -V1(3)*T4(3))+(P3(1)*(V1(4)*T4(8)+V1(5)*T4(12)+V1(6)*T4(16)
     $ -V1(3)*T4(4))+(P3(2)*(V1(4)*T4(9)+V1(5)*T4(13)+V1(6)*T4(17)
     $ -V1(3)*T4(5))+P3(3)*(V1(4)*T4(10)+V1(5)*T4(14)+V1(6)*T4(18)
     $ -V1(3)*T4(6)))))
      TMP52 = (P2(0)*(-1D0)*(V1(4)*T4(7)+V1(5)*T4(11)+V1(6)*T4(15)
     $ -V1(3)*T4(3))+(P2(1)*(V1(4)*T4(8)+V1(5)*T4(12)+V1(6)*T4(16)
     $ -V1(3)*T4(4))+(P2(2)*(V1(4)*T4(9)+V1(5)*T4(13)+V1(6)*T4(17)
     $ -V1(3)*T4(5))+P2(3)*(V1(4)*T4(10)+V1(5)*T4(14)+V1(6)*T4(18)
     $ -V1(3)*T4(6)))))
      TMP55 = (P3(0)*(-1D0)*(V1(4)*T4(4)+V1(5)*T4(5)+V1(6)*T4(6)-V1(3)
     $ *T4(3))+(P3(1)*(V1(4)*T4(8)+V1(5)*T4(9)+V1(6)*T4(10)-V1(3)*T4(7)
     $ )+(P3(2)*(V1(4)*T4(12)+V1(5)*T4(13)+V1(6)*T4(14)-V1(3)*T4(11))
     $ +P3(3)*(V1(4)*T4(16)+V1(5)*T4(17)+V1(6)*T4(18)-V1(3)*T4(15)))))
      TMP54 = (P2(0)*(-1D0)*(V1(4)*T4(4)+V1(5)*T4(5)+V1(6)*T4(6)-V1(3)
     $ *T4(3))+(P2(1)*(V1(4)*T4(8)+V1(5)*T4(9)+V1(6)*T4(10)-V1(3)*T4(7)
     $ )+(P2(2)*(V1(4)*T4(12)+V1(5)*T4(13)+V1(6)*T4(14)-V1(3)*T4(11))
     $ +P2(3)*(V1(4)*T4(16)+V1(5)*T4(17)+V1(6)*T4(18)-V1(3)*T4(15)))))
      TMP57 = (P3(0)*(-1D0)*(V2(4)*T4(7)+V2(5)*T4(11)+V2(6)*T4(15)
     $ -V2(3)*T4(3))+(P3(1)*(V2(4)*T4(8)+V2(5)*T4(12)+V2(6)*T4(16)
     $ -V2(3)*T4(4))+(P3(2)*(V2(4)*T4(9)+V2(5)*T4(13)+V2(6)*T4(17)
     $ -V2(3)*T4(5))+P3(3)*(V2(4)*T4(10)+V2(5)*T4(14)+V2(6)*T4(18)
     $ -V2(3)*T4(6)))))
      TMP56 = (P1(0)*(-1D0)*(V2(4)*T4(7)+V2(5)*T4(11)+V2(6)*T4(15)
     $ -V2(3)*T4(3))+(P1(1)*(V2(4)*T4(8)+V2(5)*T4(12)+V2(6)*T4(16)
     $ -V2(3)*T4(4))+(P1(2)*(V2(4)*T4(9)+V2(5)*T4(13)+V2(6)*T4(17)
     $ -V2(3)*T4(5))+P1(3)*(V2(4)*T4(10)+V2(5)*T4(14)+V2(6)*T4(18)
     $ -V2(3)*T4(6)))))
      TMP59 = (P1(0)*(-1D0)*(V2(4)*T4(4)+V2(5)*T4(5)+V2(6)*T4(6)-V2(3)
     $ *T4(3))+(P1(1)*(V2(4)*T4(8)+V2(5)*T4(9)+V2(6)*T4(10)-V2(3)*T4(7)
     $ )+(P1(2)*(V2(4)*T4(12)+V2(5)*T4(13)+V2(6)*T4(14)-V2(3)*T4(11))
     $ +P1(3)*(V2(4)*T4(16)+V2(5)*T4(17)+V2(6)*T4(18)-V2(3)*T4(15)))))
      TMP58 = (V1(3)*(-1D0)*(V2(4)*T4(7)+V2(5)*T4(11)+V2(6)*T4(15)
     $ -V2(3)*T4(3))+(V1(4)*(V2(4)*T4(8)+V2(5)*T4(12)+V2(6)*T4(16)
     $ -V2(3)*T4(4))+(V1(5)*(V2(4)*T4(9)+V2(5)*T4(13)+V2(6)*T4(17)
     $ -V2(3)*T4(5))+V1(6)*(V2(4)*T4(10)+V2(5)*T4(14)+V2(6)*T4(18)
     $ -V2(3)*T4(6)))))
      TMP11 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      TMP10 = (V2(3)*P3(0)-V2(4)*P3(1)-V2(5)*P3(2)-V2(6)*P3(3))
      TMP12 = (V1(3)*P3(0)-V1(4)*P3(1)-V1(5)*P3(2)-V1(6)*P3(3))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP8 = (V3(3)*P2(0)-V3(4)*P2(1)-V3(5)*P2(2)-V3(6)*P2(3))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP7 = (V3(3)*P1(0)-V3(4)*P1(1)-V3(5)*P1(2)-V3(6)*P1(3))
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      VERTEX = COUP*(TMP1*(-CI*(TMP57+TMP60)+CI*(TMP56+TMP59))+(TMP3*(
     $ -CI*(TMP62+TMP66)+CI*(TMP63+TMP67))+(TMP5*(-CI*(TMP52+TMP54)+CI
     $ *(TMP53+TMP55))+(TMP10*(-1D0)*(+CI*(TMP64+TMP68))+(TMP11*(-1D0)
     $ *(+CI*(TMP65+TMP69))+(TMP12*(+CI*(TMP65+TMP69))+(TMP58*(-CI
     $ *(TMP7)+CI*(TMP8))+(TMP61*(-CI*(TMP7)+CI*(TMP8))+TMP9*(+CI
     $ *(TMP64+TMP68))))))))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2004,2)*Metric(1,1004)*Metric(2,3) - P(2004,3)*Metric(1,1004)*M
C     etric(2,3) + P(1004,2)*Metric(1,2004)*Metric(2,3) -
C      P(1004,3)*Metric(1,2004)*Metric(2,3) - P(2004,1)*Metric(1,3)*Met
C     ric(2,1004) + P(2004,3)*Metric(1,3)*Metric(2,1004) +
C      P(3,1)*Metric(1,2004)*Metric(2,1004) - P(3,2)*Metric(1,2004)*Met
C     ric(2,1004) - P(1004,1)*Metric(1,3)*Metric(2,2004) +
C      P(1004,3)*Metric(1,3)*Metric(2,2004) + P(3,1)*Metric(1,1004)*Met
C     ric(2,2004) - P(3,2)*Metric(1,1004)*Metric(2,2004) +
C      P(2004,1)*Metric(1,2)*Metric(3,1004) - P(2004,2)*Metric(1,2)*Met
C     ric(3,1004) - P(2,1)*Metric(1,2004)*Metric(3,1004) +
C      P(2,3)*Metric(1,2004)*Metric(3,1004) + P(1,2)*Metric(2,2004)*Met
C     ric(3,1004) - P(1,3)*Metric(2,2004)*Metric(3,1004) +
C      P(1004,1)*Metric(1,2)*Metric(3,2004) - P(1004,2)*Metric(1,2)*Met
C     ric(3,2004) - P(2,1)*Metric(1,1004)*Metric(3,2004) +
C      P(2,3)*Metric(1,1004)*Metric(3,2004) + P(1,2)*Metric(2,1004)*Met
C     ric(3,2004) - P(1,3)*Metric(2,1004)*Metric(3,2004)
C     
      SUBROUTINE VVVT2P0_1(V2, V3, T4, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP10
      COMPLEX*16 T4(*)
      COMPLEX*16 TMP56
      REAL*8 P3(0:3)
      REAL*8 W1
      COMPLEX*16 TMP5
      COMPLEX*16 TMP69
      COMPLEX*16 TMP9
      COMPLEX*16 TMP60
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP57
      REAL*8 P2(0:3)
      COMPLEX*16 TMP67
      COMPLEX*16 TMP63
      COMPLEX*16 V2(*)
      REAL*8 P1(0:3)
      COMPLEX*16 TMP7
      COMPLEX*16 TMP66
      COMPLEX*16 DENOM
      COMPLEX*16 TMP62
      COMPLEX*16 V1(6)
      COMPLEX*16 TMP59
      REAL*8 M1
      COMPLEX*16 TMP65
      COMPLEX*16 COUP
      COMPLEX*16 TMP8
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      V1(1) = +V2(1)+V3(1)+T4(1)
      V1(2) = +V2(2)+V3(2)+T4(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP69 = (V2(3)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(V2(4)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(V2(5)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +V2(6)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP57 = (P3(0)*(-1D0)*(V2(4)*T4(7)+V2(5)*T4(11)+V2(6)*T4(15)
     $ -V2(3)*T4(3))+(P3(1)*(V2(4)*T4(8)+V2(5)*T4(12)+V2(6)*T4(16)
     $ -V2(3)*T4(4))+(P3(2)*(V2(4)*T4(9)+V2(5)*T4(13)+V2(6)*T4(17)
     $ -V2(3)*T4(5))+P3(3)*(V2(4)*T4(10)+V2(5)*T4(14)+V2(6)*T4(18)
     $ -V2(3)*T4(6)))))
      TMP56 = (P1(0)*(-1D0)*(V2(4)*T4(7)+V2(5)*T4(11)+V2(6)*T4(15)
     $ -V2(3)*T4(3))+(P1(1)*(V2(4)*T4(8)+V2(5)*T4(12)+V2(6)*T4(16)
     $ -V2(3)*T4(4))+(P1(2)*(V2(4)*T4(9)+V2(5)*T4(13)+V2(6)*T4(17)
     $ -V2(3)*T4(5))+P1(3)*(V2(4)*T4(10)+V2(5)*T4(14)+V2(6)*T4(18)
     $ -V2(3)*T4(6)))))
      TMP60 = (P3(0)*(-1D0)*(V2(4)*T4(4)+V2(5)*T4(5)+V2(6)*T4(6)-V2(3)
     $ *T4(3))+(P3(1)*(V2(4)*T4(8)+V2(5)*T4(9)+V2(6)*T4(10)-V2(3)*T4(7)
     $ )+(P3(2)*(V2(4)*T4(12)+V2(5)*T4(13)+V2(6)*T4(14)-V2(3)*T4(11))
     $ +P3(3)*(V2(4)*T4(16)+V2(5)*T4(17)+V2(6)*T4(18)-V2(3)*T4(15)))))
      TMP62 = (P1(0)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(P1(1)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(P1(2)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+P1(3)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP63 = (P2(0)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(P2(1)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(P2(2)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+P2(3)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP65 = (V2(3)*(-1D0)*(V3(4)*T4(7)+V3(5)*T4(11)+V3(6)*T4(15)
     $ -V3(3)*T4(3))+(V2(4)*(V3(4)*T4(8)+V3(5)*T4(12)+V3(6)*T4(16)
     $ -V3(3)*T4(4))+(V2(5)*(V3(4)*T4(9)+V3(5)*T4(13)+V3(6)*T4(17)
     $ -V3(3)*T4(5))+V2(6)*(V3(4)*T4(10)+V3(5)*T4(14)+V3(6)*T4(18)
     $ -V3(3)*T4(6)))))
      TMP66 = (P1(0)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(P1(1)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(P1(2)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +P1(3)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP67 = (P2(0)*(-1D0)*(V3(4)*T4(4)+V3(5)*T4(5)+V3(6)*T4(6)-V3(3)
     $ *T4(3))+(P2(1)*(V3(4)*T4(8)+V3(5)*T4(9)+V3(6)*T4(10)-V3(3)*T4(7)
     $ )+(P2(2)*(V3(4)*T4(12)+V3(5)*T4(13)+V3(6)*T4(14)-V3(3)*T4(11))
     $ +P2(3)*(V3(4)*T4(16)+V3(5)*T4(17)+V3(6)*T4(18)-V3(3)*T4(15)))))
      TMP59 = (P1(0)*(-1D0)*(V2(4)*T4(4)+V2(5)*T4(5)+V2(6)*T4(6)-V2(3)
     $ *T4(3))+(P1(1)*(V2(4)*T4(8)+V2(5)*T4(9)+V2(6)*T4(10)-V2(3)*T4(7)
     $ )+(P1(2)*(V2(4)*T4(12)+V2(5)*T4(13)+V2(6)*T4(14)-V2(3)*T4(11))
     $ +P1(3)*(V2(4)*T4(16)+V2(5)*T4(17)+V2(6)*T4(18)-V2(3)*T4(15)))))
      TMP9 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP8 = (V3(3)*P2(0)-V3(4)*P2(1)-V3(5)*P2(2)-V3(6)*P2(3))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP7 = (V3(3)*P1(0)-V3(4)*P1(1)-V3(5)*P1(2)-V3(6)*P1(3))
      TMP10 = (V2(3)*P3(0)-V2(4)*P3(1)-V2(5)*P3(2)-V2(6)*P3(3))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(-2D0 * CI)*(TMP5*(P2(1)*(-1D0/2D0)*(T4(4)+T4(7))
     $ +(P2(2)*(-1D0/2D0)*(T4(5)+T4(11))+(P2(3)*(-1D0/2D0)*(T4(6)
     $ +T4(15))+(P3(1)*1D0/2D0*(T4(4)+T4(7))+(P3(2)*1D0/2D0*(T4(5)
     $ +T4(11))+(P3(3)*1D0/2D0*(T4(6)+T4(15))+T4(3)*(P2(0)-P3(0))))))))
     $ +(TMP10*(V3(4)*(-1D0/2D0)*(T4(7)+T4(4))+(V3(5)*(-1D0/2D0)
     $ *(T4(11)+T4(5))+(V3(6)*(-1D0/2D0)*(T4(15)+T4(6))+V3(3)*T4(3))))
     $ +(TMP7*(V2(4)*(-1D0/2D0)*(T4(7)+T4(4))+(V2(5)*(-1D0/2D0)*(T4(11)
     $ +T4(5))+(V2(6)*(-1D0/2D0)*(T4(15)+T4(6))+V2(3)*T4(3))))+(TMP8
     $ *(V2(4)*1D0/2D0*(T4(7)+T4(4))+(V2(5)*1D0/2D0*(T4(11)+T4(5))
     $ +(V2(6)*1D0/2D0*(T4(15)+T4(6))-V2(3)*T4(3))))+(TMP9*(V3(4)*1D0
     $ /2D0*(T4(7)+T4(4))+(V3(5)*1D0/2D0*(T4(11)+T4(5))+(V3(6)*1D0/2D0
     $ *(T4(15)+T4(6))-V3(3)*T4(3))))+(V2(3)*1D0/2D0*(TMP62+TMP66
     $ -TMP63-TMP67)+(V3(3)*1D0/2D0*(TMP57+TMP60-TMP56-TMP59)+(-1D0
     $ /2D0*(P3(0)*(TMP65+TMP69))+P2(0)*1D0/2D0*(TMP65+TMP69)))))))))
      V1(4)= DENOM*(-2D0 * CI)*(TMP5*(P2(0)*1D0/2D0*(T4(7)+T4(4))
     $ +(P2(2)*(-1D0/2D0)*(T4(9)+T4(12))+(P2(3)*(-1D0/2D0)*(T4(10)
     $ +T4(16))+(P3(0)*(-1D0/2D0)*(T4(7)+T4(4))+(P3(2)*1D0/2D0*(T4(9)
     $ +T4(12))+(P3(3)*1D0/2D0*(T4(10)+T4(16))+T4(8)*(P3(1)-P2(1)))))))
     $ )+(TMP10*(V3(3)*1D0/2D0*(T4(4)+T4(7))+(V3(5)*(-1D0/2D0)*(T4(12)
     $ +T4(9))+(V3(6)*(-1D0/2D0)*(T4(16)+T4(10))-V3(4)*T4(8))))+(TMP7
     $ *(V2(3)*1D0/2D0*(T4(4)+T4(7))+(V2(5)*(-1D0/2D0)*(T4(12)+T4(9))
     $ +(V2(6)*(-1D0/2D0)*(T4(16)+T4(10))-V2(4)*T4(8))))+(TMP8*(V2(3)
     $ *(-1D0/2D0)*(T4(4)+T4(7))+(V2(5)*1D0/2D0*(T4(12)+T4(9))+(V2(6)
     $ *1D0/2D0*(T4(16)+T4(10))+V2(4)*T4(8))))+(TMP9*(V3(3)*(-1D0/2D0)
     $ *(T4(4)+T4(7))+(V3(5)*1D0/2D0*(T4(12)+T4(9))+(V3(6)*1D0/2D0
     $ *(T4(16)+T4(10))+V3(4)*T4(8))))+(V2(4)*1D0/2D0*(TMP62+TMP66
     $ -TMP63-TMP67)+(V3(4)*1D0/2D0*(TMP57+TMP60-TMP56-TMP59)+(-1D0
     $ /2D0*(P3(1)*(TMP65+TMP69))+P2(1)*1D0/2D0*(TMP65+TMP69)))))))))
      V1(5)= DENOM*(-2D0 * CI)*(TMP5*(P2(0)*1D0/2D0*(T4(11)+T4(5))
     $ +(P2(1)*(-1D0/2D0)*(T4(12)+T4(9))+(P2(3)*(-1D0/2D0)*(T4(14)
     $ +T4(17))+(P3(0)*(-1D0/2D0)*(T4(11)+T4(5))+(P3(1)*1D0/2D0*(T4(12)
     $ +T4(9))+(P3(3)*1D0/2D0*(T4(14)+T4(17))+T4(13)*(P3(2)-P2(2)))))))
     $ )+(TMP10*(V3(3)*1D0/2D0*(T4(5)+T4(11))+(V3(4)*(-1D0/2D0)*(T4(9)
     $ +T4(12))+(V3(6)*(-1D0/2D0)*(T4(17)+T4(14))-V3(5)*T4(13))))
     $ +(TMP7*(V2(3)*1D0/2D0*(T4(5)+T4(11))+(V2(4)*(-1D0/2D0)*(T4(9)
     $ +T4(12))+(V2(6)*(-1D0/2D0)*(T4(17)+T4(14))-V2(5)*T4(13))))
     $ +(TMP8*(V2(3)*(-1D0/2D0)*(T4(5)+T4(11))+(V2(4)*1D0/2D0*(T4(9)
     $ +T4(12))+(V2(6)*1D0/2D0*(T4(17)+T4(14))+V2(5)*T4(13))))+(TMP9
     $ *(V3(3)*(-1D0/2D0)*(T4(5)+T4(11))+(V3(4)*1D0/2D0*(T4(9)+T4(12))
     $ +(V3(6)*1D0/2D0*(T4(17)+T4(14))+V3(5)*T4(13))))+(V2(5)*1D0/2D0
     $ *(TMP62+TMP66-TMP63-TMP67)+(V3(5)*1D0/2D0*(TMP57+TMP60-TMP56
     $ -TMP59)+(-1D0/2D0*(P3(2)*(TMP65+TMP69))+P2(2)*1D0/2D0*(TMP65
     $ +TMP69)))))))))
      V1(6)= DENOM*(-2D0 * CI)*(TMP5*(P2(0)*1D0/2D0*(T4(15)+T4(6))
     $ +(P2(1)*(-1D0/2D0)*(T4(16)+T4(10))+(P2(2)*(-1D0/2D0)*(T4(17)
     $ +T4(14))+(P3(0)*(-1D0/2D0)*(T4(15)+T4(6))+(P3(1)*1D0/2D0*(T4(16)
     $ +T4(10))+(P3(2)*1D0/2D0*(T4(17)+T4(14))+T4(18)*(P3(3)-P2(3))))))
     $ ))+(TMP10*(V3(3)*1D0/2D0*(T4(6)+T4(15))+(V3(4)*(-1D0/2D0)
     $ *(T4(10)+T4(16))+(V3(5)*(-1D0/2D0)*(T4(14)+T4(17))-V3(6)*T4(18))
     $ ))+(TMP7*(V2(3)*1D0/2D0*(T4(6)+T4(15))+(V2(4)*(-1D0/2D0)*(T4(10)
     $ +T4(16))+(V2(5)*(-1D0/2D0)*(T4(14)+T4(17))-V2(6)*T4(18))))
     $ +(TMP8*(V2(3)*(-1D0/2D0)*(T4(6)+T4(15))+(V2(4)*1D0/2D0*(T4(10)
     $ +T4(16))+(V2(5)*1D0/2D0*(T4(14)+T4(17))+V2(6)*T4(18))))+(TMP9
     $ *(V3(3)*(-1D0/2D0)*(T4(6)+T4(15))+(V3(4)*1D0/2D0*(T4(10)+T4(16))
     $ +(V3(5)*1D0/2D0*(T4(14)+T4(17))+V3(6)*T4(18))))+(V2(6)*1D0/2D0
     $ *(TMP62+TMP66-TMP63-TMP67)+(V3(6)*1D0/2D0*(TMP57+TMP60-TMP56
     $ -TMP59)+(-1D0/2D0*(P3(3)*(TMP65+TMP69))+P2(3)*1D0/2D0*(TMP65
     $ +TMP69)))))))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)
C     
      SUBROUTINE VVVV6_0(V1, V2, V3, V4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP1
      COMPLEX*16 TMP0
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP5
      COMPLEX*16 TMP4
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP4 = (V4(3)*V1(3)-V4(4)*V1(4)-V4(5)*V1(5)-V4(6)*V1(6))
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      VERTEX = COUP*(-CI*(TMP4*TMP5)+CI*(TMP0*TMP1))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)
C     
      SUBROUTINE VVVV6P0_1(V2, V3, V4, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 DENOM
      COMPLEX*16 V3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      COMPLEX*16 TMP0
      REAL*8 W1
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP5
      COMPLEX*16 COUP
      COMPLEX*16 V1(6)
      V1(1) = +V2(1)+V3(1)+V4(1)
      V1(2) = +V2(2)+V3(2)+V4(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(-CI*(V4(3)*TMP5)+CI*(V3(3)*TMP0))
      V1(4)= DENOM*(-CI*(V4(4)*TMP5)+CI*(V3(4)*TMP0))
      V1(5)= DENOM*(-CI*(V4(5)*TMP5)+CI*(V3(5)*TMP0))
      V1(6)= DENOM*(-CI*(V4(6)*TMP5)+CI*(V3(6)*TMP0))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV8_0(V1, V2, V3, V4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 V3(*)
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP5
      COMPLEX*16 TMP4
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP4 = (V4(3)*V1(3)-V4(4)*V1(4)-V4(5)*V1(5)-V4(6)*V1(6))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      VERTEX = COUP*(-CI*(TMP4*TMP5)+CI*(TMP2*TMP3))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV8P0_1(V2, V3, V4, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 DENOM
      COMPLEX*16 V3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      COMPLEX*16 TMP2
      REAL*8 W1
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP5
      COMPLEX*16 COUP
      COMPLEX*16 V1(6)
      V1(1) = +V2(1)+V3(1)+V4(1)
      V1(2) = +V2(2)+V3(2)+V4(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(-CI*(V4(3)*TMP5)+CI*(V2(3)*TMP2))
      V1(4)= DENOM*(-CI*(V4(4)*TMP5)+CI*(V2(4)*TMP2))
      V1(5)= DENOM*(-CI*(V4(5)*TMP5)+CI*(V2(5)*TMP2))
      V1(6)= DENOM*(-CI*(V4(6)*TMP5)+CI*(V2(6)*TMP2))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV9_0(V1, V2, V3, V4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP1
      COMPLEX*16 TMP0
      COMPLEX*16 V4(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      VERTEX = COUP*(-CI*(TMP0*TMP1)+CI*(TMP2*TMP3))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV9P0_1(V2, V3, V4, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 V3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      COMPLEX*16 TMP0
      REAL*8 W1
      COMPLEX*16 V4(*)
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      COMPLEX*16 V1(6)
      V1(1) = +V2(1)+V3(1)+V4(1)
      V1(2) = +V2(2)+V3(2)+V4(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*(-CI*(V3(3)*TMP0)+CI*(V2(3)*TMP2))
      V1(4)= DENOM*(-CI*(V3(4)*TMP0)+CI*(V2(4)*TMP2))
      V1(5)= DENOM*(-CI*(V3(5)*TMP0)+CI*(V2(5)*TMP2))
      V1(6)= DENOM*(-CI*(V3(6)*TMP0)+CI*(V2(6)*TMP2))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,2005)*Metric(2,4)*Metric(3,1005) - Metric(1,4)*Metric(2,
C     2005)*Metric(3,1005) + Metric(1,1005)*Metric(2,4)*Metric(3,2005)
C      - Metric(1,4)*Metric(2,1005)*Metric(3,2005) -
C      Metric(1,2005)*Metric(2,3)*Metric(4,1005) + Metric(1,3)*Metric(2
C     ,2005)*Metric(4,1005) - Metric(1,1005)*Metric(2,3)*Metric(4,2005)
C      + Metric(1,3)*Metric(2,1005)*Metric(4,2005)
C     
      SUBROUTINE VVVVT6_0(V1, V2, V3, V4, T5, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP1
      COMPLEX*16 COUP
      COMPLEX*16 TMP44
      COMPLEX*16 TMP0
      COMPLEX*16 TMP45
      COMPLEX*16 T5(*)
      COMPLEX*16 TMP46
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP47
      COMPLEX*16 TMP5
      COMPLEX*16 TMP40
      COMPLEX*16 TMP4
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP41
      COMPLEX*16 TMP42
      COMPLEX*16 TMP43
      COMPLEX*16 V1(*)
      TMP42 = (V1(3)*(-1D0)*(V3(4)*T5(4)+V3(5)*T5(5)+V3(6)*T5(6)-V3(3)
     $ *T5(3))+(V1(4)*(V3(4)*T5(8)+V3(5)*T5(9)+V3(6)*T5(10)-V3(3)*T5(7)
     $ )+(V1(5)*(V3(4)*T5(12)+V3(5)*T5(13)+V3(6)*T5(14)-V3(3)*T5(11))
     $ +V1(6)*(V3(4)*T5(16)+V3(5)*T5(17)+V3(6)*T5(18)-V3(3)*T5(15)))))
      TMP43 = (V2(3)*(-1D0)*(V3(4)*T5(4)+V3(5)*T5(5)+V3(6)*T5(6)-V3(3)
     $ *T5(3))+(V2(4)*(V3(4)*T5(8)+V3(5)*T5(9)+V3(6)*T5(10)-V3(3)*T5(7)
     $ )+(V2(5)*(V3(4)*T5(12)+V3(5)*T5(13)+V3(6)*T5(14)-V3(3)*T5(11))
     $ +V2(6)*(V3(4)*T5(16)+V3(5)*T5(17)+V3(6)*T5(18)-V3(3)*T5(15)))))
      TMP40 = (V1(3)*(-1D0)*(V3(4)*T5(7)+V3(5)*T5(11)+V3(6)*T5(15)
     $ -V3(3)*T5(3))+(V1(4)*(V3(4)*T5(8)+V3(5)*T5(12)+V3(6)*T5(16)
     $ -V3(3)*T5(4))+(V1(5)*(V3(4)*T5(9)+V3(5)*T5(13)+V3(6)*T5(17)
     $ -V3(3)*T5(5))+V1(6)*(V3(4)*T5(10)+V3(5)*T5(14)+V3(6)*T5(18)
     $ -V3(3)*T5(6)))))
      TMP41 = (V2(3)*(-1D0)*(V3(4)*T5(7)+V3(5)*T5(11)+V3(6)*T5(15)
     $ -V3(3)*T5(3))+(V2(4)*(V3(4)*T5(8)+V3(5)*T5(12)+V3(6)*T5(16)
     $ -V3(3)*T5(4))+(V2(5)*(V3(4)*T5(9)+V3(5)*T5(13)+V3(6)*T5(17)
     $ -V3(3)*T5(5))+V2(6)*(V3(4)*T5(10)+V3(5)*T5(14)+V3(6)*T5(18)
     $ -V3(3)*T5(6)))))
      TMP46 = (V1(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V1(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V1(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V1(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP47 = (V2(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V2(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V2(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V2(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP44 = (V1(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V1(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V1(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V1(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP45 = (V2(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V2(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V2(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V2(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP4 = (V4(3)*V1(3)-V4(4)*V1(4)-V4(5)*V1(5)-V4(6)*V1(6))
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      VERTEX = COUP*(TMP0*(-1D0)*(+CI*(TMP40+TMP42))+(TMP1*(-1D0)*(+CI
     $ *(TMP45+TMP47))+(TMP4*(+CI*(TMP41+TMP43))+TMP5*(+CI*(TMP44
     $ +TMP46)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,2005)*Metric(2,1005)*Metric(3,4) + Metric(1,1005)*Metric
C     (2,2005)*Metric(3,4) - Metric(1,4)*Metric(2,2005)*Metric(3,1005)
C      - Metric(1,4)*Metric(2,1005)*Metric(3,2005) -
C      Metric(1,2005)*Metric(2,3)*Metric(4,1005) + Metric(1,2)*Metric(3
C     ,2005)*Metric(4,1005) - Metric(1,1005)*Metric(2,3)*Metric(4,2005)
C      + Metric(1,2)*Metric(3,1005)*Metric(4,2005)
C     
      SUBROUTINE VVVVT7_0(V1, V2, V3, V4, T5, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 TMP50
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP44
      COMPLEX*16 COUP
      COMPLEX*16 T5(*)
      COMPLEX*16 TMP46
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP49
      COMPLEX*16 TMP5
      COMPLEX*16 TMP4
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP41
      COMPLEX*16 TMP51
      COMPLEX*16 TMP48
      COMPLEX*16 TMP43
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      TMP51 = (V3(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V3(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V3(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V3(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP43 = (V2(3)*(-1D0)*(V3(4)*T5(4)+V3(5)*T5(5)+V3(6)*T5(6)-V3(3)
     $ *T5(3))+(V2(4)*(V3(4)*T5(8)+V3(5)*T5(9)+V3(6)*T5(10)-V3(3)*T5(7)
     $ )+(V2(5)*(V3(4)*T5(12)+V3(5)*T5(13)+V3(6)*T5(14)-V3(3)*T5(11))
     $ +V2(6)*(V3(4)*T5(16)+V3(5)*T5(17)+V3(6)*T5(18)-V3(3)*T5(15)))))
      TMP41 = (V2(3)*(-1D0)*(V3(4)*T5(7)+V3(5)*T5(11)+V3(6)*T5(15)
     $ -V3(3)*T5(3))+(V2(4)*(V3(4)*T5(8)+V3(5)*T5(12)+V3(6)*T5(16)
     $ -V3(3)*T5(4))+(V2(5)*(V3(4)*T5(9)+V3(5)*T5(13)+V3(6)*T5(17)
     $ -V3(3)*T5(5))+V2(6)*(V3(4)*T5(10)+V3(5)*T5(14)+V3(6)*T5(18)
     $ -V3(3)*T5(6)))))
      TMP46 = (V1(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V1(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V1(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V1(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP44 = (V1(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V1(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V1(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V1(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP50 = (V3(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V3(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V3(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V3(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP48 = (V1(3)*(-1D0)*(V2(4)*T5(7)+V2(5)*T5(11)+V2(6)*T5(15)
     $ -V2(3)*T5(3))+(V1(4)*(V2(4)*T5(8)+V2(5)*T5(12)+V2(6)*T5(16)
     $ -V2(3)*T5(4))+(V1(5)*(V2(4)*T5(9)+V2(5)*T5(13)+V2(6)*T5(17)
     $ -V2(3)*T5(5))+V1(6)*(V2(4)*T5(10)+V2(5)*T5(14)+V2(6)*T5(18)
     $ -V2(3)*T5(6)))))
      TMP49 = (V1(3)*(-1D0)*(V2(4)*T5(4)+V2(5)*T5(5)+V2(6)*T5(6)-V2(3)
     $ *T5(3))+(V1(4)*(V2(4)*T5(8)+V2(5)*T5(9)+V2(6)*T5(10)-V2(3)*T5(7)
     $ )+(V1(5)*(V2(4)*T5(12)+V2(5)*T5(13)+V2(6)*T5(14)-V2(3)*T5(11))
     $ +V1(6)*(V2(4)*T5(16)+V2(5)*T5(17)+V2(6)*T5(18)-V2(3)*T5(15)))))
      TMP5 = (V2(3)*V3(3)-V2(4)*V3(4)-V2(5)*V3(5)-V2(6)*V3(6))
      TMP4 = (V4(3)*V1(3)-V4(4)*V1(4)-V4(5)*V1(5)-V4(6)*V1(6))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      VERTEX = COUP*(TMP2*(-1D0)*(+CI*(TMP48+TMP49))+(TMP3*(-1D0)*(+CI
     $ *(TMP50+TMP51))+(TMP4*(+CI*(TMP41+TMP43))+TMP5*(+CI*(TMP44
     $ +TMP46)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,2005)*Metric(2,1005)*Metric(3,4) + Metric(1,1005)*Metric
C     (2,2005)*Metric(3,4) - Metric(1,2005)*Metric(2,4)*Metric(3,1005)
C      - Metric(1,1005)*Metric(2,4)*Metric(3,2005) -
C      Metric(1,3)*Metric(2,2005)*Metric(4,1005) + Metric(1,2)*Metric(3
C     ,2005)*Metric(4,1005) - Metric(1,3)*Metric(2,1005)*Metric(4,2005)
C      + Metric(1,2)*Metric(3,1005)*Metric(4,2005)
C     
      SUBROUTINE VVVVT8_0(V1, V2, V3, V4, T5, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP1
      COMPLEX*16 TMP0
      COMPLEX*16 TMP45
      COMPLEX*16 T5(*)
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP47
      COMPLEX*16 COUP
      COMPLEX*16 TMP51
      COMPLEX*16 TMP40
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP42
      COMPLEX*16 TMP49
      COMPLEX*16 TMP48
      COMPLEX*16 TMP50
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      TMP42 = (V1(3)*(-1D0)*(V3(4)*T5(4)+V3(5)*T5(5)+V3(6)*T5(6)-V3(3)
     $ *T5(3))+(V1(4)*(V3(4)*T5(8)+V3(5)*T5(9)+V3(6)*T5(10)-V3(3)*T5(7)
     $ )+(V1(5)*(V3(4)*T5(12)+V3(5)*T5(13)+V3(6)*T5(14)-V3(3)*T5(11))
     $ +V1(6)*(V3(4)*T5(16)+V3(5)*T5(17)+V3(6)*T5(18)-V3(3)*T5(15)))))
      TMP50 = (V3(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V3(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V3(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V3(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP40 = (V1(3)*(-1D0)*(V3(4)*T5(7)+V3(5)*T5(11)+V3(6)*T5(15)
     $ -V3(3)*T5(3))+(V1(4)*(V3(4)*T5(8)+V3(5)*T5(12)+V3(6)*T5(16)
     $ -V3(3)*T5(4))+(V1(5)*(V3(4)*T5(9)+V3(5)*T5(13)+V3(6)*T5(17)
     $ -V3(3)*T5(5))+V1(6)*(V3(4)*T5(10)+V3(5)*T5(14)+V3(6)*T5(18)
     $ -V3(3)*T5(6)))))
      TMP47 = (V2(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V2(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V2(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V2(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP45 = (V2(3)*(-1D0)*(V4(4)*T5(7)+V4(5)*T5(11)+V4(6)*T5(15)
     $ -V4(3)*T5(3))+(V2(4)*(V4(4)*T5(8)+V4(5)*T5(12)+V4(6)*T5(16)
     $ -V4(3)*T5(4))+(V2(5)*(V4(4)*T5(9)+V4(5)*T5(13)+V4(6)*T5(17)
     $ -V4(3)*T5(5))+V2(6)*(V4(4)*T5(10)+V4(5)*T5(14)+V4(6)*T5(18)
     $ -V4(3)*T5(6)))))
      TMP48 = (V1(3)*(-1D0)*(V2(4)*T5(7)+V2(5)*T5(11)+V2(6)*T5(15)
     $ -V2(3)*T5(3))+(V1(4)*(V2(4)*T5(8)+V2(5)*T5(12)+V2(6)*T5(16)
     $ -V2(3)*T5(4))+(V1(5)*(V2(4)*T5(9)+V2(5)*T5(13)+V2(6)*T5(17)
     $ -V2(3)*T5(5))+V1(6)*(V2(4)*T5(10)+V2(5)*T5(14)+V2(6)*T5(18)
     $ -V2(3)*T5(6)))))
      TMP49 = (V1(3)*(-1D0)*(V2(4)*T5(4)+V2(5)*T5(5)+V2(6)*T5(6)-V2(3)
     $ *T5(3))+(V1(4)*(V2(4)*T5(8)+V2(5)*T5(9)+V2(6)*T5(10)-V2(3)*T5(7)
     $ )+(V1(5)*(V2(4)*T5(12)+V2(5)*T5(13)+V2(6)*T5(14)-V2(3)*T5(11))
     $ +V1(6)*(V2(4)*T5(16)+V2(5)*T5(17)+V2(6)*T5(18)-V2(3)*T5(15)))))
      TMP51 = (V3(3)*(-1D0)*(V4(4)*T5(4)+V4(5)*T5(5)+V4(6)*T5(6)-V4(3)
     $ *T5(3))+(V3(4)*(V4(4)*T5(8)+V4(5)*T5(9)+V4(6)*T5(10)-V4(3)*T5(7)
     $ )+(V3(5)*(V4(4)*T5(12)+V4(5)*T5(13)+V4(6)*T5(14)-V4(3)*T5(11))
     $ +V3(6)*(V4(4)*T5(16)+V4(5)*T5(17)+V4(6)*T5(18)-V4(3)*T5(15)))))
      TMP1 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP0 = (V4(3)*V2(3)-V4(4)*V2(4)-V4(5)*V2(5)-V4(6)*V2(6))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP2 = (V4(3)*V3(3)-V4(4)*V3(4)-V4(5)*V3(5)-V4(6)*V3(6))
      VERTEX = COUP*(TMP0*(+CI*(TMP40+TMP42))+(TMP1*(+CI*(TMP45+TMP47))
     $ +(TMP2*(-1D0)*(+CI*(TMP48+TMP49))-TMP3*(+CI*(TMP50+TMP51)))))
      END

C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(2003,1)*Gamma(1003,2,-1)*ProjP(-1,1) - P(2003,2)*Gamma(1003,2,-
C     1)*ProjP(-1,1) + P(1003,1)*Gamma(2003,2,-1)*ProjP(-1,1) -
C      P(1003,2)*Gamma(2003,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFT5_0(F1, F2, T3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP15
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 F1(*)
      COMPLEX*16 TMP16
      COMPLEX*16 F2(*)
      COMPLEX*16 TMP14
      COMPLEX*16 VERTEX
      COMPLEX*16 T3(*)
      COMPLEX*16 COUP
      COMPLEX*16 TMP13
      P1(0) = DBLE(F1(1))
      P1(1) = DBLE(F1(2))
      P1(2) = DIMAG(F1(2))
      P1(3) = DIMAG(F1(1))
      P2(0) = DBLE(F2(1))
      P2(1) = DBLE(F2(2))
      P2(2) = DIMAG(F2(2))
      P2(3) = DIMAG(F2(1))
      TMP15 = (F1(5)*(F2(3)*(P1(0)*(T3(3)-T3(6))+(P1(1)*(T3(10)-T3(7))
     $ +(P1(2)*(T3(14)-T3(11))+P1(3)*(T3(18)-T3(15)))))+F2(4)*(P1(0)*(
     $ -1D0)*(T3(4)+CI*(T3(5)))+(P1(1)*(T3(8)+CI*(T3(9)))+(P1(2)
     $ *(T3(12)+CI*(T3(13)))+P1(3)*(T3(16)+CI*(T3(17)))))))+F1(6)
     $ *(F2(3)*(P1(0)*(+CI*(T3(5))-T3(4))+(P1(1)*(T3(8)-CI*(T3(9)))
     $ +(P1(2)*(T3(12)-CI*(T3(13)))+P1(3)*(T3(16)-CI*(T3(17))))))+F2(4)
     $ *(P1(0)*(T3(3)+T3(6))+(P1(1)*(-1D0)*(T3(7)+T3(10))+(P1(2)*(-1D0)
     $ *(T3(11)+T3(14))-P1(3)*(T3(15)+T3(18)))))))
      TMP14 = (F1(5)*(F2(3)*(P2(0)*(T3(3)-T3(15))+(P2(1)*(T3(16)-T3(4))
     $ +(P2(2)*(T3(17)-T3(5))+P2(3)*(T3(18)-T3(6)))))+F2(4)*(P2(0)*(
     $ -1D0)*(T3(7)+CI*(T3(11)))+(P2(1)*(T3(8)+CI*(T3(12)))+(P2(2)
     $ *(T3(9)+CI*(T3(13)))+P2(3)*(T3(10)+CI*(T3(14)))))))+F1(6)*(F2(3)
     $ *(P2(0)*(+CI*(T3(11))-T3(7))+(P2(1)*(T3(8)-CI*(T3(12)))+(P2(2)
     $ *(T3(9)-CI*(T3(13)))+P2(3)*(T3(10)-CI*(T3(14))))))+F2(4)*(P2(0)
     $ *(T3(3)+T3(15))+(P2(1)*(-1D0)*(T3(4)+T3(16))+(P2(2)*(-1D0)
     $ *(T3(5)+T3(17))-P2(3)*(T3(6)+T3(18)))))))
      TMP16 = (F1(5)*(F2(3)*(P2(0)*(T3(3)-T3(6))+(P2(1)*(T3(10)-T3(7))
     $ +(P2(2)*(T3(14)-T3(11))+P2(3)*(T3(18)-T3(15)))))+F2(4)*(P2(0)*(
     $ -1D0)*(T3(4)+CI*(T3(5)))+(P2(1)*(T3(8)+CI*(T3(9)))+(P2(2)
     $ *(T3(12)+CI*(T3(13)))+P2(3)*(T3(16)+CI*(T3(17)))))))+F1(6)
     $ *(F2(3)*(P2(0)*(+CI*(T3(5))-T3(4))+(P2(1)*(T3(8)-CI*(T3(9)))
     $ +(P2(2)*(T3(12)-CI*(T3(13)))+P2(3)*(T3(16)-CI*(T3(17))))))+F2(4)
     $ *(P2(0)*(T3(3)+T3(6))+(P2(1)*(-1D0)*(T3(7)+T3(10))+(P2(2)*(-1D0)
     $ *(T3(11)+T3(14))-P2(3)*(T3(15)+T3(18)))))))
      TMP13 = (F1(5)*(F2(3)*(P1(0)*(T3(3)-T3(15))+(P1(1)*(T3(16)-T3(4))
     $ +(P1(2)*(T3(17)-T3(5))+P1(3)*(T3(18)-T3(6)))))+F2(4)*(P1(0)*(
     $ -1D0)*(T3(7)+CI*(T3(11)))+(P1(1)*(T3(8)+CI*(T3(12)))+(P1(2)
     $ *(T3(9)+CI*(T3(13)))+P1(3)*(T3(10)+CI*(T3(14)))))))+F1(6)*(F2(3)
     $ *(P1(0)*(+CI*(T3(11))-T3(7))+(P1(1)*(T3(8)-CI*(T3(12)))+(P1(2)
     $ *(T3(9)-CI*(T3(13)))+P1(3)*(T3(10)-CI*(T3(14))))))+F2(4)*(P1(0)
     $ *(T3(3)+T3(15))+(P1(1)*(-1D0)*(T3(4)+T3(16))+(P1(2)*(-1D0)
     $ *(T3(5)+T3(17))-P1(3)*(T3(6)+T3(18)))))))
      VERTEX = COUP*(-CI*(TMP13+TMP15)+CI*(TMP14+TMP16))
      END




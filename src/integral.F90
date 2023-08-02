! ____________________________________________________________________________
! __________________________ INTEGRAL ROUTINES _______________________________
! ____________________________________________________________________________

      subroutine replsn ()
      use molkst_C, ONLY: numat

      use reimers_C, only: ind, f0, x, y, z, zcore, alpha, gamma, beta, &
          fg, vnn, nrep, ifroti, mic2, mic1, mf0ss, mf0sd, mf0dd, mg1sp, &
          mf2pp, mg2sd, mg1pd, mf2pd, mg3pd, mf2dd, mf4dd, mr1sppd, &
          mr2sddd, mr2sdpp, ntrmet, itrmet, nirreg, natt, nte, natm, &
          ibf, nbf, matind, na
      USE chanel_C, only : iw

      implicit none
      double precision  enr, enucrep, fracsi, fracsj, gav, gamma1, &
                        rij, sqrt3, vnn1, vnnz, xx, zero
      integer           i, j, k, &
                        i1, i2, ia, ib, ibk, &
                        j1, j2, ja, jb, jbk,&
                        kb, mte, ndump, nham, nm

      data zero /0.D0/
      data mte  /8/
!      r(i,j)= sqrt ( (x(j)-x(i))**2 + (y(j)-y(i))**2 + (z(j)-z(i))**2 )

      na = numat

 1000 format (/' ERROR: Max number of different transition metals=',i4)
 1004 format (/' Nuclear repulsion energy=', f14.5,' eV',&
     &	       ', using ZINDO gamma(ss) formula=',f12.5)
 1007 format (1h1,10x,a80,5x,2(i2,'/'),i2/)
 2000 format (' POINT CHARGE: atoms',2i4,' at r=',f8.4,' gamma=',f8.4)

!     **** this evaluates the electron repulsion integrals and also ****
!     **** the nuclear repulsion energy ****
      nrep = 1
      vnn= zero
      vnnz= zero
      nham = 2

      do i= 1,na
        ia= natm(i)
        ib= ibf(i)

        do j= i,1,-1
          ja= natm(j)
          jb= ibf(j)
          rij= sqrt ( (x(j)-x(i))**2 + (y(j)-y(i))**2 + (z(j)-z(i))**2 )
!          write (6,*) 'dist', i, j, x(i), x(j), y(i), y(j), z(i), z(j), rij

!	  **** s-s, s-p, or p-p interactions ******
          f0(1,1)= gamma1 (fg(mf0ss,ia),fg(mf0ss,ja),ia,ja,rij)
          gav= f0(1,1)
          fracsi= 1.D0
          fracsj= 1.D0
          if (nbf(i).eq.9) then
!	    **** i-j is d-sp interaction ****
            if (i.eq.j) then
!	      *** gamma(sd) is only used if one-centre integral ****
              f0(2,1)= gamma1 (fg(mf0sd,ia),fg(mf0sd,ja),ia,ja,rij)
            else
!	      *** Zindo does this, but seems crazy to me ****
              f0(2,1)= gamma1 (fg(mf0dd,ia),fg(mf0ss,ja),ia,ja,rij)
            end if
            fracsi= (fg(mic1,ia) + fg(mic2,ia)*2) / zcore(i)

            if (nbf(j).eq.9) then
!	      **** d-d interactions ****
              f0(2,2)= gamma1 (fg(mf0dd,ia),fg(mf0dd,ja),ia,ja,rij)
              fracsj= (fg(mic1,ja) + fg(mic2,ja)*2) / zcore(j)
              gav=  fracsi*fracsj*f0(1,1)+ fracsi*(1-fracsj)*f0(1,2) +&
     &          (1-fracsi)*fracsj*f0(2,1)+ (1-fracsi)*(1-fracsj)*f0(2,2)
            else
              gav= fracsi*f0(1,1) + (1-fracsi)*f0(2,1)
            end if
          end if
          if (nbf(j).eq.9) then
!	    **** i-j is sp-d interaction ****
            if (i.eq.j) then
!	      *** gamma(sd) is only used if one-centre integral ****
              f0(1,2)= gamma1 (fg(mf0sd,ia),fg(mf0sd,ja),ia,ja,rij)
            else
!	      *** Zindo does this, but seems crazy to me ****
              f0(1,2)= gamma1 (fg(mf0ss,ia),fg(mf0dd,ja),ia,ja,rij)
            end if
            if (nbf(i).lt.9) then
              fracsj= (fg(mic1,ja) + fg(mic2,ja)*2) / zcore(j)
              gav= fracsj*f0(1,1) + (1-fracsj)*f0(1,2)
            end if
          end if
!	  **** fill in atom-atom block, f0(ss), f0(sp), .... f0(dd) ****

          do i1= 0,nbf(i)-1
            i2= ib + i1
            do j1= 0,nbf(j)-1
              j2= jb + j1
              if (i2.ge.j2) then
                gamma(i2,j2)= f0(ind(i1),ind(j1))
              end if
            end do
          end do
!	  **** nuclear repulsion energy, gav= weighted s,p,d gammas ****

          if (i.eq.j) then
!	    **** calculation of core matrix elements for this atom ****
            call uform (i,fg(1,ia),f0,zcore(i))
          else

!	    **** nuclear repulsion energy ****
!            if (ia.le.9 .and. ja.le.9) alp= alpha(ia,ja)
            enr= enucrep (rij,gav)
!            vnn1= zcore(i) * zcore(j) * enr
            vnn1= zcore(i) * zcore(j) * gav
            vnn=  vnn  + vnn1
!            write (6,*) 'enr', enr, rij, gav
!            write (6,*) 'vnn', vnn1, vnn
!            write (6,*) 'zcr', zcore(i), zcore(j)
            vnnz= vnnz + zcore(i) * zcore(j) * f0(1,1)
            if (ndump.ge.2 .and. (nbf(i).eq.0 .or. nbf(j).eq.0) .and.&
     &        vnn1.ne.0.D0) write (iw,2000) i,j,rij,enr
!	    **** add nuclear-electronic attraction to U for Fock matrix ***
            do kb= 0,nbf(i)-1
              ibk= ib+kb
              nm= matind(ibk) + ibk
              beta(nm)= beta(nm) - zcore(j) * (fracsj*f0(ind(kb),1) +&
     &                (1.D0-fracsj)*f0(ind(kb),2) )
            end do
            do kb= 0,nbf(j)-1
              jbk= jb+kb
              nm= matind(jbk) + jbk
              beta(nm)= beta(nm) - zcore(i) * (fracsi*f0(1,ind(kb)) +&
     &                (1.D0-fracsi)*f0(2,ind(kb)) )
            end do
          end if

        end do

!	**** INDO atom-atom diagonal block modifications ****
        if (nham.eq.2) then

          if (nbf(i).gt.1) then
            do j= 1,3
!	      **** s-p exchange integrals ****
              gamma(ib+0,ib+j)= fg(mg1sp,ia)
!	      **** p-p coulomb diagonal integrals ****
              gamma(ib+j,ib+j)= gamma(ib+j,ib+j) + 4*fg(mf2pp,ia)
            end do
!	    **** p-p coulomb off-diagonal integrals ****
            gamma(ib+2,ib+1)= gamma(ib+2,ib+1) - 2*fg(mf2pp,ia)
            gamma(ib+3,ib+2)= gamma(ib+2,ib+1)
            gamma(ib+3,ib+1)= gamma(ib+2,ib+1)
!	    **** p-p exchange integrals ****
            gamma(ib+1,ib+2)= 3*fg(mf2pp,ia)
            gamma(ib+2,ib+3)= gamma(ib+1,ib+2)
            gamma(ib+1,ib+3)= gamma(ib+1,ib+2)
          end if

          if (nbf(i).eq.9) then
            do j= 4,8
!	      **** s-d exchange integrals ****
              gamma(ib+0,ib+j)= fg(mg2sd,ia)
!	      **** d-d diagonal coulomb integrals ****
              gamma(ib+j,ib+j)= gamma(ib+j,ib+j) + 4*fg(mf2dd,ia) + &
     &                          36*fg(mf4dd,ia)
            end do
!	    **** p-d coulomb integrals ****
            gamma(ib+4,ib+1)= gamma(ib+4,ib+1) - 2*fg(mf2pd,ia)
            gamma(ib+4,ib+2)= gamma(ib+4,ib+1)
            gamma(ib+4,ib+3)= gamma(ib+4,ib+3) + 4*fg(mf2pd,ia)
            gamma(ib+5,ib+1)= gamma(ib+5,ib+1) + 2*fg(mf2pd,ia)
            gamma(ib+5,ib+2)= gamma(ib+5,ib+1)
            gamma(ib+6,ib+1)= gamma(ib+5,ib+1)
            gamma(ib+6,ib+2)= gamma(ib+5,ib+1)
            gamma(ib+7,ib+1)= gamma(ib+5,ib+1)
            gamma(ib+7,ib+3)= gamma(ib+5,ib+1)
            gamma(ib+8,ib+2)= gamma(ib+5,ib+1)
            gamma(ib+8,ib+3)= gamma(ib+5,ib+1)
            gamma(ib+5,ib+3)= gamma(ib+5,ib+3) - 4*fg(mf2pd,ia)
            gamma(ib+6,ib+3)= gamma(ib+5,ib+3)
            gamma(ib+7,ib+2)= gamma(ib+5,ib+3)
            gamma(ib+8,ib+1)= gamma(ib+5,ib+3)
!	    **** p-d exchange terms ****
            gamma(ib+1,ib+4)=   fg(mg1pd,ia) + 18*fg(mg3pd,ia)
            gamma(ib+2,ib+4)= gamma(ib+1,ib+4)
            gamma(ib+3,ib+4)= 4*fg(mg1pd,ia) + 27*fg(mg3pd,ia)
            gamma(ib+1,ib+5)= 3*fg(mg1pd,ia) + 24*fg(mg3pd,ia)
            gamma(ib+2,ib+5)= gamma(ib+1,ib+5)
            gamma(ib+1,ib+6)= gamma(ib+1,ib+5)
            gamma(ib+2,ib+6)= gamma(ib+1,ib+5)
            gamma(ib+1,ib+7)= gamma(ib+1,ib+5)
            gamma(ib+3,ib+7)= gamma(ib+1,ib+5)
            gamma(ib+2,ib+8)= gamma(ib+1,ib+5)
            gamma(ib+3,ib+8)= gamma(ib+1,ib+5)
            gamma(ib+3,ib+5)= 15*fg(mg3pd,ia)
            gamma(ib+3,ib+6)= gamma(ib+3,ib+5)
            gamma(ib+2,ib+7)= gamma(ib+3,ib+5)
            gamma(ib+1,ib+8)= gamma(ib+3,ib+5)
!	    **** d-d off-diagonal coulomb terms ****
            gamma(ib+5,ib+4)= gamma(ib+5,ib+4) - 4*fg(mf2dd,ia) +&
     &                                 6*fg(mf4dd,ia)
            gamma(ib+6,ib+4)= gamma(ib+5,ib+4)
            gamma(ib+6,ib+5)= gamma(ib+6,ib+5) + 4*fg(mf2dd,ia) -&
     &                                34*fg(mf4dd,ia)
            gamma(ib+7,ib+4)= gamma(ib+7,ib+4) + 2*fg(mf2dd,ia) -&
     &                                24*fg(mf4dd,ia)
            gamma(ib+8,ib+4)= gamma(ib+7,ib+4)
            gamma(ib+7,ib+5)= gamma(ib+7,ib+5) - 2*fg(mf2dd,ia) -&
     &                                 4*fg(mf4dd,ia)
            gamma(ib+7,ib+6)= gamma(ib+7,ib+5)
            gamma(ib+8,ib+5)= gamma(ib+7,ib+5)
            gamma(ib+8,ib+6)= gamma(ib+7,ib+5)
            gamma(ib+8,ib+7)= gamma(ib+7,ib+5)
!	    **** d-d exchange terms ****
            gamma(ib+4,ib+5)=  4*fg(mf2dd,ia) + 15*fg(mf4dd,ia)
            gamma(ib+4,ib+6)= gamma(ib+4,ib+5)
            gamma(ib+5,ib+6)= 35*fg(mf4dd,ia)
            gamma(ib+4,ib+7)=    fg(mf2dd,ia) + 30*fg(mf4dd,ia)
            gamma(ib+4,ib+8)= gamma(ib+4,ib+7)
            gamma(ib+5,ib+7)=  3*fg(mf2dd,ia) + 20*fg(mf4dd,ia)
            gamma(ib+6,ib+7)= gamma(ib+5,ib+7)
            gamma(ib+5,ib+8)= gamma(ib+5,ib+7)
            gamma(ib+6,ib+8)= gamma(ib+5,ib+7)
            gamma(ib+7,ib+8)= gamma(ib+5,ib+7)
          end if
        end if

      end do

!     **** calculation of irregular integrals for transition metals INDO/1S ****

      ntrmet= 0
      sqrt3= sqrt (3.D0)

      if (ifroti.eq.0) then
!	**** igrore Zerners terms, loosing rotational invariance of d orbs ****
        nirreg= 0

      else

      k= 0
      do i= 1,na
        ia= natm(i)
        natt(i)= 0
        if (nbf(i).gt.4) then
          ntrmet= ntrmet + 1
          itrmet(ntrmet)= i
!	  **** check to see if this type of atom already found ****
          j= 1
          do while (j.lt.i .and. natm(j).ne.natm(i))
            j= j + 1
          end do
          if (j.lt.i) then
!	    **** atom already done ****
            natt(i)= natt(j)
          else
            k= k + 1
            if (k.gt.mte) then
              write (iw,1000) mte
              write ( 6,1000) mte
              stop 'MTE'
            end if
            natt(i)= k
            nirreg= 0

!	    ********* G1(pd) and G3(pd) contributions *******
            xx= (  -fg(mg1pd,ia) -  3*fg(mg3pd,ia)) * sqrt3
            call irreg (k,xx,1,4,1,5)
            call irreg (k,xx,1,6,2,4)
            call irreg (k,xx,2,6,1,4)
            xx=  -3*fg(mg1pd,ia) + 21*fg(mg3pd,ia)
            call irreg (k,xx,1,6,2,5)
            xx= ( 2*fg(mg1pd,ia) -  9*fg(mg3pd,ia)) * sqrt3
            call irreg (k,xx,1,7,3,4)
            call irreg (k,xx,2,8,3,4)
            xx=                    15*fg(mg3pd,ia)
            call irreg (k,xx,1,7,3,5)
            call irreg (k,xx,1,8,2,7)
            call irreg (k,xx,1,8,3,6)
            call irreg (k,xx,2,7,3,6)
            xx= (   fg(mg1pd,ia) +  3*fg(mg3pd,ia)) * sqrt3
            call irreg (k,xx,2,5,2,4)
            xx=   3*fg(mg1pd,ia) - 21*fg(mg3pd,ia)
            call irreg (k,xx,2,6,1,5)
            xx=   3*fg(mg1pd,ia) -  6*fg(mg3pd,ia)
            call irreg (k,xx,2,8,1,7)
            call irreg (k,xx,3,7,1,5)
            call irreg (k,xx,3,7,2,6)
            call irreg (k,xx,3,8,1,6)
            xx=                  - 15*fg(mg3pd,ia)
            call irreg (k,xx,2,8,3,5)
            xx= ( - fg(mg1pd,ia) + 12*fg(mg3pd,ia)) * sqrt3
            call irreg (k,xx,3,7,1,4)
            call irreg (k,xx,3,8,2,4)
            xx= - 3*fg(mg1pd,ia) +  6*fg(mg3pd,ia)
            call irreg (k,xx,3,8,2,5)

!	    ********* F2(pd) contrinutions *************
            xx=  -2*fg(mf2pd,ia)  * sqrt3
            call irreg (k,xx,1,1,4,5)
            call irreg (k,xx,4,6,1,2)
            xx=   2*fg(mf2pd,ia)  * sqrt3
            call irreg (k,xx,2,2,4,5)
            xx=     fg(mf2pd,ia)  * sqrt3
            call irreg (k,xx,4,7,1,3)
            call irreg (k,xx,4,8,2,3)
            xx=   3*fg(mf2pd,ia)
            call irreg (k,xx,5,7,1,3)
            call irreg (k,xx,6,7,2,3)
            call irreg (k,xx,6,8,1,3)
            call irreg (k,xx,7,8,1,2)
            xx= - 3*fg(mf2pd,ia)
            call irreg (k,xx,5,8,2,3)

!	    ********* F2(dd) and F4(dd) contrinutions *************
            xx= (   fg(mf2dd,ia) -  5*fg(mf4dd,ia)) * sqrt3
            call irreg (k,xx,5,7,4,7)
            call irreg (k,xx,6,7,4,8)
            call irreg (k,xx,6,8,4,7)
            xx= (-  fg(mf2dd,ia) +  5*fg(mf4dd,ia)) * sqrt3
            call irreg (k,xx,5,8,4,8)
            xx=  -3*fg(mf2dd,ia) + 15*fg(mf4dd,ia)
            call irreg (k,xx,5,8,6,7)
            xx=   3*fg(mf2dd,ia) - 15*fg(mf4dd,ia)
            call irreg (k,xx,6,8,5,7)
            xx= (-2*fg(mf2dd,ia) + 10*fg(mf4dd,ia)) * sqrt3
            call irreg (k,xx,7,7,4,5)
            call irreg (k,xx,7,8,4,6)
            xx= ( 2*fg(mf2dd,ia) - 10*fg(mf4dd,ia)) * sqrt3
            call irreg (k,xx,8,8,4,5)
        
            if (ifroti.eq.1) then
!	      ********* R1(sppd) contrinutions *************
              xx= 0.1490712D0 * fg(mr1sppd,ia)
              call irreg (k, -xx,0,1,1,4)
              call irreg (k, -xx,0,2,2,4)
              call irreg (k,2*xx,0,3,3,4)
              xx= 0.2581989D0 * fg(mr1sppd,ia)
              call irreg (k, xx,0,1,1,5)
              call irreg (k, xx,0,1,2,6)
              call irreg (k, xx,0,1,3,7)
              call irreg (k, xx,0,2,1,6)
              call irreg (k,-xx,0,2,2,5)
              call irreg (k, xx,0,2,3,8)
              call irreg (k, xx,0,3,1,7)
              call irreg (k, xx,0,3,2,8)

!	      ********* R2(sdpp) contrinutions *************
              xx= 0.0894427D0 * fg(mr2sdpp,ia)
              call irreg (k, -xx,1,1,0,4)
              call irreg (k, -xx,2,2,0,4)
              call irreg (k,2*xx,3,3,0,4)
              xx= 0.1549193D0 * fg(mr2sdpp,ia)
              call irreg (k,  xx,1,1,0,5)
              call irreg (k, -xx,2,2,0,5)
              call irreg (k,  xx,1,2,0,6)
              call irreg (k,  xx,1,3,0,7)
              call irreg (k,  xx,2,3,0,8)

!	      ********* R2(sddd) contrinutions *************
              xx= 0.0638877 * fg(mr2sddd,ia)
              call irreg (k, 2*xx,4,4,0,4)
              call irreg (k,-2*xx,5,5,0,4)
              call irreg (k,-2*xx,6,6,0,4)
              call irreg (k,   xx,7,7,0,4)
              call irreg (k,   xx,8,8,0,4)
              call irreg (k,-2*xx,4,5,0,5)
              call irreg (k,-2*xx,4,6,0,6)
              call irreg (k,   xx,4,7,0,7)
              call irreg (k,   xx,4,8,0,8)
              xx= 0.1106567D0 * fg(mr2sddd,ia)
              call irreg (k,   xx,7,7,0,5)
              call irreg (k,  -xx,8,8,0,5)
              call irreg (k,   xx,7,8,0,6)
              call irreg (k,   xx,5,7,0,7)
              call irreg (k,   xx,6,8,0,7)
              call irreg (k,  -xx,5,8,0,8)
              call irreg (k,   xx,6,7,0,8)
            end if

          end if
        end if
      end do
      end if
      nte = k

      return
      end

!     *************************************************************************

      subroutine irreg (k,xx,i1,i2,i3,i4)

      use reimers_C, only: ig1, ig2, ig3, ig4, g, nirreg
      USE chanel_C, only : iw
      implicit none
      double precision   xx
      integer            i1, i2, i3, i4, k

!     **** stores one of Zerners rotational-invariance maintaining ****
!     **** irregular two-electron integrals into integral list      ****

      if (nirreg.ge.75) stop 'IRREGULAR INTEGRALS'
      nirreg= nirreg + 1
      g(nirreg,k)= xx
      ig1(nirreg,k)= i1
      ig2(nirreg,k)= i2
      ig3(nirreg,k)= i3
      ig4(nirreg,k)= i4

      return
      end

!     *************************************************************************

      subroutine uform (i,fg,f0,zcore)

      use reimers_C, only: matind, natm, nbf, ibf, nbt, nbeta,&
          mis2, mip2, mid2, mic2, mic1, mg1sp, mf2pp,&
          mg2sd, mg1pd, mf2pd, mg3pd, mf2dd, mf4dd, beta
      USE chanel_C, only : iw

      implicit none
      double precision  fg(24), f0(2,2), &
                        sp, uu, zcore
      integer           i, k, &
                        i0, i1, ib, ic, id, icindx, icore, ip, is, nm

!     *** This constructs the U matrix from either (I+A)/2 or I		****
!     *** using mixed electronic configurations for the uncharged atoms	****

!     *** mis2= s ioniz pot from d(n-2)s2 or from sp, mic2= its probability ****
!     *** mis1= s ioniz pot from d(n-1)s1,            mic1= its probability ****
!     *** sim for p and d ioniz pots, mid0 is d ioniz pot from d(n)s0       ****
!     *** sp is the sum of the d(n-2)s2 and d(n-1)s1 probabilities	    ****

      icore= nint (zcore)
      i0= ibf(i)
      k= natm(i)
      sp= fg(mic2) + fg(mic1)
      nbeta = 5
!     **** loop over basis functions on this atom ****

      do ib= 0,nbf(i)-1
        i1= i0 + ib
        nm= matind(i1) + i1
        beta(nm)= 0.D0

!       **** loop over the 3 configurations d(n-ic)s(ic) ****

        do ic= 2,0,-1
          icindx= 2 - ic
          uu= 0.D0
          if (nbt(i1).eq.0) then
!           ***** s function, note there is no entry for d(n)s0 ****
            if (ic.gt.0) then
              uu= fg(mis2+icindx)
              if (nbeta.eq.5) then

!               **** INDO/1S I only, all atoms ****
                if (nbf(i).eq.9) then
                  is= ic
                  ip= 0
                  id= icore - is
                else
!                 **** sp configuration only ****
                  is= min (2,icore)
                  id= 0
                  ip= icore - is
                end if
                uu= uu  - (is-1) * f0(1,1)&
     &                  - ip     * (f0(1,1)-fg(mg1sp)/2.D0)&
     &                  - id     * (f0(1,2)-fg(mg2sd)/2.D0)
              else if (nbeta.lt.4) then

!               **** (I+A)/2, 2nd row code only ****
                uu= uu - (zcore-0.5D0) * f0(1,1)
 
                if (k.eq.4) then
                  uu= uu + fg(mg1sp)/4
                else if (k.gt.4) then
                  uu= uu + fg(mg1sp)/2*(zcore-1.5D0)
                end if
 
              end if
            else
!             **** renormalization as no d(n)s(0) ****
              beta(nm)= beta(nm) / sp
            end if

          else if (nbt(i1).le.3) then
!           ***** p function, note there is no entry for d(n)s0 ****
            if (ic.gt.0) then
              uu= fg(mip2+icindx)
               if (nbeta.eq.5) then
!               **** INDO/1S I only, all atoms ****
                if (nbf(i).eq.9) then
                  ip= 1
                  if (k.ne.20) then
                    is= ic - ip
                    id= icore - ic
                  else
!                   **** ZINDO has d->p for Ca and s->p for all else ****
                    is= ic
                    id= icore - ic - 1
                  end if
                else
!                 **** sp configuration only ****
                  id= 0
                  ip= max (1,icore-2)
                  is= icore - ip
                end if
                uu= uu - (ip-1.D0) * (f0(1,1) - fg(mf2pp)*2.0D0)&
     &                 - is        * (f0(1,1) - fg(mg1sp)*0.5D0)&
     &                 - id        * (f0(1,2) - fg(mg1pd)&
     &                                        - fg(mg3pd)*10.5D0)

              else if (nbeta.lt.4) then
!               **** (I+A)/2, 2nd row code only ****
                uu= uu - (zcore-0.5D0) * f0(1,1)
                 if (k.eq.3) then
                  uu= uu + fg(mg1sp)/4
                else if (k.eq.4) then
                  uu= uu + fg(mg1sp)*0.75D0
                else if (k.gt.4) then
                  uu= uu + fg(mg1sp)&
     &                   + fg(mf2pp)*2*(zcore-2.5D0)
                end if
 
              end if

            else
!             **** renormalization as no d(n)s(0) ****
              beta(nm)= beta(nm) / sp
            end if

          else
!           **** d function, INDO/1S code only, all atoms ****
            ip= 0
            is= ic
            id= icore - is
            uu= fg(mid2+icindx)
            uu= uu - (id-1.D0)  * (f0(2,2) - fg(mf2dd)*14.D0/9.D0&
     &                                     - fg(mf4dd)*14.D0)&
     &             - is         * (f0(1,2) - fg(mg2sd)*0.5D0)&
     &             - ip         * (f0(1,2) - fg(mg1pd)&
     &                                     - fg(mg3pd)*10.5D0)
          end if
!         **** multiply contribution to core integral by its weight ****
          beta(nm)= beta(nm) + uu * fg(mic2+icindx)
        end do
      end do

      return
      end

!     *************************************************************************

      function gamma1 (gammk,gammk1,k,k1,r)

      use reimers_C, only: nham, au2ev, au2ang, au2cm, zeta, tomk
      implicit none
      double precision   d, p, r, s, t, &
                         ab, alph, ar28, ar37, beta, &
                         gamma1, gamma2, gammk, gammk1, &
                         pa, pb, rab, ri, ri2, &
                         tk, za, zb, zab
      integer            i, j, k, k1, nm

!     **** two-electron repulsion integrals ****

      nham=2
!     **** mataga integrals (with Tomono/Weiss modifications) ****
      gamma2= au2ev / (r/au2ang + tomk*2.D0*au2ev/(gammk+gammk1) )

      if (nham.ne.3) gamma2= gamma2 * tomk
      gamma1= gamma2
      return

!    **** pariser integrals ****
    2 continue
      d= 1.14925/zeta(k) - 1.14925/zeta(k1)
      s= 1.14925/zeta(k) + 1.14925/zeta(k1)
      rab= r
      if (rab .ge. 2.8) then
        ri= 1.0/rab
        ri2= ri*ri
        gamma1= 7.1975*ri*(1./sqrt(1.+ri2*d*d) + 1./sqrt(1.+ri2*s*s))
        return
      else
        ab= 0.5 * (gammk + gammk1)
        ar28= ab - 2.570533/sqrt(1.+.127551*d*d)-2.570533/&
     &        sqrt(1.+.127551*s*s)
        ar37=ab-1.945276/sqrt(1.+.073046*d*d)-1.945276/&
     &        sqrt(1.+.073046*s*s)
        alph= 0.107250 *(13.69*ar28 - 7.84 * ar37)
        beta= 0.107250 * (2.8 * ar37 - 3.7 * ar28)
        gamma1= ab - rab*(alph + rab * beta)
        return
      end if

!     ***** theoretical integrals (Roothaan) *****
    3 continue
      za= zeta(k)
      if (r.eq.0.D0) then
!	**** diagonal term ****
        if (k.gt.2) then
          gamma1= 9.882703 * za
        else
          gamma1= 17.0025 * za
        end if
        return
      end if

      rab= r / au2ang
      zb= zeta(k1)
      zab= 0.5D0 * (za + zb)
      p= zab * rab
      if (k.ne.k1) then
        t= (za - zb)/( za + zb)
        tk= 0.5D0 * (t + 1./t)
        pa= za * rab
        pb= zb * rab
      end if

      if (k .le. 2) then
        if (k1.le. 2) then
!	  **** H-H interaction ****
          gamma1= au2ev*za/p*(1.-(((.16666666*p+.75)*p+1.375)*p+1.)*&
     &        exp(-2.*p))

        else
!	  **** H-non H interaction ****
          gamma1= au2ev * zab/p * (1.&
     &     + (1.-tk)**3 * exp(-2.*pa) *&
     &          ( (0.25*tk+.3125+.125*pa)*tk - 0.0625)&
     &     - (1.+tk)**2 * exp(-2.*pb) *&
     &        ( ( (.083333333*pb+.25*(2.-tk))*pb + .375*((tk-3.)*tk+3.))&
     &        *pb - ((.25*tk-.9375)*tk+1.375)*tk-.9375) )
        end if

      else
        if (k1.le.2) then
!	  **** H - non H interaction ****
          tk= - tk
          gamma1= au2ev * zab/p * (1.&
     &     + (1.-tk)**3 * exp(-2.*pb) *& 
     &          ( (0.25*tk+.3125+.125*pb)*tk - 0.0625)&
     &     - (1.+tk)**2 * exp(-2.*pa) *&
     &        ( ( (.083333333*pa+.25*(2.-tk))*pa + .375*((tk-3.)*tk+3.))&
     &        *pa - ((.25*tk-.9375)*tk+1.375)*tk-.9375) )

        else if (k.eq.k1) then
!	  **** identical non H - non H interaction ****
          gamma1= au2ev*zab/p * (1.-(((((((.00079365*p+.0083333333)*&
     &        p +.05)*p+.20833333)*p+.61979167 )*p+1.2734375)*p&
     &        +1.63671875)*p+1.0)*exp(-2.*p))


        else
!	  **** non-identical non H - non H interaction ****
          gamma1 = au2ev * zab/p * (1.+&
     &      (1.-tk)**3 * exp(-2.*pa) * &
     &         ( ( (.041666666*tk*pa+(.25*tk+.3125)*tk-.0625)*pa +&
     &             ((.625*tk+1.375)*tk+.59375)*tk-.34375 ) * pa +&
     &           ( ((.625*tk+1.875)*tk+1.6875)*tk+.0625 ) *tk -.5 )&
     &     -(1.+tk)**3 * exp(-2.*pb) *&
     &         ( ( (.041666666*tk*pb-(.25*tk-.3125)*tk+.0625)*pb +&
     &             ((.625*tk-1.375)*tk+.59375)*tk+.34375 ) * pb -&
     &           ( ((.625*tk-1.875)*tk+1.6875)*tk-.0625 ) *tk +.5 ) )
        end if
      end if
      return

!      ****** OHNO integrals *****
    4 continue
      gamma1= au2ev / sqrt ( (r/au2ang)**2 +&
     &        (2.D0*au2ev/(gammk+gammk1))**2 )
      return

      end

!     *************************************************************************

      function enucrep (r,gamma)
         implicit none
      double precision :: r, gamma, enucrep

!     **** nuclear repulsion energy ****

      enucrep= gamma

      return
      end

!     *************************************************************************

      subroutine beta1 (s,betao,beta)

      use reimers_C,only: n
      USE molkst_C, ONLY: mpack, norbs
      USE chanel_C, only : iw
      implicit none

      double precision  betao(norbs), beta(mpack), s(mpack)
      integer           i, j, nm

 1001 format (2i4,f10.6)
 1003 format (/' beta or (I+A)/2 read' / '   i   j    beta'/)
 1004 format (2i4,f10.5)

!     **** construction of the one-electron matrix ****

      nm= 0
      do i= 1,n
        do j= 1,i
          nm= nm + 1
          if (i.ne.j) beta(nm)= (betao(i) + betao(j)) * s(nm) * 0.5D0
        end do
      end do
      return

      end

!     *************************************************************************

      subroutine dipol (x,y,z,dm)

!     **** calculates dipole moment matrix elements in AO basis   ****
!     **** units of dipole are eA				  ****
      use reimers_C, only: n, na, nb2, matind, au2ev, au2ang, au2cm, &
          natm, nbf, ibf, iat, nbt, nprn, zeta, zetad, zetawt, &
          dipsym, ndtype, fact, zcore
!      USE chanel_C, only : iw

      implicit none
      double precision  x(na),y(na),z(na),dm(nb2,3), &
                        two, three, bohrmg, &
                        sqrt3, sqrt10, srt2d5, srt6d5, srt3d2, &
                        al, aaa, am, amu, beta, bmu, delx, dely, delz, &
                        pd, r, sp, t, zj, zk
      integer           i, j, &
                        i1, ia, ib, iss, &
                        j2, j3, ja, jb, jnat, &
                        k3, ka, kb, knat, mmmm, mu, nm, nu
      

      data two,three    /2.D0,3.D0/
!     **** the Borh magneton ehbar/(2.mass_e.c) in units of e Ang ****
      data bohrmg       /0.0019308D0/

      sqrt3=  sqrt (three)
      sqrt10= sqrt (10.D0)
      srt2d5= sqrt (0.4D0)
      srt6d5= sqrt (1.2D0)
      srt3d2= sqrt (1.5D0)

!     **** start by zeroing all elements ****
      ndtype = 1
      do i= 1,nb2
        do j= 1,3
          dm(i,j)= 0.D0
        end do
      end do
      if (ndtype.le.1) then

!	**** usual length formalism of dipole integrals, REAL ****

        dipsym= 1.D0
        do i= 1,na
          ia= natm(i)
          ib= ibf(i)
          do i1= ib,ib+nbf(i)-1
!	    **** contribution from atomic charges ****
            nm= matind(i1) + i1
            dm(nm,1)= - x(i)
            dm(nm,2)= - y(i)
            dm(nm,3)= - z(i)
           end do

          if (nbf(i).gt.1 .and. ndtype.eq.1) then
!	    **** include one-centre hybridization terms ****

!	    **** sp polarization term ****
            j2= 2*nprn(ib) + 1
            amu= zeta(ia)
            sp= - au2ang * j2 / amu / two / sqrt3
            if (nbf(i).gt.4) then
!	      **** pd polarization term ****
              j3= j2 - 2
              k3= j2
              pd= 0.D0
              do iss= 1,2
                bmu= zetad(iss,ia)
                t= dsqrt ( (two*amu)**j2 * (two*bmu)**j3 /&
     &             (fact(j2)*fact(j3)*5.D0) ) * fact(k3) / (amu+bmu)**k3
                pd= pd + t * zetawt(iss,ia)
              end do
            end if

            do mu= ib+1,ib+nbf(i)-1
              do nu= ib,mu-1
                nm= matind(mu) + nu
                if (nbt(nu).eq.0) then
!		  **** store any sp terms ****
                  if (nbt(mu).eq.1) then
                    dm(nm,1)= sp
                  else if (nbt(mu).eq.2) then
                    dm(nm,2)= sp
                  else if (nbt(mu).eq.3) then
                    dm(nm,3)= sp
                  end if
                else if (nbt(nu).le.3 .and. nbt(mu).gt.3) then
!		  **** pd polarization terms ****
                  mmmm= nbt(nu) + 3*(nbt(mu)-4)
      go to (1302,1304,1306,1308,1310,1800,1312,1308,1800,1314,&
     & 1800,1308,1800,1314,1312),mmmm
 1302 dm(nm,1)=au2ang*pd/sqrt3
      go to 1800
 1304 dm(nm,2)=au2ang*pd/sqrt3
      go to 1800
 1306 dm(nm,3)=-au2ang*pd*two/sqrt3
      go to 1800
 1308 dm(nm,1)=-au2ang*pd
      go to 1800
 1310 dm(nm,2)=au2ang*pd
      go to 1800
 1312 dm(nm,2)=-au2ang*pd
      go to 1800
 1314 dm(nm,3)=-au2ang*pd
      go to 1800

 1800       continue
                end if
              end do
            end do

          end if
        end do
      else if (ndtype.eq.2) then
!	**** dipole velocity eqns of Hush and Williams CPL 8 (1971) 179 ****
!	**** but note their last two eqns for Iz must be mult by two ****
!	**** only p-p integrals hence only programmed for PPP calc ****
!	**** mx,my,mz are imaginary and have units of Angstrom**-1 ****
!	**** works only for 2nd row elements			   ****
!	**** Integrals are actually IMAGINARY			   ****

        dipsym= -1.D0
        nm= 0
        do jb= 1,n
          ja= iat(jb)
          jnat= natm(ja)
!	  **** zeta is in bohr, convert to Ang**-1 ****
          zj= zeta(jnat) / au2ang

          do kb= 1,jb
            nm= nm + 1
            if (jb.ne.kb) then
              ka= iat(kb)
              knat= natm(ka)
              zk= zeta(knat) / au2ang

              delx= x(ka) - x(ja)
              dely= y(ka) - y(ja)
              delz= z(ka) - z(ja)
              r= sqrt (delx**2 + dely**2 + delz**2)
              al= 0.5D0 * r * (zk + zj)
              if (zk.ne.zj) then
!	        **** different atom types ****
                am= 0.5D0 * r * (zk - zj)
                aaa= - 32.D0 / (r/2.D0)**4 * &
     &            (zk*zj)**3 * sqrt (zk*zj) / (zk**2-zj**2)**5 *&
     &            (3.D0+3.D0*al+al**2) *&
     &            (exp(-zk*r)*(3.D0+3.D0*am+am**2) -&
     &             exp(-zj*r)*(3.D0-3.D0*am+am**2) )
              else
!	        **** same type of atom (Chong formula) ****
                aaa= - zk * al/15.D0 * exp(-al) * (3.D0+3.D0*al+al**2)
              end if
              aaa= aaa / r
              dm(nm,1)= aaa * delx
              dm(nm,2)= aaa * dely
              dm(nm,3)= aaa * delz
            end if
          end do
        end do

      else if (ndtype.eq.3) then

!	**** Magnetic transition moments, IMAGINARY ****

        dipsym= -1.D0
        do i= 1,na
          if (nbf(i).ne.1) then
!	    **** include p-p one-centre hybridization terms ****
            ib= ibf(i)
            dm(matind(ib+3)+ib+2,1)=   bohrmg
            dm(matind(ib+3)+ib+1,2)= - bohrmg
            dm(matind(ib+2)+ib+1,3)=   bohrmg
          end if
        end do

      end if

      return
      end

!     *************************************************************************

      subroutine efmods (beta,zcore,dm)
      USE chanel_C, only : iw
      use reimers_c, only: n, na, nb2, matind, ibf, vnn, ef

      implicit none
      double precision   dm(nb2,3),zcore(na),beta(nb2), ecoref
      integer            i, mu, nm, nu


!     **** Electric field effect on core energies and 1 el operators ****
!     **** using the definition: Energy= -F * sum RiQi ****
!     **** U(i)= -beta(i,i)= - electron_charge * (-F.r)= -F.r ****
!     **** beta(i,j)= electron_charge * (-F<s|x|px>) ****

!     **** mods to one-electron integrals ****
      nm= 0
      do mu= 1,n
        do nu= 1,mu
          nm= nm + 1
          beta(nm)= beta(nm) - (ef(1)*dm(nm,1) + ef(2)*dm(nm,2) +&
     &                          ef(3)*dm(nm,3) )
        end do
      end do

!     **** nuclear interaction ****
      ecoref= 0.D0
      do i= 1,na
        mu= ibf(i)
        nm= matind(mu) + mu
        ecoref= ecoref + zcore(i) * (ef(1)*dm(nm,1) + ef(2)*dm(nm,2) +&
     &                               ef(3)*dm(nm,3))
      end do
      vnn= vnn + ecoref
      return
 1018 format (/' Electric Field in new coordinates:',3f10.6,' V/Ang'/&
     &         ' Energy of cores in Field=',f15.6,' eV')
      end

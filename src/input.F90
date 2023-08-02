      subroutine ptgrp (mt)

!     **** optionally determins molecular point grp			***
!     **** and then determins the symmetry AO transformation		***
!     **** if mt=1, then transformation is of AOs but			***
!     **** if mt=0, then uses only s functions and forms trans of gamma ***

      use reimers_C, only: n, na, natm, nbf, ibf, iat, nbt, moass, kass, &
          nshell, norb, iatsym, kind, ncisym, icisym, cisym, nmptg, &
          nmrep, stwt, nptg, nr, ns, nslwr, nsupr, nstr, istr, weight, &
          wk0, cc0, x, y, z, xz, zcore, ion, ndump
      USE molkst_C, only: norbs
      USE chanel_C, only: iw

      implicit none
      logical		op(7),isinc(na),isrep
      character*20	lam(2)
      character*15	mtype(0:1),lll
      character*3       rep
      integer           ictbl(8,7),nrep(8),jsymt(0:8,7),ioploc(7,8), &
                        irs(8), isign(8), jrels(na,7)
      integer           i, j, k, l, & 
                        iao, ier, iop, ir, is, ix, ja, jsign, jr, jrs, &
                        k1, kc, krs, mc, mt, mu, nop, nrs, nu
      double precision  sum, summ, sumx, sumy, sumz, xo, yo, zo, &
                        pmom(3), ef(3)

      data lam /'     principal      ',' centroid of charge '/
      data mtype	/' GAMMA TRANSF  ',' SYMM ORBITALS '/

!     **** effect on s, p, d_eg, and d_t2g orbitals of the 7 symm operators ****
!     **** 1= symmetric, 2= asymmetric ****
!     **** line for 7 sym ops: C2z, C2y, C2x, i, SIGxy, SIGxz, SIGyz ****
      data jsymt   / 1, -1,-1,1,  1,1,  1,-1,-1,&
     &		     1, -1,1,-1,  1,1, -1,1,-1,&
     &		     1,  1,-1,-1, 1,1, -1,-1,1,&
     &		     1, -1,-1,-1, 1,1,  1,1,1,&
     &		     1,  1,1,-1,  1,1,  1,-1,-1,&
     &		     1,  1,-1,1,  1,1, -1,1,-1,&
     &		     1, -1,1,1,   1,1, -1,-1,1  /

!     **** ictbl is the character table: the E op is not included ****
!     **** it is actually D2h, others are eq spaced rows with first cols ****
      data ictbl / 1,1,1,1,-1,-1,-1,-1, 1,1,-1,-1,1,1,-1,-1,&
     &		   1,1,-1,-1,-1,-1,1,1, 1,-1,1,-1,1,-1,1,-1,&
     &		   1,-1,1,-1,-1,1,-1,1, 1,-1,-1,1,1,-1,-1,1,&
     &		   1,-1,-1,1,-1,1,1,-1 /

!     **** ioploc is allowed operators for the point groups ****
      data ioploc / 7*0, 5,6*0, 1,6*0, 1,6,7,4*0, 1,2,3,4,5,6,7,&
     &		    4,6*0, 1,4,5,4*0, 1,2,3,4*0 /

!     **** number of irreducible reps in each point group ****
      data nrep /1,2,2,4,8,2,4,4/

 1017 format (i4,3f10.5,2x,a2)
 1021 format (/' The point group of the molecule is ',a3,a13)
 1022 format (/' This molecule belongs to one of the degenerate point',&
     & 'groups.' / ' The point group given is the highest',&
     & ' non-degenerate subgroup in the orientation used')
 1024 format (/' coordinates in',a20,'axis system'/&
     & ' atom     x',9x,'y',9x,'z     at.#'/)
 1030 format (/' ERROR: PPP on non-planar molecule; use PI_DIR option')
 1040 format (/' ERROR: SYMMETRY OPERATOR #',i2,' NOT PRESENT')
 1100 format (a15,8(1x,a3,'=',i3))
 1104 format (/' MO Electron assignment into shells by symmetry'/)
 1105 format (7x,'shell=',i2)
 1110 format (' ERROR: for irred rep',i2,', only',i4,' orbs exist but '&
     &	      ,i4,' have been assigned by symmetry!')
 1111 format (' ERROR: for shell',i2,',',i4,' orbs exist but'&
     &	      ,i4,' have been assigned by symmetry!')
 1120 format (/' TRANSFORMATION:'/)

      if (nptg.gt.0) goto 7

!  correct coordinates to center of gravity

      summ= 0.D0
      sumx= 0.D0
      sumy= 0.D0
      sumz= 0.D0
      do i= 1,na
        j= natm(i)
        sumx= sumx + x(i) * weight(j)
        sumy= sumy + y(i) * weight(j)
        sumz= sumz + z(i) * weight(j)
        summ= summ +        weight(j)
      end do
      xo= sumx / summ
      yo= sumy / summ
      zo= sumz / summ
      do i= 1,na
        x(i)= x(i) - xo
        y(i)= y(i) - yo
        z(i)= z(i) - zo
      end do

!  determine principal moments of inertia, check for degenerate group
      do k= 1,3
        do l= k,3
          sum= 0.D0
          do i= 1,na
            j= natm(i)
            sum= sum + weight(j)*xz(i,k)*xz(i,l)
          end do
          cc0(k,l)= sum
          cc0(l,k)= sum
        end do
      end do

      call tred2e (3,3,cc0,pmom,wk0,cc0)
      call tql2e  (3,3,    pmom,wk0,cc0,ier)

      if (abs(pmom(1)-pmom(2)).lt..1 .or. abs(pmom(2)-pmom(3)).lt..1)&
     &                then
!       **** degenerate, dont rotate coords ****
        write (iw,1022)
        goto 10
      end if

!     determine prin axes in an equal mass system

      do k= 1,3
        do l= k,3
          sum= 0.D0
          do i= 1,na
            j= natm(i)
            sum= sum + xz(i,k)*xz(i,l)
          end do
          cc0(k,l)= sum
          cc0(l,k)= sum
        end do
      end do

      call tred2e (3,3,cc0,pmom,wk0,cc0)
      call tql2e  (3,3,    pmom,wk0,cc0,ier)

!     rotation of coordinates 

      do i= 1,na
        do k= 1,3
          sum= 0.D0
          do l= 1,3
            sum= sum + cc0(l,k)*xz(i,l)
          end do
          pmom(k)= sum
        end do
        do k= 1,3
          xz(i,k)= pmom(k)
        end do
      end do
      do k= 1,3
        pmom(k)= cc0(1,k)*ef(1) + cc0(2,k)*ef(2) + cc0(3,k)*ef(3)
      end do
      do k= 1,3
        ef(k)= pmom(k)
      end do

!     if charged ion, translate to center of charge

   10 if (ion.eq.1) go to 7
      summ= 0.D0
      sumx= 0.D0
      sumy= 0.D0
      sumz= 0.D0
      do 8 i= 1,na
      summ= summ + zcore(i)
      sumx= sumx + zcore(i)*x(i)
      sumy= sumy + zcore(i)*y(i)
    8 sumz= sumz + zcore(i)*z(i)
      xo= sumx / summ
      yo= sumy / summ
      zo= sumz / summ
      do 9 i= 1,na
      x(i)= x(i) - xo
      y(i)= y(i) - yo
    9 z(i)= z(i) - zo

7     continue

!     ********************* search for symm operators *******************
!     **** 1= yz plane, 2= xz plane, 3=xy plane, 4=x axis, 5=y axis, ****
!     **** 6= z axis, 7= centre ****

      nop= 0
      do i= 1,7
        call serch (op(i),jrels(1,i),jsymt(0,i),xz,natm,na)
        if (op(i)) nop= nop + 1
      end do

      if (nptg.gt.0) then
!       **** verify that all required operators are present ****
        do i= 1,nrep(nptg)-1
          iop= ioploc(i,nptg)
          if (.not.op(iop)) then
            write (iw,1040) iop
            write ( 6,1040) iop
            stop 'RESTART SYMMETRY'
          end if
        end do

      else

!     ************************ determine point group ********************

      if (nop.eq.0) then
!       **** no symmetry ****
        nptg= 1

      else if (nop.eq.7) then
!       **** D2H ****
        nptg= 5

      else if (nop.eq.1) then
        if (op(7)) then
!         *** Cs symmetry, is yz plane, ensure that its the xy plane ****
          nptg= 2
          call swaap (xz,ef,op,jrels,1,3,na)
        else if (op(6)) then
!         *** Cs symmetry, is xz plane, ensure that its the xy plane ****
          nptg= 2
          call swaap (xz,ef,op,jrels,2,3,na)
        else if (op(5)) then
!         *** Cs symmetry ****
          nptg= 2
        else if (op(4)) then
!         *** Ci symmetry ****
          nptg= 6
        else if (op(3)) then
!         **** C2 symmetry, is x axis, ensure that its about z axix ****
          nptg= 3
          call swaap (xz,ef,op,jrels,1,3,na)
        else if (op(2)) then
!         **** C2 symmetry, is y axis, ensure that its about z axix ****
          nptg= 3
          call swaap (xz,ef,op,jrels,2,3,na)
        else if (op(1)) then
!         **** C2 symmetry ****
          nptg= 3
        end if

      else
        if (op(4)) then
!         **** C2h symmetry, ensure C2 axis is z axis ****
          nptg= 7
          if (op(2)) then
            call swaap (xz,ef,op,jrels,2,3,na)
          else if (op(3)) then
            call swaap (xz,ef,op,jrels,1,3,na)
          end if

        else if (op(6).or.op(7)) then
!         **** C2V symmetry, ensure that C2 axis is Z axis ****
          nptg= 4
          if (op(2)) then
            call swaap (xz,ef,op,jrels,2,3,na)
          else if (op(3)) then
            call swaap (xz,ef,op,jrels,1,3,na)
          end if

        else
!         **** D2 symmetry ****
          nptg= 8
        end if
      end if

!     ********** output of modified coordinates **************

      write (iw,1021) nmptg(nptg)
      write (iw,1024) lam(ion)
      do i= 1,na
          write (iw,1017) i,x(i),y(i),z(i),iatsym(natm(i))
      end do

      end if

!     ********************************************************************

!     ******** construct symmetry transform of atomic orbitals ******

!     **** nber of irreducible reps and step siie thru character table rows ****
      nr= nrep(nptg)
      kc= 8 / nr

!     **** loop over all irreducible reps looking for symm orbs ****
      if(.not. allocated(nstr)) allocate(nstr(norbs))
      if(.not. allocated(istr)) allocate(istr(8,norbs))

      nu= 0
      do ir= 1,nr
        ns(ir)= 0
        nslwr(ir)= nu+1

!       **** flags to indicate if an atom has been included or not ****
        do i= 1,na
          isinc(i)= .false.
        end do

!       **** scan over unique atoms ****
        mu= 0
        do i= 1,na
          if (isinc(i)) then
!           **** atom found to be related to an earlier one ****
            mu= mu + nbf(i)

          else
!           **** proceed for first atom in symmetry-related set ****

            do j= 1,nbf(i)
!              write (6,*) 'processing atom',i,' basis fn',j
              mu= mu + 1
              nrs= 1
              irs(nrs)= i
              isrep= .true.
              do jr= 2,nr
                isign(jr)= 0
              end do
              isign(1)= 1

!             **** loop over allowed symmetry ops for this point grp ****
              do k= 1,nr-1
                k1= ioploc(k,nptg)

                if (jrels(i,k1).ne.0) then
!                 **** atom is related to another by this symm op ****
                  krs= 0
                  do jrs= 1,nrs
                    if (irs(jrs).eq.jrels(i,k1)) krs= jrs
                  end do
                  if (krs.eq.0) then
!                   **** add this atom to list of symmetry-related atoms ****
                    nrs= nrs + 1
                    irs(nrs)= jrels(i,k1)
                    isinc(jrels(i,k1))= .true.
                    krs= nrs
                  end if

!                 **** calc sign of atom for this irr rep ****
                  jsign= jsymt(mt*nbt(mu),k1) * ictbl(ir*kc,k)
                  if (isign(krs).eq.0) then
!                   **** sign for this basis function not yet assigned ****
                    isign(krs)= jsign
                  else if (isign(krs).ne.jsign) then
!                   **** confused sign info means this irr rep not present ***
                    isrep= .false.
                  end if

                else
!                 **** atom maps onto itself by this symm op ****
                  if (jsymt(mt*nbt(mu),k1).ne.ictbl(ir*kc,k)) then
!                   **** this rep does not transform as this ao does ****
                    isrep= .false.
                  end if
                end if

              end do

!             **** if irr rep is found, map out ao signs ****
              if (isrep) then
                nu= nu + 1
                ns(ir)= ns(ir) + 1
                nstr(nu)= nrs
                do jrs= 1,nrs
                  if (isign(jrs).eq.0) isign(jrs)= 1
                  iao= mu - ibf(irs(1)) + ibf(irs(jrs))
                  istr(jrs,nu)= isign(jrs) * iao
                end do
              end if

            end do
          end if
        end do
        nsupr(ir)= nu
      end do

!     **** weights for transform vectors ****
      if(.not. allocated(stwt)) allocate(stwt(n))

      do i= 1,n
        stwt(i)= sqrt (1.D0/nstr(i))
      end do

!!     **** expansion of transformation matrix ****
!
!      do i= 1,n
!	do j= 1,n
!	  st(i,j)= 0.D0
!	end do
!      end do
!      do i= 1,n
!	do j1= 1,nstr(i)
!	  j= abs(istr(j1,i))
!	  st(j,i)= sign(1,istr(j1,i)) * stwt(i)
!	end do
!      end do
!      write (iw,*) 'transformation matrix'
!      do i= 1,n
!	write (iw,'(i4,21f6.2)') i,(st(i,j),j=1,min(21,n))
!      end do

!     **** output of symmetry transform ****

      if (nptg.gt.1) then
        write (iw,*)
        write (iw,1100) mtype(mt),(nmrep(ir,nptg),ns(ir),ir=1,nr)
        write ( 6,1100) mtype(mt),(nmrep(ir,nptg),ns(ir),ir=1,nr)
      end if

!     **** electron assignment based on symmetry info ****

      if (mt.eq.1 .and. moass.eq.1) then
        write (iw,1104)

        do is= 1,nshell
!          read (ix,*) (kass(ir,is),ir=1,nr)
          write (lll,1105) is
          write (iw,1100) lll,(nmrep(ir,nptg),kass(ir,is),ir=1,nr)
          i= 0
          do ir= 1,nr
            i= i + kass(ir,is)
          end do
          if (i.ne.norb(is)) then
            write (iw,1111) is,norb(is),i
            stop 'in PTGRP'
          end if
        end do

        do ir= 1,nr
          i= 0
          do is= 1,nshell
            i= i + kass(ir,is)
          end do
          is= nshell + 1
          kass(ir,is)= ns(ir) - i
          if (kass(ir,is).lt.0) then
            write (iw,1110) ir,ns(ir),i
            stop 'in PTGRP'
          end if
        end do

        lll= ' UNOCC ORBITALS'
        write (iw,1100) lll,(nmrep(ir,nptg),kass(ir,is),ir=1,nr)
      end if

!     ********* max nber of states of each symmetry to use in CI calc ********

      if (mt.eq.1) then
      do i= 1,nr
        irs(i)= 0
        if (ncisym.eq.0) irs(i)= mc
      end do
      do ir= 1,nr
        rep= nmrep(ir,nptg)
        call upcas(rep)
        do i= 1,ncisym
          if (cisym(i).eq.rep) irs(ir)= icisym(i)
        end do
      end do
      do i= 1,nr
        icisym(i)= irs(i)
      end do
      end if

      if (nptg.eq.1) return

!     ************** detailed dump of symmetry transform coefficients *********

      write (iw,1120)
      mu= 0
      do ir= 1,nr
        write (iw,*) 'irred rep ',nmrep(ir,nptg),'count=',ns(ir)
        do i= 1,ns(ir)
          mu= mu + 1
          j= abs(istr(1,mu))
          ja= iat(j)
          write (iw,'(i4,1x,a2,1x,a3,1h:,8i5)') mu, iatsym(natm(ja)),&
     &        kind(nbt(j)), (istr(k,mu),k=1,nstr(mu))
        end do
      end do

      return
      end

!     *************************************************************************

      subroutine serch (op,jrels,jsymt,xz,natm,na)

      use reimers_C, only: ef

      implicit none

      double precision  tol, xz(na,3)
      integer           i, j, k, na, jrels(na), jsymt(0:8), natm(na)

      logical           op, nm
      data tol /1.D-5/

!     **** determins if symm op is present, hence atoms related by symm ****

      do i= 1,na
        jrels(i)= 0
      end do
      op= .false.
      do k= 1,3
        if (jsymt(k).eq.-1 .and. ef(k).ne.0.0) return
      end do

      do 300 i= 1,na
        if (jrels(i).gt.0) goto 300

!	**** check for not moved by symm op ****
        nm= .true.
        do k= 1,3
          if (jsymt(k).eq.-1 .and. abs(xz(i,k)).gt.tol) nm= .false.
        end do

        if (.not.nm) then

!	  **** check for pairs of atoms related by symm ****
          do 200 j= 1,na
!	    **** check if symmetry pair ****
            if (j.eq.i) goto 200
            if (natm(i) .ne. natm(j)) goto 200
            do k= 1,3
              if (abs(xz(i,k)-xz(j,k)*jsymt(k)).gt.tol) goto 200
            end do
!	    **** matching atom found ****
            jrels(i)= j
            jrels(j)= i
            goto 300
200       continue

!	  **** only can get here if its moved and no matching atom ****
          op= .false.
          return
        end if

300   continue

      op= .true.
      return
      end

!     *************************************************************************

      subroutine swaap (xz,ef,op,jrels,i,j,na)

      implicit none
      double precision  xz(na,3), ef(3), tmp
      integer           i, j, k, &
                        itmp, jrels(na,7), na
      logical           op(7), op1

!     **** interchanges coordinate axes i and j ****

      op1= op(4-i)
      op(4-i)= op(4-j)
      op(4-j)= op1
      op1= op(8-i)
      op(8-i)= op(8-j)
      op(8-j)= op1

      tmp= ef(i)
      ef(i)= ef(j)
      ef(j)= tmp

      do k= 1,na
        tmp= xz(k,i)
        xz(k,i)= xz(k,j)
        xz(k,j)= tmp
        
        itmp= jrels(k,4-i)
        jrels(k,4-i)= jrels(k,4-j)
        jrels(k,4-j)= itmp
        itmp= jrels(k,8-i)
        jrels(k,8-i)= jrels(k,8-j)
        jrels(k,8-j)= itmp
      end do

      return
      end

!     *************************************************************************


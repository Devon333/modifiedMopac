! RMG - subroutine added to control Reimers CI in MOPAC
      subroutine rci ()

      USE chanel_C, only: iw
      USE reimers_C, only: beta, cc0, dd, n, nb2, nvh, nol, & 
          ndtype, dm, x, y, z, xz, zcore, ion, ndump, gamma,    &
          ci, evalmo, tot, nbtmo, iconf, aocc, bocc,   &
          spintr, cimatr, nstop, nconf, ic, id, dmci, ldcomp, &
          evalci, nptg, nspn, dtmp, nciout, nex, &
          nold, nvhd, nolc, nvhc, ncore, mrci, dbls
      USE molkst_C, only: norbs, keywrd, nopen, nclose
      USE overlaps_C, only: fact
      use reada_I

      implicit none

      integer :: i, j, occ, vir, occd, vird, occc, virc, noh, nvl, &
                 ndoubl, nmos, iy, nc2
      double precision, allocatable :: e(:)
      logical :: sing

      n = norbs
      nb2 = n*(n+1)/2
!     ****        transform one-electron matrix in active CI space ****
      call ao2mo  (beta,cc0,dd,norbs,n,1,n)
!     ****        transform dipole matrices in active CI space ****
      j = nol

!     ****	  get symmetry properties of gamma ****
      nptg = 1
      ion = 1
      call ptgrp(0)
      if(.not. allocated(ci)) allocate(ci(n,n))

      call stgamm (gamma,ci)

      if(.not. allocated(dm))     allocate(dm(nb2,3))

! Need to compute ground-state dipole to get proper dm matrix for CI
! dipoles 
      do i=1,nb2
        dd(i,1) = dtmp(i,1)
        dd(i,2) = dtmp(i,2)
      end do
      call dipol (x,y,z,dm)

      call gsdip  (xz,dd,zcore,dm)

      do i= 1,3
        if (ndtype.lt.4) call ao2mo (dm(1,i),cc0,dd,n,n,1,n)
      end do


! Read active space information
      noh = nopen
      nvl = nclose + 1


      if (index(keywrd,' MAXCI=').eq.0) then
        nconf = 2000
      else
        nconf = nint(reada(keywrd,index(keywrd,' MAXCI=')))
      end if
      if (index(keywrd,' WRTCI=').eq.0) then
        nciout = 500
      else
        nciout = nint(reada(keywrd,index(keywrd,' WRTCI=')))
      end if

! Single excitation active space
      occ = noh
      vir = norbs - nvl + 1

      ndoubl = 99
      if (index(keywrd,'C.I.=') /= 0) then
        i = index(keywrd,'C.I.=(')
        if (i /= 0) then
          j = index(keywrd(i:i+10),',') + i - 1
          ndoubl = nint(reada(keywrd,j))
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=(') + 5))
          ndoubl = min(ndoubl,nclose)
        else if (index(keywrd,'C.I.=') /= 0) then
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=') + 5))
        else
          nmos = nopen - nclose
        endif
        nmos = min(nmos,norbs)
        if (ndoubl == 99) then
          ndoubl = nmos/2
        endif
        occ = min(occ,ndoubl)
        vir = min(vir,nmos - ndoubl)
      end if
! Double excitation active space
      occd = 0
      vird = 0
      ndoubl = 99
      if (index(keywrd,'CISD') /= 0 .and. index(keywrd,'C.I.D.') .eq. 0) then
        ! Set default double excitation active space
        occd = min(5,occ)
        vird = min(5,vir)
      else if (index(keywrd,'C.I.D.') /= 0) then
        i = index(keywrd,'C.I.D.=(')
        if (i /= 0) then
          j = index(keywrd(i:i+10),',') + i - 1
          ndoubl = nint(reada(keywrd,j))
          nmos = nint(reada(keywrd,index(keywrd,'C.I.D.=(') + 7))
          ndoubl = min(ndoubl,nclose)
        else 
          nmos = nint(reada(keywrd,index(keywrd,'C.I.D.=') + 7))
        endif
        nmos = min(nmos,norbs)
        if (ndoubl == 99) then
          ndoubl = nmos/2
        endif
        occd = min(max(occd,ndoubl),noh)
        vird = min(max(vird,nmos - ndoubl),norbs - nvl + 1)
      end if
! CAS active space
      occc = 0
      virc = 0
      ndoubl = 99
      if (index(keywrd,'MRCI') /= 0 .and. index(keywrd,'C.A.S.') .eq. 0) then
        ! Set default CAS excitation active space
        occc = min(1,occ)
        virc = min(1,vir)
      else if (index(keywrd,'C.A.S.') /= 0) then
        i = index(keywrd,'C.A.S.=(')
        if (i /= 0) then
          j = index(keywrd(i:i+10),',') + i - 1
          ndoubl = nint(reada(keywrd,j))
          nmos = nint(reada(keywrd,index(keywrd,'C.A.S.=(') + 7))
          ndoubl = min(ndoubl,nclose)
        else 
          nmos = nint(reada(keywrd,index(keywrd,'C.A.S.=') + 7))
        endif
        nmos = min(nmos,norbs)
        if (ndoubl == 99) then
          ndoubl = nmos/2
        endif
        occc = min(max(occc,ndoubl),noh)
        virc = min(max(virc,nmos - ndoubl),norbs - nvl + 1)
      end if

      nol  = max(1,nvl-occ)
      nvh  = min(norbs,nvl+vir-1)
      nold = max(1,nvl-occd)
      nvhd = min(norbs,nvl+vird-1)
      nolc = max(1,nvl-occc)
      nvhc = min(norbs,nvl+virc-1)

      sing = .false.
      dbls = .false.
      mrci = .false.
      nex = 1

      if (index(keywrd,' CIS') /= 0 .or. index(keywrd,' C.I.') /= 0 &
          .or. index(keywrd,' MRCI') /= 0) then
        sing = .true.
        nex = (noh-nol+1)*(nvh-nvl+1)+1
      end if
! Doubles count
      if (index(keywrd,' CISD') /= 0 .or. index(keywrd,' C.I.D.') /= 0) then
        dbls = .true.
        nex = (noh-nol+1)*(nvh-nvl+1)*((noh-nold+1)*(nvhd-nvl+1)+2)+1
      end if
! MRCI count
      if (index(keywrd,' MRCI') /= 0 .or. index(keywrd,' C.A.S.') /= 0) then
        mrci = .true.
        nex = nex * fact((occc+virc)*2) /(fact(occc*2) * fact(virc*2))
        if (.not. sing) then
          nol = nolc
          nvh = nvhc
        end if
      end if
      ncore = nol-1
      if (nex .gt. 1.D7) then
        write (iw,'("Warning: Generating", i10," excitations")') nex
      end if

! Make sure nconf and nciout are not too large
      nconf = min(nex,nconf)
      nciout = min(nconf,nciout)

      allocate(cimatr(nconf*(nconf+1)/2))
      allocate(e(nconf))
      do i=1,nconf*(nconf+1)/2
        cimatr(i) = 0.D0
      end do
      if(.not. allocated(iconf))  allocate(iconf(nconf))
      do i=1,nconf
        iconf(i) = 0.D0
      end do
      call excite (dm,iconf,e,cimatr)

      if (nstop.eq.1 .or. nconf.eq.0) goto 3

!     ****	  construct dipole moment matrix in config-state basis ****
      deallocate(ci)
      allocate(dmci(nconf*(nconf+1)/2,3))
      allocate(ci(nconf,nconf))
      allocate(evalci(nconf))
      open        (ic,status='scratch',form='unformatted')
      write (ic)  ((cc0(i,j),j=1,n),i=1,n),(beta(i),i=1,nb2)
      rewind ic
      call foscil (dmci,dm,iconf,e,aocc,bocc,spintr,nspn,xz,zcore)

      if  (ndump.ge.3) then
!        write  (iy,*) 'DIPOLE MOMENTS IN SPIN-ADAPTED CONFIG BASIS'
        do i= 1,3
          call dump (dmci(1,i),nconf,0,ldcomp(i))
        end do
      end if
      nc2= nconf*(nconf+1)/2

      open        (id,status='scratch',form='unformatted')

      write (id)  ((dmci(i,j),i=1,nc2),j=1,3)
      rewind id
      read (ic)   ((cc0(i,j),j=1,n),i=1,n),(beta(i),i=1,nb2)
      rewind ic
      call cimat (cc0,gamma,beta,nbtmo, iconf,e,aocc,bocc,&
               spintr,nspn, cimatr, ndump)
      call cidiag (cimatr,ci,evalci,aocc,bocc)
!     ****	  trasfer dipole integrals to CI eigenvector basis ****
      read (id)   ((dmci(i,j),i=1,nc2),j=1,3)
      j= nciout

      if(allocated(gamma))       deallocate(gamma)
      allocate(gamma(nconf,nconf))

      do i= 1,3
        call ao2mo (dmci(1,i),ci,gamma,nconf,nconf,1,j)
      end do
      call ciout  (evalmo,ci,evalci,iconf,dmci,aocc,bocc,tot)

3     continue
      return
      end subroutine rci

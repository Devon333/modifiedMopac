      subroutine symoir(itype, vects, eigs, nvecs, r, imat) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE symmetry_C, ONLY: jndex, namo, nirred, nclass, name, jx, group
      USE meci_C, only : lab 
      use molkst_C, only : keywrd, numat
      use chanel_C, only : iw
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:03  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use charmo_I 
      use charvi_I 
      use charst_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: itype 
      integer  :: nvecs 
      integer , intent(in) :: imat  
      real(double)  :: vects(nvecs,*) 
      real(double) , intent(in) :: eigs(*) 
      real(double)  :: r(3,3) 
      real(double)  :: carmat(imat,20)  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(20) :: icount 
      integer , dimension(nvecs) :: ntype 
      integer :: korb, i, jj, j, nfind, ifound, l, ik, k 
      real(double), dimension(20) :: tchar 
      real(double) :: toler, check 
      logical :: first, r3, debug 
      character :: names*4 

      save toler 
!-----------------------------------------------
!***********************************************************************
!
!    Calculates the Irreducible Representation of an eigenvector.
!
!    GROUP(i,j) is the Character of Irreducible Representation i in
!    Class j.
!
!     On input  NVECS  = No. of vectors to be assigned.
!               VECTS  = Eigenvectors, packed
!
!     On output ICOUNT = No. of I.R.s of each type.
!               NAMO   = Irreducible Representation as a 4 character
!                        string.
!
!
!***********************************************************************
      data toler/ 0.1D0/  
!
!
!  ICOUNT will hold the number of occurances of each Irreducible
!         Representation
!
      r3 = name == 'R3' 
      korb = 0 
      do i = 1, numat 
        jj = jndex(i) 
        do j = 1, jj 
          korb = korb + 1 
          ntype(korb) = 100*i + 9 + j 
        end do 
      end do 
      if (itype /= 3) then 
        nfind = nvecs 
      else 
        nfind = lab 
      endif 
      icount(:nirred) = 0 
      names = '????' 
      if (nclass == 1) names = jx(1) 
      do i = 1, nfind 
        jndex(i) = i 
        namo(i) = names 
      end do 
      if (nclass == 1) return  
      ifound = 0 
!
!   Calculate all the characters in one pass.  This is done
!   for the convenience of CHARST.  In CHARST, the M.O. transform
!   only need be evaluated once per class.
!
      first = .TRUE. 
      do j = 1, nclass 
        l = (j - 1)*imat 
        do i = 1, nfind 
          l = l + 1 
          if (itype == 1) then 
            carmat(i,j) = charmo(vects, ntype, i, j, r, nvecs, first) 
          else if (itype == 2) then 
            carmat(i,j) = charvi(vects, i, j, r, nvecs) 
          else 
            carmat(i,j) = charst(vects, ntype, i, j, r, nvecs, first) 
          endif 
        end do 
      end do 
!
!     Reset counters in CHARST
!
      if (itype == 3) i = int(charst(vects ,ntype, -1, 0, r, nvecs, first))
      debug = index(keywrd,'SYMOIR') /= 0 
      if (debug) then 
        write (iw, '(A)') ' Characters of Transform' 
        do i = 1, nfind 
          write (iw, '(I5,6F12.6)') i, (carmat(i,j),j=1,nclass) 
        end do 
      endif 
      i = 0 
   70 continue 
      ik = i + 1 
      tchar(:nclass) = 0.D0 
   90 continue 
      i = i + 1 
      if (i > nfind) then 
!
!  Convert to lower case if one-electron eigenvector
!
        if (itype == 1) then 
          do i = 1, nfind 
            j = ichar(namo(i)(1:1)) 
            if (j>=ichar('A') .and. j<=ichar('Z')) j = j + ichar('a') - ichar(&
              'A') 
            namo(i)(1:1) = char(j) 
          end do 
        endif 
        if (debug) then 
          write (iw, '(/,A,/)') &
            ' Number of Irreducible Representations of each Class' 
          write (iw, '(7(I5,1X,A4))') (icount(i),jx(i),i=1,nirred) 
        endif 
        return  
      endif 
      tchar(:nclass) = tchar(:nclass) + carmat(i,:nclass) 
      if (tchar(1)>5.1D0 .and. .not.r3) go to 70 
      l140: do k = 1, nirred 
        do j = 1, nclass 
          check = abs(tchar(j)-group(k,j)) 
          if (check <= toler) cycle  
          cycle  l140 
        end do 
!
!   Degeneracy of Eigenvector manifold is I-IK+1
!
        icount(k) = icount(k) + 1 
!
!  The only information passed to MATOU1 is JNDEX and NAMO
!
        if (i - ik + 1 > 0) then 
          jndex(ik:i) = icount(k) 
          namo(ik:i) = jx(k) 
          ifound = i - ik + 1 + ifound 
        endif 
        go to 70 
      end do l140 
!
!   If manifold registration lost, use eigenvalues to decide
!   if start of new Irreducible Representation.
!
      if (i < nfind) then 
        if (eigs(i+1) - eigs(i) > 0.1D0) go to 70 
      endif 
      go to 90 
      end subroutine symoir 

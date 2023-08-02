      subroutine scfmat (v,shft)

!     ************************************************************************
!     **** generates the Fock matrices F for each shell plus determins the ***
!     **** variational energy using trace formula, V.			   ***
!     **** then it constracts the total fock operator in A,		   ***
!     **** using SHFT as a shift of the virtual energy levels in AII.	   ***
!     **** D= each shells density matrix (d,ff in symmetric storage mode)  ***
!     ****   D is destroyed if open shell or damped closed shell	   ***
!     **** P= atomic charges from each shell; BETA,U = 1 e matr; GAMMA= 2e ***
!     **** C= previous set of MO coefficients if first element .gt.-1.     ***
!     ************************************************************************
!     **** References:							   ***
!     ****   ROHF: Edwards and Zerner Theor Chim Acta 72 (1987) 347	   ***
!     ****	   eqn 27 for v, eqn 33 for shell Fock matrices		   ***
!     ****	   eqn 24 is solved in previous MO basis for the total	    **
!     ****	   Fock operator by putting				   ***
!     ****	     alam= n(nu)/n(mu)*LAMBDA(mu,nu) - LAMBDA(nu,mu)	   ***
!     ****	     so that cross term is alam*(F(nu)-F(mu))*n(mu)/n(nu)  ***
!     ****	     alam is arbitrary non-zero, is set to 1 herein	   ***
!     ****	   eqn 24 is extended by a VIRTUAL-VIRTUAL term designed   ***
!     ****	     to produce better virtual eigenvalues		   ***
!     ****	   Errors: eqn 24, no 1/n(mu); eqn 33 F(hat)= not F=	   ***
!     ****   INDO (Rotationally invariant for Transition Metals):	   ***
!     ****	   Bacon and Zerner Theor Chim Acta 53 (1979) 21	   ***
!     ****	   with additional R integrals included			   ***
!     ************************************************************************

      use molkst_C, only: mpack, norbs
      use reimers_C,only: n, na, nb2, matind, occfr, nel, nshell, norbl, &
          norbh, nbf, ibf, iat, ig1, ig2, ig3, ig4, g, nirreg, natt,&
          avec1, bvec1, avec2, bvec2, ppg, pg, dia, dd, ff, aa, cc0,&
          beta, gamma

      implicit none
      double precision   q, v, alamij, alamji, avec10, &
                         ga, gk, hg, irreg, off, shft, wt
      integer            i, j, k, &
                         ib, is, ishell, jb, jshell, &
                         k1, k2, k3, k4, mu, nm, nu

      real eltot
      real*8 dd1(mpack,2)

      v= 0.D0

!     **** diagonal term, sum Coulomb integrals from all other atoms ****
      if (.not. allocated(pg)) allocate(pg(na,0:2))
      if (.not. allocated(ppg)) allocate(ppg(norbs,0:2))

      do ishell= 0,nshell
        do ib= 1,n
          ppg(ib,ishell)= 0.D0
        end do
        do i= 1,na
          pg(i,ishell)= 0.D0
        end do
      end do

      do ib= 1,n
        i= iat(ib)

        do jb= 1,n
          j= iat(jb)
          ga= gamma(ib,jb)
          if (ib.lt.jb) ga= gamma(jb,ib)
          nm= matind(jb) + jb
          call veccou (nm,avec1,bvec1)
          do ishell= 0,nshell
            if (ib.eq.jb) then
!	      **** weighted charges ****
              pg(i,ishell)= pg(i,ishell) + avec1(ishell)
            else if (i.ne.j) then
!	      **** weighted sums over all other atoms ****
              ppg(ib,ishell)= ppg(ib,ishell) + avec1(ishell) * ga
            end if
          end do
        end do

      end do

!     ************ Shell Fock matrices for different Hamiltonians *************

!	**** INDO Hamiltonian ****
      do i= 1,na
        do mu= ibf(i),ibf(i)+nbf(i)-1
!	  **** terms in diagonal atom-atom block ****

!	  **** zero counters for the diagonal energy ****
          do is= 0,nshell
            dia(is)= 0.D0
          end do

!	  **** loop over all other densities on this atom ****
        do nu= ibf(i),ibf(i)+nbf(i)-1

!	    **** load the (nu,nu) elements of density matrix ****
            nm= matind(nu) + nu
            call veccou (nm,avec1,bvec1)
            if (mu.gt.nu) then
!	      **** load the (mu,nu) elements of density matrix ****
              nm= matind(mu) + nu
              call veccou (nm,avec2,bvec2)
            end if

            do is= 0,nshell
              if (mu.eq.nu) then
                dia(is)= dia(is) + gamma(mu,mu) * (avec1(is)-bvec1(is))
                avec10= avec1(0)
              else if (mu.gt.nu) then
!	        **** gamma(mu,nu)= (mu,mu|nu,nu); gamma(nu,mu)= (mu,nu|mu,nu) ***
                dia(is)= dia(is) + gamma(mu,nu)*avec1(is)&
     &                           - gamma(nu,mu)*bvec1(is)
!		**** same atom but different basis functions ****
                off= gamma(nu,mu) * (2.D0*avec2(is)-bvec2(is))&
     &             - gamma(mu,nu) * bvec2(is)
                if (is.eq.0) then
                  hg= off
                  v= v + off * avec2(0)
                else
                  ff(nm,is)= hg - off
                  v= v - off * dd(nm,is)
                end if
              else
                dia(is)= dia(is) + gamma(nu,mu)*avec1(is)&
     &                           - gamma(mu,nu)*bvec1(is)
              end if
            end do
          end do

!	  **** store fully diagonal terms ****
          nm= matind(mu) + mu
          hg= beta(nm) + ppg(mu,0) + dia(0)
          v= v + (beta(nm)+hg) * avec10 * 0.5D0
          do is= 1,nshell
            q= ppg(mu,is) + dia(is)
            ff(nm,is)= hg - q
            v= v - q * dd(nm,is) * 0.5D0
          end do

!	  **** off-diagonal atom-atom blocks and energy contribution ****
          do j= 1,i-1
            do nu= ibf(j),ibf(j)+nbf(j)-1
              nm= matind(mu) + nu
              call veccou (nm,avec1,bvec1)
              hg= beta(nm) - bvec1(0) * gamma(mu,nu)
              v= v + (beta(nm)+hg) * avec1(0)
              do is= 1,nshell
                q= - bvec1(is) * gamma(mu,nu)
                ff(nm,is)= hg - q
                v= v - q * dd(nm,is)
              end do
            end do
          end do
        end do

!	**** Zerners irregular integrals ****
        k= natt(i)
        if (k.gt.0) then
          mu= ibf(i)
          do irreg= 1,nirreg
            k1= ig1(irreg,k)
            k2= ig2(irreg,k)
            k3= ig3(irreg,k)
            k4= ig4(irreg,k)
            gk= g(irreg,k)
            call scfirr (v,mu, gk,k1,k2,k3,k4)
            call scfirr (v,mu, gk,k1,k2,k4,k3)
            call scfirr (v,mu, gk,k3,k4,k1,k2)
            call scfirr (v,mu, gk,k4,k3,k1,k2)
            if (k1.ne.k2) then
              call scfirr (v,mu, gk,k2,k1,k3,k4)
              call scfirr (v,mu, gk,k2,k1,k4,k3)
              call scfirr (v,mu, gk,k3,k4,k2,k1)
              call scfirr (v,mu, gk,k4,k3,k2,k1)
            end if
          end do
        end if

      end do

!     ********************* construct total fock matrix **********************

      if (nshell.eq.1 .and. shft.eq.0.D0) then
!	**** Closed shell, no damping so just copy single Fock matrix  ****
        do nm= 1,nb2
          aa(nm)= ff(nm,1)
        end do

      else
       dd1 = dd

!	***** first, zero all matrix elements of total fock matrix ****
        do nm= 1,nb2
          aa(nm)= 0.D0
        end do
!	**** transform Fock matrices into previous MO basis set ****
!	**** doing projections of E and Z eqn 24		****
!	**** uses first two density matrices as work space	****

!	**** these coupling coeffs are almost arbitrary, choose irrational ****
        alamij= sqrt(2.D0)
        alamji= sqrt(0.33D0)
        eltot= 0.D0
        do ishell= 1,nshell
          eltot= eltot + nel(ishell)
        end do

        do ishell= 1,nshell
          wt= 1.D0
!	  **** diagonal terms ISHELL - ISHELL ****
          call ao2mo1 (aa,ff(:,ishell),cc0,dd1, norbl(ishell),norbh(ishell),&
     &          norbl(ishell),norbh(ishell), wt)
!	  **** terms ISHELL - VIRTUAL ****
          call ao2mo1 (aa,ff(:,ishell),cc0,dd1, norbl(ishell),norbh(ishell),&
     &          norbh(nshell)+1,n, wt)
!	  **** JRRs xtra VIRTUAL-VIRTUAL term (improves order of virtuals) ****
          wt= nel(ishell) / eltot
          call ao2mo1 (aa,ff(:,ishell),cc0,dd1, norbh(nshell)+1,n,&
     &          norbh(nshell)+1,n, wt)

!	  **** off-diagonal occupied terms ISHELL - JSHELL ****
          do jshell= ishell+1,nshell
            wt= 1.D0
            call ao2mo1 (aa,ff(1,ishell),cc0,dd1, norbl(ishell),norbh(ishell),&
     &            norbl(jshell),norbh(jshell), wt)
            wt= - occfr(jshell) / occfr(ishell)
            call ao2mo1 (aa,ff(1,jshell),cc0,dd1, norbl(ishell),norbh(ishell),&
     &            norbl(jshell),norbh(jshell), wt)

          end do
        end do

!	**** shift virtual energies to dampen SCF convergence problems ****

        do is= 1,max(1,nshell-1)
          do i= norbh(is)+1,n
            nm= matind(i) + i
            aa(nm)= aa(nm) + shft
          end do
        end do

!	**** back-transform Fock matrix into AO basis ****

        call mo2ao (aa,cc0,dd1,n,n)

        if (cc0(1,1).eq.0.D0) then
          do nm= 1,nb2
            aa(nm)= ff(nm,1)
          end do
        end if
      end if
      return
      end

!     *************************************************************************

      subroutine scfirr (v,mu,g,k1,k2,k3,k4)

      use reimers_C,only: nshell, avec1, bvec1, avec2, bvec2, matind, ff, dd

      implicit none
      double precision   g, q, v, hg, v1
      integer            is, k1, k2, k3, k4, mu, nm

!     **** inserts Zerners irregular ints into atomic Fock matrix block ****
!     **** add and contribution to total energy v			 ****
!     **** mu= nber of first basis function of this atom		 ****

      if (k2.ge.k1) then
!	**** Coulomb integral, use A vec ****

!	**** mixing of dens(3,4) ****
        nm= matind(mu+max(k3,k4)) + mu+min(k3,k4)
        call veccou (nm,avec1,bvec1)

!	**** mixing of dens(1,2) ****
        nm= matind(mu+k2) + mu+k1
        call veccou (nm,avec2,bvec2)

        hg= avec1(0) * g
        v1= avec2(0) * hg
        do is= 1,nshell
          q= avec1(is) * g
          ff(nm,is)= ff(nm,is) + hg - q
          v1= v1 - q * dd(nm,is)
        end do
        if (k2.eq.k1) v1= 0.5D0 * v1
        v= v + v1
      end if

      if (k3.ge.k1) then
!	**** Exchange integral, use -B vec ****

!	**** mixing of dens(2,4) ****
        nm= matind(mu+max(k2,k4)) + mu+min(k2,k4)
        call veccou (nm,avec1,bvec1)

!	**** mixing of dens(3,1) ****
        nm= matind(mu+k3) + mu+k1
        call veccou (nm,avec2,bvec2)

        hg= - bvec1(0) * g
        v1= avec2(0) * hg
        nm= matind(mu+k3) + mu+k1
        do is= 1,nshell
          q= - bvec1(is) * g
          ff(nm,is)= ff(nm,is) + hg - q
          v1= v1 - q * dd(nm,is)
        end do
        if (k3.eq.k1) v1= 0.5D0 * v1
        v= v + v1
      end if

      return
      end

!     *************************************************************************

      subroutine veccou (nm,avec1,bvec1)

!     **** determins A and B coeffs for shells of Fock matrix using	****
!     **** E&Z eqn 30c; closed shell terms are set to zero.		****
!     **** nm= matrix element index, shell 0 is the total density	****
!     **** the B vec has the half factor of eqn 30b included		****
      use reimers_C,only: vca, vcb, nshell, dd

      implicit none
      double precision   avec1(0:4),bvec1(0:4)
      integer            i, ir, is, ishell, jshell, &
                         nm, nmoou1, nmoou2, notype

      avec1(0)= dd(nm,1)
      avec1(1)= 0.D0
      bvec1(1)= 0.D0
      do ishell= 2,nshell
        avec1(0)= avec1(0) + dd(nm,ishell)
        avec1(ishell)= 0.D0
        bvec1(ishell)= 0.D0
        do jshell= 2,nshell
          avec1(ishell)= avec1(ishell) + (1.D0-vca(ishell,jshell)) *&
     &                                dd(nm,jshell)
          bvec1(ishell)= bvec1(ishell) + (1.D0-vcb(ishell,jshell)) *&
     &                                dd(nm,jshell)
        end do
        bvec1(ishell)= bvec1(ishell) / 2.D0
      end do
      bvec1(0)= avec1(0) / 2.D0

      return
      end

!     *************************************************************************

! ____________________________________________________________________________
! ________________________ SCF OUTPUT ROUTINES _______________________________
! ____________________________________________________________________________

      subroutine output (c,aii)

      USE chanel_C, only : iw
      use reimers_C, only: n, na, nham, natm, nbf, nbt, nshell, norbl, norbh, &
          nmrep, nptg, nr, nsym, iatsym, kind, noh
      USE molkst_C, ONLY: nclose

      implicit none
      double precision  aii(n),c(n,n)
      character*1       lshell(n),ll(5)
      integer           notype(8,5),notot(8)
      integer           i, j, k, l, m, ii, ir, is, kst, &
                        nhomo, nlumo, nmoou1, nmoou2
 1000 format (1h1,10x,a80,5x,2(i2,'/'),i2/)
 1140 format (/' eigvals(Ev)',15f8.3)
 1141 format (11x,15(5x,a3))
 1142 format (' shell    ',15(7x,a1))
 1143 format (' % Frag',i3,15f8.1)
 1144 format (' Frag= ',15i8)
 1145 format (10x,15i8)
 1150 format (1x,2i4,a2,a3,15f8.4)
 1180 format (/' H(D)OMO=',f8.4,' LUMO=',f8.4,' Gap=',f8.4)
 1200 format (/' Summary of symmetry occupancy in each shell:'//&
     &	       ' SHELL  ',9(a3,1x))
 1210 format (4x,a1,1x,8i4)
 1211 format (1x,a5,8i4)

!     **** write number of each orb of each symmetry in each shell ****
!     **** get shell operator name ****

      nmoou1 = 1
      nmoou2 = n
      nr = 1
      nptg = 1

      if (nshell==1) then
        ll(1) = '1'
        ll(2) = 'v'
      else 
        ll(1) = '1'
        ll(2) = '2'
        ll(3) = 'v'
      endif

      if (.not. allocated(nsym)) allocate(nsym(n))
      do i=1,n
        nsym(i) = 1
      end do

      do is= 1,nshell+1
        do ir= 1,nr
          notype(ir,is)= 0
        end do
      end do

      do is= 1,nshell
        do i= norbl(is),norbh(is)
          lshell(i)= ll(is)
          notype(1,is)= notype(1,is) + 1
        end do
      end do
      is= nshell+1
      do i= norbh(nshell)+1,n
        lshell(i)= ll(is)
          notype(1,is)= notype(1,is) + 1
      end do
      do ir= 1,nr
        notot(ir)= 0
        do is= 1,nshell+1
          notot(ir)= notot(ir) + notype(ir,is)
        end do
      end do

      write (iw,1200) (nmrep(ir,nptg),ir=1,nr)
      do is= 1,nshell+1
        write (iw,1210) ll(is),(notype(ir,is),ir=1,nr)
      end do
      write (iw,1211) 'TOTAL',(notot(ir),ir=1,nr)

!     ***** page eject for MO dump ****

      nhomo= norbh(1)
      nlumo= norbh(nshell) + 1

      if (nmoou1.gt.nmoou2) return

!     **** first eigenvalues before main dump ****
 !     if (nmoou1.gt.1) then
 !     do 205 m= 1,nmoou1-1,15
 !       k= m+14
 !       if (k.gt.nmoou1-1) k= nmoou1-1
 !       write (iw,*)
 !       write (iw,1140) (aii(i),i=m,k)
 !       write (iw,1142) (lshell(i),i=m,k)
 !       write (iw,1141) (nmrep(nsym(i),nptg),i=m,k)
 !       write (iw,1145) (i,i=m,k)
!205     continue
!      end if

!     **** dump of molecular orbitals ****
      call matout(c, aii, n, n, n) 
  !    do 230 m= nmoou1,nmoou2,15
  !      k= m+14
  !      if (k.gt.nmoou2) k= nmoou2
  !      write (iw,*)
  !      write (iw,1140) (aii(i),i=m,k)
  !      write (iw,1142) (lshell(i),i=m,k)
  !      write (iw,1141) (nmrep(nsym(i),nptg),i=m,k)
  !      write (iw,1145) (i,i=m,k)
!        i= 0
!        do 225 j= 1,na
!          do 220 l= 1,nbf(j)
!            i= i+1
!220        write (iw,1150) i,j,iatsym(natm(j)),kind(nbt(i)+1),&
!     &                         (c(i,ii),ii=m,k)
!225       if (nham.ne.1 .and. nham.ne.3) write (iw,*)
!230     continue

!     **** last eigenvalues after main dump ****
!      if (k.lt.n) then
!      kst= k+1
!      do 250 m= kst,n,15
!        k= m+14
!        if (k.gt.n) k= n
!        write (iw,*)
!        write (iw,1140) (aii(i),i=m,k)
!        write (iw,1141) (nmrep(nsym(i),nptg),i=m,k)
!        write (iw,1145) (i,i=m,k)
!250     continue
!      end if

      return
      end

!     *************************************************************************

      subroutine gsdip (xz,d,zcore,dm)

!     **** calculates ground-state dipole moment and Mulliken charges ****
      use reimers_C, only: na, nb2, matind, debye, natm, nbf, ibf, nbt, &
          iatsym, qgs, dipgs, ndtype
      USE chanel_C, only : iw
      implicit none
      double precision  d(nb2), dm(nb2,3), xz(na,3), zcore(na), cm(3), &
                        dip0(3), pt(4), pm(3), dc, dp, dt
      integer           ind(0:8), i, j, k, ib, mu, nati, nm, nu

!     ****		identify AO as s, p, d(eg) or d(t2g) ****
      data ind          /1,2,2,2,3,3,4,4,4/

 1004 format (21x,'Dipole moment of ground state in Debye units'&
     & /10x,'contribution',9x,'moment',10x,'X',10x,'Y',10x,'Z'/13x,&
     & 'atomic',f18.6,2x,3f11.6/9x,'hybridization',f15.6,2x,3f11.6/13x,&
     & 'total',f19.6,2x,3f11.6)
 1165 format (/'   Atom  ---- Orbital Occupancy ----  Charge',&
     &         5x,'---- Local Dipole Moment ----' /&
     &         12x,'s       p  d(Eg) d(T2g)' )
 1175 format (i4,1x,a2,4f7.4,f10.6:3x,3f10.6)
 1176 format (7f10.6)

      ndtype = 1

      if(.not. allocated(qgs))   allocate(qgs(na))
      if(.not. allocated(dipgs)) allocate(dipgs(na,3))

!     **** zero counters for atomic and hybridization contributions ****
      do i= 1,3
        cm(i)= 0.D0
        pm(i)= 0.D0
      end do

      do i= 1,na
        nati= abs(natm(i))
        qgs(i)= zcore(i)
        dipgs(i,1)= 0.D0
        dipgs(i,2)= 0.D0
        dipgs(i,3)= 0.D0

!	**** break charge up into s, p, d(eg) and d(t2g) contributions ****
        do j= 1,4
          pt(j)= 0.D0
        end do
        do ib= ibf(i),ibf(i)+nbf(i)-1
          nm= matind(ib) + ib
          j= ind(nbt(ib))
          pt(j)= pt(j) + d(nm)
        end do

!	**** atomic charge and its contribution to dipole moment ****
        qgs(i)= qgs(i) - pt(1) - pt(2) - pt(3) - pt(4)
        do j= 1,3
          cm(j)= cm(j) + qgs(i) * xz(i,j)
        end do

        if (ndtype.ne.1 .or. nbf(i).le.1) then
!	  **** no hybridization contribution from this atom ****
        else
!	  **** determine hybridization contribution ****
          do mu= ibf(i),ibf(i)+nbf(i)-1
            do nu= ibf(i),mu-1
              nm= matind(mu) + nu
              do k= 1,3
                dipgs(i,k)= dipgs(i,k) + d(nm) * dm(nm,k)
              end do
            end do
          end do
          do j= 1,3
            dipgs(i,j)= dipgs(i,j) * debye * 2.D0
            pm(j)= pm(j) + dipgs(i,j)
          end do
        end if

      end do

      do j= 1,3
        cm(j)= cm(j) * debye
        dip0(j)= cm(j) + pm(j)
      end do
      dc= sqrt (cm(1)**2 + cm(2)**2 + cm(3)**2)
      dp= sqrt (pm(1)**2 + pm(2)**2 + pm(3)**2)
      dt= sqrt (dip0(1)**2 + dip0(2)**2 + dip0(3)**2)
      return
      end

! Compute and write bond orders
      subroutine boa (s,c)

      use vast_kind_param, only:  double
      use reimers_C, only: n, iat, na, matind, nb2, noh
      use common_arrays_C, only: nlast
      implicit none
      real(double), dimension(:), allocatable :: bo
      double precision   s(nb2), c(n,n), p
      integer            i, j, k, ia, ja, mind, na2, nind

      na2 = n*(n+1)/2

      allocate(bo(na2))

! Initialize bo to zero
      do i = 1,na2
        bo(i) = 0.D0
      end do


! Calculate bond orders
      do i = 1,n
        do j = i+1,n
          p = 0.D0
          do k = 1,noh
            p = p + 2.D0*c(i,k)*c(j,k)
          end do
          ia = iat(i)
          ja = iat(j)
          mind = matind(max(i,j)) + min(i,j)
          nind = matind(max(ia,ja)) + min(ia,ja)
          bo(nind) = bo(nind) + p**2
        end do
      end do

! Re-zero diagonals
      do i = 1,n
        mind = matind(i) + i
        bo(mind) = 0.D0
      end do

! Print bond orders
      write(26,*)' Bond orders'
      call vecprt(bo,na)

      end


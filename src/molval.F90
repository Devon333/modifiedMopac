      subroutine molval(c, p, rhfuhf) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE molkst_C, only : numat, norbs, mpack
      USE chanel_C, only : iw
      use common_arrays_C, only : nfirst, nlast
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:36:01  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision, dimension (norbs) :: val
      real(double) , intent(in) :: rhfuhf 
      real(double) , intent(in) :: c(norbs,norbs) 
      real(double) , intent(in) :: p(mpack) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  i, jj, jl, ju, j, kk, kl, ku, k, l1, l2, l 
      real(double) :: sum 
!----------------------------------------------- 
      do i = 1, norbs 
        sum = 0.D0 
        do jj = 1, numat 
          jl = nfirst(jj) 
          ju = nlast(jj) 
          do j = jl, ju 
            do kk = 1, numat 
              if (kk == jj) cycle  
              kl = nfirst(kk) 
              ku = nlast(kk) 
              do k = kl, ku 
                l1 = max(j,k) 
                l2 = j + k - l1 
                l = (l1*(l1 - 1))/2 + l2 
                sum = sum + c(j,i)*c(k,i)*p(l) 
              end do 
            end do 
          end do 
        end do 
        val(i) = sum*rhfuhf 
      end do 
      write (iw, '(10F8.4)') (val(i),i=1,norbs) 
      return  
      end subroutine molval 

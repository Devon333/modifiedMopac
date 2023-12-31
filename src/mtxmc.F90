      subroutine mtxmc(a, nar, b, nbr, c) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:27:20  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mxm_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nar 
      integer  :: nbr 
      real(double)  :: a(nbr,nar) 
      real(double)  :: b(nbr,nar) 
      real(double)  :: c(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, i 
!-----------------------------------------------
!     MATRIX PRODUCT C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR)
!     A AND B RECTANGULAR , PACKED,
!     C LOWER LEFT TRIANGLE ONLY, PACKED IN CANONICAL ORDER.
!  NOTE ... THIS IS THE BEST VERSION ON CRAY 1.
      l = 1

! TODO: GBR future modifications       

      do i = 1, nar 
        call mxm (a(1,i), 1, b, nbr, c(l), i) 
        l = l + i 
      end do 
      return  
      end subroutine mtxmc 

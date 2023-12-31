      subroutine set(s1, s2, na, nb, rab, ii) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE overlaps_C, only : isp, ips, sa, sb
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:35:52  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use aintgs_I 
      use bintgs_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: na 
      integer , intent(in) :: nb 
      integer , intent(in) :: ii 
      real(double) , intent(in) :: s1 
      real(double) , intent(in) :: s2 
      real(double) , intent(in) :: rab 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jcall 
      real(double) :: alpha, beta 
!-----------------------------------------------
!***********************************************************************
!
!     SET IS PART OF THE OVERLAP CALCULATION, CALLED BY OVERLP.
!         IT CALLS AINTGS AND BINTGS
!
!***********************************************************************
      if (na <= nb) then 
        isp = 1 
        ips = 2 
        sa = s1 
        sb = s2 
      else 
        isp = 2 
        ips = 1 
        sa = s2 
        sb = s1 
      endif 
      j = ii + 2 
      if (ii > 3) j = j - 1 
      alpha = 0.5D00*rab*(sa + sb) 
      beta = 0.5D00*rab*(sb - sa) 
      jcall = j - 1 
      call aintgs (alpha, jcall) 
      call bintgs (beta, jcall) 
      return  
!
      end subroutine set 
      subroutine aintgs(x, k) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE overlaps_C, only : a
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:27:09  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: k 
      real(double) , intent(in) :: x 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(double) :: c 
!-----------------------------------------------
!***********************************************************************
!
!    AINTGS FORMS THE "A" INTEGRALS FOR THE OVERLAP CALCULATION.
!
!***********************************************************************
      c = exp((-x)) 
      a(1) = c/x 
      do i = 1, k 
        a(i+1) = (a(i)*i+c)/x 
      end do 
      return  
!
      end subroutine aintgs 
      subroutine bintgs(x, k) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE overlaps_C, only : b, fact
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:20:17  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: k 
      real(double) , intent(in) :: x 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: io, last, i, m  
      real(double) :: absx, expx, expmx, y, xf 
!-----------------------------------------------
!**********************************************************************
!
!     BINTGS FORMS THE "B" INTEGRALS FOR THE OVERLAP CALCULATION.
!
!**********************************************************************
      io = 0 
      absx = abs(x) 
      if (absx > 3.D00) go to 40 
      if (absx > 2.D00) then 
        if (k <= 10) go to 40 
        last = 15 
        go to 60 
      endif 
      if (absx > 1.D00) then 
        if (k <= 7) go to 40 
        last = 12 
        go to 60 
      endif 
      if (absx > 0.5D00) then 
        if (k <= 5) go to 40 
        last = 7 
        go to 60 
      endif 
      if (absx <= 1.D-6) go to 90 
      last = 6 
      go to 60 
   40 continue 
      expx = exp(x) 
      expmx = 1.D00/expx 
      b(1) = (expx - expmx)/x 
      do i = 1, k 
        b(i+1) = (i*b(i)+(-1.D00)**i*expx-expmx)/x 
      end do 
      go to 110 
   60 continue 
      do i = io, k 
        y = 0.0D00 
        do m = io, last 
          xf = 1.0D00 
          if (m /= 0) xf = fact(m) 
          y = y + (-x)**m*(2*mod(m + i + 1,2))/(xf*(m + i + 1)) 
        end do 
        b(i+1) = y 
      end do 
      go to 110 
   90 continue 
      do i = io, k 
        b(i+1) = (2*mod(i + 1,2))/(i + 1.D0) 
      end do 
  110 continue 
      return  
!
      end subroutine bintgs 



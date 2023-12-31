      subroutine mxv(a, nar, vecx, nbr, vecy) 
      USE vast_kind_param, ONLY:  double 
      Use mod_vars_cuda, only: lgpu
#if GPU
      Use mod_vars_cuda, only: prec
      Use call_gemv_cublas
#endif      
      implicit none
      integer, parameter :: incy = 1, incx = 1
      integer  :: nar 
      integer  :: nbr 
      real(double)  :: a(nar,nbr) 
      real(double)  :: vecx(nbr) 
      real(double)  :: vecy(nar) 
!
!     RECTANGULAR MATRIX-VECTOR PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.

      if (lgpu .and. (nar >= 100 .or. nbr >= 100)) then      
#if GPU
         call gemv_cublas('N',nar,nbr,1.0_prec,a,nar,vecx,incx,0.0_prec,vecy,incy)
#endif
      else                  
         call dgemv ('N', nar, nbr, 1.0d0, a, nar, vecx, incx, 0.0d0, vecy, incy)
      endif 
      return  
      end subroutine mxv 


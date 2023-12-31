    Subroutine eigenvectors_LAPACK(eigenvecs, xmat, eigvals, ndim)  
      USE vast_kind_param, ONLY:  double 
      USE chanel_C, only : iw
      Use mod_vars_cuda, only: lgpu, ngpus
#if GPU
      Use mod_vars_cuda, only: prec
#endif
#if (MAGMA)
      Use magma
      use initMagma
#endif
      
      implicit none
      Integer :: ndim
      Real(double) :: eigenvecs(ndim,ndim), &
                    & eigvals(ndim),xmat((ndim*(ndim+1))/2)                  
      integer :: i, j
      Integer :: lwork, liwork, info 
      Real(double),dimension(1:10) :: work_tmp
      Integer, dimension(1:10) :: iwork_tmp
      Real(double), allocatable :: work(:)
      Integer, allocatable :: iwork(:)          
!==============================================================================   
! Code to find all eigenvectors and all eigenvalues for a symmetric General matrix 
! using LAPACK and MAGMA
! Gerd Bruno Rocha and Julio Carvalho Maia 11/17/2013.
!==============================================================================            
      continue    
      eigvals = 0.d0      
      eigenvecs = 0.d0        

!
! Perturb secular determinant matrix to split exact degeneracies
! (This is to get around a bug in the diagonalizer that causes eigenvectors to not be orthonormal)
!
!  The following unusual construction works.  Do NOT change it unless degeneracy tests have been done,
!  and the results are reproducible.
!
      forall (i=1:ndim)            
          xmat((i*(i + 1))/2) = xmat((i*(i + 1))/2) + i*1.d-10   ! Do NOT go much higher than 1.d-10, otherwise the geometry   
      endforall                                              ! optimization might go into an endless loop.
      call dtpttr( 'u', ndim, xmat, eigenvecs, ndim, i )
    
      if (lgpu .and. (ngpus > 1 .and. ndim > 100)) then
        call mkl_dimatcopy('C', 'T' , ndim, ndim, 1.0d0, eigenvecs, ndim, ndim)
      endif
      if (i /= 0) stop 'error in dtpttr'  
      
      j = i ! Dummy - to make FORCHECK not complain about "j"
      i = j ! Dummy  
      if (i == -999) return
      lwork = -1
      liwork = -1

! GBR_new_addition
    
#if (MAGMA)           
      if (lgpu .and. ndim > 100) then
         if (ngpus > 1) then
             call magma_dsyevd_Driver1(ngpus,'v','l',ndim,eigenvecs,ndim,eigvals,&
                    & work_tmp,lwork,iwork_tmp,liwork,info)
          else
             call magma_dsyevd_Driver1(ngpus,'v','u',ndim,eigenvecs,ndim,eigvals,&
                    & work_tmp,lwork,iwork_tmp,liwork,info) 
          endif
      else
         call dsyevd('v','u',ndim,eigenvecs,ndim,eigvals,work_tmp,&
                    & lwork,iwork_tmp,liwork,info)
      endif
#else
      call dsyevd('v','u',ndim,eigenvecs,ndim,eigvals,work_tmp, &
                    & lwork,iwork_tmp,liwork,info)              
#endif                      
      lwork = int(work_tmp(1))
      liwork = iwork_tmp(1)
      allocate (work(lwork), iwork(liwork), stat = i)      
      forall (j=1:lwork) work(j) = 0.d0      
      forall (j=1:liwork) iwork(j) = 0 
#if (MAGMA)           
      if (lgpu .and. ndim > 100) then
         if (ngpus > 1) then
             call magma_dsyevd_Driver2(ngpus,'v','l',ndim,eigenvecs,ndim,eigvals,&
                    & work,lwork,iwork,liwork,info)
          else
             call magma_dsyevd_Driver2(ngpus,'v','u',ndim,eigenvecs,ndim,eigvals,&
                    & work,lwork,iwork,liwork,info) 
          endif
      else
         call dsyevd('v','u',ndim,eigenvecs,ndim,eigvals,work,& 
                                & lwork,iwork,liwork,info)
      endif
#else
                    
      call dsyevd('v','u',ndim,eigenvecs,ndim,eigvals,work, &
                                & lwork,iwork,liwork,info)
#endif                      
                    
      deallocate (iwork,work,stat = j)                
      if (info /= 0)  write(iw,*) ' dsyevd Diagonalization error., CODE =',info                                           
      continue     
      return
    End subroutine eigenvectors_LAPACK

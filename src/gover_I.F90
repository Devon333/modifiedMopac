      MODULE gover_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.4G  10:47:18  03/09/06  
      SUBROUTINE gover (NI, NJ, XJ, R, SG, A0) 
      USE vast_kind_param,ONLY: DOUBLE 
      INTEGER, INTENT(IN) :: NI 
      INTEGER, INTENT(IN) :: NJ 
      REAL(DOUBLE), DIMENSION(3), INTENT(IN) :: XJ 
      REAL(DOUBLE), INTENT(INOUT) :: R 
      REAL(DOUBLE), DIMENSION(9,9), INTENT(INOUT) :: SG 
      REAL(DOUBLE), INTENT(IN) :: A0 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 

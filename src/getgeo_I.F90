      MODULE getgeo_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.4G  10:47:17  03/09/06  
      SUBROUTINE getgeo (IREAD, LABELS, GEO, XYZ, LOPT, NA, NB, NC, int) 
      USE vast_kind_param,ONLY: DOUBLE 
      INTEGER, INTENT(IN) :: IREAD 
      INTEGER, DIMENSION(*), INTENT(OUT) :: LABELS 
      REAL(DOUBLE), DIMENSION(3,*), INTENT(OUT) :: GEO, XYZ
      INTEGER, DIMENSION(3,*), INTENT(OUT) :: LOPT 
      INTEGER, DIMENSION(*), INTENT(OUT) :: NA 
      INTEGER, DIMENSION(*), INTENT(OUT) :: NB 
      INTEGER, DIMENSION(*), INTENT(OUT) :: NC 
      logical, intent(out) :: int
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 

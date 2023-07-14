SUBROUTINE RS_SOIL(NDIM,WG,VEG,RSOIL)
!--------------------------------------------------------------------------
!
! Computation of a soil surface resistance according
! to ECMWF formulation
! Replaces the Hu formulation from ISBA (to avoid dew flux pb)
! Minimum threshold modified (wwilt -> veg*wwilt) - Albergel et al. (2012)
!
!                                            Jean-Francois MAHFOUF (11/06)
!                                            Modified (10/21)          
!--------------------------------------------------------------------------
 USE SOIL
 IMPLICIT NONE
 INTEGER, INTENT(IN)                :: NDIM          
 REAL, INTENT(IN), DIMENSION(NDIM)  :: WG, VEG
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RSOIL
 REAL, DIMENSION(NDIM) :: ZF2
!
 ZF2 = (WG - VEG*WWILT)/(WFC - VEG*WWILT)
 ZF2 = MIN(1.,MAX(0.0001,ZF2))
 RSOIL = 50./ZF2
 RETURN
END SUBROUTINE RS_SOIL

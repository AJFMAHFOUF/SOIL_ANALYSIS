SUBROUTINE SOIL_DEF(NDIM,CLAY,SAND)
!--------------------------------------------------------
!
! Initialisation of soil thermal and hydraulic properties
! (part 1)
!
!--------------------------------------------------------
 USE SOIL
 IMPLICIT NONE
 INTEGER, INTENT (IN) :: NDIM
 REAL, INTENT(IN), DIMENSION(NDIM) :: CLAY
 REAL, INTENT(IN), DIMENSION(NDIM) :: SAND
!
 ALLOCATE (WSAT(NDIM))
 ALLOCATE (WWILT(NDIM))
 ALLOCATE (WFC(NDIM))
 ALLOCATE (B(NDIM))
 ALLOCATE (CGSAT(NDIM))
 ALLOCATE (C1SAT(NDIM))
 ALLOCATE (C2REF(NDIM))
 ALLOCATE (C3(NDIM))
 ALLOCATE (A(NDIM))
 ALLOCATE (P(NDIM))
 ALLOCATE (WL(NDIM))
 ALLOCATE (C1(NDIM))
 ALLOCATE (C2(NDIM))
 ALLOCATE (CG(NDIM))
 ALLOCATE (WGEQ(NDIM))
!
 WSAT = (-108.*SAND + 494.305)*1.E-3
 WWILT = 37.1342E-3*SQRT(CLAY*100.)
 WFC = 89.0467E-3*(CLAY*100.)**(0.3496)
 B = 13.7 * CLAY + 3.501
 CGSAT = (-1.557*SAND - 1.441*CLAY + 4.7021)*1.0E-6 
 C1SAT = 5.58*CLAY + 0.8488
 C2REF = 13.815*(CLAY*100.)**(-0.954) 
 C3 = 5.327*(CLAY*100.)**(-1.043)
 A = 0.73242*(CLAY*100.)**(-0.539)
 P = 13.4*CLAY + 3.4
 WL = 1.E-5
!
END SUBROUTINE SOIL_DEF

SUBROUTINE INIT_SURF1(NDIM)
!---------------------------------------------------------------
!
! Definition of constant surface properties 
! Addition of maximum storage of interception reservoir WRMAX
!
!                                Jean-Francois MAHFOUF (11/06)
!                                                      (10/21)
!---------------------------------------------------------------
 USE SURF1
 IMPLICIT NONE
 INTEGER, INTENT (IN) :: NDIM
 ALLOCATE (ZREF(NDIM))
 ALLOCATE (Z0(NDIM))
 ALLOCATE (Z0H(NDIM))
 ALLOCATE (EMIS(NDIM))
 ALLOCATE (ALPHA(NDIM))
 ALLOCATE (RSMIN(NDIM))
 ALLOCATE (LAI(NDIM))
 ALLOCATE (D1(NDIM))
 ALLOCATE (D2(NDIM))
 ALLOCATE (GAMMA(NDIM))
 ALLOCATE (RGL(NDIM))
 ALLOCATE (CV(NDIM))
 ALLOCATE (CT(NDIM))
 ALLOCATE (VEG(NDIM))
 ALLOCATE (WRMAX(NDIM))
 ZREF = 2.0         ! reference level of the atmospheric forcing
 Z0 = 0.01           ! surface roughness length
 Z0H = 0.001         ! surface roughness length for heat
 EMIS = 0.97         ! surface emissivity
 ALPHA = 0.2         ! surface albedo
 RSMIN = 40.         ! minimum stomatal resistance
 LAI = 1.            ! leaf area index
 D1 = 0.01           ! depth of surface soil layer
 D2 = 0.50           ! depth of deep soil layer
 GAMMA = 20.         ! dependency of RS with saturation vapor deficit
 RGL = 100.          ! dependency of RS with solar radiation
 CV = 2.E-5          ! vegetation thermal coefficient 
 VEG = 0.00          ! vegetation fractionnal cover
 WRMAX = 0.2*VEG*LAI ! Maximum capacity of interception reservoir
END SUBROUTINE INIT_SURF1

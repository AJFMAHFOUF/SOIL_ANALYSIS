SUBROUTINE BIBLE(NDIM,WG,SAND,CLAY,TA,TS,T2,TB) 
!-------------------------------------------------------------
!
! BIBLE : BIosphere Bare soil L-band Emissivity
! Computation of microwave brightness temperatures
! at the top of the atmosphere as seen by an L-band radiometer
! Mostly based on LSMEM and L-MEB descriptions
!
!                                Jean-Francois MAHFOUF (06/07)   
!-------------------------------------------------------------
 USE SURF1
 USE SOIL
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 REAL, INTENT(IN),  DIMENSION(NDIM) :: WG, SAND, CLAY, TA, TS, T2
 REAL, INTENT(OUT), DIMENSION(NDIM,2) :: TB
 REAL, DIMENSION(NDIM) :: T_AU, T_AD, TAU_AT, T_SKY, TAU_VEG, &
&                         ZALT, ZB, ZWC, ZOMEGA, ZTV, ZEPS_W0, &
&                         ZTAU_W, ZBETA, ZSIG_EFF,ZTSC, Z1, Z2, &
&                         ZGAM_AT, ZGAM_VEG, ZZZ, ZTS, ZCT, ZEPS_S, ZEPS_WINF, &
&                         ZEPSA, ZTAUWB
 REAL, DIMENSION(NDIM,2) :: ZTB_VEG, ZTB_SOIL, ZRS, ZEMIS
 COMPLEX, DIMENSION (NDIM) :: ZG, ZEPS_B, ZEPS_W
 INTEGER :: IP
 COMPLEX :: J
 REAL :: ZTHETA, ZPI, ZF, ZC, ZSIGMA, ZEPS_0, ZALF, ZCOS, ZSIN, ZSAL 
!
! Physical constants
!
 J = (0.,-1.)
 ZPI = ACOS(-1.0)
 ZC = 2.998E+8      ! speed of light
 ZEPS_0 = 8.854E-12 ! permittivity of free space (Fm-1)
 ZEPS_S = 4.7       ! solid soil dielectric constant
 ZEPS_WINF = 4.9    ! high frequency water dielectric constant
 ZSAL = 0.00        ! soil water salinity (in psu)
!
! Radiometer characteristics
!
 ZF = 1.4E+9 ! frequency in Hz
 ZTHETA = (30.0)*ZPI/180. ! incidence angle
 ZCOS = COS(ZTHETA) ; ZSIN = SIN(ZTHETA)
!
! Surface altitude in km
!
 ZALT = 0.0450 
!
! Optical properties of vegetation
!
 ZB = 0.15     ! attenuation factor at L-band
 ZWC = 0.5*LAI ! vegetation water content
 ZOMEGA = 0.05 ! single scattering albedo
!
! Roughness parameter (in cm)
!
 ZSIGMA = 1.2
!
! Factor for combination of mixed dielectric constants
! 
 ZALF = 0.65
!
! Factor for estimating "effective" soil temperature
!
 ZCT = MIN(1.0,(WG/0.3)**0.3) ! L-MEB default values
 ZTS = T2 + ZCT*(TS - T2)
!
! Empirical atmospheric contributions in L-band (Pellarin et al.,2003)
!
 T_AU = EXP(4.9274 + 0.002195*TA)
 TAU_AT = EXP(-3.9262 - 0.2211*ZALT - 0.00369*TA)
 T_AU = T_AU*(1.0 - EXP(-TAU_AT/ZCOS))
 T_AD = T_AU
 T_SKY = 2.7*EXP(-TAU_AT/ZCOS)
!
! Vegetation temperature
!
 ZTV = TS
!
! Vegetation optical depth
!
 TAU_VEG = ZB*ZWC/ZCOS
!
! Dielectric constant of pure water (Dobson et al., 1985) + salinity (Klein and Swift, 1977)
!
 ZTSC = ZTS - 273.15
 ZEPS_W0 = 87.134 - 1.949E-1*ZTSC - 1.276E-2*ZTSC*ZTSC + 2.491E-4*ZTSC*ZTSC*ZTSC 
 !ZEPS_W0 = 88.045 - 0.4147*ZTSC + 6.295E-4*ZTSC*ZTSC + 1.075E-5*ZTSC*ZTSC*ZTSC
 ZEPSA = 1.0 + 1.613E-5*ZTSC*ZSAL - 3.656E-3*ZSAL - 3.21E-5*ZSAL*ZSAL - 4.232E-7*ZSAL*ZSAL*ZSAL
 ZEPS_W0 = ZEPS_W0*ZEPSA
 ZTAU_W = 1.1109E-10 - 3.824E-12*ZTSC + 6.938E-14*ZTSC*ZTSC - 5.096E-16*ZTSC*ZTSC*ZTSC
 ZTAUWB = 1.0 + 2.282E-5*ZTSC*ZSAL - 7.638E-4*ZSAL - 7.76E-6*ZSAL*ZSAL + 1.105E-8*ZSAL*ZSAL*ZSAL
 ZTAU_W = ZTAU_W*ZTAUWB
 ZSIG_EFF = 3.493 - 5.14*WSAT - 2.013*SAND + 1.594*CLAY ! CLAY and SAND between 0 and 1
 WHERE (WG > 0.0) 
   ZEPS_W = ZEPS_WINF + (ZEPS_W0 - ZEPS_WINF)/(1 + J*ZF*ZTAU_W) - &
&           J*ZSIG_EFF/(2.*ZPI*ZF*ZEPS_W0)*WSAT/MAX(WG,0.001)
 ELSEWHERE (WG <= 0.0)  
   ZEPS_W = ZEPS_WINF + (ZEPS_W0 - ZEPS_WINF)/(1 + J*ZF*ZTAU_W) 
 END WHERE
!
! Dielectric constant of soil/water medium
! 
 ZBETA = 1.33797 - 0.606*SAND - 0.166*CLAY
 !ZBETA = 1.09 - 0.11*SAND + 0.18*CLAY
 ZEPS_B = ((1.0 - WSAT)*(ZEPS_S)**ZALF + WSAT - WG + WG**ZBETA*ZEPS_W**ZALF)**(1./ZALF) 
 Z1 = REAL (ZEPS_B)
 !ZBETA = 1.2748  - 0.519*SAND - 0.152*CLAY
 ZBETA = 1.09 - 0.11*SAND + 0.18*CLAY
 ZEPS_B = ((1.0 - WSAT)*(ZEPS_S)**ZALF + WSAT - WG +  WG**ZBETA*ZEPS_W**ZALF)**(1./ZALF) 
 Z2 = AIMAG (ZEPS_B)
 ZEPS_B = CMPLX(Z1,Z2)
!
! Fresnel reflection coefficients : (1) for H and (2) for V polarisation
!
 ZG = SQRT(ZEPS_B - ZSIN**2)
 ZRS(:,1) = ABS((ZCOS - ZG)/(ZCOS + ZG))**2 
 ZRS(:,2) = ABS((ZCOS*ZEPS_B - ZG)/(ZCOS*ZEPS_B + ZG))**2
! 
! Empirical roughness dependency (Wang and Schmugge 1980)
!
 ZRS = ZRS ! * EXP(-(4.0*ZPI*ZF*ZSIGMA*ZCOS/(ZC*100.))**2)
!
! Surface emissivity in H and V polarizations
!
 ZEMIS = 1.0 - ZRS
!
! Brightness temperatures at TOA in H and V polarisation
! (Kerr and Njoku, 1990)
!
 ZGAM_AT  = EXP(-TAU_AT)
 ZGAM_VEG = EXP(-TAU_VEG)
 ZZZ = T_AD + T_SKY*ZGAM_AT
 DO IP = 1,2
   ZTB_VEG(:,IP)  = T_AU + & 
&                   ZGAM_AT*ZZZ*(1.0 - ZEMIS(:,IP))*ZGAM_VEG**2 + &
&                   ZGAM_AT*ZEMIS(:,IP)*ZTS*ZGAM_VEG + &
&                   ZGAM_AT*ZTV*(1.0 - ZOMEGA)*(1.0 - ZGAM_VEG)* &
&                   (1.0 + (1.0 - ZEMIS(:,IP))*ZGAM_VEG)
   ZTB_SOIL(:,IP) = T_AU + ZGAM_AT*ZZZ*(1.0 - ZEMIS(:,IP)) +  &
&                   ZGAM_AT*ZEMIS(:,IP)*ZTS
   TB(:,IP) = ZTB_VEG(:,IP)*VEG + (1.0 - VEG)*ZTB_SOIL(:,IP)
 ENDDO  
END SUBROUTINE BIBLE

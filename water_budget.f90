SUBROUTINE WATER_BUDGET(NDIM,DT,WG,W2,WR,WG_n,W2_n,PR,LEG,LEV,LETR,WGN,W2N,WRN, &
&                       WGN_n,W2N_n,RO) 
!--------------------------------------------------------------------
!
! Solve (explicitely) the force-restore equations for Ws, W2 and Wr
!
!
! Modified to include model error (red noise) for EnKF option
! Modified to account the evolution of the interception reservoir
! and the corresponding evaporation fluxes and runoff components
!
!
!                                 Jean-Francois MAHFOUF (11/06)
!                                                       (10/21)
!--------------------------------------------------------------------
 USE SOIL
 USE CONST
 USE SURF1
 USE SETUP
 IMPLICIT NONE
 INTEGER, INTENT(IN)                       :: NDIM
 REAL,    INTENT(IN)                       :: DT
 REAL,    INTENT(IN),     DIMENSION(NDIM)  :: PR, LEG, LEV, LETR
 REAL,    INTENT(IN),     DIMENSION(NDIM)  :: WG, W2, WR, WG_n, W2_n
 REAL,    INTENT(OUT),    DIMENSION(NDIM)  :: WGN, W2N, WRN, WGN_n, W2N_n, RO
 REAL,                    DIMENSION(NDIM)  :: RUNOFF1, RUNOFF2, RUNOFF3, RUNOFF4
 REAL, DIMENSION (NDIM) :: ZMU
 REAL :: ZZ, ZSIGMA, ZTAU, ZALPHA, ZBETA
 INTEGER :: I
!
! Define parameters for model error
!
 ZTAU   = 3.*86400.    ! temporal correlation (sec)
 ZSIGMA = 1.E-3/86400. ! standard deviation of error (m3/m3/s)
 ZALPHA = 1./(1. + DT/ZTAU)
 IF (L_ENKF) THEN
  ZBETA = 1.0
 ELSE
  ZBETA = 0.0
 ENDIF
!
! Evolution of the interception reservoir
!
 WRN = WR + DT*(VEG*PR - LEV/LV - LETR/LV) 
!
 WHERE (WRN > WRMAX)
  RUNOFF4 = (WRN - WRMAX)/DT
  WRN = WRMAX
 ELSEWHERE (WRN <= WRMAX)
  RUNOFF4 = 0.
 END WHERE 
 WHERE (WRN < 0.)
  RUNOFF4 = (0.0 - WR)/DT ! conserve water by getting the missing amount in the soil
  WRN = 0.0
 END WHERE
!  
 WRN = MAX(0.,WRN)
!
! Evolution of the surface volumetric water content
!
 WGN = (WG +  DT*(C1*((1.-VEG)*PR + RUNOFF4 - LEG/LV)/RHOW + C2*WGEQ/TAU))/(1. + C2*DT/TAU)
!
! Add model error
!
 DO I = 1, NDIM
  CALL GASDEV(ZZ)
  ZMU(I) = ZSIGMA*SQRT(1. - ZALPHA**2)*ZZ
 ENDDO 
 WGN_n = ZALPHA*WG_n + ZMU
!
 WGN = WGN + ZBETA*WGN_n*DT
!
 WGN = MAX(WL,WGN)
 WHERE (WGN > WSAT)
  RUNOFF1 = (WGN - WSAT)*D1*1000.0
  WGN = WSAT
 ELSEWHERE (WGN < WSAT)
  RUNOFF1 = 0.
 END WHERE 
!
! Evolution of mean soil moisture content
!
 W2N = W2 + DT*((1.-VEG)*PR + RUNOFF4 - LEG/LV - LETR/LV)/(D2*RHOW) - DT*C3/TAU*MAX(0.,W2 - WFC)  
!
! Add model error
!
 DO I = 1, NDIM
  CALL GASDEV(ZZ)
  ZMU(I) = ZSIGMA*SQRT(1. - ZALPHA**2)*ZZ
 ENDDO 
 W2N_n = ZALPHA*W2_n  + ZMU
!
 W2N = W2N + ZBETA*W2N_n*DT
!---------------------------------------
!write (400,*) W2N_n(1)*86400.*1000.0
!---------------------------------------
 RUNOFF3 = DT*C3/TAU*MAX(0.,W2 - WFC)*D2*1000.0
 W2N = MAX(WL,W2N)
 WHERE (W2N > WSAT)
  RUNOFF2 = (W2N - WSAT)*D2*1000.0
  W2N = WSAT
 ELSEWHERE (W2N < WSAT)
  RUNOFF2 = 0.
 END WHERE  
 RO = RUNOFF1 + RUNOFF2 + RUNOFF3 
END SUBROUTINE WATER_BUDGET

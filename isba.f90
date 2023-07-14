SUBROUTINE ISBA(NDIM,XI,DT,T_LENGTH,T_DEB,LPRINT,XF,YF)
 USE SURF1
 USE SOIL
 USE FORC
 USE SETUP
 IMPLICIT NONE
 INTERFACE
  REAL FUNCTION QSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION QSAT 
 END INTERFACE
 INTEGER, INTENT(IN) :: NDIM
 INTEGER, PARAMETER :: NVAR=4, NOBS=4
 REAL, INTENT(IN),  DIMENSION(NDIM,NVAR) :: XI
 REAL, INTENT(IN) :: DT
 REAL, INTENT(IN) :: T_LENGTH, T_DEB
 LOGICAL, INTENT(IN) :: LPRINT
 REAL, INTENT(OUT), DIMENSION(NDIM,NVAR) :: XF 
 REAL, INTENT(OUT), DIMENSION(NDIM,NOBS) :: YF
 REAL, DIMENSION (NDIM) ::  RG, RL, PR, TA, UA, VA, PS, QA, QS,         &
&                           CLAY, SAND, WG, W2, TS, T2, RSOIL, RS, RA,  &
&                           TSN, T2N, H, LEV, LEG, LE, RN, G, WGN, W2N, &
&                           U10M, V10M, T2M, Q2M, RH2M, RO, QG, Z1, Z2, &
&                           USTAR, HU, LETR, WR, WRN, DELTA
 REAL, DIMENSION (NDIM) :: WG_n, W2_n, WGN_n, W2N_n
 REAL, DIMENSION (NDIM,2) :: TB
 REAL :: ZRG, ZPR
 INTEGER :: NSTEP, I, II
 OPEN (UNIT=50,FILE='INPUT_OI.DAT')
 ZRG = 0. ; ZPR = 0.
 CALL INIT_SURF1(NDIM)
 CLAY = 0.18 ; SAND = 0.11
 CALL SOIL_DEF(NDIM,CLAY,SAND)
 NSTEP_DEB = INT(T_DEB/DT) + 1
 NSTEP_TOT = NSTEP_DEB + INT(T_LENGTH/DT) - 1
! print *,'NSTEP_DEB=',NSTEP_DEB
! print *,'NSTEP_TOT=',NSTEP_TOT
! CALL READ_FORCING
! print *,'*** ISBA for time slot between ',NSTEP_DEB,' and ',NSTEP_TOT
 DO NSTEP = NSTEP_DEB,NSTEP_TOT
!
  IF (NSTEP == NSTEP_DEB) THEN
 !   PRINT *,'========= SOIL PROPERTIES =============='
 !   PRINT *,'WILTING POINT ',WWILT(1)
 !   PRINT *,'FIELD CAPACITY ',WFC(1)
 !   PRINT *,'SATURATION',WSAT(1)
 !   PRINT *,'========================================'
    WG = WWILT + XI(:,1)*(WFC - WWILT)
    W2 = WWILT + XI(:,2)*(WFC - WWILT)
    WHERE (WG > WSAT) 
     WG = WSAT
    ELSEWHERE 
     WG = MAX(WL,WG)
    ENDWHERE
    WHERE (W2 > WSAT) 
     W2 = WSAT
    ELSEWHERE 
     W2 = MAX(WL,W2)
    END WHERE
    TS = XI(:,3)
    T2 = XI(:,4)
    WR = 0.0
  ENDIF
!
! Set-up noise (model error)
!
  IF (NSTEP == 1) THEN
    WG_n = 0.
    W2_n = 0.
    PS(:) = PSF(1)
    DO I = 1,NDIM
      QG(I) = QSAT(PS(I),TS(I))
    ENDDO
  ENDIF
!
  if (nstep == nstep_deb .and. nstep /= 1) then
    rewind (97)
    do i = 1,ndim
      read (97,*) WG_n(i), W2_n(i)
    enddo
    rewind (98)
    do i = 1,ndim
      read (98,*) QG(i)
    enddo
  endif
!  
! Define the fraction of wet vegetation 
!
  WHERE (WRMAX > 0)
    DELTA = (WR/WRMAX)**(2./3.)
  ELSEWHERE
    DELTA = 0.0
  END WHERE  
!
  CALL SOIL_PROP(NDIM,SAND,WG,W2,TS)
  CALL INTERPOL_FORCING(NDIM,DT,NSTEP,RG,RL,PR,TA,UA,VA,PS,QA)
  CALL RS_SOIL(NDIM,WG,VEG,RSOIL)
  CALL RS_VEG(NDIM,W2,PS,QA,TA,RG,RS)
  CALL DRAG_COEFF_Z0H(NDIM,TA,TS,QA,QG,UA,VA,RA,USTAR)
  CALL ENERGY_BUDGET(NDIM,DT,TS,T2,RA,RS,RSOIL,DELTA,PS,TA,QA,RG,RL,WG,TSN,T2N,HU)
  DO I = 1, NDIM
   QS(I) = QSAT(PS(I),TSN(I))
  ENDDO
  CALL FLUXES(NDIM,TSN,T2N,TA,QS,QA,PS,RG,RL,RA,RS,RSOIL,DELTA,HU,H,LEV,LEG,LETR,LE,RN,G)
  CALL WATER_BUDGET(NDIM,DT,WG,W2,WR,WG_n,W2_n,PR,LEG,LEV,LETR,WGN,W2N,WRN,WGN_n,W2N_n,RO) 
  WG = WGN
  W2 = W2N
  WR = WRN
  TS = TSN
  T2 = T2N
  WG_n = WGN_n
  W2_n = W2N_n
  Z1 = RA/(RA + RS) ; Z2 = RA/(RA + RSOIL)
  QG = QA + VEG*DELTA*(QS - QA) + VEG*(1. - DELTA)*(QS -QA)*Z1 + (1. -VEG)*(HU*QS - QA)*Z2
  CALL VDFPPCFLS(NDIM,PS,TA,TS,UA,VA,QA,QG,U10M,V10M,T2M,Q2M,RH2M)
  CALL BIBLE(NDIM,WG,SAND,CLAY,TA,TS,T2,TB)
!
  IF (NSTEP == NSTEP_TOT) THEN
! model counterpart of observations 
   YF(:,1) = T2M
   YF(:,2) = RH2M 
!   YF(:,3) = (WG - WWILT)/(WFC - WWILT)
   YF(:,3) = TB(:,1)
   YF(:,4) = TB(:,2)
! final control vector
   XF(:,1) = (WG - WWILT) / (WFC - WWILT)
   XF(:,2) = (W2 - WWILT) / (WFC - WWILT)
   XF(:,3) = TS
   XF(:,4) = T2
! save variables for next time slot
   rewind (97)
    do i = 1,ndim
      write (97,*) WG_n(i), W2_n(i)
    enddo
    rewind (98)
    do i = 1,ndim
      write (98,*) QG(i)
    enddo
  ENDIF
!
  ZRG = ZRG + RG(1)*DT
  ZPR = ZPR + PR(1)*DT 
  IF (LPRINT) THEN
   IF (L_ENKF) THEN
     CALL PRINT_OUTPUT_ENKF(NDIM,NSTEP,DT,T2M,RH2M,TS,T2,WG,W2,TB, &
&                           H,LE,RN,G,PR,RO)
   ELSE
     CALL PRINT_OUTPUT_OI  (NDIM,NSTEP,DT,T2M,RH2M,TS,T2,WG,W2,TB, &
&                           H,LE,RN,G,PR,RO)
   ENDIF
  ENDIF
 ENDDO
 IF (L_OI) THEN
   WRITE (50,*) ZRG,ZPR,ALPHA(1),VEG(1),TA(1),SQRT(UA(1)*UA(1)+VA(1)*VA(1)), &
&               LAI(1),RSMIN(1),WWILT(1),WFC(1),D1(1),D2(1)
 ENDIF
 CLOSE (UNIT=50)
 CALL CLEAN_EXIT
RETURN
END SUBROUTINE ISBA
 

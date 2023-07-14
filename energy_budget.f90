SUBROUTINE ENERGY_BUDGET(NDIM,DT,TS,T2,RA,RS,RSOIL,DELTA,PS,TA,QA,RG,RL,WG,TSN,T2N,HU)
!-------------------------------------------------------------------------
!
! Solve (implicitely) the force-restore equations for Ts and T2
! Modified to account for Hu formulation of bare soil evaporation 
! and evaporation from interception reservoir
!
!
!                                         Jean-Francois MAHFOUF (11/06)
!                                                               (10/21)
!-------------------------------------------------------------------------
 USE SURF1
 USE CONST
 USE SOIL
 IMPLICIT NONE
 INTERFACE
  REAL FUNCTION QSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION QSAT
  REAL FUNCTION DQSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION DQSAT
 END INTERFACE
 INTEGER, INTENT(IN)                  :: NDIM
 REAL, INTENT(IN)                     :: DT
 REAL, INTENT(IN),  DIMENSION(NDIM)   :: TS, T2, RA, RS, PS, TA, QA, &
&                                        RG, RL, WG, DELTA
 REAL, INTENT(INOUT), DIMENSION(NDIM) :: RSOIL
 REAL, INTENT(OUT), DIMENSION(NDIM)   :: TSN, T2N, HU
 REAL, DIMENSION (NDIM) :: ZBETA1, ZBETA2, ZRHO, ZQS, ZDQSDT, ZA, ZB, ZC, ZTAU2, ZHU
 INTEGER :: I
 REAL :: ZPI
 ZPI = 2*ASIN(1.)
 ZTAU2 = TAU/(2.0*ZPI)
!
 DO I = 1,NDIM
  ZDQSDT(I) = DQSAT (PS(I),TS(I))
  ZQS(I) = QSAT (PS(I),TS(I))
 ENDDO 
 CT = 1./(VEG/CV + (1. - VEG)/CG)
 WHERE (WG < WFC) 
   ZHU = 0.5*(1.0 - COS(WG/WFC*ZPI))
 ELSEWHERE
   ZHU = 1.0
 END WHERE   
 WHERE (ZQS < QA) ! Dew flux at potential rate
  ZHU = 1.0
  RSOIL = 0.0
 END WHERE
 WHERE (ZHU*ZQS < QA .AND. ZQS > QA) ! Avoid dew flux on very dry and warm soils
  ZHU = QA/ZQS
 END WHERE 
 HU = ZHU
 ZBETA1 = (1. - DELTA)*VEG/(RS + RA) + DELTA*VEG/RA + (1. - VEG)*ZHU/(RSOIL + RA)
 ZBETA2 = (1. - DELTA)*VEG/(RS + RA) + DELTA*VEG/RA + (1. - VEG)/(RSOIL + RA)
 ZRHO = PS / (RD * TA*(1. + 0.608*QA)) 
!
! Solve implicitely the surface energy balance (for ts)
!
 ZA  = -CT*(4.*EMIS*STEFAN*TS**3 + ZRHO*(CP/RA + LV*ZBETA1*ZDQSDT))
 ZB  =  CT*(3.*EMIS*STEFAN*TS**3 + ZRHO*LV*ZBETA1*ZDQSDT)
 ZC  =  CT*((1.-ALPHA)*RG + EMIS*RL + ZRHO*(CP*(TA + GRAV*ZREF/CP)/RA + LV*(ZBETA2*QA - ZBETA1*ZQS)))
 TSN = ((1. + ZB*DT)*TS + DT*T2/ZTAU2 + ZC*DT)/(1. - ZA*DT + DT/ZTAU2) 
! TSN = -(ZB*TS + ZC)/ZA ! solve without ground heat flux
!
! Evolution of mean surface temperature (t2)
!
 T2N = (T2 + TSN*DT/TAU)/(1. + DT/TAU)
 RETURN
END SUBROUTINE ENERGY_BUDGET

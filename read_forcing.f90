SUBROUTINE READ_FORCING
!----------------------------------------------------------------------
!
! Read atmospheric forcing for the whole period of interest
!
!                                     Jean-Francois MAHFOUF (11/06)
!----------------------------------------------------------------------
 USE FORC
 IMPLICIT NONE
 INTEGER :: I, IK
 PRINT *,'LECTURE FORCAGE '
 !OPEN(UNIT=10,FILE='../data_in/CORALYNN_05062019_smooth_50m.dat',FORM='FORMATTED')
 OPEN(UNIT=10,FILE='../data_in/CORA_LYNN_OBS_FORCING_05062019.dat',FORM='FORMATTED')
 DO I = 1,NFORC
  READ (10,*) IK,UAF(I),VAF(I),TAF(I),QAF(I),PSF(I),RGF(I),RLF(I),PRF(I)
 ENDDO  
 PRINT *,'FIN LECTURE FORCAGE' 
 CLOSE (UNIT=10)
 RETURN
END SUBROUTINE READ_FORCING

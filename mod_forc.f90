MODULE FORC
INTEGER, PARAMETER     :: NFORC=48*36+1
REAL, PARAMETER        :: DT_FORC=1800.
REAL, DIMENSION(NFORC) :: RGF, RLF, PRF, TAF, UAF, VAF, PSF, QAF
INTEGER                :: NSTEP_TOT, NSTEP_DEB
END MODULE FORC

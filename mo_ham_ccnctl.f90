MODULE mo_ham_ccnctl

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  !--- Prescribed supersaturations from CCN database provided by Dominick Spracklen (pers. comm., 2009)
  !    Included all supersaturations [%] with more than one measurement (31). 

  PRIVATE

  INTEGER,  PUBLIC, PARAMETER :: nsat =31
  REAL(dp), PUBLIC, PARAMETER :: zsat(nsat) = (/                                         &
        0.0002_dp,0.0004_dp,0.0006_dp,0.0007_dp,0.0008_dp,0.0010_dp,0.0016_dp,0.0020_dp, &
        0.0023_dp,0.0025_dp,0.0027_dp,0.0030_dp,0.0032_dp,0.0033_dp,0.0035_dp,0.0038_dp, &
        0.0040_dp,0.0050_dp,0.0055_dp,0.0060_dp,0.0065_dp,0.0066_dp,0.0070_dp,0.0075_dp, &
        0.0080_dp,0.0084_dp,0.0085_dp,0.0100_dp,0.0112_dp,0.0120_dp,0.0150_dp            /)

END MODULE mo_ham_ccnctl

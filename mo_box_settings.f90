MODULE mo_box_settings

  IMPLICIT NONE

  INTEGER, PARAMETER :: kproma = 1	! geographic block number of locations -> kbdim
  INTEGER, PARAMETER :: kbdim = 1	! maximum number of locations in block 
  INTEGER, PARAMETER :: klev = 1	! number of levels
  INTEGER, PARAMETER :: krow = 1	! geographic block number
  INTEGER, PARAMETER :: ktdia = 1	! highest vertical level for "diagnostics" -> klev 
  INTEGER, PARAMETER :: ktrac = 1	! numer of tracers -> ntrac


END MODULE mo_box_settings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! [mo_ham_m7.f90
!!
!! \brief
!! Module to provide interface to the M7 aerosol microphysics scheme
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# The original code is from J. Feichter, J. Wilson and E. Vignatti, JRC Ispra 
!! and was adapted for ECHAM by P. Stier, Oxford. Other contributions include 
!! D. O'Donnell, K. Zhang and others
!!   -# M.G. Schultz (FZ Juelich) - new code structure for integration into echam6-hammoz (2009-09-24)
!!
!! \limitations
!! None
!!
!! \details
!! This module contains the m7_interface routine and all individual routines
!! which make up M7. Parameter lists and flags are defined in mo_ham_m7ctl.
!! This module contains the following subroutines which used to be individual files.
!!       m7_interface
!!       m7_cumulative_normal      (renamed from m7_cumnor)
!!       m7
!!       m7_averageproperties
!!       m7_kappa
!!       m7_equiz
!!       m7_equimix
!!       m7_equil
!!       m7_h2so4_cs
!!       m7_prod_cond
!!       m7_nuck
!!       m7_dnum
!!       m7_dconc
!!       m7_dconc_soa
!!       m7_coaset
!!       m7_concoag
!!       m7_delcoa 
!!       m7_mass_sum
!!
!! \bibliographic_references
!!   - Vignati E., Wilson J. and Stier P., M7: a size resolved aerosol mixture module 
!!     for the use in global aerosol models, JGR 109, D22 202, doi:10.1029/2003JD004 485, 2004.
!!   - Stier P. et al, The aerosol-climate model ECHAM5-HAM, ACP, 5, 1125–1156, 2005
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_box_ham_m7

IMPLICIT NONE

PRIVATE

PUBLIC      :: m7_cumulative_normal


CONTAINS


! ---------------------------------------------------------------------------

SUBROUTINE m7_cumulative_normal ( arg, presult, ccum )
  !
  !*******************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is double precision
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is double precision
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is double precision
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  USE mo_kind, ONLY: dp
  !
  IMPLICIT NONE
  !
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287e00_dp, &
       1.6102823106855587881e02_dp, &
       1.0676894854603709582e03_dp, &
       1.8154981253343561249e04_dp, &
       6.5682337918207449113e-2_dp /)
  REAL(dp) :: arg
  REAL(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870e01_dp, &
       9.7609855173777669322e02_dp, &
       1.0260932208618978205e04_dp, &
       4.5507789335026729956e04_dp /)
  REAL(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764e-1_dp, &
       8.8831497943883759412e00_dp, &
       9.3506656132177855979e01_dp, &
       5.9727027639480026226e02_dp, &
       2.4945375852903726711e03_dp, &
       6.8481904505362823326e03_dp, &
       1.1602651437647350124e04_dp, &
       9.8427148383839780218e03_dp, &
       1.0765576773720192317e-8_dp /)
  REAL(dp) :: ccum
  REAL(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691e01_dp, &
       2.3538790178262499861e02_dp, &
       1.5193775994075548050e03_dp, &
       6.4855582982667607550e03_dp, &
       1.8615571640885098091e04_dp, &
       3.4900952721145977266e04_dp, &
       3.8912003286093271411e04_dp, &
       1.9685429676859990727e04_dp /)
  REAL(dp) :: del
!@@@ REAL(dp) :: dpmpar
  REAL(dp) :: eps
  INTEGER :: i
  REAL(dp) :: zmin
  REAL(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699e-1_dp, &
       1.274011611602473639e-1_dp, &
       2.2235277870649807e-2_dp, &
       1.421619193227893466e-3_dp, &
       2.9112874951168792e-5_dp, &
       2.307344176494017303e-2_dp /)
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121e00_dp, &
       4.68238212480865118e-1_dp, &
       6.59881378689285515e-2_dp, &
       3.78239633202758244e-3_dp, &
       7.29751555083966205e-5_dp /)
  REAL(dp) :: presult
  REAL(dp), PARAMETER :: root32 = 5.656854248_dp
  REAL(dp), PARAMETER :: sixten = 16.0_dp
  REAL(dp) :: temp
  REAL(dp), PARAMETER :: sqrpi = 3.9894228040143267794e-1_dp
  REAL(dp), PARAMETER :: thrsh = 0.66291_dp
  REAL(dp) :: x
  REAL(dp) :: xden
  REAL(dp) :: xnum
  REAL(dp) :: y
  REAL(dp) :: xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0_dp ) * 0.5_dp
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  zmin = EPSILON ( 1.0_dp )

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
    !
    !  Evaluate  anorm  for  |X| <= 0.66291
    !
    IF ( y > eps ) THEN
       xsq = x * x
    ELSE
       xsq = 0.0_dp
    END IF

    xnum = a(5) * xsq
    xden = xsq
    DO i = 1, 3
       xnum = ( xnum + a(i) ) * xsq
       xden = ( xden + b(i) ) * xsq
    END DO
    presult = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = presult
    presult = 0.5_dp + temp
    ccum = 0.5_dp - temp
    !
    !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
    !
  ELSEIF ( y <= root32 ) THEN

    xnum = c(9) * y
    xden = y
!CDIR UNROLL=7
    DO i = 1, 7
       xnum = ( xnum + c(i) ) * y
       xden = ( xden + d(i) ) * y
    END DO
    presult = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = AINT ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    presult = EXP(-xsq*xsq*0.5_dp) * EXP(-del*0.5_dp) * presult
    ccum = 1.0_dp - presult

    IF ( x > 0.0_dp ) THEN
       temp = presult
       presult = ccum
       ccum = temp
    END IF
    !
    !  Evaluate  anorm  for |X| > sqrt(32).
    !
  ELSE

    presult = 0.0_dp
    xsq = 1.0_dp / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    DO i = 1, 4
       xnum = ( xnum + p(i) ) * xsq
       xden = ( xden + q(i) ) * xsq
    END DO

    presult = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    presult = ( sqrpi - presult ) / y
    xsq = AINT ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    presult = EXP ( - xsq * xsq * 0.5_dp ) * EXP ( - del * 0.5_dp ) * presult
    ccum = 1.0_dp - presult  

    IF ( x > 0.0_dp ) THEN
       temp = presult
       presult = ccum
       ccum = temp
    END IF

  END IF

  IF ( presult < zmin ) THEN
    presult = 0.0_dp
  END IF

  IF ( ccum < zmin ) THEN
    ccum = 0.0_dp
  END IF

END SUBROUTINE m7_cumulative_normal


END MODULE mo_box_ham_m7

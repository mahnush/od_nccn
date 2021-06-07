MODULE mo_kind

  ! L. Kornblueh, MPI, August 2001, added working precision and comments 

  IMPLICIT NONE

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.

  ! Floating point section: 

  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: pi4 = 9
  INTEGER, PARAMETER :: pi8 = 14

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)

  ! Floating point working precision

  INTEGER, PARAMETER :: wp = dp   

  ! Integer section

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)

  ! Working precision for index variables
  !
  ! predefined preprocessor macros:
  !
  ! xlf         __64BIT__   checked with P6 and AIX
  ! gfortran    __LP64__    checked with Darwin and Linux
  ! Intel, PGI  __x86_64__  checked with Linux
  ! Sun         __x86_64    checked with Linux 

!K #if defined (__64BIT__) || defined (__LP64__) || defined (__x86_64__) || defined (__x86_64)
!K   INTEGER, PARAMETER :: widx = i8
!K #else
!K   INTEGER, PARAMETER :: widx = i4
!K #endif
!K 
!K CONTAINS
!K 
!K   SUBROUTINE print_kinds  
!K 
!K     PRINT *, 'single precision : ', sp
!K     PRINT *, 'double precision : ', dp
!K     PRINT *, 'working precision: ', wp
!K     PRINT *, '4 byte integer   : ', i4
!K     PRINT *, '8 byte integer   : ', i8
!K     PRINT *, 'index integer    : ', widx
!K 
!K #if defined (__64BIT__) 
!K     ! xlf
!K     PRINT *, '__64BIT__' 
!K #endif
!K #if defined (__LP64__)
!K     ! gfortran
!K     PRINT *, '__LP64__'
!K #endif
!K #if defined (__x86_64__)
!K     ! Intel, PGI
!K     PRINT *, '__x86_64__'
!K #endif
!K #if defined (__x86_64)
!K     ! Sun
!K     PRINT *, '__x86_64'
!K #endif
!K   ! NAG has no predefined macro but speed is no issue anyhow
!K 
!K   END SUBROUTINE print_kinds

END MODULE mo_kind

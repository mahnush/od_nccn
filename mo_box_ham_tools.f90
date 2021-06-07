!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_tools.f90
!!
!! \brief
!! mo_ham_tools hold auxiliary routines for the 
!! HAM aerosol model
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Philip Stier (MPI-Met) - original code (2002)
!!   -# Philip Stier (MPI-Met) - added ham_logtail (2004)
!!   -# Betty Croft (Dalhousie University) - added scavenging coefficient bilinear interpolation (2008)
!!   -# Sylvaine Ferrachat (ETH Zurich) - cleanup and security (2011)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
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

MODULE mo_box_ham_tools

  ! *mo_ham_tools* hold auxiliary routines for the 
  !                 HAM aerosol model

  USE mo_kind,               ONLY: dp
!K   USE mo_exception,          ONLY: finish

  IMPLICIT NONE

CONTAINS


!---------------------------------------------------------------------------------------------

  SUBROUTINE ham_logtail(kproma, kbdim,  klev,   krow,  kmod, &
                         ld_wetrad, ld_numb, pr,              &          
                         pfrac)

    ! *ham_logtail* calculates mass- or number-fraction larger than
    !               the radius pr for one given mode of a superposition 
    !               of nmod log-normal aerosol distributions
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !
    ! Revision:
    ! ---------
    ! Sylvaine Ferrachat, ETH Zurich,             2013
    !    --> this routine computes now one mode at a time, which allows to reduce the computational load
    !        for all irrelevant modes
    !
    ! Method:
    ! -------

    !
    ! The calculation of the activated number fraction and mass fraction
    ! from the radius of activation is done by a transformation of the 
    ! log-normal distribution to the error function which is then computed
    ! using the routine m7_cumulative_normal:
    !
    !                        / x                              _
    !                 N      |       1           1   ln(R)-ln(R)  2
    !    N(0,x) = ---------  |   --------  exp(- - ( ----------- )   ) d ln(R) 
    !             ln(sigma)  |   sqrt(2PI)       2    ln(sigma)
    !                        / 0 
    !                         
    !                         /tx                   2
    !                        |        1            t
    !           =     N      |     --------  exp(- - ) d t 
    !                        |     sqrt(2PI)       2 
    !                        /-inf
    ! 
    !    where:                   
    !
    !                        _
    !               ln(R)-ln(R)
    !    t      =   -----------
    !                ln(sigma)
    !
    !    and:
    !                        _
    !               ln(x)-ln(R)
    !    tx     =   -----------
    !                ln(sigma)


!K     USE mo_ham_m7ctl,    ONLY: nmod, sigmaln, cmedr2mmedr
    USE mo_box_aerosols, ONLY: aersigma, aercmedr2mmedr, aercmr
!K     USE mo_ham_m7,       ONLY: m7_cumulative_normal    
    USE mo_box_ham_m7,       ONLY: m7_cumulative_normal !K  
!K     USE mo_ham_streams,  ONLY: rwet, rdry 

    IMPLICIT NONE

    !--- Subroutine parameters:

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow, kmod

    LOGICAL, INTENT(in) :: ld_wetrad !wet vs dry radius switch
    LOGICAL, INTENT(in) :: ld_numb   !number vs mass switch

    REAL(dp), INTENT(in) :: pr(kbdim,klev) !lower bound radius

    REAL(dp), INTENT(out) :: pfrac(kbdim,klev)

    !--- Local variables:

    INTEGER :: jl, jk

    REAL(dp) :: zt, zdummy, zeps, zfact

!K     REAL(dp), POINTER :: r_p(:,:)

!K     REAL(dp) :: r_p(kbdim,klev)

    REAL(dp) :: r_p(kbdim,klev)		!K count median radius [m]
   
    !--- 0) 
    zeps=EPSILON(1.0_dp)

    !--- 1) 

!K     IF (ld_wetrad) THEN
!K        r_p => rwet(kmod)%ptr(1:kproma,:,krow)
!K     ELSE
!K        r_p => rdry(kmod)%ptr(1:kproma,:,krow)
!K     ENDIF

    IF (ld_wetrad) THEN
!K       r_p = rwet
      r_p(kbdim,klev) = aercmr(kmod)!K  we don't have any rwet!
    ELSE
!K        r_p = rdry
      r_p(kbdim,klev) = aercmr(kmod)*10.**(-6.) !K  count median radius: microns to m
    ENDIF




    !>>SF
    IF (ld_numb) THEN !number calculation
       zfact = 1._dp
    ELSE !mass calculation
!K        zfact = cmedr2mmedr(kmod)
       zfact = aercmedr2mmedr(kmod) !K is not used since ld_numb = True
    ENDIF
    !<<SF

    DO jk=1, klev
       DO jl=1, kproma

          IF (pr(jl,jk)>zeps .AND. r_p(jl,jk)>zeps) THEN

             !--- Transform number distribution to error function:

!K              zt=(LOG(pr(jl,jk))-LOG(r_p(jl,jk)*zfact))/sigmaln(kmod)
             zt=(LOG(pr(jl,jk))-LOG(r_p(jl,jk)*zfact))/aersigma(kmod)

             !--- Calculate the cumulative of the log-normal number distribution:

             CALL m7_cumulative_normal(zt,zdummy,pfrac(jl,jk))

             !--- Calculate the cumulative of the log-normal mass distribution:

          ELSEIF (pr(jl,jk)<zeps .AND. r_p(jl,jk)>zeps) THEN

             pfrac(jl,jk)=1.0_dp

          ELSE

             pfrac(jl,jk)=0.0_dp

          END IF

       END DO !kproma
    END DO !klev

  END SUBROUTINE ham_logtail


END MODULE mo_box_ham_tools

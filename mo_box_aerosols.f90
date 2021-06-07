MODULE mo_box_aerosols

! aerosol properties and size distributions


  USE mo_kind, 		ONLY: dp
  USE mo_box_settings,	ONLY: kproma, kbdim, klev, krow
  USE mo_constants,	ONLY: rd, api

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER, PARAMETER  :: nmod = 11
  INTEGER, PARAMETER  :: naerocomp = nmod
  INTEGER, PARAMETER  :: ntrac = nmod

  LOGICAL, PARAMETER  :: AEROSOLMASS = .TRUE.

!-------------------------------------------------------------------------------
  ! aerosol species (n = non-hygroscopic, h = hygroscopic)
  CHARACTER(LEN=3), PARAMETER :: aertype(nmod)=(/'SS1','SS2','SS3', &
  						 'DU1','DU2','DU3', &
						 'OMn','OMh','BCn','BCh','SO4'/)

  ! molecular weights [g mol-1]
  REAL(dp), PARAMETER :: mw_so4 = 96.0631_dp, &		! molecular weight SO4
                         mw_bc = 12.010_dp,   &		! molecular weight BC
			 mw_om = 180._dp,     &		! molecular weight OM (primarily OC)
			 mw_ss = 58.443_dp,   &		! molecular weight SS
			 mw_du = 250._dp		! molecular weight DU

  REAL(dp), PARAMETER :: aermw(nmod) = (/mw_ss, mw_ss, mw_ss, &
  					 mw_du, mw_du, mw_du, &
					 mw_om, mw_om, mw_bc, mw_bc, mw_so4/)


 !-------------------------------------------------------------------------------
 ! new size distributions and aerosol densities (adjusted by Karoline Block)
 
   ! densities [g cm-3]
  REAL(dp), PARAMETER :: dens_so4 = 1.84_dp, &        ! density SO4 (ECHAM-HAM - between IFS & Pinty)
		       dens_bc = 1.80_dp,    &	      ! density BC (ECHAM-HAM & Pinty)
		       dens_om = 1.76_dp,  &	      ! density OM (IFS & Pinty)
		       dens_ss = 2.16_dp,  &	      ! density SS (ECHAM-HAM & IFS & Pinty)
		       dens_du = 2.60_dp	      ! density DU (ECHAM-HAM & IFS)
 
  REAL(dp), PARAMETER :: aerdens(nmod) = (/dens_ss, dens_ss, dens_ss, &
  					   dens_du, dens_du, dens_du, &
					   dens_om, dens_om, dens_bc, dens_bc, dens_so4/)
 
  ! count median radius [microns]
  REAL(dp), PARAMETER :: cmr_so4 = 0.0355_dp, &       
		       cmr_bcn = 0.0118_dp,  &	     
		       cmr_bch = 0.0118_dp,  &
		       cmr_omn = 0.0355_dp,  &
		       cmr_omh = 0.0355_dp,  &	     
		       cmr_ss1 = 0.125_dp,   &
		       cmr_ss2 = 1.60_dp,	&
		       cmr_ss3 = 10.0_dp,	&
		       cmr_du1 = 0.135_dp,	&
		       cmr_du2 = 0.704_dp,	&	      
		       cmr_du3 = 4.4_dp 	              
 
  REAL(dp), PARAMETER :: aercmr(nmod) = (/cmr_ss1, cmr_ss2, cmr_ss3, &
  					  cmr_du1, cmr_du2, cmr_du3, &
					  cmr_omn, cmr_omh, cmr_bcn, cmr_bch, cmr_so4/)
 
  ! geometric standard deviation
   REAL(dp), PARAMETER :: gsd_so4 = 2.0_dp, &	      
			gsd_bc = 2.0_dp, &          
			gsd_om = 2.0_dp, &          
			gsd_ss1  = 2.0_dp, &
			gsd_ss2  = 2.0_dp, &
			gsd_ss3  = 2.0_dp, &
			gsd_du1  = 2.0_dp, &
			gsd_du2  = 2.0_dp, &
			gsd_du3  = 2.0_dp           

   REAL(dp), PARAMETER :: aergsd(nmod) = (/gsd_ss1, gsd_ss2, gsd_ss3, &
   					   gsd_du1, gsd_du2, gsd_du3, &
					   gsd_om, gsd_om, gsd_bc, gsd_bc, gsd_so4/)
					   
   REAL(dp), PARAMETER :: aersigma(nmod) = LOG(aergsd(:))


   
!-------------------------------------------------------------------------------


  LOGICAL, PARAMETER  :: aersol(nmod) = (/.TRUE.,.TRUE.,.TRUE., &
  					  .FALSE.,.FALSE.,.FALSE., &
					  .FALSE.,.TRUE.,.FALSE.,.TRUE.,.TRUE./)
  
  INTEGER, PARAMETER  :: aernion(nmod) = (/2,2,2,0,0,0,0,1,0,1,2/)
  
  REAL(dp), PARAMETER :: aerosm(nmod) = (/1.0_dp,1.0_dp,1.0_dp, &
  					  0.0_dp,0.0_dp,0.0_dp, &
					  0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp/)
  
  REAL(dp), PARAMETER :: aercmedr2mmedr(nmod) = EXP(3.0_dp*(aersigma(nmod)**2))

   
!-------------------------------------------------------------------------------

  !--- Subroutines:

CONTAINS

  SUBROUTINE sizedist(nParams,mmr,prho,mmr_tot,mc,nc,ncm)
  
  IMPLICIT NONE 

    ! Input parameters:  
    INTEGER,  INTENT(in)	:: nParams  
    REAL(dp), INTENT(in)	:: mmr(kbdim,klev,nParams)
    REAL(dp), INTENT(in)	:: prho(kbdim,klev)
    
    ! Output parameters: 
    REAL(dp), INTENT(out)  	:: mmr_tot(kbdim,klev,nmod), &
    				   mc(kbdim,klev,nmod), &
    				   nc(kbdim,klev,nmod), &
				   ncm(kbdim,klev,nmod)

!    ! Local parameters: 
!    REAL(dp), PARAMETER	:: rmin = 0.0001_dp, &		!min and max radius for the size integration
!  			   	   rmax = 50.0_dp	
!  
!    INTEGER, PARAMETER	:: imax = 1000    		!number of radii intervals
!
!    INTEGER			:: ir 
!    REAL(dp), DIMENSION(imax) 	:: r, dn, dm, n, m     
!    REAL(dp)			:: dlogr, ntot, mtot
!    REAL(dp)			:: cmedr2mmedr

    INTEGER			:: ip 
    REAL(dp)			:: gsd, cmr, dens, sigmaln, cmr2ram, massp


!-------------------------------------------------------------------------------

    IF(nmod .EQ. nParams) THEN             
       DO ip = 1,nParams

	 dens = aerdens(ip)	
	 cmr = aercmr(ip)
	 gsd = aergsd(ip)		

	 ! Using Hatch-Coatch Conversion factor	
	 sigmaln = LOG(gsd)
   !      cmedr2mmedr = EXP(3.0_dp*(sigmaln**2))  ! count median radius -> mass median radius 
	 cmr2ram = EXP(1.5_dp*(sigmaln**2))      ! count median radius -> radius of average mass

	 !Mass mixing ratio to Mass concentration
	 mmr_tot(kbdim,klev,ip) = mmr(kbdim,klev,ip)
	 mc(kbdim,klev,ip) = mmr_tot(kbdim,klev,ip)*prho(kbdim,klev)		!kg/kg -> kg/m3
         !WRITE(*,*) 'prho', prho(kbdim,klev)
	 !Mass concentration to Number concentration:
	 massp = 4./3.*api*dens*10.**(-3.)*((cmr*10.**(-4.))*cmr2ram)**3.		!radius: microns -> cm; density: g/cm3 -> kg/cm3
	 nc(kbdim,klev,ip) = mc(kbdim,klev,ip)/massp			!nc: /m3

	 !Number concentration to Number conc. per unit mass:
	 ncm(kbdim,klev,ip) = nc(kbdim,klev,ip)/prho(kbdim,klev)


   !      !initialisations
   !      ntot = 0.0
   !      mtot = 0.0
   !
   !      !size integration
   !      DO ir=1, imax 
   ! 	 dlogr = 1./FLOAT(imax)*(log(rmax)-log(rmin)) 
   ! 	 r(ir) = exp(log(rmin)+FLOAT(ir)*dlogr) 	  
   !
   ! 	 dm(ir)= (1./(sqrt(2.*api)*sigmaln))*exp((-0.5*(log(r(ir))-log(cmr*cmedr2mmedr))**2.)/sigmaln**2.)
   ! 	 dn(ir)= (1./(sqrt(2.*api)*sigmaln))*exp((-0.5*(log(r(ir))-log(cmr))**2.)/sigmaln**2.) 
   !
   ! 	 m(ir) = mc(kbdim,klev,ip)*dm(ir)*dlogr
   ! 	 n(ir) = nc(kbdim,klev,ip)*dn(ir)*dlogr	 
   !
   ! 	 ntot = ntot+n(ir) 
   ! 	 mtot = mtot+m(ir)
   ! 	 !WRITE(*,*) 'i', ir, 'r', r(ir), 'dn(r)', dn(ir), 'n(r)', n(ir), 'Ntot', ntot, 'dm(r)', dm(ir), 'm(r)', m(ir), 'Mtot', mtot
   !      ENDDO !imax	   
   !
   !      !WRITE(*,*) ' '			 
   !      !WRITE(*,*) 'ip', ip
   !      !WRITE(*,*) 'mixratio [kg/kg]', mmr(:,:,ip)
   !      !WRITE(*,*) 'mass concentration[kg/m3]', mc(:,:,ip)	  
   !      !WRITE(*,*) 'check mtot', mtot	
   !      !WRITE(*,*) 'number concentration[/m3]', nc(:,:,ip)
   !      !WRITE(*,*) 'check ntot', ntot	

      ENDDO  !ip 

    ELSE
      STOP "Stopped: mismatch nmod and nParams"
    ENDIF
    
    
  END SUBROUTINE sizedist





END MODULE mo_box_aerosols

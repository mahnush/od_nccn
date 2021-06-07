!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_activ.f90
!!
!! \brief
!! This code handles all ham-dependent activation processes
!! 
!!
!! \author Philip Stier (MPI-Met)
!! \author Sylvaine Ferrachat (ETHZ)
!! \author Daniel Partridge (AOPP,Oxford)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# P. Stier (MPI-Met)  - original code 
!!   -# S. Ferrachat (ETHZ) - new code structure (separate HAM-dependent pieces 
!!                            from the general activation schemes) - (2010-03)
!!   -# D.G Partridge (AOPP-Oxford)  - Added Nenes activation scheme - (2014-02)
!!
!! \limitations
!! None
!!
!! \details
!! Implementation of :
!!   - the Abdul-Razzak & Ghan scheme which calculates the number of activated aerosol 
!!     particles from the aerosol size-distribution, composition and ambient supersaturation 
!!     (see subroutine ham_activ_abdulrazzak_ghan)
!!   - a preparatory routine for Lin & Leaitch activation scheme (HAM-specific), which computes
!!     the fractional mass and number of each mode larger than the cutoff of the instrument and add them up
!!
!! \bibliographic_references
!!    - Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
!!    - Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
!!    - Pruppbacher and Klett, Kluewer Ac. Pub., 1997.
!!    - Barahona, D. et al., ACP, 10, 2467-2473, 2010.
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

MODULE mo_ham_activ

  USE mo_kind,          ONLY: dp
!K   USE mo_ham_m7ctl,     ONLY: nmod
!K   USE mo_ham,	ONLY: m7mode
  USE mo_box_updraft,       ONLY: nw, pw 


  IMPLICIT NONE

  PUBLIC ham_activ_abdulrazzak_ghan !K , ham_activ_diag_abdulrazzak_ghan_strat, ham_activ_diag_abdulrazzak_ghan_conv
!K   PUBLIC ham_avail_activ_lin_leaitch, ham_activ_diag_lin_leaitch
  PUBLIC ham_activ_koehler_ab

  PRIVATE

  !--- Subroutines:

CONTAINS

!K     SUBROUTINE ham_activ_abdulrazzak_ghan(kproma,   kbdim,   klev,  krow,  ktdia, &
!K   				       pcdncact, pesw,	prho,		     &
!K   				       pxtm1,    ptm1,	papm1, pqm1,	     &
!K   				       pw,       pwpdf,	pa,    pb,	     &
!K   				       pnact,    prc,	psmax		     )

    SUBROUTINE ham_activ_abdulrazzak_ghan(kproma,    kbdim,    klev,  krow,  ktdia, 	&
    					  ncd_activ, pcdncact, pesw,  			&
				    	  zn,        ptm1,     papm1, pqm1,	  	&
				    	  pa,        pb, 	                 	& 
				    	  pnact,     prc,      psmax		  	)



      ! *ham_activ_abdulrazzak_ghan* calculates the number of activated aerosol 
      !		particles from the aerosol size-distribution,
      !		composition and ambient supersaturation
      !
      ! Author:
      ! -------
      ! Philip Stier, MPI-MET, Caltech, University of Oxford  2002-2009
      !
      ! Method:
      ! -------
      ! The calculation of the activation can be reduced to 4 tasks:
      ! 
      ! 0)   Calculation of Koehler A/B coefficients 
      !	(now done in ham_activ_koehler_ab)
      ! I)   Calculate the maximum supersaturation
      ! II)  Calculate the corresponding radius of activation
      !	for each mode
      ! III) Calculate the number of particles that are larger
      !	then the radius of activation for each mode.
      ! 
      ! III) Calculation of the number of activated particles:
      !	See the routine ham_logtail below.
      !
      ! The calculations are now performed separately for 
      ! stratiform and convective updraft velocities.
      !
      ! References:
      ! -----------
      ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
      ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
      ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

      !>>dod soa
!K       USE mo_ham_m7ctl, ONLY: nmod, sigmaln
      !<<dod
!K       USE mo_ham_tools, ONLY: ham_logtail
      USE mo_box_ham_tools, ONLY: ham_logtail
!K       USE mo_tracdef,	ONLY: trlist, ntrac, AEROSOLMASS
!K       USE mo_ham,	ONLY: naerocomp, aerocomp, m7mode
!K       USE mo_ham_streams, ONLY: rdry, &
!K   			 a, b, sc, rc_strat, nact_strat, nact_conv !SF
!K       USE mo_constants, ONLY: rhoh2o, argas, api, cpd, g, alv, amw, amd, tmelt
      USE mo_constants, ONLY: rhoh2o, argas, api, cpd, g, alv, amw, amd, tmelt, zsten
!K       USE mo_cloud,	ONLY: cthomi ! minimum T for mixed clouds
   !>>dod soa
!K       USE mo_ham_species, ONLY: id_so2, id_so4, id_so4g, id_bc, id_oc, id_du
   !<<dod
   !>>SF
!K       USE mo_activ,	ONLY: lwetrad, nw
!K       USE mo_conv,	ONLY: cdncact_cv
      !<<SF
    !>>DP
!K     USE mo_param_switches,    ONLY: ncd_activ
    !<<DP
      !SF note: lwetrad is false in this routine, ie dry radius is used

!---------------------------------------------------------------------------------------------
    !--added for Nenes scheme
      !USE mo_ham_nenes,	ONLY: ccnspec,      &
       !                       pdfactiv,     &
!K                              props,        &   
      !                        MODE_NENES => MODE
!End of Nenes parameter additions
!---------------------------------------------------------------------------------------------

      USE mo_box_aerosols,	ONLY: nmod, ntrac, naerocomp, AEROSOLMASS, &	!K 
      				      aersigma, aercmr 				!K 
  !    USE mo_box_updraft,	ONLY: nw, pw          				!K 



      IMPLICIT NONE

      !--- Arguments:
      !--added for Nenes scheme
      !INTEGER, PARAMETER   :: nsupersat = 2     ! Switch for the calulation of the maximum
                                                 ! supersaturation:
                                                 !            2 Abdul-Razzak and Ghan (2000)
                                                 !            3 Barahona et al. (2010)

      INTEGER, INTENT(in)   :: kproma, kbdim, klev, krow, ktdia, ncd_activ

      REAL(dp), INTENT(out) :: pcdncact(kbdim,klev,nw),	& ! number of activated particles
  			       pnact(kbdim,klev),	& ! number of activated particles per mode [m-3]
  			       prc(kbdim,klev,nmod,nw),	& ! critical radius of activation per mode [m]
  			       psmax(kbdim,klev,nw)  	  ! maximum supersaturation [% 0-1]

      REAL(dp), INTENT(in)  :: ptm1(kbdim,klev),	& ! temperature
  			       papm1(kbdim,klev),	& ! pressure 
!K   			       prho(kbdim,klev),	& ! air density
  			       pqm1(kbdim,klev),	& ! specific humidity
  			       pesw(kbdim,klev),	& ! saturation water vapour pressure
  			       !pw(kbdim,klev,nw),	& ! mean or bins of updraft velocity (>0.0) [m s-1]
!msh
                               !pw(nw),                  & !updraft velocity
  			       !pwpdf(kbdim,klev,nw),	& ! PDF of updraft velocity [s m-1]
!K   			       pxtm1(kbdim,klev,ntrac),	& 
			       zn(kbdim,klev,nmod),	& !K number concentration [/m3] instead of number per unit mass
  			       pa(kbdim,klev,nmod),	& ! curvature parameter A of the Koehler equation
  			       pb(kbdim,klev,nmod)	  ! hygroscopicity parameter B of the Koehler equation
 
      !--- Local variables:
 
!K       INTEGER :: jclass,      jk,	     &
!K 		 jt,	      jl,	     &
!K 		 jw

      INTEGER :: jclass,      jk,	     &
		 jl,	      jw, 	jt


      REAL(dp):: zkv,                        &
                 zeps,      zalpha,          &
                 zgamma,                     &
                 zeta,                       &
                 zxi,                        &
                 zgrowth,                    &
                 zdif,      zxv,             &
                 zk,        zka

      REAL(dp):: zamw,			& ! molecular weight of water [kg mol-1]
  	         zamd			  ! molecular weight of dry air [kg mol-1]

      REAL(dp):: zf(nmod),  zg(nmod)

      REAL(dp):: zcdncact_top(kbdim,klev,nmod), & ! w * dw weighted activated aerosol number
  	         zcdncact_bot(kbdim,klev,nmod)    ! weighting factor  (sum of w * dw)

!K       REAL(dp):: zn(kbdim,klev,nmod),	     & ! aerosol number concentration for each mode [m-3]
!K 	      zsm(kbdim,klev,nmod),  & ! critical supersaturation for activating particles
!K 						 !    with the mode number median radius
!K 	      zrdry(kbdim,klev,nmod)   ! dry radius for each mode 
						 !@@@ currently re-stored, check if avoidable!

      REAL(dp):: zsm(kbdim,klev,nmod),       & ! critical supersaturation for activating particles
					       !    with the mode number median radius
	         zrdry(kbdim,klev,nmod)   	! dry radius for each mode 


      REAL(dp), ALLOCATABLE :: zfracn(:,:,:,:) ! fraction of activated aerosol numbers for each mode and w bin
      REAL(dp), ALLOCATABLE :: zsum(:)

!K       REAL(dp), PARAMETER :: zsten = 75.0E-3_dp ! surface tension of H2O [J m-2] 
!K       					     !   neglecting salts and temperature
!K       					      !	(also tried P&K 5.12 - erroneous!)
      LOGICAL, PARAMETER :: lwetrad = .FALSE.   ! switch between wet/dry radii in logtail calculation 
      LOGICAL, PARAMETER :: ll_numb = .TRUE.  ! switch between number/mass in logtail calculation !SF
      REAL(dp)           :: cthomi  = tmelt-35.0_dp !K instead of using cthomi from mo_cloud
!K       REAL(dp)           :: pnactnew(kbdim,klev,nmod,nw) !K


!---------------------------------------------------------------------------------------------
!--added for Nenes scheme
!msh
     
  !--- Switches:
!Deciding whether to call PDFACTIV (if modes >0)
 !     INTEGER :: modes_active

!CHECK: REMOVED (dp)
!      REAL, PARAMETER :: ACCOM  = 1            ! Accommodation coefficient (common for all CCN)

!     DATA FOR FHH Parameters
!      REAL, PARAMETER :: A_Nenes = 2.25              ! Afhh parameter (predefined from dust CCN activity experiments)
!      REAL, PARAMETER :: B_Nenes = 1.20              ! Afhh parameter (predefined from dust CCN activity experiments)
!      REAL:: AKKI(nmod)	                             ! Aerosol hygroscopicity
!      REAL:: ei(nmod)	                             ! Aerosol soluble mass fraction
!      REAL:: sigw				     ! Run Nenes module with PDF turned (1) on (0) off	
!      REAL:: nact(kbdim,klev)                        ! Number of activated particles
!      REAL:: zsm_nenes(kbdim,klev,nmod)              ! critical supersaturation for activating particles
                                                     ! with the mode number median radius
!      REAL:: psmax_nenes(kbdim,klev,nw)              ! maximum supersaturation [% 0-1]

!      REAL(dp) :: zsigma(nmod)                       ! Geometric std dev


! ** Mode 1 Information
!      MODE_NENES (1) = 1          ! Kohler Mode
!      ei   (1) = 0.         !insoluble fraction
!      AKKI (1) = 1.1        !hygroscopicity (k parameter)
!
! ** Mode 2 Information
!      MODE_NENES (2) = 1          ! Kohler Mode
!      ei   (2) = 0.         !insoluble fraction
!      AKKI (2) = 1.1        !hygroscopicity (k parameter)
!
! ** Mode 3 Information
!      MODE_NENES (3) = 1          ! FHH Mode
!      ei   (3) = 0.9        !insoluble fraction
!      AKKI (3) = 0.6        !hygroscopicity (k parameter)
!
! ** Mode 4 Information
!      MODE_NENES (4) = 1          ! FHH Mode
!      ei   (4) = 0.9        !insoluble fraction
!      AKKI (4) = 0.6        !hygroscopicity (k parameter)
!
! ** Mode 5 Information
!      MODE_NENES (5) = 1          ! FHH Mode
!      ei   (5) = 0.9        !insoluble fraction
!      AKKI (5) = 0.6        !hygroscopicity (k parameter)
!
! ** Mode 6 Information
!      MODE_NENES (6) = 1          ! FHH Mode
!      ei   (6) = 0.9        !insoluble fraction
!      AKKI (6) = 0.6        !hygroscopicity (k parameter)
!
! ** Mode 7 Information
!      MODE_NENES (7) = 1          ! FHH Mode
!      ei   (7) = 0.9        !insoluble fraction
!      AKKI (7) = 0.6        !hygroscopicity (k parameter)

!End of Nenes parameter additions
!---------------------------------------------------------------------------------------------



      !--- 0) Initializations:

      ALLOCATE(zfracn(kbdim,klev,nmod,nw))
      ALLOCATE(zsum(nw))

      zsm(1:kproma,:,:)	    = 0._dp
      zcdncact_top(1:kproma,:,:) = 0._dp
      zcdncact_bot(1:kproma,:,:) = 0._dp
      prc(1:kproma,:,:,:)        = 1._dp ! [m] initialized with 1m, only changed if activation occurs
      pnact(1:kproma,:)        = 0._dp
      pcdncact(1:kproma,:,:)       = 0._dp
      zeps=EPSILON(1._dp)

!---------------------------------------------------------------------------------------------
!--added for Nenes scheme
!msh
!      nact(1:kproma,:)        = 0.0
!      zsm_nenes(1:kproma,:,:) = 0.0

!--- convert ln(sigmag)
!      DO jclass=1, nmod
!K         zsigma(jclass) =   EXP(sigmaln(jclass))
!	  zsigma(jclass) =   EXP(aersigma(jclass))	!K
!      END DO
!End of Nenes parameter additions
!---------------------------------------------------------------------------------------------

      DO jclass=1, nmod
!K          zrdry(1:kproma,:,jclass)=rdry(jclass)%ptr(1:kproma,:,krow)
        zrdry(1:kproma,:,jclass) = aercmr(jclass)*10.**(-6.)	!K 
      END DO


      !--- Conversions to SI units [g mol-1 to kg mol-1]:

      zamw=amw*1.E-3_dp
      zamd=amd*1.E-3_dp

!K       !--- Number per unit volume for each mode:
!K       DO jclass=1, nmod
!K          jt = m7mode(jclass)%idt_no
!K          zn(1:kproma,:,jclass)=pxtm1(1:kproma,:,jt)*prho(1:kproma,:) 
!K       END DO

!    SELECT CASE (nsupersat)
!    CASE(2)  

  IF (ncd_activ==2) THEN

      ! (7):
!K       zf(:)=0.5_dp*EXP(2.5_dp*sigmaln(:)**2._dp)
      zf(:)=0.5_dp*EXP(2.5_dp*aersigma(:)**2._dp)

      ! (8):
!K       zg(:)=1._dp+0.25_dp*sigmaln(:)
      zg(:)=1._dp+0.25_dp*aersigma(:)

      !--- 1) Calculation of Koehler A/B coefficients: 
      !	 Now done in ham_activ_koehler_ab once so that they can be used in
      !	 convective and stratiform activation 

      !--- 2) Calculate maximum supersaturation:

      !--- 2.1) Abbdul-Razzak and Ghan (2000):
      !	   (Equations numbers from this paper unless otherwise quoted)

      DO jk=ktdia, klev
  	DO jl=1, kproma

      	    !--- Water vapour pressure:

    !!$	       zew=pqm1(jl,jk)*prho(jl,jk)*rv*ptm1(jl,jk)
!msh
  	    IF( (nw>1 .OR. pw(1)>zeps) .AND. &
  		pqm1(jl,jk)>zeps  	   .AND. &
  		ptm1(jl,jk)>cthomi		 ) THEN  

      	       !--- Abdul-Razzak et al. (1998) (Eq. 11):

  	       zalpha=(g*zamw*alv)/(cpd*argas*ptm1(jl,jk)**2) - &
  		      (g*zamd)/(argas*ptm1(jl,jk))

  	       zgamma=(argas*ptm1(jl,jk))/(pesw(jl,jk)*zamw) +  &
  		      (zamw*alv**2)/(cpd*papm1(jl,jk)*zamd*ptm1(jl,jk))

      	       !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:

  	       zdif=0.211_dp * (ptm1(jl,jk)/tmelt)**1.94_dp * (101325._dp/papm1(jl,jk)) *1.E-4_dp

      	       !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:

      	       ! Mole fraction of water:

  	       zxv=pqm1(jl,jk)*(zamd/zamw)

  	       zka=(5.69_dp+0.017_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp

  	       zkv=(3.78_dp+0.020_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp
                                       
      	       ! Moist air, convert to [J m-1 s-1 K-1]:

  	       zk =zka*(1._dp-(1.17_dp-1.02_dp*zkv/zka)*zxv) * 4.1868_dp*1.E2_dp

      	       !--- Abdul-Razzak et al. (1998) (Eq. 16):

  	       zgrowth=1._dp/						      &
  			 ( (rhoh2o*argas*ptm1(jl,jk))/(pesw(jl,jk)*zdif*zamw) + &
  			 (alv*rhoh2o)/(zk*ptm1(jl,jk)) * ((alv*zamw)/(ptm1(jl,jk)*argas) -1._dp) )
               !WRITE(*,*) 'zgrowth', zgrowth 
      	       !--- Summation for equation (6):

  	       zsum(:)=0._dp !K critical Ssat

    !CDIR UNROLL=7
               DO jclass=1, nmod
                 ! WRITE(*,*) 'zn', zn(jl,jk,jclass)
  		  IF (zn(jl,jk,jclass)	> zeps     .AND. &
  		      zrdry(jl,jk,jclass) > 1.E-9_dp .AND. &
  		      pa(jl,jk,jclass)	> zeps     .AND. &
  		      pb(jl,jk,jclass)	> zeps  	 ) THEN
                      !WRITE(*,*) 'zn,zdry,pa,pb > zeps'
      		     ! (9): !K crit. Ssat for activating particles with r = cmr

  		     zsm(jl,jk,jclass)=2._dp/SQRT(pb(jl,jk,jclass)) * &
  				     (pa(jl,jk,jclass)/(3._dp*zrdry(jl,jk,jclass)))**1.5_dp
		     
  		     
		     DO jw=1, nw
      			! (10):
!msh
  			zxi=2._dp*pa(jl,jk,jclass)/3._dp * SQRT(zalpha*pw(jw)/zgrowth) 
			

      			! (11):

  			zeta=((zalpha*pw(jw)/zgrowth)**1.5_dp) / &
  				(2._dp*api*rhoh2o*zgamma*zn(jl,jk,jclass))
			
      			! (6): !K crit. Ssat depending on updraft

  			IF (pw(jw)>zeps) THEN
  			   zsum(jw)=zsum(jw) + ( 1._dp/zsm(jl,jk,jclass)**2		   &
  						 * ( zf(jclass)*(zxi/zeta)**1.5_dp	   &
  						     + zg(jclass)*( zsm(jl,jk,jclass)**2._dp &
  								    / (zeta+3._dp*zxi) )**0.75_dp ) )	
  			
			END IF 

  		     END DO ! jw

  		  ENDIF

  		  WHERE (zsum(:) > zeps)
  		     psmax(jl,jk,:)=1._dp/SQRT(zsum(:))	!K max Ssat		    
  		  ELSEWHERE
  		     psmax(jl,jk,:)=0._dp
  		  END WHERE
  	       END DO ! jclass

  	    ELSE
  	       psmax(jl,jk,:)=0._dp
  	    END IF

  	END DO ! jl
      END DO ! jk

!--------------------------------------------------------------------------------
!msh      
!--added call statements for Nenes scheme (which can be found in mo_ham_nenes.f90)   
   ELSEIF (ncd_activ>2) THEN
   write(*,*) 'ncd_active>2'
!   DO jw=1, nw
!       DO jk=ktdia, klev
!          DO jl=1, kproma

!          modes_active = 0

!             IF( (nw>1 .OR. pw(jl,jk,1)>zeps) .AND. &
!                 pqm1(jl,jk)>zeps             .AND. &
!                 ptm1(jl,jk)>cthomi                 ) THEN  

                !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:
                ! Mole fraction of water:
!                zxv=pqm1(jl,jk)*(zamd/zamw)
!                zka=(5.69_dp+0.017_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp
!                zkv=(3.78_dp+0.020_dp*(ptm1(jl,jk)-273.15_dp))*1.E-5_dp
!                ! Moist air, convert to [J m-1 s-1 K-1]:
!                zk =zka*(1._dp-(1.17_dp-1.02_dp*zkv/zka)*zxv) * 4.1868_dp*1.E2_dp
               
!                DO jclass=1, nmod
!                   IF (zn(jl,jk,jclass)    > zeps     .AND. &
!                       zrdry(jl,jk,jclass) > 1.E-9_dp .AND. &
!                       pa(jl,jk,jclass)    > zeps     .AND. &
!                       pb(jl,jk,jclass)    > zeps           ) THEN

!	               CALL CCNSPEC (REAL(zn(jl,jk,jclass)), REAL(zrdry(jl,jk,jclass)), REAL(zsigma(jclass)),	&
!	      		             REAL(ptm1(jl,jk)), REAL(papm1(jl,jk)), nmod, AKKI(jclass), ei(jclass),	&
!			             A_Nenes, B_Nenes, zsm_nenes(jl,jk,jclass), REAL(pa(jl,jk,jclass)),		&
!			             REAL(pb(jl,jk,jclass)), jclass, REAL(pesw(jl,jk)), REAL(zk)		) !Calc zsm

!                       zsm(jl,jk,jclass) = zsm_nenes(jl,jk,jclass)
!                       modes_active = modes_active +1

!                   ENDIF
!                END DO ! jclass

!                sigw = 0._dp	! 0 = calculate for single updraft (no PDF)
!        	IF (modes_active > 0) THEN

   !                   DO jw=1, nw

!                	 CALL PDFACTIV (REAL(pw(jl,jk,jw)), REAL(zn(jl,jk,:)), AKKI, ei,		&
!		   			A_Nenes, B_Nenes, ACCOM, zsm_nenes(jl,jk,:),			&
!					sigw, REAL(ptm1(jl,jk)), REAL(papm1(jl,jk)),			&
!					nact(jl,jk), psmax_nenes(jl,jk,jw), nmod,REAL(zsigma(:))	) ! Calculate Ndrop

!                	 psmax(jl,jk,jw) = psmax_nenes(jl,jk,jw)

   !                   END DO ! jw

!        	ELSE
!                      psmax(jl,jk,:)=0._dp
!        	END IF

!             ELSE
!                 psmax(jl,jk,:)=0._dp
!             ENDIF             

!          END DO ! jl
!       END DO ! jk
!   END DO ! jw

  !   END SELECT
  END IF
!End of Nenes parameter additions
!--------------------------------------------------------------------------------


!K 	 !--- Diagnostics:
!K 
!K 	 DO jclass=1,nmod
!K       	sc(jclass)%ptr(1:kproma,:,krow)=zsm(1:kproma,:,jclass)
!K 	 END DO

      !--- 3) Calculate activation:

      !---3.1) Calculate the critical radius (12):

      DO jw=1, nw
         DO jclass=1, nmod
  	 WHERE (psmax(1:kproma,:,jw)>zeps 	.AND. &
  		zsm(1:kproma,:,jclass)>zeps	.AND. &
  		zn(1:kproma,:,jclass)>zeps	.AND. &
  		zrdry(1:kproma,:,jclass)>1.E-9_dp       )

  	    prc(1:kproma,:,jclass,jw)=zrdry(1:kproma,:,jclass)*(zsm(1:kproma,:,jclass)/psmax(1:kproma,:,jw))**(2._dp/3._dp)

  	 END WHERE
         END DO ! jclass
      END DO ! jw

      !--- 3.2) Calculate the fractional number of each mode
      !	   larger than the mode critical radius:


      DO jw=1,nw
         !K WRITE(*,*) ' '     
         !K WRITE(*,*) 'jw', jw
      
         DO jclass=1, nmod

    !SF Todo: this calculation should be restricted to the relevant modes only, as implemented in the
    !	 Lin & Leaitch scheme

  	 CALL ham_logtail(kproma,    kbdim,  klev,   krow, jclass, &
  			  lwetrad, ll_numb, prc(:,:,jclass,jw),	 &
  			  zfracn(:,:,jclass,jw))

        !K WRITE(*,*) 'frac', zfracn(:,:,jclass,jw)   
	   

         END DO      
      END DO


      !--- 3.3) Sum up the total number of activated particles, integrating over updraft PDF [m-3]:
      !msh: 3.3) Sum up the total number of CCN for each prescribed updraft:
      
      DO jw=1,nw
         !K WRITE(*,*) ' '      !K
         !K WRITE(*,*) 'jw', jw !K
         DO jclass=1, nmod
            pnact(1:kproma,:) = pnact(1:kproma,:)     &
                 + zfracn(1:kproma,:,jclass,jw)         &
                     * zn(1:kproma,:,jclass)
                                   

          ! ract(1:kproma,itoplev:klev,jw,jclass) = zra(1:kproma,:,jclass) 
          ! cnfrac(1:kproma,itoplev:klev,jclass,jw) = zfracn(1:kproma,:,jclass,jw)
         END DO
         pcdncact(1:kproma,:,jw) = pnact(1:kproma,:)
      END DO     


  !msh============================================================================================          
  	 !zcdncact_top(1:kproma,:,jclass)=zcdncact_top(1:kproma,:,jclass) &
  	  ! + zfracn(1:kproma,:,jclass,jw)*zn(1:kproma,:,jclass)*pwpdf(1:kproma,:,jw)
  	 !zcdncact_bot(1:kproma,:,jclass)=zcdncact_bot(1:kproma,:,jclass) &
  	  ! + pwpdf(1:kproma,:,jw)
	   

	 !K pnactnew(1:kproma,:,jclass,jw) = zcdncact_top(1:kproma,:,jclass)/zcdncact_bot(1:kproma,:,jclass) !K
	 !K WRITE(*,*) 'top', zcdncact_top(1:kproma,:,jclass), 'bot', zcdncact_bot(1:kproma,:,jclass)   !K 	   
         !K WRITE(*,*) 'jclass', jclass, 'pnactn', pnactnew(1:kproma,:,jclass,jw)  !K  

         !END DO ! jclass
      !END DO ! jw


      !DO jclass=1, nmod
       !  pnact(1:kproma,:,jclass) = zcdncact_top(1:kproma,:,jclass)/zcdncact_bot(1:kproma,:,jclass)	 
       !  pcdncact(1:kproma,:) = pcdncact(1:kproma,:) + pnact(1:kproma,:,jclass)
	!K WRITE(*,*), 'pnact', pnact(1:kproma,:,jclass), 'pcdncact', pcdncact(1:kproma,:) !K
      !END DO

      DEALLOCATE(zsum)
      DEALLOCATE(zfracn)

    END SUBROUTINE ham_activ_abdulrazzak_ghan

!K   SUBROUTINE ham_activ_diag_abdulrazzak_ghan_strat(kproma, kbdim, klev,       &
!K 						    krow,   pnact, prc,  psmax )
!K 
!K     USE mo_activ,	ONLY: nw, swat_max_strat
!K     USE mo_ham_streams, ONLY: nact_strat, rc_strat
!K 
!K     INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
!K     REAL(dp), INTENT(IN) :: pnact(kbdim,klev,nmod),  & ! number of activated particles per mode [m-3]
!K 			     prc(kbdim,klev,nmod,nw), & ! critical radius of activation per mode and w bin [m]
!K 			     psmax(kbdim,klev,nw)       ! maximum supersaturation per w bin [% 0-1]
!K 
!K     INTEGER :: jclass, jw
!K 
!K     DO jw=1, nw
!K       swat_max_strat(jw)%ptr(1:kproma,:,krow) = psmax(1:kproma,:,jw)
!K 
!K 	DO jclass=1, nmod
!K 	   WHERE(prc(1:kproma,:,jclass,jw) /= 1._dp)
!K 	     rc_strat(jclass,jw)%ptr(1:kproma,:,krow)=prc(1:kproma,:,jclass,jw)
!K 	   ENDWHERE
!K 	END DO
!K     END DO
!K 
!K     DO jclass=1, nmod
!K 	nact_strat(jclass)%ptr(1:kproma,:,krow)=pnact(1:kproma,:,jclass)
!K     END DO
!K 
!K   END SUBROUTINE ham_activ_diag_abdulrazzak_ghan_strat
!K 
!K   SUBROUTINE ham_activ_diag_abdulrazzak_ghan_conv(kproma, kbdim, klev,       &
!K 						   krow,   pnact, prc,  psmax )
!K 
!K     USE mo_activ,	ONLY: nw, swat_max_conv
!K     USE mo_ham_streams, ONLY: nact_conv, rc_conv
!K 
!K     INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
!K     REAL(dp), INTENT(IN) :: pnact(kbdim,klev,nmod),  & ! number of activated particles per mode [m-3]
!K 			     prc(kbdim,klev,nmod,nw), & ! critical radius of activation per mode and w bin [m]
!K 			     psmax(kbdim,klev,nw)       ! maximum supersaturation per w bin [% 0-1]
!K 
!K     INTEGER :: jclass, jw
!K 
!K     DO jw=1,nw
!K 	swat_max_conv(jw)%ptr(1:kproma,:,krow) = psmax(1:kproma,:,jw)
!K 	DO jclass=1, nmod
!K 	   WHERE(prc(1:kproma,:,jclass,jw) /= 1._dp)
!K 	     rc_conv(jclass,jw)%ptr(1:kproma,:,krow)=prc(1:kproma,:,jclass,jw)
!K 	   ENDWHERE
!K 	END DO
!K     END DO
!K 
!K     DO jclass=1, nmod
!K 	nact_conv(jclass)%ptr(1:kproma,:,krow)=pnact(1:kproma,:,jclass)
!K     END DO
!K 
!K   END SUBROUTINE ham_activ_diag_abdulrazzak_ghan_conv


  SUBROUTINE ham_activ_koehler_ab(kproma,   kbdim,   klev,  krow,  ktdia, &
                                  pxtm1,    ptm1,    pa,    pb            )


    ! *ham_activ_koehler_ab* calculates the Koehler A and B coefficients
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford, 2013
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

!K     USE mo_ham,         ONLY: naerocomp, aerocomp, m7mode
!K     USE mo_tracdef,     ONLY: ntrac
    USE mo_constants,   ONLY: rhoh2o, argas, amw, zsten
!K     USE mo_ham_species, ONLY: id_so4
!K     USE mo_ham_streams, ONLY: a, b
    USE mo_box_aerosols,	ONLY: nmod, ntrac, naerocomp, aermw, aerdens, aernion, aerosm, aersol	!K 

    IMPLICIT NONE

    !--- Arguments:

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow, ktdia

    REAL(dp), INTENT(IN)  :: ptm1(kbdim,klev),        & ! temperature
                             pxtm1(kbdim,klev,ntrac)	! number concentration per unit volume (number conc. * air density)



    REAL(dp), INTENT(OUT) :: pa(kbdim,klev,nmod),     & ! curvature parameter A of the Koehler equation
                             pb(kbdim,klev,nmod)        ! hygroscopicity parameter B of the Koehler equation

    !--- Local variables:

!K     INTEGER :: jclass, jt, jl, jk
    INTEGER :: jclass, jt

    REAL(dp):: zmoleweight,                &
               znion,     zosm,            &
               zrhoaer,   zeps

    REAL(dp):: zamw                          ! molecular weight of water [kg mol-1]

    REAL(dp):: zmassfrac(kbdim,klev)

    REAL(dp):: zsumtop(kbdim,klev,nmod),   & ! temporary summation field
               zsumbot(kbdim,klev,nmod)      ! temporary summation field

    REAL(dp):: zmasssum(kbdim,klev,nmod)

!K    REAL(dp), PARAMETER :: zsten = 75.0E-3_dp ! surface tension of H2O [J m-2] 

    REAL(dp) :: zafac                        ! for calculation of curvature term in Koehler equations

!K     INTEGER :: jn, jspec
    INTEGER :: jn
    INTEGER :: nion
    
    LOGICAL :: lsoluble	!K 
    
    
    zmasssum(:,:,:) = 0._dp
    pa(:,:,:)       = 0._dp
    pb(:,:,:)       = 0._dp

    zsumtop(:,:,:)  = 0._dp
    zsumbot(:,:,:)  = 0._dp

    zeps=EPSILON(1._dp)

    !--- Conversions to SI units [g mol-1 to kg mol-1]:
    
    zamw=amw*1.E-3_dp

    !---for calculation of curvature parameter A (equation (5) in first paper)
    zafac = 2._dp*zsten*zamw / (rhoh2o*argas)

    !--- Sum mass mixing ratio for each mode:
    DO jn = 1,naerocomp
!K        jclass = aerocomp(jn)%iclass
!K        jspec = aerocomp(jn)%spid
!K        jt = aerocomp(jn)%idt

	jclass = jn	!K 
	jt = jn		!K 
	
	!K zmasssum(:,:,nmod) = pxtm1(:,:,nmod) for externally mixed species

        zmasssum(1:kproma,:,jclass)=zmasssum(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)
  
      
    END DO

    !--- 1) Calculate properties for each aerosol mode:

    !--- 1.1) Calculate the auxiliary parameters A and B of the Koehler equation:

       !--- 1) Calculate weighted properties:
       !       (Abdul-Razzak & Ghan, 2000)

    DO jn=1,naerocomp
!K        zmoleweight = aerocomp(jn)%species%moleweight*1.E-3_dp         ! [kg mol-1]
!K        nion = aerocomp(jn)%species%nion
!K        znion = REAL(nion,dp)
!K        zosm = aerocomp(jn)%species%osm
!K        zrhoaer = aerocomp(jn)%species%density
!K        jt = aerocomp(jn)%idt
!K        jclass = aerocomp(jn)%iclass

	zmoleweight = aermw(jn)*1.E-3_dp! [kg mol-1]	!K 
	nion = aernion(jn)				!K 
	znion = REAL(nion,dp)				!K 
	zosm = aerosm(jn)				!K 
	zrhoaer = aerdens(jn)*1.E3_dp	! [kg m-3]	!K 
	jt = jn						!K 
	jclass = jn					!K 
        lsoluble = aersol(jn)				!K instead of using: m7mode(jclass)%lsoluble

!K       IF (nion > 0 .AND. m7mode(jclass)%lsoluble) THEN       
       IF (nion > 0 .AND. lsoluble) THEN
          WHERE(zmasssum(1:kproma,:,jclass)>zeps)
	  
             zmassfrac(1:kproma,:)=pxtm1(1:kproma,:,jt)/zmasssum(1:kproma,:,jclass) !K zmassfrac is always 1 for each of the species
      
             zsumtop(1:kproma,:,jclass)=zsumtop(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)*znion*zosm*zmassfrac(1:kproma,:)/zmoleweight
             zsumbot(1:kproma,:,jclass)=zsumbot(1:kproma,:,jclass)+pxtm1(1:kproma,:,jt)/zrhoaer
          END WHERE
	     
!K	   WRITE(*,*) ' '	   
!K	   WRITE(*,*) 'jn', jn  	   
!K	   WRITE(*,*) 'zmasssum', zmasssum(:,:,jn)
!K	   WRITE(*,*) 'zmassfrac', zmassfrac		   
!K	   WRITE(*,*) 'zsumtop', zsumtop(:,:,jn)	   
!K	   WRITE(*,*) 'zsumbot', zsumbot(:,:,jn)		   
	     

       END IF! nion>0

    END DO !naerocomp

    DO jclass=1,nmod

       WHERE (zsumbot(1:kproma,:,jclass)>zeps)

          !--- 1.1.1) Hygroscopicity parameter B (Eq. 4) [1]:

          pb(1:kproma,:,jclass)=(zamw*zsumtop(1:kproma,:,jclass))/(rhoh2o*zsumbot(1:kproma,:,jclass))

          !--- 1.1.2) Calculate the curvature parameter A [m]:

          pa(1:kproma,:,jclass)= zafac/ptm1(1:kproma,:)

       END WHERE

!K        !--- Save in stream:

!K        a(jclass)%ptr(1:kproma,:,krow)=pa(1:kproma,:,jclass)
!K        b(jclass)%ptr(1:kproma,:,krow)=pb(1:kproma,:,jclass)

    END DO !jclass=1, nmod

  END SUBROUTINE ham_activ_koehler_ab
    

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes available particules for activation
  !! 
  !! @remarks Preparatory routine for Lin&Leaitch activation scheme (HAM-specific)
  !! Basically it computes the fractional mass and number of each mode
  !! larger than the cutoff of the instrument and add them up
  !! Derived from the former aero_activ_lin_leaitch subroutine

!K  SUBROUTINE ham_avail_activ_lin_leaitch(kproma, kbdim, klev, krow, &
!K                                       prho, pxtm1)

!K    USE mo_activ,        ONLY: na, lwetrad
!K    USE mo_conv,         ONLY: na_cv
!K    USE mo_ham_tools,    ONLY: ham_logtail
!K    USE mo_tracdef,      ONLY: ntrac
!K    USE mo_ham_streams,  ONLY: frac

!K    !SF note: lwetrad is true in this routine, ie wet radius is used

!K    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
!K    REAL(dp), INTENT(IN) :: prho(kbdim,klev)        ! air density
!K    REAL(dp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac) ! tracer mmr

!K    REAL(dp), PARAMETER :: crcut=0.03*1E-6_dp ! Assumed lower cut-off of the
!K                                              ! aerosol size distribution [m]

!K    !--- Ulrike: included for activation in convective clouds
!K    REAL(dp), PARAMETER :: crcut_cv=0.02*1E-6_dp ! Assumed lower cut-off of the
!K                                                 ! aerosol size distribution in convective clouds [m]

!K    REAL(dp), PARAMETER :: cfracn(nmod)=(/1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp/)

!K    LOGICAL, PARAMETER :: lcut=.TRUE. ! explicit calculation of cut-off crcut or
!K                                      ! usage of prescribed values cfracn

!K    LOGICAL, PARAMETER :: ll_numb = .TRUE. ! switch between number/mass in logtail calculation !SF

!K    INTEGER  :: jclass, it

!K    REAL(dp) :: zr(kbdim,klev,nmod)
!K    REAL(dp) :: zfracn(kbdim,klev,nmod)
!K    REAL(dp) :: zfracn_cv(kbdim,klev,nmod)

!K    !>>dod bugfix
!K    na(1:kproma,:,krow) = 0._dp
!K    na_cv(1:kproma,:,krow) = 0._dp
!K    !<<dod

!K    !---  Calculate the fractional number of each mode
!K    !     larger than the cutoff of the instrument:
!K    IF (lcut) THEN

!K       DO jclass=1, nmod !SF the ham_logtail calculation is now done mode per mode for better efficiency
!K          IF (cfracn(jclass) > EPSILON(1._dp)) THEN !SF #279: only performs this calculation when relevant
!K             !SF stratiform:
!K             zr(1:kproma,:,jclass)=crcut
!K      
!K             CALL ham_logtail(kproma, kbdim,  klev,  krow, jclass, &
!K                              lwetrad, ll_numb, zr(:,:,jclass), zfracn(:,:,jclass) )
!K      
!K             !SF convective:
!K             zr(1:kproma,:,jclass)=crcut_cv
!K      
!K             CALL ham_logtail(kproma, kbdim,  klev,  krow, jclass, &
!K                              lwetrad, ll_numb, zr(:,:,jclass), zfracn_cv(:,:,jclass) )
!K          ELSE !SF #279: ensures that zfracn and zfracn_cv are set to 0. when mode is not relevant for activation
!K             zfracn(1:kproma,:,jclass)    = 0._dp
!K             zfracn_cv(1:kproma,:,jclass) = 0._dp
!K          ENDIF
!K       END DO
!K    ELSE

!K       DO jclass=1, nmod
!K          zfracn(1:kproma,:,jclass)    = cfracn(jclass)
!K          zfracn_cv(1:kproma,:,jclass) = cfracn(jclass) !SF was missing in previous version. Mistake??
!K       END DO

!K    END IF

!K    !--- Sum up aerosol number concentrations and convert from [kg-1] to [m-3]:
!K    DO jclass=1, nmod

!K       !>>dod soa
!K       it = m7mode(jclass)%idt_no
!K       !<<dod
!K       na(1:kproma,:,krow) = na(1:kproma,:,krow)                   &
!K                           + pxtm1(1:kproma,:,it)*prho(1:kproma,:) &
!K                             *zfracn(1:kproma,:,jclass)*cfracn(jclass)

!K       !--- Ulrike: included for NA from convection ---
!K       na_cv(1:kproma,:,krow) = na_cv(1:kproma,:,krow)                   &
!K                              + pxtm1(1:kproma,:,it)*prho(1:kproma,:)    &
!K                                *zfracn_cv(1:kproma,:,jclass)*cfracn(jclass)
!K       !--- end included

!K       frac(jclass)%ptr(1:kproma,:,krow)=zfracn(1:kproma,:,jclass)

!K    END DO

!K  END SUBROUTINE ham_avail_activ_lin_leaitch

!K  SUBROUTINE ham_activ_diag_lin_leaitch(kproma, kbdim, klev, krow, prho, pxtm1, pcdncact)

!K    USE mo_activ,       ONLY: na
!K    USE mo_conv,        ONLY: cdncact_cv
!K    USE mo_ham_streams, ONLY: frac, nact_strat, nact_conv
!K    USE mo_tracdef,     ONLY: ntrac

!K    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
!K    REAL(dp), INTENT(IN) :: prho(kbdim,klev)        ! air density
!K    REAL(dp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac) ! tracer mmr
!K    REAL(dp), INTENT(IN) :: pcdncact(kbdim,klev) ! number of activated particles

!K    INTEGER  :: jclass, it
!K    REAL(dp) :: zeps

!K    zeps = EPSILON(1._dp)

!K    DO jclass=1, nmod
!K       !>>dod soa
!K       it = m7mode(jclass)%idt_no
!K       !<<dod
!K       WHERE(na(1:kproma,:,krow)>zeps)
!K          nact_strat(jclass)%ptr(1:kproma,:,krow) = pcdncact(1:kproma,:) &
!K                                                * pxtm1(1:kproma,:,it) * prho(1:kproma,:) &
!K                                                * frac(jclass)%ptr(1:kproma,:,krow) / na(1:kproma,:,krow)

!K          nact_conv(jclass)%ptr(1:kproma,:,krow)  = cdncact_cv(1:kproma,:,krow) &
!K                                                * pxtm1(1:kproma,:,it) * prho(1:kproma,:) &
!K                                                * frac(jclass)%ptr(1:kproma,:,krow) / na(1:kproma,:,krow)
!K       ENDWHERE
!K    END DO
!K 
!K  END SUBROUTINE ham_activ_diag_lin_leaitch

END MODULE mo_ham_activ

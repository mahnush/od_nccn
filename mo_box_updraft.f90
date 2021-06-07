MODULE mo_box_updraft

 !K --- module is based on mo_activ ---

  USE mo_kind,		ONLY: dp


  IMPLICIT NONE


  INTEGER, PUBLIC, PARAMETER    :: nw =10 ! actual number of updraft velocity (w) bins 
                                      ! (can be 1 if characteristic updraft is used)

  REAL(dp), PUBLIC, PARAMETER   :: pw(nw) = (/                                  &
       0.01_dp,0.0278_dp,0.0774_dp,0.215_dp,0.599_dp,1.67_dp,4.64_dp,12.90_dp, &
       35.90_dp,100.00_dp            /)

  
!  INTEGER,  PARAMETER	 	:: nbins=20 		! number of updraft velocity bins for pdf approach
!  REAL(dp), PARAMETER 		:: w_sigma_min = 0.1_dp ! minimum value of w standard deviation [m s-1]
!  REAL(dp), PARAMETER 		:: w_min = 0.0_dp       ! minimum characteristic w for activation [m s-1]

  !K take lwetrad from mo_activ
!  LOGICAL,  PUBLIC 		:: lwetrad = .FALSE. 	!K  switch to use wet/dry radius in the activation scheme
                                      			!K  it is set to true by default (Lin & Leaitch scheme)  



  !--- Subroutines:

!CONTAINS

!  SUBROUTINE regime_updraft(clcl,nactivpdf,nactivexp,pvervel,prho,zwlarge,pw,pwpdf)

!  USE mo_box_settings,  ONLY: kproma, kbdim, klev, krow
!  USE mo_constants,	ONLY: g
   
!   IMPLICIT NONE

    ! Input parameters
!    INTEGER, INTENT(in)		:: nactivpdf,nactivexp
!    REAL(dp), INTENT(in)	:: clcl(kbdim,klev),	&
!    				   pvervel(kbdim,klev),	&	! vertical velocity [m/s] from Zheng & Rosenfeld (2015)
!    				   prho(kbdim,klev)

    
    ! Output parameters
    !REAL(dp), INTENT(out), ALLOCATABLE, DIMENSION(:,:,:) :: pw, pwpdf 
    !REAL(dp), INTENT(out) :: zwlarge(kbdim,klev) ![m s-1]
    
    ! Local parameters    
    !REAL(dp) :: zwturb(kbdim,klev) ![m s-1]


    !--- Determine mean vertical velocity, either prescribed or from file ---

    	!Large scale vertical velocity in SI units: Convert from [Pa s-1] to [m s-1]
    	!zwlarge(1:kproma,:) = -1._dp* pvervel(kbdim,klev)/(g*prho(kbdim,klev))    

    	!Vertical velocity [m/s] from Zheng & Rosenfeld (2015)
	!zwlarge(1:kproma,:) = pvervel(kbdim,klev)
    
    	!zwlarge does not change here [m s-1]
    	!zwlarge(1:kproma,:) = w_min   

    
    !--- Determine pdf width of vertical velocity spectrum (sigma), prescribed here as regime-dependent---    

    !IF(nactivexp==1) THEN    
    !--- constant sigma
     ! zwturb(:,:) = 0.4_dp
    !ENDIF

    !IF(nactivexp==2) THEN    
    !--- regime dependence: inhomo vs homo
      !IF (clcl(kbdim,klev) >= 1.0 .AND. clcl(kbdim,klev) < 2.0) THEN  !1 homo 
      ! zwturb(:,:) = 0.4_dp	 
      !ELSEIF (clcl(kbdim,klev) >= 2.0) THEN	!2 inhomo
      ! zwturb(:,:) = 0.8_dp    
      !ELSE 
      ! zwturb(:,:) = w_sigma_min	        ! no regime found      
      !ENDIF 
    !ENDIF

    !IF(nactivexp==3) THEN    
    !--- full regime dependence     
     ! IF (clcl(kbdim,klev) == 1.1) THEN 	!1 low homo
      ! zwturb(:,:) = 0.2_dp	 
      !ELSEIF (clcl(kbdim,klev) == 1.2) THEN	!2 med homo
      ! zwturb(:,:) = 0.4_dp    
      !ELSEIF (clcl(kbdim,klev) == 1.3) THEN	!3 high homo
       !zwturb(:,:) = 0.6_dp   
      !ELSEIF (clcl(kbdim,klev) == 2.1) THEN	!4 low inhomo
       !zwturb(:,:) = 0.4_dp   
      !ELSEIF (clcl(kbdim,klev) == 2.2) THEN	!5 med inhomo
      ! zwturb(:,:) = 0.8_dp   
      !ELSEIF (clcl(kbdim,klev) == 2.3) THEN	!6 high inhomo
      ! zwturb(:,:) = 1.2_dp		      
      !ELSE				        ! no regime found
       !zwturb(:,:) = w_sigma_min	    
      !ENDIF		
    !ENDIF
       
!    WRITE(*,*) 'clcl', clcl 
!    WRITE(*,*) 'zwlarge [m/s]', zwlarge(:,:)
!    WRITE(*,*) 'zwturb', zwturb(:,:)
    

    !--- Set number of updraft bins: single characteristic, or as pdf
    !IF(nactivpdf==1) THEN
     ! nw=nbins
    !ELSE
     ! nw=1
    !ENDIF
    

    !ALLOCATE(pw(kbdim,klev,nw))
    !ALLOCATE(pwpdf(kbdim,klev,nw))
    
    !--- compute single or pdf of updraft:         
    !SELECT CASE (nactivpdf)
    !CASE (0)
     !  pwpdf(1:kproma,:,nw) = 1.0_dp

       !use mean updraft velocity
      ! pw(1:kproma,:,nw) = zwlarge(1:kproma,:) !in this case it would mean updraft = 0!! but not used here
       						!or use zwlarge as read in from file

       !or use sigma as mean updraft   
       !pw(1:kproma,:,nw) = zwturb(1:kproma,:) 
  
    !CASE (1)  
      !CALL aero_activ_updraft_pdf(kproma,  kbdim,  klev,  krow, &
      !                            zwlarge, zwturb, pw,    pwpdf )      
      
    !END SELECT

      
  !END SUBROUTINE regime_updraft




!-------------------------------------------------------------------------------

  !SUBROUTINE aero_activ_updraft_pdf(kproma,  kbdim,   klev, krow, &
    !                                pwlarge, pwsigma, pw,   pwpdf )

    ! *aero_activ_updraft_* calculates Gaussian pdf of  
    !                       updraft vertical velocity
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2013
    !
    ! References:
    ! -----------
    ! West et al., ACP, 2013. 
    !

    !USE mo_constants, ONLY: api

    !IMPLICIT NONE

    !INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    !REAL(dp), INTENT(in)  :: pwlarge(kbdim,klev), & ! large-scale vertical velocity [m s-1]
     !                        pwsigma(kbdim,klev)    ! st. dev. of vertical velocity [m s-1]

    !REAL(dp), INTENT(out) :: pw(kbdim,klev,nw),   & ! vertical velocity bins [m s-1]
     !                        pwpdf(kbdim,klev,nw)   ! vertical velocity PDF [s m-1]

    !INTEGER               :: jw  !K , jl, jk

    !REAL(dp)              :: zw_width(kbdim,klev), &
     !                        zw_min(kbdim,klev), &
      !                       zw_max(kbdim,klev)

    !zw_min(1:kproma,:) = 0.0_dp
    !zw_max(1:kproma,:) = 2.0_dp*pwsigma(1:kproma,:) ! gives flexible range, harder to compare in the end
    !zw_max(1:kproma,:) = 2.0_dp ! gives fixed range (larger than biggest sigma values), easier to compare in the end, with 20 steps of 0.1 m/s

    !zw_min(1:kproma,:) = pwlarge(1:kproma,:)-4.0_dp*pwsigma(1:kproma,:)
    !zw_max(1:kproma,:) = pwlarge(1:kproma,:)+4.0_dp*pwsigma(1:kproma,:)
    !IF(zw_min(kbdim,klev) < 0.0_dp) THEN
    !  zw_min(1:kproma,:) = 0.0_dp
    !ENDIF
    
    !zw_width(1:kproma,:) = (zw_max(1:kproma,:) - zw_min(1:kproma,:)) / DBLE(nw) 

    !DO jw=1, nw
     ! pw(1:kproma,:,jw) = zw_min(1:kproma,:) + (DBLE(jw) - 0.5_dp) * zw_width(1:kproma,:)
!K       w(jw)%ptr(1:kproma,:,krow) = pw(1:kproma,:,jw)

      !pwpdf(1:kproma,:,jw) = (1.0_dp / ((2.0_dp*api)**0.5_dp))              &
       !                      * (1.0_dp / pwsigma(1:kproma,:))                    &
        !                     * EXP( -((pw(1:kproma,:,jw) - pwlarge(1:kproma,:))**2_dp &
         !                           / (2.0_dp*pwsigma(1:kproma,:)**2.0_dp)) )

!K       w_pdf(jw)%ptr(1:kproma,:,krow) = pwpdf(1:kproma,:,jw)
   ! END DO

  !END SUBROUTINE aero_activ_updraft_pdf



END MODULE mo_box_updraft

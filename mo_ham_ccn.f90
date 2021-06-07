MODULE mo_ham_ccn

  USE mo_kind,          ONLY: dp
!K  USE mo_linked_list,   ONLY: t_stream
!K  USE mo_ham_m7ctl,     ONLY: nmod
  USE mo_box_aerosols,	ONLY: nmod	!K 
!K  USE mo_submodel_diag, ONLY: vmem3d,vmem2d
!K  USE mo_param_switches,ONLY: nccndiag
  USE mo_ham_ccnctl,    ONLY: nsat, zsat

  
  

  IMPLICIT NONE

  PUBLIC ham_ccn
 
  PRIVATE

CONTAINS

!K   SUBROUTINE ham_ccn(kproma, kbdim, klev, klevp1, krow,  ktrac, &
!K                      pxtm1,  ptvm1, papm1, paphm1		   )			

  SUBROUTINE ham_ccn(kproma, kbdim, klev, krow,  ktrac, &	!K 
  		     zn, za, zb, ract, cnfrac, ccnsum	)	!K 


    ! *ham_ccn* calculates the number Cloud Condensation Nuclei
    !           CCN at fixed supersaturation x
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford  2009
    !
    ! Method:
    ! -------
    ! The calculation of the CCNx can be reduced to 3 tasks:
    ! 
    ! I)   Use prescribed supersaturation
    ! II)  Calculate the corresponding radius of activation
    !      for each mode
    ! III) Calculate the number of particles that are larger
    !      then the radius of activation for each mode.
    ! 
    ! III) Calculation of the number of activated particles:
    !      See the routine aero_activ_tail below.
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

!K     USE mo_ham,          ONLY: m7mode
!K     USE mo_ham_m7ctl,    ONLY: nmod
!K     USE mo_ham_tools,    ONLY: ham_logtail
    USE mo_box_ham_tools,    ONLY: ham_logtail
!K     USE mo_ham_streams,  ONLY:  a, b, &
!K                                ccn_2d, ccn_burden, ccn_3d, &
!K                                cn_2d, cn_burden, cn_3d
    USE mo_constants,    ONLY: g, rd


    IMPLICIT NONE

    !--- Arguments:

!K     INTEGER, INTENT(in) :: kproma, kbdim, klev, klevp1, krow, ktrac
    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow, ktrac  
 
!K     REAL(dp), INTENT(in) :: ptvm1(kbdim,klev),  	& 	! virtual temperature
!K                             papm1(kbdim,klev),  	& 	! pressure 
!K                             paphm1(kbdim,klevp1)  		! half-level pressure 

!K     REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac)	
    REAL(dp), INTENT(in) :: zn(kbdim,klev,nmod), &	!K read in zn instead of pxtm1
    			    za(kbdim,klev,nmod), &	!K  ! A & B coefficients from ham_activ_koehler_ab
			    zb(kbdim,klev,nmod)		!K   instead of a & b from mo_ham_streams


    !--- Output variables:					!K 

    REAL(dp), INTENT(out) :: ract(kbdim,klev,nsat,nmod)		!K 
    REAL(dp), INTENT(out) :: cnfrac(kbdim,klev,nsat,nmod)   	!K  
    REAL(dp), INTENT(out) :: ccnsum(kbdim,klev,nsat)		!K 



    !--- Local variables:

!K     INTEGER :: jclass, jt, jk
    INTEGER :: jclass

    REAL(dp):: zeps, ztmp

!K     REAL(dp):: zrho(kbdim,klev)	 & ! air density
!K                zdpg(kbdim,klev)              ! auxiliary variable dp/g

    REAL(dp):: zra(kbdim,klev,nmod),       & ! radius of activation, i.e. radius of smallest particle with 
                                             ! critical supersaturation larger than prescribed supersaturation [m] 
               zfracn(kbdim,klev,nmod)!K,    & ! fraction of activated aerosol numbers for each mode - stratiform
!K               zn(kbdim,klev,nmod)           ! aerosol number concentration for each mode [m-3]

    REAL(dp):: zccn(kbdim,klev) !K ,           & ! CCN at supersaturations zsat(jsat) [m-3]
!K                zcn(kbdim,klev)               ! CN(r>zrmin) at ambient conditions [m-3]

!K     REAL(dp) :: zrmin(kbdim,klev,nmod)       ! cut-off radius for CN calculation (typical measurement limit 3nm)

    LOGICAL, PARAMETER :: lwetrad = .FALSE.   ! switch between wet/dry radii in logtail calculation

    LOGICAL, PARAMETER :: ll_numb = .TRUE.  ! switch between number/mass in logtail calculation

    INTEGER :: jsat

    INTEGER :: itoplev                        ! highest level considered: ilev=klev for 2D surface diagnostics
                                              !                           ilev=1 for full 3D diagnostics
    

    !--- 0) Initializations:

!K     zrmin(1:kproma,:,:)=3.0E-9_dp             ! 3nm cut-off radius for CN calculation

    zeps=EPSILON(1._dp)


!K     IF (nccndiag==1 .OR. nccndiag==3) THEN
       itoplev=klev
!K     ELSE
!K        itoplev=1
!K 
!K        IF (nccndiag==5 .OR. nccndiag==6) THEN
!K           !--- Calculate zdpg:
!K           zdpg(1:kproma,1)=2._dp*(paphm1(1:kproma,2)-papm1(1:kproma,1))/g
!K           zdpg(1:kproma,2:klev)=(paphm1(1:kproma,3:klev+1)-paphm1(1:kproma,2:klev))/g
!K        END IF
!K     END IF

!K     !--- Number per unit volume for each mode:
!K 
!K     zrho(1:kproma,itoplev:klev) = papm1(1:kproma,itoplev:klev)/(rd*ptvm1(1:kproma,itoplev:klev))


!K     WRITE(*,*) ' '!K 
!K     WRITE(*,*) 'CLOUD CONDENSATION NUCLEI:'!K 


!K     DO jclass=1, nmod
!K        jt = m7mode(jclass)%idt_no
!K        zn(1:kproma,itoplev:klev,jclass)=pxtm1(1:kproma,itoplev:klev,jt)*zrho(1:kproma,itoplev:klev)
!K     END DO



    !--- 1) Calculate properties for each aerosol mode:
    !--- 1.1) Calculate the auxiliary parameters A and B of the Koehler equation:
    !       Now done in ham_activ_koehler_ab once so that they can be used in
    !       convective and stratiform activation 


    !--- 2) Calculate the radius of activation [m], i.e. radius of smallest particle 
    !       with critical supersaturation larger than prescribed supersaturations [%]:
    !       AR et al. (2000) (Eq.4)

    DO jsat=1, nsat

!K        zcn(1:kproma,itoplev:klev)  = 0._dp
       zccn(1:kproma,itoplev:klev) = 0._dp

       ztmp=(2._dp/zsat(jsat))**(2._dp/3._dp)
       
!K       WRITE(*,*) ' '
!K       WRITE(*,*) 'zsat', zsat(jsat) !K 
       
       DO jclass=1, nmod

          zra(1:kproma,:,jclass) = 1.0_dp ! Set to large value of 1.0 m to set CCN count to zero

!K           WHERE(b(jclass)%ptr(1:kproma,itoplev:klev,krow)>zeps)
!K              zra(1:kproma,itoplev:klev,jclass)                                          &
!K                = a(jclass)%ptr(1:kproma,itoplev:klev,krow)                              &
!K                    / (3._dp * b(jclass)%ptr(1:kproma,itoplev:klev,krow)**(1._dp/3._dp)) &
!K                    * ztmp

	  !K compute the critical radius for a specific supersat, depending on Koehler coeff.
	  WHERE(zb(1:kproma,itoplev:klev,jclass)>zeps)						!K 
	     zra(1:kproma,itoplev:klev,jclass)  					&	!K 
	       = za(1:kproma,itoplev:klev,jclass)					&	!K 
		   / (3._dp * zb(1:kproma,itoplev:klev,jclass)**(1._dp/3._dp)) 		&	!K 
		   * ztmp									!K 
          END WHERE										!K 


         !--- 3) Calculate the fractional number of each mode
         !       larger than the mode critical radius:

         CALL ham_logtail(kproma,  kbdim,   klev,   krow,   jclass, &
                          lwetrad, ll_numb, zra(:,:,jclass),        &
                          zfracn(:,:,jclass)                        )
       
         !--- 4) Sum up the total number of CCN [m-3]:

         zccn(1:kproma,itoplev:klev) = zccn(1:kproma,itoplev:klev)              &
                                         + zfracn(1:kproma,itoplev:klev,jclass) &
                                             * zn(1:kproma,itoplev:klev,jclass)

 
!K          WRITE(*,*) 'jclass', jclass, 'zra[m]', zra(:,:,jclass), 'zfracn[%]', zfracn(:,:,jclass)*100., 'zccn[/m3]', zccn(:,:)	!K 

         ract(1:kproma,itoplev:klev,jsat,jclass) = zra(1:kproma,:,jclass)		!K 
	 cnfrac(1:kproma,itoplev:klev,jsat,jclass) = zfracn(1:kproma,:,jclass)		!K 


       END DO ! jclass
   
         ccnsum(1:kproma,itoplev:klev,jsat) = zccn(1:kproma,itoplev:klev)			!K 

!K        !--- Store diagnostics in aero stream:

!K        IF (nccndiag==1 .OR. nccndiag==3 .OR. nccndiag==5) THEN

!K           ccn_2d(jsat)%ptr(1:kproma,krow)=zccn(1:kproma,klev)

!K        ELSE IF (nccndiag==2 .OR. nccndiag==4 .OR. nccndiag==6) THEN

!K           ccn_3d(jsat)%ptr(1:kproma,itoplev:klev,krow)=zccn(1:kproma,itoplev:klev)

!K        END IF

!K        IF (nccndiag==5 .OR. nccndiag==6) THEN

!K           DO jk=itoplev, klev
!K              ccn_burden(jsat)%ptr(1:kproma,krow)=ccn_burden(jsat)%ptr(1:kproma,krow) + &
!K                                                  zccn(1:kproma,jk)/zrho(1:kproma,jk)*zdpg(1:kproma,jk)
!K           END DO ! jk

!K        END IF

    END DO ! jsat

!K     IF (nccndiag > 2) THEN 

!K        DO jclass=1, nmod
          !--- 3) Calculate the fractional number of each mode
          !       larger than the minimum radius zrmin

!K           CALL ham_logtail(kproma,  kbdim,   klev,   krow,   jclass, &
!K                            lwetrad, ll_numb, zrmin(:,:,jclass),      &
!K                            zfracn(:,:,jclass)                        )

 

          !--- 4) Sum up the total number of CCN [m-3]:

!K           zcn(1:kproma,itoplev:klev) = zcn(1:kproma,itoplev:klev)               &
!K                                          + zfracn(1:kproma,itoplev:klev,jclass) &
!K                                              * zn(1:kproma,itoplev:klev,jclass)
!K        END DO ! jclass

!K     WRITE(*,*) 'zrmin', zrmin, 'zfracn', zfracn, 'zcn', zcn	!K

!K       !--- Store diagnostics in aero stream:

!K       IF (nccndiag==3 .OR. nccndiag==5) THEN

!K          cn_2d(1:kproma,krow)=zcn(1:kproma,klev)

!K       ELSE IF (nccndiag==4 .OR. nccndiag==6) THEN

!K          cn_3d(1:kproma,itoplev:klev,krow)=zcn(1:kproma,itoplev:klev)

!K       END IF

!K 	!--- 4) Sum up the CN over levels to burdens:

!K 	IF (nccndiag==5 .OR. nccndiag==6) THEN

!K 	   DO jk=itoplev, klev
!K 	      cn_burden(1:kproma,krow)=cn_burden(1:kproma,krow) + &
!K 				       zcn(1:kproma,jk)/zrho(1:kproma,jk)*zdpg(1:kproma,jk)
!K 	   END DO ! jk

!K 	END IF

!K     END IF

  END SUBROUTINE ham_ccn

END MODULE mo_ham_ccn

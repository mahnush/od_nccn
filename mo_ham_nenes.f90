  MODULE mo_ham_nenes
! *****************************COPYRIGHT*******************************
!
! In progress
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Daniel Partridge: Converted from F77 to F90, 
!			     optimised & modified Nenes parameterisation to
!                            work in current UKCA framework. Based on
!			     DPs 1D box model of Nenes Param.	     
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
!  Module containing subroutines
!
!  Author: Daniel Partridge, Oxford, 2012
!  -------
!  
!  Question: add more detailed description
!  Purpose: Contains subroutines for Nenes droplet activation parameterisation, changes found with MOD_CHANGES/Dan
!  --------  
!           
!
!  Externals:
!  ----------
!  
!
!  References:
!  -----------
!  
!
!*** -------------------------------------------------------------***!

!K  USE mo_kind,       ONLY: dp	!K is not used anyway
!K  USE mo_ham_m7ctl,  ONLY: nmod
  USE mo_box_aerosols,	ONLY: nmod	!K 
!K  USE mo_constants,	ONLY: amd, g, argas, api	!K 

 !   IMPLICIT NONE

  PUBLIC CCNSPEC
  PUBLIC PDFACTIV
  PUBLIC PROPS

REAL, DIMENSION(nmod) :: MODE

PUBLIC MODE

PRIVATE

INTEGER, PARAMETER  :: NpGauss=10
INTEGER, PARAMETER  :: MAXIT = 100		! Max iterations for solution
INTEGER :: NMD

!REAL             :: XGS, WGS
LOGICAL          :: CRIT2, CCNSPST
CHARACTER (len=9):: VERSION

!REAL, PARAMETER :: AMA = 29.e-3			! Air molecular weight
!REAL, PARAMETER :: GRAV = 9.81                  ! g constant
!REAL, PARAMETER :: RGAS = 8.31                  ! Universal gas constant
!CHANGES

REAL, PARAMETER :: AMA = 28.97e-3		     ! Air molecular weight
REAL, PARAMETER :: GRAV = 9.80665		   ! g constant
REAL, PARAMETER :: RGAS = 8.314409		    ! Universal gas constant

!K use instead from mo_constants (dp)?	!K
!K REAL, PARAMETER :: AMA   = amd  	!K 
!K REAL, PARAMETER :: GRAV  = g		!K 
!K REAL, PARAMETER :: RGAS  = argas	!K

REAL, PARAMETER :: EPS = 1e-6			! Convergence criterion
!
!REAL, PARAMETER :: PI = 3.1415927		! Some constants
REAL, PARAMETER :: PI = 3.14159265358979323846  ! Some constants

!K use instead from mo_constants (dp)?	!K
!K REAL, PARAMETER :: PI = api	!K     

REAL, PARAMETER :: ZERO = 0.
REAL, PARAMETER :: GREAT = 1.e30
REAL, PARAMETER :: SQ2PI = 2.5066282746
!
!     CHANGES BEGIN
!     DATA FOR FHH
REAL, PARAMETER :: Dw = 2.75e-10		! Water Molecule Diameter

!PK   Data for Exponent Calculation
      ! for C1     
REAL, PARAMETER :: D11 = -0.1907
REAL, PARAMETER :: D12 = -1.6929
REAL, PARAMETER :: D13 = 1.4963
REAL, PARAMETER :: D14 = -0.5644 
REAL, PARAMETER :: D15 = 0.0711
      ! for C2
REAL, PARAMETER :: D21 = -3.9310
REAL, PARAMETER :: D22 = 7.0906
REAL, PARAMETER :: D23 = -5.3436
REAL, PARAMETER :: D24 = 1.8025 
REAL, PARAMETER :: D25 = -0.2131
      ! for C3
REAL, PARAMETER :: D31 = 8.4825
REAL, PARAMETER :: D32 = -14.9297
REAL, PARAMETER :: D33 = 11.4552
REAL, PARAMETER :: D34 = -3.9115 
REAL, PARAMETER :: D35 = 0.4647
      ! for C4
REAL, PARAMETER :: D41 = -5.1774
REAL, PARAMETER :: D42 = 8.8725
REAL, PARAMETER :: D43 = -6.8527
REAL, PARAMETER :: D44 = 2.3514 
REAL, PARAMETER :: D45 = -0.2799

!/INPUTS/
REAL, SAVE:: TEMP, PRES

!/CCNSPC/
REAL, SAVE ::DPG(nmod), SIG(nmod), VHF(nmod), AMS(nmod), DENS(nmod),&
                   & DENI(nmod), AMFS(nmod), Dpc(nmod), &
                   & TP_b(nmod), AKK_b(nmod), ei_b(nmod)

!/ACTVPR/
!ORIGINAL
REAL, SAVE :: AKOH, SSPLT, ALFA, BET1

!/THERMO/
REAL, SAVE :: AMW, DENW, CPAIR, DHV, AKA, PSAT,&
                    & DAIR, SURT, XFHH, WPARCEL_b, A_b, B_b, ACCOM_b

!/PARDBG/
REAL, SAVE :: WPDBG(NpGauss), PDDBG(NpGauss), XNADBG(NpGauss), SMDBG(NpGauss)

!/SLNPAR/
REAL, SAVE :: NITER

!/GAUSSL/
REAL, SAVE ::XGS(NpGauss), WGS(NpGauss)

!PK   Changes finish for Exponent

!
      DATA CCNSPST /.FALSE./               ! Internal consistency check

  CONTAINS
  

!=======================================================================
!
! *** SUBROUTINE CCNSPEC
! *** THIS SUBROUTINE CALCULATES THE CCN SPECTRUM OF THE AEROSOL USING
!     THE APPROPRIATE FORM OF KOHLER THEORY
!
! *** ORIGINALLY WRITTEN BY ATHANASIOS NENES FOR ONLY KOHLER PARTICLES
! *** MODIFIED BY PRASHANT KUMAR AND ATHANSIOS NENES TO INCLUDE 
! *** ACTIVATION BY FHH PARTICLES

!=======================================================================
!
      SUBROUTINE CCNSPEC(TPI,DPGI,SIGI,TPARC,PPARC,NMOD,AKK,ei,A,B,SG,za,zb,jmod,PSAT_IN,zk_IN)

!
!      INCLUDE 'parametr.inc'
!      USE mod_parameters
!      REAL DPGI(NMOD), VHFI(NMOD), AMSI(NMOD),&
!                       DENSI(NMOD),SIGI(NMOD), TPI(NMOD),&
!                       AMFSI(NMOD),DENII(NMOD),TP(nmod),&
!                       TPARC, PPARC,AKK(nmod),ei(nmod),SG(nmod),A,B,zb(nmod),za(nmod)             
!
      REAL :: TPI
      REAL :: DPGI
      REAL :: SIGI
      REAL :: TPARC
      REAL :: PPARC
      REAL :: AKK
      REAL :: ei
      REAL :: A
      REAL :: B
      REAL :: SG
      REAL :: za
      REAL :: zb
      INTEGER :: jmod	! mode number
      !REAL :: PSAT_IN
      !REAL :: zk_IN

! local variables
      REAL :: Dpcm
      REAL :: DPG
      REAL :: SIG
      REAL :: TP
      REAL :: Dpc

! Dan: remove the loop over modes
      NMD  = NMOD                ! Save aerosol params in COMMON
!      DO I=1,NMD

!PK   BEGIN CHANGES
!
!PK   Check for Kohler Modes
!      IF (MODE(jmod).EQ.1) THEN
         !DPG(I) = DPGI(I)
         DPG = 2.*DPGI
         SIG = SIGI
         TP  = TPI
!PK   Finish Checking for Kohler Modes

!PK   Check for FHH Modes 
!      ELSEIF (MODE(jmod).EQ.2) THEN
!         DPG(I) = DPGI(I)/ei(I)**(1./3.)
!         DPG = 2.*DPGI/ei**(1./3.)
!         SIG = SIGI
!    	 TP  = TPI
!PK   Finish Checking for FHH Modes

!      ENDIF

! Dan      
!      ENDDO
!
      TEMP = TPARC                 ! Save parcel props in COMMON
      PRES = PPARC
!
      CALL PROPS                   ! Calculate thermophysical properties

PSAT = PSAT_IN
AKA = zk_IN

!
! *** Calculate critical properties
!
      AKOH = 4*AMW*SURT/RGAS/TEMP/DENW   ! Curvature param


!print*, 'PSAT_IN', PSAT_IN
!print*, 'zk_IN', zk_IN
!print*, 'AKOH', AKOH
!print*, 'zb', zb
!print*, 'za',za

! Dan: remove the loop over modes
!      DO K=1,NMD
!PK   BEGIN CHANGES

!PK   Solve for Kohler Modes     
!      IF (MODE(jmod).EQ.1) THEN

!PK   KOHLER ACTIVATION
!         PAR1 = 4./27./AKK(K)/DPG(K)**3 !VAK KK
!         PAR2 = SQRT(PAR1*AKOH**3)

!Dan (want to use za/zb as get these from ARG: mo_aero_activ.f90)
PAR1 = 4./27./zb/(DPG)**3 !VAK KK
PAR2 = SQRT(PAR1*AKOH**3)
!Gives same result as:
!PAR2 = SQRT(PAR1*(2*za)**3)

!WRITE(*,*) 'za', za
!WRITE(*,*) 'zb', zb

        SG = EXP(PAR2) - 1              ! Sc of Dpg

!WRITE(*,*) 'SG', SG        

!PK   Finish solving for Kohler Modes


!PK   Solve for FHH Modes     
 !     ELSEIF (MODE(jmod).EQ.2) THEN
!PK   FHH ACTIVATION
!         CALL DpcFHH_2(DPG,TPARC,AKK,ei,A,B,Dpcm)
!         Dpc = Dpcm

!PK   Calculating Critical Super Saturation by Taylor Series Expansion 
!         SG = (AKOH/Dpc)+ &
!         (-(1.-ei)*DPG**3.*AKK/(Dpc**3.-ei*DPG**3.))+ &
!         (-A*(((Dpc-ei**(1./3.)*DPG)/(2.*Dw))**(-B)))

!PK   Finish solving for FHH Modes
!      ENDIF

!      ENDDO   
!
! *** END OF SUBROUTINE CCNSPEC ****************************************
!
      RETURN
      END SUBROUTINE CCNSPEC
!=======================================================================


!K !VK   Begin Changes
!K !=======================================================================
!K !
!K ! *** SUBROUTINE DpcFHH
!K ! *** THIS SUBROUTINE CALCULATES THE CRITICAL PARTICLE DIAMETER
!K !     ACCORDING TO THE FHH ADSOSPRTION ISOTHERM THEORY.
!K !
!K ! *** WRITTEN BY PRASHANT KUMAR AND ATHANASIOS NENES
!K ! *** UPDATED BY VLASSIS KARYDIS AND ATHANASIOS NENES
!K !
!K !=======================================================================
!K !
!K       SUBROUTINE DpcFHH_2(Ddry,TPARC,dei,dhk,A,B,Dc)
!K 
!K !      Include 'parametr.inc'
!K !      USE mod_parameters
!K       REAL:: Ddry,mu0,mu,mu1,mu2,mu3,X1,X2l,Dpcm,Dpcl,Dpcu,&
!K       X3,F1,F2,X3l,X2u,X3u,FDpcl,FDpcu,FDpcm,X2m,X3m,Dc
!K 
!K       TEMP = TPARC
!K 
!K       CALL PROPS
!K 
!K 
!K             Dpcl = Ddry         !Lower Limit
!K             Dpcu = 100.*Ddry     !Upper Limit
!K 
!K 
!K  100     FDpcl = -AKOH/Dpcl**2. + &
!K          3.*(1.-dei)*Ddry**3.*dhk*Dpcl**2./(Dpcl**3.-dei*Ddry**3.)**2. + &
!K          (A*B/(2.*Dw))*(((Dpcl-dei**(1./3.)*Ddry)/ &
!K          (2.*Dw))**(-B-1.))
!K 
!K 
!K          FDpcu = -AKOH/Dpcu**2. + &
!K          3.*(1.-dei)*Ddry**3.*dhk*Dpcu**2./(Dpcu**3.-dei*Ddry**3.)**2. + &
!K          (A*B/(2.*Dw))*(((Dpcu-dei**(1./3.)*Ddry)/ &
!K          (2.*Dw))**(-B-1.))
!K 
!K 
!K             Dpcm = (Dpcu+Dpcl)/2
!K 
!K 
!K          FDpcm = -AKOH/Dpcm**2. + &
!K          3.*(1.-dei)*Ddry**3.*dhk*Dpcm**2./(Dpcm**3.-dei*Ddry**3.)**2. + &
!K          (A*B/(2.*Dw))*(((Dpcm-dei**(1./3.)*Ddry)/ &
!K          (2.*Dw))**(-B-1.))
!K 
!K             If ((FDpcl*FDpcm).Le.0) Then
!K 
!K        	If (ABS(FDpcm).Le.10e-8) Then
!K                   Goto 200
!K        	Else
!K                    Dpcl = Dpcl
!K                    Dpcu = Dpcm
!K                    goto 100
!K        	End if
!K 
!K             Else If ((FDpcl*FDpcm).GE.0) Then
!K 
!K                 If (ABS(FDpcm).Le.10e-8) Then
!K                    Goto 200
!K                 Else
!K                     Dpcl = Dpcm
!K                     Dpcu = Dpcu
!K                     goto 100
!K                 End if
!K 
!K             Else If ((FDpcl*FDpcm).Eq.0) Then
!K                     Goto 200
!K             End if
!K 
!K 200   Dc = Dpcm
!K 
!K       RETURN
!K       END SUBROUTINE DpcFHH_2
!K ! *** END OF SUBROUTINE DpcFHH_2 ***************************************
!K 

!=======================================================================
!
! *** SUBROUTINE PDFACTIV
! *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!     and Nenes (2004). THIS ROUTINE CALCULATES FOR A PDF OF
!     UPDRAFT VELOCITIES.
!
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE PDFACTIV (WPARC,TP,AKK,ei,A,B,ACCOM,SG,SIGW,TPARC, &
!MOD_CHANGES
       PPARC,NACT,SMAX,NMOD,SIG)
!
!      INCLUDE 'parametr.inc'
!      USE mod_parameters
!      REAL:: TPART, NACT, NACTI,A,B,ACCOM, &
!       TP(nmod),AKK(nmod),ei(nmod),SG(nmod)
!MOD_CHANGES
      INTEGER NMOD
      REAL:: TPART, NACT, NACTI,A,B,ACCOM, &
       TP(NMOD),AKK(NMOD),ei(NMOD),SG(NMOD),SIG(NMOD)
      REAL::             PDF
  !
! *** Case where updraft is very small
!
      
      IF (WPARC.LE.1.e-6) THEN
         SMAX  = 0.
         NACT  = 0.
         ISEC  = 1
         DPNMX = GREAT
         RETURN
      ENDIF
!
! *** Single updraft case
! 
! MOD_CHANGES
! We always have a single updraft case here

!      IF (SIGW.LT.1e-10) THEN
         CALL ACTIVATE (WPARC,TP,AKK,ei,A,B,ACCOM,SG,NACT,SMAX,SIG)
         WPDBG(1) = WParc                      ! Save debug info
         PDDBG(1) = 1.0
         XNADBG(1) = NACT
         SMDBG(1) = SMAX
!
! *** PDF of updrafts
!
!      ELSE
!         NACT  = ZERO
!         SMAX  = ZERO
!         PLIMT = 1e-3     ! Probability of High Updraft limit
!         PROBI = SQRT(-2.0*LOG(PLIMT*SIGW*SQ2PI))
!         WHI   = WPARC + SIGW*PROBI             ! Upper updrft limit
!         WLO   = 0.05  ! WPARC - SIGW*PROBI     ! Low updrft limit
!         SCAL  = 0.5*(WHI-WLO)                  ! Scaling for updrafts
!         DO I=1,NpGauss
!            WPI  = WLO + SCAL*(1.0-XGS(i))      ! Updraft
!            CALL ACTIVATE (WPI,TP,AKK,ei,A,B,ACCOM,SG,NACTI,SMAXI,SIG)     ! # of drops
!            PDF  = (1.0/SQ2PI/SIGW)*EXP(-0.5*((WPI-WPARC)/SIGW)**2) ! Prob. of updrafts
!            NACT = NACT + WGS(i)*(PDF*NACTI)    ! Integral for drops
!            SMAX = SMAX + WGS(i)*(PDF*SMAXI)    ! Integral for Smax
!            WPDBG(I) = WPI                      ! Save debug info
!            PDDBG(I) = PDF
!            XNADBG(I) = NACTI
!            SMDBG(I) = SMAXI
!            IF (PDF.LT.PLIMT) GOTO 100
!         ENDDO
! 100     NACT = NACT*SCAL                       ! Scale Integrals
!         SMAX = SMAX*SCAL
!      ENDIF
!

      RETURN
!
! *** END OF SUBROUTINE PDFACTIV ***************************************
!
      END SUBROUTINE PDFACTIV
!=======================================================================


!K !=======================================================================
!K !
!K ! *** SUBROUTINE WACTIV
!K ! *** THIS SUBROUTINE CALCULATES THE UPDRAFT NECESSARY TO ACHIEVE A DROP
!K !     CONCENTRATION.
!K !
!K ! *** WRITTEN BY ATHANASIOS NENES
!K !
!K !=======================================================================
!K !
!K       SUBROUTINE WACTIV (NACT,TP, WPARC, SMAX)
!K !
!K !      INCLUDE 'parametr.inc'
!K !      USE mod_parameters
!K       REAL:: TPART,NACT,NACT1,NACT2,NACT3, &
!K 	TP(nmod),AKK(nmod),SG(nmod)
!K !MOD_CHANGES
!K       REAL:: ei(nmod)
!K  !
!K ! *** INITIAL VALUES FOR BISECTION **************************************
!K !
!K       X1   = 1e-3          ! Low value of updraft
!K       CALL ACTIVATE (X1, TP,AKK,ei,A,B,ACCOM,SG,NACT1,SMAX,SIG)
!K       Y1   = NACT1/NACT - 1.
!K !
!K       X2   = 20           ! High value of updraft
!K       CALL ACTIVATE (X2, TP,AKK,ei,A,B,ACCOM,SG,NACT2,SMAX,SIG)
!K       Y2   = NACT2/NACT - 1.
!K !
!K ! *** PERFORM BISECTION ************************************************
!K !
!K 20    DO 30 I=1,MAXIT
!K          X3   = 0.5*(X1+X2)
!K          CALL ACTIVATE (X3, TP,AKK,ei,A,B,ACCOM,SG,NACT3,SMAX,SIG)
!K          Y3   = NACT3/NACT - 1.0
!K !
!K          IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
!K             Y2    = Y3
!K             X2    = X3
!K          ELSE
!K             Y1    = Y3
!K             X1    = X3
!K          ENDIF
!K !
!K          IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
!K          NITER = I
!K !
!K 30    CONTINUE
!K !
!K ! *** CONVERGED ; RETURN ***********************************************
!K !
!K 40    X3    = 0.5*(X1+X2)
!K       CALL ACTIVATE (X3, TP,AKK,ei,A,B,ACCOM,SG,NACT3,SMAX,SIG)
!K       Y3    = NACT3/NACT - 1.
!K       WPARC = X3
!K !
!K       RETURN
!K !
!K ! *** END OF SUBROUTINE WACTIVE ****************************************
!K !
!K       END SUBROUTINE WACTIV
!K !=======================================================================
!K 

!=======================================================================
!
! *** SUBROUTINE ACTIVATE
! *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!     and Nenes (in preparation).
!
! *** WRITTEN BY ATHANASIOS NENES FOR KOHLER PARTICLES
! *** MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH 
!     PARTICLES 
!
!=======================================================================
!
      SUBROUTINE ACTIVATE (WPARC,TP,AKK,ei,A,B,ACCOM,SG,NDRPL,SMAX,SIG)
!      INCLUDE 'parametr.inc'
!      USE mod_parameters

      REAL:: NDRPL,NDROLNR,WPARCEL,A,B,ACCOM,BETA,BET2, &
       TP(nmod),AKK(nmod),ei(nmod),SG(nmod),SIG(NMOD)

!MOD_CHANGES
      REAL:: NDRPLNR

!
! *** Setup common block variables
!
      

      PRESA = PRES/1.013d5                  ! Pressure (Pa)
      DV    = (0.211/PRESA)*(TEMP/273)**1.94
      DV    = DV*1e-4                       ! Water vapor diffusivity in air
      DBIG  = 5.0e-6
      DLOW  = 0.207683*((ACCOM)**(-0.33048))
      DLOW  = DLOW*1e-6
!
! Dv average
!
      COEF  = ((2*PI*AMW/(RGAS*TEMP))**0.5)
!
      DV    = (DV/(DBIG-DLOW))*((DBIG-DLOW)-(2*DV/ACCOM)*COEF* &
              (LOG((DBIG+(2*DV/ACCOM)*COEF)/(DLOW+(2*DV/ACCOM)* &
              COEF))))                      ! Non-continuum effects


      WPARCEL = WPARC
      
      
!
! *** Setup constants
!
      ALFA = GRAV*AMW*DHV/CPAIR/RGAS/TEMP/TEMP - GRAV*AMA/RGAS/TEMP
      BET1 = PRES*AMA/PSAT/AMW + AMW*DHV*DHV/CPAIR/RGAS/TEMP/TEMP
      BET2 = RGAS*TEMP*DENW/PSAT/DV/AMW/4. + &
             DHV*DENW/4./AKA/TEMP*(DHV*AMW/RGAS/TEMP - 1.)
      BETA = 0.5*PI*BET1*DENW/BET2/ALFA/WPARC/DAIR
      CF1  = 0.5*(((1/BET2)/(ALFA*WPARC))**0.5)
      CF2  = AKOH/3.
!
! *** INITIAL VALUES FOR BISECTION *************************************
!     
      X1   = 1.0e-5   ! Min cloud supersaturation -> 0
      CALL SINTEGRAL (X1,NDRPL,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SINTEG1,SINTEG2,SINTEG3,SIG)
      Y1   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X1 - 1.
!     
      X2   = 1.      ! MAX cloud supersaturation = 100%
      CALL SINTEGRAL (X2,NDRPL,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SINTEG1,SINTEG2,SINTEG3,SIG)
      Y2   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X2 - 1.
!
! *** PERFORM BISECTION ************************************************
!
20    DO 30 I=1,MAXIT
         X3   = 0.5*(X1+X2)
!
        CALL SINTEGRAL (X3,NDRPL,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SINTEG1,SINTEG2,SINTEG3,SIG)
         Y3 = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1.
!
         IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
	     Y2    = Y3
	     X2    = X3
	 ELSE
	     Y1    = Y3
	     X1    = X3
         ENDIF
!
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
	    NITER = I
!
30    CONTINUE
!
! *** CONVERGED ; RETURN ***********************************************
!
40    X3   = 0.5*(X1+X2)
!
      CALL SINTEGRAL (X3,NDRPL,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SINTEG1,SINTEG2,SINTEG3,SIG)
      Y3   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1.

      SMAX = X3
      CALL SINTNR(WPARCEL,TP,AKK,ei,A,B,ACCOM,SG,NDRPLNR,X3,Y3,SIG)


 113  FORMAT (G12.5)

      RETURN
!
! *** END OF SUBROUTINE ACTIVATE ***************************************
!
      END SUBROUTINE ACTIVATE
!=======================================================================


      SUBROUTINE SINTNR (WPARCEL,TP,AKK,ei,A,B,ACCOM,SG,NDRPLNR,X3,Y3,SIG)
!      INCLUDE 'parametr.inc'
!      USE mod_parameters
      REAL:: WPARCEL,NDRPLNR,A,B,ACCOM,BETA,BET2, &
       TP(nmod),AKK(nmod),ei(nmod),SG(nmod),SIG(nmod)

!
! *** Setup constants
!

      PRESA = PRES/1.013e5                  ! Pressure (Pa)
      DV    = (0.211/PRESA)*(TEMP/273)**1.94
      DV    = DV*1e-4                       ! Water vapor diffusivity in air
      DBIG  = 5.0e-6
      DLOW  = 0.207683*((ACCOM)**(-0.33048))
      DLOW  = DLOW*1e-6
!
! Dv average
!
      COEF  = ((2*PI*AMW/(RGAS*TEMP))**0.5)
!
      DV    = (DV/(DBIG-DLOW))*((DBIG-DLOW)-(2*DV/ACCOM)*COEF* &
              (LOG((DBIG+(2*DV/ACCOM)*COEF)/(DLOW+(2*DV/ACCOM)* &
              COEF))))                      ! Non-continuum effects

      BET2 = RGAS*TEMP*DENW/PSAT/DV/AMW/4. + &
             DHV*DENW/4./AKA/TEMP*(DHV*AMW/RGAS/TEMP - 1.)
      BETA = 0.5*PI*BET1*DENW/BET2/ALFA/WPARCEL/DAIR
      CF1  = 0.5*(((1/BET2)/(ALFA*WPARCEL))**0.5)
      CF2  = AKOH/3.

      CALL SINTEGRAL (X3,NDRPLNR,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SINTEG1,SINTEG2,SINTEG3,SIG)

      Y3   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1.

      RETURN

      END SUBROUTINE SINTNR
!=======================================================================


!=======================================================================
!
! *** SUBROUTINE SINTEGRAL
! *** THIS SUBROUTINE CALCULATES THE CONDENSATION INTEGRALS, ACCORDING
!     TO THE POPULATION SPLITTING ALGORITHM OF Nenes and Seinfeld (2003)
!     Modal formulation according to Fountoukis and Nenes (2004)
!
! *** WRITTEN BY ATHANASIOS NENES for Kohler Particles
! *** MODFIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH
!     PARTICLES
!=======================================================================
!
      SUBROUTINE SINTEGRAL (SPAR,SUMMA,WPARCEL,TP,AKK,ei,A,B,BET2,SG, &
       SUM,SUMMAT,SUMFHH,SIG)
!PK   FINISH CHANGES
!      INCLUDE 'parametr.inc'
!      USE mod_parameters
      REAL:: SUM, SUMMAT, SUMMA, Nd(nmod),WPARCEL,A,B,BET2, &
                       INTEG1(nmod),INTEG2(nmod),TP(nmod),AKK(nmod), &
                       ei(nmod),SG(nmod),SIG(nmod) &
!PK   CHANGES BEGIN
                       ,SUMFHH,INTEG1F(nmod),NdF(nmod)
!PK   CHANGES FINISH
      REAL::             ERF1,ERF2,ERF3,ERF4,ERF5,ERF6,ERF4F,ERF5F,ERF66F
      REAL::             ORISM1, ORISM2, ORISM3, ORISM4, ORISM5,ORISM6
!
! Ricardo M. Some new variables
!
      REAL::             intaux1p1, intaux1p2, DLGSP1, DLGSP2, scrit

!PK   BEGIN CHANGES
      REAL::             ORISM1F, ORISM2F, ORISM3F, ORISM4F, ORISM5F, &
                       ORISM6F, ORISM7F, ORISM8F, ORISM9F, ORISM10F, &
                       ORISM11F, ORISM66F
!PK   FINISH CHANGES

!PK   DETERMINATION OF EXPONENT
      C1 = (D11)+(D12/A)+(D13/(A*A))+(D14/(A*A*A))+(D15/(A*A*A*A))
      C2 = (D21)+(D22/A)+(D23/(A*A))+(D24/(A*A*A))+(D25/(A*A*A*A))
      C3 = (D31)+(D32/A)+(D33/(A*A))+(D34/(A*A*A))+(D35/(A*A*A*A))
      C4 = (D41)+(D42/A)+(D43/(A*A))+(D44/(A*A*A))+(D45/(A*A*A*A))
!
! *** Here is where the criterion with the descriminant is put. When it
!     is < 0, then set CRIT2 = .TRUE. Otherwise, set the two values of
!     SSPLT and continue.
!
      DESCR  = 1. - (16./9.)*ALFA*WPARCEL*BET2*(AKOH/SPAR**2)**2
      SQTWO  = SQRT(2.0d0)

!C ** Population Splitting -- Modified by Ricardo Morales 2013
!
      IF (DESCR.LE.0.) THEN
         CRIT2  = .TRUE.             ! SSPLT1,SSPLT2 do not exist
         scrit  = ((16d0/9d0)*ALFA*WPARCEL*BET2*(AKOH**2))**(0.25d0)
!         RATIO  = (2.0d7/3.0)*AKOH*SPAR**(-0.3824)
         RATIO  = (2.0d7/3.0)*AKOH*(SPAR**(-0.3824)-scrit**(-0.3824))   ! Computing sp1 and sp2 (sp1 = sp2)
         RATIO  = 1/SQTWO + RATIO
         IF (RATIO.GT.1.0) RATIO = 1.0
         SSPLT2 = SPAR*RATIO
      ELSE
         CRIT2  = .FALSE.
         SSPLT1 = 0.5*(1.-SQRT(DESCR)) ! min root of both
         SSPLT2 = 0.5*(1.+SQRT(DESCR)) ! max root of both
         SSPLT1 = SQRT(SSPLT1)*SPAR       ! Multiply ratios with Smax
         SSPLT2 = SQRT(SSPLT2)*SPAR
      ENDIF
!
      SSPLT = SSPLT2  ! Store Ssplit in COMMON
!
! *** Calculate integrals
!
      SUM       = 0   !Contribution of integral 1 for Kohler 
      SUMMAT    = 0   !Contribution of integral 2 for kohler
      SUMMA     = 0   !Variable that stores all droplets
!PK   BEGIN CHANGES
      SUMFHH    = 0   !Contribution of FHH integral
!PK   FINISH CHANGES
!
!      SQTWO     = SQRT(2.)
!
!CC      DO 999 j = 1, NMD
!

!PK   BEGIN CHANGES

      DO J = 1, NMD

! MOD_CHANGES
      if (SG(J) .gt. 0.) then


!PK   Solve for Kohler
      IF (MODE(J).EQ.1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Modifications introduced by Ricardo Morales start here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Modified by Dan for KBLOCK version

!        DLGSG  = DLOG(SIG(J))                            !ln(sigmai)
!        DLGSP  = DLOG(SG(J)/SPAR)                        !ln(sg/smax)
!        DLGSP2 = DLOG(SG(J)/SSPLT2)                      !ln(sg/sp2)
        DLGSG  = LOG(SIG(J))                            !ln(sigmai)
        DLGSP  = LOG(SG(J)/SPAR)                        !ln(sg/smax)
        DLGSP2 = LOG(SG(J)/SSPLT2)                      !ln(sg/sp2)
!C
        ORISM1 = 2.d0*DLGSP2/(3.d0*SQTWO*DLGSG)          ! u(sp2)
        ORISM2 = ORISM1 - 3.d0*DLGSG/(2.d0*SQTWO)        ! u(sp2)-3ln(sigmai)/(2sqrt(2)
        ORISM5 = 2.d0*DLGSP/(3.d0*SQTWO*DLGSG)           ! u(smax)
        ORISM3 = ORISM5 - 3.d0*DLGSG/(2.d0*SQTWO)        ! u(smax)-3ln(sigmai)/(2sqrt(2)
        DEQ    = AKOH*2d0/SG(j)/3d0/SQRT(3d0)            ! Dp0 = Dpc/sqrt(3) - Equilibrium diameter

        ERF2   = ERF(ORISM2)
        ERF3   = ERF(ORISM3)
 
        INTEG2(J) = (EXP(9D0/8D0*DLGSG*DLGSG)*TP(J)/SG(J))* &
                    (ERF2 - ERF3)                          ! I2(sp2,smax)

        IF (CRIT2) THEN     

          ORISM6 = (SQTWO*DLGSP2/3d0/DLGSG)-(1.5d0*DLGSG/SQTWO)
          ERF6   = ERF(ORISM6)

          INTEG1(J) = 0.0d0
          DW3       = TP(j)*DEQ*EXP(9D0/8D0*DLGSG*DLGSG)* &   ! 'inertially' limited particles
                      (1d0-ERF6)*((BET2*ALFA*WPARCEL)**0.5d0)

        ELSE
 
          EKTH    = EXP(9D0/2d0*DLGSG*DLGSG)
!Modified by Dan for KBLOCK version          
          !DLGSP1  = DLOG(SG(J)/SSPLT1)                      ! ln(sg/sp1)
          DLGSP1  = LOG(SG(J)/SSPLT1)                      ! ln(sg/sp1)
          ORISM4  = ORISM1 + 3.d0*DLGSG/SQTWO               ! u(sp2) + 3ln(sigmai)/sqrt(2)
          ERF1    = ERF(ORISM1)
          ERF4    = ERF(ORISM4)

          intaux1p2 =  TP(J)*SPAR*((1-ERF1) - &
                       0.5d0*((SG(J)/SPAR)**2)*EKTH*(1-ERF4))  ! I1(0,sp2)

          ORISM1  = 2.d0*DLGSP1/(3.d0*SQTWO*DLGSG)          ! u(sp1)
          ORISM4  = ORISM1 + 3.d0*DLGSG/SQTWO               ! u(sp1) + 3ln(sigmai)/sqrt(2)
          ORISM6  = (SQTWO*DLGSP1/3d0/DLGSG)-(1.5d0*DLGSG/SQTWO)

          ERF1 = ERF(ORISM1)
          ERF4 = ERF(ORISM4)
          ERF6 = ERF(ORISM6)

          intaux1p1 = TP(J)*SPAR*((1-ERF1) - &
                      0.5d0*((SG(J)/SPAR)**2)*EKTH*(1-ERF4))    ! I1(0,sp1)

          INTEG1(J) = (intaux1p2-intaux1p1)                   ! I1(sp1,sp2) = I1(0,sp2) - I1(0,sp1)
!
          DW3 = TP(j)*DEQ*EXP(9D0/8D0*DLGSG*DLGSG)* &           ! 'inertially' limited particles.
                (1d0-ERF6)*((BET2*ALFA*WPARCEL)**0.5d0)
 
        ENDIF

!C *** Calculate number of Drops

        ERF5     = ERF(ORISM5)
! 
        Nd(J)    = (TP(J)/2.0)*(1.0-ERF5)
        SUM      = SUM    + INTEG1(J) + DW3           !SUM OF INTEGRAL 1 FOR KOHLER
        SUMMAT   = SUMMAT + INTEG2(J)                 !SUM OF INTEGRAL 2 FOR KOHLER
        SUMMA    = SUMMA  + Nd(J)                     !SUM OF ACTIVATED KOHLER PARTICLES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Modifications Introduced by Ricardo Morales end here.
! I am comenting out the next lines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!      DLGSG     = LOG(SIG(J))      !ln(sigmai)
!      DLGSP     = LOG(SG(J)/SPAR)  !ln(sg/smax)
!      ORISM1    = 2.*LOG(SG(J)/SSPLT2)/(3.d0*SQTWO*DLGSG) !upart
!      ORISM2    = ORISM1 - 3.*DLGSG/(2.*SQTWO) 
!      ORISM3    = 2.*DLGSP/(3.*SQTWO*DLGSG)-3.*DLGSG/(2.*SQTWO)
!      ORISM4    = ORISM1 + 3.*DLGSG/SQTWO
!      ORISM5    = 2.*DLGSP/(3*SQTWO*DLGSG)
!      EKTH      = EXP(9./2.*DLGSG*DLGSG)

! ERF1=ERF(ORISM1)
! ERF2=ERF(ORISM2)
! ERF3=ERF(ORISM3)
! ERF4=ERF(ORISM4)

 !!     CALL CALCERF(ORISM1,ERF1)
 !!     CALL CALCERF(ORISM2,ERF2)
 !!     CALL CALCERF(ORISM3,ERF3)
 !!     CALL CALCERF(ORISM4,ERF4)

!      INTEG1(J) = TP(J)*SPAR*((1-ERF1) - &
!                  0.5*((SG(J)/SPAR)**2)*EKTH*(1-ERF4))
!      INTEG2(J) = (EXP(9./8.*DLGSG*DLGSG)*TP(J)/SG(J))* &
!                  (ERF2 - ERF3)
!
!!
!!      ! Correction for kinetically limited large particles (D. Barahona et al. ACPD 2009)
!!
!      DLGSP     = LOG(SG(j)/SSPLT)
!      ORISM6  = (SQTWO*DLGSP/3./DLGSG)-(1.5*DLGSG/SQTWO)
!      DEQ= AKOH*2./SG(j)/3./SQRT(3.)
!ERF6=ERF(ORISM6)
!!      CALL CALCERF(ORISM6,ERF6)
!      DW3=TP(j)*DEQ*EXP(9./8.*DLGSG*DLGSG)* &
!          (1.-ERF6)*((BET2*ALFA*WPARCEL)**0.5)
!!
!! *** Calculate number of Drops
!!
!!ERF5=ERF(ORISM5)
!!      CALL CALCERF(ORISM5,ERF5)
!      Nd(J)     = (TP(J)/2.0)*(1.0-ERF5)
!!
!      SUM       = SUM    + INTEG1(J) + DW3 !SUM OF INTEGRAL 1 FOR KOHLER
!      SUMMAT    = SUMMAT + INTEG2(J) !SUM OF INTEGRAL 2 FOR KOHLER
!      SUMMA     = SUMMA  + Nd(J)     !SUM OF ACTIVATED KOHLER PARTICLES
      
      
!PK   THIS IS WHERE CODING FINISHES FOR KOHLER

!PK   Solve fpr FHH
      ELSEIF (MODE(J).EQ.2) THEN
!      
!VK calculate the exponent according to the unified FHH theory
      XFHH = (C1) + (C2/B) + (C3/(B*B)) + (C4/(B*B*B))
      XFHH = XFHH*exp(log(-1.5/XFHH)*(1.-ei(J))** &
      (0.1693*exp(-0.988*AKK(J)))) !unified FHH theory --- VAK


      DLGSGF  = LOG(SIG(J))                  !ln(sigma,i)
      DLGSPF  = LOG(SG(J)/SPAR)              !ln(sg/smax)
      ORISM1F = (SG(J)*SG(J))/(SPAR*SPAR)     !(sg/smax)^2
      ORISM2F = EXP(2.*XFHH*XFHH*DLGSGF*DLGSGF) !exp(term)
      ORISM3F = SQTWO*XFHH*DLGSGF             !sqrt(2).x.ln(sigma,i)
      ORISM4F = DLGSPF/(-1*ORISM3F)           !Umax
      ORISM5F = ORISM3F - ORISM4F

ERF5F=ERF(ORISM5F)
!      CALL CALCERF(ORISM5F,ERF5F)

      ORISM6F = ERF5F
      ORISM7F = ORISM6F + 1
      ORISM8F = 0.5*ORISM1F*ORISM2F*ORISM7F

ERF4F=ERF(ORISM4F)
!      CALL CALCERF(ORISM4F,ERF4F)

      ORISM9F = ORISM8F + ERF4F - 1

      INTEG1F(J) =-1*TP(J)*SPAR*ORISM9F
!
!      ! Correction for kinetically limited large particles (D. Barahona et al. ACPD 2009)
!
      DLGSPF     = LOG(SG(j)/SSPLT)
      ORISM66F  = (-DLGSPF/XFHH/SQTWO/DLGSGF)-(-XFHH*DLGSGF/SQTWO)
         DEQF=Dpc(j)

ERF66F=ERF(ORISM66F)

!      CALL CALCERF(ORISM66F,ERF66F)

      DW3F=TP(j)*DEQF*EXP(XFHH**2*DLGSGF*DLGSGF)* &
          (1.-ERF66F)*((BET2*ALFA*WPARCEL)**0.5)   
!
! Ricardo M. - Dust particles are not considered kinetically limited.
      DW3F = 0.0d0

!
! *** Calculate number of drops activated by FHH theory
!

ERF4F=ERF(ORISM4F)
!      CALL CALCERF(ORISM4F,ERF4F)

      NdF(J) = (TP(J)/2.0)*(1-ERF4F)
      SUMFHH = SUMFHH + INTEG1F(J) + DW3F  !Sum of Integral 1 for FHH + GCCN correction
      SUMMA = SUMMA + NdF(J)         !Sum of ACTIVATED Kohler + FHH particles
      ENDIF
!PK   THIS IS WHERE CODING FINISHES FOR FHH

! MOD_CHANGES
      ENDIF	! if (SG(J) .gt. 0,) then

      ENDDO
!      pause
      RETURN
      END SUBROUTINE SINTEGRAL
!=======================================================================


!=======================================================================
!
! *** SUBROUTINE PROPS
! *** THIS SUBROUTINE CALCULATES THE THERMOPHYSICAL PROPERTIES
!
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE PROPS
!      INCLUDE 'parametr.inc'
!      USE mod_parameters

!K       USE mo_constants,	ONLY: amw, rhoh2o, alv, cpd, zsten	!K 
       
!MOD_CHANGES
!      REAL::  VPRES, SFT
!
!      AMW   = 18.e-3                        ! Water molecular weight
!      DENW  = 1e3                           ! Water density
!      DHV   = 2.25e6                        ! Water enthalpy of vaporization
!      CPAIR = 1.0061e3                      ! Air Cp
!CHANGES

      AMW   = 18.0154e-3		    ! Water molecular weight
      DENW  = 1e3			    ! Water density
      DHV   = 2.5008e6  		      ! Water enthalpy of vaporization
      CPAIR = 1.00546e3 		     ! Air Cp

!K use instead from mo_constants (dp)?	!K 
!K       AMW   = amw*e-3	!K     (same as in ham_activ_koehler_ab: zamw=amw*1.E-3_dp)
!K       DENW  = rhoh2o		!K 
!K       DHV   = alv		!K 
!K       CPAIR = cpd		!K 


      DAIR  = PRES*AMA/RGAS/TEMP            ! Air density
!
!      AKA   = (4.39+0.071*TEMP)*1e-3        ! Air thermal conductivity
!CHANGES
!      AKA  = (5.69D0+0.017D0*(TEMP-273.15D0))*1D-3*4.187D0 ! Air thermal conductivity
!
!

!      PSAT  = VPRES(TEMP)*(1e5/1.0e3) ! Saturation vapor pressure
!CHANGES
!      PSAT  = 100.d0*6.1121d0*exp((18.678d0 - (TEMP-273.15d0)/ 234.5d0)* &
!              (TEMP-273.15d0)/(257.14d0 + (TEMP-273.15d0))) ! Saturation vapor pressure
!
!      SURT  = SFT(TEMP)       ! Surface Tension for water (J m-2)
!CHANGES

      SURT  = 7.5E-02       ! Surface Tension for water (J m-2)

!K use instead from mo_constants (dp)?	!K 
!K      SURT  = zsten	!K
!
      RETURN
!
! *** END OF SUBROUTINE PROPS ******************************************
!
      END SUBROUTINE PROPS
!=======================================================================


!K !=======================================================================
!K !
!K ! *** FUNCTION VPRES
!K ! *** THIS FUNCTION CALCULATES SATURATED WATER VAPOUR PRESSURE AS A
!K !     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -50 AND
!K !     50 C.
!K !
!K !========================= ARGUMENTS / USAGE ===========================
!K !
!K !  INPUT:
!K !     [T]
!K !     REAL:: variable.
!K !     Ambient temperature expressed in Kelvin.
!K !
!K !  OUTPUT:
!K !     [VPRES]
!K !     REAL:: variable.
!K !     Saturated vapor pressure expressed in mbar.
!K !
!K !=======================================================================
!K !
!K       REAL FUNCTION VPRES (T)
!K       REAL:: A(0:6), T
!K       DATA A/6.107799610E+0, 4.436518521E-1, 1.428945805E-2, &
!K              2.650648471E-4, 3.031240396E-6, 2.034080948E-8, &
!K              6.136820929E-11/
!K !
!K ! Calculate polynomial (without exponentiation).
!K !
!K       TTEMP = T-273
!K       VPRES = A(6)*TTEMP
!K       DO I=5,1,-1
!K          VPRES = (VPRES + A(I))*TTEMP
!K       ENDDO
!K       VPRES = VPRES + A(0)
!K !
!K ! End of FUNCTION VPRES
!K !
!K       RETURN
!K       END FUNCTION VPRES
!K !=======================================================================
!K 
!K 
!K !=======================================================================
!K !
!K ! *** THIS FUNCTION CALCULATES WATER SURFACE TENSION AS A
!K !     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -40 AND
!K !     40 C.
!K !
!K ! ======================== ARGUMENTS / USAGE ===========================
!K !
!K !  INPUT:
!K !     [T]
!K !     REAL:: variable.
!K !     Ambient temperature expressed in Kelvin.
!K !
!K !  OUTPUT:
!K !     [SFT]
!K !     REAL:: variable.
!K !     Surface Tension expressed in J m-2.
!K !
!K !=======================================================================
!K !
!K       REAL FUNCTION SFT (T)
!K       REAL:: T
!K !
!K       TPARS = T-273
!K       SFT   = 0.0761-1.55e-4*TPARS
!K !
!K       RETURN
!K       END FUNCTION SFT
!K !=======================================================================
!K 
!K 
!K ! ***********************************************************************
!K !
!K       SUBROUTINE GAULEG (X,W,N)
!K !
!K ! Calculation of points and weights for N point GAUSS integration
!K ! ***********************************************************************
!K       DIMENSION X(N), W(N)
!K       PARAMETER (EPS=1.E-6)
!K       PARAMETER (X1=-1.0, X2=1.0)
!K !
!K ! Calculation
!K !
!K       M=(N+1)/2
!K       XM=0.5*(X2+X1)
!K       XL=0.5*(X2-X1)
!K       DO 12 I=1,M
!K         Z=COS(3.141592654*(I-.25)/(N+.5))
!K 1       CONTINUE
!K           P1=1.
!K           P2=0.
!K           DO 11 J=1,N
!K             P3=P2
!K             P2=P1
!K             P1=((2.*J-1.)*Z*P2-(J-1.)*P3)/J
!K 11        CONTINUE
!K           PP=N*(Z*P1-P2)/(Z*Z-1.)
!K           Z1=Z
!K           Z=Z1-P1/PP
!K         IF(ABS(Z-Z1).GT.EPS)GO TO 1
!K         X(I)=XM-XL*Z
!K         X(N+1-I)=XM+XL*Z
!K         W(I)=2.*XL/((1.-Z*Z)*PP*PP)
!K         W(N+1-I)=W(I)
!K 12    CONTINUE
!K       RETURN
!K       END SUBROUTINE GAULEG
!K 
!K 
!K 
!K !=======================================================================
!K !
!K ! *** REAL FUNCTION ERF
!K ! *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION
!K !
!K ! *** OBTAINED FROM NUMERICAL RECIPIES
!K !
!K !=======================================================================
!K       SUBROUTINE CALCERF(X,ERF)
!K       IF(X.LT.0.)THEN
!K         ERF=-GAMMP(.5,X**2)
!K       ELSE
!K         ERF=GAMMP(.5,X**2)
!K       ENDIF
!K       RETURN
!K       END SUBROUTINE CALCERF
!K !
!K !=======================================================================
!K !
!K 
!K       FUNCTION GAMMLN(XX)
!K !
!K !=======================================================================
!K !
!K       REAL:: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
!K       DATA COF,STP/76.18009173,-86.50532033,24.01409822, &
!K           -1.231739516,.120858003e-2,-.536382e-5,2.50662827465/
!K       DATA HALF,ONE,FPF/0.5,1.0,5.5/
!K       X=XX-ONE
!K       TMP=X+FPF
!K       TMP=(X+HALF)*LOG(TMP)-TMP
!K       SER=ONE
!K       DO 11 J=1,6
!K         X=X+ONE
!K         SER=SER+COF(J)/X
!K 11    CONTINUE
!K       GAMMLN=TMP+LOG(STP*SER)
!K       RETURN
!K       END FUNCTION GAMMLN
!K !
!K !=======================================================================
!K !
!K 
!K       FUNCTION GAMMP(A,X)
!K !
!K !=======================================================================
!K !
!K       IF(X.LT.0..OR.A.LE.0.)PAUSE
!K       IF(X.LT.A+1.)THEN
!K 
!K         CALL GSER(GAMSER,A,X,GLN)
!K 
!K         GAMMP=GAMSER
!K       ELSE
!K 
!K         CALL GCF(GAMMCF,A,X,GLN)
!K 
!K         GAMMP=1.-GAMMCF
!K       ENDIF
!K       RETURN
!K       END FUNCTION GAMMP
!K 
!K 
!K !
!K !=======================================================================
!K !
!K 
!K       SUBROUTINE GCF(GAMMCF,A,X,GLN)
!K !
!K !=======================================================================
!K !
!K       PARAMETER (ITMAX=100,EPS=3.E-7)
!K       GLN=GAMMLN(A)
!K       GOLD=0.
!K       A0=1.
!K       A1=X
!K       B0=0.
!K       B1=1.
!K       FAC=1.
!K       DO 11 N=1,ITMAX
!K         AN=FLOAT(N)
!K         ANA=AN-A
!K         A0=(A1+A0*ANA)*FAC
!K         B0=(B1+B0*ANA)*FAC
!K         ANF=AN*FAC
!K         A1=X*A0+ANF*A1
!K         B1=X*B0+ANF*B1
!K         IF(A1.NE.0.)THEN
!K           FAC=1./A1
!K           G=B1*FAC
!K           IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
!K           GOLD=G
!K         ENDIF
!K 11    CONTINUE
!K       PAUSE 'A too large, ITMAX too small'
!K 1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
!K       RETURN
!K       END SUBROUTINE GCF
!K 
!K 
!K !
!K !=======================================================================
!K !
!K       SUBROUTINE GSER(GAMSER,A,X,GLN)
!K !
!K !=======================================================================
!K !
!K       PARAMETER (ITMAX=100,EPS=3.E-7)
!K       GLN=GAMMLN(A)
!K       IF(X.LE.0.)THEN
!K         IF(X.LT.0.)PAUSE
!K         GAMSER=0.
!K         RETURN
!K       ENDIF
!K       AP=A
!K       SUM=1./A
!K       DEL=SUM
!K       DO 11 N=1,ITMAX
!K         AP=AP+1.
!K         DEL=DEL*X/AP
!K         SUM=SUM+DEL
!K         IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
!K 11    CONTINUE
!K       PAUSE 'A too large, ITMAX too small'
!K 1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
!K       RETURN
!K       END SUBROUTINE GSER

END MODULE mo_ham_nenes

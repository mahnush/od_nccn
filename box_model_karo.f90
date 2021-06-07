PROGRAM box_model_karo

  USE mo_kind,  	 	ONLY: dp
  USE mo_constants,		ONLY: g, rd, amd, argas, amw  
  USE mo_box_settings,   	ONLY: kproma, kbdim, klev, krow, ktdia, ktrac
  USE mo_box_input,		ONLY: read_input, get_input_sub, get_input_para
  USE mo_box_netcdf_io, 	ONLY: netcdf_dimensions, netcdf_input, netcdf_output
  USE mo_box_updraft,           ONLY: nw, pw
  USE mo_box_aerosols,		ONLY: nmod, aertype, sizedist
  USE mo_ham_ccnctl,		ONLY: nsat, zsat
  USE mo_ham_ccn,		ONLY: ham_ccn  
  USE mo_ham_activ,		ONLY: ham_activ_koehler_ab, ham_activ_abdulrazzak_ghan
!-------------------------------------------------------------------------------

IMPLICIT NONE

  CHARACTER(len=150)	:: ifilename1, ifilename2
  CHARACTER(len=150)	:: ofilename1
  CHARACTER(len=:), ALLOCATABLE :: ifile1,ifile2,ofile1
   
  CHARACTER(len=4)	:: year
  CHARACTER(len=2)	:: mon
  INTEGER		:: nTimes, nLev, nLat, nLon

  
  ! parameters from input file
  REAL(dp),ALLOCATABLE,DIMENSION(:)		:: time, lev, lat, lon
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)         :: spress,sz
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: press, temp, spechum
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: mixratio_SSs, mixratio_SSm, mixratio_SSl
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: mixratio_DUs, mixratio_DUm, mixratio_DUl
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: mixratio_OMh, mixratio_OMn
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: mixratio_BCh, mixratio_BCn  
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: mixratio_SU
  

  ! parameters for output file  
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_01, VarCCN_02, VarCCN_03
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)       :: VarCCN_04, VarCCN_05, VarCCN_06
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)       :: VarCCN_07, VarCCN_08, VarCCN_09
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)       :: VarCCN_10
!--msh
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: geop
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_02_OMh, VarCCN_02_BCh, VarCCN_02_SU

! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_04_SSs, VarCCN_04_SSm, VarCCN_04_SSl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_04_OMh, VarCCN_04_BCh, VarCCN_04_SU

! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_10_SSs, VarCCN_10_SSm, VarCCN_10_SSl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarCCN_10_OMh, VarCCN_10_BCh, VarCCN_10_SU
! 
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarNC_SSs, VarNC_SSm, VarNC_SSl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarNC_DUs, VarNC_DUm, VarNC_DUl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarNC_OMh, VarNC_OMn
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarNC_BCh, VarNC_BCn  
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarNC_SU, VarNC
!   
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarMC_SSs, VarMC_SSm, VarMC_SSl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarMC_DUs, VarMC_DUm, VarMC_DUl
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarMC_OMn, VarMC_OMh
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarMC_BCn, VarMC_BCh  
! REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)	:: VarMC_SU, VarMC
  
   
  ! Local parameters
  INTEGER				:: ig, it, la, lo, jclass
  INTEGER				:: missflag, nt, ncd_activ=2
  !msh
  INTEGER                               :: vr
  REAL(dp)                              :: pratio(kbdim,klev), tavg(kbdim,klev)
!  REAL(dp)                              :: a(kbdim,klev), b(kbdim,klev)

  REAL(dp):: 	ptm1(kbdim,klev), 	&
		papm1(kbdim,klev), 	&
		pqm1(kbdim,klev), 	&		
		pesw(kbdim,klev), 	&		
		prho(kbdim,klev),	&
		prho_dry(kbdim,klev)
  
  REAL(dp)::    psp1(kbdim),            &
                pz1(kbdim)
       
  REAL(dp):: 	mmr_tot(kbdim,klev,nmod),	&
  		mconc(kbdim,klev,nmod),		&
		nconc(kbdim,klev,nmod),		&
		pxtm1(kbdim,klev,nmod)

  REAL(dp):: 	ract(kbdim,klev,nsat,nmod),	&
		cnfrac(kbdim,klev,nsat,nmod),	&
		ccnsum(kbdim,klev,nsat)

  REAL(dp)::	mmr(kbdim,klev,nmod), &
  		za(kbdim,klev,nmod),  &
		zb(kbdim,klev,nmod)
!msh
  REAL(dp)::     znact(kbdim,klev), &
                 zsmax(kbdim,klev,nw), &
                 zrc(kbdim,klev,nmod,nw)

  REAL(dp)::     pcdncact(kbdim,klev,nw)
                  

  REAL(dp):: ccntest
       
  REAL(dp), PARAMETER :: rfdeltat = 1000.0_dp
  REAL(dp), PARAMETER :: tlbound =  50.0_dp	    		! lower bound [K]
  REAL(dp), PARAMETER :: tubound = 400.0_dp	    		! upper bound [K]
  INTEGER,  PARAMETER :: jptlucu1 = NINT(rfdeltat*tlbound) 	! lookup table lower bound
  INTEGER,  PARAMETER :: jptlucu2 = NINT(rfdeltat*tubound) 	! lookup table upper bound
  REAL(dp)  :: tlucuaw(jptlucu1-1:jptlucu2+1)    		! table - Es*Rd/Rv, water phase only    

        
!-------------------------------------------------------------------------------
!--- read input statements:

  WRITE(*,*) ' '    
  WRITE(*,*) 'BOXFILE SETTINGS:'
  
!  CALL read_input(ifilename1, ifilename2, ifilename3, ifilename4, &
!  			ofilename1, ofilename2, ofilename3, year, mon)
  

  CALL read_input(ifilename1, ifilename2, &
  			ofilename1, year, mon)
   
   
  ALLOCATE(character(len=LEN_TRIM(ifilename1)) :: ifile1)
  ALLOCATE(character(len=LEN_TRIM(ifilename2)) :: ifile2)  
  !ALLOCATE(character(len=LEN_TRIM(ifilename3)) :: ifile3)   
  !ALLOCATE(character(len=LEN_TRIM(ifilename4)) :: ifile4)   
  ALLOCATE(character(len=LEN_TRIM(ofilename1)) :: ofile1)   

  ifile1 = TRIM(ifilename1)
  ifile2 = TRIM(ifilename2)  
  !ifile3 = TRIM(ifilename3)  
  !ifile4 = TRIM(ifilename4)  
  ofile1 = TRIM(ofilename1)  
   
 
!  WRITE(*,*) 'ifile1 ',  ifilename1
!  WRITE(*,*) 'ifile1 ',  ifile1  
!  WRITE(*,*) 'ifile2 ',  ifilename2    
!  WRITE(*,*) 'ifile2 ',  ifile2      
!  WRITE(*,*) 'ifile3 ',  ifilename3	   
!  WRITE(*,*) 'ifile4 ',  ifilename4	   
!  WRITE(*,*) 'ofile1 ', ofilename1  
!  WRITE(*,*) 'ofile2 ', ofilename2  
!  WRITE(*,*) 'ofile3 ', ofilename3  
!  WRITE(*,*) 'year ',  year
!  WRITE(*,*) 'mon ',  mon  
            
!-------------------------------------------------------------------------------
!--- get netCDF data dimensions from reanalysis:

  WRITE(*,*) ' '
  WRITE(*,*) 'NETCDF INPUT DIMENSIONS:' 

  CALL netcdf_dimensions(ifile2, nTimes, nLev, nLat, nLon)
  
!  WRITE(*,*) 'nTimes',nTimes
!  WRITE(*,*) 'nLev',nLev
!  WRITE(*,*) 'nLat',nLat
!  WRITE(*,*) 'nLon',nLon


!-------------------------------------------------------------------------------
!--- get netCDF data parameters from reanalysis:  
  WRITE(*,*) ' '     
  WRITE(*,*) 'READING INPUT DATA'

  ALLOCATE(time(nTimes))
  ALLOCATE(lev(nLev))
  ALLOCATE(lat(nLat))
  ALLOCATE(lon(nLon))
  ALLOCATE(spress(nLon, nLat, nTimes))
  ALLOCATE(sz(nLon, nLat, nTimes))
  ALLOCATE(press(nLon, nLat, nLev, nTimes))
  ALLOCATE(temp(nLon, nLat, nLev, nTimes))  
  ALLOCATE(spechum(nLon, nLat, nLev, nTimes))    
  ALLOCATE(mixratio_SSs(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_SSm(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_SSl(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_DUs(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_DUm(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_DUl(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_OMn(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_OMh(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_BCn(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_BCh(nLon, nLat, nLev, nTimes))
  ALLOCATE(mixratio_SU(nLon, nLat, nLev, nTimes))
  !--msh
  ALLOCATE(geop(nLon, nLat, nLev, nTimes))
  

  CALL netcdf_input(ifile1, ifile2, &
  		    nTimes, nLev, nLat, nLon, time, lev, lat, lon, &
  		    spress,sz,press, temp, spechum, mixratio_SSs, mixratio_SSm, mixratio_SSl, &
  		    mixratio_DUs, mixratio_DUm, mixratio_DUl, &
  		    mixratio_OMn, mixratio_OMh, mixratio_BCn, mixratio_BCh, mixratio_SU)


!-------------------------------------------------------------------------------
!--- Allocate most of the Output:

  WRITE(*,*) ' '     
  WRITE(*,*) 'START COMPUTATIONS IN LOOPS'

  ALLOCATE(VarCCN_01(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_02(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_03(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_04(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_05(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_06(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_07(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_08(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_09(nLon, nLat, nLev, nTimes))
  ALLOCATE(VarCCN_10(nLon, nLat, nLev, nTimes))
! ALLOCATE(VarCCN_02(nLon, nLat, nLev, nTimes))
! ALLOCATE(VarCCN_02_SSs(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_02_SSm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarCCN_02_SSl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarCCN_02_OMh(nLon, nLat, nLev, nTimes))  
! ALLOCATE(VarCCN_02_BCh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_02_SU(nLon, nLat, nLev, nTimes))     

! ALLOCATE(VarCCN_04(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_04_SSs(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_04_SSm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarCCN_04_SSl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarCCN_04_OMh(nLon, nLat, nLev, nTimes))  
! ALLOCATE(VarCCN_04_BCh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_04_SU(nLon, nLat, nLev, nTimes))     

! ALLOCATE(VarCCN_10(nLon, nLat, nLev, nTimes))
! ALLOCATE(VarCCN_10_SSs(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_10_SSm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarCCN_10_SSl(nLon, nLat, nLev, nTimes))      
! ALLOCATE(VarCCN_10_OMh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_10_BCh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarCCN_10_SU(nLon, nLat, nLev, nTimes))     
    
! ALLOCATE(VarMC(nLon, nLat, nLev, nTimes))
! ALLOCATE(VarMC_SSs(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarMC_SSm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_SSl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_DUs(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_DUm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_DUl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_OMn(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_OMh(nLon, nLat, nLev, nTimes))  
! ALLOCATE(VarMC_BCn(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarMC_BCh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarMC_SU(nLon, nLat, nLev, nTimes))     

! ALLOCATE(VarNC(nLon, nLat, nLev, nTimes))  
! ALLOCATE(VarNC_SSs(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarNC_SSm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_SSl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_DUs(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_DUm(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_DUl(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_OMn(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_OMh(nLon, nLat, nLev, nTimes))  
! ALLOCATE(VarNC_BCn(nLon, nLat, nLev, nTimes))    
! ALLOCATE(VarNC_BCh(nLon, nLat, nLev, nTimes))   
! ALLOCATE(VarNC_SU(nLon, nLat, nLev, nTimes))     


  
!--- Loop over Cloud_groups & Time/FOV: 

 nt = 1

 DO lo = 1,nLon 
  DO la = 1,nLat
   !WRITE(*,*) 'lon',lo, 'lat',la  
   DO it = 1,nTimes
    DO ig = 1,nLev
    WRITE(*,*) ' '
    WRITE(*,*) 'nt: ',nt
    WRITE(*,*) 'lon',lo, 'lat',la, 'lev',ig, 'time',it
 
!-------------------------------------------------------------------------------
!--- Do corrections on input data and get additional parameters:

      !RITE(*,*) ' '
      !RITE(*,*) 'INPUT PARA:'	   

      mmr(kbdim,klev,1) = mixratio_SSs(lo,la,ig,it)	! mass mixing ratios [kg/kg]
      mmr(kbdim,klev,2) = mixratio_SSm(lo,la,ig,it)      
      mmr(kbdim,klev,3) = mixratio_SSl(lo,la,ig,it)      
      mmr(kbdim,klev,4) = mixratio_DUs(lo,la,ig,it) 
      mmr(kbdim,klev,5) = mixratio_DUm(lo,la,ig,it) 
      mmr(kbdim,klev,6) = mixratio_DUl(lo,la,ig,it) 
      mmr(kbdim,klev,7) = mixratio_OMn(lo,la,ig,it) 
      mmr(kbdim,klev,8) = mixratio_OMh(lo,la,ig,it) 
      mmr(kbdim,klev,9) = mixratio_BCn(lo,la,ig,it) 
      mmr(kbdim,klev,10) = mixratio_BCh(lo,la,ig,it) 
      mmr(kbdim,klev,11) = mixratio_SU(lo,la,ig,it) 
      
      
      WHERE(mmr(kbdim,klev,:) .LT. 0.0_dp) mmr(kbdim,klev,:) = 0.0_dp      
      
      ptm1(kbdim,klev) = temp(lo,la,ig,it) 		! temperature [K]
      papm1(kbdim,klev) = press(lo,la,ig,it)		! pressure [Pa]
      pqm1(kbdim,klev) = spechum(lo,la,ig,it)		! specific humidity [kg/kg]
      IF(ptm1(kbdim,klev) .LT. tlbound) THEN
        ptm1(kbdim,klev) = 280.0_dp
        papm1(kbdim,klev) = 100000.0_dp      
        pqm1(kbdim,klev) = 0.016_dp
	mmr(kbdim,klev,:) = 0.0_dp
        missflag = 1
      ELSE
        missflag = 0     
      ENDIF       

      IF(it .EQ. 1 .AND. ig .EQ. 1 .AND. la .EQ. 1 .AND. lo .EQ. 1) THEN  
        CALL get_input_sub(jptlucu1,jptlucu2,tlucuaw)
      ENDIF 
      !----------------------------------geopotential------------------------------------
      vr=nlev-(ig-1)
      IF(vr .EQ. nlev) THEN
         geop(lo,la,vr,it) = sz(lo,la,it)/g
         
      ELSE
         tavg(kbdim,klev) = (temp(lo,la,vr,it)+temp(lo,la,vr+1,it))/2.0
         pratio(kbdim,klev)=press(lo,la,vr,it)/press(lo,la,vr+1,it)
        ! a(kbdim,klev)=(rd*tavg(kbdim,klev))/g
        ! b(kbdim,klev)=log(pratio(kbdim,klev))
         geop(lo,la,vr,it) = geop(lo,la,vr+1,it)-(((rd*tavg(kbdim,klev))/g)*(log(pratio(kbdim,klev))))
      ENDIF
      !--------------------------------------------------------------------------------
      CALL get_input_para(tlucuaw,jptlucu1,jptlucu2,ptm1,papm1,pqm1,pesw,prho,prho_dry)

     ! DO jclass = 1,nmod
     !  WRITE(*,*)  'mmr [kg/kg]', mmr(:,:,jclass)
     ! ENDDO
              
!      WRITE(*,*) 'ptm1 [K]', ptm1
!      WRITE(*,*) 'papm1 [Pa]', papm1   
!      WRITE(*,*) 'pesw [Pa]', pesw     
!      WRITE(*,*) 'prho [kg/m3]', prho
!      WRITE(*,*) 'prho_dry [kg/m3]', prho_dry		       
!      WRITE(*,*) 'pqm1 [kg/kg]', pqm1  

                
!-------------------------------------------------------------------------------
!--- get aerosol mass and number concentration:             

      !WRITE(*,*) ' '
      !WRITE(*,*) 'AEROSOL MASS/NUMBER:' 		      
              
      CALL sizedist(nmod, mmr, prho_dry, mmr_tot, mconc, nconc, pxtm1)

      !DO jclass = 1,nmod      
      ! WRITE(*,*) 'jclass', jclass, 'aertype  ', aertype(jclass)
      ! WRITE(*,*) 'mmr_tot [kg/kg]', mmr_tot(:,:,jclass), 'mconc [kg/m3]', mconc(:,:,jclass)
      ! WRITE(*,*) 'nconc [1/m3]', nconc(:,:,jclass), 'pxtm1 [1/kg]', pxtm1(:,:,jclass)
      !ENDDO
            
!      WRITE(*,*) ' '
!      WRITE(*,*) 'sum mconc [kg/m3]', SUM(mconc(:,:,:)), ' sum nconc [1/m3]', SUM(nconc(:,:,:))
      

!-------------------------------------------------------------------------------
!--- Koehler A & B coefficients:

      !WRITE(*,*) ' '
      !WRITE(*,*) 'KOEHLER COEFFICIENTS:'      
      
      CALL ham_activ_koehler_ab(kproma,   kbdim,   klev,  krow,  ktdia, &
                                pxtm1,    ptm1,    za,    zb            )

!      DO jclass = 1,nmod 
!       WRITE(*,*) 'jclass', jclass, 'za', za(:,:,jclass), 'zb', zb(:,:,jclass)
!      ENDDO

!-------------------------------------------------------------------------------
!--- CCn at fixed prescribed supersaturations:	
      						     	     
      !WRITE(*,*) ' '
      !WRITE(*,*) 'CCN:' 

      CALL ham_ccn(kproma, kbdim, klev, krow,  ktrac,   &	
  		   nconc, za, zb, ract, cnfrac, ccnsum	)  
     
!      DO jsat = 1,nsat 
!       WRITE(*,*) 'max. ssat [% 0-1]', zsat(jsat), 'ccn [1/m3]', ccnsum(:,:,jsat)
!      ENDDO

!      WRITE(*,*) 'ccn_02[1/m3]', ccnsum(:,:,8), ' ccn_04[1/m3]', ccnsum(:,:,17), ' ccn_10[1/m3]', ccnsum(:,:,28)
      
      ccntest = cnfrac(kbdim,klev,17,11)*nconc(kbdim,klev,11)
!      WRITE(*,*) 'ccn_04_SU[1/m3]', ccntest, ' cnfrac_04_SU[%]', cnfrac(:,:,17,11)
!msh======================================================================================
      CALL ham_activ_abdulrazzak_ghan(kproma,    kbdim,    klev,  krow, ktdia,                  &
                                      ncd_activ, pcdncact, pesw,                                &
                                      nconc,     ptm1,     papm1, pqm1,                         &
                                      za,    zb,   znact, zrc, zsmax                            )
      

      
      
!-------------------------------------------------------------------------------
!--- Collect output :	

      !WRITE(*,*) ' '
      !WRITE(*,*) 'OUTPUT VARIABLE:' 

      IF(missflag .EQ. 1) THEN
         
        VarCCN_01(lo,la,ig,it) = -9999.99_dp
        VarCCN_02(lo,la,ig,it) = -9999.99_dp
        VarCCN_03(lo,la,ig,it) = -9999.99_dp
        VarCCN_04(lo,la,ig,it) = -9999.99_dp
        VarCCN_05(lo,la,ig,it) = -9999.99_dp
        VarCCN_06(lo,la,ig,it) = -9999.99_dp
        VarCCN_07(lo,la,ig,it) = -9999.99_dp
        VarCCN_08(lo,la,ig,it) = -9999.99_dp
        VarCCN_09(lo,la,ig,it) = -9999.99_dp
        VarCCN_10(lo,la,ig,it) = -9999.99_dp
!   VarNC(lo,la,ig,it) = -9999.99_dp
!   VarMC(lo,la,ig,it) = -9999.99_dp

!   VarCCN_02_SSs(lo,la,ig,it) =  -9999.99_dp
!   VarCCN_04_SSs(lo,la,ig,it) =  -9999.99_dp        
!   VarCCN_10_SSs(lo,la,ig,it) =  -9999.99_dp    
!   VarMC_SSs(lo,la,ig,it) = -9999.99_dp
!   VarNC_SSs(lo,la,ig,it) = -9999.99_dp
!   
!   VarCCN_02_SSm(lo,la,ig,it) = -9999.99_dp   
!   VarCCN_04_SSm(lo,la,ig,it) = -9999.99_dp	  
!   VarCCN_10_SSm(lo,la,ig,it) = -9999.99_dp	 
!   VarMC_SSm(lo,la,ig,it) = -9999.99_dp
!   VarNC_SSm(lo,la,ig,it) = -9999.99_dp
!   
!   VarCCN_02_SSl(lo,la,ig,it) = -9999.99_dp 
!   VarCCN_04_SSl(lo,la,ig,it) = -9999.99_dp            
!   VarCCN_10_SSl(lo,la,ig,it) = -9999.99_dp      
!   VarMC_SSl(lo,la,ig,it) = -9999.99_dp
!   VarNC_SSl(lo,la,ig,it) = -9999.99_dp
!   
!   VarMC_DUs(lo,la,ig,it) = -9999.99_dp
!   VarNC_DUs(lo,la,ig,it) = -9999.99_dp
!     
!   VarMC_DUm(lo,la,ig,it) = -9999.99_dp
!   VarNC_DUm(lo,la,ig,it) = -9999.99_dp
!   
!   VarMC_DUl(lo,la,ig,it) = -9999.99_dp
!   VarNC_DUl(lo,la,ig,it) = -9999.99_dp
!   
!   VarMC_OMn(lo,la,ig,it) = -9999.99_dp
!   VarNC_OMn(lo,la,ig,it) = -9999.99_dp

!   VarCCN_02_OMh(lo,la,ig,it) = -9999.99_dp 
!   VarCCN_04_OMh(lo,la,ig,it) = -9999.99_dp           
!   VarCCN_10_OMh(lo,la,ig,it) = -9999.99_dp      
!   VarMC_OMh(lo,la,ig,it) = -9999.99_dp
!   VarNC_OMh(lo,la,ig,it) = -9999.99_dp
!    
!   VarMC_BCn(lo,la,ig,it) = -9999.99_dp
!   VarNC_BCn(lo,la,ig,it) = -9999.99_dp

!   VarCCN_02_BCh(lo,la,ig,it) = -9999.99_dp   
!   VarCCN_04_BCh(lo,la,ig,it) = -9999.99_dp	 	 
!   VarCCN_10_BCh(lo,la,ig,it) = -9999.99_dp	 
!   VarMC_BCh(lo,la,ig,it) = -9999.99_dp
!   VarNC_BCh(lo,la,ig,it) = -9999.99_dp

!   VarCCN_02_SU(lo,la,ig,it) = -9999.99_dp
!   VarCCN_04_SU(lo,la,ig,it) = -9999.99_dp           
!   VarCCN_10_SU(lo,la,ig,it) = -9999.99_dp       
!   VarMC_SU(lo,la,ig,it) = -9999.99_dp
!   VarNC_SU(lo,la,ig,it) = -9999.99_dp
    
   ELSE    
     VarCCN_01(lo,la,ig,it) = pcdncact(kbdim,klev,1) 
     VarCCN_02(lo,la,ig,it) = pcdncact(kbdim,klev,2)     
     VarCCN_03(lo,la,ig,it) = pcdncact(kbdim,klev,3)
     VarCCN_04(lo,la,ig,it) = pcdncact(kbdim,klev,4)
     VarCCN_05(lo,la,ig,it) = pcdncact(kbdim,klev,5)
     VarCCN_06(lo,la,ig,it) = pcdncact(kbdim,klev,6)
     VarCCN_07(lo,la,ig,it) = pcdncact(kbdim,klev,7)
     VarCCN_08(lo,la,ig,it) = pcdncact(kbdim,klev,8)
     VarCCN_09(lo,la,ig,it) = pcdncact(kbdim,klev,9)
     VarCCN_10(lo,la,ig,it) = pcdncact(kbdim,klev,10)
     
!   VarMC(lo,la,ig,it) = SUM(mconc(kbdim,klev,:))
!   VarNC(lo,la,ig,it) = SUM(nconc(kbdim,klev,:))			
!   
!   VarCCN_02_SSs(lo,la,ig,it) = cnfrac(kbdim,klev,8,1)*nconc(kbdim,klev,1)    
!   VarCCN_04_SSs(lo,la,ig,it) = cnfrac(kbdim,klev,17,1)*nconc(kbdim,klev,1)	     
!   VarCCN_10_SSs(lo,la,ig,it) = cnfrac(kbdim,klev,28,1)*nconc(kbdim,klev,1)	    
!   VarMC_SSs(lo,la,ig,it) = mconc(kbdim,klev,1)
!   VarNC_SSs(lo,la,ig,it) = nconc(kbdim,klev,1)
!   
!   VarCCN_02_SSm(lo,la,ig,it) = cnfrac(kbdim,klev,8,2)*nconc(kbdim,klev,2)    
!   VarCCN_04_SSm(lo,la,ig,it) = cnfrac(kbdim,klev,17,2)*nconc(kbdim,klev,2)	    	    
!   VarCCN_10_SSm(lo,la,ig,it) = cnfrac(kbdim,klev,28,2)*nconc(kbdim,klev,2)	    
!   VarMC_SSm(lo,la,ig,it) = mconc(kbdim,klev,2)
!   VarNC_SSm(lo,la,ig,it) = nconc(kbdim,klev,2)
!   
!   VarCCN_02_SSl(lo,la,ig,it) = cnfrac(kbdim,klev,8,3)*nconc(kbdim,klev,3)    
!   VarCCN_04_SSl(lo,la,ig,it) = cnfrac(kbdim,klev,17,3)*nconc(kbdim,klev,3)	    	    
!   VarCCN_10_SSl(lo,la,ig,it) = cnfrac(kbdim,klev,28,3)*nconc(kbdim,klev,3)	    
!   VarMC_SSl(lo,la,ig,it) = mconc(kbdim,klev,3)
!   VarNC_SSl(lo,la,ig,it) = nconc(kbdim,klev,3)
!   
!   VarMC_DUs(lo,la,ig,it) = mconc(kbdim,klev,4)
!   VarNC_DUs(lo,la,ig,it) = nconc(kbdim,klev,4)
!     
!   VarMC_DUm(lo,la,ig,it) = mconc(kbdim,klev,5)
!   VarNC_DUm(lo,la,ig,it) = nconc(kbdim,klev,5)
!   
!   VarMC_DUl(lo,la,ig,it) = mconc(kbdim,klev,6)
!   VarNC_DUl(lo,la,ig,it) = nconc(kbdim,klev,6)
!   
!   VarMC_OMn(lo,la,ig,it) = mconc(kbdim,klev,7)
!   VarNC_OMn(lo,la,ig,it) = nconc(kbdim,klev,7)

!   VarCCN_02_OMh(lo,la,ig,it) = cnfrac(kbdim,klev,8,8)*nconc(kbdim,klev,8)    
!   VarCCN_04_OMh(lo,la,ig,it) = cnfrac(kbdim,klev,17,8)*nconc(kbdim,klev,8)	    	    
!   VarCCN_10_OMh(lo,la,ig,it) = cnfrac(kbdim,klev,28,8)*nconc(kbdim,klev,8)	    
!   VarMC_OMh(lo,la,ig,it) = mconc(kbdim,klev,8)
!   VarNC_OMh(lo,la,ig,it) = nconc(kbdim,klev,8)
!    
!   VarMC_BCn(lo,la,ig,it) = mconc(kbdim,klev,9)
!   VarNC_BCn(lo,la,ig,it) = nconc(kbdim,klev,9)

!   VarCCN_02_BCh(lo,la,ig,it) = cnfrac(kbdim,klev,8,10)*nconc(kbdim,klev,10)	 
!   VarCCN_04_BCh(lo,la,ig,it) = cnfrac(kbdim,klev,17,10)*nconc(kbdim,klev,10)                
!   VarCCN_10_BCh(lo,la,ig,it) = cnfrac(kbdim,klev,28,10)*nconc(kbdim,klev,10)        
!   VarMC_BCh(lo,la,ig,it) = mconc(kbdim,klev,10)
!   VarNC_BCh(lo,la,ig,it) = nconc(kbdim,klev,10)

!   VarCCN_02_SU(lo,la,ig,it) = cnfrac(kbdim,klev,8,11)*nconc(kbdim,klev,11)	
!   VarCCN_04_SU(lo,la,ig,it) = cnfrac(kbdim,klev,17,11)*nconc(kbdim,klev,11)	     	     
!   VarCCN_10_SU(lo,la,ig,it) = cnfrac(kbdim,klev,28,11)*nconc(kbdim,klev,11)	     
!   VarMC_SU(lo,la,ig,it) = mconc(kbdim,klev,11)
!   VarNC_SU(lo,la,ig,it) = nconc(kbdim,klev,11)
  
   ENDIF


!-------------------------------------------------------------------------------
!--- end loop:

!      nt = nt+1
    ENDDO !ig
    nt = nt+1
   ENDDO !it
  ENDDO !la
 ENDDO !lo
  
 
!-------------------------------------------------------------------------------
!--- Write output file:	
 WRITE(*,*) ' '			
 WRITE(*,*) 'LOOPS FINISHED - GO TO OUTPUT'
! CALL netcdf_output(   ofilename1, ofilename2, ofilename3, year, mon, &
! 			nTimes, nLev, nLat, nLon, &
!			time, lev, lat, lon, &
!			VarMC, VarNC, VarCCN_02, VarCCN_04, VarCCN_10, &
!			VarMC_SSs, VarMC_SSm, VarMC_SSl, &
!			VarMC_DUs, VarMC_DUm, VarMC_DUl, &
!			VarMC_OMn, VarMC_OMh, VarMC_BCn, VarMC_BCh, VarMC_SU, &
!			VarNC_SSs, VarNC_SSm, VarNC_SSl, &
!			VarNC_DUs, VarNC_DUm, VarNC_DUl, &
!			VarNC_OMn, VarNC_OMh, VarNC_BCn, VarNC_BCh, VarNC_SU, &
!			VarCCN_02_SSs, VarCCN_02_SSm, VarCCN_02_SSl, &
!			VarCCN_02_OMh, VarCCN_02_BCh, VarCCN_02_SU, &
!			VarCCN_04_SSs, VarCCN_04_SSm, VarCCN_04_SSl, &
!			VarCCN_04_OMh, VarCCN_04_BCh, VarCCN_04_SU, &
!			VarCCN_10_SSs, VarCCN_10_SSm, VarCCN_10_SSl, &
!			VarCCN_10_OMh, VarCCN_10_BCh, VarCCN_10_SU)
  
 CALL netcdf_output(   ofile1, year, mon, &
 			nTimes, nLev, nLat, nLon, &
			time, lev, lat, lon, &
                        VarCCN_01, VarCCN_02, VarCCN_03, &
                        VarCCN_04, VarCCN_05, VarCCN_06, &
                        VarCCN_07, VarCCN_08, VarCCN_09, &
                        VarCCN_10,geop )
 


END PROGRAM box_model_karo






















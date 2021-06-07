MODULE mo_box_netcdf_io

  USE netcdf

  IMPLICIT NONE


  !--- Subroutines:

CONTAINS

 !--- get input dimensions from MACC-II reanalysis
  SUBROUTINE netcdf_dimensions(ifile, nT, nLe, nLa, nLo)
   
  IMPLICIT NONE
   
    ! Input parameters
    CHARACTER(len=*), INTENT(in)		:: ifile

    ! Output parameters
    INTEGER,INTENT(out)				:: nT, nLe, nLa, nLo    


    ! Local parameters   
    INTEGER					:: status, ncID, timeID, levID, latID, lonID
    LOGICAL					:: fileexists


    INQUIRE(file=ifile,exist=fileexists)
    IF (fileexists) then
 

      status = nf90_OPEN(ifile,0,ncID)  			!OPEN netCDF file
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'file ID', ncID    

      status = nf90_inq_dimid(ncID,'time',timeID) 
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'time ID', timeID	  

      status = nf90_inquire_DIMENSION(ncID, timeID, len = nT)	!get number of Time
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'ntimes', nT	 

      status = nf90_inq_dimid(ncID,'lev',levID) 
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'lev ID', levID	  

      status = nf90_inquire_DIMENSION(ncID, levID, len = nLe)	!get number of Lev
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'nLev', nLe 	 

      status = nf90_inq_dimid(ncID,'lat',latID) 
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'lat ID', latID	  

      status = nf90_inquire_DIMENSION(ncID, latID, len = nLa)	!get number of Lat   
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'nLat', nLa	 

      status = nf90_inq_dimid(ncID,'lon',lonID) 
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'lon ID', lonID	  

      status = nf90_inquire_DIMENSION(ncID, lonID, len = nLo)	!get number of Lon
       IF(status /= nf90_NoErr) CALL handle_err(status)
       WRITE(*,*) 'nLon', nLo	 


      !--- close netCDF file
      status = nf90_close(ncID)
       IF(status /= nf90_NoErr) CALL handle_err(status)  

    ELSE
      STOP "Stopped: input file does not exist"
    ENDIF
   
  
  END SUBROUTINE netcdf_dimensions




  SUBROUTINE netcdf_input(ifile1, ifile2, &
  			nT, nLe, nLa, nLo, VarTime, VarLev, VarLat, VarLon, &
 		    	VarSP, VarZ, VarP, VarT, VarQ, VarMix_SSs, VarMix_SSm, VarMix_SSl, &
		    	VarMix_DUs, VarMix_DUm, VarMix_DUl, &
		    	VarMix_OMn, VarMix_OMh, VarMix_BCn, VarMix_BCh, VarMix_SU)

  USE mo_kind,          ONLY: dp

  IMPLICIT NONE
   
    ! Input parameters
    CHARACTER(len=*), INTENT(in)	:: ifile1, ifile2
    INTEGER,INTENT(in)			:: nT, nLe, nLa, nLo

    ! Output parameters
    REAL(dp), INTENT(out)		:: VarTime(nT)
    REAL(dp), INTENT(out)		:: VarLev(nLe)    
    REAL(dp), INTENT(out)		:: VarLat(nLa)    
    REAL(dp), INTENT(out)		:: VarLon(nLo)    
    REAL(dp), INTENT(out)		:: VarP(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarT(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarQ(nLo,nLa,nLe,nT)             
    REAL(dp), INTENT(out)		:: VarMix_SSs(nLo,nLa,nLe,nT) 
    REAL(dp), INTENT(out)		:: VarMix_SSm(nLo,nLa,nLe,nT)     
    REAL(dp), INTENT(out)		:: VarMix_SSl(nLo,nLa,nLe,nT)     
    REAL(dp), INTENT(out)		:: VarMix_DUs(nLo,nLa,nLe,nT)     
    REAL(dp), INTENT(out)		:: VarMix_DUm(nLo,nLa,nLe,nT)     
    REAL(dp), INTENT(out)		:: VarMix_DUl(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarMix_OMn(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarMix_OMh(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarMix_BCn(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarMix_BCh(nLo,nLa,nLe,nT)  
    REAL(dp), INTENT(out)		:: VarMix_SU(nLo,nLa,nLe,nT)
    REAL(dp), INTENT(out)               :: VarSP(nLo,nLa,nT)
    REAL(dp), INTENT(out)               :: VarZ(nLo,nLa,nT)
    ! Local parameters 
    REAL(dp)				:: VarHyam(nLe), VarHybm(nLe)
   ! REAL(dp)				:: VarSP(nLo,nLa,nT)
          
    INTEGER				:: lo, la, it
    INTEGER				:: status, ncID1, ncID2, ncID3, ncID4
    INTEGER				:: timeID, levID, latID, lonID, varID
    INTEGER				:: qID, tID, spID, zID, hyamID, hybmID
    INTEGER				:: varID1, varID2, varID3, varID4, varID5
    INTEGER				:: varID6, varID7, varID8, varID9,varID10, varID11 
    CHARACTER(len=7)			:: Paramname(11)



    ! specify aerosol parameters to be read in
    Paramname =  (/'aermr01', 'aermr02', 'aermr03', 'aermr04', 'aermr05', 'aermr06', &
 		   'aermr07', 'aermr08', 'aermr09', 'aermr10', 'aermr11'/) 	     

    WRITE(*,*) 'ifile1 ',  ifile1
    WRITE(*,*) 'ifile2 ',  ifile2
    WRITE(*,*) 'nt'  , nT
    WRITE(*,*) 'nLe' , nLe
    WRITE(*,*) 'nLa' , nLa
    WRITE(*,*) 'nLo' , nLo
    !-----------------------------------------------------
    !--- open netCDF file1
    status = nf90_OPEN(ifile1,0,ncID1) 		
      IF(status /= nf90_NoErr) CALL handle_err(status)
            WRITE(*,*) 'file ID1', ncID1
      ! time
      !status = nf90_inq_varid(ncID1,'time',timeID)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      !      WRITE(*,*) 'time ID', timeID
      !status = nf90_get_var(ncID1,timeID,VarTime)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      ! lat
      !status = nf90_inq_varid(ncID1,'lat',latID)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      !      WRITE(*,*) 'lat ID', latID
      !status = nf90_get_var(ncID1,latID,VarLat)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      ! lon
      !status = nf90_inq_varid(ncID1,'lon',lonID)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      !      WRITE(*,*) 'lon ID', lonID
      !status = nf90_get_var(ncID1,lonID,VarLon)
      !IF(status /= nf90_NoErr) CALL handle_err(status)
      
     ! surface pressure
     status = nf90_inq_varid(ncID1,'SP',spID) 
     IF(status /= nf90_NoErr) CALL handle_err(status)
        WRITE(*,*) 'sp ID', spID	       
     status = nf90_get_var(ncID1, spID, VarSP)
     IF(status /= nf90_NoErr) CALL handle_err(status)
     ! surface geopotential
     status = nf90_inq_varid(ncID1,'Z',zID)
     IF(status /= nf90_NoErr) CALL handle_err(status)
     WRITE(*,*) 'z ID', zID
     status = nf90_get_var(ncID1, zID, VarZ)
     IF(status /= nf90_NoErr) CALL handle_err(status)
     

    !--- close netCDF file1
    status = nf90_close(ncID1)
      IF(status /= nf90_NoErr) CALL handle_err(status)  



    !-----------------------------------------------------
    !--- open netCDF file2
!    status = nf90_OPEN(ifile2,0,ncID2) 		
!      IF(status /= nf90_NoErr) CALL handle_err(status)
!      WRITE(*,*) 'file ID2', ncID2    

    ! specific humidity
!    status = nf90_inq_varid(ncID2,'q',qID) 
!      IF(status /= nf90_NoErr) CALL handle_err(status)
!      WRITE(*,*) 'q ID', qID 	      
!    status = nf90_get_var(ncID2, qID, VarQ)
!      IF(status /= nf90_NoErr) CALL handle_err(status)

    !--- close netCDF file2
!    status = nf90_close(ncID2)
!      IF(status /= nf90_NoErr) CALL handle_err(status)  


    !-----------------------------------------------------
    !--- open  netCDF file2
    status = nf90_OPEN(ifile2,0,ncID2) 		
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'file ID2', ncID2    

    ! time
    status = nf90_inq_varid(ncID2,'time',timeID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'time ID', timeID	       
    status = nf90_get_var(ncID2,timeID,VarTime)		
      IF(status /= nf90_NoErr) CALL handle_err(status)

    ! lev
    status = nf90_inq_varid(ncID2,'lev',levID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'lev ID', levID	       
    status = nf90_get_var(ncID2,levID,VarLev)		
      IF(status /= nf90_NoErr) CALL handle_err(status)
 
    ! hyam
    status = nf90_inq_varid(ncID2,'hyam',hyamID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)	       
      WRITE(*,*) 'hyam ID', hyamID
    status = nf90_get_var(ncID2,hyamID,VarHyam)		
      IF(status /= nf90_NoErr) CALL handle_err(status)

!   ! hybm
    status = nf90_inq_varid(ncID2,'hybm',hybmID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)	       
      WRITE(*,*) 'hybm ID', hybmID
    status = nf90_get_var(ncID2,hybmID,VarHybm)		
      IF(status /= nf90_NoErr) CALL handle_err(status)
 
    ! lat
    status = nf90_inq_varid(ncID2,'lat',latID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'lat ID', latID	       
    status = nf90_get_var(ncID2,latID,VarLat)		
      IF(status /= nf90_NoErr) CALL handle_err(status)

    ! lon
    status = nf90_inq_varid(ncID2,'lon',lonID) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'lon ID', lonID	       
    status = nf90_get_var(ncID2,lonID,VarLon)		
    IF(status /= nf90_NoErr) CALL handle_err(status)

    ! temperature
    status = nf90_inq_varid(ncID2,'t',tID)
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 't ID', tID
    status = nf90_get_var(ncID2, tID, VarT)
    IF(status /= nf90_NoErr) CALL handle_err(status)

    ! specific humidity
    status = nf90_inq_varid(ncID2,'q',qID)
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'q ID', qID
    status = nf90_get_var(ncID2, qID, VarQ)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    

   ! mass mixing ratio for each aerosol tracer
    status = nf90_inq_varid(ncID2,Paramname(1),varID1) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID1', varID1	     
    status = nf90_get_var(ncID2, varID1, VarMix_SSs)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(2),varID2) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID2', varID2	     
    status = nf90_get_var(ncID2, varID2, VarMix_SSm)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(3),varID3) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID3', varID3	     
    status = nf90_get_var(ncID2, varID3, VarMix_SSl)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(4),varID4) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID4', varID4	     
    status = nf90_get_var(ncID2, varID4, VarMix_DUs)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(5),varID5) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID5', varID5	     
    status = nf90_get_var(ncID2, varID5, VarMix_DUm)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(6),varID6) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID6', varID6	     
    status = nf90_get_var(ncID2, varID6, VarMix_DUl)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(7),varID7) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID7', varID7	     
    status = nf90_get_var(ncID2, varID7, VarMix_OMn)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(8),varID8) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID8', varID8	     
    status = nf90_get_var(ncID2, varID8, VarMix_OMh)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(9),varID9) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID9', varID9	     
    status = nf90_get_var(ncID2, varID9, VarMix_BCn)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(10),varID10) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID10', varID10	     
    status = nf90_get_var(ncID2, varID10, VarMix_BCh)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_inq_varid(ncID2,Paramname(11),varID11) 
      IF(status /= nf90_NoErr) CALL handle_err(status)
      WRITE(*,*) 'variable ID11', varID11	     
    status = nf90_get_var(ncID2, varID11, VarMix_SU)
      IF(status /= nf90_NoErr) CALL handle_err(status)

    !--- close netCDF file2
    status = nf90_close(ncID2)
      IF(status /= nf90_NoErr) CALL handle_err(status)  



    !-----------------------------------------------------
    !--- open netCDF file4
    !status = nf90_OPEN(ifile4,0,ncID4) 		
     ! IF(status /= nf90_NoErr) CALL handle_err(status)
!      WRITE(*,*) 'file ID4', ncID4    

    ! temperature
    !status = nf90_inq_varid(ncID4,'t',tID) 
      !IF(status /= nf90_NoErr) CALL handle_err(status)
!      WRITE(*,*) 't ID', tID 	      
    !status = nf90_get_var(ncID4, tID, VarT)
     ! IF(status /= nf90_NoErr) CALL handle_err(status)

    !--- close netCDF file4
    !status = nf90_close(ncID4)
     ! IF(status /= nf90_NoErr) CALL handle_err(status)  


    !WRITE(*,*) 'VarSP', VarSP(1,1,:), ' VarQ10', VarQ(1,1,10,:),' VarQ50', VarQ(1,1,50,:)
    !WRITE(*,*) 'VarHybm', VarHybm(:), ' VarHyam', VarHyam(:)


      DO lo = 1,nLo 
       DO la = 1,nLa
        DO it = 1,nT		  
        VarP(lo,la,:,it) = VarHybm(:)*VarSP(lo,la,it)+VarHyam(:)
        ENDDO
       ENDDO
      ENDDO

    !WRITE(*,*) 'VarP10', VarP(1,1,10,:), ' VarP50', VarP(1,1,50,:)



  END SUBROUTINE netcdf_input  



! SUBROUTINE netcdf_output(ofile1, ofile2, ofile3, year, mon, &
! 			 nT, nLe, nLa, nLo, &
!			 VarTime, VarLev, VarLat, VarLon, &
!			 VarMC, VarNC, VarCCN_02, VarCCN_04, VarCCN_10, &
!			 VarMC_SSs, VarMC_SSm, VarMC_SSl, &
!			 VarMC_DUs, VarMC_DUm, VarMC_DUl, &
!			 VarMC_OMn, VarMC_OMh, VarMC_BCn, VarMC_BCh, VarMC_SU, &
!			 VarNC_SSs, VarNC_SSm, VarNC_SSl, &
!			 VarNC_DUs, VarNC_DUm, VarNC_DUl, &
!			 VarNC_OMn, VarNC_OMh, VarNC_BCn, VarNC_BCh, VarNC_SU, &
!			 VarCCN_02_SSs, VarCCN_02_SSm, VarCCN_02_SSl, &
!			 VarCCN_02_OMh, VarCCN_02_BCh, VarCCN_02_SU, &
!			 VarCCN_04_SSs, VarCCN_04_SSm, VarCCN_04_SSl, &
!			 VarCCN_04_OMh, VarCCN_04_BCh, VarCCN_04_SU, &
!			 VarCCN_10_SSs, VarCCN_10_SSm, VarCCN_10_SSl, &
!			 VarCCN_10_OMh, VarCCN_10_BCh, VarCCN_10_SU)

  SUBROUTINE netcdf_output(ofile1, year, mon, &
  			 nT, nLe, nLa, nLo, &
 			 VarTime, VarLev, VarLat, VarLon, &
			 VarCCN_01, VarCCN_02, VarCCN_03, &
                         VarCCN_04, VarCCN_05, VarCCN_06, &
                         VarCCN_07, VarCCN_08, VarCCN_09, &
                         VarCCN_10, geop)


  USE mo_kind,          ONLY: dp
  !USE mo_box_aerosols,	ONLY: aertype
  !USE mo_ham_ccnctl,	ONLY: nsat 


  IMPLICIT NONE

  ! Input parameters
  CHARACTER(len=*), INTENT(in)       			:: ofile1
  CHARACTER(len=4), INTENT(in)       			:: year
  CHARACTER(len=2), INTENT(in)       			:: mon
  INTEGER,  INTENT(in)         				:: nT, nLe, nLa, nLo  
  REAL(dp), INTENT(in), DIMENSION(nT)			:: VarTime
  REAL(dp), INTENT(in), DIMENSION(nLe)			:: VarLev
  REAL(dp), INTENT(in), DIMENSION(nLa)			:: VarLat
  REAL(dp), INTENT(in), DIMENSION(nLo)			:: VarLon
  
  REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)      	:: VarCCN_01, VarCCN_02, VarCCN_03
  REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarCCN_04, VarCCN_05, VarCCN_06
  REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarCCN_07, VarCCN_08, VarCCN_09, VarCCN_10
  REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: geop


! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)      	:: VarMC, VarNC, VarCCN_02, VarCCN_04, VarCCN_10
! 
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarMC_SSs, VarMC_SSm, VarMC_SSl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarMC_DUs, VarMC_DUm, VarMC_DUl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarMC_OMn, VarMC_OMh, VarMC_BCn, VarMC_BCh, VarMC_SU  
! 
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarNC_SSs, VarNC_SSm, VarNC_SSl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarNC_DUs, VarNC_DUm, VarNC_DUl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarNC_OMn, VarNC_OMh, VarNC_BCn, VarNC_BCh, VarNC_SU  
! 
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT) 	:: VarCCN_02_SSs, VarCCN_02_SSm, VarCCN_02_SSl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT) 	:: VarCCN_02_OMh, VarCCN_02_BCh, VarCCN_02_SU
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT) 	:: VarCCN_04_SSs, VarCCN_04_SSm, VarCCN_04_SSl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarCCN_04_OMh, VarCCN_04_BCh, VarCCN_04_SU  
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarCCN_10_SSs, VarCCN_10_SSm, VarCCN_10_SSl
! REAL(dp), INTENT(in), DIMENSION(nLo,nLa,nLe,nT)       :: VarCCN_10_OMh, VarCCN_10_BCh, VarCCN_10_SU  
 

 
  ! Local parameters
!  INTEGER :: ncid1, ncid2, ncid3, status 
  INTEGER :: ncid1, status 
  INTEGER :: timeID, levID, latID, lonID
  INTEGER :: var1, var2, var3, var4
  INTEGER :: var7, var8, var9, var10
  INTEGER :: var11, var12, var13, var14
  INTEGER :: var15, var16, var17

!  INTEGER :: var5, var6, var7, var8, var9
!  INTEGER :: var10, var11, var12, var13, var14, var15, var16, var17, var18, var19, var20
!  INTEGER :: var21, var22, var23, var24, var25, var26, var27, var28, var29, var30, var31
!  INTEGER :: var33, var34, var35, var36, var37, var38
!  INTEGER :: var40, var41, var42, var43, var44, var45
!  INTEGER :: var47, var48, var49, var50, var51, var52
! !free Integer :: var32, var39, var46

!  INTEGER :: var102, var202, var302, var402
!  INTEGER :: var103, var203, var303, var403
  
            
  REAL, PARAMETER	:: FillAtt = -9999.99
  REAL(dp), PARAMETER	:: FillAtt2 = -9999.99
  CHARACTER(len=31)	:: timestr = "hours since XXXX-YY-01 00:00:00"      
!  CHARACTER(len=67)	:: sscribe = "CCN_02 = 0.2%SuperSat, CCN_04 = 0.4%SuperSat, CCN_10 = 1.0%SuperSat"
!  CHARACTER(len=337)	:: specscribe = "SSs = small mode sea salt, &
!  SSm = medium mode sea salt, &
!  SSl = large mode sea salt, &
!  DUs = small mode dust, &
!  DUm = medium mode dust, &
!  DUl = large mode dust, &
!  OMn = hydrophobic organic matter, &
!  OMh = hydrophilic organic matter, &
!  BCn = hydrophobic black carbon, &
!  BCh = hydrophilic black carbon, &
!  SO4 = Sulfate"

!  WRITE(*,*) 'sscribe ', sscribe
!  WRITE(*,*) 'specscribe ', specscribe
!  WRITE(*,*) 'timestr ', timestr

  !--- Define the parameters
  timestr(13:16) = year
  timestr(18:19) = mon

!  WRITE(*,*) 'timestr new ', timestr
   
   !------------------------------------------------------------------------
   !--- create output netCDF file
   status = nf90_create(path = ofile1, cmode = NF90_NETCDF4, ncid = ncid1)
    IF(status /= nf90_noerr) CALL handle_err(status)
!   WRITE(*,*) 'def nc'

   !--- Define the dimensions
   status = nf90_def_dim(ncid1, "time", nf90_unlimited, timeID)
    IF(status /= nf90_noerr) CALL handle_err(status)
   status = nf90_def_dim(ncid1, "lev", nLe, levID)
    IF(status /= nf90_noerr) CALL handle_err(status)   
   status = nf90_def_dim(ncid1, "lat", nLa, latID)
    IF(status /= nf90_noerr) CALL handle_err(status)
   status = nf90_def_dim(ncid1, "lon", nLo, lonID)
    IF(status /= nf90_noerr) CALL handle_err(status)
!   WRITE(*,*) 'def dim'

   !--- Define global attributes
!   status = nf90_put_att(ncid1, NF90_GLOBAL, "Supersaturation", sscribe)
!   WRITE(*,*) 'att'

   ! Time   
   status = nf90_def_var(ncid1, "time", nf90_double,(/ timeID /), var1)
    IF(status /= nf90_NoErr) CALL handle_err(status)
   status = nf90_put_att(ncid1, var1, "long_name", "time")
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var1, "units", timestr)
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var1, "calendar", "proleptic_gregorian")
    IF(status /= nf90_NoErr) CALL handle_err(status)	       
   status = nf90_put_att(ncid1, var1, "_FillValue", FillAtt2)
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
!   WRITE(*,*) 'def time'
   
   ! Lev 
   status = nf90_def_var(ncid1, "lev", nf90_double,(/ levID /), var2)
    IF(status /= nf90_NoErr) CALL handle_err(status)
   status = nf90_put_att(ncid1, var2, "long_name", "model level")
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var2, "units", "level")
    IF(status /= nf90_NoErr) CALL handle_err(status)	    
   status = nf90_put_att(ncid1, var2, "_FillValue", FillAtt2)
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
!   WRITE(*,*) 'def lev'
   
   ! Lat   
   status = nf90_def_var(ncid1, "lat", nf90_double,(/ latID /), var3)
    IF(status /= nf90_NoErr) CALL handle_err(status)
   status = nf90_put_att(ncid1, var3, "long_name", "latitude")
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var3, "units", "degrees_north")
    IF(status /= nf90_NoErr) CALL handle_err(status)	    
   status = nf90_put_att(ncid1, var3, "_FillValue", FillAtt2)
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var3, "axis", "Y")
    IF(status /= nf90_NoErr) CALL handle_err(status)
!   WRITE(*,*) 'def lat'

   ! Lon  
   status = nf90_def_var(ncid1, "lon", nf90_double,(/ lonID /), var4)
    IF(status /= nf90_NoErr) CALL handle_err(status)
   status = nf90_put_att(ncid1, var4, "long_name", "longitude")
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var4, "units", "degrees_east")
    IF(status /= nf90_NoErr) CALL handle_err(status)	    
   status = nf90_put_att(ncid1, var4, "_FillValue", FillAtt2)
    IF(status /= nf90_NoErr) CALL handle_err(status)	   
   status = nf90_put_att(ncid1, var4, "axis", "X")
    IF(status /= nf90_NoErr) CALL handle_err(status)
!   WRITE(*,*) 'def lon'

!  ! MCONC
!  status = nf90_def_var(ncid1, "MCONC", nf90_float,(/ lonID, latID, levID, timeID /), var5)
!   IF(status /= nf90_NoErr) CALL handle_err(status)
!  status = nf90_put_att(ncid1, var5, "long_name", "total aerosol mass concentration")
!   IF(status /= nf90_NoErr) CALL handle_err(status)	  
!  status = nf90_put_att(ncid1, var5, "units", "kg m**-3")
!   IF(status /= nf90_NoErr) CALL handle_err(status)	   
!  status = nf90_put_att(ncid1, var5, "_FillValue", FillAtt)
!   IF(status /= nf90_NoErr) CALL handle_err(status)	  
!  WRITE(*,*) 'def mc'

!  ! NCONC  
!  status = nf90_def_var(ncid1, "NCONC", nf90_float,(/ lonID, latID, levID, timeID /), var6)
!   IF(status /= nf90_NoErr) CALL handle_err(status)
!  status = nf90_put_att(ncid1, var6, "long_name", "total aerosol number concentration")
!   IF(status /= nf90_NoErr) CALL handle_err(status)	  
!  status = nf90_put_att(ncid1, var6, "units", "m**-3")
!   IF(status /= nf90_NoErr) CALL handle_err(status)	   
!  status = nf90_put_att(ncid1, var6, "_FillValue", FillAtt)
!   IF(status /= nf90_NoErr) CALL handle_err(status)	  
    !  WRITE(*,*) 'def nc'

    ! CCN_01
    status = nf90_def_var(ncid1, "CCN_01", nf90_float,(/ lonID, latID, levID, timeID /), var7)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var7, "long_name", "cloud condensation nuclei at w=0.01m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var7, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var7, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn01'
    
    ! CCN_02
    status = nf90_def_var(ncid1, "CCN_02", nf90_float,(/ lonID, latID, levID, timeID /), var8)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var8, "long_name", "cloud condensation nuclei at w=0.0278m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var8, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var8, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn02'
    
    ! CCN_03
    status = nf90_def_var(ncid1, "CCN_03", nf90_float,(/ lonID, latID, levID, timeID /), var9)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var9, "long_name", "cloud condensation nuclei at w=0.0774m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var9, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var9, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn03'
    
    ! CCN_04
    status = nf90_def_var(ncid1, "CCN_04", nf90_float,(/ lonID, latID, levID, timeID /), var10)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var10, "long_name", "cloud condensation nuclei at w=0.215m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var10, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var10, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn04'
    
    ! CCN_05
    status = nf90_def_var(ncid1, "CCN_05", nf90_float,(/ lonID, latID, levID, timeID /), var11)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var11, "long_name", "cloud condensation nuclei at w=0.599m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var11, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var11, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn05'

    ! CCN_06
    status = nf90_def_var(ncid1, "CCN_06", nf90_float,(/ lonID, latID, levID, timeID /), var12)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var12, "long_name", "cloud condensation nuclei at w=1.67m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var12, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var12, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn06'

    ! CCN_07
    status = nf90_def_var(ncid1, "CCN_07", nf90_float,(/ lonID, latID, levID, timeID /), var13)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var13, "long_name", "cloud condensation nuclei at w=4.64m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var13, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var13, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn07'
    
    ! CCN_08
    status = nf90_def_var(ncid1, "CCN_08", nf90_float,(/ lonID, latID, levID, timeID /), var14)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var14, "long_name", "cloud condensation nuclei at w=12.9m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var14, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var14, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn08'
    
    ! CCN_09
    status = nf90_def_var(ncid1, "CCN_09", nf90_float,(/ lonID, latID, levID, timeID /), var15)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var15, "long_name", "cloud condensation nuclei at w=35.9m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var15, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var15, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn09'

    ! CCN_10
    status = nf90_def_var(ncid1, "CCN_10", nf90_float,(/ lonID, latID, levID, timeID /), var16)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var16, "long_name", "cloud condensation nuclei at w=100m/s")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var16, "units", "m**-3")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var16, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def ccn10'
    ! Geop
    status = nf90_def_var(ncid1, "z", nf90_float,(/ lonID, latID, levID, timeID /), var17)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var17, "long_name", "geopothential")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var16, "units", "m")
    IF(status /= nf90_NoErr) CALL handle_err(status)
    status = nf90_put_att(ncid1, var16, "_FillValue", FillAtt)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'def geop'
    

   !--- Leave define mode
   status = nf90_enddef(ncid1)
    IF(status /= nf90_NoErr) CALL handle_err(status)  

   !--- write data
   status = nf90_put_var(ncid1, var1, VarTime)	        ! Time  
    IF(status /= nf90_NoErr) CALL handle_err(status) 
!   WRITE(*,*) 'time'
   status = nf90_put_var(ncid1, var2, VarLev)	        ! Lev  
    IF(status /= nf90_NoErr) CALL handle_err(status) 
!   WRITE(*,*) 'lev'
   status = nf90_put_var(ncid1, var3, VarLat)	        ! Lat  
    IF(status /= nf90_NoErr) CALL handle_err(status) 
!   WRITE(*,*) 'lat'
   status = nf90_put_var(ncid1, var4, VarLon)	        ! Lon  
    IF(status /= nf90_NoErr) CALL handle_err(status)
!   WRITE(*,*) 'lon'
!   status = nf90_put_var(ncid1, var5, REAL(VarMC))      ! MC	       
!    IF(status /= nf90_NoErr) CALL handle_err(status)
!   WRITE(*,*) 'mc'
!   status = nf90_put_var(ncid1, var6, REAL(VarNC))      ! NC	       
!    IF(status /= nf90_NoErr) CALL handle_err(status)
!   WRITE(*,*) 'nc'
    status = nf90_put_var(ncid1, var7, REAL(VarCCN_01))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_01'
    status = nf90_put_var(ncid1, var8, REAL(VarCCN_02))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_02'
    status = nf90_put_var(ncid1, var9, REAL(VarCCN_03))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_03'
    status = nf90_put_var(ncid1, var10, REAL(VarCCN_04))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_04'
    status = nf90_put_var(ncid1, var11, REAL(VarCCN_05))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_05'
    status = nf90_put_var(ncid1, var12, REAL(VarCCN_06))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_06'
    status = nf90_put_var(ncid1, var13, REAL(VarCCN_07))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_07'
    status = nf90_put_var(ncid1, var14, REAL(VarCCN_08))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_08'
    status = nf90_put_var(ncid1, var15, REAL(VarCCN_09))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_09'
    status = nf90_put_var(ncid1, var16, REAL(VarCCN_10))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'ccn_10'
    status = nf90_put_var(ncid1, var17, REAL(geop))
    IF(status /= nf90_NoErr) CALL handle_err(status)
    !   WRITE(*,*) 'geop'
    
   !--- close output netCDF file
   status = nf90_close(ncid1)
    IF(status /= nf90_NoErr) CALL handle_err(status)  

!  WRITE(*,*) ' '  
!  WRITE(*,*) 'FINISHED OUTPUT 1'


! !-------------------------------------------------------------------------
! !--- create output netCDF file
! status = nf90_create(path = ofile2, cmode = NF90_NETCDF4, ncid = ncid2)
!  IF(status /= nf90_noerr) CALL handle_err(status)

! !--- Define the dimensions
! status = nf90_def_dim(ncid2, "time", nf90_unlimited, timeID)
!  IF(status /= nf90_noerr) CALL handle_err(status)
! status = nf90_def_dim(ncid2, "lev", nLe, levID)
!  IF(status /= nf90_noerr) CALL handle_err(status)   
! status = nf90_def_dim(ncid2, "lat", nLa, latID)
!  IF(status /= nf90_noerr) CALL handle_err(status)
! status = nf90_def_dim(ncid2, "lon", nLo, lonID)
!  IF(status /= nf90_noerr) CALL handle_err(status)

! !--- Define global attributes
! status = nf90_put_att(ncid2, NF90_GLOBAL, "Aerosol Species", specscribe)  
! 
! ! Time   
! status = nf90_def_var(ncid2, "time", nf90_double,(/ timeID /), var102)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var102, "long_name", "time")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var102, "units", timestr)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var102, "calendar", "proleptic_gregorian")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	       
! status = nf90_put_att(ncid2, var102, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! 
! !Lev
! status = nf90_def_var(ncid2, "lev", nf90_double,(/ levID /), var202)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var202, "long_name", "model level")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var202, "units", "level")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid2, var202, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! !Lat
! status = nf90_def_var(ncid2, "lat", nf90_double,(/ latID /), var302)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var302, "long_name", "latitude")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var302, "units", "degrees_north")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid2, var302, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var302, "axis", "Y")
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! !Lon
! status = nf90_def_var(ncid2, "lon", nf90_double,(/ lonID /), var402)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var402, "long_name", "longitude")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var402, "units", "degrees_east")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid2, var402, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var402, "axis", "X")
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! ! MCONC
! status = nf90_def_var(ncid2, "MCONC_SSs", nf90_float,(/ lonID, latID, levID, timeID /), var10)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var10, "long_name", "mass concentration of small mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var10, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var10, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_SSm", nf90_float,(/ lonID, latID, levID, timeID /), var11)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var11, "long_name", "mass concentration of medium mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var11, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var11, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_SSl", nf90_float,(/ lonID, latID, levID, timeID /), var12)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var12, "long_name", "mass concentration of large mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var12, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var12, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_DUs", nf90_float,(/ lonID, latID, levID, timeID /), var13)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var13, "long_name", "mass concentration of small mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var13, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var13, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_DUm", nf90_float,(/ lonID, latID, levID, timeID /), var14)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var14, "long_name", "mass concentration of medium mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var14, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var14, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_DUl", nf90_float,(/ lonID, latID, levID, timeID /), var15)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var15, "long_name", "mass concentration of large mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var15, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var15, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_OMn", nf90_float,(/ lonID, latID, levID, timeID /), var16)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var16, "long_name", "mass concentration of non-hydrophilic organic matter")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var16, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var16, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_OMh", nf90_float,(/ lonID, latID, levID, timeID /), var17)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var17, "long_name", "mass concentration of hydrophilic organic matter")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var17, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var17, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_BCn", nf90_float,(/ lonID, latID, levID, timeID /), var18)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var18, "long_name", "mass concentration of non-hydrophilic black carbon")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var18, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var18, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "MCONC_BCh", nf90_float,(/ lonID, latID, levID, timeID /), var19)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var19, "long_name", "mass concentration of hydrophilic black carbon")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var19, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var19, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_def_var(ncid2, "MCONC_SO4", nf90_float,(/ lonID, latID, levID, timeID /), var20)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var20, "long_name", "mass concentration of sulfate")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var20, "units", "kg m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var20, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! ! NCONC  
! status = nf90_def_var(ncid2, "NCONC_SSs", nf90_float,(/ lonID, latID, levID, timeID /), var21)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var21, "long_name", "number concentration of small mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var21, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var21, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_SSm", nf90_float,(/ lonID, latID, levID, timeID /), var22)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var22, "long_name", "number concentration of medium mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var22, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var22, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_SSl", nf90_float,(/ lonID, latID, levID, timeID /), var23)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var23, "long_name", "number concentration of large mode sea salt")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var23, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var23, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_DUs", nf90_float,(/ lonID, latID, levID, timeID /), var24)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var24, "long_name", "number concentration of small mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var24, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var24, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_DUm", nf90_float,(/ lonID, latID, levID, timeID /), var25)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var25, "long_name", "number concentration of medium mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var25, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var25, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_DUl", nf90_float,(/ lonID, latID, levID, timeID /), var26)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var26, "long_name", "number concentration of large mode dust")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var26, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var26, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_OMn", nf90_float,(/ lonID, latID, levID, timeID /), var27)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var27, "long_name", "number concentration of non-hydrophilic organic matter")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var27, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var27, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_OMh", nf90_float,(/ lonID, latID, levID, timeID /), var28)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var28, "long_name", "number concentration of hydrophilic organic matter")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var28, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var28, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_BCn", nf90_float,(/ lonID, latID, levID, timeID /), var29)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var29, "long_name", "number concentration of non-hydrophilic black carbon")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var29, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var29, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  

! status = nf90_def_var(ncid2, "NCONC_BCh", nf90_float,(/ lonID, latID, levID, timeID /), var30)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var30, "long_name", "number concentration of hydrophilic black carbon")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var30, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var30, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_def_var(ncid2, "NCONC_SO4", nf90_float,(/ lonID, latID, levID, timeID /), var31)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid2, var31, "long_name", "number concentration of sulfate")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	  
! status = nf90_put_att(ncid2, var31, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid2, var31, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! !--- Leave define mode
! status = nf90_enddef(ncid2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)  

! !--- write data
! status = nf90_put_var(ncid2, var102, VarTime)	       ! Time  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid2, var202, VarLev)	       ! Lev  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid2, var302, VarLat)	       ! Lat  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid2, var402, VarLon)	       ! Lon  
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_put_var(ncid2, var10, REAL(VarMC_SSs))         ! MC	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var11, REAL(VarMC_SSm))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var12, REAL(VarMC_SSl))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var13, REAL(VarMC_DUs))  	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var14, REAL(VarMC_DUm))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var15, REAL(VarMC_DUl))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var16, REAL(VarMC_OMn))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)	
! status = nf90_put_var(ncid2, var17, REAL(VarMC_OMh))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var18, REAL(VarMC_BCn))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)   
! status = nf90_put_var(ncid2, var19, REAL(VarMC_BCh))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var20, REAL(VarMC_SU))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_put_var(ncid2, var21, REAL(VarNC_SSs))         ! NC	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var22, REAL(VarNC_SSm))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var23, REAL(VarNC_SSl))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var24, REAL(VarNC_DUs))  	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var25, REAL(VarNC_DUm))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var26, REAL(VarNC_DUl))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var27, REAL(VarNC_OMn))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)	
! status = nf90_put_var(ncid2, var28, REAL(VarNC_OMh))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var29, REAL(VarNC_BCn))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)   
! status = nf90_put_var(ncid2, var30, REAL(VarNC_BCh))         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid2, var31, REAL(VarNC_SU))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! !--- close output netCDF file
! status = nf90_close(ncid2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)  

! WRITE(*,*) ' '
! WRITE(*,*) 'FINISHED OUTPUT 2'

!!--------------------------------------------------------------------------
!!--- create output netCDF file
! status = nf90_create(path = ofile3, cmode = NF90_NETCDF4, ncid = ncid3)
!  IF(status /= nf90_noerr) CALL handle_err(status)

! !--- Define the dimensions
! status = nf90_def_dim(ncid3, "time", nf90_unlimited, timeID)
!  IF(status /= nf90_noerr) CALL handle_err(status)
! status = nf90_def_dim(ncid3, "lev", nLe, levID)
!  IF(status /= nf90_noerr) CALL handle_err(status)   
! status = nf90_def_dim(ncid3, "lat", nLa, latID)
!  IF(status /= nf90_noerr) CALL handle_err(status)
! status = nf90_def_dim(ncid3, "lon", nLo, lonID)
!  IF(status /= nf90_noerr) CALL handle_err(status)

! !--- Define global attributes 
! status = nf90_put_att(ncid3, NF90_GLOBAL, "Supersaturation", sscribe)
! status = nf90_put_att(ncid3, NF90_GLOBAL, "Aerosol Species", specscribe)
! 
! ! Time   
! status = nf90_def_var(ncid3, "time", nf90_double,(/ timeID /), var103)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var103, "long_name", "time")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var103, "units", timestr)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var103, "calendar", "proleptic_gregorian")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	       
! status = nf90_put_att(ncid3, var103, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! !Lev
! status = nf90_def_var(ncid3, "lev", nf90_double,(/ levID /), var203)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var203, "long_name", "model level")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var203, "units", "level")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var203, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! !Lat
! status = nf90_def_var(ncid3, "lat", nf90_double,(/ latID /), var303)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var303, "long_name", "latitude")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var303, "units", "degrees_north")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var303, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var303, "axis", "Y")
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! !Lon
! status = nf90_def_var(ncid3, "lon", nf90_double,(/ lonID /), var403)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var403, "long_name", "longitude")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var403, "units", "degrees_east")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var403, "_FillValue", FillAtt2)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var403, "axis", "X")
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! ! CCN_02      
! status = nf90_def_var(ncid3, "CCN_02_SSs", nf90_float,(/ lonID, latID, levID, timeID /), var33)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var33, "long_name", "cloud condensation nuclei from small mode sea salt at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var33, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var33, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_02_SSm", nf90_float,(/ lonID, latID, levID, timeID /), var34)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var34, "long_name", "cloud condensation nuclei from medium mode sea salt at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var34, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var34, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_02_SSl", nf90_float,(/ lonID, latID, levID, timeID /), var35)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var35, "long_name", "cloud condensation nuclei from large mode sea salt at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var35, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var35, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_def_var(ncid3, "CCN_02_OMh", nf90_float,(/ lonID, latID, levID, timeID /), var36)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var36, "long_name", "cloud condensation nuclei from hydrophilic organic matter at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var36, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var36, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_02_BCh", nf90_float,(/ lonID, latID, levID, timeID /), var37)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var37, "long_name", "cloud condensation nuclei from hydrophilic black carbon at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var37, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var37, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_02_SO4", nf90_float,(/ lonID, latID, levID, timeID /), var38)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var38, "long_name", "cloud condensation nuclei from sulfate at 0.2% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var38, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var38, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! ! CCN_04       
! status = nf90_def_var(ncid3, "CCN_04_SSs", nf90_float,(/ lonID, latID, levID, timeID /), var40)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var40, "long_name", "cloud condensation nuclei from small mode sea salt at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var40, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var40, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_04_SSm", nf90_float,(/ lonID, latID, levID, timeID /), var41)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var41, "long_name", "cloud condensation nuclei from medium mode sea salt at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var41, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var41, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_04_SSl", nf90_float,(/ lonID, latID, levID, timeID /), var42)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var42, "long_name", "cloud condensation nuclei from large mode sea salt at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var42, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var42, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_def_var(ncid3, "CCN_04_OMh", nf90_float,(/ lonID, latID, levID, timeID /), var43)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var43, "long_name", "cloud condensation nuclei from hydrophilic organic matter at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var43, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var43, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_04_BCh", nf90_float,(/ lonID, latID, levID, timeID /), var44)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var44, "long_name", "cloud condensation nuclei from hydrophilic black carbon at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var44, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var44, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_04_SO4", nf90_float,(/ lonID, latID, levID, timeID /), var45)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var45, "long_name", "cloud condensation nuclei from sulfate at 0.4% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var45, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var45, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! ! CCN_10 
! status = nf90_def_var(ncid3, "CCN_10_SSs", nf90_float,(/ lonID, latID, levID, timeID /), var47)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var47, "long_name", "cloud condensation nuclei from small mode sea salt at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var47, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var47, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_10_SSm", nf90_float,(/ lonID, latID, levID, timeID /), var48)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var48, "long_name", "cloud condensation nuclei from medium mode sea salt at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var48, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var48, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_10_SSl", nf90_float,(/ lonID, latID, levID, timeID /), var49)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var49, "long_name", "cloud condensation nuclei from large mode sea salt at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var49, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var49, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_def_var(ncid3, "CCN_10_OMh", nf90_float,(/ lonID, latID, levID, timeID /), var50)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var50, "long_name", "cloud condensation nuclei from hydrophilic organic matter at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var50, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var50, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_10_BCh", nf90_float,(/ lonID, latID, levID, timeID /), var51)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var51, "long_name", "cloud condensation nuclei from hydrophilic black carbon at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var51, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var51, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   

! status = nf90_def_var(ncid3, "CCN_10_SO4", nf90_float,(/ lonID, latID, levID, timeID /), var52)
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_att(ncid3, var52, "long_name", "cloud condensation nuclei from sulfate at 1.0% ssat")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	   
! status = nf90_put_att(ncid3, var52, "units", "m**-3")
!  IF(status /= nf90_NoErr) CALL handle_err(status)	    
! status = nf90_put_att(ncid3, var52, "_FillValue", FillAtt)
!  IF(status /= nf90_NoErr) CALL handle_err(status)   

! !--- Leave define mode
! status = nf90_enddef(ncid3)
!  IF(status /= nf90_NoErr) CALL handle_err(status)  

! !--- write data
! status = nf90_put_var(ncid3, var103, VarTime)	       ! Time  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid3, var203, VarLev)	       ! Lev  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid3, var303, VarLat)	       ! Lat  
!  IF(status /= nf90_NoErr) CALL handle_err(status) 
! status = nf90_put_var(ncid3, var403, VarLon)	       ! Lon  
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_put_var(ncid3, var33, REAL(VarCCN_02_SSs))     !CCN_02         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var34, REAL(VarCCN_02_SSm))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var35, REAL(VarCCN_02_SSl))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var36, REAL(VarCCN_02_OMh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var37, REAL(VarCCN_02_BCh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var38, REAL(VarCCN_02_SU))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_put_var(ncid3, var40, REAL(VarCCN_04_SSs))    !CCN_04		      
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var41, REAL(VarCCN_04_SSm))	      
!   IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var42, REAL(VarCCN_04_SSl))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var43, REAL(VarCCN_04_OMh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var44, REAL(VarCCN_04_BCh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var45, REAL(VarCCN_04_SU))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)

! status = nf90_put_var(ncid3, var47, REAL(VarCCN_10_SSs))     !CCN_10         
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var48, REAL(VarCCN_10_SSm))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var49, REAL(VarCCN_10_SSl))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var50, REAL(VarCCN_10_OMh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var51, REAL(VarCCN_10_BCh))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
! status = nf90_put_var(ncid3, var52, REAL(VarCCN_10_SU))	       
!  IF(status /= nf90_NoErr) CALL handle_err(status)
!
! !--- close output netCDF file
! status = nf90_close(ncid3)
!  IF(status /= nf90_NoErr) CALL handle_err(status)  
! 
! WRITE(*,*) ' '
! WRITE(*,*) 'FINISHED OUTPUT 3'

  END SUBROUTINE netcdf_output  



  !--- ERROR check:
  SUBROUTINE handle_err(status)
       INTEGER, INTENT (in) :: status
    
       IF(status /= nf90_noerr) then
	 PRINT *, trim(nf90_strerror(status))
	 STOP "Stopped"
       END IF
  END SUBROUTINE handle_err	 
                  
		  
		   
END MODULE mo_box_netcdf_io




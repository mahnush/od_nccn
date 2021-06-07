MODULE mo_box_input

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  !--- Subroutines:

CONTAINS

!  SUBROUTINE read_input(ifilename1, ifilename2, ifilename3, ifilename4, &
!  			ofilename1, ofilename2, ofilename3, year, mon)

  SUBROUTINE read_input(ifilename1, ifilename2, &
  			ofilename1, year, mon)

  IMPLICIT NONE
  
    ! Output parameters
    CHARACTER(len=150), INTENT(out)	:: ifilename1, ifilename2            
    CHARACTER(len=150), INTENT(out)	:: ofilename1
    CHARACTER(len=4), INTENT(out)	:: year
    CHARACTER(len=2), INTENT(out)	:: mon
    
    ! Local parameters   
    LOGICAL		:: fileexists


    INQUIRE(file='boxfile.txt',exist=fileexists)
    IF (fileexists) then
      OPEN(unit=10,file='boxfile.txt',action='READ',status='old')
      !WRITE(*,*) ' '     
      !WRITE(*,*) 'BOXFILE SETTINGS:'   
      READ(10,'(A)') ifilename1  ! read input filename to know which file to open
      !WRITE(*,*) 'ifilename1:',  ifilename1  

      READ(10,'(A)') ifilename2  ! read input filename to know which file to open
      !WRITE(*,*) 'ifilename2:',  ifilename2  

      !READ(10,'(A)') ifilename3  ! read input filename to know which file to open
      !WRITE(*,*) 'ifilename3:',  ifilename3  

      !READ(10,'(A)') ifilename4  ! read input filename to know which file to open
      !WRITE(*,*) 'ifilename4:',  ifilename4  

      READ(10,'(A)') ofilename1  ! read output filename to know which file to create
      !WRITE(*,*) 'ofilename1:', ofilename1        

!      READ(10,'(A)') ofilename2  ! read output filename to know which file to create
!      !WRITE(*,*) 'ofilename2:', ofilename2    

!      READ(10,'(A)') ofilename3  ! read output filename to know which file to create
!      !WRITE(*,*) 'ofilename3:', ofilename3    

      READ(10,'(A)') year  ! read output filename to know which file to create
      !WRITE(*,*) 'year:', year    

      READ(10,'(A)') mon  ! read output filename to know which file to create
      !WRITE(*,*) 'mon:', mon    


      CLOSE(unit=10)
    ELSE
      STOP "Stopped: boxfile.txt does not exist"  
    ENDIF

  END SUBROUTINE read_input


  SUBROUTINE get_input_sub(jpt1, jpt2, tlucuaw)

  USE mo_kind,			ONLY: dp  
  USE mo_constants,		ONLY: rd, rv

  IMPLICIT NONE	 


    ! Input parameters
    INTEGER, INTENT(in) ::	jpt1, jpt2

    ! Output parameters
    REAL(dp), INTENT(out)::	tlucuaw(jpt1-1:jpt2+1)
    
    ! Local parameters    	
    REAL(dp), PARAMETER :: cavl1 = -6096.9385_dp  
    REAL(dp), PARAMETER :: cavl2 =    21.2409642_dp
    REAL(dp), PARAMETER :: cavl3 =    -2.711193_dp
    REAL(dp), PARAMETER :: cavl4 =     1.673952_dp
    REAL(dp), PARAMETER :: cavl5 =     2.433502_dp 
    REAL(dp), PARAMETER :: fdeltat = 0.001_dp  			! distance of table knots [K]    
    INTEGER		:: jt
    REAL(dp)	      :: ztt, zlinner 

    DO jt = jpt1-1, jpt2+1
      ztt = fdeltat*jt
      zlinner  = (cavl1/ztt+cavl2+cavl3*0.01_dp*ztt+cavl4*ztt*ztt*1.e-5_dp+cavl5*LOG(ztt))
      tlucuaw(jt) = EXP(zlinner)*rd/rv
    ENDDO  
  
  END SUBROUTINE get_input_sub 
  
  
  SUBROUTINE get_input_para(tlucuaw,jpt1,jpt2,ptm1,papm1,pqm1,zesw_2d,zrho,zrho_dry)

  USE mo_kind,			ONLY: dp  
  USE mo_box_settings,   	ONLY: kbdim, klev, kproma, ktdia
  USE mo_constants,		ONLY: rd, rv, vtmpc1, amd, argas

  IMPLICIT NONE

    ! Input parameters
    INTEGER, INTENT(in) :: jpt1, jpt2
    REAL(dp), INTENT(in):: tlucuaw(jpt1-1:jpt2+1)  ! table - Es*Rd/Rv, water phase only

    REAL(dp), INTENT(in):: ptm1(kbdim,klev), 	&
			   pqm1(kbdim,klev), 	&
			   papm1(kbdim,klev)

    ! Output parameters
    REAL(dp), INTENT(out)::	zesw_2d(kbdim,klev), 	&
				zrho(kbdim,klev), 	&
				zrho_dry(kbdim,klev)
		
    ! Local parameters  
    REAL(dp), PARAMETER :: pxlm1(kbdim,klev) = 0.0_dp	!liquid water content assumption at cloud base
    REAL(dp), PARAMETER :: pxim1(kbdim,klev) = 0.0_dp	!ice water content assumption at cloud base
    INTEGER		:: jl, jk
    REAL(dp)		:: zesw
    REAL(dp) 		:: ztmp(kbdim,klev),ptvm1(kbdim,klev)
    INTEGER		:: itm1_look(kbdim,klev) !< index for temperature lookup table
 

    ztmp(1:kproma,:)      = 1000._dp*ptm1(1:kproma,:)
    itm1_look(1:kproma,:) = NINT(ztmp(1:kproma,:))

    DO jk=klev,ktdia,-1
      DO jl=1,kproma
        zesw = tlucuaw(itm1_look(jl,jk))/papm1(jl,jk)        
        zesw_2d(jl,jk) = zesw*papm1(jl,jk)*rv/rd	!Saturation water vapour pressure
      END DO !jl
    ENDDO !jk    	      

    ptvm1(1:kproma,:) = ptm1(1:kproma,:)*(1._dp+vtmpc1*pqm1(1:kproma,:)-(pxlm1(1:kproma,:)+pxim1(1:kproma,:))) !virt. Temp

    DO jk = 1,klev      
      zrho(1:kproma,jk) = papm1(1:kproma,jk)/(rd*ptvm1(1:kproma,jk)) !moist air density
      zrho_dry(1:kproma,jk) = (papm1(1:kproma,jk)*amd/(argas*ptm1(1:kproma,jk)))*10.**(-3.) !dry air density
      !WRITE(*,*) 'papm1', papm1(1:kproma,jk) 
      !WRITE(*,*) 'zrho', zrho(1:kproma,jk)
      !WRITE(*,*) 'zrho_dry', zrho_dry(1:kproma,jk)
    ENDDO !jk     

  
  END SUBROUTINE get_input_para  



END MODULE mo_box_input

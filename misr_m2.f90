program misr_m2
  use netcdf
  implicit none
  
  write(*,*) 'Opening I3RC output file...'
  
    status = nf90_open(netcdf_file,0,ncid)
    if(status .ne. nf90_noerr) stop 'File Error'
    status = nf90_inq_varid(ncid,'intensity',varid)
    if(status .ne. nf90_noerr) stop 'Variable not found'
	
    status = nf90_inquire_dimension(ncid,1,len = dimension(1))
    status = nf90_inquire_dimension(ncid,2,len = dimension(2))
    status = nf90_inquire_dimension(ncid,3,len = dimension(3))
	
    nx = dimension(2)
    ny = dimension(3)
	
!   Check for angular information
!   zeniths represented as cos_inverse(mu) and azimuths as phi
    
	status = nf90_inq_varid(ncid,'intensityMus',muid)
	if(status .ne. nf90_noerr) stop 'Zenith angles not found'
        status = nf90_inq_varid(ncid,'intensityPhis',phiid)
	if(status .ne. nf90_noerr) stop 'Azimuth angles not found'
	
!   Allocating memory for storing the angular information
!   length of angular information = dimension(3)

    allocate(mu(dimension(3)),phi(dimension(3)))
	
!   Reading the angular information off the I3RC output file

	status = nf90_get_var(ncid,muid,mu)
	if(status .ne. nf90_noerr) stop 'Zenith angles not found in file'
        status = nf90_get_var(ncid,phiid,phi)
	if(status .ne. nf90_noerr) stop 'Azimuth angles not found in file'
	
!   Now calculate reference zenith and comparison zeniths

    Ref_zenith = acos(180.0/pi*mu(Ref_File_ID))
    Comp_zenith_1 = acos(180.0/pi*mu(Comp_File_1_ID))
    Comp_zenith_2 = acos(180.0/pi*mu(Comp_File_2_ID))
    
!  Now calculate the sizes of the search area box
!  Refer to MISR Level 2 Cloud Product Algorithm Theoretical Basis (Page 17), JPL D-73327
!  Defining the minimum and maximum SOM x and y disparities
    
    dtan = tan(Comp_Zenith_1 * pi/180) - tan(Ref_zenith * pi/180)
    if (dtan.ge.0) then
    dxmin = int((hmin * dtan - abs(vmax * dt))/pixel_size)
    dxmax = hmax * dtan + abs(vmax * dt)
    else 
    dxmin = hmax * dtan - abs(vmax * dt)
    dxmax = hmin * dtan + abs(vmax * dt)
    end if
    dymin = -abs(vmax * dt)
    dymax = abs(vmax * dt)
   
!  Setting the number of steps in the along-track(y) and cross-track(x) directions
    La = int(dxmin)

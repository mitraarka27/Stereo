program misr_m2
   use netcdf
   implicit none
!--**************************************-Variable definitions-***********************************************************



!--**************************************-Reading from I3RC output-*******************************************************
    write(*,*) 'Opening I3RC output file...'
  
    status = nf90_open(netcdf_file,0,ncid)
    status = nf90_inq_varid(ncid,'intensity',varid)
    status = nf90_inquire_dimension(ncid,1,len = length(1))
    status = nf90_inquire_dimension(ncid,2,len = length(2))
    status = nf90_inquire_dimension(ncid,3,len = length(3))

!   Check for angular information
!   zeniths represented as cos_inverse(mu) and azimuths as phi
    
    status = nf90_inq_varid(ncid,'intensityMus',muid)
    status = nf90_inq_varid(ncid,'intensityPhis',phiid)
	

    nx = length(1)
    ny = length(2)

!   Allocating memory to store the required radiance data from netcdf file
!   Allocating memory for storing the angular information. Length of angular information = dimension(3)
    
    allocate(FileData(length(1),length(2),length(3)))
    allocate(RefImage(nx,ny),CompImage1(nx,ny),CompImage2(nx,ny))
    allocate(mu(length(3)),phi(length(3)))
    allocate(H(nx,ny),H1(nx,ny),H2(nx,ny))

!   Reading the radiance information from the output file

    status=nf90_get_var(ncid,varid,FileData) 
    RefImage = FileData(:,:,RefImageNum)
    CompImage1 = FileData(:,:,CompImageNum1)
    CompImage2 = FileData(:,:,CompImageNum2)
	
!   Reading the angular information off the I3RC output file

    status = nf90_get_var(ncid,muid,mu)
    status = nf90_get_var(ncid,phiid,phi)

!   Now calculate reference zenith and comparison zeniths

    Ref_zenith = acos(180.0/pi*mu(RefImageNum))
    Comp_zenith_1 = acos(180.0/pi*mu(CompImageNum1))
    Comp_zenith_2 = acos(180.0/pi*mu(CompImageNum2))

!  Now let us initialize the height vectors H1, H2 and the cumulative height vector H with garbage values(-666)
    H(:,:) = -666
    H1(:,:) = -666
    H2(:,:) = -666
   
!  Now calculate the sizes of the search area box for 1st comparison image
!  Refer to MISR Level 2 Cloud Product Algorithm Theoretical Basis (Page 17), JPL D-73327
!  Defining the minimum and maximum SOM x and y disparities
    
    dtan = tan(Comp_Zenith_1 * pi/180) - tan(Ref_zenith * pi/180)
    if (dtan.ge.0) then
    dxmin = int((hmin * dtan - abs(vmax * dt))/pixel_size)-1
    dxmax = int((hmax * dtan + abs(vmax * dt))/pixel_size)+1
    else 
    dxmin = int((hmax * dtan - abs(vmax * dt))/pixel_size)-1
    dxmax = int((hmin * dtan + abs(vmax * dt))/pixel_size)+1
    end if
    dymin = int(-abs(vmax * dt)/pixel_size)-1
    dymax = int(abs(vmax * dt)/pixel_size)+1
   
!  Setting the number of steps in the along-track(y) and cross-track(x) directions
    Nc = abs(dxmax-dxmin) + 1
    Na = abs(dcmax-dcmin) + 1
    max_patches = Nc * Na + 1
    
!  Now calculate cloud top heights using the stereo code for 1st Comparison Image
    
    H1 = stereo(RefImage(:,:), CompImage1(:,:), Ref_zenith, Comp_Zenith_1)
    
!  Now calculate the sizes of the search area box for 2nd comparison image

    dtan = tan(Comp_Zenith_2 * pi/180) - tan(Ref_zenith * pi/180)
    if (dtan.ge.0) then
    dxmin = int((hmin * dtan - abs(vmax * dt))/pixel_size)-1
    dxmax = int((hmax * dtan + abs(vmax * dt))/pixel_size)+1
    else 
    dxmin = int((hmax * dtan - abs(vmax * dt))/pixel_size)-1
    dxmax = int((hmin * dtan + abs(vmax * dt))/pixel_size)+1
    end if
    dymin = int(-abs(vmax * dt)/pixel_size)-1
    dymax = int(abs(vmax * dt)/pixel_size)+1
    
!  Setting the number of steps in the along-track(y) and cross-track(x) directions again 
    Nc = abs(dxmax-dxmin) + 1
    Na = abs(dcmax-dcmin) + 1
    max_patches = Nc * Na + 1
    
!  Now calculate cloud top heights using the stereo code for 2nd comparison image
    
    H1 = stereo(RefImage(:,:), CompImage2(:,:), Ref_zenith, Comp_Zenith_2)  

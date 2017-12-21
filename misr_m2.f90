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

!--**************************************-Stereo Height Calculations-*******************************************************

!   Now calculate reference zenith and comparison zeniths

    Ref_zenith = acos(180.0/pi*mu(RefImageNum))
    Comp_zenith_1 = acos(180.0/pi*mu(CompImageNum1))
    Comp_zenith_2 = acos(180.0/pi*mu(CompImageNum2))

!  Now let us initialize the height vectors H1, H2 and the cumulative height vector H with garbage values(-666)
    H(:,:) = -9999
    H1(:,:) = -9999
    H2(:,:) = -9999
   
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
    
    H2 = stereo(RefImage(:,:), CompImage2(:,:), Ref_zenith, Comp_Zenith_2)  
    
!  Computing the heights
    H = (H1+H2)/2
   
!  Write to output file   
   print*,'Writing heights to file...'
   open(10,file="StereoHeights", action='write')
   do no=1,ny
     write(10,*) Heights(:,no)
   end do
   close(10)
   print*,'Writing to file complete'
   
 !--**************************************-Modules/Subroutines/Functions Used-**********************************************
   contains
 
 !--**************************************-Along-track Disparity Calculation-***********************************************
    real function StereoHeight(c1,r1,c2,r2,v1,v2)
    
    real, intent(in) :: c2,r2,v1,v2
    real :: d
    integer,intent(in) :: c1,r1

    d=sqrt((r2-r1)**2.+(c2-c1)**2.)*pixelSize
    StereoHeight=d/abs(tan(v2*pi/180)-tan(v1*pi/180))
    end function StereoHeight
    
 !--**************************************-MISR M2 Stereo Matcher Code-****************************************************
    function stereo(I1,I2,v1,v2)
    real, dimension(nCols,nRows) :: m2matcher
    real, intent(in) :: v1,v2
    real, dimension(0:nCols-1,0:nRows-1), intent(in) :: I1,I2
    real, dimension(0:nCols-1,0:nRows-1) :: S,Sp,Sa
    real, dimension(0:nCols-1,0:nRows-1,2) :: disp
    real, dimension(2*mwinC+1,2*mwinR+1) :: targetArray, searchArray
    logical :: test
    integer, dimension(maxvecs,2) :: vecs
    integer :: m,n,x,y,xx,yy,vc,radj,i,j,sizeT,sizeS
    real :: pct, trange, srange, tbar, sbar, sigma, Rsigma
    real :: big
    big=10.**10
    pct=0.
    m2matcher(:,:)=-9999.
    do j=0,nRows-1
      do i=0,nCols-1
        disp(i,j,:)=(/-1000, -1000/)
      end do
    end do
        vc=1
!print*,dcmin,dcmax,drmin,drmax
    do x=min(dcmin, dcmax),max(dcmin, dcmax)
      xx=x
      do y=min(drmin, drmax),max(drmin, drmax)
        yy=y
        vecs(vc,1)=xx
        vecs(vc,2)=yy
        vc=vc+1
      end do
    end do
    vc=vc-1
!print*,vc
    S(:,:)=big
    Sa(:,:)=big

    do x=1,vc
      if (real(x)/real(vc) .ge. pct) then
        write(*,'(i3,a)') int(pct*100), '% Complete'
        pct=pct+0.25
      end if
      Sp(:,:)=big*2
      sizeT=(2*mwinR+1)*(2*mwinC+1)
      sizeS=(2*mwinR+1)*(2*mwinC+1)
      do n=0,nRows-1
        do m=0,nCols-1
          if ((m+vecs(x,1)-mwinC .ge. 0).and.(n+vecs(x,2)-mwinR .ge. 0) .and. &
                (m-mwinC .ge. 0).and.(n-mwinR .ge. 0) .and. &
                (m+mwinC .lt. nCols).and.(n+mwinR .lt. nRows) .and. &
              (m+vecs(x,1)+mwinC .lt. nCols).and.(n+vecs(x,2)+mwinR .lt. nRows)) then
              targetArray(:,:)=I1(m-mwinC:m+mwinC,n-mwinR:n+mwinR)
              searchArray(:,:)=I2(m+vecs(x,1)-mwinC:m+vecs(x,1)+mwinC,n+vecs(x,2)-mwinR:n+vecs(x,2)+mwinR)
              trange=maxval(targetArray)-minval(targetArray)
              !print*,'tr',trange
              srange=maxval(searchArray)-minval(searchArray)

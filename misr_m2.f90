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
	

!   Rows = y, Columns = x
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
    H(:,:) = -9999.
    H1(:,:) = -9999.
    H2(:,:) = -9999.
   
!  Now calculate the sizes of the search area box for 1st comparison image
!  Refer to MISR Level 2 Cloud Product Algorithm Theoretical Basis (Page 17), JPL D-73327
!  Defining the minimum and maximum SOM x and y disparities
    
    dtan = abs(tan(Comp_Zenith_1 * pi/180) - tan(Ref_zenith * pi/180))
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

    dtan = abs(tan(Comp_Zenith_2 * pi/180) - tan(Ref_zenith * pi/180))
    if (dtan.ge.0) then
    xmin = int((hmin * dtan - abs(vmax * dt))/pixel_size)-1
    xmax = int((hmax * dtan + abs(vmax * dt))/pixel_size)+1
    else 
    xmin = int((hmax * dtan - abs(vmax * dt))/pixel_size)-1
    xmax = int((hmin * dtan + abs(vmax * dt))/pixel_size)+1
    end if
    ymin = int(-abs(vmax * dt)/pixel_size)-1
    ymax = int(abs(vmax * dt)/pixel_size)+1
    
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
    real function StereoHeight(c1,r1,c2,r2,theta0,theta1)
    real, intent(in) :: x1,x2,y1,y2,theta0,theta1
    real :: d

    d=sqrt((y2-y1)**2.+(x2-x1)**2.)
    StereoHeight=d/abs(tan(theta1*pi/180)-tan(theta0*pi/180))*pixel_size
    end function StereoHeight
    
 !--**************************************-MISR M2 Stereo Matcher Code-****************************************************
    function stereo(I1,I2,v1,v2)
    real, dimension(nCols,nRows) :: m2matcher
    real, intent(in) :: v1,v2
    real, dimension(0:nx-1,0:ny-1), intent(in) :: I1,I2
    real, dimension(0:nx-1,0:ny-1) :: S,Sp,Sa
    real, dimension(0:nx-1,0:ny-1,2) :: disp
    real, dimension(5,5) :: targetPatch, searchPatch
    logical :: test
    real :: pct, trange, srange, tbar, sbar, sigma, Rsigma
    real :: big
    logical :: test
    integer, dimension(max_patches,2) :: vecs
    integer :: m,n,x,y,xx,yy,vc,radj,i,j,sizeT,sizeS
    
    big=10.**10
    pct=0.
    stereo(:,:)=-9999.
    
    do j=0,ny-1
      do i=0,nx-1
        disp(i,j,:)=(/-1000, -1000/)
      end do
    end do
    vc=1

    do x=min(xmin,xmax),max(xmin,xmax)
      xx=x
      do y=min(ymin, ymax),max(ymin, ymax)
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
    Sp(:,:)=big*2

    do x = 1,vc
     if (real(x)/real(vc) .ge. pct) then
        write(*,'(i3,a)') int(pct*100), '% Complete'
        pct=pct+0.25
     end if 
      
      
    sizeT=edge**2
    sizeS=edge**2
      
      do n=0,ny-1
        do m=0,nx-1
          if ((m+vecs(x,1)-edge .ge. 0).and.(n+vecs(x,2)-edge .ge. 0) .and. &
                (m-edge .ge. 0).and.(n-edge .ge. 0) .and. &
                (m+edge .lt. nx).and.(n+edge .lt. ny) .and. &
                (m+vecs(x,1)+edge .lt. nx).and.(n+vecs(x,2)+edge .lt. ny)) then
		
              targetPatch(:,:) = I1(m-edge:m+edge,n-edge:n+edge)
              searchPatch(:,:) = I2(m+vecs(x,1)-edge:m+vecs(x,1)+edge,n+vecs(x,2)-edge:n+vecs(x,2)+edge)
              dR=maxval(targetPatch) - minval(target)
              !print*,'dC = ',dC
              dC=maxval(searchPatch) - minval(searchPatch)
	      !print*,'dR = ',dR
              avgR = sum(targetPatch)/sizeT
              avgC = sum(searchPatch)/sizeS
              sigma = sum(abs((targetPatch - avgR)/dR))
	      
              if(contrastThreshold .eqv. .true.) then
	      
                if((sigma*SNR(m,n)*dR)/(dR*avgR) .ge. contrastThreshold) then
		
                  S_M2(m,n)=sum(abs((targetPatch-avgR)/dR-(searchPatch-avgC)/dC))/sigma
                else
		
                  S_M2(m,n)=big*2
                end if
		
              else
	      
                S_M2(m,n)=sum(abs((targetPatch-avgR)/dR-(searchPatch-avgC)/dC))/sigma
		  
              end if

          end if
        end do
      end do
!print*,minval(Sp),maxval(Sp)

      do n=0,nRows-1
        do m=0,nCols-1
          if (S(m,n) .gt. Sp(m,n)) then  !1.1*Sp? need else stmt to void match or second order of S
            Sa(m,n)=S(m,n)
            S(m,n)=Sp(m,n)
            if(metricThreshold) then
              if (S(m,n) .lt. threshold) then
                disp(m,n,:)=vecs(x,:)
              end if
            else
                disp(m,n,:)=vecs(x,:)
            end if
          else if (S(m,n) .eq. Sp(m,n)) then
            if ((disp(m,n,1)**2+disp(m,n,2)**2).gt.(vecs(x,1)**2+vecs(x,2)**2)) then
              Sa(m,n)=S(m,n)
              disp(m,n,:)=vecs(x,:)
            end if
          end if
        end do
      end do
    end do

    print*,'Finding Heights...'
    do n=0,nRows-1
      do m=0,nCols-1
        if (Heights(m+1,n+1).eq.-9999.) then
          if (disp(m,n,1) .ne. -1000) then
            m2matcher(m+1,n+1)=calcHeight(0,0,disp(m,n,1),disp(m,n,2),v1,v2)
          end if
        else
          !--if height already stored, don't overwrite it
          m2matcher(m+1,n+1)=Heights(m+1,n+1)
        end if
      end do
    end do

    if (ambiguity) then
      where(Sa<(1.1*S)) m2matcher=-9999.
    end if
  end function m2matcher


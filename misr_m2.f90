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
	
!   Check for angular information
    status = 
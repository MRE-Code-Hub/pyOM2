


module module_diag_overturning
!=======================================================================
!  diagnose meridional overturning on isopycnals and depth
!=======================================================================
 implicit none
 integer :: nitts
 integer :: nlevel 
 real*8,allocatable :: sig(:),zarea(:,:)
 character (len=80) :: over_file
 real*8 :: p_ref = 0d0 ! in dbar
 real*8,allocatable :: mean_trans(:,:)
 real*8,allocatable :: mean_vsf_iso(:,:)
 real*8,allocatable :: mean_bolus_iso(:,:)
 real*8,allocatable :: mean_vsf_depth(:,:)
 real*8,allocatable :: mean_bolus_depth(:,:)
 integer :: number_masks 
 real*8,allocatable :: mask1(:,:),mask2(:,:),zarea_1(:,:),zarea_2(:,:)
 real*8,allocatable :: mean_vsf_iso_1(:,:), mean_bolus_iso_1(:,:)
 real*8,allocatable :: mean_vsf_iso_2(:,:), mean_bolus_iso_2(:,:)
 
 real*8,allocatable :: mean_heat_tr(:),mean_salt_tr(:)
 real*8,allocatable :: mean_heat_tr_bolus(:),mean_salt_tr_bolus(:)
 real*8,allocatable :: mean_heat_tr_1(:),mean_salt_tr_1(:)
 real*8,allocatable :: mean_heat_tr_2(:),mean_salt_tr_2(:)
 real*8,allocatable :: mean_heat_tr_bolus_1(:),mean_salt_tr_bolus_1(:)
 real*8,allocatable :: mean_heat_tr_bolus_2(:),mean_salt_tr_bolus_2(:)
end module module_diag_overturning



subroutine init_diag_overturning
 use main_module
 use isoneutral_module
 use module_diag_overturning
 use rossmix2_module
 implicit none
 real*8 :: dsig,sigs,sige,get_rho
 integer :: i,j,k,n
 include "netcdf.inc"
 integer :: ncid,iret
 integer :: itimedim,itimeid,z_tdim,z_tid,z_udim,z_uid
 integer :: lat_udim,lat_uid,lat_tdim,lat_tid
 integer :: sig_dim,sig_id,id
 character (len=80) :: name,unit

 nitts = 0
 nlevel = nz*4
 over_file = 'over.cdf'

 allocate( sig(nlevel) );sig=0
 allocate( zarea(js_pe-onx:je_pe+onx,nz) ); zarea=0d0
 allocate( mean_trans(js_pe-onx:je_pe+onx,nlevel) );mean_trans=0
 allocate( mean_vsf_iso(js_pe-onx:je_pe+onx,nz) );mean_vsf_iso=0
 allocate( mean_vsf_depth(js_pe-onx:je_pe+onx,nz) );mean_vsf_depth=0
 allocate( mean_heat_tr(js_pe-onx:je_pe+onx) ); mean_heat_tr=0
 allocate( mean_salt_tr(js_pe-onx:je_pe+onx) ); mean_salt_tr=0
 allocate( mean_heat_tr_bolus(js_pe-onx:je_pe+onx) ); mean_heat_tr_bolus=0
 allocate( mean_salt_tr_bolus(js_pe-onx:je_pe+onx) ); mean_salt_tr_bolus=0
 allocate( mean_bolus_iso(js_pe-onx:je_pe+onx,nz) );mean_bolus_iso=0
 allocate( mean_bolus_depth(js_pe-onx:je_pe+onx,nz) );mean_bolus_depth=0

 
 if ( (enable_neutral_diffusion .and. enable_skew_diffusion) .and. enable_rossmix2) then
  if (my_pe==0) print*,'ERROR: skew diffusion and Rossmix2 is not possible'
  call halt_stop(' in init_diag_overturing ')
 endif
 

 ! sigma levels
 p_ref=2000.0
 sige = get_rho(35d0,-2d0,p_ref)
 sigs = get_rho(35d0,30d0,p_ref)
 dsig = (sige-sigs)/(nlevel-1.)
 if (my_pe==0) then
   print'(a)',      ' sigma ranges for overturning diagnostic:' 
   print'(a,f12.6)',' start sigma0 = ',sigs
   print'(a,f12.6)',' end sigma0   = ',sige
   print'(a,f12.6)',' Delta sigma0 = ',dsig
   if (enable_neutral_diffusion .and. enable_skew_diffusion) &
      print'(a)',      ' also calculating overturning by eddy-driven velocities' 
   if (enable_rossmix2) &
      print'(a)',      ' also calculating overturning by eddy-driven velocities from Rossmix2' 
 endif
 do k=1,nlevel
   sig(k) = sigs + dsig*(k-1)
 enddo
 
 ! precalculate area below z levels
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    zarea(j,:)=zarea(j,:)+dxt(i)*cosu(j)*maskV(i,j,:)
  enddo
 enddo
 do k=2,nz
  zarea(:,k) = zarea(:,k-1) + zarea(:,k)*dzt(k)
 enddo
 call zonal_sum_vec(zarea(js_pe:je_pe,:),nz*(je_pe-js_pe+1))

 ! read in land masks
 iret=nf_open('over_mask.nc',NF_NOWRITE,ncid)
 if (iret == 0) then
 
   if (my_pe==0) print*,' reading two regional masks from file over_mask.nc'
   allocate( mean_vsf_iso_1(js_pe-onx:je_pe+onx,nz) );mean_vsf_iso_1=0
   allocate( mean_vsf_iso_2(js_pe-onx:je_pe+onx,nz) );mean_vsf_iso_2=0
   allocate( mean_bolus_iso_1(js_pe-onx:je_pe+onx,nz) );mean_bolus_iso_1=0
   allocate( mean_bolus_iso_2(js_pe-onx:je_pe+onx,nz) );mean_bolus_iso_2=0
   allocate( mask1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );mask1=0
   allocate( mask2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) );mask2=0
   allocate( zarea_1(js_pe-onx:je_pe+onx,nz) ); zarea_1=0d0
   allocate( zarea_2(js_pe-onx:je_pe+onx,nz) ); zarea_2=0d0
   allocate( mean_heat_tr_1(js_pe-onx:je_pe+onx) ); mean_heat_tr_1=0
   allocate( mean_salt_tr_1(js_pe-onx:je_pe+onx) ); mean_salt_tr_1=0
   allocate( mean_heat_tr_2(js_pe-onx:je_pe+onx) ); mean_heat_tr_2=0
   allocate( mean_salt_tr_2(js_pe-onx:je_pe+onx) ); mean_salt_tr_2=0
   allocate( mean_heat_tr_bolus_1(js_pe-onx:je_pe+onx) ); mean_heat_tr_bolus_1=0
   allocate( mean_salt_tr_bolus_1(js_pe-onx:je_pe+onx) ); mean_salt_tr_bolus_1=0
   allocate( mean_heat_tr_bolus_2(js_pe-onx:je_pe+onx) ); mean_heat_tr_bolus_2=0
   allocate( mean_salt_tr_bolus_2(js_pe-onx:je_pe+onx) ); mean_salt_tr_bolus_2=0

   iret=nf_inq_varid(ncid,'mask1',id)
   iret= nf_get_vara_double(ncid,id,(/is_pe,js_pe/), (/ie_pe-is_pe+1,je_pe-js_pe+1/),mask1(is_pe:ie_pe,js_pe:je_pe))
   iret=nf_inq_varid(ncid,'mask2',id)
   iret= nf_get_vara_double(ncid,id,(/is_pe,js_pe/), (/ie_pe-is_pe+1,je_pe-js_pe+1/),mask2(is_pe:ie_pe,js_pe:je_pe))
   iret=nf_close(ncid)
   number_masks = 2
   call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,mask1) 
   call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,mask1)
   call border_exchg_xy(is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,mask2) 
   call setcyclic_xy   (is_pe-onx,ie_pe+onx,js_pe-onx,je_pe+onx,mask2)

   do j=js_pe,je_pe
    do i=is_pe,ie_pe
     zarea_1(j,:)=zarea_1(j,:)+dxt(i)*cosu(j)*maskV(i,j,:)*mask1(i,j)
     zarea_2(j,:)=zarea_2(j,:)+dxt(i)*cosu(j)*maskV(i,j,:)*mask2(i,j)
    enddo
   enddo
   do k=2,nz
    zarea_1(:,k) = zarea_1(:,k-1) + zarea_1(:,k)*dzt(k)
    zarea_2(:,k) = zarea_2(:,k-1) + zarea_2(:,k)*dzt(k)
   enddo
   call zonal_sum_vec(zarea_1(js_pe:je_pe,:),nz*(je_pe-js_pe+1))
   call zonal_sum_vec(zarea_2(js_pe:je_pe,:),nz*(je_pe-js_pe+1))   
 else
   if (my_pe==0) print*,'WARNING: cannot read file over_mask.nc'
   number_masks = 0
 endif

 ! prepare cdf file for output
 if (my_pe==0) then
    print'(2a)',' preparing file ',over_file(1:len_trim(over_file))
    iret = nf_create (over_file, nf_clobber, ncid)
    iret=nf_set_fill(ncid, NF_NOFILL, iret)
    iTimedim  = ncddef(ncid, 'Time', nf_unlimited, iret)
    itimeid  = ncvdef (ncid,'Time', NCFLOAT,1,(/itimedim/),iret)
    name = 'Time '; unit = 'days'
    call ncaptc(ncid, itimeid, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, itimeid, 'units',   NCCHAR, len_trim(unit), unit, iret) 
    call ncaptc(ncid, iTimeid,'time_origin',NCCHAR, 20,'01-JAN-1900 00:00:00', iret)
    sig_dim = ncddef(ncid, 'sigma', nlevel , iret)
    sig_id  = ncvdef (ncid,'sigma',NCFLOAT,1,(/sig_dim/),iret)
    name = 'Sigma axis'; unit = 'kg/m^3'
    call ncaptc(ncid, sig_id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, sig_id, 'units',     NCCHAR, len_trim(unit), unit, iret) 

    z_tdim    = ncddef(ncid, 'zt',  nz, iret)
    z_udim    = ncddef(ncid, 'zu',  nz, iret)
    z_tid  = ncvdef (ncid,'zt', NCFLOAT,1,(/z_tdim/),iret)
    z_uid  = ncvdef (ncid,'zu', NCFLOAT,1,(/z_udim/),iret)
    name = 'Height on T grid     '; unit = 'm'
    call ncaptc(ncid, z_tid, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, z_tid, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    name = 'Height on U grid     '; unit = 'm'
    call ncaptc(ncid, z_uid, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, z_uid, 'units',     NCCHAR, len_trim(unit), unit, iret) 

    Lat_udim  = ncddef(ncid,'yu', ny , iret)
    Lat_uid  = ncvdef (ncid,'yu',NCFLOAT,1,(/lat_udim/),iret)
    Lat_tdim  = ncddef(ncid,'yt', ny , iret)
    Lat_tid  = ncvdef (ncid,'yt',NCFLOAT,1,(/lat_udim/),iret)
    if (coord_degree) then
       name = 'Latitude on T grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       name = 'Latitude on U grid     '; unit = 'degrees N'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    else
       name = 'meridional axis T grid'; unit = 'km'
       call ncaptc(ncid, Lat_tid, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, Lat_tid, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       name = 'meridional axis U grid'; unit = 'km'
       call ncaptc(ncid, Lat_uid, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, Lat_uid, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    endif

    id  = ncvdef (ncid,'trans',NCFLOAT,3,(/lat_udim,sig_dim,itimedim/),iret)
    name = 'Meridional transport'; unit = 'm^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'vsf_iso',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
    name = 'Meridional transport'; unit = 'm^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'vsf_depth',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
    name = 'Meridional transport'; unit = 'm^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'heat_tr',NCFLOAT,2,(/lat_udim,itimedim/),iret)
    name = 'Meridional heat transport'; unit = 'deg C m^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'salt_tr',NCFLOAT,2,(/lat_udim,itimedim/),iret)
    name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'heat_tr_bolus',NCFLOAT,2,(/lat_udim,itimedim/),iret)
    name = 'Meridional heat transport'; unit = 'deg C m^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

    id  = ncvdef (ncid,'salt_tr_bolus',NCFLOAT,2,(/lat_udim,itimedim/),iret)
    name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
    call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
    call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
    call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
    call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)


    if (number_masks >0) then
     id  = ncvdef (ncid,'vsf_mask1',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
     name = 'Meridional transport'; unit = 'm^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)
     
     id  = ncvdef (ncid,'vsf_mask2',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
     name = 'Meridional transport'; unit = 'm^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)   
     
     id  = ncvdef (ncid,'heat_tr_1',NCFLOAT,2,(/lat_udim,itimedim/),iret)
     name = 'Meridional heat transport'; unit = 'deg C m^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

     id  = ncvdef (ncid,'heat_tr_2',NCFLOAT,2,(/lat_udim,itimedim/),iret)
     name = 'Meridional heat transport'; unit = 'deg C m^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

     id  = ncvdef (ncid,'salt_tr_1',NCFLOAT,2,(/lat_udim,itimedim/),iret)
     name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

     id  = ncvdef (ncid,'salt_tr_2',NCFLOAT,2,(/lat_udim,itimedim/),iret)
     name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
     call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
     call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
     call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
     call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)
 
    endif

    if ((enable_neutral_diffusion .and. enable_skew_diffusion).or.enable_rossmix2) then
      id  = ncvdef (ncid,'bolus_iso',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
      name = 'Meridional transport'; unit = 'm^3/s'
      call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
      call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
      call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
      call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

      id  = ncvdef (ncid,'bolus_depth',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
      name = 'Meridional transport'; unit = 'm^3/s'
      call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
      call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
      call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
      call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret) 
      
      if (number_masks >0) then
       id  = ncvdef (ncid,'bolus_mask1',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
       name = 'Meridional transport'; unit = 'm^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)
     
       id  = ncvdef (ncid,'bolus_mask2',NCFLOAT,3,(/lat_udim,z_udim,itimedim/),iret)
       name = 'Meridional transport'; unit = 'm^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)  
       
       id  = ncvdef (ncid,'heat_tr_bolus_1',NCFLOAT,2,(/lat_udim,itimedim/),iret)
       name = 'Meridional heat transport'; unit = 'deg C m^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

       id  = ncvdef (ncid,'heat_tr_bolus_2',NCFLOAT,2,(/lat_udim,itimedim/),iret)
       name = 'Meridional heat transport'; unit = 'deg C m^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

       id  = ncvdef (ncid,'salt_tr_bolus_1',NCFLOAT,2,(/lat_udim,itimedim/),iret)
       name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)

       id  = ncvdef (ncid,'salt_tr_bolus_2',NCFLOAT,2,(/lat_udim,itimedim/),iret)
       name = 'Meridional salt transport'; unit = 'g/kg m^3/s'
       call ncaptc(ncid, id, 'long_name', NCCHAR, len_trim(name), name, iret) 
       call ncaptc(ncid, id, 'units',     NCCHAR, len_trim(unit), unit, iret) 
       call ncapt (ncid,id, 'missing_value',NCFLOAT,1,-1e33,iret)
       call ncapt (ncid,id, '_FillValue', NCFLOAT, 1,-1e33, iret)
        
      endif
            
    endif

    call ncendf(ncid, iret)
    iret= nf_put_vara_double(ncid,z_tid,(/1/),(/nz/),zt)
    iret= nf_put_vara_double(ncid,z_uid,(/1/),(/nz/),zw)
    iret= nf_put_vara_double(ncid,sig_id,(/1/),(/nlevel/),sig)
    iret=nf_close(ncid)
 endif


 do n=0,n_pes-1
   call fortran_barrier
   if (my_pe==n) then
     iret=nf_open(over_file,NF_WRITE,ncid)
     iret=nf_inq_varid(ncid,'yt',lat_tid)
     iret=nf_inq_varid(ncid,'yu',lat_uid)
     if (coord_degree) then
       iret= nf_put_vara_double(ncid,lat_Tid,(/js_pe/),(/je_pe-js_pe+1/) ,yt(js_pe:je_pe))
       iret= nf_put_vara_double(ncid,lat_uid,(/js_pe/),(/je_pe-js_pe+1/) ,yu(js_pe:je_pe))
     else
       iret= nf_put_vara_double(ncid,lat_Tid,(/js_pe/),(/je_pe-js_pe+1/) ,yt(js_pe:je_pe)/1e3)
       iret= nf_put_vara_double(ncid,lat_uid,(/js_pe/),(/je_pe-js_pe+1/) ,yu(js_pe:je_pe)/1e3)
     endif
     iret=nf_close(ncid)
   endif
 enddo
end subroutine init_diag_overturning



subroutine diag_overturning
 use main_module
 use isoneutral_module
 use module_diag_overturning
 use rossmix2_module
 implicit none
 integer :: i,j,k,m,m1,m2,mm,mmp,mmm
 real*8 :: get_rho
 real*8 :: trans(js_pe-onx:je_pe+onx,nlevel),fxa
 real*8 :: z_sig(js_pe-onx:je_pe+onx,nlevel) 
 real*8 :: bolus_trans(js_pe-onx:je_pe+onx,nlevel)
 real*8 :: bolus_iso(js_pe-onx:je_pe+onx,nz)
 real*8 :: vsf_iso(js_pe-onx:je_pe+onx,nz)
 real*8 :: vsf_depth(js_pe-onx:je_pe+onx,nz)
 real*8 :: bolus_depth(js_pe-onx:je_pe+onx,nz)
 real*8 :: sig_loc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,nz)
 real*8 :: heat_tr(js_pe-onx:je_pe+onx),salt_tr(js_pe-onx:je_pe+onx)
 real*8 :: heat_tr_bolus(js_pe-onx:je_pe+onx),salt_tr_bolus(js_pe-onx:je_pe+onx)
 
 real*8 :: trans_1(js_pe-onx:je_pe+onx,nlevel)
 real*8 :: z_sig_1(js_pe-onx:je_pe+onx,nlevel) 
 real*8 :: trans_2(js_pe-onx:je_pe+onx,nlevel)
 real*8 :: z_sig_2(js_pe-onx:je_pe+onx,nlevel) 
 real*8 :: bolus_trans_1(js_pe-onx:je_pe+onx,nlevel)
 real*8 :: bolus_trans_2(js_pe-onx:je_pe+onx,nlevel)
 real*8 :: vsf_iso_1(js_pe-onx:je_pe+onx,nz)
 real*8 :: vsf_iso_2(js_pe-onx:je_pe+onx,nz)
 real*8 :: bolus_iso_1(js_pe-onx:je_pe+onx,nz)
 real*8 :: bolus_iso_2(js_pe-onx:je_pe+onx,nz)
 real*8 :: heat_tr_1(js_pe-onx:je_pe+onx),salt_tr_1(js_pe-onx:je_pe+onx)
 real*8 :: heat_tr_2(js_pe-onx:je_pe+onx),salt_tr_2(js_pe-onx:je_pe+onx)
 real*8 :: heat_tr_bolus_1(js_pe-onx:je_pe+onx),salt_tr_bolus_1(js_pe-onx:je_pe+onx)
 real*8 :: heat_tr_bolus_2(js_pe-onx:je_pe+onx),salt_tr_bolus_2(js_pe-onx:je_pe+onx)
 
 ! sigma at p_ref
 do k=1,nz
  do j=js_pe,je_pe+1
    do i=is_pe,ie_pe
      sig_loc(i,j,k) =  get_rho(salt(i,j,k,tau),temp(i,j,k,tau),p_ref)
    enddo
  enddo
 enddo

 ! transports below isopycnals and area below isopycnals
 trans=0d0; z_sig=0d0
 do j=js_pe,je_pe
  do m=1,nlevel
   do k=1,nz
     do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) ) then 
         trans(j,m) = trans(j,m) + v(i,j,k,tau)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)
         z_sig(j,m) = z_sig(j,m) + dzt(k)*dxt(i)*cosu(j)*maskV(i,j,k)
       endif
     enddo
   enddo
  enddo 
 enddo
 call zonal_sum_vec(trans(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
 call zonal_sum_vec(z_sig(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))

 if (number_masks >0) then
  trans_1=0d0; z_sig_1=0d0
  trans_2=0d0; z_sig_2=0d0
  do j=js_pe,je_pe
   do m=1,nlevel
    do k=1,nz
     do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) ) then 
         trans_1(j,m) = trans_1(j,m) + v(i,j,k,tau)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask1(i,j)
         z_sig_1(j,m) = z_sig_1(j,m) + dzt(k)*dxt(i)*cosu(j)*maskV(i,j,k)*mask1(i,j)
         trans_2(j,m) = trans_2(j,m) + v(i,j,k,tau)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask2(i,j)
         z_sig_2(j,m) = z_sig_2(j,m) + dzt(k)*dxt(i)*cosu(j)*maskV(i,j,k)*mask2(i,j)
       endif
     enddo
    enddo
   enddo 
  enddo  
  call zonal_sum_vec(trans_1(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
  call zonal_sum_vec(z_sig_1(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
  call zonal_sum_vec(trans_2(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
  call zonal_sum_vec(z_sig_2(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
 endif
 

 if (enable_neutral_diffusion .and. enable_skew_diffusion) then
   ! eddy driven transports below isopycnals
   bolus_trans=0d0;
   do j=js_pe,je_pe
    do m=1,nlevel
     k=1
     do i=is_pe,ie_pe
      fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
      if (fxa > sig(m) ) then
        bolus_trans(j,m) = bolus_trans(j,m) + b1_gm(i,j,k)*dxt(i)*cosu(j)*maskV(i,j,k)
      endif  
     enddo
     do k=2,nz
      do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) ) then
          bolus_trans(j,m) = bolus_trans(j,m) + (b1_gm(i,j,k)-b1_gm(i,j,k-1))*dxt(i)*cosu(j)*maskV(i,j,k)
       endif   
      enddo
     enddo
    enddo 
   enddo
   call zonal_sum_vec(bolus_trans(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
      
   if (number_masks >0) then
    bolus_trans_1=0d0; bolus_trans_2=0d0   
    do j=js_pe,je_pe
     do m=1,nlevel
      k=1
      do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) ) then
        bolus_trans_1(j,m) = bolus_trans_1(j,m) + b1_gm(i,j,k)*dxt(i)*cosu(j)*maskV(i,j,k)*mask1(i,j)
        bolus_trans_2(j,m) = bolus_trans_2(j,m) + b1_gm(i,j,k)*dxt(i)*cosu(j)*maskV(i,j,k)*mask2(i,j)
       endif  
      enddo
      do k=2,nz
       do i=is_pe,ie_pe
        fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
        if (fxa > sig(m) ) then
          bolus_trans_1(j,m) = bolus_trans_1(j,m) + (b1_gm(i,j,k)-b1_gm(i,j,k-1))*dxt(i)*cosu(j)*maskV(i,j,k)*mask1(i,j)
          bolus_trans_2(j,m) = bolus_trans_2(j,m) + (b1_gm(i,j,k)-b1_gm(i,j,k-1))*dxt(i)*cosu(j)*maskV(i,j,k)*mask2(i,j)
        endif   
       enddo
      enddo
     enddo 
    enddo  
    call zonal_sum_vec(bolus_trans_1(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
    call zonal_sum_vec(bolus_trans_2(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
   endif 
   
 endif
 
 if (enable_rossmix2) then
  ! eddy driven transports below isopycnals by Rossmix2
  bolus_trans=0d0;
  do j=js_pe,je_pe
   do m=1,nlevel
    do k=1,nz
      do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) )  then
         bolus_trans(j,m) = bolus_trans(j,m) + ve(i,j,k)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)
       endif  
      enddo
    enddo
   enddo 
  enddo
  call zonal_sum_vec(bolus_trans(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
  
  if (number_masks >0) then
   bolus_trans_1=0d0; bolus_trans_2=0d0
   do j=js_pe,je_pe
    do m=1,nlevel
     do k=1,nz
      do i=is_pe,ie_pe
       fxa = 0.5*( sig_loc(i,j,k) + sig_loc(i,j+1,k))
       if (fxa > sig(m) )  then
         bolus_trans_1(j,m) = bolus_trans_1(j,m) + ve(i,j,k)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask1(i,j)
         bolus_trans_2(j,m) = bolus_trans_2(j,m) + ve(i,j,k)*dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask2(i,j)
       endif  
      enddo
     enddo
    enddo 
   enddo
   call zonal_sum_vec(bolus_trans_1(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
   call zonal_sum_vec(bolus_trans_2(js_pe:je_pe,:),nlevel*(je_pe-js_pe+1))
  endif
  
 endif 
 

 ! streamfunction on geopotentials
 vsf_depth = 0d0
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    vsf_depth(j,:) = vsf_depth(j,:) + dxt(i)*cosu(j)*v(i,j,:,tau)*maskV(i,j,:)
   enddo
   do k=2,nz
     vsf_depth(j,k) = vsf_depth(j,k-1) + vsf_depth(j,k)*dzt(k)
   enddo     
 enddo     
 call zonal_sum_vec(vsf_depth(js_pe:je_pe,:),nz*(je_pe-js_pe+1))

 if (enable_neutral_diffusion .and. enable_skew_diffusion) then
   ! streamfunction for eddy driven velocity on geopotentials
   bolus_depth = 0d0
   do j=js_pe,je_pe
    do i=is_pe,ie_pe
      bolus_depth(j,:) = bolus_depth(j,:) + dxt(i)*cosu(j)*b1_gm(i,j,:)
    enddo
   enddo
   call zonal_sum_vec(bolus_depth(js_pe:je_pe,:),nz*(je_pe-js_pe+1))
 endif

 if (enable_rossmix2) then
  ! streamfunction for eddy driven velocity by Rossmix2 on geopotentials
  bolus_depth = 0d0  
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
    bolus_depth(j,:) = bolus_depth(j,:) + dxt(i)*cosu(j)*ve(i,j,:)*maskV(i,j,:)
   enddo
   do k=2,nz
     bolus_depth(j,k) = bolus_depth(j,k-1) + bolus_depth(j,k)*dzt(k)
   enddo     
  enddo  
  call zonal_sum_vec(bolus_depth(js_pe:je_pe,:),nz*(je_pe-js_pe+1)) 
 endif
 
 ! interpolate from isopcnals to depth
 if (my_blk_i==1) then
 
  vsf_iso = 0d0
  do j=js_pe,je_pe
   do k=1,nz
     mm= minloc( (zarea(j,k)-z_sig(j,:))**2,1 )
     mmp = min(mm+1,nlevel)
     mmm = max(mm-1,1)
     if     (z_sig(j,mm)>zarea(j,k) .and. z_sig(j,mmm) <= zarea(j,k) ) then
         m1=mmm; m2=mm
     elseif (z_sig(j,mm)>zarea(j,k) .and. z_sig(j,mmp) <= zarea(j,k) ) then
         m1=mmp; m2=mm
     elseif  (z_sig(j,mm)<zarea(j,k) .and. z_sig(j,mmp) >= zarea(j,k) ) then
         m1=mm; m2=mmp
     elseif  (z_sig(j,mm)<zarea(j,k) .and. z_sig(j,mmm) >= zarea(j,k) ) then
         m1=mm; m2=mmm
     else 
         m1=mm;m2=mm
     endif

     fxa =  z_sig(j,m2)-z_sig(j,m1)
     if (fxa /=0d0) then
      if (zarea(j,k)-z_sig(j,m1) > z_sig(j,m2)-zarea(j,k) ) then
       fxa = (zarea(j,k)-z_sig(j,m1))/fxa
       vsf_iso(j,k)=trans(j,m1)*(1-fxa) + trans(j,m2)*fxa 
       bolus_iso(j,k)=bolus_trans(j,m1)*(1-fxa) + bolus_trans(j,m2)*fxa  ! to save time
      else
       fxa = (z_sig(j,m2)-zarea(j,k))/fxa
       vsf_iso(j,k)=trans(j,m1)*fxa + trans(j,m2)*(1-fxa) 
       bolus_iso(j,k)=bolus_trans(j,m1)*fxa + bolus_trans(j,m2)*(1-fxa)
      endif
     else
      vsf_iso(j,k)=trans(j,m1) 
      bolus_iso(j,k)=bolus_trans(j,m1) 
     endif
   enddo
  enddo
  
  if (number_masks >0) then  
   vsf_iso_1 = 0d0
   do j=js_pe,je_pe
    do k=1,nz
     mm= minloc( (zarea_1(j,k)-z_sig_1(j,:))**2,1 )
     mmp = min(mm+1,nlevel)
     mmm = max(mm-1,1)
     if     (z_sig_1(j,mm)>zarea_1(j,k) .and. z_sig_1(j,mmm) <= zarea_1(j,k) ) then
         m1=mmm; m2=mm
     elseif (z_sig_1(j,mm)>zarea_1(j,k) .and. z_sig_1(j,mmp) <= zarea_1(j,k) ) then
         m1=mmp; m2=mm
     elseif  (z_sig_1(j,mm)<zarea_1(j,k) .and. z_sig_1(j,mmp) >= zarea_1(j,k) ) then
         m1=mm; m2=mmp
     elseif  (z_sig_1(j,mm)<zarea_1(j,k) .and. z_sig_1(j,mmm) >= zarea_1(j,k) ) then
         m1=mm; m2=mmm
     else 
         m1=mm;m2=mm
     endif

     fxa =  z_sig_1(j,m2)-z_sig_1(j,m1)
     if (fxa /=0d0) then
      if (zarea_1(j,k)-z_sig_1(j,m1) > z_sig_1(j,m2)-zarea_1(j,k) ) then
       fxa = (zarea_1(j,k)-z_sig_1(j,m1))/fxa
       vsf_iso_1(j,k)=trans_1(j,m1)*(1-fxa) + trans_1(j,m2)*fxa 
       bolus_iso_1(j,k)=bolus_trans_1(j,m1)*(1-fxa) + bolus_trans_1(j,m2)*fxa  ! to save time
      else
       fxa = (z_sig_1(j,m2)-zarea_1(j,k))/fxa
       vsf_iso_1(j,k)=trans_1(j,m1)*fxa + trans_1(j,m2)*(1-fxa) 
       bolus_iso_1(j,k)=bolus_trans_1(j,m1)*fxa + bolus_trans_1(j,m2)*(1-fxa)
      endif
     else
      vsf_iso_1(j,k)=trans_1(j,m1) 
      bolus_iso_1(j,k)=bolus_trans_1(j,m1) 
     endif
    enddo
   enddo 
 
   vsf_iso_2 = 0d0
   do j=js_pe,je_pe
    do k=1,nz
     mm= minloc( (zarea_2(j,k)-z_sig_2(j,:))**2,1 )
     mmp = min(mm+1,nlevel)
     mmm = max(mm-1,1)
     if     (z_sig_2(j,mm)>zarea_2(j,k) .and. z_sig_2(j,mmm) <= zarea_2(j,k) ) then
         m1=mmm; m2=mm
     elseif (z_sig_2(j,mm)>zarea_2(j,k) .and. z_sig_2(j,mmp) <= zarea_2(j,k) ) then
         m1=mmp; m2=mm
     elseif  (z_sig_2(j,mm)<zarea_2(j,k) .and. z_sig_2(j,mmp) >= zarea_2(j,k) ) then
         m1=mm; m2=mmp
     elseif  (z_sig_2(j,mm)<zarea_2(j,k) .and. z_sig_2(j,mmm) >= zarea_2(j,k) ) then
         m1=mm; m2=mmm
     else 
         m1=mm;m2=mm
     endif

     fxa =  z_sig_2(j,m2)-z_sig_2(j,m1)
     if (fxa /=0d0) then
      if (zarea_2(j,k)-z_sig_2(j,m1) > z_sig_2(j,m2)-zarea_2(j,k) ) then
       fxa = (zarea_2(j,k)-z_sig_2(j,m1))/fxa
       vsf_iso_2(j,k)=trans_2(j,m1)*(1-fxa) + trans_2(j,m2)*fxa 
       bolus_iso_2(j,k)=bolus_trans_2(j,m1)*(1-fxa) + bolus_trans_2(j,m2)*fxa  ! to save time
      else
       fxa = (z_sig_2(j,m2)-zarea_2(j,k))/fxa
       vsf_iso_2(j,k)=trans_2(j,m1)*fxa + trans_2(j,m2)*(1-fxa) 
       bolus_iso_2(j,k)=bolus_trans_2(j,m1)*fxa + bolus_trans_2(j,m2)*(1-fxa)
      endif
     else
      vsf_iso_2(j,k)=trans_2(j,m1) 
      bolus_iso_2(j,k)=bolus_trans_2(j,m1) 
     endif
    enddo
   enddo
  endif
 
 endif ! (my_blk_i==1)


! meridional heat and salt transport
 heat_tr = 0d0;salt_tr=0d0
 do k=1,nz
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
    fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*v(i,j,k,tau)*0.5
    heat_tr(j) = heat_tr(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
    salt_tr(j) = salt_tr(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
   enddo
  enddo
 enddo
 call zonal_sum_vec(heat_tr(js_pe:je_pe),je_pe-js_pe+1)
 call zonal_sum_vec(salt_tr(js_pe:je_pe),je_pe-js_pe+1)
 
 heat_tr_bolus = 0d0;salt_tr_bolus=0d0
 if (enable_neutral_diffusion .and. enable_skew_diffusion) then
  k=1
  do j=js_pe,je_pe 
   do i=is_pe,ie_pe
    fxa = dxt(i)*cosu(j)*maskV(i,j,k)*B1_gm(i,j,k)*0.5
    heat_tr_bolus(j) = heat_tr_bolus(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
    salt_tr_bolus(j) = salt_tr_bolus(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
   enddo
  enddo
  do k=2,nz
   do j=js_pe,je_pe 
    do i=is_pe,ie_pe
     fxa = dxt(i)*cosu(j)*maskV(i,j,k)*(B1_gm(i,j,k)-B1_gm(i,j,k-1))*0.5
     heat_tr_bolus(j) = heat_tr_bolus(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_bolus(j) = salt_tr_bolus(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
    enddo
   enddo
  enddo
 endif
 
 if (enable_rossmix2) then
  do k=1,nz
   do j=js_pe,je_pe 
    do i=is_pe,ie_pe
     fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*ve(i,j,k)*0.5
     heat_tr_bolus(j) = heat_tr_bolus(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_bolus(j) = salt_tr_bolus(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
    enddo
   enddo
  enddo
 endif 
 
 call zonal_sum_vec(heat_tr_bolus(js_pe:je_pe),je_pe-js_pe+1)
 call zonal_sum_vec(salt_tr_bolus(js_pe:je_pe),je_pe-js_pe+1)
 
 ! meridional heat and salt transport in masked regions
 if (number_masks>0) then
  heat_tr_1=0d0;heat_tr_2=0d0;salt_tr_1=0d0;salt_tr_2=0d0
  do k=1,nz
   do j=js_pe,je_pe 
    do i=is_pe,ie_pe
     fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask1(i,j)*v(i,j,k,tau)*0.5
     heat_tr_1(j) = heat_tr_1(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_1(j) = salt_tr_1(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
     fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask2(i,j)*v(i,j,k,tau)*0.5
     heat_tr_2(j) = heat_tr_2(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_2(j) = salt_tr_2(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
    enddo
   enddo
  enddo
  call zonal_sum_vec(heat_tr_1(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(salt_tr_1(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(heat_tr_2(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(salt_tr_2(js_pe:je_pe),je_pe-js_pe+1)
 
  heat_tr_bolus_1=0d0;heat_tr_bolus_2=0d0;salt_tr_bolus_1=0d0;salt_tr_bolus_2=0d0
  if (enable_neutral_diffusion .and. enable_skew_diffusion) then
   k=1
   do j=js_pe,je_pe 
    do i=is_pe,ie_pe
     fxa = dxt(i)*cosu(j)*maskV(i,j,k)*mask1(i,j)*B1_gm(i,j,k)*0.5
     heat_tr_bolus_1(j) = heat_tr_bolus_1(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_bolus_1(j) = salt_tr_bolus_1(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
     fxa = dxt(i)*cosu(j)*maskV(i,j,k)*mask2(i,j)*B1_gm(i,j,k)*0.5
     heat_tr_bolus_2(j) = heat_tr_bolus_2(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
     salt_tr_bolus_2(j) = salt_tr_bolus_2(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
    enddo
   enddo
   do k=2,nz
    do j=js_pe,je_pe 
     do i=is_pe,ie_pe
      fxa = dxt(i)*cosu(j)*maskV(i,j,k)*(B1_gm(i,j,k)-B1_gm(i,j,k-1))*mask1(i,j)*0.5
      heat_tr_bolus_1(j) = heat_tr_bolus_1(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
      salt_tr_bolus_1(j) = salt_tr_bolus_1(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
      fxa = dxt(i)*cosu(j)*maskV(i,j,k)*(B1_gm(i,j,k)-B1_gm(i,j,k-1))*mask2(i,j)*0.5
      heat_tr_bolus_2(j) = heat_tr_bolus_2(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
      salt_tr_bolus_2(j) = salt_tr_bolus_2(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
     enddo
    enddo
   enddo
  endif 
  
  if (enable_rossmix2) then
   do k=1,nz
    do j=js_pe,je_pe 
     do i=is_pe,ie_pe
      fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask1(i,j)*ve(i,j,k)*0.5
      heat_tr_bolus_1(j) = heat_tr_bolus_1(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
      salt_tr_bolus_1(j) = salt_tr_bolus_1(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
      fxa = dxt(i)*cosu(j)*dzt(k)*maskV(i,j,k)*mask2(i,j)*ve(i,j,k)*0.5
      heat_tr_bolus_2(j) = heat_tr_bolus_2(j) + (temp(i,j,k,tau) + temp(i,j+1,k,tau))*fxa
      salt_tr_bolus_2(j) = salt_tr_bolus_2(j) + (salt(i,j,k,tau) + salt(i,j+1,k,tau))*fxa
     enddo
    enddo
   enddo
  endif   
  
  call zonal_sum_vec(heat_tr_bolus_1(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(salt_tr_bolus_1(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(heat_tr_bolus_2(js_pe:je_pe),je_pe-js_pe+1)
  call zonal_sum_vec(salt_tr_bolus_2(js_pe:je_pe),je_pe-js_pe+1)

 endif
 
 ! average in time
 nitts = nitts + 1
 mean_trans = mean_trans + trans
 mean_vsf_iso = mean_vsf_iso + vsf_iso
 mean_vsf_depth = mean_vsf_depth + vsf_depth
 mean_heat_tr = mean_heat_tr + heat_tr
 mean_salt_tr = mean_salt_tr + salt_tr
 if (number_masks >0) then
  mean_vsf_iso_1 = mean_vsf_iso_1 + vsf_iso_1
  mean_vsf_iso_2 = mean_vsf_iso_2 + vsf_iso_2
  mean_heat_tr_1 = mean_heat_tr_1 + heat_tr_1
  mean_salt_tr_1 = mean_salt_tr_1 + salt_tr_1
  mean_heat_tr_2 = mean_heat_tr_2 + heat_tr_2
  mean_salt_tr_2 = mean_salt_tr_2 + salt_tr_2
 endif
 if ((enable_neutral_diffusion .and. enable_skew_diffusion).or.enable_rossmix2) then
  mean_bolus_iso = mean_bolus_iso + bolus_iso
  mean_bolus_depth = mean_bolus_depth + bolus_depth
  mean_heat_tr_bolus = mean_heat_tr_bolus + heat_tr_bolus
  mean_salt_tr_bolus = mean_salt_tr_bolus + salt_tr_bolus
  if (number_masks >0) then
   mean_bolus_iso_1 = mean_bolus_iso_1 + bolus_iso_1
   mean_bolus_iso_2 = mean_bolus_iso_2 + bolus_iso_2
   mean_heat_tr_bolus_1 = mean_heat_tr_bolus_1 + heat_tr_bolus_1
   mean_salt_tr_bolus_1 = mean_salt_tr_bolus_1 + salt_tr_bolus_1
   mean_heat_tr_bolus_2 = mean_heat_tr_bolus_2 + heat_tr_bolus_2
   mean_salt_tr_bolus_2 = mean_salt_tr_bolus_2 + salt_tr_bolus_2
  endif
 endif
end subroutine diag_overturning


subroutine write_overturning
 use main_module
 use isoneutral_module
 use module_diag_overturning
 use rossmix2_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,n
 integer :: itdimid,ilen,itimeid,id
 real*8 :: fxa

 if (my_pe==0) then
   print'(a,a)',' writing overturning diagnostics to file ',over_file(1:len_trim(over_file))
   iret=nf_open(over_file,NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   ilen=ilen+1
   fxa = itt*dt_tracer/86400.0
   iret=nf_inq_varid(ncid,'Time',itimeid)
   iret= nf_put_vara_double(ncid,itimeid,(/ilen/),(/1/),(/fxa/))
   iret=nf_close(ncid)
 endif

 if (nitts/=0) then
  mean_trans = mean_trans /nitts
  mean_vsf_iso = mean_vsf_iso /nitts
  mean_vsf_depth = mean_vsf_depth /nitts
  mean_heat_tr = mean_heat_tr /nitts
  mean_salt_tr = mean_salt_tr /nitts
  if (number_masks >0) then
   mean_vsf_iso_1 = mean_vsf_iso_1 /nitts
   mean_vsf_iso_2 = mean_vsf_iso_2 /nitts
   mean_heat_tr_1 = mean_heat_tr_1 /nitts
   mean_salt_tr_1 = mean_salt_tr_1 /nitts
   mean_heat_tr_2 = mean_heat_tr_2 /nitts
   mean_salt_tr_2 = mean_salt_tr_2 /nitts
  endif
  if ((enable_neutral_diffusion .and. enable_skew_diffusion).or.enable_rossmix2) then
   mean_bolus_iso = mean_bolus_iso /nitts
   mean_bolus_depth = mean_bolus_depth /nitts
   mean_heat_tr_bolus = mean_heat_tr_bolus /nitts
   mean_salt_tr_bolus = mean_salt_tr_bolus /nitts
   if (number_masks >0) then
    mean_bolus_iso_1 = mean_bolus_iso_1 /nitts
    mean_bolus_iso_2 = mean_bolus_iso_2 /nitts
    mean_heat_tr_bolus_1 = mean_heat_tr_bolus_1 /nitts
    mean_salt_tr_bolus_1 = mean_salt_tr_bolus_1 /nitts
    mean_heat_tr_bolus_2 = mean_heat_tr_bolus_2 /nitts
    mean_salt_tr_bolus_2 = mean_salt_tr_bolus_2 /nitts
   endif
  endif
 endif

 do n=1,n_pes_j
  call fortran_barrier
  if (my_blk_j==n .and. my_blk_i==1) then
   iret=nf_open(over_file,NF_WRITE,ncid)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid,itdimid,ilen)
   iret=nf_inq_varid(ncid,'trans',id)
   iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nlevel,1/),mean_trans(js_pe:je_pe,:))
   iret=nf_inq_varid(ncid,'vsf_iso',id)
   iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_vsf_iso(js_pe:je_pe,:))
   iret=nf_inq_varid(ncid,'vsf_depth',id)
   iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_vsf_depth(js_pe:je_pe,:))
   iret=nf_inq_varid(ncid,'heat_tr',id)
   iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr(js_pe:je_pe))
   iret=nf_inq_varid(ncid,'salt_tr',id)
   iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr(js_pe:je_pe))

   if (number_masks >0) then
    iret=nf_inq_varid(ncid,'vsf_mask1',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_vsf_iso_1(js_pe:je_pe,:))
    iret=nf_inq_varid(ncid,'vsf_mask2',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_vsf_iso_2(js_pe:je_pe,:))  
    iret=nf_inq_varid(ncid,'heat_tr_1',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr_1(js_pe:je_pe))
    iret=nf_inq_varid(ncid,'salt_tr_1',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr_1(js_pe:je_pe))
    iret=nf_inq_varid(ncid,'heat_tr_2',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr_2(js_pe:je_pe))
    iret=nf_inq_varid(ncid,'salt_tr_2',id)
    iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr_2(js_pe:je_pe))
   endif
   
   if ((enable_neutral_diffusion .and. enable_skew_diffusion).or.enable_rossmix2) then
     iret=nf_inq_varid(ncid,'bolus_iso',id)
     iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_bolus_iso(js_pe:je_pe,:))
     iret=nf_inq_varid(ncid,'bolus_depth',id)
     iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_bolus_depth(js_pe:je_pe,:))
     
     iret=nf_inq_varid(ncid,'heat_tr_bolus',id)
     iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr_bolus(js_pe:je_pe))
     iret=nf_inq_varid(ncid,'salt_tr_bolus',id)
     iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr_bolus(js_pe:je_pe))

     if (number_masks >0) then
      iret=nf_inq_varid(ncid,'bolus_mask1',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_bolus_iso_1(js_pe:je_pe,:))
      iret=nf_inq_varid(ncid,'bolus_mask2',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,1,ilen/), (/je_pe-js_pe+1,nz,1/),mean_bolus_iso_2(js_pe:je_pe,:))      
      iret=nf_inq_varid(ncid,'heat_tr_bolus_1',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr_bolus_1(js_pe:je_pe))
      iret=nf_inq_varid(ncid,'salt_tr_bolus_1',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr_bolus_1(js_pe:je_pe))
      iret=nf_inq_varid(ncid,'heat_tr_bolus_2',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_heat_tr_bolus_2(js_pe:je_pe))
      iret=nf_inq_varid(ncid,'salt_tr_bolus_2',id)
      iret= nf_put_vara_double(ncid,id,(/js_pe,ilen/), (/je_pe-js_pe+1,1/),mean_salt_tr_bolus_2(js_pe:je_pe))      
     endif  
   endif
   iret=nf_close(ncid)
  endif
  call fortran_barrier
 enddo

 nitts = 0
 mean_trans = 0d0
 mean_vsf_iso = 0d0
 mean_vsf_depth =0d0
 mean_heat_tr =0d0
 mean_salt_tr =0d0
 if (number_masks >0) then
  mean_vsf_iso_1 = 0d0
  mean_vsf_iso_2 = 0d0
  mean_heat_tr_1 =0d0
  mean_salt_tr_1 =0d0
  mean_heat_tr_2 =0d0
  mean_salt_tr_2 =0d0
 endif
 if ((enable_neutral_diffusion .and. enable_skew_diffusion).or.enable_rossmix2) then
   mean_bolus_iso = 0d0
   mean_bolus_depth = 0d0
   mean_heat_tr_bolus =0d0
   mean_salt_tr_bolus =0d0
   if (number_masks >0) then
    mean_bolus_iso_1 = 0d0
    mean_bolus_iso_2 = 0d0
    mean_heat_tr_bolus_1 =0d0
    mean_salt_tr_bolus_1 =0d0
    mean_heat_tr_bolus_2 =0d0
    mean_salt_tr_bolus_2 =0d0
   endif
 endif
end subroutine write_overturning



subroutine diag_over_read_restart
!=======================================================================
! read unfinished averages from file
!=======================================================================
 use main_module
 use isoneutral_module
 use module_diag_overturning
 use rossmix2_module
 implicit none
 character (len=80) :: filename
 logical :: file_exists
 integer :: io,ierr,ny_,nz_,nl_,js_,je_

 if (my_blk_i>1) return ! no need to read anything
 
 write(filename,'(a,i5,a)')  'unfinished_over_PE_',my_pe,'.dta'
 call replace_space_zero(filename)
 inquire ( FILE=filename, EXIST=file_exists )
 if (.not. file_exists) then
      if (my_pe==0) then
         print'(a,a,a)',' file ',filename(1:len_trim(filename)),' not present'
         print'(a)',' reading no unfinished overturning diagnostics'
      endif
      return
 endif

 if (my_pe==0) print'(2a)',' reading unfinished averages from ',filename(1:len_trim(filename))
 call get_free_iounit(io,ierr)
 if (ierr/=0) goto 10
 open(io,file=filename,form='unformatted',status='old',err=10)
 read(io,err=10) nitts,ny_,nz_,nl_
 if (ny/=ny_ .or. nz/= nz_ .or. nl_ /=nlevel) then 
       if (my_pe==0) then
        print*,' read dimensions: ',ny_,nz_,nl_
        print*,' does not match dimensions   : ',ny,nz,nlevel
       endif
       goto 10
 endif
 read(io,err=10) js_,je_
 if (js_/=js_pe.or.je_/=je_pe) then
       if (my_pe==0) then
        print*,' read PE boundaries   ',js_,je_
        print*,' which does not match ',js_pe,je_pe
       endif
       goto 10
 endif
 read(io,err=10) mean_trans,mean_vsf_iso,mean_bolus_iso,mean_vsf_depth,mean_bolus_depth
 read(io,err=10) mean_heat_tr,mean_salt_tr,mean_heat_tr_bolus,mean_salt_tr_bolus
 if (number_masks >0) then
   read(io,err=10) mean_vsf_iso_1,mean_bolus_iso_1,mean_vsf_iso_2,mean_bolus_iso_2
   read(io,err=10) mean_heat_tr_1,mean_salt_tr_1,mean_heat_tr_2,mean_salt_tr_2
   read(io,err=10) mean_heat_tr_bolus_1,mean_salt_tr_bolus_1,mean_heat_tr_bolus_2,mean_salt_tr_bolus_2
 endif
 close(io)
 return
 10 continue
 print'(a)',' Warning: error reading file'
end subroutine diag_over_read_restart




subroutine diag_over_write_restart
!=======================================================================
! write unfinished averages to restart file
!=======================================================================
 use main_module
 use isoneutral_module
 use module_diag_overturning
 use rossmix2_module
 implicit none
 character (len=80) :: filename
 integer :: io,ierr

 if (my_blk_i>1) return ! no need to write anything

 write(filename,'(a,i5,a)')  'unfinished_over_PE_',my_pe,'.dta'
 call replace_space_zero(filename)
 if (my_pe==0) print'(a,a)',' writing unfinished averages to ',filename(1:len_trim(filename))
 call get_free_iounit(io,ierr)
 if (ierr/=0) goto 10
 open(io,file=filename,form='unformatted',status='unknown')
 write(io,err=10) nitts,ny,nz,nlevel
 write(io,err=10) js_pe,je_pe
 write(io,err=10) mean_trans,mean_vsf_iso,mean_bolus_iso,mean_vsf_depth,mean_bolus_depth
 write(io,err=10) mean_heat_tr,mean_salt_tr,mean_heat_tr_bolus,mean_salt_tr_bolus
 if (number_masks >0) then
  write(io,err=10) mean_vsf_iso_1,mean_bolus_iso_1,mean_vsf_iso_2,mean_bolus_iso_2
  write(io,err=10) mean_heat_tr_1,mean_salt_tr_1,mean_heat_tr_2,mean_salt_tr_2
  write(io,err=10) mean_heat_tr_bolus_1,mean_salt_tr_bolus_1,mean_heat_tr_bolus_2,mean_salt_tr_bolus_2
 endif 
 close(io)
 return
 10 continue
 print'(a)',' Warning: error writing file'
end subroutine diag_over_write_restart



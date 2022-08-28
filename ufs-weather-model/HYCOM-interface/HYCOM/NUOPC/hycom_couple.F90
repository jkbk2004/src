!===============================================================================
! MODULE: HYCOM_COUPLE
!
! DESCRIPTION:
!   This module copies internal HYCOM variables to/from NUOPC fields.
!
! SUBROUTINES:
!   hycom_couple_init
!     Initialize decomposition blocks, coordinates, mask, and area.
!
!   set_hycom_import_flag
!     Determine field connections.
!
!   export_from_hycom_deb
!     Copy internal HYCOM variable to export fortran array.
!
!   import_to_hycom_deb
!     Copy import fortran array to internal HYCOM variable.
!
!   ocn_import_forcing
!     Process uncoupled variables.
!
!   hycom_couple_final
!     Deallocate memory.
!
!===============================================================================
module hycom_couple

!===============================================================================
! use modules
!===============================================================================
  use mod_xc ! HYCOM communication interface
  use mod_cb_arrays
  use hycom_read_latlon, only: get_coord

!===============================================================================
! settings
!===============================================================================
  implicit none

  private

!===============================================================================
! public
!===============================================================================
  public :: hycom_couple_init
  public :: set_hycom_import_flag
  public :: hycom_couple_check_deb
  public :: export_from_hycom_deb
  public :: import_to_hycom_deb
  public :: ocn_import_forcing
  public :: hycom_couple_final
  public :: cpldom_type
  public :: cpldom

!===============================================================================
! module variables
!===============================================================================
  type cpldom_type
    integer              :: idim_size
    integer              :: jdim_size
    integer, allocatable :: deBList(:,:,:)
    real, allocatable    :: lon_p(:,:)
    real, allocatable    :: lat_p(:,:)
    real, allocatable    :: mask_p(:,:)
    real, allocatable    :: area_p(:,:)
    real, allocatable    :: lon_q(:,:)
    real, allocatable    :: lat_q(:,:)
    real, allocatable    :: mask_q(:,:)
    real, allocatable    :: area_q(:,:)
  end type cpldom_type

  type(cpldom_type) :: cpldom

!===============================================================================
  contains
!===============================================================================
  subroutine hycom_couple_init(nPets, diag, rc)
!   arguments
    integer, intent(in)  :: nPets
    logical, intent(in)  :: diag
    integer, intent(out) :: rc
!   local variables
    character(*), parameter :: rname="hycom_couple_init"
    real, allocatable       :: tmx(:,:)
    real, allocatable       :: tmp_e(:,:)
    integer                 :: i, j

    rc = 0 ! success

!   grid size
    cpldom%idim_size=itdm
    cpldom%jdim_size=jtdm

#ifdef ESPC_COUPLE
!   deBlockList
!   directly from HYCOM
    if (.not.allocated(cpldom%deBList)) then
      allocate(cpldom%deBList(2,2,nPets))
    endif
    do i=1, nPets
      cpldom%deBList(1,1,i)=deBlockList(1,1,i)
      cpldom%deBList(2,1,i)=deBlockList(2,1,i)
      cpldom%deBList(1,2,i)=deBlockList(1,2,i)
      cpldom%deBList(2,2,i)=deBlockList(2,2,i)
    enddo
    if (mnproc.eq.1) then
      print *,'itdm,jtdm=',itdm,jtdm
      print *,'hycom,deBList BL11 BL21 BL12 BL22'
      do i=1, nPets
        write(*,"(I4,4I8,3x,2I8)") i,                    &
          cpldom%deBList(1,1,i),cpldom%deBList(2,1,i),   &
          cpldom%deBList(1,2,i),cpldom%deBList(2,2,i),   &
          cpldom%deBList(1,2,i)-cpldom%deBList(1,1,i)+1, &
          cpldom%deBList(2,2,i)-cpldom%deBList(2,1,i)+1
      enddo
    endif
#endif

!   allocate arrays
    if (mnproc.eq.1) then
      if (.not.allocated(cpldom%lon_p)) allocate(cpldom%lon_p(itdm,jtdm))
      if (.not.allocated(cpldom%lat_p)) allocate(cpldom%lat_p(itdm,jtdm))
      if (.not.allocated(cpldom%area_p)) allocate(cpldom%area_p(itdm,jtdm))
      if (.not.allocated(cpldom%mask_p)) allocate(cpldom%mask_p(itdm,jtdm))
      if (.not.allocated(cpldom%lon_q)) allocate(cpldom%lon_q(itdm,jtdm))
      if (.not.allocated(cpldom%lat_q)) allocate(cpldom%lat_q(itdm,jtdm))
      if (.not.allocated(cpldom%area_q)) allocate(cpldom%area_q(itdm,jtdm))
      if (.not.allocated(cpldom%mask_q)) allocate(cpldom%mask_q(itdm,jtdm))
!     read hycom regional.grid.a
      call get_coord(cpldom%lat_p, cpldom%lon_p, cpldom%lat_q, &
        cpldom%lon_q, itdm, jtdm, rc)
    else
      if (.not.allocated(cpldom%lon_p)) allocate(cpldom%lon_p(1,1))
      if (.not.allocated(cpldom%lat_p)) allocate(cpldom%lat_p(1,1))
      if (.not.allocated(cpldom%area_p)) allocate(cpldom%area_p(1,1))
      if (.not.allocated(cpldom%mask_p)) allocate(cpldom%mask_p(1,1))
      if (.not.allocated(cpldom%lon_q)) allocate(cpldom%lon_q(1,1))
      if (.not.allocated(cpldom%lat_q)) allocate(cpldom%lat_q(1,1))
      if (.not.allocated(cpldom%area_q)) allocate(cpldom%area_q(1,1))
      if (.not.allocated(cpldom%mask_q)) allocate(cpldom%mask_q(1,1))
      cpldom%lat_p(:,:)=0.0
      cpldom%lon_p(:,:)=0.0
      cpldom%area_p(:,:)=0.0
      cpldom%lat_q(:,:)=0.0
      cpldom%lon_q(:,:)=0.0
      cpldom%area_q(:,:)=0.0
    endif

!   mask information
    if (.not.allocated(tmx)) allocate(tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

!   mask_p information
    tmx(:,:)=0.0
    do j=1, jj
    do i=1, ii
      tmx(i,j)=ishlf(i,j)
    enddo
    enddo
    call xcaget(cpldom%mask_p, tmx, 1)

!   mask_q information
    tmx(:,:)=0.0
    do j=1, jj
    do i=1, ii
      tmx(i,j)=iq(i,j)
    enddo
    enddo
    call xcaget(cpldom%mask_q, tmx, 1)

!   initialize couple flags
    cpl_merge = .false.
    cpl_diag = diag
    calc_wndspd=.false.
    calc_radflx=.false.
    cpl_taux=.false.
    cpl_tauy=.false.
    cpl_u10=.false.
    cpl_v10=.false.
    cpl_wndspd=.false.
    cpl_ustara=.false.
    cpl_airtmp=.false.
    cpl_vapmix=.false.
    cpl_swflx_net=.false.
    cpl_lwflx_net=.false.
    cpl_swflx_net2down=.false.
    cpl_lwflx_net2down=.false.
    cpl_swflxd=.false.
    cpl_lwflxd=.false.
    cpl_mslprs=.false.
    cpl_precip=.false.
    cpl_surtmp=.false.
    cpl_seatmp=.false.
    cpl_sbhflx=.false.
    cpl_sensflx=.false.
    cpl_lthflx=.false.
    cpl_latflx=.false.
    cpl_sic=.false.
    cpl_sitx=.false.
    cpl_sity=.false.
    cpl_siqs=.false.
    cpl_sifh=.false.
    cpl_sifs=.false.
    cpl_sifw=.false.
    cpl_sit=.false.
    cpl_sih=.false.
    cpl_siu=.false.
    cpl_siv=.false.

  end subroutine hycom_couple_init

  !-----------------------------------------------------------------------------

  subroutine set_hycom_import_flag(fieldName, rc)
!   arguments
    character(len=30), intent(in) :: fieldName
    integer, intent(out)          :: rc
!   local variables
    character(*), parameter :: rname="set_hycom_import_flag"

    rc = 0 ! success

!   set couple flags based on fieldnames
    if (fieldName.eq.'taux10') then
      cpl_taux=.true.
    elseif (fieldName.eq.'tauy10') then
      cpl_tauy=.true.
      if (.not.cpl_taux) then
        if (mnproc.eq.1) print *,"error - tauy before taux"
        call xcstop('('//rname//')')
               stop '('//rname//')'
      endif !error
    elseif (fieldName.eq.'u10') then
      cpl_u10=.true.
    elseif (fieldName.eq.'v10') then
      cpl_v10=.true.
      if (.not.cpl_u10) then
        if (mnproc.eq.1) print *,"error - v10 before u10"
        call xcstop('('//rname//')')
               stop '('//rname//')'
      endif !error
    elseif (fieldName.eq.'wndspd10') then
      cpl_wndspd=.true.
    elseif (fieldName.eq.'ustara10') then
      cpl_ustara=.true.
    elseif (fieldName.eq.'airtmp') then
      cpl_airtmp=.true.
    elseif (fieldName.eq.'airhum') then
      cpl_vapmix=.true.
    elseif (fieldName.eq.'swflx_net') then
      cpl_swflx_net=.true.
    elseif (fieldName.eq.'lwflx_net') then
      cpl_lwflx_net=.true.
    elseif (fieldName.eq.'swflx_net2down') then
      cpl_swflx_net2down=.true.
    elseif (fieldName.eq.'lwflx_net2down') then
      cpl_lwflx_net2down=.true.
    elseif (fieldName.eq.'swflxd') then
      cpl_swflxd=.true.
    elseif (fieldName.eq.'lwflxd') then
      cpl_lwflxd=.true.
    elseif (fieldName.eq.'mslprs') then
      cpl_mslprs=.true.
    elseif (fieldName.eq.'prcp') then
      cpl_precip=.true.
    elseif (fieldName.eq.'gt') then
      cpl_surtmp=.true.
      if (sstflg.ne.3) then
        cpl_seatmp=.true.
      endif
    elseif (fieldName.eq.'sbhflx') then
      cpl_sbhflx=.true.
    elseif (fieldName.eq.'sensflx') then
      cpl_sensflx=.true.
    elseif (fieldName.eq.'lthflx') then
      cpl_lthflx=.true.
    elseif (fieldName.eq.'latflx') then
      cpl_latflx=.true.
!   import ice concentration
    elseif (fieldName.eq.'sic') then
      cpl_sic=.true.
!   import ice x-stress
    elseif (fieldName.eq.'sitx') then
      cpl_sitx=.true.
!   import ice y-stress
    elseif (fieldName.eq.'sity') then
      cpl_sity=.true.
!   import solar thru grid cell ave.
    elseif (fieldName.eq.'siqs') then
      cpl_siqs=.true.
!   import freeze, melt, h. flux
    elseif (fieldName.eq.'sifh') then
      cpl_sifh=.true.
!   import salt flux
    elseif (fieldName.eq.'sifs') then
      cpl_sifs=.true.
!   import water flux
    elseif (fieldName.eq.'sifw') then
      cpl_sifw=.true.
!   import sea ice temperature
    elseif (fieldName.eq.'sit_sfc') then
      cpl_sit=.true.
!   import sea ice thickness
    elseif (fieldName.eq.'sih') then
      cpl_sih=.true.
!   import sea ice x-velocity
    elseif (fieldName.eq.'siu') then
      cpl_siu=.true.
!   import sea ice y-velocity
    elseif (fieldName.eq.'siv') then
      cpl_siv=.true.
    else ! error
      if (mnproc.eq.1) print *, "error - fieldName unknown: "//trim(fieldName)
      call xcstop('('//rname//')')
      stop '('//rname//')'
    endif !fieldName

  end subroutine set_hycom_import_flag

  !-----------------------------------------------------------------------------

  subroutine hycom_couple_check_deb(show_minmax, rc)
!   arguments
    logical, intent(in)  :: show_minmax
    integer, intent(out) :: rc
!   local variables
    character(*), parameter :: rname="hycom_couple_check_deb"
    real, allocatable       :: data_tmp(:,:)

    rc = 0 ! success

!   diagnostic output
    if (show_minmax) then
!     allocate memory
      if (mnproc.eq.1) then
        allocate(data_tmp(itdm,jtdm))
      else
        allocate(data_tmp(1,1))
      endif

!     check pang - angle between xwards and ewards
      call xcaget(data_tmp, pang, 1)
      if (mnproc.eq.1) then
        print *,rname//' pang, min,max=', &
          minval(data_tmp), maxval(data_tmp)
      endif

!     deallocate memory
      if (allocated(data_tmp)) deallocate(data_tmp)
    endif !show_minmax

  end subroutine hycom_couple_check_deb

  !-----------------------------------------------------------------------------

  subroutine export_from_hycom_deb(tlb, tub, expData, fieldName, &
  show_minmax, rc)
!   arguments
    integer, intent(in)           :: tlb(2)
    integer, intent(in)           :: tub(2)
    real, intent(inout)           :: expData(tlb(1):tub(1),tlb(2):tub(2))
    character(len=30), intent(in) :: fieldName
    logical, intent(in)           :: show_minmax
    integer, intent(out)          :: rc
!   local variables
    character(*), parameter :: rname="export_from_hycom_deb"
    real, allocatable       :: ocn_msk(:,:)
    real, allocatable       :: field_tmp(:,:)
    real, allocatable       :: tmx(:,:)
    integer                 :: i, j, jja
!   integer                 :: k
!   real                    :: mgrid(ii,jj)

    rc = 0 ! success

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

!   print *,"idm,jdm,nbdy,ii,jj=",mnproc,idm,jdm,nbdy,ii,jj

    call export_from_hycom_tiled(util2, fieldName) !can't use util1

#if defined(ARCTIC)
!   arctic (tripole) domain, top row is replicated (ignore it)
    jja=min(jj,(jtdm-1-j0))
#else
    jja=jj
#endif

#ifndef ESPC_NOCANONICAL_CONVERT
!   export unit conversions
    if (fieldName.eq.'sst') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
!         canonical unit conversion: sst (C) -> (K)
          util2(i,j)=util2(i,j)+273.15d0
        end if
      enddo
      enddo
    endif
#endif

!   copy internal data to export array
    expData(:,:)=0.0
    do j=1, jja
    do i=1, ii
!     mgrid(i,j)=util2(i,j)
      expData(i+i0,j+j0)=util2(i,j)
    enddo
    enddo

!   diagnostic output
    if (show_minmax) then
!     allocate memory
      if (mnproc.eq.1) then
        allocate(ocn_msk(itdm,jtdm))
        allocate(field_tmp(itdm,jtdm))
      else
        allocate(ocn_msk(1,1))
        allocate(field_tmp(1,1))
      endif
      allocate(tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

!     ocn_msk: sea/land mask
      tmx(:,:)=0.0
      do j=1, jja
      do i=1, ii
        tmx(i,j)=ishlf(i,j)
!x      tmx(i,j)=ip(i,j)
      enddo
      enddo
      call xcaget(ocn_msk,tmx,1)
!     call xcsync(no_flush)

!     field_tmp: export field data
      tmx(:,:)=0.0
      do j=1, jja
      do i=1, ii
        tmx(i,j)=expData(i+i0,j+j0)
      enddo
      enddo
      call xcaget(field_tmp,tmx,1)
!     call xcsync(no_flush)

!     write minmax to stdout
      if (mnproc.eq.1) then
        write(*,992) trim(fieldName),          &
          maxval(field_tmp,mask=ocn_msk.eq.1), &
          minval(field_tmp,mask=ocn_msk.eq.1), &
          (sum(field_tmp,mask=ocn_msk.eq.1)/count(ocn_msk.eq.1))
 992    format('export_from_hycom_deb,max,min,mean=',A10,3E23.15)
      endif

!     deallocate memory
      if (allocated(ocn_msk)) deallocate(ocn_msk)
      if (allocated(field_tmp)) deallocate(field_tmp)
      if (allocated(tmx)) deallocate(tmx)
    endif !show_minmax

  end subroutine export_from_hycom_deb

  !-----------------------------------------------------------------------------

  subroutine import_to_hycom_deb(tlb, tub, impData, fill_value, fieldName, &
  show_minmax, rc)
!   arguments
    integer, intent(in)           :: tlb(2), tub(2)
    real*8, intent(in)            :: impData(tlb(1):tub(1),tlb(2):tub(2))
    real*8, intent(in)            :: fill_value
    character(len=30), intent(in) :: fieldName
    logical, intent(in)           :: show_minmax
    integer, intent(out)          :: rc
!   local variables
    character(*), parameter :: rname="import_to_hycom_deb"
    integer                 :: i, j, mcnt
    real                    :: uij, vij
    logical, allocatable    :: fld_msk(:,:)
    real                    :: fld_max
    real                    :: fld_min
    real                    :: fld_sum
    real                    :: fld_cnt
    real, parameter         :: sstmin=-1.8d0
    real, parameter         :: sstmax=35.0d0
    integer                 :: jja
    real                    :: albw,degtorad
    integer                 :: ierr
!   integer                 :: k

    rc = 0 ! success

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

!   if ((k.eq.1).and.(mnproc.eq.1)) print *, "w0,w1..=", w0, w1, w2, w3

#if defined(ARCTIC)
!   arctic (tripole) domain, top row is replicated (ignore it)
    jja=min(jj,(jtdm-1-j0))
#else
    jja=jj
#endif

!   -----------------
!   check fill_value
    if ((fill_value.gt.-1.0e10).and.(fill_value.lt.1.0e10)) then
      if (mnproc.eq.1) print *,"error - fill_value close to zero"
      call xcstop('('//rname//')')
             stop '('//rname//')'
    endif

!   -----------------
!    set import flag
!   -----------------
    if (.not.cpl_merge) then
      do j=1, jja
      do i=1, ii
        if (impData(i+i0,j+j0).eq.fill_value) then
          imp_merge(i,j)=0.d0
        else
          imp_merge(i,j)=1.d0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_merge(:,:),1,1,halo_ps)
#endif
      call xctilr(imp_merge(:,:),1,1,nbdy,nbdy,halo_ps)
      cpl_merge = .true.
    endif

!   -----------------
!    import from atm
!   -----------------
!   import xstress: Pa
    if (fieldName.eq.'taux10') then
      do j=1, jja
      do i=1, ii
!       imp_taux(i,j,1)=mgrid(i,j)
        if (ishlf(i,j).eq.1) then
          imp_taux(i,j,1)=impData(i+i0,j+j0)
        else
          imp_taux(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_taux(:,:,1),1,1,halo_pv)
#endif
      call xctilr(imp_taux(:,:,1),1,1,nbdy,nbdy,halo_pv)
!   -----------------
!   import ystress: Pa
    elseif (fieldName.eq.'tauy10') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_tauy(i,j,1)=impData(i+i0,j+j0)
          if (impData(i+i0,j+j0).ne.fill_value) then
!           rotate taux and tauy to (x,y)ward
!           assumes rotation only needed for impData
            uij=imp_taux(i,j,1)
            vij=imp_tauy(i,j,1)
            imp_taux(i,j,1)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
            imp_tauy(i,j,1)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
          endif
        else
          imp_tauy(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_taux(:,:,1),1,1,halo_pv)
      call xctila(imp_tauy(:,:,1),1,1,halo_pv)
#endif
      call xctilr(imp_taux(:,:,1),1,1,nbdy,nbdy,halo_pv)
      call xctilr(imp_tauy(:,:,1),1,1,nbdy,nbdy,halo_pv)
!   -----------------
!   import u wind at 10m height: ms-1
    elseif (fieldName.eq.'u10') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_wndspx(i,j,1)=impData(i+i0,j+j0)
        else
          imp_wndspx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_wndspx(:,:,1),1,1,halo_pv)
#endif
      call xctilr(imp_wndspx(:,:,1),1,1,nbdy,nbdy,halo_pv)
!   -----------------
!   import v wind at 10m height: ms-1
    elseif (fieldName.eq.'v10') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_wndspy(i,j,1)=impData(i+i0,j+j0)
          if (impData(i+i0,j+j0).ne.fill_value) then
!           rotate u and v to (x,y)ward
!           assumes rotation only needed for impData
            uij=imp_wndspx(i,j,1)
            vij=imp_wndspy(i,j,1)
            imp_wndspx(i,j,1)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
            imp_wndspy(i,j,1)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
          endif
        else
          imp_wndspy(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_wndspx(:,:,1),1,1,halo_pv)
      call xctila(imp_wndspy(:,:,1),1,1,halo_pv)
#endif
      call xctilr(imp_wndspx(:,:,1),1,1,nbdy,nbdy,halo_pv)
      call xctilr(imp_wndspy(:,:,1),1,1,nbdy,nbdy,halo_pv)
!   -----------------
!   import wind speed: m s-1
    elseif (fieldName.eq.'wndspd10') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_wndspd(i,j,1)=impData(i+i0,j+j0)
        else
          imp_wndspd(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_wndspd(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_wndspd(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import friction speed: m s-1
    elseif (fieldName.eq.'ustara10') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_ustara(i,j,1)=impData(i+i0,j+j0)
        else
          imp_ustara(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_ustara(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_ustara(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import air temperature
    elseif (fieldName.eq.'airtmp') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if(impData(i+i0,j+j0).ne.fill_value) then
!           canonical unit conversion: airtmp (K) -> (C)
            imp_airtmp(i,j,1)=impData(i+i0,j+j0)-273.15
          else
            imp_airtmp(i,j,1)=impData(i+i0,j+j0)
          endif
        else
          imp_airtmp(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_airtmp(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_airtmp(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import specific humidity: kg kg-1
    elseif (fieldName.eq.'airhum') then
      if (flxflg.ne.5) then !Alex flxflg.eq.4 => mixing ratio
        do j=1, jja
        do i=1, ii
          if (ishlf(i,j).eq.1) then
            if (impData(i+i0,j+j0).ne.fill_value) then
!             convert from specific humidity to mixing ratio
              imp_vapmix(i,j,1)=impData(i+i0,j+j0)/(1.-impData(i+i0,j+j0))
            else
              imp_vapmix(i,j,1)=impData(i+i0,j+j0)
            endif
          else
            imp_vapmix(i,j,1)=0.01
          endif
        enddo
        enddo
      else !Alex flxflg.eq.5 => specific humidity
        do j=1, jja
        do i=1, ii
          if (ishlf(i,j).eq.1) then
            imp_vapmix(i,j,1)=impData(i+i0,j+j0)
          else
            imp_vapmix(i,j,1)=0.01
          endif
        enddo
        enddo
      endif
#if defined(ARCTIC)
      call xctila(imp_vapmix(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_vapmix(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import sw flux: w m-2
    elseif (fieldName.eq.'swflx_net') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_swflx(i,j,1)=impData(i+i0,j+j0)
        else
          imp_swflx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_swflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_swflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import downward sw flux: w m-2
    elseif ((fieldName.eq.'swflx_net2down').or. &
            (fieldName.eq.'swflxd')) then
      if (albflg.ne.0) then !swflx is Qswdn
!       use the same method as on forfun.F
!       convert to net shortwave into the ocean to swflx
!       shortwave through sea ice is handled separately
        if (albflg.eq.1) then
          do j=1, jja
          do i=1, ii
            if (ishlf(i,j).eq.1) then
              if (impData(i+i0,j+j0).ne.fill_value) then
                imp_swflx(i,j,1)=impData(i+i0,j+j0)*(1.0-0.09) !NAVGEM albedo
              else
                imp_swflx(i,j,1)=impData(i+i0,j+j0)
              endif
            else
              imp_swflx(i,j,1)=0.0
            endif
          enddo
          enddo
        else !albflg.eq.2
          degtorad=4.d0*atan(1.d0)/180.d0
          do j=1, jja
          do i=1, ii
            if (ishlf(i,j).eq.1) then
              if (impData(i+i0,j+j0).ne.fill_value) then
!               latitudinally-varying ocean albedo (Large and Yeager, 2009)
!               5.8% at the equator and 8% at the poles
                albw=(0.069-0.011*cos(2.0*degtorad*plat(i,j)))
                imp_swflx(i,j,1)=impData(i+i0,j+j0)*(1.0-albw)
              else
                imp_swflx(i,j,1)=impData(i+i0,j+j0)
              endif
            else
              imp_swflx(i,j,1)=0.0
            endif
          enddo
          enddo
        endif
      else
        do j=1, jja
        do i=1, ii
          if (ishlf(i,j).eq.1) then
            imp_swflx(i,j,1)=impData(i+i0,j+j0)
          else
            imp_swflx(i,j,1)=0.0
          endif
        enddo
        enddo
      endif !albflg
#if defined(ARCTIC)
      call xctila(imp_swflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_swflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import lw flux: w m-2
    elseif (fieldName.eq.'lwflx_net') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if (impData(i+i0,j+j0).ne.fill_value) then
!           canonical unit conversion: lwflx_net (upward) -> (downward)
            imp_lwdflx(i,j,1)=impData(i+i0,j+j0)*(-1.)
          else
            imp_lwdflx(i,j,1)=impData(i+i0,j+j0)
          endif
        else
          imp_lwdflx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_lwdflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_lwdflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import downward lw flux: w m-2
!   +ve into ocean
    elseif ((fieldName.eq.'lwflx_net2down').or. &
            (fieldName.eq.'lwflxd')) then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_lwdflx(i,j,1)=impData(i+i0,j+j0)
        else
          imp_lwdflx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_lwdflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_lwdflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import precip: m s-1
    elseif (fieldName.eq.'prcp') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if (impData(i+i0,j+j0).ne.fill_value) then
!           canonical unit conversion: prcp (kg_m-2_s-1) -> (m_s-1)
            imp_precip(i,j,1)=impData(i+i0,j+j0)*(0.001)
          else
            imp_precip(i,j,1)=impData(i+i0,j+j0)
          endif
        else
          imp_precip(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_precip(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_precip(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import surface temperature
    elseif (fieldName.eq.'gt') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if (impData(i+i0,j+j0).ne.fill_value) then
!           canonical unit conversion: gt (K) -> (C)
            imp_surtmp(i,j,1)=impData(i+i0,j+j0)-273.15
          else
            imp_surtmp(i,j,1)=impData(i+i0,j+j0)
          endif
        else
          imp_surtmp(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_surtmp(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_surtmp(:,:,1),1,1,nbdy,nbdy,halo_ps)
      if (sstflg.ne.3) then !use atmos sst as "truth"
        do j=1, jja
        do i=1, ii
          if (imp_surtmp(i,j,1).ne.fill_value) then
            imp_seatmp(i,j,1)=max(sstmin,min(imp_surtmp(i,j,1),sstmax))
          else
            imp_seatmp(i,j,1)=imp_surtmp(i,j,1)
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(imp_seatmp(:,:,1),1,1,halo_ps)
#endif
        call xctilr(imp_seatmp(:,:,1),1,1,nbdy,nbdy,halo_ps)
      endif
!   -----------------
!   import latent heat flux: w m-2
    elseif (fieldName.eq.'latflx') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_latflx(i,j,1)=impData(i+i0,j+j0)
        else
          imp_latflx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_latflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_latflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import sensible heat flux: w m-2
    elseif (fieldName.eq.'sensflx') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          imp_sensflx(i,j,1)=impData(i+i0,j+j0)
        else
          imp_sensflx(i,j,1)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_sensflx(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_sensflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   -----------------
!   import sea level pressure anomaly: Pa
    elseif (fieldName.eq.'mslprs') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if (impData(i+i0,j+j0).ne.fill_value) then
            imp_mslprs(i,j,1)=impData(i+i0,j+j0)-prsbas
          else
            imp_mslprs(i,j,1)=impData(i+i0,j+j0)
          endif
        else
          imp_mslprs(i,j,1)=101000.0-prsbas
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_mslprs(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_mslprs(:,:,1),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!    import from sea ice
!   ---------------------
!   import ice concentration
    elseif (fieldName.eq.'sic') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sic_import(i,j)=impData(i+i0,j+j0)
        else
          sic_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sic_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sic_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import ice x-stress
    elseif (fieldName.eq.'sitx') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sitx_import(i,j)=impData(i+i0,j+j0)
        else
          sitx_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sitx_import(:,:),1,1,halo_pv)
#endif
      call xctilr(sitx_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import ice y-stress
    elseif (fieldName.eq.'sity') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sity_import(i,j)=impData(i+i0,j+j0)
        else
          sity_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sity_import(:,:),1,1,halo_pv)
#endif
      call xctilr(sity_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import solar thru grid cell ave.
    elseif (fieldName.eq.'siqs') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          siqs_import(i,j)=impData(i+i0,j+j0)
        else
          siqs_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(siqs_import(:,:),1,1,halo_ps)
#endif
      call xctilr(siqs_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import freeze, melt, H. Flux
    elseif (fieldName.eq.'sifh') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sifh_import(i,j)=impData(i+i0,j+j0)
        else
          sifh_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sifh_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sifh_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import salt flux
    elseif (fieldName.eq.'sifs') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sifs_import(i,j)=impData(i+i0,j+j0)
        else
          sifs_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sifs_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sifs_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import water flux
    elseif (fieldName.eq.'sifw') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sifw_import(i,j)=impData(i+i0,j+j0)
        else
          sifw_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sifw_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sifw_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import sea ice temperature
    elseif (fieldName.eq.'sit_sfc') then
#ifndef ESPC_NOCANONICAL_CONVERT
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          if (impData(i+i0,j+j0).ne.fill_value) then
!           canonical unit conversion: sit_sfc (K) -> (C)
            sit_import(i,j)=impData(i+i0,j+j0)-273.15
          else
            sit_import(i,j)=impData(i+i0,j+j0)
          endif
        else
          sit_import(i,j)=0.0
        endif
      enddo
      enddo
#else
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sit_import(i,j)=impData(i+i0,j+j0)
        else
          sit_import(i,j)=0.0
        endif
      enddo
      enddo
#endif
#if defined(ARCTIC)
      call xctila(sit_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sit_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import sea ice thickness
    elseif (fieldName.eq.'sih') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          sih_import(i,j)=impData(i+i0,j+j0)
        else
          sih_import(i,j)=0.0
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(sih_import(:,:),1,1,halo_ps)
#endif
      call xctilr(sih_import(:,:),1,1,nbdy,nbdy,halo_ps)
!   ---------------------
!   import sea ice x-velocity
    elseif (fieldName.eq.'siu') then
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          siu_import(i,j)=impData(i+i0,j+j0)
        else
          siu_import(i,j)=0.0
        endif
      enddo
      enddo
!   ---------------------
!   import sea ice y-velocity
    elseif (fieldName.eq.'siv') then
#ifndef ESPC_NOCANONICAL_CONVERT
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          siv_import(i,j)=impData(i+i0,j+j0)
!         rotate siu and siv to (x,y)ward
!         assumes rotation only needed for impData
          if (impData(i+i0,j+j0).ne.fill_value) then
            uij=siu_import(i,j)
            vij=siv_import(i,j)
            siu_import(i,j)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
            siv_import(i,j)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
          endif
        else
          siv_import(i,j)=0.0
        endif
      enddo
      enddo
#else
      do j=1, jja
      do i=1, ii
        if (ishlf(i,j).eq.1) then
          siv_import(i,j)=impData(i+i0,j+j0)
        else
          siv_import(i,j)=0.0
        endif
#endif
#if defined(ARCTIC)
      call xctila(siu_import(:,:),1,1,halo_pv)
      call xctila(siv_import(:,:),1,1,halo_pv)
#endif
      call xctilr(siu_import(:,:),1,1,nbdy,nbdy,halo_ps)
      call xctilr(siv_import(:,:),1,1,nbdy,nbdy,halo_ps)
    else ! field unknown error
      if (mnproc.eq.1) print *, "error - fieldName unknown: "//trim(fieldName)
      call xcstop('('//rname//')')
      stop '('//rname//')'
    endif !fieldName

!   diagnostic output
    if (show_minmax) then
!     allocate local field mask memory
      allocate(fld_msk(lbound(impData,1):ubound(impData,1), &
                       lbound(impData,2):ubound(impData,2)))
!     calculate field mask using fill value and ocean mask, ignore halos
      fld_msk(:,:)=.false.
      do j=1, jja
      do i=1, ii
        fld_msk(i+i0,j+j0)=((impData(i+i0,j+j0).ne.fill_value).and. &
                            (ishlf(i,j).eq.1))
      enddo
      enddo
!     calculate local max,min,sum,cnt
      fld_max=maxval(impData,fld_msk)
      fld_min=minval(impData,fld_msk)
      fld_sum=sum(impData,fld_msk)
      fld_cnt=real(count(fld_msk))
!     reduce max,min,sum,cnt to mnproc 1
      call xcmaxr(fld_max,1)
      call xcminr(fld_min,1)
      call xcsumr(fld_sum,1)
      call xcsumr(fld_cnt,1)
!     write max,min,mean to stdout
      if (mnproc.eq.1) then
        write(*,992) trim(fieldName),          &
          fld_max, fld_min, (fld_sum/fld_cnt)
 992    format('import_to_hycom_deb,max,min,mean=',A10,3E23.15)
      endif

!     deallocate memory
      if (allocated(fld_msk)) deallocate(fld_msk)
    endif !show_minmax

  end subroutine import_to_hycom_deb

  !-----------------------------------------------------------------------------

  subroutine ocn_import_forcing(fill_value,rc)
!   arguments
    real, intent(in)     :: fill_value
    integer, intent(out) :: rc
!   local variables
    character(*), parameter :: rname="ocn_import_forcing"
    integer                 :: i, j, m, n, jja

    rc = 0 ! success

#if defined(ARCTIC)
!   arctic (tripole) domain, top row is replicated (ignore it)
    jja=min(jj,(jtdm-1-j0))
#else
    jja=jj
#endif

!   -----------------
!   check fill_value
    if ((fill_value.gt.-1.0e10).and.(fill_value.lt.1.0e10)) then
      if (mnproc.eq.1) print *,"error - fill_value close to zero"
      call xcstop('('//rname//')')
             stop '('//rname//')'
    endif

!   -----------------
!   calculate imp_radflx
    if (lwflag.eq.0 .or. lwflag.eq.2) then
      if((cpl_lwflx_net.or.cpl_lwflx_net2down.or.cpl_lwflxd).and. &
         (cpl_swflx_net.or.cpl_swflx_net2down.or.cpl_swflxd)) then
        if (mnproc.eq.1) print *, rname//" calculating radflx..."
        do j=1, jja
        do i=1, ii
          if ((imp_lwdflx(i,j,1).ne.fill_value).and. &
              (imp_swflx(i,j,1).ne.fill_value)) then
!           imp_radflx is defined as net lwdflx+swflx, +ve into ocean
            imp_radflx(i,j,1)=imp_lwdflx(i,j,1)+imp_swflx(i,j,1)
          else
            imp_radflx(i,j,1)=fill_value
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(imp_radflx(:,:,1),1,1,halo_ps)
#endif
        call xctilr(imp_radflx(:,:,1),1,1,nbdy,nbdy,halo_ps)
        calc_radflx=.true.
      else
        calc_radflx=.false.
      endif
    else
      if (mnproc.eq.1) print *,"error - lwflag .ne. 0 or 2"
      call xcstop('('//rname//')')
             stop '('//rname//')'
    endif
!   -----------------
!   calculate imp_wndspd
    if (cpl_u10.and.cpl_v10.and.(.not.cpl_wndspd)) then
      if (mnproc.eq.1) print *, rname//" calculating wndspd..."
      calc_wndspd=.true.
      do j=1, jja
      do i=1, ii
        if ((imp_wndspx(i,j,1).ne.fill_value).and. &
            (imp_wndspy(i,j,1).ne.fill_value)) then
!         imp_wndspd based on u and v components
          imp_wndspd(i,j,1)=sqrt((imp_wndspx(i,j,1)**2)+(imp_wndspy(i,j,1)**2))
        else
          imp_wndspd(i,j,1)=fill_value
        endif
      enddo
      enddo
#if defined(ARCTIC)
      call xctila(imp_wndspd(:,:,1),1,1,halo_ps)
#endif
      call xctilr(imp_wndspd(:,:,1),1,1,nbdy,nbdy,halo_ps)
    else
      calc_wndspd=.false.
    endif

  end subroutine ocn_import_forcing

  !-----------------------------------------------------------------------------

  subroutine hycom_couple_final(rc)
!   arguments
    integer, intent(out) :: rc
!   local variables
    character(*), parameter :: rname="hycom_couple_final"

    rc = 0 ! success

!   deallocate memory
    if (allocated(cpldom%deBList)) deallocate(cpldom%deBList)
    if (allocated(cpldom%lon_p)) deallocate(cpldom%lon_p)
    if (allocated(cpldom%lat_p)) deallocate(cpldom%lat_p)
    if (allocated(cpldom%area_p)) deallocate(cpldom%area_p)
    if (allocated(cpldom%mask_p)) deallocate(cpldom%mask_p)
    if (allocated(cpldom%lon_q)) deallocate(cpldom%lon_q)
    if (allocated(cpldom%lat_q)) deallocate(cpldom%lat_q)
    if (allocated(cpldom%area_q)) deallocate(cpldom%area_q)
    if (allocated(cpldom%mask_q)) deallocate(cpldom%mask_q)

  end subroutine hycom_couple_final
!===============================================================================
end module hycom_couple
!===============================================================================


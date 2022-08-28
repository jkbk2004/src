module impexpField_cdf_mod

  !-----------------------------------------------------------------------------
  ! Field NetCDF Utilities
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  implicit none

  private

  public impexp_cdf_put_flds
  public impexp_cdf_put_latlonmsk
  public impexp_cdf_put_latlonmsk_corner

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine impexp_cdf_put_flds(cname,dtg,hours,itdmx,jtdmx, &
  numFields,fieldEnable,fieldName,standName,fieldUnit,field,status,lPet, &
  rc,label)
    character(len=*),intent(in)          :: cname          ! component name
    character(len=15),intent(in)         :: dtg            ! date and time
    real,intent(in)                      :: hours          ! hours since start
    integer,intent(in)                   :: itdmx          ! i dimension size
    integer,intent(in)                   :: jtdmx          ! j dimension size
    integer,intent(in)                   :: numFields      ! number of fields
    logical,pointer,intent(in)           :: fieldEnable(:) ! enabled fields
    character(len=30),pointer,intent(in) :: fieldName(:)   ! field names
    character(len=60),pointer,intent(in) :: standName(:)   ! standard names
    character(len=30),pointer,intent(in) :: fieldUnit(:)   ! field units
    type(ESMF_Field),intent(in)          :: field(:)       ! ESMF fields
    integer,intent(out)                  :: status         ! I/O status
    integer,intent(in)                   :: lPet           ! local PET
    integer,intent(out)                  :: rc             ! return code
    character(len=*),intent(in)          :: label          ! 'exp', 'imp'

    rc = ESMF_SUCCESS
    status = 0

  end subroutine impexp_cdf_put_flds

  !-----------------------------------------------------------------------------

  subroutine impexp_cdf_put_latlonmsk(cname,itdmx,jtdmx, &
  lat_e,lon_e,mask_e,status,rc)
    character(len=*),intent(in) :: cname       ! component name
    integer,intent(in)          :: itdmx       ! i dimension size
    integer,intent(in)          :: jtdmx       ! j dimension size
    real,intent(in)             :: lat_e(:,:)  ! latitude
    real,intent(in)             :: lon_e(:,:)  ! longitude
    real,intent(in)             :: mask_e(:,:) ! field mask
    integer,intent(out)         :: status      ! I/O status
    integer,intent(out)         :: rc          ! return code

    rc = ESMF_SUCCESS
    status = 0

  end subroutine impexp_cdf_put_latlonmsk

  !-----------------------------------------------------------------------------

  subroutine impexp_cdf_put_latlonmsk_corner(cname,itdmx,jtdmx,lat_e,lon_e,mask_e,lat_q,lon_q,status,rc)
    character(len=*),intent(in) :: cname       ! component name
    integer,intent(in)          :: itdmx       ! i dimension size
    integer,intent(in)          :: jtdmx       ! j dimension size
    real,intent(in)             :: lat_e(:,:)  ! latitude
    real,intent(in)             :: lon_e(:,:)  ! longitude
    integer,intent(in)          :: mask_e(:,:) ! field mask
    real,intent(in)             :: lat_q(:,:)  ! latitude corner
    real,intent(in)             :: lon_q(:,:)  ! longitude corner
    integer,intent(out)         :: status      ! I/O status
    integer,intent(out)         :: rc          ! return code

    rc = ESMF_SUCCESS
    status = 0

  end subroutine impexp_cdf_put_latlonmsk_corner

  !-----------------------------------------------------------------------------

end module impexpField_cdf_mod

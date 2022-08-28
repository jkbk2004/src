!=========================================================================================
! HYCOM ESMF Extensions Module
!=========================================================================================
! CONFIGURATION IDENTIFICATION $HeadURL$
! CONFIGURATION IDENTIFICATION @(#)$Id$
!=========================================================================================

! ESMF macros for logging
#define FILENAME "HYCOM_ESMF_Extensions.F90"
#define CONTEXT  line=__LINE__,file=FILENAME,method=METHOD
#define PASSTHRU msg=ESMF_LOGERR_PASSTHRU,CONTEXT

!=========================================================================================
! HYCOM ESMF Extensions Module
!=========================================================================================
module HYCOM_ESMF_Extensions

  use ESMF
  use NUOPC
  use NETCDF

  implicit none

  private

  public :: HYCOM_ESMF_GridWrite
  public :: HYCOM_ESMF_LogGrid

!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================

  contains

  !-----------------------------------------------------------------------------

#define METHOD "HYCOM_ESMF_GridWrite"
!BOP
! !IROUTINE: HYCOM_ESMF_GridWrite - Write Grid data to file
! !INTERFACE:
  subroutine HYCOM_ESMF_GridWrite(grid, fileName, overwrite, status, &
    timeslice, iofmt, relaxedflag, rc)
! ! ARGUMENTS
    type(ESMF_Grid),            intent(in)            :: grid
    character(len=*),           intent(in),  optional :: fileName
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    type(ESMF_IOFmt_Flag),      intent(in),  optional :: iofmt
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
! !DESCRIPTION:
!   Write the data in {\tt grid} to {\tt file} if supported by the
!   {\tt iofmt}.
!
!   The arguments are:
!   \begin{description}
!   \item[field]
!     The {\tt ESMF\_Field} object whose data is to be written.
!   \item[fileName]
!     The name of the file to write to. If not present then the file
!     will
!     be written to the grid's name.
!   \item[{[overwrite]}]
!      A logical flag, the default is .false., i.e., existing Field data
!      may
!      {\em not} be overwritten. If .true., the
!      data corresponding to each field's name will be
!      be overwritten. If the {\tt timeslice} option is given, only data
!      for
!      the given timeslice may be overwritten.
!      Note that it is always an error to attempt to overwrite a NetCDF
!      variable with data which has a different shape.
!   \item[{[status]}]
!      The file status. Valid options are {\tt ESMF\_FILESTATUS\_NEW},
!      {\tt ESMF\_FILESTATUS\_OLD}, {\tt ESMF\_FILESTATUS\_REPLACE}, and
!      {\tt ESMF\_FILESTATUS\_UNKNOWN} (default).
!   \item[{[timeslice]}]
!     Time slice counter. Must be positive. The behavior of this
!     option may depend on the setting of the {\tt overwrite} flag:
!     \begin{description}
!     \item[{\tt overwrite = .false.}:]\ If the timeslice value is
!     less than the maximum time already in the file, the write will
!     fail.
!     \item[{\tt overwrite = .true.}:]\ Any positive timeslice value is
!     valid.
!     \end{description}
!     By default, i.e. by omitting the {\tt timeslice} argument, no
!     provisions for time slicing are made in the output file,
!     however, if the file already contains a time axis for the
!     variable,
!     a timeslice one greater than the maximum will be written.
!   \item[{[iofmt]}]
!    The IO format.  Valid options are  {\tt ESMF\_IOFMT\_BIN} and
!    {\tt ESMF\_IOFMT\_NETCDF}. If not present, file names with a {\tt
!    .bin}
!    extension will use {\tt ESMF\_IOFMT\_BIN}, and file names with a
!    {\tt .nc}
!    extension will use {\tt ESMF\_IOFMT\_NETCDF}.  Other files default
!    to
!    {\tt ESMF\_IOFMT\_NETCDF}.
!   \item[{[relaxedflag]}]
!     If {\tt .true.}, then no error is returned even if the call cannot
!     write
!     the file due to library limitations. Default is {\tt .false.}.
!   \item[{[rc]}]
!     Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    logical                 :: ioCapable
    logical                 :: doItFlag
    character(len=64)       :: lfileName
    character(len=64)       :: gridName
    type(ESMF_Array)        :: array
    type(ESMF_ArrayBundle)  :: arraybundle
    logical                 :: isPresent
    integer                 :: dimCount
    integer                 :: dimIndex
    integer,allocatable     :: coordDimCount(:)
    integer                 :: coordDimMax
    integer                 :: stat
    logical                 :: lnclScript
    logical                 :: hasCorners

    if (present(rc)) rc = ESMF_SUCCESS

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))

    doItFlag = .true. ! default
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then

      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(grid, name=gridName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        lfileName = trim(gridName)//".nc"
      endif

      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

      ! -- centers --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      endif

      ! -- corners --

      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
          call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        endif
      endif

      ! -- mask --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      endif

      ! -- area --

      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
        call ESMF_ArrayBundleAdd(arraybundle,(/array/),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
      endif

      call ESMF_ArrayBundleWrite(arraybundle, &
        fileName=trim(lfileName),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out
    endif

  end subroutine
#undef METHOD

  !-----------------------------------------------------------------------------

#define METHOD "HYCOM_ESMF_LogGrid"
!BOP
! !IROUTINE: HYCOM_ESMF_LogGrid - Write ESMF grid information to PET logs
! !INTERFACE:
  ! call using generic interface: HYCOM_ESMF_LogGrid
  subroutine HYCOM_ESMF_LogGrid(grid,label,rc)
! ! ARGUMENTS
    type(ESMF_Grid), intent(in)            :: grid
    character(len=*), intent(in), optional :: label
    integer, intent(out), optional         :: rc
! !DESCRIPTION:
!   Write ESMF grid information to PET logs
!
!   The arguments are:
!   \begin{description}
!   \end{description}
!
!EOP
  !-----------------------------------------------------------------------------
    ! local variables
    character(len=64)           :: llabel
    character(len=64)           :: gridName
    type(ESMF_DistGrid)         :: distgrid
    character(len=64)           :: transferAction
    integer                     :: localDeCount
    integer                     :: dimCount, tileCount, deCount
    integer                     :: dimIndex, tileIndex, deIndex
    integer,allocatable         :: coordDimCount(:)
    integer                     :: coordDimMax
    integer,allocatable         :: minIndexPTile(:,:), maxIndexPTile(:,:)
    integer,allocatable         :: minIndexPDe(:,:), maxIndexPDe(:,:)
    integer                     :: stat
    character(len=ESMF_MAXSTR)  :: logMsg

    if (present(rc)) rc = ESMF_SUCCESS
    if (present(label)) then
      llabel = trim(label)
    else
      llabel = 'HYCOM_ESMF_LogGrid'
    endif

    ! access localDeCount to show this is a real Grid
    call ESMF_GridGet(grid, name=gridName, &
      localDeCount=localDeCount, distgrid=distgrid, &
      dimCount=dimCount,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

    ! allocate coordDim info accord. to dimCount and tileCount
    allocate(coordDimCount(dimCount), &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of coordinate dimensions memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

    ! get coordDim info
    call ESMF_GridGet(grid, coordDimCount=coordDimCount, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

    coordDimMax = 0
    do dimIndex=1,dimCount
      coordDimMax = MAX(coordDimMax,coordDimCount(dimIndex))
    enddo

    if (coordDimMax == 1) then
      write (logMsg,"(A,A,A)") trim(llabel)//": ", &
        trim(gridName), &
        " is a rectilinear grid with 1D coordinates in each dimension."
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    endif

    deallocate(coordDimCount, &
      stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Dellocation of coordinate dimensions memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

    write (logMsg,"(A,A,(A,I0))") trim(llabel)//": ", &
      trim(gridName), &
      " local decomposition count=",localDeCount
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

    ! get dimCount and tileCount
    call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, &
      deCount=deCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

    write (logMsg,"(A,A,(A,I0))") trim(llabel)//": ", &
      trim(gridName), &
      " dimension count=",dimCount
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,(A,I0))") trim(llabel)//": ", &
      trim(gridName), &
      " tile count=",tileCount
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A,(A,I0))") trim(llabel)//": ", &
      trim(gridName), &
      " decomp count=",deCount
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)

    ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
    allocate(minIndexPTile(dimCount, tileCount), &
      maxIndexPTile(dimCount, tileCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of index array memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

    ! get minIndex and maxIndex arrays
    call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
       maxIndexPTile=maxIndexPTile, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

    do tileIndex=1,tileCount
    do dimIndex=1,dimCount
      write (logMsg,"(A,A,A,4(I0,A))") trim(llabel)//": ", &
        trim(gridName), &
        " (tile,dim,minIndexPTile,maxIndexPTile)=(", &
        tileIndex,",",dimIndex,",", &
        minIndexPTile(dimIndex,tileIndex),",", &
        maxIndexPTile(dimIndex,tileIndex),")"
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo
    enddo

    deallocate(minIndexPTile, maxIndexPTile,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of index array memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

    ! allocate minIndexPDe and maxIndexPDe accord. to dimCount and deCount
    allocate(minIndexPDe(dimCount, deCount), &
      maxIndexPDe(dimCount, deCount),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of index array memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

    ! get minIndex and maxIndex arrays
    call ESMF_DistGridGet(distgrid, minIndexPDe=minIndexPDe, &
       maxIndexPDe=maxIndexPDe, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, PASSTHRU)) return  ! bail out

    do deIndex=1,deCount
    do dimIndex=1,dimCount
      write (logMsg,"(A,A,A,4(I0,A))") trim(llabel)//": ", &
        trim(gridName), &
        " (decomp,dim,minIndexPDe,maxIndexPDe)=(", &
        deIndex,",",dimIndex,",", &
        minIndexPDe(dimIndex,deIndex),",", &
        maxIndexPDe(dimIndex,deIndex),")"
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo
    enddo

    deallocate(minIndexPDe, maxIndexPDe,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of index array memory failed.", &
      CONTEXT, rcToReturn=rc)) return  ! bail out

  end subroutine
#undef METHOD

  !-----------------------------------------------------------------------------

end module HYCOM_ESMF_Extensions


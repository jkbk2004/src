!===============================================================================
! MODULE: HYCOM OCEAN ESMF Module
!
! DESCRIPTION:
!   This module wraps HYCOM with NUOPC interfaces.
!
! SUBROUTINES:
!   SetServices
!     Register ESMF-NUOPC entry points.
!
!   InitializeP0
!     Set initialization phase definition and read configuration attributes.
!
!   HYCOM_AttributeGet
!     Read configuration attributes from model.
!
!   HYCOM_ConfigToAttribute
!     Convert configuration object to model attributes.
!
!   InitializeP1
!     Read field configuration and advertise fields.
!
!   InitializeP2
!     Initialize HYCOM model, create ESMF grid, create ESMF fields, realize ESMF
!     fields.
!
!   ModelAdvance
!     Copy import data to HYCOM, advance HYCOM, copy HYCOM data to export.
!
!   OCEAN_Final
!     Deallocate memory and print timers.
!
!   do_export
!     Copy HYCOM data to export field.
!
!   do_import
!     Copy import field to HYCOM data.
!
!===============================================================================
#include "HYCOM_NUOPC_Macros.h"
!===============================================================================
module HYCOM_Mod
#define MODNAME "HYCOM_Mod"
!===============================================================================
! use modules
!===============================================================================
! ESMF framework modules
  use ESMF
  use NUOPC
  use HYCOM_ESMF_Extensions
  use hycom_nuopc_flags
  use NUOPC_Model, only: &
    model_routine_SS    => SetServices, &
    model_label_Advance => label_Advance, &
    model_label_Finalize => label_Finalize
! HYCOM OCEAN forecast module
  use mod_hycom, only : end_of_run_cpl, &
                        end_of_run, &
                        HYCOM_Init, &
                        HYCOM_Run,  &
                        HYCOM_Final
  use mod_archiv, only: archiv_exchange
  use hycom_couple, only: cpldom, &
                          hycom_couple_init, &
                          set_hycom_import_flag, &
                          hycom_couple_check_deb, &
                          export_from_hycom_deb, &
                          import_to_hycom_deb, &
                          ocn_import_forcing, &
                          hycom_couple_final
  use mod_import, only: hycom_imp_reset
#ifdef ESPC_COUPLE
  use read_impexp_config_mod
  use impexpField_cdf_mod
#endif

!===============================================================================
! settings
!===============================================================================
  implicit none
  private
  save

!===============================================================================
! public
!===============================================================================
  public SetServices

!===============================================================================
! module variables
!===============================================================================
  type(ESMF_VM)      :: vm
  integer            :: nPets, lPet
  integer, parameter :: localDE=0
  type(ESMF_Clock)   :: intClock
  real               :: ocean_start_dtg, ocean_end_dtg
#ifdef ESPC_COUPLE
  integer                        :: numImpFields, numExpFields
  type(ESMF_Field), allocatable  :: impField(:)
  type(ESMF_Field), allocatable  :: expField(:)
  character(len=30), pointer :: expFieldName(:), impFieldName(:) => NULL()
  character(len=60), pointer :: expStandName(:), impStandName(:) => NULL()
  character(len=30), pointer :: expFieldUnit(:), impFieldUnit(:) => NULL()
  logical, pointer           :: expFieldEnable(:), impFieldEnable(:) => NULL()
  character(len=30) :: ocn_impexp_file = "dummy_file"
  integer           :: cdf_impexp_freq
  integer           :: cpl_hour, cpl_min, cpl_sec
  real              :: cpl_time_step
  logical           :: ocn_esmf_exp_output, ocn_esmf_imp_output
  logical           :: hycom_arche_output
  character(len=15) :: base_dtg
#endif
  integer           :: end_hour, end_min, end_sec
  integer           :: start_hour, start_min, start_sec
  real              :: endtime
  integer           :: itdmx, jtdmx
  logical           :: show_minmax
  type(import_flag) :: import_setting
  logical           :: import_diagnostics
  logical           :: merge_all_import
  logical           :: skip_first_import
  integer           :: mask_field_id            = 0
  integer           :: scalar_field_id          = 0
  character(len=60) :: scalar_field_name        = 'cpl_scalars'
  integer           :: scalar_field_count       = 0
  integer           :: scalar_field_idx_grid_nx = 0
  integer           :: scalar_field_idx_grid_ny = 0
#ifdef ESPC_TIMER
  real(kind=ESMF_KIND_R8) :: timer_beg, timer_end
  real(kind=ESMF_KIND_R8) :: espc_timer(6)
! espc_timer(1): Init Phase
! espc_timer(2): Run Phase
! espc_timer(3): Final Phase
! espc_timer(4): Run Phase import
! espc_timer(5): Run Phase Core
! espc_timer(6): Run Phase export
#endif
  real(ESMF_KIND_R8),parameter :: fillValue = 9.99e20_ESMF_KIND_R8

!===============================================================================
  contains
!===============================================================================
  subroutine SetServices(model, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
!   local variables
    character(32)           :: cname
    character(*), parameter :: rname="SetServices"

    rc = ESMF_FAILURE

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=OCEAN_Final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="set final entry point failed", &
      CONTEXT)) return

    rc = ESMF_SUCCESS

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(model, importState, exportState, clock, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
!   local variables
    character(32)           :: cname
    character(*), parameter :: rname="InitializeP0"
    integer                 :: verbosity, diagnostic
    character(len=64)       :: value
    integer                 :: stat
    logical                 :: isPresent, isSet

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
!    call NUOPC_CompGet(model, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(model, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(model, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_AttributeGet(model, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"off","low","high","max","bit16","maxplus"/), &
      specialValueList=(/0,9985,32513,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv00p"/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call HYCOM_AttributeGet(rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call set_impexp_fields(cpl_scalars=(scalar_field_count.gt.0))

    contains ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine HYCOM_AttributeGet(rc)
!     arguments
      integer, intent(out) :: rc
!     local variables
      character(*), parameter :: rname="HYCOM_AttributeGet"
      logical                 :: configIsPresent
      type(ESMF_Config)       :: config
      type(NUOPC_FreeFormat)  :: attrFF
      character(len=64)       :: value
      character(ESMF_MAXSTR)  :: logMsg

      rc = ESMF_SUCCESS

      ! check model for config
      call ESMF_GridCompGet(model, configIsPresent=configIsPresent, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      ! read and ingest free format component attributes
      if (configIsPresent) then
        call ESMF_GridCompGet(model, config=config, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

#ifdef ESPC_COUPLE
        call HYCOM_ConfigToAttribute(config, &
          attrName="cdf_impexp_freq", attrDflt="9999", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        call HYCOM_ConfigToAttribute(config, &
          attrName="cpl_hour", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="cpl_min", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="cpl_sec", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        call HYCOM_ConfigToAttribute(config, &
          attrName="base_dtg", attrDflt="9999999999", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        call HYCOM_ConfigToAttribute(config, &
          attrName="hycom_arche_output", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

#if defined ESPC_OCN
        call HYCOM_ConfigToAttribute(config, &
          attrName="ocn_esmf_exp_output", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="ocn_esmf_imp_output", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="ocn_impexp_file", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
#else
        call HYCOM_ConfigToAttribute(config, &
          attrName="hyc_esmf_exp_output", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="hyc_esmf_imp_output", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="hyc_impexp_file", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
#endif
#endif

        call HYCOM_ConfigToAttribute(config, &
          attrName="espc_show_impexp_minmax", attrDflt=".true.", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

!       Start Time
        call HYCOM_ConfigToAttribute(config, &
          attrName="ocean_start_dtg", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        call HYCOM_ConfigToAttribute(config, &
          attrName="start_hour", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="start_min", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="start_sec", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

!       End Time
        call HYCOM_ConfigToAttribute(config, &
          attrName="end_hour", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="end_min", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call HYCOM_ConfigToAttribute(config, &
          attrName="end_sec", attrDflt="0", rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        attrFF = NUOPC_FreeFormatCreate(config, &
          label=trim(cname)//"_attributes::", relaxedflag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_CompAttributeIngest(model, attrFF, addFlag=.true., rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

      endif

      ! read component attributes
#ifdef ESPC_COUPLE
      call ESMF_AttributeGet(model, value=value, &
        name="cdf_impexp_freq", defaultvalue="9999", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get cdf_impexp_freq failed", CONTEXT)) return
      cdf_impexp_freq = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      call ESMF_AttributeGet(model, value=value, &
        name="cpl_hour", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get cpl_hour failed", CONTEXT)) return
      cpl_hour = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="cpl_min", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get cpl_min failed", CONTEXT)) return
      cpl_min = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="cpl_sec", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get cpl_sec failed", CONTEXT)) return
      cpl_sec = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      cpl_time_step=(cpl_hour+cpl_min/60.+cpl_sec/3600.)

      call ESMF_AttributeGet(model, value=value, &
        name="base_dtg", defaultvalue="9999999999", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get base_dtg failed", CONTEXT)) return
      base_dtg = trim(value)

      call ESMF_AttributeGet(model, value=value, &
        name="hycom_arche_output", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get hycom_arche_output failed", CONTEXT)) return
      hycom_arche_output = (trim(value)==".true.")

#if defined ESPC_OCN
      call ESMF_AttributeGet(model, value=value, &
        name="ocn_esmf_exp_output", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get ocn_esmf_output failed", CONTEXT)) return
      ocn_esmf_exp_output = (trim(value)==".true.")
      call ESMF_AttributeGet(model, value=value, &
        name="ocn_esmf_imp_output", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get ocn_esmf_output failed", CONTEXT)) return
      ocn_esmf_imp_output = (trim(value)==".true.")
      call ESMF_AttributeGet(model, value=value, &
        name="ocn_impexp_file", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get ocn_impexp_file failed", CONTEXT)) return
      ocn_impexp_file = trim(value)
#else
      call ESMF_AttributeGet(model, value=value, &
        name="hyc_esmf_exp_output", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get hyc_esmf_output failed", CONTEXT)) return
      ocn_esmf_exp_output = (trim(value)==".true.")
      call ESMF_AttributeGet(model, value=value, &
        name="hyc_esmf_imp_output", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get hyc_esmf_output failed", CONTEXT)) return
      ocn_esmf_imp_output = (trim(value)==".true.")
      call ESMF_AttributeGet(model, value=value, &
        name="hyc_impexp_file", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get hyc_impexp_file failed", CONTEXT)) return
      ocn_impexp_file = trim(value)
#endif
      if (lPet.eq.0) print *,"ocn_impexp_file=",ocn_impexp_file
#endif

      call ESMF_AttributeGet(model, value=value, &
        name="espc_show_impexp_minmax", defaultvalue=".true.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get espc_show_impexp_minmax failed", CONTEXT)) return
      show_minmax = (trim(value)==".true.")

      call ESMF_AttributeGet(model, value=value, &
        name="import_setting", defaultvalue="FLEXIBLE", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get import_setting failed", CONTEXT)) return
      import_setting=value

      call ESMF_AttributeGet(model, value=value, &
        name="import_diagnostics", defaultvalue=".false.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get import_diagnostics failed", CONTEXT)) return
      import_diagnostics = (trim(value)==".true.")

      call ESMF_AttributeGet(model, value=value, &
        name="merge_all_import", defaultvalue=".false.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get merge_all_import failed", CONTEXT)) return
      merge_all_import = (trim(value)==".true.")

      call ESMF_AttributeGet(model, value=value, &
        name="skip_first_import", defaultvalue=".false.", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get skip_first_import failed", CONTEXT)) return
      skip_first_import = (trim(value)==".true.")

!     start Time
      call ESMF_AttributeGet(model, value=value, &
        name="ocean_start_dtg", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get start_dtg failed", CONTEXT)) return
      ocean_start_dtg = ESMF_UtilString2Real(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="start_hour", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get start_hour failed", CONTEXT)) return
      start_hour = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="start_min", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get start_min failed", CONTEXT)) return
      start_min = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="start_sec", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get start_sec failed", CONTEXT)) return
      start_sec = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

!     end Time
      call ESMF_AttributeGet(model, value=value, &
        name="end_hour", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get end_hour failed", CONTEXT)) return
      end_hour = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="end_min", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get end_min failed", CONTEXT)) return
      end_min = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_AttributeGet(model, value=value, &
        name="end_sec", defaultvalue="0", &
        convention="NUOPC", purpose="Instance", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="attribute get end_sec failed", CONTEXT)) return
      end_sec = ESMF_UtilString2Int(value, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      ! get scalar field name
      call NUOPC_CompAttributeGet(model, name="ScalarFieldName", &
        isPresent=isPresent, isSet=isSet, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (isPresent .and. isSet) then
        call NUOPC_CompAttributeGet(model, name="ScalarFieldName", &
          value=value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        scalar_field_name = trim(value)
        if (len_trim(scalar_field_name).le.0) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg=rname//": ERROR ScalarFieldName cannot be blank", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
      else
        scalar_field_name = "cpl_scalars"
      end if

      ! get scalar field count
      call NUOPC_CompAttributeGet(model, name="ScalarFieldCount", &
        isPresent=isPresent, isSet=isSet, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (isPresent .and. isSet) then
        call NUOPC_CompAttributeGet(model, name="ScalarFieldCount", &
          value=value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        scalar_field_count = ESMF_UtilString2Int(value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      else
        scalar_field_count = 0
      end if

      ! get scalar field id GridNX
      call NUOPC_CompAttributeGet(model, name="ScalarFieldIdxGridNX", &
        isPresent=isPresent, isSet=isSet, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (isPresent .and. isSet) then
        call NUOPC_CompAttributeGet(model, name="ScalarFieldIdxGridNX", &
          value=value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        scalar_field_idx_grid_nx = ESMF_UtilString2Int(value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      else
        scalar_field_idx_grid_nx = 0
      end if

      ! get scalar field id GridNY
      call NUOPC_CompAttributeGet(model, name="ScalarFieldIdxGridNY", &
        isPresent=isPresent, isSet=isSet, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      if (isPresent .and. isSet) then
        call NUOPC_CompAttributeGet(model, name="ScalarFieldIdxGridNY", &
          value=value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
        scalar_field_idx_grid_ny = ESMF_UtilString2Int(value, rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      else
        scalar_field_idx_grid_ny = 0
      end if

!     log configuration settings
      if (btest(verbosity,16)) then
        call ESMF_LogWrite(trim(cname)//": Settings",ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'Verbosity               = ',verbosity
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'Diagnostic              = ',diagnostic
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'cdf_impexp_freq         = ',cdf_impexp_freq
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'cpl_hour                = ',cpl_hour
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'cpl_min                 = ',cpl_min
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'cpl_sec                 = ',cpl_sec
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'base_dtg                = ',base_dtg
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'hycom_arche_output      = ',hycom_arche_output
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'hyc_esmf_exp_output     = ',ocn_esmf_exp_output
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'hyc_esmf_imp_output     = ',ocn_esmf_imp_output
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'hyc_impexp_file         = ',ocn_impexp_file
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'espc_show_impexp_minmax = ',show_minmax
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        value=import_setting
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'import_setting          = ',trim(value)
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'import_diagnostics      = ',import_diagnostics
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'merge_all_import        = ',merge_all_import
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,L1))") trim(cname)//': ', &
          'skip_first_import       = ',skip_first_import
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,F0.1))") trim(cname)//': ', &
          'ocean_start_dtg         = ',ocean_start_dtg
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'start_hour              = ',start_hour
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'start_min               = ',start_min
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'start_sec               = ',start_sec
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'end_hour                = ',end_hour
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'end_min                 = ',end_min
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'end_sec                 = ',end_sec
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,A))") trim(cname)//': ', &
          'ScalarFieldName         = ',trim(scalar_field_name)
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'ScalarFieldCount        = ',scalar_field_count
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'ScalarFieldIdxGridNX    = ',scalar_field_idx_grid_nx
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
        write (logMsg, "(A,(A,I0))") trim(cname)//': ', &
          'ScalarFieldIdxGridNY    = ',scalar_field_idx_grid_ny
        call ESMF_LogWrite(trim(logMsg),ESMF_LOGMSG_INFO)
      endif

    end subroutine HYCOM_AttributeGet

    !---------------------------------------------------------------------------

    subroutine HYCOM_ConfigToAttribute(config, attrName, attrDflt, rc)
!     arguments
      type(ESMF_Config), intent(inout)       :: config
      character(len=*), intent(in)           :: attrName
      character(len=*), intent(in), optional :: attrDflt
      integer, intent(out)                   :: rc
!     local variables
      character(*), parameter :: rname="HYCOM_ConfigToAttribute"
      character(len=64)       :: value

      call ESMF_ConfigGetAttribute(config, value, &
        label=trim(attrName)//"=", default=attrDflt, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="config get "//trim(attrName)//" failed", CONTEXT)) return
      call NUOPC_CompAttributeAdd(model, attrList=(/trim(attrName)/), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call NUOPC_CompAttributeSet(model, value=value, &
        name=trim(attrName), rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

    end subroutine HYCOM_ConfigToAttribute

  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
!   local variables
    character(32)           :: cname
    character(*), parameter :: rname="InitializeP1"
    integer                 :: i, irc

    rc = ESMF_SUCCESS

!   Get VM info
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="GridCompGet failed", &
      CONTEXT)) return

    call ESMF_VMGet(vm, petCount=nPets,localPet=lPet,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="VMGet failed", CONTEXT)) return

    if (lPet.eq.0) print *,"hycom, InitializeP1 called,nPets=",nPets

#ifdef ESPC_TIMER
    do i=1,6
      espc_timer(i)=0.
    enddo

    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

#ifdef ESPC_COUPLE
    if ((ocn_impexp_file.eq."dummy_file").or. &
        (ocn_impexp_file.eq."none      ").or. &
        (ocn_impexp_file.eq."          ")) then
      call read_impexp_config(numExpFields,numImpFields,&
        expFieldName,impFieldName,expStandName,impStandName,expFieldUnit,&
        impFieldUnit,expFieldEnable,impFieldEnable,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    else
      call read_impexp_config(ocn_impexp_file,numExpFields,numImpFields,&
        expFieldName,impFieldName,expStandName,impStandName,expFieldUnit,&
        impFieldUnit,expFieldEnable,impFieldEnable,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    allocate(impField(numImpFields))
    allocate(expField(numExpFields))

    if (lPet.eq.0) print *,"hycom,expFieldName=", &
      (expFieldName(i),i=1,numExpFields)
    if (lPet.eq.0) print *,"hycom,expFieldEnable=", &
      (expFieldEnable(i),i=1,numExpFields)
    if (lPet.eq.0) print *,"hycom,impFieldName=", &
      (impFieldName(i),i=1,numImpFields)
    if (lPet.eq.0) print *,"hycom,impFieldEnable=", &
      (impFieldEnable(i),i=1,numImpFields)

    do i=1,numImpFields
      if(import_setting.eq.IMPORT_UNCOUPLED) then
        ! disable field if import_setting=UNCOUPLED
        impFieldEnable(i)=.false.
      else
        if (impFieldEnable(i)) then
          if (lPet.eq.0) print *,"hycom,import field advertised, name=", &
            impFieldName(i),impStandName(i)
          call NUOPC_Advertise(importState, name=impFieldName(i), &
            StandardName=impStandName(i), rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
        endif
      endif
    enddo

    do i=1,numExpFields
      if (expFieldEnable(i)) then
        if (lPet.eq.0) print *,"hycom,export field advertised, name=", &
          expFieldName(i),expStandName(i)

        call NUOPC_Advertise(exportState, name=expFieldName(i), &
        StandardName=expStandName(i), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return
      endif
    enddo
#endif

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(1)=timer_end-timer_beg
!   if (lPet.eq.0) print *,"hycom, InitializeP1, timer=",espc_timer(1)
    call print_timer_stat('hycom, Init1:',timer_end-timer_beg,lPet,nPets,vm,rc)
#endif

    if (lPet.eq.0) print *,"hycom, InitializeP1 end called..."

  end subroutine InitializeP1

  !-----------------------------------------------------------------------------

  subroutine InitializeP2(model, importState, exportState, clock, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
!   local variables
    character(32)                   :: cname
    character(*), parameter         :: rname="InitializeP2"
    integer                         :: verbosity, diagnostic
    character(len=64)               :: value
    type(ESMF_Grid)                 :: gridIn
    type(ESMF_Grid)                 :: gridOut
#ifdef ESPC_COUPLE
    integer                         :: tlb(2), tub(2)
    integer(ESMF_KIND_I4), pointer  :: iptr(:,:)
    real(ESMF_KIND_R8), pointer     :: fptr(:,:)
#endif
    integer                         :: i, j
    integer                         :: decomp(2)
    real*8                          :: h_start_dtg, h_end_dtg
    integer                         :: mpiCommunicator
    type(ESMF_DistGrid)             :: ocnDistGrid
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
#ifdef ESPC_COUPLE
    type(ESMF_Field)                :: lon_field, lat_field, mask_field
    type(ESMF_Field)                :: lon_corner_field, lat_corner_field
    real(ESMF_KIND_R8), pointer     :: lon_data(:,:), lat_data(:,:)
    real(ESMF_KIND_R8), pointer     :: mask_data(:,:)
    real(ESMF_KIND_R8), pointer     :: lon_corner_data(:,:)
    real(ESMF_KIND_R8), pointer     :: lat_corner_data(:,:)
    type(ESMF_ArraySpec)            :: arraySpec2Dr
    real(ESMF_KIND_R8), allocatable :: tmp_e(:,:)
    real(ESMF_KIND_R8), allocatable :: tmp_c(:,:)
    real(ESMF_KIND_R8), pointer     :: fldmsk_fptr(:,:)
    integer(ESMF_KIND_I4), pointer  :: msk_fptr(:,:)
    logical                         :: isConnected
    integer                         :: status
#endif

    rc = ESMF_SUCCESS

    ! Query component for name, verbosity, and diagnostic values
!    call NUOPC_CompGet(model, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(model, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(model, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(model, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"off","low","high","max","bit16","maxplus"/), & 
      specialValueList=(/0,9985,32513,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

    if (lPet.eq.0) print *,"hycom, InitializeP2 called"

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

    if (ocean_start_dtg.lt.0.0) then
!     start from rest
      ocean_start_dtg=-ocean_start_dtg
      h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg
      h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
      h_start_dtg=-h_start_dtg
    else
!     normal start
      h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg
      h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
    endif

!    h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
!    h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg

    if (lPet.eq.0) then
      print *,"HYCOM_OceanComp, HYCOM_Init called..."
      print *,"end_hour,end_min,end_sec=",end_hour,end_min,end_sec
      print *,"h_start_dtg,h_end_dtg=",h_start_dtg,h_end_dtg
    endif

    call ESMF_VMGet(vm, mpiCommunicator=mpiCommunicator,rc=rc)

!   call into ocean init
    call HYCOM_Init(mpiCommunicator,h_start_dtg,h_end_dtg)

    call hycom_couple_init(nPets,import_diagnostics,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="hycom_couple_init failed", &
      CONTEXT)) return

    call hycom_couple_check_deb(show_minmax, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="hycom_couple_check_deb failed", &
      CONTEXT)) return

    if (lPet.eq.0) print *,"hycom, InitializeP2 called 2,idim_size,jdim_size", &
       cpldom%idim_size,cpldom%jdim_size

#if defined(ARCTIC)
    itdmx=cpldom%idim_size; jtdmx=cpldom%jdim_size-1
#else
    itdmx=cpldom%idim_size; jtdmx=cpldom%jdim_size
#endif

    if (lPet.eq.0) print *,"hycom, itdmx,jtdmx..=",itdmx,jtdmx

#ifdef ESPC_COUPLE
    if (lPet.eq.0) then
      allocate(tmp_e(itdmx,jtdmx))
      allocate(tmp_c(itdmx+1,jtdmx+1))
    else
      allocate(tmp_e(1,1))
      allocate(tmp_c(1,1))
    endif
#endif

#ifdef ESPC_COUPLE
    allocate(connectionList(1)) ! one connection
    call ESMF_DistGridConnectionSet(connection=connectionList(1), &
    tileIndexA=1, tileIndexB=1, positionVector=(/itdmx, 0/), rc=rc)

!   ocnDistGrid = ESMF_DistGridCreate(minIndex=(/1,1/), &
!     maxIndex=(/itdmx,jtdmx/), indexflag=ESMF_INDEX_GLOBAL, &
!     deBlockList=cpldom%deBList,connectionList=connectionList,rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", &
!     CONTEXT)) return

    ! create new grid
    ocnDistGrid = ESMF_DistGridCreate(minIndex=(/1,1/), &
      maxIndex=(/itdmx,jtdmx/), indexflag=ESMF_INDEX_GLOBAL, &
      deBlockList=cpldom%deBList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", &
      CONTEXT)) return

    gridIn = ESMF_GridCreate(distGrid=ocnDistGrid, &
      indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, &
      name="OCEAN:grid", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", &
      CONTEXT)) return

    if (allocated(connectionList))deallocate(connectionList)
#else
    do i=int(sqrt(real(nPets))),2,-1
      if (mod(nPets,i).eq.0) exit
    enddo
    decomp=(/nPets/i,i/)
    if (lPet.eq.0) print *,"nPets,decomp=",nPets,decomp

    gridIn = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
      maxIndex=(/cpldom%idim_size,cpldom%jdim_size/), &
      indexflag=ESMF_INDEX_GLOBAL, &
      regDecomp=decomp, name="OCEAN:grid", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", &
      CONTEXT)) return
#endif

#ifdef ESPC_COUPLE
!   add ESMF grid coordinate arrays
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid add coord failed", &
      CONTEXT)) return

!   add ESMF grid corner coordinate arrays
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid add corner coord failed", &
      CONTEXT)) return

!   add ESMF mask array
    call ESMF_GridAddItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid add mask failed", &
      CONTEXT)) return

    if ((ocn_esmf_imp_output.or.ocn_esmf_exp_output).and.lPet.eq.0) then
      call impexp_cdf_put_latlonmsk('hycom',itdmx,jtdmx,cpldom%lat_p, &
        cpldom%lon_p,cpldom%mask_p,status,rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="impexp_cdf_put_latlonmsk failed", CONTEXT)) return

      call impexp_cdf_put_latlonmsk('hycom-orig',cpldom%idim_size,&
        cpldom%jdim_size,cpldom%lat_p,cpldom%lon_p,cpldom%mask_p,status,rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="impexp_cdf_put_latlonmsk failed", CONTEXT)) return

      call impexp_cdf_put_latlonmsk_corner('hycom',cpldom%idim_size, &
        cpldom%jdim_size,cpldom%lat_p,cpldom%lon_p,int(cpldom%mask_p), &
        cpldom%lat_q,cpldom%lon_q,status,rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="impexp_cdf_put_latlonmsk_corner failed", CONTEXT)) return
    endif

    call ESMF_ArraySpecSet(arraySpec2Dr, rank=2, typekind=ESMF_TYPEKIND_R8, &
      rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    lon_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
      indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      name=TRIM("lon_field"), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (lPet.eq.0) then
      do i=1,itdmx
      do j=1,jtdmx
        tmp_e(i,j)=cpldom%lon_p(i,j)
      enddo
     enddo
    endif

    call ESMF_FieldScatter(lon_field,tmp_e,rootPet=0,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    lat_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
      indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      name=TRIM("lat_field"), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (lPet.eq.0) then
      do i=1,itdmx
      do j=1,jtdmx
        tmp_e(i,j)=cpldom%lat_p(i,j)
      enddo
      enddo
    endif

    call ESMF_FieldScatter(lat_field,tmp_e,rootPet=0,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

!   corner
    lon_corner_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
      indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      name=TRIM("lon_corner_field"), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (lPet.eq.0) then
      do i=1,itdmx
      do j=1,jtdmx
        tmp_c(i,j)=cpldom%lon_q(i,j)
      enddo
      enddo

!     fill edges
      do i=1,itdmx
        tmp_c(i,jtdmx+1)=tmp_c(i,jtdmx)+(cpldom%lon_q(i,jtdmx)- &
          cpldom%lon_q(i,jtdmx-1))
      end do
      do j=1,jtdmx
        tmp_c(itdmx+1,j)=tmp_c(itdmx,j)+(cpldom%lon_q(itdmx,j)- &
          cpldom%lon_q(itdmx-1,j))
      end do
      tmp_c(itdmx+1,jtdmx+1)=tmp_c(itdmx,jtdmx)+(tmp_c(itdmx,jtdmx)- &
        tmp_c(itdmx-1,jtdmx-1))
    endif

    call ESMF_FieldScatter(lon_corner_field,tmp_c,rootPet=0,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    lat_corner_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
      indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      name=TRIM("lat_corner_field"), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (lPet.eq.0) then
      do i=1,itdmx
      do j=1,jtdmx
        tmp_c(i,j)=cpldom%lat_q(i,j)
      enddo
      enddo

!     fill edges
      do i=1,itdmx
        tmp_c(i,jtdmx+1)=tmp_c(i,jtdmx)+(cpldom%lat_q(i,jtdmx)- &
          cpldom%lat_q(i,jtdmx-1))
      end do
      do j=1,jtdmx
        tmp_c(itdmx+1,j)=tmp_c(itdmx,j)+(cpldom%lat_q(itdmx,j)- &
          cpldom%lat_q(itdmx-1,j))
      end do
      tmp_c(itdmx+1,jtdmx+1)=tmp_c(itdmx,jtdmx)+(tmp_c(itdmx,jtdmx)- &
        tmp_c(itdmx-1,jtdmx-1))
    endif

    call ESMF_FieldScatter(lat_corner_field,tmp_c,rootPet=0,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    mask_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
      indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      name=TRIM("mask_field"), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    if (lPet.eq.0) then
      do i=1,itdmx
      do j=1,jtdmx
        tmp_e(i,j)=cpldom%mask_p(i,j)
      enddo
      enddo
    endif

    call ESMF_FieldScatter(mask_field,tmp_e,rootPet=0,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(lon_field,localDE,lon_data,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(lat_field,localDE,lat_data,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(lon_corner_field,localDE,lon_corner_data,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(lat_corner_field,localDE,lat_corner_data,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(mask_field,localDE,mask_data,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

!   copy in coordinate data
    call ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=1, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid get coord 1 failed", &
      CONTEXT)) return
!   fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_e(tlb(1):tub(1),tlb(2):tub(2))
    fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_data(tlb(1):tub(1),tlb(2):tub(2))

    call ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=2, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid get coord 2 failed", &
      CONTEXT)) return
!   fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_e(tlb(1):tub(1),tlb(2):tub(2))
    fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_data(tlb(1):tub(1),tlb(2):tub(2))

!   copy in corner coordinate data
    call ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=1, &
      staggerLoc=ESMF_STAGGERLOC_CORNER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid get coord 1 failed", &
      CONTEXT)) return
!   fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_e(tlb(1):tub(1),tlb(2):tub(2))
    fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_corner_data(tlb(1):tub(1), &
      tlb(2):tub(2))

    call ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=2, &
      staggerLoc=ESMF_STAGGERLOC_CORNER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid get corner coord 2 failed", &
      CONTEXT)) return
!   fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_e(tlb(1):tub(1),tlb(2):tub(2))
    fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_corner_data(tlb(1):tub(1), &
      tlb(2):tub(2))

!   copy in land/sea mask (integer)
    call ESMF_GridGetItem(gridIn, localDE=localDE, &
      itemflag=ESMF_GRIDITEM_MASK, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      totalLBound=tlb, totalUBound=tub, farrayPtr=iptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="grid get mask failed", &
      CONTEXT)) return
!   iptr(tlb(1):tub(1),tlb(2):tub(2)) = NINT(mask_e(tlb(1):tub(1), &
!     tlb(2):tub(2)))
    iptr(tlb(1):tub(1),tlb(2):tub(2)) = NINT(mask_data(tlb(1):tub(1), &
      tlb(2):tub(2)))

    call ESMF_FieldDestroy(lon_field, rc=rc)
    call ESMF_FieldDestroy(lat_field, rc=rc)
    call ESMF_FieldDestroy(lon_corner_field, rc=rc)
    call ESMF_FieldDestroy(lat_corner_field, rc=rc)
    call ESMF_FieldDestroy(mask_field, rc=rc)
#endif

    gridOut = gridIn ! for now out same as in

    ! log grid and decomp to PET logs
    if (btest(verbosity,16)) then
      call HYCOM_ESMF_LogGrid(gridOut, trim(cname)//"_"//rname,rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    endif

    ! write grid to NetCDF file.
    if (btest(diagnostic,16)) then
      call HYCOM_ESMF_GridWrite(gridOut, "diagnostic_"//trim(cname)//"_"// &
        rname//"_grid.nc", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

#ifdef ESPC_COUPLE
    do i=1,numImpFields
      if (impFieldEnable(i)) then
        isConnected = NUOPC_IsConnected(importState, &
          fieldName=impFieldName(i),rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (isConnected) then
          if (lPet.eq.0) print *,"hycom, import field created, name=", &
            impFieldName(i)
          impField(i) = ESMF_FieldCreate(name=impFieldName(i), grid=gridIn, &
            typekind=ESMF_TYPEKIND_RX, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          call ESMF_FieldFill(impField(i), dataFillScheme="const", &
            const1=fillValue, rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
          ! realize field in import state
          call NUOPC_Realize(importState, field=impField(i), rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        else
          if(import_setting.eq.IMPORT_FLEXIBLE) then
            ! remove if field is not connected and import_setting=FLEXIBLE
            if (lPet.eq.0) print *,"hycom, import field disabled, name=", &
              impFieldName(i)
            impFieldEnable(i) = .false.
            call ESMF_StateRemove(importState, (/impFieldName(i)/), &
              relaxedflag=.true., rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          else
            ! fail if field is not connected and import_setting=REQUIRED
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
              msg="ERROR: hycom, import field required: "//trim(impFieldName(i)), &
              CONTEXT, rcToReturn=rc)
            return ! bail out
          endif
        endif
      endif
      if (impFieldEnable(i)) then
        call set_hycom_import_flag(impFieldName(i),rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, &
          msg="set_hycom_import_flag failed", CONTEXT)) return
      endif
    enddo

    do i=1,numExpFields
      if (expFieldEnable(i)) then
        isConnected = NUOPC_IsConnected(exportState, &
          fieldName=expFieldName(i), rc=rc)
        if (ESMF_STDERRORCHECK(rc)) return

        if (isConnected) then
          if (lPet.eq.0) print *,"hycom, export field created, name=", &
            expFieldName(i)

          if (expFieldName(i).eq.scalar_field_name) then
            expField(i) = CreateScalarField(name=expFieldName(i), &
              field_count=scalar_field_count, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            call SetScalarFieldValues(expField(i), &
              vals=(/real(itdmx,ESMF_KIND_R8),real(jtdmx,ESMF_KIND_R8)/), &
              idxs=(/scalar_field_idx_grid_nx,scalar_field_idx_grid_ny/), &
              rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            scalar_field_id = i
          elseif(expFieldName(i).eq.'mask') then
            expField(i) = ESMF_FieldCreate(name=expFieldName(i), grid=gridOut, &
              typekind=ESMF_TYPEKIND_RX, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            call SetMaskFieldValues(expField(i), grid=gridOut, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            mask_field_id = i
          else
            expField(i) = ESMF_FieldCreate(name=expFieldName(i), grid=gridOut, &
              typekind=ESMF_TYPEKIND_RX, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            call ESMF_FieldFill(expField(i), dataFillScheme="const", &
              const1=0.D0, rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
            call do_export(expField(i),expFieldName(i),rc=rc)
            if (ESMF_STDERRORCHECK(rc)) return
          end if
          ! realize field in export state
          call NUOPC_Realize(exportState, field=expField(i), rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        else
          if (lPet.eq.0) print *,"hycom, export field disabled, name=", &
            expFieldName(i)
          expFieldEnable(i) = .false.
          call ESMF_StateRemove(exportState,(/expFieldName(i)/), &
            relaxedflag=.true., rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif
      endif
    enddo

    endtime=0
    if (ocn_esmf_exp_output.and.(mod(endtime,float(cdf_impexp_freq)).eq.0)) then
      call impexp_cdf_put_flds('hycom', base_dtg, endtime, itdmx,jtdmx, &
        numExpFields,expFieldEnable,expFieldName,expStandName,expFieldUnit, &
        expField,status,lPet,rc,'exp')
      if (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", &
        CONTEXT)) return
    endif

    deallocate(tmp_e)
    deallocate(tmp_c)
#endif

    call hycom_couple_final(rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(1)=espc_timer(1)+timer_end-timer_beg
!   if (lPet.eq.0) print *,"hycom, InitializeP2,timer=",espc_timer(1), &
!     timer_end-timer_beg
    call print_timer_stat('hycom, Init2:',timer_end-timer_beg,lPet,nPets,vm,rc)
#endif

    if (lPet.eq.0) print *,"hycom, InitializeP2 end called..."

  end subroutine InitializeP2

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
!   local variables
    character(32)           :: cname
    character(*), parameter :: rname="ModelAdvance"
    integer                 :: verbosity, diagnostic
    character(len=64)       :: value
    type(ESMF_Clock)        :: clock
    type(ESMF_State)        :: importState, exportState
    integer                 :: i, status
    type(ESMF_Time)         :: extCurrTime
    type(ESMF_Time)         :: extRefTime
    type(ESMF_TimeInterval) :: extTimeStep
    character(ESMF_MAXSTR)  :: currtimeString, reftimeString
    type(ESMF_TimeInterval) :: extTimeSinceStart
    real(kind=ESMF_KIND_R8) :: extSecSinceStarti
    real(kind=ESMF_KIND_R8) :: extSecTimeStep
    real                    :: begtime
    real*8                  :: endtime8
    real                    :: endtimex
#ifdef ESPC_TIMER
    real(kind=ESMF_KIND_R8) :: timer_tmp_beg, timer_tmp_end
#endif

    rc = ESMF_SUCCESS

    if (lPet.eq.0) print *,"hycom, ModelAdvance called"

    ! Query component for name, verbosity, and diagnostic values
!    call NUOPC_CompGet(model, name=name, verbosity=verbosity, &
!      diagnostic=diagnostic, rc=rc)
    call ESMF_GridCompGet(model, name=cname, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(model, name="Diagnostic", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    diagnostic = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max","bit16","maxplus"/), &
      specialValueList=(/0,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    call ESMF_AttributeGet(model, name="Verbosity", value=value, &
      defaultValue="0", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out
    verbosity = ESMF_UtilString2Int(value, &
      specialStringList=(/"off","low","high","max","bit16","maxplus"/), &
      specialValueList=(/0,9985,32513,65535,65536,131071/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return  ! bail out

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

!   compute the step count from the external clock
    call ESMF_ClockGet(clock, currTime=extCurrTime, refTime=extRefTime, &
      timeStep=extTimeStep,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Get extCurrTime failed", &
      CONTEXT)) return

    call ESMF_TimeGet(extCurrTime,timeString=currtimeString)
    call ESMF_TimeGet(extRefTime,timeString=reftimeString)

!   if (lPet.eq.0) print *,"hycom,extCurrTime, extRefTime=",currtimeString, &
!     reftimeString

!   elapsed time in seconds for internal clock
    extTimeSinceStart = extCurrTime - extRefTime
    call ESMF_TimeIntervalGet(extTimeSinceStart, s_r8=extSecSinceStarti, rc=rc)
    call ESMF_TimeIntervalGet(extTimeStep, s_r8=extSecTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Get time interval failed", &
      CONTEXT)) return

    if (lPet.eq.0) print *,"hycom,extSecTimeStep=",extSecTimeStep

    begtime=extSecSinceStarti/3600.
    endtime=(extSecSinceStarti+extSecTimeStep)/3600.

    if (lPet.eq.0) print *,"hycom,begtime,endtime=",begtime,endtime

    endtime=endtime+ocean_start_dtg*24

!   run atmos forward

    if (lPet.eq.0) print *,"HYCOM_OceanCom, ModelAdvance,"// &
      " Run ocean forward...,endtime=",endtime/24

#ifdef ESPC_COUPLE
    endtimex=endtime-ocean_start_dtg*24

    if (ocn_esmf_imp_output.and.((mod(endtimex,float(cdf_impexp_freq)).eq.0 &
      .or. endtimex.eq.0.5))) then

      call impexp_cdf_put_flds('hycom', base_dtg, endtimex-cpl_time_step, &
        itdmx,jtdmx,numImpFields,impFieldEnable,impFieldName,impStandName, &
        impFieldUnit,impField,status,lPet,rc,'imp')
      if (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", &
        CONTEXT)) return
    endif

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

    if (skip_first_import) then
      call hycom_imp_reset(.false.)
      skip_first_import=.false.
    else
      call hycom_imp_reset(merge_all_import)
!     transform and copy each field to internal import arrays
!     data will be copied from import arrays to hycom after
!     forcing data is read
      do i=1,numImpFields
        if (impFieldEnable(i)) then
          call do_import(impField(i),impFieldName(i),diagnostic,rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        endif
      enddo
!     calculate radflx and wndspd
      call ocn_import_forcing(fillValue, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

!     arche file of fields exchanged with ice component
      if (hycom_arche_output .and. begtime.eq.0) &
        call archiv_exchange  !arche file of fields exchanged with ice component
    endif

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(4)=espc_timer(4)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Import):',timer_tmp_end-timer_tmp_beg, &
      lPet,nPets,vm,rc)
#endif
#endif

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in InitializeP2()
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    if (lPet.eq.0) then
      call ESMF_ClockPrint(clock, options="currTime",&
        preString="------>Advancing OCN from: ", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return

      call ESMF_ClockPrint(clock, options="stopTime",&
        preString="--------------------------------> to: ", rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
    endif

    endtime8=endtime

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

    do !until end of run
      call HYCOM_Run(endtime8/24)
      if (end_of_run .or. end_of_run_cpl ) then
        exit
      endif
    enddo

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(5)=espc_timer(5)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Core):',timer_tmp_end-timer_tmp_beg, &
      lPet,nPets,vm,rc)
#endif

#ifdef ESPC_COUPLE
#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

    do i=1,numExpFields
      if (expFieldEnable(i)) then
        if (i.eq.scalar_field_id) then
!         SetScalarFieldValues called during InitializeP2
        elseif (i.eq.mask_field_id) then
!         SetMaskFieldValues called during InitializeP2
        else
          call do_export(expField(i),expFieldName(i),rc=rc)
          if (ESMF_STDERRORCHECK(rc)) return
        end if
      endif
    enddo

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(6)=espc_timer(6)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Export):',timer_tmp_end-timer_tmp_beg, &
      lPet,nPets,vm,rc)
#endif

!move to mod_hycom.F
#ifdef ESPC_COUPLE
    if (hycom_arche_output) call archiv_exchange ! arch fields exchanged w ice
#endif

    if (ocn_esmf_exp_output.and.( (mod(endtimex,float(cdf_impexp_freq)).eq.0 &
      .or. endtimex.eq.0.5)  )) then
      call impexp_cdf_put_flds('hycom', base_dtg,endtimex, &
        itdmx,jtdmx,numExpFields,expFieldEnable,expFieldName,expStandName, &
        expFieldUnit,expField,status,lPet,rc,'exp')
      if (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", &
        CONTEXT)) return
    endif
#endif

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(2)=espc_timer(2)+timer_end-timer_beg
!   if (lPet.eq.0) print *,"hycom, ModelAdvance,timer=",espc_timer(2), &
!     timer_end-timer_beg
    call print_timer_stat('hycom, Run:',timer_end-timer_beg,lPet,nPets,vm,rc)
#endif

    if (lPet.eq.0) print *,"hycom, ModelAdvance end..."
    if (lPet.eq.0) print *,"hycom, ModelAdvance end...",begtime,endtime
  end subroutine ModelAdvance

  !-----------------------------------------------------------------------------

  subroutine OCEAN_Final(model, rc)
!   arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
!   local variables
    character(32)                   :: cname
    character(*), parameter         :: rname="OCEAN_Final"
    integer                         :: lrc, i
#ifdef ESPC_TIMER
    integer                         :: j, ij
    real(ESMF_KIND_R8), allocatable :: espc_all_timer(:)
    real, allocatable               :: all_timer(:,:)
    real                            :: timer_min(6), timer_max(6)
    real                            :: timer_mean(6), timer_stdev(6)
#endif

    rc = ESMF_FAILURE

    call ESMF_LogWrite("HYCOM finalize routine called", ESMF_LOGMSG_INFO)

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
#endif

!   finalize ocean
    if (lPet.eq.0) print *,"HYCOM Final called.."
    call HYCOM_Final

#ifdef ESPC_COUPLE
    if (associated(expFieldName)) deallocate(expFieldName)
    if (associated(impFieldName)) deallocate(impFieldName)
    if (associated(expStandName)) deallocate(expStandName)
    if (associated(impStandName)) deallocate(impStandName)
    if (associated(expFieldUnit)) deallocate(expFieldUnit)
    if (associated(impFieldUnit)) deallocate(impFieldUnit)
    !if (associated(expFieldAddOffset)) deallocate(expFieldAddOffset)
    !if (associated(impFieldAddOffset)) deallocate(impFieldAddOffset)
    !if (associated(expFieldScaleFac)) deallocate(expFieldScaleFac)
    !if (associated(impFieldScaleFac)) deallocate(impFieldScaleFac)
    if (associated(expFieldEnable)) deallocate(expFieldEnable)
    if (associated(impFieldEnable)) deallocate(impFieldEnable)

    deallocate(impField)
    deallocate(expField)
#endif

#ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_VMWtime(timer_end, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    espc_timer(3)=timer_end-timer_beg

!   if (lPet.eq.0) print *,"hycom, Final,timer=",espc_timer(3)
    call print_timer_stat('hycom, Final:',timer_end-timer_beg,lPet,nPets,vm,rc)

    call print_timer_stat('       HYCOM_Init Phase:',espc_timer(1),lPet, &
      nPets,vm,rc)
    call print_timer_stat('        HYCOM_Run Phase:',espc_timer(2),lPet, &
      nPets,vm,rc)
    call print_timer_stat('      HYCOM_Final Phase:',espc_timer(3),lPet, &
      nPets,vm,rc)
    call print_timer_stat('HYCOM_Run Phase(Import):',espc_timer(4),lPet, &
      nPets,vm,rc)
    call print_timer_stat('  HYCOM_Run Phase(Core):',espc_timer(5),lPet, &
      nPets,vm,rc)
    call print_timer_stat('HYCOM_Run Phase(Export):',espc_timer(6),lPet, &
      nPets,vm,rc)
#endif

    if (lPet.eq.0) print *,"hycom, OCEAN_Final end..."

    rc = ESMF_SUCCESS

  end subroutine OCEAN_Final

  !-----------------------------------------------------------------------------

  ! Set mask data for field
  subroutine SetMaskFieldValues(field, grid, rc)
    type(ESMF_Field),   intent(inout) :: field
    type(ESMF_Grid),    intent(in)    :: grid
    integer,            intent(inout) :: rc
    ! local variables
    character(len=*), parameter     :: rname="SetMaskFieldValues"
    real(ESMF_KIND_R8), pointer     :: fldmsk_fptr(:,:)
    integer(ESMF_KIND_I4), pointer  :: msk_fptr(:,:)

    nullify(msk_fptr)
    nullify(fldmsk_fptr)
    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=msk_fptr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_FieldGet(field,farrayPtr=fldmsk_fptr,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (size(fldmsk_fptr).eq.size(msk_fptr)) then
      fldmsk_fptr(:,:) = msk_fptr(:,:)
    else
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="ERROR: Mask array sizes do not match.", &
        CONTEXT, rcToReturn=rc)
      return ! bail out
    endif

  end subroutine SetMaskFieldValues

  !-----------------------------------------------------------------------------

  ! Set scalar data for field
  subroutine SetScalarFieldValues(field, vals, idxs, rc)
    type(ESMF_Field),   intent(inout) :: field
    real(ESMF_KIND_R8), intent(in)    :: vals(:)
    integer,            intent(in)    :: idxs(:)
    integer,            intent(inout) :: rc

    ! local variables
    integer                         :: ungriddedLBound(1)
    integer                         :: ungriddedUBound(1)
    integer                         :: i
    real(ESMF_KIND_R8), pointer     :: farrayptr(:,:)
    character(len=*), parameter     :: rname="SetScalarFieldValues"

    rc = ESMF_SUCCESS

    if (size(idxs).ne.size(vals)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=rname//": ERROR must provide scalar_field_idx for each value", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_FieldGet(field, ungriddedLBound=ungriddedLBound, &
      ungriddedUBound=ungriddedUBound, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    if (any(idxs(:).lt.ungriddedLBound(1)).or. &
        any(idxs(:).gt.ungriddedUBound(1))) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=rname//": ERROR scalar_field_idx outside scalar_field_count", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_FieldGet(field, farrayPtr=farrayptr, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    do i=lbound(idxs,1), ubound(idxs,1)
      farrayptr(idxs(i),1) = vals(i)
    enddo

  end subroutine SetScalarFieldValues

  !-----------------------------------------------------------------------------

  ! create a field for scalar data
  function CreateScalarField(name, field_count, rc)
    type(ESMF_Field)         :: CreateScalarField
    character(*), intent(in) :: name
    integer, intent(in)      :: field_count
    integer, intent(inout)   :: rc

    ! local variables
    character(len=*), parameter     :: rname="CreateScalarField"
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid

    rc = ESMF_SUCCESS

    ! create a DistGrid with a single index space element
    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    ! num of scalar values
    CreateScalarField = ESMF_FieldCreate(name=trim(name), grid=grid, &
      typekind=ESMF_TYPEKIND_R8, gridToFieldMap=(/2/), &
      ungriddedLBound=(/1/), ungriddedUBound=(/field_count/), rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

  end function CreateScalarField

  !-----------------------------------------------------------------------------

#ifdef ESPC_COUPLE
  subroutine do_export(field,fieldName,rc)
!   arguments
    type(ESMF_Field), intent(inout) :: field
    character(*), intent(in)        :: fieldName
    integer, intent(out)            :: rc
!   local variables
    character(32)               :: cname
    character(*), parameter     :: rname="do_export"
    integer                     :: i, j
    real*8, allocatable         :: expData(:,:)
    real(ESMF_KIND_RX), pointer :: field_data(:,:)
    integer                     :: tlb(2), tub(2)
    integer                     :: status

    rc = ESMF_FAILURE

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=field_data, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGetBounds(field, localDe=localDe, &
      totalLBound=tlb, totalUBound=tub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) return

!    write(*,990) lPet,tlb(1),tub(1),1+i0,ii+i0,tlb(2),tub(2),1+j0,jj+j0
!990 format('lPet,tlb...=',I4, 4I6,6x,4I6)
!    write(*,991) lPet, i0,ii,j0,jj
!991 format('lPet, i0...=',I4,4I6)

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

    field_data(:,:)=0.

    allocate(expData(tlb(1):tub(1),tlb(2):tub(2) ))
    call export_from_hycom_deb(tlb,tub,expData,fieldName,show_minmax,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) return

    do j = tlb(2),tub(2)
    do i = tlb(1),tub(1)
      field_data(i,j)=expData(i,j)
    enddo
    enddo

    if (allocated(expData)) deallocate(expData)

     rc = ESMF_SUCCESS

  end subroutine do_export

  !-----------------------------------------------------------------------------

  subroutine do_import(field,fieldName,diagnostic,rc)
!   arguments
    type(ESMF_Field), intent(in) :: field
    character(*), intent(in)     :: fieldName
    integer,intent(in)           :: diagnostic
    integer,intent(out)          :: rc
!   local variables
    character(32)               :: cname
    character(*), parameter     :: rname="do_import"
    integer                     :: i, j
    real(ESMF_KIND_RX), pointer :: field_data(:,:)
    integer                     :: tlb(2), tub(2)
    integer                     :: status
    type(ESMF_Grid)             :: grid
    type(ESMF_Field)            :: dbgField

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=field_data, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGetBounds(field, localDe=localDe,  &
      totalLBound=tlb, totalUBound=tub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) return

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

    call import_to_hycom_deb(tlb,tub,field_data,fillValue,fieldName,show_minmax,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) return

    ! Reset Import Field
    field_data = fillValue

    rc = ESMF_SUCCESS

  end subroutine do_import
#endif

!===============================================================================
end module HYCOM_Mod
!===============================================================================


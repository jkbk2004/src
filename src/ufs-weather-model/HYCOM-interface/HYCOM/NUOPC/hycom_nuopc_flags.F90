!===============================================================================
! MODULE: HYCOM NUOPC Flags Module
!
! DESCRIPTION:
!   This module provides flags for the HYCOM NUOPC cap.
!
! FLAGS:
!
!===============================================================================
#include "HYCOM_NUOPC_Macros.h"
!===============================================================================
module hycom_nuopc_flags
#define MODNAME "hycom_nuopc_flags"

  use ESMF, only: ESMF_UtilStringUpperCase, ESMF_SUCCESS

!===============================================================================
! settings
!===============================================================================
  implicit none
  private
  save

!===============================================================================
! flags
!===============================================================================
  type import_flag
    sequence
    private
      integer :: imp
  end type import_flag

  type(import_flag), parameter :: &
    IMPORT_ERROR     = import_flag(-1), & ! Import setting is invalid
    IMPORT_REQUIRED  = import_flag(0),  & ! Fails if a import is not connected
    IMPORT_UNCOUPLED = import_flag(1),  & ! Removes all import fields
    IMPORT_FLEXIBLE  = import_flag(2)     ! Remove import if not connected

!===============================================================================
! public
!===============================================================================
  public import_flag
  public IMPORT_ERROR
  public IMPORT_REQUIRED
  public IMPORT_UNCOUPLED
  public IMPORT_FLEXIBLE

  public operator(==), assignment(=)

  interface operator (==)
    module procedure import_flag_eq
  end interface

  interface assignment (=)
    module procedure import_flag_toString
    module procedure import_flag_frString
  end interface

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  function import_flag_eq(val1, val2)
    logical import_flag_eq
    type(import_flag), intent(in) :: val1, val2
    import_flag_eq = (val1%imp == val2%imp)
  end function import_flag_eq

  !-----------------------------------------------------------------------------

  subroutine import_flag_toString(string, val)
    character(len=*), intent(out) :: string
    type(import_flag), intent(in) :: val
    if (val == IMPORT_REQUIRED) then
      write(string,'(a)') 'REQUIRED'
    elseif (val == IMPORT_UNCOUPLED) then
      write(string,'(a)') 'UNCOUPLED'
    elseif (val == IMPORT_FLEXIBLE) then
      write(string,'(a)') 'FLEXIBLE'
    else
      write(string,'(a)') 'ERROR'
    endif
  end subroutine import_flag_toString

  !-----------------------------------------------------------------------------

  subroutine import_flag_frString(val, string)
    type(import_flag), intent(out) :: val
    character(len=*), intent(in) :: string
    character(len=16) :: ustring
    integer :: rc
    ustring = ESMF_UtilStringUpperCase(string, rc=rc)
    if (rc .ne. ESMF_SUCCESS) then
      val = IMPORT_ERROR
    elseif (ustring .eq. 'REQUIRED') then
      val = IMPORT_REQUIRED
    elseif (ustring .eq. 'UNCOUPLED') then
      val = IMPORT_UNCOUPLED
    elseif (ustring .eq. 'FLEXIBLE') then
      val = IMPORT_FLEXIBLE
    else
      val = IMPORT_ERROR
    endif
  end subroutine import_flag_frString

  !-----------------------------------------------------------------------------

!===============================================================================
end module hycom_nuopc_flags
!===============================================================================

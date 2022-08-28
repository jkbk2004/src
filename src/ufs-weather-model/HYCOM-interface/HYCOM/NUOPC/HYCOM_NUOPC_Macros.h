!-------------------------------------------------------------------------------
! HYCOM NUOPC CPP Macros
!-------------------------------------------------------------------------------
! ESMF macros for logging
#ifndef FILENAME
#define FILENAME __FILE__
#endif
#ifndef CONTEXT
#define CONTEXT  line=__LINE__,file=__FILE__
#endif
#define PASSTHRU msg=ESMF_LOGERR_PASSTHRU,CONTEXT
#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)
#define ESMF_MSGERRORCHECK(rc,txt) ESMF_LogFoundError(rcToCheck=rc,msg=txt,line=__LINE__,file=__FILE__)

! Define ESMF real kind to match HYCOM single/double precision
#if defined(ESPC_IMPEXP_SINGLE)
#define ESMF_KIND_RX ESMF_KIND_R4
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R4
#else
#define ESMF_KIND_RX ESMF_KIND_R8
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R8
#endif

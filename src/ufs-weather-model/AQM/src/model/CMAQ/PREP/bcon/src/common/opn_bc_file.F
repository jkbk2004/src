
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/BCON/src/common/opn_bc_file.F,v 1.2 2011/10/21 16:52:33 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE OPN_BC_FILE( LOGUNIT, SDATE, STIME, TSTEP, NSPCS_OUT, 
     &                        SPNAME_OUT, BC_FNAME, RINDX  )  

C***********************************************************************
 
C  Function: Sets most of the Models-3 file description parameters for 
C            the BC output file and opens it.
              
C  Preconditions: None
  
C  Key Subroutines/Functions Called: None   
 
C  Revision History:
C    Prototype created by Jerry Gipson, January, 1998          
C    Output format modified by JG 4/24/98   
C    Added units for aerosol surface area JG 4/22/99   
C    02/25/02 Steve Howard (Jeff Young) - dynamic allocation
C    01/05/05 J.Young: vert dyn alloc - Use VGRD_DEFN
C    13 Jul 11 J.Young: Replaced I/O API include files with M3UTILIO and
C                       Namelist for species definitions
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
 
C***********************************************************************
      USE HGRD_DEFN   ! Module to store and load the horizontal grid variables
      USE VGRD_DEFN   ! vertical layer specifications
      USE M3UTILIO    ! IOAPI module
      USE CGRID_SPCS
      USE BC_PARMS    ! BCON parameters

      IMPLICIT NONE     

C Include Files: None

C Arguments: 
      CHARACTER( 16 ) :: BC_FNAME( * )    ! Logical names of BC Output file(s) 
      CHARACTER( 16 ) :: SPNAME_OUT( * )  ! Ouput file species names

      INTEGER LOGUNIT      ! Unit number for output log
      INTEGER NSPCS_OUT    ! Number of IC species on output file(s)
      INTEGER RINDX        ! Index of root file name
      INTEGER SDATE        ! Date for IC Output
      INTEGER STIME        ! Time for IC output
      INTEGER TSTEP        ! Time step

C Parameters: None

C External Functions: Nonde

C Saved Local Variables: None

C Local Variables:
      CHARACTER(  16 ) :: BUF16          ! 16 character buffer
      CHARACTER(   3 ) :: CHRDAT         ! 3 character buffer
      CHARACTER(  16 ) :: ENV_DFLT       ! Default value for env variable
      CHARACTER(  80 ) :: ENV_DESC       ! Description of env variable 
      CHARACTER(  80 ) :: MSG            ! Log message
      CHARACTER(  16 ) :: PNAME = 'OPN_BC_FILE'  ! Procedure Name
      CHARACTER( 256 ) :: RET_VAL        ! Return value of an env variable

      INTEGER  L, N, V       ! Loop indices
      INTEGER  NUMFLS        ! Number of output IC files
      INTEGER  SPOS1, EPOS1  ! Start and end position of char string
      INTEGER  SPOS2, EPOS2  ! Start and end position of char string
      INTEGER  STATUS        ! Status code
      INTEGER  VARINDX       ! Variable index

C***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Write log information
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE( LOGUNIT, 92000 )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get the value of environment variable EXECUTION_ID
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ENV_DESC = 'Value of environment variable EXECUTION_ID' 
      ENV_DFLT = '-----------'        
      CALL ENVSTR( 'EXECUTION_ID', ENV_DESC, ENV_DFLT, RET_VAL,
     &              STATUS)
      EXECN3D = RET_VAL

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set start date, start time, and time step
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SDATE3D = SDATE
      STIME3D = STIME
      TSTEP3D = TSTEP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set start date, start time, and time step
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FTYPE3D = BNDARY3
      NVARS3D = NSPCS_OUT
      NCOLS3D = NCOLS
      NROWS3D = NROWS
      NLAYS3D = NLAYS
      NTHIK3D = NTHIK

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set coordinate data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     GDNAM3D = GDNAME_GD  ! e.g.  / 'EAST_US' /
      GDNAM3D = GRID_NAME  ! from HGRD_DEFN
      GDTYP3D = GDTYP_GD   !       / LAMGRD3   /
      P_ALP3D = P_ALP_GD   !       /  30.0 /       ! degrees
      P_BET3D = P_BET_GD   !       /  60.0 /       ! degrees
      P_GAM3D = P_GAM_GD   !       / -90.0 /       ! degrees

      XCENT3D = XCENT_GD   !       / -90.0 /       ! degrees longitude
      YCENT3D = YCENT_GD   !       /  40.0 /       ! degrees latitude
      XORIG3D = XORIG_GD   !       /  -880000.0 /  ! m, for cross grid here
      YORIG3D = YORIG_GD   !       / -1720000.0 /  ! m, for cross grid here
      XCELL3D = XCELL_GD   !       /    80000.0 /  ! m
      YCELL3D = YCELL_GD   !       /    80000.0 /  ! m

      VGTYP3D = VGTYP_GD   !       / VGSGPN3 /     ! 2: Sigma-P0
      VGTOP3D = VGTOP_GD   !       / 10000.0 /     ! Pa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set vertical levels
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO L = 1, NLAYS + 1
         VGLVS3D ( L ) = VGLVS_GD ( L )
      END DO 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set character descriptors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO N = 1, MXDESC3
         FDESC3D( N ) = ' '
         UPDSC3D( N ) = ' '
      END DO
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set-up output file header data, compute the number of output
c  files needed, and open each one
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( MOD( NSPCS_OUT, MXVARS3 ) .EQ. 0 ) THEN
         NUMFLS = NSPCS_OUT / MXVARS3
      ELSE
         NUMFLS = ( NSPCS_OUT / MXVARS3 ) + 1
      END IF
 
      WRITE( LOGUNIT, 92020 ) NUMFLS

      DO N = 1, NUMFLS

         WRITE( CHRDAT, '( I3 )' ) N
         SPOS1 = LBLANK(  CHRDAT ) + 1
         EPOS1 = LEN_TRIM( CHRDAT )

         SPOS2 = LBLANK(  BCFL_ROOT( RINDX ) ) + 1
         EPOS2 = LEN_TRIM( BCFL_ROOT( RINDX ) )

         BUF16 = BCFL_ROOT( RINDX )

         BC_FNAME( N ) = BUF16( SPOS2:EPOS2 ) // CHRDAT( SPOS1:EPOS1 )

         ENV_DESC = 'BC Output file ' // BC_FNAME( N ) 
         ENV_DFLT = ' '       
         CALL ENVSTR( BC_FNAME( N ), ENV_DESC, ENV_DFLT, RET_VAL, STATUS)

         IF ( STATUS .NE. 0 ) THEN
            MSG = 'BC output file ' // BC_FNAME( N ) //
     &            ' not assigned' 
            CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
         END IF

         FDESC3D( 1 ) = 'BCON output file ' // BC_FNAME( N )
         NVARS3D = MIN( MXVARS3 , ( NSPCS_OUT - ( N - 1 ) * 
     &                  MXVARS3 ) )
          
         DO V = 1, NVARS3D

            VARINDX = ( N - 1 ) * MXVARS3 + V
            VTYPE3D( V ) = M3REAL

            VNAME3D( V ) = SPNAME_OUT( VARINDX )
            UNITS3D( V ) = 'ppmV'
            VDESC3D( V ) = 'Variable ' // SPNAME_OUT( VARINDX )
            
            IF ( INDEX1( SPNAME_OUT( VARINDX ), N_AE_SPC, AE_SPC ) .GT. 0 )
     &            THEN

               IF ( SPNAME_OUT( V )( 1:3 ) .EQ. 'NUM' ) THEN
                  UNITS3D( V ) = '#/m**3'
               ELSE IF ( SPNAME_OUT( V )( 1:3 ) .EQ. 'SRF' ) THEN
                  UNITS3D( V ) = 'm**2/m**3'
               ELSE
                  UNITS3D( V ) = 'micrograms/m**3'
               END IF

            END IF

         END DO
         
         WRITE( LOGUNIT, 92040 ) N, BC_FNAME( N )


         WRITE( LOGUNIT, 92060 ) EXECN3D, FTYPE3D, SDATE3D, STIME3D,
     &                           TSTEP3D, NCOLS3D, NROWS3D, NLAYS3D,
     &                           NTHIK3D, NVARS3D, GDTYP3D, P_ALP3D,
     &                           P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                           XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                           VGTYP3D, VGTOP3D, GDNAM3D                           

         WRITE( LOGUNIT, 92080 ) ( VGLVS3D( L ), L = 1, NLAYS + 1 )

         WRITE( LOGUNIT, 92100 )

         DO V = 1, NVARS3D
           WRITE( LOGUNIT, 92120 ) VNAME3D( V ), VTYPE3D( V ), UNITS3D( V )
         END DO

         IF ( .NOT. OPEN3( BC_FNAME( N ), FSUNKN3, PNAME ) ) THEN
            MSG = 'Could not open nor create ' //  BC_FNAME( N ) //
     &            ' file ' 
            CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
         END IF

      END DO

      RETURN

C************************* FORMAT STATEMENTS ***************************

92000 FORMAT( // 1X, 79( '#' ) 
     &         / 1X, '#  Output File Section '
     &         / 1X, 79( '#' ) )

92020 FORMAT( // 5X, 'Total number of ouput files created: ', I1 )

92040 FORMAT( // 5X, 'IO/API Header data for BC file ', I1, ': ', A )
92060 FORMAT(  /10X, 'EXEC_ID: ', A 
     &         /10X, 'FTYPE = ', I4 
     &         /10X, 'SDATE = ', I7
     &         /10X, 'STIME = ', I7
     &         /10X, 'TSTEP = ', I7
     &         /10X, 'NCOLS = ', I4
     &         /10X, 'NROWS = ', I4
     &         /10X, 'NLAYS = ', I4
     &         /10X, 'NTHIK = ', I4
     &         /10X, 'NVARS = ', I4
     &         /10X, 'GDTYP = ', I4
     &         /10X, 'P_ALP = ', 1PE14.5
     &         /10X, 'P_BET = ', 1PE14.5
     &         /10X, 'P_GAM = ', 1PE14.5
     &         /10X, 'XCENT = ', 1PE14.5
     &         /10X, 'YCENT = ', 1PE14.5
     &         /10X, 'XORIG = ', 1PE14.5
     &         /10X, 'YORIG = ', 1PE14.5
     &         /10X, 'XCELL = ', 1PE14.5
     &         /10X, 'YCELL = ', 1PE14.5
     &         /10X, 'VGTYP = ', I4
     &         /10X, 'VGTOP = ', 1PE14.5 
     &         /10X, 'GDNAM = ', A ) 
92080 FORMAT(  /10X, 'VGLVS = ',  5 ( 1PE14.5 ) /
     &        ( 11X,  '       ', 5 ( 1PE14.5 ) ) )
92100 FORMAT( //10X, 'Output Variables:'
     &         /10X, 'Variable name     Variable type    Variable units ' )
92120 FORMAT(   10X, A16, 8X, I1, 10X, A )    

      END

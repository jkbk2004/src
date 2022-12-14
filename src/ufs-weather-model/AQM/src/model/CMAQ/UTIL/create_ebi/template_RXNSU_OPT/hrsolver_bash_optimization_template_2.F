
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
C $Header$

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE HRSOLVER( JDATE, JTIME, C, R, L )


C**********************************************************************
C
C  FUNCTION: EBI solver
C
R1  PRECONDITIONS: For SAPRC99 family of mechanisms only
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRRATES, HRG1, HRG2, HRG3
C                                     HRG4, HRPRODLOSS
R5C                                     DEGRADE
C
R2  REVISION HISTORY: Prototype created by Jerry Gipson, September, 2003
C  REVISION HISTORY: Created by EBI solver program, Dec. 2, 2014
C                     Revised hrsolver.F by Golam Sarwar, December, 2007
C                             rearranged a statement suggested by CARB
C                             to reduce the chance of code crash; it does
C                             not change results
C   21 Jun 10 J.Young: convert for Namelist redesign
C   18 Jul 14 B.Hutzell: revised 1) to use the RXNS_DATA and RXNS_FUNCTION
C                        modules instead of include files and 2) to use
C                        real(8) variables
C**********************************************************************

      USE UTILIO_DEFN
      USE RXNS_DATA
      USE HRDATA

      IMPLICIT NONE

C..INCLUDES:


C..ARGUMENTS:
      INTEGER, INTENT( IN ) :: JDATE    ! Current date (YYYYDDD)
      INTEGER, INTENT( IN ) :: JTIME    ! Current time (HHMMSS)
      INTEGER, INTENT( IN ) :: C, R, L  ! Cell col, row, lev

C..PARAMETERS:
      INTEGER, PARAMETER :: MXBKUPS = 5  ! Max no. of back-ups allowed
      INTEGER, PARAMETER :: STAT = 1     ! Status code

      REAL( 8 ), PARAMETER :: DTMIN   = 1.0D-05    ! Smallest time step allowed, min
      REAL( 8 ), PARAMETER :: EPSLON  = 1.0D-30    ! Small number
      REAL( 8 ), PARAMETER :: MAXPRED = 1.0D+03    ! Upper limit on predicted conc
      REAL( 8 ), PARAMETER :: ZERO    = 1.0D-40    ! effective zero 

C..EXTERNAL FUNCTIONS:


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ),      SAVE :: PNAME  = 'HRSOLVER' ! Program name
      LOGICAL,              SAVE :: LFIRST = .TRUE.     ! Flag for first call
      LOGICAL, ALLOCATABLE, SAVE :: LEBISPFL( : )       ! Flag for EBI species
      LOGICAL, ALLOCATABLE, SAVE :: MXFL    ( : )       ! MAXPRED flag for species

C..SCRATCH LOCAL VARIABLES:

      CHARACTER( 132 ) :: MSG           ! Message text

      INTEGER CELLNO          ! Cell no. fo debug output
      INTEGER ITER            ! Loop index for Backward Euler iterations
      INTEGER S               ! Loop index for species
      INTEGER SP              ! Imbedded loop index for species
      INTEGER NEBI            ! Loop index for time steps
      INTEGER NINR            ! No. of inner time steps
      INTEGER N               ! Loop index
R5      INTEGER M               ! species index
      INTEGER EBI             ! Loop index
      INTEGER NBKUPS          ! No. of times time step reduced
      INTEGER ERR             ! Allocate error flag


      LOGICAL LEBI_CONV          ! Flag for EBI convergence
      LOGICAL LEBI_MAXE          ! Flag if any MXFL true

      REAL( 8 ) DTC              ! Time step to take
      REAL( 8 ) FXDLOSS          ! Total loss due to negative stoichiometry
      REAL( 8 ) VARLOSS          ! Loss excluding negative stoichiometry


#ifdef hrdebug
      CHARACTER( 8 ) :: NOTE  ! Convergence fail note

      INTEGER COL             ! Column to generate deboug output for
      INTEGER ROW             ! Row to generate deboug output for
      INTEGER LEV             ! Level to generate deboug output for
      INTEGER DBGOUT          ! Output unit for debu outpt

      LOGICAL LDEBUG          ! Debug output flag
      LOGICAL, SAVE  :: LOPEN = .FALSE.
#endif

R5      INTERFACE ! for degrade algorithm
R5        SUBROUTINE DEGRADE( CBLK, DT, JDATE, JTIME, LOGDEV )
R5          USE RXNS_DATA
R5          USE DEGRADE_SETUP_TOX
R5          IMPLICIT NONE
R5          REAL( 8 ), INTENT( IN ) :: CBLK( : )      ! array holding species concentrations
R5          REAL( 8 ), INTENT( IN ) :: DT             ! time step for integrations [sec]
R5          INTEGER,   INTENT( IN ) :: JDATE          ! current model date , coded YYYYDDD
R5          INTEGER,   INTENT( IN ) :: JTIME          ! current model time , coded HHMMSS
R5          INTEGER,   INTENT( IN ) :: LOGDEV         ! unit number for the log
R5        END SUBROUTINE
R5      END INTERFACE

C**********************************************************************


       IF( LFIRST ) THEN
          LFIRST = .FALSE.
          LOGDEV = INIT3( )

         ALLOCATE ( MXFL( NUMB_MECH_SPC ), LEBISPFL( NUMB_MECH_SPC ), STAT = ERR )
         IF ( ERR .NE. 0 ) THEN
            MSG = 'Error allocating MXFL, LEBISPFL'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         END IF

       END IF


c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
      COL = 0 - PECOL_OFFSET
      ROW = 0 - PEROW_OFFSET
      LEV = 0
      IF( C .EQ. COL .AND. R .EQ. ROW .AND. L .EQ. LEV ) THEN
c      IF( JTIME .EQ. 160000 ) THEN
         LDEBUG = .TRUE.
      ELSE
         LDEBUG = .FALSE.
      END IF

      IF( LDEBUG ) THEN
           IF( .NOT. LOPEN ) THEN
              DBGOUT = JUNIT()
              OPEN( UNIT = DBGOUT, FILE = 'debug.out' )
              LOPEN = .TRUE.
           END IF

           WRITE( DBGOUT, '( A, 2I4, I3, 1X, I7, 1X, I6 ) ' )
     &             'Debug output for col/row/lev/date/time:',
     &              C, R, L, JDATE, JTIME
           WRITE( DBGOUT, '( A, F7.2) ' )
     &             'EBI_TMSTEP = ', EBI_TMSTEP
           WRITE( DBGOUT, '( A )' ) 'Starting concs and rate constants'
           DO N = 1, N_SPEC
             WRITE( DBGOUT,  '( A, I3, 1X, A, 1X, 1PD13.5 )' )
     &                     'SP ',N, CHEMISTRY_SPC( N ), YC( N )
           END DO
           DO N = 1, N_RXNS
             WRITE( DBGOUT, '( A, I3, 1X, 1PD13.5 )' )
     &                     'RKI ', N, RKI(  N )
           END DO
      END IF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++



      DTC         = EBI_TMSTEP

      N_EBI_IT    = 0
      NBKUPS      = 0
      N_INR_STEPS = 1
      NEBI        = 1

      TSTEP_EBI: DO    ! EBI time-step loop
         
         SUBSTEP_EBI: DO NINR = 1, N_INR_STEPS   ! time substeps loop 
!  first attempt is sub time-step equals EBI time-step

            FORALL( S = 1:N_SPEC )YC0( S ) = YC( S )   ! Set ICs for EBI iterations

            ITER_SUBSTEP: DO ITER = 1, NEBITER  ! iteration loop solving for sub time-step

               N_EBI_IT = N_EBI_IT + 1
R3
               CALL HRRATES

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'ITER NO ', ITER
                  WRITE( DBGOUT, '( A, F12.5 )' )
     &               ' DTC=', DTC

                  IF( ITER .EQ. 1 ) THEN
                     WRITE( DBGOUT, '( A )' ) 'Starting reaction rates'
                     DO N = 1, N_RXNS
                        WRITE( DBGOUT, '( A, I3, 1X, 1PD13.5 )' )
     &                        'RXRAT ', N, RXRAT( NCELL, N )
                     END DO
                  END IF
               END IF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Call routines to compute concentrations of groups 1-4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               CALL HRG1( DTC )

               CALL HRG2( DTC )

               CALL HRG3( DTC )

               CALL HRG4( DTC )


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do the Euler backward method
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               CALL HRPRODLOSS

               DO EBI = 1, N_EBISP
                  S = EBISP( EBI )
                  YCP( S ) =  YC( S )*( ( YC0( S ) + PROD( S ) * DTC )
     &                     / ( YC( S ) + LOSS( S ) * DTC ) )
               END DO

R4

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check for convergence
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               LEBI_CONV = .TRUE.
               LEBI_MAXE = .FALSE.
               MXFL      = .FALSE.
               
               DO S = 1, N_SPEC
                  LEBISPFL( S ) = .FALSE.
                  YCP( S ) = MAX( ZERO, YCP( S ) )
                  IF(ABS(YC(S)-YCP(S)) .GT. RTOL(S)*(YC(S)+YCP(S)))THEN
                     LEBI_CONV     = .FALSE.
                     LEBISPFL( S ) = .TRUE.
                  END IF
c..test if predictions growing too large
                  IF( YCP( S ) .GT. MAXPRED ) THEN
                     MXFL( S )  = .TRUE.
                     LEBI_MAXE  = .TRUE.      
                     CYCLE
                  END IF                  
                  YC( S ) = YCP( S )
               END DO
c..test if predictions growing too large, abort interation loop
               IF( LEBI_MAXE ) EXIT ITER_SUBSTEP

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'Concs after ITER= ', ITER
                  DO S = 1, N_SPEC

                     IF( LEBISPFL( S ) ) THEN
                        NOTE = 'CONV FAIL'
                     ELSE
                        NOTE = '         '
                     END IF

                     WRITE( DBGOUT, '( I3, 1X, A, 1PD13.5, 1X, A )' )
     &                            S, CHEMISTRY_SPC( S ), YC( S ), NOTE
                  END DO
                  IF( LEBI_CONV ) WRITE( DBGOUT, '( A )' )
     &                 '****Convergence achieved'
               END IF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++

               IF( LEBI_CONV ) THEN 
R5
R5                  DTG = 60.0D0 * DTC
R5                  DO S = 1, NUMB_MECH_SPC
R5                     M = CGRID_INDEX( S )
R5                     YCCELL( M ) = YC( S )
R5                  END DO
R5                  IF(CALL_DEG)CALL DEGRADE(YCCELL, DTG, JDATE, JTIME, LOGDEV) ! :WTH Call degradation routine
R5
                 CYCLE SUBSTEP_EBI ! solve for next sub time-step
               END IF

            END DO ITER_SUBSTEP
! interating for substep failed, attempt to increase number of sub time-steps
            NBKUPS = NBKUPS + 1

            IF ( DTC .GT. DTMIN ) THEN
! reset YC and cut sub time-step in half  
               IF ( LEBI_MAXE ) THEN
                  WRITE( LOGDEV, 92009 )(C+PECOL_OFFSET),(R+PEROW_OFFSET), L, NBKUPS
                  DO S = 1, N_SPEC 
                     IF( .NOT. MXFL( S ) )CYCLE
                     WRITE( LOGDEV, 92010 ) TRIM( CHEMISTRY_SPC( S ) ), 
     &               YCP( S ), PROD( S ),LOSS( S )
                  END DO
               ELSE
                  WRITE( LOGDEV, 92000 ) C+PECOL_OFFSET, R+PEROW_OFFSET, L, NBKUPS
               END IF

               FORALL( S = 1:N_SPEC )YC( S ) = YC0( S )   ! Set ICs for EBI time-step

               DTC = 0.5D0 * DTC

               N_INR_STEPS = 2 ** NBKUPS
               
               EXIT SUBSTEP_EBI

            ELSE

               WRITE( LOGDEV, 92040 ) C+PECOL_OFFSET, R+PECOL_OFFSET, L

               WRITE( LOGDEV, 92060 )
               DO S = 1, N_SPEC
                   WRITE( LOGDEV, 92080 ) CHEMISTRY_SPC( S ),YC0( S ),
     &             YCP( S ),PROD( S ),LOSS( S ),
     &             ABS(YC0( S )-YCP( S ))/(YC0( S )+YCP( S )) ,
     &             LEBISPFL( S ),MXFL( S )
               END DO
               DO N = 1, N_RXNS
                  WRITE( LOGDEV, 92085 ) RXLABEL( N ), 
     &            ':RKI, RXRAT = ',RKI( N ), RXRAT( N )
               END DO
               MSG = 'ERROR: Stopping because of EBI convergence failures'
               CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT2 )

            END IF

         END DO SUBSTEP_EBI

         IF( LEBI_CONV )THEN 
             NEBI        = 1 + NEBI
! test for completing final EBI time-step
             IF( NEBI .GT. N_EBI_STEPS )RETURN
! test whether backups were done
             IF( NBKUPS .EQ. 0 )CYCLE TSTEP_EBI 
! Reset NBKUPS, N_INR_STEPS and sub-time step
             NBKUPS      = 0
             N_INR_STEPS = 1
             DTC         = EBI_TMSTEP
         END IF
         
      END DO TSTEP_EBI

      RETURN


92000 FORMAT( 'WARNING: EBI Euler convergence failure' /
     &        '         Reducing EBI time step because of ',
     &         'convergence failure for ' /
     &        '         Cell (', I3, ', ', I3, ', ', I3, ')' ,
     &        '  Back-up number', I2 )

92009 FORMAT( 'WARNING: EBI Euler convergence failure' /
     &        '         Reducing EBI time step because of ',
     &         'MAXPRED convergence failure for ' /
     &        '         Cell (', I3, ', ', I3, ', ', I3, ')' ,
     &        ' and the below species at Back-up number =', I2 )

92010 FORMAT( A16,' = ', ES12.4,'ppm :PROD, LOSS = ',2(ES12.4,' ppm/min ') )

92040 FORMAT( 'ERROR: Minimum EBI time step surpassed:'
     &      / '      Convergence failure for cell (', I3, ', ', I3,
     &                ', ', I3, ')' )

92041 FORMAT( 'ERROR: Max number of EBI time step reductions exceeded'
     &      / '      Convergence failure for cell (', I3, ', ', I3,
     &                ', ', I3, ')' )
92060 FORMAT( '      Convergence failures found in the following species:'
     &      / 10X,'    SPECIES     ',' YCO,ppm      ',
     &            '   YCP,ppm   ', ' PROD,ppm/min',' LOSS,ppm/min' ,
     &            '   %ERROR    ', 'NONCONVERGENC',' EXCEED MAX  ')

92080 FORMAT( 10X, A16, 5(1X,ES12.4), 2(5X,L8) )

92085 FORMAT( A16, 1X, A, 2(1X, ES13.5 ))
      END

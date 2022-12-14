
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
C $Header: /project/yoj/arc/CCTM/src/gas/ebi_cb05tucl/hrsolver.F,v 1.4 2011/10/21 16:10:54 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HRSOLVER( JDATE, JTIME, C, R, L )

C-----------------------------------------------------------------------
C  Function: EBI solver
C
C  Preconditions: For the CB05TUCL family of mechanisms
C
C  Key Subroutines/Functions Called:  HRRATES, HRG1, HRG2, HRG3,
C                                     HRG4, HRPRODLOSS
C
C  Revision History: Created by EBI solver program, Jan. 21, 2011
C                     Revised hrsolver.F by Golam Sarwar, December, 2007
C                             rearranged a statement suggested by CARB
C                             to reduce the chance of code crash; it does
C                             not change results
C   21 Jun 10 J.Young: convert for Namelist redesign
C   29 Mar 11 S.Roselle: Replaced I/O API include files with UTILIO_DEFN
C-----------------------------------------------------------------------

      USE CGRID_SPCS, ONLY: N_GC_SPC, GC_SPC     ! CGRID mechanism species      
      USE UTILIO_DEFN
      USE HRDATA

      IMPLICIT NONE

C..Includes:
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block

C..Arguments:
      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      INTEGER C, R, L         ! Cell col, row, lev

C..Parameters:
      INTEGER, PARAMETER :: MXBKUPS = 5  ! Max no. of back-ups allowed
      INTEGER, PARAMETER :: STAT = 1     ! Status code

      REAL, PARAMETER :: EPSLON = 1.0E-30     ! Small number
      REAL, PARAMETER :: MAXPRED = 1.0E+03    ! Upper limit on predicted conc
      REAL, PARAMETER :: ZERO    = 0.0               ! zero
      REAL, PARAMETER :: DTMIN   = 1.0E-5
      REAL, PARAMETER :: ATOL    = 1.0E-5
   
C..External Functions:

C..Saved Local Variables:
      CHARACTER( 16 ), SAVE ::  PNAME = 'HRSOLVER'      ! Program name
      LOGICAL, ALLOCATABLE, SAVE :: LEBISPFL( : )  ! Flag for EBI species
      LOGICAL, ALLOCATABLE, SAVE :: MXFL( : )      ! hit MAXPRED flag
      REAL,    ALLOCATABLE, SAVE :: ERROR( : )     ! hit MAXPRED flag
      LOGICAL, SAVE :: LFIRST = .TRUE.             ! Flag for first call

C..Scratch Local Variables:

      CHARACTER( 132 ) :: MSG           ! Message text

      INTEGER CELLNO        ! Cell no. fo debug output
      INTEGER ITER          ! Loop index for Backward Euler iterations
      INTEGER S             ! Loop index for species
      INTEGER NEBI          ! Loop index for time steps
      INTEGER NINR          ! No. of inner time steps
      INTEGER N             ! Loop index
      INTEGER EBI           ! Loop index
      INTEGER NBKUPS        ! No. of times time step reduced
      INTEGER ERR           ! Allocate error flag

      LOGICAL LEBI_CONV     ! Flag for EBI convergence

      REAL DTC              ! Time step to take
      REAL DYDT, YCBAR
      REAL DTNEW
      REAL FXDLOSS          ! Total loss due to negative stoichiometry
      REAL VARLOSS          ! Loss excluding negative stoichiometry


#ifdef hrdebug
      CHARACTER( 8 ) :: NOTE  ! Convergence fail note

      INTEGER COL             ! Column to generate deboug output for
      INTEGER ROW             ! Row to generate deboug output for
      INTEGER LEV             ! Level to generate deboug output for
      INTEGER DBGOUT          ! Output unit for debu outpt

      LOGICAL LDEBUG          ! Debug output flag
      LOGICAL, SAVE :: LOPEN = .FALSE.
#endif

C-----------------------------------------------------------------------

      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         LOGDEV = INIT3( )

         ALLOCATE ( MXFL( N_GC_SPC ), LEBISPFL( N_GC_SPC ), STAT = ERR )
         IF ( ERR .NE. 0 ) THEN
            MSG = 'Error allocating MXFL, LEBISPFL'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         END IF
         ALLOCATE ( ERROR( N_GC_SPC ), STAT = ERR )
         IF ( ERR .NE. 0 ) THEN
            MSG = 'Error allocating ERROR'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         END IF
      END IF

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
      COL = 0
      ROW = 0
      LEV = 0
      IF ( C .EQ. COL .AND. R .EQ. ROW .AND. L .EQ. LEV ) THEN
c      IF ( JTIME .EQ. 160000 ) THEN
         LDEBUG = .TRUE.
      ELSE
         LDEBUG = .FALSE.
      END IF

      IF ( LDEBUG ) THEN
           IF ( .NOT. LOPEN ) THEN
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
             WRITE( DBGOUT,  '( A, I3, 1X, A, 1X, 1PE13.5 )' )
     &                     'SP ',N, GC_SPC( N ), YC( N )
           END DO
           DO N = 1, N_RXNS
             WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                     'RKI ', N, RKI( N )
           END DO
      END IF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++

      N_EBI_IT = 0
      ERROR    = 0.0

      DO 3000 NEBI = 1, N_EBI_STEPS    ! outer EBI time-tep loop

         DTC = EBI_TMSTEP
         NBKUPS = 0
         N_INR_STEPS = 1

 100     CONTINUE                        !  Restart location

         DO 2000 NINR = 1, N_INR_STEPS   ! No. of time steps for back-up

            DO S = 1, N_SPEC             ! Set ICs for EBI iterations
               YC0( S ) = YC( S )
            END DO

            DO 1000 ITER = 1, NEBITER    ! EBI iteration loop

               N_EBI_IT = N_EBI_IT + 1

               CALL HRRATES

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF ( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'ITER NO ', ITER
                  WRITE( DBGOUT, '( A, F12.5 )' )
     &               ' DTC=', DTC

                  IF ( ITER .EQ. 1 ) THEN
                     WRITE( DBGOUT, '( A )' ) 'Starting reaction rates'
                     DO N = 1, N_RXNS
                        WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                        'RXRAT ', N, RXRAT( N )
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
                  YCP( S ) = YC( S ) * ( ( YC0( S ) + PROD( S ) * DTC ) /
     &                                   ( YC ( S ) + LOSS( S ) * DTC ) )
               END DO

c..Special treatment of PAR because of negative product stoichiometry
               IF ( PNEG( PAR ) .GT. 0.0 ) THEN
                  FXDLOSS = PNEG( PAR ) * DTC
                  IF ( FXDLOSS .GE. YC0( PAR ) + PROD( PAR ) * DTC ) THEN
                     YCP( PAR ) = 0.0
                  ELSE
                     VARLOSS = MAX( LOSS( PAR ) - PNEG( PAR ) , ZERO )
                     YCP( PAR ) = ( YC0( PAR ) + PROD( PAR ) * DTC  - 
     &                   FXDLOSS ) / ( 1.0 + VARLOSS * DTC / YC( PAR ) )
                  END IF
               END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check for convergence
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               LEBI_CONV = .TRUE.
               MXFL = .FALSE.
               DO S = 1, N_SPEC
                  LEBISPFL( S ) = .FALSE.
                  YCP( S ) = MAX( EPSLON, YCP( S ) )
                  ERROR( S ) = ABS( YC( S ) - YCP( S ) ) 
     &                       / (  YC( S ) + YCP( S ) ) 
                  IF ( ABS( YC( S ) - YCP( S ) ) .GT. RTOL( S ) *
     &               ( YC( S ) + YCP( S ) ) ) THEN
!                IF ( ERROR( S ) .GT. RTOL( S ) .AND.
!    &                      ABS( YCP( S ) - YC( S ) ) .GT. ATOL ) THEN
!                IF (  ABS( YC( S ) - YCP( S ) ) .GT. RTOL( S ) * ( YC( S ) + YCP( S ) )  
!   &                      .AND.  ABS( YCP( S ) - YC( S ) ) .GT. ATOL ) THEN
                     LEBI_CONV = .FALSE.
                     LEBISPFL( S ) = .TRUE.
                  END IF
c..if predictions growing too large, treat as a convergence failure
                  IF ( YCP( S ) .GT. MAXPRED ) then
                     MXFL( S )  = .TRUE.
                     CYCLE
                  END IF
                  YC( S ) = YCP( S )
               END DO
               IF( ANY( MXFL ) ) GO TO 1010 

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF ( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'Concs after ITER= ', ITER
                  DO S = 1, N_SPEC

                     IF ( LEBISPFL( S ) ) THEN
                        NOTE = 'CONV FAIL'
                     ELSE
                        NOTE = '         '
                     END IF

                     WRITE( DBGOUT, '( I3, 1X, A, 1PE13.5, 1X, A )' )
     &                            S, GC_SPC( S ), YC( S ), NOTE
                  END DO
                  IF ( LEBI_CONV ) WRITE( DBGOUT, '( A )' )
     &                 '****Convergence achieved'
               END IF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


               IF ( LEBI_CONV ) GO TO 2000

1000        CONTINUE

c...Convergence failure section; cut the inner time step in half &
c.....start inner loop over unless max backups exceeded

1010        CONTINUE

            NBKUPS = NBKUPS + 1

!            IF ( NBKUPS .LE. MXBKUPS ) THEN
            IF ( DTC .GT. DTMIN ) THEN

               IF ( ANY( MXFL ) ) THEN
                  WRITE( LOGDEV, 92009 )(C+PECOL_OFFSET),(R+PEROW_OFFSET), L, NBKUPS
                  DO S = 1, N_SPEC 
                     IF( .NOT. MXFL( S ) )CYCLE
                     WRITE( LOGDEV, 92010 ) TRIM( GC_SPC( S ) ), 
     &               YCP( S ), PROD( S ),LOSS( S )
                  END DO
               ELSE
                  WRITE( LOGDEV, 92000 ) C+PECOL_OFFSET, R+PEROW_OFFSET, L, NBKUPS
               END IF

               DO S = 1, N_SPEC
                  YC( S ) = YC0( s )
               END DO
                
!               DTNEW = DTC
!               DO EBI = 1, N_EBISP
!                 S     = EBISP( EBI )
!                  IF( MXFL( S ) )THEN
!                     DYDT = ABS( PROD( S ) - LOSS( S )*YC0( S ) )
!                  ELSE
!                     YCBAR = 0.5*(YC0( S )+YCP( S ))
!                     DYDT = ABS( PROD( S ) - LOSS( S )*YCBAR ) 
!                  END IF
!                  IF( DYDT*DTC .LE. 1.0E-30 )CYCLE
!                  DTNEW = MIN( DTNEW,
!    &                     ( 1.95 * RTOL( S ) * YC0( S ) )
!    &                   / ( (1.0 - RTOL( S ))*DYDT ) )
!              END DO
!               IF ( MOD( DTC, DTNEW ) .EQ. 0 ) THEN
!                    N_INR_STEPS = MAX( INT( DTC / DTNEW ), 1 )
!              ELSE
!                    N_INR_STEPS = INT( DTC / DTNEW ) + 1
!                    DTNEW =  DTC / FLOAT( N_INR_STEPS )
!               END IF
!               DTC = DTNEW

               DTC = 0.5 * DTC
               N_INR_STEPS = 2 ** NBKUPS

               GO TO 100

            ELSE

               WRITE( LOGDEV, 92040 ) C+PECOL_OFFSET, R+PECOL_OFFSET, L

               WRITE( LOGDEV, 92060 )
               DO S = 1, N_SPEC
                   WRITE( LOGDEV, 92080 ) GC_SPC( S ),YC0( S ),
     &             YCP( S ),PROD( S ),LOSS( S ),ERROR( S ),
     &             LEBISPFL( S ),MXFL( S )
               END DO
               DO N = 1, N_RXNS
                  WRITE( LOGDEV, '( A16, 1X, A, 2(1X, 1PE13.5 ))' )
     &                   RXLABEL( N ) , ':RKI, RXRAT = ',RKI( N ), RXRAT( N )
               END DO

               MSG = 'ERROR: Stopping because of EBI convergence failures'
               CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )

            END IF

2000     CONTINUE

3000  CONTINUE

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

92010 FORMAT( A16,' = ', E12.4,' :PROD, LOSS = ',2(ES12.4,1X) )

92040 FORMAT( 'ERROR: Minimum EBI time step surpassed:'
     &      / '      Convergence failure for cell (', I3, ', ', I3,
     &                ', ', I3, ')' )

92041 FORMAT( 'ERROR: Max number of EBI time step reductions exceeded'
     &      / '      Convergence failure for cell (', I3, ', ', I3,
     &                ', ', I3, ')' )
92060 FORMAT( '      Convergence failures found in the following species:'
     &      / 10X,'    SPECIES     ','    YCO       ',
     &            '     YCP     ', '    PROD     ','    LOSS     ' ,
     &            '    ERROR    ', 'NONCONVERGENC',' EXCEED MAX  ')

92080 FORMAT( 10X, A16, 5(1X,ES12.4), 2(5X,L8) )

      END

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

C****************************************************************************************
C     Converts hourly data values to daily
C    
C      Using 11 different operations
C   [sum, avg, min, max, hr@min, hr@max, @maxT, maxdif, 8hrmax, w126,
C    @8hrmaxO3, hr@8hrmax, sum06]
C
C          sum - sums the 24 hour values
C          avg - sums the 24 values and divides by 24
C          min - uses the minimum hourly value
C          max - uses the maximum hourly value
C          hr@min - hour at the minimum hourly value
C          hr@max - hour at the maximum hourly value
C          @maxT - uses the hourly value at maximum temperature
C          maxdif - uses the maximum hourly change
C          8hrmax - uses the maximum 8 hour period 
C          w126 - uses secondary standard ozone calculation 
C          @8hrmaxO3 - averages the value within the 8-hr-max ozone period
C          hr@8hrmax - Starting hour of the 8-hr-max period
C          sum06 - sums all values>=0.06ppm between 8am & 8pm
C
C      A daily value is marked as missing when fewer then 18 hourly values
C      are valid.  For 8-hr-max calculations, each 8-hour period must have
C      six or more valid hourly values.
C 
C      USELOCAL using local time, else uses GMT
C
c      PROGRAM hr2day.exe
c
C*******************************************************************************
      USE M3UTILIO
      USE species_def
      USE evaluator

      IMPLICIT NONE

C External functions
      integer getTZ

C local variables
      integer status
      logical rstatus
      integer logdev
      integer c, r, n, h, i 
      integer tzadj, hr1
      integer curdate, curtime, cdate, ctime
      integer first_date, first_time, last_date, last_time, runlen
      real x,y,longitude,latitude
      integer hroffset
      integer startHr, endHr
      logical useLocal
      logical useDST  
      logical partDay 

C Array to store hourly input values
      real, allocatable :: hrValues (:, :, :)
      real, allocatable :: dayValues (:, :)
      integer, allocatable :: tzoffset(:,:)
      integer, allocatable :: offset(:,:)


C... start program
      logdev = init3 ()

C... open input file
      if( .not. open3('INFILE',fsread3, 'hr2day')) then
        Call m3err('hr2day', 0, 0, 'Could not open INPUT_FILE', .TRUE.)
        endif

C... load file description from INFILE
      if( .not. desc3('INFILE')) then
        Call m3err ('hr2day', 0, 0, 'Could not load file description from IN_FILE', .TRUE.)
        endif

C... check for 1 hour time step
      if( TSTEP3D.ne.10000 ) then
        Call m3err('hr2day', 0, 0, 'Input file must have One hour time step', .TRUE.)
        endif
 
c
c... set up the correct projection
c 
      if (GDTYP3D.eq.1) then !lat/lon, nothing to initialize
       continue
      elseif (GDTYP3D.eq.2) then !initialize Lambert projection
       if( .NOT. SETLAM( Real(P_ALP3D), Real(P_BET3D), Real(P_GAM3D),                            
     &                  Real(XCENT3D), Real(YCENT3D) )) then                                    
        Call m3err ('hr2day', 0, 0, 'Lambert projection setup error', .TRUE.)
       endif     
      elseif (GDTYP3D.eq.6) then !initialize polar stereographic projection
       if( .NOT. SETPOL( Real(P_ALP3D), Real(P_BET3D), Real(P_GAM3D),                            
     &                  Real(XCENT3D), Real(YCENT3D) )) then                                    
        Call m3err ('hr2day', 0, 0, 'Polar Stereographic projection setup error', .TRUE.)
       endif  
      else
       Call m3err ('hr2day', 0, 0, 'Grid projection not supported', .TRUE.)   
      endif

C... Compute last date and time to process                    
      first_date = SDATE3D
      first_time = STIME3D                                      
      last_date = SDATE3D
      last_time = STIME3D                                      
      runlen = SEC2TIME( (MXREC3D-1) * TIME2SEC( TSTEP3D ))
      Call NEXTIME( last_date, last_time, runlen ) 

C... get species definitions from system variables
      Call loadSpecies()

C... get compute partial day switch
      partDay = ENVYN('PARTIAL_DAY',"Compute partial day values", .FALSE., status)

C... get use local time switch
      useLocal = ENVYN('USELOCAL',"Use local time", .FALSE., status)

C... get use daylight savings time switch id uselocal
      useDST = .false.
      if( useLocal ) then
        useDST = ENVYN('USEDST',"Use Daylight Savings time", .FALSE., status)
        endif

C... get hour offset
      hroffset = ENVINT('HROFFSET',"Hour offset", 0, status)

C... get starting and ending hours for time window
      startHr = ENVINT('START_HOUR',"Starting Hour of time window", 0, status)
      endHr = ENVINT('END_HOUR',"Ending Hour of time window", 23, status)

C... create output file
      NLAYS3D = 1
      STIME3D = 000000 
      TSTEP3D = 240000
      NVARS3D = NSPECVAR
      do n = 1, NSPECVAR
        VNAME3D(n) = SPECVARS(n)%NAME
        VDESC3D(n) = SPECVARS(n)%DESCRIPTION
        UNITS3D(n) = SPECVARS(n)%UNITS
        VTYPE3D(n) = M3REAL
        enddo

       if(.not. open3('OUTFILE',3,'hr2day')) then
        if(.not. open3('OUTFILE',2,'hr2day')) then
          Call m3err('average', 0, 0, 'Could not open OUTFILE file',.TRUE.)
          endif
        endif

C... Allocate memory for data arrays
      Allocate( tzoffset(NCOLS3D, NROWS3D) )
      Allocate( offset(NCOLS3D, NROWS3D) )
      Allocate( hrValues (-36:66, NCOLS3D, NROWS3D) )
      Allocate( dayValues (NCOLS3D, NROWS3D) )

      ! set tzoffset array values 
      tzoffset = 0
      if( useLocal ) then
        write(*,'(/,''Computing timezone offsets for grid cells'',/)')
        do c=1,NCOLS3D
          do r=1,NROWS3D
            x = XORIG3D + (c-0.5) * XCELL3D
            y = YORIG3D + (r-0.5) * YCELL3D
            
            if (GDTYP3D.eq.1) then !lat/lon grid, x/y already lon/lat
             longitude=x
             latitude=y
            elseif (GDTYP3D.eq.2) then ! convert Lambert coordinates to lat/lon
             if( .NOT. LAM2LL(x, y, longitude, latitude) ) then
              Call m3err ('hr2day', 0, 0, 'Lat/Lon conversion error', .TRUE.)
             endif
            elseif (GDTYP3D.eq.6) then !convert polar stereographic coordinates to lat/lon
             if( .NOT. POL2LL(x, y, longitude, latitude) ) then
              Call m3err ('hr2day', 0, 0, 'Lat/Lon conversion error', .TRUE.)
             endif
            else
             Call m3err ('hr2day', 0, 0, 'Grid projection not supported', .TRUE.)   
            endif

            tzoffset(c,r) = getTZ(longitude, latitude)
            enddo
          enddo
        endif

C... reload file description from INFILE
      if( .not. desc3('INFILE')) then
        Call m3err ('hr2day', 0, 0, 'Could not load file description from IN_FILE', .TRUE.)
        endif

C... start loop to read and process each variable
      do n = 1, NSPECVAR
        hr1 = -36
        curDate = SDATE3D
        curTime = STIME3D
        hrValues = BADVAL3
        dayValues = BADVAL3

        ! start loops to read hourly values for each day 
        Do
          cDate = curDate
          cTime = hr1 * 10000
          offset = tzoffset
          offset = offset + hroffset 
          if(useDST .and. isDSTime(cDate)) offset = offset - 1   ! adjustment for daylight savings time
	  call NEXTIME(cDate,cTime,0) !reformat to make sure hours are between 0 and 23
          Do h=hr1,66
            rstatus = .false.
            hrValues(h,:,:) = BADVAL3
            if( (SECSDIFF(first_date, first_time, cdate, ctime) .ge. 0)
     &    .and. (SECSDIFF(cdate, ctime, last_date, last_time) .ge. 0) ) then
              Call evaluate(SPECVARS(n)%EXPRESSION, cdate, ctime, 1, NCOLS3D*NROWS3D, hrValues(h,:,:))
              endif
            Call NEXTIME(cDate, cTime, 10000)
            enddo

          ! compute daily values
          if( SPECVARS(n)%OPERATION .eq. 'SUM' ) then
            Call sumValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'AVG' ) then
            Call avgValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'MIN' ) then
            Call minValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'MAX' ) then
            Call maxValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'HR@MIN' ) then
            Call minHrValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'HR@MAX' ) then
            Call maxHrValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. '@MAXT' ) then
            Call maxTValues(hrValues, dayValues, offset, curDate, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'MAXDIF' ) then
            Call maxDifValues(hrValues, dayValues, offset, startHr, endHr, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. '8HRMAX' ) then
            Call max8hr(hrValues, dayValues, offset, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'W126' ) then
            Call w126(hrValues, dayValues, offset, SPECVARS(n)%UNITS, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. '@8HRMAXO3' ) then
            Call maxO3Values(hrValues, dayValues, offset, curDate, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'HR@8HRMAX' ) then
            Call max8hrHour(hrValues, dayValues, offset, partDay)
            endif

          if( SPECVARS(n)%OPERATION .eq. 'SUM06' ) then
            Call sum06(hrValues, dayValues, offset, SPECVARS(n)%UNITS, partDay)
            endif

          ! write daily values to output
          if(.not.write3('OUTFILE',SPECVARS(n)%NAME,curDate,000000,dayValues)) then
          Call m3err ('hr2day', curDate, 0, 'Write Error for ' // SPECVARS(n)%NAME, .TRUE.)
          endif

          ! copy current day's values to next day
          do h = -36, 42
            hrValues( h,:,:) = hrvalues(h+24,:,:) 
            enddo

          ! go to next time step
          Call NEXTIME(curdate, curtime, 240000) 

          ! if current date at noon is past last date and time, then exit loop
          if(SECSDIFF(curDate, 120000, last_date, last_time) .le. 0) EXIT

          hr1 = 43
          endDo  ! end (time-step-loop)

        enddo  ! var loop

      status = SHUT3 ()
      stop
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to sum hourly values at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine sumValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)  

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real sum
 
C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          sum = 0.0
          count = 0
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              sum = sum + hrValues(h,c,r)
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = sum
          endDo      ! end row loop                                                         
        endDo      ! end column loop 

      return
      end Subroutine sumValues


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to average hourly values at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine avgValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)  

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real sum
 
C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          sum = 0.0
          count = 0
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              sum = sum + hrValues(h,c,r)
              endif
            endDo      ! end hour loop
          if( count.ge.18 ) dayValues(c,r) = sum/count
          if(partDay .and. count.ge.1) dayValues(c,r) = sum/count
          endDo      ! end row loop                                                         
        endDo      ! end column loop 

      return
      end Subroutine avgValues


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the minimum hourly value at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine minValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real minValue

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          count = 0
          dayValues(c,r) = BADVAL3
          minValue = 1.0E32 
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              if( hrValues(h,c,r).lt.minValue ) minValue = hrValues(h,c,r)
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = minValue
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine minValues
  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the hour at the minimum hourly value at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine minHrValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real minValue
      Real minHour 

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          count = 0
          dayValues(c,r) = BADVAL3
          minValue = 1.0E32
          minHour = BADVAL3
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              if( hrValues(h,c,r).lt.minValue ) then
                minValue = hrValues(h,c,r)
                minHour = h - tzoffset(c,r)
                endif
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = minHour 
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine minHrValues



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the maximum hourly value at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine maxValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real maxValue

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          maxValue = BADVAL3
          dayValues(c,r) = BADVAL3
          count = 0
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              if( hrValues(h,c,r) .gt. maxValue ) maxValue = hrValues(h,c,r)
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = maxValue
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine maxValues


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the hour at the maximum hourly value at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine maxHrValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real maxValue
      Real maxHour 

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          maxValue = BADVAL3
          maxHour = BADVAL3
          dayValues(c,r) = BADVAL3
          count = 0
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              if( hrValues(h,c,r) .gt. maxValue ) then
                maxValue = hrValues(h,c,r)
                maxHour = h - tzoffset(c,r)
                endif
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = maxHour 
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine maxHrValues

 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the value when at maximum temperature
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine maxTValues(hrValues, dayValues, tzoffset, curDate, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer curDate
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer status
      Logical, save :: first=.true.  
      Integer, save :: first_date
      Integer, save :: first_time
      Integer, save :: last_date
      Integer, save :: last_time

      Character*(16) tempvar
      Real, save,allocatable :: tvalues(:,:,:)
      Integer runlen
      Integer c, r, h, count
      Integer cdate, ctime
      Real maxTemp, value


C... first pass
      if( first ) then
        first = .false.
        Call ENVSTR('TEMPERATURE','Temperature species used in @MAXT operation',
     &              'TEMP2', tempvar, status)  
        Allocate ( tvalues(-36:66, NCOLS3D, NROWS3D) )

        ! Compute last date and time to process                    
        first_date = SDATE3D
        first_time = STIME3D                                      
        last_date = SDATE3D
        last_time = STIME3D                                      
        runlen = SEC2TIME( (MXREC3D-1) * TIME2SEC( TSTEP3D ))
        Call NEXTIME( last_date, last_time, runlen ) 
        endif   ! first pass
    
C... read temperature values and fill tvalues array
      cdate = curDate
      ctime = 0
      tvalues = BADVAL3
      do h=-36,66
        if( (SECSDIFF(first_date, first_time, cdate, ctime) .ge. 0)
     &    .and. (SECSDIFF(cdate, ctime, last_date, last_time) .ge. 0) ) then
          if(.not.READ3('INFILE', tempvar, 1, cdate, ctime, tvalues(h,:,:)) ) then
            Write(*,'(''Cannot read temperature data for @MAXT operation'')')
            endif
          endif
        Call NEXTIME(cDate, cTime, 10000)
        enddo
         
C... start loops to find maximum temperature for each cell
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          maxTemp = BADVAL3
          value = BADVAL3
          count = 0
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr
            if( hrValues(h,c,r) .gt. BADVAL3 ) then
              count = count+1
              if( tvalues(h,c,r) .gt. maxTemp ) then
                maxTemp = tvalues(h,c,r)
                value = hrValues(h,c,r)
                endif
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.18) dayValues(c,r) = value
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine maxTValues


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the value in 8-hour-max ozone period
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine maxO3Values(hrValues, dayValues, tzoffset, curDate, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer curDate
      Logical partDay

C... local variables
      Logical, save :: first=.true.  
      Integer, save :: first_date
      Integer, save :: first_time
      Integer, save :: last_date
      Integer, save :: last_time

      Character*(16) ozonevar
      Real, save,allocatable :: O3values(:,:,:)
      Integer runlen
      Integer c, r, h, i
      Integer cdate, ctime
      Integer count, nperiods
      Real o3sum, sum, maxValue, sumValue, avg
      INTEGER STATUS                ! Status code
      INTEGER HOURS_8HRMAX          ! number of 8hr values to compute 8hr max     
      CHARACTER*16 ENV_DESC         ! message string
      CHARACTER*16 PNAME            ! Program Name
      CHARACTER*80 MSG              ! Error message
      DATA  PNAME           / 'HR2DAY'        /

C... Get the HOURS_8HRMAX environment variable (default is 24)                                                          
       ENV_DESC = 'Number of 8hr values to use when computing DM8HR'                                               
       HOURS_8HRMAX= ENVINT( 'HOURS_8HRMAX', ENV_DESC, 24, STATUS)  
	 
       if ( ( HOURS_8HRMAX .NE. 24) .AND. ( HOURS_8HRMAX .NE. 17) ) THEN                                

        MSG = '**Error** Invalid value for HOURS_8HRMAX, use 24 or 17'
        CALL M3ERR( PNAME, 0, 0, MSG, .TRUE. ) 
        Stop
       Endif


C... first pass
      if( first ) then
        first = .false.
        Call ENVSTR('OZONE','Ozone species used in @8HRMAXO3 operation',
     &              'O3', ozonevar, status)  
        Allocate ( O3values(-36:66, NCOLS3D, NROWS3D) )

        ! Compute last date and time to process                    
        first_date = SDATE3D
        first_time = STIME3D                                      
        last_date = SDATE3D
        last_time = STIME3D                                      
        runlen = SEC2TIME( (MXREC3D-1) * TIME2SEC( TSTEP3D ))
        Call NEXTIME( last_date, last_time, runlen ) 
        endif   ! first pass
    
C... read ozone values and fill tvalues array
      cdate = curDate
      ctime = -36 * 10000
      call NEXTIME(cDate,cTime,0) !reformat to make sure hours are between 0 and 23
      O3values = BADVAL3
      do h=-36,66
        if( (SECSDIFF(first_date, first_time, cdate, ctime) .ge. 0)
     &    .and. (SECSDIFF(cdate, ctime, last_date, last_time) .ge. 0) ) then
          if(.not.READ3('INFILE', ozonevar, 1, cdate, ctime, O3values(h,:,:)) ) then
            Write(*,'(''Cannot read Ozone data for @8HRMAXO3 operation'')')
            endif
          endif
        Call NEXTIME(cDate, cTime, 10000)
        enddo
        
C... start loops to find max 8-hr and save sum of hrValues
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          maxValue = BADVAL3
          sumValue = BADVAL3
          nperiods = 0
	  
          if ( HOURS_8HRMAX .eq. 24 ) then ! use 24 8hr values
	  
           Do h = tzoffset(c,r),tzoffset(c,r)+23
            o3sum = 0
            sum = 0
            count = 0
            Do i = 0,7
              if(( O3values(h+i,c,r) .gt. BADVAL3 ) .and.
     *           ( hrValues(h+i,c,r) .gt. BADVAL3 )) then
                o3sum = o3sum + O3values(h+i,c,r)
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
               nperiods = nperiods + 1
               avg = o3sum / count
               if( avg.gt.maxValue) then
                 maxvalue = avg   
                 sumvalue = sum/count
                 endif
               endif
            endDo      ! end hour loop
           if(partDay .or. nperiods.ge.18) dayValues(c,r) = sumValue !require 18/24
	   
          else !use only 17 8hr values, from 7 am to 11 pm LT

           Do h = tzoffset(c,r)+7,tzoffset(c,r)+23
            o3sum = 0
            sum = 0
            count = 0
            Do i = 0,7
              if(( O3values(h+i,c,r) .gt. BADVAL3 ) .and.
     *           ( hrValues(h+i,c,r) .gt. BADVAL3 )) then
                o3sum = o3sum + O3values(h+i,c,r)
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
               nperiods = nperiods + 1
               avg = o3sum / count
               if( avg.gt.maxValue) then
                 maxvalue = avg   
                 sumvalue = sum/count
                 endif
               endif
            endDo      ! end hour loop
           if(partDay .or. nperiods.ge.13) dayValues(c,r) = sumValue !require 13/17

          endif	   
	   
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine maxO3Values


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the maximum hourly change
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine maxDifValues(hrValues, dayValues, tzoffset, startHr, endHr, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Integer startHr, endHr
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real diff    
      Real maxValue

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          count = 0
          dayValues(c,r) = BADVAL3
          maxValue = BADVAL3
          Do h = tzoffset(c,r)+startHr ,tzoffset(c,r)+endHr-1
            if( hrValues(h,c,r).gt.BADVAL3 .and. hrValues(h+1,c,r).gt.BADVAL3) then
              count = 0
              diff = ABS( hrValues(h+1,c,r) - hrValues(h,c,r) )
              if( diff .gt. maxValue ) then
                maxValue = diff
                endif
              endif
            endDo      ! end hour loop
          if(partDay .or. count.ge.17) dayValues(c,r) = maxValue
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine maxDifValues


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the 8-hour maximum value at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine max8hr(hrValues, dayValues, tzoffset, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Logical partDay

C... local variables
      Integer c, r, h, i, count, nperiods
      Real sum, maxValue, avg
      INTEGER STATUS                ! Status code
      INTEGER HOURS_8HRMAX          ! number of 8hr values to compute 8hr max     
      CHARACTER*16 ENV_DESC         ! message string
      CHARACTER*16 PNAME            ! Program Name
      CHARACTER*80 MSG              ! Error message
      DATA  PNAME           / 'HR2DAY'        /

C... Get the HOURS_8HRMAX environment variable (default is 24)                                                          
       ENV_DESC = 'Number of 8hr values to use when computing DM8HR'                                               
       HOURS_8HRMAX= ENVINT( 'HOURS_8HRMAX', ENV_DESC, 24, STATUS)  
	 
       if ( ( HOURS_8HRMAX .NE. 24) .AND. ( HOURS_8HRMAX .NE. 17) ) THEN                                

        MSG = '**Error** Invalid value for HOURS_8HRMAX, use 24 or 17'
        CALL M3ERR( PNAME, 0, 0, MSG, .TRUE. ) 
        Stop
       Endif

C... start loops to find max 8-hr 
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          maxValue = BADVAL3
          nperiods = 0
	  
          if ( HOURS_8HRMAX .eq. 24 ) then ! use 24 8hr values
	  
           Do h = tzoffset(c,r),tzoffset(c,r)+23
            sum = 0
            count = 0
            Do i = 0,7
              if( hrValues(h+i,c,r) .gt. BADVAL3 ) then
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
              nperiods = nperiods + 1
              avg = sum/count
              if( avg .gt. maxValue ) maxValue = avg
              endif
            endDo      ! end hour loop
            if(partDay .or. nperiods.ge.18) dayValues(c,r) = maxValue !require 18/24
	   
          else !use only 17 8hr values, from 7 am to 11 pm LT

           Do h = tzoffset(c,r)+7,tzoffset(c,r)+23
            sum = 0
            count = 0
            Do i = 0,7
              if( hrValues(h+i,c,r) .gt. BADVAL3 ) then
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
              nperiods = nperiods + 1
              avg = sum/count
              if( avg .gt. maxValue ) maxValue = avg
              endif
            endDo      ! end hour loop
            if(partDay .or. nperiods.ge.13) dayValues(c,r) = maxValue !require 13/17
	    
	  endif
	   
          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine max8hr
 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to find the starting hour of the 8-hour maximum period
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine max8hrHour(hrValues, dayValues, tzoffset, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Logical partDay

C... local variables
      Integer c, r, h, i, count, nperiods
      Real sum, maxValue, avg
      Integer maxHour
      INTEGER STATUS                ! Status code
      INTEGER HOURS_8HRMAX          ! number of 8hr values to compute 8hr max     
      CHARACTER*16 ENV_DESC         ! message string
      CHARACTER*16 PNAME            ! Program Name
      CHARACTER*80 MSG              ! Error message
      DATA  PNAME           / 'HR2DAY'        /

C... Get the HOURS_8HRMAX environment variable (default is 24)                                                          
       ENV_DESC = 'Number of 8hr values to use when computing DM8HR'                                               
       HOURS_8HRMAX= ENVINT( 'HOURS_8HRMAX', ENV_DESC, 24, STATUS)  
	 
       if ( ( HOURS_8HRMAX .NE. 24) .AND. ( HOURS_8HRMAX .NE. 17) ) THEN                                

        MSG = '**Error** Invalid value for HOURS_8HRMAX, use 24 or 17'
        CALL M3ERR( PNAME, 0, 0, MSG, .TRUE. ) 
        Stop
       Endif

C... start loops to find max 8-hr 
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          dayValues(c,r) = BADVAL3
          maxValue = BADVAL3
          maxHour = 0
          nperiods = 0

          if ( HOURS_8HRMAX .eq. 24 ) then ! use 24 8hr values

           Do h = tzoffset(c,r),tzoffset(c,r)+23
            sum = 0
            count = 0
            Do i = 0,7
              if( hrValues(h+i,c,r) .gt. BADVAL3 ) then
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
              nperiods = nperiods + 1
              avg = sum/count
              if( avg .gt. maxValue ) then
                maxValue = avg
                maxHour = h - tzoffset(c,r)   !! maxHour is in local time (0-23)
                endif
              endif
            endDo      ! end hour loop

           if(partDay .or. nperiods.ge.18) dayValues(c,r) = 1.0 * maxHour !require 18/24

          else !use only 17 8hr values, from 7 am to 11 pm LT
	  
           Do h = tzoffset(c,r)+7,tzoffset(c,r)+23
            sum = 0
            count = 0
            Do i = 0,7
              if( hrValues(h+i,c,r) .gt. BADVAL3 ) then
                sum = sum + hrValues(h+i,c,r)
                count = count + 1
                endif
              enddo   ! end 8 hour loop
            if(count.ge.6) then
              nperiods = nperiods + 1
              avg = sum/count
              if( avg .gt. maxValue ) then
                maxValue = avg
                maxHour = h - tzoffset(c,r)   !! maxHour is in local time (0-23)
                endif
              endif
            endDo      ! end hour loop

           if(partDay .or. nperiods.ge.13) dayValues(c,r) = 1.0 * maxHour !require 13/17
	  
	  endif

          endDo      ! end row loop
        endDo      ! end column loop

      return
      end Subroutine max8hrHour

 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to compute to W126 values, sums the weighted concentrations in ppm
c   between 8am & 8pm at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine w126(hrValues, dayValues, tzoffset, units, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Character*(*) units
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real ozone
      Real sum
      Real factor

C... set factor value to covert to ppm
      factor = 1.0
      if( INDEX(units,'ppb') .gt.0 ) factor = 0.001 ! to convert from ppb to ppm

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          sum = 0.0
          count = 0
          Do h = tzoffset(c,r)+8,tzoffset(c,r)+19     ! go from 8am to 7pm local time
            if( hrValues(h,c,r) .ge. 0.0 )  then
              count = count+1
              ozone = factor * hrValues(h,c,r)
              sum = sum + ozone / (1.0 + 4403.0 * EXP( -126*ozone ))
              endif                                                           
            endDo      ! end hour loop                                        

          if( partDay .or. count.ge.9 ) then
            dayValues(c,r) = sum / factor
            endif
          endDo      ! end row loop                                                         
        endDo      ! end column loop                                          
                                                                              
      return                                                                  
      end Subroutine w126 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Routine to sum hourly values >= 0.06ppm between 8am & 8pm at each cell
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine sum06(hrValues, dayValues, tzoffset, units, partDay)

      USE M3UTILIO

      IMPLICIT NONE

C... arguments
      Real hrValues(-36:66, NCOLS3D, NROWS3D)
      Real dayValues(NCOLS3D, NROWS3D)
      Integer tzoffset(NCOLS3D, NROWS3D)
      Character*(*) units
      Logical partDay

C... local variables
      Integer c, r, h, count
      Real sum
      Real minValue

C... set minimum cutoff value
      minValue = 0.06     ! default is 0.06 ppm
      if( INDEX(units,'ppb') .gt.0 ) minValue = 1000.0 * minValue

C... start loops to find daily sum
      Do c = 1, NCOLS3D
        Do r = 1, NROWS3D
          sum = 0.0
          count = 0
          Do h = tzoffset(c,r)+8,tzoffset(c,r)+20
            if( hrValues(h,c,r) .ge. minValue )  then
              count = count+1
              sum = sum + hrValues(h,c,r)
              endif
            endDo      ! end hour loop
          dayValues(c,r) = sum
          endDo      ! end row loop
        endDo      ! end column loop

      return     
      end Subroutine sum06

 

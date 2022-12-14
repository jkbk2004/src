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
C $Header: /project/work/rep/TOOLS/src/combine/module_evaluator.F,v 1.1.1.1 2005/07/27 12:55:20 sjr Exp $

C***********************************************************************
C
C  MODULE:  evaluates species expressions
C             
C***********************************************************************
      MODULE evaluator
      
      Public ::  evaluate
      
      Private

      Real, Private, Allocatable :: parseBuffer(:,:)
      Integer, Parameter, Private :: EXP_LEN = 1024    

      Integer, Private :: idate
      Integer, Private :: itime  
      Integer, Private :: ilayer
      Integer, Private :: isize
      
       Logical, Private :: eflag

      Interface            
         Subroutine getFld( record, delimiter, nth, del, field, exception )
            CHARACTER*(*), Intent( In  ) :: record
            CHARACTER*(*), Intent( In  ) :: delimiter
            CHARACTER,     Intent( Out ) :: del
            Integer,       Intent( In  ) :: nth
            CHARACTER*(*), Intent( Out ) :: field
            CHARACTER*(*), Optional, Intent( In ) :: exception
         End Subroutine getFld
         INTEGER FUNCTION getFldCount(record, delimiter, exception) Result(nfields)
            CHARACTER*(*), Intent( In ) :: record
            CHARACTER*(*), Intent( In ) :: delimiter
            CHARACTER*(*), Optional, Intent( In ) :: exception
         End FUNCTION getFldCount
         Subroutine LeftTrim( STRING )
            CHARACTER*(*), INTENT( INOUT ) :: STRING
         End Subroutine LeftTrim
         Subroutine RightTrim( STRING )
            CHARACTER*(*), INTENT( INOUT ) :: STRING
         End Subroutine RightTrim
         SUBROUTINE UCASE ( STR )
            CHARACTER, INTENT( INOUT ) :: STR*( * )
         END SUBROUTINE UCASE
         Subroutine replace( string, old, new )
            Character*(*), Intent( InOut ) :: string
            Character*(1), Intent( In    ) :: old    
            Character*(1), Intent( In    ) :: new    
         End Subroutine replace 
         SUBROUTINE Remove_WhiteSpaces (text)
            CHARACTER*(*), Intent( InOut ) :: text
         END SUBROUTINE Remove_WhiteSpaces
      End Interface


      CONTAINS


C  subroutine to evaluate species expression at date
C  returns buffer array values
      Subroutine evaluate(expression,jdate,jtime,jlayer,jsize,buffer)

      IMPLICIT NONE

      ! arguments
      Character*(*) expression
      Integer jdate, jtime
      Integer jlayer
      Integer jsize
      Real    buffer(jsize)

      ! local variables
      Character*(EXP_LEN) expresscp
      Character*(EXP_LEN) express
      Integer nparen
      Integer depth, maxdepth
      Integer i, n, pos1, pos2
      Character*(5) nstring
      Logical KSWIT

      ! set module variables
      idate = jdate
      itime = jtime
      ilayer = jlayer
      isize = jsize
      eflag = .False.
    
      ! make copy of expression to modify
      expresscp = expression
      Call Remove_WhiteSpaces(expresscp)
!      print*,Trim(expresscp),' ',Trim(expression)
!      expresscp = expression

      ! check for scientific notation (E+,E-,e+,e-) and replace with 10.0^
      call rmSciNot( expresscp )

      ! replace '+' characters inside [] brackets with '!' characters
      Call replace( expresscp, '+', '!' )

      ! replace '-' characters inside [] brackets with '#' characters
      Call replace( expresscp, '-', '#' )

      ! find number of parentheses and depth
      nparen = 0 
      depth = 0
      maxdepth = 0
      Do i=1,len_trim(expresscp)
        if( expresscp(i:i).eq.'(' ) then
          nparen = nparen + 1
          depth = depth + 1
          endif

        if( expresscp(i:i).eq.')' ) then
          depth = depth - 1
          endif
        
        if( depth.gt.maxdepth ) maxdepth = depth
        enddo

      !  check for unbalanced parentheses
      if( depth.ne.0 ) then
        write(*,'(/'' unbalanced parentheses in expression''/a)') trim(expresscp)
        stop
        endif

      ! allocate memory for parseBuffer if needed
      if( nparen.gt.0 ) then
        if( Allocated(parseBuffer) .and.
     &      SIZE(parseBuffer,DIM=2).lt.nparen ) then
          deAllocate(parseBuffer)
          endif

        if( .NOT.Allocated(parseBuffer) ) then
          Allocate( parseBuffer(isize,nparen) )
          endif

        parseBuffer = 0.0
        endif

      ! find depth of parentheses
      depth = maxDepth 
      Do n=1,nparen

        ! build buffer number as string
        write(nstring, '(i5)') n
        Call leftTrim(nstring)

        ! try to find parentheses at depth
        Call findDepth( expresscp, depth, pos1, pos2 )
            
        if( pos1.eq.0 ) then
          depth = depth - 1
          Call findDepth( expresscp, depth, pos1, pos2 )
          endif

        ! if parentheses found, evaluate sub expression
        if( pos1.gt.0 ) then

          ! extract expression within parentheses and
          ! evaluate to parsebuffer(1:isize,n)
          express = expresscp(pos1+1:pos2-1)
          call eval1(express, parsebuffer(1:isize,n) )

          ! replace expression within parentheses with "buffer[n]"
          express = ''
          if( pos1.gt.1 ) express = expresscp(1:pos1-1)
          express = TRIM(express) // 'buffer[' // TRIM(nstring) //
     &              ']' // TRIM(expresscp(pos2+1:))
          expresscp = express 
          endif 
        enddo

      call eval1(expresscp, buffer)
      
      end Subroutine evaluate


C  subroutine to replace scientific notation strings
      Subroutine rmSciNot(expression)

      IMPLICIT NONE

      Character*(*) expression

      Character*(2) estring(4)
      Character*(7) pstring(4)

      Integer n, i, pos, pos1, pos2

      Data estring/'E+','e+','E-','e-'/
      Data pstring/'*(10.0^', '*(10.0^', '/(10.0^', '/(10.0^'/

      do n=1,4
        do while( index(expression, estring(n)) .gt. 0 )
          pos = index(expression, estring(n))

          ! search for start of number starting at pos-1 and working back
          pos1 = pos-1
          do i=pos-1,1,-1
            if( index('0123456789.',expression(i:i)) .eq. 0 ) then
              EXIT
              endif
            pos1 = i
            enddo

          ! search for end of number starting at pos+2
          do i=pos+2,pos+12
            if( index('0123456789.',expression(i:i)) .eq. 0 ) then
              pos2=i
              EXIT
              endif
            enddo

          if( pos1 .eq. 1 ) then
            expression = '(' // expression(1:pos-1) // pstring(n) // expression(pos+2:pos2-1)
     &                // '))' // expression(pos2:)
            endif

          if( pos1 .gt. 1 ) then
            expression = expression(1:pos1-1) // '(' // expression(pos1:pos-1) //
     &                 pstring(n) // expression(pos+2:pos2-1) // '))' // expression(pos2:)
            endif

          enddo
        enddo

      return
      end Subroutine rmSciNot



C  subroutine to find location of parentheses depth
      Subroutine findDepth(expression, depth, pos1, pos2)

      IMPLICIT NONE

      Character*(*) expression
      Integer depth, pos1, pos2

      Integer i, dep

      pos1 = 0
      pos2 = 0
      dep = 0

      ! try to find parentheses at depth
      Do i = 1, len_trim(expression)  
        if( expression(i:i).eq.'(' ) then
            dep = dep+1
            if(dep.eq.depth) pos1 = i
            endif

          if( expression(i:i).eq.')' ) then
            if(dep.eq.depth) then
              pos2 = i
              return
              endif
            dep = dep-1
            endif           
         enddo

      return
      end Subroutine findDepth


C  subroutine to return buffer array value
      Subroutine getBuffer(field, buffer)
      IMPLICIT NONE

      Character*(*), Intent( InOut ) :: field
      Real,          Intent( InOut ) :: buffer(isize)
      
      Real,  Parameter :: Max_Exponential = 82.8931
      
      Integer pos1, pos2, nbuf, status
      Character*(10) string
      Character*(10) func
      Logical KSWIT
      Logical SHUT3

      Call leftTrim(field)

      ! parse field to find buffer number
      pos1 = index(field, '[') 
      pos2 = index(field, ']',.true.) 

      if(pos1.le.0 .or. pos1.ge.pos2) then
        write(*,'(/''**ERROR**  Invalid syntax in field: '',a)') trim(field)
        KSWIT = SHUT3()
        stop
        endif

      if(field(pos2+1:) .ne. ' ') then
        write(*,'(/''**ERROR**  Invalid syntax in field: '',a)') trim(field)
        KSWIT = SHUT3()
        stop
        endif

      string = field(pos1+1:pos2-1)
      read(string,'(i10)',iostat=status) nbuf      
      if(status .ne. 0) then
        write(*,'(/''**ERROR**  Invalid syntax in field: '',a)') trim(field)
        KSWIT = SHUT3()
        stop
        endif

      buffer = parsebuffer(1:isize,nbuf)      

      ! check for function
      pos1 = index(field, 'buffer[') 
      Call UCASE(field)

      if( pos1.gt.1 ) then       
        func = field(1:pos1-1)
        print*,Trim(field)
        If( func.eq.'LOG' ) Then
          If( MinVal(buffer) .Le. 0.0 )Then
            write(*,'(/''**ERROR** evaluating LOG of '',es12.4)' )
     &      MinVal(buffer)
            KSWIT = SHUT3()
            stop
          End If
          If( MaxVal(buffer) .Gt. 1.0E36 )Then
            write(*,'(/''**WARNING** evaluating LOG of '',es12.4)' )
     &      MaxVal(buffer)
          End If
          buffer = LOG(buffer)
          return
        End If
        If( func.eq.'LOG10' ) Then
          If( MinVal(buffer) .Le. 0.0 )Then
            write(*,'(/''**ERROR** evaluating LOG10 of '',es12.4)' )
     &      MinVal(buffer)
            KSWIT = SHUT3()
            stop
          End If
          If( MaxVal(buffer) .Gt. 1.0E36 )Then
            write(*,'(/''**WARNING** evaluating LOG10 of '',es12.4)' )
     &      MaxVal(buffer)
          End If
          buffer = LOG10(buffer)
          return
        End If
        If( func.eq.'EXP' ) Then
          If( MaxVal(buffer) .Ge. Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating EXP of '',es12.4)' )
     &      MaxVal(buffer)
            KSWIT = SHUT3()
            stop
          End If
          If( MinVal(buffer) .Le. -Max_Exponential )Then
            write(*,'(/''**WARNING** evaluating EXP of '',es12.4)' )
     &      MinVal(buffer)
          End If
          Where( buffer .Le. -Max_Exponential )
             buffer = 0.0
          Else Where
             buffer = EXP(buffer)
          End Where
          return
        End If
        If( func.eq.'SIN' ) Then
          buffer = SIN(buffer)
          return
        End If
        If( func.eq.'COS' ) Then
          buffer = COS(buffer)
          return
        End If
        If( func.eq.'TAN' ) Then
          buffer = TAN(buffer)
          return
        End If
        If( func.eq.'ASIN' ) Then
          If( MinVal(buffer) .Lt. -1.0 )Then
            write(*,'(/''**ERROR** evaluating ASIN of '',es12.4)' )
     &      MinVal(buffer)
            eflag = .True.
          End If
          If( MaxVal(buffer) .Gt. 1.0 )Then
            write(*,'(/''**ERROR** evaluating ASIN of '',es12.4)' )
     &      MaxVal(buffer)
            eflag = .True.
          End If
          If( eflag )Then
            KSWIT = SHUT3()
            stop
          End If
          buffer = ASIN(buffer)
          return
        End If
        If( func.eq.'ACOS' ) Then
          If( MinVal(buffer) .Lt. -1.0 )Then
            write(*,'(/''**ERROR** evaluating ACOS of '',es12.4)' )
     &      MinVal(buffer)
            eflag = .True.
          End If
          If( MaxVal(buffer) .Gt. 1.0 )Then
            write(*,'(/''**ERROR** evaluating ACOS of '',es12.4)' )
     &      MaxVal(buffer)
            eflag = .True.
          End If
          If( eflag )Then
            KSWIT = SHUT3()
            stop
          End If
          buffer = ACOS(buffer)
          return
        End If
        If( func.eq.'ATAN' ) Then
          buffer = ATAN(buffer)
          return
        End If
        If( func.eq.'SINH' ) Then
          If( MaxVal(buffer) .Ge. Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating SINH of '',es12.4)' )
     &      MaxVal(buffer)
            eflag = .True.
          End If
          If( MinVal(buffer) .Le. -Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating SINH of '',es12.4)' )
     &      MinVal(buffer)
            eflag = .True.
          End If
          If( eflag )Then
            KSWIT = SHUT3()
            stop
          End If
          buffer = SINH(buffer)
          return
        End If
        If( func.eq.'COSH' ) Then
          If( MaxVal(buffer) .Ge. Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating COSH of '',es12.4)' )
     &      MaxVal(buffer)
            eflag = .True.
          End If
          If( MinVal(buffer) .Le. -Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating COSH of '',es12.4)' )
     &      MinVal(buffer)
            eflag = .True.
          End If
          If( eflag )Then
            KSWIT = SHUT3()
            stop
          End If
          buffer = COSH(buffer)
          return
        End If
        If( func.eq.'TANH' ) Then
          If( MaxVal(buffer) .Ge. Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating TANH of '',es12.4)' )
     &      MaxVal(buffer)
            KSWIT = SHUT3()
            stop
          End If
          If( MinVal(buffer) .Le. -Max_Exponential )Then
            write(*,'(/''**ERROR** evaluating TANH of '',es12.4)' )
     &      MinVal(buffer)
            eflag = .True.
          End If
          If( eflag )Then
            KSWIT = SHUT3()
            stop
          End If
          buffer = TANH(buffer)
          return
        End If
        If( func.eq.'SQRT' ) Then
          If( MinVal(buffer) .Lt. 0.0 )Then
            write(*,'(/''**ERROR** evaluating SQRT of '',es12.4)' )
     &      MinVal(buffer)
            KSWIT = SHUT3()
            stop
          End If
          buffer = SQRT(buffer)
          return
        End If
        If( func.eq.'AINT' .Or. func.eq.'INT') Then
          buffer = AINT(buffer)
          return
        End If
        If( func.eq.'ABS' ) Then
          buffer = ABS(buffer)
          return
        End If

        write(*,'(/''**ERROR** Invalid function name: '',a)') trim(func)
        KSWIT = SHUT3()
        stop
        endif

      return
      end Subroutine getBuffer


C  subroutine to evaluate species expression (parses conditional statment if needed)
C   X = (y[1]>10) ? 10 : y[1]
C
      Subroutine eval1(expression, buffer)

      IMPLICIT NONE

      ! arguments
      Character*(*) expression
      Real buffer(isize)

      ! functions
!      Integer getFldCount
 
      ! local variables
      Logical, Allocatable :: flags(:)
      Real, Allocatable :: value1(:)
      Real, Allocatable :: value2(:)
      Character*(EXP_LEN) field
      Character operator
      Integer nmajor
      Integer i
      Logical badopr


      ! parse major fields (?:)
      nmajor = getFldCount(expression, '?:')

      ! if conditional 
      if( nmajor.eq.3 ) then 
        Allocate( flags(isize), value1(isize), value2(isize) )
        badopr = .false.

        call getFld( expression, '?:', 1, operator, field ) 
        if(operator.ne.'?') badopr = .true.
        call eval1b( field, flags)

        call getFld( expression, '?:', 2, operator, field ) 
        if(operator.ne.'?') badopr = .true.
        call eval2( field, value1)

        call getFld( expression, '?:', 3, operator, field ) 
        if(operator.ne.':') badopr = .true.
        call eval2( field, value2)

        if( badopr ) then
          Write(*,'(/''**Error** Syntax error encountered in conditional expression: '',a)') trim(expression)
          stop
          endif

        ! set buffer values 
        do i=1,isize
          if( flags(i) ) then
            buffer(i) = value1(i)
          else
            buffer(i) = value2(i)
            endif 
          enddo 

        Deallocate (flags, value1, value2)
        return
        endif

      ! if no conditional
      if( nmajor.eq.1 ) then
        call eval2( trim(expression), buffer )
        return
        endif

      ! syntax error
      Write(*,'(/''**Error** Syntax error encountered at: '',a)') trim(expression)
      stop   
      end Subroutine eval1


C  subroutine to evaluate condition expression (called from eval1) 
      Subroutine eval1b(expression, flags)

      IMPLICIT NONE

      ! arguments
      Character*(*) expression
      Logical flags(isize)

      ! functions
!      Integer getFldCount
 
      ! local variables
      Real, Allocatable :: value1(:)
      Real, Allocatable :: value2(:)
      Character*(EXP_LEN) field
      Character operator
      Integer nflds
      Integer i


      ! verify that expression contains a parse major fields (<=>)
      nflds = getFldCount(expression, '<=>')
      if( nflds.eq.0 ) then
        Write(*,'(/''**Error** Syntax error encountered in conditional: '',a)') trim(expression)
        stop
        endif

      ! parse conditional expression
      Allocate( value1(isize), value2(isize) )

      ! determine conditional operator is <=
      if( index(expression,'<=').gt.0 ) then
        call getFld( expression, '<=', 1, operator, field ) 
        call eval2( field, value1)
        call getFld( expression, '<=', 3, operator, field ) 
        call eval2( field, value2)
        flags = ( value1 .le. value2 )
        Deallocate (value1, value2)
        return
        endif

      ! determine conditional operator is >=
      if( index(expression,'>=').gt.0 ) then
        call getFld( expression, '>=', 1, operator, field ) 
        call eval2( field, value1)
        call getFld( expression, '>=', 3, operator, field ) 
        call eval2( field, value2)
        flags = ( value1 .ge. value2 )
        Deallocate (value1, value2)
        return
        endif 

      ! determine conditional operator is >
      if( index(expression,'>').gt.0 ) then
        call getFld( expression, '>', 1, operator, field ) 
        call eval2( field, value1)
        call getFld( expression, '>', 2, operator, field ) 
        call eval2( field, value2)
        flags = ( value1 .gt. value2 )
        Deallocate (value1, value2)
        return
        endif 

      ! determine conditional operator is <
      if( index(expression,'<').gt.0 ) then
        call getFld( expression, '<', 1, operator, field ) 
        call eval2( field, value1)
        call getFld( expression, '<', 2, operator, field ) 
        call eval2( field, value2)
        flags = ( value1 .lt. value2 )
        Deallocate (value1, value2)
        return
        endif 

      ! determine conditional operator is =
      if( index(expression,'=').gt.0 ) then
        call getFld( expression, '=', 1, operator, field ) 
        call eval2( field, value1)
        call getFld( expression, '=', 2, operator, field ) 
        call eval2( field, value2)
        flags = ( value1 .eq. value2 )
        Deallocate (value1, value2)
        return
        endif 

      ! syntax error
      Write(*,'(/''**Error** Syntax error encountered: '',a)') trim(expression)
      stop
    
      end Subroutine eval1b


C  subroutine to evaluate species expression (parses major fields (+-))
      Subroutine eval2(expression, buffer)

      IMPLICIT NONE

      ! arguments
      Character*(*), Intent( In )    :: expression
      Real,          Intent( InOut ) :: buffer(isize)

      ! local variables
      Real,    Allocatable :: value(:)
      Character*(EXP_LEN)  :: field
      Character operator
      Integer nmajor
      Integer n

      buffer = 0.0D0
      Allocate ( value(isize) )

      ! parse major fields (+-)
      nmajor = getFldCount(expression, '+-', '*/^')

      ! loop thru and parse each major field and evaluate
      Do n=1,nmajor

        call getFld( expression, '+-', n, operator, field, '*/^' ) 
    
        If( field.eq.' ' ) Then
          value = 0.0D0
        else
          call eval3b( field, value)
        End If

        If( operator.eq.'+' ) Then
          buffer = buffer + value
        else
          buffer = buffer - value
        End If

        End Do

      Deallocate (value)
      return
      end Subroutine eval2


C  routine to compute a field of the expression (parses minor fields (*/^))
      Subroutine eval3(expression, value)
      
      IMPLICIT NONE

      ! arguments
      
      CHARACTER*(*), Intent( In    ) :: expression
      Real,          Intent( InOut ) :: value(isize)

      Logical SHUT3

      ! local variables
      Real,    allocatable :: specValue(:)
      Character*(EXP_LEN) field
      Character      operator   
      Integer n, m, nflds, status
      Integer pos1, pos2, fnum
      Character*(16) funcName
      Character*(16) specName
      Real    constant
      Logical KSWIT

      Allocate ( specValue(isize) )
      nflds = getFldCount(trim(expression), '*/^')
      value = 1.0
         
      Do n=1,nflds
        call getFld( trim(expression), '*/^', n, operator, field ) 

        print*,Trim(field)

        ! check for buffer array
        If( index(field,'buffer[') .gt.0 ) Then
          Call getBuffer(field, specValue)
          If( operator.eq.'*' ) value = value * specValue
          If( operator.eq.'/' ) value = value / specValue
          If( operator.eq.'^' ) value = value ** specValue
          cycle
        End If
 
        ! check for species argument (special functions)
        If( index(field,'[') .gt.0 ) Then
  
          ! switch ! and # characters within [] brackets back to + and - characters
          Call replace(field, '!', '+')
          Call replace(field, '#', '-')

          ! parse field between [ ] and check If number or species name
          pos1 = index(field, '[')
          pos2 = index(field, ']',.true.)
          specName = field(pos1+1:pos2-1)
  
          read(specName,'(i16)',iostat=status) fnum

          If( status.eq.0 ) Then    !! number found
            Call readSpecies(field, specValue)
            If( operator.eq.'*' ) value = value * specValue
            If( operator.eq.'/' ) value = value / specValue
            If( operator.eq.'^' ) value = value ** specValue
            cycle
          End If    !! contains '['
        End If
      !try to read field as number
        read(field,'(f20.0)',iostat=status) constant

        If( status.eq.0 ) Then
             If( operator.eq.'*' ) value = value * constant
             If( operator.eq.'/' ) value = value / constant
             If( operator.eq.'^' ) value = value ** constant
        Else
             Write(*,'(''**Error** Invalid field encountered:'',a)') field
             stop 
        End If
      End Do

      Deallocate (specValue)
      return
      end Subroutine eval3

C  routine to compute a field of the expression (parses minor fields (*/^))
      Subroutine eval3b(expression, value)
      
      IMPLICIT NONE

      ! arguments
      
      CHARACTER*(*), Intent( In    ) :: expression
      Real,          Intent( InOut ) :: value(isize)

      Logical SHUT3

      ! local variables
      Real,     allocatable :: specValue(:)
      Character*(EXP_LEN) field
      Character      operator   
      Integer n, nflds

      Allocate ( specValue(isize) )
!      nflds = getFldCount(trim(expression), '*/^')
      nflds = getFldCount(trim(expression), '*/')
      value = 1.0
         
      SpecValue = 0.0D0
      Do n=1,nflds
!        call getFld( trim(expression), '*/^', n, operator, field ) 
        call getFld( trim(expression), '*/', n, operator, field ) 
        call eval4(field, specValue)
!        call GetValue(field, specValue)
!        print*,operator,value,specValue

        If( operator.eq.'*' ) value = value * specValue
        If( operator.eq.'/' ) value = value / specValue
!        If( operator.eq.'^' ) value = value ** specValue

      End Do

      Deallocate (specValue)
      return
      end Subroutine eval3b
      Subroutine eval4(expression, value)
      
      IMPLICIT NONE

      ! arguments
      
      CHARACTER*(*), Intent( In    ) :: expression
      Real,           Intent( InOut ) :: value(:)

      Real,    allocatable :: specValue(:)
      Real,    allocatable :: specPower(:)
      Character*(EXP_LEN)  :: field
      Character(1)         :: operator   
      Integer              :: n, pos1, nflds
      Real    Factor

      
        

      nflds = getFldCount(trim(expression), '^')
!      value = 1.0
!...No exponents found
      Allocate ( specValue(isize) )
      If( nflds .Eq. 1 )Then
          field = expression
          call GetValue(field, specValue)
           value = specValue
          Deallocate (specValue)
          Return
      End If
!...check if correct number of exponents found      
!      If( mod(nflds,2) .Ne. 0 )Then
!          eflag = .True.
!          Write(6,*)'Incorrect number of exponents in Formula: ', Trim(formula)
!          Return
!      End If
!...compute fields with 
      Allocate ( specPower(isize) )
      Value = 1.0D0
      n = nflds 
      call getFld( trim(expression), '^', n, operator, field ) 
      n = n - 1
      call GetValue(field, specPower)
      Do 
         call getFld( trim(expression), '^', n, operator, field )
         n = n - 1
          If( field(1:1) .Eq. '-' )Then
            field  = field(2:)
            Factor = -1.0D0
          Else
            Factor = 1.0D0
          End If 
         call GetValue(field, specValue)
         Value = Factor * specValue**SpecPower
         If( n .Lt. 1)EXIT
         SpecPower = Value
      End Do

      Deallocate (specValue)
      Deallocate (specPower)

      end Subroutine eval4

C  routine to compute a field of the expression (parses minor fields (*/^))
      Subroutine GetValue(expression, value)
      
      IMPLICIT NONE

      ! arguments
      
      CHARACTER*(*), Intent( In    ) :: expression
      Real,          Intent( InOut ) :: value(isize)

      Logical SHUT3

      ! local variables
      Real,    allocatable :: specValue(:)
      Character*(EXP_LEN) field
      Character      operator   
      Integer n, m, nflds, status
      Integer pos1, pos2, fnum
      Character*(16) funcName
      Character*(16) specName
      Real    constant
      Real    Factor

      
        If( expression(1:1) .Eq. '-' )Then
            field    = expression(2:)
            Factor = -1.0D0
        Else If( expression .Eq. '+' )Then
            field    = expression(2:)
            Factor = 1.0D0
        Else            
            field = expression
            Factor = 1.0D0
        End If 

!...check for buffer array
        If( index(field,'buffer[') .gt.0 ) Then
          Call getBuffer(field, Value)
          Value = Factor * Value
          Return
        End If
        
!...check for species argument (special functions)
        If( index(field,'[') .gt.0 ) Then
  
          ! switch ! and # characters within [] brackets back to + and - characters
          Call replace(field, '!', '+')
          Call replace(field, '#', '-')

          ! parse field between [ ] and check If number or species name
          pos1 = index(field, '[')
          pos2 = index(field, ']',.true.)
          specName = field(pos1+1:pos2-1)
  
          read(specName,'(i16)',iostat=status) fnum

          If( status.eq.0 ) Then    !! number found
            Call readSpecies(field, Value)
            Value = Factor * Value
            Return
          Else
            funcName = field(1:pos1-1)
            Call UCASE(funcName)
             status = -1
             If( funcName .eq. 'FAVG' )  Call avgSpecies(  specName , Value, status )
             If( funcName .eq. 'FSDEV' ) Call sdevSpecies( specName , Value, status )
             If( funcName .eq. 'FMAX' )  Call maxSpecies(  specName , Value, status )
             If( funcName .eq. 'FMIN' )  Call minSpecies(  specName , Value, status )
             If( status.eq.0 ) Then  
                Value = Factor * Value
                Return 
             End If
          End If    !! contains '['
        End If
!...try to read field as number
        read(field,'(f20.0)',iostat=status) constant

        If( status.eq.0 ) Then
          value = Factor * constant
          Return
        Else
             Write(*,'(''**Error** Invalid field encountered:'',a)') field
             stop 
        End If

      return
      end Subroutine GetValue


C  Routine to read species value array for given date and time
      Subroutine readSpecies( field, specValue)

      USE M3FILES
      USE M3UTILIO

      IMPLICIT NONE

      ! arguments
      Character*(*) field
      Real specValue(isize)

      ! local variables
      Integer pos1, pos2, status
      Character*(16) specName
      Character*(16) fileName 
      Character*(10) numfld
      Integer fnum
      Integer kdate, ktime
      Integer edate, etime, m
      Logical KSWIT


      ! parse field into species name and file number
      pos1 = index(field, '[') 
      pos2 = index(field, ']',.true.) 
      specName = field(1:pos1-1)

      if(pos1.le.0 .and. pos1.ge.pos2) then
        Write(*,'(''**ERROR** Invalid file number for species '',a)') trim(specName)
        KSWIT = SHUT3()
        stop 
        endif

      ! parse file number
      numfld = field(pos1+1:pos2-1)

      ! read file number from numfld
      read(numfld,*,iostat=status) fnum
      if( status.ne.0 ) then
        Write(*,'(/''**ERROR** Invalid file number for species: '',a)') trim(specName) 
        Write(*,'(''   file number:'',a)') trim(numfld) 
        KSWIT = SHUT3()
        stop 
        endif

      kdate = idate
      ktime = itime

      !! Check numfld for +- sign, to read values of next or previous time step
      if( INDEX(numfld,'-').gt.0 .or. INDEX(numfld,'+').gt.0 ) then
        fnum = ABS( fnum )
        if( .NOT. getDESC( fnum ) ) then
          Write( *, '(''**Error** While running getDESC on file '',i5)') fnum
          stop  
          endif

        !! adjust date/time to read previous timestep
        if( INDEX(numfld,'-') .gt. 0 ) then
          if( SECSDIFF (kDate, ktime, SDATE3D, STIME3D) .lt. 0 ) then
            Call NEXTIME( kdate, ktime, -TSTEP )
            endif
          endif

        !! adjust date/time to read next timestep
        if( INDEX(numfld,'+') .gt. 0 ) then

          ! compute ending time of file
          edate = SDATE3D
          etime = STIME3D
          do m = 1, MXREC3D-1
            Call Nextime(edate, etime, TSTEP3D)
            enddo

          if( SECSDIFF (kDate, ktime, edate, etime) .gt. 0 ) then
            Call NEXTIME( kdate, ktime, TSTEP )
            endif
          endif
        endif   !! condition to adjust timestep for read

      !! call routine to read species values from file fnum
      status = 0
      !! write(*,'(''reading values for '',a,'' at '',2i8)') TRIM(specName), kdate,ktime
      Call ReadValues( fnum, specName, ilayer, kdate, ktime, isize,    
     &                   specValue, status)

      !! check read status
      if( status.ne.0 ) then
        if( fnum.eq.0 ) filename = 'OUTFILE'
        if( fnum.gt.0 ) filename = M3FILENAME(fnum)

        Write(*,'(/''**ERROR** Invalid syntax for field: '',a)') trim(field)     
        Write(*,'(/''**ERROR** Cannot read '',a,'' from '',a)')
     &            trim(specName), trim(fileName)
        KSWIT = SHUT3()
        stop 
        endif  

      return
      end Subroutine readSpecies  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Routine to compute average species value at each cell from all input files
      Subroutine avgSpecies( specName, specValue, status)

      USE M3FILES
      USE M3UTILIO

      IMPLICIT NONE

      ! arguments
      Character*(*) specname
      Real specValue(isize)
      Integer status

      ! local variables
      Integer n, i                 
      Real, Allocatable :: values(:,:)
      Logical KSWIT

      ! allocate arrays
      Allocate( values(isize, N_M3FILES) )

      ! read species values from all input files
      do n = 1, N_M3FILES

        status = 0  
        Call ReadValues( n, specName, ilayer, idate, itime, isize, 
     &                   values(:,n), status)

        !! check read status
        if( status.ne.0 ) then
          Write(*,'(/''**ERROR** Cannot read '',a,'' from '',a)')
     &              trim(specName), trim(M3FILENAME(n))
          KSWIT = SHUT3()
          stop
          endif

        enddo    !! read loop

      !! compute averages
      do i = 1, isize
        specValue(i) = 0.0
        do n = 1, N_M3FILES
          specValue(i) = specValue(i) + values(i,n)
          enddo
        specValue(i) = specValue(i) / N_M3FILES 
        enddo  

      ! deallocate arrays
      DeAllocate( values )

      status = 0
      return
      End Subroutine avgSpecies  


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Routine to find minimum species value at each cell from all input files
      Subroutine minSpecies( specName, specValue, status)

      USE M3FILES
      USE M3UTILIO

      IMPLICIT NONE

      ! arguments
      Character*(*) specname
      Real specValue(isize)
      Integer status

      ! local variables
      Integer n, i                 
      Real, Allocatable :: values(:,:)
      Logical KSWIT

      ! allocate arrays
      Allocate( values(isize, N_M3FILES) )

      ! read species values from all input files
      do n = 1, N_M3FILES

        status = 0  
        Call ReadValues( n, specName, ilayer, idate, itime, isize, 
     &                   values(:,n), status)

        !! check read status
        if( status.ne.0 ) then
          Write(*,'(/''**ERROR** Cannot read '',a,'' from '',a)')
     &              trim(specName), trim(M3FILENAME(n))
          KSWIT = SHUT3()
          stop
          endif

        enddo    !! read loop

      !!  find minimums
      do i = 1, isize
        specValue(i) = values(i,1)
        do n = 2, N_M3FILES
          if(values(i,n) .lt. specValue(i)) specValue(i) = values(i,n)
          enddo
        enddo  

      ! deallocate arrays
      DeAllocate( values )

      status = 0
      return
      End Subroutine minSpecies  



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Routine to find maximum species value at each cell from all input files
      Subroutine maxSpecies( specName, specValue, status)

      USE M3FILES
      USE M3UTILIO

      IMPLICIT NONE

      ! arguments
      Character*(*) specname
      Real specValue(isize)
      Integer status

      ! local variables
      Integer n, i                 
      Real, Allocatable :: values(:,:)
      Logical KSWIT

      ! allocate arrays
      Allocate( values(isize, N_M3FILES) )

      ! read species values from all input files
      do n = 1, N_M3FILES

        status = 0  
        Call ReadValues( n, specName, ilayer, idate, itime, isize, 
     &                   values(:,n), status)

        !! check read status
        if( status.ne.0 ) then
          Write(*,'(/''**ERROR** Cannot read '',a,'' from '',a)')
     &              trim(specName), trim(M3FILENAME(n))
          KSWIT = SHUT3()
          stop
          endif

        enddo    !! read loop

      !!  find minimums
      do i = 1, isize
        specValue(i) = values(i,1)
        do n = 2, N_M3FILES
          if(values(i,n) .gt. specValue(i)) specValue(i) = values(i,n)
          enddo
        enddo  

      ! deallocate arrays
      DeAllocate( values )

      status = 0
      return
      End Subroutine maxSpecies  


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Routine to compute the standard deviation value at each cell from all input files
      Subroutine sdevSpecies( specName, specValue, status)

      USE M3FILES
      USE M3UTILIO

      IMPLICIT NONE

      ! arguments
      Character*(*) specname
      Real specValue(isize)
      Integer status

      ! local variables
      Integer n, i   
      Real xtotal, x2total , var
      Real, Allocatable :: values(:,:)
      Logical KSWIT

      !! if the number of files == 1, then set standard deviation values to zero
      if( N_M3FILES .le. 0 ) then
        specValue = 0.0
        status = 0
        return
        endif

      ! allocate arrays
      Allocate( values(isize, N_M3FILES) )

      ! read species values from all input files
      do n = 1, N_M3FILES

        status = 0  
        Call ReadValues( n, specName, ilayer, idate, itime, isize, 
     &                   values(:,n), status)

        !! check read status
        if( status.ne.0 ) then
          Write(*,'(/''**ERROR** Cannot read '',a,'' from '',a)')
     &              trim(specName), trim(M3FILENAME(n))
          KSWIT = SHUT3()
          stop
          endif

        enddo    !! read loop

      !!  find minimums
      do i = 1, isize
        xtotal = 0.0
        x2total = 0.0
        do n = 1, N_M3FILES
          xtotal = xtotal + values(i,n) 
          x2total = x2total + values(i,n)**2 
          enddo

        var = (N_M3FILES*x2total - xtotal**2) / (N_M3FILES * (N_M3FILES-1))
        specValue(i) = SQRT(var)

        enddo  

      ! deallocate arrays
      DeAllocate( values )

      status = 0
      return
      End Subroutine sdevSpecies  

      END MODULE evaluator

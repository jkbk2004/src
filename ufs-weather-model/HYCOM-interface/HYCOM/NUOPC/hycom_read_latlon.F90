      module hycom_read_latlon
      implicit none

      contains

      subroutine get_coord(lat_p,lon_p,lat_q,lon_q,itdm,jtdm,rc)
      implicit none

      real, dimension(itdm,jtdm), intent(inout) :: lat_p, lon_p
      real, dimension(itdm,jtdm), intent(inout) :: lat_q, lon_q
      integer, intent(out) :: rc
      real*4, dimension(itdm,jtdm) :: tmp1, tmp2
      real*4, allocatable :: pad(:)

      integer :: npad
      integer :: ios,nrecl
      integer :: i,j
      character(len=240) :: cfilea
      integer :: itdm,jtdm

      rc = 0 ! success

      cfilea = 'regional.grid.a'

      npad = 4096 - MOD(itdm*jtdm,4096)
      if(npad.eq.4096) npad=0

      allocate(pad(npad))

      INQUIRE( IOLENGTH=nrecl) tmp1,pad

#if defined(ENDIAN_IO)
      open(unit=11,file=cfilea, form='unformatted', status='old', &
               access='direct', recl=nrecl, convert="BIG_ENDIAN", &
               iostat=ios)
#else
      open(unit=11,file=cfilea, form='unformatted', status='old', &
               access='direct', recl=nrecl, &
               iostat=ios)
#endif

      IF (ios.ne.0) THEN
        print *,"error in reading regional.grid.a"
        call exit(1)
      endif

      read(11,rec=1,iostat=ios) tmp1
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, plon"
        call exit(2)
      endif
      lon_p = tmp1

      read(11,rec=2,iostat=ios) tmp1
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, plat"
        call exit(3)
      endif
      lat_p = tmp1

      read(11,rec=3,iostat=ios) tmp1
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, qlon"
        call exit(2)
      endif
      lon_q = tmp1

      read(11,rec=4,iostat=ios) tmp1
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, qlat"
        call exit(3)
      endif
      lat_q = tmp1

      do j=1, jtdm
        do i=1, itdm
          if (lon_p(i,j).ge.360.) lon_p(i,j) = lon_p(i,j)-360.
          if (lon_q(i,j).ge.360.) lon_q(i,j) = lon_q(i,j)-360.
        enddo
      enddo

      print *,'get_coord, plat, min, max=', &
          minval(lat_p),maxval(lat_p)
      print *,'get_coord, plon, min, max=', &
          minval(lon_p),maxval(lon_p)
      print *,'get_coord, qlat, min, max=', &
          minval(lat_q),maxval(lat_q)
      print *,'get_coord, qlon, min, max=', &
          minval(lon_q),maxval(lon_q)

      ! deallocate temporary arrays
      if (allocated(pad)) deallocate(pad)

      return
      end subroutine get_coord

      end module hycom_read_latlon

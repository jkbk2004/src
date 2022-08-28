#if defined (USE_NUOPC_CESMBETA) || (ESPC_COUPLE)
#define USE_NUOPC_GENERIC 1
#endif
      module mod_import

      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays

      implicit none

      public hycom_imp_reset
      public hycom_imp_mrg_diag
      public hycom_imp_mrg
      public hycom_imp_mrg_latflx
      public hycom_imp_mrg_sensflx

      contains

! --- reset imp_merge flag
      subroutine hycom_imp_reset(reset)
      logical :: reset
      character(*),parameter :: rname = 'hycom_imp_reset'

#if defined (USE_NUOPC_GENERIC)
        if(reset) then
          imp_merge = 1.d0
        else
          imp_merge = 0.d0
        endif
        cpl_merge = reset
#endif
      end subroutine hycom_imp_reset

! --- print merge import data diagnostics
      subroutine hycom_imp_mrg_diag()
      character(*),parameter :: rname = 'hycom_imp_mrg_diag'
      integer                :: i,j,jja
      real                   :: mrg_cnt
      integer                :: tlb(2), tub(2)

#if defined (USE_NUOPC_GENERIC)
#if defined(ARCTIC)
!   arctic (tripole) domain, top row is replicated (ignore it)
        jja=min(jj,(jtdm-1-j0))
#else
        jja=jj
#endif

        ! count merge cells ignoring halos
        mrg_cnt=0
        do j=1, jja
        do i=1, ii
          if (imp_merge(i,j).ne.0) mrg_cnt=mrg_cnt+1
        enddo
        enddo
        call xcsumr(mrg_cnt,1)
        if (mnproc.eq.1) then
          write(*,'(A,L1)') rname//',cpl_merge=',cpl_merge
          write(*,'(A,E23.15)') rname//',imp_merge_cnt=',mrg_cnt
          if(cpl_taux) write(*,'(A)') rname//" merge imp_taux"
          if(cpl_tauy) write(*,'(A)') rname//" merge imp_tauy"
          if(cpl_wndspd.or.calc_wndspd) &
            write(*,'(A)') rname//" merge imp_wndspd"
          if(cpl_ustara) write(*,'(A)') rname//" merge imp_ustara"
          if(cpl_airtmp) write(*,'(A)') rname//" merge imp_airtmp"
          if(cpl_vapmix) write(*,'(A)') rname//" merge imp_vapmix"
          if(cpl_precip) write(*,'(A)') rname//" merge imp_precip"
          if(cpl_surtmp) write(*,'(A)') rname//" merge imp_surtmp"
          if(cpl_seatmp) write(*,'(A)') rname//" merge imp_seatmp"
          if(cpl_swflx_net.or.cpl_swflx_net2down.or.cpl_swflxd) &
            write(*,'(A)') rname//" merge imp_swflx"
          if(cpl_lwflx_net.or.cpl_lwflx_net2down.or.cpl_lwflxd) &
            write(*,'(A)') rname//" merge imp_lwdflx"
          if(cpl_u10) write(*,'(A)') rname//" merge imp_wndspx"
          if(cpl_v10) write(*,'(A)') rname//" merge imp_wndspy"
          if(cpl_mslprs) write(*,'(A)') rname//" merge imp_mslprs"
          if(calc_radflx) write(*,'(A)') rname//" merge imp_radflx"
          if(cpl_latflx) write(*,'(A)') rname//" merge imp_latflx"
          if(cpl_sensflx) write(*,'(A)') rname//" merge imp_sensflx"
        endif
#else
        write(lp,'(/ a,a /)') 'error - ', &
          rname//' merge requires USE_NUOPC_CESMBETA or ESPC_COUPLE'
        call flush(lp)
        call xcstop(rname)
               stop rname
#endif
      end subroutine hycom_imp_mrg_diag

! --- merge import data where imp_merge is true
      subroutine hycom_imp_mrg()
      character(*),parameter :: rname = 'hycom_imp_mrg'
      integer                :: i,j,jja
      integer                :: tlb(2), tub(2)

#if defined (USE_NUOPC_GENERIC)
#if defined(ARCTIC)
!   arctic (tripole) domain, top row is replicated (ignore it)
        jja=min(jj,(jtdm-1-j0))
#else
        jja=jj
#endif

        if(natm.eq.2) then
          tlb(1)=lbound(imp_merge,1)
          tlb(2)=lbound(imp_merge,2)
          tub(1)=ubound(imp_merge,1)
          tub(2)=ubound(imp_merge,2)
          if(cpl_taux) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                taux(i,j,:) = imp_taux(i,j,1)
            enddo
            enddo
          endif
          if(cpl_tauy) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                tauy(i,j,:) = imp_tauy(i,j,1)
            enddo
            enddo
          endif
          if(cpl_wndspd.or.calc_wndspd) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                wndspd(i,j,:) = imp_wndspd(i,j,1)
            enddo
            enddo
          endif
          if(cpl_ustara) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                ustara(i,j,:) = imp_ustara(i,j,1)
            enddo
            enddo
          endif
          if(cpl_airtmp) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                airtmp(i,j,:) = imp_airtmp(i,j,1)
            enddo
            enddo
          endif
          if(cpl_vapmix) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                vapmix(i,j,:) = imp_vapmix(i,j,1)
            enddo
            enddo
          endif
          if(cpl_precip) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                precip(i,j,:) = imp_precip(i,j,1)
            enddo
            enddo
          endif
          if(cpl_surtmp) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                surtmp(i,j,:) = imp_surtmp(i,j,1)
            enddo
            enddo
          endif
          if(cpl_seatmp) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                seatmp(i,j,:) = imp_seatmp(i,j,1)
            enddo
            enddo
          endif
          if(cpl_swflx_net.or.cpl_swflx_net2down.or.cpl_swflxd) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                swflx(i,j,:) = imp_swflx(i,j,1)
            enddo
            enddo
          endif
          if(cpl_lwflx_net.or.cpl_lwflx_net2down.or.cpl_lwflxd) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                lwflx(i,j,:) = imp_lwdflx(i,j,1)
            enddo
            enddo
          endif
          if(cpl_u10) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                wndspx(i,j,:) = imp_wndspx(i,j,1)
            enddo
            enddo
          endif
          if(cpl_v10) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                wndspy(i,j,:) = imp_wndspy(i,j,1)
            enddo
            enddo
          endif
          if(cpl_mslprs) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                mslprs(i,j,:) = imp_mslprs(i,j,1)
            enddo
            enddo
          endif
          if(calc_radflx) then
            do j=tlb(2), tub(2)
            do i=tlb(1), tub(1)
              if (imp_merge(i,j).ne.0) &
                radflx(i,j,:) = imp_radflx(i,j,1)
            enddo
            enddo
          endif
        else
          write(lp,'(/ a,a /)') 'error - ', &
             rname//' merge natm.ne.2 is not implemented'
          call flush(lp)
          call xcstop(rname)
                 stop rname
        end if

!   check CICE feedback
      if ((.not.cpl_sic) .or.(.not.cpl_sitx).or.(.not.cpl_sity).or. &
        (.not.cpl_siqs).or.(.not.cpl_sifh).or.(.not.cpl_sifs).or. &
        (.not.cpl_sifw).or.(.not.cpl_sit) .or.(.not.cpl_sih) .or. &
        (.not.cpl_siu) .or.(.not.cpl_siv)) then
        if (mnproc.eq.1) print *, &
          'warning... no feedback from CICE to HYCOM('//rname//')'
        else
!       copy import data to sea ice
        do j=1, jja
        do i=1, ii
          if (imp_merge(i,j).ne.0) then
            if (ishlf(i,j).eq.1) then !standard ocean point
              if ((iceflg.ge.2).and.(icmflg.ne.3)) then
                covice(i,j)=sic_import(i,j) !Sea Ice Concentration
                si_c(i,j)=sic_import(i,j) !Sea Ice Concentration
                if (covice(i,j).gt.0.0) then
                   si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
                   si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
                  fswice(i,j) =  siqs_import(i,j) !Solar Heat Flux thru Ice to Ocean
                  flxice(i,j) =  fswice(i,j) + &
                                 sifh_import(i,j) !Ice Freezing/Melting Heat Flux
                  sflice(i,j) =  sifs_import(i,j)*1.e3 !Ice Salt Flux
                  wflice(i,j) =  sifw_import(i,j) !Ice freshwater Flux
                  temice(i,j) =   sit_import(i,j) !Sea Ice Temperature
                    si_t(i,j) =   sit_import(i,j) !Sea Ice Temperature
                  thkice(i,j) =   sih_import(i,j) !Sea Ice Thickness
                    si_h(i,j) =   sih_import(i,j) !Sea Ice Thickness
                    si_u(i,j) =   siu_import(i,j) !Sea Ice X-Velocity
                    si_v(i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
                else
                   si_tx(i,j) = 0.0
                   si_ty(i,j) = 0.0
                  fswice(i,j) = 0.0
                  flxice(i,j) = 0.0
                  sflice(i,j) = 0.0
                  wflice(i,j) = 0.0
                  temice(i,j) = 0.0
                    si_t(i,j) = 0.0
                  thkice(i,j) = 0.0
                    si_h(i,j) = 0.0
                    si_u(i,j) = 0.0
                    si_v(i,j) = 0.0
                endif !covice
              elseif ((iceflg.ge.2).and.(icmflg.eq.3)) then
                si_c(i,j)=sic_import(i,j) !Sea Ice Concentration
                if (si_c(i,j).gt.0.0) then
                  si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
                  si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
                   si_h(i,j) =   sih_import(i,j) !Sea Ice Thickness
                   si_t(i,j) =   sit_import(i,j) !Sea Ice Temperature
                   si_u(i,j) =   siu_import(i,j) !Sea Ice X-Velocity
                   si_v(i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
                else
                  si_tx(i,j) = 0.0
                  si_ty(i,j) = 0.0
                   si_h(i,j) = 0.0
                   si_t(i,j) = 0.0
                   si_u(i,j) = 0.0
                   si_v(i,j) = 0.0
                endif !covice
              endif !iceflg>=2 (icmflg)
            endif !ishlf
          endif !imp_merge
        enddo
        enddo

#if defined(ARCTIC)
!       update last active row of array
!jcx    call xctila( sic_import,1,1,halo_ps) !Sea Ice Concentration
!jcx    call xctila(sitx_import,1,1,halo_pv) !Sea Ice X-Stress
!jcx    call xctila(sity_import,1,1,halo_pv) !Sea Ice Y-Stress
!jcx    call xctila(siqs_import,1,1,halo_ps) !Solar Heat Flux thru Ice to Ocean
!jcx    call xctila(sifh_import,1,1,halo_ps) !Ice Freezing/Melting Heat Flux
!jcx    call xctila(sifs_import,1,1,halo_ps) !Ice Freezing/Melting Salt Flux
!jcx    call xctila(sifw_import,1,1,halo_ps) !Ice Net Water Flux
!jcx    call xctila( sit_import,1,1,halo_ps) !Sea Ice Temperature
!jcx    call xctila( sih_import,1,1,halo_ps) !Sea Ice Thickness
!jcx    call xctila( siu_import,1,1,halo_pv) !Sea Ice X-Velocity
!jcx    call xctila( siv_import,1,1,halo_pv) !Sea Ice Y-Velocity
        if ((iceflg.ge.2).and.(icmflg.ne.3)) then
          call xctila(covice,1,1,halo_ps) !Sea Ice Concentration
          call xctila(  si_c,1,1,halo_ps) !Sea Ice Concentration
          call xctila( si_tx,1,1,halo_pv) !Sea Ice X-Stress into ocean
          call xctila( si_ty,1,1,halo_pv) !Sea Ice Y-Stress into ocean
          call xctila(fswice,1,1,halo_ps) !Solar Heat Flux thru Ice to Ocean
          call xctila(flxice,1,1,halo_ps) !Ice Freezing/Melting Heat Flux
          call xctila(sflice,1,1,halo_ps) !Ice Salt Flux
          call xctila(wflice,1,1,halo_ps) !Ice Freshwater Flux
          call xctila(temice,1,1,halo_ps) !Sea Ice Temperature
          call xctila(  si_t,1,1,halo_ps) !Sea Ice Temperature
          call xctila(thkice,1,1,halo_ps) !Sea Ice Thickness
          call xctila(  si_h,1,1,halo_ps) !Sea Ice Thickness
          call xctila(  si_u,1,1,halo_pv) !Sea Ice X-Velocity
          call xctila(  si_v,1,1,halo_pv) !Sea Ice Y-Velocity
        elseif ((iceflg.ge.2).and.(icmflg.eq.3)) then
          call xctila(  si_c,1,1,halo_ps) !Sea Ice Concentration
          call xctila( si_tx,1,1,halo_pv) !Sea Ice X-Stress into ocean
          call xctila( si_ty,1,1,halo_pv) !Sea Ice Y-Stress into ocean
          call xctila(  si_h,1,1,halo_ps) !Sea Ice Thickness
          call xctila(  si_t,1,1,halo_ps) !Sea Ice Temperature
          call xctila(  si_u,1,1,halo_pv) !Sea Ice X-Velocity
          call xctila(  si_v,1,1,halo_pv) !Sea Ice Y-Velocity
        endif
#endif

!       smooth sea ice velocity fields
        call psmooth(si_u,0,0,ishlf,util1)
        call psmooth(si_v,0,0,ishlf,util1)
#if defined(ARCTIC)
        call xctila(si_u,1,1,halo_pv)
        call xctila(si_v,1,1,halo_pv)
#endif
!       call xctilr(si_u,1,1, nbdy,nbdy, halo_pv)
!       call xctilr(si_v,1,1, nbdy,nbdy, halo_pv)

!       copy back from si_ to _import for archive_ice
        do j=1, jja
        do i=1, ii
          if (imp_merge(i,j).ne.0) then
            if (si_c(i,j).gt.0.0) then
              siu_import(i,j)=si_u(i,j) !Sea Ice X-Velocity
              siv_import(i,j)=si_v(i,j) !Sea Ice Y-Velocity
            endif !si_c
          endif
        enddo !i
        enddo !j
      endif !feedback from CICE to HYCOM
#else
        write(lp,'(/ a,a /)') 'error - ', &
          rname//' merge requires USE_NUOPC_CESMBETA or ESPC_COUPLE'
        call flush(lp)
        call xcstop(rname)
               stop rname
#endif
      end subroutine hycom_imp_mrg

! --- return imp_latflx if imp_merge(i,j) is true
      real function hycom_imp_mrg_latflx(i,j,fval)
      integer :: i,j
      real :: fval
      character(*),parameter :: rname = 'hycom_imp_mrg_latflx'

#if defined (USE_NUOPC_GENERIC)
        if (cpl_latflx.and.(imp_merge(i,j).ne.0)) then
          hycom_imp_mrg_latflx=imp_latflx(i,j,1)
        else
          hycom_imp_mrg_latflx=fval
        endif
#else
        hycom_imp_mrg_latflx=fval
        write(lp,'(/ a,a /)') 'error - ', &
          rname//' merge requires USE_NUOPC_CESMBETA or ESPC_COUPLE'
        call flush(lp)
        call xcstop(rname)
               stop rname
#endif
      end function hycom_imp_mrg_latflx

! --- return imp_sensflx if imp_merge(i,j) is true
      real function hycom_imp_mrg_sensflx(i,j,fval)
      integer :: i,j
      real :: fval
      character(*),parameter :: rname = 'hycom_imp_mrg_sensflx'

#if defined (USE_NUOPC_GENERIC)
        if (cpl_sensflx.and.(imp_merge(i,j).ne.0)) then
          hycom_imp_mrg_sensflx=imp_sensflx(i,j,1)
        else
          hycom_imp_mrg_sensflx=fval
        endif
#else
        hycom_imp_mrg_sensflx=fval
        write(lp,'(/ a,a /)') 'error - ', &
          rname//' merge requires USE_NUOPC_CESMBETA or ESPC_COUPLE'
        call flush(lp)
        call xcstop(rname)
               stop rname
#endif
      end function hycom_imp_mrg_sensflx

      end module mod_import

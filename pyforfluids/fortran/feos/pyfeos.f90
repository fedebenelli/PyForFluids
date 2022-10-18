   subroutine ar(&
         n, v, t, &
         model, mixrule, &
         z, pc, tc, w, kij, lij, &
         ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2 &
         )
      use between, only: pr_ar, srk_ar
      implicit none
      integer, intent(in) :: n ! Number of components
      real(8), intent(in) :: v ! Volume
      real(8), intent(in) :: t ! Temperature
      character(len=10), intent(in) :: model
      character(len=10), intent(in) :: mixrule
      
      real(8), intent(in) :: z(n) ! Molar fractions

      real(8), intent(in) :: pc(n) ! Critical pressures
      real(8), intent(in) :: tc(n) ! Critical temperatures
      real(8), intent(in) :: w(n)  ! Acentric factors

      real(8), intent(in) :: kij(n, n) ! Kij matrix
      real(8), intent(in) :: lij(n, n) ! Lij matrix

      real(8), intent(out) :: ar_val   ! Residual Hellmholtz Energy
      real(8), intent(out) :: ar_dt    
      real(8), intent(out) :: ar_dv
      real(8), intent(out) :: ar_dt2
      real(8), intent(out) :: ar_dv2
      real(8), intent(out) :: ar_dtv
      real(8), intent(out) :: ar_dn(n)
      real(8), intent(out) :: ar_dn2(n, n)

      select case(model)
      case ("PR")
         call pr_ar(&
            n, v, t, &
            z, pc, tc, w, kij, lij, &
            ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2 &
         )
      case ("SRK")
         call srk_ar(&
            n, v, t, &
            z, pc, tc, w, kij, lij, &
            ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2 &
         )
      end select
   end subroutine ar

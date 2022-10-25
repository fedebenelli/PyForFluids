module between
   use feos, only: &
      scalar_property, &
      PengRobinson, PR,&
      SoaveRedlichKwong, SRK, &
      ClassicVdW, CubicFluid

contains

   subroutine pr_ar(&
         n, v, t, &
         z, pc, tc, w, kij, lij, &
         ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2, &
         ar_dtn, ar_dvn &
   )
      ! Peng Robinson EoS calculation for Helmholtz Energy
      integer, intent(in) :: n ! Number of components
      real(8), intent(in) :: v ! Volume
      real(8), intent(in) :: t ! Temperature
      
      real(8), intent(in) :: z(n) ! Molar fractions

      real(8), intent(in) :: pc(n) ! Critical pressures
      real(8), intent(in) :: tc(n) ! Critical temperatures
      real(8), intent(in) :: w(n)  ! Acentric factors

      real(8), intent(in) :: kij(n, n) ! Kij matrix
      real(8), intent(in) :: lij(n, n) ! Lij matrix

      real(8), intent(out) :: ar_val
      real(8), intent(out) :: ar_dt
      real(8), intent(out) :: ar_dv
      real(8), intent(out) :: ar_dt2
      real(8), intent(out) :: ar_dv2
      real(8), intent(out) :: ar_dtv
      real(8), intent(out) :: ar_dn(n)
      real(8), intent(out) :: ar_dn2(n, n)
      real(8), intent(out) :: ar_dtn(n)
      real(8), intent(out) :: ar_dvn(n)

      type(ClassicVdW)   :: mixing_rule
      type(PengRobinson) :: components(n)
      type(CubicFluid)   :: fluid
      type(scalar_property) :: ar_in

      integer :: i


      do i=1, n
         components(i) = PR(name="", pc=pc(i), tc=tc(i), w=w(i))
         components(i)%moles = z(i)
      end do

      mixing_rule = ClassicVdW(kij=kij, lij=lij)
      !fluid = CubicFluid(components=components, mixing_rule=mixing_rule)

      fluid%components = components
      fluid%mixing_rule = mixing_rule

      ar_in = fluid%residual_helmholtz(v, t)

      call scalar_to_arrays(&
         ar_in, &
         ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2, &
         ar_dtn, ar_dvn &
         )
      deallocate(mixing_rule%kij, mixing_rule%lij)
      deallocate(fluid%components, fluid%mixing_rule)

   end subroutine pr_ar
   
   subroutine srk_ar(&
         n, v, t, &
         z, pc, tc, w, kij, lij, &
         ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2,&
         ar_dtn, ar_dvn &
   )
      ! Peng Robinson EoS calculation for Helmholtz Energy
      integer, intent(in) :: n ! Number of components
      real(8), intent(in) :: v ! Volume
      real(8), intent(in) :: t ! Temperature
      
      real(8), intent(in) :: z(n) ! Molar fractions

      real(8), intent(in) :: pc(n) ! Critical pressures
      real(8), intent(in) :: tc(n) ! Critical temperatures
      real(8), intent(in) :: w(n)  ! Acentric factors

      real(8), intent(in) :: kij(n, n) ! Kij matrix
      real(8), intent(in) :: lij(n, n) ! Lij matrix

      real(8), intent(out) :: ar_val
      real(8), intent(out) :: ar_dt
      real(8), intent(out) :: ar_dv
      real(8), intent(out) :: ar_dt2
      real(8), intent(out) :: ar_dv2
      real(8), intent(out) :: ar_dtv
      real(8), intent(out) :: ar_dn(n)
      real(8), intent(out) :: ar_dn2(n, n)
      real(8), intent(out) :: ar_dtn(n)
      real(8), intent(out) :: ar_dvn(n)

      type(ClassicVdW)   :: mixing_rule
      type(SoaveRedlichKwong) :: components(n)
      type(CubicFluid)   :: fluid
      type(scalar_property) :: ar_in

      integer :: i


      do i=1, n
         components(i) = SRK(name="", pc=pc(i), tc=tc(i), w=w(i))
         components(i)%moles = z(i)
      end do

      mixing_rule = ClassicVdW(kij=kij, lij=lij)

      fluid%components = components
      fluid%mixing_rule = mixing_rule

      ar_in = fluid%residual_helmholtz(v, t)

      call scalar_to_arrays(&
         ar_in, &
         ar_val, ar_dt, ar_dv, ar_dt2, ar_dv2, ar_dtv, ar_dn, ar_dn2, &
         ar_dtn, ar_dvn &
         )
      deallocate(mixing_rule%kij, mixing_rule%lij)
      deallocate(fluid%components, fluid%mixing_rule)

   end subroutine srk_ar


   subroutine scalar_to_arrays(sp, val, dt, dv, dt2, dv2, dtv, dn, dn2, dtn, dvn)
      type(scalar_property), intent(in out) :: sp
      real(8), intent(out) :: val
      real(8), intent(out) :: dt
      real(8), intent(out) :: dv
      real(8), intent(out) :: dt2
      real(8), intent(out) :: dv2
      real(8), intent(out) :: dtv
      real(8), intent(out) :: dn(size(sp%dn))
      real(8), intent(out) :: dn2(size(sp%dn), size(sp%dn))
      real(8), intent(out) :: dtn(size(sp%dn))
      real(8), intent(out) :: dvn(size(sp%dn))

      val = sp%val
      dt  = sp%dt
      dv  = sp%dv
      dt2 = sp%dt2
      dv2 = sp%dv2
      dtv = sp%dtv
      dtn = sp%dtn
      dvn = sp%dvn
      dn  = sp%dn
      dn2 = sp%dn2
      deallocate(sp%dtn)
      deallocate(sp%dvn)
      deallocate(sp%dn)
      deallocate(sp%dn2)
   end subroutine scalar_to_arrays
end module between

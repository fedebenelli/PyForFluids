! ------------------------
! Thermodynamic Properties
! ------------------------
Module thermo_props  !
   Implicit None

Contains
        Subroutine mean_molecular_weight(X, M, MM)
                real(8), intent(in) :: X(21), M(21)
                real(8), intent(out) :: MM
                integer :: i

                MM = 0
                do i = 1, size(X)
                MM = MM + X(i) * M(i)
                end do
                MM = MM / 1000 ! Translate to Kilograms
        End Subroutine mean_molecular_weight

   Subroutine zeta(delta, Ar, z)
      real(8), intent(in) :: delta, Ar(3, 3)
      real(8), intent(out) :: z
      z = 1.d0 + delta * Ar(2, 1)
   End Subroutine zeta

   Subroutine isochoric_heat(tau, R, Ao, Ar, cv)
      real(8), intent(in) :: tau, R, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: cv

      cv = - tau ** 2.d0 * (Ao(3, 2) + Ar(3, 2))
      cv = cv * R
   End Subroutine isochoric_heat

   Subroutine isobaric_heat(delta, tau, R, Ao, Ar, cp)
      real(8), intent(in) :: delta, tau, R, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: cp
      real(8) :: up, down

      up = (1.d0 + delta * Ar(2, 1) - delta * tau * Ar(3, 3)) ** 2.d0
      down = 1.d0 + 2.d0 * delta * Ar(2, 1) + delta ** 2.d0 * Ar(3, 1)

      cp = - tau ** 2.d0 * (Ao(3, 2) + Ar(3, 2)) + up / down
      cp = cp * R
   End Subroutine isobaric_heat

   Subroutine sound_speed(delta, tau, R, T, M, Ao, Ar, w)
      real(8), intent(in) :: delta, tau, R, T, M, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: w
      real(8) :: up, down

      up = (1.d0 + delta * Ar(2, 1) - delta * tau * Ar(3, 3)) ** 2
      down = (tau ** 2 * (Ao(3, 2) + Ar(3, 2)))

      w = 1.d0 + 2.d0 * delta * Ar(2, 1) + delta ** 2 * Ar(3, 1) - up / down
      w = sqrt(w * R * T / M)
   End Subroutine sound_speed

   Subroutine isothermal_thermal_coefficent(delta, tau, rho, Ar, delta_t)
      real(8), intent(in) :: delta, Ar(3, 3), tau, rho
      real(8), intent(out) :: delta_t
      real(8) :: up, down

      up = 1.d0 + delta * Ar(2, 1) - delta * tau * Ar(3, 3)
      down = 1.d0 + 2.d0 * delta * Ar(2, 1) + delta ** 2 * Ar(3, 1)

      delta_t = (1.d0 - up / down) / rho
   End Subroutine isothermal_thermal_coefficent

   Subroutine dp_dt(rho, delta, tau, R, Ar, dpdt)
      real(8), intent(in) :: delta, R, rho, tau
      real(8), intent(out) :: dpdt
      real(8), dimension(3, 3), intent(in) :: Ar

      dpdt = rho * R * (1.d0 + delta * Ar(2, 1) - delta * tau * Ar(3, 3))

   End Subroutine dp_dt

   Subroutine dp_drho(T, delta, R, Ar, dpdrho)
      real(8), intent(in) :: delta, R, T
      real(8), intent(out) :: dpdrho
      real(8), dimension(3, 3), intent(in) :: Ar

      dpdrho = R * T * (1.d0 + 2.d0 * delta * Ar(2, 1) + (delta ** 2) * Ar(3, 1))

   End Subroutine dp_drho

   Subroutine dp_dv(rho, delta, T, R, Ar, dpdv) ! !Recordar que acá hablamos de volumen molar!!
      real(8), intent(in) :: delta, R, rho, T
      real(8), intent(out) :: dpdv
      real(8), dimension(3, 3), intent(in) :: Ar

      dpdv = - (rho ** 2) * R * T * (1.d0 + 2.d0 * delta * Ar(2, 1) + (delta ** 2) * Ar(3, 1))

   End Subroutine dp_dv

   Subroutine pressure(delta, rho, R, T, Ar, p)
      real(8), intent(in) :: delta, rho, R, T, Ar(3, 3)
      real(8), intent(out) :: p
      real(8) :: z

      call zeta(delta, Ar, z)
      p = z * (rho * 1000.d0) * R * T
   End Subroutine pressure

   Subroutine entropy(tau, R, Ao, Ar, s)
      real(8), intent(in) :: tau, R, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: s

      s = (tau * (Ao(2, 2) + Ar(2, 2)) - Ao(1, 1) - Ar(1, 1)) * R
   End Subroutine entropy

   Subroutine internal_energy(tau, R, T, Ao, Ar, u)
      real(8), intent(in) :: tau, R, T, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: u

      u = tau * (Ao(2, 2) + Ar(2, 2)) * R * T
   End Subroutine internal_energy

   Subroutine enthalpy(delta, tau, R, T, Ao, Ar, h)
      real(8), intent(in) :: delta, tau, R, T, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: h

      h = ((1.d0 + tau * (Ao(2, 2) + Ar(2, 2)) + delta * Ar(2, 1))) * R * T
   End Subroutine enthalpy

   Subroutine gibbs_free_energy(delta, R, T, Ao, Ar, g)
      real(8), intent(in) :: delta, R, T, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: g

      g = (1.d0 + Ao(1, 1) + Ar(1, 1) + delta * Ar(2, 1)) * R * T
   End Subroutine gibbs_free_energy

   Subroutine joule_thomson_coeff(delta, tau, rho, R, Ao, Ar, JT)
      real(8), intent(in) :: delta, tau, rho, R, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: JT
      real(8) :: up, down_1, down_2

      up = - (delta * Ar(2, 1) + delta ** 2 * Ar(3, 1) + delta * tau * Ar(3, 3))
      down_1 = (1.d0 + delta * Ar(2, 1) - delta * tau * Ar (3, 3)) ** 2
      down_2 = - tau ** 2 * (Ao(3, 2) + Ar(3, 2)) * (1.d0 + 2.d0 * delta * Ar(2, 1) + delta ** 2 * Ar(3, 1))

      JT = up / (down_1 + down_2) / (R * rho * 1000)

   End Subroutine joule_thomson_coeff

   Subroutine isentropic_exponent(delta, tau, Ao, Ar, k)
      real(8), intent(in) :: delta, tau, Ao(3, 3), Ar(3, 3)
      real(8), intent(out) :: k
      real(8) :: up1, down1, up2, down2

      up1 = 1.d0 + 2 * delta * Ar(2, 1) + delta ** 2 * Ar(3, 1)
      down1 = 1.d0 + delta * Ar(2, 1)
      up2 = (1.d0 + delta * Ar(2, 1) - delta * tau * Ar(2, 3)) ** 2
      down2 = tau ** 2 * (Ao(3, 2) + Ar(2, 2)) * (1.d0 + 2.d0 * delta * Ar(2, 1) + delta ** 2 * Ar(3, 2))

      k = (up1 / down1) * (1.d0 - up2 / down2)
   End Subroutine isentropic_exponent

   Subroutine second_thermal_virial_coeff(rho_r, Ar, B)
      real(8), intent(in) :: rho_r, Ar(3, 3)
      real(8), intent(out) :: B
      real(8) :: delta

      delta = 1d-15
      B = Ar(2, 1) / rho_r
   End Subroutine second_thermal_virial_coeff

   Subroutine third_thermal_virial_coeff(rho_r, Ar, C)
      real(8), intent(in) :: rho_r, Ar(3, 3)
      real(8), intent(out) :: C
      real(8) :: delta

      delta = 1d-15
      C = Ar(3, 1) / (rho_r ** 2)
   End Subroutine third_thermal_virial_coeff

   Subroutine helmholtz_per_mol(x, delta, tau, rho_r, T_r, ar, ar_x, ar_dx, &
                                v_r_x, T_r_x, dar_dn, dadr_dn)
           real(8), intent(in) :: x(21), delta, tau, rho_r, T_r, &
                   ar(3, 3), ar_x(21), ar_dx(21), v_r_x(21), T_r_x(21)
           real(8), dimension(21), intent(out) :: dar_dn, dadr_dn
           real(8), dimension(21) :: drhor_dn, dtr_dn
           integer :: i

           dar_dn = 0.d0
           dadr_dn = 0.d0

           drhor_dn = - rho_r ** 2 * v_r_x
           dtr_dn = v_r_x

           do i = 1, size(dtr_dn)
                drhor_dn = drhor_dn - (- x(i) * (- rho_r ** 2 * v_r_x(i)))
                dtr_dn = dtr_dn - x(i) * t_r_x(i)
           end do

           dar_dn = delta * ar(2, 1) * (1.d0 - 1.d0 / rho_r * dtr_dn)
           dar_dn = dar_dn + tau * ar(2, 2) / T_r * dtr_dn + ar_x

           dadr_dn = delta * ar(3, 1) * (1.d0 - 1.d0 / rho_r * dtr_dn)
           dadr_dn = dadr_dn + tau * ar(3, 3) / T_r * dtr_dn + ar_dx

           do i = 1, size(dtr_dn)
                   dar_dn = dar_dn - X(i) * ar_x(i)
                   dadr_dn = dadr_dn - X(i) * ar_dx(i)
           end do
   End Subroutine helmholtz_per_mol

End Module thermo_props

! Author= Federico Benelli
! Gerg-2008 Equation
! Started at: 01/07/2021
! Last Modified: mié 01 sep 2021 15:48:41
!

! -------------------
! GERG 2008 Functions
! -------------------
Subroutine reducing_funcs(X, rho_r, T_r, T_r_x, V_r_x)
   ! REDUCING DENSITY AND TEMPERATURE
   ! input:
   ! - X      (array): molar fractions
   ! - Bv     (array): beta parameters for reducing density
   ! - Gv     (array): gamma parameters for reducing density
   ! - Bt     (array): beta parameters for reducing temperature
   ! - Gt     (array): gamma parameters for reducing temperature
   ! - rho_c  (array): citical densities
   ! - T_c    (array): citical temperatures
   ! output:
   ! - rho_r  (float): Reducing density
   ! - T_r    (float): reducing temperature (Tr)
   use parameters
   real(8), dimension(21), intent(in) :: X
   real(8), intent(out) :: rho_r, T_r, T_r_x(21), V_r_x(21)
   integer :: i, j, k
   ! Internal variables
   real(8) :: vki, tki, c_v, c_t, x_sum, x_pro, xki, Bv_xki, Bt_xki, &
              fv_ki, ft_ki, vik, tik, Bv_xik, Bt_xik, fv_ik, ft_ik
   call get_params()

   rho_r = 0.d0
   T_r = 0.d0
   t_r_x = 0.d0
   v_r_x = 0.d0

   do i = 1, N
   if (X(i) > eps) then
      rho_r = rho_r + X(i) ** 2 / rho_c(i)
      T_r = T_r + X(i) ** 2 * T_c(i)
   end if
   end do

   do i = 1, N - 1
   if (X(i) > eps) then
   do j = i + 1, N
   if (X(j) > eps) then
      rho_r = rho_r + &
              2.d0 * X(i) * X(j) * Bv(i, j) * Gv(i, j) &
              * (X(i) + X(j)) / (Bv(i, j) ** 2 * X(i) + X(j)) &
              * 1.d0 / 8.d0 * (rho_c(i) ** (- 1.d0 / 3.d0) &
                          + rho_c(j) ** (- 1.d0 / 3.d0)) ** 3

      T_r = T_r + &
            2.d0 * X(i) * X(j) * Bt(i, j) * Gt(i, j) &
            * (X(i) + X(j)) / (Bt(i, j) ** 2 * X(i) + X(j)) &
            * sqrt((T_c(i) * T_c(j)))
   end if
   end do
   end if
   end do
   ! Actually, the reducing volume was calculated, let's make it into density
   rho_r = 1 / rho_r

   ! Compositional derivatives
   do i = 1, N
   if (x(i) > eps) then
      T_r_x(i) = 2.d0 * X(i) * T_c(i)
      V_r_x(i) = 2.d0 * X(i) / (rho_c(i))

      do k = 1, i - 1
         vki = 1.d0 / 8.d0 * (rho_c(k) ** (- 1.d0 / 3.d0) &
                          + rho_c(i) ** (- 1.d0 / 3.d0)) ** 3
         tki = sqrt((T_c(k) * T_c(i)))
         c_v = 2 * Bv(k, i) * Gv(k, i) * vki
         c_t = 2 * Bt(k, i) * Gt(k, i) * tki

         x_sum = X(k) + X(i)
         x_pro = X(k) * X(i)
         xki = X(k) + X(i)

         Bv_xki = Bv(k, i) ** 2 * X(k) + X(i)
         Bt_xki = Bt(k, i) ** 2 * X(k) + X(i)

         fv_ki = X(k) * x_sum / Bv_xki + x_pro / Bv_xki &
                 * (1.d0 - x_sum / Bv_xki)
         ft_ki = X(k) * x_sum / Bt_xki + x_pro / Bt_xki &
                 * (1.d0 - x_sum / Bt_xki)

         V_r_x(i) = T_r_x(i) + c_v * fv_ki
         T_r_x(i) = V_r_x(i) + c_t * ft_ki

      end do
      do k = i + 1, N
         vik = 1.d0 / 8.d0 * (rho_c(i) ** (- 1.d0 / 3.d0) &
                          + rho_c(k) ** (- 1.d0 / 3.d0)) ** 3
         tik = sqrt((T_c(i) * T_c(k)))
         c_v = 2 * Bv(i, k) * Gv(i, k) * vik
         c_t = 2 * Bt(i, k) * Gt(i, k) * tik

         x_sum = X(i) + X(k)
         x_pro = X(i) * X(k)
         xki = X(i) + X(k)

         Bv_xik = Bv(i, k) ** 2 * X(i) + X(k)
         Bt_xik = Bt(i, k) ** 2 * X(i) + X(k)

         fv_ik = X(i) * x_sum / Bv_xik &
                 + x_pro / Bv_xik * (1.d0 - x_sum / Bv_xik)
         ft_ik = X(i) * x_sum / Bt_xik &
                 + x_pro / Bt_xik * (1.d0 - x_sum / Bt_xik)

         V_r_x(i) = T_r_x(i) + c_v * fv_ik
         T_r_x(i) = V_r_x(i) + c_t * ft_ik
      end do
   end if
   end do

End Subroutine reducing_funcs

! Pure Compound Helmholtz Energy (and derivatives) Calculations
! -----------------------------------------------
Subroutine a_oio(rho, T, rho_c, T_c, n, v, aoio)
   ! IDEAL GAS ENERGY
   ! ----------------
   ! Calculate the pure compound ideal Helmholtz Enegry and its derivatives
   !
   ! input:
   ! - rho (float): Density
   ! - T (float): Temperature
   ! - rho_c (float): Critical Density
   ! - T_c (float): Critical Temperature
   ! - n (dimension): n parameters
   ! - v (dimension): v parameters
   ! output:
   ! - aoio: Ideal Gas Helmholtz Energy and its derivatives,
   ! stoikuctured like:
   ! --------------+------------------+------------------|
   ! aoio          |     0            |    0             |
   ! d(aoio)/dd    |     d(aoio)/dt   |    0             |
   ! d2(aoio)/dd2  |     d2(aoio)/dt2 |    d2(aoio)/dddt |
   ! ----------------------------------------------------|
   Implicit None
   real(8), intent(in) :: rho, T, rho_c, T_c
   real(8), intent(in), dimension(7) :: n, v
   real(8), intent(out) :: aoio(3, 3)
   real(8) :: r, Tr, Dr, eps = 1d-10
   integer :: k

   r = 8.314510d0 / 8.314472d0

   aoio = 0.d0
   Tr = T_c / T
   Dr = rho / rho_c

   do k = 4, 7
      if (v(k) > eps) then
         if (k == 4 .or. k == 6) then
            aoio(1, 1) = aoio(1, 1) + r * (n(k) * log(abs(dsinh(v(k) * Tr))))
            aoio(2, 2) = aoio(2, 2) + r * (n(k) * v(k) / dtanh(v(k) * Tr))
            aoio(3, 2) = aoio(3, 2) - n(k) * (v(k) / dsinh(v(k) * Tr)) ** 2
         else
            aoio(1, 1) = aoio(1, 1) - r * (n(k) * log(dcosh(v(k) * Tr)))
            aoio(2, 2) = aoio(2, 2) - r * (n(k) * v(k) * dtanh(v(k) * Tr))
            aoio(3, 2) = aoio(3, 2) - n(k) * (v(k) / dcosh(v(k) * Tr)) ** 2
         end if
      end if
   end do

   aoio(1, 1) = aoio(1, 1) + log(Dr) + r * (n(1) + n(2) * Tr + n(3) * log(Tr))
   aoio(2, 1) = 1.d0 / Dr
   aoio(3, 1) = - 1.d0 / Dr ** 2

   aoio(2, 2) = aoio(2, 2) + r * (n(2) + n(3) / Tr)

   aoio(3, 2) = aoio(3, 2) - (n(3) * (T / T_c) ** 2)
   aoio(3, 2) = r * aoio(3, 2) 

   aoio(3, 3) = 0.d0
End Subroutine a_oio

Subroutine a_oir(delta, tau, Kpol, Kexp, n, d, t, c, aoir)
   ! RESIDUAL ENERGY
   ! ---------------
   ! Caculate the pure compound Residual Helmholtz Energy and its derivatives
   !
   ! input:
   ! - delta  (float)         : reduced density
   ! - tau    (float)         : reduced temperature
   ! - Kpol   (integer)       : number of only polynomial parameters
   ! - Kexp   (integer)       : number of only exponential parameters
   ! - n      (dimension)     : parameters n
   ! - d      (dimension)     : parameters d
   ! - t      (dimension)     : parameters t
   ! - c      (dimension)     : parameters c
   ! output:
   ! - aoir: Ideal Gas Helmholtz Energy and its derivatives,
   ! structured like:
   ! ---------------+-------------------+-----------------|
   ! aoir           |      0            |   0             |
   ! d(aoir)/dd     |      d(aoir)/dt   |   0             |
   ! d2(aoir)/dd2   |      d2(aoir)/dt2 |   d2(aoir)/dddt |
   ! -----------------------------------------------------|
   !

   real(8), intent(in) :: delta, tau
   integer, intent(in) :: Kpol, Kexp
   real(8), dimension(24), intent(in) :: n, t
   integer, dimension(24), intent(in) :: d, c
   real(8), dimension(3, 3), intent(out) :: aoir
   integer :: k

   aoir = 0.d0
   ! Residual Helmholtz Energy
   do k = 1, Kpol
      ! !write (0, *) k, n(k), d(k), t(k), c(k)
      aoir(1, 1) = aoir(1, 1) + &
                   n(k) * delta ** d(k) &
                   * tau ** t(k)
      aoir(2, 1) = aoir(2, 1) + n(k) * d(k) * delta ** (d(k) - 1) * tau ** (t(k))
      aoir(2, 2) = aoir(2, 2) + &
                   n(k) * t(k) * delta ** d(k) &
                   * tau ** (t(k) - 1.d0)
      aoir(3, 1) = aoir(3, 1) + &
                   n(k) * d(k) * (d(k) - 1.d0) * delta ** (d(k) - 2.d0) * tau ** t(k)
      aoir(3, 2) = aoir(3, 2) + &
                   n(k) * t(k) * (t(k) - 1.d0) &
                   * delta ** d(k) * tau ** (t(k) - 2.d0)
      aoir(3, 3) = aoir(3, 3) + &
                   n(k) * d(k) * t(k) &
                   * delta ** (d(k) - 1) * tau ** (t(k) - 1.d0)
   end do
   do k = Kpol + 1, Kpol + Kexp
      aoir(1, 1) = aoir(1, 1) + &
                   n(k) * delta ** d(k) &
                   * tau ** t(k) &
                   * exp(- delta ** c(k))
      ! First Derivative with reduced density
      aoir(2, 1) = aoir(2, 1) + &
                   n(k) * delta ** (d(k) - 1) &
                   * (d(k) - c(k) * delta ** c(k)) &
                   * tau ** t(k) * exp(- delta ** c(k))
      ! First Derivative with reduced temperature
      aoir(2, 2) = aoir(2, 2) + &
                   n(k) * t(k) * delta ** d(k) &
                   * tau ** (t(k) - 1) &
                   * exp(- delta ** c(k))
      ! Second Derivative with reduced density
      aoir(3, 1) = aoir(3, 1) + &
                   n(k) * delta ** (d(k) - 2.d0) * ( &
                   (d(k) - c(k) * delta ** c(k)) &
                   * (d(k) - 1.d0 - c(k) * delta ** c(k)) &
                   - c(k) ** 2.d0 * delta ** c(k)) &
                   * tau ** t(k) * exp(- delta ** c(k))
      ! Second Derivative with reduced temperature
      aoir(3, 2) = aoir(3, 2) + &
                   n(k) * t(k) * (t(k) - 1.d0) &
                   * delta ** d(k) * tau ** (t(k) - 2.d0) &
                   * exp(- delta ** c(k))
      ! Second Derivative with reduced temperature and density
      aoir(3, 3) = aoir(3, 3) + &
                   n(k) * t(k) * delta ** (d(k) - 1) &
                   * (d(k) - c(k) * delta ** c(k)) &
                   * tau ** (t(k) - 1.d0) * exp(- delta ** c(k))
   end do

End Subroutine a_oir

Subroutine a_ijr(delta, tau, Kpolij, Kexpij, &
                 n, d, t, eta, eps, gamm, beta, aijr)
   ! DEPARTURE FUNCTION ENERGY
   ! ---------------
   ! Caculate binary departure Helmholtz Energy and its derivatives
   !
   ! input:
   ! - delta  (float)         : reduced density
   ! - tau    (float)         : reduced temperature
   ! - Kpolij (integer)       : number of only polynomial parameters
   ! - Kexpij (integer)       : number of only exponential parameters
   ! - n      (dimension)     : parameters n
   ! - d      (dimension)     : parameters d
   ! - t      (dimension)     : parameters t
   ! - eta    (dimension)     : parameters eta
   ! - eps    (dimension)     : parameters epsilon
   ! - gamm   (dimension)     : parameters gamma
   ! - beta   (dimension)     : parameters beta
   ! output:
   ! - aijr: Binary departure Helmholtz Energy and its derivatives,
   ! stoikuctured like:
   ! ---------------+------------------+------------------|
   ! aijr           |      0            |   0             |
   ! d(aijr)/dd     |      d(aijr)/dt   |   0             |
   ! d2(aijr)/dd2   |      d2(aijr)/dt2 |   d2(aijr)/dddt |
   ! -----------------------------------------------------|
   Implicit None
   real(8), intent(in) :: delta, tau
   integer, intent(in) :: Kpolij, Kexpij
   real(8), dimension(24), intent(in) :: n, t, eta, eps, gamm, beta
   integer, dimension(24), intent(in) :: d
   real(8), dimension(3, 3), intent(out) :: aijr
   integer :: k

   aijr = 0.d0
   do k = 1, Kpolij
      aijr(1, 1) = aijr(1, 1) + &
                   n(k) * delta ** d(k) * tau ** t(k)
      aijr(2, 1) = aijr(2, 1) + &
                   n(k) * d(k) * delta ** (d(k) - 1.d0) * tau ** t(k)
      aijr(3, 1) = aijr(3, 1) + &
                   n(k) * d(k) * (d(k) - 1.d0) * delta ** (d(k) - 2.d0) * tau ** t(k)
      aijr(2, 2) = aijr(2, 2) + &
                   n(k) * t(k) * delta ** d(k) * tau ** (t(k) - 1.d0)
      aijr(3, 2) = aijr(3, 2) + &
                   n(k) * t(k) * (t(k) - 1.d0) * delta ** d(k) * tau ** (t(k) - 2.d0)
      aijr(3, 3) = aijr(3, 3) + &
                   n(k) * d(k) * t(k) &
                   * delta ** (d(k) - 1.d0) * tau ** (t(k) - 1.d0)
   end do
   do k = Kpolij + 1, Kpolij + Kexpij
      aijr(1, 1) = aijr(1, 1) + &
                   n(k) * delta ** d(k) * tau ** t(k) &
                   * exp( &
                   - eta(k) * (delta - eps(k)) ** 2.d0 &
                   - beta(k) * (delta - gamm(k)) &
                   )
      aijr(2, 1) = aijr(2, 1) + &
                   n(k) * delta ** d(k) * tau ** t(k) &
                   * exp( &
                   - eta(k) * (delta - eps(k)) ** 2.d0 &
                   - beta(k) * (delta - gamm(k)) &
                   ) &
                   * (d(k) / delta - 2.d0 * eta(k) * (delta - eps(k)) - beta(k))
      aijr(3, 1) = aijr(3, 1) + &
                   n(k) * delta ** d(k) * tau ** t(k) * exp( &
                   - eta(k) * (delta - eps(k)) ** 2 - beta(k) * (delta - gamm(k))) &
                   * ((d(k) / delta - 2.d0 * eta(k) &
                      * (delta - eps(k)) - beta(k)) ** 2.d0 - d(k) / delta ** 2.d0 - 2.d0 * eta(k))
      aijr(2, 2) = aijr(2, 2) + &
                   n(k) * t(k) * delta ** d(k) * tau ** (t(k) - 1.d0) &
                   * exp( &
                   - eta(k) * (delta - eps(k)) ** 2.d0 &
                   - beta(k) * (delta - gamm(k)) &
                   )
      aijr(3, 2) = aijr(3, 2) + &
                   n(k) * t(k) * (t(k) - 1.d0) * delta ** d(k) * tau ** (t(k) - 2.d0) &
                   * exp( &
                   - eta(k) * (delta - eps(k)) ** 2.d0 - beta(k) * (delta - gamm(k)) &
                   )
      aijr(3, 3) = aijr(3, 3) + &
                   n(k) * t(k) * delta ** d(k) * tau ** (t(k) - 1) &
                   * exp( &
                   - eta(k) * (delta - eps(k)) ** 2.d0 &
                   - beta(k) * (delta - gamm(k)) &
                   ) &
                   * (d(k) / delta - 2.d0 * eta(k) * (delta - eps(k)) - beta(k))
   end do

End Subroutine a_ijr

Subroutine ideal_term(X, rho, T, rho_r, T_r, ao)
   use parameters
   real(8), intent(in) :: X(21), rho, T, rho_r, T_r
   real(8), intent(out) :: ao(3, 3)
   real(8) ::  aoio(3, 3), xi

   call get_params()

   ao = 0.d0
   do i = 1, N
      if (X(i) > eps) then
         aoio = 0
         call a_oio(rho, T, rho_c(i), T_c(i), n0i(i, :), th0i(i, :), aoio)

         xi = X(i)
         ao(1, 1) = ao(1, 1) + xi * (aoio(1, 1) + log(xi))

         ao(2, 1) = ao(2, 1) + xi * rho_r / rho_c(i) * aoio(2, 1)
         ao(3, 1) = ao(3, 1) + xi * (rho_r / rho_c(i)) ** 2 * aoio(3, 1)

         ao(2, 2) = ao(2, 2) + xi * (T_c(i) / T_r) * aoio(2, 2)
         ao(3, 2) = ao(3, 2) + xi * (T_c(i) / T_r) ** 2 * aoio(3, 2)

         ao(3, 3) = 0
      end if
   end do
End Subroutine ideal_term

Subroutine residual_term(X, delta, tau, ar, ar_x, ar_dx, ar_tx, ar_xx)
   use parameters
   Implicit None
   real(8), intent(in) :: delta, tau, X(21)
   real(8), dimension(3, 3), intent(out) :: ar
   real(8), dimension(21), intent(out) :: ar_x, ar_dx, ar_tx, ar_xx
   real(8), dimension(3, 3) :: aoir, aijr
   real(8) :: F
   integer :: i, j, k, temp_i, temp_k

   call get_params()

   ar = 0
   ar_x = 0
   ar_dx = 0
   ar_tx = 0
   ar_xx = 0

   do i = 1, N
      if (X(i) > eps) then
         aoir = 0.0

         call a_oir(delta, tau, Kpol(i), Kexp(i), &
                    noik(i, :), doik(i, :), toik(i, :), coik(i, :), aoir)

         ar(1, 1) = ar(1, 1) + X(i) * aoir(1, 1)
         ar(2, 1) = ar(2, 1) + X(i) * aoir(2, 1)
         ar(2, 2) = ar(2, 2) + X(i) * aoir(2, 2)
         ar(3, 1) = ar(3, 1) + X(i) * aoir(3, 1)
         ar(3, 2) = ar(3, 2) + X(i) * aoir(3, 2)
         ar(3, 3) = ar(3, 3) + X(i) * aoir(3, 3)

         ar_x(i) = ar(1, 1)
         ar_dx(i) = ar(2, 1)
         ar_tx(i) = ar(2, 2)
         ar_xx(i) = 0

         ! Compositional derivatives
         do k = 1, N
            temp_i = i
            temp_k = k
            F = Fij(temp_i, temp_k)
            if (F < eps .and. Fij(k, i) > eps) then
               ! Invert variables for cases like i=2, k=1
               ! TODO: NOT SURE IF THIS SHOULD BE DONE!
               temp_i = k
               temp_k = i
               F = Fij(temp_i, temp_k)
            end if
            if (k /= i .and. X(k) > eps .and. F > eps) then
               aijr = 0.d0
               call a_ijr(delta, tau, Kpolij(temp_i, temp_k), Kexpij(i, temp_k), &
                          nij(temp_i, temp_k, :), dij(temp_i, temp_k, :), &
                          tij(temp_i, temp_k, :), ethaij(temp_i, temp_k, :), &
                          epsij(temp_i, temp_k, :), gammaij(temp_i, temp_k, :), &
                          betaij(temp_i, temp_k, :), aijr)

               ar_x(i) = ar_x(i) + X(k) * Fij(i, k) * aijr(1, 1)
               ar_dx(i) = ar_dx(i) + X(k) * Fij(i, k) * aijr(2, 1)
               ar_tx(i) = ar_tx(i) + X(k) * Fij(i, k) * aijr(2, 2)
            end if
         end do
      end if
   end do

   do i = 1, N - 1
      do j = i + 1, N
      if (Fij(i, j) > eps .and. X(i) > eps .and. X(j) > eps) then
         call a_ijr(delta, tau, Kpolij(i, j), Kexpij(i, j), &
                    nij(i, j, :), dij(i, j, :), tij(i, j, :), ethaij(i, j, :), &
                    epsij(i, j, :), gammaij(i, j, :), betaij(i, j, :), &
                    aijr)
         ar(1, 1) = ar(1, 1) + X(i) * X(j) * Fij(i, j) * aijr(1, 1)
         ar(2, 1) = ar(2, 1) + X(i) * X(j) * Fij(i, j) * aijr(2, 1)
         ar(2, 2) = ar(2, 2) + X(i) * X(j) * Fij(i, j) * aijr(2, 2)
         ar(3, 1) = ar(3, 1) + X(i) * X(j) * Fij(i, j) * aijr(3, 1)
         ar(3, 2) = ar(3, 2) + X(i) * X(j) * Fij(i, j) * aijr(3, 2)
         ar(3, 3) = ar(3, 3) + X(i) * X(j) * Fij(i, j) * aijr(3, 3)
      end if
      end do
   end do

End Subroutine residual_term

Program main
        print *, "This only exist to compile with gfortran"
End Program main

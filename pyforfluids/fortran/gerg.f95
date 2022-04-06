! -----------------------------------------------------------------------------
! Author= Federico Benelli
! Gerg-2008 Equation
! Started at: 01/07/2021
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Mixing rules of critical properties
! -----------------------------------------------------------------------------
Subroutine reducing_funcs(X, &
                          rho_r, T_r, &
                          dvr_dx, dtr_dx, &
                          dvr2_dx2, dtr2_dx2, dvr2_dxx, dtr2_dxx)
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
   ! - T_r    (float): Reducing temperature (Tr)
   use parameters
   Implicit None
   real(8), dimension(N), intent(in) :: X
   real(8), intent(out) :: rho_r, T_r, dvr_dx(N), dtr_dx(N), dvr2_dx2(N), dtr2_dx2(N), &
                           dvr2_dxx(N, N), dtr2_dxx(N, N)
   integer :: i, j, k
   ! Internal variables
   real(8) :: c_v, c_t, xki, xik, Bv_xki, Bt_xki, &
              fv_ki, ft_ki, Bv_xik, Bt_xik, fv_ik, ft_ik, &
              f2v_ki, f2t_ki, f2v_ik, f2t_ik, xij, Bv_xij, Bt_xij
   call get_params()

   rho_r = sum(X**2/rho_c)
   T_r = sum(X**2*T_c)

   do i = 1, N - 1
   if (X(i) > eps) then
   do j = i + 1, N
   if (X(j) > eps) then
      rho_r = rho_r + &
              2*X(i)*X(j)*Bv(i, j)*Gv(i, j) &
              *(X(i) + X(j))/(Bv(i, j)**2*X(i) + X(j)) &
              *1.d0/8*(rho_c(i)**(-1.d0/3) &
                       + rho_c(j)**(-1.d0/3))**3

      T_r = T_r + &
            2*X(i)*X(j)*Bt(i, j)*Gt(i, j) &
            *(X(i) + X(j))/(Bt(i, j)**2*X(i) + X(j)) &
            *sqrt((T_c(i)*T_c(j)))
   end if
   end do
   end if
   end do
   ! Actually, the reducing volume was calculated, let's make it into density
   rho_r = 1.d0/rho_r

   ! Compositional derivatives
   dvr_dx = 2*X/rho_c
   dtr_dx = 2*X*T_c

   dvr2_dx2 = 2/rho_c
   dtr2_dx2 = 2*T_c

   do i = 1, N
   if (X(i) > eps) then
   do k = 1, i - 1
   if (X(k) > eps) then
      xki = X(k) + X(i)
      c_t = 2*Bt(k, i)*Gt(k, i)*sqrt(T_c(k)*T_c(i))
      c_v = 2*Bv(k, i)*Gv(k, i)*1.d0/8 &
            *(rho_c(k)**(-1.d0/3) + rho_c(i)**(-1.d0/3))**3

      Bv_xki = Bv(k, i)**2*X(k) + X(i)
      Bt_xki = Bt(k, i)**2*X(k) + X(i)

      fv_ki = (X(k)*xki + X(k)*X(i)*(1 - xki/Bv_xki))/Bv_xki
      ft_ki = (X(k)*xki + X(k)*X(i)*(1 - xki/Bt_xki))/Bt_xki

      f2v_ki = (1 - xki/Bv_xki)*(2*X(k) - X(k)*X(i)*2/Bv_xki)/Bv_xki
      f2t_ki = (1 - xki/Bt_xki)*(2*X(k) - X(k)*X(i)*2/Bt_xki)/Bt_xki

      dvr_dx(i) = dvr_dx(i) + c_v*fv_ki
      dtr_dx(i) = dtr_dx(i) + c_t*ft_ki
      dvr2_dx2(i) = dvr2_dx2(i) + c_v*f2v_ki
      dtr2_dx2(i) = dtr2_dx2(i) + c_t*f2t_ki

      ! Changes to the crossed derivatives
      ! ----------------------------------
      Bv_xik = Bv(i, k)**2*X(i) + X(k)
      Bt_xik = Bt(i, k)**2*X(i) + X(k)

      c_t = 2*Bt(i, k)*Gt(i, k)*sqrt(T_c(i)*T_c(k))
      c_v = 2*Bv(i, k)*Gv(i, k)*1.d0/8 &
            *(rho_c(i)**(-1.d0/3) + rho_c(k)**(-1.d0/3))**3

      ! Crossed derivatives
      dvr2_dxx(i, k) = c_v*( &
                       xki + X(k)*(1 - xki/Bv_xik) &
                       + X(k)*(1 - xki/Bv_xik) &
                       + X(i)*(1 - Bv(i, k)**2*xki/Bv_xik) &
                       - X(i)*X(k)/Bv_xik*(1 + Bv(i, k) - 2*Bv(i, k)**2*xki/Bv_xik) &
                       )/Bv_xik

      dtr2_dxx(i, k) = c_t*( &
                       xki + X(k)*(1 - xki/Bt_xik) &
                       + X(k)*(1 - xki/Bt_xik) &
                       + X(i)*(1 - Bt(i, k)**2*xki/Bt_xik) &
                       - X(i)*X(k)/Bt_xik*(1 + Bt(i, k) - 2*Bt(i, k)**2*xki/Bt_xik) &
                       )/Bt_xik
   end if
   end do
   do k = i + 1, N
   if (X(k) > eps) then
      xik = X(i) + X(k)
      c_t = 2*Bt(i, k)*Gt(i, k)*sqrt(T_c(i)*T_c(k))
      c_v = 2*Bv(i, k)*Gv(i, k)*1.d0/8 &
            *(rho_c(i)**(-1.d0/3) + rho_c(k)**(-1.d0/3))**3

      Bv_xik = Bv(i, k)**2*X(i) + X(k)
      Bt_xik = Bt(i, k)**2*X(i) + X(k)

      fv_ik = (X(k)*xik + X(i)*X(k)*(1 - Bv(i, k)**2*xik/Bv_xik))/Bv_xik
      ft_ik = (X(k)*xik + X(i)*X(k)*(1 - Bt(i, k)**2*xik/Bt_xik))/Bt_xik

      f2v_ik = (1 - Bv(i, k)**2*xik/Bv_xik)*(2*X(k) - X(i)*X(k)*2*Bv(i, k)**2/Bv_xik)/Bv_xik
      f2t_ik = (1 - Bt(i, k)**2*xik/Bt_xik)*(2*X(k) - X(i)*X(k)*2*Bt(i, k)**2/Bt_xik)/Bt_xik

      dvr_dx(i) = dvr_dx(i) + c_v*fv_ik
      dtr_dx(i) = dtr_dx(i) + c_t*ft_ik
      dvr2_dx2(i) = dvr2_dx2(i) + c_v*f2v_ki
      dtr2_dx2(i) = dtr2_dx2(i) + c_t*f2t_ki

      dvr2_dxx(i, k) = c_v*( &
                       xki + X(k)*(1 - xki/Bv_xik) &
                       + X(k)*(1 - xki/Bv_xik) &
                       + X(i)*(1 - Bv(i, k)**2*xki/Bv_xik) &
                       - X(i)*X(k)/Bv_xik*(1 + Bv(i, k) - 2*Bv(i, k)**2*xki/Bv_xik) &
                       )/Bv_xik

      dtr2_dxx(i, k) = c_t*( &
                       xki + X(k)*(1 - xki/Bt_xik) &
                       + X(k)*(1 - xki/Bt_xik) &
                       + X(i)*(1 - Bt(i, k)**2*xki/Bt_xik) &
                       - X(i)*X(k)/Bt_xik*(1 + Bt(i, k) - 2*Bt(i, k)**2*xki/Bt_xik) &
                       )/Bt_xik
   end if
   end do
   end if
   end do

End Subroutine reducing_funcs

! -----------------------------------------------------------------------------
! Pure Compound Helmholtz Energy (and derivatives) Calculations
! -----------------------------------------------------------------------------
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

   r = 8.314510d0/8.314472d0

   aoio = 0
   Tr = T_c/T
   Dr = rho/rho_c

   do k = 4, 7
      if (v(k) > eps) then
         if (k == 4 .or. k == 6) then
            aoio(1, 1) = aoio(1, 1) + r*(n(k)*log(abs(dsinh(v(k)*Tr))))
            aoio(2, 2) = aoio(2, 2) + r*(n(k)*v(k)/dtanh(v(k)*Tr))
            aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/dsinh(v(k)*Tr)**2)
         else
            aoio(1, 1) = aoio(1, 1) - r*(n(k)*log(dcosh(v(k)*Tr)))
            aoio(2, 2) = aoio(2, 2) - r*(n(k)*v(k)*dtanh(v(k)*Tr))
            aoio(3, 2) = aoio(3, 2) - n(k)*(v(k)/dcosh(v(k)*Tr))**2
         end if
      end if
   end do

   aoio(1, 1) = aoio(1, 1) + log(Dr) + r*(n(1) + n(2)*Tr + n(3)*log(Tr))
   aoio(2, 1) = 1.d0/Dr
   aoio(3, 1) = -1.d0/Dr**2

   aoio(2, 2) = aoio(2, 2) + r*(n(2) + n(3)/Tr)

   aoio(3, 2) = aoio(3, 2) - (n(3)*(T/T_c)**2)
   aoio(3, 2) = r*aoio(3, 2)

   aoio(3, 3) = 0
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
   real(8), intent(in) :: delta, tau
   integer, intent(in) :: Kpol, Kexp
   real(8), dimension(24), intent(in) :: n, t
   integer, dimension(24), intent(in) :: d, c
   real(8), dimension(3, 3), intent(out) :: aoir
   integer :: k

   aoir = 0
   ! Residual Helmholtz Energy
   do k = 1, Kpol
      aoir(1, 1) = aoir(1, 1) + &
                   n(k)*delta**d(k) &
                   *tau**t(k)
      aoir(2, 1) = aoir(2, 1) + n(k)*d(k)*delta**(d(k) - 1)*tau**(t(k))
      aoir(2, 2) = aoir(2, 2) + &
                   n(k)*t(k)*delta**d(k) &
                   *tau**(t(k) - 1)
      aoir(3, 1) = aoir(3, 1) + &
                   n(k)*d(k)*(d(k) - 1)*delta**(d(k) - 2)*tau**t(k)
      aoir(3, 2) = aoir(3, 2) + &
                   n(k)*t(k)*(t(k) - 1) &
                   *delta**d(k)*tau**(t(k) - 2)
      aoir(3, 3) = aoir(3, 3) + &
                   n(k)*d(k)*t(k) &
                   *delta**(d(k) - 1)*tau**(t(k) - 1)
   end do
   do k = Kpol + 1, Kpol + Kexp
      aoir(1, 1) = aoir(1, 1) + &
                   n(k)*delta**d(k) &
                   *tau**t(k) &
                   *exp(-delta**c(k))
      ! First Derivative with reduced density
      aoir(2, 1) = aoir(2, 1) + &
                   n(k)*delta**(d(k) - 1) &
                   *(d(k) - c(k)*delta**c(k)) &
                   *tau**t(k)*exp(-delta**c(k))
      ! First Derivative with reduced temperature
      aoir(2, 2) = aoir(2, 2) + &
                   n(k)*t(k)*delta**d(k) &
                   *tau**(t(k) - 1) &
                   *exp(-delta**c(k))
      ! Second Derivative with reduced density
      aoir(3, 1) = aoir(3, 1) + &
                   n(k)*delta**(d(k) - 2)*( &
                   (d(k) - c(k)*delta**c(k)) &
                   *(d(k) - 1 - c(k)*delta**c(k)) &
                   - c(k)**2*delta**c(k)) &
                   *tau**t(k)*exp(-delta**c(k))
      ! Second Derivative with reduced temperature
      aoir(3, 2) = aoir(3, 2) + &
                   n(k)*t(k)*(t(k) - 1) &
                   *delta**d(k)*tau**(t(k) - 2) &
                   *exp(-delta**c(k))
      ! Second Derivative with reduced temperature and density
      aoir(3, 3) = aoir(3, 3) + &
                   n(k)*t(k)*delta**(d(k) - 1) &
                   *(d(k) - c(k)*delta**c(k)) &
                   *tau**(t(k) - 1)*exp(-delta**c(k))
   end do
End Subroutine a_oir

! -----------------------------------------------------------------------------
! Binary interaction contribution
! -----------------------------------------------------------------------------
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

   aijr = 0
   do k = 1, Kpolij
      aijr(1, 1) = aijr(1, 1) + &
                   n(k)*delta**d(k)*tau**t(k)
      aijr(2, 1) = aijr(2, 1) + &
                   n(k)*d(k)*delta**(d(k) - 1)*tau**t(k)
      aijr(3, 1) = aijr(3, 1) + &
                   n(k)*d(k)*(d(k) - 1)*delta**(d(k) - 2)*tau**t(k)
      aijr(2, 2) = aijr(2, 2) + &
                   n(k)*t(k)*delta**d(k)*tau**(t(k) - 1)
      aijr(3, 2) = aijr(3, 2) + &
                   n(k)*t(k)*(t(k) - 1)*delta**d(k)*tau**(t(k) - 2)
      aijr(3, 3) = aijr(3, 3) + &
                   n(k)*d(k)*t(k) &
                   *delta**(d(k) - 1)*tau**(t(k) - 1)
   end do
   do k = Kpolij + 1, Kpolij + Kexpij
      aijr(1, 1) = aijr(1, 1) + &
                   n(k)*delta**d(k)*tau**t(k) &
                   *exp( &
                   -eta(k)*(delta - eps(k))**2 &
                   - beta(k)*(delta - gamm(k)) &
                   )
      aijr(2, 1) = aijr(2, 1) + &
                   n(k)*delta**d(k)*tau**t(k) &
                   *exp( &
                   -eta(k)*(delta - eps(k))**2 &
                   - beta(k)*(delta - gamm(k)) &
                   ) &
                   *(d(k)/delta - 2*eta(k)*(delta - eps(k)) - beta(k))
      aijr(3, 1) = aijr(3, 1) + &
                   n(k)*delta**d(k)*tau**t(k)*exp( &
                   -eta(k)*(delta - eps(k))**2 - beta(k)*(delta - gamm(k))) &
                   *((d(k)/delta - 2*eta(k) &
                      *(delta - eps(k)) - beta(k))**2 - d(k)/delta**2 - 2*eta(k))
      aijr(2, 2) = aijr(2, 2) + &
                   n(k)*t(k)*delta**d(k)*tau**(t(k) - 1) &
                   *exp( &
                   -eta(k)*(delta - eps(k))**2 &
                   - beta(k)*(delta - gamm(k)) &
                   )
      aijr(3, 2) = aijr(3, 2) + &
                   n(k)*t(k)*(t(k) - 1)*delta**d(k)*tau**(t(k) - 2) &
                   *exp( &
                   -eta(k)*(delta - eps(k))**2 - beta(k)*(delta - gamm(k)) &
                   )
      aijr(3, 3) = aijr(3, 3) + &
                   n(k)*t(k)*delta**d(k)*tau**(t(k) - 1) &
                   *exp( &
                   -eta(k)*(delta - eps(k))**2 &
                   - beta(k)*(delta - gamm(k)) &
                   ) &
                   *(d(k)/delta - 2*eta(k)*(delta - eps(k)) - beta(k))
   end do
End Subroutine a_ijr

! -----------------------------------------------------------------------------
! GERG-2008 TERMS
! -----------------------------------------------------------------------------
Subroutine ideal_term(X, rho, T, rho_r, T_r, ao)
   use parameters
   real(8), intent(in) :: X(N), rho, T, rho_r, T_r
   real(8), intent(out) :: ao(3, 3)
   real(8) ::  aoio(3, 3), xi

   call get_params()

   ao = 0
   do i = 1, N
      if (X(i) > eps) then
         aoio = 0
         call a_oio(rho, T, rho_c(i), T_c(i), n0i(i, :), th0i(i, :), aoio)

         xi = X(i)
         ao(1, 1) = ao(1, 1) + xi*(aoio(1, 1) + log(xi))

         ao(2, 1) = ao(2, 1) + xi*rho_r/rho_c(i)*aoio(2, 1)
         ao(3, 1) = ao(3, 1) + xi*(rho_r/rho_c(i))**2*aoio(3, 1)

         ao(2, 2) = ao(2, 2) + xi*(T_c(i)/T_r)*aoio(2, 2)
         ao(3, 2) = ao(3, 2) + xi*(T_c(i)/T_r)**2*aoio(3, 2)

         ao(3, 3) = 0
      end if
   end do
End Subroutine ideal_term

Subroutine residual_term(X, delta, tau, ar, ar_x, ar_dx, ar_tx, ar_xx)
   use parameters
   Implicit None
   real(8), intent(in) :: delta, tau, X(N)
   real(8), dimension(3, 3), intent(out) :: ar
   real(8), dimension(N), intent(out) :: ar_x, ar_dx, ar_tx
   real(8), dimension(N, N), intent(out) :: ar_xx
   real(8), dimension(3, 3) :: aoir, aijr
   integer :: i, j, k

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

         ar(1, 1) = ar(1, 1) + X(i)*aoir(1, 1)
         ar(2, 1) = ar(2, 1) + X(i)*aoir(2, 1)
         ar(2, 2) = ar(2, 2) + X(i)*aoir(2, 2)
         ar(3, 1) = ar(3, 1) + X(i)*aoir(3, 1)
         ar(3, 2) = ar(3, 2) + X(i)*aoir(3, 2)
         ar(3, 3) = ar(3, 3) + X(i)*aoir(3, 3)

         ar_x(i) = aoir(1, 1)
         ar_dx(i) = aoir(2, 1)
         ar_tx(i) = aoir(2, 2)

         ! Compositional derivatives
         do k = 1, N
            if (i /= k .and. Fij(i, k) > eps) then
               call a_ijr(delta, tau, Kpolij(i, k), Kexpij(i, k), &
                          nij(i, k, :), dij(i, k, :), &
                          tij(i, k, :), ethaij(i, k, :), &
                          epsij(i, k, :), gammaij(i, k, :), &
                          betaij(i, k, :), aijr)

               ar_x(i) = ar_x(i) + X(k)*Fij(i, k)*aijr(1, 1)
               ar_dx(i) = ar_dx(i) + X(k)*Fij(i, k)*aijr(2, 1)
               ar_tx(i) = ar_tx(i) + X(k)*Fij(i, k)*aijr(2, 2)
               ar_xx(i, k) = Fij(i, k)*aijr(1, 1)
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
         ar(1, 1) = ar(1, 1) + X(i)*X(j)*Fij(i, j)*aijr(1, 1)
         ar(2, 1) = ar(2, 1) + X(i)*X(j)*Fij(i, j)*aijr(2, 1)
         ar(2, 2) = ar(2, 2) + X(i)*X(j)*Fij(i, j)*aijr(2, 2)
         ar(3, 1) = ar(3, 1) + X(i)*X(j)*Fij(i, j)*aijr(3, 1)
         ar(3, 2) = ar(3, 2) + X(i)*X(j)*Fij(i, j)*aijr(3, 2)
         ar(3, 3) = ar(3, 3) + X(i)*X(j)*Fij(i, j)*aijr(3, 3)
      end if
      end do
   end do
End Subroutine residual_term

! -----------------------------------------------------------------------------
! MAIN PROGRAM
! -----------------------------------------------------------------------------
Program main
   print *, "This only exists to test compile with gfortran"
End Program main

! Author= Federico Benelli
! Gerg-2008 Equation
! Started at: 01/07/2021
! Last Modified: jue 26 ago 2021 10:56:11
!

! -------------------
! GERG 2008 Functions
! -------------------
Subroutine reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)
    ! REDUCING DENSITY AND TEMPERATURE
    ! input:
    !  - X      (dimension): molar fractions
    !  - Bv     (dimension): beta parameters for reducing density
    !  - Gv     (dimension): gamma parameters for reducing density
    !  - Bt     (dimension): beta parameters for reducing temperature
    !  - Gt     (dimension): gamma parameters for reducing temperature
    !  - rho_c  (dimension): citical densities
    !  - T_c    (dimension): citical temperatures
    ! output:
    !  - rho_r  (float): Reducing density (1/rhor)
    !  - T_r    (float): reducing temperature (Tr)
    real*8, dimension(21), intent(in):: X, rho_c, T_c
    real*8, dimension(21, 21), intent(in)::  Bv, Gv, Bt, Gt
    real*8, intent(out):: rho_r, T_r
    integer:: N = 21, i, j
    real*8:: eps = 1d-10

    rho_r = 0.d0
    T_r = 0.d0

    do i = 1, N
    if (X(i) > eps) then
        rho_r = rho_r + X(i)**2/rho_c(i)
        T_r = T_r + X(i)**2*T_c(i)
    end if
    end do

    do i = 1, N - 1
    if (X(i) > eps) then
    do j = i + 1, N
    if (X(j) > eps) then
        rho_r = rho_r + &
                2.d0*X(i)*X(j)*Bv(i, j)*Gv(i, j) &
                *(X(i) + X(j))/(Bv(i, j)**2.d0*X(i) + X(j)) &
                *1.d0/8.d0*(rho_c(i)**(-1.d0/3.d0) + rho_c(j)**(-1.d0/3.d0))**3

        T_r = T_r + &
              2.d0*X(i)*X(j)*Bt(i, j)*Gt(i, j) &
              *(X(i) + X(j))/(Bt(i, j)**2.d0*X(i) + X(j)) &
              *sqrt((T_c(i)*T_c(j)))
    end if
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
    !  - rho (float): Density
    !  - T (float): Temperature
    !  - rho_c (float): Critical Density
    !  - T_c (float): Critical Temperature
    !  - n (dimension): n parameters
    !  - v (dimension): v parameters
    ! output:
    ! - aoio: Ideal Gas Helmholtz Energy and its derivatives,
    !           stoikuctured like:
    ! --------------+------------------+------------------|
    ! aoio          |     0            |    0             |
    ! d(aoio)/dd    |     d(aoio)/dt   |    0             |
    ! d2(aoio)/dd2  |     d2(aoio)/dt2 |    d2(aoio)/dddt |
    ! ----------------------------------------------------|
    integer, parameter:: dp = 8
    real(dp), intent(in):: rho, T, rho_c, T_c
    real(dp), intent(in), dimension(7):: n, v
    real(dp), intent(out):: aoio(3, 3)
    real(dp):: r, Tr, Dr, eps = 1d-10

    r = 8.314510d0/8.314472d0

    aoio = 0.d0
    Tr = T_c/T
    Dr = rho/rho_c

    do k = 4, 7
        !write (0, *) 'Ideal Terms: ', k, n(k), v(k)
        if (v(k) > eps) then
        if (k == 4) then
            aoio(1, 1) = aoio(1, 1) + r*(n(k)*log(abs(sinh(v(k)*Tr))))
            aoio(2, 2) = aoio(2, 2) + r*(n(k)*v(k)/tanh(v(k)*Tr))
            aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/sinh(v(k)*Tr)**2)
        else if (k == 5 .or. k == 6) then
            aoio(1, 1) = aoio(1, 1) + r*(n(k)*log(abs(sinh(v(k)*Tr))))
            aoio(1, 1) = aoio(1, 1) - r*(n(k)*log(cosh(v(k)*Tr)))

            aoio(2, 2) = aoio(2, 2) + r*(n(k)*v(k)/tanh(v(k)*Tr))
            aoio(2, 2) = aoio(2, 2) - r*(n(k)*v(k)*tanh(v(k)*Tr))

            aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/sinh(v(k)*Tr)**2)
            aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/cosh(v(k)*Tr)**2)
        else
            aoio(1, 1) = aoio(1, 1) - r*(n(k)*log(cosh(v(k)*Tr)))
            aoio(2, 2) = aoio(2, 2) - r*(n(k)*v(k)*tanh(v(k)*Tr))
            aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/cosh(v(k)*Tr)**2)
        end if
        end if
    end do

    aoio(1, 1) = aoio(1, 1) + log(Dr) + r*(n(1) + n(2)*Tr + n(3)*log(Tr))
    aoio(2, 1) = 1.d0/Dr
    aoio(3, 1) = -1.d0/Dr**2
    aoio(2, 2) = aoio(2, 2) + r*(n(2) + n(3)/Tr)
    aoio(3, 2) = aoio(3, 2) + r*(-n(3)/Tr**2)
    aoio(3, 3) = 0.d0

End Subroutine a_oio

Subroutine a_oir(delta, tau, Kpol, Kexp, n, d, t, c, aoir)
    ! RESIDUAL ENERGY
    ! ---------------
    ! Caculate the pure compound Residual Helmholtz Energy and its derivatives
    !
    ! input:
    !  - delta  (float)         : reduced density
    !  - tau    (float)         : reduced temperature
    !  - Kpol   (integer)       : number of only polynomial parameters
    !  - Kexp   (integer)       : number of only exponential parameters
    !  - n      (dimension)     : parameters n
    !  - d      (dimension)     : parameters d
    !  - t      (dimension)     : parameters t
    !  - c      (dimension)     : parameters c
    ! output:
    ! - aoir: Ideal Gas Helmholtz Energy and its derivatives,
    !           structured like:
    ! ---------------+-------------------+-----------------|
    ! aoir           |      0            |   0             |
    ! d(aoir)/dd     |      d(aoir)/dt   |   0             |
    ! d2(aoir)/dd2   |      d2(aoir)/dt2 |   d2(aoir)/dddt |
    ! -----------------------------------------------------|
    !

    real*8, intent(in):: delta, tau
    integer, intent(in):: Kpol, Kexp
    real*8, dimension(24), intent(in):: n, t
    integer, dimension(24), intent(in):: d, c
    real*8, dimension(3, 3), intent(out):: aoir
    integer:: k

    aoir = 0.d0
    ! Residual Helmholtz Energy
    do k = 1, Kpol
        !!write (0, *) k, n(k), d(k), t(k), c(k)
        aoir(1, 1) = aoir(1, 1) + &
                     n(k)*delta**d(k) &
                     *tau**t(k)
        aoir(2, 1) = aoir(2, 1) + n(k)*d(k)*delta**(d(k) - 1)*tau**(t(k))
        aoir(2, 2) = aoir(2, 2) + &
                     n(k)*t(k)*delta**d(k) &
                     *tau**(t(k) - 1.d0)
        aoir(3, 1) = aoir(3, 1) + &
                     n(k)*d(k)*(d(k) - 1)*delta**(d(k) - 2)*tau**t(k)
        aoir(3, 2) = aoir(3, 2) + &
                     n(k)*t(k)*(t(k) - 1) &
                     *delta**d(k)*tau**(t(k) - 2.d0)
        aoir(3, 3) = aoir(3, 3) + &
                     n(k)*d(k)*t(k) &
                     *delta**(d(k) - 1)*tau**(t(k) - 1.d0)
    end do
    do k = Kpol + 1, Kpol + Kexp
        !!write (0, *) k, n(k), d(k), t(k), c(k)
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
                     n(k)*t(k)*(t(k) - 1.d0)*delta**d(k)*tau**(t(k) - 2.d0)*exp(-delta**c(k))
        ! Second Derivative with reduced temperature and density
        aoir(3, 3) = aoir(3, 3) + &
                     n(k)*t(k)*delta**(d(k) - 1) &
                     *(d(k) - c(k)*delta**c(k)) &
                     *tau**(t(k) - 1.d0)*exp(-delta**c(k))
    end do

End Subroutine a_oir

Subroutine a_ijr(delta, tau, Kpolij, Kexpij, &
                 n, d, t, eta, eps, gamm, beta, aijr)
    ! DEPARTURE FUNCTION ENERGY
    ! ---------------
    ! Caculate binary departure Helmholtz Energy and its derivatives
    !
    ! input:
    !  - delta  (float)         : reduced density
    !  - tau    (float)         : reduced temperature
    !  - Kpolij (integer)       : number of only polynomial parameters
    !  - Kexpij (integer)       : number of only exponential parameters
    !  - n      (dimension)     : parameters n
    !  - d      (dimension)     : parameters d
    !  - t      (dimension)     : parameters t
    !  - eta    (dimension)     : parameters eta
    !  - eps    (dimension)     : parameters epsilon
    !  - gamm   (dimension)     : parameters gamma
    !  - beta   (dimension)     : parameters beta
    ! output:
    ! - aijr: Binary departure Helmholtz Energy and its derivatives,
    !           stoikuctured like:
    ! ---------------+------------------+------------------|
    ! aijr           |      0            |   0             |
    ! d(aijr)/dd     |      d(aijr)/dt   |   0             |
    ! d2(aijr)/dd2   |      d2(aijr)/dt2 |   d2(aijr)/dddt |
    ! -----------------------------------------------------|
    real*8, intent(in):: delta, tau
    integer, intent(in):: Kpolij, Kexpij
    real*8, dimension(24), intent(in):: n, t, eta, eps, gamm, beta
    integer, dimension(24), intent(in):: d
    real*8, dimension(3, 3), intent(out):: aijr
    integer:: k

    aijr = 0.0
    ! Departure function
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
                     *((d(k)/delta - 2*eta(k)*(delta - eps(k)) - beta(k))**2 - d(k)/delta**2 - 2*eta(k))
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

Subroutine ideal_term(X, rho, T, rho_c, T_c, rho_r, T_r, &
                      n0i, th0i, ao)
    real*8:: X(21), rho, T, rho_r, T_r
    real*8, dimension(21, 7):: n0i, th0i
    real*8, dimension(3, 3):: ao, aoio
    real*8, dimension(21):: rho_c, T_c
    real*8::  xi, eps = 1d-10
    integer:: N = 21

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

            ! ao(3,3) should always be zero, based on the paper
        end if
    end do

End Subroutine ideal_term

Subroutine residual_term(X, delta, tau, &
                         Kpol, Kexp, noik, doik, toik, coik, &
                         Fij, Kpolij, Kexpij, dij, nij, tij, betaij, epsij, ethaij, gammaij, &
                         ar)
    real*8, dimension(21), intent(in):: X
    integer, dimension(21), intent(in):: Kpol, Kexp
    real*8, dimension(21, 24), intent(in):: noik, toik
    integer, dimension(21, 24), intent(in):: doik, coik
    real*8, dimension(21, 21, 12), intent(in):: betaij, epsij, ethaij, gammaij
    real*8, dimension(21, 21, 12), intent(in):: nij, tij
    integer, dimension(21, 21, 12), intent(in):: dij
    real*8, dimension(21, 21), intent(in):: Fij
    integer, dimension(21, 21), intent(in):: Kpolij, Kexpij
    real*8, intent(in):: delta, tau
    real*8, dimension(3, 3), intent(out):: ar
    real*8, dimension(3, 3):: aoir, aijr
    real*8:: eps = 1d-10

    integer:: i, j

    ar = 0.0

    do i = 1, 21
        if (X(i) > eps) then
            aoir = 0.0
            aijr = 0.0
            call a_oir(delta, tau, Kpol(i), Kexp(i), &
                       noik(i, :), doik(i, :), toik(i, :), coik(i, :), aoir)
            ar(1, 1) = ar(1, 1) + X(i)*aoir(1, 1)
            ar(2, 1) = ar(2, 1) + X(i)*aoir(2, 1)
            ar(2, 2) = ar(2, 2) + X(i)*aoir(2, 2)
            ar(3, 1) = ar(3, 1) + X(i)*aoir(3, 1)
            ar(3, 2) = ar(3, 2) + X(i)*aoir(3, 2)
            ar(3, 3) = ar(3, 3) + X(i)*aoir(3, 3)
        end if
    end do

    !write (0, *) "Normal:", ar(2, 1), ar(3, 1), ar(2, 2), ar(3, 2), ar(3, 3)

    do i = 1, 20
    do j = i + 1, 21
    if (Fij(i, j) > eps .and. X(i) > eps .and. X(j) > eps) then
        !!write (0, *) "_Departure Function_"
        !!write (0, *) i, j, Fij(i, j), Kpolij(i, j), Kexpij(i, j)
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
    !write (0, *) "Departure:", aijr(2, 1), aijr(3, 1), aijr(2, 2), aijr(3, 2), aijr(3, 3)

End Subroutine residual_term

Program gerg
    use parameters
    use thermo_props

    Implicit None
    ! Constants
    real*8:: R
    ! Variables
    real*8:: T, rho, tau, delta, X(21)
    ! Calculated variables
    real*8, dimension(3, 3):: ar, ao
    real*8:: P, Z, w, cp, cv, mean_M
    real*8:: T_r, rho_r
    character(len=100):: compound

    integer:: io, counter, stdin = 5, stdout = 6, stderr = 0
    character(len=100):: arg

    real*8:: test_P, test_cv, test_cp, test_w, up, down

    ! Constantes
    R = 8.314472d0

    call get_params()

    ! Propiedades de operación
    ! ------------------------
    call getarg(1, arg)
    rho = 1.d0
    T = 1.d0
    X = 0.d0

    counter = 1

    select case (arg)
    case ('bintest')
        open (1, file='GERG_test/binary/test_data')
        io = 0
        do while (io == 0)
            !write (0, *) "___________________________________________________________"
            X = 0
            !write (stderr, *) "Iteracion: ", counter
            counter = counter + 1

            read (1, *, iostat=io) i, X(i), j, X(j), T, rho, test_P, test_cv, test_cp, test_w

            if (io .ne. 0) then
                exit
            end if

            mean_M = 0
            do k = 1, size(X)
                mean_M = mean_M + X(k)*M(k)
            end do
            mean_M = mean_M/1000.d0

            ! Set delta and tau
            call reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)
            delta = rho*rho_r
            tau = T_r/T

            call residual_term(X, delta, tau, &
                               Kpol, Kexp, noik, doik, toik, coik, &
                               Fij, Kpolij, Kexpij, dij, nij, tij, betaij, epsij, ethaij, gammaij, &
                               ar)

            call ideal_term(X, rho, T, rho_c, T_c, rho_r, T_r, &
                            n0i, th0i, &
                            ao)

            ! Calculate Properties
            call zeta(delta, ar, Z)
            call isobaric_heat(delta, tau, R, ao, ar, cp)
            call isocoric_heat(tau, R, Ao, Ar, cv)
            call sound_speed(delta, tau, R, T, mean_M, Ao, Ar, w)

            ! adim * [mol/L] * [K] * [J/(mol*K)] * [L/m3] * [MPa/Pa]
            P = Z*rho*T*R*1000.d0*1d-6

            up = (1.d0 + delta*Ar(2, 1) - delta*tau*Ar(3, 3))**2
            down = (1.d0 + 2.d0*delta*Ar(2, 1) - delta**2*Ar(3, 1))

            print *, i, j, T, rho, P, cv, cp, w, test_P, test_cv, test_cp, test_w
            write (stderr, *) i, j, -tau**2*(Ao(3, 2) + Ar(3, 2)), +up/down
        end do
        close (1)

    case ('isotherm')
        open (1, file='debug_data')

        ! Get Isothermic data of mixture at given T
        read (1, *) T
        io = 0

        ! Get concentration data
        do while (io == 0)
            read (1, *, iostat=io) i, compound, X(i)
            !print *, i, X(i)
        end do
        close (1)

        ! Calculate mean molecular weight
        mean_M = 0
        do i = 1, size(X)
            mean_M = mean_M + X(i)*M(i)
        end do

        call reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)

        tau = T_r/T

        do j = 0, 5000, 10
            rho = float(j)/10
            delta = rho*rho_r
            tau = T_r/T

            call residual_term(X, delta, tau, &
                               Kpol, Kexp, noik, doik, toik, coik, &
                               Fij, Kpolij, Kexpij, dij, nij, tij, betaij, &
                               epsij, ethaij, gammaij, &
                               ar)
            call ideal_term(X, rho, T, rho_c, T_c, rho_r, T_r, &
                            n0i, th0i, &
                            ao)
            call zeta(delta, ar, Z)
            call sound_speed(R, T, mean_M/1000, delta, tau, ar, ao, w)

            ! adim * [mol/L] * [K] * [J/(mol*K)] * [L/m3]
            P = Z*rho*T*R*1000

            !print *, T, rho, P/1e6, Z, w
        end do
    end select
End Program gerg

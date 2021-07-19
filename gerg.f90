! Author= Federico Benelli
! Started at: 01/07/2021
! Last Modified: dom 18 jul 2021 15:33:49
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
    !  - rho_c  (dimension): critical densities
    !  - T_c    (dimension): critical temperatures
    ! output:
    !  - rho_r  (float): Reducing density (1/rhor)
    !  - T_r    (float): reducing temperature (Tr)
    double precision, dimension(21), intent(in):: X, rho_c, T_c
    double precision, dimension(21, 21), intent(in)::  Bv, Gv, Bt, Gt
    double precision, intent(inout):: rho_r, T_r
    double precision:: eps = epsilon(rho_r)
    integer:: N = 21, i, j
    rho_r = 0
    T_r = 0
    do i = 1, N
    if (X(i) > eps) then
        rho_r = rho_r + X(i)**2/rho_c(i)
        T_r = T_r + X(i)**2*T_c(i)
    end if
    end do

    do i = 1, N - 1
    do j = i + 1, N
    if (X(i) > eps .and. X(j) > eps) then
        rho_r = rho_r + &
                2*X(i)*X(j)*Bv(i, j)*Gv(i, j) &
                *(X(i) + X(j))/(Bv(i, j)**2*X(i) + X(j)) &
                *1/8*(rho_c(i)**(-1.0/3.0) + rho_c(j)**(-1.0/3.0))**3

        T_r = T_r + &
              2*X(i)*X(j)*Bt(i, j)*Gt(i, j) &
              *(X(i) + X(j))/(Bt(i, j)**2*X(i) + X(j)) &
              *sqrt((T_c(i)*T_c(j)))
    end if
    end do
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
    !           structured like:
    ! --------------+------------------+------------------|
    ! aoio          |     0            |    0             |
    ! d(aoio)/dd    |     d(aoio)/dt   |    0             |
    ! d2(aoio)/dd2  |     d2(aoio)/dt2 |    d2(aoio)/dddt |
    ! ----------------------------------------------------|

    double precision, intent(in):: rho, T, rho_c, T_c
    double precision, intent(in), dimension(7):: n, v
    double precision, intent(inout):: aoio(3, 3)
    double precision:: delta, tau
    double precision:: r

    r = 8.314510/8.314472

    aoio = 0
    delta = rho/rho_c
    tau = T_c/T

    aoio(1, 1) = log(delta) + r*(n(1) + n(2)*tau + n(3)*log(tau))
    aoio(2, 1) = 1/delta
    aoio(2, 2) = n(2) + n(3)/tau
    aoio(3, 1) = -(1/delta)**2.0
    aoio(3, 2) = n(3)*(T/T_c)**2
    aoio(3, 3) = 0

    do k = 4, 6
        aoio(1, 1) = aoio(1, 1) + r*n(k)*log(abs(sinh(v(k)*tau)))
        aoio(2, 2) = aoio(2, 2) + n(k)*v(k)/tanh(v(k)*tau)
        aoio(3, 2) = aoio(3, 2) + n(k)*v(k)**2/sinh(v(k)*tau)**2
    end do
    do k = 5, 7
        aoio(1, 1) = aoio(1, 1) - r*n(k)*log(cosh(v(k)*tau))
        aoio(2, 2) = aoio(2, 2) - n(k)*v(k)*tanh(v(k)*tau)
        aoio(3, 2) = aoio(3, 2) + n(k)*v(k)**2/cosh(v(k)*tau)**2
    end do

    aoio(2, 2) = r*aoio(2, 2)
    aoio(3, 2) = -r*aoio(3, 2)

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
    ! ---------------+------------------+------------------|
    ! aoir           |      0            |   0             |
    ! d(aoir)/dd     |      d(aoir)/dt   |   0             |
    ! d2(aoir)/dd2   |      d2(aoir)/dt2 |   d2(aoir)/dddt |
    ! -----------------------------------------------------|

    double precision, intent(in):: delta, tau
    integer, intent(in):: Kpol, Kexp
    double precision, dimension(24), intent(in):: n, d, t, c
    double precision, dimension(3, 3), intent(inout):: aoir
    integer:: k

    aoir = 0

    ! Residual Helmholtz Energy
    do k = 1, Kpol
        aoir(1, 1) = aoir(1, 1) + &
                     n(k)*delta**d(k) &
                     *tau**t(k)
        aoir(2, 1) = aoir(2, 1) + &
                     n(k)*d(k) &
                     *delta**(d(k) - 1) &
                     *tau**t(k)
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
                     n(k) &
                     *delta**(d(k) - 1) &
                     *(d(k) - c(k)*delta**c(k)) &
                     *tau**t(k) &
                     *exp(-delta**c(k))
        ! First Derivative with reduced temperature
        aoir(2, 2) = aoir(2, 2) + &
                     n(k)*t(k)*delta**d(k) &
                     *tau**(t(k) - 1) &
                     *exp(-delta**c(k))
        ! Second Derivative with reduced density
        aoir(3, 1) = aoir(3, 1) + &
                     n(k)*delta**(d(k) - 2) &
                     *( &
                     (d(k) - c(k)*delta**c(k)) &
                     *(d(k) - 1 - c(k)*delta**c(k)) - c(k)**2*delta**c(k)) &
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
    !           structured like:
    ! ---------------+------------------+------------------|
    ! aijr           |      0            |   0             |
    ! d(aijr)/dd     |      d(aijr)/dt   |   0             |
    ! d2(aijr)/dd2   |      d2(aijr)/dt2 |   d2(aijr)/dddt |
    ! -----------------------------------------------------|
    double precision, intent(in):: delta, tau
    integer, intent(in):: Kpolij, Kexpij
    double precision, dimension(24), intent(in):: n, d, t, eta, eps, gamm, beta
    double precision, dimension(3, 3), intent(inout):: aijr
    integer:: k

    aijr = 0
    ! Departure function
    do k = 1, Kpolij
        aijr(1, 1) = aijr(1, 1) + &
                     n(k)*delta**d(k)*tau**t(k)
        aijr(2, 1) = aijr(2, 1) + &
                     n(k)*d(k)*delta**(d(k) - 1)*tau**t(k)
        aijr(3, 1) = aijr(3, 1) + &
                     n(k)*d(k)*(d(k) - 1) &
                     *delta**(d(k) - 2)*tau**t(k)
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
                     n(k)*delta**d(k)*tau**t(k) &
                     *exp( &
                     -eta(k)*(delta - eps(k))**2 &
                     - beta(k)*(delta - gamm(k)) &
                     ) &
                     *((d(k)/delta - 2*eta(k)*(delta - eps(k)) - beta(k))**2 &
                       - d(k)/delta**2 - 2*eta(k))
        aijr(2, 2) = aijr(2, 2) + &
                     n(k)*t(k)*delta**d(k)*tau**(t(k) - 1) &
                     *exp( &
                     -eta(k)*(delta - eps(k))**2 &
                     - beta(k)*(delta - gamm(k)) &
                     )
        aijr(3, 2) = aijr(3, 2) + &
                     n(k)*t(k)*(t(k) - 1)*delta**d(k)*tau**(t(k) - 2) &
                     *exp( &
                     -eta(k)*(delta - eps(k))**2 &
                     - beta(k)*(delta - gamm(k)) &
                     )
        ! TODO esta derivada estan bien en el paper??
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
                      no, vo, ao)
    double precision:: X(21), rho, T, rho_r, T_r
    double precision, dimension(21, 7):: no, vo
    double precision, dimension(3, 3):: ao, aoio
    double precision, dimension(21):: rho_c, T_c
    double precision:: eps = epsilon(rho)
    integer:: N = 21

    ao = 0
    do i = 1, N
        if (X(i) > eps) then
            call a_oio(rho, T, rho_c(i), T_c(i), no(i, :), vo(i, :), aoio)
            ao(1, 1) = ao(1, 1) + X(i)*(aoio(1, 1) + log(X(i)))
            ao(2, 1) = ao(2, 1) + X(i)/(rho_c(i)*rho_r)*aoio(2, 1)
            ao(2, 2) = ao(2, 2) + X(i)*(T_c(i)/T_r)*aoio(2, 2)
            ao(3, 1) = ao(3, 1) + X(i)*(1/(rho_c(i)/rho_r)**2)*aoio(3, 1)
            ao(3, 2) = ao(3, 2) + X(i)*(T_c(i)/T_r)**2*aoio(3, 2)
            ao(3, 3) = ao(3, 3) + X(i)*(T_c(i)/T_r)*(1/(rho_c(i)/rho_r))*aoio(3, 3)
        end if
    end do

End Subroutine ideal_term

Subroutine residual_term(X, delta, tau, &
                         Kpol, Kexp, nr, dr, tr, cr, &
                         Fij, Kpolij, Kexpij, dij, nij, tij, betaij, epsij, etaij, gammij, &
                         ar)
    double precision, dimension(21):: X
    integer, dimension(21):: Kpol, Kexp
    double precision, dimension(21, 24):: nr, dr, tr, cr
    double precision, dimension(21, 21, 12):: betaij, epsij, etaij, gammij
    double precision, dimension(21, 21, 12):: dij, nij, tij
    double precision, dimension(21, 21):: Fij
    integer, dimension(21, 21):: Kpolij, Kexpij
    double precision::delta, tau
    double precision, dimension(3, 3):: ar, aoir, aijr
    double precision:: eps = epsilon(X(1))

    integer:: i, j

    ar = 0

    do i = 1, size(X)
        if (X(i) > eps) then
            call a_oir(delta, tau, Kpol(i), Kexp(i), &
                       nr(i, :), dr(i, :), tr(i, :), cr(i, :), aoir)
            ar(1, 1) = ar(1, 1) + X(i)*aoir(1, 1)
            ar(2, 1) = ar(2, 1) + X(i)*aoir(2, 1)
            ar(2, 2) = ar(2, 2) + X(i)*aoir(2, 2)
            ar(3, 1) = ar(3, 1) + X(i)*aoir(3, 1)
            ar(3, 2) = ar(3, 2) + X(i)*aoir(3, 2)
            ar(3, 3) = ar(3, 3) + X(i)*aoir(3, 3)
        end if
    end do

    do i = 1, size(X) - 1
    do j = i + 1, size(X)
    if (Fij(i, j) > eps .and. X(i) > eps .and. X(j) > eps) then
        call a_ijr(delta, tau, Kpolij(i, j), Kexpij(i, j), &
                   nij(i, j, :), dij(i, j, :), tij(i, j, :), etaij(i, j, :), &
                   epsij(i, j, :), gammij(i, j, :), betaij(i, j, :), &
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

Program gerg
    use get_params
    use thermo_props
    Implicit None
    ! Ideal gas parameters
    character(len=100):: no_file, vo_file
    double precision, dimension(21, 7):: no, vo
    ! Residual energy parameteres
    character(len=100)::Kr_file, nr_file, dr_file, tr_file, cr_file
    double precision, dimension(21, 24):: nr, dr, tr, cr
    integer:: Kpol(21), Kexp(21)
    ! Departure function parameters
    character(len=100):: Fij_file, Kij_file, betaij_file, dij_file, epsij_file, &
                         etaij_file, gammij_file, nij_file, tij_file
    integer:: Kpolij(21, 21), Kexpij(21, 21)
    double precision, dimension(21, 21, 12):: betaij, epsij, etaij, gammij
    double precision, dimension(21, 21, 12):: dij, nij, tij
    double precision, dimension(21, 21):: Fij
    ! Reducing function parameters
    character(len=100):: reducing_file
    double precision, dimension(21, 21):: Bv, Gv, Bt, Gt
    ! Critical properties parameters
    character(len=100):: critical_file
    double precision, dimension(21):: rho_c, T_c, M
    ! Constants
    double precision:: R
    ! Variables
    double precision:: T, rho, tau, delta, X(21)
    ! Calculated variables
    double precision, dimension(3, 3):: ar, ao
    double precision:: P, Z, w, mean_M, rho2, old_rho, B, C, cp
    double precision:: T_r, rho_r
    character(len=100):: compound
    integer:: i, j, io
    double precision:: eps = epsilon(T)

    ! Ideal Gas Parameters
    no_file = 'parameters/ideal/n'
    vo_file = 'parameters/ideal/v'
    ! Residual energy parameters
    Kr_file = 'parameters/residual/Kpol-Kexp'
    nr_file = 'parameters/residual/n'
    dr_file = 'parameters/residual/d'
    tr_file = 'parameters/residual/t'
    cr_file = 'parameters/residual/c'
    ! Departure Function parameters
    Kij_file = 'parameters/departure/Kpolij-Kexpij'
    betaij_file = 'parameters/departure/betaij'
    dij_file = 'parameters/departure/dij'
    epsij_file = 'parameters/departure/epsij'
    etaij_file = 'parameters/departure/etaij'
    gammij_file = 'parameters/departure/gammij'
    nij_file = 'parameters/departure/nij'
    tij_file = 'parameters/departure/tij'
    ! Reducing function parameters
    reducing_file = 'parameters/reducing_params'
    Fij_file = 'parameters/departure/Fij'
    ! Critical parameters
    critical_file = 'parameters/critical'

    ! Constantes
    R = 8.314472

    call get_parameters(no_file, vo_file, no, vo, &
                        Kr_file, nr_file, dr_file, tr_file, cr_file, &
                        Kpol, Kexp, nr, dr, tr, cr, &
                        Kij_file, betaij_file, dij_file, epsij_file, &
                        etaij_file, gammij_file, nij_file, tij_file, &
                        Fij_file, &
                        Kpolij, Kexpij, betaij, dij, epsij, etaij, &
                        gammij, nij, tij, Fij, &
                        reducing_file, Bv, Gv, Bt, Gt, &
                        critical_file, rho_c, T_c, M)

    ! Propiedades de operación
    ! ------------------------
    !T = 185.0
    T = 276.15
    X = 0
    open (1, file='concentrations')
    io = 0
    i = 0
    do while (io == 0)
        read (1, *, iostat=io) i, compound, X(i)
    end do

    mean_M = 0
    do i = 1, size(X)
        mean_M = mean_M + X(i)*M(i)
    end do

    X = 0
    X(1) = 1

    call reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)

    do j = 20, 50, 2
        !T = float(i)
        rho = float(j)/100

        tau = T_r/T
        delta = rho*rho_r

        call residual_term(X, delta, tau, &
                           Kpol, Kexp, nr, dr, tr, cr, &
                           Fij, Kpolij, Kexpij, dij, nij, tij, betaij, epsij, etaij, gammij, &
                           ar)
        call ideal_term(X, rho, T, rho_c, T_c, rho_r, T_r, &
                        no, vo, &
                        ao)
        call zeta(delta, ar(2, 1), Z)

        P = Z*rho*T*R*1000
        print *, T, rho, P/1e6
    end do

End Program gerg

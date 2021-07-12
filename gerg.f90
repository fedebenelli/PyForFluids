! This is file : gerg
! Author= Federico Benelli
! Started at: 01/07/2021
! Last Modified: lun 12 jul 2021 10:31:53
!

! General use subrutines
! -----------------------
Subroutine get_pure_params(f, compound, M, Tc, Pc, Dc, Zc, Kpol, Kexp, &
                           no, vo, &
                           no_r, do_r, to_r, co_r)
    ! Get the parameters of compound from file f
    integer:: i, j
    double precision, intent(inout):: M, Tc, Pc, Dc, Zc
    integer, intent(inout):: Kpol, Kexp
    double precision, dimension(7), intent(inout):: no, vo
    double precision, dimension(24), intent(inout):: no_r, do_r, to_r, co_r
    character(len=20):: label, compound, f

    open (1, file=f)
    do i = 1, 8
        select case (i)
        case (1)
            read (1, *) label, compound
        case (2)
            read (1, *) label, M, Tc, Pc, Dc, Zc, Kpol, Kexp
        case (3)
            read (1, *) label, (no(j), j=1, 7)
        case (4)
            read (1, *) label, (vo(j), j=4, 7)
        case (5)
            read (1, *) label, (no_r(j), j=1, Kpol + Kexp)
        case (6)
            read (1, *) label, (do_r(j), j=1, Kpol + Kexp)
        case (7)
            read (1, *) label, (to_r(j), j=1, Kpol + Kexp)
        case (8)
            read (1, *) label, (co_r(j), j=Kpol + 1, Kpol + Kexp)
        end select
    end do
End Subroutine get_pure_params

Subroutine get_params(rho_c, T_c, M)
    double precision, dimension(21), intent(inout):: rho_c, T_c, M
    double precision:: rc, tc, mi
    character(len=20):: compound
    integer:: N = 21

    open (1, file='data/critical')
    do i = 1, N
        read (1, *) compound, rc, tc, mi
        rho_c(i) = rc
        T_c(i) = tc
        M(i) = mi
    end do

End Subroutine get_params

Subroutine reducing_params(Bv, Gv, Bt, Gt)
    double precision, dimension(21, 21), intent(inout):: Bv, Gv, Bt, Gt
    double precision:: Bvij, Gvij, Btij, Gtij
    integer:: i, j, n

    Bv = 0
    Gv = 0
    Bt = 0
    Gt = 0
    open (1, file='./data/reducing_params')
    do n = 1, 210
        read (1, *) i, j, Bvij, Gvij, Btij, Gtij
        Bv(i, j) = Bvij
        Gv(i, j) = Gvij
        Bt(i, j) = Btij
        Gt(i, j) = Gtij
    end do

End Subroutine reducing_params

! Thermodynamic Properties
! ---------------------

Subroutine zeta(delta, dArdD, z)
    double precision:: delta, dArdD
    double precision, intent(inout):: z
    z = 1 + delta*dArdD
End Subroutine zeta

Subroutine sound_speed(R, T, M, delta, tau, Ar, Ao, w)
    double precision, intent(in):: R, T, M, delta, tau, Ar(3, 3), Ao(3, 3)
    double precision, intent(inout):: w

    w = 1 + 2*delta*Ar(2, 1) + delta**2*Ar(3, 1)
    w = w - (1 + delta*Ar(2, 1) - delta*tau*Ar(3, 3))**2 &
        /(tau**2*(Ao(3, 2) + Ar(3, 2)))
    w = sqrt(w*R*T/M)
End Subroutine sound_speed

Subroutine reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)
    ! REDUCING DENSITY AND TEMPERATURE
    double precision, dimension(21), intent(in):: X, rho_c, T_c
    double precision, dimension(21, 21), intent(in)::  Bv, Gv, Bt, Gt
    double precision, intent(inout):: rho_r, T_r
    integer:: N = 21, i, j

    rho_r = 0
    T_r = 0

    do i = 1, N
    if (X(i) > EPSILON(X(i))) then
        rho_r = rho_r + X(i)**2/rho_c(i)
        T_r = T_r + X(i)**2*T_c(i)
    end if
    end do

    do i = 1, N - 1
    do j = i + 1, N
    if (X(i) > EPSILON(X(i))) then
        if (X(j) > EPSILON(X(i))) then
            rho_r = rho_r + &
                    2*X(i)*X(j)*Bv(i, j)*Gv(i, j) &
                    *(X(i) + X(j))/(Bv(i, j)**2*X(i) + X(j)) &
                    *1/8*(rho_c(i)**(-1.0/3.0) + rho_c(j)**(-1.0/3.0))**3

            T_r = T_r + &
                  2*X(i)*X(j)*Bt(i, j)*Gt(i, j) &
                  *(X(i) + X(j))/(Bt(i, j)**2*X(i) + X(j)) &
                  *sqrt((T_c(i)*T_c(j)))
        end if
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
    ! --------------+-----------------+-----------------|
    ! aoio          |     0           |    0            |
    ! d(aoio)/dd    |     d(aor)/dt   |    0            |
    ! d2(aoio)/dd2  |     d2(aor)/dt2 |    d2(aor)/dddt |
    ! --------------------------------------------------|

    double precision, intent(in):: rho, T, rho_c, T_c
    double precision, intent(in), dimension(7):: n, v
    double precision, intent(inout):: aoio(3, 3)
    double precision:: delta, tau
    double precision:: r

    r = 8.314510/8.314472

    aoio = 0
    delta = rho/rho_c
    tau = T_c/T

    ! Ideal Gas Helmholtz Energy
    aoio(1, 1) = log(delta) + r*(n(1) + n(2)*tau + n(3)*log(tau))

    do k = 4, 6
        aoio(1, 1) = aoio(1, 1) + n(k)*log(abs(sinh(v(k)*tau)))
    end do
    do k = 5, 7
        aoio(1, 1) = aoio(1, 1) - n(k)*log(cosh(v(k)*tau))
    end do

    ! First derivative with reduced density
    aoio(2, 1) = delta

    ! Second derivative with reduced density
    aoio(3, 1) = -delta**2.0

    ! First Derivative with reduced temperature
    aoio(2, 2) = n(2) + n(3)/tau
    do k = 4, 6
        aoio(2, 2) = aoio(2, 2) + n(k)*v(k)/tanh(v(k)*tau)
    end do
    do k = 5, 7
        aoio(2, 2) = aoio(2, 2) - n(k)*v(k)*tanh(v(k)*tau)
    end do
    aoio(2, 2) = r*aoio(2, 2)

    ! Second Derivative with reduced temperature
    aoio(3, 2) = -n(3)*(1/tau)**2
    do k = 4, 6
        aoio(3, 2) = aoio(3, 2) - n(k)*(v(k)**2)/(sinh(v(k)*tau)**2)
    end do
    do k = 5, 7
        aoio(3, 2) = aoio(3, 2) - n(k)*(v(k)**2)/(cosh(v(k)*tau)**2)
    end do
    aoio(3, 2) = r*aoio(3, 2)

    ! Second Derivative with reduced density and temperature
    aoio(3, 3) = 0

End Subroutine a_oio

Subroutine a_oir(delta, tau, Kpol, Kexp, n, d, t, c, aoir)
    ! RESIDUAL ENERGY
    ! ---------------
    ! Caculate the pure compound Residual Helmholtz Energy and its derivatives
    ! aor is structured like:
    ! -------------+----------------+------------------|
    ! aor          |     0           |    0            |
    ! d(aor)/dd    |     d(aor)/dt   |    0            |
    ! d2(aor)/dd2  |     d2(aor)/dt2 |    d2(aor)/dddt |
    ! -------------------------------------------------|
    !

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
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(1, 1) = aoir(1, 1) + &
                     n(k)*delta**d(k) &
                     *tau**t(k) &
                     *exp(-delta**c(k))
    end do

    ! First Derivative with reduced density
    do k = 1, Kpol
        aoir(2, 1) = aoir(2, 1) + &
                     n(k)*d(k) &
                     *delta**(d(k) - 1) &
                     *tau**t(k)
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(2, 1) = aoir(2, 1) + &
                     n(k) &
                     *delta**(d(k) - 1) &
                     *(d(k) - c(k)*delta**c(k)) &
                     *tau**t(k) &
                     *exp(-delta**c(k))
    end do

    ! First Derivative with reduced temperature
    do k = 1, Kpol
        aoir(2, 2) = aoir(2, 2) + &
                     n(k)*t(k)*delta**d(k) &
                     *tau**(t(k) - 1)
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(2, 2) = aoir(2, 2) + &
                     n(k)*t(k)*delta**d(k) &
                     *tau**(t(k) - 1) &
                     *exp(-delta**c(k))
    end do

    ! Second Derivative with reduced density
    do k = 1, Kpol
        aoir(3, 1) = aoir(3, 1) + &
                     n(k)*d(k)*(d(k) - 1) &
                     *delta**(d(k) - 2) &
                     *tau**t(k)
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(3, 1) = aoir(3, 1) + &
                     n(k)*delta**(d(k) - 2) &
                     *( &
                     (d(k) - c(k)*delta**c(k) &
                      *(d(k) - 1 - c(k)*delta**c(k)) &
                      - c(k)**2*delta**c(k)) &
                     ) &
                     *tau**t(k)*exp(-delta**c(k))
    end do

    ! Second Derivative with reduced temperature
    do k = 1, Kpol
        aoir(3, 2) = aoir(3, 2) + &
                     n(k)*t(k)*(t(k) - 1) &
                     *delta**d(k)*tau**(t(k) - 2)
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(3, 2) = aoir(3, 2) + &
                     n(k)*t(k)*(t(k) - 1) &
                     *delta**d(k)*tau**(t(k) - 2) &
                     *exp(-delta**c(k))
    end do

    ! Second Derivative with reduced temperature and density
    do k = 1, Kpol
        aoir(3, 3) = aoir(3, 3) + &
                     n(k)*d(k)*t(k) &
                     *delta**(d(k) - 1)*tau**(t(k) - 1)
    end do
    do k = Kpol + 1, Kpol + Kexp
        aoir(3, 3) = aoir(3, 3) + &
                     n(k)*t(k)*delta**(d(k) - 1) &
                     *(d(k) - c(k)*delta**c(k)) &
                     *tau**(t(k) - 1)*exp(-delta**c(k))
    end do

End Subroutine a_oir

Program gerg
    Implicit None

    ! Constants
    double precision:: R
    double precision, dimension(21, 21):: Bv, Gv, Bt, Gt

    ! Variables
    double precision:: T, rho, tau, delta, X(21)

    ! Pure compound properties
    double precision, dimension(24):: no_r, do_r, to_r, co_r
    double precision, dimension(21):: rho_c, T_c, M
    double precision, dimension(7):: no, vo
    integer:: Kpol, Kexp, i
    character(len=20):: compound
    double precision:: Mi, Tci, Pci, Dci, Zci

    ! Calculated variables
    double precision:: aoir(3, 3), aoio(3, 3)
    double precision:: P, Z, w
    double precision:: T_r, rho_r

    ! data file
    character(len=100):: f = 'data/params'

    ! Constantes
    R = 8.314472

    call get_pure_params(f, compound, Mi, Tci, Pci, Dci, Zci, Kpol, Kexp, &
                         no, vo, &
                         no_r, do_r, to_r, co_r)
    call get_params(rho_c, T_c, M)
    call reducing_params(Bv, Gv, Bt, Gt)

    ! Propiedades de operación
    X = 0
    X(1) = 1
    call reducing_funcs(X, Bv, Gv, Bt, Gt, rho_c, rho_r, T_c, T_r)

    T = 185
    tau = T_r/T

    print *, "Isotherm of: ", compound
    print *, "T: ", T, "K"
    print *, "------------------------"

    do i = 1, 3000, 1
        rho = float(i)/100
        delta = rho*rho_r

        ! Calculate residual and ideal gas helmholtz energy
        call a_oir(delta, tau, Kpol, Kexp, no_r, do_r, to_r, co_r, aoir)
        call a_oio(rho, T, Dci, Tci, no, vo, aoio)

        call zeta(delta, aoir(2, 1), Z)
        P = delta/tau*Z/Zci*Pci

        call sound_speed(R, T, Mi, delta, tau, aoir, aoio, w)

        print *, compound, rho, P, w
    end do

End Program gerg

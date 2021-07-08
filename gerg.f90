! This is file : gerg
! Author= Federico Benelli
! Started at: 01/07/2021
! Last Modified: jue 08 jul 2021 13:05:45
!

! General use subrutines
! -----------------------
Subroutine get_params(f, compound, Tc, Pc, Dc, Zc, Kpol, Kexp, &
                      no_r, do_r, to_r, co_r)
    ! Get the parameters of compound from file f
    integer:: i, j
    double precision, intent(inout):: Tc, Pc, Dc, Zc
    integer, intent(inout):: Kpol, Kexp
    double precision, dimension(24), intent(inout):: no_r, do_r, to_r, co_r
    character(len=20):: label, compound, f

    open (1, file=f)
    do i = 1, 6
        select case (i)
        case (1)
            read (1, *) label, compound
        case (2)
            read (1, *) label, Tc, Pc, Dc, Zc, Kpol, Kexp
        case (3)
            read (1, *) label, (no_r(j), j=1, Kpol + Kexp)
        case (4)
            read (1, *) label, (do_r(j), j=1, Kpol + Kexp)
        case (5)
            read (1, *) label, (to_r(j), j=1, Kpol + Kexp)
        case (6)
            read (1, *) label, (co_r(j), j=Kpol + 1, Kpol + Kexp)
        end select
    end do
End Subroutine get_params

Subroutine zeta(delta, dArdD, z)
    double precision:: delta, dArdD
    double precision, intent(inout):: z
    z = 1 + delta*dArdD
End Subroutine zeta

! Pure Compound Helmholtz Energy (and derivatives) Calculations
! -----------------------------------------------
Subroutine a_oio(rho, T, rho_c, T_c, n, v, aoio)
    ! IDEAL GAS ENERGY
    ! ----------------
    ! Calculate the pure compound ideal Helmholtz Enegry and its derivatives

    double precision, intent(in):: rho, T, rho_c, T_c
    double precision, intent(in), dimension(7):: n, v
    double precision, intent(inout):: aoio(3, 3)
    double precision:: delta, tau
    double precision:: r

    ! <TODO> Set real value
    r = 1

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

    ! Variables
    double precision:: T, dens, tau, delta

    ! Pure compound properties
    double precision, dimension(24):: no_r, do_r, to_r, co_r
    double precision:: Tc, Pc, Dc, Zc
    integer:: Kpol, Kexp, i
    character(len=20):: compound

    ! Calculated variables
    double precision:: aoir(3, 3)
    double precision:: P, Z

    ! data file
    character(len=100):: f = 'data/params'

    ! Constantes
    R = 8.31

    call get_params(f, compound, Tc, Pc, Dc, Zc, Kpol, Kexp, &
                    no_r, do_r, to_r, co_r)

    ! Propiedades de operación
    T = 185.0
    tau = Tc/T

    do i = 0, 300, 5
        delta = i/Dc/10
        dens = float(i)/10
        ! Calculate residual helmholtz energy
        call a_oir(delta, tau, Kpol, Kexp, no_r, do_r, to_r, co_r, aoir)

        ! Calculate the compresibility factor
        call zeta(delta, aoir(2, 1), Z)

        P = delta/tau*Z/Zc*Pc

        print *, compound, dens, P
    end do

    delta = 1
    do i = 1, 5000
        T = i/10
        tau = Tc/T
        call a_oir(delta, tau, Kpol, Kexp, no_r, do_r, to_r, co_r, aoir)
        call zeta(delta, aoir(2, 1), Z)

        P = delta/tau*Z/Zc*Pc

        print *, "TP", T, P
    end do

End Program gerg

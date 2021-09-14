! Author= Federico Benelli
! Gerg-2008 Equation
! Started at: 01/07/2021
! Last Modified: mié 01 sep 2021 15:48:41
!

! -------------------
! GERG 2008 Functions
! -------------------
Subroutine reducing_funcs(X, rho_r, T_r)
    ! REDUCING DENSITY AND TEMPERATURE
    ! input:
    !  - X      (array): molar fractions
    !  - Bv     (array): beta parameters for reducing density
    !  - Gv     (array): gamma parameters for reducing density
    !  - Bt     (array): beta parameters for reducing temperature
    !  - Gt     (array): gamma parameters for reducing temperature
    !  - rho_c  (array): citical densities
    !  - T_c    (array): citical temperatures
    ! output:
    !  - rho_r  (float): Reducing density (1/rhor)
    !  - T_r    (float): reducing temperature (Tr)
    use parameters
    real*8, dimension(21), intent(in):: X
    real*8, intent(out):: rho_r, T_r
    integer:: i, j
    call get_params()

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
                *(X(i) + X(j))/(Bv(i, j)**2*X(i) + X(j)) &
                *1.d0/8.d0*(rho_c(i)**(-1.d0/3.d0) + rho_c(j)**(-1.d0/3.d0))**3

        T_r = T_r + &
              2.d0*X(i)*X(j)*Bt(i, j)*Gt(i, j) &
              *(X(i) + X(j))/(Bt(i, j)**2*X(i) + X(j)) &
              *sqrt((T_c(i)*T_c(j)))
    end if
    end do
    end if
    end do
End Subroutine reducing_funcs

Subroutine  V_calc(X, P, T, rho)
        ! Based on: https://github.com/usnistgov/AGA8/issues/12
        use parameters
        real*8, intent(in):: P, T, X(21)
        real*8, intent(out):: rho
        real*8:: vo, vi, dpdr, ar(3,3), delta, tau, rho_r, T_r, Pcalc
        integer:: it = 20
        
        rho = P/(R*T)
        rho = rho/1d3
        vo = -log(rho)
        call reducing_funcs(X, rho_r, T_r)
        do i=1,it
                delta = rho*rho_r
                tau = T_r/T
                call residual_term(X, delta, tau, ar)

                Pcalc = (1.d0 + delta*ar(2,1))
                dpdr = R*T*(1.d0 + 2.d0*delta*ar(2,1) + delta**2*ar(3,1))
                
                vi = vo + (Pcalc/(rho*1d3))/dpdr * (log(Pcalc) - log(P))

                if (abs(vi-vo) < Epsilon(vi)) then
                        rho = exp(-vi)
                        rho = rho/1d3
                return
                else
                        vo = vi
                        rho = exp(-vo)
                end if
                rho = rho / 1d3

                print *, "[VCALC]", i, rho
        end do
End Subroutine V_calc 


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
    real*8, intent(in):: rho, T, rho_c, T_c
    real*8, intent(in), dimension(7):: n, v
    real*8, intent(out):: aoio(3, 3)
    real*8:: r, Tr, Dr, eps = 1d-10

    r = 8.314510d0/8.314472d0
    
    aoio = 0.d0
    Tr = T_c/T
    Dr = rho/rho_c

    do k = 4, 7
        if (v(k) > eps) then
            if (k == 4) then
                aoio(1, 1) = aoio(1, 1) + r*(n(k)*log(abs(dsinh(v(k)*Tr))))
                aoio(2, 2) = aoio(2, 2) + r*(n(k)*v(k)/dtanh(v(k)*Tr))
                aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/dsinh(v(k)*Tr)**2)
            else if (k == 5 .or. k == 6) then
                aoio(1, 1) = aoio(1, 1) + r*(n(k)*log(abs(dsinh(v(k)*Tr))))
                aoio(1, 1) = aoio(1, 1) - r*(n(k)*log(dcosh(v(k)*Tr)))

                aoio(2, 2) = aoio(2, 2) + r*(n(k)*v(k)/dtanh(v(k)*Tr))
                aoio(2, 2) = aoio(2, 2) - r*(n(k)*v(k)*dtanh(v(k)*Tr))

                aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/dsinh(v(k)*Tr)**2)
                aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/dcosh(v(k)*Tr)**2)
            else
                aoio(1, 1) = aoio(1, 1) - r*(n(k)*log(dcosh(v(k)*Tr)))
                aoio(2, 2) = aoio(2, 2) - r*(n(k)*v(k)*dtanh(v(k)*Tr))
                aoio(3, 2) = aoio(3, 2) - r*(n(k)*v(k)**2/dcosh(v(k)*Tr)**2)
            end if
        end if
    end do

    aoio(1, 1) = aoio(1, 1) + log(Dr) + r*(n(1) + n(2)*Tr + n(3)*log(Tr))
    aoio(2, 1) = 1.d0/Dr
    aoio(3, 1) = -1.d0/Dr**2
    aoio(2, 2) = aoio(2, 2) + r*(n(2) + n(3)/Tr)
    aoio(3, 2) = aoio(3, 2) + r*(-n(3)*(T/T_c)**2)
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
                     n(k)*t(k)*(t(k) - 1.d0) &
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
                     n(k)*t(k)*(t(k) - 1.d0) &
                     *delta**d(k)*tau**(t(k) - 2.d0) &
                     *exp(-delta**c(k))
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

    aijr = 0.d0
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

Subroutine ideal_term(X, rho, T, rho_r, T_r, ao)
    use parameters
    real*8, intent(in):: X(21), rho, T, rho_r, T_r
    real*8, intent(out):: ao(3,3)
    real*8::  aoio(3,3), xi

    call get_params()

    ao = 0.d0
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

Subroutine residual_term(X, delta, tau, ar)
    use parameters
    real*8, intent(in):: delta, tau, X(21)
    real*8, dimension(3, 3), intent(out):: ar
    real*8, dimension(3, 3):: aoir, aijr
    integer:: i,j

    call get_params()

    ar = 0.0

    do i = 1, N
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

Program gerg
    use parameters
    use thermo_props
    Implicit None
    ! Variables
    real*8:: T, rho, tau, delta, X(21)
    ! Calculated variables
    real*8, dimension(3, 3):: ar, ao
    real*8:: T_r, rho_r, P, Z, w, cp, cv, mean_M
    ! bintest
    real*8:: test_P, test_cv, test_cp, test_w
    integer:: io
    ! Input args
    character(len=100):: arg
    integer:: i,j,k

    ! Get the R constant, compound properties and equations parameters
    call get_params()

    ! Propiedades de operación
    ! ------------------------
    call getarg(1, arg)
    rho = 1.d0
    T = 1.d0
    X = 0.d0

    select case (arg)
    case ('bintest')
        open (1, file='GERG_test/binary/test_data')
        io = 0
        ! Reads until the file with test data ends
        do while (io == 0)
            X = 0
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
            call reducing_funcs(X, rho_r, T_r)
            delta = rho*rho_r
            tau = T_r/T
            call residual_term(X, delta, tau, ar)
            call ideal_term(X, rho, T, rho_r, T_r, ao)

            ! Calculate Properties
            call zeta(delta, ar, Z)
            call isobaric_heat(delta, tau, R, ao, ar, cp)
            call isochoric_heat(tau, R, Ao, Ar, cv)
            call sound_speed(delta, tau, R, T, mean_M, Ao, Ar, w)

            ! adim * [mol/L] * [K] * [J/(mol*K)] * [L/m3] * [MPa/Pa]
            P = Z*rho*T*R*1000.d0*1.d-6

            print *, i, j, T, rho, P, cv, cp, w, test_P, test_cv, test_cp, test_w
            write(0,*) i,j,delta, tau, ar(3,2), ao(3,2), cv, test_cv

        end do
        close (1)
    case ('VCalc')
            rho = 12.7982864
            T = 400.0
            P = 50000.d3

            x(1 )=0.77824d0       !Methane
            x(2 )=0.02d0          !Nitrogen
            x(3 )=0.06d0          !CO2
            x(4 )=0.08d0          !Ethane
            x(5 )=0.03d0          !Propane
            x(6 )=0.003d0         !Butane
            x(7 )=0.0015d0        !Isobutane
            x(8 )=0.00165d0       !Pentane
            x(9 )=0.0005d0        !Isopentane
            x(10)=0.00215d0       !Hexane
            x(11)=0.00088d0       !Heptane
            x(12)=0.00024d0       !Octane
            x(13)=0.00015d0       !Nonane
            x(14)=0.00009d0       !Decane
            x(15)=0.004d0         !Hydrogen
            x(16)=0.005d0         !Oxygen
            x(17)=0.002d0         !CO
            x(18)=0.0001d0        !Water
            x(19)=0.0025d0        !H2S
            x(20)=0.007d0         !Helium
            x(21)=0.001d0         !Argon

            call reducing_funcs(X, rho_r, T_r)
            
            delta = rho*rho_r
            tau = T_r/T

            call residual_term(X, delta, tau, ar)
            call zeta(delta, ar, z)
            P = Z*rho*T*R*1000.d0
            print *, rho, P

            call V_calc(X, P, T, rho)
            print *, rho, P

    end select
End Program gerg

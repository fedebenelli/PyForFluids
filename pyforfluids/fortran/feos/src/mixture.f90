module multicomponent_eos
    !use constants
    !use mixing_rules
    !use cubic_eos

    implicit none

    private

    !public :: mixture

    !type :: mixture(n)
    !   integer, len :: n
    !   real(wp) :: P
    !   real(wp) :: V
    !   real(wp) :: T
    !   class(CubicEoS), allocatable :: compounds(n)
    !   real(wp), allocatable :: concentrations(:)
    !   type(ConstantBIP) :: mixing_rule(n)
    !   class(CubicEoS), allocatable :: as_a_fluid !! As a fluid representation of the mixture
    !end type mixture
end module multicomponent_eos

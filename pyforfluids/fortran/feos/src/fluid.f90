module feos
    use constants
    use properties
    use cubic_eos
    use mixing_rules

    implicit none
    private

    public :: CubicFluid
    public :: PengRobinson, PR
    public :: SoaveRedlichKwong, SRK
    public :: scalar_property, binary_property
    public :: wp, R
    public :: ClassicVdW, aij_classic, bij_classic
    public :: CubicEoS

    type :: CubicFluid
        class(CubicEoS), allocatable :: components(:)
        class(ClassicVdW), allocatable :: mixing_rule
    contains
        procedure :: residual_helmholtz
        procedure :: set_to
        procedure :: copy
    end type CubicFluid

contains

    pure subroutine set_to(self, pressure, volume,  temperature)
        class(CubicFluid), intent(in out) :: self
        real(wp), optional, intent(in) :: pressure
        real(wp), optional, intent(in) :: volume
        real(wp), optional, intent(in) :: temperature

        if (present(volume)) then
            self%components%v = volume
        end if

        if (present(temperature)) then
            self%components%t = temperature
        end if
    end subroutine

    pure function copy(self) result(new)
        class(CubicFluid), intent(in) :: self
        type(CubicFluid), allocatable :: new

        allocate(CubicFluid :: new)
        new%components = self%components
        new%mixing_rule = self%mixing_rule
    end function
    
    elemental function residual_helmholtz(self, v, t) result(ar)
        class(CubicFluid), intent(in) :: self
        real(wp), intent(in) :: v
        real(wp), intent(in) :: t
        type(scalar_property) :: ar

        real(wp) :: moles(size(self%components))
        class(CubicEoS), allocatable :: compounds(size(self%components))
        class(CubicFluid), allocatable :: fluid_tmp
        integer :: n

        type(scalar_property) :: d, b, d1

        n = size(self%components)
        ar = null_scalar_property(n)
        fluid_tmp = self%copy()
       
        call fluid_tmp%set_to(volume=v, temperature=t)

        compounds = fluid_tmp%components
        moles = fluid_tmp%components%moles
        call fluid_tmp%mixing_rule%mix(compounds, d, b, d1)

        ar = cubic_residual_helmholtz(n, moles, v, t, d, d1, b)

    end function residual_helmholtz

end module feos

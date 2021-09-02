install:
	gfortran  parameters.f90 modules.f90 gerg.f90 -Wall -Wextra -fcheck=all -o gerg
debug:
	gfortran modules.f90 parameters.f90 gerg.f90 -g -O0 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal -o gerg

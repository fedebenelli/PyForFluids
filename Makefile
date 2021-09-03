install:
	gfortran  parameters.f95 thermoprops.f95 gerg.f95 -Wall -Wextra -fcheck=all -o gerg
	f2py -c parameters.f95 thermoprops.f95 gerg.f95 -m pyforfluids
debug:
	gfortran parameters.f95 thermoprops.f95 gerg.f95 -g -O0 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal -o gerg

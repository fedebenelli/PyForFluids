install:
	gfortran modules.f90 parameters.f90 gerg.f90 -Wall -Wextra -fcheck=all -o gerg
debug:
	gfortran modules.f90 parameters.f90 gerg.f90 -g -O0 -Wall -Wextra -fcheck=all -o gerg
testing:
	gfortran modules.f90 parameters.f90 gerg.f90 -O0 -Wall -Wextra -fcheck=all -o gerg
	./run_gerg bintest


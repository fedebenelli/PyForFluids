FORTRAN_SOURCES="src/fortran"
FILES=${FORTRAN_SOURCES}/parameters.f95 ${FORTRAN_SOURCES}/thermoprops.f95 ${FORTRAN_SOURCES}/gerg.f95
COMPILE_COMMAND=gfortran ${FILES} -Wall -Wextra -fcheck=all
DEBUG_EXECUTABLE=gfortran ${FILES} -g -O0 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal

clean:
	rm *.mod ${FORTRAN_SOURCES}/*.mod

install:
	${COMPILE_COMMAND} -o gerg

debug:
	${DEBUG_EXECUTABLE} -o gerg_debug

py_module:
	f2py -c ${FORTRAN_SOURCES}/*f95 -m pyforfluids


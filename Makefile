COMPILER=gfortran
FORTRAN_SOURCES="src/fortran"
COMPILE_OPTS=-Wall -Wextra -fcheck=all
DEBUG_OPTS=-g -O0 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal
FILES=${FORTRAN_SOURCES}/parameters.f95 ${FORTRAN_SOURCES}/thermoprops.f95 ${FORTRAN_SOURCES}/gerg.f95

clean:
	@rm *.mod 

install:
	${COMPILER} ${COMPILE_OPTS} ${FILES} -o gerg

debug:
	${COMPILER} ${COMPILE_OPTS} ${DEBUG_OPTS} ${FILES} -o gerg_debug

py_module:
	f2py -c ${FORTRAN_SOURCES}/*f95 -m gerg2008

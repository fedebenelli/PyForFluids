cmake_minimum_required(VERSION 3.18) # Needed to avoid requiring embedded Python libs too

project(pyforfluids
  VERSION 1.0
  DESCRIPTION "PyForFluids Fortran Modules"
  LANGUAGES C Fortran
)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Grab Python, 3.8 or newer
find_package(Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module NumPy)
find_package(PythonExtensions REQUIRED)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

# =============================================================================
# fgerg2008 module
# =============================================================================
message(STATUS "Compiling fortran gerg2008...")
set(f2py_module_name "fgerg2008")
set(f2py_module_dir "pyforfluids/fortran")
set(
    fortran_src_file 
    "${CMAKE_SOURCE_DIR}/pyforfluids/fortran/parameters.f95"
    "${CMAKE_SOURCE_DIR}/pyforfluids/fortran/gerg.f95"
)
set(
    generated_module_file 
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX}
)

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
)

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND
   ${Python_EXECUTABLE}  -m "numpy.f2py"
    -m ${f2py_module_name}
    -c
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
)
if (SKBUILD)
    install(FILES ${generated_module_file} DESTINATION ${f2py_module_dir})
else()
  install(FILES ${generated_module_file} DESTINATION ${f2py_module_dir})
endif()
# =============================================================================

# =============================================================================
# fthermo_props module
# =============================================================================
message(STATUS "Compiling fortran thermoprops...")
set(f2py_module_name "fthermo_props")
set(f2py_module_dir "pyforfluids/fortran")
set(
    fortran_src_file 
    "${CMAKE_SOURCE_DIR}/pyforfluids/fortran/parameters.f95"
    "${CMAKE_SOURCE_DIR}/pyforfluids/fortran/thermoprops.f95"
)
set(
    generated_module_file 
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX}
)

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
)

add_custom_command(
  OUTPUT ${generated_module_file}
  COMMAND
   ${Python_EXECUTABLE}  -m "numpy.f2py"
    -m ${f2py_module_name}
    -c
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
)
if (SKBUILD)
    install(FILES ${generated_module_file} DESTINATION ${f2py_module_dir})
else()
  install(FILES ${generated_module_file} DESTINATION ${f2py_module_dir})
endif()
# =============================================================================


# =============================================================================
# pyfeos module
# =============================================================================
set(FEOSDIR "pyforfluids/fortran/feos")
set(FORTRAN_DIR "pyforfluids/fortran")
message(STATUS "Compiling fortran feos library...")

add_library(feos STATIC
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/constants.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/properties.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/cubic/base.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/cubic/cubic_eos.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/mixing_rules.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/src/fluid.f90
  ${CMAKE_SOURCE_DIR}/${FEOSDIR}/interface.f90
)
set_property(TARGET feos PROPERTY POSITION_INDEPENDENT_CODE ON)

install(
     TARGETS feos
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
 )

 message(STATUS "Compiling fortran pyfeos module...")
 
 set(f2py_module_name "pyfeos")
 set(f2py_module_dir  "${FORTRAN_DIR}")
 set(
   fortran_src_file 
   "${CMAKE_SOURCE_DIR}/${FEOSDIR}/pyfeos.f90"
  )
 set(
     generated_module_file 
     ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX}
 )
 add_custom_target(
   ${f2py_module_name} 
   ALL DEPENDS ${generated_module_file}
 )

 add_custom_command(
   OUTPUT ${generated_module_file}
   COMMAND 
   ${Python_EXECUTABLE}  -m "numpy.f2py"
   -L${CMAKE_CURRENT_BINARY_DIR} -lfeos
   -c -m ${f2py_module_name}
   ${fortran_src_file} only: ar
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   DEPENDS ${fortran_src_file} feos
 )

 install(FILES ${generated_module_file} DESTINATION ${f2py_module_dir})
# =============================================================================

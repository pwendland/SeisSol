
# Library settings:
option(HDF5 "Use HDF5 library for data output" OFF)
option(NETCDF "Use netcdf library for mesh input" OFF)
option(METIS "Use metis for partitioning" ON)
option(MPI "Use MPI parallelization" ON)


# Other settings
set(SEISSOL_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build")
set(BUILD_OPTIONS "Release" "Debug")
set_property(CACHE SEISSOL_BUILD_TYPE PROPERTY STRINGS ${BUILD_OPTIONS})


set(ORDER 6 CACHE INT "Convergence order")
set(ORDER_OPTIONS 2 3 4 5 6 7 8)
set_property(CACHE ORDER PROPERTY STRINGS ${ORDER_OPTIONS})


set(NUMBER_OF_MECHANISMS 0 CACHE INT "Number of mechanisms")


set(EQUATIONS "elastic" CACHE STRING "Equation set used")
set(EQUATIONS_OPTIONS elastic viscoelastic viscoelastic2)
set_property(CACHE EQUATIONS PROPERTY STRINGS ${EQUATIONS_OPTIONS})


set(ARCH "hsw" CACHE STRING "Type of the target architecture")
set(ARCH_OPTIONS noarch wsm snb hsw knc knl skx)
set(ARCH_ALIGNMENT   16  16  32  32  64  64  64)  # size of a vector register in bytes of a given architecture
set_property(CACHE ARCH PROPERTY STRINGS ${ARCH_OPTIONS})


set(PRECISION "double" CACHE STRING "type of floating point representation")
set(PRECISION_OPTIONS double float)
set_property(CACHE PRECISION PROPERTY STRINGS ${PRECISION_OPTIONS})


set(DYNAMIC_RUPTURE_METHOD "quadrature" CACHE STRING "Dynamic rupture method")
set(RUPTURE_OPTIONS quadrature cellaverage)
set_property(CACHE DYNAMIC_RUPTURE_METHOD PROPERTY STRINGS ${RUPTURE_OPTIONS})


option(PLASTICITY "Use plasticity")
set(PLASTICITY_METHOD "nb" CACHE STRING "Dynamic rupture method")
set(PLASTICITY_OPTIONS nb ip)
set_property(CACHE PLASTICITY_METHOD PROPERTY STRINGS ${PLASTICITY_OPTIONS})


set(NUMBER_OF_FUSED_SIMULATIONS 1 CACHE INT "A number of fused simulations")
set(MEMEORY_LAYOUT "${CMAKE_SOURCE_DIR}/../config/dense.xml" CACHE FILEPATH "A path to memory layout file")



# helpers for gpu porting
set(ACCELERATOR_TYPE "NONE" CACHE STRING "type of accelerator")
set(ACCELERATOR_TYPE_OPTIONS NONE GPU)
set_property(CACHE ACCELERATOR_TYPE PROPERTY STRINGS ${ACCELERATOR_TYPE_OPTIONS})


set(VALIDATION_MODE "COMPARE" CACHE STRING "type of accelerator")
set(VALIDATION_MODE_OPTIONS COMPARE RECORD)
set_property(CACHE VALIDATION_MODE PROPERTY STRINGS ${VALIDATION_MODE_OPTIONS})

#-------------------------------------------------------------------------------
# ------------------------------- ERROR CHECKING -------------------------------
#-------------------------------------------------------------------------------
function(check_parameter parameter_name value options)
    set(WRONG_PARAMETER -1)

    list(FIND options ${value} INDEX)
    
    if (${INDEX} EQUAL ${WRONG_PARAMETER})
        message(FATAL_ERROR "${parameter_name} is wrong. Specified \"${value}\". Allowed: ${options}")
    endif()

endfunction()



check_parameter("SEISSOL_BUILD_TYPE" ${SEISSOL_BUILD_TYPE} "${BUILD_OPTIONS}")
check_parameter("ORDER" ${ORDER} "${ORDER_OPTIONS}")
check_parameter("ARCH" ${ARCH} "${ARCH_OPTIONS}")
check_parameter("EQUATIONS" ${EQUATIONS} "${EQUATIONS_OPTIONS}")
check_parameter("PRECISION" ${PRECISION} "${PRECISION_OPTIONS}")
check_parameter("DYNAMIC_RUPTURE_METHOD" ${DYNAMIC_RUPTURE_METHOD} "${RUPTURE_OPTIONS}")
check_parameter("PLASTICITY_METHOD" ${PLASTICITY_METHOD} "${PLASTICITY_OPTIONS}")
check_parameter("ACCELERATOR_TYPE" ${ACCELERATOR_TYPE} "${ACCELERATOR_TYPE_OPTIONS}")
check_parameter("VALIDATION_MODE" ${VALIDATION_MODE} "${VALIDATION_MODE_OPTIONS}")


# check NUMBER_OF_MECHANISMS
# NOTE: NUMBER_OF_MECHANISMS must be zero if we use only elastic models
if ("${EQUATIONS}" STREQUAL "elastic" AND ${NUMBER_OF_MECHANISMS} GREATER 0)
    message(FATAL_ERROR "${EQUATIONS} does not support a NUMBER_OF_MECHANISMS > 0.")
endif()



if ("${PRECISION}" STREQUAL "double")
    set(REAL_SIZE_IN_BYTES 8)
elseif ("${PRECISION}" STREQUAL "float")
    set(REAL_SIZE_IN_BYTES 4)
endif()


# compute alignment
list(FIND ARCH_OPTIONS ${ARCH} INDEX)
list(GET ARCH_ALIGNMENT ${INDEX} ALIGNMENT)


# check NUMBER_OF_FUSED_SIMULATIONS
math(EXPR IS_ALIGNED_MULT_SIMULATIONS 
        "${NUMBER_OF_FUSED_SIMULATIONS} % (${ALIGNMENT} / ${REAL_SIZE_IN_BYTES})")

if (NOT ${NUMBER_OF_FUSED_SIMULATIONS} EQUAL 1 AND NOT ${IS_ALIGNED_MULT_SIMULATIONS} EQUAL 0)
    math(EXPR FACTOR "${ALIGNMENT} / ${REAL_SIZE_IN_BYTES}")
    message(FATAL_ERROR "a number of fused must be multiple of ${FACTOR}")
endif()

#-------------------------------------------------------------------------------
# ------------------------ COMPUTE ADDITIONAL PARAMETERS -----------------------
#-------------------------------------------------------------------------------
if ("${EQUATIONS}" STREQUAL "elastic")
    set(NUMBER_OF_QUANTITIES 6)

elseif("${EQUATIONS}" STREQUAL "viscoelastic" OR "${EQUATIONS}" STREQUAL "viscoelastic2")
    math(EXPR NUMBER_OF_QUANTITIES "9 + 6 * ${NUMBER_OF_MECHANISMS}")

endif()


# generate an internal representation of an architecture type which is used in seissol
string(SUBSTRING ${PRECISION} 0 1 PRECISION_PREFIX)
if (${PRECISION} STREQUAL "double")
    set(ARCH_STRING "d${ARCH}")
elseif(${PRECISION} STREQUAL "float")
    set(ARCH_STRING "s${ARCH}")
endif()



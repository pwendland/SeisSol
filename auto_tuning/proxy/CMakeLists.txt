cmake_minimum_required(VERSION 3.5)
project(seissol_proxy LANGUAGES CXX)



# process_users_input returs the following:
# variables: SEISSOL_BUILD_TYPE, ORDER, EQUATIONS, ARCH_STRING, 
#            PRECISION, PLASTICITY, NUMBER_OF_FUSED_SIMULATIONS
#            NUMBER_OF_QUANTITIES, REAL_SIZE_IN_BYTES, ALIGNMENT
# switches: HDF5, NETCDF, METIS, MPI
include(cmake/process_users_input.cmake)

# set hardware specific definition needed for seissol during compile time
set(HARDWARE_DEFINITIONS "ALIGNMENT=${ALIGNMENT}"
                         "REAL_SIZE=${REAL_SIZE_IN_BYTES}")


# get helper functions: get_specific_cpu_flags, 
include(cmake/auxiliary.cmake)


# sets CPU_ARCH_FLAGS at the current i.e parent, scope
get_specific_cpu_flags(${ARCH} ${CMAKE_CXX_COMPILER_ID})


# set c++ standarts
# NOTE: the commnads below doesn't support modern intel compiler
#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS ON)

# add definitions to all targets
find_package(OpenMP REQUIRED)
add_compile_options(-std=c++11 
                    #-pg  # <-- for gprof TODO: remove
                    )



# adjust common compiler flags
# TODO: check intel compiler
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CPU_FLAGS "-Werror=vla")
endif()


# create numerical specific difinitions needed during compile time
if ("${EQUATIONS}" STREQUAL "elastic" OR "${EQUATIONS}" STREQUAL "viscoelastic2")
    set(MODEL_DEFINITIONS "ENABLE_MATRIX_PREFETCH")
endif()

string(TOUPPER ${EQUATIONS} EQUATIONS_IN_UPPER_CASE)
list(APPEND MODEL_DEFINITIONS "GENERATEDKERNELS"
                              "CONVERGENCE_ORDER=${ORDER}"
                              "NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS}"
                              "NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}"
                              "EQUATIONS_${EQUATIONS_IN_UPPER_CASE}")


# adjust flags and definitions based on the building type
if (${SEISSOL_BUILD_TYPE} STREQUAL "Release")
    set(BUILD_FLAGS "-O3")
    set(BUILD_DEFINITIONS "NDEBUG")
elseif(${SEISSOL_BUILD_TYPE} STREQUAL "Debug")
    set(BUILD_FLAGS "-O0" "-g")
endif()

#----------------------------------- CODEGEN ------------------------------------

set(CODEGEN_OUTPUT_DIR ${PROJECT_BINARY_DIR}/src/generated_code)

set(GENERATED_SOURCES ${CODEGEN_OUTPUT_DIR}/tensor.cpp
                      ${CODEGEN_OUTPUT_DIR}/subroutine.cpp
                      ${CODEGEN_OUTPUT_DIR}/init.cpp
                      ${CODEGEN_OUTPUT_DIR}/kernel.cpp)

set(GENERATED_INCLUDES ${CODEGEN_OUTPUT_DIR}/subroutine.h
                       ${CODEGEN_OUTPUT_DIR}/tensor.h
                       ${CODEGEN_OUTPUT_DIR}/init.h
                       ${CODEGEN_OUTPUT_DIR}/kernel.h)

add_custom_command(
    OUTPUT
       ${GENERATED_SOURCES}
       ${GENERATED_INCLUDES}

    COMMAND 
        ${CMAKE_COMMAND} -E make_directory "${CODEGEN_OUTPUT_DIR}"

    COMMAND
        python ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
        "--equations" ${EQUATIONS}
        "--matricesDir" ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/matrices
        "--outputDir" ${CODEGEN_OUTPUT_DIR}
        "--arch" ${ARCH_STRING}
        "--order" ${ORDER}
        "--numberOfMechanisms" ${NUMBER_OF_MECHANISMS}
        "--memLayout" ${MEMEORY_LAYOUT}
        "--multipleSimulations" ${NUMBER_OF_FUSED_SIMULATIONS}
        "--dynamicRuptureMethod" ${DYNAMIC_RUPTURE_METHOD}
        "--PlasticityMethod" ${PLASTICITY_METHOD}
    
    WORKING_DIRECTORY 
        ${CMAKE_SOURCE_DIR}/generated_code
    
    DEPENDS
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/viscoelastic.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/DynamicRupture.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/Plasticity.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/multSim.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/aderdg.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/generate.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/Point.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/elastic.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/viscoelastic2.py
       ${CMAKE_CURRENT_SOURCE_DIR}/generated_code/SurfaceDisplacement.py
    
    COMMENT "RUNNING PYTHON CODEGEN")



add_library(codegen ${GENERATED_SOURCES})
target_compile_options(codegen PRIVATE ${CPU_FLAGS} 
                                       ${CPU_ARCH_FLAGS}
                                       ${BUILD_FLAGS} 
                                       "-m64"  # <-- for intel mkl
                                       ${OpenMP_CXX_FLAGS})

target_compile_definitions(codegen PRIVATE ${BUILD_DEFINITIONS}
                                           ${HARDWARE_DEFINITIONS}
                                           MKL_ILP64  # <-- for intel mkl 
                                           )

# add all necessary libs for intel mkl
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    target_link_libraries(codegen PRIVATE -L${MKLROOT}/lib/intel64 
                                          -lmkl_intel_ilp64 
                                          -lmkl_sequential 
                                          -lmkl_core 
                                          -lpthread  # TODO: disable for score-p
                                          -lm 
                                          -ldl)

elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_link_libraries(codegen PRIVATE -L${MKLROOT}/lib/intel64 
                                          -lmkl_intel_ilp64 
                                          -lmkl_sequential 
                                          -lmkl_core 
                                          -lpthread  # TODO: disable for score-p
                                          -lm 
                                          -ldl)
endif()

target_include_directories(codegen PRIVATE  ${CODEGEN_OUTPUT_DIR}
                                            ${MKLROOT}/include
                                            "../../submodules/yateto/include"
                                            "../../src")

# TODO: add other matrix-matrix product libs

#----------------------------------- PROXY -----------------------------------

# necessary_sources.cmake defines all necessary source files to build the project
include(src/necessary_sources.cmake)


add_executable(seissol_proxy ${SOURCES})
add_dependencies(seissol_proxy codegen)


# add properties to the proxy_main target
target_compile_options(seissol_proxy PRIVATE ${CPU_FLAGS} 
                                             ${CPU_ARCH_FLAGS}
                                             ${BUILD_FLAGS}
                                             ${OpenMP_CXX_FLAGS})


target_compile_definitions(seissol_proxy PRIVATE ${BUILD_DEFINITIONS}
                                                 ${HARDWARE_DEFINITIONS}
                                                 ${MODEL_DEFINITIONS})


target_include_directories(seissol_proxy PRIVATE "../../submodules"
                                                 "../../submodules/yateto/include"
                                                 "../../submodules/easi/include"
                                                 "../../src/Equations/${EQUATIONS}"
                                                 "../../src"
                                                 ${PROJECT_BINARY_DIR}/src  # <- codegen
                                                 )

# assemble the entire applications
set_target_properties(seissol_proxy PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
target_link_libraries(seissol_proxy PRIVATE codegen ${OpenMP_CXX_LIBRARIES}
                                                    #-pg  # <-- for grof TODO: remove 
                                                    )


#------------------------------- CUSTOM BLAS ---------------------------------
message(STATUS "Hardware Definitions for CUDA: ${HARDWARE_DEFINITIONS}")
set(ACCELERATOR_LIB "customBLAS/cuda")
target_include_directories(seissol_proxy PRIVATE ${ACCELERATOR_LIB}/src)

add_subdirectory(${ACCELERATOR_LIB}) # compiles 'custom_blas' library
target_link_libraries(seissol_proxy PRIVATE custom_blas)



# ----------------------------------- HDF5 -------------------------------------
find_package(HDF5 REQUIRED)

target_compile_definitions(seissol_proxy PRIVATE ${HDF5_DEFINITIONS})
target_include_directories(seissol_proxy PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(seissol_proxy PRIVATE ${HDF5_LIBRARIES})

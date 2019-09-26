
# specify all necessary sorce files needed to build proxy
set(SEISSOL_SOURCES "src/seissol_src/Initializer/GlobalData.cpp"
                    "src/seissol_src/Initializer/MemoryAllocator.cpp"
                    "src/seissol_src/Kernels/TimeCommon.cpp"
                    "src/seissol_src/Kernels/DynamicRupture.cpp"
                    "src/seissol_src/Kernels/Plasticity.cpp")


set(SEISSOL_EQUATIONS_SOURCES "src/seissol_src/Equations/${EQUATIONS}/Kernels/Time.cpp"
                              "src/seissol_src/Equations/${EQUATIONS}/Kernels/Neighbor.cpp"
                              "src/seissol_src/Equations/${EQUATIONS}/Kernels/Local.cpp")

set(SEISSOL_INITIALIZER_SOURCES "src/seissol_src/Initializer/LTS.cpp"
                                "src/seissol_src/Initializer/tree/DeviceVarInfo.cpp")

# specify all necessary proxy files needed to build the application
set(SOURCES "src/flop_counter.cpp"
            "src/proxy_seissol.cpp")

# specify additionall files
# TODO: chech whether we build another executable
#set(PROFILING_COUNTERS_SOURCES "src/flops_per_cell.cpp")


list(APPEND SOURCES ${SEISSOL_SOURCES} 
                    ${SEISSOL_EQUATIONS_SOURCES} 
                    ${PROFILING_COUNTERS_SOURCES}
                    ${SEISSOL_INITIALIZER_SOURCES})
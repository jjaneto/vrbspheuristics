if (APPLE)
    cmake_minimum_required(VERSION 3.15)
endif ()

project(VRBSP-Heuristics)

OPTION(USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

if (APPLE AND USE_OpenMP)
    set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/10.0.0_3/bin/clang")
    set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/10.0.0_3/bin/clang++")
    set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/10.0.0_3/lib")
    set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/10.0.0_3/include")

    add_compile_definitions(DEBUG_CLION)

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif ()

    if (USE_OpenMP)
        find_package(OpenMP REQUIRED)
    endif (USE_OpenMP)

    if (OPENMP_FOUND)
        include_directories("${OPENMP_INCLUDES}")
        link_directories("${OPENMP_LIBRARIES}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif (OPENMP_FOUND)

    add_executable(main
#            HeuristicBRKGA.cpp
#            BRKGA.h
#            Population.h
#            Utility.cpp
#            Utility.h
#            DifferentialEvolution.h
#            DifferentialEvolution.cpp
            VNS.cpp
            MTRand.h
            HeuristicDecoder.cpp
            HeuristicDecoder.h)
else ()
    find_package(OpenMP REQUIRED)

    if (OpenMP_CXX_FOUND)
        add_executable(main
                HeuristicBRKGA.cpp
                BRKGA.h
                Population.h
                Utility.cpp
                Utility.h
#               VNS.cpp
                MTRand.h
                HeuristicDecoder.cpp
                HeuristicDecoder.h)
    endif ()

    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif ()

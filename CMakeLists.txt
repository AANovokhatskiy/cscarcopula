cmake_minimum_required(VERSION 3.12)
project(main)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(EIGEN_INCLUDE_DIR "/usr/include/eigen3/")
include_directories(${EIGEN_INCLUDE_DIR})

set(NLOPT_DIR "/usr/include/")
find_library(NLOPT_LIB nlopt PATHS ${NLOPT_DIR})
include_directories(${NLOPT_DIR})


# Add the source files
set(SOURCES
    main.cpp
    src/io/data.cpp
    src/pdf/generated_pdf.cpp
    src/sampler/scar-ou.cpp
    src/sampler/mle.cpp
    src/auxiliary.cpp
    src/copula.cpp
    src/marginals.cpp
    src/risk_metrics.cpp
    )

# Add the executable
add_executable(main ${SOURCES})

#target_link_libraries(main Eigen3::Eigen)
target_link_libraries(main ${NLOPT_LIB})
#target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)

target_compile_options(main PRIVATE -fopenmp)
target_compile_options(main PRIVATE -Ofast)
target_compile_options(main PRIVATE -lnlopt)
target_compile_options(main PRIVATE -lm)


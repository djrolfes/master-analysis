# Preloaded CMake configuration for klft build

# General build options
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
set(BUILD_TESTING ON CACHE BOOL "Build the testing tree")

# Compiler settings
set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx" CACHE FILEPATH "CXX compiler")
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-12.9/bin/nvcc" CACHE FILEPATH "CUDA compiler")

# Compiler flags
set(CMAKE_CXX_FLAGS "" CACHE STRING "Flags used by the CXX compiler during all build types")
# Ensure maximum debug info and no optimization for the host code
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "Flags used by the CXX compiler during DEBUG builds") # <-- KEY CHANGE: Added -O0 for host
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Flags used by the CXX compiler during RELEASE builds")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "Flags used by the CXX compiler during RELWITHDEBINFO builds")

set(CMAKE_CUDA_FLAGS "" CACHE STRING "Flags used by the CUDA compiler during all build types")
# Ensure maximum debug info (via -G) and no optimization for the device code
set(CMAKE_CUDA_FLAGS_DEBUG "-G -O0" CACHE STRING "Flags used by the CUDA compiler during DEBUG builds") # <-- KEY CHANGE: Used -G for device and added -O0
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Flags used by the CUDA compiler during RELEASE builds")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "Flags used by the CUDA compiler during RELWITHDEBINFO builds")

# Kokkos-specific flags
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Kokkos turns off CXX extensions")
set(Kokkos_ARCH_PASCAL60 ON CACHE BOOL "Enable Kokkos architecture for NVIDIA Pascal GPUs")
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP in Kokkos")
set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable CUDA in Kokkos")
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable Serial in Kokkos")

# YAML-cpp-specific flags
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Disable yaml-cpp tests")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Disable yaml-cpp tools")
set(YAML_CPP_INSTALL ON CACHE BOOL "Enable yaml-cpp installation targets")

# Installation paths
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Install path prefix, prepended onto install directories")
set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "User executables (bin)")
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Object code libraries (lib)")
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "C header files (include)")
# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.9)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(PARABOLOID_TEST VERSION 0.1
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Inform CMake where to look for CMake files.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
                      "${PROJECT_SOURCE_DIR}/cmake/Utils")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

set(SOURCE "${PROJECT_SOURCE_DIR}/src")

set(JSON "${PROJECT_SOURCE_DIR}/external/json/single_include")

# Find external packages needed by IRL
find_package(Eigen) # Provide -D EIGEN_PATH=/path/to/Eigen
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# Add AMR executable
add_executable(amr_generate_result)
set_target_properties(amr_generate_result
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
    INTERPROCEDURAL_OPTIMIZATION True
    CXX_STANDARD 17)

target_include_directories(amr_generate_result PRIVATE "./")
target_link_libraries(amr_generate_result PUBLIC "${IRL_INSTALL_LOCATION}/lib/libirl.a")
target_link_libraries(amr_generate_result PUBLIC "${IRL_INSTALL_LOCATION}/lib/libabsl_all.a")
target_link_libraries(amr_generate_result PUBLIC "-lpthread")
target_link_libraries(amr_generate_result PUBLIC MPI::MPI_CXX)
if(APPLE)
  target_link_libraries(amr_generate_result PUBLIC "-framework CoreFoundation")
endif()
target_include_directories(amr_generate_result PUBLIC "${IRL_INSTALL_LOCATION}/include")
target_include_directories(amr_generate_result PUBLIC "${IRL_ROOT_LOCATION}")
target_include_directories(amr_generate_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/abseil-cpp/")
target_include_directories(amr_generate_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/triangle/")
target_include_directories(amr_generate_result SYSTEM PUBLIC "${JSON}")
target_link_libraries(amr_generate_result INTERFACE PUBLIC Eigen::Eigen)
target_include_directories(amr_generate_result SYSTEM PUBLIC "")
target_include_directories(amr_generate_result SYSTEM PUBLIC "${MPI_INCLUDE_PATH}")


# Add IRL executable
add_executable(irl_confirm_result)
set_target_properties(irl_confirm_result
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
    INTERPROCEDURAL_OPTIMIZATION True
    CXX_STANDARD 17)

target_include_directories(irl_confirm_result PRIVATE "./")
target_link_libraries(irl_confirm_result PUBLIC "${IRL_INSTALL_LOCATION}/lib/libirl.a")
target_link_libraries(irl_confirm_result PUBLIC "${IRL_INSTALL_LOCATION}/lib/libabsl_all.a")
target_link_libraries(irl_confirm_result PUBLIC "-lpthread")
target_link_libraries(irl_confirm_result PUBLIC OpenMP::OpenMP_CXX)

if(APPLE)
  target_link_libraries(irl_confirm_result PUBLIC "-framework CoreFoundation")
endif()
target_include_directories(irl_confirm_result PUBLIC "${IRL_INSTALL_LOCATION}/include")
target_include_directories(irl_confirm_result PUBLIC "${IRL_ROOT_LOCATION}")
target_include_directories(irl_confirm_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/abseil-cpp/")
target_include_directories(irl_confirm_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/triangle/")
target_include_directories(irl_confirm_result SYSTEM PUBLIC "${JSON}")
target_link_libraries(irl_confirm_result INTERFACE PUBLIC Eigen::Eigen)
target_include_directories(irl_confirm_result SYSTEM PUBLIC "")

add_subdirectory("${SOURCE}")


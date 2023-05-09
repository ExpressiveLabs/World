# CMake build system design considerations:
# - Include directories:
#   + Do not define include directories globally using the include_directories
#     command but rather at the target level using the
#     target_include_directories command. That way, it is easier to guarantee
#     that targets are built using the proper list of include directories.
#   + Use the PUBLIC and PRIVATE keywords to specify the scope of include
#     directories. That way, a target linking to a library (using the
#     target_link_librairies command) inherits from the library PUBLIC include
#     directories and not from the PRIVATE ones.

cmake_minimum_required(VERSION 3.0)

# Build Options
option(WORLD_BUILD_TESTS "Set to ON to build tests" OFF)
option(WORLD_BUILD_EXAMPLES "Set to ON to build examples" ON)

# Adding Torch to the project
set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../libtorch-cpu")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# use git version as library version
find_package(Git QUIET)
if (Git_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --count HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE git_version
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else ()
    set(git_version 0)
endif ()

project(WORLD_TCH LANGUAGES CXX VERSION 1.0.${git_version})

add_library(world_tch STATIC
        src/TorchWorld/cheaptrick.h
        src/TorchWorld/codec.h
        src/TorchWorld/common.h
        src/TorchWorld/constantnumbers.h
        src/TorchWorld/d4c.h
        src/TorchWorld/dio.h
        src/TorchWorld/fft.h
        src/TorchWorld/harvest.h
        src/TorchWorld/macrodefinitions.h
        src/TorchWorld/matlabfunctions.h
        src/TorchWorld/stonemask.h
        src/TorchWorld/synthesis.h
        src/TorchWorld/synthesisrealtime.h
        src/ext/tqdm.h
        src/ext/utils.h
        src/cheaptrick.cpp
        src/codec.cpp
        src/common.cpp
        src/d4c.cpp
        src/dio.cpp
        src/fft.cpp
        src/harvest.cpp
        src/matlabfunctions.cpp
        src/stonemask.cpp
        src/synthesis.cpp
        src/synthesisrealtime.cpp
        )

add_library(world_tool_tch STATIC
        tools/audioio.h
        tools/parameterio.h
        tools/audioio.cpp
        tools/parameterio.cpp
        )

target_compile_features(world_tch PRIVATE cxx_std_20)
target_link_libraries(world_tch PUBLIC "${TORCH_LIBRARIES}")

add_library(TorchWorld::core ALIAS world_tch)
add_library(TorchWorld::tool ALIAS world_tool_tch)

target_link_libraries(world_tool_tch PUBLIC world_tch)

foreach (lib world_tch world_tool_tch)
    target_include_directories(${lib} PUBLIC $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> PRIVATE src)
    set_target_properties(${lib}
            PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )
endforeach ()

include(GNUInstallDirs)
install(TARGETS world_tch world_tool_tch
        EXPORT world-export
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
install(DIRECTORY src/TorchWorld DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
configure_package_config_file(world-config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/world-config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/world/cmake
        PATH_VARS CMAKE_INSTALL_INCLUDEDIR
        )
write_basic_package_version_file(
        world-config-version.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/world-config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/world
        )
install(EXPORT world-export
        FILE world-config-version.cmake
        NAMESPACE world::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/world
        )
export(TARGETS world_tch world_tool_tch NAMESPACE world::
        FILE ${PROJECT_BINARY_DIR}/world-targets.cmake)

if (WORLD_BUILD_TESTS)
    add_executable(tests test/test.cpp)
    target_link_libraries(tests
            PRIVATE world_tch
            PRIVATE world_tool_tch
            )
    target_include_directories(tests
            PRIVATE src
            PRIVATE tools
            )
endif ()

if(WORLD_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
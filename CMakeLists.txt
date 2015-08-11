cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." OFF) # Makes boolean 'test' available.

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )

# Make PROJECT_SOURCE_DIR, PROJECT_BINARY_DIR, and PROJECT_NAME available.
project(GLnexus)

ExternalProject_Add(htslib
    URL https://github.com/dnanexus/htslib/archive/3c8ae031.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j8
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(htslib source_dir)
set(HTSLIB_INCLUDE_DIR ${source_dir}/htslib)
include_directories(${HTSLIB_INCLUDE_DIR})
ExternalProject_Get_Property(htslib binary_dir)
add_library(libhts STATIC IMPORTED)
set(LIBHTS_A ${binary_dir}/libhts.a)
set_property(TARGET libhts PROPERTY IMPORTED_LOCATION ${LIBHTS_A})

ExternalProject_Add(rocksdb
    URL https://github.com/facebook/rocksdb/archive/rocksdb-3.10.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j8 OPT="-DNDEBUG -O2" static_lib
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(rocksdb source_dir)
set(ROCKSDB_INCLUDE_DIR ${source_dir}/include)
include_directories(${ROCKSDB_INCLUDE_DIR})
ExternalProject_Get_Property(rocksdb binary_dir)
add_library(librocksdb STATIC IMPORTED)
set(LIBROCKSDB_A ${binary_dir}/librocksdb.a)
set_property(TARGET librocksdb PROPERTY IMPORTED_LOCATION ${LIBROCKSDB_A})

ExternalProject_Add(yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/release-0.5.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND cmake -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF .
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j8 OPT="-DNDEBUG -O2"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(yaml-cpp source_dir)
set(YAML_CPP_INCLUDE_DIR ${source_dir}/include)
include_directories(${YAML_CPP_INCLUDE_DIR})
ExternalProject_Get_Property(yaml-cpp binary_dir)
add_library(libyaml-cpp STATIC IMPORTED)
set(LIBYAML_CPP_A ${binary_dir}/libyaml-cpp.a)
set_property(TARGET libyaml-cpp PROPERTY IMPORTED_LOCATION ${LIBYAML_CPP_A})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")

################################
# Normal Libraries & Executables
################################
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES} src)
include_directories(${COMMON_INCLUDES} src/rocksdb/include)

add_library(glnexus
            include/types.h src/types.cc
            include/data.h src/data.cc
            include/diploid.h src/diploid.cc
            include/service.h src/service.cc
            include/alleles.h src/alleles.cc
            include/genotyper.h src/genotyper.cc
            include/KeyValue.h src/KeyValue.cc
            include/BCFSerialize.h src/BCFSerialize.cc
            include/BCFKeyValueData.h src/BCFKeyValueData.cc
            include/RocksKeyValue.h src/RocksKeyValue.cc)
add_dependencies(glnexus htslib)
add_dependencies(glnexus rocksdb)
add_dependencies(glnexus yaml-cpp)
add_library(libglnexus ALIAS glnexus)

add_executable(playground src/playground.cc)
add_dependencies(playground libglnexus)
target_link_libraries(playground -pthread jemalloc glnexus libhts z librocksdb libyaml-cpp)

install(TARGETS playground DESTINATION bin)

################################
# Testing
################################
if (test)

  ExternalProject_Add(catch
    URL https://github.com/philsquared/Catch/archive/master.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )

  ExternalProject_Get_Property(catch source_dir)
  include_directories(${source_dir}/single_include)

  ##############
  # Unit Tests
  ##############
  add_executable(unit_tests
                 src/diploid.cc test/diploid.cc
                 test/unit_tests.cc
                 test/htslib_behaviors.cc
                 test/rocks_behaviors.cc
                 test/types.cc
                 test/service.cc
                 test/BCFKeyValueData.cc
                 test/rocks_integration.cc)
  add_dependencies(unit_tests libglnexus)
  target_link_libraries(unit_tests glnexus librocksdb -pthread jemalloc libhts rt snappy z bz2 lz4 rt libyaml-cpp)

  add_test(unit_tests unit_tests)
endif()

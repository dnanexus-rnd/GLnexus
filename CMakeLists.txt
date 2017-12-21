cmake_minimum_required(VERSION 3.2)
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

INCLUDE(ExternalProject)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR}/external/include)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/lib)
SET(ENV{PATH} ${CMAKE_BINARY_DIR}/external/bin:$ENV{PATH})
SET(ENV{PKG_CONFIG_PATH} ${CMAKE_BINARY_DIR}/external/lib/pkgconfig)

# After build, double-check external/src/htslib-stamp/htslib-build-out.log
# to make sure CFLAGS were set as expected
ExternalProject_Add(htslib
    URL https://github.com/samtools/htslib/archive/1.4.1.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    PATCH_COMMAND sed -i "s/^CFLAGS .*$/CFLAGS = -gdwarf -O3 -DNDEBUG -march=ivybridge/" Makefile
    BUILD_IN_SOURCE 1
    BUILD_COMMAND bash -c "make -n && make -j8"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_BUILD ON
  )
ExternalProject_Get_Property(htslib source_dir)
set(HTSLIB_INCLUDE_DIR ${source_dir}/htslib)
include_directories(${HTSLIB_INCLUDE_DIR})
ExternalProject_Get_Property(htslib binary_dir)
add_library(libhts STATIC IMPORTED)
set(LIBHTS_A ${binary_dir}/libhts.a)
set_property(TARGET libhts PROPERTY IMPORTED_LOCATION ${LIBHTS_A})

# RocksDB
#
# We set DISABLE_JEMALLOC=1 to avoid link-time dependency to jemalloc; we
# DO however wish to use jemalloc, but expect user to activate it at
# runtime with LD_PRELOAD instead.
#
# Warning: it's dangerously easy to mess with the RocksDB external project
# configuration below and silently get a librocksdb.a compiled without
# optimizations, or targeting an instruction set that's too old or too new.
#
# Here are some checks to perform after any adjustment:
# 1. Verify AVX is used (nonzero results):
#      objdump -d external/src/rocksdb/librocksdb.a | grep vzeroupper | wc -l
# 2. Verify Haswell-and-newer instructions are NOT used (zero results):
#      objdump -d external/src/rocksdb/librocksdb.a | grep shlx
# 3. Verify accelerated CRC32 is used (nonzero results):
#      objdump -d external/src/rocksdb/librocksdb.a | grep crc32q
# 4. Verify compiler optimizations performed:
#      (TODO better heuristics needed)
#      a) check external/src/rocksdb-stamp/rocksdb-build-out.log to see
#         compiler invocations with -O3
#      b) go through perf/FlameGraph profiling flow and check that no
#         rocksdb function has a big chunk like 5% or more
ExternalProject_Add(rocksdb
    URL https://github.com/facebook/rocksdb/archive/v5.9.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND bash -c "export PORTABLE=1 && export DISABLE_JEMALLOC=1 && export USE_SSE=1 && export OPT='-DNDEBUG -O3 -march=ivybridge -msse4.2 -DHAVE_SSE42' && make -n static_lib && make -j8 static_lib"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_BUILD ON
  )
ExternalProject_Get_Property(rocksdb source_dir)
set(ROCKSDB_INCLUDE_DIR ${source_dir}/include)
include_directories(${ROCKSDB_INCLUDE_DIR})
ExternalProject_Get_Property(rocksdb binary_dir)
add_library(librocksdb STATIC IMPORTED)
set(LIBROCKSDB_A ${binary_dir}/librocksdb.a)
set_property(TARGET librocksdb PROPERTY IMPORTED_LOCATION ${LIBROCKSDB_A})

# After build, double-check external/src/yaml-cpp-stamp/yaml-cpp-build-out.log
# to make sure CFLAGS were set as expected
ExternalProject_Add(yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/005a6a19ee0af3e284637abe9ecea5501993fc4c.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND cmake -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-march=ivybridge .
    BUILD_IN_SOURCE 1
    BUILD_COMMAND bash -c "make -n && make -j8"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_BUILD ON
  )
ExternalProject_Get_Property(yaml-cpp source_dir)
set(YAML_CPP_INCLUDE_DIR ${source_dir}/include)
include_directories(${YAML_CPP_INCLUDE_DIR})
ExternalProject_Get_Property(yaml-cpp binary_dir)
add_library(libyaml-cpp STATIC IMPORTED)
set(LIBYAML_CPP_A ${binary_dir}/libyaml-cpp.a)
set_property(TARGET libyaml-cpp PROPERTY IMPORTED_LOCATION ${LIBYAML_CPP_A})

ExternalProject_Add(CTPL
    URL https://github.com/vit-vit/CTPL/archive/v.0.0.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(CTPL source_dir)
set(CTPL_INCLUDE_DIR ${source_dir})
include_directories(${CTPL_INCLUDE_DIR})

ExternalProject_Add(fcmm
    URL https://github.com/giacomodrago/fcmm/archive/v1.0.1.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(fcmm source_dir)
set(FCMM_INCLUDE_DIR ${source_dir})
include_directories(${FCMM_INCLUDE_DIR})

ExternalProject_Add(spdlog
    URL https://github.com/gabime/spdlog/archive/c900b25dfd.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(spdlog source_dir)
set(SPDLOG_INCLUDE_DIR ${source_dir}/include)
include_directories(${SPDLOG_INCLUDE_DIR})

ExternalProject_Add(capnp
    URL https://capnproto.org/capnproto-c++-0.6.1.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/external
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j4 check
    INSTALL_COMMAND make install
    LOG_DOWNLOAD ON
  )

execute_process(COMMAND git describe --tags --long --dirty --always
                OUTPUT_VARIABLE GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_REVISION=\"\\\"${GIT_REVISION}\\\"\" -pthread -std=c++14 -Wall -Werror=return-type -Werror=unused-result -Wno-sign-compare -fdiagnostics-color=auto -march=ivybridge")
set(CMAKE_CXX_FLAGS_RELEASE "-gdwarf -DNDEBUG -O3")

################################
# Normal Libraries & Executables
################################
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES} src)
include_directories(${COMMON_INCLUDES} src/rocksdb/include)
include_directories(${COMMON_INCLUDES} capnp/serialize)

# Invoke cap'n proto IDL compiler, to generate C++ serialization code from
# definitions.
add_custom_command(
  OUTPUT  capnp/serialize/defs.capnp.c++ capnp/serialize/defs.capnp.h
  COMMAND ${CMAKE_BINARY_DIR}/external/bin/capnp compile -I${CMAKE_CURRENT_BINARY_DIR}/external/include -oexternal/bin/capnpc-c++ capnp/serialize/defs.capnp
  DEPENDS capnp/serialize/defs.capnp
  COMMENT "GLnexus: generating cap'n proto C++ and H files")

add_library(glnexus
            capnp/serialize/defs.capnp.h capnp/serialize/defs.capnp.c++
            include/types.h src/types.cc
            include/data.h src/data.cc
            include/compare_queries.h src/compare_queries.cc
            include/diploid.h src/diploid.cc
            include/service.h src/service.cc
            include/discovery.h src/discovery.cc
            include/unifier.h src/unifier.cc
            include/genotyper.h src/genotyper.cc
            include/genotyper_utils.h
            include/residuals.h src/residuals.cc
            include/KeyValue.h src/KeyValue.cc
            include/BCFSerialize.h src/BCFSerialize.cc
            include/BCFKeyValueData.h src/BCFKeyValueData.cc
            include/RocksKeyValue.h src/RocksKeyValue.cc
            include/cli_utils.h src/cli_utils.cc
            test/utils.cc)
add_dependencies(glnexus htslib)
add_dependencies(glnexus rocksdb)
add_dependencies(glnexus yaml-cpp)
add_dependencies(glnexus CTPL)
add_dependencies(glnexus fcmm)
add_dependencies(glnexus spdlog)
add_dependencies(glnexus capnp)
add_library(libglnexus ALIAS glnexus)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -march=ivybridge -static-libstdc++ -pthread")

add_executable(playground src/playground.cc)
add_dependencies(playground libglnexus)
target_link_libraries(playground glnexus libhts librocksdb libyaml-cpp libz.a libsnappy.a libbz2.a libzstd.a liblzma.a librt.a)

add_executable(glnexus_cli cli/glnexus_cli.cc)
add_dependencies(glnexus_cli libglnexus)
target_link_libraries(glnexus_cli glnexus libhts librocksdb libyaml-cpp libz.a libsnappy.a libbz2.a libzstd.a liblzma.a librt.a libcapnp.a libkj.a)

install(TARGETS playground DESTINATION bin)

################################
# Testing
################################
if (test)

  ExternalProject_Add(catch
    URL https://github.com/philsquared/Catch/archive/v1.8.2.zip
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
                 test/genotyper.cc
                 test/service.cc
                 test/gvcf_test_cases.cc
                 test/BCFKeyValueData.cc
                 test/rocks_integration.cc
                 test/unifier.cc
                 test/cli_utils.cc
                 test/cli.cc)
  add_dependencies(unit_tests catch)
  add_dependencies(unit_tests libglnexus)
  target_link_libraries(unit_tests glnexus libhts librocksdb libyaml-cpp libz.a libsnappy.a libbz2.a libzstd.a liblzma.a librt.a libcapnp.a libkj.a)

  add_test(unit_tests unit_tests)
endif()

# ======================================================================
#  dunesim main build file
# ======================================================================


cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)


find_package(cetmodules REQUIRED)
project(dunesim LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.04.06d00)
set(${PROJECT_NAME}_FW_DIR "config_data")

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package( art REQUIRED EXPORT )
find_package( art_root_io REQUIRED EXPORT )
find_package( canvas_root_io REQUIRED EXPORT )
find_package( Boost REQUIRED )
find_package(ROOT REQUIRED)
find_package( larevt REQUIRED EXPORT )
find_package( larsim REQUIRED EXPORT )
find_package( lardataalg REQUIRED EXPORT )
find_package( lardataobj REQUIRED EXPORT )
find_package( larcore REQUIRED EXPORT )
find_package( larcorealg REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT )
find_package( ifdhc REQUIRED )
find_package( nusimdata REQUIRED EXPORT )
find_package( nurandom REQUIRED EXPORT )
find_package( dunecore REQUIRED EXPORT )
find_package( dunepdlegacy REQUIRED EXPORT )
find_package( lardata REQUIRED EXPORT )
find_package( CLHEP REQUIRED EXPORT )
find_package( nugen REQUIRED EXPORT )
find_package( Geant4 REQUIRED EXPORT )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunesim)
add_subdirectory(macros)

cet_cmake_config()

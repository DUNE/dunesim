# ======================================================================
#  dunesim main build file
# ======================================================================


cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)


find_package(cetmodules REQUIRED)
project(dunesim LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.83.01d00)
set(${PROJECT_NAME}_FW_DIR "config_data")

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package( art REQUIRED )
find_package( art_root_io REQUIRED )
find_package( canvas_root_io REQUIRED )
find_package( Boost REQUIRED )
find_package(ROOT REQUIRED)
find_package( larevt REQUIRED )
find_package( larsim REQUIRED )
find_package( lardataalg REQUIRED )
find_package( lardataobj REQUIRED )
find_package( larcore REQUIRED )
find_package( larcorealg REQUIRED )
find_package( larcoreobj REQUIRED )
find_package( nusimdata REQUIRED )
find_package( nurandom REQUIRED )
find_package( dunecore REQUIRED )
find_package( dunepdlegacy REQUIRED )
find_package( lardata REQUIRED )
find_package( CLHEP REQUIRED )
find_package( nugen REQUIRED )
find_package( Geant4 REQUIRED )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunesim)
add_subdirectory(macros)

cet_cmake_config()

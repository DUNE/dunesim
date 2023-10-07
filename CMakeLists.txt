# ======================================================================
#  dunesim main build file
# ======================================================================


cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(dunesim LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.79.00d00)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_boost( )
find_ups_root()
find_ups_product( cetbuildtools )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( lardataalg )
find_ups_product( lardataobj )
find_ups_product( larcore )
find_ups_product( larcorealg )
find_ups_product( larcoreobj )
find_ups_product( nusimdata )
find_ups_product( nurandom )
find_ups_product( dunecore )
find_ups_product( lardata )
find_ups_product( clhep )
find_ups_product( nugen )
find_ups_geant4( )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunesim)
add_subdirectory(ups)
add_subdirectory(macros)

include(UseCPack)

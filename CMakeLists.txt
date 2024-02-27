# ======================================================================
#  dunesim main build file
# ======================================================================


cmake_minimum_required (VERSION 3.20 FATAL_ERROR) ### MIGRATE-CHANGED (migrate-3.22.02) - VERSION 3.19 -> 3.20
# cmake_policy(VERSION 3.19) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.22.02) - Uncomment to preserve compatibility with older CMake versions


find_package(cetmodules REQUIRED) ### MIGRATE-ADDED (migrate-3.22.02) - find_package(cetmodules) must precede project()
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

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - add args: REPORT_THRESHOLD VERBOSE

find_package( art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art ...) -> find_package(art REQUIRED ...)
find_package( art_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art_root_io ...) -> find_package(art_root_io REQUIRED ...)
find_package( canvas_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(canvas_root_io ...) -> find_package(canvas_root_io REQUIRED ...)
find_package( Boost REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(ROOT REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package( larevt REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larevt ...) -> find_package(larevt REQUIRED ...)
find_package( larsim REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larsim ...) -> find_package(larsim REQUIRED ...)
find_package( lardataalg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataalg ...) -> find_package(lardataalg REQUIRED ...)
find_package( lardataobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataobj ...) -> find_package(lardataobj REQUIRED ...)
find_package( larcore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcore ...) -> find_package(larcore REQUIRED ...)
find_package( larcorealg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcorealg ...) -> find_package(larcorealg REQUIRED ...)
find_package( larcoreobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcoreobj ...) -> find_package(larcoreobj REQUIRED ...)
find_package( nusimdata REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nusimdata ...) -> find_package(nusimdata REQUIRED ...)
find_package( nurandom REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nurandom ...) -> find_package(nurandom REQUIRED ...)
find_package( dunecore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunecore ...) -> find_package(dunecore REQUIRED ...)
find_package( lardata REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardata ...) -> find_package(lardata REQUIRED ...)
find_package( CLHEP REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(clhep ...) -> find_package(clhep REQUIRED ...)
find_package( nugen REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nugen ...) -> find_package(nugen REQUIRED ...)
find_package( Geant4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunesim)
add_subdirectory(macros)

cet_cmake_config() ### MIGRATE-ADDED (migrate-3.22.02) - required cmd

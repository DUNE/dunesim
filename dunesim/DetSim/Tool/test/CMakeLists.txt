# dunetpc/dune/DetSim/Tool/test/CMakeLists.txt

# David Adams
# April 2017

# Build test for each service.

include(CetTest)

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})

cet_enable_asserts()

cet_test(test_IdealAdcSimulator SOURCES test_IdealAdcSimulator.cxx
  LIBRARIES
    art::Utilities
    dunecore::DuneInterface_Data
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    SQLite::SQLite3
    Boost::filesystem
)

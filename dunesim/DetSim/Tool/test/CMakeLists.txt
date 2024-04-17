# dunetpc/dune/DetSim/Tool/test/CMakeLists.txt

# David Adams
# April 2017

# Build test for each service.

include(CetTest)

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

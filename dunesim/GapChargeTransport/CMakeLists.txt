cet_make_library(
  LIBRARY_NAME dunesim_GapInfo
  SOURCE GapInfo.cc
  LIBRARIES
)

cet_build_plugin(
    GapChargeTransport
    art::tool
    LIBRARIES PRIVATE
        dunesim_GapInfo
        art::Utilities
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        larcoreobj::geo_vectors
        lardataalg::DetectorInfo
        larsim::IonizationScintillation   # <-- added dependency
)

install_headers()
install_source()
install_fhicl()

# dunetpc/dune/DetSim/Module/CMakeLists.txt
#
# David Adams
# January 2017

art_make(
  BASENAME_ONLY
  MODULE_LIBRARIES
    lardataalg::DetectorInfo
    dunecore::ArtSupport
    dunesim::DetSim
    dunecore::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
    lardata::Utilities
    larevt::Filters
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::SignalShapingServiceDUNE_service
    dunecore::SignalShapingServiceDUNEDPhase_service
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    ROOT::Core
    art_root_io::TFileService_service
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities 
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib 
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    ROOT::Geom
    ROOT::XMLIO
    ROOT::Gdml
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_source()

# dunetpc/dune/DetSim/Module/CMakeLists.txt
#
# David Adams
# January 2017

art_make(
  BASENAME_ONLY
  MODULE_LIBRARIES
    lardataalg_DetectorInfo
    dune_ArtSupport
    dune_DetSim
    dune_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    larcorealg_Geometry
    lardataobj_RawData
    SignalShapingServiceDUNE_service
    SignalShapingServiceDUNEDPhase_service
    nurandom::RandomUtils_NuRandomService_service
    ART_FRAMEWORK_CORE
    ART_FRAMEWORK_PRINCIPAL
    ART_FRAMEWORK_SERVICES_REGISTRY
    ART_ROOT_IO_TFILE_SUPPORT
    ROOT::Core
    ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
    ART_ROOT_IO_TFILESERVICE_SERVICE
    art::Persistency_Common canvas
    art::Persistency_Provenance canvas
    art::Utilities canvas
    messagefacility::MF_MessageLogger
    
    fhiclcpp::fhiclcpp
    cetlib::cetlib cetlib_except
    CLHEP
    ROOT_GEOM
    ROOT_XMLIO
    ROOT_GDML
    ROOT_BASIC_LIB_LIST
)

install_source()

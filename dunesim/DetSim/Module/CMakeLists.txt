# dunetpc/dune/DetSim/Module/CMakeLists.txt
#
# David Adams
# January 2017

art_make(
  BASENAME_ONLY
  MODULE_LIBRARIES
    dune_ArtSupport
    dune_DetSim
    dune_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    larcorealg_Geometry
    lardataobj_RawData
    SignalShapingServiceDUNE_service
    SignalShapingServiceDUNEDPhase_service
    nutools_RandomUtils_NuRandomService_service
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art_Persistency_Common canvas
    art_Persistency_Provenance canvas
    art_Utilities canvas
    ${MF_MESSAGELOGGER}
    
    ${FHICLCPP}
    cetlib cetlib_except
    ${CLHEP}
    ${ROOT_GEOM}
    ${ROOT_XMLIO}
    ${ROOT_GDML}
    ${ROOT_BASIC_LIB_LIST}
)

install_source()

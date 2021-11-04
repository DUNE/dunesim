include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
                        larcoreobj_SummaryData
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        nurandom_RandomUtils_NuRandomService_service
                        larsim_Simulation lardataobj_Simulation
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities canvas
                        messagefacility::MF_MessageLogger
                        hep_concurrency
                        FHICLCPP
                        cetlib cetlib_except
                        CLHEP
                        ROOT_BASIC_LIB_LIST
                        SQLITE3
                        ROOT_EG
                        ifdh::ifdh
                        dune_EventGenerator_ProtoDUNEbeamDataProducts
 )


install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(Radiological)
add_subdirectory(CORSIKAprotodunedp)
add_subdirectory(NEUT)
add_subdirectory(MUSUN)
add_subdirectory(NuE)
add_subdirectory(ProtoDUNEbeamDataProducts)

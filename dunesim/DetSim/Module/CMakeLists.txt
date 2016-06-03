# dunetpc/dune//larsim_DetSim/Module/CMakeLists.txt
#
# David Adams
# December 2015

art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcore_Geometry
                        dune_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardata_RawData
                        larsim_DetSim
                        SignalShapingServiceDUNE_service
			SignalShapingServiceDUNEDPhase_service
                        larsim_RandomUtils_LArSeedService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
#install_fhicl()
install_source()


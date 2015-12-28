# dunetpc/dune//DetSim/Module/CMakeLists.txt
#
# David Adams
# December 2015

art_make( BASENAME_ONLY
          MODULE_LIBRARIES Geometry
                        dune_Geometry
                        Geometry_service
                        Simulation
                        Utilities
                        Filters
                        RawData
                        DetSim
                        TimeService_service
                        SignalShapingServiceDUNE10kt_service
                        SignalShapingServiceDUNE35t_service
                        ${SEEDSERVICE_SERVICE} # artextensions
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
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

#install_headers()
#install_fhicl()
install_source()


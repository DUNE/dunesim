# dunetpc/dune/DetSim/Service/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for DetSim services.

# Add plugin for each service.

simple_plugin(Dune35tSimChannelExtractService   "service"
                DetSim
                Geometry
                Simulation
                Utilities
                SignalShapingServiceDUNE35t_service
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GenericSimChannelExtractService   "service"
                DetSim
                Simulation
                Utilities
                SignalShapingServiceDUNE35t_service
                ${ART_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DPhaseSimChannelExtractService   "service"
                DetSim
                Simulation
                Utilities
                SignalShapingServiceDUNEDPhase_service
                ${ART_UTILITIES}
		${ART_FRAMEWORK_CORE}
                ${ROOT_BASIC_LIB_LIST}
		${SEEDSERVICE_SERVICE}
		${CLHEP}
             )
simple_plugin(ExponentialChannelNoiseService   "service"
                Geometry
                SignalShapingServiceDUNE35t_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                ${ART_UTILITIES}
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(WhiteChannelNoiseService   "service"
                Geometry
                SignalShapingServiceDUNE35t_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                ${ART_UTILITIES}
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProvidedPedestalAdditionService "service"
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                ${ART_UTILITIES}
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tZeroSuppressService        "service"
                DetSim
)

simple_plugin(Legacy35tZeroSuppressService  "service"
                DetSim
)

simple_plugin(FixedZeroSuppressService      "service"
)

simple_plugin(ReplaceCompressService        "service"
)

simple_plugin(LarsoftHuffmanCompressService "service"
                RawData
                dune_DetSim_Service_ReplaceCompressService_service
)

simple_plugin(StuckBitAdcDistortionService   "service"
                Geometry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                ${ART_UTILITIES}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
                BASENAME_ONLY
             )

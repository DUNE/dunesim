# dunetpc/dune/larsim_DetSim/Service/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for larsim_DetSim services.

# Add plugin for each service.

simple_plugin(Dune35tSimChannelExtractService   "service"
                larsim_DetSim
                larcore_Geometry
                larsim_Simulation larsimobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNE_service
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GenericSimChannelExtractService   "service"
                larsim_DetSim
                larsim_Simulation larsimobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNE_service
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DPhaseSimChannelExtractService   "service"
                larsim_DetSim
                larsim_Simulation larsimobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNEDPhase_service
                art_Utilities canvas_Utilities
		${ART_FRAMEWORK_CORE}
                ${ROOT_BASIC_LIB_LIST}
		larsim_RandomUtils_LArSeedService_service
		${CLHEP}
             )
simple_plugin(ExponentialChannelNoiseService   "service"
                larcore_Geometry
                SignalShapingServiceDUNE_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                larsim_RandomUtils_LArSeedService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(WhiteChannelNoiseService   "service"
                larcore_Geometry
                SignalShapingServiceDUNE_service
                larsim_RandomUtils_LArSeedService_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProvidedPedestalAdditionService "service"
                larsim_RandomUtils_LArSeedService_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                MF_MessageLogger
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tZeroSuppressService        "service"
                larsim_DetSim
)

simple_plugin(Legacy35tZeroSuppressService  "service"
                larsim_DetSim
)

simple_plugin(FixedZeroSuppressService      "service"
)

simple_plugin(ReplaceCompressService        "service"
)

simple_plugin(LarsoftHuffmanCompressService "service"
                lardataobj_RawData
                dune_DetSim_Service_ReplaceCompressService_service
)

simple_plugin(StuckBitAdcDistortionService   "service"
                larcore_Geometry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                larsim_RandomUtils_LArSeedService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FixedDetPedestalService "service"
                larcore_Geometry
                art_Utilities canvas_Utilities
                ${ROOT_BASIC_LIB_LIST}
)

add_subdirectory(test)

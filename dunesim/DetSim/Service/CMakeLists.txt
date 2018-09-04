# dunetpc/dune/DetSim/Service/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for DetSim services.

# Add plugin for each service.

simple_plugin(Dune35tSimChannelExtractService   "service"
                dune_DetSim
                larcorealg_Geometry
                larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNE_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GenericSimChannelExtractService   "service"
                dune_DetSim
                larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNE_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DPhaseSimChannelExtractService   "service"
                dune_DetSim
                larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNEDPhase_service
                art_Utilities canvas
		${ART_FRAMEWORK_CORE}
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
		nutools_RandomUtils_NuRandomService_service
		${CLHEP}
             )

simple_plugin(ExponentialChannelNoiseService   "service"
                larcorealg_Geometry
                SignalShapingServiceDUNE_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DPhaseRealisticNoiseService   "service"
               larcorealg_Geometry
               SignalShapingServiceDUNEDPhase_service
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               nutools_RandomUtils_NuRandomService_service
               ${ART_FRAMEWORK_CORE}
               art_Utilities canvas
               cetlib cetlib_except
               ${CLHEP}
               ${ROOT_BASIC_LIB_LIST}
             )

 simple_plugin(DPhaseCoherentNoiseService   "service"
                larcorealg_Geometry
                SignalShapingServiceDUNEDPhase_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
              )

simple_plugin(ProtoDUNEChannelNoiseService   "service"
                larcorealg_Geometry
                SignalShapingServiceDUNE_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )


simple_plugin(WhiteChannelNoiseService   "service"
                larcorealg_Geometry
                SignalShapingServiceDUNE_service
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ProvidedPedestalAdditionService "service"
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                cetlib cetlib_except
                MF_MessageLogger
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(Dune35tZeroSuppressService        "service"
                dune_DetSim
                cetlib cetlib_except
)

simple_plugin(Legacy35tZeroSuppressService  "service"
                dune_DetSim
                cetlib cetlib_except
)

simple_plugin(FixedZeroSuppressService      "service"
                cetlib cetlib_except
)

simple_plugin(ReplaceCompressService        "service"
                cetlib cetlib_except
)

simple_plugin(LarsoftHuffmanCompressService "service"
                lardataobj_RawData
                dune_DetSim_Service_ReplaceCompressService_service
                cetlib cetlib_except
)

simple_plugin(StuckBitAdcDistortionService   "service"
                larcorealg_Geometry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                ${MF_MESSAGELOGGER}

                cetlib cetlib_except
                ${CLHEP}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FixedDetPedestalService "service"
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(DPhaseFixedDetPedestalService "service"
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
)

add_subdirectory(test)

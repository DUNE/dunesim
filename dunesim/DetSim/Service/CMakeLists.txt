# dunetpc/dune/DetSim/Service/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for DetSim services.

# Add plugin for each service.

cet_build_plugin(Dune35tSimChannelExtractService   art::service
                dunesim_DetSim
                larcorealg_Geometry
                larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                dunecore::SignalShapingServiceDUNE_service
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(GenericSimChannelExtractService   art::service
                dunesim_DetSim
                larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                dunecore::SignalShapingServiceDUNE_service
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseSimChannelExtractService   art::service
                dunesim_DetSim
	        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                dunecore::SignalShapingServiceDUNEDPhase_service
                dunecore::CrpGainService_service
                art::Utilities canvas
		ART_FRAMEWORK_CORE
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
		nurandom::RandomUtils_NuRandomService_service
		CLHEP
             )

cet_build_plugin(ExponentialChannelNoiseService   art::service
                larcorealg_Geometry
                dunecore::SignalShapingServiceDUNE_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )
              
cet_build_plugin(ShapedCohProtoDUNENoiseService   art::service
                dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                larcorealg_Geometry
                dunecore::SignalShapingServiceDUNE_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseRealisticNoiseService   art::service
               lardata_Utilities_LArFFT_service
	       lardataalg_DetectorInfo
               larcorealg_Geometry
               dunecore::SignalShapingServiceDUNEDPhase_service
               ART_ROOT_IO_TFILE_SUPPORT
               ROOT::Core
               ART_ROOT_IO_TFILESERVICE_SERVICE
               nurandom::RandomUtils_NuRandomService_service
               ART_FRAMEWORK_CORE
               art::Utilities canvas
               cetlib::cetlib cetlib_except
               CLHEP
               ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseCoherentNoiseService   art::service
                larcorealg_Geometry
 	        lardataalg_DetectorInfo
                dunecore::SignalShapingServiceDUNEDPhase_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
              )

cet_build_plugin(ProtoDUNEChannelNoiseService   art::service
                larcorealg_Geometry
                dunecore::SignalShapingServiceDUNE_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(SPhaseChannelNoiseService   art::service
                larcorealg_Geometry
                dunecore::SignalShapingServiceDUNE_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )


cet_build_plugin(WhiteChannelNoiseService   art::service
                larcorealg_Geometry
                dunecore::SignalShapingServiceDUNE_service
                nurandom::RandomUtils_NuRandomService_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ProvidedPedestalAdditionService art::service
                nurandom::RandomUtils_NuRandomService_service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                messagefacility::MF_MessageLogger
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(Dune35tZeroSuppressService        art::service
                dunesim_DetSim
                cetlib::cetlib cetlib_except
)

cet_build_plugin(Legacy35tZeroSuppressService  art::service
                dunesim_DetSim
                cetlib::cetlib cetlib_except
)

cet_build_plugin(FixedZeroSuppressService      art::service
                cetlib::cetlib cetlib_except
)

cet_build_plugin(ReplaceCompressService        art::service
                cetlib::cetlib cetlib_except
)

cet_build_plugin(LarsoftHuffmanCompressService art::service
                lardataobj_RawData
                dunesim_DetSim_Service_ReplaceCompressService_service
                cetlib::cetlib cetlib_except
)

cet_build_plugin(StuckBitAdcDistortionService   art::service
                larcorealg_Geometry
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                nurandom::RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                art::Utilities canvas
                messagefacility::MF_MessageLogger

                cetlib::cetlib cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FixedDetPedestalService art::service
                larcorealg_Geometry
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
)

cet_build_plugin(DPhaseFixedDetPedestalService art::service
                larcorealg_Geometry
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
)

add_subdirectory(test)

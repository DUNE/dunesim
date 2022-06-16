# dunetpc/dune/DetSim/Service/CMakeLists.txt
#
# David Adams
# December 2015
#
# CMake instructions for DetSim services.

# Add plugin for each service.

cet_build_plugin(Dune35tSimChannelExtractService   
                art::service
                dunesim::DetSim
                larcorealg::Geometry
                larsim::Simulation 
                nug4::ParticleNavigation 
                lardataobj::Simulation
                lardata::Utilities
                dunecore::SignalShapingServiceDUNE_service
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(GenericSimChannelExtractService   
                art::service
                dunesim::DetSim
                larsim::Simulation 
		nug4::ParticleNavigation 
		lardataobj::Simulation
                lardata::Utilities
                dunecore::SignalShapingServiceDUNE_service
                art::Utilities 
		canvas::canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseSimChannelExtractService   
                art::service
                dunesim::DetSim
	        larsim::Simulation 
                nug4::ParticleNavigation 
		lardataobj::Simulation
                lardata::Utilities
                dunecore::SignalShapingServiceDUNEDPhase_service
                dunecore::CrpGainService_service
                art::Utilities canvas
                art::Framework_Core
                cetlib::cetlib 
		cetlib_except::cetlib_except
                ROOT_BASIC_LIB_LIST
		nurandom::RandomUtils_NuRandomService_service
		CLHEP
             )

cet_build_plugin(ExponentialChannelNoiseService   
                art::service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )
              
cet_build_plugin(ShapedCohProtoDUNENoiseService   art::service
                dunepdlegacy_Services_ChannelMap_PdspChannelMapService_service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities 
		canvas::canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseRealisticNoiseService   
               art::service
               lardata::Utilities_LArFFT_service
	       lardataalg::DetectorInfo
               larcorealg::Geometry
               dunecore::SignalShapingServiceDUNEDPhase_service
               art_root_io::tfile_support
               ROOT::Core
               art_root_io::TFileService_service
               nurandom::RandomUtils_NuRandomService_service
               art::Framework_Core
               art::Utilities 
	       canvas::canvas
               cetlib::cetlib 
	       cetlib_except::cetlib_except
               CLHEP
               ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DPhaseCoherentNoiseService   
                art::service
                larcorealg::Geometry
 	        lardataalg::DetectorInfo
                dunecore::SignalShapingServiceDUNEDPhase_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities 
		canvas::canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
              )

cet_build_plugin(ProtoDUNEChannelNoiseService   
                art::service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(SPhaseChannelNoiseService   
                art::service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(SPhaseUniqChannelNoiseService   
                art::service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(WhiteChannelNoiseService   
                art::service
                larcorealg::Geometry
                dunecore::SignalShapingServiceDUNE_service
                nurandom::RandomUtils_NuRandomService_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                art::Framework_Core
                art::Utilities 
		canvas::canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ProvidedPedestalAdditionService 
                art::service
                nurandom::RandomUtils_NuRandomService_service
                art_root_io::tfile_support
                ROOT::Core
                art_root_io::TFileService_service
                art::Framework_Core
                art::Utilities 
		canvas::canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                messagefacility::MF_MessageLogger
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(Dune35tZeroSuppressService        
                art::service
                dunesim_DetSim
                cetlib::cetlib 
		cetlib_except::cetlib_except
)

cet_build_plugin(Legacy35tZeroSuppressService  
                art::service
                dunesim_DetSim
                cetlib::cetlib 
		cetlib_except::cetlib_except
)

cet_build_plugin(FixedZeroSuppressService      
                art::service
                cetlib::cetlib 
		cetlib_except::cetlib_except
)

cet_build_plugin(ReplaceCompressService        
                art::service
                cetlib::cetlib 
		cetlib_except::cetlib_except
)

cet_build_plugin(LarsoftHuffmanCompressService 
                art::service
                lardataobj::RawData
                dunesim::DetSim_Service_ReplaceCompressService_service
                cetlib::cetlib 
		cetlib_except::cetlib_except
)

cet_build_plugin(StuckBitAdcDistortionService   
                art::service
                larcorealg::Geometry
                art_root_io::tfile_support
                ROOT::Core
                nurandom::RandomUtils_NuRandomService_service
                art::Framework_Core
                art::Utilities 
		canvas::canvas
                messagefacility::MF_MessageLogger
                cetlib::cetlib 
		cetlib_except::cetlib_except
                CLHEP
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FixedDetPedestalService 
                art::service
                larcorealg::Geometry
                art::Utilities canvas
                cetlib::cetlib 
		cetlib_except::cetlib_except
                ROOT_BASIC_LIB_LIST
)

#cet_build_plugin(DPhaseFixedDetPedestalService 
#                art::service
#                larcorealg::Geometry
#                art::Utilities canvas
#                cetlib::cetlib 
#                cetlib_except::cetlib_except
#                ROOT_BASIC_LIB_LIST
#)

add_subdirectory(test)

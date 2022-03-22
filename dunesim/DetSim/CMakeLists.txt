# dunetpc/dune/DetSim/CMakeLists.txt

# David Adams
# December 2015

# We can drop the MODULE_LIBRARIES block after the old modules
# SimWireDUNE35t_module.cc and SimWireDUNE10kt_module.cc are
# removed.

art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcorealg::Geometry
	                   dunecore::Geometry
 	                   lardataalg::DetectorInfo
			   larcore::Geometry_Geometry_service
			   larsim::Simulation 
			   nug4::ParticleNavigation 
			   lardataobj::Simulation
	                   lardata::Utilities
			   larevt::Filters
			   lardataobj::RawData
                           dunecore::SignalShapingServiceDUNE_service
			   dunecore::SignalShapingServiceDUNEDPhase_service
                           nurandom::RandomUtils_NuRandomService_service
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           ROOT::Core
                           art_root_io::TFileService_service
                           art::Framework_Services_Optional_RandomNumberGenerator_service
			   art::Persistency_Common
			   art::Persistency_Provenance
			   art::Utilities
			   canvas::canvas
			   messagefacility::MF_MessageLogger 
 	                   fhiclcpp::fhiclcpp
			   cetlib::cetlib cetlib_except
			   CLHEP
                           ROOT_GEOM
                           ROOT_XMLIO
                           ROOT_GDML
			   ROOT_BASIC_LIB_LIST
        )

add_subdirectory(Utility)
add_subdirectory(Tool)
add_subdirectory(Service)
add_subdirectory(Module)
add_subdirectory(fcl)

install_source()

# for classes_def.xml
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make(
          MODULE_LIBRARIES dunesim_EventGenerator_ProtoDUNEbeamDataProducts
          LIB_LIBRARIES lardataobj::Simulation
          larcoreobj::headers
          nusimdata::SimulationBase
          art::Framework_Services_Registry
          art_root_io::tfile_support
          ROOT::Core
	  art_root_io::TFileService_service
          art::Persistency_Common
          art::Utilities 
          art::Persistency_Common
          art::Persistency_Provenance 
	  canvas::canvas
          messagefacility::MF_MessageLogger
          cetlib::cetlib cetlib_except::cetlib_except
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
install_fhicl()   
install_scripts()


#cet_make( LIBRARY_NAME     
#	  LIBRARIES lardataobj::Simulation
#	  LIBRARIES nusimdata::SimulationBase
#                    canvas::canvas
#                    messagefacility::MF_MessageLogger
#                    cetlib::cetlib cetlib_except::cetlib_except
#          NO_DICTIONARY )

#art_dictionary( DICTIONARY_LIBRARIES lardataobj::Simulation )

#install_headers()
#install_fhicl()
#install_source()

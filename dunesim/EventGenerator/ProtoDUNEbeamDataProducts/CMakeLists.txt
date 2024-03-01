# for classes_def.xml
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

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
        )

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

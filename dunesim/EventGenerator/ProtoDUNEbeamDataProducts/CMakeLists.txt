# for classes_def.xml
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(
          MODULE_LIBRARIES dune_EventGenerator_ProtoDUNEbeamDataProducts
          LIB_LIBRARIES lardataobj_Simulation
          ART_FRAMEWORK_CORE
          ART_FRAMEWORK_PRINCIPAL
          ART_FRAMEWORK_SERVICES_REGISTRY
          ART_ROOT_IO_TFILE_SUPPORT
          ROOT::Core
          ART_ROOT_IO_TFILESERVICE_SERVICE
          art::Persistency_Common
          art::Utilities 
          art::Persistency_Common
          art::Persistency_Provenance 
	  canvas::canvas
          messagefacility::MF_MessageLogger
          cetlib::cetlib cetlib_except
          ROOT_BASIC_LIB_LIST
        )

install_headers()
install_source()
install_fhicl()   
install_scripts()


#cet_make( LIBRARY_NAME     
#	  LIBRARIES lardataobj_Simulation
#	  LIBRARIES nusimdata::SimulationBase
#                    canvas::canvas
#                    messagefacility::MF_MessageLogger
#                    
#                    cetlib::cetlib cetlib_except
#                    ROOT_BASIC_LIB_LIST
#          NO_DICTIONARY )

#art_dictionary( DICTIONARY_LIBRARIES lardataobj_Simulation )

#install_headers()
#install_fhicl()
#install_source()

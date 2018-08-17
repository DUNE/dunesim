# for classes_def.xml
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(
          MODULE_LIBRARIES dune_EventGenerator_ProtoDUNEbeamDataProducts
          LIB_LIBRARIES lardataobj_Simulation
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
          art_Persistency_Common
          art_Utilities 
          art_Persistency_Common
          art_Persistency_Provenance 
	  canvas
          ${MF_MESSAGELOGGER}

          
          cetlib cetlib_except
          ${ROOT_BASIC_LIB_LIST}
          ProtoDUNEAnaUtils
        )

install_headers()
install_source()
install_fhicl()   
install_scripts()


#cet_make( LIBRARY_NAME     
#	  LIBRARIES lardataobj_Simulation
#	  LIBRARIES nusimdata_SimulationBase
#                    canvas
#                    MF_MessageLogger
#                    
#                    cetlib cetlib_except
#                    ${ROOT_BASIC_LIB_LIST}
#          NO_DICTIONARY )

#art_dictionary( DICTIONARY_LIBRARIES lardataobj_Simulation )

#install_headers()
#install_fhicl()
#install_source()

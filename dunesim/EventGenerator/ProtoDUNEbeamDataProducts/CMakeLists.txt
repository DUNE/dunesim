# for classes_def.xml
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(
          LIB_LIBRARIES lardataobj_Simulation
          art_Persistency_Common canvas
          art_Persistency_Provenance canvas
          ${MF_MESSAGELOGGER}
          ${MF_UTILITIES}
          cetlib cetlib_except
          ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_source()
install_fhicl()   

#cet_make( LIBRARY_NAME     
#	  LIBRARIES lardataobj_Simulation
#	  LIBRARIES nusimdata_SimulationBase
#                    canvas
#                    MF_MessageLogger
#                    MF_Utilities
#                    cetlib cetlib_except
#                    ${ROOT_BASIC_LIB_LIST}
#          NO_DICTIONARY )

#art_dictionary( DICTIONARY_LIBRARIES lardataobj_Simulation )

#install_headers()
#install_fhicl()
#install_source()

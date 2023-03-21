include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


art_make_library( LIBRARY_NAME ProtoDUNETriggeredBeamUtils
		  SOURCE ProtoDUNETriggeredBeamUtils.cxx
		  LIBRARIES
                       fhiclcpp::fhiclcpp
                       ROOT_BASIC_LIB_LIST 
	        )

art_make_exec( NAME runEventCounter
	       SOURCE runEventCounter.c++
	       LIBRARIES ProtoDUNETriggeredBeamUtils 
	    )

install_scripts()
install_fhicl()

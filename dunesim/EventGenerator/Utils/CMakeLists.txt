include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


art_make_library( LIBRARY_NAME BeamMiscUtils
                  SOURCE BeamMiscUtils.cxx
		  LIBRARIES
                       fhiclcpp::fhiclcpp
                       ROOT_BASIC_LIB_LIST
)


art_make_library( LIBRARY_NAME ProtoDUNETriggeredBeamUtils
		  SOURCE ProtoDUNETriggeredBeamUtils.cxx
		  LIBRARIES
                       fhiclcpp::fhiclcpp
                       ROOT_BASIC_LIB_LIST 
	        )

art_make_exec( NAME runEventCounter
	       SOURCE runEventCounter.c++
	       LIBRARIES
                ProtoDUNETriggeredBeamUtils 
                BeamMiscUtils
	    )

art_make_exec( NAME runYieldTest
	       SOURCE runYieldTest.c++
	       LIBRARIES
                ProtoDUNETriggeredBeamUtils 
                BeamMiscUtils
	    )

install_scripts()
install_fhicl()
INSTALL(FILES h4_vle_studies.cfg
        DESTINATION ${product}/${version}/cfg_files)

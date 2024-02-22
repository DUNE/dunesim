include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible


art_make_library( LIBRARY_NAME ProtoDUNETriggeredBeamUtils
		  SOURCE ProtoDUNETriggeredBeamUtils.cxx
		  LIBRARIES
                       fhiclcpp::fhiclcpp
		       ROOT::Core
		       ROOT::Physics
	        )

cet_make_exec( NAME runEventCounter
	       SOURCE runEventCounter.c++
	       LIBRARIES ProtoDUNETriggeredBeamUtils 
	    ) ### MIGRATE-CHANGED (migrate-3.22.02) - art_make_exec() -> cet_make_exec()

install_scripts()
install_fhicl()

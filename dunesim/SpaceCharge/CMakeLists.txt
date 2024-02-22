art_make( NO_PLUGINS
          
          LIB_LIBRARIES
          				lardataobj_RecoBase
                        lardata_Utilities
                         lardataalg_DetectorInfo
              larevt_SpaceCharge
			larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        art::Framework_Services_Registry 
                        cetlib::cetlib cetlib_except
			    art_root_io::tfile_support
			    ROOT::Core
			    art::Persistency_Common canvas
			    art::Persistency_Provenance canvas
			    art::Utilities canvas
			fhiclcpp::fhiclcpp
                            ROOT::Hist
			
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types



install_headers()
install_fhicl()
install_source()
install_scripts()

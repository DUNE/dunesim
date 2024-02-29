art_make( NO_PLUGINS
          
          LIB_LIBRARIES
          				lardataobj::RecoBase
                        lardata::Utilities
                         lardataalg::DetectorInfo
              larevt::SpaceCharge
			larcorealg::Geometry
                        messagefacility::MF_MessageLogger
                        art::Framework_Services_Registry 
                        cetlib::cetlib cetlib_except::cetlib_except
			    art_root_io::tfile_support
			    ROOT::Core
			    art::Persistency_Common
			    art::Persistency_Provenance canvas::canvas
			    art::Utilities
			fhiclcpp::fhiclcpp
                            ROOT::Hist
			
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types



install_headers()
install_fhicl()
install_source()
install_scripts()

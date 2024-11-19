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
			
        )



install_headers()
install_fhicl()
install_source()
install_scripts()

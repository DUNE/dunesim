art_make( NO_PLUGINS
          
          LIB_LIBRARIES
          				lardataobj_RecoBase
                        lardata_Utilities
                         lardataalg_DetectorInfo
              larevt_SpaceCharge
			larcorealg_Geometry
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except
            ART_FRAMEWORK_CORE
			    ART_FRAMEWORK_PRINCIPAL
			    ART_FRAMEWORK_SERVICES_REGISTRY
			    ART_ROOT_IO_TFILE_SUPPORT
			    ROOT::Core
			    ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
			    ART_ROOT_IO_TFILESERVICE_SERVICE
			    art::Persistency_Common canvas
			    art::Persistency_Provenance canvas
			    art::Utilities canvas
			fhiclcpp::fhiclcpp
                        ROOT_BASIC_LIB_LIST
			
        )



install_headers()
install_fhicl()
install_source()
install_scripts()

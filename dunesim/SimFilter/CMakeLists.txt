simple_plugin(GenExtCounterFilter "module"
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dune_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(LArG4ExtCounterFilter "module"
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dune_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()

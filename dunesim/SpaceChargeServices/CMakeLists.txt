simple_plugin(SpaceChargeServiceProtoDUNE "service"
              dune_SpaceCharge
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Persistency_Common canvas_Persistency_Common
	      art_Persistency_Provenance canvas_Persistency_Provenance
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art_Utilities canvas_Utilities
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
  	      ${CETLIB}
   	      ${ROOT_BASIC_LIB_LIST}
              
)

simple_plugin(SpaceChargeServiceDUNE35t "service"
              dune_SpaceCharge
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Persistency_Common canvas_Persistency_Common
	      art_Persistency_Provenance canvas_Persistency_Provenance
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      art_Utilities canvas_Utilities
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
  	      ${CETLIB}
   	      ${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source()

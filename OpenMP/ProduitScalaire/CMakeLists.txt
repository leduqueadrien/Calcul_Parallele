
#################################### SETUP ####################################

set( PS_EXEC_NAME produit_scalaire)

set( PS_SOURCES
	produitScalaire.c
	prodScalaire_dynamic_seq.c
)

set( PS_HEADER
	prodScalaire_dynamic_seq.h
)

############################ IMPLEMENTATION/BUILD #############################


if (PRODUIT_SCALAIRE)


	add_executable(${PS_EXEC_NAME} ${PS_SOURCES} ${PS_HEADER})
	
	if(OpenMP_C_FOUND)
		target_link_libraries(${PS_EXEC_NAME} PUBLIC OpenMP::OpenMP_C)
	endif()

endif()

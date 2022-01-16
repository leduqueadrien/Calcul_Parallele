
#################################### SETUP ####################################

set( PM_EXEC_NAME produit_matriciel)

set( PM_SOURCES
	produitMatriciel.c
	prodMatriciel_dynamic_seq.c
)

set( PM_HEADER
	prodMatriciel_dynamic_seq.h
)

############################ IMPLEMENTATION/BUILD #############################


if (PRODUIT_MATRICIEL)


	add_executable(${PM_EXEC_NAME} ${PM_SOURCES} ${PM_HEADER})
	
	if(OpenMP_C_FOUND)
		target_link_libraries(${PM_EXEC_NAME} PUBLIC OpenMP::OpenMP_C)
	endif()

endif()

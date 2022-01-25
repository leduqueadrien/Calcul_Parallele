
#################################### SETUP ####################################

set( EQ_CHAL_EXEC_NAME produit_scalaire)

set( EQ_CHAL_SOURCES
	ResolutionEDO.c
)

set( EQ_CHAL_HEADER
	ResolutionEDO.h
)

############################ IMPLEMENTATION/BUILD #############################


if (EQUATION_CHALEUR)

	add_executable(${EQ_CHAL_EXEC_NAME} ${EQ_CHAL_SOURCES} ${EQ_CHAL_HEADER})
	
	if(OpenMP_C_FOUND)
		target_link_libraries(${EQ_CHAL_EXEC_NAME} PRIVATE OpenMP::OpenMP_C)
	endif()

	target_link_libraries(${EQ_CHAL_EXEC_NAME} PRIVATE m)


endif()

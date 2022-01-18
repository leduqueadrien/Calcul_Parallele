
#################################### SETUP ####################################

set( NEWTONSDL_EXEC_NAME newton)

set( NEWTONSDL_SOURCES
	newton.c
	complexe.c
	fonctions.c
	SDLgestion.c
)

set( NEWTONSDL_HEADER
	complexe.h
	fonctions.h
	SDLgestion.h
)

############################ IMPLEMENTATION/BUILD #############################


if (NEWTON_GRAPHIQUE)

	add_executable(${NEWTONSDL_EXEC_NAME} ${NEWTONSDL_SOURCES} ${NEWTONSDL_HEADER})
	
	if(OpenMP_C_FOUND)
		target_link_libraries(${NEWTONSDL_EXEC_NAME} PUBLIC OpenMP::OpenMP_C)
	endif()

	find_package(SDL2 REQUIRED)
	include_directories(${SDL2_INCLUDE_DIRS})
	target_link_libraries(${NEWTONSDL_EXEC_NAME} PUBLIC ${SDL2_LIBRARIES})


endif()

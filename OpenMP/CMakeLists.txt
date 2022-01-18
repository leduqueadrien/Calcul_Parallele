
################################### PROJECT ###################################


project ( TP_openMP )

cmake_minimum_required(VERSION 3.1)

################################### COMPILATION OPTIONS ###################################

if (MSVC)
	# warning level 4 and all warnings as errors
	# EHsc permet d'eviter un warning avec la STL
	add_compile_options(/EHsc /W4 /std:c++14)
	# ajout de open MP
	# add_compile_options(/openmp)
endif()

if ( CMAKE_COMPILER_IS_GNUCXX )
	add_compile_options(-Wall -Wextra -pedantic)
	# ajout de open MP
	# add_compile_options(-fopenmp)
endif()


#################################### SETUP ####################################

option(PRODUIT_SCALAIRE "compilation du produit scalaire" OFF)
option(PRODUIT_MATRICIEL "compilation du produit matriciel" OFF)
option(NEWTON_GRAPHIQUE "compilation du rendu graphique des fractales de Newton" ON)
option(NEWTON_NON_GRAPHIQUE "compilation de la generation d'image png des fractales de Newton" OFF)

find_package(OpenMP)

############################ IMPLEMENTATION/BUILD #############################

add_subdirectory(ProduitMatriciel)
add_subdirectory(ProduitScalaire)
add_subdirectory(NewtonSDL)

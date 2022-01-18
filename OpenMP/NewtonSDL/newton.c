
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include "SDLgestion.h"
#include "complexe.h"
#include "fonctions.h"
#include "omp.h"

void calculTabCouleur(unsigned int* tabCouleur, unsigned int width, unsigned int height) {
	unsigned int      xe, ye;        /* coordonnees de l'ecran */
	unsigned int      c;             /* couleur */
	complexe_t        x;
	double            start_time;


	start_time = omp_get_wtime();
	#pragma omp parallel for private(x, xe, ye, c) shared(width, height, tabCouleur) num_threads(2)
	for (xe=0; xe<width; xe++) {
		for (ye=0; ye<height; ye++) {
			x.re = (double)XMIN + ((double)xe*(XMAX-XMIN))/((double)(width-1));
			x.im = (double)YMAX - ((double)ye*(YMAX-YMIN))/((double)(height-1));
			
			if (x.re==0 && x.im==0)
				c=0;
			else
				c = couleur3(x);

			tabCouleur[ye*width + xe] = c;
		}
	}
	printf("temps : %f\n", omp_get_wtime() - start_time);
}

void afficher(SDL_Renderer* renderer, unsigned int* tabCouleur, unsigned int width, unsigned int height) {
	unsigned int i;
	unsigned int j;

	for(i=0; i<width; ++i) {
		for(j=0; j<height; j++) {
			setColor(renderer, tabCouleur[j*width + i]);
			SDL_RenderDrawPoint(renderer, i, j);
		}
	}

}

int main() {
	/* -------------------------------------------------------------------------- */
	/* Les variables du programme												  */
	/* -------------------------------------------------------------------------- */
	int 			failure 	= 0;
	int 			running 	= 1;
	unsigned int* 	tabCouleur 	= NULL;

	/* -------------------------------------------------------------------------- */
	/* Les variables de la fenetre SDL											  */
	/* -------------------------------------------------------------------------- */
	unsigned int 	width 		= 1000;
	unsigned int 	height 		= 1000;
	SDL_Window* 	window;
	SDL_Renderer* 	renderer;
	SDL_Event 		event;

	/* -------------------------------------------------------------------------- */
	/* Initialisation du tableau de couleur										  */
	/* -------------------------------------------------------------------------- */
	tabCouleur = (unsigned int*) malloc(width*height * sizeof(unsigned int));
	if (tabCouleur == NULL)
		return EXIT_FAILURE;

	/* -------------------------------------------------------------------------- */
	/* Initialisation de la fenetre												  */
	/* -------------------------------------------------------------------------- */
	failure = initSDL(&window, &renderer, width, height);
	if (failure == 1) {
		free(tabCouleur);
		return EXIT_FAILURE;
	}

	/* couleur de fond */
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	SDL_RenderClear(renderer);


	/* Calcul de la couleur que doit prendre chaque pixel */
	calculTabCouleur(tabCouleur, width, height);
	
	while (running) {
		/* afficher à l'ecran */
		SDL_RenderPresent(renderer);
		
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_QUIT)
				running = 0;
		}

		afficher(renderer, tabCouleur, width, height);

	}

	SDL_Delay(5000);

	quitSDL(window, renderer);
	free(tabCouleur);

}

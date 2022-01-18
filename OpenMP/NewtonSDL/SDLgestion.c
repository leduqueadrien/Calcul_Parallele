
#include <SDL2/SDL.h>
#include "SDLgestion.h"

void quitSDL(SDL_Window* window, SDL_Renderer* renderer) {
	if (renderer != NULL)
		SDL_DestroyRenderer(renderer);
	
	if (window != NULL)
		SDL_DestroyWindow(window);

	SDL_Quit();
}

int initSDL(SDL_Window** window, SDL_Renderer** renderer, unsigned int width, unsigned int height) {
	int failure = 0;

	*window = NULL;
	*renderer = NULL;


	if (SDL_Init(SDL_INIT_VIDEO) < 0)
	{
		fprintf(stderr, "Erreur d'initialisation de la SDL : %s\n", SDL_GetError());
		failure = 1;
	} else {
		*window = SDL_CreateWindow("SDL2 Programme 0.1", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
					width, height, 
					0);
			
		if (*window == 0)
		{
			fprintf(stderr, "Erreur d'initialisation de la SDL : %s\n", SDL_GetError());
			quitSDL(*window, *renderer);
		} else {

			*renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_ACCELERATED ); /*  SDL_RENDERER_SOFTWARE */
			if (renderer == 0) {
				fprintf(stderr, "Erreur d'initialisation de la SDL : %s\n", SDL_GetError());
				quitSDL(*window, *renderer);

			}

		}

	}
	return failure;
}

void setColor(SDL_Renderer* renderer, unsigned int c) {
	if (c == 1)
		SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	else if (c ==2)
		SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
	else
		SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
}

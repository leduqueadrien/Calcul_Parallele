
#ifndef SDLGESTION_HPP
#define SDLGESTION_HPP

#include <SDL2/SDL.h>

void quitSDL(SDL_Window* window, SDL_Renderer* renderer);
int initSDL(SDL_Window** window, SDL_Renderer** renderer, unsigned int width, unsigned int height);
void setColor(SDL_Renderer* renderer, unsigned int c);	

#endif
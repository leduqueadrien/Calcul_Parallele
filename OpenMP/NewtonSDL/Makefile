# Makefile pour les applications X

CC = gcc 
CFLAGS = -I/usr/X11/include -Wall -Wextra -pedantic
SDL2LIB = -lSDL2
LDOMP = -fopenmp

SRC1 = newton.c SDLgestion.c fonctions.c complexe.c
OBJ1 = $(SRC1:.c=.o)

newton : $(OBJ1)
	$(CC) -o $@ $(OBJ1) $(SDL2LIB) $(LDOMP)

%.o : %.c 
	$(CC) $(CFLAGS) $(LDOMP) -c $<

clean:	
	\rm  newton *.o


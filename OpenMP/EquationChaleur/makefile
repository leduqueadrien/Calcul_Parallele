CFLAGS = -g -Wall -Wextra -lm
OBJ = ResolutionEDO.o
EXEC = run

all: $(EXEC)

run : $(OBJ)
	gcc -o $(EXEC) $(OBJ) $(CFLAGS)

ResolutionEDO.o : ResolutionEDO.c ResolutionEDO.h
	gcc -o ResolutionEDO.o -c ResolutionEDO.c $(CFLAGS)

clean : 
	rm -f $(OBJ)
	rm -f $(EXEC)

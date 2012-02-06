# Makefile to build driver for Van Leer flux and Jacobians routine 
# --- macros
CC=gcc
CFLAGS=  -O3 -Wall
OBJECTS= flux.o
LIBS = -lm

# --- targets
all: vanleer
vanleer: $(OBJECTS) 
	$(CC)  -o vanleer  $(OBJECTS) $(LIBS)
flux.o : flux.c
	$(CC) $(CFLAGS) -c flux.c
# --- remove object and executable files
clean:
	rm -f vanleer $(OBJECTS)

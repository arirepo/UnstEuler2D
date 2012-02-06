# Makefile to build project 
# --- macros
CC=gcc
CFLAGS=  -O3 -Wall
OBJECTS= unst2d.o flux.o util2d.o
LIBS = -lm

# --- targets
all: unst2d
unst2d: $(OBJECTS) 
	$(CC)  -o unst2d  $(OBJECTS) $(LIBS)
unst2d.o : unst2d.c
	$(CC) $(CFLAGS) -c unst2d.c
flux.o : flux.c
	$(CC) $(CFLAGS) -c flux.c
util2d.o : util2d.c
	$(CC) $(CFLAGS) -c util2d.c

# --- remove object and executable files
clean:
	rm -f unst2d $(OBJECTS)

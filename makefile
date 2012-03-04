# Makefile to build unstructured 2d solver project 
# --- macros
CC=g++
CFLAGS=  -O3 -Wall
OBJECTS= unst2d.o flux.o util2d.o grid_reader.o residuals.o explicit.o maps.o implicit.o
LIBS = -lm

# --- targets
all: unst2d
unst2d: $(OBJECTS) 
	$(CC) $(CFLAGS) -o unst2d  $(OBJECTS) $(LIBS)
unst2d.o : unst2d.c
	$(CC) $(CFLAGS) -c unst2d.c
flux.o : flux.c
	$(CC) $(CFLAGS) -c flux.c
util2d.o : util2d.c
	$(CC) $(CFLAGS) -c util2d.c
grid_reader.o : grid_reader.c
	$(CC) $(CFLAGS) -c grid_reader.c
residuals.o : residuals.c
	$(CC) $(CFLAGS) -c residuals.c
explicit.o : explicit.c
	$(CC) $(CFLAGS) -c explicit.c
maps.o : maps.cpp
	$(CC) $(CFLAGS) -c maps.cpp
implicit.o : implicit.cpp
	$(CC) $(CFLAGS) -c implicit.cpp

# --- remove object and executable files
clean:
	rm -f unst2d $(OBJECTS)

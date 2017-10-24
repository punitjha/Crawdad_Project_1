CC=g++-5
CFLAGS= -std=c++11 -g -O3 -larmadillo -llapack -lblas 
OBJ =molecule.o mole_main.o
DEPS = mass.h molecule.h
all: mole.exe

mole.exe: $(OBJ)
	$(CC) -o mole.exe $^ $(CFLAGS)

%.o: %.c $(DEP)
	$(CC) -c $< $(CFLAGS)

clean:
	rm *.o *.exe

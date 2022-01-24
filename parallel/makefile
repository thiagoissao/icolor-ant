CC = gcc
#CFLAGS = -pg -g -ggdb -Wall -O3 -DNRAND 
BIN = colorant
CFLAGS = -Wall -O3 -DNRAND  
MAIN = main.c
OBJS = ./src/icolor.o ./src/icolorant.o ./src/util.o ./src/tabucol.o ./src/ant_fixed_k.o


all: $(OBJS) 
	$(CC) -o $(BIN) $(CFLAGS) $(MAIN) $(OBJS) -lm

clean:
	rm $(OBJS) $(BIN)



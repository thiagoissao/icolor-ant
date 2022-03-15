CC = gcc
#CFLAGS = -pg -g -ggdb -Wall -O3 -DNRAND 
BIN = colorant
CFLAGS = -Wall -O3 -DNRAND  
MAIN = main.c
OBJS = ./src/icolor.o ./src/icolorant/icolorant.o ./src/icolorant/icolorant.utils.o ./src/util.o ./src/tabucol.o ./src/ant_fixed_k.o ./src/helpers.o


all: $(OBJS) 
	$(CC) -o $(BIN) $(CFLAGS) $(MAIN) $(OBJS) -lm

exec-small: 
	 ./colorant ../instances/dsjc125.5.col -p 3 -k 17 -c 1000 -w 2

exec-large: 
	 ./colorant ../instances/dsjc250.5.col -p 3 -k 28 -c 1000 -w 1

clean:
	rm $(OBJS) $(BIN)



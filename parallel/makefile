CC = gcc
#CFLAGS = -pg -g -ggdb -Wall -O3 -DNRAND 
BIN = colorant
CFLAGS = -Wall -O3 -DNRAND  
OBJS = main.o icolor.o icolorant.o util.o tabucol.o ant_fixed_k.o


all: $(OBJS) 
	$(CC) -o $(BIN) $(CFLAGS) $(OBJS) -lm
	rm *.o

clean:
	rm $(OBJS) $(BIN) *~ 



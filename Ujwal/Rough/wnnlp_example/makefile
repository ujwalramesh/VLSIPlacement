CC=gcc

all: hello

hello: examples.o wnnlp.o
	$(CC) wnnlp.o examples.o -o hello

examples.o: examples.c
	$(CC) examples.c -o examples.o

clean: 
	rm -rf *o hello 

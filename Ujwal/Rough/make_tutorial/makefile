CC=g++

all: hello

hello: main.o factorial.o hello.o 
	$(CC) main.o factorial.o hello.o -o hello 

main.o: main.cpp
	$(CC) -c main.cpp 

factorial.o: factorial.cpp 
	$(CC) -c factorial.cpp 

hello.o: hello.cpp
	$(CC) -c hello.cpp 
clean:
	rm -rf *o hello

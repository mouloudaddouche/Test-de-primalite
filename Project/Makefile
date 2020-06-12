run : main
	./main
main: main.o
	gcc -o main main.o -lgmp

main.o: main.c main.h
	gcc -Wall -o main.o -c main.c -lgmp
clean : 
	rm main.o

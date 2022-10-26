build:
	gcc -o main main.c sequential/naive/TSP.h sequential/genetic/TSP.h -lpthread
clean:
	rm main

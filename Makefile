build:
	gcc -o main main.c sequential/naive/TSP.h sequential/genetic/TSP.h -lpthread

run_genetic_seq_1:
	./main sequential_genetic input/input1.in

run_genetic_seq_2:
	./main sequential_genetic input/input2.in

run_genetic_seq_3:
	./main sequential_genetic input/input3.in

run_genetic_seq_4:
	./main sequential_genetic input/input4.in

run_genetic_seq_5:
	./main sequential_genetic input/input5.in

run_naive_1:
	./main sequential_naive input/input1.in

run_naive_2:
	./main sequential_naive input/input2.in

run_naive_3:
	./main sequential_naive input/input3.in

run_naive_4:
	./main sequential_naive input/input4.in

run_naive_5:
	./main sequential_naive input/input5.in

clean:
	rm main

naive = sequential/naive/TSP.h
genetic = sequential/genetic/TSP.h
openmp = parallel/openmp/TSP.h
pthreads = parallel/pthreads/TSP.h
mpi = parallel/mpi/TSP_mpi.h
mpi_omp = parallel/mpi_omp/TSP_mpi_omp.h


build:
	mpicc -o main main.c $(naive) $(genetic) $(openmp) $(pthreads) $(mpi) $(mpi_omp) -lpthread -fopenmp -lm

build_mpi:
	mpicc main.c $(mpi) -o main -fopenmp

run_mpi:
	mpirun -np 4 main parallel_mpi input/input4.in

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

run_openmp_1:
	./main parallel_openmp input/input1.in

run_openmp_2:
	./main parallel_openmp input/input2.in

run_openmp_3:
	./main parallel_openmp input/input3.in

run_openmp_4:
	./main parallel_openmp input/input4.in

run_openmp_5:
	./main parallel_openmp input/input5.in

run_pthreads_1:
	./main parallel_pthreads input/input1.in

run_pthreads_2:
	./main parallel_pthreads input/input2.in

run_pthreads_3:
	./main parallel_pthreads input/input3.in

run_pthreads_4:
	./main parallel_pthreads input/input4.in

run_pthreads_5:
	./main parallel_pthreads input/input5.in

run_mpi_omp_4:
	mpirun -np 4 main parallel_mpi_omp input/input4.in
clean:
	rm main

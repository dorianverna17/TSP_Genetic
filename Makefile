naive = sequential/naive/TSP.h
genetic = sequential/genetic/TSP.h
openmp = parallel/openmp/TSP.h
pthreads = parallel/pthreads/TSP.h
mpi = parallel/mpi/TSP_mpi.h
mpi_omp = parallel/mpi_omp/TSP_mpi_omp.h
mpi_pthreads = parallel/mpi-pthreads/TSP_mpi_pthreads.h

build:
	mpicc -o main main.c $(naive) $(genetic) $(openmp) $(pthreads) $(mpi) $(mpi_omp) $(mpi_pthreads) -lpthread -fopenmp -lm

run_mpi_1:
	mpirun -np 4 main parallel_mpi input/input1.in

run_mpi_2:
	mpirun -np 4 main parallel_mpi input/input2.in

run_mpi_3:
	mpirun -np 4 main parallel_mpi input/input3.in

run_mpi_4:
	mpirun -np 4 main parallel_mpi input/input4.in

run_mpi_5:
	mpirun -np 4 main parallel_mpi input/input5.in

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

run_mpi_omp_1:
	mpirun -np 2 main parallel_mpi_omp input/input1.in

run_mpi_omp_2:
	mpirun -np 2 main parallel_mpi_omp input/input2.in

run_mpi_omp_3:
	mpirun -np 2 main parallel_mpi_omp input/input3.in

run_mpi_omp_4:
	mpirun -np 2 main parallel_mpi_omp input/input4.in

run_mpi_pthreads_1:
	mpirun -np 2 main parallel_mpi_pthreads input/input1.in

run_mpi_pthreads_2:
	mpirun -np 2 main parallel_mpi_pthreads input/input2.in

run_mpi_pthreads_3:
	mpirun -np 2 main parallel_mpi_pthreads input/input3.in

run_mpi_pthreads_4:
	mpirun -np 2 main parallel_mpi_pthreads input/input4.in

run_mpi_pthreads_5:
	mpirun -np 2 main parallel_mpi_pthreads input/input5.in

clean:
	rm main

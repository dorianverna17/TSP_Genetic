#!/bin/bash

for i in 1 2 3 4 5 6 7 8 9
do
	./main sequential_genetic input/input$i.in
	./main parallel_openmp input/input$i.in
	./main parallel_pthreads input/input$i.in
	mpirun -np 4 ./main parallel_mpi input/input$i.in
	mpirun -np 2 ./main parallel_mpi_omp input/input$i.in
	mpirun -np 2 ./main parallel_mpi_pthreads input/input$i.in
	echo ""
done
#!/bin/bash

# for i in 1 2 3 4 5
# do
# 	./main sequential_genetic input/input$1.in
# done

# for i in 1 2 3 4 5
# do
# 	./main parallel_openmp input/input$1.in
# done

./main sequential_genetic input/input$1.in
./main parallel_openmp input/input$1.in
./main parallel_pthreads input/input$1.in
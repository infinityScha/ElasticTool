#!/bin/bash

# Define the range of arguments
nu=(0.85 0.88 0.91 0.94 0.97 0.82 0.79 0.76)
nu_prev=(0.85 0.85 0.88 0.91 0.94 0.85 0.82 0.79)
dA=(2.1 1.9 1.7 1.5 1.3 1.1)
dA_prev=(2.2 2.0 1.8 1.6 1.4 1.2)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..8}
do
    for j in {0..5}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=3 ForPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close_f.log
    done
done
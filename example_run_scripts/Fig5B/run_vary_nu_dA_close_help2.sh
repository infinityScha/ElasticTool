#!/bin/bash

# Define the range of arguments
nu=(0.94 0.97)
nu_prev=(0.85 0.85)
dA=(2.1 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0)
dA_prev=(2.2 2.1 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..2}
do
    for j in {0..11}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=3 ForPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close_h2.log
    done
done
#!/bin/bash

# Define the range of arguments
nu=(0.76)
nu_prev=(0.79)
dA=(1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2)
dA_prev=(1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..0}
do
    for j in {0..11}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=4 ForPhaseSpace.jl open ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_open_h.log
    done
done
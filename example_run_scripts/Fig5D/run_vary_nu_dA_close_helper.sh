#!/bin/bash

# Define the range of arguments
nu=(0.79 0.76)
nu_prev=(0.82 0.79)
dA=(1.35 1.3 1.25 1.2 1.15)
dA_prev=(1.4 1.35 1.3 1.25 1.2)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..8}
do
    echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=3.0"
    nohup julia --threads=6 ForPhaseSpace.jl close ${nu_prev[$i]} ${nu[$i]} 1.4 1.4 > curr_close_h.log
done

for i in {0..8}
do
    for j in {0..5}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=6 ForPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close.log
    done
done
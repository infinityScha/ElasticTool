#!/bin/bash

# Define the range of arguments
nu=(0.91 0.94 0.97 1.0 0.88 0.85 0.82 0.79 0.76)
nu_prev=(0.91 0.91 0.94 0.97 0.91 0.88 0.85 0.82 0.79)
dA=(1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2)
dA_prev=(1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..8}
do
    echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=1.2"
    nohup julia --threads=3 ForPhaseSpace.jl open ${nu_prev[$i]} ${nu[$i]} 1.0 1.0 > curr_open.log
done

for i in {0..8}
do
    for j in {0..11}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=3 ForPhaseSpace.jl open ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_open.log
    done
done
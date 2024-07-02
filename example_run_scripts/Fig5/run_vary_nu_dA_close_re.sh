#!/bin/bash

# Define the range of arguments
nu=(0.91 0.94 0.97 0.88 0.85 0.82 0.79 0.76 0.73)
nu_prev=(0.91 0.91 0.94 0.91 0.88 0.85 0.82 0.79 0.76)
dA=(2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0 1.9 1.8 1.7 1.6 1.5)
dA_prev=(3.0 2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0 1.9 1.8 1.7 1.6)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..8}
do
    echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=3.0"
    nohup julia --threads=3 RerunPhaseSpace.jl close ${nu_prev[$i]} ${nu[$i]} 3.0 3.0 > curr_close.log
done

for i in {0..8}
do
    for j in {0..14}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=3 RerunPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close.log
    done
done
#!/bin/bash

# Define the range of arguments
nu=(0.91 0.94 0.97 1.0 0.88 0.85 0.82 0.79 0.76 0.73 0.7)
nu_prev=(0.91 0.91 0.94 0.97 0.91 0.88 0.85 0.82 0.79 0.76 0.73)
dA=(2.8 2.6 2.4 2.2 2.0 1.8 1.6)
dA_prev=(3.0 2.8 2.6 2.4 2.2 2.0 1.8)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..10}
do
    echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=3.0"
    nohup julia --threads=4 ForPhaseSpace.jl close ${nu_prev[$i]} ${nu[$i]} 3.0 3.0 > curr_close.log
done

for i in {0..10}
do
    for j in {0..6}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=4 ForPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close.log
    done
done
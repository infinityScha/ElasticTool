#!/bin/bash

# Define the range of arguments
nu=(0.94 0.97)
nu_prev=(0.91 0.94)
dA=(2.2 2.4 2.6 2.8 3.0 1.8 1.6)
dA_prev=(2.0 2.2 2.4 2.6 2.8 2.0 1.8)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl

for i in {0..4}
do
    echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=2.0"
    nohup julia --threads=8 RerunPhaseSpace.jl open ${nu_prev[$i]} ${nu[$i]} 2.0 2.0 > curr_open.log
done

for i in {0..4}
do
    for j in {0..6}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=8 RerunPhaseSpace.jl open ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_open.log
    done
done
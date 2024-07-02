#!/bin/bash

# Define the range of arguments
nu=(0.76)
nu_prev=(0.79)
dA=(1.6 1.5 1.4 1.3)
dA_prev=(1.7 1.6 1.5 1.4)

# Loop through the arguments of nu in the outer loop and dA in the inner loop and run ForPhaseSpace.jl


for i in {0..1}
do
    nohup julia --threads=4 ForPhaseSpace.jl close ${nu_prev[$i]} ${nu[$i]} 1.7 1.7 > curr_close_h.log
    for j in {0..3}
    do
        echo "Running ForPhaseSpace.jl with nu=${nu[$i]} and dA=${dA[$j]}"
        nohup julia --threads=4 ForPhaseSpace.jl close ${nu[$i]} ${nu[$i]} ${dA_prev[$j]} ${dA[$j]} > curr_close_h.log
    done
done
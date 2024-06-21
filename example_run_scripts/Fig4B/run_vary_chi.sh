#!/bin/bash

# Define the range of arguments
chi1=(0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)
chi1_prev=(0.0 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8)

# Loop through the arguments of js1 and js1_prev and run ForHeatmap.jl
for i in {6..10}
do
    echo "Running MakePath.jl with chi1=${chi1[$i]} and chi1_prev=${chi1_prev[$i]}"
    nohup julia --threads 2 MakePath.jl ${chi1_prev[$i]} ${chi1[$i]} 0 0 &>curr.log
done
